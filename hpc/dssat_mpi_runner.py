# new v5.1 (Fixed Error Logging + Header Logic)
import os
import subprocess
import glob
import re
import csv
import shutil
import argparse

# --- GLOBAL CONFIGURATION ---
OUTPUT_FIELDNAMES = [
    'point_id', 'treatment', 'crop_code', 'latitude', 'longitude',
    'weather_station_id', 'soil_profile_id', 'dssat_file_id', 'dssat_description',
    'planting_date', 'emergence_date', 'harvest_date', 'year_planting', 'year_harvest',
    'top_weight_kg_ha', 'final_grain_kg_ha', 'removed_residue_kg_ha',
    'carbon_soil_organic_matter_start_kg_ha', 'carbon_soil_organic_matter_end_kg_ha',
    'carbon_soil_organic_matter_delta_kg_ha', 'final_irrigation_applications_count',
    'final_irrigation_amount_mm', 'inorganic_n_applied_count', 'inorganic_n_applied_kg_ha',
    'nitrate_leaching_kg_ha', 'cumulative_net_co2_emissions_kg_CO2_ha',
    'cumulative_n2o_emissions_kg_N_ha',
    # Seasonal-average stress indicators from PlantGro.OUT
    'water_stress_photosynthesis_avg',
    'water_stress_development_avg',
    'nitrogen_stress_avg'
]

# Files produced by DSSAT that are commonly large / transient.
# Used for both cleanup and (optional) archiving.
DSSAT_OUTPUT_PATTERNS_BASE = [
    "*.OUT", "*.LST", "*.lst", "*.OOV", "*.oov",
    "summary.csv", "soilorg.csv", "soilni.csv", "soilwat.csv",
    "et.csv", "etphot.csv", "evaluate.csv", "mulch.csv", "n2o.csv",
    "plantgro.csv", "plantn.csv", "plantc.csv", "soiltemp.csv",
    "somlitc.csv", "somlitn.csv", "weather.csv",
]

# --- Argparse ---
def get_args():
    parser = argparse.ArgumentParser(description="Run DSSAT Simulation with MPI")
    parser.add_argument('--base_dir', type=str, required=True, help='Base directory for simulation folders')
    parser.add_argument('--summary_dir', type=str, required=True, help='Directory for summary output')
    parser.add_argument('--exe_path', type=str, default=os.environ.get('DSSAT_EXE', 'dscsm048'), help='Path to DSSAT executable (or set env DSSAT_EXE)')
    parser.add_argument('--model_code', type=str, default='CRGRO048', help='DSSAT model code')
    parser.add_argument('--run_mode', type=str, default='sequence', choices=['sequence', 'experiment'], help='Run mode')
    parser.add_argument('--trt_start', type=int, default=1, help='Treatment start')
    parser.add_argument('--trt_end', type=int, default=1, help='Treatment end')
    parser.add_argument('--seq_start', type=int, default=1, help='Sequence start')
    parser.add_argument('--seq_end', type=int, default=1, help='Sequence end')
    parser.add_argument(
        '--cleanup_mode',
        type=str,
        default='always',
        choices=['always', 'success', 'never'],
        help=("Cleanup behavior for DSSAT outputs in each point folder: "
              "'always' (default) deletes after each run; "
              "'success' deletes only when DSSAT+parsing succeeded; "
              "'never' never deletes (debugging; can accumulate files)."),
    )
    parser.add_argument(
        '--archive_outputs',
        action='store_true',
        help=("Archive DSSAT output files for each point/treatment into a per-point "
              "'_archive_outputs/' subfolder before optional cleanup. "
              "Files are moved (not copied), so you keep one copy without cluttering the point folder."),
    )
    parser.add_argument(
        '--scratch_dir',
        type=str,
        default=os.environ.get('SCRATCH_DIR', os.environ.get('TMPDIR', '')),
        help=("Node-local scratch directory for streaming per-rank part files "
              "(e.g. $TMPDIR on a SLURM node). Writing parts here avoids a "
              "small-file metadata storm on shared Lustre/GPFS; parts are staged "
              "to the shared summary dir at the end. Blank = write parts directly "
              "to the shared summary dir (legacy behavior)."),
    )
    parser.add_argument(
        '--merge_mode',
        type=str,
        default='concat',
        choices=['concat', 'none'],
        help=("Final merge policy. 'concat' (default) streams all per-rank parts "
              "into one results_<project>.csv. 'none' leaves the per-rank parts in "
              "place (each is a complete, valid CSV) — use on very large grids to "
              "avoid a rank-0 serial read of the entire dataset."),
    )

    return parser.parse_args()

args = get_args()

# Parse --help before requiring the cluster MPI runtime. This lets users and CI
# validate the command line on ordinary workstations where mpi4py is omitted by
# design. A real run still fails immediately with an actionable message.
try:
    from mpi4py import MPI
except ImportError as exc:
    raise SystemExit(
        "mpi4py is required to run the DSSAT MPI worker. Load the cluster MPI "
        "module, then install mpi4py against that MPI implementation."
    ) from exc

# Map arguments (normalize paths early to avoid surprises with trailing slashes / relative paths)
dssat_simulation_output_base = os.path.abspath(os.path.normpath(args.base_dir))
dssat_simulation_folder = os.path.basename(dssat_simulation_output_base)
dssat_simulation_summary_folder = os.path.abspath(os.path.normpath(args.summary_dir))

executable_path = args.exe_path
dssat_model_code = args.model_code
run_mode = args.run_mode
treatment_start = args.trt_start
treatment_end = args.trt_end
sequence_start = args.seq_start
sequence_end = args.seq_end
cleanup_mode = args.cleanup_mode
archive_outputs = bool(args.archive_outputs)
scratch_dir = args.scratch_dir.strip() if args.scratch_dir else ''
merge_mode = args.merge_mode
control_file_name = 'DSSBatch.V48'

# Temp rank-part files will be written under: <summary_dir>/<project_name>/
# (prevents clutter and avoids collisions if multiple projects share the same summary_dir)
project_temp_results_dir = os.path.join(dssat_simulation_summary_folder, dssat_simulation_folder)

# --- Initialize MPI ---
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# --- Helper Functions ---
def read_csv_to_dict_list(file_path):
    data = []
    try:
        with open(file_path, 'r', newline='', encoding='utf-8', errors='ignore') as f_in:
            reader = csv.DictReader(f_in)
            for row in reader:
                data.append(row)
    except FileNotFoundError:
        pass
    except Exception as e:
        print(f"Warning: Error reading CSV file {file_path}: {e}")
    return data

def get_grouped_csv_data(file_path, group_key_column='RUN'):
    grouped_data = {}
    try:
        with open(file_path, 'r', newline='', encoding='utf-8', errors='ignore') as f_in:
            reader = csv.DictReader(f_in)
            if reader.fieldnames and group_key_column not in reader.fieldnames:
                return {}
            for row in reader:
                key_value = row.get(group_key_column)
                if key_value:
                    if key_value not in grouped_data:
                        grouped_data[key_value] = []
                    grouped_data[key_value].append(row)
    except FileNotFoundError:
        pass
    except Exception as e:
        print(f"Warning: Error reading or grouping CSV {file_path}: {e}")
    return grouped_data

def write_dict_list_to_csv(data_list, file_path, fieldnames):
    try:
        with open(file_path, 'w', newline='', encoding='utf-8') as f_out:
            writer = csv.DictWriter(f_out, fieldnames=fieldnames, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(data_list)
    except Exception as e:
        print(f"Error writing CSV to {file_path}: {e}")

def parse_plantgro_stress(folder_path):
    """
    Compute seasonal-average stress indicators for a DSSAT run from either:
      - plantgro.csv (preferred, if present)
      - PlantGro.OUT (fallback)

    Stress indices can vary throughout the season, so we return the mean across
    all daily rows present in the PlantGro output for each run.

    Returns:
        {'by_run': {runno: {'wspd': mean, 'wsgd': mean, 'nstd': mean}},
         'by_trt': {trtno: {'wspd': mean, 'wsgd': mean, 'nstd': mean}}}
    """

    def _safe_int(x):
        try:
            return int(str(x).strip())
        except Exception:
            return None

    def _safe_float(x):
        try:
            v = float(str(x).strip())
        except Exception:
            return None
        # common DSSAT missing sentinels are <= -90 (e.g., -99, -999)
        if v <= -90:
            return None
        return v

    def _mean(vals):
        return (sum(vals) / len(vals)) if vals else None

    # ---- 1) Try plantgro.csv (most DSSAT workflows export this) ----
    csv_path_candidates = [
        os.path.join(folder_path, 'plantgro.csv'),
        os.path.join(folder_path, 'PlantGro.csv'),
        os.path.join(folder_path, 'PLANTGRO.CSV'),
        os.path.join(folder_path, 'PLANTGRO.csv'),
    ]
    csv_path = next((p for p in csv_path_candidates if os.path.exists(p)), None)

    if csv_path:
        by_run_acc = {}  # run -> {'wspd':[], 'wsgd':[], 'nstd':[]}
        by_trt_acc = {}  # trt -> {'wspd':[], 'wsgd':[], 'nstd':[]}

        try:
            with open(csv_path, 'r', encoding='utf-8', errors='ignore', newline='') as f:
                reader = csv.DictReader(f)
                # Normalize fieldnames to handle odd capitalization
                field_map = { (fn or '').strip().upper(): fn for fn in (reader.fieldnames or []) }

                run_col = field_map.get('RUN') or field_map.get('RUNNO')
                trt_col = field_map.get('TRTNUM') or field_map.get('TRTNO') or field_map.get('TRNO')

                wspd_col = field_map.get('WSPD')
                wsgd_col = field_map.get('WSGD')
                nstd_col = field_map.get('NSTD')

                if not (run_col and (wspd_col or wsgd_col or nstd_col)):
                    # If the CSV is present but doesn't include these columns, treat as not available.
                    return {'by_run': {}, 'by_trt': {}}

                for row in reader:
                    runno = _safe_int(row.get(run_col))
                    trtno = _safe_int(row.get(trt_col)) if trt_col else None
                    if runno is None:
                        continue

                    if runno not in by_run_acc:
                        by_run_acc[runno] = {'wspd': [], 'wsgd': [], 'nstd': []}
                    if trtno is not None and trtno not in by_trt_acc:
                        by_trt_acc[trtno] = {'wspd': [], 'wsgd': [], 'nstd': []}

                    wspd = _safe_float(row.get(wspd_col)) if wspd_col else None
                    wsgd = _safe_float(row.get(wsgd_col)) if wsgd_col else None
                    nstd = _safe_float(row.get(nstd_col)) if nstd_col else None

                    if wspd is not None:
                        by_run_acc[runno]['wspd'].append(wspd)
                        if trtno is not None:
                            by_trt_acc[trtno]['wspd'].append(wspd)
                    if wsgd is not None:
                        by_run_acc[runno]['wsgd'].append(wsgd)
                        if trtno is not None:
                            by_trt_acc[trtno]['wsgd'].append(wsgd)
                    if nstd is not None:
                        by_run_acc[runno]['nstd'].append(nstd)
                        if trtno is not None:
                            by_trt_acc[trtno]['nstd'].append(nstd)

            by_run = {
                runno: {'wspd': _mean(v['wspd']), 'wsgd': _mean(v['wsgd']), 'nstd': _mean(v['nstd'])}
                for runno, v in by_run_acc.items()
            }
            by_trt = {
                trtno: {'wspd': _mean(v['wspd']), 'wsgd': _mean(v['wsgd']), 'nstd': _mean(v['nstd'])}
                for trtno, v in by_trt_acc.items()
            }
            return {'by_run': by_run, 'by_trt': by_trt}

        except Exception as e:
            print(f"Warning: Error parsing plantgro.csv in {folder_path}: {e}", flush=True)
            return {'by_run': {}, 'by_trt': {}}

    # ---- 2) Fallback: parse PlantGro.OUT directly ----
    out_path_candidates = [
        os.path.join(folder_path, 'PlantGro.OUT'),
        os.path.join(folder_path, 'PLANTGRO.OUT'),
        os.path.join(folder_path, 'plantgro.OUT'),
        os.path.join(folder_path, 'plantgro.out'),
    ]
    out_path = next((p for p in out_path_candidates if os.path.exists(p)), None)
    if not out_path:
        return {'by_run': {}, 'by_trt': {}}

    by_run = {}
    by_trt = {}

    current_run = None
    current_trt = None
    headers = None
    col_idx = {}
    acc = {}  # {runno: {'wspd': [...], 'wsgd': [...], 'nstd': [...]}}

    def _add_value(runno, metric, val):
        if runno not in acc:
            acc[runno] = {'wspd': [], 'wsgd': [], 'nstd': []}
        acc[runno][metric].append(val)

    def _finalize(runno, trtno):
        if runno not in acc:
            return
        out = {
            'wspd': _mean(acc[runno]['wspd']),
            'wsgd': _mean(acc[runno]['wsgd']),
            'nstd': _mean(acc[runno]['nstd']),
        }
        by_run[int(runno)] = out
        if trtno is not None:
            by_trt[int(trtno)] = out

    try:
        with open(out_path, 'r', encoding='utf-8', errors='ignore') as f:
            for raw in f:
                s = raw.strip()
                if not s:
                    continue

                if s.startswith('*RUN'):
                    if current_run is not None:
                        _finalize(current_run, current_trt)

                    current_run = None
                    current_trt = None
                    headers = None
                    col_idx = {}

                    m = re.match(r"\*RUN\s+(\d+)", s)
                    if m:
                        current_run = int(m.group(1))
                    continue

                if 'TREATMENT' in s and current_trt is None:
                    m = re.search(r"\bTREATMENT\s+(\d+)", s)
                    if m:
                        current_trt = _safe_int(m.group(1))
                    continue

                if s.startswith('@'):
                    headers = [h.strip().strip("'") for h in s[1:].split()]
                    col_idx = {h: i for i, h in enumerate(headers)}
                    continue

                if headers and s[0].isdigit():
                    parts = s.split()
                    if len(parts) < len(headers):
                        continue

                    def _get_float(colname):
                        i = col_idx.get(colname)
                        if i is None:
                            return None
                        return _safe_float(parts[i])

                    wspd = _get_float('WSPD')
                    wsgd = _get_float('WSGD')
                    nstd = _get_float('NSTD')

                    if current_run is None:
                        continue

                    if wspd is not None:
                        _add_value(current_run, 'wspd', wspd)
                    if wsgd is not None:
                        _add_value(current_run, 'wsgd', wsgd)
                    if nstd is not None:
                        _add_value(current_run, 'nstd', nstd)

        if current_run is not None:
            _finalize(current_run, current_trt)

    except Exception as e:
        print(f"Warning: Error parsing PlantGro.OUT in {folder_path}: {e}", flush=True)
        return {'by_run': {}, 'by_trt': {}}

    return {'by_run': by_run, 'by_trt': by_trt}


# --- Core Functions ---
def find_sqx_file(folder_path):
    for ext in ["*.SQX", "*.sqx"]:
        sqx_files = glob.glob(os.path.join(folder_path, ext))
        if sqx_files:
            return os.path.basename(sqx_files[0])
    return None


# DSSAT seasonal "FileX" experiment files follow the pattern <name>.??X, e.g.
# MZX (maize), WHX (wheat), SBX (soybean), RIX (rice), SGX (sorghum), BAX
# (barley). SQX is the sequence/rotation variant handled by find_sqx_file().
def find_experiment_file(folder_path):
    """Return the FileX (seasonal .??X) basename in *folder_path*, or None.

    Excludes .SQX (which is the sequence-mode file). Used by experiment mode.
    """
    for pattern in ["*.???", "*.???".lower()]:
        for fp in sorted(glob.glob(os.path.join(folder_path, pattern))):
            base = os.path.basename(fp)
            ext = base.rsplit(".", 1)[-1].upper()
            if len(ext) == 3 and ext.endswith("X") and ext != "SQX":
                return base
    return None


def find_dssat_input_file(folder_path, run_mode):
    """Pick the right DSSAT input file for the run mode.

    sequence   -> .SQX  (falls back to FileX if no .SQX present)
    experiment -> FileX .??X  (falls back to .SQX if no FileX present)
    Returns (filename, mode_flag) where mode_flag is the DSSAT run-mode letter
    ('Q' for sequence batch, 'B' for experiment/seasonal batch).
    """
    if str(run_mode).lower() == "experiment":
        filex = find_experiment_file(folder_path)
        if filex:
            return filex, "B"
        sqx = find_sqx_file(folder_path)        # graceful fallback
        if sqx:
            return sqx, "Q"
        return None, None
    else:  # sequence
        sqx = find_sqx_file(folder_path)
        if sqx:
            return sqx, "Q"
        filex = find_experiment_file(folder_path)  # graceful fallback
        if filex:
            return filex, "B"
        return None, None

def write_dssbatch_file(folder_path, sqx_file_name, config, treatment):
    if not sqx_file_name: return False
    batch_file_path = os.path.join(folder_path, config['control_file_name'])
    content_lines = []
    content_lines.append(f"$BATCH")
    content_lines.append("@FILEX                                                                                        TRTNO     RP     SQ     OP     CO")

    if config['run_mode'] == "experiment":
        line = f"{sqx_file_name:<93s} {treatment:>5d} {1:>6d} {1:>6d} {1:>6d} {0:>6d}"
        content_lines.append(line)
    elif config['run_mode'] == "sequence":
        for seq_num in range(config['sequence_start'], config['sequence_end'] + 1):
            line = f"{sqx_file_name:<93s} {treatment:>5d} {1:>6d} {seq_num:>6d} {1:>6d} {0:>6d}"
            content_lines.append(line)
    
    try:
        with open(batch_file_path, 'w', encoding='utf-8') as f_batch:
            f_batch.write("\n".join(content_lines) + "\n")
        return True
    except Exception as e:
        print(f"Error writing {config['control_file_name']} in {folder_path}: {e}")
        return False

def process_dssat_outputs_no_pandas(folder_path, point_id, treatment, config):
    processed_rows_list = []
    try:
        summary_path = os.path.join(folder_path, 'summary.csv')
        if not os.path.exists(summary_path): return None # Explicit check
        
        summary_data = read_csv_to_dict_list(summary_path)
        if not summary_data: return None

        soil_org_grouped = get_grouped_csv_data(os.path.join(folder_path, 'soilorg.csv'), 'RUN')
        soil_ni_grouped = get_grouped_csv_data(os.path.join(folder_path, 'soilni.csv'), 'RUN')
        soil_wat_grouped = get_grouped_csv_data(os.path.join(folder_path, 'soilwat.csv'), 'RUN')

        stress_grouped = parse_plantgro_stress(folder_path)
        stress_by_run = stress_grouped.get('by_run', {})
        stress_by_trt = stress_grouped.get('by_trt', {})

        for summary_row in summary_data:
            run_val = summary_row.get('RUNNO')
            if not run_val:
                 run_val = str(summary_row.get('TRNO', '')).strip() + '_' + str(summary_row.get('SQ', '')).strip()

            pyear = str(summary_row.get('PDAT', ''))[:4]
            som_start, som_end, som_delta = None, None, None
            
            if run_val and run_val in soil_org_grouped and soil_org_grouped[run_val]:
                try:
                    start_val = soil_org_grouped[run_val][0].get('SOMCT')
                    end_val = soil_org_grouped[run_val][-1].get('SOMCT')
                    if start_val and end_val:
                        som_start = float(start_val)
                        som_end = float(end_val)
                        som_delta = som_end - som_start
                except ValueError: pass

            napc, nlcc, ni_m = None, None, None
            if run_val and run_val in soil_ni_grouped and soil_ni_grouped[run_val]:
                last_ni_row = soil_ni_grouped[run_val][-1]
                napc = last_ni_row.get('NAPC')
                nlcc = last_ni_row.get('NLCC')
                ni_m = last_ni_row.get('NI#M')

            ir_c, irrc = None, None
            if run_val and run_val in soil_wat_grouped and soil_wat_grouped[run_val]:
                last_wat_row = soil_wat_grouped[run_val][-1]
                ir_c = last_wat_row.get('IR#C')
                irrc = last_wat_row.get('IRRC')


            # Seasonal-average stress indicators from PlantGro.OUT (mean across daily rows)
            wspd_avg = None
            wsgd_avg = None
            nstd_avg = None
            stress_rec = None

            # Prefer matching by DSSAT run number when available; fall back to treatment number
            run_no_int = None
            try:
                run_no_int = int(str(run_val).strip())
            except Exception:
                run_no_int = None

            if run_no_int is not None and run_no_int in stress_by_run:
                stress_rec = stress_by_run.get(run_no_int)
            else:
                try:
                    trt_int = int(treatment)
                except Exception:
                    trt_int = None
                if trt_int is not None and trt_int in stress_by_trt:
                    stress_rec = stress_by_trt.get(trt_int)

            if isinstance(stress_rec, dict):
                wspd_avg = stress_rec.get('wspd')
                wsgd_avg = stress_rec.get('wsgd')
                nstd_avg = stress_rec.get('nstd')

            processed_row = {
                'point_id': point_id,
                'treatment': treatment,
                'crop_code': summary_row.get('CR', ''),
                'latitude': summary_row.get('LAT', ''),
                'longitude': summary_row.get('LONG', ''),
                'weather_station_id': summary_row.get('WSTA', ''),
                'soil_profile_id': summary_row.get('SOIL_ID', ''),
                'dssat_file_id': summary_row.get('EXNAME', ''),
                'dssat_description': summary_row.get('TNAM', ''),
                'planting_date': summary_row.get('PDAT', ''),
                'emergence_date': summary_row.get('EDAT', ''),
                'harvest_date': summary_row.get('HDAT', ''),
                'year_planting': pyear,
                'year_harvest': summary_row.get('HYEAR', '')[:4],
                'top_weight_kg_ha': summary_row.get('CWAM', ''),
                'final_grain_kg_ha': summary_row.get('HWAM', ''),
                'removed_residue_kg_ha': summary_row.get('BWAH', ''),
                'carbon_soil_organic_matter_start_kg_ha': som_start,
                'carbon_soil_organic_matter_end_kg_ha': som_end,
                'carbon_soil_organic_matter_delta_kg_ha': som_delta,
                'final_irrigation_applications_count': ir_c,
                'final_irrigation_amount_mm': irrc,
                'inorganic_n_applied_count': ni_m,
                'inorganic_n_applied_kg_ha': napc,
                'nitrate_leaching_kg_ha': nlcc,
                'cumulative_net_co2_emissions_kg_CO2_ha': summary_row.get('CO2EM', ''),
                'cumulative_n2o_emissions_kg_N_ha': summary_row.get('N2OEM', ''),
                'water_stress_photosynthesis_avg': wspd_avg,
                'water_stress_development_avg': wsgd_avg,
                'nitrogen_stress_avg': nstd_avg
            }
            processed_rows_list.append(processed_row)
        return processed_rows_list
    except Exception as e:
        print(f"Process {rank}: Error processing outputs in {folder_path}: {e}")
        return None


def get_dssat_output_patterns(config):
    """
    Return the list of DSSAT output files we consider transient / cleanable.
    Includes the DSSBatch control file for reproducibility.
    """
    patterns = list(DSSAT_OUTPUT_PATTERNS_BASE)
    ctl = config.get('control_file_name')
    if ctl:
        patterns.append(str(ctl))
    return patterns


def _unique_dest_path(dest_path):
    """
    If dest_path already exists, append _1, _2, ... before the extension.
    """
    if not os.path.exists(dest_path):
        return dest_path

    base_dir = os.path.dirname(dest_path)
    base_name = os.path.basename(dest_path)
    stem, ext = os.path.splitext(base_name)

    i = 1
    while True:
        candidate = os.path.join(base_dir, f"{stem}_{i}{ext}")
        if not os.path.exists(candidate):
            return candidate
        i += 1


def _choose_archive_dir(folder_path, treatment, config):
    """
    Create an archive directory for this treatment run, avoiding collisions across re-runs.
    Structure: <point_folder>/_archive_outputs/<run_label>[_rN]/
    """
    archive_root = os.path.join(folder_path, "_archive_outputs")

    run_label = f"trt_{int(treatment):03d}"
    if str(config.get('run_mode', '')).lower() == 'sequence':
        try:
            s0 = int(config.get('sequence_start', 0))
            s1 = int(config.get('sequence_end', 0))
            run_label += f"_seq_{s0:03d}-{s1:03d}"
        except Exception:
            pass

    base_dir = os.path.join(archive_root, run_label)

    # If the base directory already exists and contains files, avoid overwriting by creating a new suffix
    if os.path.isdir(base_dir) and os.listdir(base_dir):
        i = 1
        while True:
            candidate = f"{base_dir}_r{i}"
            if not os.path.exists(candidate):
                base_dir = candidate
                break
            i += 1

    os.makedirs(base_dir, exist_ok=True)
    return base_dir


def archive_dssat_outputs(folder_path, point_id, treatment, config):
    """
    Move DSSAT transient outputs into an archive subfolder so you can keep per-treatment
    diagnostics even when cleanup is enabled.

    This is intentionally a MOVE (not copy) to avoid doubling storage.
    """
    try:
        patterns = get_dssat_output_patterns(config)
        if not patterns:
            return

        archive_dir = _choose_archive_dir(folder_path, treatment, config)

        results_file_to_keep = f"results_{point_id}.csv"
        moved = set()

        for pattern in patterns:
            for file_path in glob.glob(os.path.join(folder_path, pattern)):
                if not os.path.isfile(file_path):
                    continue
                base = os.path.basename(file_path)

                # Never move the per-point aggregated results file
                if base == results_file_to_keep:
                    continue

                # Avoid trying to move the same file twice if multiple patterns match
                if file_path in moved:
                    continue

                dest = os.path.join(archive_dir, base)
                dest = _unique_dest_path(dest)

                try:
                    shutil.move(file_path, dest)
                    moved.add(file_path)
                except Exception as e:
                    print(f"Rank {rank}: Warning - could not archive {file_path} -> {dest}: {e}", flush=True)

    except Exception as e:
        print(f"Rank {rank}: Warning - archiving failed in {folder_path} (point {point_id}, trt {treatment}): {e}", flush=True)


def cleanup_folder_files(folder_path, point_id, config):
    results_file_to_keep = f"results_{point_id}.csv"
    patterns = get_dssat_output_patterns(config)

    for pattern in patterns:
        for file_path in glob.glob(os.path.join(folder_path, pattern)):
            if os.path.basename(file_path) == results_file_to_keep: continue
            try: os.remove(file_path)
            except: pass

def run_dssat_simulation_orchestrator(sim_folder_path, point_id, treatment, config):
    original_cwd = os.getcwd()
    if not os.path.isdir(sim_folder_path): return point_id, treatment, None, "Folder not found"

    success = False

    try:
        os.chdir(sim_folder_path)
        # Select the input file + DSSAT run-mode flag based on run_mode:
        #   sequence   -> .SQX  with mode 'Q'
        #   experiment -> FileX .??X with mode 'B'
        # (Previously this always required a .SQX and always ran mode 'Q', so
        # pure seasonal FileX experiments could not run on HPC.)
        input_file_name, mode_flag = find_dssat_input_file(
            sim_folder_path, config['run_mode'])
        if not input_file_name:
            want = (".SQX" if str(config['run_mode']).lower() == "sequence"
                    else "FileX (.MZX/.WHX/.SBX/...)")
            raise FileNotFoundError(f"No DSSAT input file found (expected {want})")

        write_dssbatch_file(sim_folder_path, input_file_name, config, treatment)
        if os.path.exists('LUN.LST'):
            try: os.remove('LUN.LST')
            except: pass

        # Prevent stale outputs from being parsed if cleanup is disabled or a previous run left files behind
        for stale_name in ('summary.csv', 'soilorg.csv', 'soilni.csv', 'soilwat.csv'):
            stale_path = os.path.join(sim_folder_path, stale_name)
            try:
                if os.path.exists(stale_path):
                    os.remove(stale_path)
            except Exception:
                pass

        cmd = [config['executable_path'], mode_flag, config['control_file_name']]
        proc = subprocess.run(
            cmd,
            check=True,
            cwd=sim_folder_path,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        result_data = process_dssat_outputs_no_pandas(sim_folder_path, point_id, treatment, config)
        if not result_data:
            # DSSAT may exit 0 even if expected summary files were not produced/parsed.
            debug_msg = 'DSSAT completed but produced no parsable outputs (missing/empty summary.csv?)'
            if proc.stdout:
                debug_msg += f" | stdout tail: {proc.stdout[-500:]}"
            if proc.stderr:
                debug_msg += f" | stderr tail: {proc.stderr[-500:]}"
            return point_id, treatment, None, debug_msg

        success = True
        return point_id, treatment, result_data, None

    except subprocess.CalledProcessError as e:
        stderr = (e.stderr or '').strip()
        stdout = (e.stdout or '').strip()
        msg = 'DSSAT Runtime Error'
        if stderr:
            msg += f" | stderr tail: {stderr[-500:]}"
        if stdout:
            msg += f" | stdout tail: {stdout[-500:]}"
        return point_id, treatment, None, msg
    except Exception as e:
        return point_id, treatment, None, str(e)
    finally:
        # Optional archiving (controlled via --archive_outputs)
        if bool(config.get('archive_outputs', False)):
            archive_dssat_outputs(sim_folder_path, point_id, treatment, config)

        # Optional cleanup policy (controlled via --cleanup_mode)
        mode = str(config.get('cleanup_mode', 'always')).strip().lower()
        do_cleanup = (mode == 'always') or (mode == 'success' and success)

        if do_cleanup:
            cleanup_folder_files(sim_folder_path, point_id, config)

        if os.path.exists(original_cwd):
            os.chdir(original_cwd)

def process_point_sequentially(sim_folder_path, point_id, config):
    original_cwd = os.getcwd()
    point_data = []
    
    # Check directory existence first
    if not os.path.isdir(sim_folder_path):
        print(f"Rank {rank}: Error - Directory {sim_folder_path} does not exist.")
        return []

    for trt in range(config['treatment_start'], config['treatment_end'] + 1):
        pid, trt, data, err = run_dssat_simulation_orchestrator(sim_folder_path, point_id, trt, config)
        if err:
            # FIX: Print error so user knows why it failed
            print(f"Rank {rank}: [Failure] Point {point_id} Trt {trt} -> {err}")
        if data: 
            point_data.extend(data)
           
    if point_data:
        #local_csv_name = f"results_{point_id}.csv"
        local_csv_name = os.path.join(sim_folder_path, f"results_{point_id}.csv")
        write_dict_list_to_csv(point_data, local_csv_name, OUTPUT_FIELDNAMES)

    return point_data

# --- MPI message tags for the manager/worker protocol ---
TAG_READY = 11   # worker -> manager: "send me work"
TAG_WORK = 12    # manager -> worker: a folder index to process
TAG_STOP = 13    # manager -> worker: no more work, shut down


def _open_part_writer(part_path):
    """Open a streaming CSV writer (header written immediately, rows flushed
    per point). Returns (file_handle, csv_writer)."""
    os.makedirs(os.path.dirname(part_path), exist_ok=True)
    fh = open(part_path, 'w', newline='', encoding='utf-8')
    writer = csv.DictWriter(fh, fieldnames=OUTPUT_FIELDNAMES, extrasaction='ignore')
    writer.writeheader()
    fh.flush()
    return fh, writer


def _process_and_stream(info, config, writer, fh):
    """Process one point and stream its rows straight to disk (no rank-wide
    in-memory accumulation). Returns number of rows written."""
    try:
        rows = process_point_sequentially(info['path'], info['id'], config)
    except Exception as e:  # noqa: BLE001
        print(f"Rank {rank}: Critical Error {info['id']}: {e}", flush=True)
        return 0
    if rows:
        writer.writerows(rows)
        fh.flush()           # durable progress; nothing lost if a rank dies late
        return len(rows)
    return 0


# --- Main Execution ---
if __name__ == "__main__":
    config_params = {
        'dssat_simulation_output_base': dssat_simulation_output_base,
        'dssat_simulation_folder': dssat_simulation_folder,
        'dssat_simulation_summary_folder': dssat_simulation_summary_folder,
        'executable_path': executable_path,
        'dssat_model_code': dssat_model_code,
        'control_file_name': control_file_name,
        'run_mode': run_mode,
        'treatment_start': treatment_start, 'treatment_end': treatment_end,
        'sequence_start': sequence_start, 'sequence_end': sequence_end,
        'cleanup_mode': cleanup_mode,
        'archive_outputs': archive_outputs,
    }

    # Per-rank part files stream to node-local scratch when provided, then are
    # staged to the shared summary dir at the end. This keeps the per-point
    # output churn off the shared parallel filesystem (Lustre/GPFS metadata).
    if scratch_dir:
        local_part_dir = os.path.join(scratch_dir, dssat_simulation_folder, "parts")
    else:
        local_part_dir = project_temp_results_dir

    # Ensure shared output directories exist before ranks start writing files
    if rank == 0:
        os.makedirs(dssat_simulation_summary_folder, exist_ok=True)
        os.makedirs(project_temp_results_dir, exist_ok=True)
    comm.Barrier()

    # 1. Discovery (rank 0) + broadcast of the full task list
    all_folders_info = []
    if rank == 0:
        print("Rank 0: Scanning folders...", flush=True)
        if os.path.exists(config_params['dssat_simulation_output_base']):
            items = os.listdir(config_params['dssat_simulation_output_base'])
            for item in items:
                p = os.path.join(config_params['dssat_simulation_output_base'], item)
                if os.path.isdir(p) and item.isdigit():
                    all_folders_info.append({'path': p, 'id': item})
        all_folders_info.sort(key=lambda x: int(x['id']))
        print(f"Rank 0: Found {len(all_folders_info)} folders.", flush=True)

    all_folders_info = comm.bcast(all_folders_info, root=0)
    n_tasks = len(all_folders_info)

    part_file = os.path.join(local_part_dir, f"temp_results_rank_{rank}.csv")

    # ---------------------------------------------------------------------
    # 2. Processing
    #    Dynamic work-stealing (manager/worker) replaces static rank::size
    #    slicing so a few slow points can't stall a whole rank's share while
    #    other ranks sit idle. Manager (rank 0) hands out one folder index at
    #    a time on demand. Falls back to single-process loop when size == 1.
    # ---------------------------------------------------------------------
    rows_written = 0

    if size == 1:
        # Single process: just stream through every task.
        fh, writer = _open_part_writer(part_file)
        for info in all_folders_info:
            rows_written += _process_and_stream(info, config_params, writer, fh)
        fh.close()

    elif rank == 0:
        # ---- Manager: dispatch indices, never processes points itself ----
        next_idx = 0
        n_workers = size - 1
        active = n_workers
        status = MPI.Status()
        print(f"Rank 0: Manager dispatching {n_tasks} tasks to {n_workers} workers.", flush=True)
        while active > 0:
            # Wait for any worker to announce readiness.
            comm.recv(source=MPI.ANY_SOURCE, tag=TAG_READY, status=status)
            wsrc = status.Get_source()
            if next_idx < n_tasks:
                comm.send(next_idx, dest=wsrc, tag=TAG_WORK)
                next_idx += 1
            else:
                comm.send(None, dest=wsrc, tag=TAG_STOP)
                active -= 1
        print("Rank 0: All tasks dispatched.", flush=True)

    else:
        # ---- Worker: request work, stream results, until told to stop ----
        fh, writer = _open_part_writer(part_file)
        while True:
            comm.send(None, dest=0, tag=TAG_READY)
            idx = comm.recv(source=0, tag=MPI.ANY_TAG, status=MPI.Status())
            if idx is None:
                break
            rows_written += _process_and_stream(
                all_folders_info[idx], config_params, writer, fh)
        fh.close()

    comm.Barrier()

    # ---------------------------------------------------------------------
    # 3. Stage parts from node-local scratch to the shared summary dir.
    #    Each rank stages its OWN part (parallel), so this is not a rank-0
    #    serial bottleneck.
    # ---------------------------------------------------------------------
    staged_part = os.path.join(project_temp_results_dir, f"temp_results_rank_{rank}.csv")
    if scratch_dir and os.path.exists(part_file):
        try:
            os.makedirs(project_temp_results_dir, exist_ok=True)
            shutil.move(part_file, staged_part)
        except Exception as e:  # noqa: BLE001
            print(f"Rank {rank}: Warning - could not stage part file: {e}", flush=True)

    comm.Barrier()

    # ---------------------------------------------------------------------
    # 4. Final merge.
    #    merge_mode='none' leaves per-rank parts as the deliverable (avoids a
    #    rank-0 serial read of the full dataset on very large grids).
    #    merge_mode='concat' streams parts into one results file.
    # ---------------------------------------------------------------------
    if rank == 0:
        parts = sorted(glob.glob(os.path.join(project_temp_results_dir, "temp_results_rank_*.csv")))
        if merge_mode == 'none':
            print(f"\n--- merge_mode=none: leaving {len(parts)} per-rank part files in "
                  f"{project_temp_results_dir} ---", flush=True)
        else:
            print("\n--- Merging Results (streaming concat) ---", flush=True)
            final_file = os.path.join(dssat_simulation_summary_folder,
                                      f"results_{dssat_simulation_folder}.csv")
            with open(final_file, 'w', newline='', encoding='utf-8') as f_out:
                writer = csv.DictWriter(f_out, fieldnames=OUTPUT_FIELDNAMES, extrasaction='ignore')
                writer.writeheader()
                for p_file in parts:
                    try:
                        with open(p_file, 'r', encoding='utf-8') as f_in:
                            reader = csv.DictReader(f_in)
                            for row in reader:
                                writer.writerow(row)
                    except Exception as e:  # noqa: BLE001
                        print(f"Error reading {p_file}: {e}")
                    try:
                        os.remove(p_file)
                    except Exception:
                        pass
            print(f"Saved: {final_file}", flush=True)

    comm.Barrier()
