#!/bin/bash
#
# Example SLURM script for running DSSAT MPI jobs using the Python runner.
# Submit with: sbatch hpc/run_dssat_python.slurm
#
# --- RESOURCE REQUESTS (edit for your cluster) ---
#SBATCH --time=6:00:00
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=110
#SBATCH -J "dssat_grid"
#SBATCH --mail-type=BEGIN,END,FAIL
##SBATCH --mail-user=<YOUR_EMAIL>
##SBATCH --account=<YOUR_ACCOUNT_OR_WCKEY>

# --- USER CONFIGURATION (edit these) ---

# One or more DSSAT run directories (each contains point folders like 00000001/, 00000002/, ...)
# Tip: keep these on fast scratch storage.
RUN_DIRS=(
  "/scratch/$USER/<RUN_NAME_1>"
  # "/scratch/$USER/<RUN_NAME_2>"
)

# Where to write the final merged CSV results_<RUN_FOLDER>.csv
SUMMARY_DIR="${SUMMARY_DIR:-$PWD/results}"

# DSSAT executable (either put dscsm048 on your PATH, set DSSAT_EXE env var, or hardcode here)
DSSAT_EXE="${DSSAT_EXE:-/path/to/dscsm048}"

# Cleanup behavior for point folders: always | success | never
CLEANUP_MODE="${CLEANUP_MODE:-always}"

# Set to 1 to move DSSAT transient outputs into per-point _archive_outputs/ folders (per treatment)
ARCHIVE_OUTPUTS="${ARCHIVE_OUTPUTS:-0}"

# Node-local scratch for per-rank part files (avoids small-file metadata storm
# on shared Lustre/GPFS). SLURM sets $TMPDIR per node; parts are staged back to
# SUMMARY_DIR at the end. Set blank to write parts directly to SUMMARY_DIR.
SCRATCH_DIR="${SCRATCH_DIR:-$TMPDIR}"

# Final merge policy: concat (single results CSV) | none (keep per-rank parts).
# Use 'none' on very large grids to skip the rank-0 serial merge.
MERGE_MODE="${MERGE_MODE:-concat}"

# DSSAT run settings
MODEL_CODE="${MODEL_CODE:-CRGRO048}"
TRT_START="${TRT_START:-1}"
TRT_END="${TRT_END:-1}"
SEQ_START="${SEQ_START:-1}"
SEQ_END="${SEQ_END:-1}"

# Path to the Python runner in this repo
PY_RUNNER="hpc/dssat_mpi_runner.py"

# --- ENVIRONMENT SETUP (cluster-specific) ---
module load python
eval "$(conda shell.bash hook)"
# Change this to your environment name/path:
conda activate <YOUR_CONDA_ENV>

export OMP_NUM_THREADS=1
export PYTHONUNBUFFERED=1

# Archive flag logic
ARCHIVE_FLAG=""
if [[ "$ARCHIVE_OUTPUTS" == "1" ]]; then
  ARCHIVE_FLAG="--archive_outputs"
fi

# --- RUN ---
echo "Date: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Host: $(hostname)"
echo "Working Directory: $(pwd)"
echo "Nodes: $SLURM_NNODES"
echo "Tasks: $SLURM_NTASKS"

mkdir -p "$SUMMARY_DIR"

# NOTE:
# Some clusters prefer srun instead of mpirun. If so, replace mpirun with:
#   srun python "$PY_RUNNER" ...
#
# Here we use mpirun and default to SLURM_NTASKS if available.
NPROCS="${SLURM_NTASKS:-550}"

for RUN_DIR in "${RUN_DIRS[@]}"; do
  echo "============================================================"
  echo "Running DSSAT MPI job for: $RUN_DIR"
  echo "============================================================"

  mpirun -n "$NPROCS" python "$PY_RUNNER" \
    --base_dir "$RUN_DIR" \
    --summary_dir "$SUMMARY_DIR" \
    --exe_path "$DSSAT_EXE" \
    --cleanup_mode "$CLEANUP_MODE" \
    --scratch_dir "$SCRATCH_DIR" \
    --merge_mode "$MERGE_MODE" \
    $ARCHIVE_FLAG \
    --run_mode sequence \
    --model_code "$MODEL_CODE" \
    --trt_start "$TRT_START" \
    --trt_end "$TRT_END" \
    --seq_start "$SEQ_START" \
    --seq_end "$SEQ_END"

  rc=$?
  echo "Run finished with exit code $rc at: $(date)"
  if [[ "$rc" != "0" ]]; then
    echo "WARNING: Non-zero exit code. Check SLURM output for details."
  fi
done
