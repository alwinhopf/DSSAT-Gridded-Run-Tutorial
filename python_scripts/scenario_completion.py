"""Complete-point resume and atomic scenario accounting, independent of caches."""
import json
import os
from pathlib import Path
import tempfile
import hashlib
import copy


def resume_identity(identity):
    value=copy.deepcopy(identity)
    resolved=value.get('resolved_inputs',{})
    for key in list(resolved):
        if key.startswith(('weather_wth','soil_sol','soil_mapping')):
            del resolved[key]
    return value


def point_input_signatures(grid, id_col, soil_dir, weather_dir):
    result={}
    wgs=grid.to_crs(4326)
    def digest(path):
        try: return hashlib.md5(Path(path).read_bytes()).hexdigest()
        except OSError: return None
    for _,row in wgs.iterrows():
        pid=str(row[id_col])
        result[pid]=dict(soil=digest(Path(soil_dir)/(pid+'.SOL')),
                        weather=digest(Path(weather_dir)/(pid+'.WTH')),
                        longitude=float(row.geometry.x),latitude=float(row.geometry.y))
    return result


def atomic_json(value, path):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=path.name + '.', dir=path.parent)
    try:
        with os.fdopen(fd, 'w', encoding='utf-8') as fh:
            json.dump(value, fh, indent=2)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp): os.unlink(tmp)


def start_completion(result_path):
    path=str(result_path)+'.status.json'
    try: value=json.loads(Path(path).read_text())
    except (OSError,ValueError): value={}
    value.update(status='running',pid=os.getpid())
    # Retain old signatures until new results are actually published.
    atomic_json(value,path)


def complete_point_ids(data, treatments):
    if data is None or data.empty or not {'point_id', 'treatment'}.issubset(data):
        return set()
    import pandas as pd
    result = set()
    required = set(int(t) for t in treatments)
    for pid, rows in data.groupby('point_id'):
        values = pd.to_numeric(rows.treatment, errors='coerce').dropna()
        got = set(values[values == values.round()].astype(int))
        if required.issubset(got): result.add(str(pid))
    return result


def restore_aggregate_points(result_path, run_dir, identity, treatments, expected_ids, input_signatures=None):
    """Recover successes after cleanup only when the scenario identity matches."""
    import pandas as pd
    path = Path(result_path)
    try:
        previous_identity = json.loads(Path(str(path) + '.manifest.json').read_text())
        reusable=None
        if previous_identity != identity:
            status=json.loads(Path(str(path)+'.status.json').read_text())
            if status.get('resume_identity') != resume_identity(identity): return set()
            old=status.get('point_input_signatures',{})
            reusable={pid for pid,signature in (input_signatures or {}).items()
                      if signature.get('soil') and signature.get('weather') and old.get(pid)==signature}
        data = pd.read_csv(path, dtype={'point_id': str, 'soil_profile_id': str,
                                       'weather_station_id': str, 'dssat_file_id': str})
    except (OSError, ValueError):
        return set()
    complete = complete_point_ids(data, treatments) & set(map(str, expected_ids))
    if reusable is not None: complete &= reusable
    for pid in complete:
        folder = Path(run_dir) / pid; folder.mkdir(parents=True, exist_ok=True)
        rows = data[data.point_id.astype(str) == pid]
        rows.to_csv(folder / f'results_{pid}.csv', index=False)
    return complete


def validate_cached_point_results(run_dir, expected_ids, treatments):
    """Retain incomplete old CSVs for diagnosis while allowing an actual retry."""
    import pandas as pd
    for pid in map(str, expected_ids):
        path = Path(run_dir) / pid / f'results_{pid}.csv'
        if not path.exists(): continue
        try:
            data = pd.read_csv(path,dtype={'point_id':str})
            valid = pid in complete_point_ids(data,treatments)
        except (OSError,ValueError): valid = False
        if not valid:
            fd, archive = tempfile.mkstemp(prefix=f'incomplete_results_{pid}_',suffix='.csv',dir=path.parent)
            os.close(fd)
            os.replace(path,archive)


def publish_completion(result_path, data, expected_ids, treatments, excluded_ids=(), identity=None, input_signatures=None):
    expected = set(map(str, expected_ids))
    complete = complete_point_ids(data, treatments) & expected
    excluded = set(map(str, excluded_ids)) & expected - complete
    failures = expected - complete - excluded
    value = dict(schema=1, status='partial' if failures else 'complete',
                 expected_ids=sorted(expected), completed_ids=sorted(complete),
                 excluded_ids=sorted(excluded), failed_ids=sorted(failures),
                 treatments=list(map(int, treatments)))
    if identity is not None:
        value['resume_identity']=resume_identity(identity)
        value['point_input_signatures']=input_signatures or {}
    atomic_json(value, str(result_path) + '.status.json')
    return value


def scenario_is_complete(result_path):
    try:
        value = json.loads(Path(str(result_path) + '.status.json').read_text())
        import pandas as pd
        data = pd.read_csv(result_path, dtype={'point_id': str})
        actual = complete_point_ids(data, value['treatments'])
        expected = set(value['expected_ids'])
        return bool(value.get('schema') == 1 and value.get('status') == 'complete'
                    and expected and not value.get('failed_ids')
                    and set(value['completed_ids']).issubset(actual)
                    and expected <= actual | set(value['excluded_ids']))
    except (OSError, ValueError, KeyError, TypeError):
        return False
