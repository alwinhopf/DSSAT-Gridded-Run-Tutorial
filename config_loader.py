"""config_loader.py

Lightweight, dependency-tolerant loader for the shared config.yml.

Imported near the top of SECTION 0 in dssat_main_pipeline.py. Exposes one
function, cfg_get(key, default), returning the value from config.yml if present,
otherwise the supplied default. If config.yml or PyYAML is missing, every call
returns the default, so the pipeline runs exactly as before.

This is the Python twin of config_loader.R: both read the same config.yml so the
R and Python SECTION-0 settings cannot drift apart.
"""
import os

_CONFIG = {}


def _find_and_load():
    here = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(here, "config.yml"),
        os.path.join(os.getcwd(), "config.yml"),
        "config.yml",
    ]
    path = next((p for p in candidates if os.path.isfile(p)), None)
    if path is None:
        print("[config_loader] No config.yml found - using in-script defaults.")
        return {}
    try:
        import yaml  # PyYAML
    except ImportError:
        print("[config_loader] PyYAML not installed - using in-script defaults.")
        return {}
    try:
        with open(path, "r") as fh:
            cfg = yaml.safe_load(fh) or {}
        print(f"[config_loader] Loaded {len(cfg)} settings from {path}")
        return cfg
    except Exception as exc:  # noqa: BLE001
        print(f"[config_loader] Failed to parse {path} ({exc}) - using defaults.")
        return {}


_CONFIG = _find_and_load()


def cfg_get(key, default):
    """Return config.yml[key], or `default` if absent or blank ("")."""
    if key not in _CONFIG or _CONFIG[key] is None:
        return default
    val = _CONFIG[key]
    if isinstance(val, str) and val == "":   # blank => default
        return default
    return val
