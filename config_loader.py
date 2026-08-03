"""Load the shared gridded-pipeline configuration.

``config.yml`` beside this module is the canonical set of defaults.  A file
named by ``DSSAT_CONFIG_FILE`` is a *partial override* merged over those
defaults.  This keeps study-specific YAML files small without falling back to
different hidden defaults in the R and Python pipelines.

This module is mirrored by :mod:`config_loader.R`; keep their precedence and
blank-value behavior aligned.
"""
import os

_CONFIG = {}


def _deep_merge(base, override):
    """Return a recursive mapping merge without mutating either input."""
    merged = dict(base)
    for key, value in override.items():
        if isinstance(merged.get(key), dict) and isinstance(value, dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def _read_yaml(path, yaml):
    with open(path, "r", encoding="utf-8-sig") as fh:
        cfg = yaml.safe_load(fh) or {}
    if not isinstance(cfg, dict):
        raise ValueError(f"Top level must be a YAML mapping: {path}")
    return {str(k).lstrip("\ufeff"): v for k, v in cfg.items()}


def _validate_override(base, override, prefix=""):
    for key, value in override.items():
        path = f"{prefix}.{key}" if prefix else key
        if key not in base:
            raise ValueError(f"Unknown configuration key: {path}")
        expected = base[key]
        if isinstance(expected, dict) and isinstance(value, dict):
            _validate_override(expected, value, path)
        elif value is not None and value != "" and expected is not None:
            if isinstance(expected, bool) and not isinstance(value, bool):
                raise TypeError(f"{path} must be boolean, got {type(value).__name__}")
            if isinstance(expected, (int, float)) and not isinstance(expected, bool) \
                    and not isinstance(value, (int, float)):
                raise TypeError(f"{path} must be numeric, got {type(value).__name__}")
            if isinstance(expected, list) and not isinstance(value, list):
                raise TypeError(f"{path} must be a list, got {type(value).__name__}")


def _find_and_load():
    here = os.path.dirname(os.path.abspath(__file__))
    base_path = os.path.join(here, "config.yml")
    env_path = os.environ.get("DSSAT_CONFIG_FILE", "")
    try:
        import yaml  # PyYAML
    except ImportError as exc:
        raise RuntimeError(
            "PyYAML is required to load the central config.yml; install the "
            "project requirements before running the pipeline."
        ) from exc
    if not os.path.isfile(base_path):
        raise FileNotFoundError(f"Canonical pipeline config not found: {base_path}")
    if env_path and not os.path.isfile(env_path):
        raise FileNotFoundError(
            f"DSSAT_CONFIG_FILE points to a missing file: {env_path}"
        )
    try:
        cfg = _read_yaml(base_path, yaml)
        loaded = [base_path]
        if env_path and os.path.abspath(env_path) != os.path.abspath(base_path):
            override = _read_yaml(env_path, yaml)
            _validate_override(cfg, override)
            cfg = _deep_merge(cfg, override)
            loaded.append(env_path)
        print(
            f"[config_loader] Loaded {len(cfg)} settings from "
            + " + ".join(loaded)
        )
        return cfg
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(f"Failed to load pipeline configuration: {exc}") from exc


_CONFIG = _find_and_load()


def cfg_get(key, default):
    """Return config.yml[key], or `default` if absent or blank ("")."""
    if key not in _CONFIG or _CONFIG[key] is None:
        return default
    val = _CONFIG[key]
    if isinstance(val, str) and val == "":   # blank => default
        return default
    return val
