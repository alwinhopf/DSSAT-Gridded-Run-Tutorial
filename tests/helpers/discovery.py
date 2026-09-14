"""Centralized executable discovery utilities for test suites.

Enforces AGENTS.md §5:
"Never invoke external tools by bare names: Centralize discovery for Rscript,
DSSAT (dscsm048), MPI (mpirun/mpiexec). When launching subprocesses, always
pass the resolved absolute Path or string, never bare 'Rscript' or 'dscsm048'."
"""
import os
import shutil
import sys
from pathlib import Path
from typing import Optional


def find_rscript() -> Optional[str]:
    """Resolve the absolute path to the Rscript executable.

    Checks:
    1. RSCRIPT and RSCRIPT_PATH environment variables.
    2. PATH lookup via shutil.which for 'Rscript' and 'Rscript.exe'.
    3. Common platform-specific installation locations (Windows, macOS, Linux).

    Returns:
        Absolute path to the executable as a string, or None if not found.
    """
    for env_var in ("RSCRIPT", "RSCRIPT_PATH"):
        env_rscript = os.environ.get(env_var)
        if env_rscript:
            candidate = Path(env_rscript).resolve()
            if candidate.is_file() and os.access(candidate, os.X_OK):
                return str(candidate)

    # Standard PATH search (handles case and extensions on Windows)
    for name in ("Rscript", "Rscript.exe", "rscript"):
        found = shutil.which(name)
        if found:
            resolved = Path(found).resolve()
            if resolved.is_file() and os.access(resolved, os.X_OK):
                return str(resolved)

    # Platform-specific standard installation locations
    candidates: list[Path] = []
    if sys.platform == "win32":
        for root_var in ("ProgramFiles", "ProgramFiles(x86)"):
            pf = os.environ.get(root_var, r"C:\Program Files")
            r_base = Path(pf) / "R"
            if r_base.is_dir():
                for r_dir in sorted(r_base.glob("R-*"), reverse=True):
                    candidates.append(r_dir / "bin" / "Rscript.exe")
                    candidates.append(r_dir / "bin" / "x64" / "Rscript.exe")
    elif sys.platform == "darwin":
        candidates.extend([
            Path("/usr/local/bin/Rscript"),
            Path("/opt/homebrew/bin/Rscript"),
            Path("/Library/Frameworks/R.framework/Resources/bin/Rscript"),
            Path("/usr/bin/Rscript"),
        ])
    else:  # Linux / Unix
        candidates.extend([
            Path("/usr/bin/Rscript"),
            Path("/usr/local/bin/Rscript"),
        ])

    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate.resolve())

    return None
