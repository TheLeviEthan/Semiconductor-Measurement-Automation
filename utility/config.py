"""
Filename: config.py
Author: Ethan Ruddell
Date: 2026-02-12
Description: Load and manage YAML-based project configuration.

This module reads the file "config.yaml" that lives in the project root folder.
That YAML file stores default settings like GPIB addresses, start/stop
frequencies, ramp rates, etc. so they don't have to be typed in every time.

If config.yaml is missing or a key isn't defined, the code simply falls back
to hard-coded defaults elsewhere in the project — nothing will crash.

How other files use this module:
    import config
    addr = config.get("gpib", "pia", "GPIB0::24::INSTR")  # section, key, fallback
"""

import os
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Locate config.yaml — look next to main.py (project root).
# sys.argv[0] is the script that was launched (main.py), so its parent
# directory is the project root where config.yaml should live.
# ---------------------------------------------------------------------------
_PROJECT_ROOT = Path(sys.argv[0]).resolve().parent
_CONFIG_PATH = _PROJECT_ROOT / "config.yaml"

_config = {}

def _load_yaml():
    """Read config.yaml from disk and store it in the module-level dictionary.

    If the file is missing, if PyYAML isn't installed, or if the file is
    malformed, we quietly fall back to an empty dictionary so that the
    rest of the application still works with its own built-in defaults.
    """
    global _config
    try:
        import yaml  # optional dependency
        with open(_CONFIG_PATH, "r", encoding="utf-8") as fh:
            _config = yaml.safe_load(fh) or {}
    except FileNotFoundError:
        _config = {}
    except ImportError:
        # PyYAML not installed — fall back to empty config
        _config = {}
    except Exception:
        _config = {}

_load_yaml()


# ---------------------------------------------------------------------------
# Public helpers — these are the functions other files call
# ---------------------------------------------------------------------------

def get(section: str, key: str, default=None):
    """Retrieve ``config[section][key]``, returning *default* on miss.

    Example::

        gpib_addr = config.get("gpib", "pia", "GPIB0::24::INSTR")
    """
    try:
        return _config[section][key]
    except (KeyError, TypeError):
        return default


def get_section(section: str) -> dict:
    """Return the full dict for *section*, or ``{}``."""
    try:
        return dict(_config[section])
    except (KeyError, TypeError):
        return {}


def reload():
    """Re-read config.yaml from disk (e.g. after the user edits it)."""
    _load_yaml()
