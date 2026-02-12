"""
Filename: config.py
Author: Ethan Ruddell
Date: 2026-02-12
Description: Load and manage YAML-based project configuration.

Falls back gracefully to hard-coded defaults when config.yaml is
missing or incomplete.
"""

import os
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Locate config.yaml — look next to main.py (project root)
# ---------------------------------------------------------------------------
_PROJECT_ROOT = Path(sys.argv[0]).resolve().parent
_CONFIG_PATH = _PROJECT_ROOT / "config.yaml"

_config = {}

def _load_yaml():
    """Attempt to parse config.yaml.  Returns {} on any failure."""
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
# Public helpers
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
