"""
Filename: logging_config.py
Author: Ethan Ruddell
Date: 2026-02-12
Description: Centralized Python logging configuration for the project.

Usage::

    import logging_config
    logging_config.setup()          # call once at startup (main.py)

    import logging
    log = logging.getLogger(__name__)
    log.info("Measurement started")
"""

import logging
import os
import sys
from pathlib import Path

try:
    import config as _cfg
except ImportError:
    _cfg = None


def setup(level=None, log_to_file=None, log_filename=None, output_dir=None):
    """Configure the root logger with console + optional file handler.

    Parameters
    ----------
    level : str, optional
        Logging level name (DEBUG / INFO / WARNING / ERROR).
        Falls back to config.yaml → ``logging.level`` → ``"INFO"``.
    log_to_file : bool, optional
        Write log to disk?  Falls back to config.yaml → ``True``.
    log_filename : str, optional
        Log file name inside *output_dir*.
    output_dir : str, optional
        Directory for the log file.  Defaults to cwd.
    """
    # --- resolve settings from config.yaml → function args → hard defaults ---
    if _cfg is not None:
        _sec = _cfg.get_section("logging")
    else:
        _sec = {}

    level = level or _sec.get("level", "INFO")
    log_to_file = log_to_file if log_to_file is not None else _sec.get("log_to_file", True)
    log_filename = log_filename or _sec.get("log_filename", "measurement.log")

    numeric_level = getattr(logging, level.upper(), logging.INFO)

    fmt = "%(asctime)s | %(name)-24s | %(levelname)-7s | %(message)s"
    datefmt = "%Y-%m-%d %H:%M:%S"

    handlers = []

    # Console handler — always present
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(numeric_level)
    console.setFormatter(logging.Formatter(fmt, datefmt=datefmt))
    handlers.append(console)

    # File handler — optional
    if log_to_file:
        log_dir = output_dir or os.getcwd()
        os.makedirs(log_dir, exist_ok=True)
        fh = logging.FileHandler(
            os.path.join(log_dir, log_filename), encoding="utf-8"
        )
        fh.setLevel(logging.DEBUG)          # file always gets full detail
        fh.setFormatter(logging.Formatter(fmt, datefmt=datefmt))
        handlers.append(fh)

    logging.basicConfig(level=numeric_level, handlers=handlers, force=True)

    # Quiet down noisy libraries
    logging.getLogger("pyvisa").setLevel(logging.WARNING)

    logging.getLogger(__name__).debug("Logging initialised (level=%s, file=%s)", level, log_to_file)
