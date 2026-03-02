"""
Filename: file_management.py
Author: Ethan Ruddell
Date: 2025-01-20
Description: File I/O helpers for the measurement automation suite.

After every measurement the software needs to save two things:
  1. A CSV data file    – columns of numbers (frequency, voltage, current, etc.)
  2. A plot image (PNG) – a graph of those numbers

This module handles both tasks, plus a few supporting chores:
  - All files are saved under an "output/" folder next to the project.
  - File names are automatically prefixed with a date-time stamp
    (e.g. "2026-03-02_14-05-12_impedance_data.csv") so that files
    are naturally sorted and no previous results are overwritten.
  - If the exact filename already exists (same second!), a counter
    (1), (2)… is appended as a fallback.
  - Plot images are created using an explicit Figure object instead
    of matplotlib's global pyplot, making them safe to generate from
    background threads (important for the GUI).
  - The output directory can be changed at runtime via the GUI
    "Browse" button or the --output-dir CLI flag.
"""

import os
import logging
import numpy as np
from datetime import datetime
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Output directory defaults
# ---------------------------------------------------------------------------
# By default, results go to a folder called "output" located one level above
# the project directory.  For example, if the project is at:
#     C:/Users/you/NRG Scripts/Semiconductor-Measurement-Automation/
# then the default output is:
#     C:/Users/you/NRG Scripts/output/
#
# We resolve this path relative to *this* file's location so it stays stable
# regardless of how the application is launched (double-click, CLI, GUI, etc.).
# ---------------------------------------------------------------------------
_UTILITY_DIR = Path(__file__).resolve().parent
_PROJECT_DIR = _UTILITY_DIR.parent            # Semiconductor-Measurement-Automation/
default_output_dir = str((_PROJECT_DIR / ".." / "output").resolve())

# Active output directory (may be changed at runtime by the GUI or CLI)
output_dir = default_output_dir


# =============================
# Path helpers
# =============================
# These small functions build file names and make sure directories exist
# before the code tries to write into them.

def _date_prefix() -> str:
    """Return a compact date-time prefix for file naming, e.g. '2026-03-02_14-05-12'."""
    return datetime.now().strftime("%Y-%m-%d_%H-%M-%S")


def uniquify(path_string: str, date_stamp: bool = True) -> str:
    """Generate a unique file path.

    When *date_stamp* is True (the default) the filename is prefixed with a
    date-time stamp so files are naturally ordered and collisions are rare.
    If the path already exists an incrementing counter ``(1)``, ``(2)`` … is
    appended as a fallback.

    Parameters
    ----------
    path_string : str
        Desired file path (may already exist).
    date_stamp : bool
        Prepend a ``YYYY-MM-DD_HH-MM-SS_`` prefix to the stem.
    """
    path = Path(path_string)
    directory = path.parent
    stem = path.stem
    ext = path.suffix

    if date_stamp:
        stem = f"{_date_prefix()}_{stem}"
        path = directory / f"{stem}{ext}"

    counter = 1
    while path.exists():
        path = directory / f"{stem} ({counter}){ext}"
        counter += 1

    return str(path)


def ensure_output_dir(path: str) -> None:
    """Create *path* (and parents) if it does not already exist."""
    os.makedirs(path, exist_ok=True)


def set_output_dir(path: str) -> None:
    """Override the active output directory with an absolute, expanded path."""
    global output_dir
    output_dir = str(Path(path).expanduser().resolve())

# =============================
# Save helpers
# =============================
# Functions that actually write data/images to disk.  They all:
#   1. Create the target subdirectory (csv/ or images/) if it doesn't exist.
#   2. Generate a unique, date-stamped filename.
#   3. Write the file.
#   4. Return the path to the created file.

def save_csv(filename: str, data, header: str) -> str:
    """Save a NumPy array to a date-stamped CSV in ``output/csv/``.

    Returns the final file path.
    """
    csv_dir = os.path.join(output_dir, "csv")
    ensure_output_dir(csv_dir)
    csv_path = uniquify(os.path.join(csv_dir, filename))
    log.info("Saving CSV to: %s", csv_path)
    np.savetxt(csv_path, data, delimiter=",", header=header, comments="")
    return csv_path


def save_image(title: str, axis1_label: str, axis1, axis2_label: str, axis2,
               APPLY_DC_BIAS: bool = False, DC_BIAS_V: float = 0.0,
               dpi: int = 300) -> str:
    """Create a semilog-x plot and save it to ``output/images/``.

    Uses the explicit *Figure* / *Axes* API so it is safe to call from
    background threads (the GUI sets ``matplotlib.use('Agg')``).

    Parameters
    ----------
    title : str
        Plot title (also used in the filename).
    axis1_label, axis1 : str, array-like
        X-axis label and data.
    axis2_label, axis2 : str, array-like
        Y-axis label and data.
    APPLY_DC_BIAS : bool
        If True, the DC bias value is appended to the title.
    DC_BIAS_V : float
        DC bias voltage shown in the title when enabled.
    dpi : int
        Image resolution.

    Returns
    -------
    str
        The final path the image was saved to.
    """
    images_dir = os.path.join(output_dir, "images")
    ensure_output_dir(images_dir)

    if APPLY_DC_BIAS:
        title_with_bias = f"{title} (DC Bias = {DC_BIAS_V:.2f} V)"
    else:
        title_with_bias = title

    fig = Figure()
    ax = fig.add_subplot(111)
    ax.semilogx(axis1, axis2, '-o', markersize=3)
    ax.set_xlabel(axis1_label)
    ax.set_ylabel(axis2_label)
    ax.set_title(title_with_bias)
    ax.grid(True, which="both")

    image_path = uniquify(os.path.join(images_dir, f"{title_with_bias}_plot.png"))
    fig.savefig(image_path, dpi=dpi, bbox_inches="tight")
    log.info("Saved image to: %s", image_path)
    return image_path


def save_plot(filename: str, fig=None, dpi: int = 300) -> str:
    """Save an already-constructed matplotlib figure to ``output/images/``.

    Use this when the caller builds a custom plot (subplots, multi-trace,
    etc.) and only needs the path / uniquify logic handled.

    Parameters
    ----------
    filename : str
        Bare filename, e.g. ``"pulsed_iv.png"``.
    fig : matplotlib Figure, optional
        If *None* the current ``plt`` figure is used.
    dpi : int
        Image resolution.

    Returns
    -------
    str
        The final path the image was saved to.
    """
    images_dir = os.path.join(output_dir, "images")
    ensure_output_dir(images_dir)
    image_path = uniquify(os.path.join(images_dir, filename))
    if fig is None:
        plt.savefig(image_path, dpi=dpi, bbox_inches="tight")
    else:
        fig.savefig(image_path, dpi=dpi, bbox_inches="tight")
    log.info("Saved image to: %s", image_path)
    return image_path


def save_cycle_plot(title: str, x_label: str, y_label: str,
                    cycles_x, cycles_y, base_filename: str,
                    dpi: int = 300) -> str:
    """Save a multi-cycle linear overlay plot to ``output/images/``.

    Thread-safe — uses the explicit Figure API.
    """
    images_dir = os.path.join(output_dir, "images")
    ensure_output_dir(images_dir)

    fig = Figure()
    ax = fig.add_subplot(111)
    for idx, (x_vals, y_vals) in enumerate(zip(cycles_x, cycles_y), start=1):
        ax.plot(x_vals, y_vals, "-o", markersize=3, label=f"Cycle {idx}")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.grid(True)
    ax.legend()

    image_path = uniquify(os.path.join(images_dir, base_filename))
    fig.savefig(image_path, dpi=dpi, bbox_inches="tight")
    log.info("Saved image to: %s", image_path)
    return image_path
