"""
Filename: File Managment.py
Author: Ethan Ruddell
Date: 2025-01-20
Description: Contains helper functions for file management tasks, ensuring 
uniqueness of file names and allowing users to easily select destinations for 
output files.
"""
import os
import sys
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt


# File name definitions

# Get the path to the current script directory
script_dir = Path(sys.argv[0]).resolve().parent

# Default output folder lives one level up from the script directory
default_output_dir = str((script_dir / ".." / "output").resolve())

# Active output directory (may be overridden via CLI)
output_dir = default_output_dir


#=============================
# Helper functions
#=============================

measurementType = ["Impedance vs Theta", "Capacitance vs Tan Loss", "Dielectric Measurements"]

def uniquify(path_string):
    # TODO: date file naming
    """
    Generates a unique file path by appending an incrementing number 
    if the original file already exists.
    """
    path = Path(path_string)
    directory = path.parent
    filename = path.stem
    extension = path.suffix
    counter = 1

    while path.exists():
        path = directory / f"{filename} ({counter}){extension}"
        counter += 1

    return str(path)

def ensure_output_dir(path):
    """
    Checks if the output directory exists, and creates it if not.
    """
    os.makedirs(path, exist_ok=True)

def set_output_dir(path):
    """Override the active output directory using an absolute, expanded path."""
    global output_dir
    output_dir = str(Path(path).expanduser().resolve())

def save_csv(filename, data, header):
    """Save data array to a uniquified CSV file with header in output/csv folder."""
    ensure_output_dir(os.path.join(output_dir, "csv"))
    csv_path = uniquify(os.path.join(output_dir, "csv", filename))
    print(f"Saving CSV to: {csv_path}")
    np.savetxt(csv_path, data, delimiter=",", header=header, comments="")
    return csv_path

def save_image(title, axis1_label, axis1, axis2_label, axis2, APPLY_DC_BIAS=False, DC_BIAS_V=0.0):
    """Save a semilog plot with optional DC bias in title and filename."""

    # TODO : logarithmic AND linear plots as well
    ensure_output_dir(os.path.join(output_dir, "images"))
    

    # Augment title with DC bias BEFORE using it in filename/plot
    if APPLY_DC_BIAS:
        title_with_bias = f"{title} (DC Bias = {DC_BIAS_V:.2f} V)"
    else:
        title_with_bias = title

    plt.figure()
    plt.semilogx(axis1, axis2, '-o', markersize=3)
    plt.xlabel(axis1_label)
    plt.ylabel(axis2_label)
    plt.title(title_with_bias)
    plt.grid(True, which="both")

    image_path = uniquify(os.path.join(output_dir, "images", f"{title_with_bias}_plot.png"))
    print(f"Saved image to: {image_path}")
    plt.savefig(image_path, dpi=300, bbox_inches="tight")
    plt.close()
    return image_path

def save_cycle_plot(title, x_label, y_label, cycles_x, cycles_y, base_filename):
    """Save multi-cycle linear plot for stacked sweeps with unique filenames."""
    ensure_output_dir(os.path.join(output_dir, "images"))

    plt.figure()
    for idx, (x_vals, y_vals) in enumerate(zip(cycles_x, cycles_y), start=1):
        plt.plot(x_vals, y_vals, "-o", markersize=3, label=f"Cycle {idx}")
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.grid(True)
    plt.legend()

    image_path = uniquify(os.path.join(output_dir, "images", base_filename))
    print(f"Saved image to: {image_path}")
    plt.savefig(image_path, dpi=300, bbox_inches="tight")
    plt.close()
    return image_path
