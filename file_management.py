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

# File name definitions

# Get the path to the current script directory
script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))

# Name a new "output" folder in the script directory
output_dir = os.path.join(script_dir, "output")


#=============================
# Helper functions
#=============================

measurementType = ["Impedance vs Theta", "Capacitance vs Tan Loss", "Dielectric Measurements"]

def uniquify(path_string):
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

def save_csv(path, data, header):
    """Save data array to a uniquified CSV file with header."""
    ensure_output_dir(output_dir)
    csv_path = uniquify(path)
    print(f"Saving CSV to: {csv_path}")
    np.savetxt(csv_path, data, delimiter=",", header=header, comments="")
    return csv_path

def save_image(title, axis1_label, axis1, axis2_label, axis2, APPLY_DC_BIAS=False, DC_BIAS_V=0.0):
    """Save a semilog plot with optional DC bias in title and filename."""
    ensure_output_dir(output_dir)
    
    import matplotlib.pyplot as plt

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

    image_path = uniquify(os.path.join(output_dir, f"{title_with_bias}_plot.png"))
    print(f"Saved image to: {image_path}")
    plt.savefig(image_path, dpi=300, bbox_inches="tight")
    plt.close()
    return image_path

