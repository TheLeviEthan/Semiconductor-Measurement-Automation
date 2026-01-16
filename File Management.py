"""
Filename: File Managment.py
Author: Ethan Ruddell
Date: 2025-01-16
Description: Contains helper functions for file management tasks, ensuring 
uniquiness of file names and allowing users to easily select desintations for 
output files.
"""
import os
import sys
from pathlib import Path

#=============================
# Helper functions
#=============================



def uniquify(path_string):
    """
    Generates a unique file path by appending an incrementing number 
    if the original file already exists.
    """
    path = Path(path_string)
    filename = path.stem
    extension = path.suffix
    counter = 1

    while path.exists():
        path_string = f"{filename} ({counter}){extension}"
        path = Path(path_string)
        counter += 1 # concatenate file name with (counter) before extension
    
    return str(path)

def ensure_output_dir(path):
        os.makedirs(path, exist_ok=True)

def save_csv(csv_name, save_path = "null"):
    # get path to current dir
    if (save_path == "null"):
        script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
    else:
        script_dir = save_path
    
        output_dir = os.path.join(script_dir, "csv_outputs")
        ensure_output_dir(output_dir) # make sure output directory exists