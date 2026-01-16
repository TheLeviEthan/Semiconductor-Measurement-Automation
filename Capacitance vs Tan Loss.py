"""
Filename: Dielectric Measurements.py
Author: Chaitanya Sharma, Ethan Ruddell
Date: 2025-01-16
Description: Contains commands for capacitance vs tan loss measurements.
"""

import os
import sys
import pyvisa
import numpy as np
import matplotlib

# Use non-interactive backend so script finishes without waiting for plot windows
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pathlib import Path


# ==============================
# User settings (EDIT THESE)
# ==============================

GPIB_ADDRESS = "GPIB0::24::INSTR"  # Example: GPIB0, address 24

# Frequency sweep settings (LOG sweep)
FREQ_START_HZ = 1e3      # start frequency, Hz (>0 for LOG sweep)
FREQ_STOP_HZ = 1e7       # stop frequency, Hz
NUM_POINTS = 201         # number of points in LOG sweep

# DC bias options
APPLY_DC_BIAS = True     # True = apply a fixed DC bias; False = no DC bias
DC_BIAS_V = 0.0          # bias voltage (V) if APPLY_DC_BIAS = True

# Get the path to the current script directory
script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))

# Name a new "output" folder in the script directory
output_dir = os.path.join(script_dir, "output")

# Create the directory
try:
    os.makedirs(output_dir, exist_ok=True)
    print(f"Directory '{output_dir}' ensured.")
except OSError as e:
    print(f"Error creating directory {output_dir}: {e}")


# ==============================
# Helper functions
# ==============================


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
        counter += 1
    
    return str(path)


def ensure_output_dir(path):
    os.makedirs(path, exist_ok=True)


def connect_4294a(resource_name=GPIB_ADDRESS):
    """Open VISA session to Agilent 4294A."""
    rm = pyvisa.ResourceManager()
    inst = rm.open_resource(resource_name)
    inst.timeout = 120000  # ms
    inst.read_termination = '\n'
    inst.write_termination = '\n'
    return inst


def initialize_4294a_for_cpd(inst):
    """
    Basic initialization for Cp–D measurement:
    - Reset & clear
    - Internal trigger
    - Cp–D mode
    """
    inst.write("*RST")
    inst.write("*CLS")

    # Internal trigger
    inst.write("TRGS INT")

    # Cp–D measurement: Cp-D (trace A = Cp, trace B = D)
    # This matches your earlier code that used "MEAS CPD"
    inst.write("MEAS CPD")

    # ASCII data transfer
    inst.write("FORM4")

    # Oscillator: voltage mode and level (e.g. 100 mV RMS)
    inst.write("POWMOD VOLT")
    inst.write("POWE 0.5")  # 0.1 V rms AC level

    # DC bias configuration (voltage mode)
    inst.write("DCMOD VOLT")
    inst.write("DCRNG M10")  # 10 mA bias range (example)


def configure_dc_bias(inst, apply_bias, dc_bias_v):
    """
    Turn DC bias ON/OFF and set level if needed.
    """
    if apply_bias:
        # Set DC bias level and enable
        inst.write(f"DCV {dc_bias_v}")
        inst.write("DCO ON")
    else:
        # Turn off DC bias
        inst.write("DCO OFF")


def set_frequency_sweep(inst, f_start, f_stop, n_points):
    """
    Configure a LOG frequency sweep (no DC bias sweep).
    Sweep parameter = FREQ.
    """
    inst.write("SWPP FREQ")   # sweep parameter: frequency
    inst.write("SWPT LOG")    # log sweep in frequency
    inst.write(f"POIN {n_points}")
    inst.write(f"STAR {f_start}")
    inst.write(f"STOP {f_stop}")
    inst.write("SWED UP")     # sweep direction up


def single_sweep_and_wait(inst):
    """
    Start a single sweep and wait for completion via *OPC?.
    """
    inst.write("SING")
    _ = inst.query("*OPC?")  # returns '1' when sweep completed


def read_sweep_axis(inst):
    """
    Read sweep parameter values for all points (here: frequency).
    Uses OUTPSWPRM?.
    """
    axis_str = inst.query("OUTPSWPRM?")
    axis = np.array([float(x) for x in axis_str.split(',') if x.strip() != ""])
    return axis


def read_trace_main(inst, trace="A"):
    """
    Read main data from a given trace (A or B).
    OUTPDTRC? returns [val, sub, val, sub, ...].
    We take every second element starting at index 0.
    """
    inst.write(f"TRAC {trace}")
    data_str = inst.query("OUTPDTRC?")
    data = np.array([float(x) for x in data_str.split(',') if x.strip() != ""])
    main_vals = data[0::2]  # main reading (skip sub-reading)
    return main_vals


def measure_cpd_vs_freq(inst):
    """
    Perform one LOG frequency sweep and measure:
        Cp(f) and D(f),
    and return frequency array and both traces.
    """
    # Configure sweep
    set_frequency_sweep(inst, FREQ_START_HZ, FREQ_STOP_HZ, NUM_POINTS)

    # Start sweep and wait until done
    single_sweep_and_wait(inst)

    # Read frequency axis
    freq_axis = read_sweep_axis(inst)

    # Read Cp from trace A and D from trace B
    cp_vals = read_trace_main(inst, trace="A")
    d_vals = read_trace_main(inst, trace="B")

    return freq_axis, cp_vals, d_vals


# ==============================
# Main
# ==============================

def main():
    ensure_output_dir(output_dir)

    inst = connect_4294a()
    try:
        initialize_4294a_for_cpd(inst)
        configure_dc_bias(inst, APPLY_DC_BIAS, DC_BIAS_V)

        print("Running Cp–D vs frequency measurement (LOG sweep)...")
        freq_axis, cp_vals, d_vals = measure_cpd_vs_freq(inst)

        # ==========================
        # PLOTS (saved as PNG; no GUI blocking)
        # ==========================

        # Cp vs frequency (log x-axis)
        fig1 = plt.figure()
        plt.semilogx(freq_axis, cp_vals, '-o', markersize=3)
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Cp (F)")
        title1 = "Cp vs Frequency"
        if APPLY_DC_BIAS:
            title1 += f" (DC Bias = {DC_BIAS_V:.2f} V)"
        plt.title(title1)
        plt.grid(True, which="both")

        fig1_path = uniquify(os.path.join(output_dir, "cp_vs_frequency.png"))
        fig1.savefig(fig1_path, dpi=300, bbox_inches="tight")
        plt.close(fig1)

        # D vs frequency (log x-axis)
        fig2 = plt.figure()
        plt.semilogx(freq_axis, d_vals, '-o', markersize=3)
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("D (loss factor)")
        title2 = "D vs Frequency"
        if APPLY_DC_BIAS:
            title2 += f" (DC Bias = {DC_BIAS_V:.2f} V)"
        plt.title(title2)
        plt.grid(True, which="both")

        fig2_path = uniquify(os.path.join(output_dir, "d_vs_frequency.png"))
        fig2.savefig(fig2_path, dpi=300, bbox_inches="tight")
        plt.close(fig2)

        print(f"Saved plots to:\n  {fig1_path}\n  {fig2_path}")

        # ==========================
        # SAVE DATA TO CSV
        # ==========================

        csv_path = uniquify(os.path.join(output_dir, "cpd_vs_frequency.csv"))
        data = np.column_stack([freq_axis, cp_vals, d_vals])
        header = "frequency_Hz, Cp_F, D"

        np.savetxt(csv_path, data, delimiter=",", header=header, comments="")
        print(f"Saved Cp–D data to: {csv_path}")

    finally:
        # Turn off DC bias and close session safely
        try:
            inst.write("DCO OFF")
        except Exception:
            pass

        inst.close()
        print("Instrument session closed. Script finished.")


if __name__ == "__main__":
    main()
