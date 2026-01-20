"""
Filename: Dielectric Measurements.py
Author: Chaitanya Sharma, Ethan Ruddell
Date: 2025-01-16
Description: Contains commands for dielectric measurements.
"""

import os
import pyvisa
import numpy as np
import matplotlib
import sys

matplotlib.use("TkAgg")  # use Agg if TkAgg fails and you only want files
import matplotlib.pyplot as plt

# ==============================
# User settings (EDIT THESE)
# ==============================

GPIB_ADDRESS = "GPIB0::24::INSTR"  # GPIB0, address 24

# C–V / dielectric butterfly settings
freq_for_cv = 25000  # Hz, measurement frequency for C-V
v_min = 0 # V, start DC bias for C-V
v_max = 5  # V, stop DC bias for C-V
num_points = 401  # number of points per sweep (up and down)
num_cycles = 1  # number of full cycles (-2.5 -> 2.5 -> -2.5)

# HZO capacitor geometry (edit to your device)
t_hzo_nm = 10.0  # thickness in nm
diam_um = 75.0  # electrode diameter in µm

# Optional: frequency sweep example (sweep freq at fixed bias)
freq_start = 0  # 1 kHz
freq_stop = 0  # 1 MHz
freq_points = 0
bias_for_f_sweep = 0  # V

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


# checks to ensure output directory exists, if not, make it
def ensure_output_dir(path):
    os.makedirs(path, exist_ok=True)


def connect_4294a(resource_name=GPIB_ADDRESS):
    """Open VISA session to Agilent 4294A."""
    rm = pyvisa.ResourceManager()
    inst = rm.open_resource(resource_name)
    inst.timeout = 120000 # ms, increased from original 10000 to acct for lower frequency
    inst.read_termination = '\n'
    inst.write_termination = '\n'
    return inst


def initialize_4294a(inst):
    """Basic initialization: reset, clear, internal trigger, Cp-D mode."""
    inst.write("*RST")
    inst.write("*CLS")

    # Internal trigger for sweeps
    inst.write("TRGS INT")

    # Cp-D measurement (trace A = Cp, trace B = D) — from Agilent sample program
    inst.write("MEAS CPD")

    # Use ASCII transfer (FORM4) – default, but set explicitly
    inst.write("FORM4")

    # Set oscillator mode and level (voltage mode, ~100 mV as example)
    inst.write("POWMOD VOLT")
    inst.write("POWE 0.1")  # 0.1 V rms AC level (adjust as needed)

    # DC bias: voltage mode, turn ON (limit/current range etc. can be added if needed)
    inst.write("DCMOD VOLT")
    inst.write("DCRNG M10")  # 10 mA bias range (example; adjust if needed)
    inst.write("DCO ON")


def set_frequency(inst, freq_hz):
    """Set CW frequency used during DC bias sweep."""
    inst.write(f"CWFREQ {freq_hz}")


def set_dc_bias_sweep(inst, v_start, v_stop, n_points, direction="UP"):
    """
    Configure a DC bias sweep (sweep parameter = DC bias level).
    direction: 'UP' or 'DOWN'
    """
    inst.write("SWPP DCB")  # sweep parameter = DC bias
    inst.write("SWPT LIN")  # linear sweep
    inst.write(f"POIN {n_points}")
    inst.write(f"STAR {v_start}")
    inst.write(f"STOP {v_stop}")
    inst.write(f"SWED {direction}")


def set_frequency_sweep(inst, f_start, f_stop, n_points):
    """Configure a frequency sweep at fixed DC bias (for εr–f measurement)."""
    inst.write("SWPP FREQ")
    inst.write("SWPT LOG")  # LOG sweep typical for frequency
    inst.write(f"POIN {n_points}")
    inst.write(f"STAR {f_start}")
    inst.write(f"STOP {f_stop}")
    inst.write("SWED UP")


def single_sweep_and_wait(inst):
    """Start a single sweep and wait for completion via *OPC?."""
    inst.write("SING")
    _ = inst.query("*OPC?")  # returns '1' when sweep completed


def read_trace_cp(inst):
    """
    Read Cp trace from active trace.
    OUTPDTRC? returns [val, sub, val, sub, ...].
    For Cp-D, Cp is trace A, D is trace B. We read Cp on active trace.
    """
    data_str = inst.query("OUTPDTRC?")
    data = np.array([float(x) for x in data_str.split(',') if x.strip() != ""])
    cp = data[0::2]  # take every 2nd element starting at index 0 (main reading)
    return cp


def read_sweep_axis(inst):
    """
    Read sweep parameter values for all points (frequency or DC bias, depending on SWPP).
    Uses OUTPSWPRM?.
    """
    axis_str = inst.query("OUTPSWPRM?")
    axis = np.array([float(x) for x in axis_str.split(',') if x.strip() != ""])
    return axis


def compute_eps_r(cp_array, thickness_nm, diameter_um):
    """
    Compute dielectric constant εr from Cp:
        εr = Cp * t / (ε0 * A)
    where:
        t = thickness (m)
        A = electrode area (m^2) = π (d/2)^2
    Cp is in F.
    """
    eps0 = 8.854e-12  # F/m
    t_m = thickness_nm * 1e-9
    d_m = diameter_um * 1e-6
    area = np.pi * (d_m / 2.0) ** 2
    eps_r = cp_array * t_m / (eps0 * area)
    return eps_r


# ==============================
# Main measurement flow
# ==============================

# TODO: add more measurements, error handling

def measure_single_cv_cycle(inst):
    """
    Measure ONE C–V butterfly loop using two sweeps:
    1) v_min -> v_max (UP)
    2) v_max -> v_min (DOWN)
    Returns concatenated bias and Cp arrays, plus εr.
    Sweep order is preserved: -2.5 → +2.5 → -2.5.
    """
    set_frequency(inst, freq_for_cv)

    # Sweep UP: v_min → v_max
    set_dc_bias_sweep(inst, v_min, v_max, num_points, direction="UP")
    single_sweep_and_wait(inst)
    inst.write("TRAC A")
    inst.write("AUTO")
    v_up = read_sweep_axis(inst)
    cp_up = read_trace_cp(inst)

    # Sweep DOWN: v_max → v_min
    set_dc_bias_sweep(inst, v_min, v_max, num_points, direction="DOWN")
    inst.write("TRAC A")
    inst.write("AUTO")
    v_down = read_sweep_axis(inst)
    cp_down = read_trace_cp(inst)

    # Combine into a single “butterfly” trace in actual sweep order
    v_full = np.concatenate([v_up, v_down])
    cp_full = np.concatenate([cp_up, cp_down])
    eps_r_full = compute_eps_r(cp_full, t_hzo_nm, diam_um)

    return v_full, cp_full, eps_r_full


def measure_eps_vs_freq(inst):
    """
    Optional: measure εr vs frequency at a fixed DC bias.
    Returns frequency array, Cp(f), and εr(f).
    """
    # inst.write(f"DCV {bias_for_f_sweep}")
    inst.write("DCO ON")

    set_frequency_sweep(inst, freq_start, freq_stop, freq_points)
    single_sweep_and_wait(inst)

    inst.write("TRAC A")
    inst.write("AUTO")

    freq_axis = read_sweep_axis(inst)
    cp_freq = read_trace_cp(inst)
    eps_r_freq = compute_eps_r(cp_freq, t_hzo_nm, diam_um)

    return freq_axis, cp_freq, eps_r_freq


def main():
    ensure_output_dir(output_dir) # run check to ensure output directory exists

    inst = connect_4294a() # connect to instrument
    try:
        initialize_4294a(inst) # initialize instrument

        # ==========================
        # 1) C–V & dielectric butterfly (multiple cycles)
        # ==========================
        all_v_cycles = []
        all_cp_cycles = []
        all_eps_cycles = []

        all_cycle_idx = []

        for cycle in range(num_cycles):
            print(f"Running C–V cycle {cycle + 1}/{num_cycles}...")
            v, cp, eps_r = measure_single_cv_cycle(inst)

            all_v_cycles.append(v)
            all_cp_cycles.append(cp)
            all_eps_cycles.append(eps_r)
            all_cycle_idx.append(np.full_like(v, cycle + 1, dtype=int))

        # Flatten for saving
        v_all = np.concatenate(all_v_cycles)
        cp_all = np.concatenate(all_cp_cycles)
        eps_all = np.concatenate(all_eps_cycles)
        cycles_all = np.concatenate(all_cycle_idx)

        # --- Plot C-V with all cycles overlaid ---
        fig1 = plt.figure()
        for i, (v, cp) in enumerate(zip(all_v_cycles, all_cp_cycles)):
            mid = len(v) // 2
            v_plot = np.insert(v, mid, np.nan)
            cp_plot = np.insert(cp, mid, np.nan)
            plt.plot(v_plot, cp_plot, label=f"Cycle {i + 1}")
        plt.xlabel("Voltage (V)")
        plt.ylabel("Capacitance Cp (F)")
        plt.title(f"C–V, {num_cycles} cycles @ {freq_for_cv / 1e3:.1f} kHz")
        plt.grid(True)
        plt.legend()

        # --- Plot εr vs V (dielectric butterfly, all cycles overlaid) ---
        fig2 = plt.figure()
        for i, (v, eps_r) in enumerate(zip(all_v_cycles, all_eps_cycles)):
            mid = len(v) // 2  # = num_points
            v_plot = np.insert(v, mid, np.nan)
            eps_plot = np.insert(eps_r, mid, np.nan)
            plt.plot(v_plot, eps_plot, label=f"Cycle {i + 1}")
        plt.xlabel("Voltage (V)")
        plt.ylabel("Dielectric constant εr")
        plt.title(f"Dielectric Butterfly, {num_cycles} cycles @ {freq_for_cv / 1e3:.1f} kHz")
        plt.grid(True)
        plt.legend()

        # ==========================
        # 2) εr vs frequency at fixed bias
        # ==========================
        print("Running εr vs frequency sweep...")
        freq_axis, cp_f, eps_f = measure_eps_vs_freq(inst)

        fig3 = plt.figure()
        plt.semilogx(freq_axis, eps_f, '.-')
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Dielectric constant εr")
        plt.title(f"εr vs Frequency at DC Bias = {bias_for_f_sweep:.2f} V")
        plt.grid(True, which="both")

        # ==========================
        # SAVE DATA AND FIGURES
        # ==========================

        # TODO: iterate filenames to prevent overwriting
        
        # C–V + εr(V) for all cycles
        cv_filename = os.path.join(output_dir, "cv_dielectric_cycles.csv")
        cv_data = np.column_stack([cycles_all, v_all, cp_all, eps_all])
        np.savetxt(
            cv_filename,
            cv_data,
            delimiter=",",
            header="cycle_index, bias_V, Cp_F, eps_r",
            comments=""
        )
        print(f"Saved C–V & εr(V) data to: {cv_filename}")

        # εr(f) data
        freq_filename = os.path.join(output_dir, "eps_vs_freq.csv")
        freq_data = np.column_stack([freq_axis, cp_f, eps_f])
        np.savetxt(
            freq_filename,
            freq_data,
            delimiter=",",
            header="frequency_Hz, Cp_F, eps_r",
            comments=""
        )
        print(f"Saved εr(f) data to: {freq_filename}")

        # PNGs for all plots
        fig1.savefig(os.path.join(output_dir, "cv_curve.png"), dpi=300, bbox_inches="tight")
        fig2.savefig(os.path.join(output_dir, "dielectric_butterfly.png"), dpi=300, bbox_inches="tight")
        fig3.savefig(os.path.join(output_dir, "eps_vs_freq.png"), dpi=300, bbox_inches="tight")
        print(f"Saved plot images to: {output_dir}")

    finally:
        try:
            inst.write("DCO OFF")
        except Exception:
            pass
        inst.close()
        print("Instrument session closed.")

        # Show plots for first visual
        plt.show()  

if __name__ == "__main__":
    main()
