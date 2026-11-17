import os
import pyvisa
import numpy as np
import matplotlib

# Use non-interactive backend so script finishes without waiting for plot windows
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ==============================
# User settings (EDIT THESE)
# ==============================

GPIB_ADDRESS = "GPIB0::24::INSTR"  # Example: GPIB0, address 24

# Frequency sweep settings (LOG sweep)
FREQ_START_HZ = 1e3     # start frequency, Hz
FREQ_STOP_HZ = 1e7      # stop frequency, Hz
NUM_POINTS = 201        # number of points in LOG sweep

# DC bias options
APPLY_DC_BIAS = True    # True = apply a fixed DC bias; False = no DC bias
DC_BIAS_V = 0.0         # bias voltage (V) if APPLY_DC_BIAS = True

# Output folder
OUTPUT_DIR = r"C:\Users\NRG\PycharmProjects\Agilent 4294A Impedance Analyzer\output_impedance_theta"


# ==============================
# Helper functions
# ==============================

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


def initialize_4294a_for_impedance(inst):
    """
    Basic initialization for impedance |Z| – θ measurement:
    - Reset & clear
    - Internal trigger
    - Z-θ mode
    """
    inst.write("*RST")
    inst.write("*CLS")

    # Internal trigger
    inst.write("TRGS INT")

    # Impedance measurement: Z-Theta (trace A = |Z|, trace B = θ)
    # (SCPI name may vary by firmware; ZTD is common on 4294A)
    inst.write("MEAS ZTD")

    # ASCII data transfer
    inst.write("FORM4")

    # Oscillator: voltage mode and level (e.g. 100 mV RMS)
    inst.write("POWMOD VOLT")
    inst.write("POWE 0.1")  # 0.1 V rms AC level

    # DC bias configuration (voltage mode)
    inst.write("DCMOD VOLT")
    inst.write("DCRNG M10")  # 10 mA bias range (example)


def configure_dc_bias(inst, apply_bias, dc_bias_v):
    """
    Turn DC bias ON/OFF and set level if needed.
    """
    if apply_bias:
        # Set DC bias level and enable
        # (Command DCV is typical for 4294A; adjust if your firmware uses another name)
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


def measure_impedance_vs_freq(inst):
    """
    Perform one LOG frequency sweep and measure:
        |Z|(f) and θ(f) (degrees),
    and return frequency array and both traces.
    """
    # Configure sweep
    set_frequency_sweep(inst, FREQ_START_HZ, FREQ_STOP_HZ, NUM_POINTS)

    # Start sweep and wait until done
    single_sweep_and_wait(inst)

    # Read frequency axis
    freq_axis = read_sweep_axis(inst)

    # Read |Z| from trace A and θ from trace B
    z_mag = read_trace_main(inst, trace="A")
    theta_deg = read_trace_main(inst, trace="B")

    return freq_axis, z_mag, theta_deg


# ==============================
# Main
# ==============================

def main():
    ensure_output_dir(OUTPUT_DIR)

    inst = connect_4294a()
    try:
        initialize_4294a_for_impedance(inst)
        configure_dc_bias(inst, APPLY_DC_BIAS, DC_BIAS_V)

        print("Running impedance vs frequency measurement (LOG sweep)...")
        freq_axis, z_mag, theta_deg = measure_impedance_vs_freq(inst)

        # Compute real and imaginary parts from |Z| and θ
        theta_rad = np.deg2rad(theta_deg)
        z_real = z_mag * np.cos(theta_rad)
        z_imag = z_mag * np.sin(theta_rad)

        # ==========================
        # PLOTS (saved as PNG; no GUI blocking)
        # ==========================

        # |Z| vs frequency (log x-axis)
        fig1 = plt.figure()
        plt.semilogx(freq_axis, z_mag, '-o', markersize=3)
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("|Z| (Ohm)")
        title1 = f"Impedance Magnitude vs Frequency"
        if APPLY_DC_BIAS:
            title1 += f" (DC Bias = {DC_BIAS_V:.2f} V)"
        plt.title(title1)
        plt.grid(True, which="both")

        fig1_path = os.path.join(OUTPUT_DIR, "impedance_magnitude_vs_frequency.png")
        fig1.savefig(fig1_path, dpi=300, bbox_inches="tight")
        plt.close(fig1)

        # θ vs frequency (log x-axis)
        fig2 = plt.figure()
        plt.semilogx(freq_axis, theta_deg, '-o', markersize=3)
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Phase (deg)")
        title2 = f"Impedance Phase vs Frequency"
        if APPLY_DC_BIAS:
            title2 += f" (DC Bias = {DC_BIAS_V:.2f} V)"
        plt.title(title2)
        plt.grid(True, which="both")

        fig2_path = os.path.join(OUTPUT_DIR, "impedance_phase_vs_frequency.png")
        fig2.savefig(fig2_path, dpi=300, bbox_inches="tight")
        plt.close(fig2)

        print(f"Saved plots to:\n  {fig1_path}\n  {fig2_path}")

        # ==========================
        # SAVE DATA TO CSV
        # ==========================

        csv_path = os.path.join(OUTPUT_DIR, "impedance_vs_frequency.csv")
        data = np.column_stack([freq_axis, z_mag, theta_deg, z_real, z_imag])
        header = "frequency_Hz, Z_mag_ohm, theta_deg, Z_real_ohm, Z_imag_ohm"

        np.savetxt(csv_path, data, delimiter=",", header=header, comments="")
        print(f"Saved impedance data to: {csv_path}")

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
