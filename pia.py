"""
Filename: pia.py
Author: Ethan Ruddell
Date: 2025-01-23
Description: Contains all constants and functions for PIA measurements.
"""

import pyvisa
import numpy as np

# TODO: check actual needed inputs for each measurement, eg C-V does not need frequency

# =============================
# User settings and constants
# =============================
GPIB_ADDRESS = "GPIB0::24::INSTR"   # GPIB address for PIA

# Default Frequency sweep settings (LOG sweep)
FREQ_START_HZ = 4e1      # start frequency, Hz (>0 for LOG sweep)
FREQ_STOP_HZ = 5e6       # stop frequency, Hz
NUM_POINTS = 201         # number of points in LOG sweep

# Default DC bias options
APPLY_DC_BIAS = True     # True = apply a fixed DC bias; False = no DC bias
DC_BIAS_V = 0         # bias voltage (V) if APPLY_DC_BIAS = True

# =============================
# Measurement Parameters
# =============================
# This dictionary stores the current measurement parameters
# allowing for duplication across measurement types
current_parameters = {
    "freq_start_hz": FREQ_START_HZ,
    "freq_stop_hz": FREQ_STOP_HZ,
    "num_points": NUM_POINTS,
    "apply_dc_bias": APPLY_DC_BIAS,
    "dc_bias_v": DC_BIAS_V,
}

# =============================
# Relevant VARIABLES for access
# =============================
freq_axis = "" # Placeholder for measurement data storage
cap_vals = "" # Placeholder for measurement data storage
d_vals = "" # Placeholder for measurement data storage

# =============================
# Setup Functions
# =============================
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
        if( dc_bias_v < -40.0 or dc_bias_v > 40.0 ):
            print("Warning: DC Bias voltage out of range (-40V to 40V). Setting to 0V.")
            dc_bias_v = 0.0
            return
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


def measure_cpd_vs_freq(inst, freq_start=None, freq_stop=None, num_points=None):
    """
    Perform one LOG frequency sweep and measure:
        Cp(f) and D(f),
    and return frequency array and both traces.
    
    Args:
        inst: VISA instrument instance
        freq_start: Start frequency in Hz (uses FREQ_START_HZ if None)
        freq_stop: Stop frequency in Hz (uses FREQ_STOP_HZ if None)
        num_points: Number of points (uses NUM_POINTS if None)
    """
    # Use provided parameters or fall back to module constants
    if freq_start is None:
        freq_start = FREQ_START_HZ
    if freq_stop is None:
        freq_stop = FREQ_STOP_HZ
    if num_points is None:
        num_points = NUM_POINTS
    
    # Configure sweep
    set_frequency_sweep(inst, freq_start, freq_stop, num_points)

    # Start sweep and wait until done
    single_sweep_and_wait(inst)

    # Read frequency axis
    freq_axis = read_sweep_axis(inst)

    # Read Cp from trace A and D from trace B
    cp_vals = read_trace_main(inst, trace="A")
    d_vals = read_trace_main(inst, trace="B")

    return freq_axis, cp_vals, d_vals


def measure_impedance_vs_freq(inst, freq_start=None, freq_stop=None, num_points=None):
    """
    Perform one LOG frequency sweep and measure:
        |Z|(f) and θ(f) (degrees),
    and return frequency array and both traces.
    
    Args:
        inst: VISA instrument instance
        freq_start: Start frequency in Hz (uses FREQ_START_HZ if None)
        freq_stop: Stop frequency in Hz (uses FREQ_STOP_HZ if None)
        num_points: Number of points (uses NUM_POINTS if None)
    """
    # Use provided parameters or fall back to module constants
    if freq_start is None:
        freq_start = FREQ_START_HZ
    if freq_stop is None:
        freq_stop = FREQ_STOP_HZ
    if num_points is None:
        num_points = NUM_POINTS
    
    # Configure sweep
    set_frequency_sweep(inst, freq_start, freq_stop, num_points)

    # Start sweep and wait until done
    single_sweep_and_wait(inst)

    # Read frequency axis
    freq_axis = read_sweep_axis(inst)

    # Read |Z| from trace A and θ from trace B
    z_mag = read_trace_main(inst, trace="A")
    theta_deg = read_trace_main(inst, trace="B")

    return freq_axis, z_mag, theta_deg

def setup():
    """
    Connect to and initialize the 4294A for measurements.
    Returns the instrument handle.
    """
    inst = connect_4294a()
    return inst


# =============================
# Parameter Configuration and Duplication
# =============================

def get_frequency_parameters(use_last=None, freq_start=None, freq_stop=None, num_points=None):
    """
    Get and optionally update frequency sweep parameters.
    
    Args:
        use_last: If True, use previously stored parameters. If None, prompt user.
        freq_start: Start frequency in Hz (optional)
        freq_stop: Stop frequency in Hz (optional)
        num_points: Number of points (optional)
    
    Returns:
        tuple: (freq_start, freq_stop, num_points)
    """
    global current_parameters
    
    if use_last:
        freq_start = current_parameters["freq_start_hz"]
        freq_stop = current_parameters["freq_stop_hz"]
        num_points = current_parameters["num_points"]
    else:
        # Prompt user for parameters with defaults from current_parameters
        freq_start_input = input(f"Enter start frequency (Hz) [default {current_parameters['freq_start_hz']:.2e}]: ").strip()
        freq_start = float(freq_start_input) if freq_start_input else current_parameters["freq_start_hz"]
        
        freq_stop_input = input(f"Enter stop frequency (Hz) [default {current_parameters['freq_stop_hz']:.2e}]: ").strip()
        freq_stop = float(freq_stop_input) if freq_stop_input else current_parameters["freq_stop_hz"]
        
        num_points_input = input(f"Enter number of points [default {current_parameters['num_points']}]: ").strip()
        num_points = int(num_points_input) if num_points_input else current_parameters["num_points"]
        
        # Update current parameters for future duplication
        current_parameters["freq_start_hz"] = freq_start
        current_parameters["freq_stop_hz"] = freq_stop
        current_parameters["num_points"] = num_points
    
    return freq_start, freq_stop, num_points


def get_dc_bias_parameters(use_last=None, apply_bias=None, bias_voltage=None):
    """
    Get and optionally update DC bias parameters.
    
    Args:
        use_last: If True, use previously stored parameters. If None, prompt user.
        apply_bias: Boolean to apply DC bias (optional)
        bias_voltage: DC bias voltage in V (optional)
    
    Returns:
        tuple: (apply_bias, bias_voltage)
    """
    global current_parameters
    
    if use_last:
        apply_bias = current_parameters["apply_dc_bias"]
        bias_voltage = current_parameters["dc_bias_v"]
    else:
        # Prompt user for DC bias parameters
        bias_input = input(f"Apply DC bias? [y/n] [default {'y' if current_parameters['apply_dc_bias'] else 'n'}]: ").strip().lower()
        apply_bias = bias_input == 'y' if bias_input else current_parameters["apply_dc_bias"]
        
        if apply_bias:
            bias_v_input = input(f"Enter DC bias voltage (V) [default {current_parameters['dc_bias_v']}]: ").strip()
            bias_voltage = float(bias_v_input) if bias_v_input else current_parameters["dc_bias_v"]
        else:
            bias_voltage = 0.0
        
        # Update current parameters for future duplication
        current_parameters["apply_dc_bias"] = apply_bias
        current_parameters["dc_bias_v"] = bias_voltage
    
    return apply_bias, bias_voltage


def prompt_for_parameter_duplication():
    """
    Ask user if they want to duplicate the last measurement parameters.
    
    Returns:
        bool: True if user wants to duplicate, False otherwise
    """
    dup_input = input("Use last measurement parameters? (y/n): ").strip().lower()
    return dup_input == 'y'


def set_frequency_sweep_with_params(inst, freq_start, freq_stop, num_points):
    """
    Configure a LOG frequency sweep using provided parameters.
    Wrapper around set_frequency_sweep with explicit parameters.
    """
    set_frequency_sweep(inst, freq_start, freq_stop, num_points)


# =============================
# Dielectric C–V Measurement Helpers
# =============================

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


def set_frequency(inst, freq_hz):
    """Set CW frequency used during DC bias sweep."""
    inst.write(f"CWFREQ {freq_hz}")


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


def measure_single_cv_cycle(inst, freq_hz, v_min, v_max, n_points):
    """
    Measure ONE C–V butterfly loop using two sweeps:
    1) v_min -> v_max (UP)
    2) v_max -> v_min (DOWN)
    Returns concatenated bias and Cp arrays, plus εr.
    Sweep order is preserved: v_min → v_max → v_min.
    """
    set_frequency(inst, freq_hz)

    # Sweep UP: v_min → v_max
    set_dc_bias_sweep(inst, v_min, v_max, n_points, direction="UP")
    single_sweep_and_wait(inst)
    inst.write("TRAC A")
    inst.write("AUTO")
    v_up = read_sweep_axis(inst)
    cp_up = read_trace_cp(inst)

    # Sweep DOWN: v_max → v_min
    set_dc_bias_sweep(inst, v_min, v_max, n_points, direction="DOWN")
    inst.write("TRAC A")
    inst.write("AUTO")
    v_down = read_sweep_axis(inst)
    cp_down = read_trace_cp(inst)

    # Combine into a single "butterfly" trace in actual sweep order
    v_full = np.concatenate([v_up, v_down])
    cp_full = np.concatenate([cp_up, cp_down])

    return v_full, cp_full


def measure_eps_vs_freq(inst, freq_start, freq_stop, freq_points, bias_voltage=0.0):
    """
    Optional: measure εr vs frequency at a fixed DC bias.
    Returns frequency array, Cp(f), and εr(f).
    """
    inst.write(f"DCV {bias_voltage}")
    inst.write("DCO ON")

    set_frequency_sweep(inst, freq_start, freq_stop, freq_points)
    single_sweep_and_wait(inst)

    inst.write("TRAC A")
    inst.write("AUTO")

    freq_axis = read_sweep_axis(inst)
    cp_freq = read_trace_cp(inst)

    return freq_axis, cp_freq