"""
Filename: pia.py
Author: Ethan Ruddell
Date: 2025-01-23
Description: Driver for the Agilent 4294A Precision Impedance Analyzer (PIA).

The 4294A measures how a device (capacitor, dielectric film, etc.) responds to
an AC signal across a range of frequencies or DC bias voltages.  It can report
many different quantities:
  - Cp and D  (parallel capacitance and dissipation factor)
  - |Z| and θ (impedance magnitude and phase angle)
  - R and X   (resistance and reactance)
  - G and B   (conductance and susceptance)
  - |Y| and θ (admittance magnitude and phase angle)
  - εr        (relative permittivity / dielectric constant)

Communication uses plain GPIB commands (not SCPI).  The most important ones:
    - MEAS CPD / ZTD / RX / GB / YTD  – select measurement mode
    - SWPP FREQ / DCB                 – sweep frequency or DC bias
    - SING                             – start a single sweep
    - OUTPSWPRM? / OUTPDTRC?          – read back the data
    - POWMOD VOLT / POWE               – set the AC oscillator level

Typical workflow:
  1. connect_4294a()              – open a GPIB session
  2. initialize_4294a_for_*()     – reset, pick mode, set oscillator level
  3. configure_dc_bias()          – optionally apply a DC bias
  4. measure_*_vs_freq()          – run the sweep and read data
  5. close the session (handled by InstrumentSession in the caller)
"""

import pyvisa
import numpy as np
import logging
import config
from dielectric_utils import compute_eps_r_from_area
from gpib_utils import prompt_bool, safe_float_input, safe_int_input

log = logging.getLogger(__name__)

# =============================
# User settings and constants
# =============================
# These defaults are used when a measurement function is called without
# explicit parameters.  They can be overridden at runtime.
DEFAULT_GPIB_ADDRESS = "GPIB0::24::INSTR"
GPIB_ADDRESS = str(config.get("gpib", "pia", DEFAULT_GPIB_ADDRESS))

# Default Frequency sweep settings (LOG sweep)
FREQ_START_HZ = 4e1      # start frequency, Hz (must be >0 for LOG sweep)
FREQ_STOP_HZ = 5e6       # stop frequency, Hz (5 MHz, near 4294A max)
NUM_POINTS = 201         # number of data points across the sweep

# Default DC bias options
# DC bias applies a constant voltage to the DUT during the AC measurement,
# which is needed for C-V (capacitance vs voltage) measurements.
APPLY_DC_BIAS = True     # True = apply a fixed DC bias; False = no DC bias
DC_BIAS_V = 0            # bias voltage (V) if APPLY_DC_BIAS = True

# Default AC oscillator level for most PIA measurements.
DEFAULT_OSC_VOLTAGE_VRMS = float(config.get("pia", "osc_voltage_v", 0.5))
DEFAULT_IMPEDANCE_OSC_VRMS = 0.1

# =============================
# Measurement Parameters
# =============================
# This dictionary stores the parameters from the most recent measurement
# so the CLI can offer a "use last parameters?" shortcut.
current_parameters = {
    "freq_start_hz": FREQ_START_HZ,
    "freq_stop_hz": FREQ_STOP_HZ,
    "num_points": NUM_POINTS,
    "apply_dc_bias": APPLY_DC_BIAS,
    "dc_bias_v": DC_BIAS_V,
}


# =============================
# Setup Functions
# =============================
# Functions to connect to the instrument and configure it for a specific
# measurement mode.  The 4294A is configured by sending plain text commands
# over the GPIB bus.
def connect_4294a(resource_name=None):
    """Open VISA session to Agilent 4294A."""
    target = str(resource_name or config.get("gpib", "pia", GPIB_ADDRESS)).strip()
    if not target:
        raise RuntimeError("PIA GPIB address is empty. Set gpib.pia in config.yaml.")

    rm = pyvisa.ResourceManager()
    try:
        inst = rm.open_resource(target)
    except Exception as exc:
        try:
            resources = rm.list_resources()
        except Exception:
            resources = ()
        raise RuntimeError(
            "Failed to open PIA VISA resource "
            f"'{target}'. Available resources: {resources}. "
            "Check instrument power/cable, NI-MAX visibility, and gpib.pia in config.yaml."
        ) from exc

    inst.timeout = 120000  # ms
    inst.read_termination = '\n'
    inst.write_termination = '\n'
    return inst

def _initialize_4294a(inst, meas_mode="CPD", osc_level=None):
    """
    Generic 4294A initialization:
    - Reset & clear
    - Internal trigger
    - Set measurement mode
    - ASCII data transfer
    - Set oscillator level
    - Set DC bias configuration

    Args:
        inst: VISA instrument instance
        meas_mode: Measurement type string for the MEAS command
                   "CPD" = Cp-D, "ZTD" = Z-θ, "RX" = R-X,
                   "GB" = G-B, "YTD" = Y-θ
        osc_level: AC oscillator level in Vrms. If omitted, the default PIA
            oscillator level from config.yaml is used.
    """
    inst.write("*RST")
    inst.write("*CLS")
    inst.write("TRGS INT")
    inst.write(f"MEAS {meas_mode}")
    inst.write("FORM4")
    inst.write("POWMOD VOLT")
    if osc_level is None:
        osc_level = DEFAULT_OSC_VOLTAGE_VRMS
    inst.write(f"POWE {osc_level}")
    inst.write("DCMOD VOLT")
    inst.write("DCRNG M10")


# Convenience wrappers — thin calls to _initialize_4294a
def initialize_4294a_for_cpd(inst, osc_level=None):
    """Initialize for Cp-D measurement (trace A = Cp, trace B = D).

    The oscillator level is the AC drive amplitude applied to the device under
    test. Higher levels usually produce a stronger signal; lower levels are
    gentler on sensitive devices.
    """
    _initialize_4294a(inst, "CPD", osc_level=osc_level)


def initialize_4294a_for_impedance(inst, osc_level=None):
    """Initialize for Z-θ measurement (trace A = |Z|, trace B = θ).

    Impedance measurements historically use a lower AC drive level, so the
    default remains conservative unless the caller overrides it.
    """
    if osc_level is None:
        osc_level = DEFAULT_IMPEDANCE_OSC_VRMS
    _initialize_4294a(inst, "ZTD", osc_level=osc_level)


def initialize_4294a_for_rx(inst, osc_level=None):
    """Initialize for R-X measurement (trace A = R, trace B = X)."""
    _initialize_4294a(inst, "RX", osc_level=osc_level)


def initialize_4294a_for_gb(inst, osc_level=None):
    """Initialize for G-B measurement (trace A = G, trace B = B)."""
    _initialize_4294a(inst, "GB", osc_level=osc_level)


def initialize_4294a_for_ytd(inst, osc_level=None):
    """Initialize for Y-θ measurement (trace A = |Y|, trace B = θ)."""
    _initialize_4294a(inst, "YTD", osc_level=osc_level)


def configure_dc_bias(inst, apply_bias, dc_bias_v):
    """
    Turn DC bias ON/OFF and set level if needed.
    """
    if apply_bias:
        # Set DC bias level and enable
        if( dc_bias_v < -40.0 or dc_bias_v > 40.0 ):
            print("Warning: DC Bias voltage out of range (-40V to 40V). Setting to 0V.")
            dc_bias_v = 0.0
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
# These helpers support the CLI's "use last measurement parameters?" feature.
# When the user says yes, the software re-uses the values from
# current_parameters instead of prompting for new ones.

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
        return (
            current_parameters["freq_start_hz"],
            current_parameters["freq_stop_hz"],
            current_parameters["num_points"],
        )

    freq_start = safe_float_input(
        f"Enter start frequency (Hz) [default {current_parameters['freq_start_hz']:.2e}]: ",
        current_parameters["freq_start_hz"],
    )
    freq_stop = safe_float_input(
        f"Enter stop frequency (Hz) [default {current_parameters['freq_stop_hz']:.2e}]: ",
        current_parameters["freq_stop_hz"],
    )
    num_points = safe_int_input(
        f"Enter number of points [default {current_parameters['num_points']}]: ",
        current_parameters["num_points"],
    )

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
        return current_parameters["apply_dc_bias"], current_parameters["dc_bias_v"]

    apply_bias = prompt_bool("Apply DC bias?", current_parameters["apply_dc_bias"])
    if apply_bias:
        bias_voltage = safe_float_input(
            f"Enter DC bias voltage (V) [default {current_parameters['dc_bias_v']}]: ",
            current_parameters["dc_bias_v"],
        )
    else:
        bias_voltage = 0.0

    current_parameters["apply_dc_bias"] = apply_bias
    current_parameters["dc_bias_v"] = bias_voltage
    
    return apply_bias, bias_voltage


def prompt_for_parameter_duplication():
    """
    Ask user if they want to duplicate the last measurement parameters.
    
    Returns:
        bool: True if user wants to duplicate, False otherwise
    """
    return prompt_bool("Use last measurement parameters?", False)



# =============================
# Dielectric C–V Measurement Helpers
# =============================
# C-V (Capacitance vs Voltage) measurements sweep a DC bias across the DUT
# at a fixed AC frequency and measure how capacitance changes.  A "butterfly"
# cycle sweeps up then down to reveal hysteresis (memory effects).
# These are essential for characterising ferroelectric thin films (e.g. HZO).

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




def compute_eps_r(cp_array, thickness_nm, area_um2):
    """
    Compute dielectric constant εr from Cp:
        εr = Cp * t / (ε0 * A)
    where:
        t = thickness (m)
        A = electrode area (m^2)

    Args:
        cp_array: Capacitance array in F
        thickness_nm: Dielectric thickness in nm
        area_um2: Electrode area in µm²
    """
    return compute_eps_r_from_area(cp_array, thickness_nm, area_um2)


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
    cp_up = read_trace_main(inst, trace="A")

    # Sweep DOWN: v_max → v_min
    set_dc_bias_sweep(inst, v_min, v_max, n_points, direction="DOWN")
    inst.write("TRAC A")
    inst.write("AUTO")
    v_down = read_sweep_axis(inst)
    cp_down = read_trace_main(inst, trace="A")

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
    cp_freq = read_trace_main(inst, trace="A")

    return freq_axis, cp_freq


# =============================
# Additional Measurement Modes
# =============================
# R-X, G-B, and Y-θ are alternative ways of expressing impedance data.
# They all use the same generic two-trace sweep under the hood.


def _measure_two_trace_vs_freq(inst, freq_start=None, freq_stop=None, num_points=None):
    """
    Generic LOG frequency sweep reading both traces A and B.

    Returns:
        tuple: (freq_axis, trace_a_vals, trace_b_vals)
    """
    if freq_start is None:
        freq_start = FREQ_START_HZ
    if freq_stop is None:
        freq_stop = FREQ_STOP_HZ
    if num_points is None:
        num_points = NUM_POINTS

    set_frequency_sweep(inst, freq_start, freq_stop, num_points)
    single_sweep_and_wait(inst)

    freq_axis = read_sweep_axis(inst)
    a_vals = read_trace_main(inst, trace="A")
    b_vals = read_trace_main(inst, trace="B")
    return freq_axis, a_vals, b_vals


def measure_rx_vs_freq(inst, freq_start=None, freq_stop=None, num_points=None):
    """
    Perform LOG frequency sweep and measure R(f) and X(f).

    Returns:
        tuple: (freq_axis, r_vals, x_vals)
    """
    return _measure_two_trace_vs_freq(inst, freq_start, freq_stop, num_points)


def measure_gb_vs_freq(inst, freq_start=None, freq_stop=None, num_points=None):
    """
    Perform LOG frequency sweep and measure G(f) and B(f).

    Returns:
        tuple: (freq_axis, g_vals, b_vals)
    """
    return _measure_two_trace_vs_freq(inst, freq_start, freq_stop, num_points)


def measure_ytd_vs_freq(inst, freq_start=None, freq_stop=None, num_points=None):
    """
    Perform LOG frequency sweep and measure |Y|(f) and θ(f).

    Returns:
        tuple: (freq_axis, y_mag, y_theta_deg)
    """
    return _measure_two_trace_vs_freq(inst, freq_start, freq_stop, num_points)