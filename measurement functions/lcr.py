"""
Filename: lcr.py
Author: Ethan Ruddell
Date: 2026-02-12
Description: Contains all constants and functions for LCR measurements.
            For use with E4980A LCR from Agilent (Keysight).
"""

import pyvisa
import numpy as np
import time
import logging

log = logging.getLogger(__name__)

# =============================
# User settings and constants
# =============================
GPIB_ADDRESS = "GPIB0::_____::INSTR"   # GPIB address for E4980A LCR

# Default Frequency settings
FREQ_DEFAULT_HZ = 1e3    # default frequency, Hz
FREQ_MIN_HZ = 20         # minimum frequency, Hz
FREQ_MAX_HZ = 2e6        # maximum frequency, Hz (E4980A: 20 Hz to 2 MHz)

# Default AC signal level
AC_LEVEL_V = 1.0         # AC signal level in V (E4980A: 0 to 20V)
AC_LEVEL_MIN_V = 0.0
AC_LEVEL_MAX_V = 20.0

# Default DC bias settings
APPLY_DC_BIAS = False    # True = apply DC bias; False = no DC bias
DC_BIAS_V = 0.0          # DC bias voltage (V) (E4980A: -40V to +40V)
DC_BIAS_MIN_V = -40.0
DC_BIAS_MAX_V = 40.0

# Frequency sweep settings
SWEEP_MODE = "LIN"       # "LIN" or "LOG"
SWEEP_POINTS = 201       # number of sweep points

# Measurement speed/integration time
MEAS_TIME = "MED"        # "SHOR", "MED", "LONG", "LONG2" (E4980A specific)

# =============================
# Measurement Parameters Storage
# =============================
current_parameters = {
    "frequency": FREQ_DEFAULT_HZ,
    "ac_level_v": AC_LEVEL_V,
    "apply_dc_bias": APPLY_DC_BIAS,
    "dc_bias_v": DC_BIAS_V,
    "sweep_points": SWEEP_POINTS,
    "meas_time": MEAS_TIME,
}

# =============================
# Connection and Initialization
# =============================

def connect_e4980a(resource_name=None):
    """
    Open VISA session to E4980A LCR meter.
    
    Args:
        resource_name: GPIB address string (uses GPIB_ADDRESS if None)
    
    Returns:
        VISA instrument instance
    """
    if resource_name is None:
        resource_name = GPIB_ADDRESS
    
    try:
        rm = pyvisa.ResourceManager()
        print(f"Attempting to connect to: {resource_name}")
        inst = rm.open_resource(resource_name)
        inst.timeout = 30000  # 30 second timeout
        inst.read_termination = '\n'
        inst.write_termination = '\n'
        
        # Verify connection
        idn = inst.query("*IDN?")
        print(f"Connected to: {idn.strip()}")
        
        return inst
    except Exception as e:
        print(f"Error connecting to E4980A: {e}")
        raise

def disconnect_e4980a(inst):
    """Safely disconnect from E4980A."""
    try:
        # Turn off DC bias
        inst.write(":BIAS:STAT OFF")
        # Return to manual trigger (safe idle state)
        inst.write(":TRIG:SOUR HOLD")
    except Exception:
        pass
    inst.close()
    print("E4980A disconnected")

def initialize_e4980a(inst, measurement_function="CPD"):
    """
    Basic initialization for E4980A:
    - Reset & clear
    - Set measurement function
    - Configure trigger and format
    
    Args:
        inst: VISA instrument instance
        measurement_function: Measurement type
            "CPD" = Cp-D (capacitance-dissipation)
            "CPRP" = Cp-Rp (capacitance-parallel resistance)
            "CPQ" = Cp-Q (capacitance-quality factor)
            "CPG" = Cp-G (capacitance-conductance)
            "CSRS" = Cs-Rs (series capacitance-resistance)
            "CSD" = Cs-D (series capacitance-dissipation)
            "ZTD" = Z-θ (impedance magnitude-phase)
            "ZTR" = Z-θr (impedance-phase in radians)
            "GB" = G-B (conductance-susceptance)
            "YTD" = Y-θ (admittance-phase)
            "RX" = R-X (resistance-reactance)
            "LPQ" = Lp-Q (parallel inductance-quality)
            "LSD" = Ls-D (series inductance-dissipation)
    """
    inst.write("*RST")
    inst.write("*CLS")
    
    # Set measurement function
    inst.write(f":FUNC:IMP:TYPE {measurement_function}")
    
    # Set trigger mode to internal (immediate)
    inst.write(":TRIG:SOUR INT")
    
    # Set data format to ASCII
    inst.write("FORM ASC")
    
    # Set aperture (integration time) to medium, 1 average
    inst.write(f":APER {MEAS_TIME},1")
    
    print(f"E4980A initialized for {measurement_function} measurement")

def check_errors(inst):
    """Check for instrument errors."""
    errors = []
    try:
        error = inst.query(":SYST:ERR?")
        if not error.startswith('+0') and "No error" not in error:
            print(f"Instrument error: {error.strip()}")
            errors.append(error.strip())
    except Exception as e:
        print(f"Error checking errors: {e}")
    return errors

# =============================
# Measurement Configuration Functions
# =============================

def set_frequency(inst, freq_hz):
    """
    Set measurement frequency.
    
    Args:
        freq_hz: Frequency in Hz (20 Hz to 2 MHz for E4980A)
    """
    if freq_hz < FREQ_MIN_HZ or freq_hz > FREQ_MAX_HZ:
        print(f"Warning: Frequency {freq_hz} Hz out of range ({FREQ_MIN_HZ}-{FREQ_MAX_HZ} Hz)")
        freq_hz = np.clip(freq_hz, FREQ_MIN_HZ, FREQ_MAX_HZ)
    
    inst.write(f":FREQ:CW {freq_hz}")
    current_parameters["frequency"] = freq_hz

def set_ac_level(inst, ac_level_v):
    """
    Set AC signal level.
    
    Args:
        ac_level_v: AC voltage level in V (0 to 20V for E4980A)
    """
    if ac_level_v < AC_LEVEL_MIN_V or ac_level_v > AC_LEVEL_MAX_V:
        print(f"Warning: AC level {ac_level_v} V out of range ({AC_LEVEL_MIN_V}-{AC_LEVEL_MAX_V} V)")
        ac_level_v = np.clip(ac_level_v, AC_LEVEL_MIN_V, AC_LEVEL_MAX_V)
    
    inst.write(f":VOLT:LEV {ac_level_v}")
    current_parameters["ac_level_v"] = ac_level_v

def configure_dc_bias(inst, apply_bias, dc_bias_v=0.0):
    """
    Configure DC bias.
    
    Args:
        apply_bias: Boolean to enable/disable DC bias
        dc_bias_v: DC bias voltage in V (-40V to +40V)
    """
    if apply_bias:
        if dc_bias_v < DC_BIAS_MIN_V or dc_bias_v > DC_BIAS_MAX_V:
            print(f"Warning: DC bias {dc_bias_v} V out of range ({DC_BIAS_MIN_V}-{DC_BIAS_MAX_V} V)")
            dc_bias_v = np.clip(dc_bias_v, DC_BIAS_MIN_V, DC_BIAS_MAX_V)
        
        inst.write(f":BIAS:VOLT {dc_bias_v}")
        inst.write(":BIAS:STAT ON")
        current_parameters["apply_dc_bias"] = True
        current_parameters["dc_bias_v"] = dc_bias_v
    else:
        inst.write(":BIAS:STAT OFF")
        current_parameters["apply_dc_bias"] = False
        current_parameters["dc_bias_v"] = 0.0

def set_integration_time(inst, time_setting="MED"):
    """
    Set measurement integration time (aperture).
    
    Args:
        time_setting: "SHOR" (short), "MED" (medium), "LONG", or "LONG2" (E4980A specific)
    """
    valid_settings = ["SHOR", "MED", "LONG", "LONG2"]
    if time_setting not in valid_settings:
        print(f"Warning: Invalid integration time '{time_setting}'. Using 'MED'.")
        time_setting = "MED"
    
    inst.write(f":APER {time_setting},1")
    current_parameters["meas_time"] = time_setting

def set_measurement_range(inst, auto=True, range_val=None):
    """
    Set measurement range.
    
    Args:
        auto: If True, use auto-ranging
        range_val: Fixed range value if auto=False (instrument-specific)
    """
    if auto:
        inst.write(":FUNC:IMP:RANG:AUTO ON")
    else:
        inst.write(":FUNC:IMP:RANG:AUTO OFF")
        if range_val is not None:
            inst.write(f":FUNC:IMP:RANG {range_val}")

# =============================
# Sweep Configuration Functions
# =============================

def configure_frequency_sweep(inst, start_freq, stop_freq, num_points, sweep_type="LIN"):
    """
    Configure frequency sweep using LIST mode.
    
    The E4980A uses :LIST:FREQ with explicit comma-separated frequency values,
    not start/stop/points commands. We generate the frequency list and send it.
    
    Args:
        start_freq: Start frequency in Hz
        stop_freq: Stop frequency in Hz
        num_points: Number of sweep points
        sweep_type: "LIN" for linear or "LOG" for logarithmic
    """
    # Generate frequency list
    if sweep_type.upper() == "LOG":
        freq_list = np.logspace(np.log10(start_freq), np.log10(stop_freq), num_points)
    else:
        freq_list = np.linspace(start_freq, stop_freq, num_points)
    
    # E4980A LIST mode: send comma-separated frequency values
    inst.write(":DISP:PAGE LIST")
    inst.write(":LIST:MODE SEQ")
    freq_str = ",".join([f"{f:.6e}" for f in freq_list])
    inst.write(f":LIST:FREQ {freq_str}")
    
    current_parameters["sweep_points"] = num_points
    
    # Store for later retrieval
    return freq_list

def configure_voltage_sweep(inst, start_v, stop_v, num_points):
    """
    Configure DC bias voltage sweep using LIST mode.
    
    Args:
        start_v: Start voltage in V
        stop_v: Stop voltage in V
        num_points: Number of sweep points
    """
    # Generate voltage list
    voltage_list = np.linspace(start_v, stop_v, num_points)
    
    # E4980A LIST mode for DC bias sweep
    inst.write(":DISP:PAGE LIST")
    inst.write(":LIST:MODE SEQ")
    bias_str = ",".join([f"{v:.6e}" for v in voltage_list])
    inst.write(f":LIST:BIAS:VOLT {bias_str}")
    inst.write(":BIAS:STAT ON")
    
    return voltage_list

# =============================
# Single Point Measurement Functions
# =============================

def measure_single_point(inst):
    """
    Perform a single-point measurement and return primary and secondary values.
    
    The E4980A :FETCH? returns 4 values: [A, B, status_A, status_B].
    Values A and B correspond to the selected measurement function.
    
    Returns:
        tuple: (primary_value, secondary_value) based on current measurement function
               e.g., for CPD: (Cp, D)
    """
    # Ensure trigger source is set for bus triggering
    inst.write(":TRIG:SOUR BUS")
    inst.write(":TRIG:IMM")
    inst.query("*OPC?")  # Wait for measurement complete
    
    # Fetch both primary and secondary parameters
    result_str = inst.query(":FETCH?")
    values = result_str.strip().split(',')
    
    if len(values) >= 2:
        primary = float(values[0])
        secondary = float(values[1])
        return primary, secondary
    else:
        raise ValueError(f"Unexpected measurement result format: {result_str}")

def measure_impedance(inst, freq_hz=None):
    """
    Measure impedance at specified frequency.
    
    Args:
        freq_hz: Frequency in Hz (uses current setting if None)
    
    Returns:
        tuple: (Z_magnitude, theta_degrees)
    """
    if freq_hz is not None:
        set_frequency(inst, freq_hz)
    
    # Ensure we're in Z-theta mode
    inst.write(":FUNC:IMP:TYPE ZTD")
    
    z_mag, theta_deg = measure_single_point(inst)
    return z_mag, theta_deg

def measure_capacitance(inst, freq_hz=None, mode="CPD"):
    """
    Measure capacitance at specified frequency.
    
    Args:
        freq_hz: Frequency in Hz (uses current setting if None)
        mode: "CPD" (Cp-D), "CPRP" (Cp-Rp), "CPG" (Cp-G), "CSRS" (Cs-Rs), "CSD" (Cs-D)
    
    Returns:
        tuple: (capacitance, secondary_parameter)
               For CPD: (Cp, D)
               For CPRP: (Cp, Rp)
    """
    if freq_hz is not None:
        set_frequency(inst, freq_hz)
    
    inst.write(f":FUNC:IMP:TYPE {mode}")
    
    cap, secondary = measure_single_point(inst)
    return cap, secondary

def measure_inductance(inst, freq_hz=None, mode="LPQ"):
    """
    Measure inductance at specified frequency.
    
    Args:
        freq_hz: Frequency in Hz (uses current setting if None)
        mode: "LPQ" (Lp-Q), "LSD" (Ls-D), "LPRP" (Lp-Rp), "LSRS" (Ls-Rs)
    
    Returns:
        tuple: (inductance, secondary_parameter)
    """
    if freq_hz is not None:
        set_frequency(inst, freq_hz)
    
    inst.write(f":FUNC:IMP:TYPE {mode}")
    
    ind, secondary = measure_single_point(inst)
    return ind, secondary

# =============================
# Sweep Measurement Functions
# =============================

def measure_frequency_sweep(inst, start_freq, stop_freq, num_points, sweep_type="LOG", 
                            measurement_function="CPD"):
    """
    Perform frequency sweep measurement.
    
    Args:
        start_freq: Start frequency in Hz
        stop_freq: Stop frequency in Hz
        num_points: Number of points
        sweep_type: "LIN" or "LOG"
        measurement_function: Measurement type (e.g., "CPD", "ZTD")
    
    Returns:
        tuple: (freq_array, primary_array, secondary_array)
    """
    # Set measurement function
    inst.write(f":FUNC:IMP:TYPE {measurement_function}")
    
    # Configure sweep — returns the generated frequency array
    freq_array = configure_frequency_sweep(inst, start_freq, stop_freq, num_points, sweep_type)
    
    # Trigger list sweep using bus trigger
    inst.write(":TRIG:SOUR BUS")
    inst.write(":INIT:CONT ON")
    inst.write(":TRIG:IMM")
    
    # Wait for sweep completion (robust for large point counts)
    inst.query("*OPC?")
    
    # :FETCh:IMPedance:FORMatted? returns 4 values per point: [A, B, status_A, status_B]
    data_str = inst.query(":FETCh:IMPedance:FORMatted?")
    data = [float(x) for x in data_str.strip().split(',') if x.strip()]
    
    # Return to safe trigger state
    inst.write(":TRIG:SOUR HOLD")
    
    # Extract primary (A) and secondary (B) from 4-tuples
    # Data format: [A1, B1, status1a, status1b, A2, B2, status2a, status2b, ...]
    primary_vals = np.array([data[i] for i in range(0, len(data), 4)])
    secondary_vals = np.array([data[i] for i in range(1, len(data), 4)])
    
    return freq_array, primary_vals, secondary_vals

def measure_impedance_vs_frequency(inst, start_freq, stop_freq, num_points, sweep_type="LOG"):
    """
    Measure impedance (Z-θ) vs frequency.
    
    Returns:
        tuple: (freq_array, z_magnitude, theta_degrees)
    """
    return measure_frequency_sweep(inst, start_freq, stop_freq, num_points, 
                                   sweep_type, measurement_function="ZTD")

def measure_capacitance_vs_frequency(inst, start_freq, stop_freq, num_points, 
                                     sweep_type="LOG", mode="CPD"):
    """
    Measure capacitance vs frequency.
    
    Args:
        mode: "CPD" for Cp-D, "CSRS" for Cs-Rs, etc.
    
    Returns:
        tuple: (freq_array, capacitance_vals, secondary_vals)
               For CPD: (freq, Cp, D)
    """
    return measure_frequency_sweep(inst, start_freq, stop_freq, num_points, 
                                   sweep_type, measurement_function=mode)

def measure_cv_sweep(inst, freq_hz, v_min, v_max, num_points):
    """
    Perform C-V (capacitance vs voltage) sweep at fixed frequency.
    
    Args:
        freq_hz: Measurement frequency
        v_min: Minimum DC bias voltage
        v_max: Maximum DC bias voltage
        num_points: Number of voltage points
    
    Returns:
        tuple: (voltage_array, capacitance_array, dissipation_array)
    """
    # Set CW frequency (used outside list mode as well)
    set_frequency(inst, freq_hz)
    
    # Set measurement function to CPD
    inst.write(":FUNC:IMP:TYPE CPD")
    
    # Configure voltage sweep — returns the voltage array
    voltage_array = configure_voltage_sweep(inst, v_min, v_max, num_points)
    
    # Also set a constant-frequency list so LIST mode uses the chosen freq
    inst.write(f":LIST:FREQ {freq_hz:.6e}")
    
    # Trigger list sweep
    inst.write(":TRIG:SOUR BUS")
    inst.write(":INIT:CONT ON")
    inst.write(":TRIG:IMM")
    
    # Wait for sweep completion
    inst.query("*OPC?")
    data_str = inst.query(":FETCh:IMPedance:FORMatted?")
    data = [float(x) for x in data_str.strip().split(',') if x.strip()]
    
    # Return to safe trigger state
    inst.write(":TRIG:SOUR HOLD")
    
    # Extract Cp and D from 4-tuples
    cp_vals = np.array([data[i] for i in range(0, len(data), 4)])
    d_vals = np.array([data[i] for i in range(1, len(data), 4)])
    
    return voltage_array, cp_vals, d_vals

def measure_cv_butterfly(inst, freq_hz, v_min, v_max, num_points):
    """
    Perform bidirectional C-V butterfly measurement (up and down sweep).
    
    Args:
        freq_hz: Measurement frequency
        v_min: Minimum DC bias voltage
        v_max: Maximum DC bias voltage
        num_points: Number of voltage points per direction
    
    Returns:
        tuple: (voltage_full, capacitance_full, dissipation_full)
               Arrays contain concatenated up and down sweep
    """
    # Up sweep: v_min → v_max
    v_up, cp_up, d_up = measure_cv_sweep(inst, freq_hz, v_min, v_max, num_points)
    
    # Down sweep: v_max → v_min
    v_down, cp_down, d_down = measure_cv_sweep(inst, freq_hz, v_max, v_min, num_points)
    
    # Concatenate
    v_full = np.concatenate([v_up, v_down])
    cp_full = np.concatenate([cp_up, cp_down])
    d_full = np.concatenate([d_up, d_down])
    
    return v_full, cp_full, d_full

# =============================
# Dielectric Parameter Calculations
# =============================

def compute_eps_r(capacitance, thickness_nm, diameter_um):
    """
    Compute relative permittivity (dielectric constant) from capacitance.
    
    εr = C * t / (ε0 * A)
    
    Args:
        capacitance: Capacitance in F
        thickness_nm: Dielectric thickness in nm
        diameter_um: Electrode diameter in µm
    
    Returns:
        Relative permittivity (dimensionless)
    """
    eps0 = 8.854e-12  # F/m (vacuum permittivity)
    t_m = thickness_nm * 1e-9
    d_m = diameter_um * 1e-6
    area_m2 = np.pi * (d_m / 2.0) ** 2
    
    eps_r = capacitance * t_m / (eps0 * area_m2)
    return eps_r

def compute_impedance_components(z_mag, theta_deg):
    """
    Compute real and imaginary impedance from magnitude and phase.
    
    Args:
        z_mag: Impedance magnitude (Ω)
        theta_deg: Phase angle (degrees)
    
    Returns:
        tuple: (Z_real, Z_imaginary)
    """
    theta_rad = np.deg2rad(theta_deg)
    z_real = z_mag * np.cos(theta_rad)
    z_imag = z_mag * np.sin(theta_rad)
    return z_real, z_imag

# =============================
# Parameter Management
# =============================

def get_measurement_parameters(use_last=False):
    """
    Get measurement parameters from user or use last values.
    
    Args:
        use_last: If True, return stored parameters without prompting
    
    Returns:
        dict: Dictionary of measurement parameters
    """
    global current_parameters
    
    if use_last:
        return current_parameters.copy()
    
    # Prompt user for parameters
    freq_input = input(f"Enter frequency (Hz) [default {current_parameters['frequency']:.2e}]: ").strip()
    freq = float(freq_input) if freq_input else current_parameters["frequency"]
    
    ac_input = input(f"Enter AC level (V) [default {current_parameters['ac_level_v']}]: ").strip()
    ac_level = float(ac_input) if ac_input else current_parameters["ac_level_v"]
    
    bias_input = input(f"Apply DC bias? (y/n) [default {'y' if current_parameters['apply_dc_bias'] else 'n'}]: ").strip().lower()
    apply_bias = bias_input == 'y' if bias_input else current_parameters["apply_dc_bias"]
    
    if apply_bias:
        bias_v_input = input(f"Enter DC bias (V) [default {current_parameters['dc_bias_v']}]: ").strip()
        bias_v = float(bias_v_input) if bias_v_input else current_parameters["dc_bias_v"]
    else:
        bias_v = 0.0
    
    # Update current parameters
    current_parameters["frequency"] = freq
    current_parameters["ac_level_v"] = ac_level
    current_parameters["apply_dc_bias"] = apply_bias
    current_parameters["dc_bias_v"] = bias_v
    
    return current_parameters.copy()

def prompt_for_parameter_duplication():
    """Ask user if they want to use last measurement parameters."""
    dup_input = input("Use last measurement parameters? (y/n): ").strip().lower()
    return dup_input == 'y'

# =============================
# Utility Functions
# =============================

def setup(resource_name=None, measurement_function="CPD"):
    """
    Connect to and initialize E4980A.
    
    Args:
        resource_name: GPIB address (uses default if None)
        measurement_function: Initial measurement type
    
    Returns:
        VISA instrument instance
    """
    inst = connect_e4980a(resource_name)
    initialize_e4980a(inst, measurement_function)
    return inst

def reset_instrument(inst):
    """Reset instrument to default state."""
    inst.write("*RST")
    inst.write("*CLS")
    time.sleep(1)
    print("E4980A reset to default state")


# =============================
# Open / Short / Load Correction
# =============================

def perform_open_correction(inst):
    """
    Execute open-circuit correction.

    The user should remove the DUT (leave probes/fixture open) before
    calling this function.  The E4980A measures stray admittance and
    stores it for later subtraction.
    """
    print("Performing OPEN correction — ensure DUT is removed / probes are open...")
    inst.write(":CORR:OPEN")
    time.sleep(3)          # typical correction takes 1-3 s
    inst.query("*OPC?")
    inst.write(":CORR:OPEN:STAT ON")
    print("Open correction complete and enabled.")


def perform_short_correction(inst):
    """
    Execute short-circuit correction.

    The user should short the fixture/probes before calling this.
    """
    print("Performing SHORT correction — ensure probes are shorted...")
    inst.write(":CORR:SHOR")
    time.sleep(3)
    inst.query("*OPC?")
    inst.write(":CORR:SHOR:STAT ON")
    print("Short correction complete and enabled.")


def perform_load_correction(inst, ref_r, ref_x=0.0):
    """
    Execute load correction with a known reference impedance.

    Args:
        ref_r: Reference resistance (Ω)
        ref_x: Reference reactance (Ω), default 0
    """
    print(f"Performing LOAD correction (R={ref_r} Ω, X={ref_x} Ω)...")
    inst.write(f":CORR:LOAD:TYPE CKIT")
    inst.write(f":CORR:LOAD:STAN:REF:PRIM {ref_r}")
    inst.write(f":CORR:LOAD:STAN:REF:SEC {ref_x}")
    inst.write(":CORR:LOAD")
    time.sleep(3)
    inst.query("*OPC?")
    inst.write(":CORR:LOAD:STAT ON")
    print("Load correction complete and enabled.")


def disable_corrections(inst):
    """Disable all fixture corrections."""
    inst.write(":CORR:OPEN:STAT OFF")
    inst.write(":CORR:SHOR:STAT OFF")
    inst.write(":CORR:LOAD:STAT OFF")
    print("All corrections disabled.")


# =============================
# Additional Sweep Measurements
# =============================

def measure_quality_factor_vs_frequency(inst, start_freq, stop_freq, num_points,
                                         sweep_type="LOG", mode="CPQ"):
    """
    Measure quality factor Q vs frequency.

    Uses CPQ (Cp-Q) or LPQ (Lp-Q) mode depending on the DUT.

    Args:
        mode: "CPQ" for capacitor Q, "LPQ" for inductor Q

    Returns:
        tuple: (freq_array, primary_vals, q_vals)
    """
    return measure_frequency_sweep(inst, start_freq, stop_freq, num_points,
                                   sweep_type, measurement_function=mode)


def measure_rx_vs_frequency(inst, start_freq, stop_freq, num_points,
                            sweep_type="LOG"):
    """
    Measure R-X (resistance–reactance) vs frequency.

    Returns:
        tuple: (freq_array, r_vals, x_vals)
    """
    return measure_frequency_sweep(inst, start_freq, stop_freq, num_points,
                                   sweep_type, measurement_function="RX")


def measure_gb_vs_frequency(inst, start_freq, stop_freq, num_points,
                            sweep_type="LOG"):
    """
    Measure G-B (conductance–susceptance) vs frequency.

    Returns:
        tuple: (freq_array, g_vals, b_vals)
    """
    return measure_frequency_sweep(inst, start_freq, stop_freq, num_points,
                                   sweep_type, measurement_function="GB")
