"""
Filename: pspa.py
Author: Ethan Ruddell
Date: 2026-03-02
Description: Driver for the Agilent 4155C / 4156C Semiconductor Parameter Analyzer.

The 4155C/4156C has up to four SMU (Source/Measure Unit) channels.  Each SMU
can force a voltage and measure current, or force a current and measure
voltage.  This makes it ideal for I-V curves, transistor characterisation,
diode testing, breakdown measurements, and more.

IMPORTANT: This module uses the FLEX command set, NOT SCPI.
FLEX (also called the "US" command mode) is the older native language of
the instrument.  Key commands:
  US       – switch the instrument into FLEX mode
  FMT 1    – set output format to ASCII with header
  CN / CL  – enable / disable SMU channels
  DV / DI  – force DC voltage / current
  TI? / TV? – measure current / voltage (high-speed spot)
  PV       – set pulse voltage parameters
  MM       – select measurement mode (spot, sweep, pulsed, etc.)
  XE       – execute measurement
  RMD?     – read measurement data

Typical workflow:
  1. connect_pspa()             – open GPIB session, enter FLEX mode
  2. configure_smu()            – set each channel as voltage or current source
  3. measure_*()                – perform the desired measurement (sweep / pulse)
  4. disconnect_pspa()          – disable outputs, return to SCPI, close session
"""

import pyvisa
import numpy as np
import re
import time
import logging

from utility.gpib_utils import create_visa_resource_manager

log = logging.getLogger(__name__)

# =============================
# User settings and constants
# =============================
GPIB_ADDRESS = "GPIB0::16::INSTR"   # Default GPIB address for the PSPA
DEFAULT_TIMEOUT_MS = 120000

# Global state to track how each SMU channel is configured.
# Keys: Channel number (1-4)
# Values: dict with 'mode' (VOLT or CURR), 'compliance' (safety limit), 'range' (0=auto).
SMU_CONFIG = {}

# Global state for pulse timing configuration.
# Keys: 'width' (float, seconds), 'period' (float, seconds), 'hold' (float, seconds).
PULSE_CONFIG = {}

# Global state for measurement integration time.
# Stored as FLEX/SCPI-compatible keywords: SHOR, MED, LONG.
INTEGRATION_TIME = "MED"

# =============================
# Connection and Initialization
# =============================
# FLEX mode replies embed status letters before the number (e.g. '128AI+5.33E-11').
# The regex below extracts the numeric part from any FLEX reply string.

_NUM_RE = re.compile(r'([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)')

def parse_flex_number(s: str) -> float:
    """
    Extract the *last* floating-point number from a FLEX reply.
    Example: '128AI+5.331500E-11' -> 5.331500e-11
    """
    s = s.strip()
    nums = _NUM_RE.findall(s)
    if not nums:
        raise ValueError(f"No numeric value found in reply: {s!r}")
    return float(nums[-1])

def list_visa_resources():
    """List all available VISA resources."""
    try:
        rm = create_visa_resource_manager()
        resources = rm.list_resources()
        return list(resources)
    except Exception as e:
        print(f"Error listing VISA resources: {e}")
        return []

def connect_pspa(gpib_address=None):
    """
    Connect to PSPA via GPIB and initialize FLEX mode.
    """
    if gpib_address is None:
        gpib_address = GPIB_ADDRESS
    
    if not gpib_address:
        print("\nError: GPIB address not set!")
        list_visa_resources()
        raise ValueError("GPIB address not configured")
    
    try:
        rm = create_visa_resource_manager()
        print(f"Attempting to connect to: {gpib_address}")
        pspa = rm.open_resource(gpib_address)
        pspa.timeout = DEFAULT_TIMEOUT_MS
        
        # Clear interface
        pspa.clear()
        
        # Initialize 4155C/4156C into FLEX mode
        # US: User Switch (switches to FLEX command set)
        pspa.write("US")
        
        # Set Data Format to ASCII with header [cite: 1302]
        pspa.write("FMT 1")
        
        # Disable all channels initially 
        pspa.write("CL")
        
        print("Connected to PSPA and initialized FLEX mode.")
        return pspa
        
    except Exception as e:
        print(f"Error connecting: {e}")
        raise

def check_errors(pspa):
    """Check for instrument errors using standard SCPI command (valid in FLEX)."""
    # Note: :SYST:ERR? is one of the few SCPI commands valid in FLEX 
    errors = []
    try:
        err_query = pspa.query(":SYST:ERR?")
        if "No error" not in err_query and not err_query.startswith('+0'):
            print(f"Instrument error: {err_query.strip()}")
            errors.append(err_query.strip())
    except Exception:
        pass
    return errors

def reset_instrument(pspa):
    """
    Perform a soft reset and reinitialize the instrument.
    Call this if the instrument becomes unresponsive or enters an error state.
    """
    print("Attempting instrument reset...")
    try:
        # Soft reset
        pspa.write("*RST")
        time.sleep(1.0)
        
        # Clear errors
        pspa.write("*CLS")
        time.sleep(0.3)
        
        # Re-initialize FLEX mode
        pspa.write("US")
        time.sleep(0.3)
        pspa.write("FMT 1")
        time.sleep(0.3)
        pspa.write("CL")
        
        print("Instrument reset complete.")
        return True
    except Exception as e:
        print(f"Reset failed: {e}")
        return False

def diagnose_communication(pspa):
    """
    Diagnostic function to test instrument communication and report status.
    """
    print("\n" + "=" * 70)
    print("INSTRUMENT DIAGNOSTICS")
    print("=" * 70)
    
    try:
        # Test basic connectivity
        print("\n1. Testing connectivity...")
        idn = pspa.query("*IDN?")
        print(f"   ID: {idn.strip()}")
        
        # Check timeout setting
        print(f"\n2. VISA Timeout: {pspa.timeout} ms")
        
        # Check for errors
        print("\n3. Checking instrument errors...")
        errors = check_errors(pspa)
        if errors:
            print(f"   Found {len(errors)} error(s):")
            for err in errors:
                print(f"     - {err}")
        else:
            print("   No errors reported.")
        
        # Test channel communication
        print("\n4. Testing channel measurements...")
        for ch in [1, 2]:
            try:
                current = measure_channel_current(pspa, ch, retry_count=1)
                print(f"   Channel {ch}: {current:.4e} A (OK)")
            except Exception as e:
                print(f"   Channel {ch}: FAILED - {e}")
        
        print("\n" + "=" * 70)
        
    except Exception as e:
        print(f"Diagnostics error: {e}")

def disconnect_pspa(pspa):
    """Safely disconnect from PSPA."""
    try:
        # Disable all channels 
        pspa.write("CL")
    except Exception:
        pass
    pspa.close()
    print("PSPA disconnected")

def measure_channel_current(pspa, channel):
    """
    Perform high-speed spot measurement of current using FLEX command TI?
    """
    # TI? chnum, range 
    # Range 0 = Auto Range [cite: 1303, 1409]
    try:
        result_str = pspa.query(f"TI? {channel},0")
        return parse_flex_number(result_str)
    except Exception as e:
        print(f"Error measuring current on CH{channel}: {e}")
        return 0.0


def measure_channel_voltage(pspa, channel):
    """Perform a high-speed spot measurement of voltage using FLEX command TV?."""
    try:
        result_str = pspa.query(f"TV? {channel},0")
        return parse_flex_number(result_str)
    except Exception as e:
        print(f"Error measuring voltage on CH{channel}: {e}")
        return 0.0


# =============================
# Plot Display Functions
# =============================
# Control the 4155C/4156C built-in plot display on the instrument screen.
# Uses SCPI PAGE mode commands to navigate and refresh the graph display.

def goto_plot_page(pspa):
    """
    Switch to PAGE (SCPI) mode and navigate to the plot/graph display page.
    The instrument will show real-time I-V curves as data is collected.
    
    Note: After this, use FLEX commands by sending 'US' again, or use SCPI ':PAGE:...' commands.
    """
    try:
        # Switch to SCPI PAGE mode (native 4155C/4156C SCPI subset)
        pspa.write(":PAGE")
        time.sleep(0.5)
        print("Switched to PAGE (SCPI) mode for plot display.")
    except Exception as e:
        print(f"Warning: Could not switch to plot page: {e}")

def configure_plot_display(pspa, y_scale="LINEAR", title="I-V Curve"):
    """
    Configure plot display axes and appearance.
    
    Args:
        y_scale: "LINEAR" or "LOG" for Y-axis (current) scale
        title: Display title on instrument (currently not settable, for reference)
    """
    try:
        # Set Y-axis scale: LOG for better range visibility, LINEAR for direct I values
        scale_cmd = "LOG" if y_scale.upper() == "LOG" else "LIN"
        pspa.write(f":PAGE:DISP:GRAP:Y:SCAL {scale_cmd}")
        
        # Enable automatic scaling on graph insert (ZCAN = Zero-Center Auto-scale)
        pspa.write(":PAGE:MEAS:MSET:ZCAN ON")
        
        print(f"Plot display configured: Y-axis={scale_cmd}, auto-scale enabled")
    except Exception as e:
        print(f"Warning: Could not configure plot display: {e}")

def refresh_plot(pspa):
    """
    Refresh and auto-scale the plot display on the instrument screen.
    Call this after adding a batch of data points to update the visible graph.
    """
    try:
        # Auto-scale graph once (GLIS:SCAL = Graph List Scaling)
        pspa.write(":PAGE:GLIS:SCAL:AUTO ONCE")
        time.sleep(0.2)
    except Exception as e:
        print(f"Warning: Could not refresh plot: {e}")

def switch_back_to_flex(pspa):
    """
    Switch back from PAGE (SCPI) mode to FLEX mode for measurements.
    Must be called before attempting more FLEX commands.
    """
    try:
        pspa.write("US")
        pspa.write("FMT 1")
        time.sleep(0.3)
        print("Switched back to FLEX mode.")
    except Exception as e:
        print(f"Warning: Could not switch back to FLEX: {e}")


# =============================
# Source Configuration Functions
# =============================
# These functions set up what each SMU channel does (force voltage or
# force current) and at what safety limit (compliance).

def configure_smu(pspa, channel, mode='VOLT', compliance=0.1):
    """
    Configure local config state for SMU. 
    FLEX forces configuration at the moment of forcing voltage/current.
    """
    global SMU_CONFIG
    SMU_CONFIG[channel] = {
        'mode': mode,
        'compliance': compliance,
        'range': 0  # 0 = Auto Range
    }
    # Match the established PSPA transfer-characteristics flow: stage the
    # selected channel immediately so subsequent DV/DI writes take effect.
    pspa.write(f"CN {channel}")
    print(f"Channel {channel} configured (State updated).")

def set_voltage(pspa, channel, voltage):
    """
    Force DC voltage using FLEX 'DV' command.
    Requires compliance from global config.
    """
    config = SMU_CONFIG.get(channel, {'compliance': 0.1, 'range': 0})
    comp = config['compliance']
    rng = config['range']
    
    # DV chnum, range, output, compliance 
    cmd = f"DV {channel},{rng},{voltage},{comp}"
    pspa.write(cmd)

def set_current(pspa, channel, current):
    """
    Force DC current using FLEX 'DI' command.
    """
    config = SMU_CONFIG.get(channel, {'compliance': 2.0, 'range': 0})
    # For Current source, compliance is Voltage
    comp = config['compliance'] 
    rng = config['range']
    
    # DI chnum, range, output, compliance 
    cmd = f"DI {channel},{rng},{current},{comp}"
    pspa.write(cmd)

def output_on(pspa, channel):
    """Enable output for a channel."""
    # CN [chnum] 
    pspa.write(f"CN {channel}")

def output_off(pspa, channel):
    """Disable output for a channel."""
    # CL [chnum] 
    pspa.write(f"CL {channel}")


def configure_constant_voltage_source(pspa, channel, voltage, compliance=0.1):
    """Compatibility wrapper for configuring a constant-voltage source channel."""
    channel = validate_channel(channel)
    validate_compliance(compliance)
    configure_smu(pspa, channel, mode='VOLT', compliance=compliance)
    set_voltage(pspa, channel, voltage)
    output_on(pspa, channel)


def start_constant_source(pspa):
    """Compatibility wrapper for enabling the already-configured constant source."""
    # The source channel is enabled by configure_constant_voltage_source().
    return None


def stop_constant_source(pspa):
    """Compatibility wrapper for disabling the constant source setup."""
    pspa.write("CL")


def read_constant_source_current(pspa, channel):
    """Compatibility wrapper for reading current from a configured source channel."""
    return measure_channel_current(pspa, channel)

# =============================
# Helper Functions
# =============================
# Basic sanity checks and a sweep-value generator used by all measurements.

def validate_channel(channel, allow_none=False):
    """Validate channel number is within valid range (1-4 for 4156C)."""
    if channel is None:
        if allow_none:
            return None
        raise ValueError("Channel number is required.")

    # Accept GUI/CLI numeric forms like "2" or 2.0.
    if not isinstance(channel, int):
        try:
            channel = int(float(channel))
        except Exception as exc:
            raise ValueError(f"Invalid channel number: {channel}. Must be 1-4.") from exc

    if channel < 1 or channel > 4:
        raise ValueError(f"Invalid channel number: {channel}. Must be 1-4.")
    return channel

def validate_compliance(compliance):
    """Validate compliance value is positive."""
    if compliance <= 0:
        raise ValueError(f"Compliance must be positive, got {compliance}")
    return compliance

def validate_step(step):
    """Validate step value is non-zero."""
    if step == 0:
        raise ValueError("Step must be non-zero")
    return step

def normalize_integration_time(time_setting):
    """Normalize user-friendly integration time inputs to SHOR/MED/LONG."""
    if time_setting is None:
        return "MED"

    token = str(time_setting).strip().upper()
    mapping = {
        "S": "SHOR",
        "SHOR": "SHOR",
        "SHORT": "SHOR",
        "M": "MED",
        "MED": "MED",
        "MEDIUM": "MED",
        "L": "LONG",
        "LONG": "LONG",
    }
    return mapping.get(token, "MED")

def set_integration_time(pspa, time_setting="MED"):
    """
    Set PSPA integration time for spot measurements.

    Args:
        time_setting: Accepts SHOR/MED/LONG and aliases S/M/L, SHORT/MEDIUM.
    """
    global INTEGRATION_TIME
    setting = normalize_integration_time(time_setting)

    # MATLAB reference uses :PAGE:MEAS:MSET:ITIM {SHOR|MED|LONG}.
    # Enter PAGE mode first, then return to FLEX because this driver operates
    # primarily in FLEX for measurement commands.
    try:
        pspa.write(":PAGE")
        time.sleep(0.2)
        pspa.write(f":PAGE:MEAS:MSET:ITIM {setting}")
        time.sleep(0.1)
        pspa.write("US")
        pspa.write("FMT 1")
    except Exception as e:
        print(f"Warning: failed to set integration time on instrument ({e}).")

    INTEGRATION_TIME = setting
    return setting

def sweep_values(start: float, stop: float, step: float) -> np.ndarray:
    """Generate sweep values from start to stop with given step magnitude."""
    if step == 0:
        raise ValueError("step must be non-zero")

    if start == stop:
        return np.array([float(start)], dtype=float)

    # Always sweep from start toward stop, regardless of step sign provided by user.
    signed_step = abs(step) if stop > start else -abs(step)

    # Include end point with a small margin for floating-point rounding.
    vals = np.arange(start, stop + (0.5 * signed_step), signed_step, dtype=float)

    # Ensure first/last values match requested endpoints exactly.
    if len(vals) == 0:
        vals = np.array([float(start), float(stop)], dtype=float)
    else:
        vals[0] = float(start)
        if not np.isclose(vals[-1], stop):
            vals = np.append(vals, float(stop))
        else:
            vals[-1] = float(stop)

    return vals


def configure_transfer_bias(pspa, drain_ch, source_ch, vds_constant,
                            compliance=0.1, integration_time="MED",
                            gate_ch=None):
    """Apply the standard PSPA transfer-style drain/source setup."""
    drain_ch = validate_channel(drain_ch)
    source_ch = validate_channel(source_ch)
    gate_ch = validate_channel(gate_ch, allow_none=True)
    validate_compliance(compliance)

    pspa.write("US")
    pspa.write("FMT 1")
    pspa.write("CL")

    configure_smu(pspa, drain_ch, mode='VOLT', compliance=compliance)
    if gate_ch is not None:
        configure_smu(pspa, gate_ch, mode='VOLT', compliance=compliance)
    configure_smu(pspa, source_ch, mode='VOLT', compliance=compliance)

    pspa.write(f"MM 1,{drain_ch}")
    set_integration_time(pspa, integration_time)

    set_voltage(pspa, source_ch, 0.0)
    set_voltage(pspa, drain_ch, vds_constant)

    if gate_ch is None:
        pspa.write(f"CN {drain_ch},{source_ch}")
    else:
        pspa.write(f"CN {drain_ch},{gate_ch},{source_ch}")

    time.sleep(0.5)
    return drain_ch, gate_ch, source_ch

# =============================
# Transistor I-V Measurements
# =============================
# The two standard transistor characterisations:
#   Output characteristics  – sweep drain voltage (Vds) at several gate voltages (Vgs)
#   Transfer characteristics – sweep gate voltage (bidirectional) at fixed drain voltage
# Both require three SMU channels: drain, gate, and source (grounded).

def measure_transistor_output_characteristics(pspa, vds_start, vds_stop, vds_step, 
                                               vgs_start, vgs_stop, vgs_step,
                                               drain_ch=1, gate_ch=2, source_ch=3,
                                               compliance=0.1):
    # Validate inputs
    drain_ch = validate_channel(drain_ch)
    gate_ch = validate_channel(gate_ch)
    source_ch = validate_channel(source_ch)
    validate_compliance(compliance)
    validate_step(vds_step)
    validate_step(vgs_step)
    
    # 1. Ensure FLEX mode is active and clear outputs from any prior run
    pspa.write("US")          # Switch to FLEX mode
    pspa.write("FMT 1")       # Set ASCII format with header
    pspa.write("CL")          # Force all SMUs off before reconfiguration
    
    # 2. Configure channels (Updates global state and enables channels)
    configure_smu(pspa, drain_ch, mode='VOLT', compliance=compliance)
    configure_smu(pspa, gate_ch, mode='VOLT', compliance=compliance)
    configure_smu(pspa, source_ch, mode='VOLT', compliance=compliance)
    
    # 3. Force 0V on source
    set_voltage(pspa, source_ch, 0)
    
    # 4. Enable all relevant channels at once 
    pspa.write(f"CN {drain_ch},{gate_ch},{source_ch}")

    # Small dwell to allow relays/output state to settle before first read.
    settle_s = 0.05
    time.sleep(settle_s)
    
    vds_array = []
    vgs_array = []
    id_array = []

    vgs_values = sweep_values(vgs_start, vgs_stop, vgs_step)
    vds_values = sweep_values(vds_start, vds_stop, vds_step)

    print(f"\n{'='*70}")
    print("TRANSISTOR OUTPUT CHARACTERISTICS MEASUREMENT")
    print(f"{'='*70}")
    print(f"SMU Assignment: Drain=CH{drain_ch}, Gate=CH{gate_ch}, Source=CH{source_ch}")
    print(f"Measurement Config: FLEX spot reads (TI?), settle={settle_s:.3f}s")
    print(f"Vds sweep: {vds_start}V -> {vds_stop}V (step {vds_step}V)")
    print(f"Vgs range: {vgs_start}V -> {vgs_stop}V (step {vgs_step}V)")
    print(f"Compliance: {compliance} A")
    print(f"Starting Output Characteristics Sweep...")
    print(f"{'='*70}\n")
    
    for vgs in vgs_values:
        set_voltage(pspa, gate_ch, vgs)
        set_voltage(pspa, source_ch, 0)
        time.sleep(settle_s)
        print(f"Vgs = {vgs:.3f}V : ", end="", flush=True)
        
        for vds in vds_values:
            set_voltage(pspa, drain_ch, vds)
            time.sleep(settle_s)
            
            # TI? executes measurement and returns data 
            id_val = measure_channel_current(pspa, drain_ch)
            
            vds_array.append(vds)
            vgs_array.append(vgs)
            id_array.append(id_val)
            print(".", end="", flush=True)
        
        print(f" [{len(vds_values)} points]")
    
    # Disable outputs 
    pspa.write("CL")

    vds_np = np.array(vds_array)
    vgs_np = np.array(vgs_array)
    id_np = np.array(id_array)

    # 4155/4156 current sign depends on SMU current direction convention.
    # Normalize Id so transistor plots and CSVs follow positive drain current.
    if id_np.size > 0 and np.nanmedian(id_np) < 0:
        id_np = -id_np
        print("Note: Id polarity normalized (sign inverted to positive convention).")
    
    return {
        'Vds': vds_np,
        'Vgs': vgs_np,
        'Id': id_np
    }

def measure_transistor_transfer_characteristics(pspa, vgs_start, vgs_stop, vgs_step,
                                                 vds_constant,
                                                 drain_ch=2, gate_ch=3, source_ch=1,
                                                 compliance=0.1, integration_time="MED"):
    """
    Bidirectional Vgs sweep for transfer characteristics.
    Sweeps Vgs from vgs_start -> vgs_stop (forward), then vgs_stop -> vgs_start
    (reverse) to capture switching hysteresis.

    Returns dict with 'Vgs', 'Id', and 'Sweep_Direction' arrays.
    'Sweep_Direction' contains 'forward' or 'reverse' for each data point.
    """
    # Validate inputs
    drain_ch = validate_channel(drain_ch)
    gate_ch = validate_channel(gate_ch)
    source_ch = validate_channel(source_ch)
    validate_compliance(compliance)
    validate_step(vgs_step)
    
    # Use the same drain/source setup sequence as the Keithley+PSPA transfer path.
    configure_transfer_bias(pspa, drain_ch, source_ch, vds_constant,
                            compliance=compliance, integration_time=integration_time,
                            gate_ch=gate_ch)

    vgs_array = []
    id_array = []
    direction_array = []

    # Forward sweep: vgs_start -> vgs_stop
    forward_values = sweep_values(vgs_start, vgs_stop, vgs_step)
    # Reverse sweep: vgs_stop -> vgs_start (step sign flipped)
    reverse_values = sweep_values(vgs_stop, vgs_start, -vgs_step)

    print(f"\n{'='*70}")
    print("TRANSISTOR TRANSFER CHARACTERISTICS MEASUREMENT")
    print(f"{'='*70}")
    print(f"SMU Assignment: Drain=CH{drain_ch}, Gate=CH{gate_ch}, Source=CH{source_ch}")
    print(f"Measurement Config: Spot Mode (MM 1), Integration={integration_time}")
    print(f"Vds (constant): {vds_constant}V")
    print(f"Vgs sweep: {vgs_start}V -> {vgs_stop}V (step {vgs_step}V) - Bidirectional")
    print(f"Compliance: {compliance} A")
    print(f"{'='*70}\n")

    print(f"Forward sweep: {vgs_start}V -> {vgs_stop}V : ", end="", flush=True)
    for vgs in forward_values:
        set_voltage(pspa, gate_ch, vgs)

        # Id must be read from the drain SMU while Vgs is being swept.
        id_val = measure_channel_current(pspa, drain_ch)

        vgs_array.append(vgs)
        id_array.append(id_val)
        direction_array.append('forward')
        print(".", end="", flush=True)
    
    print(f" [{len(forward_values)} points]")

    print(f"Reverse sweep: {vgs_stop}V -> {vgs_start}V : ", end="", flush=True)
    for vgs in reverse_values:
        set_voltage(pspa, gate_ch, vgs)

        # Id must be read from the drain SMU while Vgs is being swept.
        id_val = measure_channel_current(pspa, drain_ch)

        vgs_array.append(vgs)
        id_array.append(id_val)
        direction_array.append('reverse')
        print(".", end="", flush=True)
    
    print(f" [{len(reverse_values)} points]")

    pspa.write("CL")

    vgs_np = np.array(vgs_array)
    id_np = np.array(id_array)
    dir_np = np.array(direction_array)

    # Normalize Id polarity to positive drain current convention.
    if id_np.size > 0 and np.nanmedian(id_np) < 0:
        id_np = -id_np
        print("Note: Id polarity normalized (sign inverted to positive convention).")

    return {
        'Vgs': vgs_np,
        'Id': id_np,
        'Sweep_Direction': dir_np
    }

# =============================
# General I-V Measurements
# =============================
# Single-channel voltage sweeps for generic two-terminal devices.
#   Unidirectional: 0 → V_max  (simple ramp)
#   Bidirectional:  0 → +V_max → -V_max → 0  (detects hysteresis)

def measure_iv_curve(pspa, v_start, v_stop, v_step, channel=1, compliance=0.1,
                     ground_ch=None, integration_time="MED", enable_plot=True):
    # Validate inputs
    channel = validate_channel(channel)
    ground_ch = validate_channel(ground_ch, allow_none=True)
    validate_compliance(compliance)
    validate_step(v_step)

    set_integration_time(pspa, integration_time)
    
    # Setup plot display ONCE at start (stay in PAGE mode briefly, then return to FLEX)
    if enable_plot:
        goto_plot_page(pspa)
        configure_plot_display(pspa, y_scale="LINEAR")
        switch_back_to_flex(pspa)
    
    configure_smu(pspa, channel, mode='VOLT', compliance=compliance)
    if ground_ch is not None:
        configure_smu(pspa, ground_ch, mode='VOLT', compliance=compliance)
        set_voltage(pspa, ground_ch, 0)
        pspa.write(f"CN {channel},{ground_ch}")
    output_on(pspa, channel)
    
    voltage_array = []
    current_array = []

    voltages = sweep_values(v_start, v_stop, v_step)

    # Collect all data in FLEX mode (no mode switching during sweep)
    for v in voltages:
        set_voltage(pspa, channel, v)
        i_val = measure_channel_current(pspa, channel)
        voltage_array.append(v)
        current_array.append(i_val)
    
    if ground_ch is not None:
        pspa.write("CL")
    else:
        output_off(pspa, channel)
    
    # Final plot refresh after all data collected (one time in PAGE mode)
    if enable_plot:
        try:
            goto_plot_page(pspa)
            refresh_plot(pspa)
            switch_back_to_flex(pspa)
        except Exception as e:
            print(f"Final plot update warning: {e}")
    
    # Convert to numpy arrays
    voltage_array = np.array(voltage_array)
    current_array = np.array(current_array)
    
    # Remove first 0V datapoint and any near-zero points from startup
    # Keep only points where we've actually applied meaningful voltage
    # mask = np.abs(voltage_array) >= 1e-6  # Filter out near-zero voltages
    # voltage_array = voltage_array[mask]
    # current_array = current_array[mask]
    
    return {'Voltage': voltage_array, 'Current': current_array}

def measure_iv_bidirectional(pspa, v_max, v_step, channel=1, compliance=0.1,
                             ground_ch=None, integration_time="MED", enable_plot=True):
    # Validate inputs
    channel = validate_channel(channel)
    ground_ch = validate_channel(ground_ch, allow_none=True)
    validate_compliance(compliance)
    validate_step(v_step)

    set_integration_time(pspa, integration_time)
    
    # Setup plot display ONCE at start (stay in PAGE mode briefly, then return to FLEX)
    if enable_plot:
        goto_plot_page(pspa)
        configure_plot_display(pspa, y_scale="LINEAR")
        switch_back_to_flex(pspa)
    
    configure_smu(pspa, channel, mode='VOLT', compliance=compliance)
    if ground_ch is not None:
        configure_smu(pspa, ground_ch, mode='VOLT', compliance=compliance)
        set_voltage(pspa, ground_ch, 0)
        pspa.write(f"CN {channel},{ground_ch}")
    output_on(pspa, channel)
    
    voltage_array = []
    current_array = []

    up = sweep_values(0, v_max, v_step)
    down = sweep_values(v_max - v_step, -v_max, -v_step)
    back = sweep_values(-v_max + v_step, 0, v_step)
    sweep = np.concatenate([up, down, back])
    
    # Collect all data in FLEX mode (no mode switching during sweep)
    for v in sweep:
        set_voltage(pspa, channel, v)
        i_val = measure_channel_current(pspa, channel)
        voltage_array.append(v)
        current_array.append(i_val)
    
    if ground_ch is not None:
        pspa.write("CL")
    else:
        output_off(pspa, channel)
    
    # Final plot refresh after all data collected (one time in PAGE mode)
    if enable_plot:
        try:
            goto_plot_page(pspa)
            refresh_plot(pspa)
            switch_back_to_flex(pspa)
        except Exception as e:
            print(f"Final plot update warning: {e}")
    
    # Convert to numpy arrays
    voltage_array = np.array(voltage_array)
    current_array = np.array(current_array)
    
    # Remove first 0V datapoint and any near-zero points from startup
    # Keep only points where we've actually applied meaningful voltage
    # mask = np.abs(voltage_array) >= 1e-6  # Filter out near-zero voltages
    # voltage_array = voltage_array[mask]
    # current_array = current_array[mask]
    
    return {'Voltage': voltage_array, 'Current': current_array}

# =============================
# Pulse Measurements
# =============================
# Pulsed measurements apply voltage in short bursts instead of continuously.
# This reduces self-heating in the device under test, giving more accurate
# results for high-power transistors and other heat-sensitive devices.
# The 4155C/4156C FLEX mode supports pulsed spot measurements (MM 3).

def configure_pulse_mode(pspa, channel, pulse_width, pulse_period):
    """
    Setup Pulse Timing. Note: In FLEX, PT sets timing globally for pulse sources.
    This function stores parameters to be used in the measurement function.
    """
    global PULSE_CONFIG
    PULSE_CONFIG = {
        'width': pulse_width,
        'period': pulse_period,
        'hold': 0.0
    }
    # We don't send PT here because PT requires hold/delay params which might vary
    print(f"Pulse config stored: {pulse_width}s width, {pulse_period}s period")

def measure_pulsed_iv(pspa, v_base, v_pulse, pulse_width, pulse_period, num_pulses,
                      channel=1, compliance=0.1):
    """
    Measure pulsed IV using '1 Channel Pulsed Spot Measurement' mode (MM 3).
    """
    # Validate inputs
    channel = validate_channel(channel)
    validate_compliance(compliance)
    if pulse_width <= 0 or pulse_period <= 0:
        raise ValueError("Pulse width and period must be positive")
    if pulse_width >= pulse_period:
        raise ValueError("Pulse width must be less than pulse period")
    if num_pulses <= 0:
        raise ValueError("Number of pulses must be positive")
    
    # 1. Setup
    pspa.write("US") # Ensure FLEX
    pspa.write("FMT 1") 
    
    # 2. Enable Channel
    pspa.write(f"CN {channel}")
    
    # 3. Setup Pulse Timing (PT hold, width, period) [cite: 1346, 1359]
    pspa.write(f"PT 0,{pulse_width},{pulse_period}")
    
    # 4. Setup Measurement Mode (MM 3 = 1ch pulsed spot) [cite: 1344]
    pspa.write(f"MM 3,{channel}")
    
    time_array = []
    voltage_array = []
    current_array = []
    
    print(f"Starting Pulsed IV (FLEX MM 3) on Ch {channel}...")
    
    # 5. Execution Loop
    for i in range(num_pulses):
        # Setup Pulse Voltage parameters: PV ch, range, base, pulse, compliance 
        # Range 0 = Auto
        pspa.write(f"PV {channel},0,{v_base},{v_pulse},{compliance}")
        
        # Execute Measurement [cite: 1344]
        pspa.write("XE")
        
        # Wait for completion (*OPC?) [cite: 1360]
        pspa.query("*OPC?")
        
        # Check errors [cite: 1360]
        err_res = pspa.query(":SYST:ERR?")
        if "+0" not in err_res and "No error" not in err_res:
            print(f"Pulse Error: {err_res}")
        
        # Read Measurement Data (RMD?) [cite: 1361]
        # MM 3 returns measured data at pulse peak
        try:
            val_str = pspa.query("RMD? 1")
            # Parse result (Format: status header + value)
            current_val = parse_flex_number(val_str)
        except Exception:
            current_val = 0.0
            
        # Store Data
        time_array.append(i * pulse_period)
        voltage_array.append(v_pulse)
        current_array.append(current_val)
        
    # 6. Cleanup
    pspa.write("CL")
    
    return {
        'Time': np.array(time_array),
        'Voltage': np.array(voltage_array),
        'Current': np.array(current_array)
    }

def measure_pulsed_transistor(pspa, vds_pulse, vgs_pulse, vds_base, vgs_base,
                               pulse_width, pulse_period, num_pulses,
                               drain_ch=1, gate_ch=2, source_ch=3, compliance=0.1):
    """
    Pulsed transistor measurement using FLEX.
    
    LIMITATION: This implementation pulses the Drain (Vds) while the Gate (Vgs) 
    is held at a constant DC level (vgs_pulse). The 4155C requires complex 
    'Staircase Sweep with Pulse' setup for synchronized dual-channel pulsing,
    which is not implemented here.
    
    The measurement pulses Vds between vds_base and vds_pulse while measuring
    drain current. Gate is held constant at vgs_pulse (vgs_base is not used
    in this simplified implementation).
    """
    # Validate inputs
    drain_ch = validate_channel(drain_ch)
    gate_ch = validate_channel(gate_ch)
    source_ch = validate_channel(source_ch)
    validate_compliance(compliance)
    if pulse_width <= 0 or pulse_period <= 0:
        raise ValueError("Pulse width and period must be positive")
    if pulse_width >= pulse_period:
        raise ValueError("Pulse width must be less than pulse period")
    if num_pulses <= 0:
        raise ValueError("Number of pulses must be positive")
    
    # Setup FLEX
    pspa.write("US")
    pspa.write("FMT 1")
    pspa.write(f"CN {drain_ch},{gate_ch},{source_ch}")
    pspa.write(f"PT 0,{pulse_width},{pulse_period}")
    
    # We use MM 3 (Pulsed Spot) on Drain. Gate is DC. [cite: 1346]
    pspa.write(f"MM 3,{drain_ch}")
    
    # Ground Source
    set_voltage(pspa, source_ch, 0)
    
    # Apply pulse voltage to Gate (DC - held constant during measurement)
    # NOTE: vgs_base is ignored in this implementation
    set_voltage(pspa, gate_ch, vgs_pulse) 
    
    time_array = []
    id_array = []
    
    print("Starting Pulsed Transistor (Pulsed Vds, DC Vgs)...")
    
    for i in range(num_pulses):
        # Pulse Drain: PV ch, range, base, pulse, compliance
        pspa.write(f"PV {drain_ch},0,{vds_base},{vds_pulse},{compliance}")
        
        pspa.write("XE")
        pspa.query("*OPC?")
        
        try:
            val = parse_flex_number(pspa.query("RMD? 1"))
        except Exception:
            val = 0.0
            
        time_array.append(i * pulse_period)
        id_array.append(val)
        
    pspa.write("CL")
    
    return {
        'Time': np.array(time_array),
        'Id': np.array(id_array)
    }


# =============================
# Additional Measurement Functions
# =============================
# Specialised measurements for diodes, gate leakage, breakdown voltage,
# and resistance (force-I / measure-V).

def measure_diode_iv(pspa, v_start, v_stop, v_step, anode_ch=1, cathode_ch=2,
                     compliance=0.1):
    """
    Measure diode I-V characteristics.

    Sweeps voltage on the anode while grounding the cathode, measures
    current through the anode at each bias point.

    Returns:
        dict with keys 'Voltage' (V) and 'Current' (A) as numpy arrays.
    """
    configure_smu(pspa, anode_ch, mode='VOLT', compliance=compliance)
    configure_smu(pspa, cathode_ch, mode='VOLT', compliance=compliance)

    set_voltage(pspa, cathode_ch, 0)
    pspa.write(f"CN {anode_ch},{cathode_ch}")

    voltages_arr = sweep_values(v_start, v_stop, v_step)
    v_out = []
    i_out = []

    print("Starting Diode I-V sweep...")
    for v in voltages_arr:
        set_voltage(pspa, anode_ch, v)
        i_val = measure_channel_current(pspa, anode_ch)
        v_out.append(v)
        i_out.append(i_val)

    pspa.write("CL")
    
    # Convert to numpy arrays
    v_out = np.array(v_out)
    i_out = np.array(i_out)
    
    # Remove first 0V datapoint and any near-zero points from startup
    # Keep only points where we've actually applied meaningful voltage
    mask = np.abs(v_out) >= 1e-6  # Filter out near-zero voltages
    v_out = v_out[mask]
    i_out = i_out[mask]
    
    return {'Voltage': v_out, 'Current': i_out}


def measure_gate_leakage(pspa, vgs_start, vgs_stop, vgs_step,
                         gate_ch=2, source_ch=3, compliance=1e-3):
    """
    Measure gate leakage current (Ig) as a function of Vgs.

    Forces voltage on the gate channel while grounding the source,
    and reads gate current at each step.

    Returns:
        dict with keys 'Vgs' (V) and 'Ig' (A) as numpy arrays.
    """
    configure_smu(pspa, gate_ch, mode='VOLT', compliance=compliance)
    configure_smu(pspa, source_ch, mode='VOLT', compliance=compliance)

    set_voltage(pspa, source_ch, 0)
    pspa.write(f"CN {gate_ch},{source_ch}")

    vgs_arr = sweep_values(vgs_start, vgs_stop, vgs_step)
    vgs_out = []
    ig_out = []

    print("Starting Gate Leakage sweep...")
    for vgs in vgs_arr:
        set_voltage(pspa, gate_ch, vgs)
        ig = measure_channel_current(pspa, gate_ch)
        vgs_out.append(vgs)
        ig_out.append(ig)

    pspa.write("CL")
    return {'Vgs': np.array(vgs_out), 'Ig': np.array(ig_out)}


def measure_breakdown_voltage(pspa, v_start, v_stop, v_step, channel=1,
                              compliance=1e-3, threshold_a=1e-4):
    """
    Sweep voltage upward until current exceeds *threshold_a*, indicating
    breakdown.

    Returns:
        dict with keys 'Voltage', 'Current' (full sweep data) and
        'breakdown_v' (float or None if threshold not reached).
    """
    configure_smu(pspa, channel, mode='VOLT', compliance=compliance)
    output_on(pspa, channel)

    voltages_arr = sweep_values(v_start, v_stop, v_step)
    v_out = []
    i_out = []
    breakdown_v = None

    print("Starting Breakdown Voltage sweep...")
    for v in voltages_arr:
        set_voltage(pspa, channel, v)
        i_val = measure_channel_current(pspa, channel)
        v_out.append(v)
        i_out.append(i_val)
        if abs(i_val) >= threshold_a and breakdown_v is None:
            breakdown_v = v
            print(f"  Breakdown detected at {v:.3f} V (I = {i_val:.4e} A)")
            break  # stop to protect DUT

    output_off(pspa, channel)
    
    # Convert to numpy arrays
    v_out = np.array(v_out)
    i_out = np.array(i_out)
    
    # Remove first 0V datapoint and any near-zero points from startup
    # Keep only points where we've actually applied meaningful voltage
    mask = np.abs(v_out) >= 1e-6  # Filter out near-zero voltages
    v_out = v_out[mask]
    i_out = i_out[mask]
    # Recalculate breakdown_v if it was filtered out
    if breakdown_v is not None and abs(breakdown_v) < 1e-6:
        breakdown_v = v_out[0] if len(v_out) > 0 else None
    
    return {
        'Voltage': v_out,
        'Current': i_out,
        'breakdown_v': breakdown_v,
    }


def measure_resistance(pspa, i_start, i_stop, i_step, channel=1,
                       compliance=10.0, sense_channel=None):
    """
    Four-point Kelvin resistance measurement: force current, measure voltage on separate channel.

    Forces current through *channel*, measures voltage on *sense_channel* (or same channel if None).
    This provides more accurate resistance measurement by avoiding lead resistance effects.
    
    Args:
        pspa: PyVISA instrument connection
        i_start: Starting current (A)
        i_stop: Stopping current (A)
        i_step: Current step size (A)
        channel: Channel to force current through (default: 1)
        compliance: Voltage compliance limit (default: 10.0 V)
        sense_channel: Channel to measure voltage on (default: None = use same channel as force)
                      For true 4-point Kelvin: use different channel (e.g., channel=1, sense_channel=2)

    Returns:
        dict with keys 'Current' (A), 'Voltage' (V), and
        'resistance_ohm' (float, slope of least-squares linear fit).
    """
    # Default to same channel if sense_channel not specified
    if sense_channel is None:
        sense_channel = channel
    
    config_backup_force = SMU_CONFIG.get(channel, {}).copy()
    config_backup_sense = SMU_CONFIG.get(sense_channel, {}).copy()
    
    # Configure force channel as current source
    configure_smu(pspa, channel, mode='CURR', compliance=compliance)
    # Configure sense channel as voltage meter
    configure_smu(pspa, sense_channel, mode='VOLT', compliance=compliance)
    
    # Enable both channels
    output_on(pspa, channel)
    output_on(pspa, sense_channel)

    currents_arr = sweep_values(i_start, i_stop, i_step)
    i_out = []
    v_out = []

    print(f"Starting Resistance measurement (Kelvin 4-point: force I on ch{channel}, measure V on ch{sense_channel})...")
    for i_val in currents_arr:
        set_current(pspa, channel, i_val)
        # Measure voltage on the sense channel
        v_str = pspa.query(f"TV? {sense_channel},0")
        v_val = parse_flex_number(v_str)
        i_out.append(i_val)
        v_out.append(v_val)

    # Turn off both channels
    output_off(pspa, channel)
    output_off(pspa, sense_channel)
    
    # Restore prior config
    if config_backup_force:
        SMU_CONFIG[channel] = config_backup_force
    if config_backup_sense:
        SMU_CONFIG[sense_channel] = config_backup_sense

    i_arr = np.array(i_out)
    v_arr = np.array(v_out)
    
    # Remove first 0A datapoint and any near-zero points
    mask = np.abs(i_arr) >= 1e-9
    i_arr = i_arr[mask]
    v_arr = v_arr[mask]
    
    # Least-squares linear fit
    if len(i_arr) > 1:
        coeffs = np.polyfit(i_arr, v_arr, 1)
        resistance = coeffs[0]
    else:
        resistance = v_arr[0] / i_arr[0] if i_arr[0] != 0 else float('inf')

    return {
        'Current': i_arr,
        'Voltage': v_arr,
        'resistance_ohm': resistance,
    }
