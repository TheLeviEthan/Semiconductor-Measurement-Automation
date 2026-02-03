"""
Filename: pspa.py
Author: Ethan Ruddell
Date: 2025-01-23
Description: Contains all constants and functions for PSPA measurements.
             Updated for Agilent 4155C/4156C FLEX Syntax.
"""

import pyvisa
import numpy as np
import time

# =============================
# User settings and constants
# =============================
GPIB_ADDRESS = "GPIB0::16::INSTR"   # Example: "GPIB0::17::INSTR"

# Global state to track SMU configuration (compliance, mode)
# Keys: Channel int (1-4), Values: {'mode': 'VOLT', 'compliance': 0.1, 'range': 0}
SMU_CONFIG = {}

# Global state for pulse configuration
PULSE_CONFIG = {}

# =============================
# Connection and Initialization
# =============================

import re

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
        rm = pyvisa.ResourceManager()
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
        rm = pyvisa.ResourceManager()
        print(f"Attempting to connect to: {gpib_address}")
        pspa = rm.open_resource(gpib_address)
        pspa.timeout = 30000
        
        # Clear interface
        pspa.clear()
        
        # Initialize 4155C/4156C into FLEX mode [cite: 1290, 1690]
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
    # Note: :SYST:ERR? is one of the few SCPI commands valid in FLEX [cite: 1711]
    errors = []
    try:
        err_query = pspa.query(":SYST:ERR?")
        if "No error" not in err_query and not err_query.startswith('+0'):
            print(f"Instrument error: {err_query.strip()}")
            errors.append(err_query.strip())
    except:
        pass
    return errors

def disconnect_pspa(pspa):
    """Safely disconnect from PSPA."""
    try:
        # Disable all channels 
        pspa.write("CL")
        # Return to SCPI mode (Page mode) [cite: 1316]
        pspa.write(":PAGE")
    except:
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


# =============================
# Source Configuration Functions
# =============================

def configure_smu(pspa, channel, mode='VOLT', compliance=0.1):
    """
    Configure local config state for SMU. 
    FLEX forces configuration at the moment of forcing voltage/current.
    """
    global SMU_CONFIG
    SMU_CONFIG[channel] = {
        'mode': mode,
        'compliance': compliance,
        'range': 0  # 0 = Auto Range [cite: 1409]
    }
    # Enable the channel 
    # Note: FLEX 'CN' command enables specified channels.
    # We should ensure we don't accidentally disable others, but for simplicity
    # this function just updates our state. The actual 'CN' is sent usually 
    # as a group or we can send it here.
    # Sending 'CN channel' enables that channel and leaves others as is (usually).
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

# =============================
# Helper Functions
# =============================

def validate_channel(channel):
    """Validate channel number is within valid range (1-4 for 4156C)."""
    if not isinstance(channel, int) or channel < 1 or channel > 4:
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

def sweep_values(start: float, stop: float, step: float) -> np.ndarray:
    """Generate sweep values from start to stop with given step."""
    if step == 0:
        raise ValueError("step must be non-zero")

    n = int(round((stop - start) / step)) + 1
    # Guard against sign mistakes
    if n <= 0:
        return np.array([], dtype=float)

    vals = start + step * np.arange(n, dtype=float)

    # Optional: hard-clip last value to exactly stop (prevents 1e-16 drift)
    if len(vals) > 0:
        vals[-1] = stop
    return vals

# =============================
# Transistor I-V Measurements
# =============================

def measure_transistor_output_characteristics(pspa, vds_start, vds_stop, vds_step, 
                                               vgs_start, vgs_stop, vgs_step,
                                               drain_ch=1, gate_ch=2, source_ch=3,
                                               compliance=0.1):
    # Validate inputs
    validate_channel(drain_ch)
    validate_channel(gate_ch)
    validate_channel(source_ch)
    validate_compliance(compliance)
    validate_step(vds_step)
    validate_step(vgs_step)
    
    # Configure channels (Updates global state and enables channels)
    configure_smu(pspa, drain_ch, mode='VOLT', compliance=compliance)
    configure_smu(pspa, gate_ch, mode='VOLT', compliance=compliance)
    configure_smu(pspa, source_ch, mode='VOLT', compliance=compliance)
    
    # Force 0V on source
    set_voltage(pspa, source_ch, 0)
    
    # Enable all relevant channels at once 
    pspa.write(f"CN {drain_ch},{gate_ch},{source_ch}")
    
    vds_array = []
    vgs_array = []
    id_array = []

    vgs_values = sweep_values(vgs_start, vgs_stop, vgs_step)
    vds_values = sweep_values(vds_start, vds_stop, vds_step)

    print("Starting Output Characteristics Sweep...")
    
    for vgs in vgs_values:
        set_voltage(pspa, gate_ch, vgs)
        
        for vds in vds_values:
            set_voltage(pspa, drain_ch, vds)
            
            # TI? executes measurement and returns data 
            id_val = measure_channel_current(pspa, drain_ch)
            
            vds_array.append(vds)
            vgs_array.append(vgs)
            id_array.append(id_val)
    
    # Disable outputs 
    pspa.write("CL")
    
    return {
        'Vds': np.array(vds_array),
        'Vgs': np.array(vgs_array),
        'Id': np.array(id_array)
    }

def measure_transistor_transfer_characteristics(pspa, vgs_start, vgs_stop, vgs_step,
                                                 vds_constant,
                                                 drain_ch=1, gate_ch=2, source_ch=3,
                                                 compliance=0.1):
    # Validate inputs
    validate_channel(drain_ch)
    validate_channel(gate_ch)
    validate_channel(source_ch)
    validate_compliance(compliance)
    validate_step(vgs_step)
    
    configure_smu(pspa, drain_ch, mode='VOLT', compliance=compliance)
    configure_smu(pspa, gate_ch, mode='VOLT', compliance=compliance)
    configure_smu(pspa, source_ch, mode='VOLT', compliance=compliance)
    
    set_voltage(pspa, source_ch, 0)
    set_voltage(pspa, drain_ch, vds_constant)
    
    pspa.write(f"CN {drain_ch},{gate_ch},{source_ch}")
    
    vgs_array = []
    id_array = []
    ig_array = []
    
    vgs_values = np.arange(vgs_start, vgs_stop + vgs_step, vgs_step)
    
    print("Starting Transfer Characteristics Sweep...")
    
    for vgs in vgs_values:
        set_voltage(pspa, gate_ch, vgs)
        
        id_val = measure_channel_current(pspa, drain_ch)
        ig_val = measure_channel_current(pspa, gate_ch)
        
        vgs_array.append(vgs)
        id_array.append(id_val)
        ig_array.append(ig_val)
    
    pspa.write("CL")
    
    return {
        'Vgs': np.array(vgs_array),
        'Id': np.array(id_array),
        'Ig': np.array(ig_array)
    }

# =============================
# General I-V Measurements
# =============================

def measure_iv_curve(pspa, v_start, v_stop, v_step, channel=1, compliance=0.1):
    # Validate inputs
    validate_channel(channel)
    validate_compliance(compliance)
    validate_step(v_step)
    
    configure_smu(pspa, channel, mode='VOLT', compliance=compliance)
    output_on(pspa, channel)
    
    voltage_array = []
    current_array = []

    voltages = sweep_values(v_start, v_stop, v_step)

    for v in voltages:
        set_voltage(pspa, channel, v)
        i_val = measure_channel_current(pspa, channel)
        voltage_array.append(v)
        current_array.append(i_val)
    
    output_off(pspa, channel)
    
    return {'Voltage': np.array(voltage_array), 'Current': np.array(current_array)}

def measure_iv_bidirectional(pspa, v_max, v_step, channel=1, compliance=0.1):
    # Validate inputs
    validate_channel(channel)
    validate_compliance(compliance)
    validate_step(v_step)
    
    configure_smu(pspa, channel, mode='VOLT', compliance=compliance)
    output_on(pspa, channel)
    
    voltage_array = []
    current_array = []

    up = sweep_values(0, v_max, v_step)
    down = sweep_values(v_max - v_step, -v_max, -v_step)
    back = sweep_values(-v_max + v_step, 0, v_step)
    sweep = np.concatenate([up, down, back])
    
    for v in sweep:
        set_voltage(pspa, channel, v)
        i_val = measure_channel_current(pspa, channel)
        voltage_array.append(v)
        current_array.append(i_val)
    
    output_off(pspa, channel)
    
    return {'Voltage': np.array(voltage_array), 'Current': np.array(current_array)}

# =============================
# Pulse Measurements
# =============================

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
    [cite: 1338, 1346]
    """
    # Validate inputs
    validate_channel(channel)
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
        except:
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
    validate_channel(drain_ch)
    validate_channel(gate_ch)
    validate_channel(source_ch)
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
            val = float(pspa.query("RMD? 1").strip().split()[-1])
        except:
            val = 0.0
            
        time_array.append(i * pulse_period)
        id_array.append(val)
        
    pspa.write("CL")
    
    return {
        'Time': np.array(time_array),
        'Id': np.array(id_array)
    }   