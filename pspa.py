"""
Filename: main.py
Author: Ethan Ruddell
Date: 2025-01-23
Description: Contains all constants and functions for PIA measurements.
"""

import pyvisa
import numpy as np


# =============================
# User settings and constants
# =============================
GPIB_ADDRESS = ""   # FIND AND RECORD GPIB ADDRESS FOR PSPA

# Default Frequency sweep settings (LOG sweep)
FREQ_START_HZ = 1e4      # start frequency, Hz (>0 for LOG sweep)
FREQ_STOP_HZ = 1e6       # stop frequency, Hz
NUM_POINTS = 201         # number of points in LOG sweep

# FIND MANUAL FOR PSPA TO GET SYNTAX (keysight)
# UTILIZES STANDARD SCPI COMMANDS
# https://en.wikipedia.org/wiki/Standard_Commands_for_Programmable_Instruments


# =============================
# Connection and Initialization
# =============================

def connect_pspa(gpib_address=None):
    """
    Connect to PSPA via GPIB.
    
    Args:
        gpib_address (str): GPIB address (e.g., "GPIB0::17::INSTR")
    
    Returns:
        pyvisa.Resource: PSPA instrument object
    """
    if gpib_address is None:
        gpib_address = GPIB_ADDRESS
    
    rm = pyvisa.ResourceManager()
    pspa = rm.open_resource(gpib_address)
    pspa.timeout = 30000  # 30 second timeout
    
    # Reset and clear
    pspa.write("*RST")
    pspa.write("*CLS")
    
    # Query instrument ID
    idn = pspa.query("*IDN?")
    print(f"Connected to: {idn.strip()}")
    
    return pspa


def disconnect_pspa(pspa):
    """
    Safely disconnect from PSPA.
    
    Args:
        pspa (pyvisa.Resource): PSPA instrument object
    """
    # Turn off all outputs
    pspa.write(":OUTP:STAT OFF")
    pspa.close()
    print("PSPA disconnected")


# =============================
# Source Configuration Functions
# =============================

def configure_smu(pspa, channel, mode='VOLT', compliance=0.1):
    """
    Configure a Source Measurement Unit (SMU) channel.
    
    Args:
        pspa (pyvisa.Resource): PSPA instrument object
        channel (int): Channel number (1-4)
        mode (str): 'VOLT' for voltage source or 'CURR' for current source
        compliance (float): Compliance limit (current for voltage source, voltage for current source)
    """
    pspa.write(f":SOUR{channel}:FUNC:MODE {mode}")
    
    if mode == 'VOLT':
        pspa.write(f":SENS{channel}:CURR:PROT {compliance}")  # Current compliance
    else:
        pspa.write(f":SENS{channel}:VOLT:PROT {compliance}")  # Voltage compliance
    
    print(f"Channel {channel} configured as {mode} source with compliance {compliance}")


def set_voltage(pspa, channel, voltage):
    """
    Set output voltage for a channel.
    
    Args:
        pspa (pyvisa.Resource): PSPA instrument object
        channel (int): Channel number
        voltage (float): Voltage in volts
    """
    pspa.write(f":SOUR{channel}:VOLT {voltage}")


def set_current(pspa, channel, current):
    """
    Set output current for a channel.
    
    Args:
        pspa (pyvisa.Resource): PSPA instrument object
        channel (int): Channel number
        current (float): Current in amps
    """
    pspa.write(f":SOUR{channel}:CURR {current}")


def output_on(pspa, channel):
    """Enable output for a channel."""
    pspa.write(f":OUTP{channel}:STAT ON")


def output_off(pspa, channel):
    """Disable output for a channel."""
    pspa.write(f":OUTP{channel}:STAT OFF")


# =============================
# Transistor I-V Measurements
# =============================

def measure_transistor_output_characteristics(pspa, vds_start, vds_stop, vds_step, 
                                               vgs_start, vgs_stop, vgs_step,
                                               drain_ch=1, gate_ch=2, source_ch=3,
                                               compliance=0.1):
    """
    Measure transistor output characteristics (Id vs Vds for various Vgs).
    
    Args:
        pspa (pyvisa.Resource): PSPA instrument object
        vds_start (float): Start drain-source voltage (V)
        vds_stop (float): Stop drain-source voltage (V)
        vds_step (float): Drain-source voltage step (V)
        vgs_start (float): Start gate-source voltage (V)
        vgs_stop (float): Stop gate-source voltage (V)
        vgs_step (float): Gate-source voltage step (V)
        drain_ch (int): Drain channel number
        gate_ch (int): Gate channel number
        source_ch (int): Source channel number (typically grounded)
        compliance (float): Current compliance (A)
    
    Returns:
        dict: Dictionary containing Vds, Vgs, and Id arrays
    """
    # Configure channels
    configure_smu(pspa, drain_ch, mode='VOLT', compliance=compliance)
    configure_smu(pspa, gate_ch, mode='VOLT', compliance=compliance)
    configure_smu(pspa, source_ch, mode='VOLT', compliance=compliance)
    
    # Ground source
    set_voltage(pspa, source_ch, 0)
    output_on(pspa, source_ch)
    
    # Prepare data storage
    vds_array = []
    vgs_array = []
    id_array = []
    
    # Generate sweep values
    vgs_values = np.arange(vgs_start, vgs_stop + vgs_step, vgs_step)
    vds_values = np.arange(vds_start, vds_stop + vds_step, vds_step)
    
    print(f"Starting transistor output characteristics sweep...")
    print(f"Vgs: {vgs_start} to {vgs_stop} V in {vgs_step} V steps")
    print(f"Vds: {vds_start} to {vds_stop} V in {vds_step} V steps")
    
    # Enable outputs
    output_on(pspa, gate_ch)
    output_on(pspa, drain_ch)
    
    # Sweep through Vgs values
    for vgs in vgs_values:
        set_voltage(pspa, gate_ch, vgs)
        
        # Sweep through Vds values
        for vds in vds_values:
            set_voltage(pspa, drain_ch, vds)
            
            # Wait for settling
            pspa.write(":TRIG:ACQ")
            
            # Measure drain current
            id_val = float(pspa.query(f":MEAS:CURR? (@{drain_ch})"))
            
            vds_array.append(vds)
            vgs_array.append(vgs)
            id_array.append(id_val)
    
    # Turn off outputs
    output_off(pspa, drain_ch)
    output_off(pspa, gate_ch)
    output_off(pspa, source_ch)
    
    print("Output characteristics measurement complete")
    
    return {
        'Vds': np.array(vds_array),
        'Vgs': np.array(vgs_array),
        'Id': np.array(id_array)
    }


def measure_transistor_transfer_characteristics(pspa, vgs_start, vgs_stop, vgs_step,
                                                 vds_constant,
                                                 drain_ch=1, gate_ch=2, source_ch=3,
                                                 compliance=0.1):
    """
    Measure transistor transfer characteristics (Id vs Vgs at constant Vds).
    
    Args:
        pspa (pyvisa.Resource): PSPA instrument object
        vgs_start (float): Start gate-source voltage (V)
        vgs_stop (float): Stop gate-source voltage (V)
        vgs_step (float): Gate-source voltage step (V)
        vds_constant (float): Constant drain-source voltage (V)
        drain_ch (int): Drain channel number
        gate_ch (int): Gate channel number
        source_ch (int): Source channel number
        compliance (float): Current compliance (A)
    
    Returns:
        dict: Dictionary containing Vgs and Id arrays
    """
    # Configure channels
    configure_smu(pspa, drain_ch, mode='VOLT', compliance=compliance)
    configure_smu(pspa, gate_ch, mode='VOLT', compliance=compliance)
    configure_smu(pspa, source_ch, mode='VOLT', compliance=compliance)
    
    # Ground source
    set_voltage(pspa, source_ch, 0)
    output_on(pspa, source_ch)
    
    # Set constant Vds
    set_voltage(pspa, drain_ch, vds_constant)
    
    # Prepare data storage
    vgs_array = []
    id_array = []
    ig_array = []
    
    # Generate sweep values
    vgs_values = np.arange(vgs_start, vgs_stop + vgs_step, vgs_step)
    
    print(f"Starting transistor transfer characteristics sweep...")
    print(f"Vds = {vds_constant} V (constant)")
    print(f"Vgs: {vgs_start} to {vgs_stop} V in {vgs_step} V steps")
    
    # Enable outputs
    output_on(pspa, gate_ch)
    output_on(pspa, drain_ch)
    
    # Sweep through Vgs values
    for vgs in vgs_values:
        set_voltage(pspa, gate_ch, vgs)
        
        # Wait for settling
        pspa.write(":TRIG:ACQ")
        
        # Measure drain current and gate current
        id_val = float(pspa.query(f":MEAS:CURR? (@{drain_ch})"))
        ig_val = float(pspa.query(f":MEAS:CURR? (@{gate_ch})"))
        
        vgs_array.append(vgs)
        id_array.append(id_val)
        ig_array.append(ig_val)
    
    # Turn off outputs
    output_off(pspa, drain_ch)
    output_off(pspa, gate_ch)
    output_off(pspa, source_ch)
    
    print("Transfer characteristics measurement complete")
    
    return {
        'Vgs': np.array(vgs_array),
        'Id': np.array(id_array),
        'Ig': np.array(ig_array)
    }


# =============================
# General I-V Measurements
# =============================

def measure_iv_curve(pspa, v_start, v_stop, v_step, channel=1, compliance=0.1):
    """
    Measure simple I-V curve for a two-terminal device.
    
    Args:
        pspa (pyvisa.Resource): PSPA instrument object
        v_start (float): Start voltage (V)
        v_stop (float): Stop voltage (V)
        v_step (float): Voltage step (V)
        channel (int): Channel number
        compliance (float): Current compliance (A)
    
    Returns:
        dict: Dictionary containing voltage and current arrays
    """
    # Configure channel
    configure_smu(pspa, channel, mode='VOLT', compliance=compliance)
    
    # Prepare data storage
    voltage_array = []
    current_array = []
    
    # Generate sweep values
    voltages = np.arange(v_start, v_stop + v_step, v_step)
    
    print(f"Starting I-V sweep on channel {channel}...")
    print(f"Voltage: {v_start} to {v_stop} V in {v_step} V steps")
    
    # Enable output
    output_on(pspa, channel)
    
    # Sweep through voltage values
    for v in voltages:
        set_voltage(pspa, channel, v)
        
        # Wait for settling
        pspa.write(":TRIG:ACQ")
        
        # Measure current
        i_val = float(pspa.query(f":MEAS:CURR? (@{channel})"))
        
        voltage_array.append(v)
        current_array.append(i_val)
    
    # Turn off output
    output_off(pspa, channel)
    
    print("I-V measurement complete")
    
    return {
        'Voltage': np.array(voltage_array),
        'Current': np.array(current_array)
    }


def measure_iv_bidirectional(pspa, v_max, v_step, channel=1, compliance=0.1):
    """
    Measure bidirectional I-V curve (0 -> +V -> 0 -> -V -> 0).
    
    Args:
        pspa (pyvisa.Resource): PSPA instrument object
        v_max (float): Maximum voltage magnitude (V)
        v_step (float): Voltage step (V)
        channel (int): Channel number
        compliance (float): Current compliance (A)
    
    Returns:
        dict: Dictionary containing voltage and current arrays
    """
    # Configure channel
    configure_smu(pspa, channel, mode='VOLT', compliance=compliance)
    
    # Prepare data storage
    voltage_array = []
    current_array = []
    
    # Generate sweep values: 0 -> +V -> 0 -> -V -> 0
    sweep = np.concatenate([
        np.arange(0, v_max + v_step, v_step),
        np.arange(v_max - v_step, -v_max - v_step, -v_step),
        np.arange(-v_max + v_step, 0 + v_step, v_step)
    ])
    
    print(f"Starting bidirectional I-V sweep on channel {channel}...")
    print(f"Voltage: 0 -> {v_max} -> 0 -> {-v_max} -> 0 V in {v_step} V steps")
    
    # Enable output
    output_on(pspa, channel)
    
    # Sweep through voltage values
    for v in sweep:
        set_voltage(pspa, channel, v)
        
        # Wait for settling
        pspa.write(":TRIG:ACQ")
        
        # Measure current
        i_val = float(pspa.query(f":MEAS:CURR? (@{channel})"))
        
        voltage_array.append(v)
        current_array.append(i_val)
    
    # Turn off output
    output_off(pspa, channel)
    
    print("Bidirectional I-V measurement complete")
    
    return {
        'Voltage': np.array(voltage_array),
        'Current': np.array(current_array)
    }


# =============================
# Pulse Measurements
# =============================

def configure_pulse_mode(pspa, channel, pulse_width, pulse_period):
    """
    Configure a channel for pulse mode operation.
    
    Args:
        pspa (pyvisa.Resource): PSPA instrument object
        channel (int): Channel number
        pulse_width (float): Pulse width in seconds
        pulse_period (float): Pulse period in seconds
    """
    pspa.write(f":SOUR{channel}:FUNC:MODE PULS")
    pspa.write(f":SOUR{channel}:PULS:WIDT {pulse_width}")
    pspa.write(f":SOUR{channel}:PULS:PER {pulse_period}")
    print(f"Channel {channel} configured for pulse mode: {pulse_width*1e6} µs width, {pulse_period*1e3} ms period")


def measure_pulsed_iv(pspa, v_base, v_pulse, pulse_width, pulse_period, num_pulses,
                      channel=1, compliance=0.1):
    """
    Measure current response to voltage pulses.
    
    Args:
        pspa (pyvisa.Resource): PSPA instrument object
        v_base (float): Base voltage level (V)
        v_pulse (float): Pulse voltage level (V)
        pulse_width (float): Pulse width in seconds
        pulse_period (float): Pulse period in seconds
        num_pulses (int): Number of pulses to apply
        channel (int): Channel number
        compliance (float): Current compliance (A)
    
    Returns:
        dict: Dictionary containing time, voltage, and current arrays
    """
    # Configure channel for pulsed operation
    configure_smu(pspa, channel, mode='VOLT', compliance=compliance)
    configure_pulse_mode(pspa, channel, pulse_width, pulse_period)
    
    # Set base and pulse levels
    pspa.write(f":SOUR{channel}:VOLT:LEV:IMM:AMPL {v_pulse}")
    pspa.write(f":SOUR{channel}:VOLT:LEV:IMM:OFFS {v_base}")
    
    # Prepare data storage
    time_array = []
    voltage_array = []
    current_array = []
    
    print(f"Starting pulsed I-V measurement on channel {channel}...")
    print(f"Base: {v_base} V, Pulse: {v_pulse} V")
    print(f"Pulse width: {pulse_width*1e6} µs, Period: {pulse_period*1e3} ms")
    print(f"Number of pulses: {num_pulses}")
    
    # Enable output
    output_on(pspa, channel)
    
    # Trigger pulse measurements
    for pulse_num in range(num_pulses):
        pspa.write(f":TRIG:SOUR{channel}:SING")
        
        # Measure at base and pulse levels
        # Base measurement
        i_base = float(pspa.query(f":MEAS:CURR? (@{channel})"))
        time_array.append(pulse_num * pulse_period)
        voltage_array.append(v_base)
        current_array.append(i_base)
        
        # Pulse measurement (during pulse)
        pspa.write(":TRIG:ACQ")
        i_pulse = float(pspa.query(f":MEAS:CURR? (@{channel})"))
        time_array.append(pulse_num * pulse_period + pulse_width/2)
        voltage_array.append(v_pulse)
        current_array.append(i_pulse)
    
    # Turn off output
    output_off(pspa, channel)
    
    print("Pulsed I-V measurement complete")
    
    return {
        'Time': np.array(time_array),
        'Voltage': np.array(voltage_array),
        'Current': np.array(current_array)
    }


def measure_pulsed_transistor(pspa, vds_pulse, vgs_pulse, vds_base, vgs_base,
                               pulse_width, pulse_period, num_pulses,
                               drain_ch=1, gate_ch=2, source_ch=3, compliance=0.1):
    """
    Measure pulsed transistor characteristics with synchronized gate and drain pulses.
    
    Args:
        pspa (pyvisa.Resource): PSPA instrument object
        vds_pulse (float): Drain-source pulse voltage (V)
        vgs_pulse (float): Gate-source pulse voltage (V)
        vds_base (float): Drain-source base voltage (V)
        vgs_base (float): Gate-source base voltage (V)
        pulse_width (float): Pulse width in seconds
        pulse_period (float): Pulse period in seconds
        num_pulses (int): Number of pulses to apply
        drain_ch (int): Drain channel number
        gate_ch (int): Gate channel number
        source_ch (int): Source channel number
        compliance (float): Current compliance (A)
    
    Returns:
        dict: Dictionary containing time, Vds, Vgs, Id, and Ig arrays
    """
    # Configure channels
    configure_smu(pspa, drain_ch, mode='VOLT', compliance=compliance)
    configure_smu(pspa, gate_ch, mode='VOLT', compliance=compliance)
    configure_smu(pspa, source_ch, mode='VOLT', compliance=compliance)
    
    # Configure pulse mode for drain and gate
    configure_pulse_mode(pspa, drain_ch, pulse_width, pulse_period)
    configure_pulse_mode(pspa, gate_ch, pulse_width, pulse_period)
    
    # Set pulse and base levels
    pspa.write(f":SOUR{drain_ch}:VOLT:LEV:IMM:AMPL {vds_pulse}")
    pspa.write(f":SOUR{drain_ch}:VOLT:LEV:IMM:OFFS {vds_base}")
    pspa.write(f":SOUR{gate_ch}:VOLT:LEV:IMM:AMPL {vgs_pulse}")
    pspa.write(f":SOUR{gate_ch}:VOLT:LEV:IMM:OFFS {vgs_base}")
    
    # Ground source
    set_voltage(pspa, source_ch, 0)
    
    # Prepare data storage
    time_array = []
    vds_array = []
    vgs_array = []
    id_array = []
    ig_array = []
    
    print(f"Starting pulsed transistor measurement...")
    print(f"Vds: {vds_base} V (base) -> {vds_pulse} V (pulse)")
    print(f"Vgs: {vgs_base} V (base) -> {vgs_pulse} V (pulse)")
    print(f"Pulse width: {pulse_width*1e6} µs, Period: {pulse_period*1e3} ms")
    print(f"Number of pulses: {num_pulses}")
    
    # Enable outputs
    output_on(pspa, source_ch)
    output_on(pspa, gate_ch)
    output_on(pspa, drain_ch)
    
    # Synchronize triggers for simultaneous pulsing
    pspa.write(":TRIG:SOUR:SYNC ON")
    
    # Trigger pulse measurements
    for pulse_num in range(num_pulses):
        # Trigger synchronized pulses
        pspa.write(":TRIG:ALL")
        
        # Measure during pulse
        pspa.write(":TRIG:ACQ")
        id_pulse = float(pspa.query(f":MEAS:CURR? (@{drain_ch})"))
        ig_pulse = float(pspa.query(f":MEAS:CURR? (@{gate_ch})"))
        
        time_array.append(pulse_num * pulse_period + pulse_width/2)
        vds_array.append(vds_pulse)
        vgs_array.append(vgs_pulse)
        id_array.append(id_pulse)
        ig_array.append(ig_pulse)
    
    # Turn off outputs
    output_off(pspa, drain_ch)
    output_off(pspa, gate_ch)
    output_off(pspa, source_ch)
    
    print("Pulsed transistor measurement complete")
    
    return {
        'Time': np.array(time_array),
        'Vds': np.array(vds_array),
        'Vgs': np.array(vgs_array),
        'Id': np.array(id_array),
        'Ig': np.array(ig_array)
    }
