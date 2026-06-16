"""
Filename: keithley.py
Author: Ethan Ruddell
Date: 2026-05-25
Description: Driver helpers for Keithley sourcemeter-based measurements.

This module contains the Keithley-specific configuration and measurement
functions used by the PSPA menu item for transistor-style voltage sweeps
with current readback from a separate Keithley instrument.
"""

import logging
import time

import numpy as np
import pyvisa

import pspa
from utility.gpib_utils import create_visa_resource_manager

log = logging.getLogger(__name__)

DEFAULT_GPIB_ADDRESS = "GPIB0::17::INSTR"


def connect_keithley_sourcemeter(resource_name=None):
    """Open a VISA session to the Keithley sourcemeter."""
    target = str(resource_name or DEFAULT_GPIB_ADDRESS).strip()
    if not target:
        raise ValueError("Keithley GPIB address not configured")

    rm = create_visa_resource_manager()
    print(f"Attempting to connect Keithley sourcemeter to: {target}")
    keithley = rm.open_resource(target)
    keithley.timeout = pspa.DEFAULT_TIMEOUT_MS
    return keithley


def disconnect_keithley_sourcemeter(keithley):
    """Safely disable and close the Keithley sourcemeter session."""
    try:
        keithley.write(":OUTP OFF")
    except Exception:
        pass
    try:
        keithley.close()
    except Exception:
        pass


def configure_current_readback(keithley, source_voltage=0.0, compliance=0.1):
    """Configure the Keithley to source 0 V and measure current."""
    keithley.write("*RST")
    keithley.write("*CLS")
    keithley.write(":SOUR:FUNC VOLT")
    keithley.write(":SOUR:VOLT:MODE FIXED")
    keithley.write(":SENS:FUNC 'CURR:DC'")
    keithley.write(":SENS:CURR:RANG:AUTO ON")
    keithley.write(":FORM:ELEM CURR")
    keithley.write(f":SOUR:VOLT:ILIM {compliance}")
    keithley.write(f":SOUR:VOLT {source_voltage}")
    keithley.write(":OUTP ON")


def measure_voltage_sweep_current(pspa_inst, keithley_address, v_start, v_stop,
                                  v_step, sweep_channel=1, compliance=0.1,
                                  settle_s=0.05):
    """
    Sweep voltage on a PSPA SMU and measure current with an external Keithley.

    The Keithley stays at 0 V and reports current directly, so no switchbox
    changes are required.
    """
    pspa.validate_channel(sweep_channel)
    pspa.validate_compliance(compliance)
    pspa.validate_step(v_step)

    voltages = pspa.sweep_values(v_start, v_stop, v_step)
    voltage_array = []
    current_array = []
    keithley = None

    pspa_inst.write("US")
    pspa_inst.write("FMT 1")
    pspa_inst.write("CL")
    pspa.configure_smu(pspa_inst, sweep_channel, mode='VOLT', compliance=compliance)
    pspa.output_on(pspa_inst, sweep_channel)

    try:
        keithley = connect_keithley_sourcemeter(keithley_address)
        configure_current_readback(keithley, source_voltage=0.0, compliance=compliance)

        print(f"Starting PSPA/Keithley sweep on CH{sweep_channel} with Keithley at {keithley_address}...")
        for voltage in voltages:
            pspa.set_voltage(pspa_inst, sweep_channel, voltage)
            time.sleep(settle_s)
            current_value = pspa.parse_flex_number(keithley.query("READ?"))
            voltage_array.append(voltage)
            current_array.append(current_value)

    finally:
        try:
            pspa_inst.write("CL")
        except Exception:
            pass
        if keithley is not None:
            disconnect_keithley_sourcemeter(keithley)

    voltage_np = np.array(voltage_array)
    current_np = np.array(current_array)

    if current_np.size > 0 and np.nanmedian(current_np) < 0:
        current_np = -current_np
        print("Note: Keithley current polarity normalized (sign inverted to positive convention).")

    return {
        'Voltage': voltage_np,
        'Current': current_np,
    }