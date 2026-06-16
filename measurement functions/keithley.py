"""
Filename: keithley.py
Author: Ethan Ruddell
Date: 2026-05-28
Description: Remote measurement functions for Keithley 2400 and 4156C multi-probe transistor characterization.

This module implements 3-probe transistor measurements with:
  - 4156C SMU2: Source probe (establishes ground reference, always 0V)
  - 4156C SMU1: Drain probe (applies Vds, measures Id)
  - 2400 (Keithley): Gate probe (steps through Vg values)

CRITICAL SEQUENCE RULES:
  Enable order:  SMU2 (source) → SMU1 (drain) → 2400 (gate)
                 This ensures proper biasing with ground reference first
  
  Disable order: 2400 (gate) → SMU1 (drain) → SMU2 (source)
                 Reverse of enable to avoid floating junctions
"""

import atexit
import logging
import time

import numpy as np
import pyvisa

import pspa
from utility.gpib_utils import create_visa_resource_manager

log = logging.getLogger(__name__)

DEFAULT_KEITHLEY_ADDRESS = "GPIB0::17::INSTR"
DEFAULT_4156C_ADDRESS = "GPIB0::16::INSTR"
# Backward-compatible alias used by the GUI/CLI call sites.
DEFAULT_GPIB_ADDRESS = DEFAULT_KEITHLEY_ADDRESS
_OPEN_KEITHLEY_SESSIONS = set()

# =====================================================================
# CONNECTION FUNCTIONS
# =====================================================================

def connect_4156c(resource_name=None):
    """Open a VISA session to the 4156C Semiconductor Parameter Analyzer."""
    target = str(resource_name or DEFAULT_4156C_ADDRESS).strip()
    if not target:
        raise ValueError("4156C GPIB address not configured")
    
    rm = pyvisa.ResourceManager()
    print(f"Connecting to 4156C at: {target}")
    inst = rm.open_resource(target)
    inst.timeout = pspa.DEFAULT_TIMEOUT_MS
    return inst


def connect_keithley_2400(resource_name=None):
    """Open a VISA session to the Keithley 2400 SourceMeter."""
    target = str(resource_name or DEFAULT_KEITHLEY_ADDRESS).strip()
    if not target:
        raise ValueError("Keithley 2400 GPIB address not configured")
    
    rm = pyvisa.ResourceManager()
    print(f"Connecting to Keithley 2400 at: {target}")
    keithley = rm.open_resource(target)
    keithley.timeout = pspa.DEFAULT_TIMEOUT_MS
    _OPEN_KEITHLEY_SESSIONS.add(keithley)
    return keithley


def connect_keithley_sourcemeter(resource_name=None):
    """Backward-compatible alias for the Keithley 2400 connection helper."""
    return connect_keithley_2400(resource_name)


def disconnect_keithley_sourcemeter(keithley):
    """Backward-compatible alias for the Keithley 2400 disconnect helper."""
    return disconnect_instrument(keithley, "Keithley 2400")


def configure_voltage_sweep_source(keithley, compliance=0.1):
    """Configure the Keithley 2400 to source voltage for a gate sweep."""
    # Match the known-good command sequence from LEGACY/Need ReRAM Mark 4.py.
    keithley.write("*RST")
    keithley.write("*CLS")
    keithley.write("ROUT:TERM REAR")
    keithley.write("SOUR:FUNC VOLT")
    keithley.write("SOUR:VOLT:MODE FIXED")
    keithley.write("SENS:FUNC 'CURR:DC'")
    keithley.write("SENS:CURR:RANG:AUTO ON")
    keithley.write("FORM:ELEM CURR")
    keithley.write(f"SOUR:VOLT:ILIM {compliance}")
    keithley.write("SOUR:VOLT 0")
    keithley.write("OUTP ON")


def configure_current_readback(keithley, source_voltage=0.0, compliance=0.1):
    """Compatibility wrapper for the older Keithley current readback setup."""
    configure_voltage_sweep_source(keithley, compliance=compliance)


def disconnect_instrument(inst, name="Instrument"):
    """Safely disable outputs and close session."""
    try:
        # Use PSPA helper for the 4156C to ensure proper FLEX shutdown
        if name and "4156" in str(name):
            try:
                pspa.disconnect_pspa(inst)
                return
            except Exception as e:
                log.warning(f"pspa.disconnect_pspa failed for {name}: {e}")

        # Default behavior for SCPI instruments (e.g., Keithley 2400)
        try:
            inst.write(":OUTP OFF")
        except Exception as e:
            log.warning(f"Error disabling {name} output: {e}")
    except Exception:
        # fall through to close attempt
        pass
    try:
        inst.close()
    except Exception as e:
        log.warning(f"Error closing {name} session: {e}")
    finally:
        _OPEN_KEITHLEY_SESSIONS.discard(inst)


def _close_open_keithley_sessions():
    """Best-effort cleanup so Keithley sessions are closed on process exit."""
    for inst in list(_OPEN_KEITHLEY_SESSIONS):
        try:
            disconnect_instrument(inst, "Keithley 2400")
        except Exception:
            pass


atexit.register(_close_open_keithley_sessions)


# =====================================================================
# INITIALIZATION & CONFIGURATION FUNCTIONS
# =====================================================================

def initialize_4156c(inst_4156c, vds_target=5.0, i_compliance=0.1):
    """
    Initialize 4156C to a known state and pre-configure both SMU channels.
    
    SMU2 (source): 0V, compliance set
    SMU1 (drain):  target Vds, compliance set
    
    Outputs remain OFF until enable sequence.
    """
    print("\n[4156C] Initializing using pspa helpers...")
    # Use the PSPA module to configure SMU channels correctly and safely.
    # Configure SMU2 (source) to 0V and set compliance
    pspa.configure_smu(inst_4156c, 2, mode='VOLT', compliance=i_compliance)
    pspa.set_voltage(inst_4156c, 2, 0.0)

    # Configure SMU1 (drain) to target Vds and set compliance
    pspa.configure_smu(inst_4156c, 1, mode='VOLT', compliance=i_compliance)
    pspa.set_voltage(inst_4156c, 1, vds_target)
    
    log.info("4156C initialization complete (via pspa).")


def initialize_keithley_2400(inst_keithley, vg_start=0.0, i_compliance=0.1):
    """
    Initialize Keithley 2400 (gate probe) with SCPI commands.
    
    Configured for voltage source mode at starting Vg value.
    Output remains OFF until enable sequence.
    """
    print(f"\n[Keithley] Initializing...")
    
    # Reset to known state
    inst_keithley.write("*RST")
    #inst_keithley.write("*CLS")
    
    # Configure for voltage source, current measurement
    print(f"  Configuring gate sweep: start {vg_start}V, compliance {i_compliance} A")
    inst_keithley.write(":ROUT:TERM REAR")
    inst_keithley.write(":SOUR:FUNC VOLT")
    inst_keithley.write(":SOUR:VOLT:MODE FIXED")
    inst_keithley.write(":SENS:FUNC 'CURR:DC'")
    inst_keithley.write(":SENS:CURR:RANG:AUTO ON")
    inst_keithley.write(":FORM:ELEM CURR")
    inst_keithley.write(f":SOUR:VOLT:ILIM {i_compliance}")
    inst_keithley.write(f":SOUR:VOLT {vg_start}")
    # DO NOT enable output yet (no :OUTP ON) - wait for enable sequence
    
    log.info("Keithley 2400 initialization complete. Output is OFF.")


# =====================================================================
# ENABLE SEQUENCE (CRITICAL ORDER)
# =====================================================================

def enable_4156c_smu2(inst_4156c):
    """
    Enable 4156C SMU2 (source probe) FIRST.
    
    This establishes the 0V ground reference at the device source terminal.
    Must be done before enabling SMU1.
    """
    print(f"\n[Enable] Step 1: Enabling 4156C SMU2 (source at 0V)")
    pspa.output_on(inst_4156c, 2)
    time.sleep(0.1)
    log.info("SMU2 output enabled - ground reference established")


def enable_4156c_smu1(inst_4156c, settle_time_s=1.0):
    """
    Enable 4156C SMU1 (drain probe) SECOND.
    
    With source reference already present, now apply Vds to drain.
    """
    print(f"\n[Enable] Step 2: Enabling 4156C SMU1 (drain applies Vds)")
    pspa.output_on(inst_4156c, 1)
    time.sleep(settle_time_s)
    log.info(f"SMU1 output enabled, settling for {settle_time_s}s")


def enable_keithley_2400(inst_keithley):
    """
    Enable Keithley 2400 (gate probe) LAST.
    
    Gate bias is always applied after drain/source are stable and referenced.
    """
    print(f"\n[Enable] Step 3: Enabling Keithley 2400 (gate probe)")
    inst_keithley.write(":OUTP ON")
    time.sleep(0.1)
    log.info("Gate output enabled - measurement ready")


def enable_all_3probe(inst_4156c, inst_keithley, vds_settle_s=1.0):
    """
    Execute the complete enable sequence in the correct order.
    
    Order: SMU2 (source) → SMU1 (drain) → 2400 (gate)
    """
    print("\n" + "="*70)
    print("BEGINNING ENABLE SEQUENCE (DO NOT INTERRUPT)")
    print("="*70)
    
    enable_4156c_smu2(inst_4156c)
    enable_4156c_smu1(inst_4156c, settle_time_s=vds_settle_s)
    enable_keithley_2400(inst_keithley)
    
    print("\n" + "="*70)
    print("ENABLE SEQUENCE COMPLETE - Device is now biased")
    print("="*70)


# =====================================================================
# MEASUREMENT LOOP
# =====================================================================

def verify_gate_voltage_settled(inst_keithley, target_vg, tolerance_v=0.01, max_tries=10):
    """
    Verify that gate voltage has settled to commanded value.
    
    Parameters:
    -----------
    inst_keithley : PyVISA instrument
        Keithley 2400 instrument handle
    target_vg : float
        Commanded gate voltage (V)
    tolerance_v : float
        Acceptable voltage error (V). Default 0.01V (10mV)
    max_tries : int
        Maximum readback attempts before timeout
    
    Returns:
    --------
    bool : True if voltage settled within tolerance, False if timeout
    """
    for attempt in range(max_tries):
        try:
            actual_vg = float(inst_keithley.query(":SOUR:VOLT?"))
            error_v = abs(actual_vg - target_vg)
            
            if error_v <= tolerance_v:
                if attempt > 0:
                    print(f"    [Synced after {attempt} readbacks] Vg={actual_vg:.4f}V (δ={error_v:.4f}V)", end="")
                return True
            
            # Voltage not settled yet, wait and retry
            time.sleep(0.01)
        except Exception as e:
            log.warning(f"Voltage readback attempt {attempt+1} failed: {e}")
            time.sleep(0.01)
    
    actual_vg = float(inst_keithley.query(":SOUR:VOLT?"))
    error_v = abs(actual_vg - target_vg)
    log.warning(f"Gate voltage did not settle: target={target_vg}V, actual={actual_vg}V, error={error_v}V")
    return False


def measure_transistor_id(inst_4156c):
    """
    Trigger a measurement on 4156C SMU1 (drain current) with status validation.
    Returns the measured drain current in Amps.
    """
    # Use pspa helper which wraps the correct TI? channel query
    return pspa.measure_channel_current(inst_4156c, 1)


def measure_voltage_sweep_current(pspa_inst, keithley_address, vgs_start, vgs_stop,
                                  vgs_step, vds_constant, drain_channel=1,
                                  source_channel=2, compliance=0.1, settle_s=0.05,
                                  integration_time="MED"):
    """
    Sweep Vgs with the Keithley while PSPA holds constant Vds and measures Id.

    This preserves the GUI/CLI API used by the merged measurement flow.
    """
    drain_channel = pspa.validate_channel(drain_channel)
    source_channel = pspa.validate_channel(source_channel)
    pspa.validate_compliance(compliance)
    pspa.validate_step(vgs_step)

    vgs_forward = pspa.sweep_values(vgs_start, vgs_stop, vgs_step)
    vgs_reverse = pspa.sweep_values(vgs_stop, vgs_start, -vgs_step)
    vgs_array = []
    id_array = []
    direction_array = []
    keithley = None

    try:
        # Mirror the normal PSPA transfer setup exactly for drain/source biasing.
        drain_channel, _, source_channel = pspa.configure_transfer_bias(
            pspa_inst,
            drain_channel,
            source_channel,
            vds_constant,
            compliance=compliance,
            integration_time=integration_time,
        )
        try:
            drain_readback = pspa.measure_channel_voltage(pspa_inst, drain_channel)
            source_readback = pspa.measure_channel_voltage(pspa_inst, source_channel)
            print(f"PSPA commanded: drain CH{drain_channel} = {vds_constant} V")
            print(f"PSPA readback drain CH{drain_channel}: {drain_readback}")
            print(f"PSPA readback source CH{source_channel}: {source_readback}")
            print(f"PSPA errors: {pspa.check_errors(pspa_inst)}")
            log.info(
                "PSPA readback before sweep: drain CH%s=%.6f V, source CH%s=%.6f V",
                drain_channel,
                drain_readback,
                source_channel,
                source_readback,
            )
        except Exception as exc:
            print(f"PSPA readback failed before sweep: {exc}")
            log.warning("PSPA voltage readback failed before sweep: %s", exc)

        # Match the legacy Keithley sweep flow exactly.
        keithley = connect_keithley_sourcemeter(keithley_address)
        configure_voltage_sweep_source(keithley, compliance=compliance)
        keithley.write(f"SOUR:VOLT:RANG {max(abs(vgs_start), abs(vgs_stop))}")
        keithley.write("OUTP ON")

        def _run_sweep(values, direction):
            for vgs in values:
                keithley.write(f"SOUR:VOLT {vgs}")
                time.sleep(settle_s)
                # Legacy pattern triggers Keithley READ? at each point.
                try:
                    keithley.query("READ?")
                except Exception:
                    pass
                #pspa.set_voltage(pspa_inst, drain_channel, vds_constant);
                current_value = pspa.measure_channel_current(pspa_inst, drain_channel)
                vgs_array.append(vgs)
                id_array.append(current_value)
                direction_array.append(direction)

        _run_sweep(vgs_forward, "forward")
        _run_sweep(vgs_reverse, "reverse")

    finally:
        try:
            if keithley is not None:
                keithley.write("OUTP OFF")
        except Exception:
            pass
        if keithley is not None:
            disconnect_keithley_sourcemeter(keithley)
        try:
            pspa.set_voltage(pspa_inst, drain_channel, 0)
            pspa.set_voltage(pspa_inst, source_channel, 0)
            pspa_inst.write("CL")
        except Exception:
            pass

    vgs_np = np.array(vgs_array)
    id_np = np.array(id_array)
    dir_np = np.array(direction_array)

    if id_np.size > 0 and np.nanmedian(id_np) < 0:
        id_np = -id_np
        log.info("Keithley current polarity normalized to positive convention")

    return {
        'Vgs': vgs_np,
        'Id': id_np,
        'Sweep_Direction': dir_np,
    }


# =====================================================================
# SHUTDOWN SEQUENCE (REVERSE OF ENABLE - CRITICAL ORDER)
# =====================================================================

def disable_keithley_2400(inst_keithley, ramp_to_zero=True, ramp_steps=10):
    """
    Disable Keithley 2400 (gate probe) FIRST in shutdown.
    
    Optionally ramp Vg back to 0V before disabling to avoid abrupt transitions.
    """
    print(f"\n[Shutdown] Step 1: Disabling Keithley 2400 (gate probe)")
    
    if ramp_to_zero:
        print(f"  Ramping gate voltage to 0V across {ramp_steps} steps...")
        current_vg = inst_keithley.query(":SOUR:VOLT?")
        current_vg = float(current_vg)
        
        for step_idx in range(1, ramp_steps + 1):
            fraction = 1.0 - (step_idx / ramp_steps)
            target_v = current_vg * fraction
            inst_keithley.write(f":SOUR:VOLT {target_v}")
            time.sleep(0.05)
    
    inst_keithley.write(":OUTP OFF")
    time.sleep(0.1)
    log.info("Gate output disabled")


def disable_4156c_smu1(inst_4156c):
    """
    Disable 4156C SMU1 (drain probe) SECOND in shutdown.
    
    Removes drain bias while source reference is still active.
    """
    print(f"\n[Shutdown] Step 2: Disabling 4156C SMU1 (drain)")
    pspa.output_off(inst_4156c, 1)
    time.sleep(0.1)
    log.info("SMU1 output disabled")


def disable_4156c_smu2(inst_4156c):
    """
    Disable 4156C SMU2 (source probe) LAST in shutdown.
    
    Remove the ground reference last to avoid leaving device in undefined state.
    """
    print(f"\n[Shutdown] Step 3: Disabling 4156C SMU2 (source)")
    pspa.output_off(inst_4156c, 2)
    time.sleep(0.1)
    log.info("SMU2 output disabled - ground reference removed")


def disable_all_3probe(inst_keithley, inst_4156c):
    """
    Execute the complete disable sequence in the REVERSE order of enable.
    
    Order: 2400 (gate) → SMU1 (drain) → SMU2 (source)
    """
    print("\n" + "="*70)
    print("BEGINNING SHUTDOWN SEQUENCE (REVERSE OF ENABLE)")
    print("="*70)
    
    disable_keithley_2400(inst_keithley, ramp_to_zero=True, ramp_steps=10)
    disable_4156c_smu1(inst_4156c)
    disable_4156c_smu2(inst_4156c)
    
    print("\n" + "="*70)
    print("SHUTDOWN SEQUENCE COMPLETE - Device is safe")
    print("="*70)


# =====================================================================
# MAIN 3-PROBE MEASUREMENT FUNCTION
# =====================================================================

def measure_transistor_vg_sweep(
    gpib_4156c=None,
    gpib_keithley=None,
    vds=5.0,
    vg_start=0.0,
    vg_stop=5.0,
    vg_step=0.1,
    i_compliance=0.1,
    vds_settle_s=0.05,
    vg_settle_s=0.05
):
    """
    Perform a complete 3-probe transistor measurement with Vg sweep.
    
    Parameters:
    -----------
    gpib_4156c : str, optional
        GPIB address of 4156C (default: "GPIB0::16::INSTR")
    gpib_keithley : str, optional
        GPIB address of Keithley 2400 (default: "GPIB0::17::INSTR")
    vds : float
        Drain-source voltage to apply (V)
    vg_start : float
        Starting gate voltage (V)
    vg_stop : float
        Ending gate voltage (V)
    vg_step : float
        Gate voltage step size (V)
    i_compliance : float
        Current compliance limit for all sources (A)
    vds_settle_s : float
        Settling time after enabling SMU1 (drain) (s)
    vg_settle_s : float
        Settling time after each gate step (s)
    
    Returns:
    --------
    dict with keys:
        'vg': array of gate voltages (V)
        'id': array of corresponding drain currents (A)
        'status': 'success' or 'failed'
    """
    
    inst_4156c = None
    inst_keithley = None
    result = {'vg': np.array([]), 'id': np.array([]), 'status': 'failed'}
    
    try:
        # Connect to instruments
        print("\n" + "="*70)
        print("CONNECTING TO INSTRUMENTS")
        print("="*70)
        inst_4156c = connect_4156c(gpib_4156c)
        inst_keithley = connect_keithley_2400(gpib_keithley)
        
        # Initialize both instruments
        print("\n" + "="*70)
        print("INITIALIZATION PHASE")
        print("="*70)
        initialize_4156c(inst_4156c, vds_target=vds, i_compliance=i_compliance)
        initialize_keithley_2400(inst_keithley, vg_start=vg_start, i_compliance=i_compliance)
        
        # Enable in correct order
        enable_all_3probe(inst_4156c, inst_keithley, vds_settle_s=vds_settle_s)
        
        # Perform measurement sweep
        print("\n" + "="*70)
        print("MEASUREMENT PHASE")
        print("="*70)
        
        vg_values = np.arange(vg_start, vg_stop + vg_step/2, vg_step)
        id_values = []
        
        for idx, vg in enumerate(vg_values, 1):
            print(f"\n[Measurement {idx}/{len(vg_values)}] Gate: {vg:.3f}V", end=" ")
            
            # Step gate to next voltage
            inst_keithley.write(f":SOUR:VOLT {vg}")
            
            # Wait for voltage to settle and verify
            if not verify_gate_voltage_settled(inst_keithley, vg, tolerance_v=0.01, max_tries=10):
                print(f"WARNING: Gate voltage did not settle at {vg}V")
            
            # Perform integration time delay as specified
            time.sleep(vg_settle_s)
            
            # Measure drain current
            i_d = measure_transistor_id(inst_4156c)
            id_values.append(i_d)
            
            print(f"→ Id: {i_d:.6e}A")
        
        # Shutdown in reverse order
        disable_all_3probe(inst_keithley, inst_4156c)
        
        # Package results
        result['vg'] = np.array(vg_values)
        result['id'] = np.array(id_values)
        result['status'] = 'success'
        
        print("\n" + "="*70)
        print("MEASUREMENT COMPLETE")
        print(f"  Vg points: {len(vg_values)}")
        print(f"  Id range: {np.min(id_values):.3e} to {np.max(id_values):.3e} A")
        print("="*70 + "\n")
        
    except Exception as e:
        log.error(f"Measurement failed: {e}", exc_info=True)
        print(f"\nERROR: {e}")
        
        # Emergency shutdown
        try:
            if inst_keithley:
                print("Emergency shutdown: Disabling gate...")
                inst_keithley.write(":OUTP OFF")
            if inst_4156c:
                print("Emergency shutdown: Disabling 4156C channels via pspa.disconnect_pspa...")
                try:
                    pspa.disconnect_pspa(inst_4156c)
                except Exception as e3:
                    log.warning(f"pspa.disconnect_pspa failed during emergency shutdown: {e3}")
        except Exception as e2:
            log.error(f"Emergency shutdown failed: {e2}")
    
    finally:
        # Close connections
        if inst_4156c:
            disconnect_instrument(inst_4156c, "4156C")
        if inst_keithley:
            disconnect_instrument(inst_keithley, "Keithley 2400")
    
    return result