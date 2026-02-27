"""
Filename: cryo.py
Author: Ethan Ruddell (updated by assistant)
Date: 2026-02-16
Description:
    Cryo-Con Model 32 / 32B temperature controller module (GPIB / IEEE-488.2 / SCPI).

Key manual-aligned behaviors:
    - Read temperature via INPUT <chan>:TEMPER? (short: INP <chan>:TEMP?)
    - Set loop setpoint via LOOP <no>:SETPT <temp>
    - Set ramp rate via LOOP <no>:RATE <rate> (Units/min)
    - Loop control type via LOOP <no>:TYPE <type> (PID, RAMPP, OFF, etc.)
    - Do NOT spam *RST (manual warns it may lock up early IEEE-488 systems)
    - Do not rely on text termination chars on GPIB; use read_raw (EOI)
"""


#TODO: heater stuff?

from __future__ import annotations

import time
import logging
import numpy as np
import pyvisa
from pyvisa import errors as visa_errors

from gpib_utils import safe_float_input

log = logging.getLogger(__name__)

# =============================
# User settings and constants
# =============================
GPIB_ADDRESS = "GPIB0::12::INSTR"

# Input channels on Model 32/32B are A and B
DEFAULT_INPUT_CHANNEL = "A"   # used for get_current_temperature(loop=2) mapping

# Temperature limits and defaults (your application constraints)
TEMP_MIN_K = 1.4
TEMP_MAX_K = 325.0
TEMP_DEFAULT_START_K = 100.0
TEMP_DEFAULT_END_K = 200.0

# Ramp rate limits and defaults (your application constraints)
RAMP_MIN_K_PER_MIN = 0.1
RAMP_MAX_K_PER_MIN = 50.0
RAMP_DEFAULT_K_PER_MIN = 5.0

MEAS_INTERVAL_MIN_K = 1.0
MEAS_INTERVAL_DEFAULT_K = 10.0

TEMP_TOLERANCE_K = 0.5
TEMP_STABILITY_TIME_S = 30.0

# VISA / comm tuning
DEFAULT_TIMEOUT_MS = 60000
QUERY_DELAY_S = 0.05
WRITE_RETRIES = 3
QUERY_RETRIES = 2

# Initialization behavior:
# Manual warns against sending reset commands before each query.
# Keep reset disabled by default.
ENABLE_RESET_ON_INIT = False
POST_RESET_WAIT_S = 5.0

current_parameters = {
    "temp_start_k": TEMP_DEFAULT_START_K,
    "temp_end_k": TEMP_DEFAULT_END_K,
    "ramp_rate_k_per_min": RAMP_DEFAULT_K_PER_MIN,
    "meas_interval_k": MEAS_INTERVAL_DEFAULT_K,
    "enable_cryo": False,
}

# =============================
# Low-level VISA helpers
# =============================

def _strip_response_bytes(b: bytes) -> str:
    # responses may include CR/LF; strip whitespace and nulls
    return b.decode(errors="ignore").strip().strip("\x00").strip()

def write_safe(inst, cmd: str, *, delay_s: float = 0.0, retries: int = WRITE_RETRIES) -> None:
    """
    Robust write with timeout retries (common if instrument is briefly busy).
    """
    last_exc = None
    backoff = 0.2
    for _ in range(retries + 1):
        try:
            inst.write(cmd)
            if delay_s:
                time.sleep(delay_s)
            return
        except visa_errors.VisaIOError as e:
            last_exc = e
            # Retry only timeouts
            if e.error_code != visa_errors.VI_ERROR_TMO:
                raise
        except Exception as e:
            last_exc = e
        time.sleep(backoff)
        backoff *= 1.7
    raise RuntimeError(f"write_safe failed for {cmd!r}: {last_exc!r}")

def query_raw(inst, cmd: str, *, delay_s: float = QUERY_DELAY_S, retries: int = QUERY_RETRIES) -> str:
    """
    Query that does NOT depend on read_termination characters (GPIB uses EOI).
    """
    last_exc = None
    for _ in range(retries + 1):
        try:
            inst.write(cmd)
            if delay_s:
                time.sleep(delay_s)
            resp = inst.read_raw()
            txt = _strip_response_bytes(resp)
            if txt != "":
                return txt
        except visa_errors.VisaIOError as e:
            last_exc = e
        except Exception as e:
            last_exc = e
        time.sleep(delay_s)
    raise RuntimeError(f"query_raw failed for {cmd!r}: {last_exc!r}")

def _open_resource(rm: pyvisa.ResourceManager, resource_name: str):
    inst = rm.open_resource(resource_name)
    inst.timeout = DEFAULT_TIMEOUT_MS

    # IMPORTANT for GPIB/EOI: do not rely on termination chars
    inst.read_termination = None
    inst.write_termination = ""  # send EOI at end of write; do not append newline

    # Some backends expose send_end; enable if available
    try:
        inst.send_end = True
    except Exception:
        pass

    # Some backends support query_delay; safe to attempt
    try:
        inst.query_delay = QUERY_DELAY_S
    except Exception:
        pass

    return inst

def _reopen_same_resource(inst):
    rm = getattr(inst, "_cryo_rm", None)
    rn = getattr(inst, "_cryo_resource_name", None)
    if rm is None or rn is None:
        return inst
    try:
        inst.close()
    except Exception:
        pass
    new_inst = _open_resource(rm, rn)
    new_inst._cryo_rm = rm
    new_inst._cryo_resource_name = rn
    return new_inst

# =============================
# Connection and Initialization
# =============================

def connect_cryocon(resource_name: str | None = None):
    """
    Open VISA session to Cryo-Con 32B temperature controller.
    """
    if resource_name is None:
        resource_name = GPIB_ADDRESS

    rm = pyvisa.ResourceManager()
    log.info(f"Attempting to connect to Cryo-Con 32B at: {resource_name}")
    print(f"Attempting to connect to Cryo-Con 32B at: {resource_name}")

    cryo = _open_resource(rm, resource_name)
    cryo._cryo_rm = rm
    cryo._cryo_resource_name = resource_name

    # ID query (supported on many IEEE-488.2 instruments; tolerate failure)
    try:
        idn = query_raw(cryo, "*IDN?")
        log.info(f"Connected to: {idn}")
        print(f"Connected to: {idn}")
    except Exception as e:
        log.warning(f"*IDN? failed/unsupported. Continuing. Error: {e}")
        print("Warning: *IDN? failed/unsupported. Continuing.")

    return cryo

def disconnect_cryocon(cryo):
    """Safely disconnect from Cryo-Con."""
    try:
        # Put loop 1 in a safe state: turn off control
        try:
            set_loop_type(cryo, "OFF", loop=2)
        except Exception:
            pass
        time.sleep(0.2)
        cryo.close()
        log.info("Cryo-Con disconnected")
        print("Cryo-Con disconnected")
    except Exception as e:
        log.warning(f"Error during disconnect: {e}")

def initialize_cryocon(cryo):
    """
    Basic initialization for Cryo-Con 32/32B.

    NOTE: Manual warns not to send reset commands repeatedly on IEEE-488 systems,
    so *RST is disabled by default. If you enable it, we wait longer afterward.
    """
    # Clear status (safe on IEEE-488.2; tolerate timeouts by reopening once)
    try:
        write_safe(cryo, "*CLS", delay_s=0.05, retries=2)
    except visa_errors.VisaIOError as e:
        if e.error_code == visa_errors.VI_ERROR_TMO:
            log.warning("*CLS timed out; reopening VISA session once.")
            cryo = _reopen_same_resource(cryo)
        else:
            raise

    if ENABLE_RESET_ON_INIT:
        try:
            write_safe(cryo, "*RST", delay_s=0.0, retries=1)
            time.sleep(POST_RESET_WAIT_S)
            # Clear after reset
            try:
                write_safe(cryo, "*CLS", delay_s=0.05, retries=1)
            except Exception:
                pass
        except Exception as e:
            log.warning(f"*RST failed/timed out; continuing without reset. {e}")

    log.info("Cryo-Con initialized")
    print("Cryo-Con initialized")
    return cryo

def check_errors(cryo):
    """
    IEEE-488.2 error queue query.
    Many instruments support SYST:ERR?; if unsupported, returns [].
    """
    errors = []
    try:
        err = query_raw(cryo, "SYST:ERR?")
        # Typical "0,No error" or "+0,\"No error\""
        if "no error" not in err.lower() and not err.lstrip().startswith(("+0", "0")):
            errors.append(err.strip())
    except Exception:
        pass
    return errors

# =============================
# Cryo-Con SCPI wrappers
# =============================

def _norm_channel(ch: str) -> str:
    ch = str(ch).strip().upper()
    if ch not in ("A", "B"):
        raise ValueError("Input channel must be 'A' or 'B'")
    return ch

def get_input_temperature(cryo, channel: str = "A") -> float | None:
    """
    Read input temperature in DISPLAY units using:
        INPUT <chan>:TEMPER?   (short form often works: INP <chan>:TEMP?)
    """
    try:
        ch = _norm_channel(channel)
        resp = query_raw(cryo, f"INPUT {ch}:TEMPER?")
        return float(resp)
    except Exception as e:
        log.error(f"Error reading temperature on input {channel}: {e}")
        return None

def get_loop_source(cryo, loop: int = 1) -> str | None:
    """
    Query which input channel controls a loop:
        LOOP <no>:SOURCE?
    Response is typically CHA or CHB.
    """
    try:
        resp = query_raw(cryo, f"LOOP {int(loop)}:SOURCE?")
        resp = resp.strip().upper()
        if "A" in resp:
            return "A"
        if "B" in resp:
            return "B"
        return resp
    except Exception as e:
        log.debug(f"Error reading loop {loop} source: {e}")
        return None

def set_loop_source(cryo, channel: str, loop: int = 1) -> None:
    """
    Set loop controlling input channel:
        LOOP <no>:SOURCE CH<chan>
    """
    ch = _norm_channel(channel)
    write_safe(cryo, f"LOOP {int(loop)}:SOURCE CH{ch}", delay_s=0.02)

def set_loop_type(cryo, loop_type: str, loop: int = 1) -> None:
    """
    Set loop control type:
        LOOP <no>:TYPE <type>
    Common types per manual include: OFF, PID, MAN, TABLE, RAMPP
    """
    lt = str(loop_type).strip()
    write_safe(cryo, f"LOOP {int(loop)}:TYPE {lt}", delay_s=0.02)

def get_loop_type(cryo, loop: int = 1) -> str | None:
    try:
        return query_raw(cryo, f"LOOP {int(loop)}:TYPE?").strip()
    except Exception:
        return None

def get_loop_ramp_active(cryo, loop: int = 1) -> str | None:
    """
    Query ramp status:
        LOOP <no>:RAMP?
    Returns ON/OFF
    """
    try:
        return query_raw(cryo, f"LOOP {int(loop)}:RAMP?").strip().upper()
    except Exception:
        return None

# =============================
# API compatible functions (your app calls these)
# =============================

def get_current_temperature(cryo, loop: int = 1):
    """
    Backwards-friendly wrapper:
      - We interpret `loop` as "control loop number"
      - We read the controlling input channel temperature.
    """
    try:
        src = get_loop_source(cryo, loop=loop) or DEFAULT_INPUT_CHANNEL
        return get_input_temperature(cryo, channel=src)
    except Exception as e:
        log.error(f"Error reading current temperature for loop {loop}: {e}")
        return None

def set_temperature_setpoint(cryo, temp_k: float, loop: int = 1):
    """
    Set loop setpoint:
        LOOP <no>:SETPT <temp>
    """
    if temp_k < TEMP_MIN_K or temp_k > TEMP_MAX_K:
        log.warning(f"Temperature {temp_k} out of range ({TEMP_MIN_K}-{TEMP_MAX_K}); clipping")
        temp_k = float(np.clip(temp_k, TEMP_MIN_K, TEMP_MAX_K))

    write_safe(cryo, f"LOOP {int(loop)}:SETPT {temp_k:.6f}", delay_s=0.02)
    log.info(f"Loop {loop} setpoint -> {temp_k:.3f}")
    print(f"Temperature setpoint set to {temp_k:.3f} (Loop {loop})")
    current_parameters["temp_setpoint"] = temp_k

def set_ramp_rate(cryo, ramp_k_per_min: float, loop: int = 1):
    """
    Set ramp rate:
        LOOP <no>:RATE <rate>
    Manual: rate is Units per minute, range 0..100 (instrument side).
    """
    if ramp_k_per_min < RAMP_MIN_K_PER_MIN or ramp_k_per_min > RAMP_MAX_K_PER_MIN:
        log.warning(
            f"Ramp rate {ramp_k_per_min} out of range ({RAMP_MIN_K_PER_MIN}-{RAMP_MAX_K_PER_MIN}); clipping"
        )
        ramp_k_per_min = float(np.clip(ramp_k_per_min, RAMP_MIN_K_PER_MIN, RAMP_MAX_K_PER_MIN))

    write_safe(cryo, f"LOOP {int(loop)}:RATE {ramp_k_per_min:.6f}", delay_s=0.02)
    log.info(f"Loop {loop} ramp rate -> {ramp_k_per_min:.3f} units/min")
    print(f"Ramp rate set to {ramp_k_per_min:.3f} units/min (Loop {loop})")
    current_parameters["ramp_rate_k_per_min"] = ramp_k_per_min

def enable_ramp_mode(cryo, loop: int = 1):
    """
    Put the loop into ramp control type.
    Manual lists RampP as a loop TYPE value for ramp mode.
    """
    set_loop_type(cryo, "RAMPP", loop=loop)

def disable_ramp(cryo, loop: int = 1):
    """
    Disable ramping by returning to PID control (common choice).
    If you prefer OFF, change to set_loop_type(..., "OFF").
    """
    try:
        set_loop_type(cryo, "PID", loop=loop)
        log.info(f"Ramp disabled on loop {loop} (TYPE -> PID)")
    except Exception as e:
        log.error(f"Error disabling ramp: {e}")

def hold_temperature(cryo, loop: int = 1):
    """
    'Hold' = stop ramping by leaving setpoint where it is and returning to PID.
    """
    disable_ramp(cryo, loop=loop)
    log.info("Temperature hold active (ramp paused)")
    print("  Temperature hold active (ramp paused)")

def resume_ramp_to_next(cryo, target_temp_k: float, ramp_rate: float, loop: int = 1):
    """
    Resume a ramp:
      - set ramp rate
      - set setpoint
      - enable ramp mode
    """
    set_ramp_rate(cryo, ramp_rate, loop=loop)
    set_temperature_setpoint(cryo, target_temp_k, loop=loop)
    enable_ramp_mode(cryo, loop=loop)

def wait_for_temperature(
    cryo,
    target_temp_k,
    tolerance_k: float = TEMP_TOLERANCE_K,
    max_wait_s: float = 3600,
    stability_time_s: float = TEMP_STABILITY_TIME_S,
    loop: int = 1,
):
    """
    Wait for controlling-input temperature to be within tolerance for stability_time_s.
    """
    start_time = time.time()
    stability_start = None

    while True:
        if (time.time() - start_time) > max_wait_s:
            log.warning(f"Temperature control timeout (max wait {max_wait_s}s)")
            print(f"Temperature control timeout. Target: {target_temp_k:.2f}, Max wait: {max_wait_s}s")
            return False

        current_temp = get_current_temperature(cryo, loop=loop)
        if current_temp is None:
            time.sleep(1)
            continue

        if abs(current_temp - target_temp_k) <= tolerance_k:
            if stability_start is None:
                stability_start = time.time()
                log.info(f"Within tolerance ({tolerance_k}) at {current_temp:.3f}")
                print(f"Temperature within tolerance at {current_temp:.3f}")

            if (time.time() - stability_start) >= stability_time_s:
                log.info(f"Stable at {current_temp:.3f} for {stability_time_s}s")
                print(f"Temperature STABLE at {current_temp:.3f}")
                return True
        else:
            stability_start = None
            print(f"Temperature: {current_temp:.3f} / Target: {target_temp_k:.3f}")

        time.sleep(2)

def generate_temperature_sweep_points(start_k, end_k, interval_k):
    if interval_k <= 0:
        raise ValueError("Measurement interval must be > 0 K")

    if end_k >= start_k:
        temps = list(np.arange(start_k, end_k + interval_k / 2, interval_k))
    else:
        temps = list(np.arange(start_k, end_k - interval_k / 2, -interval_k))

    if not temps or abs(temps[-1] - end_k) > 1e-6:
        temps.append(end_k)

    return np.array(temps, dtype=float)

# =============================
# Parameter Management
# =============================

def get_cryogenic_parameters():
    print("\n" + "=" * 60)
    print("CRYOGENIC MEASUREMENT PARAMETERS")
    print("=" * 60)

    end_temp = safe_float_input(
        f"Enter end temperature (K) [default {current_parameters['temp_end_k']}]: ",
        current_parameters["temp_end_k"],
    )
    ramp_rate = safe_float_input(
        f"Enter ramp rate (K/min) [default {current_parameters['ramp_rate_k_per_min']}]: ",
        current_parameters["ramp_rate_k_per_min"],
    )
    meas_interval = safe_float_input(
        f"Enter measurement interval (K) [default {current_parameters['meas_interval_k']}]: ",
        current_parameters["meas_interval_k"],
    )

    if end_temp < TEMP_MIN_K or end_temp > TEMP_MAX_K:
        log.warning(f"End temperature {end_temp}K out of range; clipping")
        end_temp = float(np.clip(end_temp, TEMP_MIN_K, TEMP_MAX_K))
    if ramp_rate < RAMP_MIN_K_PER_MIN or ramp_rate > RAMP_MAX_K_PER_MIN:
        log.warning(f"Ramp rate {ramp_rate}K/min out of range; clipping")
        ramp_rate = float(np.clip(ramp_rate, RAMP_MIN_K_PER_MIN, RAMP_MAX_K_PER_MIN))
    if meas_interval < MEAS_INTERVAL_MIN_K:
        log.warning(f"Measurement interval {meas_interval}K too small; using {MEAS_INTERVAL_MIN_K}K")
        meas_interval = MEAS_INTERVAL_MIN_K

    current_parameters["temp_end_k"] = end_temp
    current_parameters["ramp_rate_k_per_min"] = ramp_rate
    current_parameters["meas_interval_k"] = meas_interval
    current_parameters["enable_cryo"] = True

    print(f"\nTemperature sweep: current temperature → {end_temp:.1f}K at {ramp_rate:.2f}K/min")
    print(f"Measurement interval: {meas_interval:.1f}K")
    print("Start temperature will be read from the controller at sweep start.\n")

    return {
        "end_k": end_temp,
        "ramp_rate_k_per_min": ramp_rate,
        "meas_interval_k": meas_interval,
    }

# =============================
# Utility
# =============================

def setup(resource_name=None):
    """
    Connect to and initialize Cryo-Con.
    """
    cryo = connect_cryocon(resource_name)
    cryo = initialize_cryocon(cryo)

    set_loop_source(cryo, "B", loop=2)  # or "B" if loop 2 uses input B
    set_loop_type(cryo, "PID", loop=2)  # ensure loop 2 is active
    return cryo
