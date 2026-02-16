"""
Filename: cryo.py
Author: Ethan Ruddell
Date: 2026-02-13
Description: Contains all constants and functions for Cryocon 32B temperature controller.
             Integration module for temperature-dependent measurements.
"""

import pyvisa
import numpy as np
import time
import logging
from gpib_utils import safe_float_input

log = logging.getLogger(__name__)

# =============================
# User settings and constants
# =============================
GPIB_ADDRESS = "GPIB0::___::INSTR"   # GPIB address for Cryocon 32B (to be configured)

# Temperature limits and defaults
TEMP_MIN_K = 1.4                    # Minimum temperature (K)
TEMP_MAX_K = 325.0                  # Maximum temperature (K)
TEMP_DEFAULT_START_K = 100.0        # Default start temperature (K)
TEMP_DEFAULT_END_K = 200.0          # Default end temperature (K)

# Ramp rate limits and defaults
RAMP_MIN_K_PER_MIN = 0.1            # Minimum ramp rate (K/min)
RAMP_MAX_K_PER_MIN = 50.0           # Maximum ramp rate (K/min, depends on cooling system)
RAMP_DEFAULT_K_PER_MIN = 5.0        # Default ramp rate (K/min)

# Temperature measurement interval defaults
MEAS_INTERVAL_MIN_K = 1.0           # Minimum measurement interval (K)
MEAS_INTERVAL_DEFAULT_K = 10.0      # Default measurement interval (K)

# Temperature stability thresholds
TEMP_TOLERANCE_K = 0.5              # Temperature tolerance for "reached" state (K)
TEMP_STABILITY_TIME_S = 30.0        # Stability wait time (seconds)

# =============================
# Measurement Parameters Storage
# =============================
current_parameters = {
    "temp_start_k": TEMP_DEFAULT_START_K,
    "temp_end_k": TEMP_DEFAULT_END_K,
    "ramp_rate_k_per_min": RAMP_DEFAULT_K_PER_MIN,
    "meas_interval_k": MEAS_INTERVAL_DEFAULT_K,
    "enable_cryo": False,
}

# =============================
# Connection and Initialization
# =============================


def connect_cryocon(resource_name=None):
    """
    Open VISA session to Cryocon 32B temperature controller.
    
    Args:
        resource_name: GPIB address string (uses GPIB_ADDRESS if None)
    
    Returns:
        VISA instrument instance
    """
    if resource_name is None:
        resource_name = GPIB_ADDRESS
    
    try:
        rm = pyvisa.ResourceManager()
        log.info(f"Attempting to connect to Cryocon 32B at: {resource_name}")
        print(f"Attempting to connect to Cryocon 32B at: {resource_name}")
        
        cryo = rm.open_resource(resource_name)
        cryo.timeout = 30000  # 30 second timeout
        cryo.read_termination = '\n'
        cryo.write_termination = '\n'
        
        # Verify connection
        idn = cryo.query("*IDN?")
        log.info(f"Connected to: {idn.strip()}")
        print(f"Connected to: {idn.strip()}")
        
        return cryo
    except Exception as e:
        log.error(f"Error connecting to Cryocon 32B: {e}")
        print(f"Error connecting to Cryocon 32B: {e}")
        raise


def disconnect_cryocon(cryo):
    """Safely disconnect from Cryocon 32B."""
    try:
        # Set ramp to safe state (slow)
        cryo.write("RAMP 1,0")  # Disable ramp on loop 1
        time.sleep(0.5)
        cryo.close()
        log.info("Cryocon 32B disconnected")
        print("Cryocon 32B disconnected")
    except Exception as e:
        log.warning(f"Error during disconnect: {e}")
        pass


def initialize_cryocon(cryo):
    """
    Basic initialization for Cryocon 32B.
    """
    try:
        # Query control loop configuration
        cryo.write("*RST")
        time.sleep(1)
        
        # Get authorized state
        cryo.write("*CLS")
        
        log.info("Cryocon 32B initialized")
        print("Cryocon 32B initialized")
    except Exception as e:
        log.error(f"Initialization error: {e}")
        raise


def check_errors(cryo):
    """Check for instrument errors."""
    errors = []
    try:
        error = cryo.query(":SYST:ERR?")
        if not error.startswith('+0') and "No error" not in error:
            log.warning(f"Cryocon error: {error.strip()}")
            errors.append(error.strip())
    except Exception as e:
        log.debug(f"Error checking Cryocon errors: {e}")
    return errors


# =============================
# Temperature Control Functions
# =============================

def get_current_temperature(cryo, loop=1):
    """
    Get current temperature reading.
    
    Args:
        cryo: VISA instrument instance
        loop: Control loop number (1 or 2 for Cryocon 32B)
    
    Returns:
        float: Current temperature in Kelvin
    """
    try:
        # Query sensor temperature: KRDG? [loop]
        # Returns temperature in Kelvin
        response = cryo.query(f"KRDG? {loop}")
        temp_k = float(response.strip())
        return temp_k
    except Exception as e:
        log.error(f"Error reading temperature on loop {loop}: {e}")
        return None


def set_temperature_setpoint(cryo, temp_k, loop=1):
    """
    Set temperature setpoint.
    
    Args:
        cryo: VISA instrument instance
        temp_k: Target temperature in Kelvin
        loop: Control loop number (1 or 2)
    """
    # Validate temperature range
    if temp_k < TEMP_MIN_K or temp_k > TEMP_MAX_K:
        log.warning(f"Warning: Temperature {temp_k} K out of range ({TEMP_MIN_K}-{TEMP_MAX_K} K)")
        temp_k = np.clip(temp_k, TEMP_MIN_K, TEMP_MAX_K)
    
    try:
        # SETP loop,setpoint
        cryo.write(f"SETP {loop},{temp_k:.3f}")
        log.info(f"Temperature setpoint set to {temp_k:.3f} K on loop {loop}")
        print(f"Temperature setpoint set to {temp_k:.3f} K")
        current_parameters["temp_setpoint"] = temp_k
    except Exception as e:
        log.error(f"Error setting temperature: {e}")
        raise


def set_ramp_rate(cryo, ramp_k_per_min, loop=1):
    """
    Set temperature ramp rate.
    
    Args:
        cryo: VISA instrument instance
        ramp_k_per_min: Ramp rate in K/min
        loop: Control loop number
    """
    # Validate ramp rate
    if ramp_k_per_min < RAMP_MIN_K_PER_MIN or ramp_k_per_min > RAMP_MAX_K_PER_MIN:
        log.warning(f"Warning: Ramp rate {ramp_k_per_min} K/min out of range ({RAMP_MIN_K_PER_MIN}-{RAMP_MAX_K_PER_MIN} K/min)")
        ramp_k_per_min = np.clip(ramp_k_per_min, RAMP_MIN_K_PER_MIN, RAMP_MAX_K_PER_MIN)
    
    try:
        # RAMP loop,on/off,rate
        # on/off: 0 = OFF, 1 = ON
        cryo.write(f"RAMP {loop},1,{ramp_k_per_min:.3f}")
        log.info(f"Ramp rate set to {ramp_k_per_min:.3f} K/min on loop {loop}")
        print(f"Ramp rate set to {ramp_k_per_min:.3f} K/min")
        current_parameters["ramp_rate_k_per_min"] = ramp_k_per_min
    except Exception as e:
        log.error(f"Error setting ramp rate: {e}")
        raise


def disable_ramp(cryo, loop=1):
    """
    Disable temperature ramp (go to manual control).
    
    Args:
        cryo: VISA instrument instance
        loop: Control loop number
    """
    try:
        cryo.write(f"RAMP {loop},0")
        log.info(f"Ramp disabled on loop {loop}")
    except Exception as e:
        log.error(f"Error disabling ramp: {e}")


def hold_temperature(cryo, loop=1):
    """
    Hold temperature constant by pausing the ramp.

    The controller maintains the current setpoint but stops ramping
    toward any new setpoint.  Call resume_ramp_to_next() to continue.

    Args:
        cryo: VISA instrument instance
        loop: Control loop number
    """
    disable_ramp(cryo, loop)
    log.info("Temperature hold: ramp paused, maintaining current setpoint")
    print("  Temperature hold active (ramp paused)")


def resume_ramp_to_next(cryo, target_temp_k, ramp_rate, loop=1):
    """
    Resume ramping to a new temperature setpoint.

    Args:
        cryo: VISA instrument instance
        target_temp_k: Next target temperature (K)
        ramp_rate: Ramp rate in K/min
        loop: Control loop number
    """
    set_temperature_setpoint(cryo, target_temp_k, loop)
    set_ramp_rate(cryo, ramp_rate, loop)


def wait_for_temperature(cryo, target_temp_k, tolerance_k=TEMP_TOLERANCE_K, 
                         max_wait_s=3600, stability_time_s=TEMP_STABILITY_TIME_S, loop=1):
    """
    Wait for temperature to reach target within tolerance.
    
    Args:
        cryo: VISA instrument instance
        target_temp_k: Target temperature (K)
        tolerance_k: Temperature tolerance (K)
        max_wait_s: Maximum wait time (seconds)
        stability_time_s: Time to maintain stable temperature before returning
        loop: Control loop number
    
    Returns:
        bool: True if temperature reached, False if timeout
    """
    start_time = time.time()
    stability_start = None
    
    while True:
        elapsed_time = time.time() - start_time
        
        if elapsed_time > max_wait_s:
            log.warning(f"Temperature control timeout (max wait {max_wait_s}s)")
            print(f"Temperature control timeout. Target: {target_temp_k:.2f}K, Max wait: {max_wait_s}s")
            return False
        
        current_temp = get_current_temperature(cryo, loop)
        if current_temp is None:
            time.sleep(1)
            continue
        
        # Check if temperature is within tolerance
        if abs(current_temp - target_temp_k) <= tolerance_k:
            if stability_start is None:
                stability_start = time.time()
                log.info(f"Temperature within tolerance ({tolerance_k}K) at {current_temp:.2f}K")
                print(f"Temperature within tolerance at {current_temp:.2f}K")
            
            # Check if stable for required time
            stability_elapsed = time.time() - stability_start
            if stability_elapsed >= stability_time_s:
                log.info(f"Temperature stable at {current_temp:.2f}K for {stability_time_s}s")
                print(f"Temperature STABLE at {current_temp:.2f}K")
                return True
        else:
            # Temperature went out of tolerance
            stability_start = None
            print(f"Temperature: {current_temp:.2f}K / Target: {target_temp_k:.2f}K")
        
        time.sleep(2)  # Check every 2 seconds


def generate_temperature_sweep_points(start_k, end_k, interval_k):
    """
    Generate array of temperature measurement points.
    
    Args:
        start_k: Start temperature (K)
        end_k: End temperature (K)
        interval_k: Interval between measurements (K)
    
    Returns:
        numpy array of temperature points
    """
    # Generate points
    if end_k > start_k:
        temps = np.arange(start_k, end_k + interval_k/2, interval_k)
        if temps[-1] > end_k:
            temps = temps[:-1]
        temps = np.append(temps, end_k)  # Ensure end point is included
    else:
        temps = np.arange(start_k, end_k - interval_k/2, -interval_k)
        if temps[-1] < end_k:
            temps = temps[:-1]
        temps = np.append(temps, end_k)
    
    # Remove duplicates and sort
    temps = np.unique(np.sort(temps))
    return temps


# =============================
# Parameter Management
# =============================

def get_cryogenic_parameters():
    """
    Get cryogenic measurement parameters from user input.
    
    Returns:
        dict: Temperature parameters
    """
    print("\n" + "="*60)
    print("CRYOGENIC MEASUREMENT PARAMETERS")
    print("="*60)
    
    start_temp = safe_float_input(f"Enter start temperature (K) [default {current_parameters['temp_start_k']}]: ", current_parameters['temp_start_k'])
    end_temp = safe_float_input(f"Enter end temperature (K) [default {current_parameters['temp_end_k']}]: ", current_parameters['temp_end_k'])
    ramp_rate = safe_float_input(f"Enter ramp rate (K/min) [default {current_parameters['ramp_rate_k_per_min']}]: ", current_parameters['ramp_rate_k_per_min'])
    meas_interval = safe_float_input(f"Enter measurement interval (K) [default {current_parameters['meas_interval_k']}]: ", current_parameters['meas_interval_k'])
    
    # Validate inputs
    if start_temp < TEMP_MIN_K or start_temp > TEMP_MAX_K:
        log.warning(f"Start temperature {start_temp}K out of range, clipping to {TEMP_MIN_K}-{TEMP_MAX_K}K")
        start_temp = np.clip(start_temp, TEMP_MIN_K, TEMP_MAX_K)
    
    if end_temp < TEMP_MIN_K or end_temp > TEMP_MAX_K:
        log.warning(f"End temperature {end_temp}K out of range, clipping to {TEMP_MIN_K}-{TEMP_MAX_K}K")
        end_temp = np.clip(end_temp, TEMP_MIN_K, TEMP_MAX_K)
    
    if ramp_rate < RAMP_MIN_K_PER_MIN or ramp_rate > RAMP_MAX_K_PER_MIN:
        log.warning(f"Ramp rate {ramp_rate}K/min out of range, clipping to {RAMP_MIN_K_PER_MIN}-{RAMP_MAX_K_PER_MIN}K/min")
        ramp_rate = np.clip(ramp_rate, RAMP_MIN_K_PER_MIN, RAMP_MAX_K_PER_MIN)
    
    if meas_interval < MEAS_INTERVAL_MIN_K:
        log.warning(f"Measurement interval {meas_interval}K too small, using {MEAS_INTERVAL_MIN_K}K")
        meas_interval = MEAS_INTERVAL_MIN_K
    
    # Update stored parameters
    current_parameters["temp_start_k"] = start_temp
    current_parameters["temp_end_k"] = end_temp
    current_parameters["ramp_rate_k_per_min"] = ramp_rate
    current_parameters["meas_interval_k"] = meas_interval
    current_parameters["enable_cryo"] = True
    
    # Generate temperature points
    temp_points = generate_temperature_sweep_points(start_temp, end_temp, meas_interval)
    
    print(f"\nTemperature sweep: {start_temp:.1f}K → {end_temp:.1f}K at {ramp_rate:.2f}K/min")
    print(f"Measurement interval: {meas_interval:.1f}K")
    print(f"Total measurement points: {len(temp_points)}")
    print(f"Estimated time: {abs(end_temp - start_temp) / ramp_rate:.1f} minutes\n")
    
    return {
        "start_k": start_temp,
        "end_k": end_temp,
        "ramp_rate_k_per_min": ramp_rate,
        "meas_interval_k": meas_interval,
        "temp_points": temp_points,
    }


# =============================
# Utility Functions
# =============================

def setup(resource_name=None):
    """
    Connect to and initialize Cryocon 32B.
    
    Args:
        resource_name: GPIB address (uses default if None)
    
    Returns:
        VISA instrument instance
    """
    cryo = connect_cryocon(resource_name)
    initialize_cryocon(cryo)
    return cryo
