"""
Filename: parameter_validators.py
Author: Ethan Ruddell
Date: 2026-02-13
Description: Centralized parameter validation for all measurement types.
             Consolidates validation logic scattered across instrument modules.
"""

import numpy as np
import logging

log = logging.getLogger(__name__)


class ValidationError(ValueError):
    """Raised when parameter validation fails."""
    pass


# =============================
# Frequency Validation
# =============================

def validate_frequency(freq_hz: float, min_hz: float = 20, max_hz: float = 2e6) -> float:
    """Validate and clip frequency to acceptable range."""
    if freq_hz < min_hz or freq_hz > max_hz:
        log.warning(
            f"Frequency {freq_hz} Hz out of range [{min_hz}, {max_hz}], clipping"
        )
        return float(np.clip(freq_hz, min_hz, max_hz))
    return float(freq_hz)


def validate_frequency_range(
    freq_start: float,
    freq_stop: float,
    min_hz: float = 20,
    max_hz: float = 2e6
) -> tuple[float, float]:
    """Validate frequency sweep range."""
    freq_start = validate_frequency(freq_start, min_hz, max_hz)
    freq_stop = validate_frequency(freq_stop, min_hz, max_hz)
    
    if freq_start >= freq_stop:
        raise ValidationError(
            f"Start frequency ({freq_start}) must be less than stop ({freq_stop})"
        )
    return freq_start, freq_stop


# =============================
# Temperature Validation
# =============================

def validate_temperature(
    temp_k: float,
    min_k: float = 1.4,
    max_k: float = 325.0
) -> float:
    """Validate and clip temperature to acceptable range."""
    if temp_k < min_k or temp_k > max_k:
        log.warning(
            f"Temperature {temp_k} K out of range [{min_k}, {max_k}], clipping"
        )
        return float(np.clip(temp_k, min_k, max_k))
    return float(temp_k)


def validate_temperature_range(
    temp_start: float,
    temp_end: float,
    min_k: float = 1.4,
    max_k: float = 325.0
) -> tuple[float, float]:
    """Validate temperature sweep range."""
    temp_start = validate_temperature(temp_start, min_k, max_k)
    temp_end = validate_temperature(temp_end, min_k, max_k)
    
    if temp_start == temp_end:
        raise ValidationError(
            f"Start temperature ({temp_start}) must differ from end ({temp_end})"
        )
    return temp_start, temp_end


# =============================
# Ramp Rate Validation
# =============================

def validate_ramp_rate(
    ramp_k_per_min: float,
    min_rate: float = 0.1,
    max_rate: float = 50.0
) -> float:
    """Validate and clip ramp rate to acceptable range."""
    if ramp_k_per_min < min_rate or ramp_k_per_min > max_rate:
        log.warning(
            f"Ramp rate {ramp_k_per_min} K/min out of range [{min_rate}, {max_rate}], clipping"
        )
        return float(np.clip(ramp_k_per_min, min_rate, max_rate))
    return float(ramp_k_per_min)


# =============================
# Voltage Validation
# =============================

def validate_voltage(
    voltage_v: float,
    min_v: float = -100,
    max_v: float = 100
) -> float:
    """Validate and clip voltage to acceptable range."""
    if voltage_v < min_v or voltage_v > max_v:
        log.warning(
            f"Voltage {voltage_v} V out of range [{min_v}, {max_v}], clipping"
        )
        return float(np.clip(voltage_v, min_v, max_v))
    return float(voltage_v)


def validate_voltage_range(
    v_start: float,
    v_stop: float,
    min_v: float = -100,
    max_v: float = 100
) -> tuple[float, float]:
    """Validate voltage sweep range."""
    v_start = validate_voltage(v_start, min_v, max_v)
    v_stop = validate_voltage(v_stop, min_v, max_v)
    
    if v_start == v_stop:
        raise ValidationError(
            f"Start voltage ({v_start}) must differ from stop ({v_stop})"
        )
    return v_start, v_stop


# =============================
# Current Validation
# =============================

def validate_compliance(
    compliance_a: float,
    min_a: float = 1e-6,
    max_a: float = 10.0
) -> float:
    """Validate current compliance limit."""
    if compliance_a < min_a or compliance_a > max_a:
        log.warning(
            f"Compliance {compliance_a} A out of range [{min_a}, {max_a}], clipping"
        )
        return float(np.clip(compliance_a, min_a, max_a))
    return float(compliance_a)


# =============================
# Point Count Validation
# =============================

def validate_point_count(
    num_points: int,
    min_points: int = 2,
    max_points: int = 10000
) -> int:
    """Validate number of measurement points."""
    if num_points < min_points or num_points > max_points:
        log.warning(
            f"Point count {num_points} out of range [{min_points}, {max_points}], clipping"
        )
        return int(np.clip(num_points, min_points, max_points))
    return int(num_points)


def validate_step_size(
    step: float,
    min_step: float = 1e-6,
    max_step: float = 1e6
) -> float:
    """Validate measurement step/interval."""
    if step < min_step or step > max_step:
        log.warning(
            f"Step size {step} out of range [{min_step}, {max_step}], clipping"
        )
        return float(np.clip(step, min_step, max_step))
    return float(step)


# =============================
# AC Level/Amplitude Validation
# =============================

def validate_ac_level(
    ac_level_v: float,
    min_v: float = 0.001,
    max_v: float = 10.0
) -> float:
    """Validate AC oscillator level (Vrms)."""
    if ac_level_v < min_v or ac_level_v > max_v:
        log.warning(
            f"AC level {ac_level_v} V out of range [{min_v}, {max_v}], clipping"
        )
        return float(np.clip(ac_level_v, min_v, max_v))
    return float(ac_level_v)


# =============================
# Choice Validation
# =============================

def validate_choice(value: str, valid_options: list[str], default: str = None) -> str:
    """Validate that choice is in list of valid options."""
    value_upper = value.upper()
    valid_upper = [opt.upper() for opt in valid_options]
    
    if value_upper not in valid_upper:
        if default:
            log.warning(
                f"Invalid choice '{value}', using default '{default}'"
            )
            return default.upper()
        else:
            raise ValidationError(
                f"'{value}' not in valid options: {valid_options}"
            )
    return value_upper


# =============================
# Multiple Parameter Validation
# =============================

class ParameterSet:
    """Helper class for validating related parameters together."""
    
    def __init__(self):
        self.errors = []
        self.warnings = []
        self.values = {}
    
    def add_frequency_sweep(
        self,
        freq_start: float,
        freq_stop: float,
        num_points: int
    ):
        """Validate frequency sweep parameters."""
        try:
            self.values['freq_start'] = validate_frequency(freq_start)
            self.values['freq_stop'] = validate_frequency(freq_stop)
            self.values['num_points'] = validate_point_count(num_points)
            
            if self.values['freq_start'] >= self.values['freq_stop']:
                self.errors.append("Start frequency must be less than stop frequency")
                
        except ValidationError as e:
            self.errors.append(str(e))
    
    def add_voltage_sweep(
        self,
        v_start: float,
        v_stop: float,
        num_points: int
    ):
        """Validate voltage sweep parameters."""
        try:
            self.values['v_start'] = validate_voltage(v_start)
            self.values['v_stop'] = validate_voltage(v_stop)
            self.values['num_points'] = validate_point_count(num_points)
            
            if self.values['v_start'] >= self.values['v_stop']:
                self.errors.append("Start voltage must be less than stop voltage")
                
        except ValidationError as e:
            self.errors.append(str(e))
    
    def is_valid(self) -> bool:
        """Check if all validations passed."""
        return len(self.errors) == 0
    
    def get_validated_values(self) -> dict:
        """Return validated parameter dictionary."""
        if not self.is_valid():
            raise ValidationError(f"Validation errors: {'; '.join(self.errors)}")
        return self.values


# =============================
# Bulk Validation Helper
# =============================

def validate_measurement_params(
    measurement_type: str,
    params: dict
) -> dict:
    """
    One-stop validation for common measurement parameter sets.
    
    Args:
        measurement_type: "frequency_sweep", "voltage_sweep", "cryo_sweep", etc.
        params: Dictionary of parameter values
    
    Returns:
        Dictionary of validated and clipped parameters
    
    Raises:
        ValidationError: If critical validation fails
    """
    validated = dict(params)  # Copy input
    
    if measurement_type == "frequency_sweep":
        try:
            validated['freq_start'] = validate_frequency(params.get('freq_start', 1000))
            validated['freq_stop'] = validate_frequency(params.get('freq_stop', 1e6))
            validated['num_points'] = validate_point_count(params.get('num_points', 201))
            validated['ac_level'] = validate_ac_level(params.get('ac_level', 1.0))
            if 'dc_bias' in params:
                validated['dc_bias'] = validate_voltage(params['dc_bias'])
        except Exception as e:
            raise ValidationError(f"Frequency sweep validation failed: {e}")
    
    elif measurement_type == "voltage_sweep":
        try:
            validated['v_start'] = validate_voltage(params.get('v_start', 0))
            validated['v_stop'] = validate_voltage(params.get('v_stop', 5))
            validated['num_points'] = validate_point_count(params.get('num_points', 201))
            if 'compliance' in params:
                validated['compliance'] = validate_compliance(params['compliance'])
        except Exception as e:
            raise ValidationError(f"Voltage sweep validation failed: {e}")
    
    elif measurement_type == "cryo_sweep":
        try:
            validated['temp_start'] = validate_temperature(params.get('temp_start', 10))
            validated['temp_end'] = validate_temperature(params.get('temp_end', 300))
            validated['ramp_rate'] = validate_ramp_rate(params.get('ramp_rate', 5))
            validated['num_points'] = validate_point_count(params.get('num_points', 10))
        except Exception as e:
            raise ValidationError(f"Cryo sweep validation failed: {e}")
    
    return validated
