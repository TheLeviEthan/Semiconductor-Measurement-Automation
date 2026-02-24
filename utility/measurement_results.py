"""
Filename: measurement_results.py
Author: Ethan Ruddell
Date: 2026-02-13
Description: Standardized result classes for different measurement types.
             Provides consistent interface for all measurements across instruments.
"""

from dataclasses import dataclass, field
from typing import Optional, Any
import numpy as np
from datetime import datetime


@dataclass
class MeasurementMetadata:
    """Standard metadata for any measurement."""
    timestamp: datetime = field(default_factory=datetime.now)
    instrument: str = ""
    measurement_type: str = ""
    temperature_k: Optional[float] = None
    dc_bias_v: Optional[float] = None
    ac_level_v: Optional[float] = None
    custom_params: dict = field(default_factory=dict)
    
    def to_dict(self) -> dict:
        """Convert to dictionary for CSV/JSON serialization."""
        return {
            'timestamp': self.timestamp.isoformat(),
            'instrument': self.instrument,
            'measurement_type': self.measurement_type,
            'temperature_k': self.temperature_k,
            'dc_bias_v': self.dc_bias_v,
            'ac_level_v': self.ac_level_v,
            **self.custom_params
        }


@dataclass
class FrequencySweepResult:
    """Result from a frequency sweep measurement."""
    frequency: np.ndarray
    primary: np.ndarray              # Main measurement (Z, C, L, etc.)
    secondary: Optional[np.ndarray] = None  # Phase, D, Q, etc.
    metadata: MeasurementMetadata = field(default_factory=MeasurementMetadata)
    
    def __post_init__(self):
        """Validate dimensions match."""
        if len(self.frequency) != len(self.primary):
            raise ValueError("Frequency and primary data arrays must have same length")
        if self.secondary is not None and len(self.frequency) != len(self.secondary):
            raise ValueError("Secondary data must match frequency array length")
    
    def to_csv_array(self) -> np.ndarray:
        """Convert to CSV-friendly matrix format."""
        columns = [self.frequency, self.primary]
        if self.secondary is not None:
            columns.append(self.secondary)
        return np.column_stack(columns)
    
    def get_csv_header(self, freq_label="frequency_Hz", prim_label="primary", 
                       sec_label="secondary") -> str:
        """Generate CSV header for this result."""
        header = f"{freq_label}, {prim_label}"
        if self.secondary is not None:
            header += f", {sec_label}"
        return header


@dataclass
class VoltageSweepResult:
    """Result from voltage-dependent measurement."""
    voltage: np.ndarray
    primary: np.ndarray
    secondary: Optional[np.ndarray] = None
    metadata: MeasurementMetadata = field(default_factory=MeasurementMetadata)
    
    def __post_init__(self):
        if len(self.voltage) != len(self.primary):
            raise ValueError("Voltage and primary data arrays must have same length")
        if self.secondary is not None and len(self.voltage) != len(self.secondary):
            raise ValueError("Secondary data must match voltage array length")
    
    def to_csv_array(self) -> np.ndarray:
        columns = [self.voltage, self.primary]
        if self.secondary is not None:
            columns.append(self.secondary)
        return np.column_stack(columns)


@dataclass
class IVCurveResult:
    """Result from I-V curve measurement."""
    voltage: np.ndarray
    current: np.ndarray
    resistance_ohm: Optional[float] = None
    metadata: MeasurementMetadata = field(default_factory=MeasurementMetadata)
    
    def __post_init__(self):
        if len(self.voltage) != len(self.current):
            raise ValueError("Voltage and current arrays must have same length")
    
    def calculate_resistance(self) -> float:
        """Calculate average slope resistance from I-V curve."""
        if len(self.voltage) < 2:
            return np.nan
        # Simple linear fit
        z = np.polyfit(self.current, self.voltage, 1)
        self.resistance_ohm = z[0]
        return self.resistance_ohm
    
    def to_csv_array(self) -> np.ndarray:
        return np.column_stack([self.voltage, self.current])


@dataclass
class CycleResult:
    """Result from repeating cycle measurement (e.g., C-V butterflies)."""
    cycles: list[np.ndarray]  # List of arrays, one per cycle
    primary_cycles: list[np.ndarray]
    secondary_cycles: Optional[list[np.ndarray]] = None
    metadata: MeasurementMetadata = field(default_factory=MeasurementMetadata)
    
    @property
    def num_cycles(self) -> int:
        return len(self.cycles)
    
    @property
    def total_points(self) -> int:
        return sum(len(c) for c in self.cycles)
    
    def flatten_for_csv(self) -> tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
        """Flatten cycles into single arrays for CSV, with cycle indices."""
        cycle_indices = np.concatenate([
            np.full(len(c), i + 1, dtype=int)
            for i, c in enumerate(self.cycles)
        ])
        coord_flat = np.concatenate(self.cycles)
        primary_flat = np.concatenate(self.primary_cycles)
        
        if self.secondary_cycles:
            secondary_flat = np.concatenate(self.secondary_cycles)
        else:
            secondary_flat = None
        
        return cycle_indices, coord_flat, primary_flat, secondary_flat
    
    def to_csv_array(self) -> np.ndarray:
        """Convert cycles to CSV matrix with cycle index column."""
        cycle_idx, coord, primary, secondary = self.flatten_for_csv()
        
        if secondary is not None:
            return np.column_stack([cycle_idx, coord, primary, secondary])
        else:
            return np.column_stack([cycle_idx, coord, primary])


@dataclass
class TransistorCharacteristicsResult:
    """Result from transistor measurement (Ids vs Vds/Vgs)."""
    vgs: np.ndarray
    vds: np.ndarray
    ids: np.ndarray
    igs: Optional[np.ndarray] = None
    metadata: MeasurementMetadata = field(default_factory=MeasurementMetadata)
    
    @property
    def unique_vgs_steps(self) -> np.ndarray:
        """Get unique gate voltages applied."""
        return np.unique(self.vgs)
    
    @property  
    def unique_vds_steps(self) -> np.ndarray:
        """Get unique drain voltages applied."""
        return np.unique(self.vds)
    
    def get_output_curve_at_vgs(self, vgs_target: float, tolerance: float = 0.01) -> np.ndarray:
        """Extract output curve (Ids vs Vds) at a target gate voltage."""
        mask = np.abs(self.vgs - vgs_target) < tolerance
        return np.column_stack([self.vds[mask], self.ids[mask]])
    
    def to_csv_array(self) -> np.ndarray:
        columns = [self.vgs, self.vds, self.ids]
        if self.igs is not None:
            columns.append(self.igs)
        return np.column_stack(columns)


@dataclass  
class SinglePointResult:
    """Result from single-point measurement (no sweep)."""
    value: float
    uncertainty: Optional[float] = None
    unit: str = ""
    metadata: MeasurementMetadata = field(default_factory=MeasurementMetadata)
    
    def __str__(self) -> str:
        if self.uncertainty:
            return f"{self.value:.6e} ± {self.uncertainty:.6e} {self.unit}"
        else:
            return f"{self.value:.6e} {self.unit}"


@dataclass
class TemperatureSweepResult:
    """Result from cryogenic temperature sweep."""
    temperatures: np.ndarray
    measurement_results: list[Any]  # List of any measurement results
    metadata: MeasurementMetadata = field(default_factory=MeasurementMetadata)
    
    @property
    def num_temperatures(self) -> int:
        return len(self.temperatures)
    
    def get_results_at_temperature(self, temp_k: float, tolerance: float = 0.1) -> Any:
        """Get measurement results for a specific temperature."""
        idx = np.argmin(np.abs(self.temperatures - temp_k))
        if abs(self.temperatures[idx] - temp_k) > tolerance:
            return None  # Temperature not in sweep
        return self.measurement_results[idx]


# =============================
# Factory Functions
# =============================

def create_frequency_sweep_result(
    frequency: np.ndarray,
    primary: np.ndarray,
    secondary: Optional[np.ndarray] = None,
    **metadata_kwargs
) -> FrequencySweepResult:
    """Convenience factory for frequency sweep results."""
    metadata = MeasurementMetadata(**metadata_kwargs)
    return FrequencySweepResult(frequency, primary, secondary, metadata)


def create_voltage_sweep_result(
    voltage: np.ndarray,
    primary: np.ndarray,
    secondary: Optional[np.ndarray] = None,
    **metadata_kwargs
) -> VoltageSweepResult:
    """Convenience factory for voltage sweep results."""
    metadata = MeasurementMetadata(**metadata_kwargs)
    return VoltageSweepResult(voltage, primary, secondary, metadata)


def create_iv_curve_result(
    voltage: np.ndarray,
    current: np.ndarray,
    **metadata_kwargs
) -> IVCurveResult:
    """Convenience factory for I-V curve results."""
    metadata = MeasurementMetadata(**metadata_kwargs)
    return IVCurveResult(voltage, current, metadata=metadata)
