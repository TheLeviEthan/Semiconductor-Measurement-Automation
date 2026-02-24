"""
Filename: measurements_config.py
Author: Ethan Ruddell
Date: 2026-02-13
Description: Centralized configuration for all measurement types and parameters.
             Used by both CLI and GUI to ensure consistency.
"""

# =============================
# Measurement Lists
# =============================

PIA_MEASUREMENTS = [
    "Impedance Magnitude vs Frequency",
    "Impedance Phase vs Frequency",
    "Capacitance vs Frequency",
    "Tan Loss vs Frequency",
    "C-V Butterfly Cycles: Cp vs Voltage",
    "C-V Butterfly Cycles: Permittivity vs Voltage",
    "Permittivity vs Frequency",
    "R-X (Resistance-Reactance) vs Frequency",
    "G-B (Conductance-Susceptance) vs Frequency",
    "Y-θ (Admittance) vs Frequency",
]

PSPA_MEASUREMENTS = [
    "Transistor Output Characteristics",
    "Transistor Transfer Characteristics",
    "I-V Curve (Unidirectional)",
    "I-V Curve (Bidirectional)",
    "Pulsed I-V (Single Device)",
    "Pulsed Transistor",
    "Diode I-V Characteristics",
    "Gate Leakage Current",
    "Breakdown Voltage",
    "Resistance Measurement",
]

LCR_MEASUREMENTS = [
    "Impedance vs Frequency",
    "Capacitance vs Frequency",
    "Inductance vs Frequency",
    "C-V Sweep (Single)",
    "C-V Butterfly Cycles",
    "Single Point Impedance",
    "Single Point Capacitance",
    "Quality Factor (Q) vs Frequency",
    "R-X (Resistance-Reactance) vs Frequency",
    "G-B (Conductance-Susceptance) vs Frequency",
    "Open/Short Correction",
]

# =============================
# GUI Parameter Definitions
# =============================
# Maps (instrument, measurement_idx) to [(label, key, default_value), ...]

MEASUREMENT_PARAMS = {
    # PIA Measurements (1-4, 7-10: freq sweep with optional DC bias)
    ("PIA", 1): [("Start Frequency (Hz)", "freq_start", "1000"),
                  ("Stop Frequency (Hz)", "freq_stop", "1000000"),
                  ("Number of Points", "num_points", "201"),
                  ("DC Bias (V)", "dc_bias_v", "0")],
    ("PIA", 2): [("Start Frequency (Hz)", "freq_start", "1000"),
                  ("Stop Frequency (Hz)", "freq_stop", "1000000"),
                  ("Number of Points", "num_points", "201"),
                  ("DC Bias (V)", "dc_bias_v", "0")],
    ("PIA", 3): [("Start Frequency (Hz)", "freq_start", "1000"),
                  ("Stop Frequency (Hz)", "freq_stop", "1000000"),
                  ("Number of Points", "num_points", "201"),
                  ("DC Bias (V)", "dc_bias_v", "0")],
    ("PIA", 4): [("Start Frequency (Hz)", "freq_start", "1000"),
                  ("Stop Frequency (Hz)", "freq_stop", "1000000"),
                  ("Number of Points", "num_points", "201"),
                  ("DC Bias (V)", "dc_bias_v", "0")],
    # PIA C-V Butterfly (5-6)
    ("PIA", 5): [("Measurement Frequency (Hz)", "freq_cv", "25000"),
                  ("Min Voltage (V)", "v_min", "0"),
                  ("Max Voltage (V)", "v_max", "5"),
                  ("Points per Sweep", "n_points", "401"),
                  ("Number of Cycles", "n_cycles", "1"),
                  ("HZO Thickness (nm)", "thickness_nm", "10.0"),
                  ("Electrode Diameter (µm)", "diam_um", "75.0")],
    ("PIA", 6): [("Measurement Frequency (Hz)", "freq_cv", "25000"),
                  ("Min Voltage (V)", "v_min", "0"),
                  ("Max Voltage (V)", "v_max", "5"),
                  ("Points per Sweep", "n_points", "401"),
                  ("Number of Cycles", "n_cycles", "1"),
                  ("HZO Thickness (nm)", "thickness_nm", "10.0"),
                  ("Electrode Diameter (µm)", "diam_um", "75.0")],
    # PIA Permittivity vs Frequency (7)
    ("PIA", 7): [("Start Frequency (Hz)", "freq_start", "1000"),
                  ("Stop Frequency (Hz)", "freq_stop", "1000000"),
                  ("Number of Points", "num_points", "201"),
                  ("DC Bias (V)", "dc_bias_v", "0"),
                  ("HZO Thickness (nm)", "thickness_nm", "10.0"),
                  ("Electrode Diameter (µm)", "diam_um", "75.0")],
    # PIA additional frequency sweeps (8-10)
    ("PIA", 8): [("Start Frequency (Hz)", "freq_start", "1000"),
                  ("Stop Frequency (Hz)", "freq_stop", "1000000"),
                  ("Number of Points", "num_points", "201"),
                  ("DC Bias (V)", "dc_bias_v", "0")],
    ("PIA", 9): [("Start Frequency (Hz)", "freq_start", "1000"),
                  ("Stop Frequency (Hz)", "freq_stop", "1000000"),
                  ("Number of Points", "num_points", "201"),
                  ("DC Bias (V)", "dc_bias_v", "0")],
    ("PIA", 10): [("Start Frequency (Hz)", "freq_start", "1000"),
                   ("Stop Frequency (Hz)", "freq_stop", "1000000"),
                   ("Number of Points", "num_points", "201"),
                   ("DC Bias (V)", "dc_bias_v", "0")],
    
    # PSPA Measurements
    ("PSPA", 1): [("Start Vds (V)", "vds_start", "0"),
                   ("Stop Vds (V)", "vds_stop", "5"),
                   ("Vds Step (V)", "vds_step", "0.1"),
                   ("Start Vgs (V)", "vgs_start", "0"),
                   ("Stop Vgs (V)", "vgs_stop", "3"),
                   ("Vgs Step (V)", "vgs_step", "0.5"),
                   ("Current Compliance (A)", "compliance", "0.1"),
                   ("Drain Channel", "drain_ch", "1"),
                   ("Gate Channel", "gate_ch", "2"),
                   ("Source Channel", "source_ch", "3")],
    ("PSPA", 2): [("Start Vgs (V)", "vgs_start", "-1"),
                   ("Stop Vgs (V)", "vgs_stop", "3"),
                   ("Vgs Step (V)", "vgs_step", "0.05"),
                   ("Constant Vds (V)", "vds_constant", "5"),
                   ("Current Compliance (A)", "compliance", "0.1"),
                   ("Drain Channel", "drain_ch", "1"),
                   ("Gate Channel", "gate_ch", "2"),
                   ("Source Channel", "source_ch", "3")],
    ("PSPA", 3): [("Start Voltage (V)", "v_start", "0"),
                   ("Stop Voltage (V)", "v_stop", "5"),
                   ("Voltage Step (V)", "v_step", "0.1"),
                   ("Current Compliance (A)", "compliance", "0.1"),
                   ("Channel Number", "channel", "1")],
    ("PSPA", 4): [("Max Voltage Magnitude (V)", "v_max", "5"),
                   ("Voltage Step (V)", "v_step", "0.1"),
                   ("Current Compliance (A)", "compliance", "0.1"),
                   ("Channel Number", "channel", "1")],
    ("PSPA", 5): [("Base Voltage (V)", "v_base", "0"),
                   ("Pulse Voltage (V)", "v_pulse", "5"),
                   ("Pulse Width (µs)", "pulse_width_us", "100"),
                   ("Pulse Period (ms)", "pulse_period_ms", "10"),
                   ("Number of Pulses", "num_pulses", "10"),
                   ("Current Compliance (A)", "compliance", "0.1"),
                   ("Channel Number", "channel", "1")],
    ("PSPA", 6): [("Base Vds (V)", "vds_base", "0"),
                   ("Pulse Vds (V)", "vds_pulse", "5"),
                   ("Base Vgs (V)", "vgs_base", "0"),
                   ("Pulse Vgs (V)", "vgs_pulse", "3"),
                   ("Pulse Width (µs)", "pulse_width_us", "100"),
                   ("Pulse Period (ms)", "pulse_period_ms", "10"),
                   ("Number of Pulses", "num_pulses", "10"),
                   ("Current Compliance (A)", "compliance", "0.1"),
                   ("Drain Channel", "drain_ch", "1"),
                   ("Gate Channel", "gate_ch", "2"),
                   ("Source Channel", "source_ch", "3")],
    ("PSPA", 7): [("Start Voltage (V)", "v_start", "-1"),
                   ("Stop Voltage (V)", "v_stop", "1"),
                   ("Voltage Step (V)", "v_step", "0.02"),
                   ("Current Compliance (A)", "compliance", "0.1"),
                   ("Anode Channel", "anode_ch", "1"),
                   ("Cathode Channel", "cathode_ch", "2")],
    ("PSPA", 8): [("Start Vgs (V)", "vgs_start", "0"),
                   ("Stop Vgs (V)", "vgs_stop", "5"),
                   ("Vgs Step (V)", "vgs_step", "0.1"),
                   ("Current Compliance (A)", "compliance", "0.001"),
                   ("Gate Channel", "gate_ch", "2"),
                   ("Source Channel", "source_ch", "3")],
    ("PSPA", 9): [("Start Voltage (V)", "v_start", "0"),
                   ("Stop Voltage (V)", "v_stop", "40"),
                   ("Voltage Step (V)", "v_step", "0.5"),
                   ("Current Compliance (A)", "compliance", "0.001"),
                   ("Breakdown Threshold (A)", "threshold_a", "0.0001"),
                   ("Channel Number", "channel", "1")],
    ("PSPA", 10): [("Start Current (A)", "i_start", "0"),
                    ("Stop Current (A)", "i_stop", "0.001"),
                    ("Current Step (A)", "i_step", "0.0001"),
                    ("Voltage Compliance (V)", "compliance", "10"),
                    ("Channel Number", "channel", "1")],
    
    # LCR Measurements
    ("LCR", 1): [("Start Frequency (Hz)", "freq_start", "20"),
                  ("Stop Frequency (Hz)", "freq_stop", "2000000"),
                  ("Number of Points", "num_points", "201"),
                  ("Sweep Type (LIN/LOG)", "sweep_type", "LOG"),
                  ("AC Level (V)", "ac_level", "1.0"),
                  ("DC Bias (V)", "dc_bias_v", "0")],
    ("LCR", 2): [("Start Frequency (Hz)", "freq_start", "20"),
                  ("Stop Frequency (Hz)", "freq_stop", "2000000"),
                  ("Number of Points", "num_points", "201"),
                  ("Sweep Type (LIN/LOG)", "sweep_type", "LOG"),
                  ("Mode (CPD/CSRS/CPRP)", "mode", "CPD"),
                  ("AC Level (V)", "ac_level", "1.0"),
                  ("DC Bias (V)", "dc_bias_v", "0")],
    ("LCR", 3): [("Start Frequency (Hz)", "freq_start", "20"),
                  ("Stop Frequency (Hz)", "freq_stop", "2000000"),
                  ("Number of Points", "num_points", "201"),
                  ("Sweep Type (LIN/LOG)", "sweep_type", "LOG"),
                  ("Mode (LPQ/LSD)", "mode", "LPQ"),
                  ("AC Level (V)", "ac_level", "1.0")],
    ("LCR", 4): [("Measurement Frequency (Hz)", "freq", "1000"),
                  ("Min Voltage (V)", "v_min", "-5"),
                  ("Max Voltage (V)", "v_max", "5"),
                  ("Number of Points", "num_points", "201"),
                  ("AC Level (V)", "ac_level", "0.1")],
    ("LCR", 5): [("Measurement Frequency (Hz)", "freq", "1000"),
                  ("Min Voltage (V)", "v_min", "-5"),
                  ("Max Voltage (V)", "v_max", "5"),
                  ("Points per Direction", "num_points", "201"),
                  ("Number of Cycles", "num_cycles", "1"),
                  ("AC Level (V)", "ac_level", "0.1"),
                  ("Dielectric Thickness (nm)", "thickness_nm", "10"),
                  ("Electrode Diameter (µm)", "diameter_um", "75")],
    ("LCR", 6): [("Frequency (Hz)", "freq", "1000"),
                  ("AC Level (V)", "ac_level", "1.0"),
                  ("DC Bias (V)", "dc_bias_v", "0")],
    ("LCR", 7): [("Frequency (Hz)", "freq", "1000"),
                  ("Mode (CPD/CSRS/CPRP)", "mode", "CPD"),
                  ("AC Level (V)", "ac_level", "1.0"),
                  ("DC Bias (V)", "dc_bias_v", "0"),
                  ("Dielectric Thickness (nm)", "thickness_nm", "10"),
                  ("Electrode Diameter (µm)", "diameter_um", "75")],
    ("LCR", 8): [("Start Frequency (Hz)", "freq_start", "20"),
                  ("Stop Frequency (Hz)", "freq_stop", "2000000"),
                  ("Number of Points", "num_points", "201"),
                  ("Sweep Type (LIN/LOG)", "sweep_type", "LOG"),
                  ("Mode (CPQ/LPQ)", "mode", "CPQ"),
                  ("AC Level (V)", "ac_level", "1.0")],
    ("LCR", 9): [("Start Frequency (Hz)", "freq_start", "20"),
                  ("Stop Frequency (Hz)", "freq_stop", "2000000"),
                  ("Number of Points", "num_points", "201"),
                  ("Sweep Type (LIN/LOG)", "sweep_type", "LOG"),
                  ("AC Level (V)", "ac_level", "1.0")],
    ("LCR", 10): [("Start Frequency (Hz)", "freq_start", "20"),
                   ("Stop Frequency (Hz)", "freq_stop", "2000000"),
                   ("Number of Points", "num_points", "201"),
                   ("Sweep Type (LIN/LOG)", "sweep_type", "LOG"),
                   ("AC Level (V)", "ac_level", "1.0")],
    ("LCR", 11): [],  # Open/Short has no parameters
}

# =============================
# Measurement Aliases
# =============================
# Get measurement name by instrument and index
def get_measurement_name(instrument: str, idx: int) -> str:
    """Get measurement name by instrument and 1-indexed measurement number."""
    if instrument == "PIA":
        return PIA_MEASUREMENTS[idx - 1] if 1 <= idx <= len(PIA_MEASUREMENTS) else ""
    elif instrument == "PSPA":
        return PSPA_MEASUREMENTS[idx - 1] if 1 <= idx <= len(PSPA_MEASUREMENTS) else ""
    elif instrument == "LCR":
        return LCR_MEASUREMENTS[idx - 1] if 1 <= idx <= len(LCR_MEASUREMENTS) else ""
    return ""

def get_measurements_list(instrument: str) -> list[str]:
    """Get list of measurements for an instrument."""
    if instrument == "PIA":
        return PIA_MEASUREMENTS
    elif instrument == "PSPA":
        return PSPA_MEASUREMENTS
    elif instrument == "LCR":
        return LCR_MEASUREMENTS
    return []
