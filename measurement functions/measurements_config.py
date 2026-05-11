"""
Filename: measurements_config.py
Author: Ethan Ruddell
Date: 2026-02-13
Description: Central catalogue of every measurement the software can perform.

This file is the SINGLE SOURCE OF TRUTH for:
  - The ordered list of measurement names for each instrument
    (PIA, PSPA, LCR).  Both the GUI drop-down and the CLI menu
    read from these lists, so they always stay in sync.
  - The parameter definitions that the GUI uses to build its
    dynamic entry fields (label text, dictionary key, default value).

When you add a new measurement to an instrument driver you should:
  1. Append its name to the relevant list (PIA / PSPA / LCR).
  2. Add a MEASUREMENT_PARAMS entry for (instrument, index) with
     the required parameters.
  3. Implement the actual measurement logic in the driver module
     and the GUI/CLI execution logic.
"""

# =============================
# Measurement Lists
# =============================
# Each list below enumerates the measurements available for one instrument.
# The order matters: index 1 in the list corresponds to choice 1 in the menu.
# Names should be short and descriptive so they fit nicely in the GUI listbox.

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
    "Transistor Transfer Characteristics (Linear)",
    "Transistor Transfer Characteristics (Log)",
    "I-V Curve (Unidirectional)",
    "I-V Curve (Bidirectional)",
    "Pulsed I-V (Single Device)",
    "Pulsed Transistor",
    "Diode I-V Characteristics",
    "Gate Leakage Current",
    "Breakdown Voltage",
    "Resistance Measurement",
]

# The PSPA menu order shown to users differs from historical execution
# indices used in the underlying measurement logic.
PSPA_CHOICE_MAP = {
    1: 1,
    2: 2,
    3: 11,
    4: 3,
    5: 4,
    6: 5,
    7: 6,
    8: 7,
    9: 8,
    10: 9,
    11: 10,
}

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
# This dictionary tells the GUI which input fields to show for each
# measurement.  The key is a tuple (instrument_name, measurement_index)
# and the value is a list of (display_label, dict_key, default_value)
# tuples.  When the user selects a measurement, the GUI reads this
# dictionary and dynamically creates text entry boxes for every parameter.
#
# Maps (instrument, measurement_idx) to [(label, key, default_value), ...]

PIA_OSCILLATOR_FIELD = ("Oscillator Level (Vrms)", "osc_voltage_v", "0.5")
PIA_IMPEDANCE_OSCILLATOR_FIELD = ("Oscillator Level (Vrms)", "osc_voltage_v", "0.1")

MEASUREMENT_HELP = {
    ("PIA", 1): "Measures impedance magnitude and phase versus frequency to show the device's impedance response.",
    ("PIA", 2): "Measures impedance phase versus frequency; useful for seeing whether the device behaves more capacitively or inductively.",
    ("PIA", 3): "Measures capacitance and dissipation factor versus frequency to show capacitance dispersion and loss.",
    ("PIA", 4): "Measures tan loss versus frequency, highlighting dielectric loss behavior as the frequency changes.",
    ("PIA", 5): "Runs a C-V butterfly sweep using capacitance. The up/down bias sweep reveals hysteresis and dielectric response.",
    ("PIA", 6): "Runs a C-V butterfly sweep but converts capacitance to relative permittivity so the dielectric response is easier to compare.",
    ("PIA", 7): "Measures relative permittivity versus frequency at a fixed DC bias to show dielectric dispersion.",
    ("PIA", 8): "Measures resistance and reactance versus frequency to separate resistive and reactive behavior.",
    ("PIA", 9): "Measures conductance and susceptance versus frequency to show leakage and reactive admittance components.",
    ("PIA", 10): "Measures admittance magnitude and phase versus frequency as another view of the device's AC response.",
    ("PSPA", 1): "Sweeps drain voltage at several gate voltages to produce transistor output characteristics.",
    ("PSPA", 2): "Sweeps gate voltage at a constant drain voltage to show transfer characteristics on a linear scale.",
    ("PSPA", 3): "Sweeps gate voltage at a constant drain voltage to show transfer characteristics on a logarithmic scale.",
    ("PSPA", 4): "Measures a unidirectional I-V curve over the requested voltage range.",
    ("PSPA", 5): "Measures a bidirectional I-V curve so forward and reverse sweeps can be compared.",
    ("PSPA", 6): "Runs a pulsed I-V test on a single device to reduce heating during the measurement.",
    ("PSPA", 7): "Runs a pulsed transistor measurement with drain and gate pulsing.",
    ("PSPA", 8): "Measures diode I-V behavior across a voltage sweep.",
    ("PSPA", 9): "Measures gate leakage current to see how much current flows through the gate stack.",
    ("PSPA", 10): "Sweeps voltage until the current reaches the breakdown threshold.",
    ("PSPA", 11): "Measures resistance by forcing current and reading the resulting voltage drop.",
    ("LCR", 1): "Measures impedance magnitude and phase versus frequency with the LCR meter.",
    ("LCR", 2): "Measures capacitance and dissipation factor versus frequency with the LCR meter.",
    ("LCR", 3): "Measures inductance and quality factor versus frequency with the LCR meter.",
    ("LCR", 4): "Runs a single C-V sweep to show capacitance as the bias voltage changes.",
    ("LCR", 5): "Runs repeated C-V sweeps so hysteresis and cycle-to-cycle changes are visible.",
    ("LCR", 6): "Measures impedance at one frequency and optional DC bias point.",
    ("LCR", 7): "Measures capacitance at one frequency and optional DC bias point, with optional permittivity calculation.",
    ("LCR", 8): "Measures quality factor versus frequency to show how the sample losses change with frequency.",
    ("LCR", 9): "Measures resistance and reactance versus frequency with the LCR meter.",
    ("LCR", 10): "Measures conductance and susceptance versus frequency with the LCR meter.",
    ("LCR", 11): "Runs open/short correction to remove fixture parasitics before measurement.",
}

PARAM_HELP_OVERRIDES = {
    "freq_start": "First frequency in the sweep. Lower values capture low-frequency behavior; higher values capture high-frequency response.",
    "freq_stop": "Last frequency in the sweep.",
    "num_points": "Number of points in the sweep. More points give smoother curves but take longer.",
    "sweep_type": "Chooses linear or logarithmic spacing between sweep points.",
    "dc_bias_v": "Constant DC voltage applied during the measurement. Bias can change the measured response.",
    "osc_voltage_v": "AC test-signal amplitude for the PIA. Higher values usually strengthen the signal; lower values are gentler on the device.",
    "freq_cv": "Fixed frequency used during the C-V butterfly measurement.",
    "freq": "Measurement frequency used for a single-point or fixed-frequency sweep.",
    "v_min": "Starting voltage for the bias sweep.",
    "v_max": "Ending voltage for the bias sweep.",
    "v_step": "Voltage increment between sweep points.",
    "n_points": "Number of points in each sweep direction.",
    "n_cycles": "Number of complete up/down cycles to repeat.",
    "thickness_nm": "Dielectric film thickness used when converting capacitance to permittivity.",
    "area_um2": "Electrode area used when converting capacitance to permittivity.",
    "ac_level": "AC stimulus amplitude for the LCR meter. Larger values can improve signal strength, but may disturb sensitive devices.",
    "vds_start": "Starting drain-source voltage.",
    "vds_stop": "Ending drain-source voltage.",
    "vds_step": "Drain-source voltage increment.",
    "vgs_start": "Starting gate-source voltage.",
    "vgs_stop": "Ending gate-source voltage.",
    "vgs_step": "Gate-source voltage increment.",
    "vds_constant": "Drain-source voltage held constant during a transfer measurement.",
    "compliance": "Safety limit that prevents the instrument from sourcing more current or voltage than allowed.",
    "drain_ch": "Instrument channel assigned to the drain terminal.",
    "gate_ch": "Instrument channel assigned to the gate terminal.",
    "source_ch": "Instrument channel assigned to the source terminal.",
    "channel": "Instrument channel used for the measurement.",
    "anode_ch": "Instrument channel assigned to the diode anode.",
    "cathode_ch": "Instrument channel assigned to the diode cathode.",
    "sense_channel": "Sense channel used for the measurement. In 2-wire mode this is usually the same as the force channel.",
    "mode": "Selects the measurement function or mode used by the instrument.",
    "integration_time": "Sets the measurement aperture. Longer times are slower but usually less noisy.",
    "threshold_a": "Current threshold used to identify breakdown.",
    "i_start": "Starting current for a source-current sweep.",
    "i_stop": "Ending current for a source-current sweep.",
    "i_step": "Current increment for a source-current sweep.",
    "v_base": "Baseline voltage used before the pulse.",
    "v_pulse": "Voltage level during the pulse.",
    "pulse_width_us": "Duration of each pulse in microseconds.",
    "pulse_period_ms": "Time between pulses in milliseconds.",
    "num_pulses": "Number of pulses to apply.",
}


def get_measurement_help(instrument: str, idx: int) -> str:
    """Return a short help description for a measurement option."""
    return MEASUREMENT_HELP.get((instrument, idx), "No additional help text is available for this measurement.")


def get_parameter_help(key: str, label: str = "") -> str:
    """Return a short help description for a parameter field."""
    if key in PARAM_HELP_OVERRIDES:
        return PARAM_HELP_OVERRIDES[key]

    label_lower = label.lower()
    if "frequency" in label_lower and "start" in label_lower:
        return "First frequency in the sweep."
    if "frequency" in label_lower and "stop" in label_lower:
        return "Last frequency in the sweep."
    if "points" in label_lower:
        return "Number of points used for the sweep."
    if "voltage" in label_lower and "start" in label_lower:
        return "Starting voltage for the sweep."
    if "voltage" in label_lower and "stop" in label_lower:
        return "Ending voltage for the sweep."
    if "bias" in label_lower:
        return "DC bias applied during the measurement."
    if "area" in label_lower:
        return "Electrode area used in the permittivity calculation."
    if "thickness" in label_lower:
        return "Dielectric thickness used in the permittivity calculation."
    return "Hover help for this parameter."

MEASUREMENT_PARAMS = {
    # PIA measurements always expose the oscillator level because it controls
    # the AC test-signal amplitude used by the analyzer.
    ("PIA", 1): [("Start Frequency (Hz)", "freq_start", "1000"),
                  ("Stop Frequency (Hz)", "freq_stop", "1000000"),
                  ("Number of Points", "num_points", "201"),
                  ("DC Bias (V)", "dc_bias_v", "0"),
                  PIA_OSCILLATOR_FIELD],
    ("PIA", 2): [("Start Frequency (Hz)", "freq_start", "1000"),
                  ("Stop Frequency (Hz)", "freq_stop", "1000000"),
                  ("Number of Points", "num_points", "201"),
                  ("DC Bias (V)", "dc_bias_v", "0"),
                  PIA_IMPEDANCE_OSCILLATOR_FIELD],
    ("PIA", 3): [("Start Frequency (Hz)", "freq_start", "1000"),
                  ("Stop Frequency (Hz)", "freq_stop", "1000000"),
                  ("Number of Points", "num_points", "201"),
                  ("DC Bias (V)", "dc_bias_v", "0"),
                  PIA_OSCILLATOR_FIELD],
    ("PIA", 4): [("Start Frequency (Hz)", "freq_start", "1000"),
                  ("Stop Frequency (Hz)", "freq_stop", "1000000"),
                  ("Number of Points", "num_points", "201"),
                  ("DC Bias (V)", "dc_bias_v", "0"),
                  PIA_OSCILLATOR_FIELD],
    # PIA C-V Butterfly (5-6)
    ("PIA", 5): [("Measurement Frequency (Hz)", "freq_cv", "25000"),
                  ("Min Voltage (V)", "v_min", "0"),
                  ("Max Voltage (V)", "v_max", "5"),
                  ("Points per Sweep", "n_points", "401"),
                  ("Number of Cycles", "n_cycles", "1"),
                  ("HZO Thickness (nm)", "thickness_nm", "10.0"),
                  ("Electrode Area (µm²)", "area_um2", "4418"),
                  PIA_OSCILLATOR_FIELD],
    ("PIA", 6): [("Measurement Frequency (Hz)", "freq_cv", "25000"),
                  ("Min Voltage (V)", "v_min", "0"),
                  ("Max Voltage (V)", "v_max", "5"),
                  ("Points per Sweep", "n_points", "401"),
                  ("Number of Cycles", "n_cycles", "1"),
                  ("HZO Thickness (nm)", "thickness_nm", "10.0"),
                  ("Electrode Area (µm²)", "area_um2", "4418"),
                  PIA_OSCILLATOR_FIELD],
    # PIA Permittivity vs Frequency (7)
    ("PIA", 7): [("Start Frequency (Hz)", "freq_start", "1000"),
                  ("Stop Frequency (Hz)", "freq_stop", "1000000"),
                  ("Number of Points", "num_points", "201"),
                  ("DC Bias (V)", "dc_bias_v", "0"),
                  ("HZO Thickness (nm)", "thickness_nm", "10.0"),
                  ("Electrode Area (µm²)", "area_um2", "4418"),
                  PIA_OSCILLATOR_FIELD],
    # PIA additional frequency sweeps (8-10)
    ("PIA", 8): [("Start Frequency (Hz)", "freq_start", "1000"),
                  ("Stop Frequency (Hz)", "freq_stop", "1000000"),
                  ("Number of Points", "num_points", "201"),
                  ("DC Bias (V)", "dc_bias_v", "0"),
                  PIA_OSCILLATOR_FIELD],
    ("PIA", 9): [("Start Frequency (Hz)", "freq_start", "1000"),
                  ("Stop Frequency (Hz)", "freq_stop", "1000000"),
                  ("Number of Points", "num_points", "201"),
                  ("DC Bias (V)", "dc_bias_v", "0"),
                  PIA_OSCILLATOR_FIELD],
    ("PIA", 10): [("Start Frequency (Hz)", "freq_start", "1000"),
                   ("Stop Frequency (Hz)", "freq_stop", "1000000"),
                   ("Number of Points", "num_points", "201"),
                   ("DC Bias (V)", "dc_bias_v", "0"),
                   PIA_OSCILLATOR_FIELD],
    
    # PSPA Measurements
    ("PSPA", 1): [("Start Vds (V)", "vds_start", "0"),
                   ("Stop Vds (V)", "vds_stop", "5"),
                   ("Vds Step (V)", "vds_step", "0.1"),
                   ("Start Vgs (V)", "vgs_start", "0"),
                   ("Stop Vgs (V)", "vgs_stop", "3"),
                   ("Vgs Step (V)", "vgs_step", "0.5"),
                   ("Current Compliance (A)", "compliance", "0.1"),
                   ("Integration Time (SHOR/MED/LONG)", "integration_time", "MED"),
                   ("Drain SMU Channel", "drain_ch", "1"),
                   ("Gate SMU Channel", "gate_ch", "2"),
                   ("Source SMU Channel", "source_ch", "3")],
    ("PSPA", 2): [("Start Vgs (V)", "vgs_start", "-1"),
                   ("Stop Vgs (V)", "vgs_stop", "3"),
                   ("Vgs Step (V)", "vgs_step", "0.05"),
                   ("Constant Vds (V)", "vds_constant", "5"),
                   ("Current Compliance (A)", "compliance", "0.1"),
                   ("Integration Time (SHOR/MED/LONG)", "integration_time", "MED"),
                   ("Drain SMU Channel", "drain_ch", "2"),
                   ("Gate SMU Channel", "gate_ch", "3"),
                   ("Source SMU Channel", "source_ch", "1")],
    ("PSPA", 3): [("Start Voltage (V)", "v_start", "0"),
                   ("Stop Voltage (V)", "v_stop", "5"),
                   ("Voltage Step (V)", "v_step", "0.1"),
                   ("Current Compliance (A)", "compliance", "0.1"),
                   ("Channel Number", "channel", "1"),
                   ("Integration Time (SHOR/MED/LONG)", "integration_time", "MED")],
    ("PSPA", 4): [("Max Voltage Magnitude (V)", "v_max", "5"),
                   ("Voltage Step (V)", "v_step", "0.1"),
                   ("Current Compliance (A)", "compliance", "0.1"),
                   ("Channel Number", "channel", "1"),
                   ("Integration Time (SHOR/MED/LONG)", "integration_time", "MED")],
    ("PSPA", 5): [("Base Voltage (V)", "v_base", "0"),
                   ("Pulse Voltage (V)", "v_pulse", "5"),
                   ("Pulse Width (µs)", "pulse_width_us", "100"),
                   ("Pulse Period (ms)", "pulse_period_ms", "10"),
                   ("Number of Pulses", "num_pulses", "10"),
                   ("Current Compliance (A)", "compliance", "0.1"),
                   ("Integration Time (SHOR/MED/LONG)", "integration_time", "MED"),
                   ("Channel Number", "channel", "1")],
    ("PSPA", 6): [("Base Vds (V)", "vds_base", "0"),
                   ("Pulse Vds (V)", "vds_pulse", "5"),
                   ("Base Vgs (V)", "vgs_base", "0"),
                   ("Pulse Vgs (V)", "vgs_pulse", "3"),
                   ("Pulse Width (µs)", "pulse_width_us", "100"),
                   ("Pulse Period (ms)", "pulse_period_ms", "10"),
                   ("Number of Pulses", "num_pulses", "10"),
                   ("Current Compliance (A)", "compliance", "0.1"),
                   ("Integration Time (SHOR/MED/LONG)", "integration_time", "MED"),
                   ("Drain SMU Channel", "drain_ch", "1"),
                   ("Gate SMU Channel", "gate_ch", "2"),
                   ("Source SMU Channel", "source_ch", "3")],
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
                   ("Gate SMU Channel", "gate_ch", "2"),
                   ("Source SMU Channel", "source_ch", "3")],
    ("PSPA", 9): [("Start Voltage (V)", "v_start", "0"),
                   ("Stop Voltage (V)", "v_stop", "40"),
                   ("Voltage Step (V)", "v_step", "0.5"),
                   ("Current Compliance (A)", "compliance", "0.001"),
                   ("Breakdown Threshold (A)", "threshold_a", "0.0001"),
                   ("Channel Number", "channel", "1")],
    ("PSPA", 10): [("Start Current (A)", "i_start", "0"),
                    ("Stop Current (A)", "i_stop", "0.0001"),
                    ("Current Step (A)", "i_step", "0.00001"),
                    ("Voltage Compliance (V)", "compliance", "10"),
                    ("Channel Number", "channel", "1"),
                    ("Sense Channel (same as force for 2-wire)", "sense_channel", "1")],
    ("PSPA", 11): [("Start Vgs (V)", "vgs_start", "-1"),
                    ("Stop Vgs (V)", "vgs_stop", "3"),
                    ("Vgs Step (V)", "vgs_step", "0.05"),
                    ("Constant Vds (V)", "vds_constant", "5"),
                    ("Current Compliance (A)", "compliance", "0.1"),
                    ("Integration Time (SHOR/MED/LONG)", "integration_time", "MED"),
                    ("Drain SMU Channel", "drain_ch", "2"),
                    ("Gate SMU Channel", "gate_ch", "3"),
                    ("Source SMU Channel", "source_ch", "1")],
    
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
                  ("Electrode Area (µm²)", "area_um2", "4418")],
    ("LCR", 6): [("Frequency (Hz)", "freq", "1000"),
                  ("AC Level (V)", "ac_level", "1.0"),
                  ("DC Bias (V)", "dc_bias_v", "0")],
    ("LCR", 7): [("Frequency (Hz)", "freq", "1000"),
                  ("Mode (CPD/CSRS/CPRP)", "mode", "CPD"),
                  ("AC Level (V)", "ac_level", "1.0"),
                  ("DC Bias (V)", "dc_bias_v", "0"),
                  ("Dielectric Thickness (nm)", "thickness_nm", "10"),
                  ("Electrode Area (µm²)", "area_um2", "4418")],
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
# Helper Functions
# =============================
# Convenience functions so the rest of the code doesn't have to know the
# internal structure of the lists above.

# Get measurement name by instrument and index
def get_measurement_name(instrument: str, idx: int) -> str:
    """Return the human-readable name of a measurement.

    Args:
        instrument: "PIA", "PSPA", or "LCR"
        idx:        1-indexed measurement number (matches the menu / listbox)
    """
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


def normalize_pspa_choice(choice: int) -> int:
    """Translate displayed PSPA menu index to execution index."""
    return PSPA_CHOICE_MAP.get(choice, choice)
