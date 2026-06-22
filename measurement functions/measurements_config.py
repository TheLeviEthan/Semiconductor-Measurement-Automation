"""
Filename: measurements_config.py
Author: Ethan Ruddell
Date: 2026-5-11
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

Parameter tuple formats
-----------------------
Regular text field (3-element tuple):
    (display_label, dict_key, default_value)

Multi-select checkbox field (4-element tuple, 3rd element is a list):
    (display_label, dict_key, [(key, display_label), ...], [default_keys])
"""

# =============================
# Measurement Lists
# =============================
# Each list below enumerates the measurements available for one instrument.
# The order matters: index 1 in the list corresponds to choice 1 in the menu.
# Names should be short and descriptive so they fit nicely in the GUI listbox.
#
# PIA note: the Agilent 4294A always measures TWO parameters simultaneously
# in a single sweep (Trace A and Trace B).  Menu items that previously ran the
# same physical sweep twice just to plot different traces have been merged:
#
#   Old 1 + 2  (both Z-θ sweeps)        → 1: Impedance (|Z| & Phase)
#   Old 3-6    (CPD freq sweep)          → 2: Dielectrics
#   Old 7      (CPD C-V butterfly)       → 3: C-V Butterfly Cycles
#   Old 8      R-X                       → 4
#   Old 9      G-B                       → 5
#   Old 10     Y-θ                       → 6

PIA_MEASUREMENTS = [
    "Impedance (|Z| & Phase) vs Frequency",
    "Dielectrics",
    "C-V Butterfly Cycles",
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
    "Keithley Sourcemeter Transistors",
]

# The PSPA menu order shown to users maps directly to execution indices.
PSPA_CHOICE_MAP = {i: i for i in range(1, len(PSPA_MEASUREMENTS) + 1)}

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
# and the value is a list of parameter tuples (see module docstring for formats).

PIA_OSCILLATOR_FIELD = ("Oscillation Voltage (V)", "osc_voltage_v", "0.5")
PIA_IMPEDANCE_OSCILLATOR_FIELD = ("Impedance Oscillation Voltage (V)", "osc_voltage_v", "0.1")

MEASUREMENT_HELP = {
    # ------------------------------------------------------------------
    # PIA measurements
    # ------------------------------------------------------------------
    ("PIA", 1): (
        "Impedance sweep — measures |Z| (magnitude) and θ (phase angle) simultaneously "
        "in a single Z-θ sweep. Both are always saved to the CSV; use the checkboxes to "
        "choose which graphs to display. Z_real and Z_imag are derived from |Z| and θ."
    ),
    ("PIA", 2): (
        "Dielectric frequency characterization using CPD mode. The Agilent 4294A "
        "measures Cp (parallel capacitance) and D (dissipation factor / tan δ) "
        "SIMULTANEOUSLY in a single frequency sweep.\n\n"
        "Saves Cp(f), D(f), εr(f), εr″(f), and Rp(f) to CSV.\n"
        "Select which graphs to display with the checkboxes."
    ),
    ("PIA", 3): (
        "C-V butterfly measurement using CPD mode. The DC bias is swept "
        "V_min→V_max→V_min at a fixed AC frequency for the requested number of cycles. "
        "Cp(V), D(V), and εr(V) are recorded for each cycle, revealing ferroelectric "
        "hysteresis and polarization switching in HZO FeCAPs."
    ),
    ("PIA", 4): (
        "R-X sweep — measures resistance R and reactance X simultaneously in a single "
        "RX-mode sweep. Both are saved to the CSV and plotted."
    ),
    ("PIA", 5): (
        "G-B sweep — measures conductance G and susceptance B simultaneously in a single "
        "GB-mode sweep. Both are saved to the CSV and plotted."
    ),
    ("PIA", 6): (
        "Admittance sweep — measures admittance magnitude |Y| and phase θ simultaneously "
        "in a single YTD-mode sweep. Both are saved to the CSV and plotted."
    ),
    # ------------------------------------------------------------------
    # PSPA measurements
    # ------------------------------------------------------------------
    ("PSPA", 1): "Sweeps drain voltage at several gate voltages to produce transistor output characteristics.",
    ("PSPA", 2): "Sweeps gate voltage at a constant drain voltage to show transfer characteristics; choose Linear, Log, or Both using the 'Measurement Scale' option.",
    ("PSPA", 3): "Measures a unidirectional I-V curve over the requested voltage range.",
    ("PSPA", 4): "Measures a bidirectional I-V curve so forward and reverse sweeps can be compared.",
    ("PSPA", 5): "Runs a pulsed I-V test on a single device to reduce heating during the measurement.",
    ("PSPA", 6): "Runs a pulsed transistor measurement with drain and gate pulsing.",
    ("PSPA", 7): "Measures diode I-V behavior across a voltage sweep.",
    ("PSPA", 8): "Measures gate leakage current to see how much current flows through the gate stack.",
    ("PSPA", 9): "Sweeps voltage until the current reaches the breakdown threshold.",
    ("PSPA", 10): "Measures resistance by forcing current and reading the resulting voltage drop.",
    ("PSPA", 11): "Keithley sourcemeter transfer measurement with user-selected drain and source SMU channels and both linear and logarithmic plots.",
    # ------------------------------------------------------------------
    # LCR measurements
    # ------------------------------------------------------------------
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
    "freq_cv": "Fixed AC frequency used during the C-V butterfly measurement.",
    "freq": "Measurement frequency used for a single-point or fixed-frequency sweep.",
    "v_min": "Starting voltage for the bias sweep.",
    "v_max": "Ending voltage for the bias sweep.",
    "v_step": "Voltage increment between sweep points.",
    "n_points": "Number of points in each C-V sweep direction.",
    "n_cycles": "Number of complete up/down C-V cycles to repeat.",
    "thickness_nm": "Dielectric film thickness used when converting capacitance to permittivity. Set to 0 to skip permittivity calculations.",
    "area_um2": "Electrode area used when converting capacitance to permittivity. Set to 0 to skip permittivity calculations.",
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
    "sweep_ch": "PSPA channel used to force the voltage sweep.",
    "settle_s": "Delay after each voltage step to allow the current reading to settle.",
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
    "display_graphs": (
        "Select which graphs to generate and display. All collected data is always "
        "saved to the CSV file regardless of this selection."
    ),
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

# ---------------------------------------------------------------------------
# Multi-select graph options
# ---------------------------------------------------------------------------
_PIA_IMPEDANCE_GRAPHS = [
    ("z_mag",   "|Z| (Impedance Magnitude) vs Frequency"),
    ("z_phase", "Phase (θ) vs Frequency"),
    ("z_real",  "Z_real (Resistance) vs Frequency"),
    ("z_imag",  "Z_imag (Reactance) vs Frequency"),
]
_PIA_IMPEDANCE_GRAPH_DEFAULTS = ["z_mag", "z_phase"]

_PIA_DIELECTRICS_GRAPHS = [
    ("eps_vs_freq",      "εr (Permittivity) vs Frequency"),
    ("cp_vs_freq",       "Cp (Capacitance) vs Frequency"),
    ("d_vs_freq",        "D (Tan Loss) vs Frequency"),
    ("eps_imag_vs_freq", "εr″ (Imaginary Permittivity) vs Frequency"),
]
_PIA_DIELECTRICS_GRAPH_DEFAULTS = ["eps_vs_freq", "d_vs_freq"]

_PIA_CV_BUTTERFLY_GRAPHS = [
    ("cv_butterfly_cp",  "C-V Butterfly: Cp vs Voltage"),
    ("cv_butterfly_eps", "C-V Butterfly: εr vs Voltage"),
    ("cv_butterfly_d",   "C-V Butterfly: D vs Voltage"),
]
_PIA_CV_BUTTERFLY_GRAPH_DEFAULTS = ["cv_butterfly_cp", "cv_butterfly_eps"]

_PIA_FREQ_FIELDS = [
    ("Start Frequency (Hz)", "freq_start", "1000"),
    ("Stop Frequency (Hz)",  "freq_stop",  "1000000"),
    ("Number of Points",     "num_points", "201"),
    ("DC Bias (V)",          "dc_bias_v",  "0"),
]

MEASUREMENT_PARAMS = {
    # ------------------------------------------------------------------
    # PIA 1: Impedance — |Z| & θ measured simultaneously (Z-θ mode)
    # ------------------------------------------------------------------
    ("PIA", 1): [
        *_PIA_FREQ_FIELDS,
        PIA_IMPEDANCE_OSCILLATOR_FIELD,
        ("Graphs to Display", "display_graphs",
         _PIA_IMPEDANCE_GRAPHS, _PIA_IMPEDANCE_GRAPH_DEFAULTS),
    ],

    # ------------------------------------------------------------------
    # PIA 2: Dielectrics — CPD frequency sweep only
    # ------------------------------------------------------------------
    ("PIA", 2): [
        ("Start Frequency (Hz)",         "freq_start",   "1000"),
        ("Stop Frequency (Hz)",          "freq_stop",    "1000000"),
        ("Number of Points (freq sweep)","num_points",   "201"),
        ("DC Bias During Freq Sweep (V)","dc_bias_v",    "0"),
        ("Dielectric Thickness (nm)",    "thickness_nm", "10.0"),
        ("Electrode Area (µm²)",         "area_um2",     "4418"),
        PIA_OSCILLATOR_FIELD,
        ("Graphs to Display", "display_graphs",
         _PIA_DIELECTRICS_GRAPHS, _PIA_DIELECTRICS_GRAPH_DEFAULTS),
    ],

    # ------------------------------------------------------------------
    # PIA 3: C-V Butterfly Cycles
    # ------------------------------------------------------------------
    ("PIA", 3): [
        ("C-V Measurement Frequency (Hz)","freq_cv",     "25000"),
        ("C-V Min Voltage (V)",           "v_min",        "0"),
        ("C-V Max Voltage (V)",           "v_max",        "5"),
        ("C-V Points per Sweep",          "n_points",     "401"),
        ("Number of C-V Cycles",          "n_cycles",     "1"),
        ("Dielectric Thickness (nm)",     "thickness_nm", "10.0"),
        ("Electrode Area (µm²)",          "area_um2",     "4418"),
        PIA_OSCILLATOR_FIELD,
        ("Graphs to Display", "display_graphs",
         _PIA_CV_BUTTERFLY_GRAPHS, _PIA_CV_BUTTERFLY_GRAPH_DEFAULTS),
    ],

    # ------------------------------------------------------------------
    # PIA 4-6: R-X, G-B, Y-θ
    # ------------------------------------------------------------------
    ("PIA", 4): [*_PIA_FREQ_FIELDS, PIA_OSCILLATOR_FIELD],
    ("PIA", 5): [*_PIA_FREQ_FIELDS, PIA_OSCILLATOR_FIELD],
    ("PIA", 6): [*_PIA_FREQ_FIELDS, PIA_OSCILLATOR_FIELD],

    # ------------------------------------------------------------------
    # PSPA Measurements (unchanged)
    # ------------------------------------------------------------------
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
                   ("Measurement Scale (Linear/Log/Both)", "plot_scale", "Both"),
                   ("Drain SMU Channel", "drain_ch", "2"),
                   ("Gate SMU Channel", "gate_ch", "3"),
                   ("Source SMU Channel", "source_ch", "1")],
    ("PSPA", 3): [("Start Voltage (V)", "v_start", "0"),
                   ("Stop Voltage (V)", "v_stop", "5"),
                   ("Voltage Step (V)", "v_step", "0.1"),
                   ("Current Compliance (A)", "compliance", "0.1"),
                   ("Force Channel (applied voltage)", "channel", "1"),
                   ("Ground/Common Channel (optional)", "ground_ch", ""),
                   ("Integration Time (SHOR/MED/LONG)", "integration_time", "MED")],
    ("PSPA", 4): [("Max Voltage Magnitude (V)", "v_max", "5"),
                   ("Voltage Step (V)", "v_step", "0.1"),
                   ("Current Compliance (A)", "compliance", "0.1"),
                   ("Force Channel (applied voltage)", "channel", "1"),
                   ("Ground/Common Channel (optional)", "ground_ch", ""),
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
                    ("Measurement Scale (Linear/Log/Both)", "plot_scale", "Both"),
                    ("Drain SMU Channel", "drain_ch", "1"),
                    ("Source SMU Channel", "source_ch", "2")],

    # ------------------------------------------------------------------
    # LCR Measurements (unchanged)
    # ------------------------------------------------------------------
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

def get_measurement_name(instrument: str, idx: int) -> str:
    """Return the human-readable name of a measurement."""
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
