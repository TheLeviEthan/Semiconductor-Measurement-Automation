"""
Filename: main.py
Author: Ethan Ruddell
Date: 2026-2-12
Description: Implements CLI for automating semiconductor measurements using
the Keysight 4294A Precision Impedance Analyzer, the Agilent 4155C/4156C
Semiconductor Parameter Analyzer, and the Keysight E4980A LCR Meter.
"""

import os
import sys
import argparse
import logging
import numpy as np
import matplotlib.pyplot as plt

# Add directories to path
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'utility'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'measurement functions'))

# Import modules
import pia
import pspa
import lcr
import file_management
import config
import logging_config
from gpib_utils import (
    InstrumentSession, prompt_float, prompt_int, prompt_choice, prompt_bool,
)

log = logging.getLogger(__name__)

# =============================
# User settings and constants
# =============================

pia_measurements = [
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

pspa_measurements = [
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

lcr_measurements = [
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


def parse_args():
    """Parse CLI arguments for the measurement automation script."""
    parser = argparse.ArgumentParser(description="Automate semiconductor measurements.")
    parser.add_argument(
        "-o",
        "--output-dir",
        default=None,
        help=f"Output directory for results (default: {file_management.default_output_dir})",
    )
    return parser.parse_args()

def main():
    args = parse_args()

    # --- Logging ----------------------------------------------------------
    logging_config.setup(output_dir=file_management.default_output_dir)
    log.info("Welcome to the NRG Semiconductor Measurement Automation")

    print("Welcome to the NRG Semiconductor Measurement Automation...\n")
    
    # Prompt for output directory once at the start
    default_dir = args.output_dir or config.get("output", "directory", "") or file_management.default_output_dir
    output_path = input(f"Enter output directory [default: {default_dir}]: ").strip() or default_dir
    file_management.set_output_dir(output_path)
    file_management.ensure_output_dir(file_management.output_dir)
    log.info("Output directory set to: %s", file_management.output_dir)
    print(f"Output directory set to: {file_management.output_dir}\n")
    
    while True:
        # Prompt user to select instrument or quit
        print("="*60)
        print("SELECT INSTRUMENT")
        print("="*60)
        print("1. PIA (Precision Impedance Analyzer)")
        print("2. PSPA (Parameter/Source Analyzer)")
        print("3. LCR (E4980A LCR Meter)")
        print("q. Quit")
        
        instrument_choice = input("\nEnter your choice (1, 2, 3, or q): ").strip().lower()
        
        if instrument_choice == 'q':
            log.info("User exited program")
            print("Exiting the program. Goodbye!")
            break
        elif instrument_choice == '1':
            run_pia_measurements()
        elif instrument_choice == '2':
            run_pspa_measurements()
        elif instrument_choice == '3':
            run_lcr_measurements()
        else:
            print("Invalid choice. Please try again.")


def run_pia_measurements():
    """Run PIA measurement menu."""
    while True:
        print("\n" + "="*60)
        print("PIA MEASUREMENTS")
        print("="*60)
        for i, meas in enumerate(pia_measurements):
            print(f"{i + 1}. {meas}")
        print("b. Back to instrument selection")
        
        choice = input(f"\nEnter your choice (1-{len(pia_measurements)} or b): ").strip().lower()
        
        if choice == 'b':
            break
        
        try:
            choice_num = int(choice)
            if 1 <= choice_num <= len(pia_measurements):
                execute_pia_measurement(choice_num)
            else:
                print("Invalid choice. Please try again.")
        except ValueError:
            print("Invalid input. Please enter a number or 'b'.")


def _pia_get_params():
    """Shared prompt logic for PIA frequency + DC bias parameters."""
    if (pia.current_parameters["freq_start_hz"] != pia.FREQ_START_HZ
            or pia.current_parameters["apply_dc_bias"] != pia.APPLY_DC_BIAS):
        use_last_params = pia.prompt_for_parameter_duplication()
    else:
        use_last_params = False
    freq_start, freq_stop, num_points = pia.get_frequency_parameters(use_last=use_last_params)
    apply_dc_bias, dc_bias_v = pia.get_dc_bias_parameters(use_last=use_last_params)
    return freq_start, freq_stop, num_points, apply_dc_bias, dc_bias_v


def _pia_safe_close(inst):
    """Turn off DC bias and close PIA instrument."""
    try:
        inst.write("DCO OFF")
    except Exception:
        pass
    inst.close()


def _pia_cv_butterfly(plot_mode="cp"):
    """
    Shared C-V butterfly measurement for PIA choices 5 and 6.

    Args:
        plot_mode: "cp" to save Cp plot, "eps" to save εr plot.
    """
    label = "Cp vs Voltage" if plot_mode == "cp" else "Permittivity vs Voltage"
    print(f"C–V Butterfly Cycle: {label}")

    freq_cv = prompt_float("Enter measurement frequency (Hz)", 25000)
    v_min = prompt_float("Enter minimum voltage (V)", 0)
    v_max = prompt_float("Enter maximum voltage (V)", 5)
    n_points = prompt_int("Enter number of points per sweep", 401)
    n_cycles = prompt_int("Enter number of cycles", 1)
    thickness_nm = prompt_float("Enter HZO thickness (nm)", 10.0)
    diam_um = prompt_float("Enter electrode diameter (µm)", 75.0)

    with InstrumentSession(pia.setup, _pia_safe_close) as inst:
        pia.initialize_4294a_for_cpd(inst)

        all_v_cycles = []
        all_cp_cycles = []
        all_eps_cycles = []
        all_cycle_idx = []

        for cycle in range(n_cycles):
            print(f"Running C–V cycle {cycle + 1}/{n_cycles}...")
            v, cp = pia.measure_single_cv_cycle(inst, freq_cv, v_min, v_max, n_points)
            eps_r = pia.compute_eps_r(cp, thickness_nm, diam_um)

            all_v_cycles.append(v)
            all_cp_cycles.append(cp)
            all_eps_cycles.append(eps_r)
            all_cycle_idx.append(np.full_like(v, cycle + 1, dtype=int))

        # Flatten for CSV
        v_all = np.concatenate(all_v_cycles)
        cp_all = np.concatenate(all_cp_cycles)
        eps_all = np.concatenate(all_eps_cycles)
        cycles_all = np.concatenate(all_cycle_idx)

        csv_data = np.column_stack([cycles_all, v_all, cp_all, eps_all])
        file_management.save_csv(
            "cv_dielectric_cycles.csv", csv_data,
            "cycle_index, bias_V, Cp_F, eps_r"
        )

        if plot_mode == "cp":
            file_management.save_cycle_plot(
                "C-V Butterfly: Cp vs Voltage", "Bias Voltage (V)", "Cp (F)",
                all_v_cycles, all_cp_cycles, "cv_butterfly_cp.png",
            )
        else:
            file_management.save_cycle_plot(
                "C-V Butterfly: εr vs Voltage", "Bias Voltage (V)",
                "Dielectric constant εr",
                all_v_cycles, all_eps_cycles, "cv_butterfly_eps_r.png",
            )
        log.info("C–V %s measurement complete", label)
        print(f"C–V {label} measurement complete. Data saved to output folder.")


def execute_pia_measurement(choice):
    """Execute the selected PIA measurement."""
    repeating = True
    while repeating:
        if choice == 1:
            # Impedance Magnitude vs Frequency
            print("You selected Impedance Magnitude vs Frequency measurement.")
            freq_start, freq_stop, num_points, apply_dc_bias, dc_bias_v = _pia_get_params()

            with InstrumentSession(pia.setup, _pia_safe_close) as inst:
                pia.initialize_4294a_for_impedance(inst)
                pia.configure_dc_bias(inst, apply_dc_bias, dc_bias_v)
                freq_axis, z_mag, theta_deg = pia.measure_impedance_vs_freq(inst, freq_start, freq_stop, num_points)

                file_management.save_image(
                    "Impedance Magnitude vs Frequency", "Frequency (Hz)", freq_axis, "|Z| (Ohm)", z_mag,
                    APPLY_DC_BIAS=apply_dc_bias, DC_BIAS_V=dc_bias_v
                )

                z_real = z_mag * np.cos(np.deg2rad(theta_deg))
                z_imag = z_mag * np.sin(np.deg2rad(theta_deg))
                csv_data = np.column_stack([freq_axis, z_mag, theta_deg, z_real, z_imag])
                file_management.save_csv(
                    "impedance_data.csv", csv_data,
                    "frequency_Hz, Z_mag_ohm, theta_deg, Z_real_ohm, Z_imag_ohm"
                )

        elif choice == 2:
            # Impedance Phase vs Frequency
            print("You selected Impedance Phase vs Frequency measurement.")
            freq_start, freq_stop, num_points, apply_dc_bias, dc_bias_v = _pia_get_params()

            with InstrumentSession(pia.setup, _pia_safe_close) as inst:
                pia.initialize_4294a_for_impedance(inst)
                pia.configure_dc_bias(inst, apply_dc_bias, dc_bias_v)
                freq_axis, z_mag, theta_deg = pia.measure_impedance_vs_freq(inst, freq_start, freq_stop, num_points)

                file_management.save_image(
                    "Impedance Phase vs Frequency", "Frequency (Hz)", freq_axis, "Phase (Degrees)", theta_deg,
                    APPLY_DC_BIAS=apply_dc_bias, DC_BIAS_V=dc_bias_v
                )

                z_real = z_mag * np.cos(np.deg2rad(theta_deg))
                z_imag = z_mag * np.sin(np.deg2rad(theta_deg))
                csv_data = np.column_stack([freq_axis, z_mag, theta_deg, z_real, z_imag])
                file_management.save_csv(
                    "impedance_data.csv", csv_data,
                    "frequency_Hz, Z_mag_ohm, theta_deg, Z_real_ohm, Z_imag_ohm"
                )

        elif choice == 3:
            # Capacitance vs Frequency
            print("You selected Capacitance vs Frequency measurement.")
            freq_start, freq_stop, num_points, apply_dc_bias, dc_bias_v = _pia_get_params()

            with InstrumentSession(pia.setup, _pia_safe_close) as inst:
                pia.initialize_4294a_for_cpd(inst)
                pia.configure_dc_bias(inst, apply_dc_bias, dc_bias_v)
                freq_axis, cp_vals, d_vals = pia.measure_cpd_vs_freq(inst, freq_start, freq_stop, num_points)

                file_management.save_image(
                    "Capacitance vs Frequency", "Frequency (Hz)", freq_axis, "Cp (F)", cp_vals,
                    APPLY_DC_BIAS=apply_dc_bias, DC_BIAS_V=dc_bias_v
                )

                csv_data = np.column_stack([freq_axis, cp_vals, d_vals])
                file_management.save_csv(
                    "cpd_data.csv", csv_data,
                    "frequency_Hz, Cp_F, D"
                )

        elif choice == 4:
            # Tan Loss vs Frequency
            print("You selected Tan Loss vs Frequency measurement.")
            freq_start, freq_stop, num_points, apply_dc_bias, dc_bias_v = _pia_get_params()

            with InstrumentSession(pia.setup, _pia_safe_close) as inst:
                pia.initialize_4294a_for_cpd(inst)
                pia.configure_dc_bias(inst, apply_dc_bias, dc_bias_v)
                freq_axis, cp_vals, d_vals = pia.measure_cpd_vs_freq(inst, freq_start, freq_stop, num_points)

                file_management.save_image(
                    "Tan Loss vs Frequency", "Frequency (Hz)", freq_axis, "D (loss factor)", d_vals,
                    APPLY_DC_BIAS=apply_dc_bias, DC_BIAS_V=dc_bias_v
                )

                csv_data = np.column_stack([freq_axis, cp_vals, d_vals])
                file_management.save_csv(
                    "cpd_data.csv", csv_data,
                    "frequency_Hz, Cp_F, D"
                )

        elif choice == 5:
            _pia_cv_butterfly(plot_mode="cp")

        elif choice == 6:
            _pia_cv_butterfly(plot_mode="eps")

        elif choice == 7:
            # εr vs frequency measurement
            print("Permittivity vs Frequency Measurement")

            if (pia.current_parameters["freq_start_hz"] != pia.FREQ_START_HZ
                    or pia.current_parameters["apply_dc_bias"] != pia.APPLY_DC_BIAS):
                use_last_params = pia.prompt_for_parameter_duplication()
            else:
                use_last_params = False

            if use_last_params:
                freq_start = pia.current_parameters["freq_start_hz"]
                freq_stop = pia.current_parameters["freq_stop_hz"]
                freq_points = pia.current_parameters["num_points"]
                bias_voltage = pia.current_parameters["dc_bias_v"]
            else:
                freq_start = prompt_float("Enter start frequency (Hz)", 1000)
                freq_stop = prompt_float("Enter stop frequency (Hz)", 1e6)
                freq_points = prompt_int("Enter number of frequency points", 201)
                bias_voltage = prompt_float("Enter DC bias voltage (V)", 0)

                pia.current_parameters["freq_start_hz"] = freq_start
                pia.current_parameters["freq_stop_hz"] = freq_stop
                pia.current_parameters["num_points"] = freq_points
                pia.current_parameters["dc_bias_v"] = bias_voltage
                pia.current_parameters["apply_dc_bias"] = (bias_voltage != 0)

            thickness_nm = prompt_float("Enter HZO thickness (nm)", 10.0)
            diam_um = prompt_float("Enter electrode diameter (µm)", 75.0)

            with InstrumentSession(pia.setup, _pia_safe_close) as inst:
                pia.initialize_4294a_for_cpd(inst)

                print("Running εr vs frequency sweep...")
                freq_axis, cp_f = pia.measure_eps_vs_freq(inst, freq_start, freq_stop, freq_points, bias_voltage)
                eps_f = pia.compute_eps_r(cp_f, thickness_nm, diam_um)

                csv_data = np.column_stack([freq_axis, cp_f, eps_f])
                file_management.save_csv("eps_vs_freq.csv", csv_data, "frequency_Hz, Cp_F, eps_r")

                file_management.save_image(
                    "Permittivity vs Frequency", "Frequency (Hz)", freq_axis, "Dielectric constant εr", eps_f,
                    APPLY_DC_BIAS=(bias_voltage != 0), DC_BIAS_V=bias_voltage
                )
                print("Frequency sweep complete. Data and plot saved to output folder.")

        elif choice == 8:
            # R-X vs Frequency
            print("R-X (Resistance-Reactance) vs Frequency")
            freq_start, freq_stop, num_points, apply_dc_bias, dc_bias_v = _pia_get_params()

            with InstrumentSession(pia.setup, _pia_safe_close) as inst:
                pia.initialize_4294a_for_rx(inst)
                pia.configure_dc_bias(inst, apply_dc_bias, dc_bias_v)
                freq_axis, r_vals, x_vals = pia.measure_rx_vs_freq(inst, freq_start, freq_stop, num_points)

                file_management.save_image(
                    "Resistance vs Frequency", "Frequency (Hz)", freq_axis, "R (Ohm)", r_vals,
                    APPLY_DC_BIAS=apply_dc_bias, DC_BIAS_V=dc_bias_v
                )
                file_management.save_image(
                    "Reactance vs Frequency", "Frequency (Hz)", freq_axis, "X (Ohm)", x_vals,
                    APPLY_DC_BIAS=apply_dc_bias, DC_BIAS_V=dc_bias_v
                )

                csv_data = np.column_stack([freq_axis, r_vals, x_vals])
                file_management.save_csv("rx_data.csv", csv_data, "frequency_Hz, R_ohm, X_ohm")
                print("R-X measurement complete. Data saved to output folder.")

        elif choice == 9:
            # G-B vs Frequency
            print("G-B (Conductance-Susceptance) vs Frequency")
            freq_start, freq_stop, num_points, apply_dc_bias, dc_bias_v = _pia_get_params()

            with InstrumentSession(pia.setup, _pia_safe_close) as inst:
                pia.initialize_4294a_for_gb(inst)
                pia.configure_dc_bias(inst, apply_dc_bias, dc_bias_v)
                freq_axis, g_vals, b_vals = pia.measure_gb_vs_freq(inst, freq_start, freq_stop, num_points)

                file_management.save_image(
                    "Conductance vs Frequency", "Frequency (Hz)", freq_axis, "G (S)", g_vals,
                    APPLY_DC_BIAS=apply_dc_bias, DC_BIAS_V=dc_bias_v
                )
                file_management.save_image(
                    "Susceptance vs Frequency", "Frequency (Hz)", freq_axis, "B (S)", b_vals,
                    APPLY_DC_BIAS=apply_dc_bias, DC_BIAS_V=dc_bias_v
                )

                csv_data = np.column_stack([freq_axis, g_vals, b_vals])
                file_management.save_csv("gb_data.csv", csv_data, "frequency_Hz, G_S, B_S")
                print("G-B measurement complete. Data saved to output folder.")

        elif choice == 10:
            # Y-θ vs Frequency
            print("Y-θ (Admittance) vs Frequency")
            freq_start, freq_stop, num_points, apply_dc_bias, dc_bias_v = _pia_get_params()

            with InstrumentSession(pia.setup, _pia_safe_close) as inst:
                pia.initialize_4294a_for_ytd(inst)
                pia.configure_dc_bias(inst, apply_dc_bias, dc_bias_v)
                freq_axis, y_mag, y_theta = pia.measure_ytd_vs_freq(inst, freq_start, freq_stop, num_points)

                file_management.save_image(
                    "Admittance Magnitude vs Frequency", "Frequency (Hz)", freq_axis, "|Y| (S)", y_mag,
                    APPLY_DC_BIAS=apply_dc_bias, DC_BIAS_V=dc_bias_v
                )
                file_management.save_image(
                    "Admittance Phase vs Frequency", "Frequency (Hz)", freq_axis, "θ (Degrees)", y_theta,
                    APPLY_DC_BIAS=apply_dc_bias, DC_BIAS_V=dc_bias_v
                )

                csv_data = np.column_stack([freq_axis, y_mag, y_theta])
                file_management.save_csv("ytd_data.csv", csv_data, "frequency_Hz, Y_mag_S, theta_deg")
                print("Y-θ measurement complete. Data saved to output folder.")

        repeat_input = input("\nDo you want to perform another measurement? (y/n): ").strip().lower()
        if repeat_input != 'y':
            repeating = False


def run_pspa_measurements():
    """Run PSPA measurement menu."""
    while True:
        print("\n" + "="*60)
        print("PSPA MEASUREMENTS")
        print("="*60)
        for i, meas in enumerate(pspa_measurements):
            print(f"{i + 1}. {meas}")
        print("b. Back to instrument selection")
        
        choice = input(f"\nEnter your choice (1-{len(pspa_measurements)} or b): ").strip().lower()
        
        if choice == 'b':
            break
        
        try:
            choice_num = int(choice)
            if 1 <= choice_num <= len(pspa_measurements):
                execute_pspa_measurement(choice_num)
            else:
                print("Invalid choice. Please try again.")
        except ValueError:
            print("Invalid input. Please enter a number or 'b'.")


def execute_pspa_measurement(choice):
    """Execute the selected PSPA measurement."""
    repeating = True
    while repeating:
        if choice == 1:
            # PSPA Transistor Output Characteristics
            print("PSPA: Transistor Output Characteristics (Id-Vds curves)")

            vds_start = prompt_float("Enter start Vds (V)", 0)
            vds_stop = prompt_float("Enter stop Vds (V)", 5)
            vds_step = prompt_float("Enter Vds step (V)", 0.1)
            vgs_start = prompt_float("Enter start Vgs (V)", 0)
            vgs_stop = prompt_float("Enter stop Vgs (V)", 3)
            vgs_step = prompt_float("Enter Vgs step (V)", 0.5)
            compliance = prompt_float("Enter current compliance (A)", 0.1)
            drain_ch = prompt_int("Enter drain channel", 1)
            gate_ch = prompt_int("Enter gate channel", 2)
            source_ch = prompt_int("Enter source channel", 3)

            try:
                with InstrumentSession(pspa.connect_pspa, pspa.disconnect_pspa) as pspa_inst:
                    print("\nRunning transistor output characteristics sweep...")
                    data = pspa.measure_transistor_output_characteristics(
                        pspa_inst, vds_start, vds_stop, vds_step,
                        vgs_start, vgs_stop, vgs_step,
                        drain_ch, gate_ch, source_ch, compliance
                    )

                    csv_data = np.column_stack([data['Vds'], data['Vgs'], data['Id']])
                    file_management.save_csv("transistor_output_chars.csv", csv_data, "Vds_V, Vgs_V, Id_A")

                    plt.figure(figsize=(10, 6))
                    for vgs in np.unique(data['Vgs']):
                        mask = data['Vgs'] == vgs
                        plt.plot(data['Vds'][mask], data['Id'][mask] * 1e3, marker='o', label=f"Vgs = {vgs:.2f} V")
                    plt.xlabel('Vds (V)')
                    plt.ylabel('Id (mA)')
                    plt.title('Transistor Output Characteristics')
                    plt.legend()
                    plt.grid(True)
                    plt.tight_layout()
                    file_management.save_plot("transistor_output_chars.png")
                    plt.close()
                    print("Measurement complete. Data saved to output folder.")
            except Exception as e:
                log.error("PSPA output chars failed: %s", e)
                print(f"Error during measurement: {e}")

        elif choice == 2:
            # PSPA Transistor Transfer Characteristics
            print("PSPA: Transistor Transfer Characteristics (Id-Vgs curve)")

            vgs_start = prompt_float("Enter start Vgs (V)", -1)
            vgs_stop = prompt_float("Enter stop Vgs (V)", 3)
            vgs_step = prompt_float("Enter Vgs step (V)", 0.05)
            vds_constant = prompt_float("Enter constant Vds (V)", 5)
            compliance = prompt_float("Enter current compliance (A)", 0.1)
            drain_ch = prompt_int("Enter drain channel", 1)
            gate_ch = prompt_int("Enter gate channel", 2)
            source_ch = prompt_int("Enter source channel", 3)

            try:
                with InstrumentSession(pspa.connect_pspa, pspa.disconnect_pspa) as pspa_inst:
                    print("\nRunning transistor transfer characteristics sweep...")
                    data = pspa.measure_transistor_transfer_characteristics(
                        pspa_inst, vgs_start, vgs_stop, vgs_step, vds_constant,
                        drain_ch, gate_ch, source_ch, compliance
                    )

                    csv_data = np.column_stack([data['Vgs'], data['Id'], data['Ig']])
                    file_management.save_csv("transistor_transfer_chars.csv", csv_data, "Vgs_V, Id_A, Ig_A")

                    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
                    ax1.plot(data['Vgs'], data['Id'] * 1e3, marker='o', label='Id')
                    ax1.set_xlabel('Vgs (V)')
                    ax1.set_ylabel('Id (mA)')
                    ax1.set_title(f'Transfer Characteristics (Vds = {vds_constant} V) - Linear')
                    ax1.grid(True)
                    ax1.legend()

                    ax2.semilogy(data['Vgs'], np.abs(data['Id']), marker='o', label='|Id|')
                    ax2.semilogy(data['Vgs'], np.abs(data['Ig']), marker='s', label='|Ig|')
                    ax2.set_xlabel('Vgs (V)')
                    ax2.set_ylabel('Current (A)')
                    ax2.set_title(f'Transfer Characteristics (Vds = {vds_constant} V) - Log')
                    ax2.grid(True)
                    ax2.legend()

                    plt.tight_layout()
                    file_management.save_plot("transistor_transfer_chars.png")
                    plt.close()
                    print("Measurement complete. Data saved to output folder.")
            except Exception as e:
                log.error("PSPA transfer chars failed: %s", e)
                print(f"Error during measurement: {e}")

        elif choice == 3:
            # PSPA I-V Curve (Unidirectional)
            print("PSPA: I-V Curve (Unidirectional)")

            v_start = prompt_float("Enter start voltage (V)", 0)
            v_stop = prompt_float("Enter stop voltage (V)", 5)
            v_step = prompt_float("Enter voltage step (V)", 0.1)
            compliance = prompt_float("Enter current compliance (A)", 0.1)
            channel = prompt_int("Enter channel number", 1)

            try:
                with InstrumentSession(pspa.connect_pspa, pspa.disconnect_pspa) as pspa_inst:
                    print("\nRunning I-V measurement...")
                    data = pspa.measure_iv_curve(pspa_inst, v_start, v_stop, v_step, channel, compliance)

                    csv_data = np.column_stack([data['Voltage'], data['Current']])
                    file_management.save_csv("iv_curve.csv", csv_data, "Voltage_V, Current_A")

                    plt.figure(figsize=(10, 6))
                    plt.plot(data['Voltage'], data['Current'] * 1e3, marker='o')
                    plt.xlabel('Voltage (V)')
                    plt.ylabel('Current (mA)')
                    plt.title('I-V Curve (Unidirectional)')
                    plt.grid(True)
                    plt.tight_layout()
                    file_management.save_plot("iv_curve.png")
                    plt.close()
                    print("Measurement complete. Data saved to output folder.")
            except Exception as e:
                log.error("PSPA IV curve failed: %s", e)
                print(f"Error during measurement: {e}")

        elif choice == 4:
            # PSPA I-V Curve (Bidirectional)
            print("PSPA: I-V Curve (Bidirectional)")

            v_max = prompt_float("Enter maximum voltage magnitude (V)", 5)
            v_step = prompt_float("Enter voltage step (V)", 0.1)
            compliance = prompt_float("Enter current compliance (A)", 0.1)
            channel = prompt_int("Enter channel number", 1)

            try:
                with InstrumentSession(pspa.connect_pspa, pspa.disconnect_pspa) as pspa_inst:
                    print("\nRunning bidirectional I-V measurement...")
                    data = pspa.measure_iv_bidirectional(pspa_inst, v_max, v_step, channel, compliance)

                    csv_data = np.column_stack([data['Voltage'], data['Current']])
                    file_management.save_csv("iv_curve_bidirectional.csv", csv_data, "Voltage_V, Current_A")

                    plt.figure(figsize=(10, 6))
                    plt.plot(data['Voltage'], data['Current'] * 1e3, marker='o')
                    plt.xlabel('Voltage (V)')
                    plt.ylabel('Current (mA)')
                    plt.title('I-V Curve (Bidirectional)')
                    plt.grid(True)
                    plt.tight_layout()
                    file_management.save_plot("iv_curve_bidirectional.png")
                    plt.close()
                    print("Measurement complete. Data saved to output folder.")
            except Exception as e:
                log.error("PSPA bidirectional IV failed: %s", e)
                print(f"Error during measurement: {e}")

        elif choice == 5:
            # Pulsed I-V (Single Device)
            print("PSPA: Pulsed I-V (Single Device)")

            v_base = prompt_float("Enter base voltage (V)", 0)
            v_pulse = prompt_float("Enter pulse voltage (V)", 5)
            pulse_width = prompt_float("Enter pulse width (µs)", 100) * 1e-6
            pulse_period = prompt_float("Enter pulse period (ms)", 10) * 1e-3
            num_pulses = prompt_int("Enter number of pulses", 10)
            compliance = prompt_float("Enter current compliance (A)", 0.1)
            channel = prompt_int("Enter channel number", 1)

            try:
                with InstrumentSession(pspa.connect_pspa, pspa.disconnect_pspa) as pspa_inst:
                    print("\nRunning pulsed I-V measurement...")
                    data = pspa.measure_pulsed_iv(
                        pspa_inst, v_base, v_pulse, pulse_width, pulse_period, num_pulses,
                        channel, compliance
                    )

                    csv_data = np.column_stack([data['Time'], data['Voltage'], data['Current']])
                    file_management.save_csv("pulsed_iv.csv", csv_data, "Time_s, Voltage_V, Current_A")

                    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
                    ax1.plot(data['Time'] * 1e3, data['Current'] * 1e3, marker='o')
                    ax1.set_xlabel('Time (ms)')
                    ax1.set_ylabel('Current (mA)')
                    ax1.set_title('Pulsed I-V: Current Response')
                    ax1.grid(True)

                    ax2.plot(data['Time'] * 1e3, data['Voltage'], marker='o')
                    ax2.set_xlabel('Time (ms)')
                    ax2.set_ylabel('Voltage (V)')
                    ax2.set_title('Pulsed I-V: Applied Voltage')
                    ax2.grid(True)

                    plt.tight_layout()
                    file_management.save_plot("pulsed_iv.png")
                    plt.close()
                    print("Measurement complete. Data saved to output folder.")
            except Exception as e:
                log.error("PSPA pulsed IV failed: %s", e)
                print(f"Error during measurement: {e}")

        elif choice == 6:
            # Pulsed Transistor
            print("PSPA: Pulsed Transistor Measurement")

            vds_base = prompt_float("Enter base Vds (V)", 0)
            vds_pulse = prompt_float("Enter pulse Vds (V)", 5)
            vgs_base = prompt_float("Enter base Vgs (V)", 0)
            vgs_pulse = prompt_float("Enter pulse Vgs (V)", 3)
            pulse_width = prompt_float("Enter pulse width (µs)", 100) * 1e-6
            pulse_period = prompt_float("Enter pulse period (ms)", 10) * 1e-3
            num_pulses = prompt_int("Enter number of pulses", 10)
            compliance = prompt_float("Enter current compliance (A)", 0.1)
            drain_ch = prompt_int("Enter drain channel", 1)
            gate_ch = prompt_int("Enter gate channel", 2)
            source_ch = prompt_int("Enter source channel", 3)

            try:
                with InstrumentSession(pspa.connect_pspa, pspa.disconnect_pspa) as pspa_inst:
                    print("\nRunning pulsed transistor measurement...")
                    data = pspa.measure_pulsed_transistor(
                        pspa_inst, vds_pulse, vgs_pulse, vds_base, vgs_base,
                        pulse_width, pulse_period, num_pulses,
                        drain_ch, gate_ch, source_ch, compliance
                    )

                    csv_data = np.column_stack([data['Time'], data['Id']])
                    file_management.save_csv("pulsed_transistor.csv", csv_data, "Time_s, Id_A")

                    plt.figure(figsize=(10, 6))
                    plt.plot(data['Time'] * 1e3, data['Id'] * 1e3, marker='o', label='Id')
                    plt.xlabel('Time (ms)')
                    plt.ylabel('Id (mA)')
                    plt.title(f'Pulsed Transistor: Drain Current (Vds={vds_pulse}V, Vgs={vgs_pulse}V)')
                    plt.grid(True)
                    plt.legend()
                    plt.tight_layout()
                    file_management.save_plot("pulsed_transistor.png")
                    plt.close()
                    print("Measurement complete. Data saved to output folder.")
            except Exception as e:
                log.error("PSPA pulsed transistor failed: %s", e)
                print(f"Error during measurement: {e}")

        elif choice == 7:
            # Diode I-V Characteristics
            print("PSPA: Diode I-V Characteristics")

            v_start = prompt_float("Enter start voltage (V)", -1)
            v_stop = prompt_float("Enter stop voltage (V)", 1)
            v_step = prompt_float("Enter voltage step (V)", 0.02)
            compliance = prompt_float("Enter current compliance (A)", 0.1)
            anode_ch = prompt_int("Enter anode channel", 1)
            cathode_ch = prompt_int("Enter cathode channel", 2)

            try:
                with InstrumentSession(pspa.connect_pspa, pspa.disconnect_pspa) as pspa_inst:
                    print("\nRunning diode I-V measurement...")
                    data = pspa.measure_diode_iv(
                        pspa_inst, v_start, v_stop, v_step, anode_ch, cathode_ch, compliance
                    )

                    csv_data = np.column_stack([data['Voltage'], data['Current']])
                    file_management.save_csv("diode_iv.csv", csv_data, "Voltage_V, Current_A")

                    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
                    ax1.plot(data['Voltage'], data['Current'] * 1e3, marker='o')
                    ax1.set_xlabel('Voltage (V)')
                    ax1.set_ylabel('Current (mA)')
                    ax1.set_title('Diode I-V (Linear)')
                    ax1.grid(True)

                    ax2.semilogy(data['Voltage'], np.abs(data['Current']), marker='o')
                    ax2.set_xlabel('Voltage (V)')
                    ax2.set_ylabel('|Current| (A)')
                    ax2.set_title('Diode I-V (Log)')
                    ax2.grid(True)

                    plt.tight_layout()
                    file_management.save_plot("diode_iv.png")
                    plt.close()
                    print("Measurement complete. Data saved to output folder.")
            except Exception as e:
                log.error("PSPA diode IV failed: %s", e)
                print(f"Error during measurement: {e}")

        elif choice == 8:
            # Gate Leakage Current
            print("PSPA: Gate Leakage Current")

            vgs_start = prompt_float("Enter start Vgs (V)", 0)
            vgs_stop = prompt_float("Enter stop Vgs (V)", 5)
            vgs_step = prompt_float("Enter Vgs step (V)", 0.1)
            compliance = prompt_float("Enter current compliance (A)", 1e-3)
            gate_ch = prompt_int("Enter gate channel", 2)
            source_ch = prompt_int("Enter source channel", 3)

            try:
                with InstrumentSession(pspa.connect_pspa, pspa.disconnect_pspa) as pspa_inst:
                    print("\nRunning gate leakage measurement...")
                    data = pspa.measure_gate_leakage(
                        pspa_inst, vgs_start, vgs_stop, vgs_step, gate_ch, source_ch, compliance
                    )

                    csv_data = np.column_stack([data['Vgs'], data['Ig']])
                    file_management.save_csv("gate_leakage.csv", csv_data, "Vgs_V, Ig_A")

                    plt.figure(figsize=(10, 6))
                    plt.semilogy(data['Vgs'], np.abs(data['Ig']), marker='o')
                    plt.xlabel('Vgs (V)')
                    plt.ylabel('|Ig| (A)')
                    plt.title('Gate Leakage Current')
                    plt.grid(True)
                    plt.tight_layout()
                    file_management.save_plot("gate_leakage.png")
                    plt.close()
                    print("Measurement complete. Data saved to output folder.")
            except Exception as e:
                log.error("PSPA gate leakage failed: %s", e)
                print(f"Error during measurement: {e}")

        elif choice == 9:
            # Breakdown Voltage
            print("PSPA: Breakdown Voltage Measurement")

            v_start = prompt_float("Enter start voltage (V)", 0)
            v_stop = prompt_float("Enter stop voltage (V)", 40)
            v_step = prompt_float("Enter voltage step (V)", 0.5)
            compliance = prompt_float("Enter current compliance (A)", 1e-3)
            threshold = prompt_float("Enter breakdown current threshold (A)", 1e-4)
            channel = prompt_int("Enter channel number", 1)

            try:
                with InstrumentSession(pspa.connect_pspa, pspa.disconnect_pspa) as pspa_inst:
                    print("\nRunning breakdown voltage measurement...")
                    data = pspa.measure_breakdown_voltage(
                        pspa_inst, v_start, v_stop, v_step, channel, compliance, threshold
                    )

                    csv_data = np.column_stack([data['Voltage'], data['Current']])
                    file_management.save_csv("breakdown_voltage.csv", csv_data, "Voltage_V, Current_A")

                    plt.figure(figsize=(10, 6))
                    plt.semilogy(data['Voltage'], np.abs(data['Current']), marker='o')
                    plt.xlabel('Voltage (V)')
                    plt.ylabel('|Current| (A)')
                    if data['breakdown_v'] is not None:
                        plt.axvline(x=data['breakdown_v'], color='r', linestyle='--',
                                    label=f"Breakdown = {data['breakdown_v']:.2f} V")
                        plt.legend()
                    plt.title('Breakdown Voltage')
                    plt.grid(True)
                    plt.tight_layout()
                    file_management.save_plot("breakdown_voltage.png")
                    plt.close()

                    if data['breakdown_v'] is not None:
                        print(f"Breakdown voltage: {data['breakdown_v']:.3f} V")
                    else:
                        print("Breakdown threshold not reached in sweep range.")
                    print("Measurement complete. Data saved to output folder.")
            except Exception as e:
                log.error("PSPA breakdown failed: %s", e)
                print(f"Error during measurement: {e}")

        elif choice == 10:
            # Resistance Measurement
            print("PSPA: Resistance Measurement (force I, measure V)")

            i_start = prompt_float("Enter start current (A)", 0)
            i_stop = prompt_float("Enter stop current (A)", 1e-3)
            i_step = prompt_float("Enter current step (A)", 1e-4)
            compliance = prompt_float("Enter voltage compliance (V)", 10.0)
            channel = prompt_int("Enter channel number", 1)

            try:
                with InstrumentSession(pspa.connect_pspa, pspa.disconnect_pspa) as pspa_inst:
                    print("\nRunning resistance measurement...")
                    data = pspa.measure_resistance(
                        pspa_inst, i_start, i_stop, i_step, channel, compliance
                    )

                    csv_data = np.column_stack([data['Current'], data['Voltage']])
                    file_management.save_csv("resistance.csv", csv_data, "Current_A, Voltage_V")

                    plt.figure(figsize=(10, 6))
                    plt.plot(data['Current'] * 1e3, data['Voltage'], marker='o')
                    plt.xlabel('Current (mA)')
                    plt.ylabel('Voltage (V)')
                    plt.title(f"Resistance Measurement (R ≈ {data['resistance_ohm']:.4e} Ω)")
                    plt.grid(True)
                    plt.tight_layout()
                    file_management.save_plot("resistance.png")
                    plt.close()

                    print(f"\nMeasured Resistance: {data['resistance_ohm']:.4e} Ω")
                    print("Measurement complete. Data saved to output folder.")
            except Exception as e:
                log.error("PSPA resistance failed: %s", e)
                print(f"Error during measurement: {e}")

        repeat_input = input("\nDo you want to perform another measurement? (y/n): ").strip().lower()
        if repeat_input != 'y':
            repeating = False


def run_lcr_measurements():
    """Run LCR measurement menu."""
    while True:
        print("\n" + "="*60)
        print("LCR MEASUREMENTS (E4980A)")
        print("="*60)
        for i, meas in enumerate(lcr_measurements):
            print(f"{i + 1}. {meas}")
        print("b. Back to instrument selection")
        
        choice = input(f"\nEnter your choice (1-{len(lcr_measurements)} or b): ").strip().lower()
        
        if choice == 'b':
            break
        
        try:
            choice_num = int(choice)
            if 1 <= choice_num <= len(lcr_measurements):
                execute_lcr_measurement(choice_num)
            else:
                print("Invalid choice. Please try again.")
        except ValueError:
            print("Invalid input. Please enter a number or 'b'.")


def execute_lcr_measurement(choice):
    """Execute the selected LCR measurement."""
    repeating = True
    while repeating:
        if choice == 1:
            # Impedance vs Frequency
            print("LCR: Impedance vs Frequency")

            freq_start = prompt_float("Enter start frequency (Hz)", 20)
            freq_stop = prompt_float("Enter stop frequency (Hz)", 2e6)
            num_points = prompt_int("Enter number of points", 201)
            sweep_type = prompt_choice("Sweep type", ["LIN", "LOG"], "LOG")
            ac_level = prompt_float("Enter AC level (V)", 1.0)
            apply_bias = prompt_bool("Apply DC bias?", False)
            bias_v = 0.0
            if apply_bias:
                bias_v = prompt_float("Enter DC bias (V)", 0)

            try:
                with InstrumentSession(
                    lambda: lcr.setup(measurement_function="ZTD"),
                    lcr.disconnect_e4980a
                ) as lcr_inst:
                    lcr.set_ac_level(lcr_inst, ac_level)
                    lcr.configure_dc_bias(lcr_inst, apply_bias, bias_v)

                    print("\nRunning impedance vs frequency sweep...")
                    freq_axis, z_mag, theta_deg = lcr.measure_impedance_vs_frequency(
                        lcr_inst, freq_start, freq_stop, num_points, sweep_type
                    )

                    file_management.save_image(
                        "Impedance Magnitude vs Frequency (LCR)", "Frequency (Hz)", freq_axis,
                        "|Z| (Ohm)", z_mag, APPLY_DC_BIAS=apply_bias, DC_BIAS_V=bias_v
                    )
                    file_management.save_image(
                        "Impedance Phase vs Frequency (LCR)", "Frequency (Hz)", freq_axis,
                        "Phase (Degrees)", theta_deg, APPLY_DC_BIAS=apply_bias, DC_BIAS_V=bias_v
                    )

                    z_real, z_imag = lcr.compute_impedance_components(z_mag, theta_deg)
                    csv_data = np.column_stack([freq_axis, z_mag, theta_deg, z_real, z_imag])
                    file_management.save_csv(
                        "lcr_impedance_vs_freq.csv", csv_data,
                        "frequency_Hz, Z_mag_ohm, theta_deg, Z_real_ohm, Z_imag_ohm"
                    )
                    print("Measurement complete. Data saved to output folder.")
            except Exception as e:
                log.error("LCR impedance sweep failed: %s", e)
                print(f"Error during measurement: {e}")

        elif choice == 2:
            # Capacitance vs Frequency
            print("LCR: Capacitance vs Frequency")

            freq_start = prompt_float("Enter start frequency (Hz)", 20)
            freq_stop = prompt_float("Enter stop frequency (Hz)", 2e6)
            num_points = prompt_int("Enter number of points", 201)
            sweep_type = prompt_choice("Sweep type", ["LIN", "LOG"], "LOG")
            mode = prompt_choice("Measurement mode", ["CPD", "CSRS", "CPRP"], "CPD")
            ac_level = prompt_float("Enter AC level (V)", 1.0)
            apply_bias = prompt_bool("Apply DC bias?", False)
            bias_v = 0.0
            if apply_bias:
                bias_v = prompt_float("Enter DC bias (V)", 0)

            try:
                with InstrumentSession(
                    lambda: lcr.setup(measurement_function=mode),
                    lcr.disconnect_e4980a
                ) as lcr_inst:
                    lcr.set_ac_level(lcr_inst, ac_level)
                    lcr.configure_dc_bias(lcr_inst, apply_bias, bias_v)

                    print(f"\nRunning capacitance vs frequency sweep ({mode})...")
                    freq_axis, cap_vals, secondary_vals = lcr.measure_capacitance_vs_frequency(
                        lcr_inst, freq_start, freq_stop, num_points, sweep_type, mode
                    )

                    secondary_label = "D (loss factor)" if mode == "CPD" else "Secondary"

                    file_management.save_image(
                        f"Capacitance vs Frequency ({mode}) (LCR)", "Frequency (Hz)", freq_axis,
                        "Capacitance (F)", cap_vals, APPLY_DC_BIAS=apply_bias, DC_BIAS_V=bias_v
                    )
                    file_management.save_image(
                        f"{secondary_label} vs Frequency (LCR)", "Frequency (Hz)", freq_axis,
                        secondary_label, secondary_vals, APPLY_DC_BIAS=apply_bias, DC_BIAS_V=bias_v
                    )

                    csv_data = np.column_stack([freq_axis, cap_vals, secondary_vals])
                    file_management.save_csv(
                        f"lcr_capacitance_vs_freq_{mode.lower()}.csv", csv_data,
                        f"frequency_Hz, capacitance_F, {secondary_label.replace(' ', '_')}"
                    )
                    print("Measurement complete. Data saved to output folder.")
            except Exception as e:
                log.error("LCR capacitance sweep failed: %s", e)
                print(f"Error during measurement: {e}")

        elif choice == 3:
            # Inductance vs Frequency
            print("LCR: Inductance vs Frequency")

            freq_start = prompt_float("Enter start frequency (Hz)", 20)
            freq_stop = prompt_float("Enter stop frequency (Hz)", 2e6)
            num_points = prompt_int("Enter number of points", 201)
            sweep_type = prompt_choice("Sweep type", ["LIN", "LOG"], "LOG")
            mode = prompt_choice("Measurement mode", ["LPQ", "LSD"], "LPQ")
            ac_level = prompt_float("Enter AC level (V)", 1.0)

            try:
                with InstrumentSession(
                    lambda: lcr.setup(measurement_function=mode),
                    lcr.disconnect_e4980a
                ) as lcr_inst:
                    lcr.set_ac_level(lcr_inst, ac_level)

                    print(f"\nRunning inductance vs frequency sweep ({mode})...")
                    freq_axis, ind_vals, secondary_vals = lcr.measure_frequency_sweep(
                        lcr_inst, freq_start, freq_stop, num_points, sweep_type, mode
                    )

                    secondary_label = "Q (quality factor)" if mode == "LPQ" else "D (loss factor)"

                    file_management.save_image(
                        f"Inductance vs Frequency ({mode}) (LCR)", "Frequency (Hz)", freq_axis,
                        "Inductance (H)", ind_vals, APPLY_DC_BIAS=False, DC_BIAS_V=0.0
                    )
                    file_management.save_image(
                        f"{secondary_label} vs Frequency (LCR)", "Frequency (Hz)", freq_axis,
                        secondary_label, secondary_vals, APPLY_DC_BIAS=False, DC_BIAS_V=0.0
                    )

                    csv_data = np.column_stack([freq_axis, ind_vals, secondary_vals])
                    file_management.save_csv(
                        f"lcr_inductance_vs_freq_{mode.lower()}.csv", csv_data,
                        f"frequency_Hz, inductance_H, {secondary_label.replace(' ', '_')}"
                    )
                    print("Measurement complete. Data saved to output folder.")
            except Exception as e:
                log.error("LCR inductance sweep failed: %s", e)
                print(f"Error during measurement: {e}")

        elif choice == 4:
            # C-V Sweep (Single)
            print("LCR: C-V Sweep (Single Direction)")

            freq = prompt_float("Enter measurement frequency (Hz)", 1000)
            v_min = prompt_float("Enter minimum voltage (V)", -5)
            v_max = prompt_float("Enter maximum voltage (V)", 5)
            num_points = prompt_int("Enter number of points", 201)
            ac_level = prompt_float("Enter AC level (V)", 0.1)

            try:
                with InstrumentSession(
                    lambda: lcr.setup(measurement_function="CPD"),
                    lcr.disconnect_e4980a
                ) as lcr_inst:
                    lcr.set_ac_level(lcr_inst, ac_level)

                    print("\nRunning C-V sweep...")
                    voltage, capacitance, dissipation = lcr.measure_cv_sweep(
                        lcr_inst, freq, v_min, v_max, num_points
                    )

                    plt.figure(figsize=(10, 6))
                    plt.plot(voltage, capacitance, '-o', markersize=3)
                    plt.xlabel('Voltage (V)')
                    plt.ylabel('Capacitance (F)')
                    plt.title(f'C-V Sweep at {freq:.0f} Hz (LCR)')
                    plt.grid(True)
                    file_management.save_plot("lcr_cv_sweep.png")
                    plt.close()

                    csv_data = np.column_stack([voltage, capacitance, dissipation])
                    file_management.save_csv(
                        "lcr_cv_sweep.csv", csv_data,
                        "voltage_V, capacitance_F, dissipation_D"
                    )
                    print("Measurement complete. Data saved to output folder.")
            except Exception as e:
                log.error("LCR C-V sweep failed: %s", e)
                print(f"Error during measurement: {e}")

        elif choice == 5:
            # C-V Butterfly Cycles
            print("LCR: C-V Butterfly Cycles")

            freq = prompt_float("Enter measurement frequency (Hz)", 1000)
            v_min = prompt_float("Enter minimum voltage (V)", -5)
            v_max = prompt_float("Enter maximum voltage (V)", 5)
            num_points = prompt_int("Enter number of points per direction", 201)
            num_cycles = prompt_int("Enter number of cycles", 1)
            ac_level = prompt_float("Enter AC level (V)", 0.1)
            calc_eps = prompt_bool("Calculate permittivity?", False)
            thickness_nm = 0.0
            diameter_um = 0.0
            if calc_eps:
                thickness_nm = prompt_float("Enter dielectric thickness (nm)", 10)
                diameter_um = prompt_float("Enter electrode diameter (µm)", 75)

            try:
                with InstrumentSession(
                    lambda: lcr.setup(measurement_function="CPD"),
                    lcr.disconnect_e4980a
                ) as lcr_inst:
                    lcr.set_ac_level(lcr_inst, ac_level)

                    all_v_cycles = []
                    all_cp_cycles = []
                    all_d_cycles = []
                    all_eps_cycles = []

                    for cycle in range(num_cycles):
                        print(f"\nRunning C-V butterfly cycle {cycle + 1}/{num_cycles}...")
                        voltage, capacitance, dissipation = lcr.measure_cv_butterfly(
                            lcr_inst, freq, v_min, v_max, num_points
                        )
                        all_v_cycles.append(voltage)
                        all_cp_cycles.append(capacitance)
                        all_d_cycles.append(dissipation)
                        if calc_eps:
                            all_eps_cycles.append(lcr.compute_eps_r(capacitance, thickness_nm, diameter_um))

                    file_management.save_cycle_plot(
                        f"C-V Butterfly at {freq:.0f} Hz (LCR)", "Voltage (V)", "Capacitance (F)",
                        all_v_cycles, all_cp_cycles, "lcr_cv_butterfly.png"
                    )

                    if calc_eps:
                        file_management.save_cycle_plot(
                            f"Permittivity Butterfly at {freq:.0f} Hz (LCR)", "Voltage (V)",
                            "Relative Permittivity εr",
                            all_v_cycles, all_eps_cycles, "lcr_eps_butterfly.png"
                        )

                    v_all = np.concatenate(all_v_cycles)
                    cp_all = np.concatenate(all_cp_cycles)
                    d_all = np.concatenate(all_d_cycles)
                    cycle_idx = np.concatenate([np.full(len(v), i+1) for i, v in enumerate(all_v_cycles)])

                    if calc_eps:
                        eps_all = np.concatenate(all_eps_cycles)
                        csv_data = np.column_stack([cycle_idx, v_all, cp_all, d_all, eps_all])
                        header = "cycle, voltage_V, capacitance_F, dissipation_D, eps_r"
                    else:
                        csv_data = np.column_stack([cycle_idx, v_all, cp_all, d_all])
                        header = "cycle, voltage_V, capacitance_F, dissipation_D"

                    file_management.save_csv("lcr_cv_butterfly.csv", csv_data, header)
                    print("Measurement complete. Data saved to output folder.")
            except Exception as e:
                log.error("LCR C-V butterfly failed: %s", e)
                print(f"Error during measurement: {e}")

        elif choice == 6:
            # Single Point Impedance
            print("LCR: Single Point Impedance Measurement")

            freq = prompt_float("Enter frequency (Hz)", 1000)
            ac_level = prompt_float("Enter AC level (V)", 1.0)
            apply_bias = prompt_bool("Apply DC bias?", False)
            bias_v = 0.0
            if apply_bias:
                bias_v = prompt_float("Enter DC bias (V)", 0)

            try:
                with InstrumentSession(
                    lambda: lcr.setup(measurement_function="ZTD"),
                    lcr.disconnect_e4980a
                ) as lcr_inst:
                    lcr.set_ac_level(lcr_inst, ac_level)
                    lcr.configure_dc_bias(lcr_inst, apply_bias, bias_v)

                    print("\nMeasuring impedance...")
                    z_mag, theta_deg = lcr.measure_impedance(lcr_inst, freq)
                    z_real, z_imag = lcr.compute_impedance_components(z_mag, theta_deg)

                    print(f"\n--- Results ---")
                    print(f"Frequency: {freq:.2e} Hz")
                    print(f"|Z|: {z_mag:.6e} Ω")
                    print(f"θ: {theta_deg:.2f}°")
                    print(f"Z_real: {z_real:.6e} Ω")
                    print(f"Z_imag: {z_imag:.6e} Ω")
                    if apply_bias:
                        print(f"DC Bias: {bias_v:.2f} V")
            except Exception as e:
                log.error("LCR single impedance failed: %s", e)
                print(f"Error during measurement: {e}")

        elif choice == 7:
            # Single Point Capacitance
            print("LCR: Single Point Capacitance Measurement")

            freq = prompt_float("Enter frequency (Hz)", 1000)
            mode = prompt_choice("Measurement mode", ["CPD", "CSRS", "CPRP"], "CPD")
            ac_level = prompt_float("Enter AC level (V)", 1.0)
            apply_bias = prompt_bool("Apply DC bias?", False)
            bias_v = 0.0
            if apply_bias:
                bias_v = prompt_float("Enter DC bias (V)", 0)
            calc_eps = prompt_bool("Calculate permittivity?", False)
            thickness_nm = 0.0
            diameter_um = 0.0
            if calc_eps:
                thickness_nm = prompt_float("Enter dielectric thickness (nm)", 10)
                diameter_um = prompt_float("Enter electrode diameter (µm)", 75)

            try:
                with InstrumentSession(
                    lambda: lcr.setup(measurement_function=mode),
                    lcr.disconnect_e4980a
                ) as lcr_inst:
                    lcr.set_ac_level(lcr_inst, ac_level)
                    lcr.configure_dc_bias(lcr_inst, apply_bias, bias_v)

                    print("\nMeasuring capacitance...")
                    cap, secondary = lcr.measure_capacitance(lcr_inst, freq, mode)
                    secondary_label = "D" if mode == "CPD" else "Secondary"

                    print(f"\n--- Results ---")
                    print(f"Frequency: {freq:.2e} Hz")
                    print(f"Mode: {mode}")
                    print(f"Capacitance: {cap:.6e} F")
                    print(f"{secondary_label}: {secondary:.6e}")
                    if apply_bias:
                        print(f"DC Bias: {bias_v:.2f} V")
                    if calc_eps:
                        eps_r = lcr.compute_eps_r(cap, thickness_nm, diameter_um)
                        print(f"Relative Permittivity εr: {eps_r:.2f}")
            except Exception as e:
                log.error("LCR single capacitance failed: %s", e)
                print(f"Error during measurement: {e}")

        elif choice == 8:
            # Quality Factor (Q) vs Frequency
            print("LCR: Quality Factor (Q) vs Frequency")

            freq_start = prompt_float("Enter start frequency (Hz)", 20)
            freq_stop = prompt_float("Enter stop frequency (Hz)", 2e6)
            num_points = prompt_int("Enter number of points", 201)
            sweep_type = prompt_choice("Sweep type", ["LIN", "LOG"], "LOG")
            mode = prompt_choice("Measurement mode", ["CPQ", "LPQ"], "CPQ")
            ac_level = prompt_float("Enter AC level (V)", 1.0)

            try:
                with InstrumentSession(
                    lambda: lcr.setup(measurement_function=mode),
                    lcr.disconnect_e4980a
                ) as lcr_inst:
                    lcr.set_ac_level(lcr_inst, ac_level)

                    print(f"\nRunning Q vs frequency sweep ({mode})...")
                    freq_axis, primary_vals, q_vals = lcr.measure_quality_factor_vs_frequency(
                        lcr_inst, freq_start, freq_stop, num_points, sweep_type, mode
                    )

                    primary_label = "Cp (F)" if mode == "CPQ" else "Lp (H)"
                    file_management.save_image(
                        f"Quality Factor vs Frequency ({mode}) (LCR)", "Frequency (Hz)", freq_axis,
                        "Q (quality factor)", q_vals, APPLY_DC_BIAS=False, DC_BIAS_V=0.0
                    )

                    csv_data = np.column_stack([freq_axis, primary_vals, q_vals])
                    file_management.save_csv(
                        f"lcr_q_vs_freq_{mode.lower()}.csv", csv_data,
                        f"frequency_Hz, {primary_label.replace(' ', '_')}, Q"
                    )
                    print("Measurement complete. Data saved to output folder.")
            except Exception as e:
                log.error("LCR Q sweep failed: %s", e)
                print(f"Error during measurement: {e}")

        elif choice == 9:
            # R-X vs Frequency
            print("LCR: R-X (Resistance-Reactance) vs Frequency")

            freq_start = prompt_float("Enter start frequency (Hz)", 20)
            freq_stop = prompt_float("Enter stop frequency (Hz)", 2e6)
            num_points = prompt_int("Enter number of points", 201)
            sweep_type = prompt_choice("Sweep type", ["LIN", "LOG"], "LOG")
            ac_level = prompt_float("Enter AC level (V)", 1.0)

            try:
                with InstrumentSession(
                    lambda: lcr.setup(measurement_function="RX"),
                    lcr.disconnect_e4980a
                ) as lcr_inst:
                    lcr.set_ac_level(lcr_inst, ac_level)

                    print("\nRunning R-X vs frequency sweep...")
                    freq_axis, r_vals, x_vals = lcr.measure_rx_vs_frequency(
                        lcr_inst, freq_start, freq_stop, num_points, sweep_type
                    )

                    file_management.save_image(
                        "Resistance vs Frequency (LCR)", "Frequency (Hz)", freq_axis,
                        "R (Ohm)", r_vals, APPLY_DC_BIAS=False, DC_BIAS_V=0.0
                    )
                    file_management.save_image(
                        "Reactance vs Frequency (LCR)", "Frequency (Hz)", freq_axis,
                        "X (Ohm)", x_vals, APPLY_DC_BIAS=False, DC_BIAS_V=0.0
                    )

                    csv_data = np.column_stack([freq_axis, r_vals, x_vals])
                    file_management.save_csv("lcr_rx_vs_freq.csv", csv_data, "frequency_Hz, R_ohm, X_ohm")
                    print("Measurement complete. Data saved to output folder.")
            except Exception as e:
                log.error("LCR R-X sweep failed: %s", e)
                print(f"Error during measurement: {e}")

        elif choice == 10:
            # G-B vs Frequency
            print("LCR: G-B (Conductance-Susceptance) vs Frequency")

            freq_start = prompt_float("Enter start frequency (Hz)", 20)
            freq_stop = prompt_float("Enter stop frequency (Hz)", 2e6)
            num_points = prompt_int("Enter number of points", 201)
            sweep_type = prompt_choice("Sweep type", ["LIN", "LOG"], "LOG")
            ac_level = prompt_float("Enter AC level (V)", 1.0)

            try:
                with InstrumentSession(
                    lambda: lcr.setup(measurement_function="GB"),
                    lcr.disconnect_e4980a
                ) as lcr_inst:
                    lcr.set_ac_level(lcr_inst, ac_level)

                    print("\nRunning G-B vs frequency sweep...")
                    freq_axis, g_vals, b_vals = lcr.measure_gb_vs_frequency(
                        lcr_inst, freq_start, freq_stop, num_points, sweep_type
                    )

                    file_management.save_image(
                        "Conductance vs Frequency (LCR)", "Frequency (Hz)", freq_axis,
                        "G (S)", g_vals, APPLY_DC_BIAS=False, DC_BIAS_V=0.0
                    )
                    file_management.save_image(
                        "Susceptance vs Frequency (LCR)", "Frequency (Hz)", freq_axis,
                        "B (S)", b_vals, APPLY_DC_BIAS=False, DC_BIAS_V=0.0
                    )

                    csv_data = np.column_stack([freq_axis, g_vals, b_vals])
                    file_management.save_csv("lcr_gb_vs_freq.csv", csv_data, "frequency_Hz, G_S, B_S")
                    print("Measurement complete. Data saved to output folder.")
            except Exception as e:
                log.error("LCR G-B sweep failed: %s", e)
                print(f"Error during measurement: {e}")

        elif choice == 11:
            # Open/Short Correction
            print("LCR: Open/Short Correction")
            print("  1. Open correction")
            print("  2. Short correction")
            print("  3. Open + Short (recommended)")
            print("  4. Disable all corrections")

            corr_choice = prompt_int("Select correction type", 3)

            try:
                with InstrumentSession(
                    lambda: lcr.setup(measurement_function="CPD"),
                    lcr.disconnect_e4980a
                ) as lcr_inst:
                    if corr_choice in (1, 3):
                        input("Remove DUT / open probes, then press Enter...")
                        lcr.perform_open_correction(lcr_inst)

                    if corr_choice in (2, 3):
                        input("Short probes together, then press Enter...")
                        lcr.perform_short_correction(lcr_inst)

                    if corr_choice == 4:
                        lcr.disable_corrections(lcr_inst)

                    print("Correction procedure complete.")
            except Exception as e:
                log.error("LCR correction failed: %s", e)
                print(f"Error during correction: {e}")

        repeat_input = input("\nDo you want to perform another measurement? (y/n): ").strip().lower()
        if repeat_input != 'y':
            repeating = False


if __name__ == "__main__":
    main()
