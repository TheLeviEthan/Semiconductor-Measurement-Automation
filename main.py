"""
Filename: main.py
Author: Ethan Ruddell
Date: 2025-01-23
Description: Implements CLI for automating semiconductor measurements using
the Keysight 4294A Precision Impedance Analyzer and a DC bias source.
"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import pia
import pspa
import file_management

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
    "Permittivity vs Frequency"
]

pspa_measurements = [
    "Transistor Output Characteristics",
    "Transistor Transfer Characteristics",
    "I-V Curve (Unidirectional)",
    "I-V Curve (Bidirectional)",
    "Pulsed I-V (Single Device)",
    "Pulsed Transistor"
]


def safe_float_input(prompt, default):
    """Safely get float input from user with default value and error handling."""
    while True:
        try:
            user_input = input(prompt).strip()
            if not user_input:
                return float(default)
            return float(user_input)
        except ValueError:
            print(f"Invalid input. Please enter a number.")


def safe_int_input(prompt, default):
    """Safely get integer input from user with default value and error handling."""
    while True:
        try:
            user_input = input(prompt).strip()
            if not user_input:
                return int(default)
            return int(user_input)
        except ValueError:
            print(f"Invalid input. Please enter an integer.")


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
    # TODO: add save path in CLI, first prompt should select tool

    args = parse_args()
    
    print("Welcome to the NRG Semiconductor Measurement Automation...\n")
    
    # Prompt for output directory once at the start
    default_dir = args.output_dir or file_management.default_output_dir
    output_path = input(f"Enter output directory [default: {default_dir}]: ").strip() or default_dir
    file_management.set_output_dir(output_path)
    file_management.ensure_output_dir(file_management.output_dir)
    print(f"Output directory set to: {file_management.output_dir}\n")

    while True:
        # Prompt user to select instrument or quit
        print("\n" + "="*60)
        print("SELECT INSTRUMENT")
        print("="*60)
        print("1. PIA (Precision Impedance Analyzer)")
        print("2. PSPA (Parameter/Source Analyzer)")
        print("q. Quit")
        
        instrument_choice = input("\nEnter your choice (1, 2, or q): ").strip().lower()
        
        if instrument_choice == 'q':
            print("Exiting the program. Goodbye!")
            break
        elif instrument_choice == '1':
            # PIA measurements menu
            run_pia_measurements()
        elif instrument_choice == '2':
            # PSPA measurements menu
            run_pspa_measurements()
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


def execute_pia_measurement(choice):
    """Execute the selected PIA measurement."""
    repeating = True
    while repeating:
        if choice == 1:
            # Impedance Magnitude vs Frequency
            print("You selected Impedance Magnitude vs Frequency measurement.")
            
            # Prompt for parameter duplication
            if pia.current_parameters["freq_start_hz"] != pia.FREQ_START_HZ or pia.current_parameters["apply_dc_bias"] != pia.APPLY_DC_BIAS:
                use_last_params = pia.prompt_for_parameter_duplication()
            else:
                use_last_params = False
            
            # Get frequency and DC bias parameters
            freq_start, freq_stop, num_points = pia.get_frequency_parameters(use_last=use_last_params)
            apply_dc_bias, dc_bias_v = pia.get_dc_bias_parameters(use_last=use_last_params)
            
            # Call the measurement
            inst = pia.setup()
            pia.initialize_4294a_for_impedance(inst)
            pia.configure_dc_bias(inst, apply_dc_bias, dc_bias_v)
            freq_axis, z_mag, theta_deg = pia.measure_impedance_vs_freq(inst, freq_start, freq_stop, num_points)
            
            # Save magnitude plot only
            file_management.save_image(
                "Impedance Magnitude vs Frequency", "Frequency (Hz)", freq_axis, "|Z| (Ohm)", z_mag,
                APPLY_DC_BIAS=apply_dc_bias, DC_BIAS_V=dc_bias_v
            )
            
            # Save CSV with all impedance data
            z_real = z_mag * np.cos(np.deg2rad(theta_deg))
            z_imag = z_mag * np.sin(np.deg2rad(theta_deg))
            csv_data = np.column_stack([freq_axis, z_mag, theta_deg, z_real, z_imag])
            file_management.save_csv(
                "impedance_data.csv",
                csv_data,
                "frequency_Hz, Z_mag_ohm, theta_deg, Z_real_ohm, Z_imag_ohm"
            )
            inst.close()
            
        elif choice == 2:
            # Impedance Phase vs Frequency
            print("You selected Impedance Phase vs Frequency measurement.")
            
            # Prompt for parameter duplication
            if pia.current_parameters["freq_start_hz"] != pia.FREQ_START_HZ or pia.current_parameters["apply_dc_bias"] != pia.APPLY_DC_BIAS:
                use_last_params = pia.prompt_for_parameter_duplication()
            else:
                use_last_params = False
            
            # Get frequency and DC bias parameters
            freq_start, freq_stop, num_points = pia.get_frequency_parameters(use_last=use_last_params)
            apply_dc_bias, dc_bias_v = pia.get_dc_bias_parameters(use_last=use_last_params)
            
            # Call the measurement
            inst = pia.setup()
            pia.initialize_4294a_for_impedance(inst)
            pia.configure_dc_bias(inst, apply_dc_bias, dc_bias_v)
            freq_axis, z_mag, theta_deg = pia.measure_impedance_vs_freq(inst, freq_start, freq_stop, num_points)
            
            # Save phase plot only
            file_management.save_image(
                "Impedance Phase vs Frequency", "Frequency (Hz)", freq_axis, "Phase (Degrees)", theta_deg,
                APPLY_DC_BIAS=apply_dc_bias, DC_BIAS_V=dc_bias_v
            )
            
            # Save CSV with all impedance data
            z_real = z_mag * np.cos(np.deg2rad(theta_deg))
            z_imag = z_mag * np.sin(np.deg2rad(theta_deg))
            csv_data = np.column_stack([freq_axis, z_mag, theta_deg, z_real, z_imag])
            file_management.save_csv(
                "impedance_data.csv",
                csv_data,
                "frequency_Hz, Z_mag_ohm, theta_deg, Z_real_ohm, Z_imag_ohm"
            )
            inst.close()
            
        elif choice == 3:
            # Capacitance vs Frequency
            print("You selected Capacitance vs Frequency measurement.")
            
            # Prompt for parameter duplication
            if pia.current_parameters["freq_start_hz"] != pia.FREQ_START_HZ or pia.current_parameters["apply_dc_bias"] != pia.APPLY_DC_BIAS:
                use_last_params = pia.prompt_for_parameter_duplication()
            else:
                use_last_params = False
            
            # Get frequency and DC bias parameters
            freq_start, freq_stop, num_points = pia.get_frequency_parameters(use_last=use_last_params)
            apply_dc_bias, dc_bias_v = pia.get_dc_bias_parameters(use_last=use_last_params)
            
            # Call the measurement
            inst = pia.setup()
            pia.initialize_4294a_for_cpd(inst)
            pia.configure_dc_bias(inst, apply_dc_bias, dc_bias_v)
            freq_axis, cp_vals, d_vals = pia.measure_cpd_vs_freq(inst, freq_start, freq_stop, num_points)
            
            # Save Cp plot only
            file_management.save_image(
                "Capacitance vs Frequency", "Frequency (Hz)", freq_axis, "Cp (F)", cp_vals,
                APPLY_DC_BIAS=apply_dc_bias, DC_BIAS_V=dc_bias_v
            )
            
            # Save Cp-D data to CSV
            csv_data = np.column_stack([freq_axis, cp_vals, d_vals])
            file_management.save_csv(
                "cpd_data.csv",
                csv_data,
                "frequency_Hz, Cp_F, D"
            )
            inst.close()
            
        elif choice == 4:
            # Tan Loss vs Frequency
            print("You selected Tan Loss vs Frequency measurement.")
            
            # Prompt for parameter duplication
            if pia.current_parameters["freq_start_hz"] != pia.FREQ_START_HZ or pia.current_parameters["apply_dc_bias"] != pia.APPLY_DC_BIAS:
                use_last_params = pia.prompt_for_parameter_duplication()
            else:
                use_last_params = False
            
            # Get frequency and DC bias parameters
            freq_start, freq_stop, num_points = pia.get_frequency_parameters(use_last=use_last_params)
            apply_dc_bias, dc_bias_v = pia.get_dc_bias_parameters(use_last=use_last_params)
            
            # Call the measurement
            inst = pia.setup()
            pia.initialize_4294a_for_cpd(inst)
            pia.configure_dc_bias(inst, apply_dc_bias, dc_bias_v)
            freq_axis, cp_vals, d_vals = pia.measure_cpd_vs_freq(inst, freq_start, freq_stop, num_points)
            
            # Save D plot only
            file_management.save_image(
                "Tan Loss vs Frequency", "Frequency (Hz)", freq_axis, "D (loss factor)", d_vals,
                APPLY_DC_BIAS=apply_dc_bias, DC_BIAS_V=dc_bias_v
            )
            
            # Save Cp-D data to CSV
            csv_data = np.column_stack([freq_axis, cp_vals, d_vals])
            file_management.save_csv(
                "cpd_data.csv",
                csv_data,
                "frequency_Hz, Cp_F, D"
            )
            inst.close()
            
        elif choice == 5:
            # C–V Butterfly Cycle: Cp vs Voltage
            print("C–V Butterfly Cycle: Cp vs Voltage")
            freq_cv = float(input("Enter measurement frequency (Hz) [default 25000]: ") or "25000")
            v_min = float(input("Enter minimum voltage (V) [default 0]: ") or "0")
            v_max = float(input("Enter maximum voltage (V) [default 5]: ") or "5")
            n_points = int(input("Enter number of points per sweep [default 401]: ") or "401")
            n_cycles = int(input("Enter number of cycles [default 1]: ") or "1")
            thickness_nm = float(input("Enter HZO thickness (nm) [default 10.0]: ") or "10.0")
            diam_um = float(input("Enter electrode diameter (µm) [default 75.0]: ") or "75.0")
            
            inst = pia.setup()
            try:
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
                
                # Flatten for saving
                v_all = np.concatenate(all_v_cycles)
                cp_all = np.concatenate(all_cp_cycles)
                eps_all = np.concatenate(all_eps_cycles)
                cycles_all = np.concatenate(all_cycle_idx)
                
                # Save C–V & εr(V) data
                csv_data = np.column_stack([cycles_all, v_all, cp_all, eps_all])
                file_management.save_csv(
                    "cv_dielectric_cycles.csv",
                    csv_data,
                    "cycle_index, bias_V, Cp_F, eps_r"
                )

                # Save Cp butterfly plot only
                file_management.save_cycle_plot(
                    "C-V Butterfly: Cp vs Voltage",
                    "Bias Voltage (V)",
                    "Cp (F)",
                    all_v_cycles,
                    all_cp_cycles,
                    "cv_butterfly_cp.png",
                )
                print("C–V Cp measurement complete. Data saved to output folder.")
            finally:
                try:
                    inst.write("DCO OFF")
                except Exception:
                    pass
                inst.close()
                
        elif choice == 6:
            # C–V Butterfly Cycle: Permittivity vs Voltage
            print("C–V Butterfly Cycle: Permittivity vs Voltage")
            freq_cv = float(input("Enter measurement frequency (Hz) [default 25000]: ") or "25000")
            v_min = float(input("Enter minimum voltage (V) [default 0]: ") or "0")
            v_max = float(input("Enter maximum voltage (V) [default 5]: ") or "5")
            n_points = int(input("Enter number of points per sweep [default 401]: ") or "401")
            n_cycles = int(input("Enter number of cycles [default 1]: ") or "1")
            thickness_nm = float(input("Enter HZO thickness (nm) [default 10.0]: ") or "10.0")
            diam_um = float(input("Enter electrode diameter (µm) [default 75.0]: ") or "75.0")
            
            inst = pia.setup()
            try:
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
                
                # Flatten for saving
                v_all = np.concatenate(all_v_cycles)
                cp_all = np.concatenate(all_cp_cycles)
                eps_all = np.concatenate(all_eps_cycles)
                cycles_all = np.concatenate(all_cycle_idx)
                
                # Save C–V & εr(V) data
                csv_data = np.column_stack([cycles_all, v_all, cp_all, eps_all])
                file_management.save_csv(
                    "cv_dielectric_cycles.csv",
                    csv_data,
                    "cycle_index, bias_V, Cp_F, eps_r"
                )

                # Save eps_r butterfly plot only
                file_management.save_cycle_plot(
                    "C-V Butterfly: eps_r vs Voltage",
                    "Bias Voltage (V)",
                    "Dielectric constant eps_r",
                    all_v_cycles,
                    all_eps_cycles,
                    "cv_butterfly_eps_r.png",
                )
                print("C–V permittivity measurement complete. Data saved to output folder.")
            finally:
                try:
                    inst.write("DCO OFF")
                except Exception:
                    pass
                inst.close()
                
        elif choice == 7:
            # εr vs frequency measurement
            print("Permittivity vs Frequency Measurement")
            
            # Prompt for parameter duplication
            if pia.current_parameters["freq_start_hz"] != pia.FREQ_START_HZ or pia.current_parameters["apply_dc_bias"] != pia.APPLY_DC_BIAS:
                use_last_params = pia.prompt_for_parameter_duplication()
            else:
                use_last_params = False
            
            if use_last_params:
                freq_start = pia.current_parameters["freq_start_hz"]
                freq_stop = pia.current_parameters["freq_stop_hz"]
                freq_points = pia.current_parameters["num_points"]
                bias_voltage = pia.current_parameters["dc_bias_v"]
            else:
                freq_start = float(input("Enter start frequency (Hz) [default 1000]: ") or "1000")
                freq_stop = float(input("Enter stop frequency (Hz) [default 1e6]: ") or "1e6")
                freq_points = int(input("Enter number of frequency points [default 201]: ") or "201")
                bias_voltage = float(input("Enter DC bias voltage (V) [default 0]: ") or "0")
                
                # Update current parameters for future duplication
                pia.current_parameters["freq_start_hz"] = freq_start
                pia.current_parameters["freq_stop_hz"] = freq_stop
                pia.current_parameters["num_points"] = freq_points
                pia.current_parameters["dc_bias_v"] = bias_voltage
                pia.current_parameters["apply_dc_bias"] = (bias_voltage != 0)
            
            thickness_nm = float(input("Enter HZO thickness (nm) [default 10.0]: ") or "10.0")
            diam_um = float(input("Enter electrode diameter (µm) [default 75.0]: ") or "75.0")
            
            inst = pia.setup()
            try:
                pia.initialize_4294a_for_cpd(inst)
                
                print("Running εr vs frequency sweep...")
                freq_axis, cp_f = pia.measure_eps_vs_freq(inst, freq_start, freq_stop, freq_points, bias_voltage)
                eps_f = pia.compute_eps_r(cp_f, thickness_nm, diam_um)
                
                # Save εr(f) data
                csv_data = np.column_stack([freq_axis, cp_f, eps_f])
                file_management.save_csv(
                    "eps_vs_freq.csv",
                    csv_data,
                    "frequency_Hz, Cp_F, eps_r"
                )
                
                # Save plot
                file_management.save_image(
                    "Permittivity vs Frequency", "Frequency (Hz)", freq_axis, "Dielectric constant εr", eps_f,
                    APPLY_DC_BIAS=(bias_voltage != 0), DC_BIAS_V=bias_voltage
                )
                
                print("Frequency sweep complete. Data and plot saved to output folder.")
            finally:
                try:
                    inst.write("DCO OFF")
                except Exception:
                    pass
                inst.close()
        
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
            
            vds_start = safe_float_input("Enter start Vds (V) [default 0]: ", 0)
            vds_stop = safe_float_input("Enter stop Vds (V) [default 5]: ", 5)
            vds_step = safe_float_input("Enter Vds step (V) [default 0.1]: ", 0.1)
            vgs_start = safe_float_input("Enter start Vgs (V) [default 0]: ", 0)
            vgs_stop = safe_float_input("Enter stop Vgs (V) [default 3]: ", 3)
            vgs_step = safe_float_input("Enter Vgs step (V) [default 0.5]: ", 0.5)
            compliance = safe_float_input("Enter current compliance (A) [default 0.1]: ", 0.1)
            
            drain_ch = safe_int_input("Enter drain channel [default 1]: ", 1)
            gate_ch = safe_int_input("Enter gate channel [default 2]: ", 2)
            source_ch = safe_int_input("Enter source channel [default 3]: ", 3)
            
            try:
                pspa_inst = pspa.connect_pspa()
                
                print("\nRunning transistor output characteristics sweep...")
                data = pspa.measure_transistor_output_characteristics(
                    pspa_inst, vds_start, vds_stop, vds_step,
                    vgs_start, vgs_stop, vgs_step,
                    drain_ch, gate_ch, source_ch, compliance
                )
                
                # Save data to CSV
                csv_data = np.column_stack([data['Vds'], data['Vgs'], data['Id']])
                file_management.save_csv(
                    "transistor_output_chars.csv",
                    csv_data,
                    "Vds_V, Vgs_V, Id_A"
                )
                
                # Plot Id vs Vds for each Vgs
                plt.figure(figsize=(10, 6))
                vgs_values = np.unique(data['Vgs'])
                for vgs in vgs_values:
                    mask = data['Vgs'] == vgs
                    plt.plot(data['Vds'][mask], data['Id'][mask] * 1e3, marker='o', label=f"Vgs = {vgs:.2f} V")
                plt.xlabel('Vds (V)')
                plt.ylabel('Id (mA)')
                plt.title('Transistor Output Characteristics')
                plt.legend()
                plt.grid(True)
                plt.tight_layout()
                file_management.ensure_output_dir(os.path.join(file_management.output_dir, "images"))
                plot_path = file_management.uniquify(os.path.join(file_management.output_dir, "images", "transistor_output_chars.png"))
                plt.savefig(plot_path, dpi=300)
                plt.close()
                print(f"Plot saved: {plot_path}")
                
                print("Measurement complete. Data saved to output folder.")
            except Exception as e:
                print(f"Error during measurement: {e}")
            finally:
                try:
                    pspa.disconnect_pspa(pspa_inst)
                except Exception:
                    pass
        
        elif choice == 2:
            # PSPA Transistor Transfer Characteristics
            print("PSPA: Transistor Transfer Characteristics (Id-Vgs curve)")
            
            vgs_start = safe_float_input("Enter start Vgs (V) [default -1]: ", -1)
            vgs_stop = safe_float_input("Enter stop Vgs (V) [default 3]: ", 3)
            vgs_step = safe_float_input("Enter Vgs step (V) [default 0.05]: ", 0.05)
            vds_constant = safe_float_input("Enter constant Vds (V) [default 5]: ", 5)
            compliance = safe_float_input("Enter current compliance (A) [default 0.1]: ", 0.1)
            
            drain_ch = safe_int_input("Enter drain channel [default 1]: ", 1)
            gate_ch = safe_int_input("Enter gate channel [default 2]: ", 2)
            source_ch = safe_int_input("Enter source channel [default 3]: ", 3)
            
            try:
                pspa_inst = pspa.connect_pspa()
                
                print("\nRunning transistor transfer characteristics sweep...")
                data = pspa.measure_transistor_transfer_characteristics(
                    pspa_inst, vgs_start, vgs_stop, vgs_step, vds_constant,
                    drain_ch, gate_ch, source_ch, compliance
                )
                
                # Save data to CSV
                csv_data = np.column_stack([data['Vgs'], data['Id'], data['Ig']])
                file_management.save_csv(
                    "transistor_transfer_chars.csv",
                    csv_data,
                    "Vgs_V, Id_A, Ig_A"
                )
                
                # Plot Id and Ig vs Vgs
                fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
                
                # Linear scale
                ax1.plot(data['Vgs'], data['Id'] * 1e3, marker='o', label='Id')
                ax1.set_xlabel('Vgs (V)')
                ax1.set_ylabel('Id (mA)')
                ax1.set_title(f'Transfer Characteristics (Vds = {vds_constant} V) - Linear')
                ax1.grid(True)
                ax1.legend()
                
                # Log scale
                ax2.semilogy(data['Vgs'], np.abs(data['Id']), marker='o', label='|Id|')
                ax2.semilogy(data['Vgs'], np.abs(data['Ig']), marker='s', label='|Ig|')
                ax2.set_xlabel('Vgs (V)')
                ax2.set_ylabel('Current (A)')
                ax2.set_title(f'Transfer Characteristics (Vds = {vds_constant} V) - Log')
                ax2.grid(True)
                ax2.legend()
                
                plt.tight_layout()
                file_management.ensure_output_dir(os.path.join(file_management.output_dir, "images"))
                plot_path = file_management.uniquify(os.path.join(file_management.output_dir, "images", "transistor_transfer_chars.png"))
                plt.savefig(plot_path, dpi=300)
                plt.close()
                print(f"Plot saved: {plot_path}")
                
                print("Measurement complete. Data saved to output folder.")
            except Exception as e:
                print(f"Error during measurement: {e}")
            finally:
                try:
                    pspa.disconnect_pspa(pspa_inst)
                except Exception:
                    pass
        
        elif choice == 3:
            # PSPA I-V Curve (Unidirectional)
            print("PSPA: I-V Curve (Unidirectional)")
            
            v_start = safe_float_input("Enter start voltage (V) [default 0]: ", 0)
            v_stop = safe_float_input("Enter stop voltage (V) [default 5]: ", 5)
            v_step = safe_float_input("Enter voltage step (V) [default 0.1]: ", 0.1)
            compliance = safe_float_input("Enter current compliance (A) [default 0.1]: ", 0.1)
            channel = safe_int_input("Enter channel number [default 1]: ", 1)
            
            try:
                pspa_inst = pspa.connect_pspa()
                
                print("\nRunning I-V measurement...")
                data = pspa.measure_iv_curve(pspa_inst, v_start, v_stop, v_step, channel, compliance)
                
                # Save data to CSV
                csv_data = np.column_stack([data['Voltage'], data['Current']])
                file_management.save_csv(
                    "iv_curve.csv",
                    csv_data,
                    "Voltage_V, Current_A"
                )
                
                # Plot I-V curve
                plt.figure(figsize=(10, 6))
                plt.plot(data['Voltage'], data['Current'] * 1e3, marker='o')
                plt.xlabel('Voltage (V)')
                plt.ylabel('Current (mA)')
                plt.title('I-V Curve (Unidirectional)')
                plt.grid(True)
                plt.tight_layout()
                file_management.ensure_output_dir(os.path.join(file_management.output_dir, "images"))
                plot_path = file_management.uniquify(os.path.join(file_management.output_dir, "images", "iv_curve.png"))
                plt.savefig(plot_path, dpi=300)
                plt.close()
                print(f"Plot saved: {plot_path}")
                
                print("Measurement complete. Data saved to output folder.")
            except Exception as e:
                print(f"Error during measurement: {e}")
            finally:
                try:
                    pspa.disconnect_pspa(pspa_inst)
                except Exception:
                    pass
        
        elif choice == 4:
            # PSPA I-V Curve (Bidirectional)
            print("PSPA: I-V Curve (Bidirectional)")
            
            v_max = safe_float_input("Enter maximum voltage magnitude (V) [default 5]: ", 5)
            v_step = safe_float_input("Enter voltage step (V) [default 0.1]: ", 0.1)
            compliance = safe_float_input("Enter current compliance (A) [default 0.1]: ", 0.1)
            channel = safe_int_input("Enter channel number [default 1]: ", 1)
            
            try:
                pspa_inst = pspa.connect_pspa()
                
                print("\nRunning bidirectional I-V measurement...")
                data = pspa.measure_iv_bidirectional(pspa_inst, v_max, v_step, channel, compliance)
                
                # Save data to CSV
                csv_data = np.column_stack([data['Voltage'], data['Current']])
                file_management.save_csv(
                    "iv_curve_bidirectional.csv",
                    csv_data,
                    "Voltage_V, Current_A"
                )
                
                # Plot I-V curve
                plt.figure(figsize=(10, 6))
                plt.plot(data['Voltage'], data['Current'] * 1e3, marker='o')
                plt.xlabel('Voltage (V)')
                plt.ylabel('Current (mA)')
                plt.title('I-V Curve (Bidirectional)')
                plt.grid(True)
                plt.tight_layout()
                file_management.ensure_output_dir(os.path.join(file_management.output_dir, "images"))
                plot_path = file_management.uniquify(os.path.join(file_management.output_dir, "images", "iv_curve_bidirectional.png"))
                plt.savefig(plot_path, dpi=300)
                plt.close()
                print(f"Plot saved: {plot_path}")
                
                print("Measurement complete. Data saved to output folder.")
            except Exception as e:
                print(f"Error during measurement: {e}")
            finally:
                try:
                    pspa.disconnect_pspa(pspa_inst)
                except Exception:
                    pass
        
        elif choice == 5:
            # Pulsed I-V (Single Device)
            print("PSPA: Pulsed I-V (Single Device)")
            
            v_base = safe_float_input("Enter base voltage (V) [default 0]: ", 0)
            v_pulse = safe_float_input("Enter pulse voltage (V) [default 5]: ", 5)
            pulse_width = safe_float_input("Enter pulse width (µs) [default 100]: ", 100) * 1e-6
            pulse_period = safe_float_input("Enter pulse period (ms) [default 10]: ", 10) * 1e-3
            num_pulses = safe_int_input("Enter number of pulses [default 10]: ", 10)
            compliance = safe_float_input("Enter current compliance (A) [default 0.1]: ", 0.1)
            channel = safe_int_input("Enter channel number [default 1]: ", 1)
            
            try:
                pspa_inst = pspa.connect_pspa()
                
                print("\nRunning pulsed I-V measurement...")
                data = pspa.measure_pulsed_iv(
                    pspa_inst, v_base, v_pulse, pulse_width, pulse_period, num_pulses,
                    channel, compliance
                )
                
                # Save data to CSV
                csv_data = np.column_stack([data['Time'], data['Voltage'], data['Current']])
                file_management.save_csv(
                    "pulsed_iv.csv",
                    csv_data,
                    "Time_s, Voltage_V, Current_A"
                )
                
                # Plot pulsed measurements
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
                file_management.ensure_output_dir(os.path.join(file_management.output_dir, "images"))
                plot_path = file_management.uniquify(os.path.join(file_management.output_dir, "images", "pulsed_iv.png"))
                plt.savefig(plot_path, dpi=300)
                plt.close()
                print(f"Plot saved: {plot_path}")
                
                print("Measurement complete. Data saved to output folder.")
            except Exception as e:
                print(f"Error during measurement: {e}")
            finally:
                try:
                    pspa.disconnect_pspa(pspa_inst)
                except Exception:
                    pass
        
        elif choice == 6:
            # Pulsed Transistor
            print("PSPA: Pulsed Transistor Measurement")
            print("NOTE: Gate voltage is held constant (DC) while drain is pulsed.")
            
            vds_base = safe_float_input("Enter base Vds (V) [default 0]: ", 0)
            vds_pulse = safe_float_input("Enter pulse Vds (V) [default 5]: ", 5)
            vgs_base = safe_float_input("Enter base Vgs (V) [default 0]: ", 0)
            vgs_pulse = safe_float_input("Enter pulse Vgs (V) [default 3]: ", 3)
            pulse_width = safe_float_input("Enter pulse width (µs) [default 100]: ", 100) * 1e-6
            pulse_period = safe_float_input("Enter pulse period (ms) [default 10]: ", 10) * 1e-3
            num_pulses = safe_int_input("Enter number of pulses [default 10]: ", 10)
            compliance = safe_float_input("Enter current compliance (A) [default 0.1]: ", 0.1)
            
            drain_ch = safe_int_input("Enter drain channel [default 1]: ", 1)
            gate_ch = safe_int_input("Enter gate channel [default 2]: ", 2)
            source_ch = safe_int_input("Enter source channel [default 3]: ", 3)
            
            try:
                pspa_inst = pspa.connect_pspa()
                
                print("\nRunning pulsed transistor measurement...")
                data = pspa.measure_pulsed_transistor(
                    pspa_inst, vds_pulse, vgs_pulse, vds_base, vgs_base,
                    pulse_width, pulse_period, num_pulses,
                    drain_ch, gate_ch, source_ch, compliance
                )
                
                # Save data to CSV
                csv_data = np.column_stack([data['Time'], data['Id']])
                file_management.save_csv(
                    "pulsed_transistor.csv",
                    csv_data,
                    "Time_s, Id_A"
                )
                
                # Plot pulsed measurements
                plt.figure(figsize=(10, 6))
                
                plt.plot(data['Time'] * 1e3, data['Id'] * 1e3, marker='o', label='Id')
                plt.xlabel('Time (ms)')
                plt.ylabel('Id (mA)')
                plt.title(f'Pulsed Transistor: Drain Current (Vds={vds_pulse}V, Vgs={vgs_pulse}V)')
                plt.grid(True)
                plt.legend()
                
                plt.tight_layout()
                file_management.ensure_output_dir(os.path.join(file_management.output_dir, "images"))
                plot_path = file_management.uniquify(os.path.join(file_management.output_dir, "images", "pulsed_transistor.png"))
                plt.savefig(plot_path, dpi=300)
                plt.close()
                print(f"Plot saved: {plot_path}")
                
                print("Measurement complete. Data saved to output folder.")
            except Exception as e:
                print(f"Error during measurement: {e}")
            finally:
                try:
                    pspa.disconnect_pspa(pspa_inst)
                except Exception:
                    pass
        
        repeat_input = input("\nDo you want to perform another measurement? (y/n): ").strip().lower()
        if repeat_input != 'y':
            repeating = False


if __name__ == "__main__":
    main()
