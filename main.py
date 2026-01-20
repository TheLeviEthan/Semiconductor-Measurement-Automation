import os
import sys
import pyvisa
import numpy as np
import pia
import file_management

# =============================
# User settings and constants
# =============================

measurementType = ["Impedance vs Theta", "Capacitance vs Tan Loss", "Dielectric Measurements"]

def main():
    print("Welcome to the NRG Semiconductor Measurement Automation...\n")

    while True:
        # Prompt user for measurement type or to quit
        user_input = input("Enter 'q' to quit or press Enter to continue: ").strip().lower()
        if user_input == 'q':
            print("Exiting the program. Goodbye!")
            break
        print("Please type a number to select the measurement you will be performing:\n")
        for i, meas in enumerate(measurementType):
            print(f"{i + 1}. {meas}")
        choice = int(input("\nEnter your choice (1-3): "))

        if choice == 1:
            repeating = True
            while repeating:
                print("You selected Impedance/Theta measurement.")
                # Call the relevant function or module for Impedance vs Theta
                inst = pia.setup()
                pia.initialize_4294a_for_impedance(inst)
                pia.configure_dc_bias(inst, pia.APPLY_DC_BIAS, pia.DC_BIAS_V)
                freq_axis, z_mag, theta_deg = pia.measure_impedance_vs_freq(inst)
                file_management.ensure_output_dir(file_management.output_dir)
                
                # Save impedance data to CSV
                z_real = z_mag * np.cos(np.deg2rad(theta_deg))
                z_imag = z_mag * np.sin(np.deg2rad(theta_deg))
                csv_data = np.column_stack([freq_axis, z_mag, theta_deg, z_real, z_imag])
                file_management.save_csv(
                    os.path.join(file_management.output_dir, "impedance_data.csv"),
                    csv_data,
                    "frequency_Hz, Z_mag_ohm, theta_deg, Z_real_ohm, Z_imag_ohm"
                )
                
                # Save magnitude and phase plots
                file_management.save_image(
                    "Impedance Magnitude vs Frequency", "Frequency (Hz)", freq_axis, "|Z| (Ohm)", z_mag,
                    APPLY_DC_BIAS=pia.APPLY_DC_BIAS, DC_BIAS_V=pia.DC_BIAS_V
                )
                file_management.save_image(
                    "Impedance Phase vs Frequency", "Frequency (Hz)", freq_axis, "Phase (Degrees)", theta_deg,
                    APPLY_DC_BIAS=pia.APPLY_DC_BIAS, DC_BIAS_V=pia.DC_BIAS_V
                )
                inst.close()
                repeat_input = input("Do you want to perform another Impedance/Theta measurement? (y/n): ").strip().lower()
                if repeat_input != 'y':
                    repeating = False
        elif choice == 2:
            repeating = True
            while repeating:
                print("You selected Capacitance/Tan Loss measurement.")
                # Cp-D frequency sweep
                inst = pia.setup()
                pia.initialize_4294a_for_cpd(inst)
                pia.configure_dc_bias(inst, pia.APPLY_DC_BIAS, pia.DC_BIAS_V)
                freq_axis, cp_vals, d_vals = pia.measure_cpd_vs_freq(inst)
                file_management.ensure_output_dir(file_management.output_dir)
                
                # Save Cp-D data to CSV
                csv_data = np.column_stack([freq_axis, cp_vals, d_vals])
                file_management.save_csv(
                    os.path.join(file_management.output_dir, "cpd_data.csv"),
                    csv_data,
                    "frequency_Hz, Cp_F, D"
                )
                
                # Save Cp and D plots
                file_management.save_image(
                    "Capacitance vs Frequency", "Frequency (Hz)", freq_axis, "Cp (F)", cp_vals,
                    APPLY_DC_BIAS=pia.APPLY_DC_BIAS, DC_BIAS_V=pia.DC_BIAS_V
                )
                file_management.save_image(
                    "Tan Loss vs Frequency", "Frequency (Hz)", freq_axis, "D (loss factor)", d_vals,
                    APPLY_DC_BIAS=pia.APPLY_DC_BIAS, DC_BIAS_V=pia.DC_BIAS_V
                )
                inst.close()
                repeat_input = input("Do you want to perform another Capacitance/Tan Loss measurement? (y/n): ").strip().lower()
                if repeat_input != 'y':
                    repeating = False
        elif choice == 3:
            print("You selected Dielectric Measurements.")
            dielectric_submenu = True
            while dielectric_submenu:
                print("\nSelect dielectric measurement type:")
                print("1. C–V Butterfly Cycles")
                print("2. Permittivity vs Frequency")
                print("3. Back to main menu")
                sub_choice = int(input("Enter your choice (1-3): "))
                
                if sub_choice == 1:
                    # C–V butterfly cycle measurement
                    print("C–V Butterfly Cycle Measurement")
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
                        file_management.ensure_output_dir(file_management.output_dir)
                        
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
                            os.path.join(file_management.output_dir, "cv_dielectric_cycles.csv"),
                            csv_data,
                            "cycle_index, bias_V, Cp_F, eps_r"
                        )
                        
                        print("C–V measurement complete. Data saved to output folder.")
                    finally:
                        try:
                            inst.write("DCO OFF")
                        except Exception:
                            pass
                        inst.close()
                
                elif sub_choice == 2:
                    # εr vs frequency measurement
                    print("Permittivity vs Frequency Measurement")
                    freq_start = float(input("Enter start frequency (Hz) [default 1000]: ") or "1000")
                    freq_stop = float(input("Enter stop frequency (Hz) [default 1e6]: ") or "1e6")
                    freq_points = int(input("Enter number of frequency points [default 201]: ") or "201")
                    bias_voltage = float(input("Enter DC bias voltage (V) [default 0]: ") or "0")
                    thickness_nm = float(input("Enter HZO thickness (nm) [default 10.0]: ") or "10.0")
                    diam_um = float(input("Enter electrode diameter (µm) [default 75.0]: ") or "75.0")
                    
                    inst = pia.setup()
                    try:
                        pia.initialize_4294a_for_cpd(inst)
                        file_management.ensure_output_dir(file_management.output_dir)
                        
                        print("Running εr vs frequency sweep...")
                        freq_axis, cp_f = pia.measure_eps_vs_freq(inst, freq_start, freq_stop, freq_points, bias_voltage)
                        eps_f = pia.compute_eps_r(cp_f, thickness_nm, diam_um)
                        
                        # Save εr(f) data
                        csv_data = np.column_stack([freq_axis, cp_f, eps_f])
                        file_management.save_csv(
                            os.path.join(file_management.output_dir, "eps_vs_freq.csv"),
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
                
                elif sub_choice == 3:
                    dielectric_submenu = False
                else:
                    print("Invalid choice. Please try again.")



if __name__ == "__main__":
    main()

