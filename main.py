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
import cryo
import file_management
import config
import logging_config
from gpib_utils import (
    InstrumentSession, prompt_choice, prompt_bool,
    safe_float_input, safe_int_input,
)

log = logging.getLogger(__name__)

# =============================
# Global Cryo Parameters
# =============================
cryo_enabled = False
cryo_params = None  # Will hold temperature sweep parameters

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


def run_cryo_sweep(measurement_queue):
    """
    Execute multiple measurements at each temperature point in a cryo sweep.

    The temperature ramp is paused (held) while measurements are running at
    each point, then resumed to move to the next setpoint.

    Data for each temperature point is saved into a subdirectory named by
    temperature (e.g. ``output/cryo_sweep/100.0K/``).

    Args:
        measurement_queue: list of (name, callable) tuples returned by
                           ``_prepare_*_cryo()`` functions.
    """
    global cryo_params

    if not cryo_enabled or cryo_params is None:
        print("Error: Cryogenic parameters are not configured.")
        return
    if not measurement_queue:
        print("Error: No measurements queued.")
        return

    temp_points = cryo_params['temp_points']
    ramp_rate = cryo_params['ramp_rate_k_per_min']
    n_temps = len(temp_points)
    n_meas = len(measurement_queue)
    base_output_dir = file_management.output_dir  # remember for restore

    print(f"\n{'='*60}")
    print("STARTING CRYOGENIC MEASUREMENT SWEEP")
    print(f"{'='*60}")
    print(f"Temperature sweep : {temp_points[0]:.1f}K → {temp_points[-1]:.1f}K")
    print(f"Ramp rate         : {ramp_rate:.2f} K/min")
    print(f"Temperature points: {n_temps}")
    print(f"Measurements/point: {n_meas}")
    for idx, (mname, _) in enumerate(measurement_queue, 1):
        print(f"  {idx}. {mname}")
    print()

    cryocon = None
    try:
        # Connect to Cryocon 32B
        cryocon = cryo.setup()

        for temp_idx, target_temp in enumerate(temp_points):
            print(f"\n{'─'*60}")
            print(f"Temperature point {temp_idx + 1}/{n_temps}: {target_temp:.2f} K")
            print(f"{'─'*60}")

            # Ramp to target temperature
            cryo.resume_ramp_to_next(cryocon, target_temp, ramp_rate, loop=2)

            # Wait for temperature to stabilise
            success = cryo.wait_for_temperature(
                cryocon, target_temp,
                tolerance_k=cryo.TEMP_TOLERANCE_K,
                stability_time_s=cryo.TEMP_STABILITY_TIME_S,
                loop=2)

            if not success:
                log.warning("Temperature stabilisation failed at point %d", temp_idx + 1)
                print(f"Warning: stabilisation failed at {target_temp:.2f} K – continuing anyway…")

            # HOLD temperature — pause ramp while measurements run
            cryo.hold_temperature(cryocon, loop=2)

            current_temp = cryo.get_current_temperature(cryocon, loop=2)
            if current_temp is not None:
                print(f"  Confirmed temperature: {current_temp:.2f} K")

            # Direct output into a temperature-specific subdirectory
            temp_subdir = os.path.join(base_output_dir, "cryo_sweep", f"{target_temp:.1f}K")
            file_management.set_output_dir(temp_subdir)
            file_management.ensure_output_dir(temp_subdir)

            # Execute every queued measurement at this temperature
            for m_idx, (mname, executor) in enumerate(measurement_queue, 1):
                print(f"\n  [{m_idx}/{n_meas}] {mname}")
                try:
                    executor()
                    print(f"  [{m_idx}/{n_meas}] Complete")
                except Exception as exc:
                    log.error("Measurement '%s' failed at %.1f K: %s", mname, target_temp, exc)
                    print(f"  [{m_idx}/{n_meas}] FAILED: {exc}")

            print(f"\nAll measurements done at {target_temp:.2f} K")

        # Finished sweep
        cryo.disable_ramp(cryocon, loop=2)
        cryo.disconnect_cryocon(cryocon)
        cryocon = None

        log.info("Cryogenic sweep complete – %d temps × %d measurements", n_temps, n_meas)
        print(f"\n{'='*60}")
        print("CRYOGENIC SWEEP COMPLETE")
        print(f"{'='*60}")
        print(f"Data saved under: {os.path.join(base_output_dir, 'cryo_sweep')}")

    except Exception as e:
        log.error("Cryogenic sweep error: %s", e)
        print(f"\nCritical error during sweep: {e}")
        raise
    finally:
        # Always restore the original output directory
        file_management.set_output_dir(base_output_dir)
        if cryocon is not None:
            try:
                cryo.disable_ramp(cryocon, loop=2)
                cryo.disconnect_cryocon(cryocon)
            except Exception:
                pass


def main():
    """Main entry point for the measurement automation system."""
    global cryo_enabled, cryo_params
    
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
        print("0. CRYO - Cryogenic Temperature Controller")
        print("1. PIA (Precision Impedance Analyzer)")
        print("2. PSPA (Parameter/Source Analyzer)")
        print("3. LCR (E4980A LCR Meter)")
        print("q. Quit")
        
        instrument_choice = input("\nEnter your choice (0, 1, 2, 3, or q): ").strip().lower()
        
        if instrument_choice == 'q':
            log.info("User exited program")
            print("Exiting the program. Goodbye!")
            break
        elif instrument_choice == '0':
            # ── Cryogenic temperature-controlled measurements ──
            cryo_params = cryo.get_cryogenic_parameters()
            cryo_enabled = True
            log.info("Cryogenic measurements enabled")

            # Build measurement queue
            measurement_queue = []   # list of (name, callable)

            _TOOL_MENU = {
                '1': ("PIA", pia_measurements, _prepare_pia_cryo),
                '2': ("PSPA", pspa_measurements, _prepare_pspa_cryo),
                '3': ("LCR", lcr_measurements, _prepare_lcr_cryo),
            }

            while True:
                print(f"\n{'='*60}")
                print("CRYO SWEEP – QUEUE MEASUREMENTS")
                print(f"{'='*60}")
                if measurement_queue:
                    print("Current queue:")
                    for qi, (qn, _) in enumerate(measurement_queue, 1):
                        print(f"  {qi}. {qn}")
                else:
                    print("Current queue: (empty)")
                print()
                print("1. Add PIA measurement(s)")
                print("2. Add PSPA measurement(s)")
                print("3. Add LCR measurement(s)")
                print("v. View / edit queue")
                print("s. Start sweep")
                print("b. Cancel and return to main menu")

                qchoice = input("\nEnter choice: ").strip().lower()

                if qchoice == 'b':
                    cryo_enabled = False
                    cryo_params = None
                    measurement_queue.clear()
                    break

                elif qchoice == 's':
                    if not measurement_queue:
                        print("Queue is empty – add at least one measurement first.")
                        continue
                    # Confirm and start
                    print(f"\nReady to sweep {len(cryo_params['temp_points'])} temperature "
                          f"points × {len(measurement_queue)} measurement(s).")
                    go = input("Press Enter to start (or 'b' to go back): ").strip().lower()
                    if go == 'b':
                        continue
                    run_cryo_sweep(measurement_queue)
                    # After sweep, reset cryo state
                    cryo_enabled = False
                    cryo_params = None
                    measurement_queue.clear()
                    break

                elif qchoice == 'v':
                    if not measurement_queue:
                        print("Queue is empty.")
                        continue
                    print("\nCurrent queue:")
                    for qi, (qn, _) in enumerate(measurement_queue, 1):
                        print(f"  {qi}. {qn}")
                    edit = input("\nEnter number to remove, 'c' to clear all, or 'b' to go back: ").strip().lower()
                    if edit == 'c':
                        measurement_queue.clear()
                        print("Queue cleared.")
                    elif edit == 'b':
                        pass
                    else:
                        try:
                            rm_idx = int(edit) - 1
                            if 0 <= rm_idx < len(measurement_queue):
                                removed = measurement_queue.pop(rm_idx)
                                print(f"Removed: {removed[0]}")
                            else:
                                print("Invalid index.")
                        except ValueError:
                            print("Invalid input.")

                elif qchoice in _TOOL_MENU:
                    tool_label, meas_list, prepare_fn = _TOOL_MENU[qchoice]
                    print(f"\n{tool_label} MEASUREMENTS")
                    for mi, mname in enumerate(meas_list, 1):
                        print(f"  {mi}. {mname}")
                    sel = input(f"\nSelect measurement(s) (comma-separated, e.g. 1,3) or 'b': ").strip().lower()
                    if sel == 'b':
                        continue
                    try:
                        choices = [int(x.strip()) for x in sel.split(',')]
                    except ValueError:
                        print("Invalid input – use numbers separated by commas.")
                        continue
                    for ch in choices:
                        if ch < 1 or ch > len(meas_list):
                            print(f"  Skipping invalid choice {ch}")
                            continue
                        name, executor = prepare_fn(ch)
                        if executor is None:
                            print(f"  '{name}' is not available in cryo mode – skipping.")
                        else:
                            measurement_queue.append((f"[{tool_label}] {name}", executor))
                            print(f"  Queued: [{tool_label}] {name}")
                else:
                    print("Invalid choice.")
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

    freq_cv = safe_float_input("Enter measurement frequency (Hz) [default 25000]: ", 25000)
    v_min = safe_float_input("Enter minimum voltage (V) [default 0]: ", 0)
    v_max = safe_float_input("Enter maximum voltage (V) [default 5]: ", 5)
    n_points = safe_int_input("Enter number of points per sweep [default 401]: ", 401)
    n_cycles = safe_int_input("Enter number of cycles [default 1]: ", 1)
    thickness_nm = safe_float_input("Enter HZO thickness (nm) [default 10.0]: ", 10.0)
    diam_um = safe_float_input("Enter electrode diameter (µm) [default 75.0]: ", 75.0)

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


# =============================
# Cryo Prepare Functions
# =============================
# Each _prepare_*_cryo(choice) function collects measurement parameters
# upfront and returns (name, callable) where the callable executes the
# measurement with the pre-bound params.  Used by run_cryo_sweep().

def _prepare_pia_cryo(choice):
    """Collect PIA parameters and return (name, executor) for cryo queue."""
    name = pia_measurements[choice - 1]
    print(f"\n  Configuring: {name}")

    if choice in (1, 2):
        # Impedance Magnitude / Phase vs Frequency
        freq_start, freq_stop, num_points, apply_dc_bias, dc_bias_v = _pia_get_params()

        def run():
            with InstrumentSession(pia.setup, _pia_safe_close) as inst:
                pia.initialize_4294a_for_impedance(inst)
                pia.configure_dc_bias(inst, apply_dc_bias, dc_bias_v)
                freq_axis, z_mag, theta_deg = pia.measure_impedance_vs_freq(
                    inst, freq_start, freq_stop, num_points)
                y_label = "|Z| (Ohm)" if choice == 1 else "Phase (Degrees)"
                y_data = z_mag if choice == 1 else theta_deg
                file_management.save_image(
                    name, "Frequency (Hz)", freq_axis, y_label, y_data,
                    APPLY_DC_BIAS=apply_dc_bias, DC_BIAS_V=dc_bias_v)
                z_real = z_mag * np.cos(np.deg2rad(theta_deg))
                z_imag = z_mag * np.sin(np.deg2rad(theta_deg))
                file_management.save_csv("impedance_data.csv",
                    np.column_stack([freq_axis, z_mag, theta_deg, z_real, z_imag]),
                    "frequency_Hz, Z_mag_ohm, theta_deg, Z_real_ohm, Z_imag_ohm")
        return name, run

    elif choice in (3, 4):
        # Capacitance / Tan Loss vs Frequency
        freq_start, freq_stop, num_points, apply_dc_bias, dc_bias_v = _pia_get_params()

        def run():
            with InstrumentSession(pia.setup, _pia_safe_close) as inst:
                pia.initialize_4294a_for_cpd(inst)
                pia.configure_dc_bias(inst, apply_dc_bias, dc_bias_v)
                freq_axis, cp_vals, d_vals = pia.measure_cpd_vs_freq(
                    inst, freq_start, freq_stop, num_points)
                y_label = "Cp (F)" if choice == 3 else "D (loss factor)"
                y_data = cp_vals if choice == 3 else d_vals
                file_management.save_image(
                    name, "Frequency (Hz)", freq_axis, y_label, y_data,
                    APPLY_DC_BIAS=apply_dc_bias, DC_BIAS_V=dc_bias_v)
                file_management.save_csv("cpd_data.csv",
                    np.column_stack([freq_axis, cp_vals, d_vals]),
                    "frequency_Hz, Cp_F, D")
        return name, run

    elif choice in (5, 6):
        # C–V Butterfly: Cp / εr vs Voltage
        plot_mode = "cp" if choice == 5 else "eps"
        freq_cv = safe_float_input("Enter measurement frequency (Hz) [default 25000]: ", 25000)
        v_min = safe_float_input("Enter minimum voltage (V) [default 0]: ", 0)
        v_max = safe_float_input("Enter maximum voltage (V) [default 5]: ", 5)
        n_points = safe_int_input("Enter number of points per sweep [default 401]: ", 401)
        n_cycles = safe_int_input("Enter number of cycles [default 1]: ", 1)
        thickness_nm = safe_float_input("Enter HZO thickness (nm) [default 10.0]: ", 10.0)
        diam_um = safe_float_input("Enter electrode diameter (µm) [default 75.0]: ", 75.0)

        def run():
            with InstrumentSession(pia.setup, _pia_safe_close) as inst:
                pia.initialize_4294a_for_cpd(inst)
                all_v, all_cp, all_eps, all_idx = [], [], [], []
                for cyc in range(n_cycles):
                    v, cp = pia.measure_single_cv_cycle(inst, freq_cv, v_min, v_max, n_points)
                    eps_r = pia.compute_eps_r(cp, thickness_nm, diam_um)
                    all_v.append(v); all_cp.append(cp); all_eps.append(eps_r)
                    all_idx.append(np.full_like(v, cyc + 1, dtype=int))
                file_management.save_csv("cv_dielectric_cycles.csv",
                    np.column_stack([np.concatenate(all_idx), np.concatenate(all_v),
                                     np.concatenate(all_cp), np.concatenate(all_eps)]),
                    "cycle_index, bias_V, Cp_F, eps_r")
                if plot_mode == "cp":
                    file_management.save_cycle_plot("C-V Butterfly: Cp vs Voltage",
                        "Bias Voltage (V)", "Cp (F)", all_v, all_cp, "cv_butterfly_cp.png")
                else:
                    file_management.save_cycle_plot("C-V Butterfly: εr vs Voltage",
                        "Bias Voltage (V)", "Dielectric constant εr",
                        all_v, all_eps, "cv_butterfly_eps_r.png")
        return name, run

    elif choice == 7:
        # Permittivity vs Frequency
        if (pia.current_parameters["freq_start_hz"] != pia.FREQ_START_HZ
                or pia.current_parameters["apply_dc_bias"] != pia.APPLY_DC_BIAS):
            use_last = pia.prompt_for_parameter_duplication()
        else:
            use_last = False
        if use_last:
            freq_start = pia.current_parameters["freq_start_hz"]
            freq_stop = pia.current_parameters["freq_stop_hz"]
            freq_points = pia.current_parameters["num_points"]
            bias_voltage = pia.current_parameters["dc_bias_v"]
        else:
            freq_start = safe_float_input("Enter start frequency (Hz) [default 1000]: ", 1000)
            freq_stop = safe_float_input("Enter stop frequency (Hz) [default 1e6]: ", 1e6)
            freq_points = safe_int_input("Enter number of frequency points [default 201]: ", 201)
            bias_voltage = safe_float_input("Enter DC bias voltage (V) [default 0]: ", 0)
            pia.current_parameters.update({
                "freq_start_hz": freq_start, "freq_stop_hz": freq_stop,
                "num_points": freq_points, "dc_bias_v": bias_voltage,
                "apply_dc_bias": (bias_voltage != 0)})
        thickness_nm = safe_float_input("Enter HZO thickness (nm) [default 10.0]: ", 10.0)
        diam_um = safe_float_input("Enter electrode diameter (µm) [default 75.0]: ", 75.0)

        def run():
            with InstrumentSession(pia.setup, _pia_safe_close) as inst:
                pia.initialize_4294a_for_cpd(inst)
                freq_axis, cp_f = pia.measure_eps_vs_freq(
                    inst, freq_start, freq_stop, freq_points, bias_voltage)
                eps_f = pia.compute_eps_r(cp_f, thickness_nm, diam_um)
                file_management.save_csv("eps_vs_freq.csv",
                    np.column_stack([freq_axis, cp_f, eps_f]),
                    "frequency_Hz, Cp_F, eps_r")
                file_management.save_image("Permittivity vs Frequency",
                    "Frequency (Hz)", freq_axis, "Dielectric constant εr", eps_f,
                    APPLY_DC_BIAS=(bias_voltage != 0), DC_BIAS_V=bias_voltage)
        return name, run

    elif choice == 8:
        # R-X vs Frequency
        freq_start, freq_stop, num_points, apply_dc_bias, dc_bias_v = _pia_get_params()

        def run():
            with InstrumentSession(pia.setup, _pia_safe_close) as inst:
                pia.initialize_4294a_for_rx(inst)
                pia.configure_dc_bias(inst, apply_dc_bias, dc_bias_v)
                freq_axis, r_vals, x_vals = pia.measure_rx_vs_freq(
                    inst, freq_start, freq_stop, num_points)
                file_management.save_image("Resistance vs Frequency",
                    "Frequency (Hz)", freq_axis, "R (Ohm)", r_vals,
                    APPLY_DC_BIAS=apply_dc_bias, DC_BIAS_V=dc_bias_v)
                file_management.save_image("Reactance vs Frequency",
                    "Frequency (Hz)", freq_axis, "X (Ohm)", x_vals,
                    APPLY_DC_BIAS=apply_dc_bias, DC_BIAS_V=dc_bias_v)
                file_management.save_csv("rx_data.csv",
                    np.column_stack([freq_axis, r_vals, x_vals]),
                    "frequency_Hz, R_ohm, X_ohm")
        return name, run

    elif choice == 9:
        # G-B vs Frequency
        freq_start, freq_stop, num_points, apply_dc_bias, dc_bias_v = _pia_get_params()

        def run():
            with InstrumentSession(pia.setup, _pia_safe_close) as inst:
                pia.initialize_4294a_for_gb(inst)
                pia.configure_dc_bias(inst, apply_dc_bias, dc_bias_v)
                freq_axis, g_vals, b_vals = pia.measure_gb_vs_freq(
                    inst, freq_start, freq_stop, num_points)
                file_management.save_image("Conductance vs Frequency",
                    "Frequency (Hz)", freq_axis, "G (S)", g_vals,
                    APPLY_DC_BIAS=apply_dc_bias, DC_BIAS_V=dc_bias_v)
                file_management.save_image("Susceptance vs Frequency",
                    "Frequency (Hz)", freq_axis, "B (S)", b_vals,
                    APPLY_DC_BIAS=apply_dc_bias, DC_BIAS_V=dc_bias_v)
                file_management.save_csv("gb_data.csv",
                    np.column_stack([freq_axis, g_vals, b_vals]),
                    "frequency_Hz, G_S, B_S")
        return name, run

    elif choice == 10:
        # Y-θ vs Frequency
        freq_start, freq_stop, num_points, apply_dc_bias, dc_bias_v = _pia_get_params()

        def run():
            with InstrumentSession(pia.setup, _pia_safe_close) as inst:
                pia.initialize_4294a_for_ytd(inst)
                pia.configure_dc_bias(inst, apply_dc_bias, dc_bias_v)
                freq_axis, y_mag, y_theta = pia.measure_ytd_vs_freq(
                    inst, freq_start, freq_stop, num_points)
                file_management.save_image("Admittance Magnitude vs Frequency",
                    "Frequency (Hz)", freq_axis, "|Y| (S)", y_mag,
                    APPLY_DC_BIAS=apply_dc_bias, DC_BIAS_V=dc_bias_v)
                file_management.save_image("Admittance Phase vs Frequency",
                    "Frequency (Hz)", freq_axis, "θ (Degrees)", y_theta,
                    APPLY_DC_BIAS=apply_dc_bias, DC_BIAS_V=dc_bias_v)
                file_management.save_csv("ytd_data.csv",
                    np.column_stack([freq_axis, y_mag, y_theta]),
                    "frequency_Hz, Y_mag_S, theta_deg")
        return name, run

    else:
        log.warning("PIA choice %d not supported in cryo mode", choice)
        return name, None


def _prepare_pspa_cryo(choice):
    """Collect PSPA parameters and return (name, executor) for cryo queue."""
    name = pspa_measurements[choice - 1]
    print(f"\n  Configuring: {name}")

    if choice == 1:
        # Transistor Output Characteristics
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

        def run():
            with InstrumentSession(pspa.connect_pspa, pspa.disconnect_pspa) as inst:
                data = pspa.measure_transistor_output_characteristics(
                    inst, vds_start, vds_stop, vds_step,
                    vgs_start, vgs_stop, vgs_step,
                    drain_ch, gate_ch, source_ch, compliance)
                file_management.save_csv("transistor_output_chars.csv",
                    np.column_stack([data['Vds'], data['Vgs'], data['Id']]),
                    "Vds_V, Vgs_V, Id_A")
                plt.figure(figsize=(10, 6))
                for vgs in np.unique(data['Vgs']):
                    m = data['Vgs'] == vgs
                    plt.plot(data['Vds'][m], data['Id'][m] * 1e3,
                             marker='o', label=f"Vgs = {vgs:.2f} V")
                plt.xlabel('Vds (V)'); plt.ylabel('Id (mA)')
                plt.title('Transistor Output Characteristics')
                plt.legend(); plt.grid(True); plt.tight_layout()
                file_management.save_plot("transistor_output_chars.png"); plt.close()
        return name, run

    elif choice == 2:
        # Transistor Transfer Characteristics
        vgs_start = safe_float_input("Enter start Vgs (V) [default -1]: ", -1)
        vgs_stop = safe_float_input("Enter stop Vgs (V) [default 3]: ", 3)
        vgs_step = safe_float_input("Enter Vgs step (V) [default 0.05]: ", 0.05)
        vds_constant = safe_float_input("Enter constant Vds (V) [default 5]: ", 5)
        compliance = safe_float_input("Enter current compliance (A) [default 0.1]: ", 0.1)
        drain_ch = safe_int_input("Enter drain channel [default 1]: ", 1)
        gate_ch = safe_int_input("Enter gate channel [default 2]: ", 2)
        source_ch = safe_int_input("Enter source channel [default 3]: ", 3)

        def run():
            with InstrumentSession(pspa.connect_pspa, pspa.disconnect_pspa) as inst:
                data = pspa.measure_transistor_transfer_characteristics(
                    inst, vgs_start, vgs_stop, vgs_step, vds_constant,
                    drain_ch, gate_ch, source_ch, compliance)
                file_management.save_csv("transistor_transfer_chars.csv",
                    np.column_stack([data['Vgs'], data['Id'], data['Ig']]),
                    "Vgs_V, Id_A, Ig_A")
                fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
                ax1.plot(data['Vgs'], data['Id'] * 1e3, marker='o', label='Id')
                ax1.set_xlabel('Vgs (V)'); ax1.set_ylabel('Id (mA)')
                ax1.set_title(f'Transfer Chars (Vds = {vds_constant} V) - Linear')
                ax1.grid(True); ax1.legend()
                ax2.semilogy(data['Vgs'], np.abs(data['Id']), marker='o', label='|Id|')
                ax2.semilogy(data['Vgs'], np.abs(data['Ig']), marker='s', label='|Ig|')
                ax2.set_xlabel('Vgs (V)'); ax2.set_ylabel('Current (A)')
                ax2.set_title(f'Transfer Chars (Vds = {vds_constant} V) - Log')
                ax2.grid(True); ax2.legend()
                plt.tight_layout()
                file_management.save_plot("transistor_transfer_chars.png", fig); plt.close()
        return name, run

    elif choice == 3:
        # I-V Curve (Unidirectional)
        v_start = safe_float_input("Enter start voltage (V) [default 0]: ", 0)
        v_stop = safe_float_input("Enter stop voltage (V) [default 5]: ", 5)
        v_step = safe_float_input("Enter voltage step (V) [default 0.1]: ", 0.1)
        compliance = safe_float_input("Enter current compliance (A) [default 0.1]: ", 0.1)
        channel = safe_int_input("Enter channel number [default 1]: ", 1)

        def run():
            with InstrumentSession(pspa.connect_pspa, pspa.disconnect_pspa) as inst:
                data = pspa.measure_iv_curve(inst, v_start, v_stop, v_step, channel, compliance)
                file_management.save_csv("iv_curve.csv",
                    np.column_stack([data['Voltage'], data['Current']]),
                    "Voltage_V, Current_A")
                plt.figure(figsize=(10, 6))
                plt.plot(data['Voltage'], data['Current'] * 1e3, marker='o')
                plt.xlabel('Voltage (V)'); plt.ylabel('Current (mA)')
                plt.title('I-V Curve (Unidirectional)')
                plt.grid(True); plt.tight_layout()
                file_management.save_plot("iv_curve.png"); plt.close()
        return name, run

    elif choice == 4:
        # I-V Curve (Bidirectional)
        v_max = safe_float_input("Enter maximum voltage magnitude (V) [default 5]: ", 5)
        v_step = safe_float_input("Enter voltage step (V) [default 0.1]: ", 0.1)
        compliance = safe_float_input("Enter current compliance (A) [default 0.1]: ", 0.1)
        channel = safe_int_input("Enter channel number [default 1]: ", 1)

        def run():
            with InstrumentSession(pspa.connect_pspa, pspa.disconnect_pspa) as inst:
                data = pspa.measure_iv_bidirectional(inst, v_max, v_step, channel, compliance)
                file_management.save_csv("iv_curve_bidirectional.csv",
                    np.column_stack([data['Voltage'], data['Current']]),
                    "Voltage_V, Current_A")
                plt.figure(figsize=(10, 6))
                plt.plot(data['Voltage'], data['Current'] * 1e3, marker='o')
                plt.xlabel('Voltage (V)'); plt.ylabel('Current (mA)')
                plt.title('I-V Curve (Bidirectional)')
                plt.grid(True); plt.tight_layout()
                file_management.save_plot("iv_curve_bidirectional.png"); plt.close()
        return name, run

    elif choice == 5:
        # Pulsed I-V (Single Device)
        v_base = safe_float_input("Enter base voltage (V) [default 0]: ", 0)
        v_pulse = safe_float_input("Enter pulse voltage (V) [default 5]: ", 5)
        pulse_width = safe_float_input("Enter pulse width (µs) [default 100]: ", 100) * 1e-6
        pulse_period = safe_float_input("Enter pulse period (ms) [default 10]: ", 10) * 1e-3
        num_pulses = safe_int_input("Enter number of pulses [default 10]: ", 10)
        compliance = safe_float_input("Enter current compliance (A) [default 0.1]: ", 0.1)
        channel = safe_int_input("Enter channel number [default 1]: ", 1)

        def run():
            with InstrumentSession(pspa.connect_pspa, pspa.disconnect_pspa) as inst:
                data = pspa.measure_pulsed_iv(
                    inst, v_base, v_pulse, pulse_width, pulse_period,
                    num_pulses, channel, compliance)
                file_management.save_csv("pulsed_iv.csv",
                    np.column_stack([data['Time'], data['Voltage'], data['Current']]),
                    "Time_s, Voltage_V, Current_A")
                fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
                ax1.plot(data['Time'] * 1e3, data['Current'] * 1e3, marker='o')
                ax1.set_xlabel('Time (ms)'); ax1.set_ylabel('Current (mA)')
                ax1.set_title('Pulsed I-V: Current Response'); ax1.grid(True)
                ax2.plot(data['Time'] * 1e3, data['Voltage'], marker='o')
                ax2.set_xlabel('Time (ms)'); ax2.set_ylabel('Voltage (V)')
                ax2.set_title('Pulsed I-V: Applied Voltage'); ax2.grid(True)
                plt.tight_layout()
                file_management.save_plot("pulsed_iv.png", fig); plt.close()
        return name, run

    elif choice == 6:
        # Pulsed Transistor
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

        def run():
            with InstrumentSession(pspa.connect_pspa, pspa.disconnect_pspa) as inst:
                data = pspa.measure_pulsed_transistor(
                    inst, vds_pulse, vgs_pulse, vds_base, vgs_base,
                    pulse_width, pulse_period, num_pulses,
                    drain_ch, gate_ch, source_ch, compliance)
                file_management.save_csv("pulsed_transistor.csv",
                    np.column_stack([data['Time'], data['Id']]),
                    "Time_s, Id_A")
                plt.figure(figsize=(10, 6))
                plt.plot(data['Time'] * 1e3, data['Id'] * 1e3, marker='o', label='Id')
                plt.xlabel('Time (ms)'); plt.ylabel('Id (mA)')
                plt.title(f'Pulsed Transistor (Vds={vds_pulse}V, Vgs={vgs_pulse}V)')
                plt.grid(True); plt.legend(); plt.tight_layout()
                file_management.save_plot("pulsed_transistor.png"); plt.close()
        return name, run

    elif choice == 7:
        # Diode I-V Characteristics
        v_start = safe_float_input("Enter start voltage (V) [default -1]: ", -1)
        v_stop = safe_float_input("Enter stop voltage (V) [default 1]: ", 1)
        v_step = safe_float_input("Enter voltage step (V) [default 0.02]: ", 0.02)
        compliance = safe_float_input("Enter current compliance (A) [default 0.1]: ", 0.1)
        anode_ch = safe_int_input("Enter anode channel [default 1]: ", 1)
        cathode_ch = safe_int_input("Enter cathode channel [default 2]: ", 2)

        def run():
            with InstrumentSession(pspa.connect_pspa, pspa.disconnect_pspa) as inst:
                data = pspa.measure_diode_iv(
                    inst, v_start, v_stop, v_step, anode_ch, cathode_ch, compliance)
                file_management.save_csv("diode_iv.csv",
                    np.column_stack([data['Voltage'], data['Current']]),
                    "Voltage_V, Current_A")
                fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
                ax1.plot(data['Voltage'], data['Current'] * 1e3, marker='o')
                ax1.set_xlabel('Voltage (V)'); ax1.set_ylabel('Current (mA)')
                ax1.set_title('Diode I-V (Linear)'); ax1.grid(True)
                ax2.semilogy(data['Voltage'], np.abs(data['Current']), marker='o')
                ax2.set_xlabel('Voltage (V)'); ax2.set_ylabel('|Current| (A)')
                ax2.set_title('Diode I-V (Log)'); ax2.grid(True)
                plt.tight_layout()
                file_management.save_plot("diode_iv.png"); plt.close()
        return name, run

    elif choice == 8:
        # Gate Leakage Current
        vgs_start = safe_float_input("Enter start Vgs (V) [default 0]: ", 0)
        vgs_stop = safe_float_input("Enter stop Vgs (V) [default 5]: ", 5)
        vgs_step = safe_float_input("Enter Vgs step (V) [default 0.1]: ", 0.1)
        compliance = safe_float_input("Enter current compliance (A) [default 1e-3]: ", 1e-3)
        gate_ch = safe_int_input("Enter gate channel [default 2]: ", 2)
        source_ch = safe_int_input("Enter source channel [default 3]: ", 3)

        def run():
            with InstrumentSession(pspa.connect_pspa, pspa.disconnect_pspa) as inst:
                data = pspa.measure_gate_leakage(
                    inst, vgs_start, vgs_stop, vgs_step, gate_ch, source_ch, compliance)
                file_management.save_csv("gate_leakage.csv",
                    np.column_stack([data['Vgs'], data['Ig']]),
                    "Vgs_V, Ig_A")
                plt.figure(figsize=(10, 6))
                plt.semilogy(data['Vgs'], np.abs(data['Ig']), marker='o')
                plt.xlabel('Vgs (V)'); plt.ylabel('|Ig| (A)')
                plt.title('Gate Leakage Current')
                plt.grid(True); plt.tight_layout()
                file_management.save_plot("gate_leakage.png"); plt.close()
        return name, run

    elif choice == 9:
        # Breakdown Voltage
        v_start = safe_float_input("Enter start voltage (V) [default 0]: ", 0)
        v_stop = safe_float_input("Enter stop voltage (V) [default 40]: ", 40)
        v_step = safe_float_input("Enter voltage step (V) [default 0.5]: ", 0.5)
        compliance = safe_float_input("Enter current compliance (A) [default 1e-3]: ", 1e-3)
        threshold = safe_float_input("Enter breakdown current threshold (A) [default 1e-4]: ", 1e-4)
        channel = safe_int_input("Enter channel number [default 1]: ", 1)

        def run():
            with InstrumentSession(pspa.connect_pspa, pspa.disconnect_pspa) as inst:
                data = pspa.measure_breakdown_voltage(
                    inst, v_start, v_stop, v_step, channel, compliance, threshold)
                file_management.save_csv("breakdown_voltage.csv",
                    np.column_stack([data['Voltage'], data['Current']]),
                    "Voltage_V, Current_A")
                plt.figure(figsize=(10, 6))
                plt.semilogy(data['Voltage'], np.abs(data['Current']), marker='o')
                plt.xlabel('Voltage (V)'); plt.ylabel('|Current| (A)')
                if data['breakdown_v'] is not None:
                    plt.axvline(x=data['breakdown_v'], color='r', linestyle='--',
                                label=f"Breakdown = {data['breakdown_v']:.2f} V")
                    plt.legend()
                plt.title('Breakdown Voltage')
                plt.grid(True); plt.tight_layout()
                file_management.save_plot("breakdown_voltage.png"); plt.close()
                if data['breakdown_v'] is not None:
                    print(f"  Breakdown voltage: {data['breakdown_v']:.3f} V")
        return name, run

    elif choice == 10:
        # Resistance Measurement
        i_start = safe_float_input("Enter start current (A) [default 0]: ", 0)
        i_stop = safe_float_input("Enter stop current (A) [default 1e-3]: ", 1e-3)
        i_step = safe_float_input("Enter current step (A) [default 1e-4]: ", 1e-4)
        compliance = safe_float_input("Enter voltage compliance (V) [default 10.0]: ", 10.0)
        channel = safe_int_input("Enter channel number [default 1]: ", 1)

        def run():
            with InstrumentSession(pspa.connect_pspa, pspa.disconnect_pspa) as inst:
                data = pspa.measure_resistance(
                    inst, i_start, i_stop, i_step, channel, compliance)
                file_management.save_csv("resistance.csv",
                    np.column_stack([data['Current'], data['Voltage']]),
                    "Current_A, Voltage_V")
                plt.figure(figsize=(10, 6))
                plt.plot(data['Current'] * 1e3, data['Voltage'], marker='o')
                plt.xlabel('Current (mA)'); plt.ylabel('Voltage (V)')
                plt.title(f"Resistance (R ≈ {data['resistance_ohm']:.4e} Ω)")
                plt.grid(True); plt.tight_layout()
                file_management.save_plot("resistance.png"); plt.close()
                print(f"  Measured Resistance: {data['resistance_ohm']:.4e} Ω")
        return name, run

    else:
        log.warning("PSPA choice %d not supported in cryo mode", choice)
        return name, None


def _prepare_lcr_cryo(choice):
    """Collect LCR parameters and return (name, executor) for cryo queue."""
    name = lcr_measurements[choice - 1]
    print(f"\n  Configuring: {name}")

    if choice == 1:
        # Impedance vs Frequency
        freq_start = safe_float_input("Enter start frequency (Hz) [default 20]: ", 20)
        freq_stop = safe_float_input("Enter stop frequency (Hz) [default 2e6]: ", 2e6)
        num_points = safe_int_input("Enter number of points [default 201]: ", 201)
        sweep_type = prompt_choice("Sweep type", ["LIN", "LOG"], "LOG")
        ac_level = safe_float_input("Enter AC level (V) [default 1.0]: ", 1.0)
        apply_bias = prompt_bool("Apply DC bias?", False)
        bias_v = safe_float_input("Enter DC bias (V) [default 0]: ", 0) if apply_bias else 0.0

        def run():
            with InstrumentSession(
                lambda: lcr.setup(measurement_function="ZTD"),
                lcr.disconnect_e4980a) as inst:
                lcr.set_ac_level(inst, ac_level)
                lcr.configure_dc_bias(inst, apply_bias, bias_v)
                freq_axis, z_mag, theta_deg = lcr.measure_impedance_vs_frequency(
                    inst, freq_start, freq_stop, num_points, sweep_type)
                file_management.save_image(
                    "Impedance Magnitude vs Frequency (LCR)", "Frequency (Hz)",
                    freq_axis, "|Z| (Ohm)", z_mag,
                    APPLY_DC_BIAS=apply_bias, DC_BIAS_V=bias_v)
                file_management.save_image(
                    "Impedance Phase vs Frequency (LCR)", "Frequency (Hz)",
                    freq_axis, "Phase (Degrees)", theta_deg,
                    APPLY_DC_BIAS=apply_bias, DC_BIAS_V=bias_v)
                z_real, z_imag = lcr.compute_impedance_components(z_mag, theta_deg)
                file_management.save_csv("lcr_impedance_vs_freq.csv",
                    np.column_stack([freq_axis, z_mag, theta_deg, z_real, z_imag]),
                    "frequency_Hz, Z_mag_ohm, theta_deg, Z_real_ohm, Z_imag_ohm")
        return name, run

    elif choice == 2:
        # Capacitance vs Frequency
        freq_start = safe_float_input("Enter start frequency (Hz) [default 20]: ", 20)
        freq_stop = safe_float_input("Enter stop frequency (Hz) [default 2e6]: ", 2e6)
        num_points = safe_int_input("Enter number of points [default 201]: ", 201)
        sweep_type = prompt_choice("Sweep type", ["LIN", "LOG"], "LOG")
        mode = prompt_choice("Measurement mode", ["CPD", "CSRS", "CPRP"], "CPD")
        ac_level = safe_float_input("Enter AC level (V) [default 1.0]: ", 1.0)
        apply_bias = prompt_bool("Apply DC bias?", False)
        bias_v = safe_float_input("Enter DC bias (V) [default 0]: ", 0) if apply_bias else 0.0

        def run():
            with InstrumentSession(
                lambda: lcr.setup(measurement_function=mode),
                lcr.disconnect_e4980a) as inst:
                lcr.set_ac_level(inst, ac_level)
                lcr.configure_dc_bias(inst, apply_bias, bias_v)
                freq_axis, cap_vals, sec_vals = lcr.measure_capacitance_vs_frequency(
                    inst, freq_start, freq_stop, num_points, sweep_type, mode)
                sec_label = "D (loss factor)" if mode == "CPD" else "Secondary"
                file_management.save_image(
                    f"Capacitance vs Frequency ({mode}) (LCR)", "Frequency (Hz)",
                    freq_axis, "Capacitance (F)", cap_vals,
                    APPLY_DC_BIAS=apply_bias, DC_BIAS_V=bias_v)
                file_management.save_image(
                    f"{sec_label} vs Frequency (LCR)", "Frequency (Hz)",
                    freq_axis, sec_label, sec_vals,
                    APPLY_DC_BIAS=apply_bias, DC_BIAS_V=bias_v)
                file_management.save_csv(
                    f"lcr_capacitance_vs_freq_{mode.lower()}.csv",
                    np.column_stack([freq_axis, cap_vals, sec_vals]),
                    f"frequency_Hz, capacitance_F, {sec_label.replace(' ', '_')}")
        return name, run

    elif choice == 3:
        # Inductance vs Frequency
        freq_start = safe_float_input("Enter start frequency (Hz) [default 20]: ", 20)
        freq_stop = safe_float_input("Enter stop frequency (Hz) [default 2e6]: ", 2e6)
        num_points = safe_int_input("Enter number of points [default 201]: ", 201)
        sweep_type = prompt_choice("Sweep type", ["LIN", "LOG"], "LOG")
        mode = prompt_choice("Measurement mode", ["LPQ", "LSD"], "LPQ")
        ac_level = safe_float_input("Enter AC level (V) [default 1.0]: ", 1.0)

        def run():
            with InstrumentSession(
                lambda: lcr.setup(measurement_function=mode),
                lcr.disconnect_e4980a) as inst:
                lcr.set_ac_level(inst, ac_level)
                freq_axis, ind_vals, sec_vals = lcr.measure_frequency_sweep(
                    inst, freq_start, freq_stop, num_points, sweep_type, mode)
                sec_label = "Q (quality factor)" if mode == "LPQ" else "D (loss factor)"
                file_management.save_image(
                    f"Inductance vs Frequency ({mode}) (LCR)", "Frequency (Hz)",
                    freq_axis, "Inductance (H)", ind_vals,
                    APPLY_DC_BIAS=False, DC_BIAS_V=0.0)
                file_management.save_image(
                    f"{sec_label} vs Frequency (LCR)", "Frequency (Hz)",
                    freq_axis, sec_label, sec_vals,
                    APPLY_DC_BIAS=False, DC_BIAS_V=0.0)
                file_management.save_csv(
                    f"lcr_inductance_vs_freq_{mode.lower()}.csv",
                    np.column_stack([freq_axis, ind_vals, sec_vals]),
                    f"frequency_Hz, inductance_H, {sec_label.replace(' ', '_')}")
        return name, run

    elif choice == 4:
        # C-V Sweep (Single)
        freq = safe_float_input("Enter measurement frequency (Hz) [default 1000]: ", 1000)
        v_min = safe_float_input("Enter minimum voltage (V) [default -5]: ", -5)
        v_max = safe_float_input("Enter maximum voltage (V) [default 5]: ", 5)
        num_points = safe_int_input("Enter number of points [default 201]: ", 201)
        ac_level = safe_float_input("Enter AC level (V) [default 0.1]: ", 0.1)

        def run():
            with InstrumentSession(
                lambda: lcr.setup(measurement_function="CPD"),
                lcr.disconnect_e4980a) as inst:
                lcr.set_ac_level(inst, ac_level)
                voltage, capacitance, dissipation = lcr.measure_cv_sweep(
                    inst, freq, v_min, v_max, num_points)
                plt.figure(figsize=(10, 6))
                plt.plot(voltage, capacitance, '-o', markersize=3)
                plt.xlabel('Voltage (V)'); plt.ylabel('Capacitance (F)')
                plt.title(f'C-V Sweep at {freq:.0f} Hz (LCR)')
                plt.grid(True)
                file_management.save_plot("lcr_cv_sweep.png"); plt.close()
                file_management.save_csv("lcr_cv_sweep.csv",
                    np.column_stack([voltage, capacitance, dissipation]),
                    "voltage_V, capacitance_F, dissipation_D")
        return name, run

    elif choice == 5:
        # C-V Butterfly Cycles
        freq = safe_float_input("Enter measurement frequency (Hz) [default 1000]: ", 1000)
        v_min = safe_float_input("Enter minimum voltage (V) [default -5]: ", -5)
        v_max = safe_float_input("Enter maximum voltage (V) [default 5]: ", 5)
        num_points = safe_int_input("Enter number of points per direction [default 201]: ", 201)
        num_cycles = safe_int_input("Enter number of cycles [default 1]: ", 1)
        ac_level = safe_float_input("Enter AC level (V) [default 0.1]: ", 0.1)
        calc_eps = prompt_bool("Calculate permittivity?", False)
        thickness_nm = 0.0
        diameter_um = 0.0
        if calc_eps:
            thickness_nm = safe_float_input("Enter dielectric thickness (nm) [default 10]: ", 10)
            diameter_um = safe_float_input("Enter electrode diameter (µm) [default 75]: ", 75)

        def run():
            with InstrumentSession(
                lambda: lcr.setup(measurement_function="CPD"),
                lcr.disconnect_e4980a) as inst:
                lcr.set_ac_level(inst, ac_level)
                all_v, all_cp, all_d, all_eps_list = [], [], [], []
                for cyc in range(num_cycles):
                    voltage, capacitance, dissipation = lcr.measure_cv_butterfly(
                        inst, freq, v_min, v_max, num_points)
                    all_v.append(voltage); all_cp.append(capacitance); all_d.append(dissipation)
                    if calc_eps:
                        all_eps_list.append(lcr.compute_eps_r(capacitance, thickness_nm, diameter_um))
                file_management.save_cycle_plot(
                    f"C-V Butterfly at {freq:.0f} Hz (LCR)", "Voltage (V)",
                    "Capacitance (F)", all_v, all_cp, "lcr_cv_butterfly.png")
                if calc_eps:
                    file_management.save_cycle_plot(
                        f"Permittivity Butterfly at {freq:.0f} Hz (LCR)", "Voltage (V)",
                        "Relative Permittivity εr", all_v, all_eps_list, "lcr_eps_butterfly.png")
                v_all = np.concatenate(all_v)
                cp_all = np.concatenate(all_cp)
                d_all = np.concatenate(all_d)
                cyc_idx = np.concatenate([np.full(len(v), i+1) for i, v in enumerate(all_v)])
                if calc_eps:
                    eps_all = np.concatenate(all_eps_list)
                    csv_data = np.column_stack([cyc_idx, v_all, cp_all, d_all, eps_all])
                    header = "cycle, voltage_V, capacitance_F, dissipation_D, eps_r"
                else:
                    csv_data = np.column_stack([cyc_idx, v_all, cp_all, d_all])
                    header = "cycle, voltage_V, capacitance_F, dissipation_D"
                file_management.save_csv("lcr_cv_butterfly.csv", csv_data, header)
        return name, run

    elif choice == 6:
        # Single Point Impedance
        freq = safe_float_input("Enter frequency (Hz) [default 1000]: ", 1000)
        ac_level = safe_float_input("Enter AC level (V) [default 1.0]: ", 1.0)
        apply_bias = prompt_bool("Apply DC bias?", False)
        bias_v = safe_float_input("Enter DC bias (V) [default 0]: ", 0) if apply_bias else 0.0

        def run():
            with InstrumentSession(
                lambda: lcr.setup(measurement_function="ZTD"),
                lcr.disconnect_e4980a) as inst:
                lcr.set_ac_level(inst, ac_level)
                lcr.configure_dc_bias(inst, apply_bias, bias_v)
                z_mag, theta_deg = lcr.measure_impedance(inst, freq)
                z_real, z_imag = lcr.compute_impedance_components(z_mag, theta_deg)
                print(f"  |Z|={z_mag:.6e} Ω  θ={theta_deg:.2f}°  "
                      f"Z_real={z_real:.6e} Ω  Z_imag={z_imag:.6e} Ω")
                file_management.save_csv("lcr_single_impedance.csv",
                    np.array([[freq, z_mag, theta_deg, z_real, z_imag]]),
                    "frequency_Hz, Z_mag_ohm, theta_deg, Z_real_ohm, Z_imag_ohm")
        return name, run

    elif choice == 7:
        # Single Point Capacitance
        freq = safe_float_input("Enter frequency (Hz) [default 1000]: ", 1000)
        mode = prompt_choice("Measurement mode", ["CPD", "CSRS", "CPRP"], "CPD")
        ac_level = safe_float_input("Enter AC level (V) [default 1.0]: ", 1.0)
        apply_bias = prompt_bool("Apply DC bias?", False)
        bias_v = safe_float_input("Enter DC bias (V) [default 0]: ", 0) if apply_bias else 0.0
        calc_eps = prompt_bool("Calculate permittivity?", False)
        thickness_nm = 0.0
        diameter_um = 0.0
        if calc_eps:
            thickness_nm = safe_float_input("Enter dielectric thickness (nm) [default 10]: ", 10)
            diameter_um = safe_float_input("Enter electrode diameter (µm) [default 75]: ", 75)

        def run():
            with InstrumentSession(
                lambda: lcr.setup(measurement_function=mode),
                lcr.disconnect_e4980a) as inst:
                lcr.set_ac_level(inst, ac_level)
                lcr.configure_dc_bias(inst, apply_bias, bias_v)
                cap, secondary = lcr.measure_capacitance(inst, freq, mode)
                sec_label = "D" if mode == "CPD" else "Secondary"
                result_str = f"  Cap={cap:.6e} F  {sec_label}={secondary:.6e}"
                if calc_eps:
                    eps_r = lcr.compute_eps_r(cap, thickness_nm, diameter_um)
                    result_str += f"  εr={eps_r:.2f}"
                    file_management.save_csv("lcr_single_capacitance.csv",
                        np.array([[freq, cap, secondary, eps_r]]),
                        f"frequency_Hz, capacitance_F, {sec_label}, eps_r")
                else:
                    file_management.save_csv("lcr_single_capacitance.csv",
                        np.array([[freq, cap, secondary]]),
                        f"frequency_Hz, capacitance_F, {sec_label}")
                print(result_str)
        return name, run

    elif choice == 8:
        # Quality Factor (Q) vs Frequency
        freq_start = safe_float_input("Enter start frequency (Hz) [default 20]: ", 20)
        freq_stop = safe_float_input("Enter stop frequency (Hz) [default 2e6]: ", 2e6)
        num_points = safe_int_input("Enter number of points [default 201]: ", 201)
        sweep_type = prompt_choice("Sweep type", ["LIN", "LOG"], "LOG")
        mode = prompt_choice("Measurement mode", ["CPQ", "LPQ"], "CPQ")
        ac_level = safe_float_input("Enter AC level (V) [default 1.0]: ", 1.0)

        def run():
            with InstrumentSession(
                lambda: lcr.setup(measurement_function=mode),
                lcr.disconnect_e4980a) as inst:
                lcr.set_ac_level(inst, ac_level)
                freq_axis, primary_vals, q_vals = lcr.measure_quality_factor_vs_frequency(
                    inst, freq_start, freq_stop, num_points, sweep_type, mode)
                primary_label = "Cp (F)" if mode == "CPQ" else "Lp (H)"
                file_management.save_image(
                    f"Quality Factor vs Frequency ({mode}) (LCR)", "Frequency (Hz)",
                    freq_axis, "Q (quality factor)", q_vals,
                    APPLY_DC_BIAS=False, DC_BIAS_V=0.0)
                file_management.save_csv(
                    f"lcr_q_vs_freq_{mode.lower()}.csv",
                    np.column_stack([freq_axis, primary_vals, q_vals]),
                    f"frequency_Hz, {primary_label.replace(' ', '_')}, Q")
        return name, run

    elif choice == 9:
        # R-X vs Frequency
        freq_start = safe_float_input("Enter start frequency (Hz) [default 20]: ", 20)
        freq_stop = safe_float_input("Enter stop frequency (Hz) [default 2e6]: ", 2e6)
        num_points = safe_int_input("Enter number of points [default 201]: ", 201)
        sweep_type = prompt_choice("Sweep type", ["LIN", "LOG"], "LOG")
        ac_level = safe_float_input("Enter AC level (V) [default 1.0]: ", 1.0)

        def run():
            with InstrumentSession(
                lambda: lcr.setup(measurement_function="RX"),
                lcr.disconnect_e4980a) as inst:
                lcr.set_ac_level(inst, ac_level)
                freq_axis, r_vals, x_vals = lcr.measure_rx_vs_frequency(
                    inst, freq_start, freq_stop, num_points, sweep_type)
                file_management.save_image(
                    "Resistance vs Frequency (LCR)", "Frequency (Hz)",
                    freq_axis, "R (Ohm)", r_vals,
                    APPLY_DC_BIAS=False, DC_BIAS_V=0.0)
                file_management.save_image(
                    "Reactance vs Frequency (LCR)", "Frequency (Hz)",
                    freq_axis, "X (Ohm)", x_vals,
                    APPLY_DC_BIAS=False, DC_BIAS_V=0.0)
                file_management.save_csv("lcr_rx_vs_freq.csv",
                    np.column_stack([freq_axis, r_vals, x_vals]),
                    "frequency_Hz, R_ohm, X_ohm")
        return name, run

    elif choice == 10:
        # G-B vs Frequency
        freq_start = safe_float_input("Enter start frequency (Hz) [default 20]: ", 20)
        freq_stop = safe_float_input("Enter stop frequency (Hz) [default 2e6]: ", 2e6)
        num_points = safe_int_input("Enter number of points [default 201]: ", 201)
        sweep_type = prompt_choice("Sweep type", ["LIN", "LOG"], "LOG")
        ac_level = safe_float_input("Enter AC level (V) [default 1.0]: ", 1.0)

        def run():
            with InstrumentSession(
                lambda: lcr.setup(measurement_function="GB"),
                lcr.disconnect_e4980a) as inst:
                lcr.set_ac_level(inst, ac_level)
                freq_axis, g_vals, b_vals = lcr.measure_gb_vs_frequency(
                    inst, freq_start, freq_stop, num_points, sweep_type)
                file_management.save_image(
                    "Conductance vs Frequency (LCR)", "Frequency (Hz)",
                    freq_axis, "G (S)", g_vals,
                    APPLY_DC_BIAS=False, DC_BIAS_V=0.0)
                file_management.save_image(
                    "Susceptance vs Frequency (LCR)", "Frequency (Hz)",
                    freq_axis, "B (S)", b_vals,
                    APPLY_DC_BIAS=False, DC_BIAS_V=0.0)
                file_management.save_csv("lcr_gb_vs_freq.csv",
                    np.column_stack([freq_axis, g_vals, b_vals]),
                    "frequency_Hz, G_S, B_S")
        return name, run

    elif choice == 11:
        # Open/Short Correction — not supported in cryo mode
        print("  Open/Short Correction is not supported in cryo sweep mode.")
        print("  Run it from the normal LCR menu before starting a cryo sweep.")
        return name, None

    else:
        log.warning("LCR choice %d not supported in cryo mode", choice)
        return name, None


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
                file_management.save_csv("cpd_data.csv", csv_data, "frequency_Hz, Cp_F, D")
                print("Measurement complete. Data saved to output folder.")
            
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
                file_management.save_csv("cpd_data.csv", csv_data, "frequency_Hz, Cp_F, D")
                print("Measurement complete. Data saved to output folder.")
            
        elif choice == 5:
            # C–V Butterfly Cycle: Cp vs Voltage
            _pia_cv_butterfly(plot_mode="cp")

        elif choice == 6:
            # C–V Butterfly Cycle: Permittivity vs Voltage
            _pia_cv_butterfly(plot_mode="eps")

        elif choice == 7:
            # εr vs frequency measurement
            print("Permittivity vs Frequency Measurement")

            try:
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
                    freq_start = safe_float_input("Enter start frequency (Hz) [default 1000]: ", 1000)
                    freq_stop = safe_float_input("Enter stop frequency (Hz) [default 1e6]: ", 1e6)
                    freq_points = safe_int_input("Enter number of frequency points [default 201]: ", 201)
                    bias_voltage = safe_float_input("Enter DC bias voltage (V) [default 0]: ", 0)

                    pia.current_parameters["freq_start_hz"] = freq_start
                    pia.current_parameters["freq_stop_hz"] = freq_stop
                    pia.current_parameters["num_points"] = freq_points
                    pia.current_parameters["dc_bias_v"] = bias_voltage
                    pia.current_parameters["apply_dc_bias"] = (bias_voltage != 0)

                thickness_nm = safe_float_input("Enter HZO thickness (nm) [default 10.0]: ", 10.0)
                diam_um = safe_float_input("Enter electrode diameter (µm) [default 75.0]: ", 75.0)

                with InstrumentSession(pia.setup, _pia_safe_close) as inst:
                    pia.initialize_4294a_for_cpd(inst)

                    print("Running εr vs frequency sweep...")
                    freq_axis, cp_f = pia.measure_eps_vs_freq(inst, freq_start, freq_stop, freq_points, bias_voltage)
                    eps_f = pia.compute_eps_r(cp_f, thickness_nm, diam_um)

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
            except Exception as e:
                log.error("PIA εr vs frequency measurement failed: %s", e)
                print(f"Error during measurement: {e}")

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
                    file_management.save_plot("transistor_output_chars.png")
                    plt.close()

                    print("Measurement complete. Data saved to output folder.")
            except Exception as e:
                log.error("PSPA output chars failed: %s", e)
                print(f"Error during measurement: {e}")

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
                    file_management.save_plot("transistor_transfer_chars.png", fig)
                    plt.close()

                    print("Measurement complete. Data saved to output folder.")
            except Exception as e:
                log.error("PSPA transfer chars failed: %s", e)
                print(f"Error during measurement: {e}")

        elif choice == 3:
            # PSPA I-V Curve (Unidirectional)
            print("PSPA: I-V Curve (Unidirectional)")
            
            v_start = safe_float_input("Enter start voltage (V) [default 0]: ", 0)
            v_stop = safe_float_input("Enter stop voltage (V) [default 5]: ", 5)
            v_step = safe_float_input("Enter voltage step (V) [default 0.1]: ", 0.1)
            compliance = safe_float_input("Enter current compliance (A) [default 0.1]: ", 0.1)
            channel = safe_int_input("Enter channel number [default 1]: ", 1)
            
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
            
            v_max = safe_float_input("Enter maximum voltage magnitude (V) [default 5]: ", 5)
            v_step = safe_float_input("Enter voltage step (V) [default 0.1]: ", 0.1)
            compliance = safe_float_input("Enter current compliance (A) [default 0.1]: ", 0.1)
            channel = safe_int_input("Enter channel number [default 1]: ", 1)
            
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
            
            v_base = safe_float_input("Enter base voltage (V) [default 0]: ", 0)
            v_pulse = safe_float_input("Enter pulse voltage (V) [default 5]: ", 5)
            pulse_width = safe_float_input("Enter pulse width (µs) [default 100]: ", 100) * 1e-6
            pulse_period = safe_float_input("Enter pulse period (ms) [default 10]: ", 10) * 1e-3
            num_pulses = safe_int_input("Enter number of pulses [default 10]: ", 10)
            compliance = safe_float_input("Enter current compliance (A) [default 0.1]: ", 0.1)
            channel = safe_int_input("Enter channel number [default 1]: ", 1)
            
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
                    file_management.save_plot("pulsed_iv.png", fig)
                    plt.close()

                    print("Measurement complete. Data saved to output folder.")
            except Exception as e:
                log.error("PSPA pulsed IV failed: %s", e)
                print(f"Error during measurement: {e}")

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

            v_start = safe_float_input("Enter start voltage (V) [default -1]: ", -1)
            v_stop = safe_float_input("Enter stop voltage (V) [default 1]: ", 1)
            v_step = safe_float_input("Enter voltage step (V) [default 0.02]: ", 0.02)
            compliance = safe_float_input("Enter current compliance (A) [default 0.1]: ", 0.1)
            anode_ch = safe_int_input("Enter anode channel [default 1]: ", 1)
            cathode_ch = safe_int_input("Enter cathode channel [default 2]: ", 2)

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

            vgs_start = safe_float_input("Enter start Vgs (V) [default 0]: ", 0)
            vgs_stop = safe_float_input("Enter stop Vgs (V) [default 5]: ", 5)
            vgs_step = safe_float_input("Enter Vgs step (V) [default 0.1]: ", 0.1)
            compliance = safe_float_input("Enter current compliance (A) [default 1e-3]: ", 1e-3)
            gate_ch = safe_int_input("Enter gate channel [default 2]: ", 2)
            source_ch = safe_int_input("Enter source channel [default 3]: ", 3)

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

            v_start = safe_float_input("Enter start voltage (V) [default 0]: ", 0)
            v_stop = safe_float_input("Enter stop voltage (V) [default 40]: ", 40)
            v_step = safe_float_input("Enter voltage step (V) [default 0.5]: ", 0.5)
            compliance = safe_float_input("Enter current compliance (A) [default 1e-3]: ", 1e-3)
            threshold = safe_float_input("Enter breakdown current threshold (A) [default 1e-4]: ", 1e-4)
            channel = safe_int_input("Enter channel number [default 1]: ", 1)

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

            i_start = safe_float_input("Enter start current (A) [default 0]: ", 0)
            i_stop = safe_float_input("Enter stop current (A) [default 1e-3]: ", 1e-3)
            i_step = safe_float_input("Enter current step (A) [default 1e-4]: ", 1e-4)
            compliance = safe_float_input("Enter voltage compliance (V) [default 10.0]: ", 10.0)
            channel = safe_int_input("Enter channel number [default 1]: ", 1)

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

            freq_start = safe_float_input("Enter start frequency (Hz) [default 20]: ", 20)
            freq_stop = safe_float_input("Enter stop frequency (Hz) [default 2e6]: ", 2e6)
            num_points = safe_int_input("Enter number of points [default 201]: ", 201)
            sweep_type = prompt_choice("Sweep type", ["LIN", "LOG"], "LOG")
            ac_level = safe_float_input("Enter AC level (V) [default 1.0]: ", 1.0)
            apply_bias = prompt_bool("Apply DC bias?", False)
            bias_v = 0.0
            if apply_bias:
                bias_v = safe_float_input("Enter DC bias (V) [default 0]: ", 0)

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

            freq_start = safe_float_input("Enter start frequency (Hz) [default 20]: ", 20)
            freq_stop = safe_float_input("Enter stop frequency (Hz) [default 2e6]: ", 2e6)
            num_points = safe_int_input("Enter number of points [default 201]: ", 201)
            sweep_type = prompt_choice("Sweep type", ["LIN", "LOG"], "LOG")
            mode = prompt_choice("Measurement mode", ["CPD", "CSRS", "CPRP"], "CPD")
            ac_level = safe_float_input("Enter AC level (V) [default 1.0]: ", 1.0)
            apply_bias = prompt_bool("Apply DC bias?", False)
            bias_v = 0.0
            if apply_bias:
                bias_v = safe_float_input("Enter DC bias (V) [default 0]: ", 0)

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

            freq_start = safe_float_input("Enter start frequency (Hz) [default 20]: ", 20)
            freq_stop = safe_float_input("Enter stop frequency (Hz) [default 2e6]: ", 2e6)
            num_points = safe_int_input("Enter number of points [default 201]: ", 201)
            sweep_type = prompt_choice("Sweep type", ["LIN", "LOG"], "LOG")
            mode = prompt_choice("Measurement mode", ["LPQ", "LSD"], "LPQ")
            ac_level = safe_float_input("Enter AC level (V) [default 1.0]: ", 1.0)

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

            freq = safe_float_input("Enter measurement frequency (Hz) [default 1000]: ", 1000)
            v_min = safe_float_input("Enter minimum voltage (V) [default -5]: ", -5)
            v_max = safe_float_input("Enter maximum voltage (V) [default 5]: ", 5)
            num_points = safe_int_input("Enter number of points [default 201]: ", 201)
            ac_level = safe_float_input("Enter AC level (V) [default 0.1]: ", 0.1)

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

            freq = safe_float_input("Enter measurement frequency (Hz) [default 1000]: ", 1000)
            v_min = safe_float_input("Enter minimum voltage (V) [default -5]: ", -5)
            v_max = safe_float_input("Enter maximum voltage (V) [default 5]: ", 5)
            num_points = safe_int_input("Enter number of points per direction [default 201]: ", 201)
            num_cycles = safe_int_input("Enter number of cycles [default 1]: ", 1)
            ac_level = safe_float_input("Enter AC level (V) [default 0.1]: ", 0.1)
            calc_eps = prompt_bool("Calculate permittivity?", False)
            thickness_nm = 0.0
            diameter_um = 0.0
            if calc_eps:
                thickness_nm = safe_float_input("Enter dielectric thickness (nm) [default 10]: ", 10)
                diameter_um = safe_float_input("Enter electrode diameter (µm) [default 75]: ", 75)

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

            freq = safe_float_input("Enter frequency (Hz) [default 1000]: ", 1000)
            ac_level = safe_float_input("Enter AC level (V) [default 1.0]: ", 1.0)
            apply_bias = prompt_bool("Apply DC bias?", False)
            bias_v = 0.0
            if apply_bias:
                bias_v = safe_float_input("Enter DC bias (V) [default 0]: ", 0)

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

            freq = safe_float_input("Enter frequency (Hz) [default 1000]: ", 1000)
            mode = prompt_choice("Measurement mode", ["CPD", "CSRS", "CPRP"], "CPD")
            ac_level = safe_float_input("Enter AC level (V) [default 1.0]: ", 1.0)
            apply_bias = prompt_bool("Apply DC bias?", False)
            bias_v = 0.0
            if apply_bias:
                bias_v = safe_float_input("Enter DC bias (V) [default 0]: ", 0)
            calc_eps = prompt_bool("Calculate permittivity?", False)
            thickness_nm = 0.0
            diameter_um = 0.0
            if calc_eps:
                thickness_nm = safe_float_input("Enter dielectric thickness (nm) [default 10]: ", 10)
                diameter_um = safe_float_input("Enter electrode diameter (µm) [default 75]: ", 75)

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

            freq_start = safe_float_input("Enter start frequency (Hz) [default 20]: ", 20)
            freq_stop = safe_float_input("Enter stop frequency (Hz) [default 2e6]: ", 2e6)
            num_points = safe_int_input("Enter number of points [default 201]: ", 201)
            sweep_type = prompt_choice("Sweep type", ["LIN", "LOG"], "LOG")
            mode = prompt_choice("Measurement mode", ["CPQ", "LPQ"], "CPQ")
            ac_level = safe_float_input("Enter AC level (V) [default 1.0]: ", 1.0)

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

            freq_start = safe_float_input("Enter start frequency (Hz) [default 20]: ", 20)
            freq_stop = safe_float_input("Enter stop frequency (Hz) [default 2e6]: ", 2e6)
            num_points = safe_int_input("Enter number of points [default 201]: ", 201)
            sweep_type = prompt_choice("Sweep type", ["LIN", "LOG"], "LOG")
            ac_level = safe_float_input("Enter AC level (V) [default 1.0]: ", 1.0)

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

            freq_start = safe_float_input("Enter start frequency (Hz) [default 20]: ", 20)
            freq_stop = safe_float_input("Enter stop frequency (Hz) [default 2e6]: ", 2e6)
            num_points = safe_int_input("Enter number of points [default 201]: ", 201)
            sweep_type = prompt_choice("Sweep type", ["LIN", "LOG"], "LOG")
            ac_level = safe_float_input("Enter AC level (V) [default 1.0]: ", 1.0)

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

            corr_choice = safe_int_input("Select correction type [default 3]: ", 3)

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
