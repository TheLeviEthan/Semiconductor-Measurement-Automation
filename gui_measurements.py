"""
Filename: gui_measurements.py
Author: Ethan Ruddell
Date: 2026-2-27
Description: GUI-specific measurement execution functions.
These take a params dict from the GUI and call instrument functions directly,
without any interactive input() prompts.
"""

import os
import sys
import logging
import numpy as np
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for thread safety
import matplotlib.pyplot as plt

# Add directories to path
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'utility'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'measurement functions'))

import pia
import pspa
import lcr
import file_management
from gpib_utils import InstrumentSession

log = logging.getLogger(__name__)


# =============================
# Helper
# =============================

def _pia_safe_close(inst):
    """Turn off DC bias and close PIA instrument."""
    try:
        inst.write("DCO OFF")
    except Exception:
        pass
    inst.close()


# =============================
# PIA GUI Measurements
# =============================

def execute_pia_gui(choice, params):
    """Execute a PIA measurement using GUI-provided params (no input() calls)."""

    if choice in (1, 2):
        # Impedance Magnitude (1) or Phase (2) vs Frequency
        freq_start = float(params.get("freq_start", 1000))
        freq_stop = float(params.get("freq_stop", 1e6))
        num_points = int(float(params.get("num_points", 201)))
        dc_bias_v = float(params.get("dc_bias_v", 0))
        apply_dc_bias = dc_bias_v != 0

        with InstrumentSession(pia.setup, _pia_safe_close) as inst:
            pia.initialize_4294a_for_impedance(inst)
            pia.configure_dc_bias(inst, apply_dc_bias, dc_bias_v)
            freq_axis, z_mag, theta_deg = pia.measure_impedance_vs_freq(
                inst, freq_start, freq_stop, num_points)

            y_label = "|Z| (Ohm)" if choice == 1 else "Phase (Degrees)"
            y_data = z_mag if choice == 1 else theta_deg
            title = "Impedance Magnitude vs Frequency" if choice == 1 else "Impedance Phase vs Frequency"

            file_management.save_image(
                title, "Frequency (Hz)", freq_axis, y_label, y_data,
                APPLY_DC_BIAS=apply_dc_bias, DC_BIAS_V=dc_bias_v)

            z_real = z_mag * np.cos(np.deg2rad(theta_deg))
            z_imag = z_mag * np.sin(np.deg2rad(theta_deg))
            file_management.save_csv("impedance_data.csv",
                np.column_stack([freq_axis, z_mag, theta_deg, z_real, z_imag]),
                "frequency_Hz, Z_mag_ohm, theta_deg, Z_real_ohm, Z_imag_ohm")

    elif choice in (3, 4):
        # Capacitance (3) or Tan Loss (4) vs Frequency
        freq_start = float(params.get("freq_start", 1000))
        freq_stop = float(params.get("freq_stop", 1e6))
        num_points = int(float(params.get("num_points", 201)))
        dc_bias_v = float(params.get("dc_bias_v", 0))
        apply_dc_bias = dc_bias_v != 0

        with InstrumentSession(pia.setup, _pia_safe_close) as inst:
            pia.initialize_4294a_for_cpd(inst)
            pia.configure_dc_bias(inst, apply_dc_bias, dc_bias_v)
            freq_axis, cp_vals, d_vals = pia.measure_cpd_vs_freq(
                inst, freq_start, freq_stop, num_points)

            y_label = "Cp (F)" if choice == 3 else "D (loss factor)"
            y_data = cp_vals if choice == 3 else d_vals
            title = "Capacitance vs Frequency" if choice == 3 else "Tan Loss vs Frequency"

            file_management.save_image(
                title, "Frequency (Hz)", freq_axis, y_label, y_data,
                APPLY_DC_BIAS=apply_dc_bias, DC_BIAS_V=dc_bias_v)
            file_management.save_csv("cpd_data.csv",
                np.column_stack([freq_axis, cp_vals, d_vals]),
                "frequency_Hz, Cp_F, D")

    elif choice in (5, 6):
        # C-V Butterfly: Cp (5) or Permittivity (6) vs Voltage
        plot_mode = "cp" if choice == 5 else "eps"
        freq_cv = float(params.get("freq_cv", 25000))
        v_min = float(params.get("v_min", 0))
        v_max = float(params.get("v_max", 5))
        n_points = int(float(params.get("n_points", 401)))
        n_cycles = int(float(params.get("n_cycles", 1)))
        thickness_nm = float(params.get("thickness_nm", 10.0))
        diam_um = float(params.get("diam_um", 75.0))

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

    elif choice == 7:
        # Permittivity vs Frequency
        freq_start = float(params.get("freq_start", 1000))
        freq_stop = float(params.get("freq_stop", 1e6))
        num_points = int(float(params.get("num_points", 201)))
        dc_bias_v = float(params.get("dc_bias_v", 0))
        thickness_nm = float(params.get("thickness_nm", 10.0))
        diam_um = float(params.get("diam_um", 75.0))

        with InstrumentSession(pia.setup, _pia_safe_close) as inst:
            pia.initialize_4294a_for_cpd(inst)
            freq_axis, cp_f = pia.measure_eps_vs_freq(
                inst, freq_start, freq_stop, num_points, dc_bias_v)
            eps_f = pia.compute_eps_r(cp_f, thickness_nm, diam_um)

            file_management.save_csv("eps_vs_freq.csv",
                np.column_stack([freq_axis, cp_f, eps_f]),
                "frequency_Hz, Cp_F, eps_r")
            file_management.save_image("Permittivity vs Frequency",
                "Frequency (Hz)", freq_axis, "Dielectric constant εr", eps_f,
                APPLY_DC_BIAS=(dc_bias_v != 0), DC_BIAS_V=dc_bias_v)

    elif choice == 8:
        # R-X vs Frequency
        freq_start = float(params.get("freq_start", 1000))
        freq_stop = float(params.get("freq_stop", 1e6))
        num_points = int(float(params.get("num_points", 201)))
        dc_bias_v = float(params.get("dc_bias_v", 0))
        apply_dc_bias = dc_bias_v != 0

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

    elif choice == 9:
        # G-B vs Frequency
        freq_start = float(params.get("freq_start", 1000))
        freq_stop = float(params.get("freq_stop", 1e6))
        num_points = int(float(params.get("num_points", 201)))
        dc_bias_v = float(params.get("dc_bias_v", 0))
        apply_dc_bias = dc_bias_v != 0

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

    elif choice == 10:
        # Y-θ vs Frequency
        freq_start = float(params.get("freq_start", 1000))
        freq_stop = float(params.get("freq_stop", 1e6))
        num_points = int(float(params.get("num_points", 201)))
        dc_bias_v = float(params.get("dc_bias_v", 0))
        apply_dc_bias = dc_bias_v != 0

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

    else:
        raise ValueError(f"Unknown PIA measurement index: {choice}")


# =============================
# PSPA GUI Measurements
# =============================

def execute_pspa_gui(choice, params):
    """Execute a PSPA measurement using GUI-provided params (no input() calls)."""

    if choice == 1:
        # Transistor Output Characteristics
        vds_start = float(params.get("vds_start", 0))
        vds_stop = float(params.get("vds_stop", 5))
        vds_step = float(params.get("vds_step", 0.1))
        vgs_start = float(params.get("vgs_start", 0))
        vgs_stop = float(params.get("vgs_stop", 3))
        vgs_step = float(params.get("vgs_step", 0.5))
        compliance = float(params.get("compliance", 0.1))
        drain_ch = int(float(params.get("drain_ch", 1)))
        gate_ch = int(float(params.get("gate_ch", 2)))
        source_ch = int(float(params.get("source_ch", 3)))

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

    elif choice == 2:
        # Transistor Transfer Characteristics
        vgs_start = float(params.get("vgs_start", -1))
        vgs_stop = float(params.get("vgs_stop", 3))
        vgs_step = float(params.get("vgs_step", 0.05))
        vds_constant = float(params.get("vds_constant", 5))
        compliance = float(params.get("compliance", 0.1))
        drain_ch = int(float(params.get("drain_ch", 1)))
        gate_ch = int(float(params.get("gate_ch", 2)))
        source_ch = int(float(params.get("source_ch", 3)))

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

    elif choice == 3:
        # I-V Curve (Unidirectional)
        v_start = float(params.get("v_start", 0))
        v_stop = float(params.get("v_stop", 5))
        v_step = float(params.get("v_step", 0.1))
        compliance = float(params.get("compliance", 0.1))
        channel = int(float(params.get("channel", 1)))

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

    elif choice == 4:
        # I-V Curve (Bidirectional)
        v_max = float(params.get("v_max", 5))
        v_step = float(params.get("v_step", 0.1))
        compliance = float(params.get("compliance", 0.1))
        channel = int(float(params.get("channel", 1)))

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

    elif choice == 5:
        # Pulsed I-V (Single Device)
        v_base = float(params.get("v_base", 0))
        v_pulse = float(params.get("v_pulse", 5))
        pulse_width = float(params.get("pulse_width_us", 100)) * 1e-6
        pulse_period = float(params.get("pulse_period_ms", 10)) * 1e-3
        num_pulses = int(float(params.get("num_pulses", 10)))
        compliance = float(params.get("compliance", 0.1))
        channel = int(float(params.get("channel", 1)))

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

    elif choice == 6:
        # Pulsed Transistor
        vds_base = float(params.get("vds_base", 0))
        vds_pulse = float(params.get("vds_pulse", 5))
        vgs_base = float(params.get("vgs_base", 0))
        vgs_pulse = float(params.get("vgs_pulse", 3))
        pulse_width = float(params.get("pulse_width_us", 100)) * 1e-6
        pulse_period = float(params.get("pulse_period_ms", 10)) * 1e-3
        num_pulses = int(float(params.get("num_pulses", 10)))
        compliance = float(params.get("compliance", 0.1))
        drain_ch = int(float(params.get("drain_ch", 1)))
        gate_ch = int(float(params.get("gate_ch", 2)))
        source_ch = int(float(params.get("source_ch", 3)))

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

    elif choice == 7:
        # Diode I-V Characteristics
        v_start = float(params.get("v_start", -1))
        v_stop = float(params.get("v_stop", 1))
        v_step = float(params.get("v_step", 0.02))
        compliance = float(params.get("compliance", 0.1))
        anode_ch = int(float(params.get("anode_ch", 1)))
        cathode_ch = int(float(params.get("cathode_ch", 2)))

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

    elif choice == 8:
        # Gate Leakage Current
        vgs_start = float(params.get("vgs_start", 0))
        vgs_stop = float(params.get("vgs_stop", 5))
        vgs_step = float(params.get("vgs_step", 0.1))
        compliance = float(params.get("compliance", 1e-3))
        gate_ch = int(float(params.get("gate_ch", 2)))
        source_ch = int(float(params.get("source_ch", 3)))

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

    elif choice == 9:
        # Breakdown Voltage
        v_start = float(params.get("v_start", 0))
        v_stop = float(params.get("v_stop", 40))
        v_step = float(params.get("v_step", 0.5))
        compliance = float(params.get("compliance", 1e-3))
        threshold = float(params.get("threshold_a", 1e-4))
        channel = int(float(params.get("channel", 1)))

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

    elif choice == 10:
        # Resistance Measurement
        i_start = float(params.get("i_start", 0))
        i_stop = float(params.get("i_stop", 1e-3))
        i_step = float(params.get("i_step", 1e-4))
        compliance = float(params.get("compliance", 10.0))
        channel = int(float(params.get("channel", 1)))

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

    else:
        raise ValueError(f"Unknown PSPA measurement index: {choice}")


# =============================
# LCR GUI Measurements
# =============================

def execute_lcr_gui(choice, params):
    """Execute an LCR measurement using GUI-provided params (no input() calls)."""

    if choice == 1:
        # Impedance vs Frequency
        freq_start = float(params.get("freq_start", 20))
        freq_stop = float(params.get("freq_stop", 2e6))
        num_points = int(float(params.get("num_points", 201)))
        sweep_type = str(params.get("sweep_type", "LOG")).upper()
        ac_level = float(params.get("ac_level", 1.0))
        dc_bias_v = float(params.get("dc_bias_v", 0))
        apply_bias = dc_bias_v != 0

        with InstrumentSession(
            lambda: lcr.setup(measurement_function="ZTD"),
            lcr.disconnect_e4980a) as inst:
            lcr.set_ac_level(inst, ac_level)
            lcr.configure_dc_bias(inst, apply_bias, dc_bias_v)
            freq_axis, z_mag, theta_deg = lcr.measure_impedance_vs_frequency(
                inst, freq_start, freq_stop, num_points, sweep_type)

            file_management.save_image(
                "Impedance Magnitude vs Frequency (LCR)", "Frequency (Hz)",
                freq_axis, "|Z| (Ohm)", z_mag,
                APPLY_DC_BIAS=apply_bias, DC_BIAS_V=dc_bias_v)
            file_management.save_image(
                "Impedance Phase vs Frequency (LCR)", "Frequency (Hz)",
                freq_axis, "Phase (Degrees)", theta_deg,
                APPLY_DC_BIAS=apply_bias, DC_BIAS_V=dc_bias_v)
            z_real, z_imag = lcr.compute_impedance_components(z_mag, theta_deg)
            file_management.save_csv("lcr_impedance_vs_freq.csv",
                np.column_stack([freq_axis, z_mag, theta_deg, z_real, z_imag]),
                "frequency_Hz, Z_mag_ohm, theta_deg, Z_real_ohm, Z_imag_ohm")

    elif choice == 2:
        # Capacitance vs Frequency
        freq_start = float(params.get("freq_start", 20))
        freq_stop = float(params.get("freq_stop", 2e6))
        num_points = int(float(params.get("num_points", 201)))
        sweep_type = str(params.get("sweep_type", "LOG")).upper()
        mode = str(params.get("mode", "CPD")).upper()
        ac_level = float(params.get("ac_level", 1.0))
        dc_bias_v = float(params.get("dc_bias_v", 0))
        apply_bias = dc_bias_v != 0

        with InstrumentSession(
            lambda: lcr.setup(measurement_function=mode),
            lcr.disconnect_e4980a) as inst:
            lcr.set_ac_level(inst, ac_level)
            lcr.configure_dc_bias(inst, apply_bias, dc_bias_v)
            freq_axis, cap_vals, sec_vals = lcr.measure_capacitance_vs_frequency(
                inst, freq_start, freq_stop, num_points, sweep_type, mode)
            sec_label = "D (loss factor)" if mode == "CPD" else "Secondary"
            file_management.save_image(
                f"Capacitance vs Frequency ({mode}) (LCR)", "Frequency (Hz)",
                freq_axis, "Capacitance (F)", cap_vals,
                APPLY_DC_BIAS=apply_bias, DC_BIAS_V=dc_bias_v)
            file_management.save_image(
                f"{sec_label} vs Frequency (LCR)", "Frequency (Hz)",
                freq_axis, sec_label, sec_vals,
                APPLY_DC_BIAS=apply_bias, DC_BIAS_V=dc_bias_v)
            file_management.save_csv(
                f"lcr_capacitance_vs_freq_{mode.lower()}.csv",
                np.column_stack([freq_axis, cap_vals, sec_vals]),
                f"frequency_Hz, capacitance_F, {sec_label.replace(' ', '_')}")

    elif choice == 3:
        # Inductance vs Frequency
        freq_start = float(params.get("freq_start", 20))
        freq_stop = float(params.get("freq_stop", 2e6))
        num_points = int(float(params.get("num_points", 201)))
        sweep_type = str(params.get("sweep_type", "LOG")).upper()
        mode = str(params.get("mode", "LPQ")).upper()
        ac_level = float(params.get("ac_level", 1.0))

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

    elif choice == 4:
        # C-V Sweep (Single)
        freq = float(params.get("freq", 1000))
        v_min = float(params.get("v_min", -5))
        v_max = float(params.get("v_max", 5))
        num_points = int(float(params.get("num_points", 201)))
        ac_level = float(params.get("ac_level", 0.1))

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

    elif choice == 5:
        # C-V Butterfly Cycles
        freq = float(params.get("freq", 1000))
        v_min = float(params.get("v_min", -5))
        v_max = float(params.get("v_max", 5))
        num_points = int(float(params.get("num_points", 201)))
        num_cycles = int(float(params.get("num_cycles", 1)))
        ac_level = float(params.get("ac_level", 0.1))
        thickness_nm = float(params.get("thickness_nm", 0))
        diameter_um = float(params.get("diameter_um", 0))
        calc_eps = thickness_nm > 0 and diameter_um > 0

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

    elif choice == 6:
        # Single Point Impedance
        freq = float(params.get("freq", 1000))
        ac_level = float(params.get("ac_level", 1.0))
        dc_bias_v = float(params.get("dc_bias_v", 0))
        apply_bias = dc_bias_v != 0

        with InstrumentSession(
            lambda: lcr.setup(measurement_function="ZTD"),
            lcr.disconnect_e4980a) as inst:
            lcr.set_ac_level(inst, ac_level)
            lcr.configure_dc_bias(inst, apply_bias, dc_bias_v)
            z_mag, theta_deg = lcr.measure_impedance(inst, freq)
            z_real, z_imag = lcr.compute_impedance_components(z_mag, theta_deg)
            log.info("|Z|=%.6e Ω  θ=%.2f°  Z_real=%.6e Ω  Z_imag=%.6e Ω",
                     z_mag, theta_deg, z_real, z_imag)
            file_management.save_csv("lcr_single_impedance.csv",
                np.array([[freq, z_mag, theta_deg, z_real, z_imag]]),
                "frequency_Hz, Z_mag_ohm, theta_deg, Z_real_ohm, Z_imag_ohm")

    elif choice == 7:
        # Single Point Capacitance
        freq = float(params.get("freq", 1000))
        mode = str(params.get("mode", "CPD")).upper()
        ac_level = float(params.get("ac_level", 1.0))
        dc_bias_v = float(params.get("dc_bias_v", 0))
        apply_bias = dc_bias_v != 0
        thickness_nm = float(params.get("thickness_nm", 0))
        diameter_um = float(params.get("diameter_um", 0))
        calc_eps = thickness_nm > 0 and diameter_um > 0

        with InstrumentSession(
            lambda: lcr.setup(measurement_function=mode),
            lcr.disconnect_e4980a) as inst:
            lcr.set_ac_level(inst, ac_level)
            lcr.configure_dc_bias(inst, apply_bias, dc_bias_v)
            cap, secondary = lcr.measure_capacitance(inst, freq, mode)
            sec_label = "D" if mode == "CPD" else "Secondary"
            if calc_eps:
                eps_r = lcr.compute_eps_r(cap, thickness_nm, diameter_um)
                file_management.save_csv("lcr_single_capacitance.csv",
                    np.array([[freq, cap, secondary, eps_r]]),
                    f"frequency_Hz, capacitance_F, {sec_label}, eps_r")
            else:
                file_management.save_csv("lcr_single_capacitance.csv",
                    np.array([[freq, cap, secondary]]),
                    f"frequency_Hz, capacitance_F, {sec_label}")

    elif choice == 8:
        # Quality Factor (Q) vs Frequency
        freq_start = float(params.get("freq_start", 20))
        freq_stop = float(params.get("freq_stop", 2e6))
        num_points = int(float(params.get("num_points", 201)))
        sweep_type = str(params.get("sweep_type", "LOG")).upper()
        mode = str(params.get("mode", "CPQ")).upper()
        ac_level = float(params.get("ac_level", 1.0))

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

    elif choice == 9:
        # R-X vs Frequency
        freq_start = float(params.get("freq_start", 20))
        freq_stop = float(params.get("freq_stop", 2e6))
        num_points = int(float(params.get("num_points", 201)))
        sweep_type = str(params.get("sweep_type", "LOG")).upper()
        ac_level = float(params.get("ac_level", 1.0))

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

    elif choice == 10:
        # G-B vs Frequency
        freq_start = float(params.get("freq_start", 20))
        freq_stop = float(params.get("freq_stop", 2e6))
        num_points = int(float(params.get("num_points", 201)))
        sweep_type = str(params.get("sweep_type", "LOG")).upper()
        ac_level = float(params.get("ac_level", 1.0))

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

    elif choice == 11:
        # Open/Short Correction - not supported from GUI without interactive prompts
        raise ValueError(
            "Open/Short Correction requires interactive prompts "
            "(connect/disconnect DUT). Please use CLI mode (--cli) for this measurement.")

    else:
        raise ValueError(f"Unknown LCR measurement index: {choice}")
