# AI Coding Agent Instructions for Semiconductor-Measurement-Automation

These instructions make AI agents immediately productive in this repo by explaining the architecture, workflows, and project-specific patterns.

## Big Picture
- Purpose: Automate measurements using an Agilent/Keysight 4294A (PIA) via GPIB using `pyvisa`, and generate plots/CSVs.
- Components:
  - `main.py`: CLI entrypoint. Selects measurement type and orchestrates instrument setup, sweeps, and saving results.
  - `pia.py`: Instrument constants and SCPI helpers for the 4294A. Provides initialization and sweep functions for impedance (|Z|–θ) and Cp–D.
  - `file_management.py`: Centralized output handling: ensures an `output/` folder, uniquifies filenames, and saves plots/CSVs.
  - Legacy scripts: [`impedence_theta.py`](Semiconductor-Measurement-Automation/impedence_theta.py), [`capacitance_vs_tan_loss.py`](Semiconductor-Measurement-Automation/capacitance_vs_tan_loss.py), [`dielectric_meas.py`](Semiconductor-Measurement-Automation/dielectric_meas.py) contain earlier standalone flows now being consolidated.
  - `measurement_functions.py`: Draft/placeholder for additional flows (C–V cycles, εr computations) – many symbols are not yet wired.

## How Things Flow
- CLI: `main.py` prompts selection, then uses `pia.py` to connect and configure the 4294A, runs a sweep, and hands data to `file_management.py` for saving.
- Instrument session: `pia.connect_4294a()` → `pia.initialize_4294a_for_impedance()` or `pia.initialize_4294a_for_cpd()` → configure DC bias via `pia.configure_dc_bias()` → set sweep (`pia.set_frequency_sweep()`) → run + block on completion (`pia.single_sweep_and_wait()`).
- Data acquisition: `pia.read_sweep_axis()` returns frequency; `pia.read_trace_main(trace)` returns main values for trace A/B (Cp, D, |Z|, θ).
- Output: `file_management.ensure_output_dir(file_management.output_dir)` → `file_management.save_image(...)`; a CSV saver exists but needs parameterization to accept arrays and headers.

## Conventions and Patterns
- Paths: Output root is `Semiconductor-Measurement-Automation/output/` resolved from `sys.argv[0]` in `file_management.py`.
- Uniquify filenames: Use `file_management.uniquify(path)` to avoid clobbering.
- Plotting: Use `matplotlib` with `semilogx` for frequency sweeps. Titles optionally include DC bias.
- DC Bias: `pia.APPLY_DC_BIAS` (bool) and `pia.DC_BIAS_V` (float). When saving plots, pass both so titles/files reflect the actual bias.
- Sweep types: Frequency sweeps are LOG (`SWPT LOG`), DC bias sweeps are LIN (`SWPT LIN`). Use `SWPP FREQ` vs `SWPP DCB` appropriately.
- Traces: For Cp–D mode, trace A = Cp, trace B = D. For Z–θ mode, trace A = |Z|, trace B = θ (degrees).

## Developer Workflows
- Run CLI:
  ```powershell
  & C:\Python313\python.exe "c:/Users/rudde/NRG Scripts Local/Semiconductor-Measurement-Automation/main.py"
  ```
- Quick impedance sweep (programmatic example):
  - Connect and init: `inst = pia.connect_4294a(); pia.initialize_4294a_for_impedance(inst)`
  - Bias: `pia.configure_dc_bias(inst, pia.APPLY_DC_BIAS, pia.DC_BIAS_V)`
  - Measure: `freq, z_mag, theta = pia.measure_impedance_vs_freq(inst)`
  - Save: `file_management.save_image("Impedance Magnitude vs Frequency", "Frequency (Hz)", freq, "|Z| (Ohm)", z_mag, APPLY_DC_BIAS=pia.APPLY_DC_BIAS, DC_BIAS_V=pia.DC_BIAS_V)`

## Project-Specific Gotchas
- Naming mismatches: The legacy files `impedence_theta.py` and `dielectric_meas.py` differ from the names mentioned in consolidation notes; ensure imports reference actual filenames.
- `main.py` uses `repeating` without initialization and calls `file_management.save_csv(path)` but the function signature currently takes no arguments. Update `save_csv` to accept data/headers and a target path.
- `file_management.save_image` sets the plot title before appending the DC bias string; pass `DC_BIAS_V` and update ordering so titles and filenames match.
- `measurement_functions.py` references undefined variables (`freq_for_cv`, `v_min`, etc.); treat it as a future integration area or move those constants into `pia.py`.

## Extending Measurements
- Add menu handlers in [`main.py`](Semiconductor-Measurement-Automation/main.py) for options 2 (Cp–D) and 3 (Dielectric C–V):
  - Cp–D: call `pia.initialize_4294a_for_cpd(inst)` then `pia.measure_cpd_vs_freq(inst)`; save Cp and D plots and CSV via `file_management`.
  - Dielectric C–V cycles: port `measure_single_cv_cycle` and `compute_eps_r` from [`dielectric_meas.py`](Semiconductor-Measurement-Automation/dielectric_meas.py) into a cohesive API (e.g., `pia_cv.py` or expand `pia.py`).

## External Dependencies
- Requires `pyvisa`, `numpy`, `matplotlib`, and a working VISA/GPIB setup. GPIB address defaults to `GPIB0::24::INSTR` in `pia.py`.

## Style and Structure
- Keep instrument SCPI in `pia.py`; keep plotting/CSV in `file_management.py`; keep orchestration in `main.py`.
- Prefer passing arrays and labels into `file_management` rather than duplicating plotting code in measurement modules.

Feedback: If any of the above is unclear or incomplete (e.g., exact DC bias handling or C–V cycle integration), tell me which parts to refine and I’ll iterate.