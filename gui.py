"""
Filename: gui.py
Author: Ethan Ruddell
Date: 2026-2-12
Description: Implements a simple GUI for semiconductor measurement automation
using tkinter. Provides an interface to configure and run measurements on
the Keysight 4294A, Agilent 4155C/4156C, and Keysight E4980A instruments.
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import sys
import os
import logging
import threading

# Add directories to path
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'utility'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'measurement functions'))

import file_management
import logging_config
from measurements_config import (
    PIA_MEASUREMENTS, PSPA_MEASUREMENTS, LCR_MEASUREMENTS,
    MEASUREMENT_PARAMS, get_measurements_list,
)

log = logging.getLogger(__name__)


class MeasurementGUI:
    """Main GUI class for semiconductor measurement automation."""

    PIA_MEASUREMENTS = PIA_MEASUREMENTS
    PSPA_MEASUREMENTS = PSPA_MEASUREMENTS
    LCR_MEASUREMENTS = LCR_MEASUREMENTS
    MEASUREMENT_PARAMS = MEASUREMENT_PARAMS

    def __init__(self, root, measurement_executor):
        """
        Initialize the GUI.
        
        Args:
            root: Tkinter root window
            measurement_executor: Callback function to execute measurements
        """
        self.root = root
        self.measurement_executor = measurement_executor
        self.root.title("NRG Semiconductor Measurement Automation")
        self.root.geometry("900x700")
        
        # Setup logging
        logging_config.setup(output_dir=file_management.default_output_dir)
        
        # Track cryo state and measurement queue
        self.cryo_enabled = False
        self.measurement_queue = []
        
        self.create_widgets()
        self.toggle_cryo_params()  # Hide cryo UI initially
        self.update_measurement_list()

    def create_widgets(self):
        """Create GUI widgets."""
        # --- Main Frame ---
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        # Configure rows to expand/contract proportionally
        main_frame.rowconfigure(4, weight=1)  # Measurement selection
        main_frame.rowconfigure(5, weight=2)  # Parameters (largest)
        main_frame.rowconfigure(7, weight=1)  # Status

        # --- Output Directory ---
        ttk.Label(main_frame, text="Output Directory:", font=("Arial", 10, "bold")).grid(
            row=0, column=0, sticky=tk.W, pady=(0, 5))
        output_frame = ttk.Frame(main_frame)
        output_frame.grid(row=0, column=1, sticky=(tk.W, tk.E), pady=(0, 5))
        output_frame.columnconfigure(0, weight=1)
        
        self.output_dir_var = tk.StringVar(value=file_management.default_output_dir)
        output_entry = ttk.Entry(output_frame, textvariable=self.output_dir_var, state='readonly')
        output_entry.grid(row=0, column=0, sticky=(tk.W, tk.E))
        
        browse_btn = ttk.Button(output_frame, text="Browse", command=self.browse_output_dir)
        browse_btn.grid(row=0, column=1, padx=(5, 0))

        # --- Cryo Integration ---
        self.cryo_var = tk.BooleanVar(value=False)
        cryo_check = ttk.Checkbutton(main_frame, text="Enable Cryogenic Temperature Sweep",
                                     variable=self.cryo_var, command=self.toggle_cryo_params)
        cryo_check.grid(row=1, column=0, columnspan=2, sticky=tk.W, pady=(10, 5))

        # --- Cryo Parameters Frame (initially hidden) ---
        cryo_params_frame = ttk.LabelFrame(main_frame, text="Cryogenic Parameters", padding="10")
        cryo_params_frame.grid(row=2, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(0, 10))
        cryo_params_frame.columnconfigure(1, weight=1)
        
        self.cryo_params_widgets = []
        
        # Ending temperature (target)
        lbl = ttk.Label(cryo_params_frame, text="Target Temperature (K):")
        lbl.grid(row=0, column=0, sticky=tk.W, pady=2)
        self.cryo_temp_end = ttk.Entry(cryo_params_frame, width=15)
        self.cryo_temp_end.insert(0, "10")
        self.cryo_temp_end.grid(row=0, column=1, sticky=tk.W, padx=(10, 0), pady=2)
        self.cryo_params_widgets.append(self.cryo_temp_end)

        # Measurement interval
        lbl = ttk.Label(cryo_params_frame, text="Measurement Interval (K):")
        lbl.grid(row=1, column=0, sticky=tk.W, pady=2)
        self.cryo_meas_interval = ttk.Entry(cryo_params_frame, width=15)
        self.cryo_meas_interval.insert(0, "10")
        self.cryo_meas_interval.grid(row=1, column=1, sticky=tk.W, padx=(10, 0), pady=2)
        self.cryo_params_widgets.append(self.cryo_meas_interval)

        # Ramp rate
        lbl = ttk.Label(cryo_params_frame, text="Ramp Rate (K/min):")
        lbl.grid(row=2, column=0, sticky=tk.W, pady=2)
        self.cryo_ramp_rate = ttk.Entry(cryo_params_frame, width=15)
        self.cryo_ramp_rate.insert(0, "5.0")
        self.cryo_ramp_rate.grid(row=2, column=1, sticky=tk.W, padx=(10, 0), pady=2)
        self.cryo_params_widgets.append(self.cryo_ramp_rate)
        
        self.cryo_params_frame = cryo_params_frame

        # --- Instrument Selection ---
        ttk.Label(main_frame, text="Instrument:", font=("Arial", 10, "bold")).grid(
            row=3, column=0, sticky=tk.W, pady=(10, 5))
        
        self.instrument_var = tk.StringVar(value="PIA")
        instrument_frame = ttk.Frame(main_frame)
        instrument_frame.grid(row=3, column=1, sticky=tk.W, pady=(10, 5))
        
        instruments = [
            ("PIA (Precision Impedance Analyzer)", "PIA"),
            ("PSPA (Parameter/Source Analyzer)", "PSPA"),
            ("LCR (E4980A LCR Meter)", "LCR"),
        ]
        
        for label, value in instruments:
            rb = ttk.Radiobutton(
                instrument_frame, text=label, variable=self.instrument_var,
                value=value, command=self.update_measurement_list
            )
            rb.pack(anchor=tk.W)

        # --- Measurement Selection ---
        ttk.Label(main_frame, text="Measurement:", font=("Arial", 10, "bold")).grid(
            row=4, column=0, sticky=(tk.W, tk.N), pady=(10, 5))
        
        measurement_frame = ttk.Frame(main_frame)
        measurement_frame.grid(row=4, column=1, sticky=(tk.W, tk.E, tk.N, tk.S), pady=(10, 5))
        measurement_frame.columnconfigure(0, weight=1)
        measurement_frame.rowconfigure(0, weight=1)
        
        scrollbar = ttk.Scrollbar(measurement_frame)
        scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))
        
        self.measurement_listbox = tk.Listbox(
            measurement_frame, height=6, yscrollcommand=scrollbar.set
        )
        self.measurement_listbox.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        self.measurement_listbox.bind('<<ListboxSelect>>', self.on_measurement_selected)
        scrollbar.config(command=self.measurement_listbox.yview)

        # --- Parameters Frame (Scrollable) ---
        ttk.Label(main_frame, text="Parameters:", font=("Arial", 10, "bold")).grid(
            row=5, column=0, sticky=(tk.W, tk.N), pady=(10, 5))
        
        params_outer = ttk.Frame(main_frame)
        params_outer.grid(row=5, column=1, sticky=(tk.W, tk.E, tk.N, tk.S), pady=(10, 5))
        params_outer.columnconfigure(0, weight=1)
        params_outer.rowconfigure(0, weight=1)
        
        params_canvas = tk.Canvas(params_outer, highlightthickness=0)
        params_scrollbar = ttk.Scrollbar(params_outer, orient="vertical", command=params_canvas.yview)
        params_scrollable_frame = ttk.Frame(params_canvas)
        
        params_scrollable_frame.bind(
            "<Configure>",
            lambda e: params_canvas.configure(scrollregion=params_canvas.bbox("all"))
        )
        
        params_canvas.create_window((0, 0), window=params_scrollable_frame, anchor="nw")
        params_canvas.configure(yscrollcommand=params_scrollbar.set)
        
        params_canvas.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        params_scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))
        
        self.params_frame = params_scrollable_frame
        self.param_entries = {}

        # --- Measurement Queue Frame (for cryo) ---
        queue_frame = ttk.LabelFrame(main_frame, text="Measurement Queue (Cryo)", padding="10")
        queue_frame.grid(row=6, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(10, 5))
        queue_frame.columnconfigure(0, weight=1)
        
        queue_scrollbar = ttk.Scrollbar(queue_frame)
        queue_scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))
        
        self.queue_listbox = tk.Listbox(
            queue_frame, height=3, yscrollcommand=queue_scrollbar.set
        )
        self.queue_listbox.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        queue_scrollbar.config(command=self.queue_listbox.yview)
        
        self.queue_frame = queue_frame
        self.queue_frame.grid_remove()  # Hide initially
        ttk.Label(main_frame, text="Status:", font=("Arial", 10, "bold")).grid(
            row=7, column=0, sticky=(tk.W, tk.N), pady=(10, 5))
        
        log_frame = ttk.Frame(main_frame)
        log_frame.grid(row=7, column=1, sticky=(tk.W, tk.E, tk.N, tk.S), pady=(10, 5))
        log_frame.columnconfigure(0, weight=1)
        log_frame.rowconfigure(0, weight=1)
        
        log_scrollbar = ttk.Scrollbar(log_frame)
        log_scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))
        
        self.status_text = tk.Text(
            log_frame, height=4, state='disabled', yscrollcommand=log_scrollbar.set
        )
        self.status_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        log_scrollbar.config(command=self.status_text.yview)

        # --- Buttons Frame ---
        button_frame = ttk.Frame(main_frame)
        button_frame.grid(row=8, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(10, 0))
        button_frame.columnconfigure(1, weight=1)
        
        self.run_btn = ttk.Button(button_frame, text="Run Measurement", command=self.run_measurement)
        self.run_btn.grid(row=0, column=0, padx=(0, 5))
        
        self.queue_btn = ttk.Button(button_frame, text="Queue Measurement", 
                                    command=self.queue_measurement)
        self.queue_btn.grid(row=0, column=0, padx=(0, 5))
        self.queue_btn.grid_remove()  # Hide initially
        
        self.start_sweep_btn = ttk.Button(button_frame, text="Start Cryo Sweep",
                                         command=self.start_cryo_sweep)
        self.start_sweep_btn.grid(row=0, column=1, padx=(0, 5))
        self.start_sweep_btn.grid_remove()  # Hide initially
        
        self.clear_queue_btn = ttk.Button(button_frame, text="Clear Queue",
                                         command=self.clear_queue)
        self.clear_queue_btn.grid(row=0, column=2, padx=(0, 5))
        self.clear_queue_btn.grid_remove()  # Hide initially
        
        self.progress_var = tk.DoubleVar()
        progress = ttk.Progressbar(button_frame, variable=self.progress_var, maximum=100)
        progress.grid(row=0, column=3, sticky=(tk.W, tk.E), padx=(0, 5))
        
        exit_btn = ttk.Button(button_frame, text="Exit", command=self.root.quit)
        exit_btn.grid(row=0, column=4)

    def populate_default_params(self):
        """Populate parameter frame with default parameters based on selected measurement."""
        # Clear existing parameters
        for widget in self.params_frame.winfo_children():
            widget.destroy()
        self.param_entries.clear()

        # Get selected measurement
        try:
            instrument = self.instrument_var.get()
            selection = self.measurement_listbox.curselection()
            if not selection:
                return
            
            measurement_idx = selection[0] + 1  # 1-indexed
            params = self.MEASUREMENT_PARAMS.get((instrument, measurement_idx), [])
        except (tk.TclError, IndexError):
            params = []

        # Add parameter fields
        if not params:
            label = ttk.Label(self.params_frame, text="(No parameters for this measurement)", 
                            font=("Arial", 9, "italic"), foreground="gray")
            label.grid(row=0, column=0, columnspan=2, sticky=tk.W, pady=5)
        else:
            for label_text, key, default in params:
                self.add_param_field(label_text, key, default)

    def add_param_field(self, label, key, default=""):
        """Add a single parameter input field."""
        row = len(self.param_entries)
        ttk.Label(self.params_frame, text=label, font=("Arial", 9)).grid(
            row=row, column=0, sticky=tk.W, pady=2)
        entry = ttk.Entry(self.params_frame, width=20)
        entry.insert(0, default)
        entry.grid(row=row, column=1, sticky=tk.W, padx=(10, 0), pady=2)
        self.param_entries[key] = entry

    def update_measurement_list(self):
        """Update measurement listbox based on selected instrument."""
        self.measurement_listbox.delete(0, tk.END)
        instrument = self.instrument_var.get()
        
        if instrument == "PIA":
            measurements = self.PIA_MEASUREMENTS
        elif instrument == "PSPA":
            measurements = self.PSPA_MEASUREMENTS
        else:  # LCR
            measurements = self.LCR_MEASUREMENTS
        
        for meas in measurements:
            self.measurement_listbox.insert(tk.END, meas)
        
        if measurements:
            self.measurement_listbox.selection_set(0)
            self.populate_default_params()

    def on_measurement_selected(self, event):
        """Handle measurement selection change."""
        self.populate_default_params()

    def browse_output_dir(self):
        """Browse for output directory."""
        folder = filedialog.askdirectory(title="Select Output Directory")
        if folder:
            self.output_dir_var.set(folder)

    def toggle_cryo_params(self):
        """Show/hide cryo parameters based on checkbox state."""
        if self.cryo_var.get():
            self.cryo_params_frame.grid()
            self.queue_frame.grid()
            self.run_btn.grid_remove()
            self.queue_btn.grid()
            self.start_sweep_btn.grid()
            self.clear_queue_btn.grid()
            self.cryo_enabled = True
            self.update_status("Cryo mode enabled - measurements will be queued")
        else:
            self.cryo_params_frame.grid_remove()
            self.queue_frame.grid_remove()
            self.run_btn.grid()
            self.queue_btn.grid_remove()
            self.start_sweep_btn.grid_remove()
            self.clear_queue_btn.grid_remove()
            self.measurement_queue.clear()
            self.queue_listbox.delete(0, tk.END)
            self.cryo_enabled = False
            self.update_status("Cryo mode disabled - running single measurements")

    def queue_measurement(self):
        """Queue the selected measurement (with current params) for cryo sweep."""
        try:
            instrument = self.instrument_var.get()
            measurement_idx, measurement_name = self.get_selected_measurement()
            params = self.get_params_dict()  # Capture params at queue time
            
            self.measurement_queue.append((instrument, measurement_idx, measurement_name, params))
            self.queue_listbox.insert(tk.END, f"[{instrument}] {measurement_name}")
            self.update_status(f"Queued: [{instrument}] {measurement_name}")
        except ValueError as e:
            messagebox.showerror("Input Error", str(e))

    def clear_queue(self):
        """Clear all queued measurements."""
        if not self.measurement_queue:
            messagebox.showinfo("Queue Empty", "No measurements queued.")
            return
        
        if messagebox.askyesno("Clear Queue", "Clear all queued measurements?"):
            self.measurement_queue.clear()
            self.queue_listbox.delete(0, tk.END)
            self.update_status("Queue cleared")

    def start_cryo_sweep(self):
        """Start the cryogenic temperature sweep."""
        if not self.measurement_queue:
            messagebox.showerror("Empty Queue", "Please queue at least one measurement before starting sweep.")
            return
        
        try:
            import cryo
            # Get cryo parameters
            temp_end = float(self.cryo_temp_end.get())
            meas_interval = float(self.cryo_meas_interval.get())
            ramp_rate = float(self.cryo_ramp_rate.get())
            
            # Validate
            if meas_interval <= 0:
                raise ValueError("Measurement interval must be > 0")
            if ramp_rate <= 0:
                raise ValueError("Ramp rate must be > 0")

            if temp_end < cryo.TEMP_MIN_K or temp_end > cryo.TEMP_MAX_K:
                raise ValueError("Target temperature is out of range")
            
            output_dir = self.output_dir_var.get()
            file_management.set_output_dir(output_dir)
            file_management.ensure_output_dir(output_dir)
            
            self.update_status(f"Starting cryo sweep: current temperature → {temp_end}K")
            self.update_status(f"Ramp rate: {ramp_rate} K/min")
            self.update_status(f"Measurement interval: {meas_interval} K")
            self.update_status(f"Queued measurements: {len(self.measurement_queue)}")
            
            # Disable buttons during sweep
            self.queue_btn.config(state='disabled')
            self.start_sweep_btn.config(state='disabled')
            self.clear_queue_btn.config(state='disabled')
            self.progress_var.set(50)
            
            # Run sweep in background thread
            thread = threading.Thread(
                target=self._run_cryo_sweep_thread,
                args=(temp_end, meas_interval, ramp_rate)
            )
            thread.daemon = True
            thread.start()
        
        except ValueError as e:
            messagebox.showerror("Input Error", str(e))
        except Exception as e:
            messagebox.showerror("Error", f"Error: {str(e)}")
            self.queue_btn.config(state='normal')
            self.start_sweep_btn.config(state='normal')
            self.clear_queue_btn.config(state='normal')
            self.progress_var.set(0)

    def _run_cryo_sweep_thread(self, temp_end, meas_interval, ramp_rate):
        """Execute cryo sweep in background thread."""
        def _restore_buttons():
            self.queue_btn.config(state='normal')
            self.start_sweep_btn.config(state='normal')
            self.clear_queue_btn.config(state='normal')

        try:
            import cryo
            # Import execute functions from cli
            from gui_measurements import execute_pia_gui, execute_pspa_gui, execute_lcr_gui
            
            self.update_status("Connecting to cryogenic controller...")
            cryocon = cryo.setup()

            start_temp = cryo.get_current_temperature(cryocon, loop=1)
            if start_temp is None:
                raise RuntimeError("Unable to read current temperature from controller")
            if start_temp < temp_end - cryo.TEMP_TOLERANCE_K:
                raise RuntimeError(
                    f"Current temperature ({start_temp:.2f} K) is below target ({temp_end:.2f} K). "
                    "Cryo sweep only supports ramping down."
                )

            temp_points = cryo.generate_temperature_sweep_points(start_temp, temp_end, meas_interval)
            n_temps = len(temp_points)
            n_meas = len(self.measurement_queue)

            self.update_status(f"Temperature sweep: {temp_points[0]:.1f}K → {temp_points[-1]:.1f}K")
            self.update_status(f"Temperature points: {n_temps}")
            
            for temp_idx, target_temp in enumerate(temp_points):
                self.update_status(f"\n--- Temperature {temp_idx + 1}/{n_temps}: {target_temp:.2f}K ---")
                
                try:
                    # Ramp to target temperature
                    self.update_status(f"Ramping to {target_temp:.2f}K...")
                    cryo.resume_ramp_to_next(cryocon, target_temp, ramp_rate, loop=1)
                    
                    # Wait for stabilization
                    success = cryo.wait_for_temperature(
                        cryocon, target_temp,
                        tolerance_k=cryo.TEMP_TOLERANCE_K,
                        stability_time_s=cryo.TEMP_STABILITY_TIME_S,
                        loop=1
                    )
                    
                    if not success:
                        self.update_status(f"Warning: Temperature stabilization uncertain at {target_temp:.2f}K")
                    
                    # Hold temperature during measurements
                    cryo.hold_temperature(cryocon, loop=1)
                    current_temp = cryo.get_current_temperature(cryocon, loop=1)
                    if current_temp:
                        self.update_status(f"Confirmed temperature: {current_temp:.2f}K")
                    
                    # Create subdirectory for this temperature
                    base_dir = file_management.output_dir
                    temp_subdir = os.path.join(base_dir, "cryo_sweep", f"{target_temp:.1f}K")
                    file_management.set_output_dir(temp_subdir)
                    file_management.ensure_output_dir(temp_subdir)
                    
                    # Run queued measurements
                    for m_idx, (instrument, meas_idx, meas_name, meas_params) in enumerate(self.measurement_queue, 1):
                        self.update_status(f"  [{m_idx}/{n_meas}] {meas_name}...")
                        try:
                            if instrument == "PIA":
                                execute_pia_gui(meas_idx, meas_params)
                            elif instrument == "PSPA":
                                execute_pspa_gui(meas_idx, meas_params)
                            elif instrument == "LCR":
                                execute_lcr_gui(meas_idx, meas_params)
                            self.update_status(f"  [{m_idx}/{n_meas}] Complete")
                        except Exception as e:
                            log.error("Measurement failed: %s", e)
                            self.update_status(f"  [{m_idx}/{n_meas}] FAILED: {e}")
                    
                    progress = ((temp_idx + 1) / n_temps) * 100
                    self.root.after(0, lambda p=progress: self.progress_var.set(p))
                
                except Exception as e:
                    log.error("Temperature point failed: %s", e)
                    self.update_status(f"Error at {target_temp:.2f}K: {e}")
            
            # Cleanup
            try:
                cryo.disable_ramp(cryocon, loop=1)
                cryo.disconnect_cryocon(cryocon)
            except:
                pass
            
            self.update_status("\n✓ Cryogenic sweep complete!")
            self.root.after(0, lambda: self.progress_var.set(100))
            
        except Exception as e:
            log.error("Cryo sweep error: %s", e)
            error_msg = str(e)
            self.update_status(f"✗ Cryo sweep failed: {error_msg}")
            self.root.after(0, lambda: [
                self.progress_var.set(0),
                messagebox.showerror("Cryo Sweep Error",
                    f"Cryogenic sweep failed:\n\n{error_msg}")
            ])
        finally:
            # Re-enable buttons on main thread
            self.root.after(0, _restore_buttons)

    def get_selected_measurement(self):
        """Get the selected measurement index and name."""
        selection = self.measurement_listbox.curselection()
        if not selection:
            raise ValueError("No measurement selected")
        idx = selection[0]
        name = self.measurement_listbox.get(idx)
        return idx + 1, name  # 1-indexed for compatibility

    def get_params_dict(self):
        """Get all parameter values as a dictionary."""
        params = {}
        for key, entry in self.param_entries.items():
            value = entry.get()
            # Try to convert to float, falls back to string
            try:
                params[key] = float(value)
            except ValueError:
                params[key] = value
        return params

    def update_status(self, message):
        """Update the status text area (thread-safe)."""
        def _update():
            self.status_text.config(state='normal')
            self.status_text.insert(tk.END, message + "\n")
            self.status_text.see(tk.END)
            self.status_text.config(state='disabled')
        # Schedule on main thread if called from a background thread
        self.root.after(0, _update)

    def run_measurement(self):
        """Run the selected measurement in a separate thread."""
        try:
            # Validate selections
            instrument = self.instrument_var.get()
            measurement_idx, measurement_name = self.get_selected_measurement()
            output_dir = self.output_dir_var.get()

            # Set output directory
            file_management.set_output_dir(output_dir)
            file_management.ensure_output_dir(output_dir)

            self.update_status(f"Starting measurement: {measurement_name}")
            self.update_status(f"Output directory: {output_dir}")
            
            # Run measurement in background thread
            self.run_btn.config(state='disabled')
            self.progress_var.set(50)
            
            thread = threading.Thread(
                target=self._run_measurement_thread,
                args=(instrument, measurement_idx, measurement_name)
            )
            thread.daemon = True
            thread.start()

        except ValueError as e:
            messagebox.showerror("Input Error", str(e))
        except Exception as e:
            messagebox.showerror("Error", f"Error: {str(e)}")
            self.run_btn.config(state='normal')
            self.progress_var.set(0)

    def _run_measurement_thread(self, instrument, measurement_idx, measurement_name):
        """Execute measurement in background thread."""
        try:
            params = self.get_params_dict()
            
            # Call the measurement executor
            self.measurement_executor(instrument, measurement_idx, params)
            
            self.update_status(f"✓ Measurement completed: {measurement_name}")
            self.root.after(0, lambda: self.progress_var.set(100))
            
            # Reset after a delay
            self.root.after(1000, lambda: [
                self.progress_var.set(0),
                self.run_btn.config(state='normal')
            ])

        except Exception as e:
            log.error("Measurement execution failed: %s", e)
            error_msg = str(e)
            self.update_status(f"✗ Error: {error_msg}")
            self.root.after(0, lambda: [
                self.progress_var.set(0),
                self.run_btn.config(state='normal'),
                messagebox.showerror("Measurement Error",
                    f"Measurement '{measurement_name}' failed:\n\n{error_msg}")
            ])


def run_gui_mode(measurement_executor):
    """Launch the GUI application."""
    root = tk.Tk()
    gui = MeasurementGUI(root, measurement_executor)
    root.mainloop()
