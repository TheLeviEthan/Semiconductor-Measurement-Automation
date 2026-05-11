"""
Filename: gui.py
Author: Ethan Ruddell
Date: 2026-5-11
Description: Graphical User Interface (GUI) for the measurement automation.

This file builds a window using Python's built-in tkinter library.  The window
lets you:
  1. Choose an instrument (PIA, PSPA, or LCR).
  2. Pick a measurement from a scrollable list.
  3. Fill in parameters (frequency, voltage, oscillator level, etc.) in dynamically generated
     text fields.
  4. Optionally enable a cryogenic temperature sweep, which queues multiple
     measurements to be run at each temperature point.
  5. Click "Run Measurement" (or "Start Cryo Sweep") to begin.

Measurements run in a background thread so the window stays responsive.
Status updates and errors are displayed in real-time in the status area
at the bottom of the window.

The GUI does NOT contain any instrument communication code.  It delegates
all actual measurements to gui_measurements.py, which in turn calls the
instrument driver modules (pia.py, pspa.py, lcr.py, cryo.py).
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from contextlib import contextmanager, redirect_stderr, redirect_stdout
import sys
import os
import logging
import threading

# Add the 'utility' and 'measurement functions' folders to Python's search
# path so that module imports work without package prefixes.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'utility'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'measurement functions'))

import file_management
import logging_config
import usb_switchbox
from measurements_config import (
    PIA_MEASUREMENTS, PSPA_MEASUREMENTS, LCR_MEASUREMENTS,
    MEASUREMENT_PARAMS, get_measurements_list,
    normalize_pspa_choice, get_measurement_help, get_parameter_help,
)

log = logging.getLogger(__name__)

class HoverTooltip:
    """Small semi-transparent tooltip that appears after a short hover delay."""

    def __init__(self, widget, text, delay_ms=500, wraplength=320):
        self.widget = widget
        self.text = text
        self.delay_ms = delay_ms
        self.wraplength = wraplength
        self._after_id = None
        self._tip_window = None

        widget.bind("<Enter>", self._schedule_show, add="+")
        widget.bind("<Leave>", self._hide, add="+")
        widget.bind("<ButtonPress>", self._hide, add="+")

    def _schedule_show(self, event=None):
        self._cancel_scheduled_show()
        self._after_id = self.widget.after(self.delay_ms, self._show)

    def _cancel_scheduled_show(self):
        if self._after_id is not None:
            try:
                self.widget.after_cancel(self._after_id)
            except Exception:
                pass
            self._after_id = None

    def _show(self):
        self._after_id = None
        if self._tip_window is not None or not self.widget.winfo_exists():
            return

        x, y = self.widget.winfo_pointerxy()
        tip = tk.Toplevel(self.widget)
        tip.wm_overrideredirect(True)
        tip.wm_attributes("-topmost", True)
        try:
            tip.wm_attributes("-alpha", 0.92)
        except tk.TclError:
            pass
        tip.geometry(f"+{x + 14}+{y + 14}")

        frame = tk.Frame(tip, bg="#f5f5f5", bd=1, relief="solid")
        frame.pack(fill="both", expand=True)
        label = tk.Label(
            frame,
            text=self.text,
            bg="#f5f5f5",
            fg="#222222",
            justify="left",
            anchor="w",
            wraplength=self.wraplength,
            padx=8,
            pady=6,
        )
        label.pack(fill="both", expand=True)
        self._tip_window = tip

    def _hide(self, event=None):
        self._cancel_scheduled_show()
        if self._tip_window is not None:
            try:
                self._tip_window.destroy()
            except Exception:
                pass
            self._tip_window = None


@contextmanager
def _suppress_console_output():
    """Silence stdout/stderr while a background measurement is running."""
    with open(os.devnull, "w") as sink:
        with redirect_stdout(sink), redirect_stderr(sink):
            yield


class MeasurementGUI:
    """
    Main GUI class for semiconductor measurement automation.

    This class manages the entire window: instrument selection, measurement
    selection, parameter entry, optional cryogenic sweep, and status display.
    All long-running operations (measurements, cryo sweeps) are executed in
    background threads so the window never freezes.
    """

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
        self.root.title("NRG Semiconductor Measurement Automation - Ethan Ruddell - University of Florida")
        self.root.geometry("900x700")
        
        # Setup logging
        logging_config.setup(output_dir=file_management.default_output_dir)
        
        # Track cryo state and measurement queue
        self.cryo_enabled = False
        self.measurement_queue = []

        # USB switchbox (multiplexer) control
        self.switchbox = usb_switchbox.create_switchbox_from_config()
        self.switchbox_enabled_var = tk.BooleanVar(value=self.switchbox.enabled)
        self.switchbox_status_var = tk.StringVar()
        self._instrument_buttons_layout_in_progress = False
        self._instrument_buttons_last_width = None
        self.graph_window = None
        self.graph_canvas = None
        self.graph_photo = None
        self._tooltips = []
        self.measurement_choice_var = tk.IntVar(value=1)
        self.measurement_rows = []
        
        self.create_widgets()
        self.toggle_cryo_params()  # Hide cryo UI initially
        self.update_measurement_list()
        self.update_switchbox_status()
        if self.switchbox_enabled_var.get():
            self.switch_now_from_selection(raise_on_error=False)
        self.root.protocol("WM_DELETE_WINDOW", self.on_close)

    def create_widgets(self):
        """Build all visual elements of the window (labels, buttons, text fields, etc.)."""
        # --- Scrollable Main Area ---
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_container = ttk.Frame(self.root)
        main_container.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        main_container.columnconfigure(0, weight=1)
        main_container.rowconfigure(0, weight=1)

        self.main_canvas = tk.Canvas(main_container, highlightthickness=0)
        self.main_scrollbar = ttk.Scrollbar(main_container, orient="vertical", command=self.main_canvas.yview)
        self.main_canvas.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        self.main_scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))
        self.main_canvas.configure(yscrollcommand=self.main_scrollbar.set)

        self.main_frame = ttk.Frame(self.main_canvas, padding="10")
        self.main_canvas_window = self.main_canvas.create_window((0, 0), window=self.main_frame, anchor="nw")
        self.main_frame.columnconfigure(1, weight=1)
        # Configure rows to expand/contract proportionally
        self.main_frame.rowconfigure(4, weight=1)  # Measurement selection
        self.main_frame.rowconfigure(5, weight=2)  # Parameters (largest)
        self.main_frame.rowconfigure(7, weight=2)  # Status

        self.main_frame.bind(
            "<Configure>",
            lambda event: self.main_canvas.configure(scrollregion=self.main_canvas.bbox("all")),
        )
        self.main_canvas.bind(
            "<Configure>",
            lambda event: self.main_canvas.itemconfigure(self.main_canvas_window, width=event.width),
        )
        self._bind_mousewheel(self.main_canvas)
        self.root.bind_all("<MouseWheel>", self._on_main_mousewheel, add="+")
        self.root.bind_all("<Shift-MouseWheel>", self._on_main_mousewheel, add="+")
        self.root.bind_all("<Button-4>", self._on_main_mousewheel, add="+")
        self.root.bind_all("<Button-5>", self._on_main_mousewheel, add="+")

        # --- Output Directory ---
        ttk.Label(self.main_frame, text="Output Directory:", font=("Arial", 10, "bold")).grid(
            row=0, column=0, sticky=tk.W, pady=(0, 5))
        output_frame = ttk.Frame(self.main_frame)
        output_frame.grid(row=0, column=1, sticky=(tk.W, tk.E), pady=(0, 5))
        output_frame.columnconfigure(0, weight=1)
        
        self.output_dir_var = tk.StringVar(value=file_management.default_output_dir)
        output_entry = ttk.Entry(output_frame, textvariable=self.output_dir_var, state='readonly')
        output_entry.grid(row=0, column=0, sticky=(tk.W, tk.E))
        
        browse_btn = ttk.Button(output_frame, text="Browse", command=self.browse_output_dir)
        browse_btn.grid(row=0, column=1, padx=(5, 0))

        # --- Cryo Integration ---
        self.cryo_var = tk.BooleanVar(value=False)
        cryo_check = ttk.Checkbutton(self.main_frame, text="Enable Cryogenic Temperature Sweep",
                                     variable=self.cryo_var, command=self.toggle_cryo_params)
        cryo_check.grid(row=1, column=0, columnspan=2, sticky=tk.W, pady=(10, 5))

        # --- Cryo Parameters Frame (initially hidden) ---
        cryo_params_frame = ttk.LabelFrame(self.main_frame, text="Cryogenic Parameters", padding="10")
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
        ttk.Label(self.main_frame, text="Instrument:", font=("Arial", 10, "bold")).grid(
            row=3, column=0, sticky=tk.W, pady=(10, 5))
        
        self.instrument_var = tk.StringVar(value="PIA")
        instrument_frame = ttk.Frame(self.main_frame)
        instrument_frame.grid(row=3, column=1, sticky=(tk.W, tk.E), pady=(10, 5))
        instrument_frame.columnconfigure(0, weight=1)

        tools_frame = ttk.Frame(instrument_frame)
        tools_frame.grid(row=0, column=0, sticky=(tk.W, tk.E))

        self.instrument_tools_frame = tools_frame
        self.instrument_buttons = []
        
        instruments = [
            ("PIA (Precision Impedance Analyzer)", "PIA"),
            ("PSPA (Parameter/Source Analyzer)", "PSPA"),
            ("LCR (E4980A LCR Meter)", "LCR"),
            ("Ferro (Switchbox Position D)", "FERRO"),
        ]
        
        for label, value in instruments:
            rb = ttk.Radiobutton(
                tools_frame, text=label, variable=self.instrument_var,
                value=value, command=self.on_instrument_changed
            )
            self.instrument_buttons.append(rb)

        switchbox_frame = ttk.LabelFrame(instrument_frame, text="USB Switchbox", padding="8")
        switchbox_frame.grid(row=1, column=0, sticky=(tk.W, tk.E), pady=(8, 0))
        self.switchbox_frame = switchbox_frame

        switchbox_enable = ttk.Checkbutton(
            switchbox_frame,
            text="Enable automatic routing",
            variable=self.switchbox_enabled_var,
            command=self.on_switchbox_toggle,
        )
        switchbox_enable.grid(row=0, column=0, sticky=tk.W)

        switch_now_btn = ttk.Button(
            switchbox_frame,
            text="Switch Now",
            command=self.switch_now_from_selection,
        )
        switch_now_btn.grid(row=0, column=1, sticky=tk.E, padx=(10, 0))

        self.switchbox_status_label = ttk.Label(
            switchbox_frame,
            textvariable=self.switchbox_status_var,
        )
        self.switchbox_status_label.grid(row=1, column=0, columnspan=2, sticky=tk.W, pady=(5, 0))

        self.instrument_tools_frame.bind("<Configure>", self._layout_instrument_buttons)
        self.root.after_idle(self._layout_instrument_buttons)

        # --- Measurement Selection ---
        ttk.Label(self.main_frame, text="Measurement:", font=("Arial", 10, "bold")).grid(
            row=4, column=0, sticky=(tk.W, tk.N), pady=(10, 5))

        measurement_outer = ttk.Frame(self.main_frame)
        measurement_outer.grid(row=4, column=1, sticky=(tk.W, tk.E, tk.N, tk.S), pady=(10, 5))
        measurement_outer.columnconfigure(0, weight=1)
        measurement_outer.rowconfigure(0, weight=1)

        self.measurement_canvas = tk.Canvas(measurement_outer, highlightthickness=0)
        measurement_scrollbar = ttk.Scrollbar(measurement_outer, orient="vertical", command=self.measurement_canvas.yview)
        self.measurement_canvas.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        measurement_scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))
        self.measurement_canvas.configure(yscrollcommand=measurement_scrollbar.set)

        self.measurement_rows_frame = ttk.Frame(self.measurement_canvas)
        self.measurement_canvas.create_window((0, 0), window=self.measurement_rows_frame, anchor="nw")
        self.measurement_rows_frame.bind(
            "<Configure>",
            lambda event: self.measurement_canvas.configure(scrollregion=self.measurement_canvas.bbox("all")),
        )

        # --- Parameters Frame (Scrollable) ---
        ttk.Label(self.main_frame, text="Parameters:", font=("Arial", 10, "bold")).grid(
            row=5, column=0, sticky=(tk.W, tk.N), pady=(10, 5))
        
        params_outer = ttk.Frame(self.main_frame)
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
        self.params_canvas = params_canvas
        self.param_entries = {}
        self.param_entries_list = []  # Track entry order for tab traversal
        
        # Prevent canvas clicks from deselecting focused entry
        params_canvas.bind('<Button-1>', self._on_canvas_click)

        # --- Measurement Queue Frame (for cryo) ---
        queue_frame = ttk.LabelFrame(self.main_frame, text="Measurement Queue (Cryo)", padding="10")
        queue_frame.grid(row=6, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(10, 5))
        queue_frame.columnconfigure(0, weight=1)
        
        queue_scrollbar = ttk.Scrollbar(queue_frame)
        queue_scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))
        
        self.queue_listbox = tk.Listbox(
            queue_frame, height=3, yscrollcommand=queue_scrollbar.set,
            exportselection=False
        )
        self.queue_listbox.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        queue_scrollbar.config(command=self.queue_listbox.yview)
        
        self.queue_frame = queue_frame
        self.queue_frame.grid_remove()  # Hide initially
        ttk.Label(self.main_frame, text="Status:", font=("Arial", 10, "bold")).grid(
            row=7, column=0, sticky=(tk.W, tk.N), pady=(10, 5))
        
        log_frame = ttk.Frame(self.main_frame)
        log_frame.grid(row=7, column=1, sticky=(tk.W, tk.E, tk.N, tk.S), pady=(10, 5))
        log_frame.columnconfigure(0, weight=1)
        log_frame.rowconfigure(0, weight=1)
        
        log_scrollbar = ttk.Scrollbar(log_frame)
        log_scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))
        
        self.status_text = tk.Text(
            log_frame, height=10, state='disabled', yscrollcommand=log_scrollbar.set
        )
        self.status_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        log_scrollbar.config(command=self.status_text.yview)

        # --- Buttons Frame ---
        button_frame = ttk.Frame(self.main_frame)
        button_frame.grid(row=9, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(10, 0))
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
        
        exit_btn = ttk.Button(button_frame, text="Exit", command=self.on_close)
        exit_btn.grid(row=0, column=4)

    def update_switchbox_status(self, extra_text=""):
        """Refresh switchbox status text in GUI."""
        state = "Enabled" if self.switchbox_enabled_var.get() else "Disabled"
        ports = self.switchbox.settings.ports or [self.switchbox.settings.port]
        base = f"Status: {state} | Ports: {', '.join(ports)}"
        if extra_text:
            base = f"{base} | {extra_text}"
        self.switchbox_status_var.set(base)

    def _layout_instrument_buttons(self, event=None):
        """Lay out instrument buttons in rows that wrap with the frame width."""
        if not getattr(self, "instrument_buttons", None):
            return

        if self._instrument_buttons_layout_in_progress:
            return

        frame = getattr(self, "instrument_tools_frame", None)
        if frame is None or not frame.winfo_exists():
            return

        frame_width = frame.winfo_width()
        if frame_width <= 1:
            frame_width = frame.winfo_reqwidth()
        if frame_width <= 1:
            frame_width = self.root.winfo_width()
        if frame_width <= 1:
            frame_width = 800

        if self._instrument_buttons_last_width == frame_width:
            return

        self._instrument_buttons_layout_in_progress = True
        try:
            for button in self.instrument_buttons:
                button.grid_forget()

            horizontal_pad = 12
            vertical_pad = 4
            current_row = 0
            current_col = 0
            remaining_width = frame_width

            for button in self.instrument_buttons:
                button.update_idletasks()
                button_width = button.winfo_reqwidth()
                needed_width = button_width if current_col == 0 else button_width + horizontal_pad

                if current_col > 0 and needed_width > remaining_width:
                    current_row += 1
                    current_col = 0
                    remaining_width = frame_width
                    needed_width = button_width

                button.grid(
                    row=current_row,
                    column=current_col,
                    sticky=tk.W,
                    padx=(0, horizontal_pad if current_col > 0 else 0),
                    pady=(0, vertical_pad),
                )
                remaining_width -= needed_width
                current_col += 1
        finally:
            self._instrument_buttons_last_width = frame_width
            self._instrument_buttons_layout_in_progress = False

    def _create_help_icon(self, parent, tooltip_text):
        """Create a small grey question-mark icon with hover help text."""
        canvas = tk.Canvas(parent, width=16, height=16, highlightthickness=0, bd=0, relief="flat")
        canvas.configure(cursor="question_arrow")
        canvas.create_oval(2, 2, 14, 14, fill="#8a8a8a", outline="#8a8a8a")
        canvas.create_text(8, 8, text="?", fill="white", font=("Arial", 9, "bold"))
        self._tooltips.append(HoverTooltip(canvas, tooltip_text))
        return canvas

    def _bind_tab_skip(self, widget):
        """Make Tab/Shift-Tab jump to the parameter area from a measurement control."""
        widget.bind('<Tab>', self._on_measurement_tab_pressed)
        widget.bind('<Shift-Tab>', self._on_measurement_shift_tab_pressed)

    def _bind_mousewheel(self, widget):
        """Bind mousewheel scrolling for the main window scroll region."""
        widget.bind("<MouseWheel>", self._on_main_mousewheel)
        widget.bind("<Shift-MouseWheel>", self._on_main_mousewheel)
        widget.bind("<Button-4>", self._on_main_mousewheel)
        widget.bind("<Button-5>", self._on_main_mousewheel)

    def _on_main_mousewheel(self, event):
        """Scroll the main GUI canvas when the mouse wheel is used."""
        delta = 0
        if hasattr(event, "delta") and event.delta:
            delta = -1 if event.delta > 0 else 1
        elif getattr(event, "num", None) == 4:
            delta = -1
        elif getattr(event, "num", None) == 5:
            delta = 1

        if delta:
            self.main_canvas.yview_scroll(delta, "units")
        return "break"

    def _focus_first_param_entry(self):
        """Move focus to the first parameter field, if any exist."""
        if self.param_entries_list:
            self.param_entries_list[0].focus_set()
            self.param_entries_list[0].select_range(0, tk.END)
            return True
        try:
            self.run_btn.focus_set()
        except Exception:
            pass
        return False

    def on_switchbox_toggle(self):
        """Enable/disable switchbox use from GUI checkbox."""
        enabled = self.switchbox_enabled_var.get()
        self.switchbox.set_enabled(enabled)
        if not enabled:
            self.switchbox.close()
            self.update_switchbox_status("Routing inactive")
            self.update_status("USB switchbox disabled")
            return

        self.update_switchbox_status("Routing active")
        self.update_status("USB switchbox enabled")
        self.switch_now_from_selection(raise_on_error=False)

    def on_instrument_changed(self):
        """Refresh measurements and auto-route switchbox when instrument changes."""
        self.update_measurement_list()
        self.switch_now_from_selection(raise_on_error=False)

    def switch_now_from_selection(self, raise_on_error=False):
        """Route switchbox to currently selected instrument."""
        instrument = self.instrument_var.get()
        self.route_switchbox(instrument, raise_on_error=raise_on_error)

    def route_switchbox(self, instrument, raise_on_error=True):
        """Route USB switchbox to the requested instrument if enabled."""
        if not self.switchbox_enabled_var.get():
            return

        self.switchbox.set_enabled(True)
        try:
            message = self.switchbox.switch_to_instrument(instrument)
            self.update_status(message)
            self.update_switchbox_status(f"Routed to {instrument}")
        except Exception as e:
            error_msg = f"Switchbox routing failed for {instrument}: {e}"
            self.update_status(error_msg)
            self.update_switchbox_status("Routing error")
            if raise_on_error:
                raise RuntimeError(error_msg)

    def on_close(self):
        """Cleanup resources before closing the GUI."""
        try:
            self.switchbox.close()
        except Exception:
            pass
        try:
            if self.graph_window is not None and self.graph_window.winfo_exists():
                self.graph_window.destroy()
        except Exception:
            pass
        self.root.destroy()

    def _on_canvas_click(self, event):
        """Handle canvas click to focus on entry fields if clicked on empty space."""
        # Only focus the first entry if clicking in empty space (not on a widget)
        widget = event.widget.winfo_containing(event.x_root, event.y_root)
        if widget == event.widget and self.param_entries_list:
            # Clicked on empty canvas area, focus first entry
            self.param_entries_list[0].focus_set()
            return 'break'  # Prevent default behavior
        return 'break'  # Prevent canvas from stealing focus

    def _on_measurement_tab_pressed(self, event):
        """Move from measurement selection to the first parameter field."""
        self._focus_first_param_entry()
        return 'break'

    def _on_measurement_shift_tab_pressed(self, event):
        """Move backward from measurement selection to the instrument area."""
        try:
            selection = int(self.measurement_choice_var.get())
            if self.measurement_rows and 1 <= selection <= len(self.measurement_rows):
                self.measurement_rows[selection - 1].focus_set()
            elif self.measurement_rows:
                self.measurement_rows[0].focus_set()
            else:
                self.instrument_buttons[0].focus_set()
        except Exception:
            try:
                self.root.focus_set()
            except Exception:
                pass
        return 'break'
    
    def _on_tab_pressed(self, event):
        """Handle Tab key to advance to next entry field."""
        if not self.param_entries_list:
            return

        current_widget = event.widget
        # If widget is not in our tracked entries (e.g. focus came from listbox),
        # jump to the first parameter entry.
        if current_widget not in self.param_entries_list:
            self.param_entries_list[0].focus_set()
            self.param_entries_list[0].select_range(0, tk.END)
            return 'break'

        try:
            current_idx = self.param_entries_list.index(current_widget)
            # Move to next entry; if at end, move to next focusable widget instead
            next_idx = current_idx + 1
            if next_idx < len(self.param_entries_list):
                self.param_entries_list[next_idx].focus_set()
                self.param_entries_list[next_idx].select_range(0, tk.END)
            else:
                # Move focus to the run button as the next logical target
                try:
                    self.run_btn.focus_set()
                except Exception:
                    pass
            return 'break'  # Prevent default tab behavior
        except (ValueError, IndexError):
            pass
    
    def _on_shift_tab_pressed(self, event):
        """Handle Shift+Tab key to go to previous entry field."""
        if not self.param_entries_list:
            return
        
        current_widget = event.widget
        # If shift-tab from an entry that is the first in the list, return focus
        # to the measurement listbox and preserve its selection.
        try:
            current_idx = self.param_entries_list.index(current_widget)
            if current_idx == 0:
                try:
                    # Focus the measurement listbox and keep the current selection
                    if self.measurement_rows:
                        selection = int(self.measurement_choice_var.get())
                        if 1 <= selection <= len(self.measurement_rows):
                            self.measurement_rows[selection - 1].focus_set()
                        else:
                            self.measurement_rows[0].focus_set()
                except Exception:
                    pass
                return 'break'

            prev_idx = current_idx - 1
            if prev_idx >= 0:
                self.param_entries_list[prev_idx].focus_set()
                self.param_entries_list[prev_idx].select_range(0, tk.END)
                return 'break'
        except (ValueError, IndexError):
            pass

        return None

    def _on_listbox_tab_pressed(self, event):
        """Backward-compatible alias for measurement selection tab handling."""
        return self._on_measurement_tab_pressed(event)

    def _on_listbox_shift_tab_pressed(self, event):
        """Backward-compatible alias for reverse measurement selection tab handling."""
        return self._on_measurement_shift_tab_pressed(event)

    def populate_default_params(self):
        """Populate parameter frame with default parameters based on selected measurement."""
        # Clear existing parameters
        for widget in self.params_frame.winfo_children():
            widget.destroy()
        self.param_entries.clear()
        self.param_entries_list.clear()

        # Get selected measurement
        try:
            instrument = self.instrument_var.get()
            measurement_idx = int(self.measurement_choice_var.get())
            if measurement_idx <= 0:
                return
            if instrument == "PSPA":
                measurement_idx = normalize_pspa_choice(measurement_idx)
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
        label_frame = ttk.Frame(self.params_frame)
        label_frame.grid(row=row, column=0, sticky=tk.W, pady=2)
        ttk.Label(label_frame, text=label, font=("Arial", 9)).grid(row=0, column=0, sticky=tk.W)
        help_icon = self._create_help_icon(label_frame, get_parameter_help(key, label))
        help_icon.grid(row=0, column=1, sticky=tk.W, padx=(5, 0))
        entry = ttk.Entry(self.params_frame, width=20)
        entry.insert(0, default)
        entry.grid(row=row, column=1, sticky=tk.W, padx=(10, 0), pady=2)
        
        # Bind tab/shift+tab for custom focus traversal
        entry.bind('<Tab>', self._on_tab_pressed)
        entry.bind('<Shift-Tab>', self._on_shift_tab_pressed)
        
        self.param_entries[key] = entry
        self.param_entries_list.append(entry)

    def update_measurement_list(self):
        """Update measurement options based on selected instrument."""
        for widget in self.measurement_rows_frame.winfo_children():
            widget.destroy()
        self.measurement_rows.clear()

        instrument = self.instrument_var.get()
        measurements = get_measurements_list(instrument)

        if not measurements:
            self.measurement_choice_var.set(0)
            empty_label = ttk.Label(self.measurement_rows_frame, text="(No measurements available)", foreground="gray")
            empty_label.grid(row=0, column=0, sticky=tk.W)
            return

        self.measurement_choice_var.set(1)
        for idx, meas in enumerate(measurements, 1):
            row_frame = ttk.Frame(self.measurement_rows_frame)
            row_frame.grid(row=idx - 1, column=0, sticky=(tk.W, tk.E), pady=2)
            row_frame.columnconfigure(0, weight=1)

            rb = ttk.Radiobutton(
                row_frame,
                text=meas,
                variable=self.measurement_choice_var,
                value=idx,
                command=self.populate_default_params,
            )
            rb.grid(row=0, column=0, sticky=tk.W)
            self._bind_tab_skip(rb)

            help_icon = self._create_help_icon(row_frame, get_measurement_help(instrument, idx))
            help_icon.grid(row=0, column=1, sticky=tk.W, padx=(6, 0))

            self.measurement_rows.append(rb)

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
        """Show or hide the cryogenic parameter fields when the checkbox is toggled."""
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
        """Add the currently selected measurement (with its parameters) to the cryo queue."""
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
        """Validate cryo parameters and launch the temperature sweep in a background thread."""
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
        """Background thread: ramp temperature down and run queued measurements at each point."""
        def _restore_buttons():
            self.queue_btn.config(state='normal')
            self.start_sweep_btn.config(state='normal')
            self.clear_queue_btn.config(state='normal')

        try:
            import cryo
            # Import execute functions from cli
            from gui_measurements import execute_pia_gui, execute_pspa_gui, execute_lcr_gui
            
            self.update_status("Connecting to cryogenic controller...")
            with _suppress_console_output():
                cryocon = cryo.setup()

            with _suppress_console_output():
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
                    with _suppress_console_output():
                        cryo.resume_ramp_to_next(cryocon, target_temp, ramp_rate, loop=1)
                    
                    # Wait for stabilization
                    with _suppress_console_output():
                        success = cryo.wait_for_temperature(
                            cryocon, target_temp,
                            tolerance_k=cryo.TEMP_TOLERANCE_K,
                            stability_time_s=cryo.TEMP_STABILITY_TIME_S,
                            loop=1
                        )
                    
                    if not success:
                        self.update_status(f"Warning: Temperature stabilization uncertain at {target_temp:.2f}K")
                    
                    # Hold temperature during measurements
                    with _suppress_console_output():
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
                            with _suppress_console_output():
                                self.route_switchbox(instrument, raise_on_error=True)
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
                with _suppress_console_output():
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
        idx = int(self.measurement_choice_var.get())
        if idx <= 0:
            raise ValueError("No measurement selected")
        measurements = get_measurements_list(self.instrument_var.get())
        if idx > len(measurements):
            raise ValueError("Selected measurement is unavailable")
        return idx, measurements[idx - 1]

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
        """Append a line to the status box at the bottom of the window (thread-safe)."""
        def _update():
            self.status_text.config(state='normal')
            self.status_text.insert(tk.END, message + "\n")
            self.status_text.see(tk.END)
            self.status_text.config(state='disabled')
        # Schedule on main thread if called from a background thread
        self.root.after(0, _update)

    def display_image(self, image_path):
        """Display an image from a file path in a separate popup window (thread-safe)."""
        def _display():
            try:
                from PIL import Image, ImageTk
                
                # Load image
                img = Image.open(image_path)

                if self.graph_window is None or not self.graph_window.winfo_exists():
                    self.graph_window = tk.Toplevel(self.root)
                    self.graph_window.title("Measurement Graph")
                    self.graph_window.geometry("920x620")
                    self.graph_window.minsize(640, 420)
                    self.graph_window.transient(self.root)

                    container = ttk.Frame(self.graph_window, padding=10)
                    container.pack(fill="both", expand=True)
                    container.rowconfigure(0, weight=1)
                    container.columnconfigure(0, weight=1)

                    self.graph_canvas = tk.Canvas(container, bg="white", highlightthickness=0)
                    self.graph_canvas.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

                    self.graph_window.protocol("WM_DELETE_WINDOW", self.graph_window.withdraw)

                canvas_width = self.graph_canvas.winfo_width()
                canvas_height = self.graph_canvas.winfo_height()
                if canvas_width <= 1:
                    canvas_width = 900
                if canvas_height <= 1:
                    canvas_height = 600

                img.thumbnail((canvas_width - 10, canvas_height - 10), Image.Resampling.LANCZOS)

                self.graph_photo = ImageTk.PhotoImage(img)

                self.graph_canvas.delete("all")
                self.graph_canvas.create_image(
                    canvas_width // 2, canvas_height // 2,
                    image=self.graph_photo
                )
                self.graph_window.deiconify()
                self.graph_window.lift()
                self.graph_window.focus_force()
            except Exception as e:
                log.error(f"Failed to display image: {e}")
                if self.graph_canvas is not None and self.graph_canvas.winfo_exists():
                    self.graph_canvas.delete("all")
                    self.graph_canvas.create_text(
                        300, 125, text=f"Failed to load image:\n{str(e)}", fill='red'
                    )
        
        # Schedule on main thread if called from a background thread
        self.root.after(0, _display)

    def run_measurement(self):
        """Run the selected measurement in a separate thread."""
        try:
            # Validate selections
            instrument = self.instrument_var.get()
            measurement_idx, measurement_name = self.get_selected_measurement()
            output_dir = self.output_dir_var.get()

            with _suppress_console_output():
                self.route_switchbox(instrument, raise_on_error=True)

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
        """Background thread: execute a single measurement and update the status bar."""
        try:
            params = self.get_params_dict()
            
            # Call the measurement executor and capture returned image path
            with _suppress_console_output():
                image_path = self.measurement_executor(instrument, measurement_idx, params)
            
            self.update_status(f"✓ Measurement completed: {measurement_name}")
            
            # Display the image if one was generated
            if image_path:
                self.display_image(image_path)
                self.update_status(f"  Image saved to: {image_path}")
            
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
