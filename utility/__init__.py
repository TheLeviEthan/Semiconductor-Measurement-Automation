"""
Utility package for Semiconductor Measurement Automation.

This folder holds helper modules shared across the entire project:
  - config.py           – Loads settings from config.yaml (GPIB addresses, defaults)
  - logging_config.py   – Sets up the logging system (console + optional file)
  - gpib_utils.py       – Safe connect/disconnect for instruments, user-input helpers
  - file_management.py  – Saves CSV data and plot images, manages the output folder
  - parameter_validators.py – Checks that measurement parameters are within safe limits
  - measurement_results.py  – Data classes that hold measurement results in a standard format
  - find_instruments.py – Discovers VISA/GPIB instruments connected to the computer
  - gpib_identification.py – Quick one-liner to list GPIB instruments

None of these files talk to instruments directly; they provide building-block
functions used by the measurement modules (pia.py, pspa.py, lcr.py, cryo.py)
and the user interfaces (gui.py, cli.py).
"""
