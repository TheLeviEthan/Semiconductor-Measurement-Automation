"""
Filename: main.py
Author: Ethan Ruddell
Date: 2025-01-23
Description: Contains all constants and functions for PIA measurements.
"""

import pyvisa
import numpy as np


# =============================
# User settings and constants
# =============================
GPIB_ADDRESS = ""   # FIND AND RECORD GPIB ADDRESS FOR PSPA

# Default Frequency sweep settings (LOG sweep)
FREQ_START_HZ = 1e4      # start frequency, Hz (>0 for LOG sweep)
FREQ_STOP_HZ = 1e6       # stop frequency, Hz
NUM_POINTS = 201         # number of points in LOG sweep

# FIND MANUAL FOR PSPA TO GET SYNTAX (keysight)
# UTILIZES STANDARD SCPI COMMANDS
# https://en.wikipedia.org/wiki/Standard_Commands_for_Programmable_Instruments


