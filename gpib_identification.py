"""
Quick-check script: prints every VISA/GPIB instrument the computer can see.

Run this from a terminal to verify that your GPIB interface card and cables
are working before attempting any measurements.  The output will look like:
    ('GPIB0::16::INSTR', 'GPIB0::24::INSTR', ...)
Each string is the address of one instrument.
"""

import pyvisa

rm = pyvisa.ResourceManager()
print(rm.list_resources())