"""
Filename: find_instruments.py
Author: Ethan Ruddell
Description: Utility script to discover all VISA/GPIB instruments connected
             to this computer.

Run this file directly from a terminal:
    python find_instruments.py

It will print the GPIB address and identity (manufacturer, model) of every
instrument it can find.  Copy the address string (e.g. 'GPIB0::16::INSTR')
into config.yaml or the relevant module constant to use that instrument.

This is especially useful the first time you set up the lab computer or when
an instrument has been moved to a different GPIB address.
"""

import pyvisa

def main():
    print("="*60)
    print("VISA Resource Finder")
    print("="*60)
    
    try:
        rm = pyvisa.ResourceManager()
        print(f"\nVISA Library: {rm}")
        
        resources = rm.list_resources()
        
        if not resources:
            print("\nNo VISA resources found!")
            print("\nPossible issues:")
            print("  - GPIB interface not installed")
            print("  - Instrument not powered on")
            print("  - GPIB cable not connected")
            print("  - NI-VISA or similar driver not installed")
            return
        
        print(f"\nFound {len(resources)} resource(s):\n")
        
        for i, resource in enumerate(resources):
            print(f"{i+1}. {resource}")
            
            # Try to open and query each resource
            try:
                inst = rm.open_resource(resource)
                inst.timeout = 5000  # 5 second timeout
                
                # Try to get instrument ID
                try:
                    idn = inst.query("*IDN?")
                    print(f"   ID: {idn.strip()}")
                except Exception:
                    print(f"   (Could not query *IDN?)")
                
                inst.close()
                
            except Exception as e:
                print(f"   (Could not open: {e})")
            
            print()
        
        print("="*60)
        print("\nTo use an instrument, copy its address (e.g., 'GPIB0::17::INSTR')")
        print("and set GPIB_ADDRESS in pspa.py")
        print("="*60)
        
    except Exception as e:
        print(f"\nError: {e}")
        print("\nMake sure NI-VISA or PyVISA-py is installed:")
        print("  pip install pyvisa")
        print("  pip install pyvisa-py  # if NI-VISA is not installed")

if __name__ == "__main__":
    main()
