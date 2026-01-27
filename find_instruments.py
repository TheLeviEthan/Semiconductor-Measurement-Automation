"""
Utility script to find and test VISA instruments.
Run this to discover the GPIB address of your PSPA.
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
                except:
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
