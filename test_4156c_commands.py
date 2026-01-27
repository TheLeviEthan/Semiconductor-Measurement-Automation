"""
Test script for Agilent 4156C SCPI commands.
This script verifies the correct command syntax for your specific instrument.
"""

import pyvisa
import time

# UPDATE THIS with your actual GPIB address
GPIB_ADDRESS = ""  # e.g., "GPIB0::17::INSTR"

def test_4156c_commands():
    """Test SCPI commands for Agilent 4156C Precision Semiconductor Parameter Analyzer."""
    
    if not GPIB_ADDRESS:
        print("ERROR: Please set GPIB_ADDRESS in this script first!")
        print("Run find_instruments.py to find your 4156C's address.")
        return
    
    try:
        rm = pyvisa.ResourceManager('@py')
        print(f"Connecting to {GPIB_ADDRESS}...")
        pspa = rm.open_resource(GPIB_ADDRESS)
        pspa.timeout = 30000
        
        # Identify instrument
        idn = pspa.query("*IDN?")
        print(f"Connected to: {idn}")
        
        # Clear any existing errors
        pspa.write("*CLS")
        time.sleep(0.1)
        
        print("\n" + "="*60)
        print("Testing Agilent 4156C SCPI Commands")
        print("="*60)
        
        # Test 1: Configure SMU channel
        print("\n--- Test 1: Configure SMU Channel ---")
        test_commands = [
            (":PAGE:CHAN:SMU1:FUNC VOLT", "Set SMU1 to voltage source"),
            (":PAGE:CHAN:SMU1:ICOMP 0.1", "Set current compliance to 0.1A"),
        ]
        
        for cmd, desc in test_commands:
            print(f"\nCommand: {cmd}")
            print(f"Description: {desc}")
            pspa.write(cmd)
            time.sleep(0.1)
            error = pspa.query(":SYST:ERR?")
            print(f"Error response: {error}")
            
        # Test 2: Set voltage
        print("\n--- Test 2: Set Voltage ---")
        test_commands = [
            (":PAGE:CHAN:SMU1:VOLT 0", "Set SMU1 to 0V"),
        ]
        
        for cmd, desc in test_commands:
            print(f"\nCommand: {cmd}")
            print(f"Description: {desc}")
            pspa.write(cmd)
            time.sleep(0.1)
            error = pspa.query(":SYST:ERR?")
            print(f"Error response: {error}")
            
        # Test 3: Enable output
        print("\n--- Test 3: Enable/Disable Output ---")
        test_commands = [
            (":PAGE:CHAN:SMU1:STATE ON", "Enable SMU1 output"),
            (":PAGE:CHAN:SMU1:STATE OFF", "Disable SMU1 output"),
        ]
        
        for cmd, desc in test_commands:
            print(f"\nCommand: {cmd}")
            print(f"Description: {desc}")
            pspa.write(cmd)
            time.sleep(0.1)
            error = pspa.query(":SYST:ERR?")
            print(f"Error response: {error}")
            
        # Test 4: Measurement commands
        print("\n--- Test 4: Measurement Commands ---")
        test_commands = [
            (":MEAS:CURR? SMU1", "Measure current on SMU1"),
            (":MEAS:VOLT? SMU1", "Measure voltage on SMU1"),
        ]
        
        for cmd, desc in test_commands:
            print(f"\nCommand: {cmd}")
            print(f"Description: {desc}")
            try:
                result = pspa.query(cmd)
                print(f"Result: {result}")
            except Exception as e:
                print(f"Exception: {e}")
            time.sleep(0.1)
            error = pspa.query(":SYST:ERR?")
            print(f"Error response: {error}")
            
        print("\n" + "="*60)
        print("Testing Complete!")
        print("="*60)
        print("\nIf all error responses show '0,\"No error\"', the commands are working correctly.")
        
        # Clean up
        pspa.close()
        print("\nDisconnected from instrument.")
        
    except Exception as e:
        print(f"\nERROR: {e}")
        print("\nMake sure:")
        print("1. The 4156C is powered on and connected via GPIB")
        print("2. NI-VISA or compatible GPIB driver is installed")
        print("3. GPIB address is correct (check instrument settings)")
        
if __name__ == "__main__":
    test_4156c_commands()
