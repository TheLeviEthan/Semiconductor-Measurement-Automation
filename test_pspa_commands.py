"""
Test script to verify PSPA SCPI commands.
Use this to test which command syntax works with your PSPA.
"""

import pyvisa
import time

# Update this with your PSPA's GPIB address
GPIB_ADDRESS = "GPIB0::16::INSTR"  # Change this to match your instrument

def test_commands():
    print("="*60)
    print("PSPA Command Test")
    print("="*60)
    
    try:
        rm = pyvisa.ResourceManager()
        print(f"\nConnecting to {GPIB_ADDRESS}...")
        pspa = rm.open_resource(GPIB_ADDRESS)
        pspa.timeout = 10000  # 10 second timeout
        
        # Test 1: Get instrument ID
        print("\n1. Testing *IDN? ...")
        try:
            idn = pspa.query("*IDN?")
            print(f"   Success: {idn.strip()}")
        except Exception as e:
            print(f"   Failed: {e}")
            return
        
        # Test 2: Clear errors
        print("\n2. Clearing errors with *CLS...")
        try:
            pspa.write("*CLS")
            print("   Success")
        except Exception as e:
            print(f"   Failed: {e}")
        
        # Test 3: Check for errors
        print("\n3. Checking for errors...")
        try:
            error = pspa.query(":SYST:ERR?")
            print(f"   Error status: {error.strip()}")
        except Exception as e:
            print(f"   Failed: {e}")
        
        # Test 4: Try different measurement command formats
        channel = 1
        print(f"\n4. Testing measurement commands on channel {channel}...")
        
        # First configure the channel
        print("   Configuring channel...")
        try:
            pspa.write(f":SOUR{channel}:FUNC:MODE VOLT")
            pspa.write(f":SENS{channel}:CURR:PROT 0.1")
            pspa.write(f":SOUR{channel}:VOLT 0")
            print("   Configuration success")
        except Exception as e:
            print(f"   Configuration failed: {e}")
        
        # Test different read formats
        test_formats = [
            f":READ? (@{channel})",
            f":READ{channel}?",
            f":MEAS:CURR? (@{channel})",
            f":MEAS{channel}:CURR?",
            f":FETC? (@{channel})",
        ]
        
        for cmd in test_formats:
            print(f"\n   Testing: {cmd}")
            try:
                result = pspa.query(cmd)
                print(f"   Success: {result.strip()}")
            except Exception as e:
                print(f"   Failed: {e}")
            
            # Check for errors after each command
            try:
                error = pspa.query(":SYST:ERR?")
                if not error.startswith('0,'):
                    print(f"   Instrument error: {error.strip()}")
                    pspa.write("*CLS")  # Clear the error
            except:
                pass
            
            time.sleep(0.1)
        
        print("\n" + "="*60)
        print("Test complete!")
        print("="*60)
        
        pspa.close()
        
    except Exception as e:
        print(f"\nError: {e}")

if __name__ == "__main__":
    print("\nIMPORTANT: Update GPIB_ADDRESS at the top of this file first!")
    input("Press Enter to continue or Ctrl+C to exit...")
    test_commands()
