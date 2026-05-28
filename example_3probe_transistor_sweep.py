"""
Example: 3-Probe Transistor Measurement with Keithley 2400 and 4156C

This demonstrates the complete workflow for measuring transistor Id-Vg characteristics
using the newly implemented keithley module with proper enable/disable sequencing.

HARDWARE SETUP:
===============
  Keithley 2400 (at GPIB0::17):
    - Output HI → Device Gate (G)
    - Output LO → 4156C circuit common (CRITICAL!)
  
  4156C (at GPIB0::16):
    - SMU2 (CH2) → Device Source (S) at 0V
    - SMU1 (CH1) → Device Drain (D) [measures Id with Vds applied]
    - Circuit common → Device substrate/body (if used)

MEASUREMENT PARAMETERS:
=======================
  Vds (drain-source):   5.0 V  (applied by 4156C SMU1)
  Vg sweep:             0 V to +5 V in 0.1 V steps
  Current compliance:   0.1 A (all channels)
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Add the measurement functions module to path
meas_func_dir = Path(__file__).parent / "measurement functions"
sys.path.insert(0, str(meas_func_dir))

import keithley


def main():
    """Run the 3-probe transistor measurement."""
    
    print("\n" + "="*70)
    print("3-PROBE TRANSISTOR MEASUREMENT")
    print("Keithley 2400 Gate + 4156C Drain/Source")
    print("="*70)
    
    # Measurement parameters
    VDS = 5.0          # Drain-source voltage (V)
    VG_START = 0.0     # Starting gate voltage (V)
    VG_STOP = 5.0      # Ending gate voltage (V)
    VG_STEP = 0.1      # Gate voltage step (V)
    I_COMPLIANCE = 0.1 # Current limit (A)
    
    # Connection parameters (update if your instruments are at different addresses)
    GPIB_4156C = "GPIB0::16::INSTR"
    GPIB_KEITHLEY = "GPIB0::17::INSTR"
    
    # Settling times
    VDS_SETTLE_S = 1.0   # Vds settling after enable (s)
    VG_SETTLE_S = 0.05   # Per-step gate settling (s)
    
    # Run the measurement
    result = keithley.measure_transistor_vg_sweep(
        gpib_4156c=GPIB_4156C,
        gpib_keithley=GPIB_KEITHLEY,
        vds=VDS,
        vg_start=VG_START,
        vg_stop=VG_STOP,
        vg_step=VG_STEP,
        i_compliance=I_COMPLIANCE,
        vds_settle_s=VDS_SETTLE_S,
        vg_settle_s=VG_SETTLE_S
    )
    
    # Check if measurement succeeded
    if result['status'] != 'success':
        print("\nMeasurement failed!")
        return
    
    # Extract data
    vg_data = result['vg']
    id_data = result['id']
    
    # Display summary
    print(f"\nMeasurement Summary:")
    print(f"  Points collected: {len(vg_data)}")
    print(f"  Vg range: {vg_data[0]:.2f}V to {vg_data[-1]:.2f}V")
    print(f"  Id range: {np.min(id_data):.3e}A to {np.max(id_data):.3e}A")
    print(f"  Vds applied: {VDS:.1f}V")
    
    # Create summary plots
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    
    # Linear I-V plot
    axes[0].plot(vg_data, id_data, 'b-o', linewidth=1.5, markersize=4)
    axes[0].set_xlabel('Gate Voltage (V)')
    axes[0].set_ylabel('Drain Current (A)')
    axes[0].set_title(f'Id-Vg (Linear) @ Vds = {VDS}V')
    axes[0].grid(True, alpha=0.3)
    
    # Log plot
    axes[1].semilogy(vg_data, np.abs(id_data), 'r-o', linewidth=1.5, markersize=4)
    axes[1].set_xlabel('Gate Voltage (V)')
    axes[1].set_ylabel('|Drain Current| (A)')
    axes[1].set_title(f'|Id|-Vg (Log) @ Vds = {VDS}V')
    axes[1].grid(True, alpha=0.3, which='both')
    
    plt.tight_layout()
    plt.show()
    
    print("\nMeasurement complete!")


if __name__ == "__main__":
    main()
