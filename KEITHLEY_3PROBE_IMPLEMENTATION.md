# 3-Probe Transistor Measurement Implementation

## Overview

This implementation provides remote measurement functions for Keithley 2400 and 4156C semiconductor parameter analyzers, specifically designed for **3-probe transistor characterization** with proper biasing order and shutdown sequence.

### Instruments Used

1. **4156C Semiconductor Parameter Analyzer** (HP/Agilent)
   - SMU1 (Channel 1): Drain probe - applies `Vds`, measures `Id`
   - SMU2 (Channel 2): Source probe - always at 0V (ground reference)

2. **Keithley 2400 SourceMeter**
   - Gate probe - steps through `Vg` values

### Key Innovation: Proper Enable/Disable Sequencing

The implementation enforces a **critical sequence** based on semiconductor physics:

#### Enable Sequence (Powers up device safely)
1. **SMU2 (Source) First** → Establishes 0V ground reference at device source
2. **SMU1 (Drain) Second** → Applies Vds with source reference already present
3. **2400 (Gate) Last** → Applies Vg bias after drain/source are stable

**Why this order?**
- Forward biasing prevention: Source at 0V first prevents source-bulk diodes from forward biasing
- Reference stability: Establishes measurement reference before loading circuit
- Junction protection: Prevents undefined states between SMUs

#### Disable Sequence (REVERSE order - essential!)
1. **2400 (Gate) First** → Remove gate bias before drain
2. **SMU1 (Drain) Second** → Remove drain bias before source reference
3. **SMU2 (Source) Last** → Remove ground reference last

**Why reverse order?**
- Avoids floating drain/gate when source is still referenced
- Prevents forward-biasing parasitic junctions
- Ensures clean shutdown and device is left in defined state
- Critical for back-to-back measurements without device damage

---

## Module Structure

### Connection Functions

```python
connect_4156c(resource_name=None)
connect_keithley_2400(resource_name=None)
disconnect_instrument(inst, name="Instrument")
```

Open VISA sessions to instruments and safely close them.

### Initialization Functions

```python
verify_ground_connection()
initialize_4156c(inst_4156c, vds_target=5.0, i_compliance=0.1)
initialize_keithley_2400(inst_keithley, vg_start=0.0, i_compliance=0.1)
```

Initialize instruments to known states with matching compliance settings.
**Note:** Outputs remain OFF until enable sequence is called.

### Enable Sequence Functions

```python
enable_4156c_smu2(inst_4156c)              # Step 1: Source reference
enable_4156c_smu1(inst_4156c, settle_time_s=1.0)  # Step 2: Apply Vds
enable_keithley_2400(inst_keithley)        # Step 3: Apply Vg
enable_all_3probe(inst_4156c, inst_keithley, vds_settle_s=1.0)
```

Execute enable sequence in correct order with settling times.

### Measurement Function

```python
measure_transistor_id(inst_4156c)
```

Trigger spot measurement on 4156C SMU1 and return drain current.

### Disable Sequence Functions

```python
disable_keithley_2400(inst_keithley, ramp_to_zero=True, ramp_steps=10)
disable_4156c_smu1(inst_4156c)
disable_4156c_smu2(inst_4156c)
disable_all_3probe(inst_keithley, inst_4156c)
```

Execute shutdown in reverse order with optional gate ramp-down.

### Main Measurement Function

```python
result = measure_transistor_vg_sweep(
    gpib_4156c="GPIB0::16::INSTR",
    gpib_keithley="GPIB0::17::INSTR",
    vds=5.0,
    vg_start=0.0,
    vg_stop=5.0,
    vg_step=0.1,
    i_compliance=0.1,
    vds_settle_s=1.0,
    vg_settle_s=0.05
)
```

**Parameters:**
- `gpib_4156c`: GPIB address of 4156C (default: "GPIB0::16::INSTR")
- `gpib_keithley`: GPIB address of Keithley 2400 (default: "GPIB0::17::INSTR")
- `vds`: Drain-source voltage (V)
- `vg_start`, `vg_stop`, `vg_step`: Gate voltage sweep range
- `i_compliance`: Current limit for all channels (A)
- `vds_settle_s`: Settling time after enabling drain bias (s)
- `vg_settle_s`: Settling time per gate step (s)

**Returns:**
```python
{
    'vg': np.array([...]),      # Gate voltages
    'id': np.array([...]),      # Drain currents
    'status': 'success'|'failed'
}
```

---

## Hardware Setup

### Connections

```
Keithley 2400 (GPIB0::17):
  Output HI  ─────→ Device Gate (G)
  Output LO  ─────→ 4156C Circuit Common (CRITICAL!)

4156C (GPIB0::16):
  SMU2 (CH2) ─────→ Device Source (S)
  SMU1 (CH1) ─────→ Device Drain (D)
  Common     ←─────┘ (shared with Keithley LO)
  
DUT (semiconductor device):
  Gate ←────── Keithley output
  Drain ←───── 4156C SMU1
  Source ←──── 4156C SMU2
  Body/Sub ← (optional, to common)
```

### Critical Hardware Check

Before running measurements:

1. **Verify Common Ground Connection**
   - Keithley 2400 LO must be connected to 4156C circuit common
   - This is YOUR responsibility - the software will prompt for confirmation
   - Incorrect grounding can damage both instruments and DUT

2. **Check Device Connections**
   - Gate → Keithley (not 4156C!)
   - Drain → 4156C SMU1
   - Source → 4156C SMU2
   - Body/Substrate → Circuit common (if used)

---

## Safety Features

### Built-in Protections

1. **Ground Connection Verification**
   - User must confirm proper connections before measurement starts
   - Prevents accidental data collection with incorrect wiring

2. **Proper Enable/Disable Sequence**
   - Enforces correct order to prevent junction forward-biasing
   - Automatic settling delays between steps

3. **Current Compliance**
   - All channels limited to prevent device damage
   - Uniform compliance across all measurement points

4. **Emergency Shutdown**
   - If measurement fails, automatic emergency disable sequence
   - Attempts to safely disable all outputs in reverse order

5. **Settling Times**
   - Configurable delays allow device transient effects to settle
   - Default values suitable for typical MOSFETs and BJTs

### Compliance Settings

- **Typical devices**: 0.1 A (safe for most devices)
- **Low-power devices**: 1 mA to 10 mA (adjust as needed)
- **Power devices**: 1 A+ (verify device specs first)

All channels must use the SAME compliance value.

---

## Usage Example

```python
import keithley

# Run measurement
result = keithley.measure_transistor_vg_sweep(
    gpib_4156c="GPIB0::16::INSTR",
    gpib_keithley="GPIB0::17::INSTR",
    vds=5.0,           # 5V drain-source
    vg_start=0.0,      # Gate sweep 0V to 5V
    vg_stop=5.0,
    vg_step=0.1,       # 0.1V steps
    i_compliance=0.1,  # 100mA limit
    vds_settle_s=1.0,  # 1s VDS settling
    vg_settle_s=0.05   # 50ms per step
)

# Check results
if result['status'] == 'success':
    vg = result['vg']
    id = result['id']
    print(f"Measured {len(vg)} points")
    print(f"Id range: {min(id):.3e} to {max(id):.3e} A")
```

---

## Technical Details

### 4156C FLEX Commands Used

- `US` – Switch to FLEX command mode
- `FMT 1` – ASCII format with header
- `CL` – Clear (disable) all channels
- `*RST` – Full reset
- `*CLS` – Clear error queue
- `CH1/CH2` – Select channel
- `FUNC VOLT` – Set voltage source mode
- `LIMI {value}` – Set current compliance
- `DV {value}` – Set/apply voltage
- `CN` – Enable (connect) channel
- `CL` – Disable (disconnect) channel
- `TI?` – Trigger and read current (spot measurement)

### Keithley 2400 SCPI Commands Used

- `:SOUR:FUNC VOLT` – Voltage source mode
- `:SOUR:VOLT:MODE FIXED` – Fixed voltage output
- `:SENS:FUNC 'CURR:DC'` – Current measurement
- `:SENS:CURR:RANG:AUTO ON` – Auto current range
- `:FORM:ELEM CURR` – Output current only
- `:SOUR:VOLT:ILIM` – Current compliance
- `:SOUR:VOLT` – Set voltage value
- `:OUTP ON/OFF` – Enable/disable output
- `:*RST` – Reset
- `:*CLS` – Clear errors

---

## Troubleshooting

### "Keithley GPIB address not configured"
- Check GPIB address in measurement call
- Verify instrument is connected and powered on
- Run `keithley.list_visa_resources()` to see available devices

### Measurement doesn't start after ground verification
- Ensure both instruments show no error condition
- Try manual reset: `*RST` on both instruments
- Check compliance values aren't too low (try 0.1 A)

### Zero current reading
- Verify DUT is correctly connected
- Check device isn't in off-state or broken
- Ensure compliance limit can supply needed current

### Erratic measurements
- Increase settling times (`vds_settle_s`, `vg_settle_s`)
- Reduce gate step size (`vg_step`)
- Check for electromagnetic noise near instruments

### Device appears damaged after measurement
- **Likely cause**: Improper enable/disable sequence or reversed connections
- **Prevention**: Always follow the documented wiring and use provided functions
- Use example script as template for new measurements

---

## Performance Considerations

### Measurement Speed
- Each gate step: ~50 ms (default) + instrument overhead
- 51-point sweep: ~5-10 seconds typical
- Adjust `vg_settle_s` down for faster sweeps on stable devices

### Accuracy
- Current resolution limited by instrument specs (typically nA-pA range)
- Use compliance settings matched to expected current
- Auto-range on Keithley ensures best measurement accuracy

### Noise Reduction
- Increase settling times on noisy supplies
- Wrap cables in shielding where possible
- Ensure proper grounding (star point at circuit common)

---

## Reference: Manual Entry for "Need ReRAM Mark 4.py"

The legacy file `LEGACY/Need ReRAM Mark 4.py` shows a single-channel 2400 sweep approach.
The new implementation extends this to proper 3-probe transistor characterization with:

✅ Correct enable/disable sequence (prevents damage)
✅ Automatic ground verification (safety)
✅ Proper 4156C configuration (dual-SMU support)
✅ Emergency shutdown (fault tolerance)
✅ Comprehensive error handling (reliability)

---

## Version History

- **2026-05-28**: Initial 3-probe implementation with proper sequencing
  - Implements exact user-specified enable/disable sequence
  - Adds comprehensive error handling and safety checks
  - Full documentation and example script
