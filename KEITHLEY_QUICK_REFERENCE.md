# 3-Probe Measurement - Quick Reference Card

## 🎯 System Components

| Component | Role | GPIB Address | Key Channel |
|-----------|------|--------|----------|
| 4156C | Dual SMU parameter analyzer | GPIB0::16 | SMU1 (drain), SMU2 (source) |
| Keithley 2400 | Gate bias SourceMeter | GPIB0::17 | Output HI (gate bias) |

---

## ⚡ ENABLE SEQUENCE (Device Power-Up Order)

### ✅ CORRECT ORDER - DO THIS

```
Step 1: Enable 4156C SMU2 (Source)
        → Applies 0V at device source
        → Establishes ground reference
        
Step 2: Enable 4156C SMU1 (Drain)
        → Applies Vds at device drain
        → With source reference already present
        
Step 3: Enable Keithley 2400 (Gate)
        → Applies Vg gate bias
        → After drain/source stable
```

### ❌ WRONG ORDER - DON'T DO THIS

```
❌ Gate first (damage via source forward-bias)
❌ Random order  (undefined device state)
❌ All at once   (transient overshoots)
```

### 🛡️ WHY This Order?

**Without source reference first:**
- Source floats → undefined voltage
- Drain bias applied → forward biases source-bulk junction
- Device damaged or behaves unpredictably

**With source at 0V first:**
- Sets known reference point
- Source-bulk junction reverse-biased
- Drain bias applied safely with reference present
- Gate bias applied to stable configuration
- Device behavior is predictable and safe

---

## 🛑 SHUTDOWN SEQUENCE (Device Power-Down Order)

### ✅ CORRECT ORDER - DO THIS (REVERSE OF ENABLE)

```
Step 1: Disable Keithley 2400 (Gate)
        → Remove gate bias first
        → Optional: ramp down to 0V
        
Step 2: Disable 4156C SMU1 (Drain)
        → Remove drain bias second
        
Step 3: Disable 4156C SMU2 (Source)
        → Remove ground reference last
```

### ❌ WRONG ORDER - DON'T DO THIS

```
❌ Source first    (leaves drain/gate floating - damage!)
❌ Random order    (undefined shutdown state)
❌ All at once     (parasitic transient currents)
```

### 🛡️ WHY Reverse Order?

**Without draining gate first:**
- Gate still biased when drain removed
- Drain node floats
- Parasitic current paths activate
- Source becomes undefined
- Device enters invalid state

**With proper reverse sequence:**
- Gate bias removed → no charge injection
- Drain bias removed → clean disconnect
- Source reference removed last
- Device left in safe, defined state
- Ready for next measurement cycle

---

## 🔌 Hardware Checklist

Before every measurement:

```
□ Verify common ground connection
  - Keithley 2400 LO tied to 4156C common
  - NOT left floating
  
□ Check device connections
  - Gate → Keithley (NOT 4156C!)
  - Drain → 4156C SMU1 (CH1)
  - Source → 4156C SMU2 (CH2)
  - Body/Sub → 4156C common
  
□ Both instruments powered on
  - 4156C: GPIB0::16
  - Keithley: GPIB0::17
  
□ Compliance settings reasonable
  - All channels: same limit
  - Matches expected current
```

---

## 📊 Typical Measurement Flow

```
1. INITIALIZATION
   ├─ verify_ground_connection()
   ├─ connect_4156c()
   ├─ connect_keithley_2400()
   ├─ initialize_4156c(Vds=5V, I_comp=0.1A)
   ├─ initialize_keithley_2400(Vg_start=0V, I_comp=0.1A)
   └─ [All outputs OFF - safe state]

2. ENABLE SEQUENCE (in order!)
   ├─ enable_4156c_smu2()          [source to 0V]
   ├─ enable_4156c_smu1(settle=1s) [apply Vds]
   └─ enable_keithley_2400()       [apply Vg]
   
3. MEASUREMENT LOOP
   ├─ For each gate voltage step:
   │  ├─ Set Keithley to next Vg
   │  ├─ Wait for settling
   │  └─ Measure drain current
   └─ Store (Vg, Id) pairs

4. SHUTDOWN SEQUENCE (in reverse order!)
   ├─ disable_keithley_2400()   [remove gate]
   ├─ disable_4156c_smu1()      [remove drain]
   └─ disable_4156c_smu2()      [remove source]
```

---

## ⚙️ Code Usage Example

```python
import keithley

# Single function call - handles everything!
result = keithley.measure_transistor_vg_sweep(
    gpib_4156c="GPIB0::16::INSTR",
    gpib_keithley="GPIB0::17::INSTR",
    vds=5.0,            # Apply 5V to drain
    vg_start=0.0,       # Gate sweep 0V to +5V
    vg_stop=5.0,
    vg_step=0.1,        # 0.1V steps (51 points)
    i_compliance=0.1,   # 100mA limit (all channels)
    vds_settle_s=1.0,   # Wait 1s after enabling drain
    vg_settle_s=0.05    # Wait 50ms per gate step
)

# Extract results
if result['status'] == 'success':
    vg = result['vg']   # Gate voltages (array)
    id = result['id']   # Drain currents (array)
    print(f"Measured {len(vg)} points successfully")
```

---

## 🚨 Emergency Procedures

### If measurement fails mid-sweep:

1. **DO NOT** attempt manual commands - may worsen state
2. Let the code execute its **emergency shutdown sequence**:
   - Automatically disables gate first
   - Then drain
   - Then source
3. Check hardware connections
4. Verify no smoke or burning smell
5. Wait 30 seconds before retry

### If device appears damaged:

1. Power down both instruments immediately
2. Visually inspect device for:
   - Burn marks
   - Discoloration
   - Component damage
3. Check connections for shorts/reversal
4. If device is sacrificial (test die): replace and proceed
5. If device is valuable: debug connection first

---

## 📋 Key Parameters Explained

| Parameter | Units | Purpose | Typical Range | Notes |
|-----------|-------|---------|----------------|-------|
| `Vds` | V | Drain-source bias | 0.1 - 50 | Set by application |
| `Vg_start` | V | Gate sweep minimum | -5 to +5 | Device dependent |
| `Vg_stop` | V | Gate sweep maximum | -5 to +5 | Device dependent |
| `Vg_step` | V | Gate increment | 0.01 - 1.0 | Smaller = slower but smoother |
| `i_compliance` | A | Current limit (safety) | 1mA - 1A | Protects device |
| `vds_settle_s` | s | After drain enable | 0.5 - 5.0 | Longer for reactive loads |
| `vg_settle_s` | s | After gate step | 0.01 - 0.5 | Longer for parasitic RC |

---

## 🔍 Troubleshooting Quick Links

| Symptom | Likely Cause | Fix |
|---------|-------------|-----|
| "Ground not verified" | User declined confirmation | Run again, verify connections |
| No current (zero reading) | Device off / disconnected | Check wiring, retry |
| Noisy data | Electromagnetic interference | Increase settling times |
| Measurement stops | VISA bus error | Restart instruments, reconnect |
| Zero Vds applied | Compliance too low | Increase both to 0.1-1.0 A |
| Device damaged | Reversed connections | Check gate NOT on 4156C |

---

## 📚 Further Reading

- Full implementation: `keithley.py`
- Complete documentation: `KEITHLEY_3PROBE_IMPLEMENTATION.md`
- Example usage: `example_3probe_transistor_sweep.py`
- Legacy reference: `LEGACY/Need ReRAM Mark 4.py`

---

**Last Updated:** May 28, 2026  
**Status:** Production Ready ✅
