import pyvisa
import time
import numpy as np
import matplotlib.pyplot as plt
import csv
import os
from datetime import datetime

# =======================
# VISA SETUP
# =======================
rm = pyvisa.ResourceManager()
inst = rm.open_resource("GPIB0::17::INSTR")
inst.timeout = 10000

print(inst.query("*IDN?"))

# =======================
# USER PARAMETERS
# =======================
V_form   = 0
I_comp   = 0.001
Vpos     = 5
Vneg     = -5
step     = 0.1
cycles   = 1
settle  = 0.05

# =======================
# OUTPUT
# =======================
outdir = "output"
os.makedirs(outdir, exist_ok=True)
ts = datetime.now().strftime("%Y%m%d_%H%M%S")

# =======================
# CONFIGURE SMU (CRITICAL)
# =======================
inst.write("*RST")
inst.write("ROUT:TERM REAR")
inst.write("SOUR:FUNC VOLT")
inst.write("SOUR:VOLT:MODE FIXED")
inst.write("SENS:FUNC 'CURR:DC'")
inst.write("SENS:CURR:RANG:AUTO ON")
inst.write("FORM:ELEM CURR")           # <<< THIS FIXES EVERYTHING
inst.write(f"SOUR:VOLT:ILIM {I_comp}")

# =======================
# MEASUREMENT FUNCTION
# =======================
def sweep(voltage_list):
    V, I = [], []
    inst.write("OUTP ON")
    for v in voltage_list:
        inst.write(f"SOUR:VOLT {v}")
        time.sleep(settle)
        i = float(inst.query("READ?"))
        V.append(v)
        I.append(i)
    inst.write("OUTP OFF")
    return np.array(V), np.array(I)

# =======================
# FORMING
# =======================
forming_v = np.arange(0, V_form + step, step)
Vf, If = sweep(forming_v)

# Save CSV
with open(f"{outdir}/Forming_{ts}.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["Voltage (V)", "Current (A)"])
    w.writerows(zip(Vf, If))

# =======================
# NORMAL CYCLES
# =======================
all_cycles = []

inst.write(f"SOUR:VOLT:RANG {max(abs(Vpos), abs(Vneg))}")

for c in range(cycles):
    sweep_v = np.concatenate([
        np.arange(0, Vpos + step, step),
        np.arange(Vpos, 0, -step),
        np.arange(0, Vneg - step, -step),
        np.arange(Vneg, step, step)
    ])
    Vc, Ic = sweep(sweep_v)
    all_cycles.append((Vc, Ic))

    with open(f"{outdir}/Cycle_{c+1}_{ts}.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Voltage (V)", "Current (A)"])
        w.writerows(zip(Vc, Ic))

# =======================
# PLOTS — LINEAR
# =======================
plt.figure(figsize=(8,6))
plt.plot(Vf, If, '--', label="Forming")
for i,(V,I) in enumerate(all_cycles):
    plt.plot(V, I, label=f"Cycle {i+1}")
plt.xlabel("Voltage (V)")
plt.ylabel("Current (A)")
plt.title("I–V (Linear)")
plt.legend()
plt.grid(True)
plt.savefig(f"{outdir}/IV_Linear_{ts}.png", dpi=300)
plt.show()

# =======================
# PLOTS — LOG |I|
# =======================
plt.figure(figsize=(8,6))
plt.semilogy(Vf, np.abs(If), '--', label="Forming")
for i,(V,I) in enumerate(all_cycles):
    plt.semilogy(V, np.abs(I), label=f"Cycle {i+1}")
plt.xlabel("Voltage (V)")
plt.ylabel("|Current| (A)")
plt.title("|I|–V (Log)")
plt.legend()
plt.grid(True, which="both")
plt.savefig(f"{outdir}/IV_Log_{ts}.png", dpi=300)
plt.show()

# =======================
# CLEANUP
# =======================
inst.write("OUTP OFF")
inst.close()
rm.close()
