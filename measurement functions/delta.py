import re
import time
import math
from collections import deque
import pyvisa


CHAMBER_ADDR = "GPIB0::27::INSTR"

HIGH_TEMP_C = 150.0
STEP_C = 5.0

TEMP_TOL_C = 0.5
POLL_INTERVAL_S = 5

# Stabilization detection
TREND_WINDOW_POINTS = 6          # ~30 s window at 5 s polling
MAX_ABS_SLOPE_C_PER_MIN = 0.05   # "flat enough"


class DeltaChamber:
    def __init__(self, address=CHAMBER_ADDR):
        self.rm = pyvisa.ResourceManager()
        self.inst = self.rm.open_resource(address)
        self.inst.timeout = 5000
        self.inst.write_termination = "\r"
        self.inst.read_termination = "\n"

    @staticmethod
    def _parse_float(resp: str) -> float:
        m = re.search(r"[-+]?\d+(?:\.\d+)?", resp)
        if not m:
            raise ValueError(f"Could not parse float from response: {resp!r}")
        return float(m.group(0))

    def get_temperature(self) -> float:
        return self._parse_float(self.inst.query("TEMP?"))

    def get_setpoint(self) -> float:
        return self._parse_float(self.inst.query("SETP?"))

    def set_setpoint(self, temp_c: float):
        self.inst.write(f"SETP {temp_c:.1f}")

    def close(self):
        try:
            self.inst.close()
        finally:
            self.rm.close()


def frange_inclusive(start, stop, step):
    vals = []
    x = start
    while x <= stop + 1e-9:
        vals.append(round(x, 1))
        x += step
    return vals


def next_step_at_or_above(temp_c, step_c=5.0):
    return math.ceil(temp_c / step_c) * step_c


def compute_slope_c_per_min(samples):
    if len(samples) < 2:
        return None

    t0, y0 = samples[0]
    t1, y1 = samples[-1]
    dt = t1 - t0
    if dt <= 0:
        return None

    return (y1 - y0) / dt * 60.0


def wait_until_stable(chamber, target_c):
    samples = deque(maxlen=TREND_WINDOW_POINTS)

    while True:
        now = time.time()
        temp = chamber.get_temperature()
        err = abs(temp - target_c)

        samples.append((now, temp))
        slope = compute_slope_c_per_min(samples)

        slope_str = "N/A" if slope is None else f"{slope:+.3f} C/min"
        print(f"Temp={temp:.2f} C | Target={target_c:.2f} C | Err={err:.2f} C | Slope={slope_str}")

        in_band = err <= TEMP_TOL_C
        flat = slope is not None and abs(slope) <= MAX_ABS_SLOPE_C_PER_MIN
        enough = len(samples) >= TREND_WINDOW_POINTS

        if in_band and flat and enough:
            print(">>> Temperature stabilized <<<")
            return temp

        time.sleep(POLL_INTERVAL_S)


def measure_with_agilent(target_c, actual_c):
    print(f">>> MEASURE at target={target_c:.1f} C, actual={actual_c:.2f} C <<<")
    # Replace later with real Agilent code


def main():
    chamber = DeltaChamber(CHAMBER_ADDR)

    try:
        current_temp = chamber.get_temperature()
        current_setpoint = chamber.get_setpoint()

        print("Connected to chamber")
        print(f"Current temp    : {current_temp:.2f} C")
        print(f"Current setpoint: {current_setpoint:.2f} C")

        start_temp = next_step_at_or_above(current_temp, STEP_C)

        if start_temp > HIGH_TEMP_C:
            print("Current temperature is already above the requested high temperature.")
            return

        temps = frange_inclusive(start_temp, HIGH_TEMP_C, STEP_C)

        print("\nHeating sequence:")
        print(temps)

        for i, target in enumerate(temps, start=1):
            print("\n" + "=" * 60)
            print(f"Step {i}/{len(temps)} -> {target:.1f} C")

            print(f"Setting chamber to {target:.1f} C")
            chamber.set_setpoint(target)

            # Confirm the setpoint actually changed
            time.sleep(1)
            confirmed_setpoint = chamber.get_setpoint()
            print(f"Controller setpoint now: {confirmed_setpoint:.1f} C")

            stabilized_temp = wait_until_stable(chamber, target)
            measure_with_agilent(target, stabilized_temp)

        print("\nDONE")

    finally:
        chamber.close()


if __name__ == "__main__":
    main()