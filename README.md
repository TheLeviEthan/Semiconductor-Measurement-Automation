# Semiconductor-Measurement-Automation

Automated measurements for PIA (4294A), PSPA (4155C/4156C), LCR (E4980A), and a Cryo-Con Model 32/32B temperature controller.

## Cryogenic Sweep Behavior

- The sweep starts at the current controller temperature and ramps down to the target temperature.
- The controller cools quickly on its own; the warmer is used to slow the cool-down to the requested ramp rate.
- Ramping up is not supported in cryo mode.
- The system holds temperature at each point while measurements run, then resumes the ramp to the next point.
- Temperature points are generated from the current temperature to the target using the measurement interval.

## CLI Usage

1. Run `python main.py`.
2. Select `0. CRYO`.
3. Enter target temperature, ramp rate, and measurement interval.
4. Queue measurements and start the sweep.

## GUI Usage

1. Run `python main.py --gui`.
2. Enable the cryo sweep checkbox.
3. Enter target temperature, ramp rate, and measurement interval.
4. Queue measurements and start the sweep.