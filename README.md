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

## USB Switchover Box (Multiplexer)

- The GUI includes a **USB Switchbox** section under instrument selection.
- Enable **automatic routing** to switch the DUT path when you select `PIA`, `PSPA`, or `LCR`.
- The selected route is also enforced before each measurement (including each queued cryo measurement).

### Configuration

Edit `config.yaml` under `usb_switchbox`:

- `port`: your USB serial port (placeholder by default)
- `baudrate`, `timeout_s`, `write_terminator`
- `command_map`: per-instrument route commands from your switchbox manual

For Model 7197 (`M7197_5073-01RS.pdf`), default protocol assumptions are now:

- RS-422 serial through your USB-to-RS422 adapter (`COMx`)
- 9600 baud, 8 data bits, no parity, 1 stop bit, no flow control
- Commands are ASCII control characters (for example `CTRL_A`, `CTRL_B`, etc.)
- No Enter/newline terminator at end of command

Supported command value formats in `command_map`:

- `CTRL_A` ... `CTRL_Z`
- `HEX:01` style byte value
- `TEXT:...` plain ASCII command text

### Dependency

Install pyserial for USB serial communication:

`pip install pyserial`