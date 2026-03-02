"""
Filename: main.py
Author: Ethan Ruddell
Date: 2026-3-2
Description: Entry point for the NRG Semiconductor Measurement Automation.

This is the file you run to start the software:
    python main.py            → launches the graphical user interface (GUI)
    python main.py --cli      → launches the command-line interface (CLI)
    python main.py -o <dir>   → set a custom output directory for data files

The GUI is the default because it is easier for day-to-day use.  The CLI
is available for scripting, remote sessions, or when a graphical display
is not available.

This file is intentionally very short.  All it does is:
  1. Parse command-line arguments (--cli, --output-dir).
  2. Set up the logging system.
  3. Hand off control to gui.py (GUI mode) or cli.py (CLI mode).
"""



# TODO: add check to confirm measurements before proceding, DO YOU KNOW WHAT YOU'RE DOING
import os
import sys
import argparse
import logging

# Add the 'utility' and 'measurement functions' folders to Python's search
# path so that "import file_management" (etc.) works without package prefixes.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'utility'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'measurement functions'))

import file_management
import logging_config

log = logging.getLogger(__name__)


def parse_args():
    """Parse CLI arguments for the measurement automation script."""
    parser = argparse.ArgumentParser(description="Automate semiconductor measurements.")
    parser.add_argument(
        "-o",
        "--output-dir",
        default=None,
        help=f"Output directory for results (default: {file_management.default_output_dir})",
    )
    parser.add_argument(
        "--cli",
        action="store_true",
        help="Launch the CLI instead of the GUI (default is GUI)",
    )
    return parser.parse_args()


def gui_execute_measurement(instrument, measurement_idx, params):
    """
    Execute a measurement from the GUI (non-interactive).

    Args:
        instrument: "PIA", "PSPA", or "LCR"
        measurement_idx: 1-indexed measurement number
        params: Dictionary of parameters from GUI entry fields
    
    Returns:
        String path to the most recently saved image file, or None if no image was generated
    """
    from gui_measurements import execute_pia_gui, execute_pspa_gui, execute_lcr_gui

    try:
        if instrument == "PIA":
            return execute_pia_gui(measurement_idx, params)
        elif instrument == "PSPA":
            return execute_pspa_gui(measurement_idx, params)
        elif instrument == "LCR":
            return execute_lcr_gui(measurement_idx, params)
        else:
            raise ValueError(f"Unknown instrument: {instrument}")
    except Exception as e:
        log.error("Measurement failed: %s", e)
        raise



def main():
    """Main entry point for the measurement automation system."""
    args = parse_args()

    # --- Logging ----------------------------------------------------------
    logging_config.setup(output_dir=file_management.default_output_dir)
    log.info("Welcome to the NRG Semiconductor Measurement Automation")

    if args.cli:
        log.info("Launching CLI mode")
        from cli import cli_main
        cli_main(output_dir=args.output_dir)
    else:
        log.info("Launching GUI mode")
        import gui
        gui.run_gui_mode(gui_execute_measurement)


if __name__ == "__main__":
    main()
