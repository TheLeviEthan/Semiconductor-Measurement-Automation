"""
Filename: main.py
Author: Ethan Ruddell
Date: 2026-2-27
Description: Entry point for the NRG Semiconductor Measurement Automation.
Launches the GUI by default. Pass --cli to use the command-line interface.
"""

import os
import sys
import argparse
import logging

# Add directories to path
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
    """
    from gui_measurements import execute_pia_gui, execute_pspa_gui, execute_lcr_gui

    try:
        if instrument == "PIA":
            execute_pia_gui(measurement_idx, params)
        elif instrument == "PSPA":
            execute_pspa_gui(measurement_idx, params)
        elif instrument == "LCR":
            execute_lcr_gui(measurement_idx, params)
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
