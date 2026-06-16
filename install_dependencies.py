"""Install the Python dependencies required by Semiconductor-Measurement-Automation.

Run this script from the project root:

    python install_dependencies.py

It installs the packages listed in requirements.txt using the current Python
interpreter.  For GPIB instruments on Windows, NI-VISA (or a compatible VISA
implementation) is still required separately; this script installs pyvisa and
pyvisa-py so the project can use either backend when available.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent
REQUIREMENTS_FILE = PROJECT_ROOT / "requirements.txt"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Install the Python dependencies used by this project."
    )
    parser.add_argument(
        "--upgrade-pip",
        action="store_true",
        help="Upgrade pip before installing project dependencies.",
    )
    parser.add_argument(
        "--no-input",
        action="store_true",
        help="Pass --no-input to pip for unattended installs.",
    )
    return parser.parse_args()


def run_command(command: list[str]) -> None:
    print("Running:", " ".join(command))
    subprocess.run(command, check=True)


def main() -> int:
    args = parse_args()

    if not REQUIREMENTS_FILE.exists():
        print(f"Missing requirements file: {REQUIREMENTS_FILE}")
        return 1

    pip_base = [sys.executable, "-m", "pip"]

    if args.upgrade_pip:
        run_command(pip_base + ["install", "--upgrade", "pip"])

    install_command = pip_base + ["install", "-r", str(REQUIREMENTS_FILE)]
    if args.no_input:
        install_command.append("--no-input")

    run_command(install_command)

    print()
    print("Dependency installation complete.")
    print(
        "If you plan to use GPIB instruments on Windows, install NI-VISA or a"
        " compatible VISA driver separately."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())