"""
Filename: usb_switchbox.py
Author: Ethan Ruddell
Date: 2026-04-02
Description: USB switchover box (multiplexer) helper for automatic instrument routing.

This module provides a small abstraction around a USB serial-connected
switchbox so GUI and CLI flows can route DUT connections to the selected
instrument before a measurement runs.

Hardware/protocol values are read from config.yaml under `usb_switchbox`.
All command strings are placeholders by default and should be replaced with
real commands from your switchbox manual.
"""

from __future__ import annotations

import logging
import threading
import importlib
import re
from dataclasses import dataclass

import config

log = logging.getLogger(__name__)


class SwitchBoxError(RuntimeError):
    """Raised when USB switchbox operations fail."""


@dataclass
class SwitchBoxSettings:
    """Configuration for the USB switchbox serial interface."""

    enabled: bool = False
    port: str = "COM__PLACEHOLDER__"
    baudrate: int = 9600
    timeout_s: float = 2.0
    write_terminator: str = ""
    read_response: bool = True
    command_map: dict[str, str] | None = None


DEFAULT_COMMAND_MAP = {
    "PIA": "CTRL_A",
    "PSPA": "CTRL_B",
    "LCR": "CTRL_C",
}


def _decode_command_to_bytes(command: str, terminator: str = "") -> bytes:
    """Convert configured command string to bytes for serial transmission.

    Supported forms:
      - CTRL_A / CTRL-A ... CTRL_Z / CTRL-Z
      - HEX:01 (single byte, hex)
      - TEXT:... (ASCII text command)
      - Raw text (ASCII)
    """
    command = str(command).strip()
    normalized = command.upper().replace("-", "_")

    m = re.fullmatch(r"CTRL_([A-Z])", normalized)
    if m:
        ctrl_val = ord(m.group(1)) - 64
        payload = bytes([ctrl_val])
        return payload

    if normalized.startswith("HEX:"):
        hex_part = command.split(":", 1)[1].strip().replace("0X", "")
        if len(hex_part) % 2 != 0:
            hex_part = f"0{hex_part}"
        return bytes.fromhex(hex_part)

    if normalized.startswith("TEXT:"):
        text_part = command.split(":", 1)[1]
        return f"{text_part}{terminator}".encode("ascii", errors="strict")

    return f"{command}{terminator}".encode("ascii", errors="strict")


def load_switchbox_settings() -> SwitchBoxSettings:
    """Load switchbox settings from config.yaml with safe fallbacks."""
    section = config.get_section("usb_switchbox")
    command_map = section.get("command_map", {}) if isinstance(section, dict) else {}

    merged_map = dict(DEFAULT_COMMAND_MAP)
    if isinstance(command_map, dict):
        for key, value in command_map.items():
            merged_map[str(key).upper()] = str(value)

    return SwitchBoxSettings(
        enabled=bool(section.get("enabled", False)),
        port=str(section.get("port", "COM__PLACEHOLDER__")),
        baudrate=int(section.get("baudrate", 9600)),
        timeout_s=float(section.get("timeout_s", 2.0)),
        write_terminator=str(section.get("write_terminator", "")),
        read_response=bool(section.get("read_response", True)),
        command_map=merged_map,
    )


class UsbSwitchBox:
    """USB serial switchbox controller."""

    def __init__(self, settings: SwitchBoxSettings | None = None):
        self.settings = settings or load_switchbox_settings()
        self._serial = None
        self._lock = threading.Lock()
        self.last_routed_instrument = None

    @property
    def enabled(self) -> bool:
        """Whether switchbox control is enabled."""
        return self.settings.enabled

    def set_enabled(self, enabled: bool) -> None:
        """Enable/disable switchbox control at runtime."""
        self.settings.enabled = bool(enabled)

    def close(self) -> None:
        """Close serial connection if open."""
        with self._lock:
            if self._serial is not None:
                try:
                    self._serial.close()
                except Exception:
                    pass
                self._serial = None

    def _ensure_connection(self):
        """Open serial connection if needed."""
        if self._serial is not None and getattr(self._serial, "is_open", False):
            return self._serial

        try:
            serial_module = importlib.import_module("serial")
        except ImportError as exc:
            raise SwitchBoxError(
                "pyserial is required for USB switchbox control. Install with: pip install pyserial"
            ) from exc

        try:
            self._serial = serial_module.Serial(
                port=self.settings.port,
                baudrate=self.settings.baudrate,
                timeout=self.settings.timeout_s,
                write_timeout=self.settings.timeout_s,
                bytesize=serial_module.EIGHTBITS,
                parity=serial_module.PARITY_NONE,
                stopbits=serial_module.STOPBITS_ONE,
                xonxoff=False,
                rtscts=False,
                dsrdtr=False,
            )
            return self._serial
        except Exception as exc:
            raise SwitchBoxError(
                f"Unable to open switchbox serial port '{self.settings.port}': {exc}"
            ) from exc

    def switch_to_instrument(self, instrument: str) -> str:
        """Route the switchbox to the given instrument and return a status message."""
        instrument = str(instrument).upper().strip()

        if not self.enabled:
            return "Switchbox disabled"

        command_map = self.settings.command_map or DEFAULT_COMMAND_MAP
        command = command_map.get(instrument)
        if not command:
            raise SwitchBoxError(f"No switch command configured for instrument '{instrument}'")

        with self._lock:
            ser = self._ensure_connection()
            payload = _decode_command_to_bytes(command, self.settings.write_terminator)

            try:
                ser.write(payload)
                ser.flush()
            except Exception as exc:
                raise SwitchBoxError(f"Failed to write switch command for {instrument}: {exc}") from exc

            response = ""
            if self.settings.read_response:
                try:
                    response = ser.readline().decode("utf-8", errors="replace").strip()
                except Exception as exc:
                    raise SwitchBoxError(f"Failed to read switchbox response: {exc}") from exc

            self.last_routed_instrument = instrument
            if response:
                msg = f"Switchbox routed to {instrument} (response: {response})"
            else:
                msg = f"Switchbox routed to {instrument}"
            log.info(msg)
            return msg


def create_switchbox_from_config() -> UsbSwitchBox:
    """Factory helper for app code."""
    return UsbSwitchBox(load_switchbox_settings())
