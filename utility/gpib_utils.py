"""
Filename: gpib_utils.py
Author: Ethan Ruddell
Date: 2026-02-12
Description: Shared GPIB / VISA helpers used across all instruments.

Provides:
    - ``gpib_retry``   — decorator that retries flaky VISA queries
    - ``InstrumentSession`` — context-manager for safe connect / disconnect
    - ``prompt_float``, ``prompt_int``, ``prompt_choice`` — typed CLI helpers
"""

import functools
import logging
import time

log = logging.getLogger(__name__)

# ============================================================
# Retry decorator for VISA calls
# ============================================================

def gpib_retry(max_retries: int = 3, delay: float = 0.5, backoff: float = 2.0):
    """Decorator: retry a function up to *max_retries* times on Exception.

    Parameters
    ----------
    max_retries : int
        Total attempts before re-raising.
    delay : float
        Initial wait between retries (seconds).
    backoff : float
        Multiplicative factor applied to *delay* after each retry.

    Example::

        @gpib_retry(max_retries=3)
        def read_trace(inst):
            return inst.query("OUTPDTRC?")
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            current_delay = delay
            last_exc = None
            for attempt in range(1, max_retries + 1):
                try:
                    return func(*args, **kwargs)
                except Exception as exc:
                    last_exc = exc
                    log.warning(
                        "%s attempt %d/%d failed: %s",
                        func.__name__, attempt, max_retries, exc,
                    )
                    if attempt < max_retries:
                        time.sleep(current_delay)
                        current_delay *= backoff
            # All retries exhausted — re-raise the last exception
            raise last_exc  # type: ignore[misc]
        return wrapper
    return decorator


# ============================================================
# Instrument context manager
# ============================================================

class InstrumentSession:
    """Context manager that guarantees instrument cleanup.

    Parameters
    ----------
    connect_fn : callable
        Zero-argument callable that returns an instrument handle.
    disconnect_fn : callable
        One-argument callable that disconnects the handle.

    Example::

        with InstrumentSession(lcr.setup, lcr.disconnect_e4980a) as inst:
            lcr.measure_impedance(inst, 1000)
        # disconnect_e4980a(inst) is guaranteed even on exception
    """

    def __init__(self, connect_fn, disconnect_fn):
        self._connect = connect_fn
        self._disconnect = disconnect_fn
        self._inst = None

    def __enter__(self):
        self._inst = self._connect()
        return self._inst

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._inst is not None:
            try:
                self._disconnect(self._inst)
            except Exception:
                log.warning("Error during instrument disconnect", exc_info=True)
        return False  # do not suppress exceptions


# ============================================================
# CLI prompt helpers
# ============================================================

def safe_float_input(prompt: str, default: float) -> float:
    """Safely get float input from user with default value and retry on error.

    Loops until valid input is received.  Pressing Enter returns *default*.

    Example::

        freq = safe_float_input("Enter frequency (Hz) [default 1000]: ", 1000)
    """
    while True:
        try:
            raw = input(prompt).strip()
            if not raw:
                return float(default)
            return float(raw)
        except ValueError:
            print("Invalid input. Please enter a number.")


def safe_int_input(prompt: str, default: int) -> int:
    """Safely get integer input from user with default value and retry on error.

    Loops until valid input is received.  Pressing Enter returns *default*.
    """
    while True:
        try:
            raw = input(prompt).strip()
            if not raw:
                return int(default)
            return int(raw)
        except ValueError:
            print("Invalid input. Please enter an integer.")


def prompt_bool(label: str, default: bool = False) -> bool:
    """Prompt user for a yes/no answer.

    Returns *default* when the user presses Enter without typing.

    Example::

        apply = prompt_bool("Apply DC bias?", False)
    """
    hint = "Y/n" if default else "y/N"
    raw = input(f"{label} ({hint}): ").strip().lower()
    if not raw:
        return default
    return raw.startswith("y")


def prompt_choice(label: str, options: list[str], default: str) -> str:
    """Prompt user to pick from a list of string options.

    Returns the chosen value in upper-case for convenience.

    Example::

        mode = prompt_choice("Sweep type", ["LIN", "LOG"], "LOG")
    """
    opts_str = "/".join(options)
    raw = input(f"{label} ({opts_str}) [default {default}]: ").strip().upper()
    if not raw:
        return default.upper()
    if raw in [o.upper() for o in options]:
        return raw
    log.warning("Invalid choice '%s', using default '%s'", raw, default)
    return default.upper()

