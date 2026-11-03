"""Progress reporting shared by the engine, the CLI and the desktop app.

Every stage reports ``(stage, done, total)``. The default reporter writes machine-readable lines
``@progress <stage> <done> <total>`` to the ``spot.progress`` logger at most every ``interval`` seconds
(and always at 0 and at the end), which the app reads from the engine's output.
"""

from __future__ import annotations

import logging
import threading
import time
from collections.abc import Callable

Progress = Callable[[str, int, int], None]
PREFIX = "@progress"
STAGES = ("locate", "align", "test", "validate", "export", "dotplots", "plots")

_log = logging.getLogger("spot.progress")


class Reporter:
    def __init__(self, interval: float = 0.5) -> None:
        self.interval = interval
        self._last: dict[str, float] = {}
        self._lock = threading.Lock()

    def __call__(self, stage: str, done: int, total: int) -> None:
        now = time.monotonic()
        with self._lock:
            if 0 < done < total and now - self._last.get(stage, 0.0) < self.interval:
                return
            self._last[stage] = now
        _log.info("%s %s %d %d", PREFIX, stage, done, total)


def parse(line: str) -> tuple[str, int, int] | None:
    """(stage, done, total) from a progress line, or None."""
    i = line.find(PREFIX)
    if i < 0:
        return None
    parts = line[i + len(PREFIX) :].split()
    if len(parts) != 3 or not parts[1].isdigit() or not parts[2].isdigit():
        return None
    return parts[0], int(parts[1]), int(parts[2])


def silent(stage: str, done: int, total: int) -> None:
    del stage, done, total
