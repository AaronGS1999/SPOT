"""CPU budget: how many worker threads SPOT uses when ``threads: auto``."""

from __future__ import annotations

import os
import subprocess
import sys


def performance_cores() -> int:
    """Cores worth filling with CPU-bound work.

    On Apple Silicon only the performance cores are counted (efficiency cores run BLAST several times
    slower and would hold back the last chunks); elsewhere the CPUs this process may run on.
    """
    if sys.platform == "darwin":
        for key in ("hw.perflevel0.physicalcpu", "hw.physicalcpu"):
            try:
                out = subprocess.run(["sysctl", "-n", key], capture_output=True, text=True, timeout=5)
                n = int(out.stdout.strip())
                if n > 0:
                    return n
            except (OSError, ValueError, subprocess.TimeoutExpired):
                continue
    if hasattr(os, "sched_getaffinity"):
        return max(1, len(os.sched_getaffinity(0)))
    return max(1, os.cpu_count() or 1)


def resolve(threads: int | str) -> int:
    if isinstance(threads, int):
        return threads
    return performance_cores()
