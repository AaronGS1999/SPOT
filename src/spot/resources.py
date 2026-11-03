"""CPU budget: how many worker threads SPOT uses when ``threads: auto``."""

from __future__ import annotations

import os
import subprocess
import sys


def performance_cores() -> int:
    """Cores worth filling with CPU-bound work.

    On Apple Silicon the cores of every performance level except the slowest (the efficiency cores, which
    run BLAST several times slower and would hold back the last chunks); elsewhere the CPUs this process
    may run on.
    """
    if sys.platform == "darwin":
        levels = _sysctl("hw.nperflevels")
        if levels and levels > 1:
            cores = [_sysctl(f"hw.perflevel{i}.physicalcpu") or 0 for i in range(levels - 1)]
            if sum(cores) > 0:
                return sum(cores)
        physical = _sysctl("hw.physicalcpu")
        if physical:
            return physical
    if hasattr(os, "sched_getaffinity"):
        return max(1, len(os.sched_getaffinity(0)))
    return max(1, os.cpu_count() or 1)


def _sysctl(key: str) -> int | None:
    # full path: apps started from the Finder may not have /usr/sbin on PATH
    for exe in ("/usr/sbin/sysctl", "sysctl"):
        try:
            out = subprocess.run([exe, "-n", key], capture_output=True, text=True, timeout=5)
            return int(out.stdout.strip())
        except (OSError, ValueError, subprocess.TimeoutExpired):
            continue
    return None


def resolve(threads: int | str) -> int:
    if isinstance(threads, int):
        return threads
    return performance_cores()
