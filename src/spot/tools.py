from __future__ import annotations

import logging
import os
import re
import shlex
import shutil
import subprocess
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from importlib import metadata
from pathlib import Path
from typing import Any

from .config import Config
from .errors import ToolError

log = logging.getLogger("spot")

_VERSION = re.compile(r"(\d+\.\d+(?:[.\-][0-9A-Za-z]+)*)")
PYTHON_PACKAGES = ("pydantic", "PyYAML", "numpy", "pyfaidx", "matplotlib")


def as_args(options: Mapping[str, Any], *, single_dash: bool = False) -> list[str]:
    args = []
    for key, value in options.items():
        if value is None or value is False:
            continue
        args.append(("-" if single_dash or len(key) == 1 else "--") + key)
        if value is not True:
            args.append(str(value))
    return args


def run(cmd: Sequence[str | Path], *, stdout: Path | None = None) -> str:
    args = [str(c) for c in cmd]
    log.debug("$ %s", shlex.join(args))
    try:
        if stdout is None:
            proc = subprocess.run(args, capture_output=True, text=True)
        else:
            with open(stdout, "w") as fh:
                proc = subprocess.run(args, stdout=fh, stderr=subprocess.PIPE, text=True)
    except FileNotFoundError:
        raise ToolError(f"{args[0]}: command not found") from None
    if proc.returncode != 0:
        tail = "\n".join(proc.stderr.strip().splitlines()[-15:])
        raise ToolError(f"{args[0]} exited with status {proc.returncode}\n{tail}")
    return proc.stdout or ""


@dataclass(frozen=True)
class Tool:
    name: str
    package: str
    version_cmd: tuple[str, ...]
    pattern: re.Pattern[str] = _VERSION


TOOLS = (
    Tool("minimap2", "minimap2", ("minimap2", "--version")),
    Tool("syri", "syri", ("syri", "--version")),
    Tool("blastn", "blast", ("blastn", "-version")),
    Tool("java", "openjdk", ("java", "-version"), re.compile(r'version "([^"]+)"')),
)


@dataclass(frozen=True)
class Status:
    name: str
    package: str
    required: bool
    path: str | None
    version: str | None
    problem: str = ""

    @property
    def ok(self) -> bool:
        return self.path is not None and not self.problem


@dataclass(frozen=True)
class Gepard:
    jar: Path

    @property
    def version(self) -> str | None:
        m = re.search(r"Gepard-([\d.]+)\.jar$", self.jar.name)
        return m.group(1) if m else None

    @property
    def matrices(self) -> Path:
        return self.jar.parent.parent / "resources" / "matrices"

    def matrix(self, name: str) -> Path:
        path = Path(name)
        return path if path.is_absolute() else self.matrices / name


def java_executable() -> str:
    """The Java of the active environment when there is one (conda's openjdk lives in lib/jvm), then
    $JAVA_HOME, then whatever ``java`` is on PATH (on macOS /usr/bin/java is only a launcher stub)."""
    candidates = [Path(sys.prefix, "lib", "jvm", "bin", "java")]
    if home := os.environ.get("JAVA_HOME"):
        candidates.append(Path(home, "bin", "java"))
    for path in candidates:
        if path.is_file() and os.access(path, os.X_OK):
            return str(path)
    return shutil.which("java") or "java"


def find_gepard() -> Gepard | None:
    candidates: list[Path] = []
    if env := os.environ.get("SPOT_GEPARD_JAR"):
        candidates.append(Path(env))
    if script := shutil.which("gepardcmd"):
        try:
            text = Path(script).read_text(errors="replace")
        except OSError:
            text = ""
        candidates += [Path(p) for p in re.findall(r"(\S+\.jar)", text)]
    candidates += sorted(Path(sys.prefix, "share", "gepard", "dist").glob("Gepard-*.jar"), reverse=True)
    return next((Gepard(jar) for jar in candidates if jar.is_file()), None)


def version_of(cmd: Sequence[str], pattern: re.Pattern[str] = _VERSION) -> str | None:
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.TimeoutExpired):
        return None
    m = pattern.search(proc.stdout + proc.stderr)
    return m.group(1) if m else None


def required(config: Config | None) -> set[str]:
    if config is None:
        return {"minimap2", "syri", "blastn", "java", "gepard"}
    need = set()
    if config.syri_out is None:
        need |= {"minimap2", "syri"}
    if config.validation.enabled:
        need.add("blastn")
    if config.dotplots.pairs:
        need |= {"java", "gepard"}
    return need


def check(config: Config | None = None) -> list[Status]:
    need = required(config)
    report = []
    for tool in TOOLS:
        if tool.name == "java":
            java = java_executable()
            path = java if Path(java).is_file() else shutil.which(java)
            cmd: tuple[str, ...] = (java, *tool.version_cmd[1:])
        else:
            path = shutil.which(tool.version_cmd[0])
            cmd = tool.version_cmd
        version = version_of(cmd, tool.pattern) if path else None
        report.append(Status(tool.name, tool.package, tool.name in need, path, version))
    gepard = find_gepard()
    if gepard is None:
        report.append(Status("gepard", "gepard", "gepard" in need, None, None))
    else:
        matrix = gepard.matrix(config.gepard.matrix if config else "edna.mat")
        problem = "" if matrix.is_file() else f"matrix not found: {matrix}"
        report.append(Status("gepard", "gepard", "gepard" in need, str(gepard.jar), gepard.version, problem))
    return report


def python_packages() -> list[tuple[str, str | None]]:
    versions: list[tuple[str, str | None]] = []
    for name in PYTHON_PACKAGES:
        try:
            versions.append((name, metadata.version(name)))
        except metadata.PackageNotFoundError:
            versions.append((name, None))
    return versions
