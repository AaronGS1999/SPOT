from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import yaml

from . import __version__, config, orthofinder, pipeline, tools
from .config import Config
from .errors import SpotError
from .log import setup

log = logging.getLogger("spot")


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    setup(args.verbose, args.progress)
    try:
        code: int = args.func(args)
        return code
    except SpotError as e:
        log.error("%s", e)
        return 2
    except KeyboardInterrupt:
        log.error("interrupted")
        return 130
    except BrokenPipeError:  # e.g. `spot config | less` closed early
        sys.stderr.close()
        return 0


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="spot", description="Synteny-based Positional Orthology Test")
    parser.add_argument("--version", action="version", version=f"spot {__version__}")
    parser.add_argument("-v", "--verbose", action="store_true", help="debug messages")
    parser.add_argument(
        "--progress", action="store_true", help="machine-readable progress lines (for the app)"
    )
    sub = parser.add_subparsers(dest="command", required=True, metavar="COMMAND")

    settings = argparse.ArgumentParser(add_help=False)
    settings.add_argument("-c", "--config", type=Path, required=True, help="YAML configuration")
    settings.add_argument("-p", "--profile", help="built-in profile name or a profile YAML")
    settings.add_argument("-x", dest="preset", choices=["asm5", "asm10", "asm20"], help="minimap2 preset")
    settings.add_argument("--anchors", help="SyRI classes used as anchors, e.g. SYN,INV")
    settings.add_argument("-t", "--threads", type=_threads, help="worker threads, or auto (default)")
    settings.add_argument("-o", "--outdir", type=Path)
    settings.add_argument(
        "--set",
        dest="assignments",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="override any setting, e.g. validation.min_frac=0.3",
    )

    p = sub.add_parser("check", help="check external tools")
    p.add_argument("-c", "--config", type=Path, help="only flag the tools this configuration needs")
    p.set_defaults(func=_check)

    p = sub.add_parser("profiles", help="list profiles or show one")
    p.add_argument("name", nargs="?")
    p.set_defaults(func=_profiles)

    p = sub.add_parser("config", parents=[settings], help="print the resolved configuration")
    p.set_defaults(func=_config)

    p = sub.add_parser("inspect", parents=[settings], help="summarise the OrthoFinder runs")
    p.set_defaults(func=_inspect)

    p = sub.add_parser("run", parents=[settings], help="run the whole analysis")
    p.set_defaults(func=_run)

    p = sub.add_parser("gui", help="open the desktop app")
    p.add_argument("config", nargs="?", help="configuration to open")
    p.set_defaults(func=_gui)
    return parser


def _threads(text: str) -> int | str:
    if text == "auto":
        return text
    try:
        return int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected a number or auto, got {text!r}") from None


def _load(args: argparse.Namespace) -> Config:
    overrides: dict[str, Any] = {}
    for text in args.assignments:
        overrides = config.merge(overrides, config.assignment(text))
    if args.preset:
        overrides = config.merge(overrides, {"minimap2": {"x": args.preset}})
    if args.anchors:
        overrides = config.merge(overrides, {"anchors": {"types": args.anchors}})
    if args.threads is not None:
        overrides["threads"] = args.threads
    if args.outdir is not None:
        overrides["outdir"] = str(args.outdir.expanduser().resolve())
    return config.load(args.config, profile=args.profile, overrides=overrides)


def _check(args: argparse.Namespace) -> int:
    cfg = config.load(args.config) if args.config else None
    report = tools.check(cfg)
    rows = [["tool", "package", "needed", "version", "location"]]
    for s in report:
        where = s.path or "not found"
        if s.problem:
            where += f" ({s.problem})"
        rows.append([s.name, s.package, "yes" if s.required else "no", s.version or "-", where])
    print(_table(rows, numeric=False))
    packages = ", ".join(f"{n} {v or 'missing'}" for n, v in tools.python_packages())
    print(f"\npython {sys.version.split()[0]}; {packages}")
    missing = [s.name for s in report if s.required and not s.ok]
    if missing:
        log.error("missing or broken: %s", ", ".join(missing))
        return 1
    return 0


def _profiles(args: argparse.Namespace) -> int:
    if args.name:
        print(yaml.safe_dump(config.profile_settings(args.name), sort_keys=False), end="")
        return 0
    rows = [[name, config.load_profile(name)[0]] for name in config.builtin_profiles()]
    print(_table(rows, numeric=False))
    return 0


def _config(args: argparse.Namespace) -> int:
    print(config.dump(_load(args)), end="")
    return 0


def _inspect(args: argparse.Namespace) -> int:
    cfg = _load(args)
    runs = orthofinder.load_all(cfg)
    rows: list[list[str]] = []
    for run in runs:
        mapping = ", ".join(f"{sp} <- {col}" for sp, col in run.columns.items())
        kind = {True: "DNA", False: "protein", None: "input type unknown"}[run.dna]
        log.info(
            "%s: %s (OrthoFinder %s, %s, %s)", run.name, mapping, run.version or "?", kind, run.search or "?"
        )
        counts = run.summary()
        if not rows:
            rows.append(["run", *counts])
        rows.append([run.name, *(str(v) for v in counts.values())])
    print(_table(rows))
    return 0


def _run(args: argparse.Namespace) -> int:
    pipeline.run(_load(args))
    return 0


def _gui(args: argparse.Namespace) -> int:
    from .gui import main as gui_main

    return gui_main([args.config] if args.config else [])


def _table(rows: list[list[str]], numeric: bool = True) -> str:
    widths = [max(len(r[i]) for r in rows) for i in range(len(rows[0]))]
    lines = []
    for row in rows:
        cells = [
            cell.rjust(width) if numeric and i else cell.ljust(width)
            for i, (cell, width) in enumerate(zip(row, widths, strict=True))
        ]
        lines.append("  ".join(cells).rstrip())
    return "\n".join(lines)
