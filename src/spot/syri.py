from __future__ import annotations

import bisect
import json
import logging
from collections import defaultdict
from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path

from .cache import digest
from .config import SYRI_CLASSES, Config
from .errors import InputError
from .progress import Progress, silent
from .tools import as_args, run
from .util import natural_key

log = logging.getLogger("spot")


@dataclass(frozen=True)
class Alignment:
    ref_chr: str
    ref_start: int
    ref_end: int
    qry_chr: str
    qry_start: int
    qry_end: int
    kind: str
    block: str

    @property
    def length(self) -> int:
        return self.ref_end - self.ref_start + 1


def read(path: Path) -> tuple[list[Alignment], list[Alignment]]:
    """Alignments (SYNAL, INVAL, ...) and the blocks they belong to (SYN, INV, ...) from syri.out."""
    alignments: list[Alignment] = []
    blocks: list[Alignment] = []
    for n, line in enumerate(_lines(path), start=1):
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        is_alignment = f[10].endswith("AL")
        kind = f[10].removesuffix("AL")
        if kind not in SYRI_CLASSES:
            continue
        block = f[9] if is_alignment else f[8]
        try:
            a = Alignment(f[0], int(f[1]), int(f[2]), f[5], int(f[6]), int(f[7]), kind, block)
        except ValueError:
            raise InputError(f"{path}:{n}: unexpected coordinates") from None
        (alignments if is_alignment else blocks).append(a)
    if not alignments:
        raise InputError(f"{path}: no SyRI alignments found")
    return alignments, blocks


def _lines(path: Path) -> Iterator[str]:
    try:
        with open(path) as fh:
            yield from fh
    except OSError as e:
        raise InputError(f"cannot read {path}: {e.strerror}") from None


def homologous_pairs(alignments: Iterable[Alignment], min_frac: float = 0.1) -> set[tuple[str, str]]:
    """Reference-query chromosome pairs joined by syntenic (SYN) alignments.

    A pair is kept when its aligned SYN length is at least ``min_frac`` of the best pair of either of its
    two chromosomes, so fusions, fissions and duplicated chromosomes give N:M pairs. ``min_frac = 1``
    keeps only the best partner of each chromosome.
    """
    aligned: dict[tuple[str, str], int] = defaultdict(int)
    for a in alignments:
        if a.kind == "SYN":
            aligned[(a.ref_chr, a.qry_chr)] += a.length
    best_ref: dict[str, int] = defaultdict(int)
    best_qry: dict[str, int] = defaultdict(int)
    for (ref, qry), n in aligned.items():
        best_ref[ref] = max(best_ref[ref], n)
        best_qry[qry] = max(best_qry[qry], n)
    if min_frac >= 1:
        return {(r, q) for (r, q), n in aligned.items() if n == best_ref[r]}
    return {
        (r, q) for (r, q), n in aligned.items() if n >= min_frac * best_ref[r] or n >= min_frac * best_qry[q]
    }


def homology_table(
    alignments: Iterable[Alignment], pairs: set[tuple[str, str]]
) -> list[tuple[str, str, int, float]]:
    """Rows (reference chr, query chr, SYN bp, share of the reference chromosome's SYN) per chosen pair."""
    aligned: dict[tuple[str, str], int] = defaultdict(int)
    for a in alignments:
        if a.kind == "SYN":
            aligned[(a.ref_chr, a.qry_chr)] += a.length
    total: dict[str, int] = defaultdict(int)
    for (ref, _), n in aligned.items():
        total[ref] += n
    rows = [
        (r, q, aligned.get((r, q), 0), aligned.get((r, q), 0) / total[r] if total[r] else 0.0)
        for r, q in pairs
    ]
    return sorted(rows, key=lambda row: (natural_key(row[0]), natural_key(row[1])))


class AlignedIndex:
    """Merged aligned intervals of one genome (any SyRI alignment class), for coverage queries."""

    def __init__(self, alignments: Iterable[Alignment], side: str) -> None:
        spans: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for a in alignments:
            if side == "reference":
                spans[a.ref_chr].append((a.ref_start, a.ref_end))
            else:
                spans[a.qry_chr].append((min(a.qry_start, a.qry_end), max(a.qry_start, a.qry_end)))
        self._starts: dict[str, list[int]] = {}
        self._ends: dict[str, list[int]] = {}
        for chrom, items in spans.items():
            merged: list[list[int]] = []
            for s, e in sorted(items):
                if merged and s <= merged[-1][1] + 1:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])
            self._starts[chrom] = [m[0] for m in merged]
            self._ends[chrom] = [m[1] for m in merged]

    def covered(self, chrom: str, start: int, end: int) -> int:
        """Bases of [start, end] (1-based, inclusive) inside aligned intervals."""
        starts, ends = self._starts.get(chrom), self._ends.get(chrom)
        if not starts or ends is None or end < start:
            return 0
        i = max(bisect.bisect_right(starts, start) - 1, 0)
        total = 0
        while i < len(starts) and starts[i] <= end:
            lo, hi = max(starts[i], start), min(ends[i], end)
            if hi >= lo:
                total += hi - lo + 1
            i += 1
        return total


def align(config: Config, workdir: Path, progress: Progress = silent) -> Path:
    """minimap2 + SyRI. Each step is reused when its inputs and settings are unchanged, so a run that
    stopped in SyRI does not repeat the (long) alignment."""
    ref = config.species[config.reference].genome
    qry = config.species[config.query].genome
    if ref is None or qry is None:
        raise InputError("both genomes are required to align them")
    minimap2 = config.minimap2.model_dump()
    preset = minimap2.pop("x")
    genomes = {"reference": _fingerprint(ref), "query": _fingerprint(qry)}
    paf_stamp = {**genomes, "minimap2": {"x": preset, **minimap2}}
    syri_stamp = {**paf_stamp, "syri": config.syri.model_dump()}
    workdir.mkdir(parents=True, exist_ok=True)
    paf, out = workdir / "genomes.paf", workdir / "syri.out"
    threads = str(config.n_threads)
    progress("align", 0, 2)
    if _fresh(paf, workdir / "genomes.stamp.json", paf_stamp):
        log.info("reusing %s", paf)
    else:
        log.info("aligning %s (reference) and %s (query) with minimap2 -x %s", ref.name, qry.name, preset)
        run(
            ["minimap2", "-c", "--eqx", "-x", preset, "-t", threads, *as_args(minimap2), ref, qry],
            stdout=paf,
            cwd=workdir,
        )
        (workdir / "genomes.stamp.json").write_text(json.dumps(paf_stamp, indent=2))
    progress("align", 1, 2)
    if _fresh(out, workdir / "syri.stamp.json", syri_stamp):
        log.info("reusing %s", out)
    else:
        log.info("running SyRI")
        run(
            [
                "syri",
                "-c",
                paf,
                "-r",
                ref,
                "-q",
                qry,
                "-F",
                "P",
                "--dir",
                workdir,
                "--nc",
                threads,
                *as_args(config.syri.model_dump()),
            ],
            cwd=workdir,  # SyRI opens its log (syri.log) in the current folder before reading --dir
        )
        if not out.is_file():
            raise InputError(f"SyRI finished without writing {out}")
        (workdir / "syri.stamp.json").write_text(json.dumps(syri_stamp, indent=2))
    progress("align", 2, 2)
    return out


def _fresh(path: Path, stamp_file: Path, stamp: Mapping[str, object]) -> bool:
    if not (path.is_file() and path.stat().st_size > 0 and stamp_file.is_file()):
        return False
    try:
        return bool(json.loads(stamp_file.read_text()) == stamp)
    except (OSError, ValueError):
        return False


def _fingerprint(path: Path) -> dict[str, object]:
    return {"path": str(path), "sha256": digest(path)}
