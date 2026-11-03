"""Quick checks of how a species' sequences will be placed, without running the analysis.

Used by the app to tell, while a project is being set up, whether the coordinates in the sequence names
are found and whether they fall on sequences of the chosen genome (or need lift files).
"""

from __future__ import annotations

import re
from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path

from .config import HEADER_REGEX, column_pattern
from .lift import Lift
from .orthofinder import _sequence_ids
from .util import natural_key

HEADER_CANDIDATES = (
    HEADER_REGEX,  # chr1:100-200
    r"(?P<chr>[^\s:]+):(?P<start>\d+)\.\.(?P<end>\d+)",  # chr1:100..200
    r"(?P<chr>[^\s_|]+)[_|](?P<start>\d+)[_|-](?P<end>\d+)",  # chr1_100_200, chr1|100-200
)


@dataclass
class SpeciesReport:
    sequences: int = 0
    with_coordinates: int = 0
    on_genome: int | None = None
    chromosomes: Counter[str] = field(default_factory=Counter)
    missing: Counter[str] = field(default_factory=Counter)
    genome_names: list[str] = field(default_factory=list)
    lift_failed: int = 0
    example: str = ""

    def text(self) -> str:
        if self.sequences == 0:
            return "no sequences of this species in the OrthoFinder runs"
        parts = [f"{self.sequences:,} sequences"]
        if self.with_coordinates < self.sequences:
            parts.append(f"{self.with_coordinates:,} with coordinates in their names")
        else:
            parts.append(f"coordinates in every name (e.g. {self.example})")
        if self.on_genome is not None:
            parts.append(f"{self.on_genome:,} on sequences of the genome")
        lines = [" · ".join(parts)]
        if self.lift_failed:
            lines.append(f"{self.lift_failed:,} fall outside the lift files")
        if self.on_genome is not None and self.missing:
            names = ", ".join(n for n, _ in sorted(self.missing.items(), key=lambda x: natural_key(x[0]))[:5])
            more = "…" if len(self.missing) > 5 else ""
            lines.append(f"{sum(self.missing.values()):,} on sequences not in the genome ({names}{more})")
            if self.on_genome == 0 and self.genome_names:
                shown = ", ".join(self.genome_names[:3])
                lines.append(
                    f"the genome has {shown}…: if the coordinates refer to another assembly, add its lift "
                    "files (AGP or split table)"
                )
        return "\n".join(lines)


def headers(run_paths: Iterable[Path], match: str) -> list[str]:
    """Descriptions (FASTA header after the ID) of the sequences of the species matching ``match``."""
    pattern = re.compile(column_pattern(match)) if match else None
    out: list[str] = []
    for path in run_paths:
        descriptions, origin, _ = _sequence_ids(path)
        for gene, desc in descriptions.items():
            if pattern is None or pattern.search(origin.get(gene, "")):
                out.append(desc)
    return out


def detect_header_regex(descriptions: list[str]) -> str | None:
    """The first known coordinate pattern that matches most descriptions."""
    best, score = None, 0
    for candidate in HEADER_CANDIDATES:
        rx = re.compile(candidate)
        n = sum(1 for d in descriptions[:500] if rx.search(d))
        if n > score:
            best, score = candidate, n
    return best


def genome_names(path: Path) -> list[str]:
    """Sequence names of a FASTA, from its .fai index when present (fast), else from its headers."""
    fai = Path(str(path) + ".fai")
    if fai.is_file() and fai.stat().st_mtime >= path.stat().st_mtime:
        return [line.split("\t", 1)[0] for line in fai.read_text().splitlines() if line]
    names = []
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                names.append(line[1:].split()[0].decode())
    return names


def species_report(
    descriptions: list[str],
    header_regex: str,
    genome: Path | None = None,
    lift: list[Path] | None = None,
    strip_suffix: str | None = None,
) -> SpeciesReport:
    report = SpeciesReport(sequences=len(descriptions))
    rx = re.compile(header_regex)
    lifter = Lift.read(lift or [], strip_suffix) if (lift or strip_suffix) else None
    coords = []
    for d in descriptions:
        m = rx.search(d)
        if m is None:
            continue
        report.with_coordinates += 1
        if not report.example:
            report.example = m.group(0)
        chrom = m["chr"]
        if lifter is not None:
            lifted = lifter(chrom, int(m["start"]), int(m["end"]))
            if lifted is None:
                report.lift_failed += 1
                continue
            chrom = lifted[0]
        coords.append(chrom)
    report.chromosomes = Counter(coords)
    if genome is not None and genome.is_file():
        names = genome_names(genome)
        report.genome_names = sorted(names, key=natural_key)
        known = set(names)
        report.on_genome = sum(n for c, n in report.chromosomes.items() if c in known)
        report.missing = Counter({c: n for c, n in report.chromosomes.items() if c not in known})
    return report


def column_prefixes(run_paths: Iterable[Path]) -> list[str]:
    """Species prefixes of the OrthoFinder columns (text before the first _ or -), in order of appearance."""
    seen: list[str] = []
    for path in run_paths:
        table = path / "Orthogroups" / "Orthogroups.tsv"
        if not table.is_file():
            continue
        with open(table, encoding="utf-8") as fh:
            header = fh.readline().rstrip("\r\n").split("\t")[1:]
        for column in header:
            prefix = re.split(r"[_\-.]", column, maxsplit=1)[0]
            if prefix and prefix not in seen:
                seen.append(prefix)
    return seen
