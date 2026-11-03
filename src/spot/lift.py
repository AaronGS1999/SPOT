from __future__ import annotations

import bisect
from dataclasses import dataclass
from pathlib import Path

from .errors import InputError

AGP_COMPONENTS = set("ADFGOPW")


@dataclass(frozen=True)
class Segment:
    obj: str
    obj_start: int
    comp: str
    comp_start: int
    comp_end: int
    strand: str


class LiftStep:
    """One assembly step: sequences of the previous assembly (components) placed into objects.

    Components the step does not mention are passed through unchanged, so a table that only
    splits one chromosome can be chained with an AGP that covers everything else.
    """

    def __init__(self, source: Path, segments: list[Segment]) -> None:
        self.source = source
        self._segments: dict[str, list[Segment]] = {}
        for seg in sorted(segments, key=lambda s: (s.comp, s.comp_start)):
            self._segments.setdefault(seg.comp, []).append(seg)
        self._starts = {comp: [s.comp_start for s in segs] for comp, segs in self._segments.items()}

    @classmethod
    def read(cls, path: Path) -> LiftStep:
        try:
            lines = path.read_text().splitlines()
        except OSError as e:
            raise InputError(f"cannot read {path}: {e.strerror}") from None
        agp = path.suffix.lower() == ".agp" or any(line.startswith("##agp") for line in lines[:5])
        segments = _agp(lines, path) if agp else _table(lines, path)
        if not segments:
            raise InputError(f"{path}: no segments found")
        return cls(path, segments)

    def map(self, chrom: str, start: int, end: int) -> tuple[str, int, int] | None:
        segs = self._segments.get(chrom)
        if segs is None:
            return chrom, start, end
        i = bisect.bisect_right(self._starts[chrom], start) - 1
        if i < 0 or end > segs[i].comp_end:
            return None
        seg = segs[i]
        if seg.strand == "-":
            return seg.obj, seg.obj_start + seg.comp_end - end, seg.obj_start + seg.comp_end - start
        return seg.obj, seg.obj_start + start - seg.comp_start, seg.obj_start + end - seg.comp_start


class Lift:
    def __init__(self, steps: list[LiftStep], strip_suffix: str | None = None) -> None:
        self.steps = steps
        self.strip_suffix = strip_suffix

    @classmethod
    def read(cls, paths: list[Path], strip_suffix: str | None = None) -> Lift:
        return cls([LiftStep.read(p) for p in paths], strip_suffix)

    def __call__(self, chrom: str, start: int, end: int) -> tuple[str, int, int] | None:
        loc: tuple[str, int, int] = (chrom, start, end)
        for step in self.steps:
            mapped = step.map(*loc)
            if mapped is None:
                return None
            loc = mapped
        chrom, start, end = loc
        if self.strip_suffix and chrom.endswith(self.strip_suffix):
            chrom = chrom[: -len(self.strip_suffix)]
        return chrom, start, end


def _agp(lines: list[str], path: Path) -> list[Segment]:
    segments = []
    for n, line in enumerate(lines, start=1):
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9:
            raise InputError(f"{path}:{n}: expected 9 AGP columns")
        if f[4] not in AGP_COMPONENTS:
            continue
        try:
            segments.append(Segment(f[0], int(f[1]), f[5], int(f[6]), int(f[7]), "-" if f[8] == "-" else "+"))
        except ValueError:
            raise InputError(f"{path}:{n}: non-numeric coordinates") from None
    return segments


def _table(lines: list[str], path: Path) -> list[Segment]:
    segments: list[Segment] = []
    for n, line in enumerate(lines, start=1):
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split()
        if len(f) < 4:
            raise InputError(f"{path}:{n}: expected name, parent, start, end[, strand]")
        try:
            start, end = int(f[2]), int(f[3])
        except ValueError:
            if not segments:
                continue
            raise InputError(f"{path}:{n}: non-numeric coordinates") from None
        segments.append(Segment(f[0], 1, f[1], start, end, f[4] if len(f) > 4 and f[4] == "-" else "+"))
    return segments
