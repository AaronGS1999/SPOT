from __future__ import annotations

import bisect
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass

from .syri import Alignment

COLLINEAR = {"SYN", "INV"}


@dataclass(frozen=True)
class Anchor:
    start: int
    end: int
    other_chr: str
    other_at_start: int
    other_at_end: int
    kind: str
    block: str

    def project(self, x: int) -> float:
        if self.end == self.start:
            return float(self.other_at_start)
        slope = (self.other_at_end - self.other_at_start) / (self.end - self.start)
        return self.other_at_start + (x - self.start) * slope


@dataclass(frozen=True)
class Bracket:
    left: Anchor
    right: Anchor
    gap_left: int
    gap_right: int
    proj_left: float
    proj_right: float
    span: int

    @property
    def gap(self) -> int:
        return max(self.gap_left, self.gap_right)


class AnchorIndex:
    def __init__(self, anchors: Iterable[tuple[str, Anchor]]) -> None:
        by_chr: dict[str, list[Anchor]] = defaultdict(list)
        for chrom, anchor in anchors:
            by_chr[chrom].append(anchor)
        self._anchors = {c: sorted(v, key=lambda a: (a.start, a.end)) for c, v in by_chr.items()}
        self._starts = {c: [a.start for a in v] for c, v in self._anchors.items()}
        self._maxlen = max((a.end - a.start for v in self._anchors.values() for a in v), default=0)

    def __len__(self) -> int:
        return sum(len(v) for v in self._anchors.values())

    @classmethod
    def build(
        cls,
        alignments: Iterable[Alignment],
        side: str,
        kinds: Iterable[str],
        homologous: set[tuple[str, str]] | None = None,
    ) -> AnchorIndex:
        kinds = set(kinds)
        anchors = []
        for a in alignments:
            if a.kind not in kinds or (homologous is not None and (a.ref_chr, a.qry_chr) not in homologous):
                continue
            if side == "reference":
                anchors.append(
                    (
                        a.ref_chr,
                        Anchor(a.ref_start, a.ref_end, a.qry_chr, a.qry_start, a.qry_end, a.kind, a.block),
                    )
                )
            elif a.qry_start <= a.qry_end:
                anchors.append(
                    (
                        a.qry_chr,
                        Anchor(a.qry_start, a.qry_end, a.ref_chr, a.ref_start, a.ref_end, a.kind, a.block),
                    )
                )
            else:
                anchors.append(
                    (
                        a.qry_chr,
                        Anchor(a.qry_end, a.qry_start, a.ref_chr, a.ref_end, a.ref_start, a.kind, a.block),
                    )
                )
        return cls(anchors)

    def bracket(self, chrom: str, start: int, end: int, gapmax: int) -> Bracket | None:
        anchors = self._anchors.get(chrom)
        if not anchors:
            return None
        starts = self._starts[chrom]
        lo = bisect.bisect_left(starts, start - gapmax - self._maxlen)
        hi = bisect.bisect_right(starts, end + gapmax)
        near = anchors[lo:hi]
        lefts = [a for a in near if a.start < start]
        rights = [a for a in near if a.end > end]
        if not lefts or not rights:
            return None
        left = max(lefts, key=lambda a: min(a.end, start - 1))
        right = min(rights, key=lambda a: max(a.start, end + 1))
        xl, xr = min(left.end, start - 1), max(right.start, end + 1)
        return Bracket(
            left=left,
            right=right,
            gap_left=start - 1 - xl,
            gap_right=xr - (end + 1),
            proj_left=left.project(xl),
            proj_right=right.project(xr),
            span=xr - xl,
        )
