from __future__ import annotations

import bisect
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass, field

from .anchors import COLLINEAR, AnchorIndex
from .config import Anchors
from .locate import Location
from .syri import AlignedIndex

SYNTENIC_ORTHOLOG = "syntenic_ortholog"
FLANKS_NOT_CONFIRMED = "syntenic_locus_flanks_not_confirmed"
NOT_SIGNIFICANT = "syntenic_locus_not_significant"  # a member is in the window, but no more than by chance
OTHER = "syntenic_locus_other"
EMPTY_ALIGNED = "empty_locus_aligned"  # the orthologous locus is aligned and holds no sequence of this set
EMPTY_GAP = "empty_locus_gap"  # the projected window is mostly unaligned, so absence cannot be told
EVALUABLE = (SYNTENIC_ORTHOLOG, OTHER, EMPTY_ALIGNED, FLANKS_NOT_CONFIRMED, NOT_SIGNIFICANT)
NOT_EVALUABLE = (
    "no_anchor",
    "anchors_too_far",
    "anchors_split",
    "anchors_inconsistent",
    EMPTY_GAP,
    "unplaced",
)

# projected anchor ends further apart than this are not one locus
SPREAD_FACTOR = 3
SPREAD_SLACK = 10_000


@dataclass
class Test:
    gene: str
    species: str
    klass: str
    evidence: str = ""
    block_type: str = ""
    blocks: str = ""
    gap_left: int | None = None
    gap_right: int | None = None
    other_chr: str = ""
    window_start: int | None = None
    window_end: int | None = None
    window_aligned: float | None = None
    orthologs: list[str] = field(default_factory=list)
    others: list[str] = field(default_factory=list)
    unconfirmed: list[str] = field(default_factory=list)
    flank_left: int | None = None
    flank_right: int | None = None
    flank_back_left: int | None = None
    flank_back_right: int | None = None
    inverted: bool | None = None
    p_value: float | None = None
    q_value: float | None = None


class GeneIndex:
    def __init__(self, locations: Iterable[Location]) -> None:
        by_chr: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
        for loc in locations:
            if loc.located:
                by_chr[loc.chrom].append((loc.start, loc.end, loc.gene))
        self._genes = {c: sorted(v) for c, v in by_chr.items()}
        self._starts = {c: [g[0] for g in v] for c, v in self._genes.items()}
        self._maxlen = max((e - s for v in self._genes.values() for s, e, _ in v), default=0)

    def overlapping(self, chrom: str, lo: int, hi: int) -> list[str]:
        genes = self._genes.get(chrom, [])
        i = bisect.bisect_left(self._starts.get(chrom, []), lo - self._maxlen)
        out = []
        for start, end, gene in genes[i:]:
            if start > hi:
                break
            if end >= lo:
                out.append(gene)
        return out


def window_tolerance(gap: int, settings: Anchors) -> int:
    """Slack added to both ends of the projected window: a fixed minimum, or a fraction of the distance
    to the farthest anchor when that is larger (indels between anchor and locus shift the projection)."""
    return max(settings.tol, int(settings.tol_frac * gap))


def test(
    loc: Location,
    same_group: set[str],
    anchors: AnchorIndex,
    others: GeneIndex,
    settings: Anchors,
    aligned: AlignedIndex | None = None,
) -> Test:
    if not loc.located:
        return Test(loc.gene, loc.species, "unplaced")
    b = anchors.bracket(loc.chrom, loc.start, loc.end, settings.gapmax)
    if b is None:
        return Test(loc.gene, loc.species, "no_anchor")
    t = Test(loc.gene, loc.species, "", gap_left=b.gap_left, gap_right=b.gap_right)
    if b.gap > settings.gapmax:
        t.klass = "anchors_too_far"
        return t
    if b.left.other_chr != b.right.other_chr:
        t.klass = "anchors_split"
        return t
    if abs(b.proj_right - b.proj_left) > SPREAD_FACTOR * b.span + SPREAD_SLACK:
        t.klass = "anchors_inconsistent"
        return t
    blocks = [b.left.block] if b.left.block == b.right.block else [b.left.block, b.right.block]
    kinds = {b.left.kind, b.right.kind}
    t.blocks = "|".join(blocks)
    t.block_type = ("collinear" if kinds <= COLLINEAR else "rearranged") + (
        "_two_blocks" if len(blocks) > 1 else ""
    )
    t.evidence = "flank_anchors" if b.gap <= settings.flank else "bracketing_anchors"
    t.other_chr = b.left.other_chr
    tol = window_tolerance(b.gap, settings)
    t.window_start = max(int(min(b.proj_left, b.proj_right)) - tol, 1)
    t.window_end = int(max(b.proj_left, b.proj_right)) + tol
    found = others.overlapping(t.other_chr, t.window_start, t.window_end)
    t.orthologs = [g for g in found if g in same_group]
    t.others = [g for g in found if g not in same_group]
    if t.orthologs:
        t.klass = SYNTENIC_ORTHOLOG
    elif t.others:
        t.klass = OTHER
    else:
        t.klass = EMPTY_ALIGNED
        if aligned is not None:
            span = t.window_end - t.window_start + 1
            t.window_aligned = aligned.covered(t.other_chr, t.window_start, t.window_end) / span
            if t.window_aligned < settings.min_aligned_frac:
                t.klass = EMPTY_GAP
    return t
