"""Flank validation of candidate positional orthologs with BLASTN.

The two loci are extracted with ``window`` bp on each side and aligned with BLASTN. In ``colinear`` mode
(the default) an HSP only counts for a flank when it joins the matching flank of the other locus: left
flank to left flank and right to right on the plus strand, or left to right and right to left on the
minus strand (a pair inside an inversion). The orientation is chosen from the HSPs, not from the strands
of the two sequences, because NUXTs of two species may come from organelle genomes assembled in opposite
orientations. ``coverage`` mode reproduces the v3 script: any HSP anywhere counts. With ``reciprocal``
the same HSPs are also read from the candidate's side, so both flanks of both loci must pass.
"""

from __future__ import annotations

import logging
import tempfile
from collections.abc import Iterable, Mapping
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

from .config import Blastn, Validation
from .locate import Genome, Location
from .tools import as_args, run

log = logging.getLogger("spot")


@dataclass(frozen=True)
class Region:
    """A locus with its flanks, as plus-strand genome sequence."""

    name: str
    chrom: str
    start: int
    end: int
    seq: str
    before: int
    core: int

    @property
    def after(self) -> int:
        return len(self.seq) - self.before - self.core


@dataclass(frozen=True)
class Hsp:
    """One BLASTN HSP (1-based, inclusive). qstart < qend; sstart > send on the minus strand."""

    pident: float
    length: int
    qstart: int
    qend: int
    sstart: int
    send: int

    def swapped(self) -> Hsp:
        """The same HSP read with subject and query exchanged."""
        if self.sstart <= self.send:
            return Hsp(self.pident, self.length, self.sstart, self.send, self.qstart, self.qend)
        return Hsp(self.pident, self.length, self.send, self.sstart, self.qend, self.qstart)


@dataclass(frozen=True)
class FlankCheck:
    ok: bool
    left: int
    right: int
    back_left: int | None = None
    back_right: int | None = None
    inverted: bool | None = None


def region(genome: Genome, loc: Location, window: int) -> Region:
    start = max(1, loc.start - window)
    end = min(genome.lengths[loc.chrom], loc.end + window)
    return Region(
        loc.gene,
        loc.chrom,
        start,
        end,
        genome.fetch(loc.chrom, start, end),
        loc.start - start,
        loc.end - loc.start + 1,
    )


def required(available: int, min_frac: float, min_bp: int) -> float:
    """Aligned bases a flank needs: the larger of a fraction and a minimum, capped at half the flank."""
    return min(max(min_bp, min_frac * available), 0.5 * available)


def _union(spans: Iterable[tuple[int, int]]) -> int:
    """Total length of half-open spans [s, e) after merging overlaps."""
    total, cur_s, cur_e = 0, 0, 0
    for s, e in sorted(spans):
        if e <= s:
            continue
        if s > cur_e:
            total += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    return total + cur_e - cur_s


def covered(hsps: Iterable[tuple[int, int]], lo: int, hi: int) -> int:
    """Bases of [lo, hi) (0-based) covered by HSPs given as 1-based inclusive query coordinates."""
    return _union((max(s - 1, lo), min(e, hi)) for s, e in hsps if e > lo and s - 1 < hi)


def colinear(
    hsps: Iterable[Hsp], query: Region, subject: Region, max_offset: int, inverted: bool
) -> tuple[int, int]:
    """Bases of the query's left and right flank aligned colinearly to the matching flank of the subject.

    With ``inverted`` the subject locus is read on its reverse strand, so the query's left flank must join
    the subject's right flank on the minus strand (a pair inside an inversion). The part of an HSP that
    falls in a query flank must land on the same side of the subject locus, overshooting it by at most
    ``max_offset`` bp. Distance to the locus is not constrained, so indels next to a locus are tolerated.
    """
    bq, cq, lq, ls = query.before, query.core, len(query.seq), len(subject.seq)
    bs = subject.after if inverted else subject.before
    cs = subject.core
    left: list[tuple[int, int]] = []
    right: list[tuple[int, int]] = []
    for h in hsps:
        sa, sb = (ls - h.sstart + 1, ls - h.send + 1) if inverted else (h.sstart, h.send)
        if sb <= sa and h.qend > h.qstart:
            continue  # joined in the other orientation
        q1, q2 = h.qstart, h.qend
        slope = (sb - sa) / (q2 - q1) if q2 > q1 else 1.0

        def at(x: int, q1: int = q1, sa: int = sa, slope: float = slope) -> float:
            return sa + (x - q1) * slope

        a, b = max(q1, 1), min(q2, bq)
        if a <= b and all(at(x) <= bs + max_offset for x in (a, b)):
            left.append((a - 1, b))
        a, b = max(q1, bq + cq + 1), min(q2, lq)
        if a <= b and all(at(x) > bs + cs - max_offset for x in (a, b)):
            right.append((a - 1, b))
    return _union(left), _union(right)


def best_colinear(hsps: list[Hsp], query: Region, subject: Region, max_offset: int) -> tuple[int, int, bool]:
    """Colinear flank bases in the better of the two relative orientations (direct or inverted)."""
    direct = colinear(hsps, query, subject, max_offset, inverted=False)
    inverse = colinear(hsps, query, subject, max_offset, inverted=True)
    if min(inverse) > min(direct) or (min(inverse) == min(direct) and sum(inverse) > sum(direct)):
        return inverse[0], inverse[1], True
    return direct[0], direct[1], False


def _passes(left: int, right: int, before: int, after: int, s: Validation) -> bool:
    return (
        before > 0
        and after > 0
        and left >= required(before, s.min_frac, s.min_bp)
        and right >= required(after, s.min_frac, s.min_bp)
    )


def judge(hsps: list[Hsp], tested: Region, candidate: Region, settings: Validation) -> FlankCheck:
    """Decide whether a candidate pair is confirmed by its flanks, from the HSPs of tested vs candidate."""
    s = settings
    back = [h.swapped() for h in hsps]
    if s.mode == "coverage":
        left = covered(((h.qstart, h.qend) for h in hsps), 0, tested.before)
        right = covered(((h.qstart, h.qend) for h in hsps), tested.before + tested.core, len(tested.seq))
        ok = _passes(left, right, tested.before, tested.after, s)
        if not s.reciprocal:
            return FlankCheck(ok, left, right)
        bl = covered(((h.qstart, h.qend) for h in back), 0, candidate.before)
        br = covered(
            ((h.qstart, h.qend) for h in back), candidate.before + candidate.core, len(candidate.seq)
        )
        return FlankCheck(ok and _passes(bl, br, candidate.before, candidate.after, s), left, right, bl, br)
    left, right, inverted = best_colinear(hsps, tested, candidate, s.max_offset)
    ok = _passes(left, right, tested.before, tested.after, s)
    if not s.reciprocal:
        return FlankCheck(ok, left, right, inverted=inverted)
    bl, br = colinear(back, candidate, tested, s.max_offset, inverted)
    back_ok = _passes(bl, br, candidate.before, candidate.after, s)
    return FlankCheck(ok and back_ok, left, right, bl, br, inverted)


def parse_hsps(text: str, min_identity: float, min_length: int) -> list[Hsp]:
    hsps = []
    for line in text.splitlines():
        fields = line.split()
        if len(fields) < 6:
            continue
        qs, qe, ss, se = (int(x) for x in fields[2:6])
        h = Hsp(float(fields[0]), int(fields[1]), qs, qe, ss, se)
        if h.pident >= min_identity and h.length >= min_length:
            hsps.append(h)
    return hsps


class Validator:
    def __init__(
        self, settings: Validation, blastn: Blastn, genomes: Mapping[str, Genome], threads: int
    ) -> None:
        self.settings = settings
        self.identity = blastn.perc_identity
        self.options = as_args(blastn.model_dump(), single_dash=True)
        self.genomes = genomes
        self.threads = threads

    def check(self, jobs: Iterable[tuple[Location, Location]]) -> dict[tuple[str, str], FlankCheck]:
        regions: dict[str, Region] = {}
        pairs: list[tuple[str, str]] = []
        for tested, candidate in jobs:
            for loc in (tested, candidate):
                if loc.gene not in regions:
                    regions[loc.gene] = region(self.genomes[loc.species], loc, self.settings.window)
            pairs.append((tested.gene, candidate.gene))
        pairs = list(dict.fromkeys(pairs))
        if not pairs:
            return {}
        log.info(
            "validating %d candidate pairs with BLASTN (+-%d bp, %s%s)",
            len(pairs),
            self.settings.window,
            self.settings.mode,
            ", reciprocal" if self.settings.reciprocal else "",
        )
        with tempfile.TemporaryDirectory(prefix="spot_blastn_") as tmp:
            files: dict[str, Path] = {}
            for i, (name, reg) in enumerate(regions.items()):
                files[name] = Path(tmp, f"{i}.fa")
                files[name].write_text(f">{name}\n{reg.seq}\n")

            def one(pair: tuple[str, str]) -> FlankCheck:
                a, b = pair
                return self._one(regions[a], regions[b], files[a], files[b])

            with ThreadPoolExecutor(max_workers=self.threads) as pool:
                return dict(zip(pairs, pool.map(one, pairs), strict=True))

    def _one(self, tested: Region, candidate: Region, query: Path, subject: Path) -> FlankCheck:
        out = run(
            [
                "blastn",
                "-query",
                query,
                "-subject",
                subject,
                "-outfmt",
                "6 pident length qstart qend sstart send",
                *self.options,
            ]
        )
        return judge(parse_hsps(out, self.identity, self.settings.min_hsp), tested, candidate, self.settings)
