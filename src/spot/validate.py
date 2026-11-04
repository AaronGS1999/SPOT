"""Flank validation of candidate positional orthologs with BLASTN.

The two loci are extracted with ``window`` bp on each side and aligned with BLASTN. An HSP only counts
for a flank when it joins the matching flank of the other locus: left flank to left flank and right to
right on the plus strand, or left to right and right to left on the minus strand (a pair inside an
inversion). The orientation is chosen from the HSPs, not from the strands
of the two sequences, because NUXTs of two species may come from organelle genomes assembled in opposite
orientations. With ``reciprocal`` the same HSPs are also read from the candidate's side, so both flanks
of both loci must pass.
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
from .progress import Progress, silent
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
    left, right, inverted = best_colinear(hsps, tested, candidate, s.max_offset)
    ok = _passes(left, right, tested.before, tested.after, s)
    if not s.reciprocal:
        return FlankCheck(ok, left, right, inverted=inverted)
    bl, br = colinear(back, candidate, tested, s.max_offset, inverted)
    back_ok = _passes(bl, br, candidate.before, candidate.after, s)
    return FlankCheck(ok and back_ok, left, right, bl, br, inverted)


class Validator:
    """Runs BLASTN for all candidate pairs at once.

    All flanked loci go into one BLAST database; the tested loci are split into chunks searched in
    parallel, one single-threaded ``blastn`` per worker, each restricted (``-seqidlist``) to the
    candidates of its own loci, and only the HSPs of the requested pairs are kept. ``-dbsize`` is set
    to the length of one flanked locus, so E-values are those of a pairwise comparison and do not grow
    with the number of pairs.
    """

    def __init__(
        self,
        settings: Validation,
        blastn: Blastn,
        genomes: Mapping[str, Genome],
        threads: int,
        progress: Progress = silent,
    ) -> None:
        self.settings = settings
        self.identity = blastn.perc_identity
        self.options = as_args(blastn.model_dump(), single_dash=True)
        self.genomes = genomes
        self.threads = max(1, threads)
        self.progress = progress

    def regions(
        self, jobs: Iterable[tuple[Location, Location]]
    ) -> tuple[dict[str, Region], list[tuple[str, str]]]:
        regions: dict[str, Region] = {}
        pairs: list[tuple[str, str]] = []
        for tested, candidate in jobs:
            for loc in (tested, candidate):
                if loc.gene not in regions:
                    regions[loc.gene] = region(self.genomes[loc.species], loc, self.settings.window)
            pairs.append((tested.gene, candidate.gene))
        return regions, list(dict.fromkeys(pairs))

    def check(self, jobs: Iterable[tuple[Location, Location]]) -> dict[tuple[str, str], FlankCheck]:
        regions, pairs = self.regions(jobs)
        if not pairs:
            return {}
        log.info(
            "validating %d candidate pairs with BLASTN (+-%d bp%s, %d threads)",
            len(pairs),
            self.settings.window,
            ", reciprocal" if self.settings.reciprocal else "",
            self.threads,
        )
        hsps = self.search(regions, pairs)
        return {(a, b): judge(hsps.get((a, b), []), regions[a], regions[b], self.settings) for a, b in pairs}

    def search(
        self, regions: Mapping[str, Region], pairs: list[tuple[str, str]]
    ) -> dict[tuple[str, str], list[Hsp]]:
        names = list(regions)
        index = {name: i for i, name in enumerate(names)}
        wanted: dict[int, set[int]] = {}
        for a, b in pairs:
            wanted.setdefault(index[a], set()).add(index[b])
        queries = sorted(wanted)
        n_chunks = min(len(queries), self.threads * 4)
        chunks = [queries[i::n_chunks] for i in range(n_chunks)]
        dbsize = max(len(r.seq) for r in regions.values())
        found: dict[tuple[str, str], list[Hsp]] = {}
        done = 0
        self.progress("validate", 0, len(queries))
        with tempfile.TemporaryDirectory(prefix="spot_blastn_") as tmp:
            db = Path(tmp, "regions")
            with open(Path(tmp, "regions.fa"), "w") as fh:
                for name in names:
                    fh.write(f">r{index[name]}\n{regions[name].seq}\n")
            run(
                [
                    "makeblastdb",
                    "-dbtype",
                    "nucl",
                    "-parse_seqids",
                    "-in",
                    Path(tmp, "regions.fa"),
                    "-out",
                    db,
                ]
            )
            for k, chunk in enumerate(chunks):
                with open(Path(tmp, f"q{k}.fa"), "w") as fh:
                    for i in chunk:
                        fh.write(f">r{i}\n{regions[names[i]].seq}\n")
                # each chunk is only searched against the candidates of its own loci
                subjects = sorted({s for i in chunk for s in wanted[i]})
                Path(tmp, f"s{k}.txt").write_text("".join(f"r{s}\n" for s in subjects))

            def one(k: int) -> tuple[int, str]:
                out = run(
                    [
                        "blastn",
                        "-query",
                        Path(tmp, f"q{k}.fa"),
                        "-db",
                        db,
                        "-seqidlist",
                        Path(tmp, f"s{k}.txt"),
                        "-dbsize",
                        str(dbsize),
                        "-max_target_seqs",
                        str(len(names)),
                        "-num_threads",
                        "1",
                        "-outfmt",
                        "6 qseqid sseqid pident length qstart qend sstart send",
                        *self.options,
                    ]
                )
                return len(chunks[k]), out

            with ThreadPoolExecutor(max_workers=self.threads) as pool:
                for n, out in pool.map(one, range(len(chunks))):
                    for line in out.splitlines():
                        f = line.split("\t")
                        if len(f) < 8:
                            continue
                        q, s = int(f[0][1:]), int(f[1][1:])
                        if s not in wanted.get(q, ()):
                            continue
                        h = Hsp(float(f[2]), int(f[3]), int(f[4]), int(f[5]), int(f[6]), int(f[7]))
                        if h.pident >= self.identity and h.length >= self.settings.min_hsp:
                            found.setdefault((names[q], names[s]), []).append(h)
                    done += n
                    self.progress("validate", done, len(queries))
        return found
