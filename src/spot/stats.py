"""Statistics of the positional test.

P-value. A member of the same orthogroup could fall in the projected window by chance, more easily for
large orthogroups and long windows. Under the null hypothesis the window is placed uniformly at random on
the aligned part of the chromosomes homologous to the tested sequence; it then contains a member with
probability equal to the length of the positions from which it would overlap one (each member widened by
the window length, overlaps merged) divided by the aligned length. This is the exact value of a
permutation test that relocates the window, computed without sampling. The p-values of a run are adjusted
with the Benjamini-Hochberg procedure (SciPy) and a test is kept when q <= alpha.

Calibration. The orthogroup labels of each species are shuffled within the run, keeping the orthogroup
sizes and every position, and the positional test is repeated. The mean number of hits over the
permutations, divided by the observed number, estimates the share of hits expected by chance.
"""

from __future__ import annotations

import zlib
from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass

import numpy as np
from scipy.stats import false_discovery_control

from .locate import Location
from .positional import Test
from .syri import AlignedIndex


def union_length(spans: Iterable[tuple[int, int]]) -> int:
    """Total length of 1-based inclusive spans after merging overlaps."""
    total, cur_s, cur_e = 0, 0, -1
    for s, e in sorted(spans):
        if s > cur_e + 1:
            total += cur_e - cur_s + 1 if cur_e >= cur_s else 0
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    return total + (cur_e - cur_s + 1 if cur_e >= cur_s else 0)


def window_pvalue(width: int, members: Iterable[tuple[int, int]], space: int) -> float:
    """Probability that a ``width`` bp window placed uniformly in ``space`` aligned bp overlaps a member."""
    if space <= 0 or width >= space:
        return 1.0
    hit = union_length((max(1, s - width + 1), e) for s, e in members)
    return min(1.0, hit / space)


class NullSpace:
    """Where projected windows can land in one species: aligned bases of the chromosomes homologous to
    each chromosome of the other species."""

    def __init__(self, aligned: AlignedIndex, partners: Mapping[str, set[str]] | None) -> None:
        """``partners`` None: windows may land on any chromosome (``anchors.same_chromosome`` off)."""
        self.aligned = aligned
        self.partners = partners
        self._all = set(aligned.chromosomes())

    def targets(self, chrom: str) -> set[str]:
        return self._all if self.partners is None else self.partners.get(chrom, set())

    def length(self, chrom: str) -> int:
        return sum(self.aligned.total(c) for c in self.targets(chrom))


def partners(homologous: set[tuple[str, str]], side: str) -> dict[str, set[str]]:
    """Chromosomes of the other species homologous to each chromosome of ``side``."""
    out: dict[str, set[str]] = defaultdict(set)
    for ref, qry in homologous:
        if side == "reference":
            out[ref].add(qry)
        else:
            out[qry].add(ref)
    return out


def pvalue(test: Test, loc: Location, same_group: Iterable[Location], space: NullSpace) -> float:
    targets = space.targets(loc.chrom)
    members = [(m.start, m.end) for m in same_group if m.located and m.chrom in targets]
    width = (test.window_end or 0) - (test.window_start or 0) + 1
    return window_pvalue(width, members, space.length(loc.chrom))


def adjust(pvalues: list[float]) -> list[float]:
    """Benjamini-Hochberg q-values."""
    if not pvalues:
        return []
    return [float(q) for q in false_discovery_control(np.asarray(pvalues), method="bh")]


@dataclass(frozen=True)
class Calibration:
    run: str
    observed: int
    permutations: int
    null_mean: float
    null_sd: float
    null_p95: float

    @property
    def expected_false(self) -> float:
        """Share of the observed hits expected by chance (null mean / observed, at most 1)."""
        return min(1.0, self.null_mean / self.observed) if self.observed else 0.0


def calibrate(
    run_name: str,
    labels: Mapping[str, str],
    windows: Mapping[str, list[str]],
    permutations: int,
    seed: int,
) -> Calibration:
    """Shuffle orthogroup labels within each species and count positional hits.

    ``labels`` maps every assigned sequence to ``species/orthogroup``; ``windows`` maps each tested
    sequence to the sequences found in its projected window.
    """
    genes = list(labels)
    species = np.array([labels[g].split("/", 1)[0] for g in genes])
    groups = np.array([labels[g] for g in genes], dtype=object)
    index = {g: i for i, g in enumerate(genes)}
    tested = [
        (index[g], [index[o] for o in found if o in index]) for g, found in windows.items() if g in index
    ]

    def hits(assigned: np.ndarray) -> int:
        og = np.array([a.split("/", 1)[1] for a in assigned], dtype=object)
        return sum(1 for i, found in tested if any(og[j] == og[i] for j in found))

    observed = hits(groups)
    rng = np.random.default_rng(seed + zlib.crc32(run_name.encode()))
    null = []
    for _ in range(permutations):
        shuffled = groups.copy()
        for sp in np.unique(species):
            idx = np.flatnonzero(species == sp)
            shuffled[idx] = groups[rng.permutation(idx)]
        null.append(hits(shuffled))
    values = np.asarray(null, dtype=float) if null else np.zeros(1)
    return Calibration(
        run_name,
        observed,
        permutations,
        float(values.mean()),
        float(values.std(ddof=1)) if len(values) > 1 else 0.0,
        float(np.percentile(values, 95)),
    )
