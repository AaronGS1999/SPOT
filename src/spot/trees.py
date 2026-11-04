"""Concordance with the orthologues that OrthoFinder infers from its gene trees.

OrthoFinder writes, for every pair of species, the orthologues found by reconciling each gene tree with the
species tree (``Orthologues/Orthologues_<A>/<A>__v__<B>.tsv``). A pair that SPOT confirms by position and
that the tree also calls orthologous is supported by two independent kinds of evidence; a pair the tree
calls paralogous points to a duplication before the split of the species, or to an error in the tree.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

from .orthofinder import OrthoFinderRun

log = logging.getLogger("spot")

ORTHOLOG = "ortholog"
NOT_ORTHOLOG = "not_ortholog"


def orthologue_table(run: OrthoFinderRun, reference: str, query: str) -> Path:
    a, b = run.columns[reference], run.columns[query]
    return run.path / "Orthologues" / f"Orthologues_{a}" / f"{a}__v__{b}.tsv"


def orthologues(run: OrthoFinderRun, reference: str, query: str) -> set[tuple[str, str]] | None:
    """(reference, query) gene pairs called orthologous by the gene trees, or None if the table is missing."""
    path = orthologue_table(run, reference, query)
    if not path.is_file():
        log.warning("%s: %s not found, tree concordance is not reported", run.name, path)
        return None
    pairs: set[tuple[str, str]] = set()
    with path.open(newline="") as fh:
        rows = csv.reader(fh, delimiter="\t")
        next(rows, None)
        for row in rows:
            if len(row) < 3:
                continue
            refs = [g.strip() for g in row[1].split(",") if g.strip()]
            qrys = [g.strip() for g in row[2].split(",") if g.strip()]
            pairs.update((r, q) for r in refs for q in qrys)
    return pairs


def concordance(pair: tuple[str, str], called: set[tuple[str, str]] | None) -> str:
    if called is None:
        return ""
    return ORTHOLOG if pair in called else NOT_ORTHOLOG
