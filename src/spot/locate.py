from __future__ import annotations

import csv
import logging
import re
from collections import Counter
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path

import pyfaidx

from .config import Species
from .errors import InputError
from .fasta import revcomp
from .lift import Lift
from .orthofinder import OrthoFinderRun

log = logging.getLogger("spot")

LOCATED = "located"
OFF_ASSEMBLY = (
    "off_assembly"  # the chromosome or scaffold is not in the genome FASTA (e.g. unplaced scaffolds)
)
OUT_OF_RANGE = "out_of_range"  # the chromosome exists but the coordinates fall outside it
MISMATCH = "sequence_mismatch"


@dataclass(frozen=True)
class Location:
    gene: str
    species: str
    status: str
    chrom: str = ""
    start: int = 0
    end: int = 0
    strand: str = ""
    identity: float | None = None

    @property
    def located(self) -> bool:
        return self.status == LOCATED


class Genome:
    def __init__(self, path: Path) -> None:
        if not path.is_file():
            raise InputError(f"genome not found: {path}")
        try:
            self._fasta = pyfaidx.Fasta(str(path), as_raw=True, sequence_always_upper=True)
        except (OSError, pyfaidx.FastaIndexingError) as e:
            raise InputError(f"cannot index {path}: {e}") from None
        self.path = path
        self.lengths = {name: len(record) for name, record in self._fasta.items()}

    def fetch(self, chrom: str, start: int, end: int) -> str:
        return str(self._fasta[chrom][start - 1 : end])


def coordinates(
    species: Species, genes: Iterable[tuple[str, str]], table: Mapping[str, tuple[str, int, int]] | None
) -> dict[str, tuple[str, int, int]]:
    found = {}
    if table is not None:
        for gene, _ in genes:
            if gene in table:
                found[gene] = table[gene]
        return found
    pattern = re.compile(species.header_regex)
    for gene, description in genes:
        m = pattern.search(description)
        if m:
            found[gene] = (m["chr"], int(m["start"]), int(m["end"]))
    return found


def read_table(path: Path) -> dict[str, tuple[str, int, int]]:
    try:
        with open(path, newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
    except OSError as e:
        raise InputError(f"cannot read {path}: {e.strerror}") from None
    if rows and not {"id", "chr", "start", "end"} <= set(rows[0]):
        raise InputError(f"{path}: expected columns id, chr, start, end")
    try:
        return {r["id"]: (r["chr"], int(r["start"]), int(r["end"])) for r in rows}
    except ValueError:
        raise InputError(f"{path}: non-numeric coordinates") from None


def locate(
    runs: list[OrthoFinderRun],
    species: Mapping[str, Species],
    sequences: Mapping[str, Mapping[str, str]],
    genomes: Mapping[str, Genome],
) -> dict[str, Location]:
    tables = {k: read_table(sp.coords) for k, sp in species.items() if isinstance(sp.coords, Path)}
    lifts = {
        k: Lift.read(sp.lift, sp.strip_suffix) for k, sp in species.items() if sp.lift or sp.strip_suffix
    }
    locations: dict[str, Location] = {}
    for run in runs:
        for key, sp in species.items():
            genes = [
                (g.id, g.description)
                for g in run.genes.values()
                if g.species == key and g.id not in locations
            ]
            coords = coordinates(sp, genes, tables.get(key))
            for gene, _ in genes:
                locations[gene] = _place(
                    gene,
                    key,
                    coords.get(gene),
                    lifts.get(key),
                    genomes[key],
                    sequences[run.name][gene],
                    sp.min_identity,
                )
    for key in species:
        counts = Counter(loc.status for loc in locations.values() if loc.species == key)
        log.info("%s: %s", key, ", ".join(f"{status} {n}" for status, n in sorted(counts.items())))
    return locations


def _place(
    gene: str,
    species: str,
    coords: tuple[str, int, int] | None,
    lift: Lift | None,
    genome: Genome,
    seq: str,
    min_identity: float = 1.0,
) -> Location:
    if coords is None:
        return Location(gene, species, "no_coordinates")
    if lift is not None:
        lifted = lift(*coords)
        if lifted is None:
            return Location(gene, species, "lift_failed")
        coords = lifted
    chrom, start, end = coords
    if chrom not in genome.lengths:
        return Location(gene, species, OFF_ASSEMBLY, chrom, start, end)
    if start < 1 or end > genome.lengths[chrom] or end < start:
        return Location(gene, species, OUT_OF_RANGE, chrom, start, end)
    ref = genome.fetch(chrom, start, end)
    seq = seq.upper()
    if ref == seq:
        return Location(gene, species, LOCATED, chrom, start, end, "+", 1.0)
    rc = revcomp(seq)
    if ref == rc:
        return Location(gene, species, LOCATED, chrom, start, end, "-", 1.0)
    plus, minus = identity(ref, seq), identity(ref, rc)
    strand, best = ("+", plus) if plus >= minus else ("-", minus)
    if best >= min_identity:
        return Location(gene, species, LOCATED, chrom, start, end, strand, best)
    return Location(gene, species, MISMATCH, chrom, start, end, identity=best)


def identity(a: str, b: str) -> float:
    """Fraction of identical positions of two sequences of the same length (0 when the lengths differ).

    Ambiguous bases (N) in either sequence count as mismatches, so a soft-masked or gap-filled
    copy can only lower the identity.
    """
    if len(a) != len(b) or not a:
        return 0.0
    same = sum(x == y and x != "N" for x, y in zip(a, b, strict=True))
    return same / len(a)
