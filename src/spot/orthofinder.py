from __future__ import annotations

import logging
import os
import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from .config import Config, OrthoFinderInput, Species
from .errors import InputError
from .fasta import read_fasta

log = logging.getLogger("spot")

ORTHOGROUPS = Path("Orthogroups", "Orthogroups.tsv")
UNASSIGNED = Path("Orthogroups", "Orthogroups_UnassignedGenes.tsv")
MAX_DEPTH = 4
PROTEIN_SEARCH = {"blast", "diamond", "diamond_ultra_sens", "mmseqs"}


@dataclass(frozen=True)
class Gene:
    id: str
    species: str
    orthogroup: str
    assigned: bool
    description: str = ""


@dataclass(frozen=True)
class Orthogroup:
    id: str
    members: Mapping[str, tuple[str, ...]]

    @property
    def size(self) -> int:
        return sum(len(ids) for ids in self.members.values())


@dataclass
class OrthoFinderRun:
    name: str
    path: Path
    columns: dict[str, str]
    version: str | None
    dna: bool | None
    search: str | None
    orthogroups: dict[str, Orthogroup]
    genes: dict[str, Gene]
    trees: frozenset[str]
    resolved_trees: frozenset[str]

    @property
    def species(self) -> tuple[str, ...]:
        return tuple(self.columns)

    def shared(self) -> list[Orthogroup]:
        return [og for og in self.orthogroups.values() if all(og.members[s] for s in self.species)]

    def only(self, species: str) -> list[Orthogroup]:
        return [
            og
            for og in self.orthogroups.values()
            if og.members[species] and not any(og.members[s] for s in self.species if s != species)
        ]

    def unassigned(self, species: str) -> list[Gene]:
        return [g for g in self.genes.values() if g.species == species and not g.assigned]

    def sequences(self, og: str) -> Path:
        return self.path / "Orthogroup_Sequences" / f"{og}.fa"

    def msa(self, og: str) -> Path:
        return self.path / "MultipleSequenceAlignments" / f"{og}.fa"

    def tree(self, og: str, resolved: bool = False) -> Path | None:
        if og not in (self.resolved_trees if resolved else self.trees):
            return None
        return self.path / ("Resolved_Gene_Trees" if resolved else "Gene_Trees") / f"{og}_tree.txt"

    def summary(self) -> dict[str, int]:
        shared = self.shared()
        counts = {"orthogroups": len(self.orthogroups), "shared": len(shared)}
        for sp in self.species:
            counts[f"{sp}_only"] = len(self.only(sp))
        for sp in self.species:
            counts[f"{sp}_unassigned"] = len(self.unassigned(sp))
        counts["trees"] = len(self.trees)
        counts["shared_with_tree"] = sum(og.id in self.trees for og in shared)
        return counts


def discover(root: Path) -> dict[str, Path]:
    if not root.is_dir():
        raise InputError(f"OrthoFinder directory not found: {root}")
    if (root / ORTHOGROUPS).is_file():
        return {root.name: root}
    found: dict[str, list[Path]] = {}
    for dirpath, dirnames, _ in os.walk(root):
        here = Path(dirpath)
        rel = here.relative_to(root)
        if (here / ORTHOGROUPS).is_file():
            found.setdefault(rel.parts[0], []).append(here)
            dirnames[:] = []
        elif len(rel.parts) >= MAX_DEPTH:
            dirnames[:] = []
        else:
            dirnames[:] = sorted(d for d in dirnames if not d.startswith("."))
    if not found:
        raise InputError(f"no OrthoFinder results (Orthogroups/Orthogroups.tsv) under {root}")
    runs = {}
    for name, paths in sorted(found.items()):
        latest = max(paths, key=lambda p: (p / ORTHOGROUPS).stat().st_mtime)
        if len(paths) > 1:
            log.warning("%s: %d OrthoFinder results, using %s", name, len(paths), latest.relative_to(root))
        runs[name] = latest
    return runs


def select(source: OrthoFinderInput) -> dict[str, Path]:
    if not source.runs:
        return discover(source.root)
    found = discover(source.root) if any(r.path is None for r in source.runs.values()) else {}
    runs = {}
    for name, run in source.runs.items():
        path = run.path or found.get(name)
        if path is None:
            available = ", ".join(found) or "none"
            raise InputError(f"run {name!r} not found under {source.root} (found: {available})")
        if not (path / ORTHOGROUPS).is_file():
            raise InputError(f"{path} is not an OrthoFinder results directory")
        runs[name] = path
    return runs


def load_all(config: Config) -> list[OrthoFinderRun]:
    return [load(name, path, config.species) for name, path in select(config.orthofinder).items()]


def load(name: str, path: Path, species: Mapping[str, Species]) -> OrthoFinderRun:
    table = path / ORTHOGROUPS
    header, rows = _read_table(table)
    cols = _columns(header, species, table)
    columns = {sp: header[i] for sp, i in cols.items()}

    entries: dict[str, tuple[str, str, bool]] = {}
    orthogroups: dict[str, Orthogroup] = {}

    def add(gene: str, sp: str, og: str, assigned: bool) -> None:
        if gene in entries:
            raise InputError(f"{name}: {gene} is in {entries[gene][1]} and {og}")
        entries[gene] = (sp, og, assigned)

    for row in rows:
        members = {sp: _split(row[i]) for sp, i in cols.items()}
        orthogroups[row[0]] = Orthogroup(row[0], members)
        for sp, ids in members.items():
            for gene in ids:
                add(gene, sp, row[0], True)

    if (path / UNASSIGNED).is_file():
        uheader, urows = _read_table(path / UNASSIGNED)
        for sp, i in _columns(uheader, species, path / UNASSIGNED).items():
            for row in urows:
                for gene in _split(row[i]):
                    add(gene, sp, row[0], False)
    else:
        log.warning("%s: %s not found, unassigned genes are ignored", name, UNASSIGNED)

    descriptions, origin, _ = _sequence_ids(path)
    if descriptions:
        missing = [g for g in entries if g not in descriptions]
        if missing:
            raise InputError(
                f"{name}: {len(missing)} gene IDs are not in WorkingDirectory/SequenceIDs.txt "
                f"(e.g. {', '.join(missing[:3])})"
            )
        wrong = [g for g, (sp, _, _) in entries.items() if (origin.get(g) or columns[sp]) != columns[sp]]
        if wrong:
            examples = ", ".join(wrong[:3])
            raise InputError(f"{name}: genes listed under the wrong species column (e.g. {examples})")
    else:
        log.warning("%s: WorkingDirectory/SequenceIDs.txt not found, FASTA headers unavailable", name)

    genes = {
        g: Gene(g, sp, og, assigned, descriptions.get(g, "")) for g, (sp, og, assigned) in entries.items()
    }
    version, dna, search = _log_info(path)
    if version and not version.startswith("2."):
        log.warning("%s: OrthoFinder %s; only 2.5 output has been tested", name, version)
    if dna and search in PROTEIN_SEARCH:
        log.warning(
            "%s: DNA input was searched with %s, a protein aligner; use -S blast_nucl (or no -S) with -d",
            name,
            search,
        )
    return OrthoFinderRun(
        name=name,
        path=path,
        columns=columns,
        version=version,
        dna=dna,
        search=search,
        orthogroups=orthogroups,
        genes=genes,
        trees=_listing(path / "Gene_Trees", "_tree.txt"),
        resolved_trees=_listing(path / "Resolved_Gene_Trees", "_tree.txt"),
    )


def _read_table(path: Path) -> tuple[list[str], list[list[str]]]:
    with path.open() as fh:
        header = fh.readline().rstrip("\r\n").split("\t")
        rows = []
        for line in fh:
            if not line.strip():
                continue
            row = line.rstrip("\r\n").split("\t")
            rows.append(row + [""] * (len(header) - len(row)))
    return header, rows


def _columns(header: list[str], species: Mapping[str, Species], source: Path) -> dict[str, int]:
    found: dict[str, int] = {}
    for i, column in enumerate(header[1:], start=1):
        hits = [key for key, sp in species.items() if sp.column_regex.search(column)]
        if not hits:
            patterns = ", ".join(f"{k}: {s.match}" for k, s in species.items())
            raise InputError(f"{source}: column {column!r} matches no species ({patterns})")
        if len(hits) > 1:
            raise InputError(f"{source}: column {column!r} matches species {', '.join(hits)}")
        if hits[0] in found:
            raise InputError(f"{source}: species {hits[0]} matches more than one column")
        found[hits[0]] = i
    missing = [k for k in species if k not in found]
    if missing:
        raise InputError(f"{source}: no column for species {', '.join(missing)}")
    return {k: found[k] for k in species}


def _split(cell: str) -> tuple[str, ...]:
    return tuple(gene for gene in (x.strip() for x in cell.split(",")) if gene)


def load_sequences(run: OrthoFinderRun) -> dict[str, str]:
    wd = run.path / "WorkingDirectory"
    _, _, internal = _sequence_ids(run.path)
    seqs: dict[str, str] = {}
    if internal:
        for f in sorted(wd.glob("Species*.fa")):
            if re.fullmatch(r"Species\d+\.fa", f.name):
                for name, seq in read_fasta(f):
                    gene = internal.get(name.split()[0])
                    if gene in run.genes:
                        seqs[gene] = seq
    if len(seqs) < len(run.genes):
        folder = run.path / "Orthogroup_Sequences"
        if folder.is_dir():
            for f in sorted(folder.glob("*.fa")):
                if not f.name.startswith("."):
                    for name, seq in read_fasta(f):
                        gene = name.split()[0]
                        if gene in run.genes:
                            seqs.setdefault(gene, seq)
    missing = [g for g in run.genes if g not in seqs]
    if missing:
        raise InputError(f"{run.name}: no sequence for {len(missing)} genes (e.g. {', '.join(missing[:3])})")
    return seqs


def _sequence_ids(path: Path) -> tuple[dict[str, str], dict[str, str], dict[str, str]]:
    wd = path / "WorkingDirectory"
    seq_ids, species_ids = wd / "SequenceIDs.txt", wd / "SpeciesIDs.txt"
    if not (seq_ids.is_file() and species_ids.is_file()):
        return {}, {}, {}
    files = {}
    for line in species_ids.read_text().splitlines():
        key, sep, name = line.partition(": ")
        if sep and not key.startswith("#"):
            files[key.strip()] = Path(name.strip()).stem
    descriptions, origin, internal = {}, {}, {}
    with seq_ids.open() as fh:
        for line in fh:
            key, sep, header = line.rstrip("\r\n").partition(": ")
            if not sep:
                continue
            parts = header.split(None, 1)
            if not parts:
                continue
            descriptions[parts[0]] = parts[1] if len(parts) > 1 else ""
            origin[parts[0]] = files.get(key.split("_")[0], "")
            internal[key.strip()] = parts[0]
    return descriptions, origin, internal


def _listing(directory: Path, suffix: str) -> frozenset[str]:
    if not directory.is_dir():
        return frozenset()
    with os.scandir(directory) as it:
        return frozenset(
            e.name.removesuffix(suffix) for e in it if e.name.endswith(suffix) and not e.name.startswith(".")
        )


def _log_info(path: Path) -> tuple[str | None, bool | None, str | None]:
    log_file = path / "Log.txt"
    if not log_file.is_file():
        return None, None, None
    text = log_file.read_text(errors="replace")
    version = re.search(r"OrthoFinder version (\S+)", text)
    command = re.search(r"^Command Line: (.+)$", text, re.MULTILINE)
    if not command:
        return (version.group(1) if version else None), None, None
    args = command.group(1).split()
    dna = "-d" in args or "--dna" in args
    search = "blast_nucl" if dna else "diamond"
    for flag, value in zip(args, args[1:], strict=False):
        if flag in ("-S", "--search"):
            search = value
    return (version.group(1) if version else None), dna, search
