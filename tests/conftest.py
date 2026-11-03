from pathlib import Path

import pytest

ORTHOGROUPS = """Orthogroup\tBB_set\tAA_set
OG0000000\tBB_1\tAA_1, AA_2
OG0000001\t\tAA_3, AA_4
OG0000002\tBB_2\tAA_5
OG0000003\tBB_3, BB_4\t
"""

UNASSIGNED = """Orthogroup\tBB_set\tAA_set
OG0000004\t\tAA_6
OG0000005\tBB_5\t
"""

SPECIES_IDS = "0: AA_set.fasta\n1: BB_set.fasta\n"

LOG = """2026-10-27 18:30:15 : Started OrthoFinder version 2.5.5
Command Line: /usr/local/bin/orthofinder -f input -d -S diamond -M msa -A mafft -T fasttree -X -o results
"""


def write_run(results: Path) -> Path:
    (results / "Orthogroups").mkdir(parents=True)
    (results / "Orthogroups" / "Orthogroups.tsv").write_text(ORTHOGROUPS)
    (results / "Orthogroups" / "Orthogroups_UnassignedGenes.tsv").write_text(UNASSIGNED)
    (results / "Log.txt").write_text(LOG)
    wd = results / "WorkingDirectory"
    wd.mkdir()
    (wd / "SpeciesIDs.txt").write_text(SPECIES_IDS)
    lines = [f"0_{i}: AA_{i + 1} chr1:{100 * i + 1}-{100 * i + 80} K=0.1" for i in range(6)]
    lines += [f"1_{i}: BB_{i + 1} chr2:{100 * i + 1}-{100 * i + 90}" for i in range(5)]
    (wd / "SequenceIDs.txt").write_text("\n".join(lines) + "\n")
    trees = results / "Gene_Trees"
    trees.mkdir()
    for og in ("OG0000000", "OG0000003"):
        (trees / f"{og}_tree.txt").write_text("(a,b);\n")
    (trees / "._OG0000002_tree.txt").write_bytes(b"\x00\x05\x16\x07")
    seqs = results / "Orthogroup_Sequences"
    seqs.mkdir()
    for i in range(6):
        (seqs / f"OG000000{i}.fa").write_text(">x\nACGT\n")
    return results


@pytest.fixture
def of_root(tmp_path: Path) -> Path:
    root = tmp_path / "orthofinder"
    write_run(root / "runA" / "results" / "Results_Oct27")
    write_run(root / "runB" / "results" / "Results_Oct27")
    (root / "runA" / "input").mkdir()
    return root


@pytest.fixture
def config_file(tmp_path: Path, of_root: Path) -> Path:
    path = tmp_path / "spot.yaml"
    path.write_text(
        f"""orthofinder:
  root: {of_root.name}
species:
  AA:
    label: Species A
    match: "^AA_"
    genome: genomes/a.fasta
  BB:
    label: Species B
    match: "^BB_"
    genome: /abs/b.fasta
reference: AA
query: BB
"""
    )
    return path
