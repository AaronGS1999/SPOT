import random
import shutil
from pathlib import Path

import pytest

from spot import config, pipeline
from spot.fasta import revcomp, write_fasta

pytestmark = pytest.mark.skipif(shutil.which("blastn") is None, reason="blastn not on PATH")

GENES = {
    "AA_1": ("chr1", 20_001, 20_400),
    "AA_2": ("chr1", 40_001, 40_300),
    "AA_3": ("chr1", 50_001, 50_200),
    "BB_1": ("chr1", 20_001, 20_400),
    "BB_2": ("chr1", 8_001, 8_300),
}
ORTHOGROUPS = (
    "Orthogroup\tAA_set\tBB_set\r\nOG0000000\tAA_1\tBB_1\r\nOG0000001\tAA_2\tBB_2\r\nOG0000002\tAA_3\t\r\n"
)
SYRI = (
    "chr1\t1\t600000\t-\t-\tchr1\t1\t600000\tSYN1\t-\tSYN\t-\n"
    "chr1\t1\t600000\t-\t-\tchr1\t1\t600000\tSYNAL1\tSYN1\tSYNAL\t-\n"
)


@pytest.fixture
def project(tmp_path: Path) -> Path:
    rng = random.Random(7)
    chrom = "".join(rng.choice("ACGT") for _ in range(600_000))
    (tmp_path / "genomes").mkdir()
    write_fasta(tmp_path / "genomes" / "a.fa", [("chr1", chrom)])
    write_fasta(tmp_path / "genomes" / "b.fa", [("chr1", chrom)])

    res = tmp_path / "of" / "run1" / "Results_Nov01"
    (res / "Orthogroups").mkdir(parents=True)
    (res / "Orthogroups" / "Orthogroups.tsv").write_text(ORTHOGROUPS)
    (res / "Orthogroups" / "Orthogroups_UnassignedGenes.tsv").write_text("Orthogroup\tAA_set\tBB_set\r\n")
    (res / "Log.txt").write_text("Started OrthoFinder version 2.5.5\nCommand Line: orthofinder -f x -d\n")
    wd = res / "WorkingDirectory"
    wd.mkdir()
    (wd / "SpeciesIDs.txt").write_text("0: AA_set.fa\n1: BB_set.fa\n")
    ids, species = [], {"0": [], "1": []}
    for i, (gene, (chr_, start, end)) in enumerate(GENES.items()):
        sp = "0" if gene.startswith("AA") else "1"
        internal = f"{sp}_{i}"
        ids.append(f"{internal}: {gene} {chr_}:{start}-{end}")
        seq = chrom[start - 1 : end]
        species[sp].append((internal, revcomp(seq) if gene == "BB_1" else seq))
    (wd / "SequenceIDs.txt").write_text("\n".join(ids) + "\n")
    for sp, records in species.items():
        write_fasta(wd / f"Species{sp}.fa", records)
    (tmp_path / "syri.out").write_text(SYRI)

    (tmp_path / "spot.yaml").write_text(
        """orthofinder:
  root: of
  runs:
    run1: {label: Run 1, unit: genes}
species:
  AA: {label: Species A, match: "^AA_", genome: genomes/a.fa}
  BB: {label: Species B, match: "^BB_", genome: genomes/b.fa}
reference: AA
query: BB
syri_out: syri.out
outdir: out
threads: 2
dotplots: {pairs: false}
"""
    )
    return tmp_path


def test_run(project: Path) -> None:
    analysis = pipeline.run(config.load(project / "spot.yaml"))
    (result,) = analysis.results
    assert result.status == {
        "OG0000000": "syntenic_ortholog",
        "OG0000001": "evaluable_not_syntenic",
    }
    assert result.pairs["OG0000000"] == [("AA_1", "BB_1")]
    assert analysis.locations["BB_1"].strand == "-"
    hit = result.tests["AA_1"]
    assert hit.p_value is not None and hit.q_value is not None and hit.q_value <= 0.05
    assert result.calibration is not None and result.calibration.observed == 2

    out = project / "out" / "run1"
    assert (out / "syntenic_ortholog" / "OG0000000" / "OG0000000.pairs.fa").is_file()
    assert (out / "AA_only" / "OG0000002" / "OG0000002.fa").is_file()
    rows = (out / "orthogroups.tsv").read_text().splitlines()
    assert rows[0].split("\t")[:3] == ["orthogroup", "zone", "status"]
    assert len(rows) == 4
    for name in ("spot.log", "run.yaml", "versions.tsv", "genome/dotplot.png"):
        assert (project / "out" / name).is_file(), name
