import os
import random
import shutil
from pathlib import Path

import pytest

from spot import config, syri
from spot.fasta import revcomp, write_fasta

pytestmark = pytest.mark.skipif(
    shutil.which("syri") is None or shutil.which("minimap2") is None, reason="syri/minimap2 not on PATH"
)


def mutate(seq: str, rate: float, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") if rng.random() < rate else b for b in seq)


def test_align_writes_only_in_its_folder(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    rng = random.Random(11)
    ref = "".join(rng.choice("ACGT") for _ in range(300_000))
    qry = mutate(ref[:100_000] + revcomp(ref[100_000:150_000]) + ref[150_000:], 0.005, rng)
    write_fasta(tmp_path / "r.fa", [("chr1", ref)])
    write_fasta(tmp_path / "q.fa", [("chr1", qry)])
    (tmp_path / "c.yaml").write_text(
        "orthofinder: {root: of}\n"
        "species:\n  R: {label: r, match: R, genome: r.fa}\n  Q: {label: q, match: Q, genome: q.fa}\n"
        "reference: R\nquery: Q\nthreads: 2\n"
    )
    cfg = config.load(tmp_path / "c.yaml")
    elsewhere = tmp_path / "cwd"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    monkeypatch.setenv("SPOT_CACHE", str(tmp_path / "cache"))
    out = syri.align(cfg, tmp_path / "genome")
    assert out.is_file() and (tmp_path / "genome" / "syri.log").is_file()
    assert os.listdir(elsewhere) == []
    alignments, _ = syri.read(out)
    assert {"SYN", "INV"} <= {a.kind for a in alignments}
