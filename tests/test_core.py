from pathlib import Path

import pytest

from spot.anchors import AnchorIndex
from spot.cache import Checkpoints, digest, key_of
from spot.config import Anchors
from spot.locate import Genome, _place, identity
from spot.positional import GeneIndex, window_tolerance
from spot.positional import test as positional_test
from spot.syri import AlignedIndex, Alignment, homologous_pairs


def aln(r: str, rs: int, re: int, q: str, qs: int, qe: int, kind: str = "SYN") -> Alignment:
    return Alignment(r, rs, re, q, qs, qe, kind, f"{kind}1")


def test_homology_keeps_split_chromosomes() -> None:
    alignments = [
        aln("chr1", 1, 1_000_000, "A", 1, 1_000_000),
        aln("chr1", 1_000_001, 1_600_000, "B", 1, 600_000),
        aln("chr1", 1_600_001, 1_620_000, "C", 1, 20_000),
        aln("chr2", 1, 900_000, "D", 1, 900_000),
    ]
    assert homologous_pairs(alignments, 0.1) == {("chr1", "A"), ("chr1", "B"), ("chr2", "D"), ("chr1", "C")}
    assert homologous_pairs(alignments, 1.0) == {("chr1", "A"), ("chr2", "D")}


def test_aligned_index_merges() -> None:
    idx = AlignedIndex([aln("c", 1, 100, "q", 1, 100), aln("c", 50, 200, "q", 300, 150)], "query")
    assert idx.covered("q", 1, 300) == 251
    ref = AlignedIndex([aln("c", 1, 100, "q", 1, 100), aln("c", 50, 200, "q", 300, 150)], "reference")
    assert ref.covered("c", 150, 400) == 51
    assert ref.covered("zz", 1, 10) == 0


def test_empty_locus_in_unaligned_gap() -> None:
    alignments = [aln("c", 1, 10_000, "c", 1, 10_000), aln("c", 20_001, 30_000, "c", 40_001, 50_000)]
    anchors = AnchorIndex.build(alignments, "reference", ["SYN"], None)
    others = GeneIndex([])
    loc = _loc("x", 12_000, 12_100)
    settings = Anchors(tol=1000, tol_frac=0.0)
    t = positional_test(loc, set(), anchors, others, settings, AlignedIndex(alignments, "query"))
    assert t.klass == "empty_locus_gap" and t.window_aligned is not None and t.window_aligned < 0.5
    assert not t.evaluable


def test_window_tolerance() -> None:
    assert window_tolerance(5_000, Anchors(tol=2000, tol_frac=0.1)) == 2000
    assert window_tolerance(40_000, Anchors(tol=2000, tol_frac=0.1)) == 4000


def _loc(gene: str, start: int, end: int):  # type: ignore[no-untyped-def]
    from spot.locate import Location

    return Location(gene, "A", "located", "c", start, end, "+")


def test_identity_and_statuses(tmp_path: Path) -> None:
    (tmp_path / "g.fa").write_text(">chr1\n" + "ACGT" * 50 + "\n")
    g = Genome(tmp_path / "g.fa")
    seq = ("ACGT" * 25)[:100]
    assert _place("a", "A", ("chr1", 1, 100), None, g, seq).strand == "+"
    one_off = "T" + seq[1:]
    assert _place("b", "A", ("chr1", 1, 100), None, g, one_off).status == "sequence_mismatch"
    placed = _place("b", "A", ("chr1", 1, 100), None, g, one_off, 0.99)
    assert placed.located and placed.identity == pytest.approx(0.99)
    assert _place("c", "A", ("scaf9", 1, 100), None, g, seq).status == "off_assembly"
    assert _place("d", "A", ("chr1", 150, 300), None, g, seq).status == "out_of_range"
    assert identity("ACGN", "ACGN") == 0.75
    assert identity("AC", "ACG") == 0.0


def test_checkpoints_and_digest(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("SPOT_CACHE", str(tmp_path / "cache"))
    f = tmp_path / "x.txt"
    f.write_text("hello")
    assert digest(f) == digest(f) == "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824"
    calls: list[int] = []

    def compute() -> int:
        calls.append(1)
        return 42

    cp = Checkpoints(tmp_path / "out")
    assert cp.stage("s", key_of("a"), compute) == 42
    assert cp.stage("s", key_of("a"), compute) == 42
    assert cp.stage("s", key_of("b"), compute) == 42
    assert len(calls) == 2
    assert Checkpoints(tmp_path / "out", enabled=False).stage("s", key_of("a"), compute) == 42
    assert len(calls) == 3
