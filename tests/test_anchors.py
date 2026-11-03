from pathlib import Path

from spot import syri
from spot.anchors import AnchorIndex
from spot.config import Anchors
from spot.locate import Location
from spot.positional import GeneIndex
from spot.positional import test as positional_test

SYRI = """chr1\t1\t999\t-\t-\t-\t-\t-\tNOTAL1\t-\tNOTAL\t-
chr1\t1000\t50000\t-\t-\tchr1\t1000\t50000\tSYN1\t-\tSYN\t-
chr1\t1000\t19000\t-\t-\tchr1\t1000\t19000\tSYNAL1\tSYN1\tSYNAL\t-
chr1\t21000\t50000\t-\t-\tchr1\t21000\t50000\tSYNAL2\tSYN1\tSYNAL\t-
chr1\t60000\t70000\t-\t-\tchr1\t80000\t70000\tINV1\t-\tINV\t-
chr1\t60000\t70000\t-\t-\tchr1\t80000\t70000\tINVAL1\tINV1\tINVAL\t-
chr1\t71000\t72000\t-\t-\tchr2\t5000\t6000\tTRANSAL1\tTRANS1\tTRANSAL\t-
chr1\t500\t501\tA\tG\tchr1\t500\t501\tSNP1\tSYN1\tSNP\t-
"""


def read(tmp_path: Path) -> tuple[list[syri.Alignment], list[syri.Alignment]]:
    path = tmp_path / "syri.out"
    path.write_text(SYRI)
    return syri.read(path)


def test_read_and_homologs(tmp_path: Path) -> None:
    alignments, blocks = read(tmp_path)
    assert [a.kind for a in alignments] == ["SYN", "SYN", "INV", "TRANS"]
    assert [b.block for b in blocks] == ["SYN1", "INV1"]
    assert alignments[0].block == "SYN1"
    assert syri.homologous_pairs(alignments) == {("chr1", "chr1")}


def test_index_filters_and_query_side(tmp_path: Path) -> None:
    alignments, _ = read(tmp_path)
    homologous = syri.homologous_pairs(alignments)
    assert len(AnchorIndex.build(alignments, "reference", ["SYN", "INV", "TRANS"], homologous)) == 3
    assert len(AnchorIndex.build(alignments, "reference", ["SYN", "INV", "TRANS"], None)) == 4
    qry = AnchorIndex.build(alignments, "query", ["INV"], None)
    b = qry.bracket("chr1", 74000, 74100, 50_000)
    assert b is not None
    assert b.left.start == 70000 and b.left.other_at_start == 70000
    assert b.proj_left == 66001 and b.proj_right == 65899


def loc(gene: str, species: str, start: int, end: int, chrom: str = "chr1") -> Location:
    return Location(gene, species, "located", chrom, start, end, "+")


def test_positional_classes(tmp_path: Path) -> None:
    alignments, _ = read(tmp_path)
    anchors = AnchorIndex.build(alignments, "reference", ["SYN", "INV"], syri.homologous_pairs(alignments))
    others = GeneIndex([loc("q1", "B", 19500, 20400), loc("q2", "B", 45000, 45100)])
    settings = Anchors()
    t = positional_test(loc("r1", "A", 19200, 20800), {"q1"}, anchors, others, settings)
    assert t.klass == "syntenic_ortholog" and t.orthologs == ["q1"]
    assert t.evidence == "flank_anchors" and t.block_type == "collinear" and t.blocks == "SYN1"
    assert (
        positional_test(loc("r1", "A", 19200, 20800), {"zz"}, anchors, others, settings).klass
        == "syntenic_locus_other"
    )
    assert (
        positional_test(loc("r2", "A", 30000, 30100), set(), anchors, others, settings).klass
        == "empty_locus_aligned"
    )
    assert positional_test(loc("r3", "A", 100, 200), set(), anchors, others, settings).klass == "no_anchor"
    wide = AnchorIndex.build(alignments, "reference", ["SYN", "INV", "TRANS"], None)
    r4 = loc("r4", "A", 70500, 70600)
    assert positional_test(r4, set(), wide, others, Anchors(gapmax=450)).klass == "anchors_too_far"
    assert positional_test(r4, set(), wide, others, settings).klass == "anchors_split"
    assert (
        positional_test(Location("r5", "A", "off_assembly"), set(), anchors, others, settings).klass
        == "unplaced"
    )
