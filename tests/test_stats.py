from pathlib import Path

import numpy as np
import pytest

from spot import stats, trees
from spot.locate import Location
from spot.orthofinder import OrthoFinderRun
from spot.positional import Test as PositionalTest
from spot.syri import AlignedIndex, Alignment


def aln(r: str, rs: int, re: int, q: str, qs: int, qe: int) -> Alignment:
    return Alignment(r, rs, re, q, qs, qe, "SYN", "SYN1")


def test_union_length_merges_overlaps() -> None:
    assert stats.union_length([]) == 0
    assert stats.union_length([(1, 10), (5, 20), (30, 30)]) == 21
    assert stats.union_length([(1, 10), (11, 20)]) == 20


def test_window_pvalue_matches_placement_by_simulation() -> None:
    space, width, members = 100_000, 2_000, [(10_000, 10_150), (11_000, 11_200), (60_000, 60_100)]
    p = stats.window_pvalue(width, members, space)
    rng = np.random.default_rng(0)
    starts = rng.integers(1, space - width + 2, 200_000)
    ends = starts + width - 1
    hit = np.zeros(starts.size, bool)
    for s, e in members:
        hit |= (starts <= e) & (ends >= s)
    assert p == pytest.approx(hit.mean(), abs=0.003)
    assert stats.window_pvalue(width, members, 0) == 1.0
    assert stats.window_pvalue(10**6, members, space) == 1.0


def test_pvalue_uses_homologous_chromosomes_only() -> None:
    aligned = AlignedIndex(
        [aln("a1", 1, 50_000, "b1", 1, 50_000), aln("a2", 1, 50_000, "b2", 1, 50_000)], "query"
    )
    space = stats.NullSpace(aligned, stats.partners({("a1", "b1"), ("a2", "b2")}, "reference"))
    test = PositionalTest("x", "A", "syntenic_ortholog", window_start=1_000, window_end=1_999)
    loc = Location("x", "A", "located", "a1", 1_000, 1_100, "+")
    members = [
        Location("y", "B", "located", "b1", 1_500, 1_600, "+"),
        Location("z", "B", "located", "b2", 10, 20, "+"),
    ]
    assert stats.pvalue(test, loc, members, space) == pytest.approx(1_100 / 50_000)
    anywhere = stats.NullSpace(aligned, None)
    assert anywhere.length("a1") == 100_000


def test_benjamini_hochberg() -> None:
    q = stats.adjust([0.01, 0.04, 0.03, 1.0])
    assert q == pytest.approx([0.04, 0.04 * 4 / 3, 0.04 * 4 / 3, 1.0])
    assert stats.adjust([]) == []


def test_calibration_is_reproducible_and_low_for_real_signal() -> None:
    labels = {f"a{i}": f"A/OG{i}" for i in range(50)} | {f"b{i}": f"B/OG{i}" for i in range(50)}
    windows = {f"a{i}": [f"b{i}", f"b{(i + 1) % 50}"] for i in range(50)}
    one = stats.calibrate("run", labels, windows, 50, seed=3)
    assert one == stats.calibrate("run", labels, windows, 50, seed=3)
    assert one.observed == 50
    assert one.null_mean < 5 and one.expected_false < 0.1


def test_tree_concordance(tmp_path: Path) -> None:
    folder = tmp_path / "Orthologues" / "Orthologues_A_x"
    folder.mkdir(parents=True)
    (folder / "A_x__v__B_x.tsv").write_text("Orthogroup\tA_x\tB_x\nOG1\ta1, a2\tb1\nOG2\ta3\tb3\n")
    run = OrthoFinderRun(
        "r", tmp_path, {"A": "A_x", "B": "B_x"}, None, None, None, {}, {}, frozenset(), frozenset()
    )
    called = trees.orthologues(run, "A", "B")
    assert called == {("a1", "b1"), ("a2", "b1"), ("a3", "b3")}
    assert trees.concordance(("a2", "b1"), called) == trees.ORTHOLOG
    assert trees.concordance(("a1", "b3"), called) == trees.NOT_ORTHOLOG
    run.columns = {"A": "A_x", "B": "C_x"}
    assert trees.orthologues(run, "A", "B") is None
    assert trees.concordance(("a1", "b1"), None) == ""
