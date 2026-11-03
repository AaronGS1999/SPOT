from spot.config import Validation
from spot.tools import as_args
from spot.validate import Hsp, Region, colinear, covered, judge, parse_hsps, required

# two loci with 1000 bp flanks and a 100 bp core: [1, 1000] left, [1001, 1100] core, [1101, 2100] right
A = Region("a", "chr1", 1, 2100, "A" * 2100, 1000, 100)
B = Region("b", "chr1", 1, 2100, "C" * 2100, 1000, 100)
STRICT = Validation(min_frac=0.5, min_bp=500, max_offset=200)


def test_required_follows_fraction_minimum_and_half_cap() -> None:
    assert required(10_000, 0.2, 1000) == 2000
    assert required(10_000, 0.0, 1000) == 1000
    assert required(1_500, 0.2, 1000) == 750
    assert required(0, 0.2, 1000) == 0


def test_covered_merges_and_clips() -> None:
    hsps = [(1, 100), (50, 150), (300, 400)]
    assert covered(hsps, 0, 1000) == 251
    assert covered(hsps, 0, 120) == 120
    assert covered(hsps, 350, 1000) == 50
    assert covered([], 0, 10) == 0


def test_direct_flanks_are_colinear() -> None:
    hsps = [Hsp(95, 900, 50, 950, 60, 960), Hsp(95, 900, 1150, 2050, 1140, 2040)]
    assert colinear(hsps, A, B, 200, inverted=False) == (901, 901)
    assert judge(hsps, A, B, STRICT).ok


def test_one_hsp_across_the_core_counts_for_both_flanks() -> None:
    assert colinear([Hsp(99, 2100, 1, 2100, 1, 2100)], A, B, 200, inverted=False) == (1000, 1000)


def test_crossed_flanks_are_rejected() -> None:
    # the left flank of A aligns to the right flank of B in the same orientation: not one locus
    hsps = [Hsp(95, 900, 50, 950, 1150, 2050), Hsp(95, 900, 1150, 2050, 50, 950)]
    assert colinear(hsps, A, B, 200, inverted=False) == (0, 0)
    assert colinear(hsps, A, B, 200, inverted=True) == (0, 0)
    assert not judge(hsps, A, B, STRICT).ok
    # the old coverage mode accepts them
    assert judge(hsps, A, B, Validation(mode="coverage", min_frac=0.5, min_bp=500)).ok


def test_inverted_pair_is_accepted() -> None:
    # B lies in an inversion: A's left flank joins B's right flank on the minus strand
    hsps = [Hsp(95, 900, 50, 950, 2050, 1150), Hsp(95, 900, 1150, 2050, 950, 50)]
    assert colinear(hsps, A, B, 200, inverted=True) == (901, 901)
    check = judge(hsps, A, B, STRICT)
    assert check.ok and check.inverted
    assert check.back_left == 901 and check.back_right == 901


def test_indel_next_to_the_locus_is_tolerated() -> None:
    # a 3 kb insertion in B's left flank: the flanks still join on the right sides, only further away
    far = Region("b", "chr1", 1, 5100, "C" * 5100, 4000, 100)
    hsps = [Hsp(95, 900, 50, 950, 60, 960), Hsp(95, 900, 1150, 2050, 4140, 5040)]
    assert colinear(hsps, A, far, 200, inverted=False) == (901, 901)


def test_piece_on_the_wrong_side_is_rejected() -> None:
    # A's left flank aligns forward, but into B's right flank
    hsps = [Hsp(95, 900, 50, 950, 1150, 2050)]
    assert colinear(hsps, A, B, 200, inverted=False) == (0, 0)


def test_reciprocal_needs_both_loci() -> None:
    # the candidate has a 3 kb left flank of which only 900 bp align: enough for A, not for the candidate
    long = Region("c", "chr1", 1, 4100, "C" * 4100, 3000, 100)
    hsps = [Hsp(95, 900, 50, 950, 2050, 2950), Hsp(95, 900, 1150, 2050, 3150, 4050)]
    one_way = Validation(min_frac=0.5, min_bp=500, max_offset=200, reciprocal=False)
    assert judge(hsps, A, long, one_way).ok
    check = judge(hsps, A, long, STRICT)
    assert not check.ok and check.back_left == 901


def test_parse_hsps_filters() -> None:
    text = "90.0\t300\t1\t300\t5\t304\n70.0\t500\t1\t500\t1\t500\n95.0\t100\t1\t100\t1\t100\n"
    assert parse_hsps(text, 80, 200) == [Hsp(90.0, 300, 1, 300, 5, 304)]
    assert Hsp(90, 300, 1, 300, 400, 101).swapped() == Hsp(90, 300, 101, 400, 300, 1)


def test_as_args() -> None:
    assert as_args({"x": "asm10", "eqx": True, "secondary": False, "N": None}) == ["-x", "asm10", "--eqx"]
    assert as_args({"evalue": 1e-10, "perc_identity": 80}, single_dash=True) == [
        "-evalue",
        "1e-10",
        "-perc_identity",
        "80",
    ]
