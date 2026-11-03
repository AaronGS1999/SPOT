from pathlib import Path

import pytest

from spot.errors import InputError
from spot.lift import Lift, LiftStep

AGP = """## agp-version 2.1
# AGP created by RagTag v2.1.1
chr1_RagTag\t1\t100\t1\tW\tctgA\t1\t100\t+
chr1_RagTag\t101\t200\t2\tU\t100\tscaffold\tyes\talign_genus
chr1_RagTag\t201\t250\t3\tW\tctgB\t1\t50\t-
ctgC\t1\t80\t1\tW\tctgC\t1\t80\t+
"""


@pytest.fixture
def steps(tmp_path: Path) -> tuple[Path, Path]:
    split = tmp_path / "split.tsv"
    split.write_text("ctgA\told1\t1\t100\nctgB\told1\t1101\t1150\n")
    agp = tmp_path / "ragtag.scaffold.agp"
    agp.write_text(AGP)
    return split, agp


def test_agp_forward_and_reverse(steps: tuple[Path, Path]) -> None:
    step = LiftStep.read(steps[1])
    assert step.map("ctgA", 10, 20) == ("chr1_RagTag", 10, 20)
    assert step.map("ctgB", 1, 10) == ("chr1_RagTag", 241, 250)
    assert step.map("ctgB", 41, 50) == ("chr1_RagTag", 201, 210)
    assert step.map("ctgC", 5, 6) == ("ctgC", 5, 6)
    assert step.map("other", 5, 6) == ("other", 5, 6)
    assert step.map("ctgA", 95, 105) is None


def test_chain(steps: tuple[Path, Path]) -> None:
    lift = Lift.read(list(steps), strip_suffix="_RagTag")
    assert lift("old1", 10, 20) == ("chr1", 10, 20)
    assert lift("old1", 1101, 1110) == ("chr1", 241, 250)
    assert lift("old1", 500, 510) is None
    assert lift("old2", 1, 5) == ("old2", 1, 5)


def test_table_header_and_errors(tmp_path: Path) -> None:
    table = tmp_path / "split.tsv"
    table.write_text("name\tparent\tstart\tend\nA\tX\t11\t20\n")
    assert LiftStep.read(table).map("X", 11, 12) == ("A", 1, 2)
    table.write_text("A\tX\t11\n")
    with pytest.raises(InputError, match="expected name, parent, start, end"):
        LiftStep.read(table)
