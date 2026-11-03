from pathlib import Path

from spot import diagnose
from spot.config import HEADER_REGEX


def test_detects_species_and_coordinates(of_root: Path, tmp_path: Path) -> None:
    run = of_root / "runA" / "results" / "Results_Oct27"
    assert diagnose.column_prefixes([run]) == ["BB", "AA"]
    names = diagnose.headers([run], "AA")
    assert len(names) == 6 and names[0].startswith("chr1:1-80")
    assert diagnose.detect_header_regex(names) == HEADER_REGEX
    genome = tmp_path / "g.fa"
    genome.write_text(">chr2\nACGT\n>chr3\nACGT\n")
    report = diagnose.species_report(names, HEADER_REGEX, genome)
    assert report.on_genome == 0 and report.missing == {"chr1": 6}
    assert "add its lift" in report.text()
    lift = tmp_path / "split.tsv"
    lift.write_text("chr2\tchr1\t1\t10000\n")
    lifted = diagnose.species_report(names, HEADER_REGEX, genome, [lift])
    assert lifted.on_genome == 6


def test_other_header_formats() -> None:
    assert diagnose.detect_header_regex(["Chr1_100_200 x", "Chr2_5_9"]) != HEADER_REGEX


def test_initials() -> None:
    assert diagnose.initials("M. stenopetala") == "MS"
    assert diagnose.initials("") == ""
