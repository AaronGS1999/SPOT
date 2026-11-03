import logging
from pathlib import Path

import pytest

from spot import config, orthofinder
from spot.config import OrthoFinderInput, Species
from spot.errors import InputError

SPECIES = {
    "AA": Species(label="A", match="^AA_"),
    "BB": Species(label="B", match="^BB_"),
}


def test_discover(of_root: Path) -> None:
    runs = orthofinder.discover(of_root)
    assert list(runs) == ["runA", "runB"]
    assert runs["runA"] == of_root / "runA" / "results" / "Results_Oct27"
    assert orthofinder.discover(runs["runA"]) == {"Results_Oct27": runs["runA"]}


def test_select_keeps_config_order(of_root: Path) -> None:
    source = OrthoFinderInput(root=of_root, runs={"runB": {}, "runA": {}})
    assert list(orthofinder.select(source)) == ["runB", "runA"]
    with pytest.raises(InputError, match="'runC' not found"):
        orthofinder.select(OrthoFinderInput(root=of_root, runs=["runC"]))


def test_load(of_root: Path) -> None:
    run = orthofinder.load("runA", of_root / "runA" / "results" / "Results_Oct27", SPECIES)
    assert run.columns == {"AA": "AA_set", "BB": "BB_set"}
    assert run.version == "2.5.5"
    assert run.dna is True
    assert [og.id for og in run.shared()] == ["OG0000000", "OG0000002"]
    assert [og.id for og in run.only("AA")] == ["OG0000001"]
    assert [og.id for og in run.only("BB")] == ["OG0000003"]
    assert [g.id for g in run.unassigned("AA")] == ["AA_6"]
    assert run.genes["AA_2"].orthogroup == "OG0000000"
    assert run.genes["AA_2"].description == "chr1:101-180 K=0.1"
    assert run.genes["BB_5"].assigned is False
    assert run.trees == {"OG0000000", "OG0000003"}
    assert run.tree("OG0000002") is None
    assert run.tree("OG0000000").name == "OG0000000_tree.txt"
    assert run.summary() == {
        "orthogroups": 4,
        "shared": 2,
        "AA_only": 1,
        "BB_only": 1,
        "AA_unassigned": 1,
        "BB_unassigned": 1,
        "trees": 2,
        "shared_with_tree": 1,
    }


def test_column_matching_two_species(of_root: Path) -> None:
    species = {"AA": Species(label="A", match="_set"), "BB": Species(label="B", match="^BB_")}
    with pytest.raises(InputError, match="matches species AA, BB"):
        orthofinder.load("runA", of_root / "runA" / "results" / "Results_Oct27", species)


def test_column_without_species(of_root: Path) -> None:
    species = {"AA": Species(label="A", match="^AA_"), "BB": Species(label="B", match="^CC_")}
    with pytest.raises(InputError, match="'BB_set' matches no species"):
        orthofinder.load("runA", of_root / "runA" / "results" / "Results_Oct27", species)


def test_gene_in_two_orthogroups(of_root: Path) -> None:
    results = of_root / "runA" / "results" / "Results_Oct27"
    table = results / "Orthogroups" / "Orthogroups.tsv"
    table.write_text(table.read_text().replace("BB_2\tAA_5", "BB_1\tAA_5"))
    with pytest.raises(InputError, match="BB_1 is in OG0000000 and OG0000002"):
        orthofinder.load("runA", results, SPECIES)


def test_ids_missing_from_sequence_ids(of_root: Path) -> None:
    results = of_root / "runA" / "results" / "Results_Oct27"
    ids = results / "WorkingDirectory" / "SequenceIDs.txt"
    ids.write_text(ids.read_text().replace(": AA_3 ", ": AA_3x "))
    with pytest.raises(InputError, match="AA_3"):
        orthofinder.load("runA", results, SPECIES)


@pytest.mark.parametrize(
    ("command", "search", "warned"),
    [
        ("orthofinder -f in -d -S diamond -X", "diamond", True),
        ("orthofinder -f in -S diamond -d -X", "diamond", True),
        ("orthofinder -f in -d -X", "blast_nucl", False),
        ("orthofinder -f in -d -S blast_nucl", "blast_nucl", False),
        ("orthofinder -f in -S mmseqs", "mmseqs", False),
    ],
)
def test_search_program(
    of_root: Path, caplog: pytest.LogCaptureFixture, command: str, search: str, warned: bool
) -> None:
    results = of_root / "runA" / "results" / "Results_Oct27"
    (results / "Log.txt").write_text(f"Started OrthoFinder version 2.5.5\nCommand Line: {command}\n")
    logging.getLogger("spot").addHandler(caplog.handler)
    try:
        run = orthofinder.load("runA", results, SPECIES)
    finally:
        logging.getLogger("spot").removeHandler(caplog.handler)
    assert run.search == search
    assert ("protein aligner" in caplog.text) is warned


def test_load_all(config_file: Path) -> None:
    runs = orthofinder.load_all(config.load(config_file))
    assert [r.name for r in runs] == ["runA", "runB"]
