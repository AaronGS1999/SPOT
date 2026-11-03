from pathlib import Path

import pytest

from spot import tools
from spot.cli import main


def test_inspect(config_file: Path, capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["inspect", "-c", str(config_file)]) == 0
    out = capsys.readouterr().out.splitlines()
    assert out[0].split() == [
        "run",
        "orthogroups",
        "shared",
        "AA_only",
        "BB_only",
        "AA_unassigned",
        "BB_unassigned",
        "trees",
        "shared_with_tree",
    ]
    assert out[1].split() == ["runA", "4", "2", "1", "1", "1", "1", "2", "1"]


def test_config_overrides(config_file: Path, capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["config", "-c", str(config_file), "-p", "relaxed", "-x", "asm20", "--set", "threads=2"]) == 0
    out = capsys.readouterr().out
    assert "profile: relaxed" in out
    assert "x: asm20" in out
    assert "threads: 2" in out


def test_bad_config_exits_2(config_file: Path) -> None:
    assert main(["config", "-c", str(config_file), "--anchors", "SYN,NOPE"]) == 2


def test_profiles(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["profiles"]) == 0
    assert [line.split()[0] for line in capsys.readouterr().out.splitlines()] == ["default", "relaxed"]


def test_gepard_jar_from_wrapper(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    share = tmp_path / "share" / "gepard"
    (share / "dist").mkdir(parents=True)
    (share / "resources" / "matrices").mkdir(parents=True)
    jar = share / "dist" / "Gepard-2.1.jar"
    jar.touch()
    (share / "resources" / "matrices" / "edna.mat").touch()
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    wrapper = bin_dir / "gepardcmd"
    wrapper.write_text(f"#!/bin/bash\n\njava -cp {jar} org.gepard.client.cmdline.CommandLine\n")
    wrapper.chmod(0o755)
    monkeypatch.setenv("PATH", str(bin_dir))
    monkeypatch.delenv("SPOT_GEPARD_JAR", raising=False)
    gepard = tools.find_gepard()
    assert gepard is not None
    assert gepard.jar == jar
    assert gepard.version == "2.1"
    assert gepard.matrix("edna.mat").is_file()


def test_run_reports_missing_inputs(config_file: Path, capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["run", "-c", str(config_file), "--set", "species.AA.centromeres=cen.bed"]) == 2
    err = capsys.readouterr().err
    assert "species AA: genome: not found" in err
    assert "species AA: centromeres: not found" in err
    assert "species BB: genome: not found: /abs/b.fasta" in err
