import os
from pathlib import Path

import pytest
import yaml

pytest.importorskip("PySide6")
pytest.importorskip("pytestqt")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from spot import config  # noqa: E402
from spot.gui.main import MainWindow, _absolute, _relative  # noqa: E402
from spot.progress import parse  # noqa: E402


def test_paths_round_trip(tmp_path: Path) -> None:
    data = {
        "orthofinder": {"root": "of"},
        "species": {"A": {"genome": "g/a.fa", "lift": ["l/x.agp"], "coords": "header"}},
        "outdir": "out",
    }
    absolute = _absolute(data, tmp_path)
    assert absolute["species"]["A"]["genome"] == str(tmp_path / "g" / "a.fa")
    assert absolute["species"]["A"]["coords"] == "header"
    assert _relative(absolute, tmp_path) == data


def test_progress_lines() -> None:
    assert parse("09:00:00 INFO    @progress validate 3 10") == ("validate", 3, 10)
    assert parse("nothing here") is None


def test_window_loads_edits_and_saves(qtbot, config_file: Path, monkeypatch: pytest.MonkeyPatch) -> None:  # type: ignore[no-untyped-def]
    window = MainWindow()
    qtbot.addWidget(window)
    window.open_config(config_file)
    assert window.project.reference.currentText() == "AA"
    assert window.project.runs.rowCount() == 2
    row = window.parameters.forms["validation"].rows["window"]
    row.editor.set_value(5000)
    row._edited()
    window.save()
    saved = yaml.safe_load(config_file.read_text())
    assert saved["validation"] == {"window": 5000}
    assert saved["species"]["AA"]["genome"] == "genomes/a.fasta"
    assert saved["threads"] == "auto"
    loaded = config.load(config_file)
    assert loaded.validation.window == 5000 and loaded.reference == "AA"
    row._reset()
    assert "validation" not in window.data()
    window.dirty = False
