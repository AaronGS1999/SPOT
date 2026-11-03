"""Project page: OrthoFinder input, the two species and where the results go."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from PySide6.QtCore import QObject, QRunnable, Qt, QThreadPool, QTimer, Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QPushButton,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from .. import diagnose
from ..config import HEADER_REGEX
from ..errors import SpotError
from ..orthofinder import discover
from ..resources import performance_cores
from .widgets import ColorButton, PathEdit, PathList

FASTA = "FASTA (*.fa *.fasta *.fna *.fas);;All files (*)"


class _Signals(QObject):
    done = Signal(int, str)


class _Diagnose(QRunnable):
    """Checks a species' coordinates in a worker thread (reading a genome index can take seconds)."""

    def __init__(
        self,
        generation: int,
        runs: list[Path],
        match: str,
        regex: str | None,
        genome: Path | None,
        lift: list[Path],
        strip: str | None,
    ) -> None:
        super().__init__()
        self.signals = _Signals()
        self.args = (generation, runs, match, regex, genome, lift, strip)

    def run(self) -> None:
        generation, runs, match, regex, genome, lift, strip = self.args
        try:
            descriptions = diagnose.headers(runs, match)
            pattern = regex or diagnose.detect_header_regex(descriptions) or HEADER_REGEX
            text = diagnose.species_report(descriptions, pattern, genome, lift, strip).text()
        except Exception as e:  # noqa: BLE001 - shown to the user, never fatal
            text = f"cannot check: {e}"
        self.signals.done.emit(generation, text)


class SpeciesEditor(QGroupBox):
    changed = Signal()

    def __init__(self, title: str) -> None:
        super().__init__(title)
        self.runs: list[Path] = []
        self.key = QLineEdit()
        self.key.setMaximumWidth(100)
        self.key.setToolTip("Short name used in tables and file names (e.g. MS)")
        self.label = QLineEdit()
        self.label.setPlaceholderText("species name for figures, e.g. M. stenopetala")
        self.match = QLineEdit()
        self.match.setPlaceholderText("MS")
        self.match.setToolTip(
            "Start of this species' column names in Orthogroups.tsv (the input FASTA names), e.g. MS for "
            "MS_NUMT_evento2. A regular expression is also accepted."
        )
        self.genome = PathEdit("file", "Genome FASTA", FASTA)
        self.genome.setPlaceholderText("the assembly used in the whole-genome alignment")
        self.centromeres = PathEdit("file", "Centromeres BED", "BED (*.bed);;All files (*)")
        self.centromeres.setPlaceholderText("optional, only drawn in the synteny figure")
        self.coords_table = PathEdit("file", "Coordinate table", "TSV (*.tsv *.txt);;All files (*)")
        self.coords_table.setPlaceholderText(
            "optional: id chr start end, when the names carry no coordinates"
        )
        self.custom_regex = QCheckBox("custom")
        self.header_regex = QLineEdit(HEADER_REGEX)
        self.header_regex.setVisible(False)
        self.custom_regex.toggled.connect(self.header_regex.setVisible)
        self.detected = ""
        self.lift = PathList("Lift files (AGP or split tables, in order)")
        self.lift.setToolTip(
            "Only when the coordinates in the sequence names refer to another assembly than the genome, "
            "e.g. a scaffolding (RagTag AGP) or a split chromosome (split table), applied in order"
        )
        self.strip_suffix = QLineEdit()
        self.strip_suffix.setPlaceholderText("removed from lifted names, e.g. _RagTag")
        self.min_identity = QDoubleSpinBox()
        self.min_identity.setRange(0.5, 1.0)
        self.min_identity.setSingleStep(0.005)
        self.min_identity.setDecimals(3)
        self.min_identity.setValue(0.99)
        self.min_identity.setToolTip("Minimum identity between a sequence and the genome at its coordinates")
        self.status = QLabel()
        self.status.setWordWrap(True)
        self.status.setStyleSheet("color: #8a8a8a")
        form = QFormLayout(self)
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
        top = QHBoxLayout()
        top.addWidget(QLabel("key"))
        top.addWidget(self.key)
        top.addWidget(QLabel("label"))
        top.addWidget(self.label, 1)
        form.addRow(top)
        form.addRow("column prefix", self.match)
        form.addRow("genome", self.genome)
        regex_row = QHBoxLayout()
        self.regex_note = QLabel("detected from the names")
        self.regex_note.setStyleSheet("color: #8a8a8a")
        regex_row.addWidget(self.regex_note, 1)
        regex_row.addWidget(self.custom_regex)
        form.addRow("coordinates", regex_row)
        form.addRow("", self.header_regex)
        form.addRow("coordinate table", self.coords_table)
        form.addRow("lift", self.lift)
        form.addRow("strip suffix", self.strip_suffix)
        form.addRow("min. identity", self.min_identity)
        form.addRow("centromeres", self.centromeres)
        form.addRow(self.status)
        for w in (self.key, self.label, self.match, self.header_regex, self.strip_suffix):
            w.textChanged.connect(self.changed)
        for p in (self.genome, self.centromeres, self.coords_table):
            p.changed.connect(self.changed)
        self.lift.changed.connect(self.changed)
        self.min_identity.valueChanged.connect(self.changed)
        self.custom_regex.toggled.connect(self.changed)
        self._generation = 0
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(500)
        self._timer.timeout.connect(self._diagnose)
        for w in (self.match, self.header_regex, self.strip_suffix):
            w.textChanged.connect(self._timer.start)
        self.genome.changed.connect(self._timer.start)
        self.lift.changed.connect(self._timer.start)
        self.custom_regex.toggled.connect(self._timer.start)

    def set_runs(self, runs: list[Path]) -> None:
        self.runs = runs
        self._timer.start()

    def _diagnose(self) -> None:
        if not self.runs or not self.match.text().strip():
            self.status.setText("")
            return
        self._generation += 1
        regex = self.header_regex.text().strip() if self.custom_regex.isChecked() else None
        job = _Diagnose(
            self._generation,
            list(self.runs),
            self.match.text().strip(),
            regex,
            self.genome.path(),
            self.lift.paths(),
            self.strip_suffix.text().strip() or None,
        )
        job.signals.done.connect(self._diagnosed)
        self.status.setText("checking the coordinates…")
        QThreadPool.globalInstance().start(job)

    def _diagnosed(self, generation: int, text: str) -> None:
        if generation == self._generation:
            self.status.setText(text)

    def load(self, key: str, data: dict[str, Any]) -> None:
        self.key.setText(key)
        self.label.setText(str(data.get("label", "")))
        self.match.setText(str(data.get("match", "")))
        self.genome.set_path(data.get("genome"))
        self.centromeres.set_path(data.get("centromeres"))
        coords = data.get("coords", "header")
        self.coords_table.set_path(None if coords in (None, "header") else coords)
        regex = data.get("header_regex")
        self.custom_regex.setChecked(bool(regex))
        self.header_regex.setText(str(regex or HEADER_REGEX))
        self.lift.set_paths(data.get("lift", []) or [])
        self.strip_suffix.setText(str(data.get("strip_suffix") or ""))
        self.min_identity.setValue(float(data.get("min_identity", 0.99)))

    def store(self) -> tuple[str, dict[str, Any]]:
        data: dict[str, Any] = {"label": self.label.text().strip(), "match": self.match.text().strip()}
        if self.genome.path():
            data["genome"] = str(self.genome.path())
        if self.coords_table.path():
            data["coords"] = str(self.coords_table.path())
        regex = self.header_regex.text().strip()
        if self.custom_regex.isChecked() and regex and regex != HEADER_REGEX:
            data["header_regex"] = regex
        elif not self.custom_regex.isChecked() and self.runs and self.match.text().strip():
            detected = diagnose.detect_header_regex(diagnose.headers(self.runs, self.match.text().strip()))
            if detected and detected != HEADER_REGEX:
                data["header_regex"] = detected
        if self.lift.paths():
            data["lift"] = [str(p) for p in self.lift.paths()]
        if self.strip_suffix.text().strip():
            data["strip_suffix"] = self.strip_suffix.text().strip()
        if abs(self.min_identity.value() - 0.99) > 1e-9:
            data["min_identity"] = round(self.min_identity.value(), 4)
        if self.centromeres.path():
            data["centromeres"] = str(self.centromeres.path())
        return self.key.text().strip(), data


class RunsTable(QTableWidget):
    """OrthoFinder runs found under the root folder: include, name, label, colour, unit."""

    changed = Signal()
    COLUMNS = ("use", "run", "label", "colour", "unit", "results folder")

    def __init__(self) -> None:
        super().__init__(0, len(self.COLUMNS))
        self.setHorizontalHeaderLabels(self.COLUMNS)
        self.verticalHeader().setVisible(False)
        header = self.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(5, QHeaderView.ResizeMode.Stretch)
        self.setMinimumHeight(120)
        self.itemChanged.connect(lambda _: self.changed.emit())

    def fill(self, found: dict[str, Path], settings: dict[str, dict[str, Any]], select_all: bool) -> None:
        self.blockSignals(True)
        self.setRowCount(0)
        names = list(settings) + [n for n in found if n not in settings]  # keep the configured order
        for name in names:
            row = self.rowCount()
            self.insertRow(row)
            use = QCheckBox()
            use.setChecked(select_all or name in settings)
            use.toggled.connect(self.changed)
            self.setCellWidget(row, 0, use)
            item = QTableWidgetItem(name)
            item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.setItem(row, 1, item)
            s = settings.get(name, {}) or {}
            self.setItem(row, 2, QTableWidgetItem(str(s.get("label") or "")))
            color = ColorButton(s.get("color"))
            color.changed.connect(self.changed)
            self.setCellWidget(row, 3, color)
            self.setItem(row, 4, QTableWidgetItem(str(s.get("unit") or "")))
            where = found.get(name)
            path = QTableWidgetItem(str(where) if where else "not found under the root folder")
            path.setFlags(path.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.setItem(row, 5, path)
        self.blockSignals(False)

    def _text(self, row: int, column: int) -> str:
        item = self.item(row, column)
        return item.text().strip() if item is not None else ""

    def settings(self) -> dict[str, dict[str, Any]]:
        runs: dict[str, dict[str, Any]] = {}
        for row in range(self.rowCount()):
            use = self.cellWidget(row, 0)
            if not isinstance(use, QCheckBox) or not use.isChecked():
                continue
            name = self._text(row, 1)
            entry: dict[str, Any] = {}
            label, unit = self._text(row, 2), self._text(row, 4)
            color = self.cellWidget(row, 3)
            if label:
                entry["label"] = label
            if isinstance(color, ColorButton) and color.color():
                entry["color"] = color.color()
            if unit:
                entry["unit"] = unit
            runs[name] = entry
        return runs


class ProjectPage(QWidget):
    changed = Signal()

    def __init__(self) -> None:
        super().__init__()
        # OrthoFinder
        self.root = PathEdit("dir", "OrthoFinder results folder")
        self.root.setPlaceholderText("folder with one OrthoFinder run per subfolder")
        detect = QPushButton("Find runs")
        detect.clicked.connect(lambda: self.detect(select_all=False))
        self.runs = RunsTable()
        self.runs_note = QLabel()
        self.runs_note.setStyleSheet("color: #888")
        of = QGroupBox("OrthoFinder")
        of_layout = QGridLayout(of)
        of_layout.addWidget(QLabel("root"), 0, 0)
        of_layout.addWidget(self.root, 0, 1)
        of_layout.addWidget(detect, 0, 2)
        of_layout.addWidget(self.runs, 1, 0, 1, 3)
        of_layout.addWidget(self.runs_note, 2, 0, 1, 3)
        # species
        self.species = [SpeciesEditor("Species 1"), SpeciesEditor("Species 2")]
        sp_row = QHBoxLayout()
        for s in self.species:
            sp_row.addWidget(s)
        # analysis
        self.reference = QComboBox()
        self.query = QComboBox()
        self.outdir = PathEdit("dir", "Output folder")
        self.threads = QSpinBox()
        self.threads.setRange(0, 512)
        self.threads.setSpecialValueText(f"auto ({performance_cores()})")
        self.threads.setToolTip("Worker threads; auto = performance cores of this computer")
        self.use_syri = QCheckBox("use an existing syri.out")
        self.syri_out = PathEdit("file", "syri.out", "SyRI output (*.out);;All files (*)")
        self.syri_out.setEnabled(False)
        self.use_syri.toggled.connect(self.syri_out.setEnabled)
        self.resume = QCheckBox("reuse finished stages of a previous run in the same folder")
        self.resume.setChecked(True)
        run_box = QGroupBox("Analysis")
        form = QFormLayout(run_box)
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
        roles = QHBoxLayout()
        roles.addWidget(QLabel("reference"))
        roles.addWidget(self.reference)
        roles.addWidget(QLabel("query"))
        roles.addWidget(self.query)
        roles.addStretch()
        form.addRow(roles)
        form.addRow("output folder", self.outdir)
        form.addRow("threads", self.threads)
        form.addRow(self.use_syri, self.syri_out)
        form.addRow(self.resume)
        layout = QVBoxLayout(self)
        layout.addWidget(of)
        layout.addLayout(sp_row)
        layout.addWidget(run_box)
        layout.addStretch()
        for s in self.species:
            s.changed.connect(self._keys)
            s.changed.connect(self.changed)
        self.root.changed.connect(self.changed)
        self.runs.changed.connect(self.changed)
        for w in (self.reference, self.query):
            w.currentTextChanged.connect(self.changed)
        self.outdir.changed.connect(self.changed)
        self.threads.valueChanged.connect(self.changed)
        self.use_syri.toggled.connect(self.changed)
        self.syri_out.changed.connect(self.changed)
        self.resume.toggled.connect(self.changed)

    def detect(self, select_all: bool) -> None:
        root = self.root.path()
        current = self.runs.settings()
        found: dict[str, Path] = {}
        if root and root.is_dir():
            try:
                found = discover(root)
                self.runs_note.setText(f"{len(found)} OrthoFinder runs found")
            except SpotError as e:
                self.runs_note.setText(str(e))
        else:
            self.runs_note.setText("choose the folder that holds the OrthoFinder runs")
        self.runs.fill(found, current, select_all or not current)
        self._species_from_runs(list(found.values()))
        self.changed.emit()

    def _species_from_runs(self, paths: list[Path]) -> None:
        """Hand the runs to the species editors and fill empty species from the column prefixes."""
        for editor in self.species:
            editor.set_runs(paths)
        prefixes = diagnose.column_prefixes(paths)
        empty = [e for e in self.species if not e.match.text().strip()]
        free = [p for p in prefixes if p not in {e.match.text().strip() for e in self.species}]
        # a species whose name is already typed gets the prefix made of its initials (M. stenopetala -> MS)
        for editor in list(empty):
            guess = diagnose.initials(editor.label.text())
            if guess in free:
                self._assign(editor, guess)
                empty.remove(editor)
                free.remove(guess)
        for editor, prefix in zip(empty, free, strict=False):
            self._assign(editor, prefix)
        self._keys()

    @staticmethod
    def _assign(editor: SpeciesEditor, prefix: str) -> None:
        if editor.key.text().strip() in ("", "SP1", "SP2"):
            editor.key.setText(prefix)
        editor.match.setText(prefix)

    def _keys(self) -> None:
        keys = [k for k in (s.key.text().strip() for s in self.species) if k]
        chosen = [self.reference.currentText(), self.query.currentText()]
        for i, box in enumerate((self.reference, self.query)):
            box.blockSignals(True)
            box.clear()
            box.addItems(keys)
            if chosen[i] in keys:
                box.setCurrentText(chosen[i])
            elif len(keys) > i:
                box.setCurrentIndex(i)
            box.blockSignals(False)
        if len(keys) == 2 and self.reference.currentText() == self.query.currentText():
            self.query.setCurrentText(next(k for k in keys if k != self.reference.currentText()))

    def load(self, data: dict[str, Any]) -> None:
        of = data.get("orthofinder", {}) or {}
        self.root.set_path(of.get("root"))
        runs = of.get("runs", {}) or {}
        if isinstance(runs, list):
            runs = {name: {} for name in runs}
        found: dict[str, Path] = {}
        root = self.root.path()
        if root and root.is_dir():
            try:
                found = discover(root)
            except SpotError:
                found = {}
        self.runs.fill(found, runs, not runs)
        species = list((data.get("species") or {}).items())
        for editor, (key, sp) in zip(self.species, species + [("", {})] * (2 - len(species)), strict=False):
            editor.load(key, sp or {})
        self._species_from_runs(list(found.values()))
        if data.get("reference"):
            self.reference.setCurrentText(str(data["reference"]))
        if data.get("query"):
            self.query.setCurrentText(str(data["query"]))
        self.outdir.set_path(data.get("outdir"))
        threads = data.get("threads", "auto")
        self.threads.setValue(0 if threads in (None, "auto") else int(threads))
        self.use_syri.setChecked(bool(data.get("syri_out")))
        self.syri_out.set_path(data.get("syri_out"))
        self.resume.setChecked(bool(data.get("resume", True)))

    def store(self, data: dict[str, Any]) -> None:
        of: dict[str, Any] = {"root": str(self.root.path()) if self.root.path() else ""}
        runs = self.runs.settings()
        if runs:
            of["runs"] = runs
        data["orthofinder"] = of
        data["species"] = dict(s.store() for s in self.species)
        data["reference"] = self.reference.currentText()
        data["query"] = self.query.currentText()
        if self.outdir.path():
            data["outdir"] = str(self.outdir.path())
        data["threads"] = "auto" if self.threads.value() == 0 else self.threads.value()
        if self.use_syri.isChecked() and self.syri_out.path():
            data["syri_out"] = str(self.syri_out.path())
        else:
            data.pop("syri_out", None)
        if not self.resume.isChecked():
            data["resume"] = False
        else:
            data.pop("resume", None)
