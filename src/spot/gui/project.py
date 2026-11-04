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
        self.label = QLineEdit()
        self.match = QLineEdit()
        self.genome = PathEdit("file", "Genome FASTA", FASTA)
        self.centromeres = PathEdit("file", "Centromeres BED", "BED (*.bed);;All files (*)")
        self.use_table = QCheckBox("from a table")
        self.coords_table = PathEdit("file", "Coordinate table", "TSV (*.tsv *.txt);;All files (*)")
        self.custom_regex = QCheckBox("custom pattern")
        self.header_regex = QLineEdit(HEADER_REGEX)
        self.use_lift = QCheckBox("they refer to another assembly")
        self.lift = PathList("Lift files (AGP or split tables, in order)")
        self.strip_suffix = QLineEdit()
        self.min_identity = QDoubleSpinBox()
        self.min_identity.setRange(0.5, 1.0)
        self.min_identity.setSingleStep(0.005)
        self.min_identity.setDecimals(3)
        self.min_identity.setValue(0.99)
        self.status = QLabel()
        self.status.setWordWrap(True)
        self.status.setStyleSheet("color: #8a8a8a")
        form = QFormLayout(self)
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
        self.form = form

        def row(name: str, widget: QWidget | QHBoxLayout, help_text: str) -> QLabel:
            label = QLabel(name)
            label.setToolTip(help_text)
            if isinstance(widget, QWidget):
                widget.setToolTip(help_text)
            form.addRow(label, widget)
            return label

        names = QHBoxLayout()
        names.addWidget(self.key)
        key_help = "Short name of the species in tables, file names and the confirmed-pair columns, e.g. SP1."
        self.key.setToolTip(key_help)
        label_title = QLabel("name")
        label_help = "Species name used in the figures, e.g. A. thaliana (shown in italics)."
        label_title.setToolTip(label_help)
        self.label.setToolTip(label_help)
        names.addWidget(label_title)
        names.addWidget(self.label, 1)
        row("key", names, key_help)
        row(
            "column prefix",
            self.match,
            "Start of this species' column in Orthogroups.tsv, i.e. of its input FASTA name: SP1 for "
            "SP1_genes.fasta. Filled in by Find runs. A regular expression is also accepted.",
        )
        row(
            "genome",
            self.genome,
            "Genome FASTA of this species: the assembly the whole-genome alignment (SyRI) was made on. "
            "Every sequence is checked against it at its coordinates.",
        )
        coords = QHBoxLayout()
        coords.addWidget(QLabel("from the sequence names"), 1)
        coords.addWidget(self.custom_regex)
        coords.addWidget(self.use_table)
        row(
            "coordinates",
            coords,
            "Where each sequence lies in the genome. By default SPOT reads chr:start-end from the FASTA "
            "header of each sequence (the format is detected). Use a table only if the headers have no "
            "coordinates.",
        )
        self.regex_label = row(
            "pattern",
            self.header_regex,
            "Regular expression with the named groups chr, start and end, applied to the FASTA header.",
        )
        self.table_label = row(
            "table",
            self.coords_table,
            "Tab-separated table with the columns id, chr, start and end (1-based, inclusive), one row "
            "per sequence. Only needed when the FASTA headers carry no coordinates.",
        )
        lift_row = QHBoxLayout()
        lift_row.addWidget(self.use_lift)
        lift_row.addStretch()
        row(
            "lift",
            lift_row,
            "Tick when the coordinates in the headers were taken on another assembly than the genome above "
            "(for example the original scaffolds before RagTag, or before a chromosome was split).",
        )
        self.lift_label = row(
            "lift files",
            self.lift,
            "AGP files (e.g. ragtag.scaffold.agp) or split tables (piece, chromosome, start, end) that turn "
            "the old coordinates into coordinates of the genome, applied from top to bottom.",
        )
        self.strip_label = row(
            "strip suffix",
            self.strip_suffix,
            "Text removed from the end of the lifted chromosome names, e.g. _RagTag (chr1_RagTag -> chr1).",
        )
        row(
            "min. identity",
            self.min_identity,
            "Each sequence is compared with the genome at its coordinates (both strands). It is accepted "
            "when it is identical or at least this similar (0.99 = 99 %); otherwise it is reported as "
            "sequence_mismatch and left out, because its coordinates cannot be trusted.",
        )
        row(
            "centromeres",
            self.centromeres,
            "Optional BED file of centromeres, only drawn in the synteny figure.",
        )
        form.addRow(self.status)
        for checkbox, rows in (
            (self.custom_regex, (self.regex_label,)),
            (self.use_table, (self.table_label,)),
            (self.use_lift, (self.lift_label, self.strip_label)),
        ):
            checkbox.toggled.connect(lambda on, rows=rows: self._show_rows(rows, on))
            self._show_rows(rows, False)
        self.use_lift.toggled.connect(self._timer_start)
        self.use_table.toggled.connect(self.changed)
        self.use_lift.toggled.connect(self.changed)
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

    def _show_rows(self, labels: tuple[QLabel, ...], on: bool) -> None:
        for label in labels:
            self.form.setRowVisible(label, on)

    def _timer_start(self) -> None:
        if hasattr(self, "_timer"):
            self._timer.start()

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
            self.lift.paths() if self.use_lift.isChecked() else [],
            (self.strip_suffix.text().strip() or None) if self.use_lift.isChecked() else None,
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
        self.use_table.setChecked(coords not in (None, "header"))
        self.use_lift.setChecked(bool(data.get("lift") or data.get("strip_suffix")))
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
        if self.use_table.isChecked() and self.coords_table.path():
            data["coords"] = str(self.coords_table.path())
        regex = self.header_regex.text().strip()
        if self.custom_regex.isChecked() and regex and regex != HEADER_REGEX:
            data["header_regex"] = regex
        elif not self.custom_regex.isChecked() and self.runs and self.match.text().strip():
            detected = diagnose.detect_header_regex(diagnose.headers(self.runs, self.match.text().strip()))
            if detected and detected != HEADER_REGEX:
                data["header_regex"] = detected
        if self.use_lift.isChecked() and self.lift.paths():
            data["lift"] = [str(p) for p in self.lift.paths()]
        if self.use_lift.isChecked() and self.strip_suffix.text().strip():
            data["strip_suffix"] = self.strip_suffix.text().strip()
        if abs(self.min_identity.value() - 0.99) > 1e-9:
            data["min_identity"] = round(self.min_identity.value(), 4)
        if self.centromeres.path():
            data["centromeres"] = str(self.centromeres.path())
        return self.key.text().strip(), data


class RunsTable(QTableWidget):
    """OrthoFinder runs found under the root folder: include, name, label, colour, unit."""

    changed = Signal()
    COLUMNS = ("use", "run", "name in figures", "colour", "sequence type", "results folder")
    HELP = (
        "Include this OrthoFinder run in the analysis.",
        "Folder name of the run.",
        "Title of this run in the figures (Venn diagram, legends), e.g. Genes. Empty: the run name.",
        "Colour of this run in all figures (syntenic orthologs in the Venn diagram, links and dots in the "
        "synteny figure and the dotplot). auto: a default palette.",
        "What the sequences are called in the Venn diagram counts, e.g. genes or TEs. Empty: sequences.",
        "OrthoFinder results folder used (the latest Results_* of the run).",
    )

    def __init__(self) -> None:
        super().__init__(0, len(self.COLUMNS))
        self.setHorizontalHeaderLabels(self.COLUMNS)
        for i, text in enumerate(self.HELP):
            header_item = self.horizontalHeaderItem(i)
            if header_item is not None:
                header_item.setToolTip(text)
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
        self.root.setToolTip(
            "Folder that holds the OrthoFinder runs, one per subfolder (each with its Results_* folder). "
            "Find runs lists them; the latest Results_* of each run is used."
        )
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
        self.threads.setToolTip(
            "Programs run in parallel (BLASTN searches, Gepard dotplots, minimap2, SyRI). auto uses the "
            "performance cores of this computer."
        )
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
        for widget, text in (
            (self.reference, "Species whose genome is the reference of the whole-genome alignment."),
            (self.query, "Species aligned to the reference."),
            (self.outdir, "Folder where SPOT writes its results; a run can be resumed in the same folder."),
            (
                self.use_syri,
                "Use the syri.out of an alignment already made between these two genomes (reference first) "
                "instead of aligning them again with minimap2 and SyRI.",
            ),
            (self.syri_out, "syri.out made with this reference and query genome."),
            (
                self.resume,
                "Keep the stages a previous run in this folder already finished with the same inputs and "
                "settings (locating sequences, BLASTN validation, the alignment).",
            ),
        ):
            widget.setToolTip(text)
            label = form.labelForField(widget)
            if label is not None:
                label.setToolTip(text)
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
        # a species whose name is already typed gets the prefix made of its initials (A. thaliana -> AT)
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
