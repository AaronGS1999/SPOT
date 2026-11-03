"""Results page: summary, filterable tables, figures and the dotplot of any confirmed pair."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

from PySide6.QtCore import (
    QAbstractTableModel,
    QModelIndex,
    QPersistentModelIndex,
    QSortFilterProxyModel,
    Qt,
    QUrl,
)
from PySide6.QtGui import QDesktopServices, QPixmap
from PySide6.QtWidgets import (
    QComboBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPlainTextEdit,
    QPushButton,
    QScrollArea,
    QSplitter,
    QTableView,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

Index = QModelIndex | QPersistentModelIndex


class TsvModel(QAbstractTableModel):
    def __init__(self) -> None:
        super().__init__()
        self.header: list[str] = []
        self.rows: list[list[str]] = []

    def load(self, path: Path) -> None:
        self.beginResetModel()
        self.header, self.rows = [], []
        if path.is_file():
            with open(path, newline="", encoding="utf-8") as fh:
                reader = csv.reader(fh, delimiter="\t")
                self.header = next(reader, [])
                self.rows = [r + [""] * (len(self.header) - len(r)) for r in reader]
        self.endResetModel()

    def rowCount(self, parent: Index = QModelIndex()) -> int:  # noqa: N802, B008
        return 0 if parent.isValid() else len(self.rows)

    def columnCount(self, parent: Index = QModelIndex()) -> int:  # noqa: N802, B008
        return 0 if parent.isValid() else len(self.header)

    def data(self, index: Index, role: int = Qt.ItemDataRole.DisplayRole) -> Any:
        if not index.isValid():
            return None
        value = self.rows[index.row()][index.column()]
        if role == Qt.ItemDataRole.DisplayRole:
            return value
        if role == Qt.ItemDataRole.UserRole:  # numeric sort
            try:
                return float(value)
            except ValueError:
                return value
        return None

    def headerData(  # noqa: N802
        self, section: int, orientation: Qt.Orientation, role: int = Qt.ItemDataRole.DisplayRole
    ) -> Any:
        if role == Qt.ItemDataRole.DisplayRole and orientation == Qt.Orientation.Horizontal:
            return self.header[section]
        return None

    def column(self, name: str) -> int:
        return self.header.index(name) if name in self.header else -1


class TablePanel(QWidget):
    def __init__(self) -> None:
        super().__init__()
        self.model = TsvModel()
        self.proxy = QSortFilterProxyModel()
        self.proxy.setSourceModel(self.model)
        self.proxy.setFilterCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
        self.proxy.setFilterKeyColumn(-1)
        self.proxy.setSortRole(Qt.ItemDataRole.UserRole)
        self.filter = QLineEdit()
        self.filter.setPlaceholderText("filter (any column)")
        self.filter.textChanged.connect(self.proxy.setFilterFixedString)
        self.count = QLabel()
        self.proxy.rowsInserted.connect(self._count)
        self.proxy.modelReset.connect(self._count)
        self.filter.textChanged.connect(self._count)
        self.view = QTableView()
        self.view.setModel(self.proxy)
        self.view.setSortingEnabled(True)
        self.view.setSelectionBehavior(QTableView.SelectionBehavior.SelectRows)
        self.view.setAlternatingRowColors(True)
        top = QHBoxLayout()
        top.addWidget(self.filter, 1)
        top.addWidget(self.count)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addLayout(top)
        layout.addWidget(self.view)

    def load(self, path: Path) -> None:
        self.model.load(path)
        self.view.resizeColumnsToContents()
        self._count()

    def _count(self) -> None:
        self.count.setText(f"{self.proxy.rowCount():,} / {self.model.rowCount():,} rows")

    def selected_row(self) -> dict[str, str] | None:
        rows = self.view.selectionModel().selectedRows()
        if not rows:
            return None
        source = self.proxy.mapToSource(rows[0])
        return dict(zip(self.model.header, self.model.rows[source.row()], strict=True))


class ImagePanel(QScrollArea):
    def __init__(self) -> None:
        super().__init__()
        self.label = QLabel("no image")
        self.label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.setWidget(self.label)
        self.setWidgetResizable(True)
        self.pixmap = QPixmap()

    def show_image(self, path: Path | None) -> None:
        if path is not None and path.is_file():
            self.pixmap = QPixmap(str(path))
            self._fit()
        else:
            self.pixmap = QPixmap()
            self.label.setText("no image" if path is None else f"not found: {path.name}")

    def _fit(self) -> None:
        if self.pixmap.isNull():
            return
        size = self.viewport().size()
        self.label.setPixmap(
            self.pixmap.scaled(
                size, Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation
            )
        )

    def resizeEvent(self, event: Any) -> None:  # noqa: N802
        super().resizeEvent(event)
        self._fit()


class ResultsPage(QWidget):
    def __init__(self) -> None:
        super().__init__()
        self.outdir: Path | None = None
        self.run = QComboBox()
        self.run.currentTextChanged.connect(self._show_run)
        reload = QPushButton("Reload")
        reload.clicked.connect(lambda: self.open(self.outdir))
        folder = QPushButton("Open folder")
        folder.clicked.connect(self._open_folder)
        self.where = QLabel()
        self.where.setStyleSheet("color: #888")
        top = QHBoxLayout()
        top.addWidget(QLabel("run"))
        top.addWidget(self.run)
        top.addWidget(reload)
        top.addWidget(folder)
        top.addWidget(self.where, 1)
        self.summary = QPlainTextEdit()
        self.summary.setReadOnly(True)
        self.orthogroups, self.pairs, self.sequences = TablePanel(), TablePanel(), TablePanel()
        self.figure = ImagePanel()
        self.figure_choice = QComboBox()
        self.figure_choice.currentTextChanged.connect(self._show_figure)
        figures = QWidget()
        fl = QVBoxLayout(figures)
        fl.addWidget(self.figure_choice)
        fl.addWidget(self.figure, 1)
        self.dotplot = ImagePanel()
        self.tree = QPlainTextEdit()
        self.tree.setReadOnly(True)
        self.tree.setMaximumHeight(110)
        side = QWidget()
        sl = QVBoxLayout(side)
        sl.addWidget(QLabel("Dotplot of the selected pair"))
        sl.addWidget(self.dotplot, 1)
        sl.addWidget(QLabel("Gene tree of its orthogroup"))
        sl.addWidget(self.tree)
        pairs_split = QSplitter()
        pairs_split.addWidget(self.pairs)
        pairs_split.addWidget(side)
        pairs_split.setSizes([700, 450])
        self.tabs = QTabWidget()
        self.tabs.addTab(self.summary, "Summary")
        self.tabs.addTab(self.orthogroups, "Orthogroups")
        self.tabs.addTab(pairs_split, "Confirmed pairs")
        self.tabs.addTab(self.sequences, "Sequences")
        self.tabs.addTab(figures, "Figures")
        layout = QVBoxLayout(self)
        layout.addLayout(top)
        layout.addWidget(self.tabs, 1)
        self.pairs.view.selectionModel().selectionChanged.connect(self._pair_selected)

    def open(self, outdir: Path | None) -> None:
        self.outdir = outdir
        self.run.blockSignals(True)
        self.run.clear()
        runs = []
        if outdir and outdir.is_dir():
            runs = sorted(p.name for p in outdir.iterdir() if (p / "orthogroups.tsv").is_file())
        self.run.addItems(runs)
        self.run.blockSignals(False)
        self.where.setText(str(outdir) if outdir else "")
        if runs:
            self._show_run(runs[0])
        else:
            self.summary.setPlainText(
                "No results yet: run SPOT, or choose its output folder on the Project page."
            )
        self._figures()

    def _show_run(self, name: str) -> None:
        if not self.outdir or not name:
            return
        folder = self.outdir / name
        summary = folder / "summary.txt"
        self.summary.setPlainText(summary.read_text(encoding="utf-8") if summary.is_file() else "")
        self.orthogroups.load(folder / "orthogroups.tsv")
        self.pairs.load(folder / "pairs.tsv")
        self.pairs.view.selectionModel().selectionChanged.connect(self._pair_selected)
        self.sequences.load(folder / "sequences.tsv")
        self._figures()

    def _figures(self) -> None:
        self.figure_choice.blockSignals(True)
        self.figure_choice.clear()
        if self.outdir:
            names = ["venn_panel.png", "synteny.png", "genome/dotplot.png"]
            if self.run.currentText():
                names.insert(0, f"{self.run.currentText()}/venn.png")
            self.figure_choice.addItems([n for n in names if (self.outdir / n).is_file()])
        self.figure_choice.blockSignals(False)
        self._show_figure(self.figure_choice.currentText())

    def _show_figure(self, name: str) -> None:
        self.figure.show_image(self.outdir / name if self.outdir and name else None)

    def _pair_selected(self) -> None:
        row = self.pairs.selected_row()
        if row is None or not self.outdir:
            return
        og = row.get("orthogroup", "")
        ids = [k for k in row if k.endswith("_id")]
        folder = self.outdir / self.run.currentText() / "syntenic_ortholog" / og
        if len(ids) == 2:
            self.dotplot.show_image(folder / "dotplots" / f"{row[ids[0]]}__{row[ids[1]]}.png")
        tree = folder / f"{og}_tree.txt"
        self.tree.setPlainText(tree.read_text(encoding="utf-8") if tree.is_file() else "no gene tree")

    def _open_folder(self) -> None:
        if self.outdir:
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.outdir)))
