"""Small reusable widgets: path pickers, path lists and colour buttons."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Signal
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QColorDialog,
    QFileDialog,
    QHBoxLayout,
    QLineEdit,
    QListWidget,
    QPushButton,
    QToolButton,
    QVBoxLayout,
    QWidget,
)


class PathEdit(QWidget):
    """A path field with a Browse button. ``mode`` is ``file``, ``dir`` or ``save``."""

    changed = Signal()

    def __init__(
        self, mode: str = "file", caption: str = "Choose", filters: str = "", parent: QWidget | None = None
    ):
        super().__init__(parent)
        self.mode, self.caption, self.filters = mode, caption, filters
        self.edit = QLineEdit()
        self.edit.setClearButtonEnabled(True)
        self.edit.textChanged.connect(self.changed)
        button = QToolButton()
        button.setText("…")
        button.setToolTip("Browse")
        button.clicked.connect(self._browse)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.edit)
        layout.addWidget(button)

    def _browse(self) -> None:
        start = self.edit.text() or str(Path.home())
        if self.mode == "dir":
            path = QFileDialog.getExistingDirectory(self, self.caption, start)
        elif self.mode == "save":
            path, _ = QFileDialog.getSaveFileName(self, self.caption, start, self.filters)
        else:
            path, _ = QFileDialog.getOpenFileName(self, self.caption, start, self.filters)
        if path:
            self.edit.setText(path)

    def path(self) -> Path | None:
        text = self.edit.text().strip()
        return Path(text).expanduser() if text else None

    def set_path(self, path: str | Path | None) -> None:
        self.edit.setText("" if path is None else str(path))

    def setPlaceholderText(self, text: str) -> None:  # noqa: N802 (Qt naming)
        self.edit.setPlaceholderText(text)


class PathList(QWidget):
    """An ordered list of files (e.g. lift steps) with add, remove and reorder buttons."""

    changed = Signal()

    def __init__(self, caption: str = "Add files", parent: QWidget | None = None):
        super().__init__(parent)
        self.caption = caption
        self.list = QListWidget()
        self.list.setMaximumHeight(72)
        buttons = QVBoxLayout()
        for text, slot in (("+", self._add), ("−", self._remove), ("↑", self._up), ("↓", self._down)):
            b = QToolButton()
            b.setText(text)
            b.clicked.connect(slot)
            buttons.addWidget(b)
        buttons.addStretch()
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.list)
        layout.addLayout(buttons)

    def paths(self) -> list[Path]:
        return [Path(self.list.item(i).text()) for i in range(self.list.count())]

    def set_paths(self, paths: list[str] | list[Path]) -> None:
        self.list.clear()
        self.list.addItems([str(p) for p in paths])

    def _add(self) -> None:
        files, _ = QFileDialog.getOpenFileNames(self, self.caption, str(Path.home()))
        self.list.addItems(files)
        if files:
            self.changed.emit()

    def _remove(self) -> None:
        for item in self.list.selectedItems():
            self.list.takeItem(self.list.row(item))
        self.changed.emit()

    def _move(self, step: int) -> None:
        row = self.list.currentRow()
        if row < 0 or not 0 <= row + step < self.list.count():
            return
        item = self.list.takeItem(row)
        self.list.insertItem(row + step, item)
        self.list.setCurrentRow(row + step)
        self.changed.emit()

    def _up(self) -> None:
        self._move(-1)

    def _down(self) -> None:
        self._move(1)


class ColorButton(QPushButton):
    """A button showing a colour; clicking it opens a colour picker. Empty means 'default'."""

    changed = Signal()

    def __init__(self, color: str | None = None, parent: QWidget | None = None):
        super().__init__(parent)
        self._color: str | None = None
        self.setFixedWidth(64)
        self.clicked.connect(self._pick)
        self.set_color(color)

    def color(self) -> str | None:
        return self._color

    def set_color(self, color: str | None) -> None:
        self._color = color or None
        if self._color:
            self.setText("")
            self.setStyleSheet(
                f"background-color: {self._color}; border: 1px solid #888; border-radius: 3px;"
            )
        else:
            self.setText("auto")
            self.setStyleSheet("")

    def _pick(self) -> None:
        chosen = QColorDialog.getColor(QColor(self._color or "#C1121F"), self, "Run colour")
        if chosen.isValid():
            self.set_color(chosen.name().upper())
            self.changed.emit()
