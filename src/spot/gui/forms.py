"""Forms built from the configuration models, so every setting of the engine is editable in the app.

Each field shows the value in force and where it comes from: the profile (grey tag) or a value set in
this project (bold tag, with a reset button that returns it to the profile's value). Only the values set
in the project are written to the configuration file, exactly like a hand-written YAML.
"""

from __future__ import annotations

import math
import types
import typing
from typing import Any, Literal

import yaml
from annotated_types import Ge, Gt, Le
from pydantic import BaseModel
from pydantic.fields import FieldInfo
from PySide6.QtCore import Signal
from PySide6.QtGui import QDoubleValidator
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QPlainTextEdit,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from ..config import SYRI_CLASSES
from .widgets import ColorButton

INT_MAX = 2_000_000_000


class Editor(QWidget):
    """One field: get/set its value and say when the user changes it."""

    edited = Signal()

    def value(self) -> Any:
        raise NotImplementedError

    def set_value(self, value: Any) -> None:
        raise NotImplementedError


def _bounds(field: FieldInfo) -> tuple[float | None, float | None]:
    lo = hi = None
    for m in field.metadata:
        if isinstance(m, Ge | Gt):
            lo = float(m.ge if isinstance(m, Ge) else m.gt)  # type: ignore[arg-type]
        if isinstance(m, Le):
            hi = float(m.le)  # type: ignore[arg-type]
    return lo, hi


class BoolEditor(Editor):
    def __init__(self) -> None:
        super().__init__()
        self.box = QCheckBox()
        self.box.toggled.connect(self.edited)
        _wrap(self, self.box)

    def value(self) -> bool:
        return self.box.isChecked()

    def set_value(self, value: Any) -> None:
        self.box.blockSignals(True)
        self.box.setChecked(bool(value))
        self.box.blockSignals(False)


class IntEditor(Editor):
    def __init__(self, field: FieldInfo) -> None:
        super().__init__()
        lo, hi = _bounds(field)
        self.box = QSpinBox()
        self.box.setRange(int(lo) if lo is not None else -INT_MAX, int(hi) if hi is not None else INT_MAX)
        self.box.setGroupSeparatorShown(True)
        self.box.setMinimumWidth(140)
        self.box.valueChanged.connect(self.edited)
        _wrap(self, self.box)

    def value(self) -> int:
        return int(self.box.value())

    def set_value(self, value: Any) -> None:
        self.box.blockSignals(True)
        self.box.setValue(int(value))
        self.box.blockSignals(False)


class FloatEditor(Editor):
    """A float typed as text, so values such as 1e-10 keep their notation."""

    def __init__(self, field: FieldInfo) -> None:
        super().__init__()
        lo, hi = _bounds(field)
        self.edit = QLineEdit()
        self.edit.setMaximumWidth(160)
        validator = QDoubleValidator(
            lo if lo is not None else -math.inf, hi if hi is not None else math.inf, 12
        )
        validator.setNotation(QDoubleValidator.Notation.ScientificNotation)
        self.edit.setValidator(validator)
        self.edit.textEdited.connect(self.edited)
        _wrap(self, self.edit)

    def value(self) -> float:
        try:
            return float(self.edit.text())
        except ValueError:
            return math.nan

    def set_value(self, value: Any) -> None:
        self.edit.setText(f"{float(value):g}")


class ChoiceEditor(Editor):
    def __init__(self, choices: list[str]) -> None:
        super().__init__()
        self.box = QComboBox()
        self.box.addItems(choices)
        self.box.currentTextChanged.connect(self.edited)
        _wrap(self, self.box)

    def value(self) -> str:
        return str(self.box.currentText())

    def set_value(self, value: Any) -> None:
        self.box.blockSignals(True)
        self.box.setCurrentText(str(value))
        self.box.blockSignals(False)


class TextEditor(Editor):
    def __init__(self) -> None:
        super().__init__()
        self.edit = QLineEdit()
        self.edit.textEdited.connect(self.edited)
        _wrap(self, self.edit)

    def value(self) -> str:
        return str(self.edit.text())

    def set_value(self, value: Any) -> None:
        self.edit.setText("" if value is None else str(value))


class ClassesEditor(Editor):
    """SyRI classes as a row of check boxes."""

    def __init__(self) -> None:
        super().__init__()
        self.boxes = {name: QCheckBox(name) for name in SYRI_CLASSES}
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        for box in self.boxes.values():
            box.toggled.connect(self.edited)
            layout.addWidget(box)
        layout.addStretch()

    def value(self) -> list[str]:
        return [name for name, box in self.boxes.items() if box.isChecked()]

    def set_value(self, value: Any) -> None:
        chosen = set(value or [])
        for name, box in self.boxes.items():
            box.blockSignals(True)
            box.setChecked(name in chosen)
            box.blockSignals(False)


class HomologyEditor(Editor):
    """auto, or one 'reference query' pair per line."""

    def __init__(self) -> None:
        super().__init__()
        self.mode = QComboBox()
        self.mode.addItems(["auto", "manual"])
        self.pairs = QPlainTextEdit()
        self.pairs.setToolTip(
            "One homologous pair per line: reference chromosome, a space, query chromosome."
        )
        self.pairs.setMaximumHeight(80)
        self.mode.currentTextChanged.connect(self._mode)
        self.mode.currentTextChanged.connect(self.edited)
        self.pairs.textChanged.connect(self.edited)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.mode)
        layout.addWidget(self.pairs, 1)
        layout.addStretch()
        self._mode(self.mode.currentText())

    def _mode(self, mode: str) -> None:
        self.pairs.setVisible(mode == "manual")

    def value(self) -> str | list[list[str]]:
        if self.mode.currentText() == "auto":
            return "auto"
        return [line.split()[:2] for line in self.pairs.toPlainText().splitlines() if len(line.split()) >= 2]

    def set_value(self, value: Any) -> None:
        self.blockSignals(True)
        if value == "auto" or not value:
            self.mode.setCurrentText("auto")
        else:
            self.mode.setCurrentText("manual")
            self.pairs.setPlainText("\n".join(f"{a} {b}" for a, b in value))
        self.blockSignals(False)


class ColorEditor(Editor):
    def __init__(self) -> None:
        super().__init__()
        self.button = ColorButton()
        self.button.changed.connect(self.edited)
        _wrap(self, self.button)

    def value(self) -> str:
        return self.button.color() or "#000000"

    def set_value(self, value: Any) -> None:
        self.button.set_color(str(value) if value else None)


class OptionsEditor(QWidget):
    """Extra command-line options of a tool as a two-column table (option, value)."""

    edited = Signal()

    def __init__(self, tool: str) -> None:
        super().__init__()
        self.table = QTableWidget(0, 2)
        self.table.setHorizontalHeaderLabels(["option", "value"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.verticalHeader().setVisible(False)
        self.table.setMaximumHeight(110)
        self.table.itemChanged.connect(lambda _: self.edited.emit())
        self.setToolTip(
            f"Any other {tool} option, written as in its manual without the dash (e.g. max_hsps with "
            "value 50 gives -max_hsps 50). Leave the value empty for a flag."
        )
        add, remove = QToolButton(), QToolButton()
        add.setText("+")
        remove.setText("−")
        add.clicked.connect(self._add)
        remove.clicked.connect(self._remove)
        buttons = QVBoxLayout()
        buttons.addWidget(add)
        buttons.addWidget(remove)
        buttons.addStretch()
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.table, 1)
        layout.addLayout(buttons)

    def _add(self) -> None:
        self.table.insertRow(self.table.rowCount())
        self.table.editItem(
            self.table.item(self.table.rowCount() - 1, 0) or self._blank(self.table.rowCount() - 1)
        )

    def _blank(self, row: int) -> QTableWidgetItem:
        item = QTableWidgetItem("")
        self.table.setItem(row, 0, item)
        return item

    def _remove(self) -> None:
        rows = sorted({i.row() for i in self.table.selectedIndexes()}, reverse=True)
        for r in rows or [self.table.rowCount() - 1]:
            if r >= 0:
                self.table.removeRow(r)
        self.edited.emit()

    def options(self) -> dict[str, Any]:
        out: dict[str, Any] = {}
        for r in range(self.table.rowCount()):
            key_item, value_item = self.table.item(r, 0), self.table.item(r, 1)
            key = key_item.text().strip().lstrip("-") if key_item else ""
            if not key:
                continue
            text = value_item.text().strip() if value_item else ""
            out[key] = True if text == "" else yaml.safe_load(text)
        return out

    def set_options(self, options: dict[str, Any]) -> None:
        self.table.blockSignals(True)
        self.table.setRowCount(0)
        for key, value in options.items():
            r = self.table.rowCount()
            self.table.insertRow(r)
            self.table.setItem(r, 0, QTableWidgetItem(str(key)))
            self.table.setItem(r, 1, QTableWidgetItem("" if value is True else str(value)))
        self.table.blockSignals(False)


def _wrap(owner: QWidget, child: QWidget) -> None:
    layout = QHBoxLayout(owner)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.addWidget(child)
    layout.addStretch()


def editor_for(name: str, field: FieldInfo) -> Editor:
    annotation = field.annotation
    if name == "types":
        return ClassesEditor()
    if name == "homology":
        return HomologyEditor()
    extra = field.json_schema_extra
    if isinstance(extra, dict) and extra.get("color"):
        return ColorEditor()
    if annotation is bool:
        return BoolEditor()
    if annotation is int:
        return IntEditor(field)
    if annotation is float:
        return FloatEditor(field)
    if typing.get_origin(annotation) is Literal:
        return ChoiceEditor([str(a) for a in typing.get_args(annotation)])
    if isinstance(annotation, types.UnionType) or typing.get_origin(annotation) is typing.Union:
        return TextEditor()
    return TextEditor()


class FieldRow(QWidget):
    edited = Signal()

    def __init__(self, name: str, field: FieldInfo) -> None:
        super().__init__()
        self.name = name
        self.editor = editor_for(name, field)
        self.tag = QLabel()
        self.tag.setMinimumWidth(56)
        self.reset = QToolButton()
        self.reset.setText("↺")
        self.reset.setToolTip("Back to the profile value")
        self.reset.clicked.connect(self._reset)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.editor, 1)
        layout.addWidget(self.tag)
        layout.addWidget(self.reset)
        self.editor.edited.connect(self._edited)
        self.default: Any = None
        self.user = False
        if field.description:
            self.setToolTip(field.description)

    def show_value(self, default: Any, user: Any = None, user_set: bool = False) -> None:
        self.default = default
        self.user = user_set
        self.editor.set_value(user if user_set else default)
        self._style()

    def _edited(self) -> None:
        self.user = True
        self._style()
        self.edited.emit()

    def _reset(self) -> None:
        self.user = False
        self.editor.set_value(self.default)
        self._style()
        self.edited.emit()

    def _style(self) -> None:
        if self.user:
            self.tag.setText("<b>project</b>")
            self.tag.setToolTip("Set in this project")
        else:
            self.tag.setText("<span style='color:#888'>profile</span>")
            self.tag.setToolTip("Value of the profile")
        self.reset.setEnabled(self.user)


class SectionForm(QGroupBox):
    """All fields of one configuration section, plus free options for the external tools."""

    edited = Signal()

    def __init__(self, section: str, model: type[BaseModel], title: str, extra_options: bool = False) -> None:
        super().__init__(title)
        self.section, self.model = section, model
        layout = QFormLayout(self)
        layout.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        self.rows: dict[str, FieldRow] = {}
        for name, field in model.model_fields.items():
            row = FieldRow(name, field)
            row.edited.connect(self.edited)
            label = QLabel(name)
            if field.description:
                label.setToolTip(field.description)
            layout.addRow(label, row)
            self.rows[name] = row
        self.layout_ = layout
        self.base_title = title
        self.extra: OptionsEditor | None = None
        if extra_options:
            self.extra = OptionsEditor(section)
            self.extra.edited.connect(self.edited)
            label = QLabel("other options")
            label.setToolTip(self.extra.toolTip())
            layout.addRow(label, self.extra)

    def show(self, profile: dict[str, Any], user: dict[str, Any]) -> None:  # type: ignore[override]
        for name, row in self.rows.items():
            row.show_value(profile.get(name), user.get(name), name in user)
        if self.extra is not None:
            known = set(self.rows)
            self.extra.set_options({k: v for k, v in user.items() if k not in known})
        super().show()

    def set_active(self, active: bool, reason: str = "") -> None:
        """Grey out a section whose settings are not used with the current choices."""
        self.setEnabled(active)
        self.setTitle(self.base_title if active else f"{self.base_title}  (not used: {reason})")

    def set_rows_enabled(self, enabled: bool, keep: tuple[str, ...] = ()) -> None:
        for name, row in self.rows.items():
            if name not in keep:
                row.setEnabled(enabled)
                label = self.layout_.labelForField(row)
                if label is not None:
                    label.setEnabled(enabled)

    def overrides(self) -> dict[str, Any]:
        out = {name: row.editor.value() for name, row in self.rows.items() if row.user}
        if self.extra is not None:
            out.update(self.extra.options())
        return out


SECTION_TITLES: dict[str, tuple[str, bool]] = {
    "anchors": ("Anchors and positional test", False),
    "stats": ("Statistics", False),
    "validation": ("Flank validation (BLASTN)", False),
    "blastn": ("BLASTN", True),
    "dotplots": ("Dotplots", False),
    "gepard": ("Gepard (dotplots of the confirmed pairs)", True),
    "minimap2": ("minimap2 (whole-genome alignment)", True),
    "syri": ("SyRI (whole-genome alignment)", True),
    "figures": ("Figure colours", False),
}
