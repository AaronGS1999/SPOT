"""Parameters page: profile and every engine setting, generated from the configuration models."""

from __future__ import annotations

from typing import Any

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QComboBox, QHBoxLayout, QLabel, QScrollArea, QVBoxLayout, QWidget

from .. import config
from ..errors import SpotError
from .forms import SECTION_TITLES, SectionForm


class ParametersPage(QWidget):
    changed = Signal()

    def __init__(self) -> None:
        super().__init__()
        self.profile = QComboBox()
        self.profile.addItems(list(config.builtin_profiles()))
        self.description = QLabel()
        self.description.setWordWrap(True)
        self.description.setStyleSheet("color: #666")
        top = QHBoxLayout()
        top.addWidget(QLabel("profile"))
        top.addWidget(self.profile)
        top.addWidget(self.description, 1)
        self.forms: dict[str, SectionForm] = {}
        body = QWidget()
        column = QVBoxLayout(body)
        for section, (title, extra) in SECTION_TITLES.items():
            form = SectionForm(section, config.SECTIONS[section], title, extra)
            form.edited.connect(self.changed)
            self.forms[section] = form
            column.addWidget(form)
        column.addStretch()
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(body)
        layout = QVBoxLayout(self)
        layout.addLayout(top)
        layout.addWidget(scroll)
        self._user: dict[str, dict[str, Any]] = {}
        self._external_alignment = False
        self.profile.currentTextChanged.connect(self._profile_changed)
        self.changed.connect(self.update_dependencies)

    def set_external_alignment(self, external: bool) -> None:
        """minimap2 and SyRI settings only matter when SPOT aligns the genomes itself."""
        self._external_alignment = external
        self.update_dependencies()

    def update_dependencies(self) -> None:
        validation = self.forms["validation"].rows["enabled"].editor.value()
        self.forms["validation"].set_rows_enabled(bool(validation), keep=("enabled",))
        self.forms["blastn"].set_active(bool(validation), "flank validation is off")
        pairs = self.forms["dotplots"].rows["pairs"].editor.value()
        self.forms["gepard"].set_active(bool(pairs), "no dotplots of the pairs")
        for tool in ("minimap2", "syri"):
            self.forms[tool].set_active(not self._external_alignment, "an existing syri.out is used")

    def _profile_changed(self, name: str) -> None:
        self._user = self.overrides()
        self._show(name)
        self.changed.emit()

    def _show(self, name: str) -> None:
        try:
            description, _ = config.load_profile(name)
            settings = config.profile_settings(name)
        except SpotError as e:
            self.description.setText(str(e))
            return
        self.description.setText(description)
        for section, form in self.forms.items():
            form.show(settings.get(section, {}), self._user.get(section, {}))
        self.update_dependencies()

    def load(self, data: dict[str, Any]) -> None:
        self._user = {s: dict(data.get(s) or {}) for s in self.forms}
        name = str(data.get("profile") or "default")
        self.profile.blockSignals(True)
        if self.profile.findText(name) < 0:
            self.profile.addItem(name)
        self.profile.setCurrentText(name)
        self.profile.blockSignals(False)
        self._show(name)

    def overrides(self) -> dict[str, dict[str, Any]]:
        return {section: form.overrides() for section, form in self.forms.items()}

    def store(self, data: dict[str, Any]) -> None:
        data["profile"] = self.profile.currentText()
        for section, values in self.overrides().items():
            if values:
                data[section] = values
            else:
                data.pop(section, None)
