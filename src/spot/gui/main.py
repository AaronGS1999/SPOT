"""Main window: a project is one SPOT configuration file, edited on four pages."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import yaml
from PySide6.QtCore import QSettings
from PySide6.QtGui import QAction, QCloseEvent, QKeySequence
from PySide6.QtWidgets import QFileDialog, QMainWindow, QMessageBox, QTabWidget

from .. import __version__, config, diagnose, orthofinder
from ..errors import SpotError
from .parameters import ParametersPage
from .project import ProjectPage
from .results import ResultsPage
from .runner import RunPage

ICON = Path(__file__).with_name("assets") / "spot.png"
YAML = "SPOT configuration (*.yaml *.yml);;All files (*)"
PATH_KEYS = ("genome", "centromeres", "coords")


def _absolute(data: dict[str, Any], base: Path) -> dict[str, Any]:
    """Paths of a configuration made absolute (they are relative to the file in YAML)."""

    def fix(value: Any) -> Any:
        if not value or value == "header":
            return value
        p = Path(str(value)).expanduser()
        return str(p if p.is_absolute() else Path(os.path.normpath(base / p)))

    out = dict(data)
    of = dict(out.get("orthofinder") or {})
    if of.get("root"):
        of["root"] = fix(of["root"])
    for run in (of.get("runs") or {}).values() if isinstance(of.get("runs"), dict) else []:
        if isinstance(run, dict) and run.get("path"):
            run["path"] = fix(run["path"])
    out["orthofinder"] = of
    species = {}
    for key, sp in (out.get("species") or {}).items():
        sp = dict(sp or {})
        for k in PATH_KEYS:
            if sp.get(k):
                sp[k] = fix(sp[k])
        if sp.get("lift"):
            sp["lift"] = [fix(p) for p in sp["lift"]]
        species[key] = sp
    out["species"] = species
    for k in ("outdir", "syri_out"):
        if out.get(k):
            out[k] = fix(out[k])
    if str(out.get("profile", "")).endswith((".yaml", ".yml")):
        out["profile"] = fix(out["profile"])
    return out


def _relative(data: dict[str, Any], base: Path) -> dict[str, Any]:
    """Paths inside ``base`` written relative to it, so a project folder can be moved or shared."""

    def rel(value: Any) -> Any:
        if not isinstance(value, str) or not value or value == "header":
            return value
        p = Path(value)
        try:
            return str(p.relative_to(base)) if p.is_absolute() else value
        except ValueError:
            return value

    out = yaml.safe_load(yaml.safe_dump(data))  # deep copy
    of = out.get("orthofinder") or {}
    if of.get("root"):
        of["root"] = rel(of["root"])
    for sp in (out.get("species") or {}).values():
        for k in PATH_KEYS:
            if sp.get(k):
                sp[k] = rel(sp[k])
        if sp.get("lift"):
            sp["lift"] = [rel(p) for p in sp["lift"]]
    for k in ("outdir", "syri_out", "profile"):
        if out.get(k):
            out[k] = rel(out[k])
    return dict(out)


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.settings = QSettings()
        self.path: Path | None = None
        self.dirty = False
        self.project = ProjectPage()
        self.parameters = ParametersPage()
        self.runner = RunPage()
        self.results = ResultsPage()
        self.tabs = QTabWidget()
        self.tabs.addTab(self.project, "Project")
        self.tabs.addTab(self.parameters, "Parameters")
        self.tabs.addTab(self.runner, "Run")
        self.tabs.addTab(self.results, "Results")
        self.setCentralWidget(self.tabs)
        self._menus()
        self.project.changed.connect(self._touched)
        self.project.use_syri.toggled.connect(self.parameters.set_external_alignment)
        self.parameters.changed.connect(self._touched)
        self.runner.check_button.clicked.connect(self.check)
        self.runner.run_button.clicked.connect(self.run)
        self.runner.finished.connect(self._finished)
        self.tabs.currentChanged.connect(self._tab)
        self.resize(1200, 860)
        self.new()

    def _menus(self) -> None:
        menu = self.menuBar().addMenu("&File")
        bar = self.addToolBar("File")
        bar.setMovable(False)
        for text, shortcut, slot in (
            ("New", QKeySequence.StandardKey.New, self.new),
            ("Open…", QKeySequence.StandardKey.Open, self.open_dialog),
            ("Save", QKeySequence.StandardKey.Save, self.save),
            ("Save as…", QKeySequence.StandardKey.SaveAs, self.save_as),
        ):
            action = QAction(text, self)
            action.setShortcut(shortcut)
            action.triggered.connect(slot)
            menu.addAction(action)
            bar.addAction(action)
        bar.addSeparator()
        run = QAction("Run", self)
        run.setShortcut(QKeySequence("Ctrl+R"))
        run.triggered.connect(self.run)
        bar.addAction(run)
        help_menu = self.menuBar().addMenu("&Help")
        about = QAction("About SPOT", self)
        about.triggered.connect(
            lambda: QMessageBox.about(
                self, "SPOT", f"SPOT {__version__}\nSynteny-based Positional Orthology Test"
            )
        )
        help_menu.addAction(about)
        tools_action = QAction("External tools…", self)
        tools_action.triggered.connect(self.show_tools)
        help_menu.addAction(tools_action)

    def data(self) -> dict[str, Any]:
        data: dict[str, Any] = {}
        self.project.store(data)
        self.parameters.store(data)
        return data

    def load(self, data: dict[str, Any]) -> None:
        self.project.load(data)
        self.parameters.load(data)
        self.parameters.set_external_alignment(self.project.use_syri.isChecked())
        self.results.open(Path(data["outdir"]) if data.get("outdir") else None)
        self.dirty = False
        self._title()

    def new(self) -> None:
        if not self._may_discard():
            return
        self.path = None
        self.load({"profile": "default", "threads": "auto", "species": {"": {}, " ": {}}})

    def open_dialog(self) -> None:
        if not self._may_discard():
            return
        start = str(self.settings.value("last_dir", str(Path.home())))
        path, _ = QFileDialog.getOpenFileName(self, "Open SPOT configuration", start, YAML)
        if path:
            self.open_config(path)

    def open_config(self, path: str | Path) -> None:
        path = Path(path).expanduser().resolve()
        try:
            raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        except (OSError, yaml.YAMLError) as e:
            QMessageBox.critical(self, "SPOT", f"Cannot read {path}:\n{e}")
            return
        if not isinstance(raw, dict):
            QMessageBox.critical(self, "SPOT", f"{path} is not a SPOT configuration")
            return
        self.path = path
        self.settings.setValue("last_dir", str(path.parent))
        self.load(_absolute(raw, path.parent))

    def save(self) -> bool:
        if self.path is None:
            return self.save_as()
        text = yaml.safe_dump(_relative(self.data(), self.path.parent), sort_keys=False, allow_unicode=True)
        self.path.write_text(text, encoding="utf-8")
        self.dirty = False
        self._title()
        self.statusBar().showMessage(f"Saved {self.path}", 4000)
        return True

    def save_as(self) -> bool:
        start = str(self.path or Path(str(self.settings.value("last_dir", str(Path.home())))) / "spot.yaml")
        path, _ = QFileDialog.getSaveFileName(self, "Save SPOT configuration", start, YAML)
        if not path:
            return False
        self.path = Path(path)
        self.settings.setValue("last_dir", str(self.path.parent))
        return self.save()

    def _touched(self) -> None:
        if not self.dirty:
            self.dirty = True
            self._title()

    def _title(self) -> None:
        name = self.path.name if self.path else "untitled"
        self.setWindowTitle(f"SPOT {__version__} - {name}{' *' if self.dirty else ''}")

    def _may_discard(self) -> bool:
        if not self.dirty:
            return True
        answer = QMessageBox.question(
            self,
            "SPOT",
            "Save the changes to this project?",
            QMessageBox.StandardButton.Save
            | QMessageBox.StandardButton.Discard
            | QMessageBox.StandardButton.Cancel,
        )
        if answer == QMessageBox.StandardButton.Save:
            return self.save()
        return answer == QMessageBox.StandardButton.Discard

    def closeEvent(self, event: QCloseEvent) -> None:  # noqa: N802
        if self.runner.process is not None:
            answer = QMessageBox.question(self, "SPOT", "SPOT is still running. Stop it and quit?")
            if answer != QMessageBox.StandardButton.Yes:
                event.ignore()
                return
            self.runner.cancel()
        if self._may_discard():
            event.accept()
        else:
            event.ignore()

    def show_tools(self) -> None:
        from .. import tools

        lines = []
        for s in tools.check(None):
            where = s.path or "not found"
            lines.append(
                f"{s.name:<9} {s.version or '-':<14} {where}{f'  ({s.problem})' if s.problem else ''}"
            )
        box = QMessageBox(self)
        box.setWindowTitle("External tools")
        box.setText("Tools SPOT uses, as found by this app:")
        box.setDetailedText("\n".join(lines))
        box.setInformativeText("\n".join(lines))
        box.exec()

    def _snapshot(self) -> Path | None:
        """Write the configuration the engine will run (absolute paths) into the output folder."""
        data = self.data()
        if not data.get("outdir"):
            QMessageBox.warning(self, "SPOT", "Choose an output folder on the Project page.")
            self.tabs.setCurrentWidget(self.project)
            return None
        outdir = Path(data["outdir"])
        outdir.mkdir(parents=True, exist_ok=True)
        path = outdir / "spot_app.yaml"
        path.write_text(yaml.safe_dump(data, sort_keys=False, allow_unicode=True), encoding="utf-8")
        return path

    def check(self) -> bool:
        self.tabs.setCurrentWidget(self.runner)
        self.runner.log.clear()
        path = self._snapshot()
        if path is None:
            return False
        try:
            cfg = config.load(path)
            runs = orthofinder.load_all(cfg)
        except SpotError as e:
            self.runner.message(f"Configuration problem:\n{e}")
            self.runner.status.setText("Fix the configuration")
            return False
        self.runner.message(f"Configuration OK ({path})\n")
        for run in runs:
            counts = ", ".join(f"{k} {v:,}" for k, v in run.summary().items())
            columns = ", ".join(f"{sp} <- {col}" for sp, col in run.columns.items())
            self.runner.message(f"{run.name}: {columns}\n    {counts}")
        blocking = False
        paths = [run.path for run in runs]
        for key, sp in cfg.species.items():
            if isinstance(sp.coords, Path):
                continue
            names = diagnose.headers(paths, sp.match)
            report = diagnose.species_report(names, sp.header_regex, sp.genome, sp.lift, sp.strip_suffix)
            self.runner.message(f"\n{key} ({sp.label}):\n    " + report.text().replace("\n", "\n    "))
            if report.on_genome == 0:
                blocking = True
        if blocking:
            self.runner.message(
                "\nNo sequence of a species falls on its genome: check the genome and lift files."
            )
            self.runner.status.setText("Fix the configuration")
            return False
        self.runner.message(f"\nthreads: {cfg.n_threads}, profile: {cfg.profile}")
        self.runner.status.setText("Ready to run")
        return True

    def run(self) -> None:
        if self.runner.process is not None:
            return
        if not self.check():
            return
        path = self._snapshot()
        if path is not None:
            self.runner.start(path)

    def _finished(self, ok: bool) -> None:
        outdir = self.data().get("outdir")
        self.results.open(Path(outdir) if outdir else None)
        if ok:
            self.tabs.setCurrentWidget(self.results)

    def _tab(self, index: int) -> None:
        if self.tabs.widget(index) is self.results and self.results.outdir is None:
            outdir = self.data().get("outdir")
            self.results.open(Path(outdir) if outdir else None)
