"""Run page: checks the configuration, launches the engine in its own process and follows it.

The engine runs as ``python -m spot --progress run -c <config>`` in a separate process group, so the
window never freezes, a crash cannot take the app down and Cancel stops BLAST, Gepard and every other
child process at once.
"""

from __future__ import annotations

import os
import signal
import sys
import time
from pathlib import Path

from PySide6.QtCore import QProcess, QProcessEnvironment, QTimer, Signal
from PySide6.QtGui import QFont, QTextCursor
from PySide6.QtWidgets import (
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from ..progress import STAGES, parse

STAGE_LABELS = {
    "locate": "Locate sequences",
    "align": "Genome alignment (minimap2)",
    "syri": "Synteny (SyRI)",
    "test": "Positional test",
    "validate": "Flank validation (BLASTN)",
    "export": "Write results",
    "dotplots": "Dotplots (Gepard)",
    "plots": "Figures",
}


class RunPage(QWidget):
    finished = Signal(bool)  # True when the run succeeded
    started = Signal()

    def __init__(self) -> None:
        super().__init__()
        self.check_button = QPushButton("Check configuration")
        self.run_button = QPushButton("Run SPOT")
        self.run_button.setDefault(True)
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.setEnabled(False)
        self.cancel_button.clicked.connect(self.cancel)
        self.status = QLabel("Ready")
        buttons = QHBoxLayout()
        for b in (self.check_button, self.run_button, self.cancel_button):
            buttons.addWidget(b)
        buttons.addWidget(self.status, 1)
        grid = QGridLayout()
        self.bars: dict[str, QProgressBar] = {}
        for i, stage in enumerate(STAGES):
            bar = QProgressBar()
            bar.setRange(0, 1)
            bar.setValue(0)
            bar.setFormat("")
            grid.addWidget(QLabel(STAGE_LABELS[stage]), i, 0)
            grid.addWidget(bar, i, 1)
            self.bars[stage] = bar
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumBlockCount(20000)
        self.log.setFont(QFont("Menlo" if sys.platform == "darwin" else "Monospace", 10))
        layout = QVBoxLayout(self)
        layout.addLayout(buttons)
        layout.addLayout(grid)
        layout.addWidget(self.log, 1)
        self.process: QProcess | None = None
        self._buffer = ""
        self._killer = QTimer(self)
        self._killer.setSingleShot(True)
        self._killer.timeout.connect(self._kill)
        self._stage, self._since = "", 0.0
        self._clock = QTimer(self)
        self._clock.setInterval(1000)
        self._clock.timeout.connect(self._tick)

    def message(self, text: str) -> None:
        self.log.appendPlainText(text)
        self.log.moveCursor(QTextCursor.MoveOperation.End)

    def reset(self) -> None:
        for bar in self.bars.values():
            bar.setRange(0, 1)
            bar.setValue(0)
            bar.setFormat("")

    def start(self, config_path: Path) -> None:
        if self.process is not None:
            return
        self.reset()
        self.log.clear()
        self.message(f"$ spot run -c {config_path}")
        process = QProcess(self)
        env = QProcessEnvironment.systemEnvironment()
        env.insert("PYTHONUNBUFFERED", "1")
        process.setProcessEnvironment(env)
        process.setWorkingDirectory(str(config_path.parent))  # never the app's read-only folder
        process.setProgram(sys.executable)
        process.setArguments(["-m", "spot", "--progress", "run", "-c", str(config_path)])
        process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        if hasattr(process, "setChildProcessModifier") and os.name == "posix":
            process.setChildProcessModifier(os.setsid)  # own process group: Cancel stops the children too
        process.readyReadStandardOutput.connect(self._read)
        process.finished.connect(self._done)
        process.errorOccurred.connect(self._error)
        self.process = process
        self.run_button.setEnabled(False)
        self.check_button.setEnabled(False)
        self.cancel_button.setEnabled(True)
        self.status.setText("Running…")
        self._clock.start()
        process.start()
        self.started.emit()

    def _read(self) -> None:
        if self.process is None:
            return
        self._buffer += bytes(self.process.readAllStandardOutput().data()).decode(errors="replace")
        *lines, self._buffer = self._buffer.split("\n")
        for line in lines:
            update = parse(line)
            if update is None:
                if line.strip():
                    self.message(line.rstrip())
                continue
            self._progress(*update)

    def _progress(self, stage: str, done: int, total: int) -> None:
        bar = self.bars.get(stage)
        if bar is None:
            return
        finished = done >= total
        if finished:
            bar.setRange(0, 1)
            bar.setValue(1)
            bar.setFormat("done")
        elif done == 0 or total == 1:
            bar.setRange(0, 0)  # busy: moving bar until the step reports progress
        else:
            bar.setRange(0, total)
            bar.setValue(done)
            bar.setFormat(f"{done:,} / {total:,}")
        if not finished:
            if stage != self._stage:
                self._stage, self._since = stage, time.monotonic()
            self._tick()

    def _tick(self) -> None:
        if self.process is None or not self._stage:
            return
        seconds = int(time.monotonic() - self._since)
        self.status.setText(
            f"{STAGE_LABELS.get(self._stage, self._stage)}… {seconds // 60}:{seconds % 60:02d}"
        )

    def _done(self, code: int, status: QProcess.ExitStatus) -> None:
        self._read()
        if self._buffer.strip():
            self.message(self._buffer.rstrip())
        self._buffer = ""
        self._killer.stop()
        self._clock.stop()
        self._stage = ""
        for bar in self.bars.values():
            if bar.maximum() == 0:  # a stage still busy when the run stopped
                bar.setRange(0, 1)
                bar.setValue(0)
        ok = status == QProcess.ExitStatus.NormalExit and code == 0
        self.status.setText("Finished" if ok else f"Stopped (exit code {code})")
        self.process = None
        self.run_button.setEnabled(True)
        self.check_button.setEnabled(True)
        self.cancel_button.setEnabled(False)
        self.finished.emit(ok)

    def _error(self, error: QProcess.ProcessError) -> None:
        if error == QProcess.ProcessError.FailedToStart:
            self.message(f"could not start {sys.executable}")

    def cancel(self) -> None:
        if self.process is None:
            return
        self.status.setText("Cancelling…")
        self._signal(signal.SIGTERM)
        self._killer.start(5000)

    def _kill(self) -> None:
        if self.process is not None:
            self._signal(signal.SIGKILL if hasattr(signal, "SIGKILL") else signal.SIGTERM)

    def _signal(self, sig: int) -> None:
        if self.process is None:
            return
        pid = int(self.process.processId())
        if pid and os.name == "posix":
            try:
                os.killpg(pid, sig)
                return
            except (ProcessLookupError, PermissionError):
                pass
        self.process.kill() if sig != signal.SIGTERM else self.process.terminate()
