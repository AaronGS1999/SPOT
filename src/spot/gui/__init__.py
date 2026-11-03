"""SPOT desktop app (Qt). Start it with ``spot gui`` or ``spot-gui``."""

from __future__ import annotations

import sys


def main(argv: list[str] | None = None) -> int:
    try:
        from PySide6.QtWidgets import QApplication
    except ImportError:
        print(
            "The SPOT app needs PySide6: conda install -c conda-forge pyside6 (or pip install 'spot[gui]')",
            file=sys.stderr,
        )
        return 1
    from .main import MainWindow

    args = [a for a in (sys.argv[1:] if argv is None else argv) if not a.startswith("-psn_")]
    existing = QApplication.instance()
    app = existing if isinstance(existing, QApplication) else QApplication(sys.argv[:1])
    app.setApplicationName("SPOT")
    app.setApplicationDisplayName("SPOT")
    app.setOrganizationName("SPOT")
    from PySide6.QtGui import QIcon

    from .main import ICON

    if ICON.is_file():
        app.setWindowIcon(QIcon(str(ICON)))
    window = MainWindow()
    if args:
        window.open_config(args[0])
    window.show()
    return int(app.exec())
