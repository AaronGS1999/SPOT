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

    args = sys.argv[1:] if argv is None else argv
    app = QApplication.instance() or QApplication(sys.argv[:1])
    app.setApplicationName("SPOT")
    app.setOrganizationName("SPOT")
    window = MainWindow()
    if args:
        window.open_config(args[0])
    window.show()
    return int(app.exec())
