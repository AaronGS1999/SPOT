import logging
import sys
from pathlib import Path

NAME = "spot"
PROGRESS = "spot.progress"


class _NoProgress(logging.Filter):
    def filter(self, record: logging.LogRecord) -> bool:
        return record.name != PROGRESS


def add_file(path: Path) -> logging.Handler:
    handler = logging.FileHandler(path, mode="w", encoding="utf-8")
    handler.setLevel(logging.DEBUG)
    handler.addFilter(_NoProgress())
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)-7s %(message)s", "%Y-%m-%d %H:%M:%S"))
    logging.getLogger(NAME).addHandler(handler)
    return handler


def setup(verbose: bool = False, progress: bool = False) -> logging.Logger:
    """Console logging. Progress lines (``@progress stage done total``) are only shown with ``progress``,
    which is how the desktop app follows a run."""
    log = logging.getLogger(NAME)
    log.setLevel(logging.DEBUG)
    log.handlers.clear()
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    if not progress:
        handler.addFilter(_NoProgress())
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)-7s %(message)s", "%H:%M:%S"))
    log.addHandler(handler)
    log.propagate = False
    return log
