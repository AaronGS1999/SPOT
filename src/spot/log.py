import logging
import sys
from pathlib import Path

NAME = "spot"


def add_file(path: Path) -> logging.Handler:
    handler = logging.FileHandler(path, mode="w")
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)-7s %(message)s", "%Y-%m-%d %H:%M:%S"))
    logging.getLogger(NAME).addHandler(handler)
    return handler


def setup(verbose: bool = False) -> logging.Logger:
    log = logging.getLogger(NAME)
    log.setLevel(logging.DEBUG)
    log.handlers.clear()
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)-7s %(message)s", "%H:%M:%S"))
    log.addHandler(handler)
    log.propagate = False
    return log
