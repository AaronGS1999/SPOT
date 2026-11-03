"""Content digests and per-stage checkpoints.

Digests of large inputs (genomes, syri.out) are computed once and remembered by path, size and
modification time in ``~/.cache/spot/digests.json`` (or ``$SPOT_CACHE``), so a re-run does not hash
gigabytes again. Checkpoints are pickles in ``<outdir>/.spot/`` named after the stage and keyed by a
SHA-256 of everything the stage depends on; a different key means the stage runs again.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import pickle
import tempfile
from collections.abc import Callable
from pathlib import Path
from typing import Any, TypeVar

log = logging.getLogger("spot")

T = TypeVar("T")
CHUNK = 1 << 20
CHECKPOINT_VERSION = 1


def cache_dir() -> Path:
    env = os.environ.get("SPOT_CACHE")
    return Path(env) if env else Path.home() / ".cache" / "spot"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(CHUNK):
            h.update(chunk)
    return h.hexdigest()


def digest(path: Path) -> str:
    """SHA-256 of a file, remembered across runs while its size and modification time stay the same."""
    path = path.resolve()
    st = path.stat()
    key = f"{path}|{st.st_size}|{st.st_mtime_ns}"
    store = cache_dir() / "digests.json"
    known: dict[str, str] = {}
    try:
        known = json.loads(store.read_text())
    except (OSError, ValueError):
        known = {}
    if key in known:
        return known[key]
    value = _sha256(path)
    known[key] = value
    try:
        store.parent.mkdir(parents=True, exist_ok=True)
        _atomic_write(store, json.dumps(known, indent=0).encode())
    except OSError as e:
        log.debug("cannot write %s: %s", store, e)
    return value


def key_of(*parts: Any) -> str:
    """Stable SHA-256 of JSON-serialisable parts (paths are turned into strings)."""
    text = json.dumps(parts, sort_keys=True, default=str)
    return hashlib.sha256(text.encode()).hexdigest()


def _atomic_write(path: Path, data: bytes) -> None:
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.")
    try:
        with os.fdopen(fd, "wb") as fh:
            fh.write(data)
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


class Checkpoints:
    def __init__(self, outdir: Path, enabled: bool = True) -> None:
        self.folder = outdir / ".spot"
        self.enabled = enabled

    def stage(self, name: str, key: str, compute: Callable[[], T]) -> T:
        """Return the stored result of ``name`` if its key matches, otherwise compute and store it."""
        path = self.folder / f"{name}.pkl"
        if self.enabled and path.is_file():
            try:
                with open(path, "rb") as fh:
                    stored = pickle.load(fh)
                if stored.get("version") == CHECKPOINT_VERSION and stored.get("key") == key:
                    log.info("%s: reusing checkpoint", name)
                    return stored["value"]  # type: ignore[no-any-return]
            except (OSError, pickle.UnpicklingError, EOFError, AttributeError, KeyError) as e:
                log.debug("ignoring checkpoint %s: %s", path, e)
        value = compute()
        if self.enabled:
            self.folder.mkdir(parents=True, exist_ok=True)
            payload = {"version": CHECKPOINT_VERSION, "key": key, "value": value}
            _atomic_write(path, pickle.dumps(payload, protocol=pickle.HIGHEST_PROTOCOL))
        return value
