from __future__ import annotations

import logging
import shutil
import tempfile
from collections.abc import Iterable
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from .config import Gepard as GepardOptions
from .errors import ToolError
from .tools import Gepard, as_args, run
from .validate import Region

log = logging.getLogger("spot")

MAIN_CLASS = "org.gepard.client.cmdline.CommandLine"


class GepardPlotter:
    def __init__(self, install: Gepard, options: GepardOptions, threads: int) -> None:
        settings = options.model_dump()
        self.matrix = install.matrix(settings.pop("matrix"))
        if not self.matrix.is_file():
            raise ToolError(f"Gepard matrix not found: {self.matrix}")
        self.jar = install.jar
        self.options = as_args(settings, single_dash=True)
        self.threads = threads

    def plot(self, jobs: Iterable[tuple[Path, Region, Region]]) -> int:
        jobs = list(jobs)
        if not jobs:
            return 0
        log.info("drawing %d Gepard dotplots", len(jobs))
        with ThreadPoolExecutor(max_workers=self.threads) as pool:
            return sum(pool.map(lambda job: self._one(*job), jobs))

    def _one(self, out: Path, a: Region, b: Region) -> int:
        # Gepard splits -seq on spaces and appends the pair number to -outfile, so it works in a temp dir
        with tempfile.TemporaryDirectory(prefix="spot_gepard_") as tmp:
            fa, fb, png = Path(tmp, "a.fa"), Path(tmp, "b.fa"), Path(tmp, "dotplot.png")
            for path, reg in ((fa, a), (fb, b)):
                path.write_text(f">{reg.name} {reg.chrom}:{reg.start}-{reg.end}\n{reg.seq}\n")
            run(
                [
                    "java",
                    "-Djava.awt.headless=true",
                    "-cp",
                    self.jar,
                    MAIN_CLASS,
                    "-seq",
                    fa,
                    fb,
                    "-matrix",
                    self.matrix,
                    "-outfile",
                    png,
                    "-silent",
                    *self.options,
                ]
            )
            written = next((p for p in (Path(tmp, "dotplot0.png"), png) if p.is_file()), None)
            if written is None:
                raise ToolError(f"Gepard did not write a dotplot for {a.name} vs {b.name}")
            out.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(written, out)
        return 1
