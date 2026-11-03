from __future__ import annotations

import logging
import shutil
import tempfile
from collections.abc import Iterable
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from .config import Gepard as GepardOptions
from .errors import ToolError
from .progress import Progress, silent
from .tools import Gepard, as_args, java_executable, run
from .validate import Region

log = logging.getLogger("spot")

MAIN_CLASS = "org.gepard.client.cmdline.CommandLine"


class GepardPlotter:
    """Draws Gepard dotplots in a few long-lived JVMs.

    Starting a JVM costs about a second, more than a dotplot of two 20 kb regions, so the jobs are split
    into one chunk per worker and each chunk is drawn by a single Gepard call (``-seq a1 b1 a2 b2 ...``).
    Each JVM is told it has one processor so that ``threads`` JVMs fill ``threads`` cores without the
    JIT and garbage-collector threads competing for them.
    """

    def __init__(
        self,
        install: Gepard,
        options: GepardOptions,
        threads: int,
        progress: Progress = silent,
        memory_mb: int = 1024,
    ) -> None:
        settings = options.model_dump()
        self.matrix = install.matrix(settings.pop("matrix"))
        if not self.matrix.is_file():
            raise ToolError(f"Gepard matrix not found: {self.matrix}")
        self.jar = install.jar
        self.java = java_executable()
        self.options = as_args(settings, single_dash=True)
        self.threads = max(1, threads)
        self.progress = progress
        self.memory_mb = memory_mb

    def plot(self, jobs: Iterable[tuple[Path, Region, Region]]) -> int:
        jobs = list(jobs)
        if not jobs:
            return 0
        n_chunks = min(len(jobs), self.threads)
        chunks = [jobs[i::n_chunks] for i in range(n_chunks)]
        log.info("drawing %d Gepard dotplots in %d JVMs", len(jobs), n_chunks)
        done = 0
        self.progress("dotplots", 0, len(jobs))
        with ThreadPoolExecutor(max_workers=n_chunks) as pool:
            for n in pool.map(self._chunk, chunks):
                done += n
                self.progress("dotplots", done, len(jobs))
        return done

    def _chunk(self, jobs: list[tuple[Path, Region, Region]]) -> int:
        # Gepard splits -seq on spaces and names pair k <outfile stem><2k>.png, so it works in a temp dir
        with tempfile.TemporaryDirectory(prefix="spot_gepard_") as tmp:
            seqs: list[Path] = []
            for k, (_, a, b) in enumerate(jobs):
                for tag, reg in (("a", a), ("b", b)):
                    path = Path(tmp, f"{tag}{k}.fa")
                    path.write_text(f">{reg.name} {reg.chrom}:{reg.start}-{reg.end}\n{reg.seq}\n")
                    seqs.append(path)
            run(
                [
                    self.java,
                    "-Djava.awt.headless=true",
                    "-XX:ActiveProcessorCount=1",
                    f"-Xmx{self.memory_mb}m",
                    "-cp",
                    self.jar,
                    MAIN_CLASS,
                    "-seq",
                    *seqs,
                    "-matrix",
                    self.matrix,
                    "-outfile",
                    Path(tmp, "dp.png"),
                    "-silent",
                    *self.options,
                ]
            )
            for k, (out, a, b) in enumerate(jobs):
                written = Path(tmp, f"dp{2 * k}.png")
                if not written.is_file():
                    raise ToolError(f"Gepard did not write a dotplot for {a.name} vs {b.name}")
                out.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(written, out)
        return len(jobs)
