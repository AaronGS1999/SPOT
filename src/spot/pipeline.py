from __future__ import annotations

import logging
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

from . import __version__, orthofinder, plots, positional, syri, tools
from .anchors import AnchorIndex
from .cache import Checkpoints, digest, key_of
from .classify import SYNTENIC, RunResult, classify
from .config import Config, dump
from .dotplot import GepardPlotter
from .errors import ConfigError, ToolError
from .export import export_run
from .locate import Genome, Location, locate
from .log import add_file
from .orthofinder import OrthoFinderRun
from .positional import FLANKS_NOT_CONFIRMED, SYNTENIC_ORTHOLOG, GeneIndex, Test
from .progress import Progress, Reporter, silent
from .syri import AlignedIndex, Alignment
from .util import natural_key
from .validate import Region, Validator, region

log = logging.getLogger("spot")


@dataclass
class Analysis:
    config: Config
    runs: list[OrthoFinderRun]
    sequences: dict[str, dict[str, str]]
    genomes: dict[str, Genome]
    locations: dict[str, Location]
    syri_out: Path
    alignments: list[Alignment]
    blocks: list[Alignment]
    homologous: set[tuple[str, str]]
    results: list[RunResult]


def check_inputs(config: Config) -> None:
    problems = [f"species {k}: genome is required" for k, sp in config.species.items() if sp.genome is None]
    files: list[tuple[str, Path | None]] = [("syri_out", config.syri_out)]
    for k, sp in config.species.items():
        files += [(f"species {k}: genome", sp.genome), (f"species {k}: centromeres", sp.centromeres)]
        files += [(f"species {k}: lift", p) for p in sp.lift]
        if isinstance(sp.coords, Path):
            files.append((f"species {k}: coords", sp.coords))
    problems += [f"{what}: not found: {path}" for what, path in files if path and not path.is_file()]
    if problems:
        raise ConfigError("\n".join(problems))


def homology(config: Config, alignments: list[Alignment]) -> set[tuple[str, str]]:
    """Homologous chromosome pairs: given in the configuration, or inferred from the SYN alignments."""
    if config.anchors.homology != "auto":
        return {(r, q) for r, q in config.anchors.homology}
    return syri.homologous_pairs(alignments, config.anchors.min_syn_frac)


def analyse(config: Config, checkpoints: Checkpoints | None = None, progress: Progress = silent) -> Analysis:
    check_inputs(config)
    cp = checkpoints or Checkpoints(config.outdir, enabled=False)
    ref, qry = config.reference, config.query
    other = {ref: qry, qry: ref}

    runs = orthofinder.load_all(config)
    sequences = {run.name: orthofinder.load_sequences(run) for run in runs}
    genomes = {k: Genome(sp.genome) for k, sp in config.species.items() if sp.genome is not None}
    progress("locate", 0, 1)
    with ThreadPoolExecutor(max_workers=len(genomes) or 1) as pool:  # hashing is I/O bound
        digests = dict(zip(genomes, pool.map(lambda g: digest(g.path), genomes.values()), strict=True))
    species_key = key_of(
        {k: sp.model_dump(mode="json") for k, sp in config.species.items()},
        digests,
        {run.name: key_of(sorted(sequences[run.name].items())) for run in runs},
    )
    locations = cp.stage("locate", species_key, lambda: locate(runs, config.species, sequences, genomes))
    progress("locate", 1, 1)
    progress("align", 0, 1)

    syri_out = config.syri_out or syri.align(config, config.outdir / "genome")
    alignments, blocks = syri.read(syri_out)
    progress("align", 1, 1)
    homologous = homology(config, alignments)
    order = sorted(homologous, key=lambda p: (natural_key(p[0]), natural_key(p[1])))
    log.info("homologous chromosomes: %s", ", ".join(f"{r}-{q}" for r, q in order))
    kinds = config.anchors.types
    pairs = homologous if config.anchors.same_chromosome else None
    anchors = {
        ref: AnchorIndex.build(alignments, "reference", kinds, pairs),
        qry: AnchorIndex.build(alignments, "query", kinds, pairs),
    }
    aligned = {ref: AlignedIndex(alignments, "reference"), qry: AlignedIndex(alignments, "query")}
    log.info("anchors (%s): %d alignments", ", ".join(kinds), len(anchors[ref]))
    genes = {sp: GeneIndex(loc for loc in locations.values() if loc.species == sp) for sp in (ref, qry)}

    tests: dict[str, dict[str, Test]] = {}
    jobs: list[tuple[Location, Location]] = []
    progress("test", 0, len(runs))
    for n_run, run in enumerate(runs, start=1):
        tests[run.name] = {}
        for g in run.genes.values():
            same = set(run.orthogroups[g.orthogroup].members[other[g.species]]) if g.assigned else set()
            target = other[g.species]
            t = positional.test(
                locations[g.id], same, anchors[g.species], genes[target], config.anchors, aligned[target]
            )
            tests[run.name][g.id] = t
            if t.klass == SYNTENIC_ORTHOLOG:
                jobs.extend((locations[g.id], locations[o]) for o in t.orthologs)
        progress("test", n_run, len(runs))

    if config.validation.enabled:
        validator = Validator(config.validation, config.blastn, genomes, config.n_threads, progress)
        job_key = key_of(
            species_key,
            config.validation.model_dump(mode="json"),
            config.blastn.model_dump(mode="json"),
            sorted({(a.gene, b.gene) for a, b in jobs}),
        )
        checks = cp.stage("validation", job_key, lambda: validator.check(jobs))
        for run_tests in tests.values():
            for t in run_tests.values():
                if t.klass != SYNTENIC_ORTHOLOG:
                    continue
                found = {o: checks[(t.gene, o)] for o in t.orthologs}
                best = max(found.values(), key=lambda c: (c.ok, min(c.left, c.right)))
                t.flank_left, t.flank_right = best.left, best.right
                t.flank_back_left, t.flank_back_right = best.back_left, best.back_right
                t.inverted = best.inverted
                t.unconfirmed = [o for o, c in found.items() if not c.ok]
                t.orthologs = [o for o, c in found.items() if c.ok]
                if not t.orthologs:
                    t.klass = FLANKS_NOT_CONFIRMED

    results = [classify(run, tests[run.name], ref, qry) for run in runs]
    return Analysis(
        config, runs, sequences, genomes, locations, syri_out, alignments, blocks, homologous, results
    )


def run(config: Config, progress: Progress | None = None) -> Analysis:
    progress = progress or Reporter()
    out = config.outdir
    out.mkdir(parents=True, exist_ok=True)
    handler = add_file(out / "spot.log")
    try:
        log.info("spot %s, profile %s, output %s", __version__, config.profile, out)
        check_inputs(config)
        (out / "run.yaml").write_text(dump(config))
        report = tools.check(config)
        _write_versions(out / "versions.tsv", report)
        broken = [s.name for s in report if s.required and not s.ok]
        if broken:
            raise ToolError(f"missing or broken: {', '.join(broken)} (see spot check)")

        log.info("threads: %d", config.n_threads)
        analysis = analyse(config, Checkpoints(out, enabled=config.resume), progress)
        _write_homology(out / "homology.tsv", analysis)
        species = [plots.SpeciesStyle(k, config.species[k].label) for k in (config.reference, config.query)]
        styles = _run_styles(config, analysis.runs)
        palette = plots.Palette.from_config(config.figures)
        lengths = {k: g.lengths for k, g in analysis.genomes.items()}
        dotplots: list[tuple[Path, Region, Region]] = []
        progress("export", 0, len(analysis.results))

        def export(result: RunResult) -> Path:
            return export_run(result, analysis.locations, analysis.sequences[result.run.name], out)

        # writing thousands of small files is I/O bound (slow on network drives): one thread per run
        with ThreadPoolExecutor(max_workers=max(1, min(len(analysis.results), config.n_threads))) as pool:
            folders = list(pool.map(export, analysis.results))
        progress("export", len(analysis.results), len(analysis.results))
        for result, style, folder in zip(analysis.results, styles, folders, strict=True):
            plots.venn(result, species, style, folder / "venn", palette)
            if config.dotplots.pairs:
                dotplots += _dotplot_jobs(result, analysis, folder)
        if dotplots:
            gepard = tools.find_gepard()
            if gepard is None:
                raise ToolError("Gepard not found")
            GepardPlotter(gepard, config.gepard, config.n_threads, progress).plot(dotplots)
        progress("plots", 0, 1)

        plots.venn_panel(analysis.results, species, styles, out / "venn_panel", palette)
        centromeres = {k: plots.read_bed(sp.centromeres) for k, sp in config.species.items()}
        plots.synteny(
            analysis.results,
            styles,
            species,
            analysis.locations,
            analysis.blocks,
            analysis.homologous,
            lengths,
            centromeres,
            out / "synteny",
            palette,
        )
        if config.dotplots.genome:
            links = [
                (style, _links(result, analysis.locations))
                for result, style in zip(analysis.results, styles, strict=True)
            ]
            plots.genome_dotplot(
                analysis.alignments,
                species,
                lengths,
                analysis.homologous,
                links,
                out / "genome" / "dotplot",
                palette,
            )
        progress("plots", 1, 1)
        _report(analysis)
        return analysis
    finally:
        log.removeHandler(handler)
        handler.close()


def _run_styles(config: Config, runs: list[OrthoFinderRun]) -> list[plots.RunStyle]:
    styles = []
    for i, run in enumerate(runs):
        settings = config.orthofinder.runs.get(run.name)
        styles.append(
            plots.RunStyle(
                label=(settings.label if settings and settings.label else run.name),
                color=(
                    settings.color
                    if settings and settings.color
                    else plots.DEFAULT_COLORS[i % len(plots.DEFAULT_COLORS)]
                ),
                unit=(settings.unit if settings and settings.unit else "sequences"),
            )
        )
    return styles


def _dotplot_jobs(result: RunResult, analysis: Analysis, folder: Path) -> list[tuple[Path, Region, Region]]:
    window = analysis.config.validation.window
    jobs: list[tuple[Path, Region, Region]] = []
    for og, pairs in result.pairs.items():
        for a, b in pairs:
            la, lb = analysis.locations[a], analysis.locations[b]
            jobs.append(
                (
                    folder / SYNTENIC / og / "dotplots" / f"{a}__{b}.png",
                    region(analysis.genomes[la.species], la, window),
                    region(analysis.genomes[lb.species], lb, window),
                )
            )
    return jobs


def _links(result: RunResult, locations: dict[str, Location]) -> list[tuple[str, float, str, float]]:
    links = []
    for pairs in result.pairs.values():
        for a, b in pairs:
            la, lb = locations[a], locations[b]
            links.append((la.chrom, (la.start + la.end) / 2, lb.chrom, (lb.start + lb.end) / 2))
    return links


def _write_homology(path: Path, analysis: Analysis) -> None:
    rows = syri.homology_table(analysis.alignments, analysis.homologous)
    lines = [f"{r}\t{q}\t{n}\t{frac:.4f}\n" for r, q, n, frac in rows]
    path.write_text("reference_chr\tquery_chr\tsyn_bp\tfraction_of_reference_syn\n" + "".join(lines))


def _write_versions(path: Path, report: list[tools.Status]) -> None:
    rows = [("spot", __version__, ""), ("python", sys.version.split()[0], sys.executable)]
    rows += [(s.name, s.version or "", s.path or "") for s in report]
    rows += [(name, version or "", "") for name, version in tools.python_packages()]
    path.write_text("tool\tversion\tpath\n" + "".join("\t".join(r) + "\n" for r in rows))


def _report(analysis: Analysis) -> None:
    header = ["run", "orthogroups", "shared", SYNTENIC, "evaluable_not_syntenic", "not_evaluable", "pairs"]
    rows = [header]
    for result in analysis.results:
        status = Counter(result.status.values())
        rows.append(
            [
                result.run.name,
                str(len(result.run.orthogroups)),
                str(len(result.status)),
                str(status[SYNTENIC]),
                str(status["evaluable_not_syntenic"]),
                str(status["not_evaluable"]),
                str(sum(len(p) for p in result.pairs.values())),
            ]
        )
    widths = [max(len(r[i]) for r in rows) for i in range(len(header))]
    for row in rows:
        log.info(
            "%s",
            "  ".join(
                c.ljust(w) if i == 0 else c.rjust(w) for i, (c, w) in enumerate(zip(row, widths, strict=True))
            ),
        )
