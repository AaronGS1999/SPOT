from __future__ import annotations

import csv
import shutil
from collections import Counter
from pathlib import Path

from .classify import SYNTENIC, RunResult
from .fasta import write_fasta
from .locate import Location
from .orthofinder import OrthoFinderRun
from .positional import SYNTENIC_ORTHOLOG

SEQUENCE_COLUMNS = [
    "id",
    "species",
    "orthogroup",
    "zone",
    "og_status",
    "og_direction",
    "class",
    "evidence",
    "block_type",
    "syri_blocks",
    "chr",
    "start",
    "end",
    "strand",
    "location",
    "identity",
    "anchor_gap_left",
    "anchor_gap_right",
    "other_chr",
    "window_start",
    "window_end",
    "window_aligned_frac",
    "orthologs_in_window",
    "others_in_window",
    "unconfirmed",
    "flank_left_bp",
    "flank_right_bp",
    "flank_back_left_bp",
    "flank_back_right_bp",
    "flank_orientation",
    "description",
]


def pair_directions(result: RunResult, og: str) -> dict[tuple[str, str], str]:
    found: dict[tuple[str, str], list[str]] = {}
    ref, qry = result.reference, result.query
    for sp, other in ((ref, qry), (qry, ref)):
        for g in result.run.orthogroups[og].members[sp]:
            for o in result.tests[g].orthologs:
                key = (g, o) if sp == ref else (o, g)
                found.setdefault(key, []).append(f"{sp}->{other}")
    return {k: "+".join(v) for k, v in found.items()}


def export_run(
    result: RunResult, locations: dict[str, Location], sequences: dict[str, str], outdir: Path
) -> Path:
    run = result.run
    folder = outdir / run.name
    if folder.exists():
        shutil.rmtree(folder)
    folder.mkdir(parents=True)
    ref, qry = result.reference, result.query

    with open(folder / "sequences.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(SEQUENCE_COLUMNS)
        for sp in (ref, qry):
            for g in (g for g in run.genes.values() if g.species == sp):
                t, loc = result.tests[g.id], locations[g.id]
                zone = result.zone(g.id)
                w.writerow(
                    [
                        g.id,
                        sp,
                        g.orthogroup if g.assigned else "",
                        zone,
                        result.status.get(g.orthogroup, "") if zone == "shared" else "",
                        result.direction.get(g.orthogroup, "") if zone == "shared" else "",
                        t.klass,
                        t.evidence,
                        t.block_type,
                        t.blocks,
                        loc.chrom,
                        loc.start or "",
                        loc.end or "",
                        loc.strand,
                        loc.status,
                        "" if loc.identity is None else f"{loc.identity:.4f}",
                        _blank(t.gap_left),
                        _blank(t.gap_right),
                        t.other_chr,
                        _blank(t.window_start),
                        _blank(t.window_end),
                        "" if t.window_aligned is None else f"{t.window_aligned:.3f}",
                        ",".join(t.orthologs),
                        ",".join(t.others),
                        ",".join(t.unconfirmed),
                        _blank(t.flank_left),
                        _blank(t.flank_right),
                        _blank(t.flank_back_left),
                        _blank(t.flank_back_right),
                        {None: "", True: "inverted", False: "direct"}[t.inverted],
                        g.description,
                    ]
                )

    with open(folder / "orthogroups.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["orthogroup", "zone", "status", "direction", f"n_{ref}", f"n_{qry}", "tree", "pairs"])
        for og in run.orthogroups.values():
            pairs = ";".join(f"{a}|{b}" for a, b in result.pairs.get(og.id, []))
            w.writerow(
                [
                    og.id,
                    result.og_zone(og.id),
                    result.status.get(og.id, ""),
                    result.direction.get(og.id, ""),
                    len(og.members[ref]),
                    len(og.members[qry]),
                    "yes" if og.id in run.trees else "no",
                    pairs,
                ]
            )

    with open(folder / "pairs.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(
            [
                "orthogroup",
                f"{ref}_id",
                f"{ref}_chr",
                f"{ref}_start",
                f"{ref}_end",
                f"{ref}_strand",
                f"{qry}_id",
                f"{qry}_chr",
                f"{qry}_start",
                f"{qry}_end",
                f"{qry}_strand",
                "found_from",
            ]
        )
        for og_id, og_pairs in result.pairs.items():
            directions = pair_directions(result, og_id)
            for a, b in og_pairs:
                la, lb = locations[a], locations[b]
                w.writerow(
                    [
                        og_id,
                        a,
                        la.chrom,
                        la.start,
                        la.end,
                        la.strand,
                        b,
                        lb.chrom,
                        lb.start,
                        lb.end,
                        lb.strand,
                        directions[(a, b)],
                    ]
                )

    for og in run.orthogroups.values():
        target = folder / result.folder(og.id) / og.id
        target.mkdir(parents=True)
        members = [g for sp in (ref, qry) for g in og.members[sp]]
        write_fasta(target / f"{og.id}.fa", ((_header(run, g), sequences[g]) for g in members))
        paired = sorted({g for pair in result.pairs.get(og.id, []) for g in pair}, key=members.index)
        if paired:
            write_fasta(target / f"{og.id}.pairs.fa", ((_header(run, g), sequences[g]) for g in paired))
        tree = run.tree(og.id)
        if tree is not None:
            shutil.copyfile(tree, target / f"{og.id}_tree.txt")

    for sp in (ref, qry):
        genes = run.unassigned(sp)
        if genes:
            (folder / "unassigned").mkdir(exist_ok=True)
            write_fasta(
                folder / "unassigned" / f"{sp}.fa", ((_header(run, g.id), sequences[g.id]) for g in genes)
            )

    (folder / "summary.txt").write_text(summary(result))
    return folder


def summary(result: RunResult) -> str:
    run, ref, qry = result.run, result.reference, result.query
    zones = Counter(result.og_zone(og) for og in run.orthogroups)
    status = Counter(result.status.values())
    syn = [og for og, s in result.status.items() if s == SYNTENIC]
    directions = Counter(result.direction[og] for og in syn).most_common()
    lines = [
        f"{run.name}: {len(run.orthogroups)} orthogroups, {len(run.genes)} sequences",
        f"{ref}-only: {zones[f'{ref}_only']} orthogroups + {len(run.unassigned(ref))} unassigned",
        f"{qry}-only: {zones[f'{qry}_only']} orthogroups + {len(run.unassigned(qry))} unassigned",
        f"shared: {zones['shared']} orthogroups; syntenic orthologs {status[SYNTENIC]}, "
        f"evaluable not syntenic {status['evaluable_not_syntenic']}, not evaluable {status['not_evaluable']}",
        "syntenic orthogroups by direction: " + (", ".join(f"{d} {n}" for d, n in directions) or "-"),
        f"confirmed pairs: {sum(len(p) for p in result.pairs.values())}",
        "",
        "sequences by zone and class:",
    ]
    table = Counter((result.zone(g), t.klass) for g, t in result.tests.items())
    for zone in sorted({z for z, _ in table}):
        parts = ", ".join(f"{k}={n}" for (z, k), n in sorted(table.items()) if z == zone)
        lines.append(f"  {zone:24s} {parts}")
    evidence = Counter(
        (t.species, t.evidence, t.block_type) for t in result.tests.values() if t.klass == SYNTENIC_ORTHOLOG
    )
    lines += ["", "syntenic_ortholog tests by evidence:"]
    lines += [f"  {sp} {ev} {bt}: {n}" for (sp, ev, bt), n in sorted(evidence.items())]
    return "\n".join(lines) + "\n"


def _header(run: OrthoFinderRun, gene: str) -> str:
    description = run.genes[gene].description
    return f"{gene} {description}" if description else gene


def _blank(value: int | None) -> object:
    return "" if value is None else value
