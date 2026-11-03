from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import matplotlib
import numpy as np
from matplotlib.axes import Axes
from matplotlib.collections import LineCollection
from matplotlib.colors import to_hex, to_rgb
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Patch, PathPatch, Polygon, Rectangle
from matplotlib.path import Path as MPath

from .classify import SYNTENIC, RunResult
from .locate import Location
from .positional import SYNTENIC_ORTHOLOG
from .syri import Alignment
from .util import natural_key

matplotlib.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    }
)

DEFAULT_COLORS = ["#C1121F", "#1B7F3B", "#7CC78E", "#2B59C3", "#E09F3E", "#7B2CBF"]
BLOCK_COLORS = {
    "SYN": "#DCDCDC",
    "INV": "#F4B97A",
    "TRANS": "#B9B3E0",
    "INVTR": "#B9B3E0",
    "DUP": "#9CCBE8",
    "INVDP": "#9CCBE8",
}
BLOCK_LABELS = {"SYN": "Syntenic", "INV": "Inversion", "TRANS": "Translocation", "DUP": "Duplication"}
DOT_COLORS = {
    "SYN": "#4D4D4D",
    "INV": "#E07B00",
    "TRANS": "#6A5ACD",
    "INVTR": "#6A5ACD",
    "DUP": "#2E86C1",
    "INVDP": "#2E86C1",
}
REF_GREY, QRY_GREY, CREAM = "#3A3A3A", "#8C8C8C", "#FCEFC4"


@dataclass(frozen=True)
class RunStyle:
    label: str
    color: str
    unit: str


@dataclass(frozen=True)
class SpeciesStyle:
    key: str
    label: str


def _mix(color: str, other: str, t: float) -> str:
    a, b = to_rgb(color), to_rgb(other)
    mixed = tuple(float(x * (1 - t) + y * t) for x, y in zip(a, b, strict=True))
    return to_hex((mixed[0], mixed[1], mixed[2]))


def _save(fig: Figure, prefix: Path, formats: Sequence[str], dpi: int = 400) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for ext in formats:
        fig.savefig(f"{prefix}.{ext}", dpi=dpi, bbox_inches="tight", facecolor="white")


def _venn_counts(result: RunResult) -> tuple[dict[str, list[int]], Counter[str], Counter[str]]:
    ref, qry = result.reference, result.query
    zones: dict[str, list[int]] = defaultdict(lambda: [0, 0, 0])
    for og in result.run.orthogroups.values():
        zone = result.og_zone(og.id)
        if zone == "shared":
            zone = "syn" if result.status[og.id] == SYNTENIC else "non"
        zones[zone][0] += 1
        zones[zone][1] += len(og.members[ref])
        zones[zone][2] += len(og.members[qry])
    unassigned = Counter(g.species for g in result.run.genes.values() if not g.assigned)
    total = Counter(g.species for g in result.run.genes.values())
    return zones, unassigned, total


def _venn(ax: Axes, result: RunResult, species: Sequence[SpeciesStyle], style: RunStyle) -> None:
    ref, qry = species
    pal = {
        ref.key: _mix(style.color, "#000000", 0.1),
        qry.key: _mix(style.color, "#FFFFFF", 0.35),
        "syn": _mix(style.color, "#FFFFFF", 0.45),
        "non": _mix(style.color, "#FFFFFF", 0.82),
        "txt": _mix(style.color, "#000000", 0.6),
    }
    zones, unassigned, total = _venn_counts(result)
    r0, dx = 1.5, 0.95
    yy = np.sqrt(r0**2 - dx**2)

    def lens(sign: int) -> np.ndarray:
        t = np.linspace(0, np.arcsin(yy / r0), 200)
        outline: np.ndarray = np.r_[
            np.c_[-dx + r0 * np.cos(t), sign * r0 * np.sin(t)],
            np.c_[dx - r0 * np.cos(t[::-1]), sign * r0 * np.sin(t[::-1])],
        ]
        return outline

    def block(x: float, y: float, big: str, small: str, color: str = "#333333", size: float = 10.5) -> None:
        ax.text(
            x,
            y + 0.13,
            big,
            ha="center",
            va="center",
            fontsize=size,
            color=color,
            fontweight="bold",
            zorder=5,
        )
        ax.text(
            x,
            y - 0.17,
            small,
            ha="center",
            va="center",
            fontsize=7.6,
            color=color,
            linespacing=1.15,
            zorder=5,
        )

    ax.set_aspect("equal")
    ax.axis("off")
    for cx, sp in ((-dx, ref.key), (dx, qry.key)):
        ax.add_patch(Circle((cx, 0), r0, fc=pal[sp], alpha=0.13, ec="none", zorder=1))
    ax.add_patch(Polygon(lens(1), closed=True, fc=pal["syn"], ec="none", zorder=2))
    ax.add_patch(Polygon(lens(-1), closed=True, fc=pal["non"], ec="none", zorder=2))
    ax.plot([-(r0 - dx), r0 - dx], [0, 0], color="white", lw=1.4, zorder=3)
    for cx, sp in ((-dx, ref.key), (dx, qry.key)):
        ax.add_patch(Circle((cx, 0), r0, fc="none", ec=pal[sp], lw=1.6, zorder=4))
    xo = dx + 0.62
    only_ref, only_qry = zones[f"{ref.key}_only"], zones[f"{qry.key}_only"]
    block(-xo, 0, f"{only_ref[0]} OGs", f"{only_ref[1]:,} {style.unit}\n+ {unassigned[ref.key]:,} unassigned")
    block(xo, 0, f"{only_qry[0]} OGs", f"{only_qry[2]:,} {style.unit}\n+ {unassigned[qry.key]:,} unassigned")
    syn, non = zones["syn"], zones["non"]
    block(0, yy * 0.4, f"{syn[0]} OGs", f"{syn[1]} {ref.key}\n{syn[2]} {qry.key}", color=pal["txt"], size=9.3)
    block(0, -yy * 0.4, f"{non[0]} OGs", f"{non[1]} {ref.key}\n{non[2]} {qry.key}", size=9.3)
    for cx, spec in ((-dx, ref), (dx, qry)):
        label_color = _mix(pal[spec.key], "#000000", 0.3)
        ax.text(
            cx,
            r0 + 0.24,
            spec.label,
            ha="center",
            va="bottom",
            fontsize=10.5,
            color=label_color,
            style="italic",
            fontweight="bold",
        )
        ax.text(
            cx, r0 + 0.08, f"n = {total[spec.key]:,}", ha="center", va="bottom", fontsize=8, color=label_color
        )
    ax.text(0, r0 + 0.72, style.label, ha="center", va="bottom", fontsize=11.5, fontweight="bold")
    handles = [
        Patch(fc=pal["syn"], ec="#BBBBBB", lw=0.5, label="Syntenic orthologs"),
        Patch(fc=pal["non"], ec="#BBBBBB", lw=0.5, label="Non-syntenic homologs"),
    ]
    ax.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.02),
        ncol=2,
        frameon=False,
        fontsize=8,
        handlelength=1.1,
        handleheight=1.1,
        columnspacing=1.4,
    )
    ax.set_xlim(-dx - r0 - 0.05, dx + r0 + 0.05)
    ax.set_ylim(-r0 - 0.05, r0 + 1.05)


def venn(result: RunResult, species: Sequence[SpeciesStyle], style: RunStyle, prefix: Path) -> None:
    fig = Figure(figsize=(5.6, 4.6))
    ax = fig.subplots()
    _venn(ax, result, species, style)
    _save(fig, prefix, ("pdf", "svg", "png"), dpi=600)


def venn_panel(
    results: Sequence[RunResult], species: Sequence[SpeciesStyle], styles: Sequence[RunStyle], prefix: Path
) -> None:
    fig = Figure(figsize=(5.0 * len(results), 4.6))
    axs = fig.subplots(1, len(results), squeeze=False)
    for ax, result, style in zip(axs[0], results, styles, strict=True):
        _venn(ax, result, species, style)
    fig.tight_layout(w_pad=1.5)
    _save(fig, prefix, ("pdf", "svg", "png"), dpi=600)


def _ribbon(ax: Axes, a: float, b: float, c: float, d: float, y1: float, y0: float, color: str) -> None:
    ym = (y1 + y0) / 2
    verts = [(a, y1), (b, y1), (b, ym), (d, ym), (d, y0), (c, y0), (c, ym), (a, ym), (a, y1), (a, y1)]
    codes = [
        MPath.MOVETO,
        MPath.LINETO,
        MPath.CURVE4,
        MPath.CURVE4,
        MPath.CURVE4,
        MPath.LINETO,
        MPath.CURVE4,
        MPath.CURVE4,
        MPath.CURVE4,
        MPath.CLOSEPOLY,
    ]
    ax.add_patch(PathPatch(MPath(verts, codes), fc=color, ec="none", alpha=0.88, zorder=1))


def _curve(ax: Axes, x1: float, y1: float, x0: float, y0: float, color: str) -> None:
    ym = (y1 + y0) / 2
    path = MPath([(x1, y1), (x1, ym), (x0, ym), (x0, y0)], [MPath.MOVETO] + [MPath.CURVE4] * 3)
    ax.add_patch(PathPatch(path, fc="none", ec=color, lw=0.75, alpha=0.95, zorder=3))


def read_bed(path: Path | None) -> dict[str, list[tuple[int, int]]]:
    regions: dict[str, list[tuple[int, int]]] = defaultdict(list)
    if path is None:
        return regions
    for line in path.read_text().splitlines():
        f = line.split("\t")
        if len(f) >= 3 and not line.startswith(("#", "track", "browser")):
            regions[f[0]].append((int(f[1]), int(f[2])))
    return regions


def synteny(
    results: Sequence[RunResult],
    styles: Sequence[RunStyle],
    species: Sequence[SpeciesStyle],
    locations: Mapping[str, Location],
    blocks: Sequence[Alignment],
    homologous: set[tuple[str, str]],
    lengths: Mapping[str, Mapping[str, int]],
    centromeres: Mapping[str, Mapping[str, list[tuple[int, int]]]],
    prefix: Path,
    bin_size: int = 250_000,
) -> None:
    ref, qry = species
    chroms = sorted(
        (p for p in homologous if p[0] in lengths[ref.key] and p[1] in lengths[qry.key]),
        key=lambda p: natural_key(p[0]),
    )
    if not chroms:
        return
    tracks = []
    for result, style in zip(results, styles, strict=True):
        density: dict[str, dict[str, np.ndarray]] = {ref.key: {}, qry.key: {}}
        positions: dict[str, dict[str, list[float]]] = {
            ref.key: defaultdict(list),
            qry.key: defaultdict(list),
        }
        dots: dict[str, dict[str, list[float]]] = {ref.key: defaultdict(list), qry.key: defaultdict(list)}
        for gene, t in result.tests.items():
            loc = locations[gene]
            if loc.located:
                mid = (loc.start + loc.end) / 2
                positions[loc.species][loc.chrom].append(mid)
                if t.klass == SYNTENIC_ORTHOLOG:
                    dots[loc.species][loc.chrom].append(mid)
        for sp in (ref.key, qry.key):
            for chrom, values in positions[sp].items():
                if chrom in lengths[sp]:
                    edges = np.arange(0, lengths[sp][chrom] + bin_size, bin_size)
                    density[sp][chrom] = np.histogram(values, bins=edges)[0]
        top = max([h.max() for sp in density for h in density[sp].values()] + [1])
        links = []
        for pairs in result.pairs.values():
            for a, b in pairs:
                la, lb = locations[a], locations[b]
                links.append((la.chrom, (la.start + la.end) / 2, lb.chrom, (lb.start + lb.end) / 2))
        tracks.append((style, density, top, dots, links))

    y_ref, y_qry, track_h, track_gap = 1.0, 0.0, 0.11, 0.03
    upper = y_ref + 0.1 + len(tracks) * (track_h + track_gap)
    lower = y_qry - 0.1 - len(tracks) * (track_h + track_gap)
    height = 1.9 * len(chroms) + 1.1
    xmax = max(max(lengths[ref.key][r], lengths[qry.key][q]) for r, q in chroms) / 1e6
    step, pad = _grid_step(xmax), 0.004 * xmax
    fig = Figure(figsize=(10, height))
    axs = fig.subplots(len(chroms), 1, squeeze=False)
    by_pair: dict[tuple[str, str], list[Alignment]] = defaultdict(list)
    for blk in blocks:
        by_pair[(blk.ref_chr, blk.qry_chr)].append(blk)
    for ax, (rc, qc) in zip(axs[:, 0], chroms, strict=True):
        for sp_key, y, chrom in ((ref.key, y_ref, rc), (qry.key, y_qry, qc)):
            for c0, c1 in centromeres.get(sp_key, {}).get(chrom, []):
                ax.add_patch(
                    Rectangle((c0 / 1e6, y - 0.09), (c1 - c0) / 1e6, 0.18, fc=CREAM, ec="none", zorder=0.5)
                )
        for blk in by_pair[(rc, qc)]:
            _ribbon(
                ax,
                blk.ref_start / 1e6,
                blk.ref_end / 1e6,
                blk.qry_start / 1e6,
                blk.qry_end / 1e6,
                y_ref,
                y_qry,
                BLOCK_COLORS[blk.kind],
            )
        for style, _, _, _, links in tracks:
            for a_chr, a_mid, b_chr, b_mid in links:
                if a_chr == rc and b_chr == qc:
                    _curve(ax, a_mid / 1e6, y_ref, b_mid / 1e6, y_qry, style.color)
        ax.plot(
            [0, lengths[ref.key][rc] / 1e6],
            [y_ref, y_ref],
            color=REF_GREY,
            lw=2.2,
            solid_capstyle="round",
            zorder=4,
        )
        ax.plot(
            [0, lengths[qry.key][qc] / 1e6],
            [y_qry, y_qry],
            color=QRY_GREY,
            lw=2.2,
            solid_capstyle="round",
            zorder=4,
        )
        for i, (style, density, top, dots, _) in enumerate(tracks):
            for sp_key, y, chrom in ((ref.key, y_ref, rc), (qry.key, y_qry, qc)):
                xs = np.array(dots[sp_key].get(chrom, [])) / 1e6
                ax.scatter(
                    xs, np.full(len(xs), y), s=9, c=style.color, edgecolors="white", linewidths=0.3, zorder=5
                )
            for sp_key, chrom, base, sign in (
                (ref.key, rc, y_ref + 0.12 + i * (track_h + track_gap), 1),
                (qry.key, qc, y_qry - 0.12 - i * (track_h + track_gap), -1),
            ):
                n_bins = int(np.ceil(lengths[sp_key][chrom] / bin_size))
                xb = (np.arange(n_bins) * bin_size + bin_size / 2) / 1e6
                h = density[sp_key].get(chrom, np.zeros(n_bins))[:n_bins]
                ax.fill_between(
                    xb, base, base + sign * track_h * h / top, color=style.color, alpha=0.55, lw=0, zorder=2
                )
                ax.plot(
                    [0, lengths[sp_key][chrom] / 1e6],
                    [base, base],
                    color=style.color,
                    lw=0.4,
                    alpha=0.6,
                    zorder=2,
                )
        for g in np.arange(0, xmax + step, step):
            ax.axvline(g, color="#B0B0B0", lw=0.5, ls=(0, (4, 3)), zorder=0)
        ax.set_ylim(lower - 0.05, upper + 0.05)
        ax.set_yticks([])
        ax.set_xlim(0, xmax)
        for side in ("left", "right", "top", "bottom"):
            ax.spines[side].set_visible(False)
        ax.tick_params(axis="x", colors="#555555", labelsize=7, length=0)
        if (rc, qc) != chroms[-1]:
            ax.set_xticklabels([])
        ax.text(
            -0.015,
            0.5,
            rc if rc == qc else f"{rc}/{qc}",
            transform=ax.transAxes,
            ha="right",
            va="center",
            fontsize=9,
        )
        ax.text(
            lengths[ref.key][rc] / 1e6 + pad,
            y_ref,
            ref.label,
            va="center",
            fontsize=7,
            style="italic",
            color=REF_GREY,
        )
        ax.text(
            lengths[qry.key][qc] / 1e6 + pad,
            y_qry,
            qry.label,
            va="center",
            fontsize=7,
            style="italic",
            color=QRY_GREY,
        )
    axs[-1, 0].set_xlabel("Chromosome position (Mbp)", fontsize=9)
    genomes = [
        Line2D([], [], color=REF_GREY, lw=2, label=ref.label),
        Line2D([], [], color=QRY_GREY, lw=2, label=qry.label),
    ]
    annotations = [Patch(fc=BLOCK_COLORS[k], label=v) for k, v in BLOCK_LABELS.items()]
    if any(centromeres.values()):
        annotations.append(Patch(fc=CREAM, label="Centromere"))
    runs = [
        Line2D([], [], color=s.color, lw=1.2, marker="o", ms=4, mec="white", mew=0.4, label=s.label)
        for s, *_ in tracks
    ]
    for i, (handles, title, x) in enumerate(
        ((genomes, "Genomes", 0.10), (annotations, "Annotations", 0.38), (runs, "Syntenic orthologs", 0.74))
    ):
        fig.legend(
            handles=handles,
            title=title,
            loc="lower center",
            bbox_to_anchor=(x, 1 - 0.3 / height),
            ncol=1 if i == 0 else 2,
            frameon=False,
            fontsize=7.5,
            title_fontsize=8,
            alignment="left",
            handlelength=1.4,
            columnspacing=1.0,
        )
    fig.tight_layout(rect=(0.02, 0, 0.97, 1 - 0.9 / height), h_pad=0.1)
    _save(fig, prefix, ("pdf", "png"))


def _grid_step(span: float, target: int = 10) -> float:
    raw = span / target
    magnitude = float(10 ** np.floor(np.log10(raw)))
    return next(m * magnitude for m in (1, 2, 2.5, 5, 10) if m * magnitude >= raw)


def genome_dotplot(
    alignments: Sequence[Alignment],
    species: Sequence[SpeciesStyle],
    lengths: Mapping[str, Mapping[str, int]],
    homologous: set[tuple[str, str]],
    links: Sequence[tuple[RunStyle, list[tuple[str, float, str, float]]]],
    prefix: Path,
) -> None:
    ref, qry = species
    ref_chroms = sorted(lengths[ref.key], key=natural_key)
    partner = dict(homologous)
    qry_chroms = [partner[c] for c in ref_chroms if partner.get(c) in lengths[qry.key]]
    qry_chroms += sorted((c for c in lengths[qry.key] if c not in qry_chroms), key=natural_key)
    xo, yo, x, y = {}, {}, 0, 0
    for c in ref_chroms:
        xo[c], x = x, x + lengths[ref.key][c]
    for c in qry_chroms:
        yo[c], y = y, y + lengths[qry.key][c]
    fig = Figure(figsize=(9, 9))
    ax = fig.subplots()
    segments: dict[str, list[tuple[tuple[float, float], tuple[float, float]]]] = defaultdict(list)
    for a in alignments:
        if a.ref_chr in xo and a.qry_chr in yo:
            start = ((xo[a.ref_chr] + a.ref_start) / 1e6, (yo[a.qry_chr] + a.qry_start) / 1e6)
            end = ((xo[a.ref_chr] + a.ref_end) / 1e6, (yo[a.qry_chr] + a.qry_end) / 1e6)
            segments[a.kind].append((start, end))
    for kind, segs in segments.items():
        ax.add_collection(
            LineCollection(segs, colors=DOT_COLORS[kind], linewidths=1.8, zorder=2, capstyle="round")
        )
    for style, pairs in links:
        pts = [((xo[a] + am) / 1e6, (yo[b] + bm) / 1e6) for a, am, b, bm in pairs if a in xo and b in yo]
        if pts:
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            ax.scatter(
                xs,
                ys,
                s=5,
                c=style.color,
                alpha=0.85,
                zorder=3,
                label=style.label,
                edgecolors="none",
            )
    for pos in list(xo.values())[1:]:
        ax.axvline(pos / 1e6, color="#CCCCCC", lw=0.5, zorder=0)
    for pos in list(yo.values())[1:]:
        ax.axhline(pos / 1e6, color="#CCCCCC", lw=0.5, zorder=0)
    ax.set_xlim(0, x / 1e6)
    ax.set_ylim(0, y / 1e6)
    ax.set_xticks(
        [(xo[c] + lengths[ref.key][c] / 2) / 1e6 for c in ref_chroms], ref_chroms, rotation=90, fontsize=7
    )
    ax.set_yticks([(yo[c] + lengths[qry.key][c] / 2) / 1e6 for c in qry_chroms], qry_chroms, fontsize=7)
    ax.set_xlabel(f"{ref.label} (reference)", style="italic")
    ax.set_ylabel(f"{qry.label} (query)", style="italic")
    handles = [Line2D([], [], color=DOT_COLORS[k], lw=1.5, label=v) for k, v in BLOCK_LABELS.items()]
    handles += [Line2D([], [], color=s.color, marker="o", lw=0, ms=4, label=s.label) for s, _ in links]
    ax.legend(handles=handles, loc="upper left", fontsize=7, frameon=False)
    _save(fig, prefix, ("png", "pdf"), dpi=300)
