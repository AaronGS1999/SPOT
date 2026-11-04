"""Gene tree viewer: draws a Newick tree as a rectangular phylogram."""

from __future__ import annotations

from io import StringIO

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PySide6.QtWidgets import QLabel, QScrollArea, QVBoxLayout, QWidget


def read_newick(text: str) -> Clade:
    """Root clade of a Newick string (Biopython)."""
    tree = Phylo.read(StringIO(text.strip()), "newick")  # type: ignore[attr-defined, no-untyped-call]
    root: Clade = tree.root
    return root


def layout(root: Clade) -> dict[Clade, tuple[float, float]]:
    """(x, y) of every clade: x = distance from the root, y = leaf order, parents centred on children."""
    pos: dict[Clade, tuple[float, float]] = {}
    leaves = 0

    def walk(clade: Clade, depth: float) -> float:
        nonlocal leaves
        if not clade.clades:
            y = float(leaves)
            leaves += 1
        else:
            ys = [walk(c, depth + max(c.branch_length or 0.0, 0.0)) for c in clade.clades]
            y = (ys[0] + ys[-1]) / 2
        pos[clade] = (depth, y)
        return y

    walk(root, 0.0)
    return pos


class TreeView(QWidget):
    """Draws a tree; leaves are coloured by species and the selected pair is highlighted."""

    def __init__(self) -> None:
        super().__init__()
        self.message = QLabel("Select a confirmed pair to see the gene tree of its orthogroup.")
        self.message.setWordWrap(True)
        self.area = QScrollArea()
        self.area.setWidgetResizable(True)
        self.canvas: FigureCanvasQTAgg | None = None
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.message)
        layout.addWidget(self.area, 1)

    def clear(self, text: str) -> None:
        self.message.setText(text)
        self.message.show()
        self.area.takeWidget()
        self.canvas = None

    def show_tree(
        self, newick: str, species_of: dict[str, str], colors: dict[str, str], highlight: set[str], title: str
    ) -> None:
        try:
            root = read_newick(newick)
        except (ValueError, IndexError, Phylo.NewickIO.NewickError) as e:
            self.clear(f"Cannot read this tree: {e}")
            return
        pos = layout(root)
        leaves: list[Clade] = root.get_terminals()  # type: ignore[no-untyped-call]
        width = max(pos[n][0] for n in leaves) or 1.0
        height_in = max(2.5, 0.17 * len(leaves) + 0.8)
        fig = Figure(figsize=(5.5, height_in), dpi=100)
        ax = fig.add_axes((0.02, 0.4 / height_in, 0.62, 1 - 0.8 / height_in))

        def draw(n: Clade) -> None:
            x, y = pos[n]
            for child in n.clades:
                cx, cy = pos[child]
                ax.plot([x, x], [y, cy], color="#555555", lw=0.8, solid_capstyle="butt")
                ax.plot([x, cx], [cy, cy], color="#555555", lw=0.8, solid_capstyle="butt")
                draw(child)
            support = n.confidence if n.confidence is not None else n.name
            if n.clades and support is not None and _is_number(str(support)):
                ax.text(
                    x,
                    y,
                    f" {support:g}" if isinstance(support, float) else f" {support}",
                    fontsize=6,
                    color="#888888",
                    va="center",
                )

        draw(root)
        for leaf in leaves:
            lx, ly = pos[leaf]
            sp = species_of.get(leaf.name, "")
            color = colors.get(sp, "#333333")
            picked = leaf.name in highlight
            ax.plot(lx, ly, "o", ms=4 if picked else 2.5, color=color, zorder=3)
            ax.text(
                lx + width * 0.02,
                ly,
                leaf.name,
                fontsize=7,
                va="center",
                color=color,
                fontweight="bold" if picked else "normal",
                bbox={"fc": "#FFF3B0", "ec": "none", "pad": 0.6} if picked else None,
                clip_on=False,
            )
        ax.set_ylim(len(leaves) - 0.5, -0.5)
        ax.set_xlim(-width * 0.01, width * 1.01)
        ax.axis("off")
        ax.plot([0, width * 0.1], [len(leaves) - 0.1] * 2, color="#555555", lw=0.8, clip_on=False)
        ax.text(
            width * 0.05,
            len(leaves) + 0.1,
            f"{width * 0.1:.2g}",
            fontsize=6,
            ha="center",
            va="top",
            color="#555555",
        )
        ax.set_title(title, fontsize=8, loc="left")
        self.canvas = FigureCanvasQTAgg(fig)  # type: ignore[no-untyped-call]
        self.canvas.setMinimumHeight(int(height_in * 100))
        self.area.setWidget(self.canvas)
        self.message.hide()


def _is_number(text: str) -> bool:
    try:
        float(text)
    except ValueError:
        return False
    return True
