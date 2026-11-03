"""Gene tree viewer: draws a Newick tree as a rectangular phylogram."""

from __future__ import annotations

from dataclasses import dataclass, field

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PySide6.QtWidgets import QLabel, QScrollArea, QVBoxLayout, QWidget


@dataclass
class Node:
    name: str = ""
    length: float = 0.0
    children: list[Node] = field(default_factory=list)
    x: float = 0.0
    y: float = 0.0

    def leaves(self) -> list[Node]:
        if not self.children:
            return [self]
        return [leaf for child in self.children for leaf in child.leaves()]


def parse_newick(text: str) -> Node:
    """Parse a Newick string (names, branch lengths, internal labels, quoted names)."""
    text = text.strip().rstrip(";")
    pos = 0

    def name() -> str:
        nonlocal pos
        if pos < len(text) and text[pos] == "'":
            end = text.index("'", pos + 1)
            value, pos = text[pos + 1 : end], end + 1
            return value
        start = pos
        while pos < len(text) and text[pos] not in ",():;":
            pos += 1
        return text[start:pos].strip()

    def node() -> Node:
        nonlocal pos
        current = Node()
        if pos < len(text) and text[pos] == "(":
            pos += 1
            current.children.append(node())
            while text[pos] == ",":
                pos += 1
                current.children.append(node())
            if text[pos] != ")":
                raise ValueError(f"expected ')' at {pos}")
            pos += 1
        current.name = name()
        if pos < len(text) and text[pos] == ":":
            pos += 1
            start = pos
            while pos < len(text) and text[pos] not in ",();":
                pos += 1
            current.length = float(text[start:pos] or 0)
        return current

    root = node()
    return root


def layout(root: Node) -> None:
    """x = distance from the root; y = leaf order, parents centred on their children."""
    leaves = 0

    def walk(n: Node, depth: float) -> None:
        nonlocal leaves
        n.x = depth
        if not n.children:
            n.y = leaves
            leaves += 1
            return
        for child in n.children:
            walk(child, depth + max(child.length, 0.0))
        n.y = (n.children[0].y + n.children[-1].y) / 2

    walk(root, 0.0)


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
            root = parse_newick(newick)
        except (ValueError, IndexError) as e:
            self.clear(f"Cannot read this tree: {e}")
            return
        layout(root)
        leaves = root.leaves()
        width = max(n.x for n in leaves) or 1.0
        height_in = max(2.5, 0.17 * len(leaves) + 0.8)
        fig = Figure(figsize=(5.5, height_in), dpi=100)
        ax = fig.add_axes((0.02, 0.4 / height_in, 0.62, 1 - 0.8 / height_in))

        def draw(n: Node) -> None:
            for child in n.children:
                ax.plot([n.x, n.x], [n.y, child.y], color="#555555", lw=0.8, solid_capstyle="butt")
                ax.plot([n.x, child.x], [child.y, child.y], color="#555555", lw=0.8, solid_capstyle="butt")
                draw(child)
            if n.children and n.name and _is_number(n.name):
                ax.text(n.x, n.y, f" {n.name}", fontsize=6, color="#888888", va="center")

        draw(root)
        for leaf in leaves:
            sp = species_of.get(leaf.name, "")
            color = colors.get(sp, "#333333")
            picked = leaf.name in highlight
            ax.plot(leaf.x, leaf.y, "o", ms=4 if picked else 2.5, color=color, zorder=3)
            ax.text(
                leaf.x + width * 0.02,
                leaf.y,
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
