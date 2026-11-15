"""Draw the seven-vertex recursive tree and the two pairs used in Figure 1."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


PARENTS = {2: 1, 3: 1, 4: 2, 5: 3, 6: 4, 7: 4}
POSITIONS = {
    1: (0.0, 2.15), 2: (-0.7, 1.45), 3: (0.7, 1.45),
    4: (-0.7, 0.77), 5: (0.7, 0.77),
    6: (-1.10, 0.08), 7: (-0.30, 0.08),
}
ROOT = Path(__file__).resolve().parents[1]


def lineage(vertex):
    result = [vertex]
    while vertex in PARENTS:
        vertex = PARENTS[vertex]
        result.append(vertex)
    return result


def draw_pair(ax, u, v, panel):
    first, second = lineage(u), lineage(v)
    common = next(label for label in first if label in second)
    path = first[:first.index(common) + 1] + second[:second.index(common)][::-1]
    assert len(path) - 1 == 2
    selected_edges = {frozenset(edge) for edge in zip(path, path[1:])}
    ink, highlight, inactive = "#243541", "#125d76", "#bbc1c5"
    for child, parent in PARENTS.items():
        active = frozenset((parent, child)) in selected_edges
        ax.plot(*zip(POSITIONS[parent], POSITIONS[child]),
                color=highlight if active else inactive,
                lw=2.6 if active else 1.1, solid_capstyle="round", zorder=1)
    for label, (x, y) in POSITIONS.items():
        endpoint = label in (u, v)
        ax.scatter(x, y, s=178, facecolor=highlight if endpoint else "white",
                   edgecolor=highlight if label in path else inactive, linewidth=1.2, zorder=3)
        if label == common:
            ax.scatter(x, y, s=270, facecolor="none", edgecolor=highlight,
                       linewidth=1.15, zorder=4)
        ax.text(x, y, str(label), ha="center", va="center", fontsize=8.5,
                color="white" if endpoint else ink, zorder=5)
    ax.set_title(f"({panel})   Pair {u} and {v} meets at {common}", fontsize=9, loc="center", pad=2)
    ax.set(xlim=(-1.55, 1.25), ylim=(-0.20, 2.46), aspect="equal")
    ax.axis("off")


def main():
    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(6.0, 1.95))
    fig.subplots_adjust(left=0.04, right=0.96, bottom=0.03, top=0.86, wspace=0.06)
    draw_pair(axes[0], 2, 3, "a")
    draw_pair(axes[1], 6, 7, "b")
    output = ROOT / "figures" / "two_equal_paths.pdf"
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, metadata={"Title": "Equal distances and different ancestors",
                                  "CreationDate": None, "ModDate": None})
    plt.close(fig)


if __name__ == "__main__":
    main()
