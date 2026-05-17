"""
Figure 1 — Proof cascade.

A vertical chain of boxes summarising the main deductive route of the
paper:

    Petersen K(5,2)
        |  Theorem 1 (Bose-Mesner)
    Primitive idempotent E_1  +  10-line ETF in R^5
        |  Theorem 5.6 (five-tetrahedra hinge bridge)
    A_5-equivariant Veronese lift
        |  Theorem 4.3 (Schur descent)
    Unique RP^2 configuration up to O(3) x S_5
        |  Theorem 3 (corollary)
    Signed Gram rigidity

Reads-in-one-glance: rigid all the way down.
"""

from __future__ import annotations

import os

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "..", "paper", "figures")


def main():
    fig, ax = plt.subplots(figsize=(7.6, 8.2))
    ax.set_xlim(0, 11)
    ax.set_ylim(0.8, 11)
    ax.axis("off")

    # Each entry: (centre_y, title, subtitle, fillcolor, box_h_override)
    boxes = [
        (10.0, "Petersen graph $K(5,2)$",
         "$10$ vertices, adjacency $A$, all-ones $J$", "#d9e7f5", None),
        ( 8.2, "Petersen primitive idempotent $E_1$ — 10-line ETF in $\\mathbb{R}^5$",
         "$K = \\frac{3}{2}(Q - J/3) = 2\\,E_1$ (Theorem 1)", "#d9f0d9", None),
        ( 6.4, "$A_5$-equivariant Veronese lift",
         "five-tetrahedra hinge bridge (Theorem 5.6)", "#f5e9c8", None),
        ( 4.6, "Unique $\\mathbb{RP}^2$ configuration up to $O(3) \\times S_5$",
         "Schur descent (Theorem 4.3): comparison is $\\pm I$; $-I$ kills the Veronese surface",
         "#f3d9d9", None),
        ( 2.5, "Signed Gram rigidity (Theorem 3)",
         "$H = P_\\sigma\\,D\\,H_{\\mathrm{ico}}\\,D\\,P_\\sigma^T$",
         "#e6d9f3", None),
    ]

    arrow_labels = [
        "Theorem 1: $K = 2E_1$\n(Bose–Mesner)",
        "Theorem 5.6: five-tetrahedra\nhinge bridge",
        "Theorem 4.3: Schur descent\n+ $-I$ spectral obstruction",
        "Theorem 3:\nrank-$3$ Gram factorisation",
    ]

    box_w = 7.8
    box_h = 1.05
    cx = 5.5

    drawn_boxes = []
    for y, title, sub, color, _ in boxes:
        box = FancyBboxPatch(
            (cx - box_w / 2, y - box_h / 2),
            box_w, box_h,
            boxstyle="round,pad=0.04,rounding_size=0.18",
            linewidth=1.0,
            edgecolor="#333",
            facecolor=color,
        )
        ax.add_patch(box)
        ax.text(cx, y + 0.18, title, ha="center", va="center",
                fontsize=9.5, fontweight="bold")
        ax.text(cx, y - 0.27, sub, ha="center", va="center",
                fontsize=7.8, color="#333")
        drawn_boxes.append(y)

    # Arrows + labels
    for i, label in enumerate(arrow_labels):
        y_top = drawn_boxes[i] - box_h / 2
        y_bot = drawn_boxes[i + 1] + box_h / 2
        arrow = FancyArrowPatch(
            (cx, y_top - 0.03),
            (cx, y_bot + 0.03),
            arrowstyle="-|>",
            mutation_scale=14,
            linewidth=1.2,
            color="#444",
        )
        ax.add_patch(arrow)
        ax.text(cx + box_w / 2 + 0.25, (y_top + y_bot) / 2, label,
                ha="left", va="center", fontsize=7.8, color="#444",
                style="italic")

    ax.set_title("Proof cascade", fontsize=12, pad=12)

    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, "fig1_proof_cascade.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "fig1_proof_cascade.png"), bbox_inches="tight", dpi=200)
    print("Wrote paper/figures/fig1_proof_cascade.pdf")
    print("Wrote paper/figures/fig1_proof_cascade.png")


if __name__ == "__main__":
    main()
