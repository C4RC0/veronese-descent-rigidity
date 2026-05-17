"""
Figure 4 — Five-tetrahedra hinge bridge.

Three panels + a certificate box that summarise the technical heart of
§5.

(A) T_0 standard cube-inscribed tetrahedron in R^3:
        a_1 = (+1, +1, +1) / sqrt(3),
        a_2 = (+1, -1, -1) / sqrt(3),
        a_3 = (-1, +1, -1) / sqrt(3),
        a_4 = (-1, -1, +1) / sqrt(3).
    3D projection (oblique) with the 4 axes drawn from the origin.

(B) K_5-incidence diagram of the 5 projective tetrahedra T_0..T_4:
    every pair (T_k, T_l) shares exactly one of the 10 = C(5,2)
    projective axes, labelled by {k, l}.

(C) Hinge parameterisation: T_k = T_k(eps_k, phi_k), with 6 hinge
    constraints x_{kl} = +/- x_{lk}.

(D) Exact certificate summary box.
"""

from __future__ import annotations

import itertools
import math
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "..", "paper", "figures")

SQ3 = math.sqrt(3.0)
A_VERTS = np.array([
    [+1, +1, +1],
    [+1, -1, -1],
    [-1, +1, -1],
    [-1, -1, +1],
], dtype=float) / SQ3


# ----- Panel A: T_0 in R^3 -----
def draw_T0(ax):
    # Draw 4 axes from origin to a_k.
    for k, a in enumerate(A_VERTS):
        ax.quiver(0, 0, 0, a[0], a[1], a[2],
                  color="#444", linewidth=1.4,
                  arrow_length_ratio=0.10)
        ax.text(a[0] * 1.18, a[1] * 1.18, a[2] * 1.18,
                f"$a_{{{k+1}}}$", fontsize=10, fontweight="bold")
    # Tetrahedron faces (4 of them) for visual context.
    faces = [
        [A_VERTS[0], A_VERTS[1], A_VERTS[2]],
        [A_VERTS[0], A_VERTS[1], A_VERTS[3]],
        [A_VERTS[0], A_VERTS[2], A_VERTS[3]],
        [A_VERTS[1], A_VERTS[2], A_VERTS[3]],
    ]
    poly = Poly3DCollection(faces, alpha=0.13, facecolor="#88aaff",
                            edgecolor="#5577cc", linewidth=0.6)
    ax.add_collection3d(poly)

    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_zlim(-1, 1)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.view_init(elev=18, azim=30)
    ax.set_title("(A) Standard tetrahedron $T_0$\n"
                 "$a_k = \\frac{1}{\\sqrt{3}}(\\pm 1, \\pm 1, \\pm 1)$, "
                 "$\\prod \\mathrm{signs} = +1$", fontsize=9)


# ----- Panel B: K_5-incidence -----
def draw_K5_incidence(ax):
    # 5 tetrahedra around a pentagon.
    n = 5
    R = 1.0
    angles = [-math.pi / 2 + k * 2 * math.pi / n for k in range(n)]
    pos = [(R * math.cos(a), R * math.sin(a)) for a in angles]

    # Edges: every pair (k, l) carries the Kneser label {k, l}.
    for i in range(n):
        for j in range(i + 1, n):
            p1, p2 = pos[i], pos[j]
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]],
                    color="#888", lw=0.8, zorder=1)
            mx = (p1[0] + p2[0]) / 2
            my = (p1[1] + p2[1]) / 2
            # Push the label slightly off the line so it stays readable.
            ax.text(mx, my, f"$\\{{{i},{j}\\}}$",
                    fontsize=7.5, ha="center", va="center",
                    bbox=dict(boxstyle="round,pad=0.15", facecolor="#fff",
                              edgecolor="#aaa", linewidth=0.4),
                    zorder=2)

    # 5 vertices labelled T_0..T_4.
    for k, (x, y) in enumerate(pos):
        ax.scatter([x], [y], s=380, c="#ffe5a8", edgecolors="#aa8844",
                   linewidths=1.4, zorder=3)
        ax.text(x, y, f"$T_{k}$", fontsize=11, fontweight="bold",
                ha="center", va="center", zorder=4)

    ax.set_xlim(-1.45, 1.45)
    ax.set_ylim(-1.35, 1.35)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("(B) $K_5$-incidence of the 5 projective tetrahedra\n"
                 "every pair $(T_k, T_l)$ shares one axis labelled "
                 "$\\{k, l\\}$;  total $\\binom{5}{2} = 10$ axes",
                 fontsize=9)


# ----- Panel C: Hinge parameterisation -----
def draw_hinge_param(ax):
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    ax.set_title("(C) Hinge parameterisation\n"
                 "$T_k = T_k(\\varepsilon_k, \\varphi_k)$, "
                 "6 hinge constraints", fontsize=9)

    # T_0 box
    ax.add_patch(FancyBboxPatch(
        (3.4, 7.8), 3.2, 1.2,
        boxstyle="round,pad=0.05,rounding_size=0.18",
        facecolor="#dbe9f8", edgecolor="#446699", linewidth=1.0))
    ax.text(5.0, 8.4, "$T_0$ fixed", ha="center", va="center",
            fontsize=10, fontweight="bold")

    # T_k boxes (k = 1..4)
    tk_positions = [(0.6, 5.0), (3.0, 5.0), (5.4, 5.0), (7.8, 5.0)]
    for k, (x, y) in enumerate(tk_positions, start=1):
        ax.add_patch(FancyBboxPatch(
            (x, y), 1.6, 1.6,
            boxstyle="round,pad=0.04,rounding_size=0.14",
            facecolor="#f5e9c8", edgecolor="#aa8844", linewidth=0.9))
        ax.text(x + 0.8, y + 1.15, f"$T_{k}$", ha="center", va="center",
                fontsize=11, fontweight="bold")
        ax.text(x + 0.8, y + 0.65,
                f"$\\varepsilon_{k}, \\varphi_{k}$",
                ha="center", va="center", fontsize=9)
        ax.text(x + 0.8, y + 0.25,
                "axis $a_{" + str(k) + "}$ shared",
                ha="center", va="center", fontsize=7.5, color="#444")
        # Connector to T_0
        ax.plot([x + 0.8, 5.0], [y + 1.6, 7.8], color="#aaa", lw=0.6,
                zorder=0)

    # Hinge-pair edges T_k <-> T_l for 1<=k<l<=4
    pairs = list(itertools.combinations(range(4), 2))
    for (k, l) in pairs:
        x1 = tk_positions[k][0] + 0.8
        y1 = tk_positions[k][1]
        x2 = tk_positions[l][0] + 0.8
        y2 = tk_positions[l][1]
        # Draw the hinge edge with a slight curve down so they don't all
        # overlap the T_0 lines.
        ax.annotate("", xy=(x2, y2 - 0.15), xytext=(x1, y1 - 0.15),
                    arrowprops=dict(arrowstyle="-", color="#cc6644",
                                    lw=0.8,
                                    connectionstyle="arc3,rad=0.18"))

    ax.text(5.0, 3.0,
            "6 hinge equations:  "
            "$u_k^{(\\pi_k(l))} = \\delta_{kl}\\, u_l^{(\\pi_l(k))}$",
            ha="center", va="center", fontsize=9)
    ax.text(5.0, 2.0,
            "$\\Rightarrow$  18 scalar equations in "
            "$(c_1, s_1, ..., c_4, s_4)$",
            ha="center", va="center", fontsize=9, style="italic")
    ax.text(5.0, 1.0,
            "rank $= 7$ (Theorem 5.4)",
            ha="center", va="center", fontsize=9,
            fontweight="bold", color="#3a7d3a")


# ----- Panel D: Certificate box -----
def draw_certificate(ax):
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    ax.set_title("(D) Exact 1024-case certificate (Theorem 5.5)",
                 fontsize=9)

    # Outer box
    ax.add_patch(FancyBboxPatch(
        (0.3, 0.5), 9.4, 9.0,
        boxstyle="round,pad=0.08,rounding_size=0.25",
        facecolor="#f7f3e6", edgecolor="#aa8844", linewidth=1.0))

    # Partition rows
    rows = [
        ("1024 cases scanned", "", "#222"),
        ("", "", "#222"),
        ("inconsistent linear",        " 992  (96.9%)",  "#666"),
        ("full-rank, unit-norm invalid", "  16  (1.6%)",  "#cc6644"),
        ("rank-7 canonical quadratic",   "  16  (1.6%)",  "#3a7d3a"),
        ("", "", "#222"),
        ("partition TOTAL",            "1024",            "#222"),
        ("", "", "#222"),
        ("$16 \\times 2 = 32$ raw algebraic solutions", "", "#222"),
        ("$\\Rightarrow$ icosahedral orbit", "", "#3a7d3a"),
        ("(exact $(\\sigma, D)$ + symmetry reduction)", "", "#3a7d3a"),
    ]

    y_top = 8.7
    line_h = 0.75
    for k, (left, right, color) in enumerate(rows):
        y = y_top - k * line_h
        if left:
            weight = "bold" if k in (0, 6, 9, 10) else "normal"
            ax.text(0.8, y, left, ha="left", va="center",
                    fontsize=9.0, color=color, fontweight=weight)
        if right:
            ax.text(8.9, y, right, ha="right", va="center",
                    fontsize=9.0, color=color, family="monospace")


def main():
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.18)

    axA = fig.add_subplot(gs[0, 0], projection="3d")
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, 0])
    axD = fig.add_subplot(gs[1, 1])

    draw_T0(axA)
    draw_K5_incidence(axB)
    draw_hinge_param(axC)
    draw_certificate(axD)

    fig.suptitle("Five-tetrahedra hinge bridge (§5)", fontsize=12, y=0.995)

    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, "fig4_hinge_bridge.pdf"),
                bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "fig4_hinge_bridge.png"),
                bbox_inches="tight", dpi=200)
    print(f"Wrote {os.path.join(OUT, 'fig4_hinge_bridge.pdf')}")
    print(f"Wrote {os.path.join(OUT, 'fig4_hinge_bridge.png')}")


if __name__ == "__main__":
    main()
