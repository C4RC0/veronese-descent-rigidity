"""
Figure 2 — Petersen / Kneser / ETF layer.

Three panels visualising Theorem 1 (the Bose-Mesner identification).

(A) Petersen K(5,2): 10 vertices labelled by {a,b} in C(5,2); edge iff
    the two 2-subsets are disjoint. Classical Petersen-star drawing.

(B) Two-valued squared-Gram pattern:
        Q_ij = 5/9   on Petersen edges
        Q_ij = 1/9   on complement edges  (i != j)
    rendered as a 10x10 heatmap.

(C) Centred lift K = (3/2)(Q - J/3) = 2 E_1, with spectrum {2^5, 0^5}.
    Bar chart of the 10 eigenvalues.
"""

from __future__ import annotations

import itertools
import os

import matplotlib.pyplot as plt
import numpy as np


HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "..", "paper", "figures")


# ----- Petersen graph K(5,2) -----
LABELS = list(itertools.combinations(range(5), 2))
N = 10
ADJ = np.zeros((N, N), dtype=int)
for i in range(N):
    for j in range(N):
        if i != j and set(LABELS[i]).isdisjoint(LABELS[j]):
            ADJ[i, j] = 1
J = np.ones((N, N))
I_N = np.eye(N)

Q = (5/9) * ADJ + (1/9) * (J - I_N - ADJ) + I_N      # squared Gram
K_mat = (3/2) * (Q - J / 3)                          # K = 2 E_1


# Kneser-consistent permutation for the classical Petersen-star drawing
# (vertex 0 -> slot 0, etc.). Hard-coded by enumeration.
def find_kneser_layout_perm():
    classic_edges = {
        (0, 1), (1, 2), (2, 3), (3, 4), (0, 4),
        (0, 5), (1, 6), (2, 7), (3, 8), (4, 9),
        (5, 7), (7, 9), (6, 9), (6, 8), (5, 8),
    }

    def is_classic(p, q):
        a, b = (p, q) if p < q else (q, p)
        return (a, b) in classic_edges

    pos = [-1] * N
    used = [False] * N

    def back(i):
        if i == N:
            return True
        for p in range(N):
            if used[p]:
                continue
            ok = True
            for j in range(i):
                kneser = ADJ[i, j] == 1
                cl = is_classic(p, pos[j])
                if kneser != cl:
                    ok = False
                    break
            if not ok:
                continue
            pos[i] = p
            used[p] = True
            if back(i + 1):
                return True
            pos[i] = -1
            used[p] = False
        return False

    back(0)
    return pos


def petersen_positions():
    perm = find_kneser_layout_perm()
    slots = []
    R_out, R_in = 1.0, 0.45
    for k in range(5):
        a = -np.pi / 2 + k * 2 * np.pi / 5
        slots.append((R_out * np.cos(a), R_out * np.sin(a)))
    for k in range(5):
        a = -np.pi / 2 + (k + 0.5) * 2 * np.pi / 5
        slots.append((R_in * np.cos(a), R_in * np.sin(a)))
    return [slots[perm[i]] for i in range(N)]


def draw_petersen(ax):
    pos = petersen_positions()
    for i in range(N):
        for j in range(i + 1, N):
            if ADJ[i, j]:
                ax.plot([pos[i][0], pos[j][0]], [pos[i][1], pos[j][1]],
                        color="#cc4444", lw=1.6, zorder=1)
    for i, (x, y) in enumerate(pos):
        ax.scatter([x], [y], s=300, c="#fff5d6", edgecolors="#333",
                   linewidths=1.2, zorder=2)
        a, b = LABELS[i]
        ax.text(x, y, f"{{{a},{b}}}", ha="center", va="center",
                fontsize=8.5, zorder=3)
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.3, 1.3)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("(A) Petersen graph $K(5,2)$\n"
                 "vertex = 2-subset $\\{a,b\\} \\subset \\{0,..,4\\}$; "
                 "edge $\\Longleftrightarrow$ disjoint", fontsize=9)


def draw_Q_heatmap(ax):
    im = ax.imshow(Q, cmap="viridis", vmin=0, vmax=1,
                   interpolation="nearest")
    ax.set_xticks(range(N))
    ax.set_yticks(range(N))
    ax.set_xticklabels([f"{a}{b}" for a, b in LABELS], fontsize=7)
    ax.set_yticklabels([f"{a}{b}" for a, b in LABELS], fontsize=7)
    ax.set_title("(B) Squared-Gram $Q_{ij} = (v_i \\cdot v_j)^2$\n"
                 "$5/9$ on Petersen edges, $1/9$ on complement",
                 fontsize=9)
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.ax.tick_params(labelsize=7)
    cb.set_ticks([1/9, 5/9, 1.0])
    cb.set_ticklabels(["1/9", "5/9", "1"])


def draw_K_spectrum(ax):
    eigs = np.sort(np.linalg.eigvalsh(K_mat))[::-1]
    indices = np.arange(1, N + 1)
    colors = ["#3a7d3a" if abs(e - 2) < 1e-6 else "#aaaaaa" for e in eigs]
    ax.bar(indices, eigs, color=colors, edgecolor="#222", linewidth=0.6)
    ax.set_xticks(indices)
    ax.set_xticklabels([str(i) for i in indices], fontsize=8)
    ax.set_xlabel("eigenvalue index", fontsize=9)
    ax.set_ylabel("eigenvalue", fontsize=9)
    ax.axhline(2, color="#3a7d3a", lw=0.6, ls="--", alpha=0.5)
    ax.axhline(0, color="#222", lw=0.5)
    ax.set_ylim(-0.4, 2.5)
    ax.set_title("(C) $K = \\frac{3}{2}(Q - J/3) = 2\\,E_1$\n"
                 "$\\mathrm{spec}(K) = \\{2^{(5)}, 0^{(5)}\\}$ "
                 "(Theorem 1)", fontsize=9)


def main():
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.5),
                              gridspec_kw={"width_ratios": [1, 1.3, 1]})
    draw_petersen(axes[0])
    draw_Q_heatmap(axes[1])
    draw_K_spectrum(axes[2])

    fig.tight_layout()
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, "fig2_petersen_scheme.pdf"),
                bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "fig2_petersen_scheme.png"),
                bbox_inches="tight", dpi=200)
    print(f"Wrote {os.path.join(OUT, 'fig2_petersen_scheme.pdf')}")
    print(f"Wrote {os.path.join(OUT, 'fig2_petersen_scheme.png')}")


if __name__ == "__main__":
    main()
