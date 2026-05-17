"""
Figure 3 — Veronese descent / Schur obstruction.

Two side-by-side panels visualising Lemmas 3.1 + 4.2 + Schur (4.1).

(A) Veronese image  U = v v^T - I/3,  with spectrum
       spec(U) = (2/3, -1/3, -1/3)
    and equivalently  U + I/3  has spectrum (1, 0, 0): rank-1 PSD.

(B) Negative branch  -U,  with spectrum
       spec(-U) = (-2/3, 1/3, 1/3),
    so that  -U + I/3  has spectrum (-1/3, 2/3, 2/3): one NEGATIVE
    eigenvalue, hence not PSD, hence -U is NOT in the Veronese image.

The figure renders both spectra as bar charts plus a "PSD?" badge.
"""

from __future__ import annotations

import os

import matplotlib.pyplot as plt
import numpy as np


HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "..", "paper", "figures")


def spec_bars(ax, spec, color, title, subtitle, is_psd, shift_label):
    indices = np.arange(len(spec))
    ax.bar(indices, spec, color=color, edgecolor="#222", linewidth=0.8,
           width=0.6)
    ax.axhline(0, color="#222", lw=0.6)
    ax.set_ylim(-1.0, 1.2)
    ax.set_xticks(indices)
    ax.set_xticklabels([f"$\\lambda_{i+1}$" for i in indices], fontsize=10)
    ax.set_ylabel("eigenvalue", fontsize=9)
    ax.set_title(f"{title}\n{subtitle}", fontsize=10, pad=10)
    for i, v in enumerate(spec):
        ax.text(i, v + (0.04 if v >= 0 else -0.10),
                f"{v:+.3g}", ha="center", va="bottom" if v >= 0 else "top",
                fontsize=9)
    badge_color = "#cce8d1" if is_psd else "#f5cccc"
    badge_text = ("PSD: rank-1 projector\n" + shift_label) if is_psd \
        else ("NOT PSD\n" + shift_label)
    ax.text(0.5, -0.32, badge_text, transform=ax.transAxes, ha="center",
            fontsize=9.5,
            bbox=dict(boxstyle="round,pad=0.4", facecolor=badge_color,
                      edgecolor="#333", linewidth=0.6))


def main():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 5.0))

    spec_U = np.array([2 / 3, -1 / 3, -1 / 3])
    spec_minus_U_plus_I3 = np.array([-1 / 3, 2 / 3, 2 / 3])

    spec_bars(
        ax1, spec_U,
        color="#88c088",
        title="(A) $U = v v^T - I/3$ is a centred Veronese point",
        subtitle="$\\mathrm{spec}(U) = (2/3, -1/3, -1/3)$",
        is_psd=True,
        shift_label=("$U + I/3 = v v^T$:  $(1, 0, 0)$\n"
                     "rank-1 PSD $\\Rightarrow$ Veronese image"),
    )

    spec_bars(
        ax2, spec_minus_U_plus_I3,
        color="#c08888",
        title="(B) $-U$ is not a centred Veronese point",
        subtitle="$\\mathrm{spec}(-U + I/3) = (-1/3, 2/3, 2/3)$",
        is_psd=False,
        shift_label=("negative eigenvalue $\\Rightarrow$ not PSD\n"
                     "$-U + I/3$ cannot equal $v v^T$"),
    )

    # Schur diagram between the two panels (annotation only).
    fig.suptitle("Veronese descent / Schur obstruction\n"
                 "(Lemma 3.1, Lemma 4.2, Theorem 4.3)",
                 fontsize=11, y=1.02)
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])
    # Schur summary box at the bottom, centred.
    fig.text(0.5, -0.02,
             "Schur (Lemma 4.1): an $A_5$-equivariant isometry "
             "$T : \\mathbb{R}^5 \\to \\mathbb{R}^5$ is $\\pm I$. "
             "The case $T = -I$ violates (B): only $T = +I$ survives, "
             "giving $O(3)$-rigidity (Theorem 4.3).",
             ha="center", fontsize=9, color="#222",
             bbox=dict(boxstyle="round,pad=0.4", facecolor="#fff8e1",
                       edgecolor="#aa8844", linewidth=0.8))

    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, "fig3_schur_obstruction.pdf"),
                bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "fig3_schur_obstruction.png"),
                bbox_inches="tight", dpi=200)
    print("Wrote paper/figures/fig3_schur_obstruction.pdf")
    print("Wrote paper/figures/fig3_schur_obstruction.png")


if __name__ == "__main__":
    main()
