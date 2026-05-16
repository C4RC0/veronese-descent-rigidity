# Supplementary Material — Veronese Descent Rigidity of the Petersen ETF

This document is the auditor-facing companion to the manuscript. It records
the reproduction recipe, the script-to-theorem mapping, the expected
outputs (verbatim log files), the SHA-256 hashes of every artefact, and the
runtime table.

---

## A. Reproducibility manifest

### A.1 Software environment

- Python 3.10+ (developed under CPython 3.14)
- `numpy >= 1.24`
- `sympy >= 1.13`
- `mpmath >= 1.3`
- `matplotlib >= 3.7`

No compiled extensions, no external CAS or LP solvers, no GPU. All scripts
run on a single CPU thread.

Reproducing the environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### A.2 Script-to-theorem mapping

| Script | Role |
|--------|------|
| `verify_idempotent_kernel.py` | Independent numerical check of Theorem 1 and Lemma 5.1 |
| `p5_jordan_descent.py` | Independent numerical check of Lemma 3.1, Lemma 4.2, and the $\mathbf{1} \oplus \mathbf{4} \oplus \mathbf{5}$ character decomposition |
| `five_tet_hinge_exact.py` | Theorem 5.4 — exact derivation of the canonical quadratic on the representative branch |
| `five_tet_hinge_enum.py` | Theorem 5.5 — exhaustive 1024-case partition; exact-rank verification on all cases (`EXACT_RANK_MODE = True` by default) |
| `five_tet_orbit_exact.py` | Theorem 5.5 — exact $(\sigma, D)$ orbit certificate for the canonical chirality branch |
| `five_tet_hinge.py` | Diagnostic floating Newton/least-squares cross-check |
| `sgr_lemma_proof.py` | Cross-check of Theorem 3 (POCS-based) |
| `sgr_lemma_proof_exact.py` | Cross-check of Theorem 3 (exact sympy backtracking; heavy optional) |
| `sgr6_lemma_proof.py` | Analogous 6-vertex equiangular case (Discussion §7.3, item 1) |
| `platonic_census.py` | Platonic-group orbit census (Discussion §7.1) |
| `edge_orbit_analysis.py` | 15-edge $A_5$-orbit analysis (Discussion §7.3, item 2) |

### A.3 Reproduction recipe

**Algebraic certificates** (Theorems 5.4 and 5.5):

```bash
python3 scripts/five_tet_hinge_exact.py
python3 scripts/five_tet_hinge_enum.py
python3 scripts/five_tet_orbit_exact.py
```

**Independent numerical checks** (Theorem 1, Lemmas 3.1, 4.2, 5.1):

```bash
python3 scripts/verify_idempotent_kernel.py
python3 scripts/p5_jordan_descent.py
```

**Cross-checks and discussion-section scripts**:

```bash
python3 scripts/sgr_lemma_proof.py
python3 scripts/sgr6_lemma_proof.py
python3 scripts/platonic_census.py
python3 scripts/edge_orbit_analysis.py
# Heavy optional (several minutes):
python3 scripts/sgr_lemma_proof_exact.py
# Diagnostic floating cross-check:
python3 scripts/five_tet_hinge.py
```

### A.4 Expected outputs

Verbatim log files for every script (captured with `python3 -u` and
`tee`) are stored in `supplement/logs/`. Each filename matches the
corresponding script: e.g. `supplement/logs/five_tet_hinge_enum.log`.

Key invariants the auditor should look for:

- `verify_idempotent_kernel.py` — $\|K - 2 E_1\| < 10^{-14}$ (machine precision).
- `p5_jordan_descent.py` — Veronese spectrum $(2/3, -1/3, -1/3)$ matches at $\leq 10^{-12}$; character multiplicities $(1, 1, 1)$ for $\mathbf{1} \oplus \mathbf{4} \oplus \mathbf{5}$.
- `five_tet_hinge_exact.py` — canonical quadratic roots $s_4 = (\sqrt{3} \pm \sqrt{15})/8$; icosahedral residuals exactly zero.
- `five_tet_hinge_enum.py` — partition `992 + 16 + 16 = 1024`; exact-rank verification, zero discrepancies.
- `five_tet_orbit_exact.py` — `SUCCESS: orbit equivalence verified exactly`; explicit $\tau \in S_5$ and $D \in \{\pm 1\}^{10}$.

### A.5 SHA-256 hashes

The file `MANIFEST.sha256` in the repository root lists the SHA-256 hashes
of every script, figure, and reproduction log. A reviewer can verify
artefact integrity by:

```bash
sha256sum -c MANIFEST.sha256
```

### A.6 Runtime table

Indicative wall-clock times on a single modern CPU thread:

| Script | Wall-clock |
|--------|------------|
| `verify_idempotent_kernel.py` | < 1 s |
| `p5_jordan_descent.py` | < 1 s |
| `five_tet_hinge_exact.py` | ~ 30 s |
| `five_tet_hinge_enum.py` (exact rank) | ~ 20 s |
| `five_tet_orbit_exact.py` | ~ 1 min |
| `five_tet_hinge.py` (numerical diagnostic) | ~ 10 s |
| `sgr_lemma_proof.py` | ~ 1 s |
| `sgr_lemma_proof_exact.py` (heavy optional) | ~ several min |
| `sgr6_lemma_proof.py` | ~ 1 s |
| `platonic_census.py` | < 1 s |
| `edge_orbit_analysis.py` | < 1 s |
