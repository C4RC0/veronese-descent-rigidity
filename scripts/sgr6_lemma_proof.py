"""
SGR Lemma analogue for the 6-vertex equiangular orbit of A_5.

Claim.  Let H in R^{6x6} be real symmetric with:
    (i)   H_ii = 1
    (ii)  |H_ij| = 1/sqrt(5) for all i != j
    (iii) H is positive semidefinite
    (iv)  rank(H) = 3
Then H equals (up to switching D = diag(+/-1) and a permutation
sigma in S_6) the icosahedral 6-vertex Gram matrix.

Strategy: same Schur-complement reduction as the 10-pair case
(sgr_lemma_proof_exact.py), specialised to the simpler 1-distance setting.

The classical Lemmens-Seidel theorem (1973) establishes uniqueness of
6 equiangular lines in R^3 with angle arccos(1/sqrt 5) -- this script
gives an independent, finite-exhaustive verification by enumerating
signed Gram matrices and quotienting by switching + S_6.

Expected outcome: exactly 1 equivalence class -- the icosahedral
vertex-pair Gram matrix  H_{ico6}.
"""

from __future__ import annotations

import itertools
import math

import numpy as np
import sympy as sp

n = 6
sqrt5 = sp.sqrt(5)
mag = 1 / sqrt5
absH = sp.eye(n)
for i in range(n):
    for j in range(i + 1, n):
        absH[i, j] = absH[j, i] = mag


# ----------------------------------------------------------------------
#   Aut(K_6) = S_6 (order 720)
# ----------------------------------------------------------------------
auts = list(itertools.permutations(range(n)))
assert len(auts) == math.factorial(n) == 720


def permute_matrix(H, p):
    Hp = sp.zeros(n)
    for i in range(n):
        for j in range(n):
            Hp[p[i], p[j]] = H[i, j]
    return Hp


def switch_first_row_positive(H):
    """Switching gauge: enforce H[0, j] > 0 for j = 1..n-1."""
    D = [1] * n
    for j in range(1, n):
        D[j] = 1 if float(H[0, j].evalf()) >= 0 else -1
    H2 = sp.zeros(n)
    for i in range(n):
        for j in range(n):
            H2[i, j] = D[i] * H[i, j] * D[j]
    return H2


def _sign_vector(H):
    """Extract the 15-entry off-diagonal sign vector from a 6x6 H whose
    entries are exactly +/- 1/sqrt(5) off-diagonal and 1 on-diagonal.
    Returns a tuple in {-1, +1}^15.  Uses numeric evaluation, NOT symbolic
    simplify -- much faster.
    """
    out = []
    for i in range(n):
        for j in range(i + 1, n):
            v = float(H[i, j].evalf())
            out.append(1 if v > 0 else -1)
    return tuple(out)


def flat_key(H):
    return _sign_vector(H)


def canonical_key(H):
    """Smallest sign-vector signature over Aut(K_6) = S_6 orbit
    (modulo switching gauge that pins H[0, j] > 0)."""
    s = _sign_vector(H)                       # numpy-friendly form for permute
    # Build pair-index map: (i, j) i<j -> position in sign vector
    pair_to_pos = {}
    pos = 0
    for i in range(n):
        for j in range(i + 1, n):
            pair_to_pos[(i, j)] = pos
            pos += 1

    best = None
    for p in auts:
        # Permute sign vector by p: new pair (a, b) = sorted(p[i], p[j])
        # carries the sign of original (i, j).
        permuted = [0] * len(s)
        for i in range(n):
            for j in range(i + 1, n):
                a, b = p[i], p[j]
                if a > b:
                    a, b = b, a
                permuted[pair_to_pos[(a, b)]] = s[pair_to_pos[(i, j)]]
        # Apply switching gauge: vertex 0 row's signs all = +1.
        # Switching flips signs of all edges incident to vertices with eps=-1.
        eps = [1] * n
        for j in range(1, n):
            if permuted[pair_to_pos[(0, j)]] < 0:
                eps[j] = -1
        switched = [0] * len(permuted)
        for i in range(n):
            for j in range(i + 1, n):
                switched[pair_to_pos[(i, j)]] = \
                    eps[i] * eps[j] * permuted[pair_to_pos[(i, j)]]
        sig = tuple(switched)
        if best is None or sig < best:
            best = sig
    return best


# ----------------------------------------------------------------------
#   Schur-complement enumeration with 3 anchors  S = {0, 1, 2}.
# ----------------------------------------------------------------------
anchors = [0, 1, 2]
rest = [i for i in range(n) if i not in anchors]


def anchor_gram(anchor_signs):
    G = sp.eye(3)
    for s, (r, c) in zip(anchor_signs, [(0, 1), (0, 2), (1, 2)]):
        i, j = anchors[r], anchors[c]
        G[r, c] = G[c, r] = s * absH[i, j]
    return sp.simplify(G)


def x_vector(j, s1, s2):
    return sp.Matrix([
        absH[anchors[0], j],
        s1 * absH[anchors[1], j],
        s2 * absH[anchors[2], j],
    ])


solutions = []
NUM_TOL = 1e-9


# Numerical version of absH^2
absH2_num = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        if i != j:
            absH2_num[i, j] = 1 / 5.0


for anchor_signs in itertools.product([1, -1], repeat=3):
    G3 = anchor_gram(anchor_signs)
    G3_num = np.array(G3.evalf(), dtype=float)
    if abs(np.linalg.det(G3_num)) < 1e-12:
        continue
    ev = np.linalg.eigvalsh(G3_num)
    if ev[0] <= 1e-9:
        continue

    G3inv_num = np.linalg.inv(G3_num)

    # Build candidate x-vectors numerically for each non-anchor j.
    candidates_num: dict[int, list[np.ndarray]] = {}
    candidates_sym: dict[int, list[sp.Matrix]] = {}
    for j in rest:
        lst_num = []
        lst_sym = []
        for s1, s2 in itertools.product([1, -1], repeat=2):
            x_sym = x_vector(j, s1, s2)
            x_num = np.array(x_sym.evalf(), dtype=float).flatten()
            norm = float(x_num @ G3inv_num @ x_num)
            if abs(norm - 1) < NUM_TOL:
                lst_num.append(x_num)
                lst_sym.append(x_sym)
        candidates_num[j] = lst_num
        candidates_sym[j] = lst_sym

    if any(len(candidates_num[j]) == 0 for j in rest):
        continue

    order = sorted(rest, key=lambda j: len(candidates_num[j]))
    xs_num: dict[int, np.ndarray] = {}
    xs_sym: dict[int, sp.Matrix] = {}

    def backtrack(pos):
        if pos == len(order):
            # Build H numerically first; do symbolic check only on survivors.
            H_num = np.eye(n)
            for s, (r, c) in zip(anchor_signs, [(0, 1), (0, 2), (1, 2)]):
                i, j = anchors[r], anchors[c]
                H_num[i, j] = H_num[j, i] = s * (1 / math.sqrt(5))
            for j in rest:
                x = xs_num[j]
                for r, i in enumerate(anchors):
                    H_num[i, j] = H_num[j, i] = x[r]
            for j1, j2 in itertools.combinations(rest, 2):
                val = float(xs_num[j1] @ G3inv_num @ xs_num[j2])
                H_num[j1, j2] = H_num[j2, j1] = val
            # Magnitude check
            for i in range(n):
                for jj in range(n):
                    if i == jj:
                        continue
                    if abs(H_num[i, jj] ** 2 - 1 / 5.0) > NUM_TOL:
                        return
            evals = np.linalg.eigvalsh(H_num)
            if evals[0] < -1e-8:
                return
            if np.linalg.matrix_rank(H_num, tol=1e-8) != 3:
                return
            # Build the symbolic H for canonicalisation / display.
            H = sp.eye(n)
            for s, (r, c) in zip(anchor_signs, [(0, 1), (0, 2), (1, 2)]):
                i, j = anchors[r], anchors[c]
                H[i, j] = H[j, i] = s * absH[i, j]
            for j in rest:
                x_sym = xs_sym[j]
                for r, i in enumerate(anchors):
                    H[i, j] = H[j, i] = x_sym[r]
            for j1, j2 in itertools.combinations(rest, 2):
                # Snap the numerical cross product to its exact +/- 1/sqrt(5) form.
                v = float(xs_num[j1] @ G3inv_num @ xs_num[j2])
                s = 1 if v > 0 else -1
                H[j1, j2] = H[j2, j1] = s * absH[anchors[0], rest[0]]  # = 1/sqrt(5)
            solutions.append(H)
            return

        j = order[pos]
        for x_num, x_sym in zip(candidates_num[j], candidates_sym[j]):
            ok = True
            for k, y_num in xs_num.items():
                val = float(x_num @ G3inv_num @ y_num)
                if abs(val ** 2 - 1 / 5.0) > NUM_TOL:
                    ok = False
                    break
            if not ok:
                continue
            xs_num[j] = x_num
            xs_sym[j] = x_sym
            backtrack(pos + 1)
            del xs_num[j]
            del xs_sym[j]

    backtrack(0)


classes = {}
for H in solutions:
    classes.setdefault(canonical_key(H), H)

print("=" * 70)
print("SGR Lemma for 6-vertex equiangular orbit (single-distance K_6)")
print("=" * 70)
print(f"n = {n}, |H_ij| = 1/sqrt(5) for i != j, target rank = 3, PSD")
print(f"Aut(K_6) = S_6, order {len(auts)}")
print()
print(f"raw rank-3 PSD signed lifts:                {len(solutions)}")
print(f"classes modulo switching + S_6:             {len(classes)}")

rep = next(iter(classes.values()))
rep_num = np.array(rep.evalf(), dtype=float)
evals = np.linalg.eigvalsh(rep_num)
print(f"representative numerical eigenvalues:       "
      f"{np.round(evals, 10).tolist()}")
print(f"representative rank:                        "
      f"{np.linalg.matrix_rank(rep_num, tol=1e-8)}")
print(f"representative is PSD:                      "
      f"{bool(evals[0] > -1e-8)}")
print()
print("Representative H (icosahedral 6-vertex Gram):")
print(rep)
print()

# Seidel-matrix view: S = sqrt(5) * (H - I), off-diag in {-1, +1}
S = sqrt5 * (rep - sp.eye(n))
S_num = np.array(S.evalf(), dtype=float)
S_evals = sorted(np.round(np.linalg.eigvalsh(S_num), 8).tolist(), reverse=True)
print(f"Seidel matrix S = sqrt(5) * (H - I) spectrum: {S_evals}")
print(f"  (classical regular two-graph: {{sqrt(5)^3, -sqrt(5)^3}})")

# Tight-frame check: H decomposes as VV^T with V (6x3); V^T V should be 2*I_3
print()
w, U = np.linalg.eigh(rep_num)
idx = np.argsort(w)[::-1][:3]
V = (U[:, idx] * np.sqrt(w[idx])).reshape(n, 3)
VtV = V.T @ V
print(f"V^T V (should be 2 * I_3 since trace(H) = 6, rank 3):")
print(np.round(VtV, 8))
