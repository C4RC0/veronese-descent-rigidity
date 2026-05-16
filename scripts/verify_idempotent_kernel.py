"""
Verify the two structural claims that the literature search highlighted as
the spine of the clean proof:

   (1)  K := (3/2)(Q - J/3)  equals  2 * E_1
        where E_1 is the primitive idempotent of the Petersen association
        scheme corresponding to the eigenvalue +1 of the Petersen adjacency.

   (2)  ker K = span of the five 4-coclique indicator vectors of Petersen.
        (Petersen has exactly five 4-cocliques, one per element of {1,...,5}
        in the Kneser K(5,2) realisation.)

Both claims are immediate from association-scheme theory once stated; this
script numerically checks them and prints the explicit witness.
"""

from __future__ import annotations

import itertools
import math

import numpy as np


# ----------------------------------------------------------------------
#   Petersen graph in the Kneser K(5,2) labelling.
# ----------------------------------------------------------------------
LABELS = list(itertools.combinations(range(5), 2))
n = len(LABELS)
assert n == 10

A = np.zeros((n, n), dtype=float)
for i in range(n):
    for j in range(n):
        if i == j:
            continue
        # edge iff 2-subsets are disjoint
        if set(LABELS[i]).isdisjoint(LABELS[j]):
            A[i, j] = 1.0

# Verify Petersen invariants
deg = A.sum(axis=1)
assert all(d == 3 for d in deg), "Petersen should be 3-regular"
spec_A = np.sort(np.linalg.eigvalsh(A))[::-1]
print(f"Petersen spectrum (descending): {np.round(spec_A, 6).tolist()}")
print(f"  expected: [3, 1, 1, 1, 1, 1, -2, -2, -2, -2]")

I = np.eye(n)
J = np.ones((n, n))


# ----------------------------------------------------------------------
#   Q = (8/9) I + (4/9) A + (1/9) J   (squared Gram from Petersen)
#   K = (3/2)(Q - J/3) = (4/3) I + (2/3) A - (1/3) J
# ----------------------------------------------------------------------
Q = (8 / 9) * I + (4 / 9) * A + (1 / 9) * J
K = (3 / 2) * (Q - J / 3)
K_alt = (4 / 3) * I + (2 / 3) * A - (1 / 3) * J
print()
print(f"K vs K_alt (closed form) match: {np.allclose(K, K_alt)}")


# ----------------------------------------------------------------------
#   Primitive idempotent E_1 of the Petersen association scheme.
#   The Petersen scheme has three primitive idempotents:
#     E_0 = J/n                              (eigenvalue 3 of A, mult 1)
#     E_1                                    (eigenvalue 1 of A, mult 5)
#     E_2                                    (eigenvalue -2 of A, mult 4)
#   Explicit formulas (well-known):
#     E_1 = (1/n) * sum_{lambda in {1,1,1,1,1}} v v^T  for eigenvectors
#   Easier: build E_1 from the spectral projector.
# ----------------------------------------------------------------------
w, V = np.linalg.eigh(A)
# Find eigenvectors for eigenvalue +1 (multiplicity 5)
mask1 = np.abs(w - 1.0) < 1e-9
V1 = V[:, mask1]
E1 = V1 @ V1.T
print()
print(f"trace(E_1) = {np.trace(E1):.6f} (should equal mult = 5)")
print(f"E_1 is idempotent: ||E_1^2 - E_1|| = {np.linalg.norm(E1 @ E1 - E1):.2e}")

print()
print(f"||K - 2 * E_1|| = {np.linalg.norm(K - 2 * E1):.2e}")
print(f"Claim 1: K == 2 * E_1: {np.allclose(K, 2 * E1)}")


# ----------------------------------------------------------------------
#   Find the five 4-cocliques of Petersen.
#   A 4-coclique = 4 vertices pairwise non-adjacent = 4 2-subsets pairwise
#   intersecting (sharing 1 element).  In the Kneser model these are the
#   "stars" around a common element x in {0,...,4}.
# ----------------------------------------------------------------------
star_indicators = []
for x in range(5):
    indicator = np.zeros(n)
    for i, lab in enumerate(LABELS):
        if x in lab:
            indicator[i] = 1
    # check it's an independent set (no edges)
    rows = [i for i in range(n) if indicator[i] == 1]
    is_coclique = all(A[r1, r2] == 0 for r1 in rows for r2 in rows if r1 != r2)
    assert is_coclique, f"star at {x} is not a coclique"
    assert int(indicator.sum()) == 4, f"star at {x} has size != 4"
    star_indicators.append(indicator)

S = np.stack(star_indicators).T   # shape (10, 5): columns are indicator vectors
print()
print(f"Found {S.shape[1]} 4-cocliques (expected 5)")
print(f"sum of 5 indicator vectors: {S.sum(axis=1).tolist()}")
print(f"  (each vertex sits in exactly 2 of the 5 stars: 2-subset has 2 elements)")


# ----------------------------------------------------------------------
#   Check ker K = span(S).
# ----------------------------------------------------------------------
# (a) Each indicator should be in ker K:
for k in range(5):
    Kx = K @ S[:, k]
    print(f"  K * 1_{{star {k}}} = {np.round(Kx, 8).tolist()}")
KS = K @ S
print()
print(f"||K * S|| = {np.linalg.norm(KS):.2e}")
print(f"All indicators in ker K: {np.linalg.norm(KS) < 1e-8}")

# (b) Dimension match: rank(S) should equal nullity(K).
rk_S = np.linalg.matrix_rank(S, tol=1e-8)
rk_K = np.linalg.matrix_rank(K, tol=1e-8)
nullity_K = n - rk_K
print()
print(f"rank(S) = {rk_S},  rank(K) = {rk_K},  nullity(K) = {nullity_K}")
print(f"Claim 2 (dim match): {rk_S == nullity_K}")

# (c) Combined check: K @ S == 0 AND rank(S) == nullity(K) ⇒ ker K = col(S).
claim2 = (np.linalg.norm(K @ S) < 1e-8) and (rk_S == nullity_K)
print(f"Claim 2: ker K = span of 5 4-coclique indicators: {claim2}")


# ----------------------------------------------------------------------
#   Verify the immediate geometric consequence:
#   sum_{i in C_a} U_i = 0  for each 4-coclique C_a.
# ----------------------------------------------------------------------
# Build the 10 icosahedral face-pair vectors (V_ico) and the U_i lift.
import math as _m

PHI = (1.0 + _m.sqrt(5.0)) / 2.0
ICO_VERTS = np.array([
    [0,  1,  PHI], [0, -1, -PHI], [0, -1,  PHI], [0,  1, -PHI],
    [1,  PHI, 0], [-1, -PHI, 0], [-1,  PHI, 0], [1, -PHI, 0],
    [PHI, 0,  1], [-PHI, 0, -1], [-PHI, 0,  1], [PHI, 0, -1],
], dtype=float)

edge_len = 2.0
faces = []
for i, j, k in itertools.combinations(range(12), 3):
    d_ij = np.linalg.norm(ICO_VERTS[i] - ICO_VERTS[j])
    d_jk = np.linalg.norm(ICO_VERTS[j] - ICO_VERTS[k])
    d_ik = np.linalg.norm(ICO_VERTS[i] - ICO_VERTS[k])
    if (abs(d_ij - edge_len) < 1e-6
            and abs(d_jk - edge_len) < 1e-6
            and abs(d_ik - edge_len) < 1e-6):
        faces.append((i, j, k))
centers = np.array([ICO_VERTS[list(f)].mean(axis=0) for f in faces])

matched = [False] * 20
pair_reps = []
for i in range(20):
    if matched[i]:
        continue
    for j in range(i + 1, 20):
        if matched[j]:
            continue
        if np.linalg.norm(centers[i] + centers[j]) < 1e-6:
            pair_reps.append(centers[i])
            matched[i] = matched[j] = True
            break
V_ico = np.array(pair_reps)
V_ico = V_ico / np.linalg.norm(V_ico, axis=1, keepdims=True)

# Identify the Petersen labelling of V_ico (which face-pair sits at which
# 2-subset).  This is via the cube-inscription: each face-pair is a body
# diagonal of exactly 2 of the 5 inscribed cubes; the 2-subset names these.
# But for our verification we just need the existence of SOME labelling that
# matches the K = 2 E_1 condition.  Use the Gram-derived adjacency directly.
H_ico = V_ico @ V_ico.T
A_ico = np.zeros((10, 10), dtype=float)
for i in range(10):
    for j in range(10):
        if i != j and abs(abs(H_ico[i, j]) - _m.sqrt(5) / 3) < 1e-6:
            A_ico[i, j] = 1.0
# Permutation matching A_ico to A: brute force isomorphism (10! too big,
# use simple BFS-canonicalisation by signature).
# Simpler: A_ico and A are both Petersen, so they're isomorphic; we find a
# permutation P such that P A P^T = A_ico, by enumerating permutations
# constrained by the degree sequence + adjacency.
def find_isomorphism(A1, A2):
    """Brute-force permutation pi with A1[pi[i], pi[j]] == A2[i, j]."""
    n = A1.shape[0]
    # backtrack via BFS-tree match
    pi = [-1] * n
    used = [False] * n

    def backtrack(k):
        if k == n:
            return True
        for c in range(n):
            if used[c]:
                continue
            # check edges to previously placed vertices
            ok = True
            for prev in range(k):
                if A1[c, pi[prev]] != A2[k, prev]:
                    ok = False
                    break
            if not ok:
                continue
            pi[k] = c
            used[c] = True
            if backtrack(k + 1):
                return True
            used[c] = False
            pi[k] = -1
        return False
    backtrack(0)
    return pi

pi = find_isomorphism(A_ico, A)   # pi: index in A's labelling -> index in A_ico
P = np.zeros((10, 10), dtype=float)
for i, p in enumerate(pi):
    P[p, i] = 1.0
# Sanity: P A P^T should equal A_ico
assert np.allclose(P @ A @ P.T, A_ico), "isomorphism failed"

# Now reorder V_ico to match A's labelling
V_reordered = P.T @ V_ico          # shape (10, 3)

# Compute U_i (raw, unnormalised)
U = np.zeros((10, 3, 3))
for i in range(10):
    v = V_reordered[i]
    U[i] = np.outer(v, v) - np.eye(3) / 3

print()
print("Geometric consequence: for each 4-coclique C_a (star around x in {0..4}),")
print("    sum_{i in C_a} U_i  should equal  0  (zero matrix).")
print()
for x in range(5):
    rows = [i for i in range(10) if x in LABELS[i]]
    S_x = U[rows].sum(axis=0)
    err = np.linalg.norm(S_x)
    print(f"  star around x={x}: indices {rows}")
    print(f"    ||sum_C U_i|| = {err:.2e}")
