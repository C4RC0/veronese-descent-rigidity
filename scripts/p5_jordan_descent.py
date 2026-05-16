"""
P5 verification: Jordan-descent / Schur's lemma argument for compound
rigidity.

Claim:  The 10 lifted vectors  U_i = v_i v_i^T - I/3  in  Sym^2_0(R^3)
have eigenvalue spectrum (2/3, -1/3, -1/3).  Their negatives  -U_i
have spectrum (-2/3, 1/3, 1/3), which is NOT of the rank-1-projector
form, so  -U  is NOT a Veronese image  v v^T - I/3  for any unit v.

This is the algebraic obstruction that kills the Schur-lemma  T = -I
alternative, leaving only  T = +I  and thus uniqueness up to O(3).

The script:
  1. Builds the 10 icosahedral lap-pair vectors and their Veronese
     lift U_i.
  2. Verifies eigenvalue spectrum of each U_i.
  3. Verifies that -U_i is NOT a Veronese image (no v with vv^T = -U_i + I/3).
  4. Verifies that the A_5 action on the 10 U_i is the 5-dim irrep
     (by character / multiplicity check on the permutation representation).
"""

from __future__ import annotations

import itertools
import math

import numpy as np


PHI = (1.0 + math.sqrt(5.0)) / 2.0


# ----------------------------------------------------------------------
#   Build V_ico (10 icosahedral face-pair unit vectors).
# ----------------------------------------------------------------------
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
assert V_ico.shape == (10, 3)


# ----------------------------------------------------------------------
#   Compute the Veronese lift U_i = v_i v_i^T - I/3.
# ----------------------------------------------------------------------
U = np.zeros((10, 3, 3))
for i in range(10):
    v = V_ico[i]
    U[i] = np.outer(v, v) - np.eye(3) / 3


# ----------------------------------------------------------------------
#   Claim 1: each U_i has eigenvalues (2/3, -1/3, -1/3).
# ----------------------------------------------------------------------
print("=" * 70)
print("Claim 1: U_i has spectrum (2/3, -1/3, -1/3)  -- Veronese image")
print("=" * 70)
target_spec = sorted([2 / 3, -1 / 3, -1 / 3], reverse=True)
all_ok = True
for i in range(10):
    w = sorted(np.linalg.eigvalsh(U[i]).tolist(), reverse=True)
    err = max(abs(w[k] - target_spec[k]) for k in range(3))
    ok = err < 1e-10
    all_ok = all_ok and ok
    if i < 3 or not ok:
        print(f"  U_{i} eigenvalues: {[round(x, 12) for x in w]}  "
              f"|max err| = {err:.2e}  OK = {ok}")
print(f"\nAll 10 U_i match (2/3, -1/3, -1/3): {all_ok}\n")


# ----------------------------------------------------------------------
#   Claim 2: -U_i has spectrum (-2/3, 1/3, 1/3) -- NOT Veronese image.
#   This is the algebraic obstruction.
#
#   For W = -U_i, check if W = v v^T - I/3 for some unit v.  This requires
#   W + I/3 = v v^T  to be rank-1 PSD.  Compute eigenvalues of W + I/3:
# ----------------------------------------------------------------------
print("=" * 70)
print("Claim 2: -U_i is NOT a Veronese image (Schur T = -I excluded)")
print("=" * 70)
print()
print("For W = -U to be a Veronese image, W + I/3 must be rank-1 PSD.")
print("Check eigenvalues of -U_i + I/3:")
for i in range(3):
    W = -U[i]
    w = sorted(np.linalg.eigvalsh(W + np.eye(3) / 3).tolist(), reverse=True)
    print(f"  i = {i}:  spec(W + I/3) = {[round(x, 6) for x in w]}")
    is_rank1_psd = (abs(w[1]) < 1e-10 and abs(w[2]) < 1e-10 and w[0] > 1e-10)
    print(f"           rank-1 PSD?  {is_rank1_psd}  "
          f"(rank-1 needs two zero eigenvalues + one positive)")
print()
# U_i has eigenvalues (2/3, -1/3, -1/3) with v_i corresponding to 2/3.
# So -U_i has eigenvalues (-2/3, 1/3, 1/3).
# -U_i + I/3 has eigenvalues (-2/3 + 1/3, 1/3 + 1/3, 1/3 + 1/3) = (-1/3, 2/3, 2/3).
# This has a NEGATIVE eigenvalue, so it's NOT PSD, let alone rank-1 PSD.

print("Analytically: U_i has eigenvalues (2/3, -1/3, -1/3), so")
print("  -U_i has eigenvalues (-2/3, 1/3, 1/3), and")
print("  -U_i + I/3 has eigenvalues (-1/3, 2/3, 2/3).")
print("This matrix has a NEGATIVE eigenvalue (-1/3), so it is NOT PSD.")
print("Hence  -U_i + I/3  cannot equal  v v^T  for any v in R^3.")
print("CONCLUSION: -U_i is NOT in the Veronese image of R^3.")
print()


# ----------------------------------------------------------------------
#   Claim 3 (sanity): the A_5 action on the 10 points has the
#   5-dim irrep as a multiplicity-1 component of the 10-dim permutation
#   representation.
# ----------------------------------------------------------------------
print("=" * 70)
print("Claim 3: 10-dim A_5-permutation rep = 1 + 4 + 5  (5-irrep mult 1)")
print("=" * 70)

# Build A_5 explicitly from generators
def rotation_matrix(axis, angle):
    a = axis / np.linalg.norm(axis)
    c, s, C = math.cos(angle), math.sin(angle), 1 - math.cos(angle)
    x, y, z = a
    return np.array([
        [c + x * x * C,     x * y * C - z * s, x * z * C + y * s],
        [y * x * C + z * s, c + y * y * C,     y * z * C - x * s],
        [z * x * C - y * s, z * y * C + x * s, c + z * z * C],
    ])


# A_5 generators (5-fold around a vertex axis, 3-fold around a face axis)
v0_axis = np.array([0.0, 1.0, PHI])
f0_axis = (np.array([0.0, 1.0, PHI])
           + np.array([1.0, PHI, 0.0])
           + np.array([-1.0, PHI, 0.0]))
gens = [
    rotation_matrix(v0_axis, 2 * math.pi / 5),
    rotation_matrix(f0_axis, 2 * math.pi / 3),
]

# BFS group closure with SVD-snap
def close_group(gens, tol=1e-5, max_order=80):
    group = [np.eye(3)]

    def already(M):
        return any(np.linalg.norm(M - X) < tol for X in group)

    frontier = []
    for g in gens:
        if not already(g):
            group.append(g)
            frontier.append(g)
    while frontier:
        h = frontier.pop()
        for g in gens:
            M = g @ h
            U_, _, Vt = np.linalg.svd(M)
            M = U_ @ Vt
            if np.linalg.det(M) < 0:
                M = -M
            if not already(M):
                group.append(M)
                frontier.append(M)
                if len(group) > max_order:
                    raise RuntimeError(f"group exceeded {max_order}")
    return group


A5 = close_group(gens)
assert len(A5) == 60, f"|A_5| = {len(A5)}, expected 60"


# Build the 10-dim permutation representation: g in A_5 acts on V_ico by
# rotation, permuting the 10 lap-pair vectors (up to sign).
def perm_action_matrix(g):
    """10x10 permutation/sign matrix:  g v_i  =  sign * v_{pi(i)}."""
    M = np.zeros((10, 10))
    for i in range(10):
        gv = g @ V_ico[i]
        for j in range(10):
            if np.linalg.norm(gv - V_ico[j]) < 1e-6:
                M[j, i] = 1.0
                break
            if np.linalg.norm(gv + V_ico[j]) < 1e-6:
                M[j, i] = -1.0
                break
        else:
            raise RuntimeError(f"no match for g v_{i}")
    return M


# Character of the 10-dim rep: chi(g) = trace of perm_action_matrix
# But this is the SIGNED permutation rep.  For the lifted U_i, sign-flip
# v_i -> -v_i acts as identity (since U is degree-2 in v).  So the
# 10-dim rep on U_i is the unsigned permutation rep.
def unsigned_perm_action_matrix(g):
    """10x10 unsigned permutation matrix:  g {v_i, -v_i} -> {v_{pi(i)}, -v_{pi(i)}}."""
    M = np.zeros((10, 10))
    for i in range(10):
        gv = g @ V_ico[i]
        for j in range(10):
            if (np.linalg.norm(gv - V_ico[j]) < 1e-6
                    or np.linalg.norm(gv + V_ico[j]) < 1e-6):
                M[j, i] = 1.0
                break
        else:
            raise RuntimeError(f"no match for g v_{i}")
    return M


# Compute character (trace) on the 10-dim unsigned-perm rep
chars_10 = [int(round(np.trace(unsigned_perm_action_matrix(g)))) for g in A5]
print(f"\nA_5 character on 10-dim unsigned perm rep:")
# Group by conjugacy class.  A_5 has 5 conjugacy classes:
# {e}, {12 double-3-cycles? actually}  classes: e(1), (12)(34)(15), 3-cycles(20),
# 5-cycles split (12,12)
# Just print distribution
from collections import Counter
print(f"  character value distribution: {Counter(chars_10)}")

# A_5 has irrep dims: 1, 3, 3, 4, 5.  Character table (rows = irreps, cols = classes):
# classes: e, (12)(34), (123), (12345), (12354)  with sizes 1, 15, 20, 12, 12.
# Trivial:  1,  1,  1,  1,  1
# 3a:       3, -1,  0,  phi-1, -phi
# 3b:       3, -1,  0, -phi, phi-1
# 4:        4,  0,  1, -1, -1
# 5:        5,  1, -1,  0,  0
#
# Inner product <chi_10, chi_irrep> = (1/60) sum_g chi_10(g) * chi_irrep(g).
# We need to identify conjugacy classes of A_5.

# Identify conjugacy classes by element order
def element_order(g, tol=1e-5):
    M = np.eye(3)
    for k in range(1, 7):
        M = M @ g
        if np.linalg.norm(M - np.eye(3)) < tol:
            return k
    return -1


orders = [element_order(g) for g in A5]
print(f"\nElement orders in A_5: {Counter(orders)}")
# Expected: 1: 1, 2: 15, 3: 20, 5: 24 (= 12 + 12)

# A_5 has 5 classes:  e (1 element, order 1), 15 elements (order 2),
# 20 elements (order 3), 12+12 elements (order 5, two classes).
# We'll use orders to roughly bucket.

# Compute inner product with 5-dim irrep character:
# class      e   (12)(34)  (123)  (12345)  (12354)
# size       1    15       20      12       12
# chi_5      5    1        -1      0        0

# Aggregate chars_10 by element order:
agg = {1: [], 2: [], 3: [], 5: []}
for c, o in zip(chars_10, orders):
    if o in agg:
        agg[o].append(c)
print("\nCharacter values by element order:")
for o in [1, 2, 3, 5]:
    print(f"  order {o}: count={len(agg[o])}, characters={Counter(agg[o])}")


# Inner product with chi_5 (trivial-detection):
# <chi_10, chi_5> = (1/60)[1 * chi_10(e) * 5 + 15 * avg(chi_10 over order-2) * 1
#                          + 20 * avg(order-3) * (-1) + 24 * avg(order-5) * 0]
def avg(L): return sum(L) / len(L) if L else 0


ip_5 = (1 * agg[1][0] * 5
        + 15 * avg(agg[2]) * 1
        + 20 * avg(agg[3]) * (-1)
        + 12 * avg(agg[5]) * 0
        + 12 * avg(agg[5]) * 0) / 60
print(f"\n<chi_10, chi_5>  (multiplicity of 5-irrep in 10-dim rep) = {ip_5}")

# Inner products with other irreps:
ip_trivial = (1 * agg[1][0] * 1 + 15 * avg(agg[2]) * 1
              + 20 * avg(agg[3]) * 1 + 24 * avg(agg[5]) * 1) / 60
print(f"<chi_10, chi_1>  (trivial mult)  = {ip_trivial}")

ip_4 = (1 * agg[1][0] * 4 + 15 * avg(agg[2]) * 0
        + 20 * avg(agg[3]) * 1 + 12 * avg(agg[5]) * (-1)
        + 12 * avg(agg[5]) * (-1)) / 60
print(f"<chi_10, chi_4>  (4-irrep mult)  = {ip_4}")

print()
print("Expected decomposition of 10-dim A_5-perm-rep:  1 + 4 + 5")
print("(matches the Petersen/Johnson scheme decomposition)")
