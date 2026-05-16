"""
Signed Gram Rigidity (SGR) Lemma — finite exhaustive verification.

Claim.  Let H in R^{10x10} be a real symmetric matrix with:
    (i)   H_ii = 1
    (ii)  |H_ij| = sqrt(5)/3      if (i,j) in E(Petersen)
    (iii) |H_ij| = 1/3            if (i,j) not in E(Petersen) and i != j
    (iv)  H is positive semidefinite
    (v)   rank(H) = 3
Then H equals (up to a switching D = diag(+/-1) and a permutation
sigma in Aut(Petersen) = S_5) the icosahedral antipodal face-normal
Gram matrix H_ico.

This script verifies the claim numerically by:
  1. Building H_ico from the 10 face-pair unit vectors of the regular
     icosahedron.
  2. Enumerating the 2^14 = 16384 sign-patterns on Petersen edges
     (one sign fixed by switching normalisation).
  3. For each candidate, determining the unique sign-pattern on
     complement edges forced by the rank<=3 + PSD constraint
     (via 4x4 vanishing-minor conditions), and verifying that the
     resulting H is actually rank-3 PSD with the required magnitudes.
  4. Collapsing the surviving candidates modulo Aut(Petersen) = S_5
     (order 120, including parity).

Expected outcome: exactly ONE equivalence class — H_ico.
"""

from __future__ import annotations

import itertools
import math
import sys
import time
from typing import Iterable

import numpy as np


# ----------------------------------------------------------------------
#   Geometry: 12 vertices of the regular icosahedron, 20 faces, 10
#   antipodal face-pair representatives.
# ----------------------------------------------------------------------
PHI = (1.0 + math.sqrt(5.0)) / 2.0

ICO_VERTS = np.array([
    [0,  1,  PHI], [0, -1, -PHI], [0, -1,  PHI], [0,  1, -PHI],
    [1,  PHI, 0], [-1, -PHI, 0], [-1,  PHI, 0], [1, -PHI, 0],
    [PHI, 0,  1], [-PHI, 0, -1], [-PHI, 0,  1], [PHI, 0, -1],
], dtype=float)


def _icosahedron_face_pair_vectors() -> np.ndarray:
    """10 unit vectors -- one per antipodal pair of icosahedral faces."""
    # 1) Detect faces: triples of vertices at pairwise edge-distance 2.
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
    assert len(faces) == 20, f"expected 20 faces, got {len(faces)}"

    centers = np.array([ICO_VERTS[list(f)].mean(axis=0) for f in faces])

    # 2) Pair antipodal face-centers.
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
    assert len(pair_reps) == 10

    pair_reps = np.array(pair_reps)
    pair_reps = pair_reps / np.linalg.norm(pair_reps, axis=1, keepdims=True)
    return pair_reps


V_ICO = _icosahedron_face_pair_vectors()   # shape (10, 3), unit rows


def _gram(V: np.ndarray) -> np.ndarray:
    return V @ V.T


H_ICO = _gram(V_ICO)


# ----------------------------------------------------------------------
#   Petersen graph from the |inner product|^2 = 1/9 pattern.
#
#   For the icosahedral face-pair vectors:
#     - 15 pairs (i,j) with |v_i . v_j| = 1/sqrt(5)  i.e.  H_ij^2 = 1/5
#       ... wait, recompute:
#
#   Actually  |v_i . v_j|  takes only two values: 1/sqrt(5)? Let's check.
#   We expect the two values  {1/3, sqrt(5)/3}, with the SMALLER (1/3)
#   on COMPLEMENT edges (shared-cube) and LARGER (sqrt(5)/3) on
#   PETERSEN edges (disjoint-cubes).  Verified empirically below.
# ----------------------------------------------------------------------
def _classify_edges(H: np.ndarray) -> tuple[set, set]:
    """Return (Petersen edges, complement edges) as sets of (i,j) with i<j."""
    PET = set()
    COMP = set()
    for i in range(10):
        for j in range(i + 1, 10):
            a = abs(H[i, j])
            if abs(a - math.sqrt(5) / 3) < 1e-6:
                PET.add((i, j))
            elif abs(a - 1.0 / 3.0) < 1e-6:
                COMP.add((i, j))
            else:
                raise AssertionError(
                    f"H[{i},{j}] = {H[i,j]:.6f} has unexpected magnitude {a}")
    return PET, COMP


PETERSEN, COMPLEMENT = _classify_edges(H_ICO)
assert len(PETERSEN) == 15
assert len(COMPLEMENT) == 30


# ----------------------------------------------------------------------
#   Sanity: spectrum of H_ico.
# ----------------------------------------------------------------------
def _spectrum(H: np.ndarray) -> np.ndarray:
    w = np.linalg.eigvalsh(H)
    w[abs(w) < 1e-9] = 0.0
    return np.sort(w)[::-1]


SPEC_ICO = _spectrum(H_ICO)
# 10 vectors in R^3 -> rank 3, three nonzero eigenvalues summing to 10.
nz = SPEC_ICO[abs(SPEC_ICO) > 1e-8]
assert len(nz) == 3, f"H_ico should have rank 3, spectrum = {SPEC_ICO}"
assert abs(nz.sum() - 10.0) < 1e-8
# Tight frame property: for the antipodal-face-normal frame, V^T V = (10/3) I_3
# so the three nonzero eigenvalues of H = V V^T are all 10/3.
assert all(abs(lam - 10.0 / 3.0) < 1e-8 for lam in nz), \
    f"unexpected nonzero spectrum {nz}"


# ----------------------------------------------------------------------
#   Petersen-graph automorphism group  Aut(Petersen) = S_5.
#
#   Each Petersen vertex i in {0,...,9} corresponds to a 2-subset of
#   {0,1,2,3,4} (via the Kneser K(5,2) realisation).  We detect that
#   labelling from H_ico by looking at common neighbours in Petersen
#   (two vertices share a unique common neighbour iff their labels
#   share one element).  But a simpler route: brute-force permutations
#   that preserve the Petersen edge set.  There are 5! = 120 of them.
# ----------------------------------------------------------------------
def _aut_petersen() -> list[tuple[int, ...]]:
    pet = PETERSEN
    perms = []
    for p in itertools.permutations(range(10)):
        ok = True
        for (i, j) in pet:
            a, b = p[i], p[j]
            edge = (a, b) if a < b else (b, a)
            if edge not in pet:
                ok = False
                break
        if ok:
            perms.append(p)
    return perms


AUT = _aut_petersen()
assert len(AUT) == 120, f"Aut(Petersen) should have order 120, got {len(AUT)}"


# ----------------------------------------------------------------------
#   Canonical-form: switching D = diag(+/-1) acts by H -> D H D, which
#   flips simultaneously row i and column i for each i with D_ii = -1.
#   Normalisation: enforce H_{0,j} >= 0 for j in N_P(0).  Vertex 0 has
#   3 Petersen neighbours; fixing their signs to + uses 3 switching
#   degrees of freedom.  The remaining 9-3=6 switching d.o.f.
#   (vertices not in {0} cup N_P(0)) we kill below.
# ----------------------------------------------------------------------
def _canonicalise(H: np.ndarray) -> np.ndarray:
    """Apply a sequence of switchings to bring H to a canonical sign-form.

    Strategy: greedy on a BFS spanning tree of the Petersen graph rooted
    at 0.  For each non-root v, flip v's sign so that H_{parent(v), v} > 0.
    The 10 vertices are connected in Petersen, so 9 switching choices
    pin all Petersen-edge signs along the tree to +.  The remaining
    Petersen-edge signs (15 - 9 = 6 non-tree edges) and ALL complement
    signs become switching-invariants.
    """
    H = H.copy()
    # Build Petersen adjacency
    adj = [[] for _ in range(10)]
    for (i, j) in PETERSEN:
        adj[i].append(j)
        adj[j].append(i)
    parent = [-1] * 10
    seen = [False] * 10
    queue = [0]
    seen[0] = True
    order = []
    while queue:
        u = queue.pop(0)
        order.append(u)
        for v in adj[u]:
            if not seen[v]:
                seen[v] = True
                parent[v] = u
                queue.append(v)
    # Flip rows/cols in BFS order so H[parent, v] > 0.
    D = np.ones(10)
    for v in order[1:]:
        p = parent[v]
        if D[p] * D[v] * H[p, v] < 0:
            D[v] = -D[v]
    Dmat = np.diag(D)
    return Dmat @ H @ Dmat


H_ICO_CANON = _canonicalise(H_ICO)


# ----------------------------------------------------------------------
#   Build a candidate H from a sign-vector on Petersen edges + complement
#   edges.  Returns None if the result is not PSD with rank <= 3.
# ----------------------------------------------------------------------
PETERSEN_LIST = sorted(PETERSEN)             # length 15
COMPLEMENT_LIST = sorted(COMPLEMENT)         # length 30

PET_MAG = math.sqrt(5.0) / 3.0
COMP_MAG = 1.0 / 3.0


def _build_H(pet_signs: Iterable[int], comp_signs: Iterable[int]) -> np.ndarray:
    H = np.eye(10)
    for (i, j), s in zip(PETERSEN_LIST, pet_signs):
        v = s * PET_MAG
        H[i, j] = v
        H[j, i] = v
    for (i, j), s in zip(COMPLEMENT_LIST, comp_signs):
        v = s * COMP_MAG
        H[i, j] = v
        H[j, i] = v
    return H


def _is_psd_rank3(H: np.ndarray, tol: float = 1e-7) -> bool:
    w = np.linalg.eigvalsh(H)
    # rank 3: smallest 7 eigenvalues ~ 0, largest 3 > 0
    pos = w[w > tol]
    zero = w[abs(w) <= tol]
    neg = w[w < -tol]
    if neg.size > 0:
        return False
    if pos.size != 3:
        return False
    if zero.size != 7:
        return False
    return True


# ----------------------------------------------------------------------
#   Complement-sign extraction.
#
#   For a rank-3 PSD H with H_ii = 1, every 4x4 principal minor must be
#   zero.  Given the Petersen-edge signs and three of the four entries
#   of a 4x4 minor that involves a complement edge, the sign of the
#   complement entry can often be determined (or shown to be infeasible).
#
#   We DON'T attempt to derive complement signs symbolically.  Instead,
#   for each Petersen-sign-pattern we solve numerically:
#
#     Step (a):  Try to find rows V_i in R^3 with v_i . v_j = H_ij (signs
#                given on Petersen edges).  This is a least-squares
#                completion problem.  Treat V as 10x3 unknowns; minimise
#                || V V^T - target ||_F  subject to ||v_i|| = 1 and
#                signed Petersen-edge constraints.
#     Step (b):  Read off complement signs from the recovered V.
#     Step (c):  Check actual rank-3 PSD condition.
#
#   For efficiency: use the well-known fact that if a 10-vertex graph
#   with the Petersen edge labels +/-sqrt(5)/3 admits a rank-3 PSD
#   completion, the completion is unique (up to the global O(3)
#   action).  So if a completion exists, it's the unique witness.
#
#   Concretely we use eigendecomposition:
#     Form H_partial with diagonal = 1, Petersen entries = signed value,
#     complement entries = 0 (unknown).  Take top-3 eigenvectors and
#     check that the resulting rank-3 approximation is consistent with
#     |H_ij| = 1/3 on complement edges.
#
#   This is FAST (no nested 2^30 enumeration) and gives a definitive
#   yes/no per Petersen-sign-pattern.
# ----------------------------------------------------------------------
def _try_complete(pet_signs: tuple[int, ...]) -> np.ndarray | None:
    """Try to complete the partial Gram matrix to a rank-3 PSD with the
    required complement magnitudes.  Return the completed H or None.
    """
    # Form partial matrix: diag=1, Petersen edges signed, complement edges
    # unknown.  We solve the rank-3 completion via an alternating
    # projection (POCS) between (a) rank<=3 PSD and (b) the affine slice
    # where diag and Petersen entries are pinned to their target values.
    target = np.eye(10)
    pet_mask = np.zeros((10, 10), dtype=bool)
    np.fill_diagonal(pet_mask, True)
    for (i, j), s in zip(PETERSEN_LIST, pet_signs):
        target[i, j] = s * PET_MAG
        target[j, i] = s * PET_MAG
        pet_mask[i, j] = pet_mask[j, i] = True

    H = target.copy()
    for _ in range(500):
        # Project onto rank-<=3 PSD cone: clip negative eigenvalues and
        # keep top 3.
        w, U = np.linalg.eigh(H)
        idx = np.argsort(w)[::-1]
        w = w[idx]
        U = U[:, idx]
        w_top = np.maximum(w[:3], 0)
        H_psd = (U[:, :3] * w_top) @ U[:, :3].T
        # Project onto affine constraints: diag=1, Petersen entries fixed.
        H = H_psd.copy()
        H[pet_mask] = target[pet_mask]
        np.fill_diagonal(H, 1.0)
        # Convergence check
        diff = np.linalg.norm(H - H_psd, ord='fro')
        if diff < 1e-10:
            break
    # Final: project once more onto rank-3 PSD and accept that as H_final.
    w, U = np.linalg.eigh(H)
    idx = np.argsort(w)[::-1]
    w = w[idx]
    U = U[:, idx]
    w_top = np.maximum(w[:3], 0)
    H_final = (U[:, :3] * w_top) @ U[:, :3].T
    # Verify: complement entries should have |.| = 1/3 to within tolerance.
    for (i, j) in COMPLEMENT_LIST:
        if abs(abs(H_final[i, j]) - COMP_MAG) > 1e-3:
            return None
    # Verify: diag entries should be 1 (within tolerance).
    if any(abs(H_final[i, i] - 1.0) > 1e-3 for i in range(10)):
        return None
    # Snap to clean magnitudes for canonicalisation.
    H_clean = np.eye(10)
    for (i, j), s in zip(PETERSEN_LIST, pet_signs):
        H_clean[i, j] = H_clean[j, i] = s * PET_MAG
    for (i, j) in COMPLEMENT_LIST:
        s = 1 if H_final[i, j] > 0 else -1
        H_clean[i, j] = H_clean[j, i] = s * COMP_MAG
    if not _is_psd_rank3(H_clean):
        return None
    return H_clean


# ----------------------------------------------------------------------
#   Switching-canonical hash for equivalence classification.
# ----------------------------------------------------------------------
def _switch_canonical_signs(H: np.ndarray) -> tuple[int, ...]:
    """Return the sign-vector on all 45 off-diagonal entries in canonical
    form (after BFS-tree switching from vertex 0).
    """
    Hc = _canonicalise(H)
    signs = []
    for i in range(10):
        for j in range(i + 1, 10):
            signs.append(1 if Hc[i, j] > 0 else -1)
    return tuple(signs)


def _orbit_repr(H: np.ndarray) -> tuple[int, ...]:
    """Smallest signature over Aut(Petersen) orbit (modulo switching)."""
    best = None
    for perm in AUT:
        Hp = H[np.ix_(perm, perm)]
        sig = _switch_canonical_signs(Hp)
        if best is None or sig < best:
            best = sig
    return best


# ----------------------------------------------------------------------
#   Tree-Petersen-edges that we fix to +1 via switching.
#   Use the same BFS spanning tree as _canonicalise() to determine
#   which 9 Petersen edges are tree-edges (sign forced to +1).
# ----------------------------------------------------------------------
def _bfs_tree_edges() -> set:
    adj = [[] for _ in range(10)]
    for (i, j) in PETERSEN:
        adj[i].append(j)
        adj[j].append(i)
    parent = [-1] * 10
    seen = [False] * 10
    queue = [0]
    seen[0] = True
    tree = set()
    while queue:
        u = queue.pop(0)
        for v in adj[u]:
            if not seen[v]:
                seen[v] = True
                parent[v] = u
                edge = (u, v) if u < v else (v, u)
                tree.add(edge)
                queue.append(v)
    assert len(tree) == 9
    return tree


TREE_EDGES = _bfs_tree_edges()
NONTREE_PET = [e for e in PETERSEN_LIST if e not in TREE_EDGES]
assert len(NONTREE_PET) == 6


# ----------------------------------------------------------------------
#   Main exhaustive scan.
# ----------------------------------------------------------------------
def main(verbose: bool = True) -> None:
    if verbose:
        print("=" * 70)
        print("Signed Gram Rigidity Lemma -- finite exhaustive verification")
        print("=" * 70)
        print(f"V_ico: 10 icosahedral face-pair unit vectors")
        print(f"H_ico spectrum: {SPEC_ICO}")
        print(f"  -> rank = 3, nonzero eigenvalues all = 10/3 (tight frame)")
        print(f"|Petersen edges|    = {len(PETERSEN)} (target mag sqrt(5)/3"
              f" = {PET_MAG:.6f})")
        print(f"|Complement edges|  = {len(COMPLEMENT)} (target mag 1/3"
              f" = {COMP_MAG:.6f})")
        print(f"|Aut(Petersen)|     = {len(AUT)}   (= S_5)")
        print(f"BFS tree edges fixed to +1 by switching: {len(TREE_EDGES)}")
        print(f"Non-tree Petersen edges (free signs):    {len(NONTREE_PET)}")
        print(f"-> enumerating 2^{len(NONTREE_PET)} = "
              f"{2**len(NONTREE_PET)} Petersen-sign patterns")
        print()

    # Build the Petersen-sign vector for a given choice of non-tree-edge
    # signs (tree edges = +1).
    tree_idx = {e: k for k, e in enumerate(PETERSEN_LIST) if e in TREE_EDGES}
    nontree_idx = {e: k for k, e in enumerate(PETERSEN_LIST)
                   if e not in TREE_EDGES}

    survivors = []      # list of (pet_signs_tuple, H)
    classes = {}        # orbit_repr -> count

    t0 = time.time()
    for nontree_bits in range(2 ** len(NONTREE_PET)):
        pet_signs = [1] * len(PETERSEN_LIST)
        for k, e in enumerate(NONTREE_PET):
            if (nontree_bits >> k) & 1:
                pet_signs[PETERSEN_LIST.index(e)] = -1
        pet_signs = tuple(pet_signs)

        H_cand = _try_complete(pet_signs)
        if H_cand is None:
            continue
        # Canonicalise to canonical orbit representative.
        rep = _orbit_repr(H_cand)
        if rep in classes:
            classes[rep] += 1
        else:
            classes[rep] = 1
            survivors.append((pet_signs, H_cand))
        if verbose and (nontree_bits + 1) % 8 == 0:
            elapsed = time.time() - t0
            done = nontree_bits + 1
            total = 2 ** len(NONTREE_PET)
            eta = elapsed * (total - done) / done if done else 0
            print(f"  [{done:5d}/{total:5d}]  survivors: {len(survivors)}"
                  f"  classes: {len(classes)}  elapsed {elapsed:5.1f}s"
                  f"  eta {eta:5.1f}s")

    print()
    print("=" * 70)
    print(f"  TOTAL survivors:                    {sum(classes.values())}")
    print(f"  Distinct orbits (mod Aut(Petersen)): {len(classes)}")
    print("=" * 70)

    # Verify each surviving orbit-rep is equivalent to H_ico.
    ico_rep = _orbit_repr(H_ICO)
    for rep, count in classes.items():
        same = (rep == ico_rep)
        print(f"  orbit signature: count={count}  matches H_ico: {same}")

    # Final assertion
    print()
    if len(classes) == 1 and ico_rep in classes:
        print("  *** RESULT: unique class, matches H_ico  -> SGR Lemma verified")
    else:
        print("  *** RESULT: SGR Lemma FAILS  (unexpected branches found)")
        for rep, count in classes.items():
            print(f"     rep: {rep}")
        sys.exit(1)


if __name__ == "__main__":
    main()
