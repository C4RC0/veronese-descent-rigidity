"""
Detailed analysis of the 15-edge-pair orbit of A_5 on S^2.

The icosahedron has 30 edges, giving 15 antipodal pairs after identification.
The census found this orbit has FOUR distance classes (not two), so it is
not a simple two-distance / strongly-regular set.  This script:

  1.  Constructs the 15 unit vectors and their Gram matrix.
  2.  Classifies the four distance classes by size and value.
  3.  For each non-zero distance class, builds the graph on 15 vertices
      and reports its spectrum, regularity, triangle count, and identifies
      it against standard graphs (Johnson J(6,2), triangular T(6),
      Kneser, line graphs, etc.).
  4.  Checks whether the four classes together form an association scheme
      (i.e. whether they form a coherent algebra via intersection numbers).
  5.  Reports the Veronese-lift spectrum and confirms 5D tight frame.
  6.  Identifies which 4-class scheme this is in standard tables
      (suspect: the 4-class scheme of edges of the icosahedron, or the
      bipartite double of T(6), or related to the half-cube graph).
"""

from __future__ import annotations

import itertools
import math

import numpy as np


PHI = (1.0 + math.sqrt(5.0)) / 2.0


def _rotation_matrix(axis, angle):
    a = axis / np.linalg.norm(axis)
    c, s, C = math.cos(angle), math.sin(angle), 1 - math.cos(angle)
    x, y, z = a
    return np.array([
        [c + x * x * C,     x * y * C - z * s, x * z * C + y * s],
        [y * x * C + z * s, c + y * y * C,     y * z * C - x * s],
        [z * x * C - y * s, z * y * C + x * s, c + z * z * C],
    ])


def _close_group(generators, tol=1e-5, max_order=100):
    group = [np.eye(3)]

    def already(M):
        for X in group:
            if np.linalg.norm(M - X) < tol:
                return True
        return False

    frontier = []
    for g in generators:
        if not already(g):
            group.append(g)
            frontier.append(g)
    while frontier:
        h = frontier.pop()
        for g in generators:
            M = g @ h
            U, _, Vt = np.linalg.svd(M)
            M = U @ Vt
            if np.linalg.det(M) < 0:
                M = -M
            if not already(M):
                group.append(M)
                frontier.append(M)
                if len(group) > max_order:
                    raise RuntimeError("group blew up")
    return group


def _A5() -> list[np.ndarray]:
    v0 = np.array([0.0, 1.0, PHI])
    f0 = (np.array([0.0, 1.0, PHI])
          + np.array([1.0, PHI, 0.0])
          + np.array([-1.0, PHI, 0.0]))
    g1 = _rotation_matrix(v0, 2 * math.pi / 5)
    g2 = _rotation_matrix(f0, 2 * math.pi / 3)
    return _close_group([g1, g2])


def _orbit(seed, G, tol=1e-6):
    pts = []
    for g in G:
        v = g @ seed
        v = v / np.linalg.norm(v)
        is_new = True
        for w in pts:
            if np.linalg.norm(v - w) < tol:
                is_new = False
                break
        if is_new:
            pts.append(v)
    return np.array(pts)


def _antipodal_pair_reps(pts, tol=1e-6):
    used = [False] * len(pts)
    reps = []
    for i in range(len(pts)):
        if used[i]:
            continue
        used[i] = True
        for j in range(i + 1, len(pts)):
            if not used[j] and np.linalg.norm(pts[i] + pts[j]) < tol:
                used[j] = True
                break
        reps.append(pts[i])
    return np.array(reps)


def _classify_distances(H, tol=1e-5):
    n = H.shape[0]
    abs_vals = []
    for i in range(n):
        for j in range(i + 1, n):
            abs_vals.append(abs(H[i, j]))
    classes = {}
    abs_vals.sort()
    if not abs_vals:
        return classes
    current = [abs_vals[0]]
    for v in abs_vals[1:]:
        if abs(v - current[-1]) < tol:
            current.append(v)
        else:
            classes[round(sum(current) / len(current), 8)] = len(current)
            current = [v]
    classes[round(sum(current) / len(current), 8)] = len(current)
    return classes


def _build_adjacency(H, target_val, tol=1e-5):
    n = H.shape[0]
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i + 1, n):
            if abs(abs(H[i, j]) - target_val) < tol:
                A[i, j] = A[j, i] = 1
    return A


def _graph_diagnostics(A: np.ndarray, label: str):
    n = A.shape[0]
    deg = A.sum(axis=1)
    regular = all(d == deg[0] for d in deg)
    spec = sorted(np.round(np.linalg.eigvalsh(A.astype(float)), 6).tolist(),
                  reverse=True)
    tri = int(np.trace(A @ A @ A) // 6)
    return {
        "label": label,
        "n": n,
        "edges": int(A.sum() // 2),
        "regular_degree": deg[0] if regular else f"non-reg deg={deg.tolist()}",
        "spectrum": spec,
        "triangles": tri,
    }


def _check_association_scheme(adjs: list[np.ndarray]) -> dict:
    """Check whether {I, A_1, A_2, ..., A_d} forms a (commutative) association
    scheme:  each A_k A_l should be a non-negative integer linear combination
    of {I, A_1, ..., A_d}.
    """
    n = adjs[0].shape[0]
    I = np.eye(n, dtype=int)
    basis = [I] + adjs
    d = len(basis)
    # Vectorise each basis matrix as length-n^2 vector.
    B = np.stack([M.flatten() for M in basis])
    Btarget = B.astype(float)

    coeffs = np.zeros((d, d, d))
    is_scheme = True
    for k in range(1, d):
        for l in range(k, d):
            prod = basis[k] @ basis[l]
            prod_vec = prod.flatten().astype(float)
            # Solve prod_vec = sum_m c_m * basis_m_vec
            # Using least squares; check residual.
            c, res, _, _ = np.linalg.lstsq(Btarget.T, prod_vec, rcond=None)
            recon = (Btarget.T @ c)
            err = np.linalg.norm(recon - prod_vec)
            if err > 1e-6:
                is_scheme = False
            # Check non-negative integers
            for m, cm in enumerate(c):
                if abs(cm - round(cm)) > 1e-6 or round(cm) < 0:
                    if abs(cm) > 1e-6:
                        is_scheme = False
                coeffs[k, l, m] = round(cm)
                coeffs[l, k, m] = round(cm)
    return {
        "is_scheme": is_scheme,
        "intersection_numbers": coeffs,
    }


def main():
    G = _A5()
    assert len(G) == 60, f"|A_5| should be 60, got {len(G)}"
    # Edge midpoint seed: midpoint of an edge between two adjacent vertices
    e0 = (np.array([0.0, 1.0, PHI]) + np.array([0.0, -1.0, PHI])) / 2
    e0 = e0 / np.linalg.norm(e0)
    pts = _orbit(e0, G)
    V = _antipodal_pair_reps(pts)
    print(f"orbit size: {len(pts)}  ->  antipodal pair reps: {len(V)}")

    H = V @ V.T
    classes = _classify_distances(H)
    print()
    print("Distance classes (|H_ij| value -> count):")
    for v, c in sorted(classes.items(), key=lambda x: -x[0]):
        print(f"  |H_ij| = {v:.6f}   count = {c}")

    nonzero_vals = [v for v in classes if v > 1e-6]
    print()
    print(f"Non-zero distance classes: {len(nonzero_vals)}")
    adjs = []
    for v in sorted(nonzero_vals, reverse=True):
        A = _build_adjacency(H, v)
        diag = _graph_diagnostics(A, f"|H|={v:.4f}")
        adjs.append(A)
        print(f"  graph (|H_ij| = {v:.6f}):  "
              f"n={diag['n']}, edges={diag['edges']}, "
              f"reg_deg={diag['regular_degree']}, "
              f"triangles={diag['triangles']}")
        print(f"    spectrum: {diag['spectrum']}")

    # Add an "orthogonal" (H_ij == 0) graph as a separate scheme class
    A_ortho = _build_adjacency(H, 0.0)
    diag_ortho = _graph_diagnostics(A_ortho, "|H|=0 (orthogonal)")
    adjs.append(A_ortho)
    print(f"  graph (|H_ij| = 0, orthogonal):  "
          f"n={diag_ortho['n']}, edges={diag_ortho['edges']}, "
          f"reg_deg={diag_ortho['regular_degree']}, "
          f"triangles={diag_ortho['triangles']}")
    print(f"    spectrum: {diag_ortho['spectrum']}")

    # Verify A_1 + A_2 + A_3 + A_4 + I = J  (partition of complete graph)
    total = sum(adjs) + np.eye(len(V), dtype=int)
    assert np.all(total == 1), f"classes do not partition K_15: max = {total.max()}"
    print()
    print("All 4 class-graphs + I partition K_15: OK")

    print()
    print("Association-scheme check (intersection numbers):")
    res = _check_association_scheme(adjs)
    print(f"  forms a commutative association scheme: {res['is_scheme']}")
    if res["is_scheme"]:
        print(f"  intersection-number tensor shape: "
              f"{res['intersection_numbers'].shape}")
        # Print intersection numbers as tables
        d = len(adjs) + 1
        names = ["I", "A1", "A2", "A3", "A4"]
        print()
        print(f"  Intersection numbers p_{{kl}}^m  (rows k, cols l within each m-table):")
        c = res["intersection_numbers"]
        for m in range(d):
            print(f"\n   m = {m} ({names[m]}):  p_{{kl}}^{m}")
            print("        " + "  ".join(f"{names[l]:>4}" for l in range(d)))
            for k in range(d):
                row = "  ".join(f"{int(c[k][l][m]):>4}" for l in range(d))
                print(f"   {names[k]}: {row}")

    # Compare with Johnson J(6, 2) = T(6) intersection numbers:
    # J(6, 2): vertices = 2-subsets of [6], 3 nonzero distance classes (0, 1, 2)
    # distance = |A symm-diff B|/2.
    # That's only 3 classes (i.e. 2 nondiagonal).  Our orbit has 4 nondiagonal.
    # So this is a FINER scheme than J(6, 2).
    #
    # Candidate: the 4-class fission scheme on 15 = C(6, 2) refining J(6, 2).
    # Or: edge graph of icosahedron as a graph of its own.

    print()
    print("Veronese lift Gram (in Sym^2_0):")
    G5 = np.zeros((len(V), len(V)))
    for i in range(len(V)):
        for j in range(len(V)):
            G5[i, j] = (V[i] @ V[j]) ** 2 - 1 / 3
    spec5 = sorted(np.round(np.linalg.eigvalsh(G5), 6).tolist(), reverse=True)
    print(f"  spectrum: {spec5}")
    rank5 = np.linalg.matrix_rank(G5, tol=1e-8)
    print(f"  rank: {rank5}")

    # Tight-frame check
    s2 = math.sqrt(2.0); s6 = math.sqrt(6.0)
    BASIS = [
        np.array([[1 / s2, 0, 0], [0, -1 / s2, 0], [0, 0, 0]]),
        np.array([[1 / s6, 0, 0], [0, 1 / s6, 0], [0, 0, -2 / s6]]),
        np.array([[0, 1 / s2, 0], [1 / s2, 0, 0], [0, 0, 0]]),
        np.array([[0, 0, 1 / s2], [0, 0, 0], [1 / s2, 0, 0]]),
        np.array([[0, 0, 0], [0, 0, 1 / s2], [0, 1 / s2, 0]]),
    ]

    def lift(v):
        U = np.outer(v, v) - np.eye(3) / 3
        return np.array([np.sum(U * B) for B in BASIS])

    Y5 = np.array([lift(v) for v in V])
    Y5tY5 = Y5.T @ Y5
    print(f"  Y^T Y = ")
    print(np.round(Y5tY5, 4))
    lam = Y5tY5[0, 0]
    is_tight = np.linalg.norm(Y5tY5 - lam * np.eye(5)) < 1e-4
    print(f"  5D tight frame (Y^T Y = lambda * I_5)?  {is_tight}, "
          f"lambda = {lam:.4f}")


if __name__ == "__main__":
    main()
