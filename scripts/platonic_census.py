"""
Platonic / finite-rotation-group census of two-distance Veronese lifts.

Goal:  given finite rotation subgroups G < SO(3), enumerate their canonical
orbits on the sphere (vertex, face-centre, edge-midpoint, antipodal-pair
representatives), and for each orbit compute:

  * orbit size  n  (= antipodal pairs if we project to RP^2)
  * the off-diagonal absolute Gram values  {|v_i . v_j| : i != j}
  * the number of distinct distance classes
  * if two-distance: the "large-distance" graph  (Petersen? Johnson? K_n?)
  * Veronese lift U_i = v_i v_i^T - I/3, its Gram, rank, spectrum
  * tight-frame check  (V^T V = lambda I_5  for the lifted vectors)
  * signed-Gram rigidity (raw count of rank-3 PSD signed lifts, # orbits
    modulo switching + automorphism)  --- skipped for n > 12 to keep cost
    bounded.

The output is a printed census table.

Conventions: all vectors are unit; we identify v_i with [v_i] in RP^2 by
collecting one representative per antipodal pair.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field
from typing import Callable

import numpy as np


PHI = (1.0 + math.sqrt(5.0)) / 2.0


# ----------------------------------------------------------------------
#   Finite rotation groups in SO(3)
# ----------------------------------------------------------------------
def _rotation_matrix(axis: np.ndarray, angle: float) -> np.ndarray:
    a = axis / np.linalg.norm(axis)
    c, s, C = math.cos(angle), math.sin(angle), 1 - math.cos(angle)
    x, y, z = a
    return np.array([
        [c + x * x * C,     x * y * C - z * s, x * z * C + y * s],
        [y * x * C + z * s, c + y * y * C,     y * z * C - x * s],
        [z * x * C - y * s, z * y * C + x * s, c + z * z * C],
    ])


def _close_group(generators: list[np.ndarray], tol: float = 1e-5,
                 max_order: int = 200) -> list[np.ndarray]:
    """BFS closure of <generators> in SO(3).  Uses a slightly loose tol to
    tolerate accumulated floating-point drift in long products.
    """
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
            # snap to nearest orthogonal to fight numerical drift
            U, _, Vt = np.linalg.svd(M)
            M = U @ Vt
            if np.linalg.det(M) < 0:
                M = -M
            if not already(M):
                group.append(M)
                frontier.append(M)
                if len(group) > max_order:
                    raise RuntimeError(
                        f"group blew up beyond {max_order} elements "
                        f"(tol = {tol})")
    return group


def group_T() -> list[np.ndarray]:
    """Tetrahedral rotation group A_4, order 12."""
    g1 = _rotation_matrix(np.array([1, 1, 1]), 2 * math.pi / 3)
    g2 = _rotation_matrix(np.array([1, 0, 0]), math.pi)
    return _close_group([g1, g2])


def group_O() -> list[np.ndarray]:
    """Octahedral rotation group S_4, order 24."""
    g1 = _rotation_matrix(np.array([0, 0, 1]), math.pi / 2)
    g2 = _rotation_matrix(np.array([1, 1, 1]), 2 * math.pi / 3)
    return _close_group([g1, g2])


def group_I() -> list[np.ndarray]:
    """Icosahedral rotation group A_5, order 60.

    Generators: a 5-fold rotation around a vertex axis and a 3-fold
    rotation around a face-centre axis (where the face must contain that
    vertex, i.e. the two axes are adjacent in the icosahedron).
    """
    ico_vert0 = np.array([0.0, 1.0, PHI])
    # A valid face containing ico_vert0: {(0,1,phi), (1,phi,0), (-1,phi,0)}
    face0 = (np.array([0.0, 1.0, PHI])
             + np.array([1.0, PHI, 0.0])
             + np.array([-1.0, PHI, 0.0]))
    g1 = _rotation_matrix(ico_vert0, 2 * math.pi / 5)
    g2 = _rotation_matrix(face0, 2 * math.pi / 3)
    return _close_group([g1, g2])


# ----------------------------------------------------------------------
#   Orbit + antipodal-pair representative extraction
# ----------------------------------------------------------------------
def _orbit(seed: np.ndarray, G: list[np.ndarray], tol: float = 1e-7) -> np.ndarray:
    """Orbit of seed under G; returns array of shape (orbit_size, 3)."""
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


def _antipodal_pair_reps(pts: np.ndarray, tol: float = 1e-7) -> np.ndarray:
    """Pick one representative per antipodal pair.  If a point has no
    antipodal partner in the orbit, include it as-is (single-pair orbit).
    """
    used = [False] * len(pts)
    reps = []
    for i in range(len(pts)):
        if used[i]:
            continue
        used[i] = True
        # find antipodal partner
        for j in range(i + 1, len(pts)):
            if not used[j] and np.linalg.norm(pts[i] + pts[j]) < tol:
                used[j] = True
                break
        reps.append(pts[i])
    return np.array(reps)


# ----------------------------------------------------------------------
#   Veronese lift + diagnostics
# ----------------------------------------------------------------------
def _veronese_lift_gram(V: np.ndarray) -> np.ndarray:
    """Frobenius-Gram of the lifted matrices U_i = v_i v_i^T - I/3."""
    n = len(V)
    G = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            G[i, j] = (V[i] @ V[j]) ** 2 - 1 / 3
    return G


def _signed_gram(V: np.ndarray) -> np.ndarray:
    return V @ V.T


def _classify_distances(H: np.ndarray, tol: float = 1e-5) -> dict[float, int]:
    """Return {rounded absolute off-diagonal value: count}."""
    n = H.shape[0]
    abs_vals = []
    for i in range(n):
        for j in range(i + 1, n):
            abs_vals.append(abs(H[i, j]))
    classes: dict[float, int] = {}
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


def _graph_signature(H: np.ndarray, large_val: float, tol: float = 1e-5) -> dict:
    """Edge-graph: i ~ j iff |H_ij| close to large_val. Report basic invariants."""
    n = H.shape[0]
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i + 1, n):
            if abs(abs(H[i, j]) - large_val) < tol:
                A[i, j] = A[j, i] = 1
    deg = A.sum(axis=1)
    # eigenvalues of A
    spec_A = sorted(np.round(np.linalg.eigvalsh(A), 4).tolist(), reverse=True)
    # triangle count
    tri = int(np.trace(A @ A @ A) // 6)
    return {
        "n": n,
        "edges": int(A.sum() // 2),
        "regular_degree": deg[0] if all(d == deg[0] for d in deg) else None,
        "spectrum": spec_A,
        "triangles": tri,
    }


def _identify_graph(sig: dict) -> str:
    """Try to name the graph from invariants."""
    n = sig["n"]
    e = sig["edges"]
    deg = sig["regular_degree"]
    spec = sig["spectrum"]
    tri = sig["triangles"]
    if e == 0:
        return "empty"
    if e == n * (n - 1) // 2:
        return f"K_{n}"
    # Petersen K(5,2): n=10, deg=3, spec=[3, 1^5, -2^4], triangle-free
    if (n == 10 and deg == 3 and tri == 0
            and abs(spec[0] - 3) < 1e-3 and abs(spec[-1] + 2) < 1e-3):
        return "Petersen K(5,2)"
    # K_{3,3}: n=6, deg=3, bipartite (spec [3, 0, 0, 0, 0, -3])
    if n == 6 and deg == 3 and tri == 0:
        return "K_{3,3}"
    # K_4: n=4, deg=3, triangle count = 4
    if n == 4 and e == 6:
        return "K_4"
    # Triangle K_3: n=3, deg=2
    if n == 3 and e == 3:
        return "K_3"
    # 3K_2 perfect matching: n=6, deg=1
    if deg == 1 and n % 2 == 0:
        return f"{n // 2}K_2 (perfect matching)"
    # K(4,2) triangular: n=6, deg=4 (cocktail party 3K_2 complement)
    if n == 6 and deg == 4:
        return "K_{2,2,2} (cocktail-party / octahedron-edge)"
    # Petersen complement: triangular T(5) = J(5,2): n=10, deg=6
    if n == 10 and deg == 6:
        return "T(5) = J(5,2) (triangular / Petersen complement)"
    # Kneser K(n,k) generic flag
    if n == 15 and deg == 6:
        return "L(K_6) = T(6) (triangular line graph)?"
    # Cube graph Q_3: n=8, deg=3, bipartite (no triangles), 4-regular faces
    if n == 8 and deg == 3 and tri == 0:
        return "Q_3 (3-cube)"
    # Bipartite K_{n/2, n/2}
    if deg == n // 2 and tri == 0 and n % 2 == 0:
        return f"K_{{{n//2},{n//2}}}"
    return f"unknown (n={n}, deg={deg}, spec={spec})"


@dataclass
class OrbitReport:
    group_name: str
    seed_name: str
    orbit_size: int        # total orbit on S^2
    pair_size: int         # antipodal-pair representatives in RP^2
    distance_classes: dict[float, int] = field(default_factory=dict)
    signed_gram_spec: list[float] = field(default_factory=list)
    veronese_gram_spec: list[float] = field(default_factory=list)
    tight_frame_3d: bool = False
    tight_frame_5d: bool = False
    graph_name: str = "-"


def _analyse_orbit(group_name: str, seed_name: str, seed: np.ndarray,
                   G: list[np.ndarray]) -> OrbitReport:
    pts = _orbit(seed, G)
    V = _antipodal_pair_reps(pts)
    n = len(V)
    report = OrbitReport(
        group_name=group_name,
        seed_name=seed_name,
        orbit_size=len(pts),
        pair_size=n,
    )
    if n < 2:
        return report
    H_signed = _signed_gram(V)
    report.distance_classes = _classify_distances(H_signed)
    spec = sorted(np.round(np.linalg.eigvalsh(H_signed), 4).tolist(), reverse=True)
    report.signed_gram_spec = spec

    # 3D tight-frame check
    Y = V                                    # n x 3
    YtY = Y.T @ Y
    if np.linalg.norm(YtY - (n / 3) * np.eye(3)) < 1e-4:
        report.tight_frame_3d = True

    # Veronese lift
    G5 = _veronese_lift_gram(V)
    spec5 = sorted(np.round(np.linalg.eigvalsh(G5), 4).tolist(), reverse=True)
    report.veronese_gram_spec = spec5

    # 5D tight-frame check via lift coordinates
    # Compute U_i in Sym^2_0 (5-dim coordinate via orthonormal basis):
    s2 = math.sqrt(2.0)
    s6 = math.sqrt(6.0)
    BASIS = [
        np.array([[1 / s2, 0, 0], [0, -1 / s2, 0], [0, 0, 0]]),
        np.array([[1 / s6, 0, 0], [0, 1 / s6, 0], [0, 0, -2 / s6]]),
        np.array([[0, 1 / s2, 0], [1 / s2, 0, 0], [0, 0, 0]]),
        np.array([[0, 0, 1 / s2], [0, 0, 0], [1 / s2, 0, 0]]),
        np.array([[0, 0, 0], [0, 0, 1 / s2], [0, 1 / s2, 0]]),
    ]

    def lift_coords(v):
        U = np.outer(v, v) - np.eye(3) / 3
        return np.array([np.sum(U * B) for B in BASIS])

    Y5 = np.array([lift_coords(v) for v in V])
    Y5tY5 = Y5.T @ Y5
    lam = Y5tY5[0, 0]
    if np.linalg.norm(Y5tY5 - lam * np.eye(5)) < 1e-4 and lam > 1e-6:
        report.tight_frame_5d = True

    # Identify graph (large-distance edges) if exactly two distance classes
    if len(report.distance_classes) == 2:
        vals = sorted(report.distance_classes.keys(), reverse=True)
        large = vals[0]
        sig = _graph_signature(H_signed, large)
        report.graph_name = _identify_graph(sig)
    elif len(report.distance_classes) == 1:
        v = list(report.distance_classes.keys())[0]
        if v < 1e-5:
            report.graph_name = "no edges (orthogonal frame)"
        else:
            report.graph_name = f"K_{n} (one distance {v:.4f})"

    return report


# ----------------------------------------------------------------------
#   Census driver
# ----------------------------------------------------------------------
def _polyhedron_seeds() -> list[tuple[str, str, np.ndarray, Callable[[], list]]]:
    """Return a list of (group_name, seed_name, seed_vector, group_factory)."""
    items = []

    # Tetrahedron (T = A_4, order 12)
    # vertices: (1,1,1), (1,-1,-1), (-1,1,-1), (-1,-1,1)
    items.append(("T (A_4)", "vertex", np.array([1, 1, 1]) / math.sqrt(3), group_T))
    # face-centre = opposite vertex direction (same orbit modulo antipodal here)
    items.append(("T (A_4)", "face", -np.array([1, 1, 1]) / math.sqrt(3), group_T))
    # edge midpoint: e.g. ((1,1,1) + (1,-1,-1))/2 = (1, 0, 0) -> length 1
    items.append(("T (A_4)", "edge", np.array([1.0, 0, 0]), group_T))

    # Octahedron / Cube (O = S_4, order 24)
    items.append(("O (S_4)", "octa-vertex", np.array([1.0, 0, 0]), group_O))
    items.append(("O (S_4)", "cube-vertex", np.array([1, 1, 1]) / math.sqrt(3), group_O))
    items.append(("O (S_4)", "edge-midpoint", np.array([1, 1, 0]) / math.sqrt(2), group_O))

    # Icosahedron / Dodecahedron (I = A_5, order 60)
    vert_seed = np.array([0.0, 1.0, PHI]) / math.sqrt(1 + PHI ** 2)
    items.append(("I (A_5)", "ico-vertex", vert_seed, group_I))
    # face centre: barycentre of a valid face
    f0 = (np.array([0.0, 1.0, PHI])
          + np.array([1.0, PHI, 0.0])
          + np.array([-1.0, PHI, 0.0]))
    items.append(("I (A_5)", "ico-face", f0 / np.linalg.norm(f0), group_I))
    # edge midpoint (a valid edge connecting two adjacent vertices)
    e0 = (np.array([0.0, 1.0, PHI]) + np.array([0.0, -1.0, PHI])) / 2
    items.append(("I (A_5)", "ico-edge", e0 / np.linalg.norm(e0), group_I))

    return items


def main():
    print("=" * 100)
    print(" Platonic / finite-rotation-group census of two-distance Veronese lifts")
    print("=" * 100)
    header = (f"{'group':<10} {'orbit':<13} {'orb':<4} {'pair':<5} "
              f"{'#cls':<5} {'distances':<32} {'5D tight':<9} {'graph (large)'}")
    print(header)
    print("-" * 100)

    reports = []
    for group_name, seed_name, seed, factory in _polyhedron_seeds():
        G = factory()
        rep = _analyse_orbit(group_name, seed_name, seed, G)
        reports.append(rep)
        dist_str = ", ".join(f"{v:.4f}({c})" for v, c in
                             sorted(rep.distance_classes.items(),
                                    key=lambda x: -x[0]))
        print(f"{rep.group_name:<10} {rep.seed_name:<13} {rep.orbit_size:<4} "
              f"{rep.pair_size:<5} {len(rep.distance_classes):<5} "
              f"{dist_str:<32} {'Y' if rep.tight_frame_5d else '-':<9} "
              f"{rep.graph_name}")

    print("-" * 100)
    print("\nDetailed dump for two-distance + 5D-tight-frame orbits:")
    print()
    for rep in reports:
        if rep.tight_frame_5d and len(rep.distance_classes) == 2:
            print(f"  {rep.group_name} / {rep.seed_name}:")
            print(f"     orbit_size = {rep.orbit_size}, pair_size = {rep.pair_size}")
            print(f"     distance classes: {rep.distance_classes}")
            print(f"     graph (large dist): {rep.graph_name}")
            print(f"     signed Gram spec: {rep.signed_gram_spec}")
            print(f"     Veronese Gram spec: {rep.veronese_gram_spec}")
            print(f"     3D tight: {rep.tight_frame_3d}, 5D tight: {rep.tight_frame_5d}")
            print()


if __name__ == "__main__":
    main()
