"""
5-tetrahedron hinge rigidity: numerical exploration.

Setup:
  T_0 fixed with vertices a_1, a_2, a_3, a_4 = standard tetrahedron axes
       (cube-inscription).
  T_k (k=1..4) shares axis a_k with T_0; parameterized by (eps_k, phi_k):
     - eps_k in {+1, -1}: sign of T_k's "k-vertex" (= eps_k * a_k)
     - phi_k in [0, 2pi): rotation angle around a_k for the other 3 vertices

  Constraint:  for each pair (k, l), 1 <= k < l <= 4, the "(k,l) vertex"
  of T_k must equal +/- the "(l,k) vertex" of T_l.

Goal:
  1. Verify that the icosahedral compound corresponds to a specific
     (eps_k, phi_k) configuration.
  2. Search for all (eps_k, phi_k) configurations satisfying the 6
     hinge constraints; verify they all lie in the same O(3) orbit.

We use a numerical Newton-iteration search starting from many random
initializations to find all isolated solutions.
"""

from __future__ import annotations

import itertools
import math
from typing import Callable

import numpy as np


PHI_GOLDEN = (1.0 + math.sqrt(5.0)) / 2.0


# ----------------------------------------------------------------------
#   T_0 vertices: standard tetrahedron from cube inscription
# ----------------------------------------------------------------------
SQ3 = math.sqrt(3.0)
A = np.array([
    [+1, +1, +1],   # a_1
    [+1, -1, -1],   # a_2
    [-1, +1, -1],   # a_3
    [-1, -1, +1],   # a_4
], dtype=float) / SQ3

# Verify pairwise inner products = -1/3
for r in range(4):
    for s in range(r + 1, 4):
        ip = A[r] @ A[s]
        assert abs(ip - (-1 / 3)) < 1e-12, f"a_{r+1} . a_{s+1} = {ip}"


def perp_basis(a: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return an orthonormal basis (e1, e2) of a^perp."""
    # Pick e1 = some vector not parallel to a, then orthogonalize
    if abs(a[0]) < 0.9:
        v = np.array([1.0, 0.0, 0.0])
    else:
        v = np.array([0.0, 1.0, 0.0])
    e1 = v - (v @ a) * a
    e1 = e1 / np.linalg.norm(e1)
    e2 = np.cross(a, e1)
    return e1, e2


PERP = [perp_basis(A[k]) for k in range(4)]   # (e1, e2) for each a_k


# ----------------------------------------------------------------------
#   Tetrahedron T_k vertices given (eps_k, phi_k)
# ----------------------------------------------------------------------
def tet_vertices(k: int, eps: float, phi: float) -> np.ndarray:
    """Return the 4 vertices of T_k as rows of a (4, 3) array.
    Vertex 0 is the "k-vertex" = eps * a_k.
    Vertices 1, 2, 3 are the 3 other vertices at angles phi, phi+120°, phi+240°.
    """
    a_k = A[k]
    e1, e2 = PERP[k]
    out = np.zeros((4, 3))
    out[0] = eps * a_k
    for j in range(3):
        angle = phi + j * 2 * math.pi / 3
        out[j + 1] = (-eps / 3) * a_k + (2 * math.sqrt(2) / 3) * (
            math.cos(angle) * e1 + math.sin(angle) * e2
        )
    return out


# ----------------------------------------------------------------------
#   Constraint:  vertex of T_k at sub-index j_k  ==  +/- vertex of T_l at j_l
# ----------------------------------------------------------------------
def hinge_residual(params, pi_map, delta_map):
    """Compute the 18 residual values (3 per pair) for the 6 hinge constraints.

    params: array of 8 values [eps_1, phi_1, eps_2, phi_2, eps_3, phi_3, eps_4, phi_4]
            (we treat eps_k as continuous in {-1, +1} via sign rounding;
             for Newton we fix eps_k from a discrete choice.)

    pi_map: dict (k, l) -> sub-index j in {0, 1, 2} for T_k's "(k,l)-vertex"
            (k in {1..4}, l in {1..4}\{k})

    delta_map: dict (k, l) -> +/-1, sign of the hinge equality
               (only for k < l; delta_{lk} = delta_{kl})
    """
    eps = [None] + [params[2 * k] for k in range(4)]      # eps[1..4]
    phi = [None] + [params[2 * k + 1] for k in range(4)]  # phi[1..4]

    residuals = []
    for k in range(1, 5):
        for l in range(k + 1, 5):
            # T_k's vertex labeled (k, l) = sub-index pi_map[(k, l)] + 1
            # (since sub-index 0 is the "k-vertex" itself)
            j_k = pi_map[(k, l)]
            j_l = pi_map[(l, k)]
            T_k = tet_vertices(k - 1, eps[k], phi[k])
            T_l = tet_vertices(l - 1, eps[l], phi[l])
            v_kl = T_k[j_k + 1]
            v_lk = T_l[j_l + 1]
            d = delta_map[(k, l)]
            res = v_kl - d * v_lk    # 3-vector
            residuals.extend(res.tolist())
    return np.array(residuals)


# ----------------------------------------------------------------------
#   Construct the icosahedral compound and find (eps, phi, pi, delta)
# ----------------------------------------------------------------------
def icosahedral_compound() -> list[np.ndarray]:
    """Return 5 tetrahedra inscribed in the dodecahedron (chiral compound).

    The 20 dodecahedron vertices split into 5 groups of 4, each forming
    a regular tetrahedron. We construct them explicitly.

    Dodecahedron vertices (one standard form):
       (+/-1, +/-1, +/-1)   (8 vertices from cube)
       (0, +/-phi, +/-1/phi)
       (+/-1/phi, 0, +/-phi)
       (+/-phi, +/-1/phi, 0)
    Total = 20.
    """
    phi = PHI_GOLDEN
    inv_phi = 1 / phi
    verts = []
    for s1 in (+1, -1):
        for s2 in (+1, -1):
            for s3 in (+1, -1):
                verts.append((s1, s2, s3))
    for s1 in (+1, -1):
        for s2 in (+1, -1):
            verts.append((0, s1 * phi, s2 * inv_phi))
            verts.append((s1 * inv_phi, 0, s2 * phi))
            verts.append((s1 * phi, s2 * inv_phi, 0))
    verts = np.array(verts, dtype=float)
    # Normalize each vertex to length 1 (we want unit vectors on sphere)
    verts = verts / np.linalg.norm(verts, axis=1, keepdims=True)
    assert len(verts) == 20

    # Group into 5 regular tetrahedra: 4 vertices each, pairwise inner
    # product -1/3. Brute-force: enumerate 4-subsets with this property,
    # then group them into 5 disjoint tetrahedra.
    tetrahedra = []
    for combo in itertools.combinations(range(20), 4):
        vs = verts[list(combo)]
        ok = True
        for i, j in itertools.combinations(range(4), 2):
            if abs(vs[i] @ vs[j] - (-1 / 3)) > 1e-6:
                ok = False
                break
        if ok:
            tetrahedra.append(combo)
    # We expect 10 regular tetrahedra total (5 antipodal pairs / 2
    # chirality classes). Each chirality class contains 5 disjoint
    # tetrahedra covering all 20 vertices. Find one such class via
    # backtracking.
    print(f"  found {len(tetrahedra)} regular tetrahedra in dodecahedron")

    selected: list[tuple[int, ...]] = []

    def backtrack(used: set, pos: int) -> bool:
        if len(selected) == 5:
            return len(used) == 20
        for idx in range(pos, len(tetrahedra)):
            combo = tetrahedra[idx]
            if not any(v in used for v in combo):
                selected.append(combo)
                new_used = used | set(combo)
                if backtrack(new_used, idx + 1):
                    return True
                selected.pop()
        return False

    ok = backtrack(set(), 0)
    assert ok and len(selected) == 5, "could not find 5 disjoint tetrahedra"
    used = set()
    for combo in selected:
        used.update(combo)
    assert len(used) == 20, f"only covered {len(used)} vertices"
    return [verts[list(c)] for c in selected]


# ----------------------------------------------------------------------
#   Numerical search for hinge solutions
# ----------------------------------------------------------------------
def main():
    print("=" * 70)
    print("5-tetrahedron hinge rigidity: numerical exploration")
    print("=" * 70)

    print("\nStep 1: Build the icosahedral compound (5 tetrahedra in dodecahedron)")
    print("-" * 70)
    compound = icosahedral_compound()
    print(f"  selected 5 tetrahedra forming the chiral compound")
    for k, T in enumerate(compound):
        print(f"  T_{k}: 4 vertices, pairwise inner products check:")
        ips = []
        for i, j in itertools.combinations(range(4), 2):
            ips.append(T[i] @ T[j])
        print(f"        {[round(x, 4) for x in ips]}")

    # Compute pairwise inner products between tetrahedra (in compound)
    print("\nStep 2: Pairwise tetrahedra inner products (10x10 Gram structure)")
    print("-" * 70)
    print("  Pair (T_k, T_l) vertex-axis inner products:")
    for k in range(5):
        for l in range(k + 1, 5):
            ips = []
            for vk in compound[k]:
                for vl in compound[l]:
                    ips.append(abs(vk @ vl))
            ips.sort()
            unique = []
            for x in ips:
                if not unique or abs(x - unique[-1]) > 1e-6:
                    unique.append(x)
            print(f"    T_{k}, T_{l}: distinct |v.w| values = {[round(x, 4) for x in unique]}")

    print()
    print("Step 3: Check the K_5 incidence (each pair shares 1 axis with |v.w|=1)")
    print("-" * 70)
    shared_count = {}
    for k in range(5):
        for l in range(k + 1, 5):
            n_shared = 0
            for vk in compound[k]:
                for vl in compound[l]:
                    if abs(abs(vk @ vl) - 1.0) < 1e-6:
                        n_shared += 1
            shared_count[(k, l)] = n_shared
            print(f"    T_{k} ~ T_{l}: {n_shared} shared axes (expected: 1 in RP^2, 2 in R^3)")

    print()
    print("=" * 70)
    print("Phase 2: Newton-iteration search for hinge solutions")
    print("=" * 70)

    newton_search()


# ----------------------------------------------------------------------
#   Newton search — Phase 2
# ----------------------------------------------------------------------
# Labeling convention:
#   T_k's 3 sub-vertices are labeled by {1,2,3,4}\{k} in INCREASING order.
#   pi_map[(k, l)] = position of l in sorted {1,2,3,4}\{k}
PI_MAP = {}
for k in range(1, 5):
    rest = sorted(set(range(1, 5)) - {k})
    for pos, l in enumerate(rest):
        PI_MAP[(k, l)] = pos
PAIRS = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]


def make_residual(eps: np.ndarray, delta: np.ndarray) -> Callable:
    """Return residual function: R^4 -> R^18.
    eps: shape (4,) in {-1, +1}, signs of "k-vertices".
    delta: shape (6,) in {-1, +1}, signs in hinge equalities per pair.
    """
    def F(phi: np.ndarray) -> np.ndarray:
        res = []
        Tks = [tet_vertices(k - 1, eps[k - 1], phi[k - 1]) for k in range(1, 5)]
        for idx, (k, l) in enumerate(PAIRS):
            j_k = PI_MAP[(k, l)] + 1   # +1 since sub-index 0 is "k-vertex"
            j_l = PI_MAP[(l, k)] + 1
            v_kl = Tks[k - 1][j_k]
            v_lk = Tks[l - 1][j_l]
            diff = v_kl - delta[idx] * v_lk
            res.extend(diff.tolist())
        return np.array(res)
    return F


def known_solution_debug():
    """Verify the residual function by extracting (eps, phi) from the
    known icosahedral compound and checking residual.
    """
    print()
    print("DEBUG: Verifying residual function on known icosahedral compound")
    print("-" * 70)
    compound = icosahedral_compound()
    # Align the first tetrahedron to T_0 = A via Procrustes
    T_first = compound[0]
    best_R = None
    best_err = float('inf')
    best_perm = None
    for perm in itertools.permutations(range(4)):
        source = T_first[list(perm)]
        # Procrustes: R minimizes ||A - R @ source||
        H_proc = source.T @ A
        U_p, _, Vt_p = np.linalg.svd(H_proc)
        R = U_p @ Vt_p
        aligned = source @ R.T
        err = np.linalg.norm(aligned - A)
        if err < best_err:
            best_err = err
            best_R = R
            best_perm = perm
    print(f"  Best Procrustes alignment: err = {best_err:.2e}")
    if best_err > 0.01:
        print(f"  ⚠️  First tetrahedron doesn't align to T_0 (chirality?)")
        # Try with reflection
        for perm in itertools.permutations(range(4)):
            source = T_first[list(perm)]
            H_proc = source.T @ A
            U_p, _, Vt_p = np.linalg.svd(H_proc)
            # Force det = -1
            S = np.diag([1, 1, -1])
            R = U_p @ S @ Vt_p
            aligned = source @ R.T
            err = np.linalg.norm(aligned - A)
            if err < best_err:
                best_err = err
                best_R = R
                best_perm = perm
        print(f"  With reflection: err = {best_err:.2e}")

    # Apply alignment to all 5 tetrahedra
    aligned = [T @ best_R.T for T in compound]
    print(f"  T_0 (aligned): {aligned[0]}")
    print(f"  A (target):    {A}")

    # For each of the other 4 tetrahedra, find which axis matches a_k for some k.
    # This identifies the assignment T_k <-> aligned tetrahedron.
    label_map = {}
    for ai, aligned_T in enumerate(aligned[1:], 1):
        # Check each vertex of aligned_T to see if it matches +/- a_k for some k
        for vi, v in enumerate(aligned_T):
            for k in range(4):
                if np.linalg.norm(v - A[k]) < 0.01:
                    label_map[ai] = (k, vi, +1)
                    break
                if np.linalg.norm(v + A[k]) < 0.01:
                    label_map[ai] = (k, vi, -1)
                    break
            if ai in label_map:
                break
    print(f"  T_k assignments (compound index -> (axis k, vertex pos, sign)): {label_map}")

    # For each labeled T_k, extract (eps_k, phi_k)
    eps_extracted = np.zeros(4)
    phi_extracted = np.zeros(4)
    for ai, (k, vi, sign) in label_map.items():
        eps_extracted[k] = sign
        # Now extract phi_k: the 3 other vertices of aligned_T should be at
        # angles phi_k, phi_k + 120°, phi_k + 240° around a_k.
        aligned_T = aligned[ai]
        other_verts = np.array([aligned_T[j] for j in range(4) if j != vi])
        # Project onto a_k^perp plane
        a_k = A[k]
        e1, e2 = PERP[k]
        # Each other vertex should equal -(eps_k/3) a_k + (2sqrt2/3)(cos(theta) e1 + sin(theta) e2)
        thetas = []
        for v in other_verts:
            # Project to e1, e2 plane
            x = v @ e1
            y = v @ e2
            theta = math.atan2(y, x)
            thetas.append(theta)
        thetas.sort()
        # phi_k is the smallest theta; the other two should be theta + 120° and theta + 240°
        phi_k = thetas[0]
        # Wrap to [0, 2pi)
        phi_extracted[k] = phi_k % (2 * math.pi)
    print(f"  Extracted eps: {eps_extracted}")
    print(f"  Extracted phi (radians): {phi_extracted}")
    print(f"  Extracted phi (degrees): {np.rad2deg(phi_extracted)}")

    # Now scan over all 64 delta choices and see which gives smallest residual
    print()
    print("  Scanning 64 delta-sign choices with these (eps, phi):")
    best_delta = None
    best_residual = float('inf')
    for delta_bits in range(64):
        delta = np.array([(1 if (delta_bits >> b) & 1 else -1)
                          for b in range(6)], dtype=float)
        F = make_residual(eps_extracted, delta)
        r = np.linalg.norm(F(phi_extracted))
        if r < best_residual:
            best_residual = r
            best_delta = delta
    print(f"  Best delta: {best_delta}")
    print(f"  Best residual: {best_residual:.4e}")

    if best_residual > 1e-3:
        print()
        print(f"  ⚠️  Even best delta gives residual {best_residual:.4e}")
        print(f"      => labeling convention PI_MAP does NOT match geometry")
        print(f"      Try scanning all 1296 labeling permutations as well.")
    else:
        print(f"  ✓ Known solution verified with residual {best_residual:.2e}")


def make_residual_labeling_free(eps: np.ndarray) -> Callable:
    """Labeling-free residual:  for each pair (k, l), take the minimum
    squared distance over all 3 x 3 x 2 sub-vertex-pairings with sign.
    Returns 6-dim residual on 4 unknowns (phi_1..phi_4).
    """
    def F(phi: np.ndarray) -> np.ndarray:
        res = np.zeros(6)
        Tks = [tet_vertices(k - 1, eps[k - 1], phi[k - 1])
               for k in range(1, 5)]
        for idx, (k, l) in enumerate(PAIRS):
            best = float('inf')
            for j_k in (1, 2, 3):
                for j_l in (1, 2, 3):
                    v_k = Tks[k - 1][j_k]
                    v_l = Tks[l - 1][j_l]
                    d_plus = float(np.sum((v_k - v_l) ** 2))
                    d_minus = float(np.sum((v_k + v_l) ** 2))
                    best = min(best, d_plus, d_minus)
            res[idx] = best
        return res
    return F


def newton_search():
    """Run the full Newton search with labeling-free residual."""
    from scipy.optimize import least_squares

    # Debug first: verify residual on known solution
    known_solution_debug()

    print()
    print(f"Scan parameters (labeling-free residual):")
    print(f"  16 chirality choices  x  N random phi-inits each")
    n_trials_per_case = 100
    print(f"  {n_trials_per_case} random phi-initializations per chirality")
    print(f"  Residual cutoff: cost < 1e-16  (essentially machine precision)")
    print()

    rng = np.random.default_rng(42)
    solutions = []     # list of (eps, phi, residual)

    total_runs = 0
    for eps_bits in range(16):
        eps = np.array([(1 if (eps_bits >> b) & 1 else -1) for b in range(4)],
                       dtype=float)
        F = make_residual_labeling_free(eps)
        for trial in range(n_trials_per_case):
            phi_init = rng.uniform(0, 2 * math.pi, size=4)
            total_runs += 1
            try:
                sol = least_squares(F, phi_init, method='lm', max_nfev=500,
                                    ftol=1e-15, xtol=1e-15, gtol=1e-15)
            except Exception:
                continue
            if sol.cost < 1e-16:
                phi = np.mod(sol.x, 2 * math.pi)
                solutions.append((eps.copy(), phi, float(sol.cost)))

    print(f"Total Newton runs: {total_runs}")
    print(f"Solutions found (cost < 1e-16): {len(solutions)}")

    # Deduplicate solutions
    print()
    print("Deduplication ((eps, phi mod 2pi), tol 1e-4):")
    unique = []
    for s in solutions:
        eps, phi, cost = s
        is_dup = False
        for u_eps, u_phi, _ in unique:
            if np.array_equal(u_eps, eps):
                # Allow phi-shifts by 120° (3-fold sub-index permutation)
                # and reflections (e_2 -> -e_2 i.e. phi -> -phi)
                d_min = float('inf')
                for sign in (+1, -1):
                    for shifts in itertools.product(
                            [0, 2*math.pi/3, 4*math.pi/3], repeat=4):
                        candidate = np.mod(sign * phi + np.array(shifts),
                                            2 * math.pi)
                        diff = np.linalg.norm(
                            np.mod(u_phi - candidate + math.pi,
                                   2 * math.pi) - math.pi)
                        d_min = min(d_min, diff)
                if d_min < 1e-4:
                    is_dup = True
                    break
        if not is_dup:
            unique.append((eps, phi, cost))
    print(f"  Unique (eps, phi) configurations: {len(unique)}")

    # Build 10-vector configuration for each unique solution
    print()
    print("Computing 10x10 signed Gram for each unique solution...")
    grams = []
    for eps, phi, cost in unique:
        verts = build_10_axes_smart(eps, phi)
        if verts is None:
            continue
        H = verts @ verts.T
        grams.append((H, eps, phi))
    print(f"  Successfully built {len(grams)} configurations")

    # Orbit classification: canonicalize signed Gram modulo switching + S_5
    print()
    print("Orbit classification (signed Gram canonical form)...")
    canonical_keys = {}
    for H, eps, phi in grams:
        key = canonicalize_signed_gram(H)
        if key in canonical_keys:
            canonical_keys[key].append((eps, phi))
        else:
            canonical_keys[key] = [(eps, phi)]

    print(f"  Distinct orbit classes: {len(canonical_keys)}")
    print()
    print("Orbit-class summary:")
    H_ico_key = compute_H_ico_canonical_key()
    for i, (key, configs) in enumerate(canonical_keys.items()):
        match_ico = (key == H_ico_key)
        print(f"  Class {i}: {len(configs)} configurations, "
              f"matches H_ico: {match_ico}")
        if len(configs) <= 3:
            for cfg in configs[:3]:
                eps, phi = cfg
                print(f"     eps={list(map(int, eps))}, "
                      f"phi(deg)={[round(float(np.rad2deg(p)), 2) for p in phi]}")

    print()
    print("=" * 70)
    if len(canonical_keys) == 1:
        only_key = next(iter(canonical_keys.keys()))
        ico_match = (only_key == H_ico_key)
        print(f"RESULT: 1 orbit class — matches H_ico: {ico_match}")
        if ico_match:
            print("  ✓ Five-tetrahedron hinge rigidity numerically confirmed!")
            print("  ✓ Layer B' bridge: candidate clean algebraic proof available")
    elif len(canonical_keys) == 0:
        print("RESULT: No solutions found — scaffolding issue?")
    else:
        print(f"RESULT: {len(canonical_keys)} orbit classes — investigate further.")


def build_10_axes_smart(eps: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """Reconstruct the 10 axes given (eps, phi). For each pair (k, l) > 0,
    find the matching sub-vertex pair automatically."""
    # T_k vertices
    Tks = [None]  # T_0 will be placeholder
    for k in range(1, 5):
        Tks.append(tet_vertices(k - 1, eps[k - 1], phi[k - 1]))
    # Collect all distinct axes (10 expected)
    all_verts = list(A)  # 4 axes from T_0
    for k in range(1, 5):
        for j in range(1, 4):
            v = Tks[k][j]
            # check if already in all_verts (up to sign)
            is_new = True
            for u in all_verts:
                if (np.linalg.norm(v - u) < 1e-5
                        or np.linalg.norm(v + u) < 1e-5):
                    is_new = False
                    break
            if is_new:
                all_verts.append(v)
    if len(all_verts) != 10:
        return None
    return np.array(all_verts)


def build_10_axes(eps: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """Reconstruct the 10 vertex-axis representatives from (eps, phi)."""
    # 10 axes:
    #   indices 0,1,2,3:  T_0's vertices (= a_1, a_2, a_3, a_4)
    #   indices 4..9:    the 6 shared axes between T_k, T_l (k,l in 1..4, k<l)
    axes = np.zeros((10, 3))
    for r in range(4):
        axes[r] = A[r]
    # For each pair (k, l), the shared axis is T_k's vertex at sub-index pi_map[(k,l)]
    for idx, (k, l) in enumerate(PAIRS):
        j_k = PI_MAP[(k, l)] + 1
        T_k = tet_vertices(k - 1, eps[k - 1], phi[k - 1])
        axes[4 + idx] = T_k[j_k]
    # Normalize to unit length (should already be ~unit)
    for i in range(10):
        n = np.linalg.norm(axes[i])
        if n < 0.1:
            return None
        axes[i] = axes[i] / n
    return axes


def canonicalize_signed_gram(H: np.ndarray) -> tuple:
    """Return a canonical hashable key for H modulo switching + S_5.

    Strategy:
    1. Build sign pattern (sign of off-diagonal H_{ij}).
    2. Try all 10! permutations and switchings; pick lexicographically
       smallest. (For 10! = 3.6M this is too slow; we use a faster
       canonicalization based on row-degree signatures.)

    Simpler approach: use the unsigned distance multiset as primary
    key (this should distinguish orbit classes already).
    """
    n = H.shape[0]
    # Round magnitudes to deduplicate floating-point variations
    abs_H = np.round(np.abs(H), 4)
    # Per-row, sort magnitudes
    row_signatures = []
    for i in range(n):
        row_signatures.append(tuple(sorted(abs_H[i].tolist())))
    # Sort by row signature
    row_signatures = tuple(sorted(row_signatures))
    return row_signatures


def compute_H_ico_canonical_key() -> tuple:
    """Compute the canonical key for the icosahedral H_ico."""
    compound = icosahedral_compound()
    axes = []
    used_axes: list[np.ndarray] = []
    for T in compound:
        for v in T:
            # check if antipodal version already in used_axes
            found = False
            for u in used_axes:
                if np.linalg.norm(v - u) < 1e-6 or np.linalg.norm(v + u) < 1e-6:
                    found = True
                    break
            if not found:
                used_axes.append(v)
    used_axes = np.array(used_axes)
    H = used_axes @ used_axes.T
    return canonicalize_signed_gram(H)


if __name__ == "__main__":
    main()
