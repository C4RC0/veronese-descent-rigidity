"""
Exact orbit verification for the five-tetrahedron hinge solutions.

For each of the two real roots
    s_4 = (sqrt(3) - sqrt(15)) / 8   (the "-"-root, classical icosahedral)
    s_4 = (sqrt(3) + sqrt(15)) / 8   (the "+"-root)
of the canonical quadratic of Theorem 5.4, this script:

  1. Reconstructs the 10 axes of the five-tetrahedron compound exactly,
     using the rank-7 solution family
        c_1 = c_2 = -1/4,        c_3 = c_4 = sqrt(3)*s_4 - 1/2
        s_1 = s_2 = 2*s_4 - sqrt(3)/4,   s_3 = s_4,
     with chirality eps = (-1, -1, -1, -1).

  2. Builds the 10x10 signed Gram matrix G_{ij} = <v_i, v_j> in exact
     SymPy arithmetic, where the 10 axes v_i are indexed by the
     canonical Kneser labels (i, j) in C(5,2).

  3. Sanity checks:
       - Petersen magnitude pattern: |G_{ij}|^2 in {1, 5/9, 1/9}.
       - Bose-Mesner identity: Q = G entrywise squared equals
            8/9 I + 4/9 A + 1/9 J     (Proposition 2.1).
       - rank(G) = 3 (computed via the 3x10 matrix V of axes).

  4. Orbit-belonging certificate: takes the "-"-root signed Gram as the
     reference icosahedral signed Gram matrix H_ico, and verifies that
     the "+"-root signed Gram G equals
            G = P_sigma . D . H_ico . D . P_sigma^T
     for some sigma in S_5 (acting on the 10 2-subsets) and some
     switching D in {+/-1}^10. The search is exact (sympy equality):
     we iterate sigma over the 120 elements of S_5, fix D_0 = +1, and
     derive D_i = G_{sigma(0), sigma(i)} / H_ico_{0, i} for i > 0; the
     candidate (sigma, D) is accepted iff the derived equality holds
     for all 100 entries (i, j).

This direct exact search establishes the orbit-belonging claim of
Theorem 5.5 and Theorem 5.6 without invoking Theorem 2.3 (which only
gives uniqueness of the 10-line ETF in R^5 up to O(5), not the rigidity
of the signed Gram lift used here).

This script directly verifies the canonical (eps, delta) chirality
branch (both algebraic roots). The remaining 15 consistent (eps, delta)
cases of Theorem 5.5 yield the same canonical quadratic and the same
algebraic null-space family, and are identified with the canonical
branch by the discrete switching, chirality, and S_5-relabelling
symmetries recorded by the enumeration routine in
`five_tet_hinge_enum.py` (they are NOT re-verified explicitly here).
"""

from __future__ import annotations

import itertools

import sympy as sp


sq2 = sp.sqrt(2)
sq3 = sp.sqrt(3)
sq5 = sp.sqrt(5)
sq15 = sp.sqrt(15)
golden = (1 + sq5) / 2


# ----------------------------------------------------------------------
#   T_0 cube-inscribed tetrahedron axes (exact)
# ----------------------------------------------------------------------
A_sym = [
    sp.Matrix([+1, +1, +1]) / sq3,
    sp.Matrix([+1, -1, -1]) / sq3,
    sp.Matrix([-1, +1, -1]) / sq3,
    sp.Matrix([-1, -1, +1]) / sq3,
]


def perp_basis_sym(a: sp.Matrix) -> tuple[sp.Matrix, sp.Matrix]:
    v = sp.Matrix([1, 0, 0])
    e1_raw = v - (v.dot(a)) * a
    e1 = e1_raw / e1_raw.norm()
    cross = a.cross(e1)
    e2 = cross / cross.norm()
    return e1, e2


PERP_sym = [perp_basis_sym(a) for a in A_sym]


def tet_vertices_sym(k: int, eps_k, c_k, s_k):
    """Return the 4 vertices of T_k as sympy column vectors (exact)."""
    a_k = A_sym[k]
    e1, e2 = PERP_sym[k]
    radius = sp.Rational(2) * sq2 / 3
    out = [eps_k * a_k]
    out.append(-eps_k / 3 * a_k + radius * (c_k * e1 + s_k * e2))
    cos_b = -c_k / 2 - s_k * sq3 / 2
    sin_b = -s_k / 2 + c_k * sq3 / 2
    out.append(-eps_k / 3 * a_k + radius * (cos_b * e1 + sin_b * e2))
    cos_c = -c_k / 2 + s_k * sq3 / 2
    sin_c = -s_k / 2 - c_k * sq3 / 2
    out.append(-eps_k / 3 * a_k + radius * (cos_c * e1 + sin_c * e2))
    return out


# ----------------------------------------------------------------------
#   Closed-form (c_k, s_k) family from Theorem 5.4
# ----------------------------------------------------------------------
def canonical_family(s4):
    """Closed-form values of (c_1..c_4, s_1..s_4) given the free param s_4."""
    c1 = c2 = -sp.Rational(1, 4)
    c34 = sq3 * s4 - sp.Rational(1, 2)
    c3 = c4 = c34
    s1 = s2 = 2 * s4 - sq3 / 4
    s3 = s4
    return [c1, c2, c3, c4], [s1, s2, s3, s4]


# ----------------------------------------------------------------------
#   Petersen graph on C(5, 2) (Kneser realisation)
# ----------------------------------------------------------------------
KNESER = list(itertools.combinations(range(5), 2))     # 10 vertices
assert len(KNESER) == 10


def is_petersen_edge(i: int, j: int) -> bool:
    return set(KNESER[i]).isdisjoint(KNESER[j])


# ----------------------------------------------------------------------
#   Build the canonical labeling from the icosahedral ("-"-root) solution
# ----------------------------------------------------------------------
EPS_CANONICAL = [-1, -1, -1, -1]
PAIRS = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]


def vec_eq(v1, v2) -> bool:
    diff = sp.simplify(v1 - v2)
    return all(diff[i] == 0 for i in range(3))


def compute_labeling(c_vals, s_vals):
    """Find (j_k, j_l, delta) for each hinge pair from a concrete solution."""
    Tks = [tet_vertices_sym(k, EPS_CANONICAL[k], c_vals[k], s_vals[k])
           for k in range(4)]
    labeling = {}
    for (k, l) in PAIRS:
        found = None
        for j_k in (1, 2, 3):
            for j_l in (1, 2, 3):
                for d_sign in (+1, -1):
                    if vec_eq(Tks[k - 1][j_k], d_sign * Tks[l - 1][j_l]):
                        found = (j_k, j_l, d_sign)
                        break
                if found:
                    break
            if found:
                break
        labeling[(k, l)] = found
    return labeling, Tks


# The 10 Kneser axes:
# axis (0, k):   the k-vertex of T_0  (k = 1..4)        -> 4 axes
# axis (k, l):   the shared hinge vertex T_k <-> T_l    -> 6 axes
# Total: 10 axes, labelled by {0,1,2,3,4} 2-subsets.
def extract_ten_axes(Tks, labeling, T0_axes):
    """Return the 10 axes indexed in Kneser order C({0..4}, 2)."""
    axes = {}
    for k in range(1, 5):
        axes[(0, k)] = T0_axes[k - 1]
    for (k, l) in PAIRS:
        j_k, j_l, d = labeling[(k, l)]
        axes[(k, l)] = Tks[k - 1][j_k]
    out = []
    for pair in KNESER:
        out.append(axes[pair])
    return out


def signed_gram(axes):
    """Return the 10x10 signed Gram matrix G_{ij} = <v_i, v_j> (exact)."""
    n = len(axes)
    G = sp.zeros(n, n)
    for i in range(n):
        for j in range(n):
            G[i, j] = sp.nsimplify(
                sp.expand(axes[i].dot(axes[j])),
                [sq2, sq3, sq5, sq15],
                rational=False,
            )
    return G


def verify_petersen_pattern(G):
    """Check |G_{ij}|^2 == 1, 5/9, or 1/9 according to Petersen pattern."""
    n = G.shape[0]
    failures = []
    for i in range(n):
        for j in range(n):
            val_sq = sp.expand(G[i, j] ** 2)
            if i == j:
                expected = sp.Integer(1)
            elif is_petersen_edge(i, j):
                expected = sp.Rational(5, 9)
            else:
                expected = sp.Rational(1, 9)
            diff = sp.nsimplify(sp.expand(val_sq - expected), rational=True)
            if diff != 0:
                failures.append((i, j, val_sq, expected))
    return failures


def verify_squared_gram_bose_mesner(G):
    """Check Q = G entrywise squared equals 8/9 I + 4/9 A + 1/9 J."""
    n = G.shape[0]
    A = sp.zeros(n, n)
    for i in range(n):
        for j in range(n):
            if i != j and is_petersen_edge(i, j):
                A[i, j] = 1
    I = sp.eye(n)
    J = sp.ones(n, n)
    Q_expected = sp.Rational(8, 9) * I + sp.Rational(4, 9) * A + sp.Rational(1, 9) * J
    Q_actual = sp.zeros(n, n)
    for i in range(n):
        for j in range(n):
            Q_actual[i, j] = sp.simplify(G[i, j] ** 2)
    diff = sp.simplify(Q_actual - Q_expected)
    return all(diff[i, j] == 0 for i in range(n) for j in range(n))


def verify_rank(G):
    """Check rank(G) == 3 in exact arithmetic.

    We work via a 3x10 reconstruction: build the (3 x 10) matrix V whose
    columns are the axes (3D vectors), and check rank(V) <= 3 ; then
    G = V^T V automatically has rank <= 3 (= rank of V).
    """
    return G.rank()


def verify_rank_via_V(axes):
    """rank(V) where V is the 3x10 matrix of axes; rank(G) = rank(V)."""
    n = len(axes)
    V = sp.zeros(3, n)
    for j in range(n):
        for i in range(3):
            V[i, j] = sp.nsimplify(sp.expand(axes[j][i]), [sq2, sq3, sq5, sq15],
                                   rational=False)
    return V.rank()


# ----------------------------------------------------------------------
#   S_5 action on the 10 Kneser vertices
# ----------------------------------------------------------------------
def s5_action_on_kneser(tau: tuple[int, ...]) -> list[int]:
    """Given tau in S_5 (permutation of {0..4}), return the induced
    permutation of the 10 Kneser vertices: sigma(idx) where
    KNESER[idx] = (a, b) is sent to the index of (sorted) (tau(a), tau(b)).
    """
    idx_of = {pair: i for i, pair in enumerate(KNESER)}
    sigma = []
    for pair in KNESER:
        a, b = pair
        new_pair = tuple(sorted((tau[a], tau[b])))
        sigma.append(idx_of[new_pair])
    return sigma


def find_switching_for_sigma(G_sym: sp.Matrix,
                             H_ico_sym: sp.Matrix,
                             sigma: list[int]) -> list[int] | None:
    """Given sigma in S_5 acting on the 10 Kneser vertices, try to find
    D in {+/-1}^10 such that
        G[sigma(i), sigma(j)] == D[i] * D[j] * H_ico[i, j]
    for all (i, j).  Fix D[0] = +1 and derive D[i] from the (0, i) entry,
    then verify all (i, j).  Returns D (as list of +/-1) or None.
    """
    n = 10
    D = [None] * n
    D[0] = 1
    for i in range(1, n):
        Gv = G_sym[sigma[0], sigma[i]]
        Hv = H_ico_sym[0, i]
        # Need: Gv = D[0] * D[i] * Hv  =>  D[i] = Gv / Hv  (in {+/-1})
        if Hv == 0:
            # both should be zero
            if Gv != 0:
                return None
            D[i] = 1  # free choice; but since H is rank-3 PSD non-degenerate
                     # off-diagonal entries are never zero, so this branch
                     # should not occur for a Petersen ETF.
            continue
        ratio_sym = sp.nsimplify(sp.expand(Gv - Hv), rational=True)
        if ratio_sym == 0:
            D[i] = 1
        else:
            ratio_sym = sp.nsimplify(sp.expand(Gv + Hv), rational=True)
            if ratio_sym == 0:
                D[i] = -1
            else:
                return None
    # Verify D[i] D[j] H_ij == G[sigma(i), sigma(j)] for all (i, j)
    for i in range(n):
        for j in range(n):
            lhs = G_sym[sigma[i], sigma[j]]
            rhs = D[i] * D[j] * H_ico_sym[i, j]
            diff = sp.nsimplify(sp.expand(lhs - rhs), rational=True)
            if diff != 0:
                return None
    return D


def verify_orbit_equivalence(G_sym: sp.Matrix, H_ico_sym: sp.Matrix):
    """Exact search over S_5 x {+/-1}^10 (with D[0] = +1 fixed) for
    (sigma, D) such that
        G = P_sigma . D . H_ico . D . P_sigma^T.
    Returns (tau, D) on success (tau in S_5), else None.
    """
    for tau in itertools.permutations(range(5)):
        sigma = s5_action_on_kneser(tau)
        D = find_switching_for_sigma(G_sym, H_ico_sym, sigma)
        if D is not None:
            return tau, D
    return None


# ----------------------------------------------------------------------
#   Main verification: two roots of the canonical quadratic
# ----------------------------------------------------------------------
def verify_root(s4_value, root_name):
    print()
    print("=" * 70)
    print(f"Verifying root: s_4 = {root_name}  =  {s4_value}")
    print("=" * 70)

    s4 = sp.simplify(s4_value)
    c_vals, s_vals = canonical_family(s4)

    # Unit-norm sanity check
    print()
    print("Unit-norm check (c_k^2 + s_k^2 = 1):")
    for k in range(4):
        val = sp.simplify(c_vals[k] ** 2 + s_vals[k] ** 2)
        print(f"  k={k+1}: c^2 + s^2 = {val}")
        assert val == 1, f"unit-norm fails for k={k+1}"

    # Labeling extraction (from the icosahedral / "-" root)
    # Note: the labeling is the same for both roots because it is fixed
    # by the discrete chirality eps = (-1, -1, -1, -1).
    labeling, Tks = compute_labeling(c_vals, s_vals)
    if any(v is None for v in labeling.values()):
        print()
        print("  Direct matching failed for this root.")
        print("  Searching with relaxed labeling extraction...")
        s4_iso = (sq3 - sq15) / 8
        c_iso, s_iso = canonical_family(s4_iso)
        labeling, _ = compute_labeling(c_iso, s_iso)
        Tks = [tet_vertices_sym(k, EPS_CANONICAL[k], c_vals[k], s_vals[k])
               for k in range(4)]

    print()
    print("Canonical labeling (j_k, j_l, delta) per pair (k, l):")
    for p in PAIRS:
        print(f"  pair {p}: {labeling[p]}")

    axes = extract_ten_axes(Tks, labeling, A_sym)

    G = signed_gram(axes)
    print()
    print("Signed Gram matrix computed (10x10 exact).")

    # Verify Petersen pattern
    failures = verify_petersen_pattern(G)
    if failures:
        print(f"  Petersen-pattern FAILURES: {len(failures)}")
        for f in failures[:5]:
            print(f"    {f}")
    else:
        print("  Petersen-pattern check: PASS")

    # Verify Bose-Mesner identity for Q = G entrywise squared
    bm_ok = verify_squared_gram_bose_mesner(G)
    print(f"  Q == 8/9 I + 4/9 A + 1/9 J  (Proposition 2.1): {'PASS' if bm_ok else 'FAIL'}")

    # Verify rank 3 via V (3x10), automatically PSD as G = V^T V
    rk = verify_rank_via_V(axes)
    print(f"  rank(V) = {rk}  (so rank(G) = {rk}; expected 3): {'PASS' if rk == 3 else 'FAIL'}")

    success = (not failures) and bm_ok and rk == 3
    return success, G


def main():
    s4_minus = (sq3 - sq15) / 8
    s4_plus = (sq3 + sq15) / 8

    ok_minus, G_minus = verify_root(s4_minus, "(sqrt(3) - sqrt(15)) / 8")
    ok_plus, G_plus = verify_root(s4_plus, "(sqrt(3) + sqrt(15)) / 8")

    # Take the "-"-root signed Gram as the reference H_ico.
    # (This is the classical icosahedral compound by Theorem 5.4.)
    H_ico = G_minus

    print()
    print("=" * 70)
    print("Orbit-belonging certificate: G^+ vs H_ico (= G^-)")
    print("=" * 70)
    print("Searching for (sigma, D) in S_5 x {+/-1}^10 such that")
    print("   G^+ = P_sigma . D . H_ico . D . P_sigma^T")
    print("in exact sympy arithmetic. The search is over 120 elements of S_5")
    print("(with D[0] = +1 fixed; the remaining signs are forced).")

    result = verify_orbit_equivalence(G_plus, H_ico)
    if result is None:
        print()
        print("  FAILURE: no (sigma, D) found.  The '+'-root is NOT in the")
        print("  icosahedral orbit; this would contradict Theorem 5.4.")
        orbit_ok = False
    else:
        tau, D = result
        print()
        print(f"  SUCCESS: orbit equivalence verified exactly.")
        print(f"    tau in S_5 (permutation of {{0..4}}):  {tau}")
        print(f"    induced sigma on 10 Kneser vertices:  "
              f"{s5_action_on_kneser(tau)}")
        print(f"    switching D in {{+/-1}}^10:  {D}")
        orbit_ok = True

    print()
    print("=" * 70)
    print("Final summary")
    print("=" * 70)
    print(f"  '-'-root  Petersen ETF + rank-3 PSD check:  {ok_minus}")
    print(f"  '+'-root  Petersen ETF + rank-3 PSD check:  {ok_plus}")
    print(f"  Orbit equivalence (G^+ vs G^- under S_5 x switching):  {orbit_ok}")
    if ok_minus and ok_plus and orbit_ok:
        print()
        print("  Both algebraic branches of the canonical (eps, delta) case")
        print("  produce 10-axis configurations whose signed Gram matrices")
        print("  are switching- and S_5-equivalent in exact arithmetic.")
        print("  Hence both lie in the icosahedral")
        print("  O(3) x {+/-1}^10 x S_5  orbit (without invoking Theorem 2.3).")
        print()
        print("  The other 15 consistent (eps, delta) cases are identified")
        print("  with the canonical branch by the discrete switchings,")
        print("  chirality changes, and S_5-relabellings recorded by the")
        print("  enumeration routine in five_tet_hinge_enum.py; they are")
        print("  NOT re-verified explicitly by this script.")
        print()
        print("  EXACT orbit verification COMPLETE for the canonical branch.")


if __name__ == "__main__":
    main()
