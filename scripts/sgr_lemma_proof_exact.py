"""
Signed Gram Rigidity Lemma -- exact symbolic verification.

A second, exact-symbolic witness for the SGR Lemma, complementing the
numerical POCS verification in `sgr_lemma_proof.py`. NOTE: this is a
heavy optional cross-check (exact sympy backtracking over 2^3 anchor
sign choices, with full polynomial squared-Gram constraints); it can
run for several minutes. The certified proof of Theorem 3 does not
depend on this script.

Strategy:
 * Use the K(5,2) Kneser realisation of Petersen explicitly.
 * Anchor three vertices (0, 1, 2); enumerate the 2^3 = 8 sign choices on the
   anchor triangle.
 * For each remaining vertex j, find x in R^3 in the anchor coordinate
   system satisfying x^T G_3^{-1} x = 1 (the unit-norm constraint) with the
   correct |H_{i,j}| magnitudes for i in {0,1,2}.
 * Backtrack through the rest, enforcing all squared-Gram constraints
   exactly via sympy.
 * Canonicalise modulo switching (fix H_{0,j} > 0) and Aut(K(5,2)) = S_5.
 * Report raw count and orbit-class count.

Also computes the squared-Gram matrix Q = (5/9)A + (1/9)(J - I - A) + I and
shows its centered spectrum:  spec(Q - J/3) = {(4/3)^5, 0^5},  i.e. the
5D tight-frame property of the Veronese lift.
"""

import itertools
import numpy as np
import sympy as sp

sqrt5 = sp.sqrt(5)
a = sqrt5 / 3
b = sp.Rational(1, 3)

# Petersen = Kneser K(5, 2)
labels = list(itertools.combinations(range(5), 2))
n = 10
label_to_index = {lab: i for i, lab in enumerate(labels)}


def is_edge(i, j):
    return set(labels[i]).isdisjoint(labels[j])


absH = sp.eye(n)
A = np.zeros((n, n), dtype=int)
for i in range(n):
    for j in range(i + 1, n):
        val = a if is_edge(i, j) else b
        absH[i, j] = absH[j, i] = val
        if is_edge(i, j):
            A[i, j] = A[j, i] = 1


# Aut(K(5,2)) = S_5 acting on 2-subsets.
auts = []
for perm in itertools.permutations(range(5)):
    p = []
    for lab in labels:
        img = tuple(sorted((perm[lab[0]], perm[lab[1]])))
        p.append(label_to_index[img])
    auts.append(p)


def permute_matrix(H, p):
    Hp = sp.zeros(n)
    for i in range(n):
        for j in range(n):
            Hp[p[i], p[j]] = H[i, j]
    return Hp


def switch_first_row_positive(H):
    D = [1] * n
    for j in range(1, n):
        D[j] = 1 if float(H[0, j].evalf()) >= 0 else -1
    H2 = sp.zeros(n)
    for i in range(n):
        for j in range(n):
            H2[i, j] = D[i] * H[i, j] * D[j]
    return H2


def flat_key(H):
    return tuple(str(sp.simplify(H[i, j])) for i in range(n) for j in range(i + 1, n))


def canonical_key(H):
    keys = []
    for p in auts:
        Hp = switch_first_row_positive(permute_matrix(H, p))
        keys.append(flat_key(Hp))
    return min(keys)


anchors = [0, 1, 2]
rest = [i for i in range(n) if i not in anchors]


def anchor_gram(anchor_signs):
    G = sp.eye(3)
    for s, (r, c) in zip(anchor_signs, [(0, 1), (0, 2), (1, 2)]):
        i, j = anchors[r], anchors[c]
        G[r, c] = G[c, r] = s * absH[i, j]
    return sp.simplify(G)


def x_vector(j, s1, s2):
    # H[anchor0, j] is positive by switching gauge.
    return sp.Matrix([
        absH[anchors[0], j],
        s1 * absH[anchors[1], j],
        s2 * absH[anchors[2], j],
    ])


solutions = []

print("=== sgr_lemma_proof_exact.py ===", flush=True)
print(f"  heavy exact backtracking; this can take several minutes", flush=True)

for anchor_idx, anchor_signs in enumerate(itertools.product([1, -1], repeat=3)):
    print(f"  [progress] anchor_signs {anchor_idx + 1}/8 = {anchor_signs}", flush=True)
    G3 = anchor_gram(anchor_signs)
    if G3.det() == 0:
        continue
    ev = np.linalg.eigvalsh(np.array(G3.evalf(), dtype=float))
    if ev[0] <= 1e-9:
        continue

    G3inv = sp.simplify(G3.inv())

    candidates = {}
    for j in rest:
        lst = []
        for s1, s2 in itertools.product([1, -1], repeat=2):
            x = x_vector(j, s1, s2)
            norm = sp.simplify((x.T * G3inv * x)[0])
            if sp.simplify(norm - 1) == 0:
                lst.append(x)
        candidates[j] = lst

    if any(len(candidates[j]) == 0 for j in rest):
        continue

    order = sorted(rest, key=lambda j: len(candidates[j]))
    xs = {}

    def backtrack(pos):
        if pos == len(order):
            H = sp.eye(n)
            for s, (r, c) in zip(anchor_signs, [(0, 1), (0, 2), (1, 2)]):
                i, j = anchors[r], anchors[c]
                H[i, j] = H[j, i] = s * absH[i, j]
            for j in rest:
                x = xs[j]
                for r, i in enumerate(anchors):
                    H[i, j] = H[j, i] = x[r]
            for j1, j2 in itertools.combinations(rest, 2):
                val = sp.simplify((xs[j1].T * G3inv * xs[j2])[0])
                H[j1, j2] = H[j2, j1] = val
            for i in range(n):
                for j in range(n):
                    if sp.simplify(H[i, j] ** 2 - absH[i, j] ** 2) != 0:
                        return
            Hnum = np.array(H.evalf(), dtype=float)
            evals = np.linalg.eigvalsh(Hnum)
            if evals[0] < -1e-8:
                return
            if np.linalg.matrix_rank(Hnum, tol=1e-8) != 3:
                return
            solutions.append(H)
            return

        j = order[pos]
        for x in candidates[j]:
            ok = True
            for k, y in xs.items():
                val = sp.simplify((x.T * G3inv * y)[0])
                if sp.simplify(val ** 2 - absH[j, k] ** 2) != 0:
                    ok = False
                    break
            if not ok:
                continue
            xs[j] = x
            backtrack(pos + 1)
            del xs[j]

    backtrack(0)


classes = {}
for H in solutions:
    classes.setdefault(canonical_key(H), H)

rep = next(iter(classes.values()))
rep_num = np.array(rep.evalf(), dtype=float)
evals = np.linalg.eigvalsh(rep_num)

print("=== Signed Gram rigidity enumeration (exact symbolic) ===")
print("Petersen model:", "K(5, 2)")
print("vertex labels:", labels)
print("degrees:", A.sum(axis=1).tolist())
print("Aut(K(5,2)) size:", len(auts))
print("raw rank-3 PSD signed lifts:", len(solutions))
print("classes modulo switching + Aut(K(5,2)):", len(classes))
print("representative numerical eigenvalues:", np.round(evals, 12).tolist())
print("representative rank:", np.linalg.matrix_rank(rep_num, tol=1e-8))
print("representative is PSD:", bool(evals[0] > -1e-8))
print()
print("Representative H entries encoded:")
print("diag = 1")
print("edge abs = sqrt(5)/3, nonedge abs = 1/3")
print("H =")
print(rep)

# Squared Gram Q analysis.
Q = np.zeros((n, n), dtype=float)
for i in range(n):
    for j in range(n):
        if i == j:
            Q[i, j] = 1
        else:
            Q[i, j] = (5 / 9) if is_edge(i, j) else (1 / 9)

G_centered = Q - np.ones((n, n)) / 3
print()
print("Squared Gram Q eigenvalues:",
      np.round(np.linalg.eigvalsh(Q), 12).tolist())
print("Centered Gram Q - J/3 eigenvalues:",
      np.round(np.linalg.eigvalsh(G_centered), 12).tolist())
print("rank(Q):", np.linalg.matrix_rank(Q, tol=1e-9))
print("rank(Q - J/3):", np.linalg.matrix_rank(G_centered, tol=1e-9))
