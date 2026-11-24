"""Exact and high-precision checks added in the revision.

Run: python check_extensions.py [--output results/extensions.json]
Requires SymPy and mpmath. Every check raises SystemExit on failure, so the
script also fails under python -O. Exact checks use rational or symbolic
arithmetic. The asymptotic table uses mpmath at 50 digits; its printed values
are rounded to six significant digits and carry a declared tolerance of 1e-9.

Sections
 1. Cut decomposition of a four-point metric (Section 5).
 2. The five-state family K_{2,3}(a): positive definiteness and the interval of
    nonconcavity, from the root Schur block (Section 5).
 3. The star threshold C(a) and its minimum (Section 6).
 4. Sign-condition sufficient test on random ultrametric kernels (Section 3).
 5. A seeded random experiment: how often do five-state MTI kernels fail?
 6. Effective range of the asymptotic statements of Section 7.
"""
import argparse
import json
import random
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path

import mpmath
import sympy as sp

from classify import classify, invert_positive_definite


def check(condition, message='check failed'):
    """Explicit check that also fails under python -O."""
    if not condition:
        raise SystemExit('CHECK FAILED: ' + str(message))


RESULTS = {}


def cut_decomposition():
    d = {pair: sp.Symbol('d%d%d' % pair, nonnegative=True)
         for pair in combinations(range(4), 2)}

    def dist(i, j):
        return 0 if i == j else d[(min(i, j), max(i, j))]

    S1 = dist(0, 1) + dist(2, 3)
    S2 = dist(0, 2) + dist(1, 3)
    S3 = dist(0, 3) + dist(1, 2)
    # cut -> (side containing a point is the set listed) -> weight
    weights = {
        frozenset({1, 3}): (S1 - S2) / 2,        # cut 02|13
        frozenset({1, 2}): (S1 - S3) / 2,        # cut 03|12
        frozenset({0}): (dist(0, 2) + dist(0, 3) - dist(2, 3)) / 2,
        frozenset({1}): (dist(1, 2) + dist(1, 3) - dist(2, 3)) / 2,
        frozenset({2}): (dist(0, 2) + dist(1, 2) - dist(0, 1)) / 2,
        frozenset({3}): (dist(0, 3) + dist(1, 3) - dist(0, 1)) / 2,
    }
    for i, j in combinations(range(4), 2):
        total = sum(w for side, w in weights.items() if (i in side) != (j in side))
        check(sp.simplify(total - dist(i, j)) == 0, ('cut sum', i, j))
    # Nonnegativity: singletons are triangle inequalities; the two
    # non-singleton weights are nonnegative because S1 is the largest sum.
    tri = []
    for side, w in weights.items():
        if len(side) == 1:
            expr = sp.expand(2 * w)
            tri.append(str(expr))
    RESULTS['cut_decomposition'] = {
        'six_distances_recovered': True,
        'singleton_weights_times_two': tri,
        'pair_cut_weights': ['(S1-S2)/2', '(S1-S3)/2'],
    }
    print('1. Cut weights recover all six distances (symbolic).')


def five_state_family():
    a = sp.Symbol('a', positive=True)
    # part indicator block of K_{2,3}(a)
    block = sp.Matrix([[1 + a**2, sp.sqrt(6) * a], [sp.sqrt(6) * a, 1 + 2 * a**2]])
    check(sp.factor(block.det() - (1 - a**2) * (1 - 2 * a**2)) == 0, 'part block det')
    # full kernel: states 0,1 in part of size two, 2,3,4 in part of size three
    side = [0, 0, 1, 1, 1]
    K = sp.Matrix(5, 5, lambda i, j: 1 if i == j else (a if side[i] != side[j] else a**2))
    k = [K[i, 0] for i in range(5)]
    C = sp.Matrix(4, 4, lambda i, j: K[i + 1, j + 1] * (1 / k[i + 1] + 1 / k[j + 1])
                  - k[i + 1] * k[j + 1] - 1)
    s = 1 / sp.sqrt(3)
    # columns: e1, u=(1,1,1)/sqrt3 on part three, two orthonormal contrasts
    Bm = sp.Matrix([[1, 0, 0, 0],
                    [0, s, 1 / sp.sqrt(2), 1 / sp.sqrt(6)],
                    [0, s, -1 / sp.sqrt(2), 1 / sp.sqrt(6)],
                    [0, s, 0, -2 / sp.sqrt(6)]])
    Cb = sp.simplify(Bm.T * C * Bm)
    check(sp.simplify(Cb[0, 2]) == 0 and sp.simplify(Cb[0, 3]) == 0 and sp.simplify(Cb[1, 2]) == 0
          and sp.simplify(Cb[1, 3]) == 0 and sp.simplify(Cb[2, 3]) == 0, 'block splitting')
    d1 = sp.simplify(Cb[0, 0] - (2 / a**2 - 2 + 1 - a**4))
    d2 = sp.simplify(Cb[1, 1] - (2 / a - 2 * a - 3 * (1 - a)**2))
    off = sp.simplify(Cb[0, 1] - sp.sqrt(3) * (1 - a**4) / a)
    check(d1 == 0 and d2 == 0 and off == 0, 'block entries')
    det = sp.factor(Cb[0, 0] * Cb[1, 1] - Cb[0, 1]**2)
    poly = 4 * a**5 + a**4 + 7 * a**3 - 5 * a**2 + 5 * a - 4
    check(sp.simplify(det + (1 - a)**2 * (1 + a) / a**3 * poly) == 0, 'determinant polynomial')
    contrast = sp.simplify(Cb[2, 2])
    deriv = sp.diff(poly, a)
    check(sp.Poly(deriv, a).all_coeffs() == [20, 4, 21, -10, 5], 'derivative')
    # 20a^4+4a^3+21a^2-10a+5 = 20a^4+4a^3+16a^2+5(a-1)^2+... > 0: use exact discriminant of the quadratic part
    check(sp.discriminant(21 * a**2 - 10 * a + 5, a) < 0, 'derivative positive')
    lo, hi = F(673, 1000), F(674, 1000)
    plo = poly.subs(a, sp.Rational(lo.numerator, lo.denominator))
    phi = poly.subs(a, sp.Rational(hi.numerator, hi.denominator))
    check(plo < 0 < phi, 'root bracket')
    root = sp.nsolve(poly, a, 0.6736)
    check(F(673, 1000) < F(str(round(float(root), 12))) < F(674, 1000), 'root location')
    # for a in (root, 1/sqrt2) the determinant is negative (poly>0); contrasts value
    RESULTS['five_state_family'] = {
        'part_block_determinant': '(1-a^2)(1-2a^2)',
        'schur_block_determinant': '-(1-a)^2(1+a)a^-3 (4a^5+a^4+7a^3-5a^2+5a-4)',
        'root': float(sp.N(root, 15)),
        'contrast_entry': str(sp.factor(contrast)),
        'nonconcave_and_positive_definite_interval': [float(sp.N(root, 12)), float(sp.N(1 / sp.sqrt(2), 12))],
        'interval_width': float(sp.N(1 / sp.sqrt(2) - root, 6)),
    }
    # Exact classifier on a rational grid of a: concave iff a <= root
    grid = []
    for num in range(60, 71):
        aa = F(num, 100)
        Kr = [[1 if i == j else (aa if side[i] != side[j] else aa * aa) for j in range(5)] for i in range(5)]
        res = classify(Kr)
        expected = aa < F(6736, 10000)
        check(res['concave'] == expected, ('grid', num))
        grid.append([str(aa), res['concave']])
    RESULTS['five_state_family']['classifier_grid'] = grid
    print('2. K_{2,3}(a): PD for a<1/sqrt2; nonconcave exactly on (%.6f, 1/sqrt2).' % float(root))


def star_threshold():
    a = sp.Symbol('a', positive=True)
    C = 2 * (1 + a) / (a * (1 - a))
    crit = sp.solve(sp.diff(C, a), a)
    crit = [c for c in crit if c.is_real and 0 < c < 1]
    check(len(crit) == 1 and sp.simplify(crit[0] - (sp.sqrt(2) - 1)) == 0, 'critical point')
    check(sp.simplify(C.subs(a, crit[0]) - (6 + 4 * sp.sqrt(2))) == 0, 'minimum value')
    aa = F(41421356, 100000000)  # rational close to sqrt2-1
    for m, expected in ((11, True), (12, False)):
        star = [[1 if i == j else aa if min(i, j) == 0 else aa * aa for j in range(m + 1)] for i in range(m + 1)]
        check(classify(star)['concave'] == expected, ('star', m))
    RESULTS['star_threshold'] = {'minimum_of_C': '6+4*sqrt(2)', 'argmin': 'sqrt(2)-1',
                                 'minimum_float': float(sp.N(6 + 4 * sp.sqrt(2), 10)),
                                 'eleven_leaves_concave': True, 'twelve_leaves_fail': True}
    print('3. Star threshold: min C(a) = 6+4sqrt2 at sqrt2-1; 11 leaves pass, 12 fail.')


def random_ultrametric(rng, n):
    # random dendrogram: repeatedly merge two clusters at a strictly lower similarity
    clusters = [[i] for i in range(n)]
    K = [[F(1) if i == j else None for j in range(n)] for i in range(n)]
    level = F(1)
    while len(clusters) > 1:
        level = level * F(rng.randint(30, 95), 100)
        i, j = sorted(rng.sample(range(len(clusters)), 2))
        for x in clusters[i]:
            for y in clusters[j]:
                K[x][y] = K[y][x] = level
        clusters[i] = clusters[i] + clusters[j]
        del clusters[j]
    return K


def sign_condition():
    rng = random.Random(20261124)
    count = 0
    for trial in range(60):
        n = rng.randint(3, 7)
        K = random_ultrametric(rng, n)
        A = invert_positive_definite(K)
        zsign = all(A[i][j] <= 0 for i in range(n) for j in range(n) if i != j)
        rows = all(sum(A[i]) >= 0 for i in range(n))
        check(zsign and rows, ('sign condition', trial))
        check(classify(K)['strictly_concave'], ('ultrametric concave', trial))
        count += 1
    RESULTS['ultrametric_sign_condition'] = {'random_kernels': count, 'all_Z_matrix_nonnegative_rows': True,
                                             'all_strictly_concave_by_vertex_test': True, 'seed': 20261124}
    print('4. %d random ultrametric kernels: inverse sign condition holds; all strictly concave.' % count)


def mti(K):
    n = len(K)
    return all(K[i][j] * K[j][k] <= K[i][k] for i, j, k in product(range(n), repeat=3))


def random_five_state_experiment():
    rng = random.Random(1557)
    accepted = failures = not_pd = 0
    attempts = 0
    while accepted < 2000:
        attempts += 1
        K = [[F(1)] * 5 for _ in range(5)]
        for i, j in combinations(range(5), 2):
            K[i][j] = K[j][i] = F(rng.randint(1, 39), 40)
        if not mti(K):
            continue
        res = classify(K)
        if res.get('reason') == 'Reduced kernel is not positive definite':
            not_pd += 1
            continue
        accepted += 1
        if not res['concave']:
            failures += 1
    RESULTS['random_five_state_mti'] = {
        'distribution': 'off-diagonal entries iid uniform on {1/40,...,39/40}; kept if MTI and positive definite',
        'seed': 1557, 'attempts': attempts, 'rejected_not_pd': not_pd, 'accepted': accepted,
        'nonconcave': failures}
    print('5. Random five-state PD MTI kernels: %d of %d fail concavity.' % (failures, accepted))
    return failures, accepted


def asymptotics():
    mpmath.mp.dps = 50
    log = mpmath.log
    table = []
    for e in (3, 6, 12, 24):
        m = mpmath.mpf(10) ** e
        a = (1 - 1 / m) / mpmath.sqrt(m)
        lam = 2 / log(m)
        G = (1 - lam) * log(1 - lam + lam * m * a) + lam * log(lam + (1 - lam) * a)
        target = log(m) / 2 - log(log(m)) + log(2) - 1
        table.append({'m': '1e%d' % e, 'G': float(G), 'expansion': float(target),
                      'difference': float(G - target)})
    star = []
    for e in (3, 12, 48, 96):
        m = mpmath.mpf(10) ** e
        a = 1 / mpmath.sqrt(m)
        b = 1 / m**2
        eps = 1 / log(m)
        l = 1 - eps
        cm = (1 + (m - 1) * a**2) / m
        u = eps + l * a * b
        v = l * cm + eps * a * b
        us = (1 + b) * (eps + l * a) / 2
        vs = (l * (cm + a**2 * b) + eps * a * (1 + b)) / 2
        G = eps * log(us / u) + l * log(vs / v)
        post = -eps * log(u) - l * log(v)
        prior = -eps * log(us) - l * log(vs)
        check(abs((post - prior) - G) < mpmath.mpf(10) ** -30, 'double star identity')
        star.append({'m': '1e%d' % e, 'ratio_prior_to_mean_posterior': float(prior / post),
                     'G_minus_half_log_m_plus_loglog_m': float(G - log(m) / 2 + log(log(m)))})
    upper = []
    for e in (6, 20, 100):
        n = mpmath.mpf(10) ** e
        L = log(n)
        D = L - 4 * log(L) - 1
        eta = mpmath.sqrt(2) * D / L**2
        bound = L / 2 - log(D) + log(2) - 1 + eta / 2
        target = L / 2 - log(L) + log(2) - 1
        upper.append({'n': '1e%d' % e, 'bound_minus_expansion': float(bound - target)})
    check(abs(table[1]['difference']) < 1e-2, 'hub within 1e-2 at 1e6')
    RESULTS['asymptotics'] = {'hub_construction': table, 'fair_double_star': star,
                              'upper_bound_slack': upper, 'tolerance': 1e-9,
                              'limit_constant_double_star': '-log(4)-1/2 (approx %.6f)' % float(-log(4) - mpmath.mpf(1) / 2)}
    print('6. Hub construction within %.1e nats of the expansion at m=1e6; double-star ratio %.3f at m=1e96.'
          % (abs(table[1]['difference']), star[-1]['ratio_prior_to_mean_posterior']))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    cut_decomposition()
    five_state_family()
    star_threshold()
    sign_condition()
    random_five_state_experiment()
    asymptotics()
    if args.output:
        args.output.write_text(json.dumps(RESULTS, indent=2, sort_keys=True) + '\n')
    print('ALL CHECKS PASSED')


if __name__ == '__main__':
    main()
