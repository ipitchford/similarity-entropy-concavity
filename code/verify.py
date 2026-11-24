"""Exact four-state coefficient certificate and five-state curvature check.

Run with Python 3.9+ and the standard library: python verify.py
"""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import comb


def check(condition, message='check failed'):
    """Explicit check that also fails under python -O."""
    if not condition:
        raise SystemExit('CHECK FAILED: ' + str(message))


# A cut is encoded by the indicators at sites 0,1,2,3. Complements
# represent the same cut, so the indicator at site 0 is always zero.
CUTS = [cut for cut in product((0, 1), repeat=4)
        if cut[0] == 0 and any(cut) and cut != (0, 0, 1, 1)]
ZERO = (0,) * len(CUTS)


def add_exponents(left, right, sign=1):
    return tuple(a + sign * b for a, b in zip(left, right))


def separation(i, j):
    return tuple(abs(cut[i] - cut[j]) for cut in CUTS)


def entry(i, j):
    """Laurent coefficients of the root Schur complement in b=exp(weight)."""
    ri, rj, dij = separation(0, i), separation(0, j), separation(i, j)
    terms = Counter()
    terms[add_exponents(ri, dij, -1)] += 1
    terms[add_exponents(rj, dij, -1)] += 1
    terms[tuple(-s for s in add_exponents(ri, rj))] -= 1
    terms[ZERO] -= 1
    return terms


def determinant(indices):
    result = Counter()
    for perm in permutations(indices):
        inversions = sum(x > y for x, y in combinations(perm, 2))
        terms = {ZERO: (-1) ** inversions}
        for i, j in zip(indices, perm):
            expanded = Counter()
            for left, left_coeff in terms.items():
                for right, right_coeff in entry(i, j).items():
                    expanded[add_exponents(left, right)] += left_coeff * right_coeff
            terms = expanded
        result.update(terms)
    return {power: coeff for power, coeff in result.items() if coeff}


def translate_to_orthant(terms):
    """Remove the common Laurent monomial and substitute b_l=1+t_l."""
    minima = tuple(map(min, zip(*terms)))
    current = {add_exponents(power, minima, -1): coeff
               for power, coeff in terms.items()}
    for axis in range(len(CUTS)):
        new = Counter()
        for power, coeff in current.items():
            exponent = power[axis]
            for j in range(exponent + 1):
                key = power[:axis] + (j,) + power[axis + 1:]
                new[key] += comb(exponent, j) * coeff
        current = {power: coeff for power, coeff in new.items() if coeff}
    return current


def check_four():
    expected = {
        (1,): (255, 2, 161), (2,): (63, 2, 53), (3,): (63, 2, 53),
        (1, 2): (1557, 4, 3672), (1, 3): (1557, 4, 3672),
        (2, 3): (1245, 4, 1836), (1, 2, 3): (15172, 8, 464400),
    }
    print('Four-state Schur complement: index set, support, minimum, maximum')
    for indices, reference in expected.items():
        coeffs = translate_to_orthant(determinant(indices)).values()
        check(all(coefficient > 0 for coefficient in coeffs), 'all(coefficient > 0 for coefficient in coeffs)')
        stats = (len(coeffs), min(coeffs), max(coeffs))
        check(stats == reference, (indices, stats))
        print(indices, *stats)


def check_five():
    a = F(7, 10)
    n = 5
    side = (0, 0, 1, 1, 1)
    kernel = [[1 if i == j else a if side[i] != side[j] else a**2
               for j in range(n)] for i in range(n)]
    for i, j, k in product(range(n), repeat=3):
        check(kernel[i][j] * kernel[j][k] <= kernel[i][k], 'kernel[i][j] * kernel[j][k] <= kernel[i][k]')
    check((1 - a*a) * (1 - 2*a*a) > 0, '(1 - a*a) * (1 - 2*a*a) > 0')
    p = [F(36, 40)] + [F(1, 40)] * 4
    v = [2, 1, -1, -1, -1]

    def multiply(z):
        return [sum(x*y for x, y in zip(row, z)) for row in kernel]

    q, w = multiply(p), multiply(v)
    curvature = sum(pi*(wi/qi)**2 - 2*vi*wi/qi
                    for pi, vi, wi, qi in zip(p, v, w, q))
    check(curvature == F(2748940385160, 93148730546881), 'curvature == F(2748940385160, 93148730546881)')
    h = F(1, 1000)
    check(all(pi > h*abs(vi) for pi, vi in zip(p, v)), 'all(pi > h*abs(vi) for pi, vi in zip(p, v))')
    check(all(qi - h*abs(wi) > F(1, 2) for qi, wi in zip(q, w)), 'all(qi - h*abs(wi) > F(1, 2) for qi, wi in zip(q, w))')
    ratio_bound = F(4, 5)
    check(all(2*abs(wi) <= ratio_bound for wi in w), 'all(2*abs(wi) <= ratio_bound for wi in w)')
    third_bound = 3 * sum(map(abs, v))*ratio_bound**2 + 2*ratio_bound**3
    gap_lower = h**2*(curvature - h*third_bound)/2
    check(third_bound < 13 and gap_lower > F(8, 10**9), 'third_bound < 13 and gap_lower > F(8, 10**9)')
    print('Five-state second derivative of H_K along p+tv (a=7/10):', curvature)
    print('Two-posterior mean entropy exceeds prior entropy by > 8e-9 nats.')


if __name__ == '__main__':
    check_four()
    check_five()
