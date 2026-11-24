"""Regression and formula checks for the rational finite-kernel classifier."""
import json
from fractions import Fraction as F
from pathlib import Path
from classify import classify, positive_semidefinite


def check(condition, message='check failed'):
    """Explicit check that also fails under python -O."""
    if not condition:
        raise SystemExit('CHECK FAILED: ' + str(message))


def main():
    witness = json.loads(Path(__file__).with_name('five-state.json').read_text())
    check(not classify(witness)['concave'], "not classify(witness)['concave']")
    for removed in range(5):
        keep = [i for i in range(5) if i != removed]
        sub = [[witness[i][j] for j in keep] for i in keep]
        check(classify(sub)['strictly_concave'], "classify(sub)['strictly_concave']")
    for n in range(1, 6):
        identity = [[int(i == j) for j in range(n)] for i in range(n)]
        check(classify(identity)['strictly_concave'], "classify(identity)['strictly_concave']")
        constants = classify([[1] * n for _ in range(n)])
        check(constants['concave'] and constants['strictly_concave'] == (n == 1), "constants['concave'] and constants['strictly_concave'] == (n == 1)")
    check(not positive_semidefinite([[F(0), F(1)], [F(1), F(0)]]), 'not positive_semidefinite([[F(0), F(1)], [F(1), F(0)]])')
    check(positive_semidefinite([[F(0), F(0)], [F(0), F(1)]]), 'positive_semidefinite([[F(0), F(0)], [F(0), F(1)]])')

    singular = [[1 if i == j else F(1, 2) if (i < 3) != (j < 3) else F(1, 4)
                 for j in range(6)] for i in range(6)]
    check(not classify(singular)['concave'], "not classify(singular)['concave']")

    checks = 0
    for a in (F(1, 4), F(1, 3), F(2, 5), F(1, 2), F(2, 3)):
        threshold = 2 * (1 + a) / (a * (1 - a))
        floor = threshold.numerator // threshold.denominator
        for m in (floor - 1, floor, floor + 1):
            star = [[1 if i == j else a if min(i, j) == 0 else a*a
                     for j in range(m + 1)] for i in range(m + 1)]
            expected = m <= threshold
            result = classify(star)
            check((result['concave'], result['strictly_concave']) == (expected, expected), "(result['concave'], result['strictly_concave']) == (expected, expected)")
            checks += 1
    kernel = [[F(q) for q in row] for row in witness]
    p = [F(36, 40)] + [F(1, 40)] * 4
    z = [6, 9, -5, -5, -5]

    def third_order_value(law):
        typicality = [sum(x*y for x, y in zip(row, law)) for row in kernel]
        return (1 - sum(pi * q**2 for pi, q in zip(law, typicality))) / 2

    plus = [pi + F(zi, 1000) for pi, zi in zip(p, z)]
    minus = [pi - F(zi, 1000) for pi, zi in zip(p, z)]
    gap = (third_order_value(plus) + third_order_value(minus))/2 - third_order_value(p)
    check(gap == F(1557, 50000000000), 'gap == F(1557, 50000000000)')
    print('Classifier checks passed, including', checks, 'star threshold comparisons.')
    print('Ricotta–Szeidl order-three five-state Jensen gap:', gap)


if __name__ == '__main__':
    main()
