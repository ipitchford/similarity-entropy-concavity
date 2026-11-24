"""Decide concavity of a rational order-one similarity entropy.

Usage: python classify.py kernel.json
Input is a square symmetric JSON array with entries in [0,1] and diagonal 1.
Use integer entries or quoted fractions such as "7/10" for exact inputs.
JSON decimal numbers are also read exactly, through their literal text.
"""
import argparse
import json
from fractions import Fraction
from pathlib import Path


def positive_semidefinite(matrix, strict=False):
    work = [row[:] for row in matrix]
    for i in range(len(work)):
        pivot = work[i][i]
        if pivot < 0:
            return False
        if pivot == 0:
            if strict or any(work[i][j] for j in range(i + 1, len(work))):
                return False
            continue
        for j in range(i + 1, len(work)):
            for k in range(j, len(work)):
                work[j][k] -= work[j][i] * work[i][k] / pivot
                work[k][j] = work[j][k]
    return True


def invert_positive_definite(matrix):
    n = len(matrix)
    rows = [row[:] + [Fraction(j == i) for j in range(n)]
            for i, row in enumerate(matrix)]
    for i in range(n):
        pivot = rows[i][i]
        rows[i] = [value / pivot for value in rows[i]]
        for j in range(n):
            if j == i:
                continue
            multiple = rows[j][i]
            rows[j] = [v - multiple * u for v, u in zip(rows[j], rows[i])]
    return [row[n:] for row in rows]


def classify(kernel):
    matrix = [[Fraction(v) for v in row] for row in kernel]
    n = len(matrix)
    if n == 0 or any(len(row) != n for row in matrix):
        raise ValueError('Supply a nonempty square matrix')
    for i, row in enumerate(matrix):
        if row[i] != 1:
            raise ValueError('The diagonal must be one')
        for j, value in enumerate(row):
            if not 0 <= value <= 1 or value != matrix[j][i]:
                raise ValueError('The similarity must be symmetric and in [0,1]')
    groups = {}
    for i, row in enumerate(matrix):
        groups.setdefault(tuple(row), []).append(i)
    classes = list(groups.values())
    representatives = [members[0] for members in classes]
    reduced = [[matrix[i][j] for j in representatives] for i in representatives]
    report = {'row_classes_zero_based': classes, 'concave': False,
              'strictly_concave': False}
    if not positive_semidefinite(reduced, strict=True):
        report['reason'] = 'Reduced kernel is not positive definite'
        return report
    inverse = invert_positive_definite(reduced)
    m = len(reduced)
    for k in range(m):
        test = [[inverse[i][j] * (reduced[i][k] + reduced[j][k])
                 - (i == j == k) for j in range(m)] for i in range(m)]
        if not positive_semidefinite(test):
            report['reason'] = 'Failed vertex positive-semidefiniteness test'
            report['failed_vertex_zero_based'] = representatives[k]
            return report
    report['concave'] = True
    report['strictly_concave'] = m == n
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('kernel', type=Path)
    input_matrix = json.loads(parser.parse_args().kernel.read_text(), parse_float=str)
    print(json.dumps(classify(input_matrix), indent=2))
