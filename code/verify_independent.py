"""Independent polynomial check in the SymPy integer polynomial ring.

Compute a diagonally congruent Schur matrix directly in shifted variables.
This route has no Laurent expansions, translation passes or determinant
permutation code. Install SymPy and run python verify_independent.py.
"""
from itertools import combinations
from math import prod
from sympy.polys.rings import ring
from sympy.polys.domains import ZZ

def run():
    R, *variables = ring("u0,u1,u2,u3,u4,u5", ZZ)
    # Singleton splits followed by 02|13 and 03|12. Each selected side
    # is the side away from the root 0.
    sides = [{1,2,3}, {1}, {2}, {3}, {1,3}, {1,2}]
    factors = [1 + t for t in variables]
    idx = [1,2,3]
    def weight(i):
        return prod((b for b,s in zip(factors,sides) if i in s),start=R.one)
    def entry(i,j):
        # w_i=1/K_i0. The matrix entry is w_i w_j C_ij.
        overlap = prod((b*b for b,s in zip(factors,sides) if i in s and j in s),start=R.one)
        wi,wj=weight(i),weight(j)
        return overlap*(wi+wj)-wi*wj-1
    M = {(i,j):entry(i,j) for i in idx for j in idx}
    a,b,c = (M[i,i] for i in idx)
    d,e,f = M[1,2],M[1,3],M[2,3]
    polynomials = [a,b,c,a*b-d*d,a*c-e*e,b*c-f*f,
                   a*(b*c-f*f)-d*(d*c-e*f)+e*(d*f-b*e)]
    for subset,poly in zip([*combinations(idx,1),*combinations(idx,2),tuple(idx)],polynomials):
        values = list(poly.values())
        if not values or min(values) < 0:
            raise SystemExit("CHECK FAILED: " + f"negative or missing coefficient for {subset}")
        print(subset,len(values),min(values),max(values))

if __name__=='__main__':
    run()
