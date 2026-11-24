"""Exact checks of the dimension-five experiment and a one-bit tree example."""
from fractions import Fraction as F


def check(condition, message='check failed'):
    """Explicit check that also fails under python -O."""
    if not condition:
        raise SystemExit('CHECK FAILED: ' + str(message))


def two_parts():
    a = F(12, 17)
    A = (1 + a*a)/2
    B = (1 + 2*a*a)/3
    ratio = (A+a)*(B+a)/(4*A*B)
    check(1-2*a*a == F(1, 289) and ratio == F(999949, 999364), '1-2*a*a == F(1, 289) and ratio == F(999949, 999364)')
    print('Five-state twice-exponentiated Jensen gap:', ratio, '> 1')


def tree_example():
    m, a, e = 1024, F(1, 32), F(1, 4)
    b, l = F(1, m*m), 1-e
    c = (1+(m-1)*a*a)/m
    posterior_hub = e+l*a*b
    posterior_leaf = l*c+e*a*b
    prior_hub = (1+b)*(e+l*a)/2
    prior_leaf = (l*(c+a*a*b)+e*a*(1+b))/2
    hub_ratio = prior_hub/posterior_hub
    leaf_ratio = prior_leaf/posterior_leaf
    check(hub_ratio == F(7340039, 13421774), 'hub_ratio == F(7340039, 13421774)')
    check(leaf_ratio == F(39842851, 12576832), 'leaf_ratio == F(39842851, 12576832)')
    check(hub_ratio*leaf_ratio**3 > 16, 'hub_ratio*leaf_ratio**3 > 16')
    print('Double star on', 2*m+2, 'vertices. Exact exp(4G) > 16, hence G > log(2).')


if __name__ == '__main__':
    two_parts()
    tree_example()
