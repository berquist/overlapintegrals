from functools import reduce
from math import exp, floor, pi
from typing import Tuple, Union

from scipy.special import binom as binomial


Number = Union[int, float]


def dist2(
    xa: Number, ya: Number, za: Number, xb: Number, yb: Number, zb: Number
) -> float:
    return (xa - xb) * (xa - xb) + (ya - yb) * (ya - yb) + (za - zb) * (za - zb)


def product_center_1d(za: float, xa: float, zb: float, xb: float) -> float:
    return ((za * xa) + (zb * xb)) / (za + zb)


def fact2(n):
    """Compute the double factorial

    as in https://stackoverflow.com/a/4740229/3249688
    """
    return reduce(int.__mul__, range(n, 0, -2), 1)


def binomial_prefactor(s: int, ia: int, ib: int, xpa: float, xpb: float) -> float:
    return sum(
        binomial(ia, s - t)
        * binomial(ib, t)
        * (xpa ** (ia - s + t))
        * (xpb ** (ib - t))
        for t in range(s + 1)
        if (s - ia) <= t and t <= ib
    )


def overlap1d(l1: int, l2: int, pax: float, pbx: float, gamma: float) -> float:
    return sum(
        binomial_prefactor(2 * i, l1, l2, pax, pbx)
        * fact2(2 * i - 1)
        / ((2 * gamma) ** i)
        for i in range(1 + floor(0.5 * (l1 + l2)))
    )


def tho66(
    alpha1: float,
    alpha2: float,
    ra: Tuple[float, float, float],
    rb: Tuple[float, float, float],
    la: Tuple[float, float, float],
    lb: Tuple[float, float, float],
) -> float:
    gamma = alpha1 + alpha2
    xa, ya, za = ra
    xb, yb, zb = rb
    l1, m1, n1 = la
    l2, m2, n2 = lb
    xp = product_center_1d(alpha1, xa, alpha2, xb)
    yp = product_center_1d(alpha1, ya, alpha2, yb)
    zp = product_center_1d(alpha1, za, alpha2, zb)

    pre = exp(-alpha1 * alpha2 * dist2(xa, ya, za, xb, yb, zb) / gamma) * (pi / gamma) ** 1.5

    wx = overlap1d(l1, l2, xp - xa, xp - xb, gamma)
    wy = overlap1d(m1, m2, yp - ya, yp - yb, gamma)
    wz = overlap1d(n1, n2, zp - za, zp - zb, gamma)
    return pre * wx * wy * wz
