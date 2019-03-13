from pytest import approx, raises

import overlapintegrals


def test_dist2():
    assert overlapintegrals.dist2(0.5, 0.6, 0.7, 0.8, 0.9, 1.0) == approx(0.27, 1e-5)


def test_product_center_1d():
    assert overlapintegrals.product_center_1d(2, 3, 4, 5) == 26 / 6
    with raises(ZeroDivisionError):
        overlapintegrals.product_center_1d(0, 3, 0, 5)
    alpha1 = 1.8
    alpha2 = 2.8
    xa = 0.0
    xb = 0.5
    assert overlapintegrals.product_center_1d(alpha1, xa, alpha2, xb) == approx(0.304348, 1e-6)


def test_fact2():
    # this one is by definition
    assert overlapintegrals.fact2(0) == 1
    assert overlapintegrals.fact2(1) == 1
    assert overlapintegrals.fact2(2) == 2
    assert overlapintegrals.fact2(3) == (3 * 1)
    assert overlapintegrals.fact2(4) == (4 * 2)
    assert overlapintegrals.fact2(5) == (5 * 3 * 1)
    assert overlapintegrals.fact2(6) == (6 * 4 * 2)
    assert overlapintegrals.fact2(7) == (7 * 5 * 3 * 1)


def test_binomial_prefactor():
    assert overlapintegrals.binomial_prefactor(s=1, ia=1, ib=1, xpa=0.1, xpb=0.2) == approx(0.3, 1e-7)
    assert overlapintegrals.binomial_prefactor(s=1, ia=1, ib=1, xpa=0.3, xpb=0.4) == approx(0.7, 1e-7)
    assert overlapintegrals.binomial_prefactor(s=1, ia=3, ib=1, xpa=0.1, xpb=0.2) == approx(0.007, 1e-7)
    assert overlapintegrals.binomial_prefactor(s=2, ia=3, ib=1, xpa=0.1, xpb=0.2) == approx(0.09, 1e-7)


def test_overlap1d():
    assert overlapintegrals.overlap1d(l1=1, l2=1, pax=0.1, pbx=0.2, gamma=1.0) == approx(0.52, 1e-5)


def test_tho66():
    za = 1.8
    zb = 2.8
    ra = [0.0, 0.0, 0.0]
    rb = [0.5, 0.8, -0.2]

    integral = overlapintegrals.tho66(za, zb, ra, rb, [0, 0, 0], [0, 0, 0])
    assert integral == approx(0.20373275913014607, 1.0e-16)

    integral = overlapintegrals.tho66(za, zb, ra, rb, [1, 0, 0], [0, 0, 0])
    assert integral == approx(0.062005622343957505, 1.0e-16)

    integral = overlapintegrals.tho66(za, zb, ra, rb, [1, 1, 0], [1, 1, 0])
    assert integral == approx(-0.00043801221837779696, 1.0e-16)

    integral = overlapintegrals.tho66(za, zb, ra, rb, [2, 1, 0], [1, 1, 0])
    assert integral == approx(-0.0002385994651113168, 1.0e-16)
