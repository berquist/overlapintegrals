from pytest import approx

from sympy import integrate, symbols, oo

from overlapintegrals.sympy_exprs import pgbf_1d, pgbf


def test_pgbf_integrate():
    alpha = 0.5
    x, y, z = symbols("x y z")
    assert integrate(pgbf_1d(x, 0, 0, alpha), (x, -oo, oo)).evalf() == approx(2.50663, 1e-5)
    # Effectively do a 1D integral over a Gaussian like above.
    assert integrate(pgbf(x, 0, 0, 0, 0, 0, 0, 0, 0, alpha), (x, -oo, oo)).evalf() == approx(2.50663, 1e-5)
    pair = pgbf(x, y, z, 0, 0, 0, 0, 0, 0, 1.8) * pgbf(x, y, z, 0.5, 0.8, -0.2, 0, 0, 0, 2.8)
    assert integrate(integrate(integrate(pair, (x, -oo, oo)), (y, -oo, oo)), (z, -oo, oo)).evalf() == approx(0.20373275913014607, 1e-16)
