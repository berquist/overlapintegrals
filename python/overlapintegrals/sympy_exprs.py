from sympy import exp  #, integrate, symbols


def pgbf_1d(rx, ax, l, alpha):
    xa = rx - ax
    return (xa ** l) * exp(-alpha * xa ** 2)


def pgbf(rx, ry, rz, ax, ay, az, l, m, n, alpha):
    return pgbf_1d(rx, ax, l, alpha) * pgbf_1d(ry, ay, l, alpha) * pgbf_1d(rz, az, l, alpha)
