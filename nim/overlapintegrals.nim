import math
import ./test_utils
import unittest

func dist2(xa, ya, za, xb, yb, zb: SomeFloat): SomeFloat =
  pow(xa - xb, 2.0) + pow(ya - yb, 2.0) + pow(za - zb, 2.0)

func productCenter1d(za, xa, zb, xb: SomeNumber): SomeFloat =
   ((za * xa) + (zb * xb)) / (za + zb)

func fact2(n: SomeInteger): SomeInteger =
  result = 1
  var counter = n
  while counter > 0:
    result *= counter
    counter -= 2

func binomialPrefactor(s, ia, ib: SomeInteger, xpa, xpb: SomeFloat): SomeFloat =
  for t in 0..s:
    if (s - ia) <= t and (t <= ib):
      result += float(binom(ia, s - t) * binom(ib, t)) * pow(xpa, float(ia - s + t)) * pow(xpb, float(ib - t))

func overlap1d(l1, l2: SomeInteger, pax, pbx, gamma: SomeFloat): SomeFloat =
  for i in 0..int(floor(0.5 * float(l1 + l2))):
    result += binomialPrefactor(2 * i, l1, l2, pax, pbx) * float(fact2(2 * i - 1)) / pow(2 * gamma, float(i))

func tho66(alpha1, alpha2: SomeFloat, ra, rb: (SomeFloat, SomeFloat, SomeFloat), la, lb: (SomeInteger, SomeInteger, SomeInteger)): SomeFloat =
  let
    gamma = alpha1 + alpha2
    rab2 = dist2(ra[0], ra[1], ra[2], rb[0], rb[1], rb[2])
    pre = exp(-alpha1 * alpha2 * rab2 / gamma) * pow(PI / gamma, 1.5)
    xp = productCenter1d(alpha1, ra[0], alpha2, rb[0])
    yp = productCenter1d(alpha1, ra[1], alpha2, rb[1])
    zp = productCenter1d(alpha1, ra[2], alpha2, rb[2])
    wx = overlap1d(la[0], lb[0], xp - ra[0], xp - rb[0], gamma)
    wy = overlap1d(la[1], lb[1], yp - ra[1], yp - rb[1], gamma)
    wz = overlap1d(la[2], lb[2], zp - ra[2], zp - rb[2], gamma)
  pre * wx * wy * wz

suite "overlapintergrals":
  test "dist2":
    check: approx(dist2(0.5, 0.6, 0.7, 0.8, 0.9, 1.0), 0.27, 1e-5)

  test "productCenter1d":
    check: productCenter1d(2, 3, 4, 5) == 26 / 6

  test "fact2":
    check: fact2(0) == 1
    check: fact2(1) == 1
    check: fact2(2) == 2
    check: fact2(3) == (3 * 1)
    check: fact2(4) == (4 * 2)
    check: fact2(5) == (5 * 3 * 1)
    check: fact2(6) == (6 * 4 * 2)
    check: fact2(7) == (7 * 5 * 3 * 1)

  test "binomialPrefactor":
    check: approx(binomialPrefactor(1, 1, 1, 0.1, 0.2), 0.3, 1e-7)
    check: approx(binomialPrefactor(1, 1, 1, 0.3, 0.4), 0.7, 1e-7)
    check: approx(binomialPrefactor(1, 3, 1, 0.1, 0.2), 0.007, 1e-7)
    check: approx(binomialPrefactor(2, 3, 1, 0.1, 0.2), 0.09, 1e-7)

  test "overlap1d":
    check: approx(overlap1d(1, 1, 0.1, 0.2, 1.0), 0.52, 1e-5)
    check: approx(overlap1d(3, 1, 0.1, 0.2, 1.0), 0.7952, 1e-5)

  let
    za = 1.8
    zb = 2.8
    ra = (0.0, 0.0, 0.0)
    rb = (0.5, 0.8, -0.2)
  test "tho66":
      check: approx(tho66(za, zb, ra, rb, (0, 0, 0), (0, 0, 0)), 0.20373275913014607, 1.0e-16)
      check: approx(tho66(za, zb, ra, rb, (1, 0, 0), (0, 0, 0)), 0.062005622343957505, 1.0e-16)
      check: approx(tho66(za, zb, ra, rb, (1, 1, 0), (1, 1, 0)), -0.00043801221837779696, 1.0e-16)
      check: approx(tho66(za, zb, ra, rb, (2, 1, 0), (1, 1, 0)), -0.0002385994651113168, 1.0e-16)
