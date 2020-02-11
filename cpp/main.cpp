#include <boost/math/special_functions/math_fwd.hpp>
#include <cmath>
#define _USE_MATH_DEFINES
#include <boost/math/special_functions/binomial.hpp>

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>


double dist2(double xa, double ya, double za, double xb, double yb, double zb) {
    const double dx = xa - xb;
    const double dy = ya - yb;
    const double dz = za - zb;
    return (dx * dx) + (dy * dy) + (dz * dz);
}

double product_center_1d(double za, double xa, double zb, double xb) {
    return ((za * xa) + (zb * xb)) / (za + zb);
}

double binomial_prefactor(int s, int ia, int ib, double xpa, double xpb) {
  double sum = 0.0;
  for (size_t t = 0; t < s + 1; t++) {
    if (((s - ia) <= t) && (t <= ib)) {
        sum += boost::math::binomial_coefficient<double>(ia, s - t) * boost::math::binomial_coefficient<double>(ib, t) * pow(xpa, ia - s + t)* pow(xpb, ib - t);
      }
  }
  return sum;
}

double overlap1d(int l1, int l2, double pax, double pbx, double gamma) {
  double sum = 0.0;
  for (size_t i = 0; i < (1 + floor(0.5 * (l1 + l2))); i++) {
      sum += binomial_prefactor(2 * i, l1, l2, pax, pbx) * boost::math::double_factorial<double>(2 * i - 1) / pow(2 * gamma, i);
  }
    return sum;  
}

double tho66(double alpha1, double alpha2, double xa, double ya, double za,
             double xb, double yb, double zb, double l1, double m1, double n1,
             double l2, double m2, double n2) {
    const double gamma = alpha1 + alpha2;
    const double xp = product_center_1d(alpha1, xa, alpha2, xb);
    const double yp = product_center_1d(alpha1, ya, alpha2, yb);
    const double zp = product_center_1d(alpha1, za, alpha2, zb);
    const double pre = exp(-alpha1 * alpha2 * dist2(xa, ya, za, xb, yb, zb)) * pow(M_PI / gamma, 1.5);
    const double wx = overlap1d(l1, l2, xp - xa, xp - xb, gamma);
    const double wy = overlap1d(m1, m2, yp - ya, yp - yb, gamma);
    const double wz = overlap1d(n1, n2, zp - za, zp - zb, gamma);
    return pre * wx * wy * wz;
}

TEST_CASE("dist2", "dist2") {
    REQUIRE(dist2(0.5, 0.6, 0.7, 0.8, 0.9, 1.0) == Approx(0.27));
}

TEST_CASE("binomial_prefactor", "binomial_prefactor") {
    REQUIRE(binomial_prefactor(1, 1, 1, 0.1, 0.2) == Approx(0.3));
    REQUIRE(binomial_prefactor(1, 1, 1, 0.3, 0.4) == Approx(0.7));
    REQUIRE(binomial_prefactor(1, 3, 1, 0.1, 0.2) == Approx(0.007));
    REQUIRE(binomial_prefactor(2, 3, 1, 0.1, 0.2) == Approx(0.09));
}

TEST_CASE("overlap1d", "overlap1d") {
    REQUIRE(overlap1d(1, 1, 0.1, 0.2, 1.0) == Approx(0.52));
    REQUIRE(overlap1d(3, 1, 0.1, 0.2, 1.0) == Approx(0.7952));
}

TEST_CASE("tho66", "tho66") {
  const auto alpha1 = 1.8;
  const auto alpha2 = 2.8;
  const auto xa = 0.0;
  const auto ya = 0.0;
  const auto za = 0.0;
  const auto xb = 0.5;
  const auto yb = 0.8;
  const auto zb = -0.2;

  double integral;

  integral = tho66(za, zb, xa, ya, za, xb, yb, zb, 0, 0, 0, 0, 0, 0);
  REQUIRE(integral == Approx(0.20373275913014607));

  integral = tho66(za, zb, xa, ya, za, xb, yb, zb, 1, 0, 0, 0, 0, 0);
  REQUIRE(integral == Approx(0.062005622343957505));

  integral = tho66(za, zb, xa, ya, za, xb, yb, zb, 1, 1, 0, 1, 1, 0);
  REQUIRE(integral == Approx(-0.00043801221837779696));

  integral = tho66(za, zb, xa, ya, za, xb, yb, zb, 2, 1, 0, 1, 1, 0);
  REQUIRE(integral == Approx(-0.0002385994651113168));
}
