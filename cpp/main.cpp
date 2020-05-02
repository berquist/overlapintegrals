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

double binomial(unsigned n, unsigned k) {
    // Mimic the behavior of scipy.special.binomial
    if (k > n) {
        return 0.0;
    } else {
        return boost::math::binomial_coefficient<double>(n, k);
    }
}

double fact2(long long i) {
    // Yes, I know this is incorrect
    if (i < 0) {
        return 1.0;
    } else {
        return boost::math::double_factorial<double>(i);
    }
}

double binomial_prefactor(long long s, long long ia, long long ib, double xpa, double xpb) {
    double sum = 0.0;
    for (long long t = 0; t < s + 1; t++) {
        if (((s - ia) <= t) && (t <= ib)) {
            // std::cout << "(ia, s - t): (" << ia << ", " << s - t << ")" << std::endl;
            // std::cout << "(ib, t): (" << ib << ", " << t << ")" << std::endl;
            sum += binomial(ia, s - t) * binomial(ib, t) * pow(xpa, ia - s + t) * pow(xpb, ib - t);
        }
    }
    return sum;
}

double overlap1d(long long l1, long long l2, double pax, double pbx, double gamma) {
    double sum = 0.0;
    const long long lim = 1 + floor(0.5 * (l1 + l2));
    // std::cout << "lim: " << lim << std::endl;
    for (long long i = 0; i < lim; i++) {
        const long long arg = 2 * i - 1;
        // std::cout << "2 * i - 1: " << arg << std::endl;
        // if (arg >= 0) {
            const double p1 = binomial_prefactor(2 * i, l1, l2, pax, pbx);
            const double p2 = fact2(2 * i - 1);
            const double p3 = pow(2 * gamma, i);
            // std::cout << p1 << " " << p2 << " " << p3 << std::endl;
            const double term = p1 * p2 / p3;
            sum += term;
        // }
    }
    return sum;
}

double tho66(double alpha1, double alpha2, double xa, double ya, double za,
             double xb, double yb, double zb, long long l1, long long m1, long long n1,
             long long l2, long long m2, long long n2) {
    const double gamma = alpha1 + alpha2;
    const double xp = product_center_1d(alpha1, xa, alpha2, xb);
    const double yp = product_center_1d(alpha1, ya, alpha2, yb);
    const double zp = product_center_1d(alpha1, za, alpha2, zb);
    const double pre = exp(-alpha1 * alpha2 * dist2(xa, ya, za, xb, yb, zb) / gamma) * pow(M_PI / gamma, 1.5);
    const double wx = overlap1d(l1, l2, xp - xa, xp - xb, gamma);
    const double wy = overlap1d(m1, m2, yp - ya, yp - yb, gamma);
    const double wz = overlap1d(n1, n2, zp - za, zp - zb, gamma);
    return pre * wx * wy * wz;
}

// No need for this with the presence of tests.
//
// int main() {
//     std::cout << binomial_prefactor(1, 1, 1, 0.1, 0.2) << std::endl;
//     std::cout << binomial_prefactor(1, 1, 1, 0.3, 0.4) << std::endl;
//     std::cout << binomial_prefactor(1, 3, 1, 0.1, 0.2) << std::endl;
//     std::cout << binomial_prefactor(2, 3, 1, 0.1, 0.2) << std::endl;
//     for (unsigned long long i = 0; i <= 21; i++) {
//         // the template parameter needs to be floating-point
//         std::cout << boost::math::double_factorial<double>(i) << std::endl;
//     }
//     std::cout << overlap1d(1, 1, 0.1, 0.2, 1.0) << std::endl;
//     std::cout << overlap1d(3, 1, 0.1, 0.2, 1.0) << std::endl;
//     return 0;
// }

TEST_CASE("dist2", "dist2") {
    REQUIRE(dist2(0.5, 0.6, 0.7, 0.8, 0.9, 1.0) == Approx(0.27));
}

TEST_CASE("product_center_1d", "product_center_1d") {
    REQUIRE(product_center_1d(1.8, 0.0, 2.8, 0.5) == Approx(0.304348));
}

TEST_CASE("binomial", "binomial") {
    REQUIRE(binomial(0, 0) == 1.0);
    REQUIRE(binomial(1, 0) == 1.0);
    REQUIRE(binomial(2, 0) == 1.0);
    REQUIRE(binomial(3, 0) == 1.0);
    REQUIRE(binomial(0, 1) == 0.0);
    REQUIRE(binomial(0, 2) == 0.0);
    REQUIRE(binomial(0, 3) == 0.0);
    REQUIRE(binomial(4, 2) == 6.0);
    REQUIRE(binomial(10, 3) == 120.0);
}

TEST_CASE("fact2", "fact2") {
    REQUIRE(fact2(-1) == 1);
    REQUIRE(fact2(0) == 1);
    REQUIRE(fact2(1) == 1);
    REQUIRE(fact2(2) == 2);
    REQUIRE(fact2(3) == (3 * 1));
    REQUIRE(fact2(4) == (4 * 2));
    REQUIRE(fact2(5) == (5 * 3 * 1));
    REQUIRE(fact2(6) == (6 * 4 * 2));
    REQUIRE(fact2(7) == (7 * 5 * 3 * 1));
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
    const double alpha1 = 1.8;
    const double alpha2 = 2.8;
    const double xa = 0.0;
    const double ya = 0.0;
    const double za = 0.0;
    const double xb = 0.5;
    const double yb = 0.8;
    const double zb = -0.2;

    double integral;

    integral = tho66(alpha1, alpha2, xa, ya, za, xb, yb, zb, 0, 0, 0, 0, 0, 0);
    REQUIRE(integral == Approx(0.20373275913014607));

    integral = tho66(alpha1, alpha2, xa, ya, za, xb, yb, zb, 1, 0, 0, 0, 0, 0);
    REQUIRE(integral == Approx(0.062005622343957505));

    integral = tho66(alpha1, alpha2, xa, ya, za, xb, yb, zb, 1, 1, 0, 1, 1, 0);
    REQUIRE(integral == Approx(-0.00043801221837779696));

    integral = tho66(alpha1, alpha2, xa, ya, za, xb, yb, zb, 2, 1, 0, 1, 1, 0);
    REQUIRE(integral == Approx(-0.0002385994651113168));
}
