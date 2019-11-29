#[macro_use]
extern crate approx;
use num_integer;
use std::f64::consts;

fn dist2(xa: f64, ya: f64, za: f64, xb: f64, yb: f64, zb: f64) -> f64 {
    (xa - xb).powf(2.0) + (ya - yb).powf(2.0) + (za - zb).powf(2.0)
}

fn product_center_1d(za: f64, xa: f64, zb: f64, xb: f64) -> f64 {
    ((za * xa) + (zb * xb)) / (za + zb)
}

fn fact2(n: i64) -> i64 {
    let mut result = 1;
    let mut counter = n;
    while counter > 0 {
        result *= counter;
        if counter == 1 {
            break;
        }
        counter -= 2;
    }
    result
}

fn binomial_prefactor(s: i64, ia: i64, ib: i64, xpa: f64, xpb: f64) -> f64 {
    let mut result = 0.0;
    for t in 0..=s {
        // println!("s: {} ia: {} t: {} ib: {}", s, ia, t, ib);
        if ((s - ia) <= t) && (t <= ib) {
            result += (num_integer::binomial(ia, s - t) * num_integer::binomial(ib, t)) as f64
                * xpa.powi((ia - s + t) as i32)
                * xpb.powi((ib - t) as i32);
        }
    }
    result
}

fn overlap1d(l1: i64, l2: i64, pax: f64, pbx: f64, gamma: f64) -> f64 {
    let mut result = 0.0;
    for i in 0..=(0.5 * (l1 + l2) as f64).floor() as i64 {
        result += binomial_prefactor(2 * i, l1, l2, pax, pbx) * fact2(2 * i - 1) as f64
            / (2.0 * gamma).powi(i as i32);
    }
    result
}

fn tho66(
    alpha1: f64,
    alpha2: f64,
    ra: (f64, f64, f64),
    rb: (f64, f64, f64),
    la: (i64, i64, i64),
    lb: (i64, i64, i64),
) -> f64 {
    let gamma = alpha1 + alpha2;
    let rab2 = dist2(ra.0, ra.1, ra.2, rb.0, rb.1, rb.2);
    let pre = (-alpha1 * alpha2 * rab2 / gamma).exp() * (consts::PI / gamma).powf(1.5);
    let xp = product_center_1d(alpha1, ra.0, alpha2, rb.0);
    let yp = product_center_1d(alpha1, ra.1, alpha2, rb.1);
    let zp = product_center_1d(alpha1, ra.2, alpha2, rb.2);
    let wx = overlap1d(la.0, lb.0, xp - ra.0, xp - rb.0, gamma);
    let wy = overlap1d(la.1, lb.1, yp - ra.1, yp - rb.1, gamma);
    let wz = overlap1d(la.2, lb.2, zp - ra.2, zp - rb.2, gamma);
    pre * wx * wy * wz
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dist2() {
        assert_abs_diff_eq!(dist2(0.5, 0.6, 0.7, 0.8, 0.9, 1.0), 0.27);
    }

    #[test]
    fn test_product_center_1d() {
        assert_eq!(product_center_1d(2.0, 3.0, 4.0, 5.0), 26.0 / 6.0);
    }

    #[test]
    fn test_fact2() {
        assert_eq!(fact2(0), 1);
        assert_eq!(fact2(1), 1);
        assert_eq!(fact2(2), 2);
        assert_eq!(fact2(3), (3 * 1));
        assert_eq!(fact2(4), (4 * 2));
        assert_eq!(fact2(5), (5 * 3 * 1));
        assert_eq!(fact2(6), (6 * 4 * 2));
        assert_eq!(fact2(7), (7 * 5 * 3 * 1));
    }

    #[test]
    fn test_binomial_prefactor() {
        assert_abs_diff_eq!(binomial_prefactor(1, 1, 1, 0.1, 0.2), 0.3);
        assert_abs_diff_eq!(binomial_prefactor(1, 1, 1, 0.3, 0.4), 0.7);
        assert_abs_diff_eq!(binomial_prefactor(1, 3, 1, 0.1, 0.2), 0.007);
        assert_abs_diff_eq!(binomial_prefactor(2, 3, 1, 0.1, 0.2), 0.09);
    }

    #[test]
    fn test_overlap1d() {
        assert_abs_diff_eq!(overlap1d(1, 1, 0.1, 0.2, 1.0), 0.52, epsilon = 1.0e-40);
        assert_abs_diff_eq!(overlap1d(3, 1, 0.1, 0.2, 1.0), 0.7952, epsilon = 1.0e-40);
    }

    #[test]
    fn test_tho66() {
        let za = 1.8;
        let zb = 2.8;
        let ra = (0.0, 0.0, 0.0);
        let rb = (0.5, 0.8, -0.2);
        assert_abs_diff_eq!(
            tho66(za, zb, ra, rb, (0, 0, 0), (0, 0, 0)),
            0.20373275913014607
        );
        assert_abs_diff_eq!(
            tho66(za, zb, ra, rb, (1, 0, 0), (0, 0, 0)),
            0.062005622343957505
        );
        assert_abs_diff_eq!(
            tho66(za, zb, ra, rb, (1, 1, 0), (1, 1, 0)),
            -0.00043801221837779696
        );
        assert_abs_diff_eq!(
            tho66(za, zb, ra, rb, (2, 1, 0), (1, 1, 0)),
            -0.0002385994651113168
        );
    }
}
