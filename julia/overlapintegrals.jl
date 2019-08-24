function fact2(n)
    if n <= 0
        return 1
    else
        return n * fact2(n - 2)
    end
end

function binomial_prefactor(s::Integer, ia::Integer, ib::Integer, xpa::Number, xpb::Number)
    acc = 0
    for t = 0 : s + 1
        if (s - ia) <= t && t <= ib
            acc += binomial(ia, s - t) * binomial(ib, t) * (xpa ^ (ia - s + t)) * (xpb ^ (ib - t))
        end
    end
    return acc
end


function overlap1d(l1::Integer, l2::Integer, pax::Number, pbx::Number, gamma::Number)
    acc = 0.0
    for i = 0 : 1 + floor(0.5 * (l1 + l2))
        acc += binomial_prefactor(Integer(2 * i), l1, l2, pax, pbx) * fact2(2 * i - 1) / ((2 * gamma) ^ i)
    end
    return acc
end

function dist2(xa, ya, za, xb, yb, zb)
    return (xa - xb) * (xa - xb) + (ya - yb) * (ya - yb) + (za - zb) * (za - zb)
end

function product_center_1d(za, xa, zb, xb)
    return ((za * xa) + (zb * xb)) / (za + zb)
end

function tho66(alpha1::Number, alpha2::Number, ra, rb, la, lb)
    gamma = alpha1 + alpha2
    xa, ya, za = ra
    xb, yb, zb = rb
    rab2 = dist2(xa, ya, za, xb, yb, zb)
    l1, m1, n1 = la
    l2, m2, n2 = lb
    xp = product_center_1d(alpha1, xa, alpha2, xb)
    yp = product_center_1d(alpha1, ya, alpha2, yb)
    zp = product_center_1d(alpha1, za, alpha2, zb)

    pre = exp(-alpha1 * alpha2 * rab2 / gamma) * (pi / gamma) ^ 1.5

    wx = overlap1d(l1, l2, xp - xa, xp - xb, gamma)
    wy = overlap1d(m1, m2, yp - ya, yp - yb, gamma)
    wz = overlap1d(n1, n2, zp - za, zp - zb, gamma)
    return pre * wx * wy * wz
end
