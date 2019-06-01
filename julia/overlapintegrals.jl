function fact2(n)
    if n <= 0
        return 1
    else
        return n * fact2(n - 2)
    end
end

function test_fact2()
    @assert fact2(0) == 1
    @assert fact2(1) == 1
    @assert fact2(2) == 2
    @assert fact2(3) == (3 * 1)
    @assert fact2(4) == (4 * 2)
    @assert fact2(5) == (5 * 3 * 1)
    @assert fact2(6) == (6 * 4 * 2)
    @assert fact2(7) == (7 * 5 * 3 * 1)
end

test_fact2()

function binomial_prefactor(s::Integer, ia::Integer, ib::Integer, xpa::Number, xpb::Number)
    acc = 0
    for t = 0 : s + 1
        if (s - ia) <= t && t <= ib
            acc += binomial(ia, s - t) * binomial(ib, t) * (xpa ^ (ia - s + t)) * (xpb ^ (ib - t))
        end
    end
    return acc
end


function test_binomial_prefactor()
    @assert isapprox(binomial_prefactor(1, 1, 1, 0.1, 0.2), 0.3; atol = 1e-7)
    @assert isapprox(binomial_prefactor(1, 1, 1, 0.3, 0.4), 0.7; atol = 1e-7)
    @assert isapprox(binomial_prefactor(1, 3, 1, 0.1, 0.2), 0.007; atol = 1e-7)
    @assert isapprox(binomial_prefactor(2, 3, 1, 0.1, 0.2), 0.09; atol = 1e-7)
end

test_binomial_prefactor()

function overlap1d(l1::Integer, l2::Integer, pax::Number, pbx::Number, gamma::Number)
    acc = 0.0
    for i = 0 : 1 + floor(0.5 * (l1 + l2))
        acc += binomial_prefactor(Integer(2 * i), l1, l2, pax, pbx) * fact2(2 * i - 1) / ((2 * gamma) ^ i)
    end
    return acc
end

function test_overlap1d()
    @assert isapprox(overlap1d(1, 1, 0.1, 0.2, 1.0), 0.52; atol = 1e-5)
    @assert isapprox(overlap1d(3, 1, 0.1, 0.2, 1.0), 0.7952; atol = 1e-5)
end

test_overlap1d()

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

function test_tho66()
    za = 1.8
    zb = 2.8
    ra = [0.0, 0.0, 0.0]
    rb = [0.5, 0.8, -0.2]

    integral = tho66(za, zb, ra, rb, [0, 0, 0], [0, 0, 0])
    @assert isapprox(integral, 0.20373275913014607; atol = 1.0e-16)

    integral = tho66(za, zb, ra, rb, [1, 0, 0], [0, 0, 0])
    @assert isapprox(integral, 0.062005622343957505; atol = 1.0e-16)

    integral = tho66(za, zb, ra, rb, [1, 1, 0], [1, 1, 0])
    @assert isapprox(integral, -0.00043801221837779696; atol = 1.0e-16)

    integral = tho66(za, zb, ra, rb, [2, 1, 0], [1, 1, 0])
    @assert isapprox(integral, -0.0002385994651113168; atol = 1.0e-16)
end

test_tho66()
