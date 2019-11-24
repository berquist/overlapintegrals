using OverlapIntegrals, Test

@testset "fact2" begin
    @test fact2(0) == 1
    @test fact2(1) == 1
    @test fact2(2) == 2
    @test fact2(3) == (3 * 1)
    @test fact2(4) == (4 * 2)
    @test fact2(5) == (5 * 3 * 1)
    @test fact2(6) == (6 * 4 * 2)
    @test fact2(7) == (7 * 5 * 3 * 1)
end

@testset "binomial_prefactor" begin
    @test isapprox(binomial_prefactor(1, 1, 1, 0.1, 0.2), 0.3; atol = 1e-7)
    @test isapprox(binomial_prefactor(1, 1, 1, 0.3, 0.4), 0.7; atol = 1e-7)
    @test isapprox(binomial_prefactor(1, 3, 1, 0.1, 0.2), 0.007; atol = 1e-7)
    @test isapprox(binomial_prefactor(2, 3, 1, 0.1, 0.2), 0.09; atol = 1e-7)
end

@testset "overlap1d" begin
    @test isapprox(overlap1d(1, 1, 0.1, 0.2, 1.0), 0.52; atol = 1e-5)
    @test isapprox(overlap1d(3, 1, 0.1, 0.2, 1.0), 0.7952; atol = 1e-5)
end

@testset "tho66" begin
    za = 1.8
    zb = 2.8
    ra = [0.0, 0.0, 0.0]
    rb = [0.5, 0.8, -0.2]

    integral = tho66(za, zb, ra, rb, [0, 0, 0], [0, 0, 0])
    @test isapprox(integral, 0.20373275913014607; atol = 1.0e-16)

    integral = tho66(za, zb, ra, rb, [1, 0, 0], [0, 0, 0])
    @test isapprox(integral, 0.062005622343957505; atol = 1.0e-16)

    integral = tho66(za, zb, ra, rb, [1, 1, 0], [1, 1, 0])
    @test isapprox(integral, -0.00043801221837779696; atol = 1.0e-16)

    integral = tho66(za, zb, ra, rb, [2, 1, 0], [1, 1, 0])
    @test isapprox(integral, -0.0002385994651113168; atol = 1.0e-16)
end
