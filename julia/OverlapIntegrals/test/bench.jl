using BenchmarkTools
using OverlapIntegrals

za = 1.8
zb = 2.8
ra = [0.0, 0.0, 0.0]
rb = [0.5, 0.8, -0.2]

@btime tho66(za, zb, ra, rb, [0, 0, 0], [0, 0, 0])
@btime tho66(za, zb, ra, rb, [2, 1, 0], [1, 1, 0])
@btime tho66(za, zb, ra, rb, [3, 3, 3], [1, 1, 0])
