import overlapintegrals
import criterion

const
  za = 1.8
  zb = 2.8
  ra = (0.0, 0.0, 0.0)
  rb = (0.5, 0.8, -0.2)

var cfg = newDefaultConfig()

benchmark cfg:
  proc measure_tho66_0_0_0_0_0_0() {.measure.} =
    blackBox tho66(za, zb, ra, rb, (0, 0, 0), (0, 0, 0))
  proc measure_tho66_2_1_0_1_1_0() {.measure.} =
    blackBox tho66(za, zb, ra, rb, (2, 1, 0), (1, 1, 0))
  proc measure_tho66_3_3_3_1_1_0() {.measure.} =
    blackBox tho66(za, zb, ra, rb, (3, 3, 3), (1, 1, 0))
