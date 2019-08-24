from scipy.special import binom as binomial


def test_binomial():
    assert binomial(0, 0) == 1.0
    assert binomial(1, 0) == 1.0
    assert binomial(2, 0) == 1.0
    assert binomial(3, 0) == 1.0
    assert binomial(0, 1) == 0.0
    assert binomial(0, 2) == 0.0
    assert binomial(0, 3) == 0.0
    assert binomial(4, 2) == 6.0
    assert binomial(10, 3) == 120.0
    return None
