from bioeqpy.power import achieved_power, power_curve, sample_size


def test_achieved_power_is_probability():
    power = achieved_power(n=24, cv=0.25, gmr=0.95)

    assert 0.0 <= power <= 1.0


def test_sample_size_returns_even_n():
    n = sample_size(cv=0.25, gmr=0.95, power=0.80)

    assert n >= 4
    assert n % 2 == 0
    assert achieved_power(n=n, cv=0.25, gmr=0.95) >= 0.80


def test_power_curve_shape():
    curve = power_curve(range(12, 18, 2), cv=0.25)

    assert list(curve.columns) == ["N", "Power"]
    assert len(curve) == 3

