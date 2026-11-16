import pytest
from bmi import calculate_bmi


def test_standard_bmi():
    assert calculate_bmi(180.0, 81.0) == 25.0


def test_bmi_rounding():
    assert calculate_bmi(175.0, 70.0) == 22.86


def test_bmi_float_inputs():
    assert calculate_bmi(162.5, 58.5) == 22.15


def test_zero_height_raises_zero_division_error():
    with pytest.raises(ZeroDivisionError):
        calculate_bmi(0.0, 70.0)
