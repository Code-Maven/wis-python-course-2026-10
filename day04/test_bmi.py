import unittest
from bmi import calculate_bmi


class TestCalculateBmi(unittest.TestCase):
    def test_standard_bmi(self):
        # 180 cm, 81 kg -> 81 / (1.8^2) = 81 / 3.24 = 25.0
        self.assertEqual(calculate_bmi(180.0, 81.0), 25.0)

    def test_bmi_rounding(self):
        # 175 cm, 70 kg -> 70 / (1.75^2) = 70 / 3.0625 = 22.8571... -> 22.86
        self.assertEqual(calculate_bmi(175.0, 70.0), 22.86)

    def test_bmi_float_inputs(self):
        # 162.5 cm, 58.5 kg -> 58.5 / (1.625^2) = 58.5 / 2.640625 = 22.1538... -> 22.15
        self.assertEqual(calculate_bmi(162.5, 58.5), 22.15)

    def test_zero_height_raises_zero_division_error(self):
        with self.assertRaises(ZeroDivisionError):
            calculate_bmi(0.0, 70.0)


if __name__ == "__main__":
    unittest.main()
