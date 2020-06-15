import unittest
from numpy import allclose
from pyNastran.gui.gui_objects.utils import autotype_value


class TestSettings(unittest.TestCase):
    def test_settings_bool(self):
        value = autotype_value('true', bool)
        assert value is True

        value = autotype_value('false', bool)
        assert value is False

        value = autotype_value(True, bool)
        assert value is True

        value = autotype_value(False, bool)
        assert value is False

    def test_settings_int(self):
        value = autotype_value('1', int)
        assert allclose(value, 1)

        value = autotype_value(1, int)
        assert allclose(value, 1)

        with self.assertRaises(ValueError):
            value = autotype_value('4.2', int)

        value = autotype_value(['1', '2', '3'], int)
        assert allclose(value, [1, 2, 3])

        value = autotype_value([1, 2, 3], int)
        assert allclose(value, [1, 2, 3])

        value = autotype_value(('1', '2', '3'), int)
        assert allclose(value, (1, 2, 3))

        value = autotype_value((1, 2, 3), int)
        assert allclose(value, (1, 2, 3))

    def test_settings_float(self):
        value = autotype_value('1.1', float)
        assert allclose(value, 1.1)

        value = autotype_value(1.1, float)
        assert allclose(value, 1.1)

        value = autotype_value(['1.1', '2.2', '3.3'], float)
        assert allclose(value, [1.1, 2.2, 3.3])

        value = autotype_value([1.1, 2.2, 3.3], float)
        assert allclose(value, [1.1, 2.2, 3.3])

        value = autotype_value(('1.1', '2.2', '3.3'), float)
        assert allclose(value, (1.1, 2.2, 3.3))

        value = autotype_value((1.1, 2.2, 3.3), float)
        assert allclose(value, (1.1, 2.2, 3.3))

if __name__ == '__main__':   # pragma: no cover
    unittest.main()
