import os
import unittest
import matplotlib.pyplot as plt

#try:  # pragma: no cover
    #plt.figure()
    #plt.close()
#except:  # pragma: no cover
plt.switch_backend('Agg')

import pyNastran
from pyNastran.utils.log import get_logger2
from pyNastran.f06.utils import split_float_colons, split_int_colon
from pyNastran.f06.parse_flutter import plot_flutter_f06

PKG_PATH = pyNastran.__path__[0]

class TestF06Utils(unittest.TestCase):
    def test_split_float_colon(self):
        """tests split_float_colon"""
        a = split_float_colons('1:')
        b = split_float_colons('1:5')
        c = split_float_colons(':4')

        assert a == [1.0, None], a
        assert b == [1.0, 5.0], b
        assert c == [None, 4.0], c
        with self.assertRaises(AssertionError):
            split_float_colons('1:5:2')

    def test_split_int_colon(self):
        """tests split_int_colon"""
        a = split_int_colon('1:5')
        assert a == [1, 2, 3, 4, 5], a

        b = split_int_colon('1:5:2')
        assert b == [1, 3, 5], b

        c = split_int_colon(':4')
        assert c == [0, 1, 2, 3, 4], c

        d = split_int_colon('1:5,10:15')
        assert d == [1, 2, 3, 4, 5, 10, 11, 12, 13, 14, 15], d

        d = split_int_colon('10:15,1:5')
        assert d == [1, 2, 3, 4, 5, 10, 11, 12, 13, 14, 15], d

    def test_plot_flutter(self):
        """tests plot_flutter_f06"""
        f06_filename = os.path.join(PKG_PATH, '..', 'models', 'aero', 'bah_plane', 'bah_plane.f06')
        log = get_logger2(log=None, debug=None, encoding='utf-8')
        plot_flutter_f06(f06_filename, show=False, log=log)

    def test_plot_flutter2(self):
        """tests plot_flutter_f06"""
        f06_filename = os.path.join(PKG_PATH, '..', 'models', 'aero', '2_mode_flutter', '0012_flutter.f06')
        log = get_logger2(log=None, debug=None, encoding='utf-8')
        plot_flutter_f06(f06_filename,
                         plot_vg=True, plot_vg_vf=True, plot_root_locus=True,
                         plot_kfreq_damping=True,
                         export_zona=True,
                         export_veas=True,
                         export_f06=True,
                         show=False, log=log)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
