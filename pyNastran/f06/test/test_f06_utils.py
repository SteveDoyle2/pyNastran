import os
import unittest

import matplotlib
matplotlib.use('Agg')

import pyNastran
from pyNastran.utils.log import get_logger2
from pyNastran.f06.utils import plot_flutter_f06

pkg_path = pyNastran.__path__[0]

class TestF06Utils(unittest.TestCase):
    def test_plot_flutter(self):
        f06_filename = os.path.join(pkg_path, '..', 'models', 'aero', 'bah_plane', 'bah_plane.f06')
        log = get_logger2(log=None, debug=False, encoding='utf-8')
        plot_flutter_f06(f06_filename, show=False, log=log)

    def test_plot_flutter2(self):
        f06_filename = os.path.join(pkg_path, '..', 'models', 'aero', '2_mode_flutter', '0012_flutter.f06')
        log = get_logger2(log=None, debug=False, encoding='utf-8')
        plot_flutter_f06(f06_filename,
                         plot_vg=True, plot_vg_vf=True, plot_root_locus=True,
                         plot_kfreq_damping=True,
                         show=False, log=None)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
