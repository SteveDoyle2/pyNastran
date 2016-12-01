import os
import unittest
import matplotlib

import pyNastran
from pyNastran.f06.utils import plot_flutter_f06

matplotlib.use('Agg')
pkg_path = pyNastran.__path__[0]

class TestF06Utils(unittest.TestCase):
    def test_plot_flutter(self):
        f06_filename = os.path.join(pkg_path, '..', 'models', 'aero', 'bah_plane', 'bah_plane.f06')
        plot_flutter_f06(f06_filename, show=False)

if __name__ == '__main__':
    unittest.main()
