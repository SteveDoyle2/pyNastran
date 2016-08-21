from __future__ import print_function
import os
import unittest

import pyNastran
from pyNastran.gui.utils import load_csv, load_user_geom

pkg_path = pyNastran.__path__[0]


class GuiUtils(unittest.TestCase):
    def test_gui_utils_01(self):
        csv_filename = os.path.join(pkg_path, '..', 'models', 'solid_bending', 'solid_bending.txt')
        load_csv(csv_filename)

    def test_gui_utils_02(self):
        csv_filename = os.path.join(pkg_path, '..', 'models', 'solid_bending', 'solid_bending_multi.txt')
        load_csv(csv_filename)

    def test_gui_utils_03(self):
        csv_filename = os.path.join(pkg_path, '..', 'models', 'solid_bending', 'solid_bending_multi.csv')
        load_csv(csv_filename)

    def test_gui_utils_03(self):
        csv_filename = os.path.join(pkg_path, '..', 'models', 'custom_geom.csv')
        load_user_geom(csv_filename)

if __name__ == '__main__':
    unittest.main()
