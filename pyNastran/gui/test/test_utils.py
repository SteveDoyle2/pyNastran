from __future__ import print_function
import os
import unittest

import pyNastran
from pyNastran.gui.gui_utils.utils import load_csv, load_deflection_csv, load_user_geom

pkg_path = pyNastran.__path__[0]


class GuiUtils(unittest.TestCase):
    def test_gui_csv_01(self):
        csv_filename = os.path.join(pkg_path, '..', 'models', 'solid_bending', 'solid_bending.txt')
        with self.assertRaises(ValueError):
            load_deflection_csv(csv_filename)
        load_csv(csv_filename)

    def test_gui_csv_02(self):
        csv_filename = os.path.join(pkg_path, '..', 'models', 'solid_bending', 'solid_bending_multi_node.txt')
        with self.assertRaises(RuntimeError):
            load_deflection_csv(csv_filename)
        load_csv(csv_filename)

    def test_gui_csv_03(self):
        csv_filename = os.path.join(pkg_path, '..', 'models', 'solid_bending', 'solid_bending_multi_node.csv')
        with self.assertRaises(RuntimeError):
            load_deflection_csv(csv_filename)
        load_csv(csv_filename)

    def test_gui_deflection_csv_01(self):
        csv_filename = os.path.join(pkg_path, '..', 'models', 'solid_bending', 'solid_bending_multi_node.csv')
        with self.assertRaises(RuntimeError):
            load_deflection_csv(csv_filename)

        csv_filename = os.path.join(pkg_path, '..', 'models', 'solid_bending', 'solid_bending_multi_deflection_node.txt')
        load_deflection_csv(csv_filename)

    def test_gui_custom_geom_01(self):
        csv_filename = os.path.join(pkg_path, '..', 'models', 'custom_geom.csv')
        load_user_geom(csv_filename)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
