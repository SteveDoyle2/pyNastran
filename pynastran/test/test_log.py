# -*- coding: utf-8 -*-

import os
import unittest

import pyNastran
from pyNastran.bdf.bdf import BDF
pkg_path = pyNastran.__path__[0]

bdf_filename = os.path.join(pkg_path, '..', 'models', 'solid_bending', 'solid_bending.bdf')

class TestLog(unittest.TestCase):

    def test_log_01(self):
        print('---------------------------------------------------------------')
        model = BDF(debug=True)
        model.read_bdf(bdf_filename)

        print('---------------------------------------------------------------')
        model2 = BDF(debug=False)
        model2.read_bdf(bdf_filename)

        print('---------------------------------------------------------------')
        model3 = BDF(debug=None)
        model3.read_bdf(bdf_filename)
        print('---------------------------------------------------------------')

if __name__ == "__main__":
    import unittest
    unittest.main()

