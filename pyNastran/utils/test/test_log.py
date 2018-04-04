"""tests log.py"""

import os
import unittest

from pyNastran.utils.log import make_log
import pyNastran
from pyNastran.bdf.bdf import BDF

PKG_PATH = pyNastran.__path__[0]
BDF_FILENAME = os.path.join(PKG_PATH, '..', 'models', 'solid_bending', 'solid_bending.bdf')

class TestLog(unittest.TestCase):

    @staticmethod
    def test_log_01():
        model = BDF(debug=True)
        model.read_bdf(BDF_FILENAME)

        model2 = BDF(debug=False)
        model2.read_bdf(BDF_FILENAME)

        model3 = BDF(debug=None)
        model3.read_bdf(BDF_FILENAME)

    def test_make_log(self):
        """tests make_log"""
        make_log()

if __name__ == "__main__":  # pragma: no cover
    unittest.main()
