import unittest

import os
import pyNastran
from pyNastran.bdf.bdf import BDF

root_path = pyNastran.__path__[0]
test_path = os.path.join(root_path, 'bdf', 'test')

class TestSolids(unittest.TestCase):

    def test_dmig_1(self):
        """
        Tests DMIG reading
        """
        from numpy import array_equal
        model = BDF()
        bdf_name = os.path.join(test_path, 'cards', 'dmig.bdf')
        model.read_bdf(bdf_name, xref=False, punch=True)
        out = model.dmigs['MAAX'].getMatrix(is_sparse=False)
        MAAX_actual, rowsReversed, colsReversed = out
        print "---MAAX_actual---\n", MAAX_actual
        #print "---out---\n", out

        MAAX_expected = [
            [1.0,  0.5,  0.25],
            [0.5,  2.0,  0.75],
            [0.25, 0.75, 3.0 ],
        ]
        self.assertTrue(array_equal(MAAX_expected, MAAX_actual))

        #--------------------------------
        out = model.dmigs['BAAX'].getMatrix(is_sparse=False)
        BAAX_actual, rowsReversed, colsReversed = out
        print "---BAAX_actual---\n", BAAX_actual
        BAAX_expected = [
            [1.0, 0.5, 0.25],
            [0.0, 2.0, 0.75],
            [0.0, 0.0, 3.0 ],
        ]
        self.assertTrue(array_equal(BAAX_expected, BAAX_actual))

        #model2 = BDF()
        #bdf_name = os.path.join(test_path, 'include_dir', 'include.inc')
        #model2.read_bdf(bdf_name, xref=False, punch=True)


if __name__ == '__main__':
    unittest.main()
