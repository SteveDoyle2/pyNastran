import unittest

import os
import pyNastran
from pyNastran.bdf.bdf import BDF

root_path = pyNastran.__path__[0]
test_path = os.path.join(root_path, 'bdf', 'test', 'unit')

class TestReadWrite(unittest.TestCase):

    def test_write_1(self):
        """
        Tests 1 read method and various write methods
        """
        model = BDF()

        bdf_name = os.path.join(test_path, 'test_mass.dat')
        model.readBDF(bdf_name)
        model.write_bdf('test_mass1a.out', size=8, debug=False)
        model.write_bdf('test_mass2a.out', size=8, debug=False)

        model.writeBDF('test_mass1b.out', size=8, debug=False)
        model.writeBDFAsPatran('test_mass2b.out', size=8, debug=False)
        os.remove(os.path.join(test_path, 'test_mass1a.out'))
        os.remove(os.path.join(test_path, 'test_mass2a.out'))
        os.remove(os.path.join(test_path, 'test_mass1b.out'))
        os.remove(os.path.join(test_path, 'test_mass2b.out'))

    def test_punch_1(self):
        """
        Tests punch file reading
        """
        model = BDF()
        bdf_name = os.path.join(test_path, 'include_dir', 'include.inc')
        model.readBDF(bdf_name, xref=False, punch=True)

        model2 = BDF()
        #bdf_name = os.path.join(test_path, 'include_dir', 'include.inc')
        model2.read_bdf(bdf_name, xref=False, punch=True)

    def test_read_include_dir_1(self):
        """
        Tests various read methods using various include files
        """
        # fails
        model = BDF()
        bdf_name = os.path.join(test_path, 'test_include.bdf')
        self.assertRaises(IOError, model.read_bdf, bdf_name, include_dir=None, xref=True, punch=False)

        # passes
        full_path = os.path.join(test_path, 'include_dir')
        model2 = BDF()
        model2.read_bdf('test_include.bdf', include_dir=full_path, xref=True, punch=False)

        # repeat with slightly different calling signature - passes
        model2 = BDF()
        model2.readBDF('test_include.bdf', includeDir=full_path, xref=True, punch=False)


if __name__ == '__main__':
    # passes if you're in the local folder, fails if you aren't
    #model2 = BDF()
    #model2.read_bdf('test_include.bdf', include_dir='include_dir', xref=True, punch=False)

    unittest.main()