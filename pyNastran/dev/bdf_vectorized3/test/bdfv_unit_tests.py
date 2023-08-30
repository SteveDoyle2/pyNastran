import pathlib
import unittest
import pyNastran
from pyNastran.dev.bdf_vectorized3.bdf import read_bdf

PKG_PATH = pathlib.Path(pyNastran.__path__[0])
TEST_PATH = PKG_PATH / 'bdf' / 'test'
MODEL_PATH = PKG_PATH / '..' / 'models'

class TestBdfVectorized3(unittest.TestCase):
    def test_bwb(self):
        bdf_filename = MODEL_PATH / 'bwb' / 'bwb_saero.bdf'
        bdf_filename_out = MODEL_PATH / 'bwb' / 'bwb_saero_vectorized.bdf'
        model = read_bdf(
            bdf_filename, validate=True, xref=True, punch=False,
            save_file_structure=False, skip_cards=None, read_cards=None,
            encoding=None, log=None, debug=True, mode='msc')
        model.write_bdf(bdf_filename_out)

if __name__ == '__main__':
    unittest.main()
