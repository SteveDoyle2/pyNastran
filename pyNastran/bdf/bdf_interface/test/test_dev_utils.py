import os
import unittest
import pyNastran
from pyNastran.bdf.mesh_utils.dev.create_vectorized_numbered import create_vectorized_numbered

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, '../', 'models')

class DevUtils(unittest.TestCase):
    """tests various dev functions"""

    def test_convert_bdf(self):
        """tests create_vectorized_numbered"""
        bdf_filename_in = os.path.join(MODEL_PATH, 'elements', 'static_elements.bdf')
        bdf_filename_out = os.path.join(MODEL_PATH, 'elements', 'static_elements_convert.bdf')
        create_vectorized_numbered(bdf_filename_in, bdf_filename_out, debug=False)
        os.remove(bdf_filename_out)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

