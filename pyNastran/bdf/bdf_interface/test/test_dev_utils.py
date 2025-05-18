import os
import unittest
import pyNastran

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, '../', 'models')

class DevUtils(unittest.TestCase):
    """tests various dev functions"""

    def test_convert_bdf(self):
        """test removed"""
        pass


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
