"""tests log.py"""
import os
import unittest

from pyNastran.utils.log import make_log

class TestLog(unittest.TestCase):

    def test_make_log(self):
        """tests make_log"""
        make_log()
        os.remove('pyNastran.log')

if __name__ == "__main__":  # pragma: no cover
    unittest.main()
