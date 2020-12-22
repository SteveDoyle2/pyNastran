"""tests log.py"""
import os
import unittest

from cpylog.test_log import TestLog
from pyNastran.utils.log import make_log

class TestMakeLog(unittest.TestCase):

    def test_make_log(self):
        """tests make_log"""
        make_log()
        os.remove('pyNastran.log')

if __name__ == "__main__":  # pragma: no cover
    unittest.main()
