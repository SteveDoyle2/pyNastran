from __future__ import print_function

from pyNastran.bdf.mesh_utils.test.all_tests_no_matplotlib import *
from pyNastran.bdf.mesh_utils.test.test_cutting_plane import TestCuttingPlane


if __name__ == "__main__":  # pragma: no cover
    import os
    import unittest
    on_rtd = os.environ.get('READTHEDOCS', None)
    if on_rtd is None:
        unittest.main()

