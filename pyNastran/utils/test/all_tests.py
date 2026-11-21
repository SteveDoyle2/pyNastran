"""tests for ``pyNastran.utils``"""
try:
    import scipy
    IS_SCIPY = True
except ImportError:
    IS_SCIPY = False

from pyNastran.utils.test.test_utils import *
from pyNastran.utils.test.test_atmosphere import *
if IS_SCIPY:
    from pyNastran.utils.test.test_dict_to_h5py import *
    from pyNastran.utils.test.test_concave_hull import *


if __name__ == '__main__':  # pragma: no cover
    import unittest
    unittest.main()
