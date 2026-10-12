from pyNastran.all_tests_gui import *
from pyNastran.all_tests_no_gui import *


if __name__ == "__main__":  # pragma: no cover
    import unittest
    with np.errstate(divide='raise', over='raise', under='raise', invalid='raise'):
        unittest.main()
