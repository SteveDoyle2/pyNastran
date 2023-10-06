import unittest
try:
    import vtk
    USE_VTK = True
except ImportError:
    USE_VTK = False

from pyNastran.dev.bdf_vectorized3.test.test_models import *
if USE_VTK:
    from pyNastran.dev.bdf_vectorized3.test.test_gui_models import *

from pyNastran.dev.bdf_vectorized3.test.bdfv_unit_tests import *
from pyNastran.dev.bdf_vectorized3.test.test_numpy_utils import *

from pyNastran.dev.bdf_vectorized3.cards.test.all_tests import *
#from pyNastran.dev.bdf_vectorized3.solver.test_solver_springs import *

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
