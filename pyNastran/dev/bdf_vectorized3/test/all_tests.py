import unittest
try:
    import vtkmodules
    USE_VTK = True
except ImportError:
    USE_VTK = False

from pyNastran.dev.bdf_vectorized3.test.test_vector_models import *
from pyNastran.dev.bdf_vectorized3.mesh_utils.test.all_tests import *
if USE_VTK:
    from pyNastran.dev.bdf_vectorized3.test.test_vector_gui_models import *

from pyNastran.dev.bdf_vectorized3.test.test_vector_bdf_unit_tests import *
from pyNastran.dev.bdf_vectorized3.test.test_vector_numpy_utils import *

from pyNastran.dev.bdf_vectorized3.cards.test.all_tests import *
from pyNastran.dev.bdf_vectorized3.solver.test_vector_solver_springs import *
from pyNastran.dev.op2_vectorized3.test.test_vector_op2_unit import *


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
