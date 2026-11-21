from pyNastran.op2.test.test_op2_unit_tests import *
from pyNastran.op2.test.test_op2_no_scipy import *
from pyNastran.op2.test.matrices.test_matrices import *
from pyNastran.op2.test.examples.test_op2_in_material_coord import *
from pyNastran.op2.test.examples.test_op2_in_material_coord_panel_SOL_108 import *
from pyNastran.op2.tables.geom.test.test_geom import *
from pyNastran.op2.tables.test.test_grid_point_forces import *
from pyNastran.op2.tables.test.test_oug import *
from pyNastran.op2.writer.test_op2_writer import *
from pyNastran.op2.op2_interface.test.test_results_set import *
from pyNastran.op2.tools.test_combination import *
from pyNastran.op2.test.test_op2_bugfixes import *


if __name__ == "__main__":  # pragma: no cover
    import unittest
    unittest.main()
