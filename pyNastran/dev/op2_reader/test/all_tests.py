from pyNastran.op2.test.op2_unit_tests import TestOP2
from pyNastran.op2.test.matrices.test_matrices import TestOP2Matrix
from pyNastran.op2.test.examples.test_op2_in_material_coord import TestMaterialCoordReal
from pyNastran.op2.test.examples.test_op2_in_material_coord_panel_SOL_108 import TestMaterialCoordComplex
from pyNastran.op2.tables.geom.test.test_geom import TestOP2GeomUnit
from pyNastran.op2.tables.test.test_grid_point_forces import TestGridPointForces


if __name__ == "__main__":  # pragma: no cover
    import unittest
    unittest.main()
