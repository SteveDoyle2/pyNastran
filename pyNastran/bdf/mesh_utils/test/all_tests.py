from pyNastran.bdf.mesh_utils.test.test_equivalence import TestEquiv
from pyNastran.bdf.mesh_utils.test.test_convert import TestConvert
from pyNastran.bdf.mesh_utils.test.test_cutting_plane import TestCuttingPlane
from pyNastran.bdf.mesh_utils.test.test_mass import TestMass
from pyNastran.bdf.mesh_utils.test.test_mesh_quality import TestMeshQuality
from pyNastran.bdf.mesh_utils.test.test_mesh_utils import TestMeshUtils
from pyNastran.bdf.mesh_utils.test.test_renumber import TestRenumber
from pyNastran.bdf.mesh_utils.test.test_remove_unused import TestRemoveUnused
from pyNastran.bdf.mesh_utils.test.test_sum_loads import TestLoadSum
from pyNastran.bdf.mesh_utils.test.test_refine import TestRefine
from pyNastran.bdf.mesh_utils.test.test_flutter import TestFlutter

if __name__ == "__main__":  # pragma: no cover
    import os
    import unittest
    on_rtd = os.environ.get('READTHEDOCS', None)
    if on_rtd is None:
        unittest.main()
