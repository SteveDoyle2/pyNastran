from pyNastran.dev.bdf_vectorized3.mesh_utils.test.test_vector_equivalence import TestEquiv
#from pyNastran.dev.bdf_vectorized3.mesh_utils.test.test_convert import TestConvert
#from pyNastran.dev.bdf_vectorized3.mesh_utils.test.test_cutting_plane import TestCuttingPlane
#from pyNastran.dev.bdf_vectorized3.mesh_utils.test.test_mass import TestMass
#from pyNastran.dev.bdf_vectorized3.mesh_utils.test.test_mesh_quality import TestMeshQuality
#from pyNastran.dev.bdf_vectorized3.mesh_utils.test.test_mesh_utils import TestMeshUtils
#from pyNastran.dev.bdf_vectorized3.mesh_utils.test.test_renumber import TestRenumber
from pyNastran.dev.bdf_vectorized3.mesh_utils.test.test_vector_remove_unused import TestRemoveUnused
#from pyNastran.dev.bdf_vectorized3.mesh_utils.test.test_sum_loads import TestLoadSum
#from pyNastran.dev.bdf_vectorized3.mesh_utils.test.test_refine import TestRefine
#from pyNastran.dev.bdf_vectorized3.mesh_utils.test.test_flutter import TestFlutter

if __name__ == "__main__":  # pragma: no cover
    import os
    import unittest
    on_rtd = os.environ.get('READTHEDOCS', None)
    if on_rtd is None:
        unittest.main()
