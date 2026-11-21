from pyNastran.dev.bdf_vectorized3.mesh_utils.test.test_inertia_relief import *
from pyNastran.dev.bdf_vectorized3.mesh_utils.test.test_mass_matrix import *
from pyNastran.dev.bdf_vectorized3.mesh_utils.test.test_vector_quality import *
from pyNastran.dev.bdf_vectorized3.mesh_utils.test.test_vector_equivalence import *
#from pyNastran.dev.bdf_vectorized3.mesh_utils.test.test_convert import *
#from pyNastran.dev.bdf_vectorized3.mesh_utils.test.test_cutting_plane import *
#from pyNastran.dev.bdf_vectorized3.mesh_utils.test.test_mass import *
#from pyNastran.dev.bdf_vectorized3.mesh_utils.test.test_mesh_quality import *
#from pyNastran.dev.bdf_vectorized3.mesh_utils.test.test_mesh_utils import *
#from pyNastran.dev.bdf_vectorized3.mesh_utils.test.test_renumber import *
from pyNastran.dev.bdf_vectorized3.mesh_utils.test.test_vector_remove_unused import *
#from pyNastran.dev.bdf_vectorized3.mesh_utils.test.test_sum_loads import *
#from pyNastran.dev.bdf_vectorized3.mesh_utils.test.test_refine import *
#from pyNastran.dev.bdf_vectorized3.mesh_utils.test.test_flutter import *

if __name__ == "__main__":  # pragma: no cover
    import os
    import unittest
    on_rtd = os.environ.get('READTHEDOCS', None)
    if on_rtd is None:
        unittest.main()
