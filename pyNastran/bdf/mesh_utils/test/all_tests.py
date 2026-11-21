try:
    import scipy
    IS_SCIPY = True
except ImportError:
    IS_SCIPY = False
from pyNastran.bdf.mesh_utils.cut.test_cutting_plane import *
from pyNastran.bdf.mesh_utils.cut.test_stiffness_plot import *

from pyNastran.bdf.mesh_utils.aero.test_flutter import *
from pyNastran.bdf.mesh_utils.aero.test_aero_utils import *

if IS_SCIPY:
    from pyNastran.bdf.mesh_utils.test.test_equivalence import *
    from pyNastran.bdf.mesh_utils.test.test_mesh_utils import *
from pyNastran.bdf.mesh_utils.test.test_convert import *
from pyNastran.bdf.mesh_utils.test.test_mass import *
from pyNastran.bdf.mesh_utils.test.test_mesh_quality import *
from pyNastran.bdf.mesh_utils.test.test_refine import *
from pyNastran.bdf.mesh_utils.test.test_remove_unused import *
from pyNastran.bdf.mesh_utils.test.test_renumber import *
from pyNastran.bdf.mesh_utils.test.test_run_host_jobs import *
from pyNastran.bdf.mesh_utils.test.test_sum_loads import *
from pyNastran.bdf.mesh_utils.test.test_find_remove_duplicates import *


if __name__ == "__main__":  # pragma: no cover
    import os
    import unittest
    on_rtd = os.environ.get('READTHEDOCS', None)
    if on_rtd is None:
        unittest.main()
