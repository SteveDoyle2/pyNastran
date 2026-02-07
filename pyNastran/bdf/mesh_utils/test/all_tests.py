from pyNastran.bdf.mesh_utils.cut.test_cutting_plane import TestCuttingPlane
from pyNastran.bdf.mesh_utils.cut.test_stiffness_plot import TestStiffnessPlot

from pyNastran.bdf.mesh_utils.aero.test_flutter import TestBDFFlutter
from pyNastran.bdf.mesh_utils.aero.test_aero_utils import TestMeshUtilsAero

from pyNastran.bdf.mesh_utils.test.test_equivalence import TestEquiv
from pyNastran.bdf.mesh_utils.test.test_convert import TestConvert
from pyNastran.bdf.mesh_utils.test.test_mass import TestMass
from pyNastran.bdf.mesh_utils.test.test_mesh_quality import TestMeshQuality
from pyNastran.bdf.mesh_utils.test.test_mesh_utils import (
    TestMeshUtils, TestMeshUtilsCmdLine, TestRbeTools)
from pyNastran.bdf.mesh_utils.test.test_refine import TestRefine
from pyNastran.bdf.mesh_utils.test.test_remove_unused import TestRemoveUnused
from pyNastran.bdf.mesh_utils.test.test_renumber import TestRenumber
from pyNastran.bdf.mesh_utils.test.test_run_host_jobs import TestRunHostJobs
from pyNastran.bdf.mesh_utils.test.test_sum_loads import TestLoadSum


if __name__ == "__main__":  # pragma: no cover
    import os
    import unittest
    on_rtd = os.environ.get('READTHEDOCS', None)
    if on_rtd is None:
        unittest.main()
