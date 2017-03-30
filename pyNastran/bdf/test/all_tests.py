from __future__ import print_function

from pyNastran.bdf.test.test_utils import TestBdfUtils
from pyNastran.bdf.test.test_field_writer import Testfield_writer_8
from pyNastran.bdf.test.bdf_unit_tests import TestBDF, TestBaseCard
from pyNastran.bdf.test.test_case_control_deck import CaseControlTest

from pyNastran.bdf.cards.test.all_tests import *
from pyNastran.bdf.mesh_utils.test.test_mesh_utils import TestMeshUtils
from pyNastran.bdf.mesh_utils.test.test_convert import TestConvert
from pyNastran.bdf.mesh_utils.test.test_remove_unused import TestRemoveUnused
from pyNastran.bdf.patran_utils.test_patran import TestPatran

# unit
from pyNastran.bdf.test.unit.test_mass import TestMass
from pyNastran.bdf.test.unit.test_read_write import TestReadWrite
from pyNastran.bdf.test.unit.test_sum_loads import TestLoadSum
from pyNastran.bdf.test.test_openmdao import TestOpenMDAO

# bdf_interface
from pyNastran.bdf.bdf_interface.test.test_assign_type import TestAssignType
from pyNastran.bdf.bdf_interface.test.test_dev_utils import DevUtils


if __name__ == "__main__":  # pragma: no cover
    import os
    import unittest
    on_rtd = os.environ.get('READTHEDOCS', None)
    if on_rtd is None:
        unittest.main()

