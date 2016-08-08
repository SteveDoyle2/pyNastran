from __future__ import print_function

from pyNastran.bdf.test.test_utils import *
from pyNastran.bdf.test.test_field_writer import Testfield_writer_8
from pyNastran.bdf.test.bdf_unit_tests import TestBDF, TestBaseCard
from pyNastran.bdf.test.test_case_control_deck import CaseControlTest

from pyNastran.bdf.cards.test.all_tests import *
from pyNastran.bdf.bdf_interface.test.test_dev_utils import DevUtils
from pyNastran.bdf.mesh_utils.test_mesh_utils import TestMeshUtils
# unit
from pyNastran.bdf.test.unit.test_mass import *
from pyNastran.bdf.test.unit.test_assign_type import *
from pyNastran.bdf.test.unit.test_read_write import *
from pyNastran.bdf.test.unit.test_sum_loads import *


if __name__ == "__main__":  # pragma: no cover
    import os
    import unittest
    on_rtd = os.environ.get('READTHEDOCS', None)
    if on_rtd is None:
        unittest.main()

