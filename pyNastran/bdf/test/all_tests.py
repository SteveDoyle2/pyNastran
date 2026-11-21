# unit
from pyNastran.bdf.test.test_field_writer import *
from pyNastran.bdf.test.test_utils import *
from pyNastran.bdf.test.test_bdf_unit_tests import *
from pyNastran.bdf.test.unit.test_read_write import *
from pyNastran.bdf.test.unit.test_parsing import *
from pyNastran.bdf.test.test_openmdao import *

from pyNastran.bdf.cards.test.all_tests import *
from pyNastran.bdf.mesh_utils.test.all_tests import *

from pyNastran.bdf.patran_utils.test_patran import *

# bdf_interface
from pyNastran.bdf.bdf_interface.dev.test_mass import *
from pyNastran.bdf.bdf_interface.test.test_pybdf import *
from pyNastran.bdf.bdf_interface.test.test_assign_type import *
from pyNastran.bdf.bdf_interface.test.test_bdf_interface import *
from pyNastran.bdf.bdf_interface.test.test_dev_utils import *
from pyNastran.bdf.bdf_interface.test.test_case_control_deck import *


if __name__ == "__main__":  # pragma: no cover
    import os
    import unittest
    on_rtd = os.environ.get('READTHEDOCS', None)
    if on_rtd is None:
        unittest.main()
