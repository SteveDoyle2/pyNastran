from pyNastran.bdf.test.test_field_writer import Testfield_writer_8
from pyNastran.bdf.test.test_bdf_unit_tests import TestBDFUnit
from pyNastran.bdf.bdf_interface.dev.test_mass import TestMassGeneration

from pyNastran.bdf.cards.test.all_tests import *
from pyNastran.bdf.mesh_utils.test.all_tests import *

from pyNastran.bdf.patran_utils.test_patran import TestPatran, TestPatranSyntax

# unit
from pyNastran.bdf.test.unit.test_read_write import TestReadWrite, TestReadWriteFiles
from pyNastran.bdf.test.unit.test_parsing import TestBDFParsing
from pyNastran.bdf.test.test_openmdao import TestOpenMDAO

# bdf_interface
from pyNastran.bdf.bdf_interface.test.test_pybdf import TestPyBDF
from pyNastran.bdf.bdf_interface.test.test_assign_type import TestAssignType
from pyNastran.bdf.bdf_interface.test.test_bdf_interface import TestBDFInterface
from pyNastran.bdf.bdf_interface.test.test_dev_utils import DevUtils
from pyNastran.bdf.bdf_interface.test.test_case_control_deck import CaseControlTest


if __name__ == "__main__":  # pragma: no cover
    import os
    import unittest
    on_rtd = os.environ.get('READTHEDOCS', None)
    if on_rtd is None:
        unittest.main()

