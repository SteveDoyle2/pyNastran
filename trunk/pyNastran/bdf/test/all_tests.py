from __future__ import print_function

from pyNastran.bdf.test.test_field_writer import TestFieldWriter
from pyNastran.bdf.test.bdf_unit_tests import TestBDF, BaseCard_Test
from pyNastran.bdf.test.test_case_control_deck import CaseControlTest

from pyNastran.bdf.cards.test.all_tests import *

# unit
from pyNastran.bdf.test.unit.test_mass import *
from pyNastran.bdf.test.unit.test_assign_type import *
from pyNastran.bdf.test.unit.test_read_write import *
from pyNastran.bdf.test.unit.test_sum_loads import *


if __name__ == "__main__":  # pragma: no cover
    import unittest
    unittest.main()

