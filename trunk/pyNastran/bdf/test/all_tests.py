from __future__ import print_function

from pyNastran.bdf.test.test_field_writer import TestFieldWriter
from pyNastran.bdf.test.bdf_unit_tests import TestBDF, BaseCard_Test
from pyNastran.bdf.test.test_case_control_deck import CaseControlTest

from pyNastran.bdf.test.cards.test_coords import *
from pyNastran.bdf.test.cards.test_constraints import *
from pyNastran.bdf.test.cards.test_sets import *

# standard elements
from pyNastran.bdf.test.cards.test_rods import *
from pyNastran.bdf.test.cards.test_bars import *
from pyNastran.bdf.test.cards.test_beams import *
from pyNastran.bdf.test.cards.test_dmig import *
from pyNastran.bdf.test.cards.test_springs import *
from pyNastran.bdf.test.cards.test_solids import *
from pyNastran.bdf.test.cards.test_shells import *
from pyNastran.bdf.test.cards.test_contact import *

# other elements
from pyNastran.bdf.test.cards.test_rigid import *
from pyNastran.bdf.test.cards.test_elements import *

from pyNastran.bdf.test.cards.test_loads import *
from pyNastran.bdf.test.cards.test_materials import *

from pyNastran.bdf.test.cards.test_tables import *

# unit
from pyNastran.bdf.test.unit.test_mass import *
from pyNastran.bdf.test.unit.test_assign_type import *
from pyNastran.bdf.test.unit.test_read_write import *
from pyNastran.bdf.test.unit.test_sum_loads import *


if __name__ == "__main__":  # pragma: no cover
    import unittest
    unittest.main()

