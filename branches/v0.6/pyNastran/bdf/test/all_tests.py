from pyNastran.bdf.test.test_field_writer import TestFieldWriter
from pyNastran.bdf.test.bdf_unit_tests import TestBDF, BaseCard_Test
from pyNastran.bdf.test.test_case_control_deck import CaseControlTest
from pyNastran.bdf.test.cards.test_coords import *
from pyNastran.bdf.test.cards.test_constraints import *
from pyNastran.bdf.test.cards.test_rigid import *
from pyNastran.bdf.test.cards.test_beams import *
from pyNastran.bdf.test.cards.test_materials import *
from pyNastran.bdf.test.cards.test_shells import *

if __name__ == "__main__":
    import unittest
    unittest.main()