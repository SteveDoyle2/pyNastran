from __future__ import print_function

from pyNastran.bdf.cards.test.test_coords import *
from pyNastran.bdf.cards.test.test_constraints import *
from pyNastran.bdf.cards.test.test_sets import *

# standard elements
from pyNastran.bdf.cards.test.test_rods import *
from pyNastran.bdf.cards.test.test_bars import *
from pyNastran.bdf.cards.test.test_beams import *
from pyNastran.bdf.cards.test.test_dmig import *
from pyNastran.bdf.cards.test.test_springs import *
from pyNastran.bdf.cards.test.test_solids import *
from pyNastran.bdf.cards.test.test_shells import *
from pyNastran.bdf.cards.test.test_contact import *

# other elements
from pyNastran.bdf.cards.test.test_rigid import *
from pyNastran.bdf.cards.test.test_elements import *

from pyNastran.bdf.cards.test.test_loads import *
from pyNastran.bdf.cards.test.test_materials import *

from pyNastran.bdf.cards.test.test_tables import *


if __name__ == "__main__":  # pragma: no cover
    import unittest
    unittest.main()

