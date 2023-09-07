import unittest
from pyNastran.dev.bdf_vectorized3.test.test_models import *
from pyNastran.dev.bdf_vectorized3.test.bdfv_unit_tests import *
from pyNastran.dev.bdf_vectorized3.test.test_numpy_utils import *

#from pyNastran.dev.bdf_vectorized3.cards.test.test_vector_aero import *
#from pyNastran.dev.bdf_vectorized3.cards.test.test_vector_beams import *
#from pyNastran.dev.bdf_vectorized3.cards.test.test_vector_loads import *

# good
from pyNastran.dev.bdf_vectorized3.cards.test.test_vector_dmig import *
from pyNastran.dev.bdf_vectorized3.cards.test.test_vector_bars import *
#from pyNastran.dev.bdf_vectorized3.cards.test.test_vector_coords import *
from pyNastran.dev.bdf_vectorized3.cards.test.test_vector_shells import TestShells
from pyNastran.dev.bdf_vectorized3.cards.test.test_vector_solids import *

if __name__ == '__main__':
    unittest.main()
