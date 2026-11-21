import unittest

from pyNastran.converters.nastran.test_nastran import *

from pyNastran.converters.cart3d.test_cart3d import *
from pyNastran.converters.fast.test_fast import *
from pyNastran.converters.fld.test_fld import *
from pyNastran.converters.fluent.test_fluent import *
from pyNastran.converters.neu.test_neu import *
from pyNastran.converters.panair.test_panair import *
from pyNastran.converters.stl.test_stl import *
from pyNastran.converters.tecplot.test_tecplot import *
from pyNastran.converters.abaqus.test_unit_abaqus import *
from pyNastran.converters.su2.test_su2 import *
from pyNastran.converters.shabp.test_shabp import *
from pyNastran.converters.usm3d.test_usm3d import *

from pyNastran.converters.aflr.aflr2.test_bedge import *
from pyNastran.converters.aflr.ugrid.test_ugrid import *
from pyNastran.converters.aflr.surf.test_surf import *

#try:
from pyNastran.converters.dev.avus.test_avus import *
#except ModuleNotFoundError:  # pragma: no cover
    #pass

try:
    from pyNastran.converters.dev.vrml.test_vrml import *
except ModuleNotFoundError:  # pragma: no cover
    pass

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
