import unittest

from pyNastran.converters.nastran.test_nastran import TestNastran

from pyNastran.converters.cart3d.test_cart3d import TestCart3d
from pyNastran.converters.fast.test_fast import TestFast
from pyNastran.converters.openfoam.test_openfoam import TestOpenFOAM
from pyNastran.converters.panair.test_panair import TestPanair
from pyNastran.converters.stl.test_stl import TestSTL
from pyNastran.converters.tecplot.test_tecplot import TestTecplot
from pyNastran.converters.abaqus.test_unit_abaqus import TestAbaqus
from pyNastran.converters.su2.test_su2 import TestSU2
from pyNastran.converters.shabp.test_shabp import TestShabp
from pyNastran.converters.usm3d.test_usm3d import TestUsm3d

from pyNastran.converters.aflr.aflr2.test_bedge import TestBEdge
from pyNastran.converters.aflr.ugrid.test_ugrid import TestUgrid


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
