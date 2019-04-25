import unittest

from pyNastran_gui.converters.nastran.test_nastran import TestNastran

from pyNastran_gui.converters.cart3d.test_cart3d import TestCart3d
from pyNastran_gui.converters.fast.test_fast import TestFast
from pyNastran_gui.converters.panair.test_panair import TestPanair
from pyNastran_gui.converters.stl.test_stl import TestSTL
from pyNastran_gui.converters.tecplot.test_tecplot import TestTecplot
from pyNastran_gui.converters.abaqus.test_unit_abaqus import TestAbaqus
from pyNastran_gui.converters.usm3d.test_usm3d import TestUsm3d

from pyNastran_gui.converters.aflr.aflr2.test_bedge import TestBEdge
from pyNastran_gui.converters.aflr.ugrid.test_ugrid import TestUgrid

#try:
from pyNastran_gui.converters.dev.avus.test_avus import TestAvus
#except ImportError:  # pragma: no cover
    #pass

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
