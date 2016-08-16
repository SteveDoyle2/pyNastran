import unittest

from pyNastran.converters.nastran.test_nastran_gui import TestNastranGUI
from pyNastran.converters.cart3d.test_cart3d import TestCart3d
from pyNastran.converters.cart3d.test_cart3d_gui import TestCart3dGUI
from pyNastran.converters.panair.test_panair import TestPanair
from pyNastran.converters.panair.test_panair_gui import TestPanairGUI
from pyNastran.converters.shabp.test_shabp_gui import TestShabpGUI
from pyNastran.converters.stl.test_stl import TestSTL
from pyNastran.converters.stl.test_stl_gui import STL_GUITest
from pyNastran.converters.tecplot.test_tecplot_gui import TestTecplotGUI
from pyNastran.converters.LaWGS.test_wgs_gui import TestLawgsGUI

if __name__ == '__main__':
    unittest.main()
