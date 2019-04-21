import unittest

from pyNastran.converters.nastran.test_nastran_gui import TestNastranGUI
from pyNastran.converters.avl.test_avl_gui import TestAvlGUI
from pyNastran.converters.abaqus.test_abaqus_gui import TestAbaqusGui
from pyNastran.converters.cart3d.test_cart3d_gui import TestCart3dGUI
from pyNastran.converters.fast.test_fast_gui import TestFastGUI
from pyNastran.converters.panair.test_panair_gui import TestPanairGUI
from pyNastran.converters.shabp.test_shabp_gui import TestShabpGUI
from pyNastran.converters.stl.test_stl_gui import STL_GUITest
from pyNastran.converters.tecplot.test_tecplot_gui import TestTecplotGUI
from pyNastran.converters.lawgs.test_wgs_gui import TestLawgsGUI
from pyNastran.converters.su2.test_su2_gui import TestSU2GUI
from pyNastran.converters.tetgen.test_tetgen_gui import TestTetgenGUI
from pyNastran.converters.usm3d.test_usm3d_gui import TestUsm3dGUI
from pyNastran.converters.openfoam.test_openfoam_gui import TestOpenFoamGUI

from pyNastran.converters.aflr.aflr2.test_bedge_gui import TestBEdgeGUI
from pyNastran.converters.aflr.surf.test_surf_gui import TestSurfGui
from pyNastran.converters.aflr.ugrid.test_ugrid_gui import TestUgridGui


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
