import unittest

from pyNastran_gui.converters.nastran.test_nastran_gui import TestNastranGUI
from pyNastran_gui.converters.avl.test_avl_gui import TestAvlGUI
from pyNastran_gui.converters.abaqus.test_abaqus_gui import TestAbaqusGui
from pyNastran_gui.converters.cart3d.test_cart3d_gui import TestCart3dGUI
from pyNastran_gui.converters.fast.test_fast_gui import TestFastGUI
from pyNastran_gui.converters.panair.test_panair_gui import TestPanairGUI
from pyNastran_gui.converters.shabp.test_shabp_gui import TestShabpGUI
from pyNastran_gui.converters.stl.test_stl_gui import STL_GUITest
from pyNastran_gui.converters.tecplot.test_tecplot_gui import TestTecplotGUI
from pyNastran_gui.converters.lawgs.test_wgs_gui import TestLawgsGUI
from pyNastran_gui.converters.su2.test_su2_gui import TestSU2GUI
from pyNastran_gui.converters.tetgen.test_tetgen_gui import TestTetgenGUI
from pyNastran_gui.converters.usm3d.test_usm3d_gui import TestUsm3dGUI
from pyNastran_gui.converters.openfoam.test_openfoam_gui import TestOpenFoamGUI

from pyNastran_gui.converters.aflr.aflr2.test_bedge_gui import TestBEdgeGUI
from pyNastran_gui.converters.aflr.surf.test_surf_gui import TestSurfGui
from pyNastran_gui.converters.aflr.ugrid.test_ugrid_gui import TestUgridGui

#try:
from pyNastran_gui.converters.dev.avus.test_avus_gui import TestAvusGUI
from pyNastran_gui.converters.dev.openvsp.test_openvsp_gui import TestOpenVSP_GUI

#except ImportError:
    #pass

#try:
from pyNastran_gui.converters.dev.obj.test_obj import TestObjGUI
#except ImportError:
    #pass

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
