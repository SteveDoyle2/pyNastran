import unittest

from pyNastran.converters.nastran.test_nastran_gui import *
from pyNastran.converters.avl.test_avl_gui import *
from pyNastran.converters.abaqus.test_abaqus_gui import *
from pyNastran.converters.cart3d.test_cart3d_gui import *
from pyNastran.converters.fast.test_fast_gui import *
from pyNastran.converters.panair.test_panair_gui import *
from pyNastran.converters.shabp.test_shabp_gui import *
from pyNastran.converters.stl.test_stl_gui import *
from pyNastran.converters.tecplot.test_tecplot_gui import *
from pyNastran.converters.lawgs.test_wgs_gui import *
from pyNastran.converters.su2.test_su2_gui import *
from pyNastran.converters.tetgen.test_tetgen_gui import *
from pyNastran.converters.usm3d.test_usm3d_gui import *
from pyNastran.converters.fluent.test_fluent_gui import *
from pyNastran.converters.fld.test_fld_gui import *

from pyNastran.converters.aflr.aflr2.test_bedge_gui import *
from pyNastran.converters.aflr.surf.test_surf_gui import *
from pyNastran.converters.aflr.ugrid.test_ugrid_gui import *

#try:
from pyNastran.converters.dev.avus.test_avus_gui import *
from pyNastran.converters.dev.openvsp.test_openvsp_gui import *

#except ModuleNotFoundError:
    #pass

#try:
from pyNastran.converters.dev.obj.test_obj import *
#except ModuleNotFoundError:
    #pass

try:
    from pyNastran.converters.dev.vrml.test_vrml_gui import *
except ModuleNotFoundError:  # pragma: no cover
    pass

if __name__ == '__main__':  # pragma: no cover
    import os
    res = unittest.main(exit=False)
    if res.result.wasSuccessful():
        os._exit(0)
    else:
        os._exit(1)
