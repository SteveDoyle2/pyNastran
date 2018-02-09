"""import various codes with backup for failed imports"""
#try:

from pyNastran.converters.dev.avus.avus_io import AvusIO
from pyNastran.converters.cart3d.cart3d_io import Cart3dIO
from pyNastran.converters.panair.panair_io import PanairIO
from pyNastran.converters.stl.stl_io import STL_IO
from pyNastran.converters.usm3d.usm3d_io import Usm3dIO
from pyNastran.converters.tecplot.tecplot_io import TecplotIO
from pyNastran.converters.lawgs.wgs_io import LaWGS_IO
from pyNastran.converters.tetgen.tetgen_io import TetgenIO
from pyNastran.converters.shabp.shabp_io import ShabpIO
from pyNastran.converters.su2.su2_io import SU2_IO

CLASS_MAP = {
    'cart3d' : Cart3dIO,
    'avus' : AvusIO,
    'panair' : PanairIO,
    'stl' : STL_IO,
    'usm3d' : Usm3dIO,
    'tecplot' : TecplotIO,
    'lawgs' : LaWGS_IO,
    'tetgen' : TetgenIO,
    'shabp' : ShabpIO,
    'su2' : SU2_IO,
}


from pyNastran.converters.nastran.nastran_io import NastranIO
is_nastran = True


#from pyNastran.converters.dev.plot3d.plot3d_io import Plot3d_io



from pyNastran.converters.fast.fast_io import FastIO
is_fast = True

from pyNastran.converters.aflr.aflr2.bedge_io import BEdge_IO
is_bedge = True

from pyNastran.converters.aflr.surf.surf_io import SurfIO
is_surf = True

from pyNastran.converters.aflr.ugrid.ugrid_io import UGRID_IO
is_ugrid = True

try:
    from pyNastran.converters.dev.openvsp.adb_io import ADB_IO
    is_openvsp = True
except ImportError:
    #raise
    class ADB_IO(object):
        """dummy adb gui class"""
        def __init__(self):
            """dummy gui init"""
            pass
    is_openvsp = False

try:
    from pyNastran.converters.dev.openvsp.degen_geom_io import DegenGeomIO
    is_degen_geom = True
except ImportError:
    #raise
    class DegenGeomIO(object):
        """dummy degen_geom gui class"""
        def __init__(self):
            """dummy gui init"""
            pass
    is_degen_geom = False

try:
    from pyNastran.converters.abaqus.abaqus_io import AbaqusIO
    is_abaqus = True
except ImportError:
    #raise
    class AbaqusIO(object):
        """dummy abaqus gui class"""
        def __init__(self):
            """dummy gui init"""
            pass
    is_abaqus = False


try:
    from pyNastran.converters.openfoam.openfoam_io import OpenFoamIO
    is_openfoam = True
except ImportError:
    raise

try:
    from pyNastran.converters.dev.obj.obj_io import ObjIO
    is_obj = True
except ImportError:
    class ObjIO(object):
        """dummy SU2 gui class"""
        def __init__(self):
            """dummy gui init"""
            pass
    is_obj = False
