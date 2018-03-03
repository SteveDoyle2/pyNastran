"""import various codes with backup for failed imports"""
CLASS_MAP = {}

try:
    from pyNastran.converters.dev.avus.avus_io import AvusIO
    CLASS_MAP['avus'] = AvusIO
except ImportError:  # pragma: no cover
    pass

try:
    from pyNastran.converters.cart3d.cart3d_io import Cart3dIO
    CLASS_MAP['cart3d'] = Cart3dIO
except ImportError:  # pragma: no cover
    pass

try:
    from pyNastran.converters.panair.panair_io import PanairIO
    CLASS_MAP['panair'] = PanairIO
except ImportError:  # pragma: no cover
    pass

try:
    from pyNastran.converters.stl.stl_io import STL_IO
    CLASS_MAP['stl'] = STL_IO
except ImportError:  # pragma: no cover
    pass

try:
    from pyNastran.converters.usm3d.usm3d_io import Usm3dIO
    CLASS_MAP['usm3d'] = Usm3dIO
except ImportError:  # pragma: no cover
    pass

try:
    from pyNastran.converters.tecplot.tecplot_io import TecplotIO
    CLASS_MAP['tecplot'] = TecplotIO
except ImportError:  # pragma: no cover
    pass

try:
    from pyNastran.converters.lawgs.wgs_io import LaWGS_IO
    CLASS_MAP['lawgs'] = LaWGS_IO
except ImportError:  # pragma: no cover
    pass

try:
    from pyNastran.converters.tetgen.tetgen_io import TetgenIO
    CLASS_MAP['tetgen'] = TetgenIO
except ImportError:  # pragma: no cover
    pass

try:
    from pyNastran.converters.shabp.shabp_io import ShabpIO
    CLASS_MAP['shabp'] = ShabpIO
except ImportError:  # pragma: no cover
    pass

try:
    from pyNastran.converters.su2.su2_io import SU2_IO
    CLASS_MAP['su2'] = SU2_IO
except ImportError:
    pass


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
