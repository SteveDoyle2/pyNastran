"""import various codes with backup for failed imports"""
#try:
from pyNastran.converters.cart3d.cart3d_io import Cart3dIO
is_cart3d = True
#except ImportError:
    #raise
    #class Cart3dIO(object):
    #    """dummy cart3d gui class"""
    #    def __init__(self):
    #        """dummy gui init"""
    #        pass
    #is_cart3d = False

#try:
from pyNastran.converters.nastran.nastranIOv import NastranIO
is_nastran = True
#except ImportError:
    #raise
    #class NastranIO(object):
    #    """dummy nastran gui class"""
    #    def __init__(self):
    #        """dummy gui init"""
    #        pass
    #    def load_nastran_geometry(self, infile_name, dirname):
    #        pass
    #is_nastran = False


#try:
from pyNastran.converters.LaWGS.wgs_io import LaWGS_IO
is_lawgs = True
#except ImportError:
    #raise
    #class LaWGS_IO(object):
    #    """dummy lawgs gui class"""
    #    def __init__(self):
    #        """dummy gui init"""
    #        pass
    #    def load_lawgs_geometry(self, infile_name, dirname):
    #        pass
    #is_lawgs = False


#try:
from pyNastran.converters.panair.panair_io import PanairIO
is_panair = True
#except ImportError:
    #raise
    #class PanairIO(object):
    #    """dummy panair gui class"""
    #    def __init__(self):
    #        """dummy gui init"""
    #        pass
    #is_panair = False

try:
    from pyNastran.converters.dev.plot3d.plot3d_io import Plot3d_io
    is_plot3d = True
#except ImportError:
except:
    class Plot3d_io(object):
        """dummy plot3d gui class"""
        def __init__(self):
            """dummy gui init"""
            pass
    is_plot3d = False
    #raise

#try:
from pyNastran.converters.shabp.shabp_io import ShabpIO
is_shabp = True
#except ImportError:
    #raise
    #class ShabpIO(object):
        #"""dummy shabp gui class"""
        #def __init__(self):
            #"""dummy gui init"""
            #pass
    #is_shabp = False

#try:
from pyNastran.converters.stl.stl_io import STL_IO
is_stl = True
#except ImportError:
    #raise
    #class STL_IO(object):
        #"""dummy stl gui class"""
        #def __init__(self):
            #"""dummy gui init"""
            #pass
    #is_stl = False


#try:
from pyNastran.converters.tecplot.tecplot_io import TecplotIO
    #is_tecplot = True
#except ImportError:
    #raise
    #class TecplotIO(object):
    #    """dummy tecplot gui class"""
    #    def __init__(self):
    #        """dummy gui init"""
    #        pass
    #is_tecplot = False

try:
    from pyNastran.converters.dev.avus.avus_io import AvusIO
    is_avus = True
except ImportError:
    #raise
    class AvusIO(object):
        """dummy avus gui class"""
        def __init__(self):
            """dummy gui init"""
            pass
    is_avus = False

#try:
from pyNastran.converters.tetgen.tetgen_io import TetgenIO
is_tetgen = True
#except ImportError:
    #raise
    #class TetgenIO(object):
        #"""dummy tetgen gui class"""
        #def __init__(self):
            #"""dummy gui init"""
            #pass
    #is_tetgen = False

#try:
from pyNastran.converters.usm3d.usm3d_io import Usm3dIO
is_usm3d = True
#except ImportError:
    #raise
    #class Usm3dIO(object):
        #"""dummy usm3d gui class"""
        #def __init__(self):
            #"""dummy gui init"""
            #pass
    #is_usm3d = False


#try:
from pyNastran.converters.fast.fast_io import FastIO
is_fast = True
#except ImportError:
    #raise
    #class FastIO(object):
    #    """dummy fast gui class"""
    #    def __init__(self):
    #        """dummy gui init"""
    #        pass
    #is_fast = False

#try:
from pyNastran.converters.aflr.aflr2.bedge_io import BEdge_IO
is_bedge = True
#except ImportError:
    #raise
    #class BEdge_IO(object):
        #"""dummy bedge gui class"""
        #def __init__(self):
            #"""dummy gui init"""
            #pass
    #is_bedge = False

#try:
from pyNastran.converters.aflr.surf.surf_io import SurfIO
is_surf = True
#except ImportError:
    #raise
    #class SurfIO(object):
        #"""dummy surf gui class"""
        #def __init__(self):
            #"""dummy gui init"""
            #pass
    #is_surf = False

#try:
from pyNastran.converters.aflr.ugrid.ugrid_io import UGRID_IO
is_ugrid = True
#except ImportError:
    #raise
    #class UGRID_IO(object):
        #"""dummy ugrid gui class"""
        #def __init__(self):
            #"""dummy gui init"""
            #pass
    #is_ugrid = False

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
    from pyNastran.converters.su2.su2_io import SU2_IO
    SU2_IO = True
except ImportError:
    class SU2_IO(object):
        """dummy SU2 gui class"""
        def __init__(self):
            """dummy gui init"""
            pass
    is_su2 = False
    #raise
