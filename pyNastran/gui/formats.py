try:
    from pyNastran.converters.cart3d.cart3d_io import Cart3dIO
    is_cart3d = True
except ImportError:
    raise
    class Cart3dIO(object):
        def __init__(self):
            pass
    is_cart3d = False

try:
    from pyNastran.converters.nastran.nastranIOv import NastranIO
    is_nastran = True
except ImportError:
    raise
    class NastranIO(object):
        def __init__(self):
            pass
        def load_nastran_geometry(self, infile_name, dirname):
            pass
    is_nastran = False


try:
    from pyNastran.converters.LaWGS.wgs_io import LaWGS_IO
    is_lawgs = True
except ImportError:
    raise
    class LaWGS_IO(object):
        def __init__(self):
            pass
        def load_lawgs_geometry(self, infile_name, dirname):
            pass
    is_lawgs = False


try:
    from pyNastran.converters.panair.panair_io import PanairIO
    is_panair = True
except ImportError:
    raise
    class PanairIO(object):
        def __init__(self):
            pass
    is_panair = False

try:
    from pyNastran.converters.shabp.shabp_io import ShabpIO
    is_shabp = True
except ImportError:
    #raise
    class ShabpIO(object):
        def __init__(self):
            pass
    is_shabp = False

try:
    from pyNastran.converters.stl.stl_io import STL_IO
    is_stl = True
except ImportError:
    #raise
    class STL_IO(object):
        def __init__(self):
            pass
    is_stl = False


try:
    from pyNastran.converters.tecplot.tecplot_io import TecplotIO
    is_tecplot = True
except ImportError:
    raise
    class TecplotIO(object):
        def __init__(self):
            pass
    is_tecplot = False

try:
    from pyNastran.converters.ugrid.surf_io import SurfIO
    is_surf = True
except ImportError:
    #raise
    class SurfIO(object):
        def __init__(self):
            pass
    is_surf = False

try:
    from pyNastran.converters.ugrid.ugrid_io import UGRID_IO
    is_ugrid = True
except ImportError:
    #raise
    class UGRID_IO(object):
        def __init__(self):
            pass
    is_ugrid = False
