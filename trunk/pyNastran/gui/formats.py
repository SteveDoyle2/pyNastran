try:
    from pyNastran.converters.cart3d.cart3dIO import Cart3dIO
    is_cart3d = True
except ImportError:
    raise
    class Cart3dIO(object):
        def __init__(self):
            pass
    is_cart3d = False

try:
    from pyNastran.converters.nastran.nastranIO import NastranIO
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
    from pyNastran.converters.LaWGS.wgsIO import LaWGS_IO
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
    from pyNastran.converters.panair.panairIO import PanairIO
    is_panair = True
except ImportError:
    raise
    class PanairIO(object):
        def __init__(self):
            pass
    is_panair = False

try:
    from pyNastran.converters.plot3d.plot3d_io import Plot3d_io
    is_plot3d = True
except ImportError:
    class Plot3d_io(object):
        def __init__(self):
            pass
    is_plot3d = False
    raise

try:
    from pyNastran.converters.shabp.shabp_io import ShabpIO
    is_shabp = True
except ImportError:
    raise
    class ShabpIO(object):
        def __init__(self):
            pass
    is_shabp = False

try:
    from pyNastran.converters.stl.stl_io import STL_IO
    is_stl = True
except ImportError:
    raise
    class STL_IO(object):
        def __init__(self):
            pass
    is_stl = False


try:
    from pyNastran.converters.tecplot.tecplot_io import TecplotIO
    is_tecplot = True
except ImportError:
    #raise
    class TecplotIO(object):
        def __init__(self):
            pass
    is_tecplot = False

try:
    from pyNastran.converters.tetgen.tetgen_io import TetgenIO
    is_tetgen = True
except ImportError:
    raise
    class TetgenIO(object):
        def __init__(self):
            pass
    is_tetgen = False

try:
    from pyNastran.converters.usm3d.usm3d_io import Usm3dIO
    is_usm3d = True
except ImportError:
    raise
    class Usm3dIO(object):
        def __init__(self):
            pass
    is_usm3d = False

