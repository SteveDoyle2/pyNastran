try:
    from pyNastran.gui.nastranIO import NastranIO
    is_nastran = True
except ImportError:
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
    class LaWGS_IO(object):
        def __init__(self):
            pass
    is_lawgs = False


try:
    from pyNastran.converters.panair.panairIO import PanairIO
    is_panair = True
except ImportError:
    class PanairIO(object):
        def __init__(self):
            pass
    is_panair = False

try:
    from pyNastran.converters.cart3d.cart3dIO import Cart3dIO
    is_cart3d = True
except ImportError:
    class Cart3dIO(object):
        def __init__(self):
            pass
    is_cart3d = False

try:
    from pyNastran.converters.stl.stl_io import STL_IO
    is_stl = True
except ImportError:
    class STL_IO(object):
        def __init__(self):
            pass
    is_stl = False


try:
    from pyNastran.converters.tetgen.tetgen_io import TetgenIO
    is_tetgen = True
except ImportError:
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