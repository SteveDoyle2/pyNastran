try:
    from pyNastran.gui.nastranIO import NastranIO
    is_nastran = True
except ImportError:
    class NastranIO(object):
        def __init__(self):
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
    class PanairIO(object):
        def __init__(self):
            pass
    is_cart3d = False