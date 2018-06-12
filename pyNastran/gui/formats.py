"""import various codes with backup for failed imports"""
CLASS_MAP = {}


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
    from pyNastran.converters.fast.fast_io import FastIO
    CLASS_MAP['fast'] = FastIO
except ImportError:  # pragma: no cover
    pass


from pyNastran.converters.nastran.nastran_io import NastranIO
#CLASS_MAP['nastran'] = NastranIO



try:
    from pyNastran.converters.aflr.aflr2.bedge_io import BEdge_IO
    CLASS_MAP['bedge'] = BEdge_IO
except ImportError:  # pragma: no cover
    pass

try:
    from pyNastran.converters.aflr.surf.surf_io import SurfIO
    CLASS_MAP['surf'] = SurfIO
except ImportError:  # pragma: no cover
    pass

try:
    from pyNastran.converters.aflr.ugrid.ugrid_io import UGRID_IO
    CLASS_MAP['ugrid'] = UGRID_IO
    CLASS_MAP['ugrid3d'] = UGRID_IO
except ImportError:  # pragma: no cover
    pass

try:
    from pyNastran.converters.abaqus.abaqus_io import AbaqusIO
    CLASS_MAP['abaqus'] = AbaqusIO
except ImportError:  # pragma: no cover
    pass

try:
    from pyNastran.converters.openfoam.openfoam_io import OpenFoamIO
    CLASS_MAP['openfoam_hex'] = OpenFoamIO
    CLASS_MAP['openfoam_shell'] = OpenFoamIO
    CLASS_MAP['openfoam_faces'] = OpenFoamIO
except ImportError:  # pragma: no cover
    pass
