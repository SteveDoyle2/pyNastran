"""import various codes with backup for failed imports"""
CLASS_MAP = {}
from pyNastran import DEV

try:
    from pyNastran.converters.dev.avus.avus_io import AvusIO
    CLASS_MAP['avus'] = AvusIO
except ModuleNotFoundError:  # pragma: no cover
    pass

try:
    from pyNastran.converters.cart3d.cart3d_io import Cart3dIO
    CLASS_MAP['cart3d'] = Cart3dIO
except ModuleNotFoundError:  # pragma: no cover
    pass

try:
    from pyNastran.converters.panair.panair_io import PanairIO
    CLASS_MAP['panair'] = PanairIO
except ModuleNotFoundError:  # pragma: no cover
    pass

try:
    from pyNastran.converters.stl.stl_io import STL_IO
    CLASS_MAP['stl'] = STL_IO
except ModuleNotFoundError:  # pragma: no cover
    pass

try:
    from pyNastran.converters.usm3d.usm3d_io import Usm3dIO
    CLASS_MAP['usm3d'] = Usm3dIO
except ModuleNotFoundError:  # pragma: no cover
    pass

try:
    from pyNastran.converters.tecplot.tecplot_io import TecplotIO
    CLASS_MAP['tecplot'] = TecplotIO
except ModuleNotFoundError:  # pragma: no cover
    pass

try:
    from pyNastran.converters.lawgs.wgs_io import LaWGS_IO
    CLASS_MAP['lawgs'] = LaWGS_IO
except ModuleNotFoundError:  # pragma: no cover
    pass

try:
    from pyNastran.converters.tetgen.tetgen_io import TetgenIO
    CLASS_MAP['tetgen'] = TetgenIO
except ModuleNotFoundError:  # pragma: no cover
    pass

try:
    from pyNastran.converters.shabp.shabp_io import ShabpIO
    CLASS_MAP['shabp'] = ShabpIO
except ModuleNotFoundError:  # pragma: no cover
    pass

try:
    from pyNastran.converters.su2.su2_io import SU2_IO
    CLASS_MAP['su2'] = SU2_IO
except ModuleNotFoundError:  # pragma: no cover
    pass

try:
    from pyNastran.converters.fast.fast_io import FastIO
    CLASS_MAP['fast'] = FastIO
except ModuleNotFoundError:  # pragma: no cover
    pass

try:
    from pyNastran.converters.avl.avl_io import AVL_IO
    CLASS_MAP['avl'] = AVL_IO
except ModuleNotFoundError:  # pragma: no cover
    raise


from pyNastran.gui.gui_objects.vtk_io import VtkIO
CLASS_MAP['vtk'] = VtkIO

from pyNastran.converters.nastran.gui.nastran_io import NastranIO
#CLASS_MAP['nastran'] = NastranIO


#from pyNastran.converters.dev.plot3d.plot3d_io import Plot3d_io

try:
    from pyNastran.converters.aflr.aflr2.bedge_io import BEdge_IO
    CLASS_MAP['bedge'] = BEdge_IO
except ModuleNotFoundError:  # pragma: no cover
    pass

try:
    from pyNastran.converters.aflr.surf.surf_io import SurfIO
    CLASS_MAP['surf'] = SurfIO
except ModuleNotFoundError:  # pragma: no cover
    pass

try:
    from pyNastran.converters.aflr.ugrid.ugrid_io import UGRID_IO
    CLASS_MAP['ugrid'] = UGRID_IO
    CLASS_MAP['ugrid3d'] = UGRID_IO
except ModuleNotFoundError:  # pragma: no cover
    pass

try:
    from pyNastran.converters.abaqus.abaqus_io import AbaqusIO
    CLASS_MAP['abaqus'] = AbaqusIO
except ModuleNotFoundError:  # pragma: no cover
    pass
