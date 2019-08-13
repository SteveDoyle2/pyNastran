from pyNastran.converters.cart3d.cart3d import Cart3D, read_cart3d
from pyNastran.converters.tecplot.tecplot import Tecplot, Zone

def cart3d_to_tecplot(cart3d_filename, tecplot_filename, log=None, debug=False):
    """
    Converts Cart3d to Tecplot
    """
    if isinstance(cart3d_filename, Cart3D):
        model = cart3d_filename
    else:
        model = read_cart3d(cart3d_filename, log=log, debug=debug)

    tecplot = Tecplot()
    tecplot.log = model.log
    zone = Zone(model.log)
    zone.headers_dict['VARIABLES'] = ['X', 'Y', 'Z']
    zone.xyz = model.points
    zone.tri_elements = model.elements + 1
    tecplot.zones = [zone]

    tecplot.write_tecplot(tecplot_filename, adjust_nids=False)
    return tecplot
