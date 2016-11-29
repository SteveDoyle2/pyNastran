"""
Defines:
 - tecplot_to_cart3d_filename(tecplot_filename, cart3d_filename, debug=True)
 - tecplot_to_cart3d(tecplot_filename, cart3d_filename=None, debug=True)
"""

import numpy as np
#from pyNastran.bdf.bdf_interface.dev_utils import remove_unassociated_nodes
from pyNastran.converters.tecplot.tecplot import read_tecplot
from pyNastran.converters.cart3d.cart3d import Cart3D

def tecplot_to_cart3d_filename(tecplot_filename, cart3d_filename, debug=True):
    """Converts a Tecplot file to Cart3d."""
    return tecplot_to_cart3d(tecplot_filename, cart3d_filename, debug=debug)


def tecplot_to_cart3d(tecplot_filename, cart3d_filename=None, debug=True):
    """Converts a Tecplot file to Cart3d."""
    if isinstance(tecplot_filename, str):
        model = read_tecplot(tecplot_filename, debug=debug)
    else:
        model = tecplot_filename

    headers_no_xyz = model.variables[3:] # drop the xyz, get what's left
    iCp = headers_no_xyz.index('cp')
    cp = model.results[:, iCp]
    tris = model.quad_elements[:, :3]

    npoints = model.xyz.shape[0]
    nelements = tris.shape[0]
    assert tris.shape[1] == 3, tris.shape
    print('npoints=%s nelements=%s' % (npoints, nelements))

    #removed_nodes = False
    assert len(cp) == npoints
    regions = np.zeros(nelements, dtype='int32')
    ones_float = np.ones(npoints, dtype='float64')

    cart3d_model = Cart3D()
    cart3d_model.points = model.xyz
    cart3d_model.regions = regions
    cart3d_model.elements = tris + 1

    cart3d_model.loads = {
        'Cp' : cp, # nodal Cp
        'rho' : ones_float,
        'rhoU' : ones_float,
        'rhoV' : ones_float,
        'rhoW' : ones_float,
        'E' : ones_float,
    }
    if cart3d_filename is not None:
        cart3d_model.write_cart3d(cart3d_filename)
    return cart3d_model


#tecplot_filename2 = r'C:\NASA\m4\formats\git\pyNastran\pyNastran\converters\cart3d\PressureMapping\70degSpoiler_6degAlpha_tec_boundary.dat'
tecplot_filename2 = r'C:\NASA\m4\formats\git\pyNastran\pyNastran\converters\cart3d\PressureMapping\point_fmt.dat'
cart3d_filename2 = 'wing.tri'
tecplot_to_cart3d_filename(tecplot_filename2, cart3d_filename2, debug=True)
