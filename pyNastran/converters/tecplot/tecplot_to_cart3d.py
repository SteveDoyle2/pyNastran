"""
Defines:
 - tecplot_to_cart3d_filename(tecplot_filename, cart3d_filename, debug=True)
 - tecplot_to_cart3d(tecplot_filename, cart3d_filename=None, debug=True)
"""

import numpy as np
#from pyNastran.bdf.mesh_utils.remove_unused import remove_unassociated_nodes
from pyNastran.converters.tecplot.tecplot import read_tecplot
from pyNastran.converters.cart3d.cart3d import Cart3D


def tecplot_to_cart3d_filename(tecplot_filename, cart3d_filename, debug=True):
    """Converts a Tecplot file to Cart3d."""
    return tecplot_to_cart3d(tecplot_filename, cart3d_filename, debug=debug)


def tecplot_to_cart3d(tecplot_filename, cart3d_filename=None, debug=True):
    """
    Converts a Tecplot file to Cart3d.

    It's assumed that quads are actually degenerate triangles
    """
    if isinstance(tecplot_filename, str):
        model = read_tecplot(tecplot_filename, debug=debug)
    else:
        model = tecplot_filename

    if len(model.tri_elements) and len(model.quad_elements):
        tris = np.vstack([
            model.tri_elements,
            model.quad_elements[:, :3],
        ])
    elif len(model.tri_elements):
        tris = model.tri_elements
    elif len(model.quad_elements):
        tris = model.quad_elements[:, :3]
    else:
        raise NotImplementedError('need quads/tris')

    npoints = model.xyz.shape[0]
    nelements = tris.shape[0]
    assert tris.shape[1] == 3, tris.shape
    #print('npoints=%s nelements=%s' % (npoints, nelements))

    #removed_nodes = False
    regions = np.zeros(nelements, dtype='int32')
    ones_float = np.ones(npoints, dtype='float64')

    cart3d_model = Cart3D()
    cart3d_model.points = model.xyz
    cart3d_model.regions = regions
    cart3d_model.elements = tris + 1

    headers_no_xyz = model.variables[3:] # drop the xyz, get what's left
    if 'cp' in headers_no_xyz:
        iCp = headers_no_xyz.index('cp')
        cp = model.results[:, iCp]
        assert len(cp) == npoints
        cart3d_model.loads = {
            'Cp' : cp, # nodal Cp
            'rho' : ones_float,
            'rhoU' : ones_float,
            'rhoV' : ones_float,
            'rhoW' : ones_float,
            'E' : ones_float,
        }
    else:
        model.log.debug('skipping Cp')

    if cart3d_filename is not None:
        cart3d_model.write_cart3d(cart3d_filename)
    return cart3d_model

def main():
    """runs the test problem"""
    tecplot_filename2 = r'PressureMapping\point_fmt.dat'
    cart3d_filename2 = 'wing.tri'
    tecplot_to_cart3d_filename(tecplot_filename2, cart3d_filename2, debug=True)

if __name__ == '__main__':
    main()
