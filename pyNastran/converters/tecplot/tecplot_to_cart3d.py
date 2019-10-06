"""
Defines:
 - tecplot_to_cart3d_filename(tecplot_filename, cart3d_filename, debug=True)
 - tecplot_to_cart3d(tecplot_filename, cart3d_filename=None, debug=True)
"""

import numpy as np
#from pyNastran.bdf.mesh_utils.remove_unused import remove_unassociated_nodes
from pyNastran.converters.tecplot.tecplot import read_tecplot
from pyNastran.converters.cart3d.cart3d import Cart3D


def tecplot_to_cart3d_filename(tecplot_filename, cart3d_filename,
                               remove_degenerate_tris=True, log=None, debug=True):
    """Converts a Tecplot file to Cart3d."""
    return tecplot_to_cart3d(tecplot_filename, cart3d_filename,
                             remove_degenerate_tris=remove_degenerate_tris, log=log, debug=debug)


def tecplot_to_cart3d(tecplot_filename, cart3d_filename=None,
                      remove_degenerate_tris=True, log=None, debug=True):
    """
    Converts a Tecplot file to Cart3d.

    Parameters
    ----------
    remove_degenerate_tris : bool; default=False
        removes degenerate triangles (triangles with an area of 0.0)
    """
    if isinstance(tecplot_filename, str):
        model = read_tecplot(tecplot_filename, log=log, debug=debug)
    else:
        model = tecplot_filename

    assert len(model.zones) == 1, 'only 1 zone is supported'
    for izone, zone in enumerate(model.zones):
        #print(zone)
        tris, xyz = get_zone_tris_xyz(zone)
        if remove_degenerate_tris:
            assert tris.shape[0] > 0, tris.shape
            assert xyz.shape[0] > 0, xyz.shape
            A = _get_tri_area(tris, xyz)
            iarea = np.where(A > 0.)[0]
            tris = tris[iarea, :]

        npoints = xyz.shape[0]
        assert npoints > 0, xyz.shape
        nelements = tris.shape[0]
        assert nelements > 0, nelements
        assert tris.shape[1] == 3, tris.shape
        #print('npoints=%s nelements=%s' % (npoints, nelements))

        #removed_nodes = False
        regions = np.zeros(nelements, dtype='int32')
        ones_float = np.ones(npoints, dtype='float64')

    cart3d_model = Cart3D(log=log)
    cart3d_model.points = xyz
    cart3d_model.regions = regions
    cart3d_model.elements = tris # + 1
    assert npoints > 0, npoints

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

def _get_tri_area(tris, xyz):
    n1 = xyz[tris[:, 0], :]
    n2 = xyz[tris[:, 1], :]
    n3 = xyz[tris[:, 2], :]
    normal = np.cross(n2 - n1, n3 - n1)
    A = np.linalg.norm(normal, axis=1)
    return A

def get_zone_tris_xyz(zone):
    """gets the tris and associated xyz points"""
    ntris = len(zone.tri_elements)
    nquads = len(zone.quad_elements)
    if ntris and nquads:
        # double stack the quads to size the array
        # then overwrite the second set of quads
        tris = np.vstack([
            zone.tri_elements,
            zone.quad_elements[:, :3],
            zone.quad_elements[:, :3],
        ])
        tris[ntris+nquads:, [0, 1]] = zone.quad_elements[:, 3:]
        tris[ntris+nquads:, 2] = zone.quad_elements[:, 0]
    elif ntris:
        tris = zone.tri_elements
    elif nquads:
        # double stack the quads to size the array
        # then overwrite the second set of quads
        tris = np.vstack([
            zone.quad_elements[:, :3],
            zone.quad_elements[:, :3],
        ])
        tris[ntris+nquads:, [0, 1]] = zone.quad_elements[:, 3:]
        tris[ntris+nquads:, 2] = zone.quad_elements[:, 0]
    else:
        raise NotImplementedError('need quads/tris')

    xyz = zone.get_xyz()
    return tris, xyz

def main():  # pragma: no cover
    """runs the test problem"""
    tecplot_filename2 = r'PressureMapping\point_fmt.dat'
    cart3d_filename2 = 'wing.tri'
    tecplot_to_cart3d_filename(tecplot_filename2, cart3d_filename2, debug=True)

if __name__ == '__main__':  # pragma: no cover
    main()
