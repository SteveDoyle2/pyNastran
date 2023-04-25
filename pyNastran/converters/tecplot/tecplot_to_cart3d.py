"""
Defines:
 - tecplot_to_cart3d_filename(tecplot_filename, cart3d_filename, debug=True)
 - tecplot_to_cart3d(tecplot_filename, cart3d_filename=None, debug=True)
"""
from typing import Optional
import numpy as np
#from pyNastran.bdf.mesh_utils.remove_unused import remove_unassociated_nodes
from cpylog import SimpleLogger
from pyNastran.converters.tecplot.tecplot import read_tecplot, Zone, Tecplot
from pyNastran.converters.cart3d.cart3d import Cart3D


def tecplot_to_cart3d_filename(tecplot_filename: str, cart3d_filename: str,
                               remove_degenerate_tris: bool=True,
                               log: Optional[SimpleLogger]=None, debug:
                               bool=True) -> Cart3D:
    """Converts a Tecplot file to Cart3d."""
    return tecplot_to_cart3d(tecplot_filename, cart3d_filename,
                             remove_degenerate_tris=remove_degenerate_tris,
                             log=log, debug=debug)


def tecplot_to_cart3d(tecplot_filename: str | Tecplot,
                      cart3d_filename: Optional[str]=None,
                      remove_degenerate_tris: bool=True,
                      log=None, debug: bool=True) -> Cart3D:
    """
    Converts a Tecplot file to Cart3d.

    Parameters
    ----------
    remove_degenerate_tris : bool; default=False
        removes degenerate triangles (triangles with an area of 0.0)
    """
    if isinstance(tecplot_filename, Tecplot):
        model = tecplot_filename
    else:
        model = read_tecplot(tecplot_filename, log=log, debug=debug)

    assert len(model.zones) >= 1, 'no zones found'
    nnodes_total = 0
    xyz_list = []
    regions_list = []
    ones_float_list = []
    elements_list = []
    for izone, zone in enumerate(model.zones):
        #print(zone)
        assert len(zone.hexa_elements) == 0, zone
        assert len(zone.tet_elements) == 0, zone
        xyzi = zone.get_xyz()
        tri = get_zone_tris(zone.tri_elements, zone.quad_elements)

        assert tri.shape[0] > 0, tri.shape
        assert xyzi.shape[0] > 0, xyzi.shape
        if remove_degenerate_tris:
            area = _get_tri_area(tri, xyzi)
            iarea = np.where(area > 0.)[0]
            tri = tri[iarea, :]

        npoints = xyzi.shape[0]
        assert npoints > 0, xyzi.shape
        nelements = tri.shape[0]
        assert nelements > 0, nelements
        assert tri.shape[1] == 3, tri.shape
        #print('npoints=%s nelements=%s' % (npoints, nelements))

        #removed_nodes = False
        regionsi = np.ones(nelements, dtype='int32') * (izone + 1)
        ones_floati = np.ones(npoints, dtype='float64')
        regions_list.append(regionsi)
        ones_float_list.append(ones_floati)
        xyz_list.append(xyzi)
        elements_list.append(nnodes_total + tri)
        nnodes_total += npoints
    regions = np.hstack(regions_list)
    ones_float = np.hstack(ones_float_list)
    xyz = np.vstack(xyz_list)
    elements = np.vstack(elements_list)
    del regions_list, ones_float_list, elements_list
    nnodes = len(xyz)

    cart3d_model = Cart3D(log=log)
    cart3d_model.points = xyz
    cart3d_model.regions = regions
    cart3d_model.elements = elements # + 1
    assert nnodes > 0, nnodes

    nxyz = model.zones[0].nxyz
    headers_no_xyz = model.variables[nxyz:] # drop the xyz, get what's left
    nodal_results = model.stack_results()
    assert len(nodal_results) == nnodes, nnodes
    if 'cp' in headers_no_xyz:
        iCp = headers_no_xyz.index('cp')
        cp = nodal_results[:, iCp]
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

def _get_tri_area(tris, xyz: np.ndarray) -> np.ndarray:
    n1 = xyz[tris[:, 0], :]
    n2 = xyz[tris[:, 1], :]
    n3 = xyz[tris[:, 2], :]
    normal = np.cross(n2 - n1, n3 - n1)
    area = np.linalg.norm(normal, axis=1)
    return area

def get_zone_tris(tri_elements: np.ndarray,
                  quad_elements: np.ndarray) -> np.ndarray:
    """gets the tris and associated xyz points"""
    ntris = len(tri_elements)
    nquads = len(quad_elements)
    if ntris and nquads:
        # double stack the quads to size the array
        # then overwrite the second set of quads
        tris = np.vstack([
            zone.tri_elements,
            zone.quad_elements[:, :3],
            zone.quad_elements[:, :3],
        ])
        tris[ntris+nquads:, [0, 1]] = quad_elements[:, 3:]
        tris[ntris+nquads:, 2] = quad_elements[:, 0]
    elif ntris:
        tris = tri_elements
    elif nquads:
        # double stack the quads to size the array
        # then overwrite the second set of quads
        tris = np.vstack([
            quad_elements[:, :3],
            quad_elements[:, :3],
        ])
        tris[ntris+nquads:, [0, 1]] = quad_elements[:, 3:]
        tris[ntris+nquads:, 2] = quad_elements[:, 0]
    else:
        raise NotImplementedError('need quads/tris')

    return tris

def cmd_line_tecplot_to_cart3d():  # pragma: no cover
    """runs the test problem"""
    import sys
    tecplot_filename2 = sys.argv[1]
    cart3d_filename2 = sys.argv[2]
    tecplot_to_cart3d_filename(tecplot_filename2, cart3d_filename2, debug=True)

if __name__ == '__main__':  # pragma: no cover
    cmd_line_tecplot_to_cart3d()
