from pyNastran.converters.cart3d.cart3d_reader import Cart3DReader
from itertools import izip, count
from numpy import unique

def cart3d_to_usm3d_bc(cart3d, log=None, debug=False):
    """
    Converts a Cart3DReader object to STL format.

    :param cart3d: a Cart3DReader object
    :param log:    a logger object (or None)
    :param debug:  True/False (used if log is not defined)

    :returns stl: an STLReader object
    """
    normals = cart3d.normals()
    asdf
    #stl = STLReader(log=log, debug=debug)
    #stl.nodes = nodes
    #stl.elements = elements
    #stl.write_stl(stl_filename)
    #return stl

def cart3d_to__usm3d_bc_filename(cart3d_filename, usm3d_bc_filename, log=None, debug=False):
    """
    Converts a Cart3D file to STL format.

    :param cart3d_filename: path to the input Cart3D file
    :param usm3d_bc_filename:    path to the output BC file
    :param log:             a logger object (or None)
    :param debug:           True/False (used if log is not defined)
    """
    cart3d = Cart3DReader(log=log, debug=debug)
    (nodes, elements, regions, loads) = cart3d.read_cart3d(cart3d_filename)

    #nodes = cart3d.nodes
    #elements = cart3d.elements
    #regions = cart3d.regions

    usm3d_bc = open(usm3d_bc_filename, 'wb')
    patches = unique(regions)
    npatches = len(patches)
    nelements, three = elements.shape

    usm3d_bc.write('%-8s %-8s %-8s %s\n'  % (nelements, 'intA', npatches, 'intB'))
    usm3d_bc.write('Triangle   Patch            Nodes\n')
    for i, element, iregion in izip(count(), elements, regions):
        (n1, n2, n3) = element
        usm3d_bc.write('%-8s %-8s %-8s %-8s %s\n' % (i+1, iregion, n1, n2, n3))
    usm3d_bc.close()

if __name__ == '__main__':
    cart3d_filename = 'threePlugs.a.tri'
    usm3d_bc_filename = 'threePlugs.bc'
    cart3d_to__usm3d_bc_filename(cart3d_filename, usm3d_bc_filename)