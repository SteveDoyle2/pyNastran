from pyNastran.converters.cart3d.cart3d import Cart3D
from six.moves import zip
from itertools import count
from numpy import unique

def cart3d_to_usm3d_bc(cart3d, usm3d_bc_filename):
    """
    Creates a USM3D boundary condtion file from the Cart3d regions.

    Parameters
    ----------
    cart3d : Cart3D()
        a Cart3D object
    usm3d_bc_filename : str
        path to the output BC file
    """
    # nodes = cart3d.points
    elements = cart3d.elements
    regions = cart3d.regions
    #nodes = cart3d.nodes
    #elements = cart3d.elements
    #regions = cart3d.regions

    with open(usm3d_bc_filename, 'wb') as usm3d_bc:
        patches = unique(regions)
        npatches = len(patches)
        nelements = elements.shape[0]

        usm3d_bc.write('%-8s %-8s %-8s %s\n'  % (nelements, 'intA', npatches, 'intB'))
        usm3d_bc.write('Triangle   Patch            Nodes\n')
        for i, element, iregion in zip(count(), elements, regions):
            (n1, n2, n3) = element
            usm3d_bc.write('%-8s %-8s %-8s %-8s %s\n' % (i+1, iregion, n1, n2, n3))

def cart3d_to_usm3d_bc_filename(cart3d_filename, usm3d_bc_filename, log=None, debug=False):
    """
    Creates a USM3D boundary condtion file from the Cart3d regions.

    Parameters
    ----------
    cart3d_filename : str
        path to the input Cart3D file
    usm3d_bc_filename : str
        path to the output BC file
    log : log(); default=None
        a logger object
    debug : bool; default=False
        True/False (used if log is not defined)
    """
    cart3d = Cart3D(log=log, debug=debug)
    cart3d.read_cart3d(cart3d_filename)
    cart3d_to_usm3d_bc(cart3d, usm3d_bc_filename)


def main():  # pragma: no cover
    cart3d_filename = 'threePlugs.a.tri'
    usm3d_bc_filename = 'threePlugs.bc'
    cart3d_to_usm3d_bc_filename(cart3d_filename, usm3d_bc_filename)

if __name__ == '__main__':  # pragma: no cover
    main()
