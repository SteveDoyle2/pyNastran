from six import iteritems
from numpy import zeros, ones

from pyNastran.bdf.bdf import BDF
from pyNastran.converters.cart3d.cart3d import Cart3D


def nastran_to_cart3d(bdf, log=None, debug=False):
    """
    Converts a Nastran BDF() object to a Cart3D() object.

    Parameters
    ----------
    bdf : BDF()
        a BDF object
    log : log; default=None -> dummyLogger
        a logger object
    debug : bool; default=False
        True/False (used if log is not defined)

    Returns
    -------
    cart3d : Cart3D()
        a Cart3D object
    """
    cart3d = Cart3D(log=log, debug=debug)

    nnodes = len(bdf.nodes)
    nelements = len(bdf.elements)

    nodes = zeros((nnodes, 3), 'float64')
    elements = zeros((nelements, 3), 'int32')
    regions = zeros(nelements, 'int32')

    i = 0
    for node_id, node in sorted(iteritems(bdf.nodes)):
        nodes[i, :] = node.get_position()
        i += 1

    i = 0
    for element_id, element in sorted(iteritems(bdf.elements)):
        if element.type == 'CTRIA3':
            nids = element.node_ids
            elements[i, :] = nids
            regions[i] = element.Mid()
        else:
            raise NotImplementedError(element.type)
        i += 1


    assert elements.min() > 0, elements
    cart3d.points = nodes
    cart3d.elements = elements
    cart3d.regions = regions
    return cart3d


def nastran_to_cart3d_filename(bdf_filename, cart3d_filename, log=None, debug=False):
    """
    Creates a Nastran BDF from a Cart3D file.

    Parameters
    ----------
    bdf_filename : str
        the path to the bdf file
    cart3d_filename : str
        the path to the cart3d output file
    log : log; default=None -> dummyLogger
        a logger object
    debug : bool; default=False
        True/False (used if log is not defined)
    """
    model = BDF(log=log, debug=debug)
    model.read_bdf(bdf_filename)
    nnodes = len(model.nodes)
    nelements = len(model.elements)

    f = open(cart3d_filename, 'wb')
    f.write('%s %s\n' % (nnodes, nelements))
    node_id_shift = {}
    i = 1
    for node_id, node in sorted(iteritems(model.nodes)):
        node_id_shift[node_id] = i
        x, y, z = node.get_position()
        f.write('%s %s %s\n' % (x, y, z))
        i += 1
    mids = ''
    j = 0
    for element_id, element in sorted(iteritems(model.elements)):
        if element.type in ['CQUADR', 'CQUAD4', 'CONM2']:
            print('element type=%s is not supported' % element.type)
            continue
        assert element.type in ['CTRIA3', 'CTRIAR'], element.type


        out = element.node_ids
        try:
            n1, n2, n3 = out
        except:
            print("type =", element.type)
            raise
        #print out
        n1 = node_id_shift[n1]
        n2 = node_id_shift[n2]
        n3 = node_id_shift[n3]
        mid = element.Mid()
        f.write('%s %s %s\n' % (n1, n2, n3))
        mids += '%s ' % mid
        if j != 0 and j % 20 == 0:
            mids += '\n'
        j += 1
    f.write(mids + '\n')
    f.close()

if __name__ == '__main__':  # pragma: no cover
    bdf_filename = 'g278.bdf'
    cart3d_filename = 'g278.tri'
    nastran_to_cart3d_filename(bdf_filename, cart3d_filename)
