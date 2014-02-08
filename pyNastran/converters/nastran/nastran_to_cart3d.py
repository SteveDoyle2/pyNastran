from numpy import zeros, ones

from pyNastran.bdf.bdf import BDF
from pyNastran.converters.cart3d.cart3d_reader import Cart3DReader


def nastran_to_cart3d(bdf, log=None, debug=False):
    """
    Converts a Nastran BDF object to Cart3D format.
    
    :param bdf:    a BDF object
    :param log:    a logger object (or None)
    :param debug:  True/False (used if log is not defined)

    :returns cart3d: a Cart3D object
    """
    cart3d = Cart3DReader(log=log, debug=debug)
    
    nnodes = len(model.nodes)
    nelements = len(model.elements)
    
    nodes = zeros((nnodes, 3), 'float64')
    elements = zeros((nelements, 3), 'int32')
    regions = zeros(nelements, 'int32')

    i = 0
    for node_id, node in sorted(model.nodes.iteritems()):
        elements[i, :] = node.Position()
    for element_id, element in sorted(model.elements.iteritems()):
        if element.type == 'CTRIA3':
            elements[i, :] = element.NodeIDs()
            regions[i] = element.Mid()
        else:
            raise NotImplementedError(element.type)

    cart3d.nodes = nodes
    cart3d.elements = elements
    cart3d.regions = regions
    return cart3d


def nastran_to_cart3d_filename(bdf_filename, cart3d_filename, log=None, debug=False):
    """
    Converts a Nastran file to Cart3D format.
    
    :param bdf_filename: the path to the BDF
    :param cart3d_filename: the path to the Cart3D output file
    :param log:    a logger object (or None)
    :param debug:  True/False (used if log is not defined)
    """
    model = BDF(log=log, debug=debug)
    model.read_bdf(bdf_filename)
    nnodes = len(model.nodes)
    nelements = len(model.elements)

    f = open(cart3d_filename, 'wb')
    f.write('%s %s\n' % (nnodes, nelements))
    node_id_shift = {}
    i = 1
    for node_id, node in sorted(model.nodes.iteritems()):
        node_id_shift[node_id] = i
        x, y, z = node.Position()
        f.write('%s %s %s\n' % (x, y, z))
        i += 1
    mids = ''
    j = 0
    for element_id, element in sorted(model.elements.iteritems()):
        if element.type in ['CQUADR', 'CONM2']:
            continue
        assert element.type in ['CTRIA3', 'CTRIAR'], element.type


        out = element.nodeIDs()
        try:
            n1, n2, n3 = out
        except:
            print "type =", element.type
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

if __name__ == '__main__':
    bdf_filename = 'g278.bdf'
    cart3d_filename = 'g278.tri'
    nastran_to_cart3d_filename(bdf_filename, cart3d_filename)