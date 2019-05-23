"""
defines:
 - cart3d = nastran_to_cart3d(bdf, log=None, debug=False)
 - nastran_to_cart3d_filename(bdf_filename, cart3d_filename,
                              log=None, debug=False)

"""
from numpy import zeros, arange, array, array_equal
import numpy as np

from pyNastran.bdf.bdf import BDF
from pyNastran.converters.cart3d.cart3d import Cart3D

LINE_ELEMENTS = ['CBAR', 'CBEAM', 'CROD', 'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4']

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
    nelements = 0
    if 'CTRIA3' in bdf.card_count:
        nelements += bdf.card_count['CTRIA3']
    if 'CQUAD4' in bdf.card_count:
        nelements += bdf.card_count['CQUAD4'] * 2

    nodes = zeros((nnodes, 3), 'float64')
    elements = zeros((nelements, 3), 'int32')
    regions = zeros(nelements, 'int32')

    nids = array(list(bdf.nodes.keys()), dtype='int32')
    nids_expected = arange(1, len(nids) + 1)

    if array_equal(nids, nids_expected):
        _store_sequential_nodes(bdf, nodes, elements, regions)
    else:
        #print('not equal...')
        #i = 0
        nid_map = {}

        for node_id, node in sorted(bdf.nodes.items()):
            i = np.where(nids == node_id)[0][0]
            nodes[i, :] = node.get_position()
            nid_map[node_id] = i + 1
        #print('nid_map =', nid_map)

        i = 0
        for unused_element_id, element in sorted(bdf.elements.items()):
            if element.type == 'CTRIA3':
                nids = element.node_ids
                elements[i, :] = [nid_map[nid] for nid in nids]
                regions[i] = element.material_ids[0]
            elif element.type == 'CQUAD4':
                nids = element.node_ids
                quad = [nid_map[nid] for nid in nids]
                mid = element.material_ids[0]

                # TODO: splits on edge 1-3, not the max angle
                #       since we're just comparing the models, it doesn't matter
                elements[i, :] = [quad[0], quad[1], quad[2]]
                regions[i] = mid
                i += 1
                elements[i, :] = [quad[0], quad[2], quad[3]]
                regions[i] = mid
            elif element.type in LINE_ELEMENTS:
                continue
            else:
                raise NotImplementedError(element.type)
            i += 1

    assert elements.min() > 0, elements
    cart3d.nodes = nodes
    cart3d.elements = elements - 1
    cart3d.regions = regions
    return cart3d

def _store_sequential_nodes(bdf, nodes, elements, regions):
    # we don't need to renumber the nodes
    # so we don't need to make an nid_map
    i = 0
    for node_id, node in sorted(bdf.nodes.items()):
        nodes[i, :] = node.get_position()
        i += 1

    j = 0
    for unused_element_id, element in sorted(bdf.elements.items()):
        if element.type == 'CTRIA3':
            nids = element.node_ids
            elements[j, :] = nids
            regions[j] = element.Mid()
        elif element.type in LINE_ELEMENTS:
            pass
        #elif element.type == 'CQUAD4':
            #nids = element.node_ids
            #elements[i, :] = nids
            #regions[i] = element.Mid()
        else:
            raise NotImplementedError(element.type)
        j += 1

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

    with open(cart3d_filename, 'w', encoding='utf8') as cart3d:
        cart3d.write('%s %s\n' % (nnodes, nelements))
        node_id_shift = {}
        i = 1
        for node_id, node in sorted(model.nodes.items()):
            node_id_shift[node_id] = i
            x, y, z = node.get_position()
            cart3d.write('%s %s %s\n' % (x, y, z))
            i += 1
        mids = ''
        j = 0
        for unused_element_id, element in sorted(model.elements.items()):
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
            n1 = node_id_shift[n1]
            n2 = node_id_shift[n2]
            n3 = node_id_shift[n3]
            mid = element.Mid()
            cart3d.write('%i %i %i\n' % (n1, n2, n3))
            mids += '%i ' % mid
            if j != 0 and j % 20 == 0:
                mids += '\n'
            j += 1
        cart3d.write(mids + '\n')
