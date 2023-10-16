"""
defines:
 - cart3d = nastran_to_cart3d(bdf, log=None, debug=False)
 - nastran_to_cart3d_filename(bdf_filename, cart3d_filename,
                              log=None, debug=False)

"""
import numpy as np

from pyNastran.converters.cart3d.cart3d import Cart3D

LINE_ELEMENTS = ['CBAR', 'CBEAM', 'CROD', 'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4']
from pyNastran.bdf.bdf import BDF


def nastran_to_cart3d(bdf: BDF, log=None, debug: bool=False) -> Cart3D:
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
    nnodes, nelements = nnodes_nelements(bdf)
    nodes = np.zeros((nnodes, 3), 'float64')
    elements = np.zeros((nelements, 3), 'int32')
    regions = np.zeros(nelements, 'int32')

    nids = np.array(list(bdf.nodes.keys()), dtype='int32')
    nids_expected = np.arange(1, len(nids) + 1)

    if np.array_equal(nids, nids_expected):
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
            if element.type in {'CTRIA3', 'CTRIAR'}:
                nids = element.node_ids
                elements[i, :] = [nid_map[nid] for nid in nids]
                regions[i] = element.pid
            elif element.type in {'CQUAD4', 'CQUADR'}:
                nids = element.node_ids
                quad = [nid_map[nid] for nid in nids]
                pid = element.pid

                # TODO: splits on edge 1-3, not the max angle
                #       since we're just comparing the models, it doesn't matter
                elements[i, :] = [quad[0], quad[1], quad[2]]
                regions[i] = pid
                i += 1
                elements[i, :] = [quad[0], quad[2], quad[3]]
                regions[i] = pid
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

def _store_sequential_nodes(bdf: BDF, nodes, elements, regions) -> None:
    # we don't need to renumber the nodes
    # so we don't need to make an nid_map
    i = 0
    for node_id, node in sorted(bdf.nodes.items()):
        nodes[i, :] = node.get_position()
        i += 1

    j = 0
    for unused_element_id, element in sorted(bdf.elements.items()):
        if element.type in {'CTRIA3', 'CTRIAR'}:
            nids = element.node_ids
            elements[j, :] = nids
            regions[j] = element.Pid()
        elif element.type in LINE_ELEMENTS:
            pass
        elif element.type in {'CQUAD4', 'CQUADR'}:
            nids = element.node_ids
            pid = element.Pid()
            elements[j, :] = nids[:-1]
            regions[j] = pid

            elements[j+1, :] = [nids[0], nids[2], nids[3]]
            regions[j+1] = pid
            j += 1
        else:
            raise NotImplementedError(element.type)
        j += 1

def nastran_to_cart3d_filename(bdf_filename: str,
                               cart3d_filename: str,
                               log=None, debug: bool=False) -> None:
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
    nnodes, nelements = nnodes_nelements(model)

    log = model.log
    with open(cart3d_filename, 'w', encoding='utf8') as cart3d:
        cart3d.write('%s %s\n' % (nnodes, nelements))
        node_id_shift = {}
        i = 1
        for node_id, node in sorted(model.nodes.items()):
            node_id_shift[node_id] = i
            x, y, z = node.get_position()
            cart3d.write('%s %s %s\n' % (x, y, z))
            i += 1
        pids = ''
        j = 0
        for unused_element_id, element in sorted(model.elements.items()):
            if element.type in ['CONM2']:
                continue

            card_type = element.type
            if card_type in {'CQUADR', 'CQUAD4'}:
                out = element.node_ids
                #print('element type=%s is not supported' % element.type)
                n1, n2, n3, n4 = out

                n1 = node_id_shift[n1]
                n2 = node_id_shift[n2]
                n3 = node_id_shift[n3]
                n4 = node_id_shift[n4]
                pid = element.pid
                cart3d.write('%i %i %i\n' % (n1, n2, n3))
                cart3d.write('%i %i %i\n' % (n1, n3, n4))
                pids += '%i ' % pid
                if j != 0 and j % 20 == 0:
                    pids += '\n'
                j += 1

                pids += '%i ' % pid
                if j != 0 and j % 20 == 0:
                    pids += '\n'
                j += 1

            elif card_type in {'CTRIA3', 'CTRIAR'}:
                out = element.node_ids
                n1, n2, n3 = out

                n1 = node_id_shift[n1]
                n2 = node_id_shift[n2]
                n3 = node_id_shift[n3]
                pid = element.pid
                cart3d.write('%i %i %i\n' % (n1, n2, n3))
                pids += '%i ' % pid
                if j != 0 and j % 20 == 0:
                    pids += '\n'
                j += 1
            elif element.type in LINE_ELEMENTS:
                continue
            else:  # pragma: no cover
                log.error(f'card_type={card_type} is not supported')
                raise NotImplementedError(element.type)
                #continue

        cart3d.write(pids + '\n')

def nnodes_nelements(bdf: BDF) -> tuple[int, int]:
    nnodes = len(bdf.nodes)
    nelements = 0
    if 'CTRIA3' in bdf.card_count:
        nelements += bdf.card_count['CTRIA3']
    if 'CTRIAR' in bdf.card_count:
        nelements += bdf.card_count['CTRIAR']
    if 'CQUAD4' in bdf.card_count:
        nelements += bdf.card_count['CQUAD4'] * 2
    if 'CQUADR' in bdf.card_count:
        nelements += bdf.card_count['CQUADR'] * 2
    assert nelements > 0, bdf.card_count
    return nnodes, nelements
