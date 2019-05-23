"""
defines:
 - stl = nastran_to_stl_filename(bdf_filename, stl_filename, is_binary=False,
                                 log=None)
 - stl = nastran_to_stl(bdf_filename, stl_filename, is_binary=False,
                        log=None, stop_on_failure=False)

"""
import numpy as np
from pyNastran.bdf.bdf import read_bdf
from pyNastran.converters.stl.stl import STL

def nastran_to_stl_filename(bdf_filename, stl_filename, is_binary=False, log=None):
    """Converts a Nastran model to an STL"""
    return nastran_to_stl(bdf_filename, stl_filename, is_binary=is_binary)

def nastran_to_stl(bdf_filename, stl_filename, is_binary=False, log=None, stop_on_failure=False):
    """
    Converts a Nastran model to an STL

    Parameters
    ----------
    bdf_filename : varies
        str : the path to a BDF input file
        BDF() : a BDF() model object
    stl_filename : str
        the output STL path
    is_binary : bool; default=False
        should the output file be binary
    log : Logger()
        a Python logging object
    stop_on_failure : bool; default=False
        should the code stop if an error is encountered

    """
    if isinstance(bdf_filename, str):
        model = read_bdf(bdf_filename, log=log)
    else:
        model = bdf_filename

    #log.info('card_count = %s' % model.card_count)

    nnodes = len(model.nodes)
    nodes = np.zeros((nnodes, 3), dtype='float64')
    elements = []

    i = 0
    nodeid_to_i_map = {}
    offset = False
    if offset:
        nid = list(model.nodes.keys())[0]
        xyz0 = model.nodes[nid].get_position()
    else:
        xyz0 = np.zeros(3, dtype='float64')
    for node_id, node in sorted(model.nodes.items()):
        xyz = node.get_position()
        nodes[i, :] = xyz - xyz0
        nodeid_to_i_map[node_id] = i
        i += 1
    assert len(model.nodes) == i, 'model.nodes=%s i=%s' % (len(model.nodes), i)
    for unused_eid, element in sorted(model.elements.items()):
        if element.type in ['CQUADR']:
            continue
        elif element.type in ['CBAR', 'CBEAM', 'CONM2', 'RBE2', 'RBE3',
                              'CBUSH', 'CBUSH1D', 'CBUSH2D',
                              'CONROD', 'CROD',
                              'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
                              'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4',]:
            continue
        elif element.type in ['CQUAD4']:
            n1, n2, n3, n4 = element.node_ids
            i1, i2, i3, i4 = (nodeid_to_i_map[n1], nodeid_to_i_map[n2],
                              nodeid_to_i_map[n3], nodeid_to_i_map[n4])
            elements.append([i1, i2, i3])
            elements.append([i3, i4, i1])
        elif element.type in ['CTRIA3', 'CTRIAR']:
            nids = element.node_ids
            unids = np.unique(nids)
            if len(unids) == 2:
                continue
            n1, n2, n3 = nids
            i1, i2, i3 = nodeid_to_i_map[n1], nodeid_to_i_map[n2], nodeid_to_i_map[n3]
            elements.append([i1, i2, i3])
        else:
            model.log.warning('skipping %s' % element.type)
    elements = np.array(elements, dtype='int32')
    stl = STL(log=model.log)
    stl.nodes = nodes
    #stl.nodes -= nodes[0, :]
    stl.elements = elements
    stl.write_stl(stl_filename, is_binary=is_binary, stop_on_failure=stop_on_failure)
    return stl
