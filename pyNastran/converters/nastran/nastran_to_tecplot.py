"""
defines:
 - tecplot = nastran_to_tecplot(model)
 - tecplot = nastran_to_tecplot_filename(bdf_filename, tecplot_filename, z
                                         log=None, debug=False)

"""
import numpy as np
from pyNastran.bdf.bdf import BDF
from pyNastran.converters.tecplot.tecplot import Tecplot, Zone


def nastran_to_tecplot(model):
    """assumes sequential nodes"""
    tecplot = Tecplot(log=model.log)

    nnodes = len(model.nodes)
    inode_max = max(model.nodes)
    if nnodes == inode_max:
        xyz = np.zeros((nnodes, 3), dtype='float64')
        i = 0
        for unused_nid, node in sorted(model.nodes.items()):
            xyz[i, :] = node.get_position()
            i += 1
    else:
        msg = 'sequential node IDs required; nnodes=%s inode_max=%s' % (
            nnodes, inode_max)
        raise RuntimeError(msg)
    zone = Zone(model.log)
    zone.headers_dict['VARIABLES'] = ['X', 'Y', 'Z']
    zone.xyz = xyz

    nquads = model.card_count['CQUAD4'] if 'CQUAD4' in model.card_count else 0
    ntets = model.card_count['CTETRA'] if 'CTETRA' in model.card_count else 0
    #ntrias = model.card_count['CTRIA3'] if 'CTRIA3' in model.card_count else 0
    nhexas = model.card_count['CHEXA'] if 'CHEXA' in model.card_count else 0
    nelements = len(model.elements)
    tris = []
    quads = []
    tets = []
    hexas = []
    pentas = []
    #i = 0
    #pids = np.zeros(nelements, dtype='int32')
    #mids = np.zeros(nelements, dtype='int32')
    unhandled_types = set()
    for unused_eid, element in model.elements.items():
        if element.type in ['CTRIA3']:
            tris.append(element.node_ids)
        elif element.type in ['CQUAD4']:
            quads.append(element.node_ids)
        elif element.type == 'CTETRA':
            tets.append(element.node_ids[:4])
        elif element.type == 'CPENTA':
            pentas.append(element.node_ids[:6])
        elif element.type == 'CHEXA':
            hexas.append(element.node_ids[:8])
        else:
            unhandled_types.add(element.type)
        #pid = element.Pid()
        #mid = element.Mid()
        #pids[i] = pid
        #mids[i] = mid
        #i += 1

    for etype in unhandled_types:
        print('ignoring %s' % etype)

    # only supports nodal results
    #tecplot.nodal_results = vstack([pids, mids])#.T
    #print(tecplot.nodal_results.shape)
    #tecplot.result_names = ['PropertyID', 'MaterialID']

    ntris = len(tris)
    nquads = len(quads)
    nshells = ntris + nquads

    ntets = len(tets)
    npentas = len(pentas)
    nhexas = len(hexas)
    nsolids = ntets + npentas + nhexas
    nnot_tris = nquads
    nnot_quads = ntris
    nnot_tets = npentas + nhexas
    nnot_hexas = ntets + npentas
    if ntris and not nnot_tris and not nsolids:
        zone.tri_elements = np.array(tris, dtype='int32')
    elif nquads and not nnot_quads and not nsolids:
        zone.quad_elements = np.array(quads, dtype='int32')
    elif ntets and not nnot_tets and not nshells:
        zone.tet_elements = np.array(tets, dtype='int32')
    elif nhexas and not nnot_hexas and not nshells:
        zone.hexa_elements = np.array(hexas, dtype='int32')
    elif not nshells:
        elements = np.zeros((nelements, 8), dtype='int32')
        if ntets:
            tets = np.array(tets, dtype='int32')
            elements[:ntets, :4] = tets
            elements[:ntets, 4] = elements[:ntets, 3]
            elements[:ntets, 5] = elements[:ntets, 3]
            elements[:ntets, 6] = elements[:ntets, 3]
            elements[:ntets, 7] = elements[:ntets, 3]
        if npentas:
            # penta6
            pentas = np.array(pentas, dtype='int32')
            elements[ntets:ntets + npentas, :6] = pentas
            elements[ntets:ntets + npentas, 6] = elements[:ntets, 5]
            elements[ntets:ntets + npentas, 7] = elements[:ntets, 5]
        if nhexas:
            hexas = np.array(hexas, dtype='int32')
            elements[ntets + npentas:ntets + npentas + nhexas, :6] = pentas
            elements[ntets + npentas:ntets + npentas + nhexas, 6] = elements[:ntets, 5]
            elements[ntets + npentas:ntets + npentas + nhexas, 7] = elements[:ntets, 5]
        zone.hexa_elements = np.array(elements)
    elif not nsolids:
        elements = np.zeros((nelements, 4), dtype='int32')
        tris = np.array(tris, dtype='int32')
        elements[:ntris, :3] = tris
        elements[:ntris, 4] = elements[:ntets, 3]

        quads = np.array(quads, dtype='int32')
        elements[ntris:, :] = quads
    else:
        msg = 'Only solids or shells are allowed (not both)\n'
        msg += '  nsolids=%s nshells=%s\n' % (nsolids, nshells)
        msg += '  ntris=%s nquads=%s\n' % (ntris, nquads)
        msg += '  ntets=%s npentas=%s nhexas=%s\n' % (ntets, npentas, nhexas)
        raise NotImplementedError(msg)
    tecplot.zones = [zone]
    return tecplot

def nastran_to_tecplot_filename(bdf_filename, tecplot_filename, log=None, debug=False):
    """converts a BDF file to Tecplot format; supports solid elements"""
    model = BDF(log=log, debug=debug)
    model.read_bdf(bdf_filename)
    # tecplot = nastran_to_tecplot(model)

    #log.info('card_count = %s' % model.card_count)
    nnodes = len(model.nodes)
    nodes = np.zeros((nnodes, 3), dtype='float64')
    elements = []

    i = 0
    nodeid_to_i_map = {}
    for node_id, node in sorted(model.nodes.items()):
        xyz = node.get_position()
        nodes[i, :] = xyz
        nodeid_to_i_map[node_id] = i
        i += 1
    assert len(model.nodes) == i, 'model.nodes=%s i=%s' % (len(model.nodes), i)

    for unused_eid, element in sorted(model.elements.items()):
        if element.type in ['CTETRA']:
            n1, n2, n3, n4 = element.node_ids
            i1, i2, i3, i4 = (nodeid_to_i_map[n1], nodeid_to_i_map[n2],
                              nodeid_to_i_map[n3], nodeid_to_i_map[n4])
            elements.append([i1, i2, i3, i4,
                             i4, i4, i4, i4])
        elif element.type in ['CPENTA']:
            n1, n2, n3, n4, n5, n6 = element.node_ids
            i1, i2, i3, i4, i5, i6 = (
                nodeid_to_i_map[n1], nodeid_to_i_map[n2], nodeid_to_i_map[n3], nodeid_to_i_map[n4],
                nodeid_to_i_map[n5], nodeid_to_i_map[n6])
            elements.append([i1, i2, i3, i4,
                             i5, i6, i6, i6])
        elif element.type in ['CPYRAM']:
            n1, n2, n3, n4, n5 = element.node_ids
            i1, i2, i3, i4, i5 = (
                nodeid_to_i_map[n1], nodeid_to_i_map[n2], nodeid_to_i_map[n3], nodeid_to_i_map[n4],
                nodeid_to_i_map[n5])
            elements.append([i1, i2, i3, i4,
                             i5, i5, i5, i5])
        elif element.type in ['CHEXA']:
            n1, n2, n3, n4, n5, n6, n7, n8 = element.node_ids
            i1, i2, i3, i4, i5, i6, i7, i8 = (
                nodeid_to_i_map[n1], nodeid_to_i_map[n2], nodeid_to_i_map[n3], nodeid_to_i_map[n4],
                nodeid_to_i_map[n5], nodeid_to_i_map[n6], nodeid_to_i_map[n7], nodeid_to_i_map[n8])
            elements.append([i1, i2, i3, i4,
                             i5, i6, i7, i8])
        else:
            model.log.info('skip etype=%r' % element.type)
            model.log.info(element)
    elements = np.array(elements, dtype='int32')

    tecplot = Tecplot(log=model.log)
    zone = Zone(model.log)
    zone.headers_dict['VARIABLES'] = ['X', 'Y', 'Z']
    zone.xyz = nodes
    zone.hexa_elements = elements
    zone.nodal_results = np.array([], dtype='float32')
    tecplot.zones = [zone]
    tecplot.write_tecplot(tecplot_filename)
    return tecplot
