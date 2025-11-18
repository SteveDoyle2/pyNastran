"""
defines:
 - tecplot = nastran_to_tecplot(model)
 - tecplot = nastran_to_tecplot_filename(bdf_filename, tecplot_filename, z
                                         log=None, debug=False)

"""
import numpy as np
from pyNastran.bdf.bdf import BDF
from pyNastran.converters.tecplot.tecplot import Tecplot, Zone


def nastran_to_tecplot(model: BDF) -> Tecplot:
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
    zone.zone_data = xyz

    nquad = model.card_count['CQUAD4'] if 'CQUAD4' in model.card_count else 0
    ntet = model.card_count['CTETRA'] if 'CTETRA' in model.card_count else 0
    #ntria = model.card_count['CTRIA3'] if 'CTRIA3' in model.card_count else 0
    nhexa = model.card_count['CHEXA'] if 'CHEXA' in model.card_count else 0
    nelement = len(model.elements)
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
        print(f'ignoring {etype}')

    # only supports nodal results
    #tecplot.nodal_results = vstack([pids, mids])#.T
    #print(tecplot.nodal_results.shape)
    #tecplot.result_names = ['PropertyID', 'MaterialID']

    ntri = len(tris)
    nquad = len(quads)
    tri_elements = np.array(tris, dtype='int32') - 1
    quad_elements = np.array(quads, dtype='int32') - 1
    del tris, quads

    ntet = len(tets)
    npenta = len(pentas)
    nhexa = len(hexas)
    tet_elements = np.array(tets, dtype='int32') - 1
    penta_elements = np.array(pentas, dtype='int32') - 1
    hexa_elements = np.array(hexas, dtype='int32') - 1
    del tets, pentas, hexas

    nshell = ntri + nquad
    nsolid = ntet + npenta + nhexa
    nnot_tri = nquad
    nnot_quad = ntri
    nnot_tet = npenta + nhexa
    nnot_hexa = ntet + npenta

    if ntri and not nnot_tri and not nsolid:
        zone.tri_elements = tri_elements
        zone_type = 'FETRIANGLE'
    elif nquad and not nnot_quad and not nsolid:
        zone.quad_elements = quad_elements
        zone_type = 'FEQUADRILATERAL'
    elif ntet and not nnot_tet and not nshell:
        zone.tet_elements = tet_elements
        zone_type = 'FETETRAHEDRON'
    elif npenta and not nhexa and not nshell:
        zone.hexa_elements = penta_elements
        zone_type = 'FEBRICK'
    elif nhexa and not nnot_hexa and not nshell:
        zone.hexa_elements = hexa_elements
        zone_type = 'FEBRICK'
    elif not nshell:
        elements = np.zeros((nelement, 8), dtype='int32')
        if ntet:
            elements[:ntet, :4] = tet_elements
            elements[:ntet, 4:] = tet_elements[:, 3]
            # elements[:ntet, 4] = tet_elements[:, 3]
            # elements[:ntet, 5] = tet_elements[:, 3]
            # elements[:ntet, 6] = tet_elements[:, 3]
            # elements[:ntet, 7] = tet_elements[:, 3]
        if npenta:
            # penta6
            n0 = ntet
            n1 = ntet + npenta
            elements[n0:n1, :6] = penta_elements
            elements[n0:n1, 6:] = penta_elements[:, 5]
            # elements[n0:n1, 6] = penta_elements[:, 5]
            # elements[n0:n1, 7] = penta_elements[:, 5]
        if nhexa:
            n0 = ntet + npenta
            n1 = ntet + npenta + nhexa
            elements[n0:n1, :] = hexa_elements
        zone.hexa_elements = elements
        zone_type = 'FEBRICK'
    elif not nsolid:
        elements = np.zeros((nelement, 4), dtype='int32')
        tris = np.array(tris, dtype='int32') - 1
        elements[:ntri, :3] = tris
        elements[:ntri, 4] = elements[:ntet, 3]

        quads = np.array(quads, dtype='int32') - 1
        elements[ntri:, :] = quads
        zone_type = 'FEQUADRILATERAL'
    else:
        msg = 'Only solids or shells are allowed (not both)\n'
        msg += f'  nsolids={nsolid:d} nshells={nshell:d}\n'
        msg += f'  ntris={ntri:d} nquads={nquad:d}\n'
        msg += f'  ntets={ntet:d} npentas={npenta:d} nhexas={nhexa:d}\n'
        raise NotImplementedError(msg)
    zone.headers_dict['ZONETYPE'] = zone_type

    tecplot.zones = [zone]
    str(zone)
    return tecplot

def nastran_to_tecplot_filename(bdf_filename, tecplot_filename, log=None, debug=False):
    """converts a BDF file to Tecplot format; supports solid elements"""
    model = BDF(log=log, debug=debug)
    model.read_bdf(bdf_filename)
    # tecplot = nastran_to_tecplot(model)

    #log.info('card_count = %s' % model.card_count)
    nnodes = len(model.nodes)
    nodes = np.zeros((nnodes, 3), dtype='float64')

    i = 0
    nodeid_to_i_map = {}
    for node_id, node in sorted(model.nodes.items()):
        xyz = node.get_position()
        nodes[i, :] = xyz
        nodeid_to_i_map[node_id] = i
        i += 1
    assert len(model.nodes) == i, 'model.nodes=%s i=%s' % (len(model.nodes), i)

    elements_list = []
    zone_type = 'FEBRICK'
    for unused_eid, element in sorted(model.elements.items()):
        if element.type in ['CTETRA']:
            n1, n2, n3, n4 = element.node_ids
            i1, i2, i3, i4 = (nodeid_to_i_map[n1], nodeid_to_i_map[n2],
                              nodeid_to_i_map[n3], nodeid_to_i_map[n4])
            elements_list.append([i1, i2, i3, i4,
                                  i4, i4, i4, i4])
        elif element.type in ['CPENTA']:
            n1, n2, n3, n4, n5, n6 = element.node_ids
            i1, i2, i3, i4, i5, i6 = (
                nodeid_to_i_map[n1], nodeid_to_i_map[n2], nodeid_to_i_map[n3], nodeid_to_i_map[n4],
                nodeid_to_i_map[n5], nodeid_to_i_map[n6])
            elements_list.append([i1, i2, i3, i4,
                                  i5, i6, i6, i6])
        elif element.type in ['CPYRAM']:
            n1, n2, n3, n4, n5 = element.node_ids
            i1, i2, i3, i4, i5 = (
                nodeid_to_i_map[n1], nodeid_to_i_map[n2], nodeid_to_i_map[n3], nodeid_to_i_map[n4],
                nodeid_to_i_map[n5])
            elements_list.append([i1, i2, i3, i4,
                                  i5, i5, i5, i5])
        elif element.type in ['CHEXA']:
            n1, n2, n3, n4, n5, n6, n7, n8 = element.node_ids
            i1, i2, i3, i4, i5, i6, i7, i8 = (
                nodeid_to_i_map[n1], nodeid_to_i_map[n2], nodeid_to_i_map[n3], nodeid_to_i_map[n4],
                nodeid_to_i_map[n5], nodeid_to_i_map[n6], nodeid_to_i_map[n7], nodeid_to_i_map[n8])
            elements_list.append([i1, i2, i3, i4,
                                  i5, i6, i7, i8])
        else:
            model.log.info('skip etype=%r' % element.type)
            model.log.info(element)
    elements = np.array(elements_list, dtype='int32')

    tecplot = Tecplot(log=model.log)
    zone = Zone(model.log)
    zone.headers_dict['ZONETYPE'] = zone_type
    zone.headers_dict['VARIABLES'] = ['X', 'Y', 'Z']
    zone.zone_data = nodes
    zone.hexa_elements = elements
    tecplot.zones = [zone]
    str(zone)
    tecplot.write_tecplot(tecplot_filename)
    return tecplot
