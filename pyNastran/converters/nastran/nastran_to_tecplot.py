from six import iteritems
from numpy import zeros, array
from pyNastran.bdf.bdf import BDF
from pyNastran.converters.tecplot.tecplot import Tecplot

def nastran_to_tecplot_filename(bdf_filename, tecplot_filename, log=None):
    model = BDF(log=log)
    model.read_bdf(bdf_filename)
    #log.info('card_count = %s' % model.card_count)

    nnodes = len(model.nodes)
    nodes = zeros((nnodes, 3), dtype='float64')
    elements = []

    i = 0
    nodeid_to_i_map = {}
    for node_id, node in sorted(iteritems(model.nodes)):
        xyz = node.get_position()
        nodes[i, :] = xyz
        nodeid_to_i_map[node_id] = i
        i += 1
    assert len(model.nodes) == i, 'model.nodes=%s i=%s' % (len(model.nodes), i)

    for eid, element in sorted(iteritems(model.elements)):
        if element.type in ['CTETRA']:
            n1, n2, n3, n4 = element.node_ids
            i1, i2, i3, i4 = nodeid_to_i_map[n1], nodeid_to_i_map[n2], nodeid_to_i_map[n3], nodeid_to_i_map[n4]
            elements.append([i1, i2, i3, i4,
                             i4, i4, i4, i4])
        elif element.type in ['CPENTA']:
            n1, n2, n3, n4, n5, n6 = element.node_ids
            i1, i2, i3, i4, i5, i6 = (
                nodeid_to_i_map[n1], nodeid_to_i_map[n2], nodeid_to_i_map[n3], nodeid_to_i_map[n4]
                nodeid_to_i_map[n5], nodeid_to_i_map[n6])
            elements.append([i1, i2, i3, i4,
                             i5, i6, i6, i6])
        elif element.type in ['CYPRAM']:
            n1, n2, n3, n4, n5 = element.node_ids
            i1, i2, i3, i4, i5 = (
                nodeid_to_i_map[n1], nodeid_to_i_map[n2], nodeid_to_i_map[n3], nodeid_to_i_map[n4]
                nodeid_to_i_map[n5])
            elements.append([i1, i2, i3, i4,
                             i5, i5, i5, i5])
        elif element.type in ['CHEXA']:
            n1, n2, n3, n4, n5, n6, n7, n8 = element.node_ids
            i1, i2, i3, i4, i5, i6, i7, i8 = (
                nodeid_to_i_map[n1], nodeid_to_i_map[n2], nodeid_to_i_map[n3], nodeid_to_i_map[n4]
                nodeid_to_i_map[n5], nodeid_to_i_map[n6], nodeid_to_i_map[n7], nodeid_to_i_map[n8])
            elements.append([i1, i2, i3, i4,
                             i5, i6, i7, i8])
        else:
            print(element.type)
    elements = array(elements, dtype='int32')
    tecplot = Tecplot()
    tecplot.xyz = nodes
    tecplot.elements = elements
    tecplot.write_tecplot(tecplot_filename)
    tecplot.results = array([], dtype='float32')


if __name__ == '__main__':  # pragma: no cover
    bdf_filename = 'threePlugs.bdf'
    tecplot_filename = 'threePlugs.plt'
    nastran_to_tecplot_filename(bdf_filename, tecplot_filename)
