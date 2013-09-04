from pyNastran.bdf.bdf import BDF


def to_cart3d(fname):
    BDF()
    BDF.read_bdf(fname)
    nnodes = len(model.nodes)
    nelements = len(model.elements)

    fname = open(fname, 'wb')
    fname.write('%s %s\n' % (nnodes, nelements)
    node_id_shift = {}
    i = 1
    for node_id, node in sorted(model.nodes()):
        node_id_shift[node_id] = i
        x, y, z = node.Position()
        fname.write('%s %s %s\n' % (x, y, z))

    mids = ''
    for element_id, element in sorted(model.elements):
        assert element.type == 'CTRIA3', element.type
        n1, n2, n3 = element.node_ids()
        n1 = node_id_shift[n1]
        n2 = node_id_shift[n2]
        n3 = node_id_shift[n3]
        mid = element.material_id()
        fname.write('%s %s %s\n' % (n1, n2, n3))
        mids += '%s ' % mid
    f.write(mids)
    f.close()

if __name__ == '__main__':
    fname = 'out.bdf'
    to_cart3d(fname)