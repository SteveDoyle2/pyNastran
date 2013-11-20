from pyNastran.bdf.bdf import BDF


def nastran_to_cart3d(bdf_filename, cart3d_filename):
    model = BDF()
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
    nastran_to_cart3d(bdf_filename, cart3d_filename)