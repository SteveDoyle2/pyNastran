from pyNastran.bdf.bdf import BDF


def nastran_to_stl_filename(bdf_filename, stl_filename, log=None):
    model = BDF(log=log)
    model.read_bdf(bdf_filename)
    log.info('card_count = %s' % model.card_count)

    f = open(stl_filename, 'wb')
    positions = {}

    for node_id, node in sorted(model.nodes.iteritems()):
        positions[node_id] = node.Position()

    f.write('solid nastran_model\n')
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

#solid nastran_model
#   facet normal -6.665299e-001 6.795624e-001 3.064844e-001
#      outer loop
#         vertex 8.142845e-002 2.731541e-001 1.190024e+001
#         vertex 8.186898e-002 2.727136e-001 1.190215e+001
#         vertex 8.467505e-002 2.754588e-001 1.190215e+001
#      endloop
#   endfacet
        n1 = positions[n1]
        n2 = positions[n2]
        n3 = positions[n3]

        n = element.Normal()
        #f.write('solid nastran_model\n')
        f.write('  facet normal %f %f %f %i\n' % ( n[0], n[1], n[2], element_id ))
        f.write('    outer loop\n')
        f.write('      vertex %f %f %f\n' % ( n1[0], n1[1], n1[2] ))
        f.write('      vertex %f %f %f\n' % ( n2[0], n2[1], n2[2] ))
        f.write('      vertex %f %f %f\n' % ( n3[0], n3[1], n3[2] ))
        f.write('    endloop\n')
        f.write('  endfacet\n')
    f.write('endsolid\n')
    f.close()

if __name__ == '__main__':
    bdf_filename = 'g278.bdf'
    stl_filename = 'g278.stl'
    nastran_to_stl(bdf_filename, stl_filename)