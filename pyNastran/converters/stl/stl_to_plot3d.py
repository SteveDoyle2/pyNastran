from six.moves import range
from numpy import savetxt

from pyNastran.converters.stl.stl import STL
#from pyNastran.bdf.field_writer_8 import print_card

def stl_to_plot3d_filename(stl_filename, p3d_filename, log=None, ascii=True):
    model = STL(log=log)
    model.read_stl(stl_filename)

    #nodal_normals = model.get_normals_at_nodes(model.elements)

    with open(p3d_filename, 'wb') as p3d:
        nblocks = len(model.elements)
        #nblocks = 10
        p3d.write('%i\n' % nblocks)
        for iblock in range(nblocks):
            p3d.write('2 2 1\n')

        nodes = model.nodes
        elements = model.elements
        if 0:
            for i in [0, 1, 2]:
                for iblock in range(nblocks):
                    (n1, n2, n3) = elements[iblock]
                    p1 = nodes[n1, :]
                    p2 = nodes[n2, :]
                    p3 = nodes[n3, :]
                    p4 = p3
                    xi = [[p1[i], p2[i], p3[i], p4[i]]]
                    savetxt(p3d, xi, fmt='%f')
        else:
            for iblock in range(nblocks):
                for i in [0, 1, 2]:
                    (n1, n2, n3) = elements[iblock]
                    p1 = nodes[n1, :]
                    p2 = nodes[n2, :]
                    p3 = nodes[n3, :]
                    p4 = p3
                    xi = [[p1[i], p2[i], p3[i], p4[i]]]
                    savetxt(p3d, xi, fmt='%f')

                #p3d.write("----\n")
                #x = [[p1[0], p2[0], p3[0], p4[0]]]
                #y = [[p1[1], p2[1], p3[1], p4[1]]]
                #z = [[p1[2], p2[2], p3[2], p4[2]]]

                #savetxt(p3d, x, fmt='%f')
                #savetxt(p3d, y, fmt='%f')
                #savetxt(p3d, z, fmt='%f')

            #p3d.write('\n')
            #p3d.write(' '.join(x) + '\n')
            #p3d.write(' '.join(y) + '\n')
            #p3d.write(' '.join(z) + '\n')
            #xxxz
            #yyyz
            #zzzz

if __name__ == '__main__':  # pragma: no cover
    stl_filename = 'spw.out.STL'
    p3d_filename = 'spw.out.p3d'
    stl_to_plot3d_filename(stl_filename, p3d_filename)
