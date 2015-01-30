
import os
from six import iteritems
from six.moves import range, zip
from numpy import array

def main(su2_filename):
    su2 = SU2Reader()
    su2.read_su2(su2_filename)
    su2.to_cart3d()

class SU2Reader(object):
    def __init__(self, log=None, debug=False):
        self.log = log
        self.debug = debug

    def read_2d(self, f, nelem):
        # elements
        nelem = int(f.readline().split('=')[1])
        tris = []
        quads = []
        for ne in range(nelem):
            # what's the 0th slot?
            #Line   	3
            #Triangle	5
            #Quadrilateral	9
            #Tetrahedral	10
            #Hexahedral	12
            #Wedge  	13
            #Pyramid	14
            Type, *nodes = f.readline().split()[1:-1]
            if Type == '9':
                quads.append(nodes)
            elif Type == '5':
                tris.append(nodes)
        tris = array(tris, dtype='int32')
        quads = array(quads, dtype='int32')

        nnodes = int(f.readline().split('=')[1])
        nodes = zeros((nelem, 2), dtype='int32')
        for np in range(nnodes):
            x, y = f.readline().split()[:-1]
            nodes[ne, :] = [x, y, z]

        nmark = int(f.readline().split('=')[1])
        for imark in nmark:
            marker = f.readline().split('=')[1].strip()
            nelements_mark = int(f.readline().split('=')[1])

            if ndim == 2:
                lines = zeros(nelements_mark, 2)
                for ne in range(nelements_mark):
                    # what are the 3 slots?
                    #Line          (2D) 	3
                    #Triangle      (3D) 	5
                    #Quadrilateral (3D) 	9
                    Type, n1, n2 = f.readline().split()
                    lines[ne] = [n1, n2]

        elements = {
            5 : tris,
            9 : quads,
        }
        regions = {3 : lines}
        return nodes, elements, regions

    def read_3d(self, f, nelem):
        nelem = int(f.readline().split('=')[1])
        tets = []
        hexs = []
        wedges = []
        pyramids = []
        for ne in range(nelem):
            Type, *nodes = f.readline().split()[1:-1]
            if Type == '10':
                tets.append(nodes)
            elif Type == '12':
                hexs.append(nodes)
            elif Type == '13':
                wedges.append(nodes)
            elif Type == '14':
                pyramids.append(nodes)
        tets = array(tets, dtype='int32')
        pents = array(pents, dtype='int32')
        wedges = array(wedges, dtype='int32')
        pyramids = array(pyramids, dtype='int32')
        elements = {
            10 : tets,
            12 : hexs,
            13 : wedges,
            14 : pyramids,
        }

        nnodes = int(f.readline().split('=')[1])
        nodes = zeros((nelem, 2), dtype='int32')
        for np in range(nnodes):
            x, y, z = f.readline().split()[:-1]
            nodes[ne, :] = [x, y, z]

        nmark = int(f.readline().split('=')[1])
        for imark in nmark:
            marker = f.readline().split('=')[1].strip()
            nelements_mark = int(f.readline().split('=')[1])
            tris = []
            quads = []
            for ne in range(nelements_mark):
                Type, *nodes = f.readline().split()[1:-1]
                if Type == '9':
                    quads.append(nodes)
                elif Type == '5':
                    tris.append(nodes)
            tris = array(tris, dtype='int32')
            quads = array(quads, dtype='int32')

        regions = {
            5 : tris,
            9 : quads,
        }
        return nodes, elements, regions


    def read_su2(self, su2_filename):
        self.su2_filename = su2_filename
        f = open(su2_filename, 'r')
        ndim = int(f.readline().split('=')[1])

        if ndim == 2:
            nodes, elements, regions = read_2d(f, nelem)
        elif ndim == 3:
            nodes, elements, regions = read_3d(f, nelem)
        else:
            raise RuntimeError(ndim)
        return ndim, nodes, elements, regions

    def write_su2(self, su2_filename, nodes, elements, regions):
        nnodes, ndim = nodes.shape
        nnodes, ndim = nodes.shape
        f = open(su2_filename, 'wb')

        f.write('NDIM = %i\n' % ndim)
        self.Type_nnodes_map = {
            #Line   	3
            #Triangle	5
            #Quadrilateral	9
            #Tetrahedral	10
            #Hexahedral	12
            #Wedge  	13
            #Pyramid	14
            3 : 2,
            5 : 3,
            9 : 4,
            10 : 4,
            12 : 8,
            13 : 6,
            14 : 5,
        }
        if ndim == 2:
            for Type, elementsi in sorted(iteritems(elements)):
                n = self.Type_nnodes_map[Type]
                fmt = '%%s' + ' %%s' * (n-1) + '\n'
                for element in elementsi:
                    f.write(fmt % element)

            f.write('NPOINTS = %i\n' % nnodes)
            for inode, node in enumerate(nodes):
                f.write('%i %i %i\n' % (node[0], node[1], inode))
        elif ndim == 3:
            f.write('NPOINTS = %i\n' % nnodes)
            for inode, node in enumerate(nodes):
                f.write('%i %i %i %i\n' % (node[0], node[1], node[2], inode))

main('mesh_NACA0012_1E-4m.su2')