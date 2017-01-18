from __future__ import print_function
import os
from six import iteritems
from six.moves import range, zip
import numpy as np


class SU2Reader(object):
    def __init__(self, log=None, debug=False):
        self.log = log
        self.debug = debug

    def read_2d(self, f, ndim):
        # elements
        nelem = int(f.readline().split('=')[1])
        tris = []
        quads = []
        for ne in range(nelem):
            # what's the 0th slot?
            #Line           3
            #Triangle       5
            #Quadrilateral  9
            #Tetrahedral    10
            #Hexahedral     12
            #Wedge          13
            #Pyramid        14
            data = f.readline().split()[1:-1]
            Type = data[0]
            nodes = data[1:]
            if Type == '9':
                quads.append(nodes)
            elif Type == '5':
                tris.append(nodes)
        tris = np.array(tris, dtype='int32')
        quads = np.array(quads, dtype='int32')

        nnodes = int(f.readline().split('=')[1])
        nodes = np.zeros((nelem, 2), dtype='int32')
        for inode in range(nnodes):
            sline = f.readline().split()
            assert len(sline) == 3, sline
            x, y, z = sline
            #print(x, y, z)
            nodes[inode, :] = [float(x), float(y)]

        nmark = int(f.readline().split('=')[1])
        for imark in range(nmark):
            marker = f.readline().split('=')[1].strip()
            nelements_mark = int(f.readline().split('=')[1])

            if ndim == 2:
                lines = np.zeros((nelements_mark, 2), dtype='int32')
                for ne in range(nelements_mark):
                    # what are the 3 slots?
                    #Line          (2D)     3
                    #Triangle      (3D)     5
                    #Quadrilateral (3D)     9
                    Type, n1, n2 = f.readline().split()
                    lines[ne] = [n1, n2]

        elements = {
            5 : tris,
            9 : quads,
        }
        regions = {3 : lines}
        return nodes, elements, regions

    def read_3d(self, f, ndim):
        nelem = int(f.readline().split('=')[1])
        tets = []
        hexs = []
        wedges = []
        pyramids = []
        for ne in range(nelem):
            data = f.readline().split()[1:-1]
            Type = data[0]
            nodes = data[1:]
            if Type == '10':
                tets.append(nodes)
            elif Type == '12':
                hexs.append(nodes)
            elif Type == '13':
                wedges.append(nodes)
            elif Type == '14':
                pyramids.append(nodes)
            else:
                raise NotImplementedError(Type)
        tets = np.array(tets, dtype='int32')
        pents = np.array(pents, dtype='int32')
        wedges = np.array(wedges, dtype='int32')
        pyramids = np.array(pyramids, dtype='int32')
        elements = {
            10 : tets,
            12 : hexs,
            13 : wedges,
            14 : pyramids,
        }

        nnodes = int(f.readline().split('=')[1])
        nodes = zeros((nelem, 2), dtype='int32')
        for inode in range(nnodes):
            x, y, z = f.readline().split()[:-1]
            nodes[inode, :] = [x, y, z]

        nmark = int(f.readline().split('=')[1])
        for imark in nmark:
            marker = f.readline().split('=')[1].strip()
            nelements_mark = int(f.readline().split('=')[1])
            tris = []
            quads = []
            for ne in range(nelements_mark):
                data = f.readline().split()[1:-1]
                Type = data[0]
                nodes = data[1:]
                if Type == '9':
                    quads.append(nodes)
                elif Type == '5':
                    tris.append(nodes)
            tris = np.array(tris, dtype='int32')
            quads = np.array(quads, dtype='int32')

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
            nodes, elements, regions = self.read_2d(f, ndim)
        elif ndim == 3:
            nodes, elements, regions = self.read_3d(f, ndim)
        else:
            raise RuntimeError(ndim)
        return ndim, nodes, elements, regions

    def write_su2(self, su2_filename, nodes, elements, regions):
        nnodes, ndim = nodes.shape
        nnodes, ndim = nodes.shape
        f = open(su2_filename, 'wb')

        f.write('NDIM = %i\n' % ndim)
        self.Type_nnodes_map = {
            #Line       3
            #Triangle   5
            #Quadrilateral  9
            #Tetrahedral    10
            #Hexahedral 12
            #Wedge      13
            #Pyramid    14
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


def main(su2_filename):
    su2 = SU2Reader()
    su2.read_su2(su2_filename)
    #su2.to_cart3d()

if __name__ == '__main__':
    main('mesh_naca0012_inv.su2')

