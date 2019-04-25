from __future__ import print_function
import numpy as np


def read_su2(su2_filename, log=None, debug=False):
    model = SU2Reader(log=log, debug=debug)
    unused_ndim, nodes, elements, regions = model.read_su2(su2_filename)
    #model.to_cart3d()
    return model, nodes, elements, regions

class SU2Reader(object):
    etype_nnodes_map = {
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
    def __init__(self, log=None, debug=False):
        self.log = log
        self.debug = debug
        self.su2_filename = None

    def read_2d(self, su2_file, ndim):
        # elements
        nelem = int(su2_file.readline().split('=')[1])
        tris = []
        quads = []
        for unused_ne in range(nelem):
            # what's the 0th slot?
            #Line           3
            #Triangle       5
            #Quadrilateral  9
            #Tetrahedral    10
            #Hexahedral     12
            #Wedge          13
            #Pyramid        14
            data = su2_file.readline().split()[:-1]
            etype = data[0]
            nodes = data[1:]
            if etype == '9':
                assert len(nodes) == 4, nodes
                quads.append(nodes)
            elif etype == '5':
                assert len(nodes) == 3, nodes
                tris.append(nodes)
            else:
                raise NotImplementedError(etype)

        tris = np.asarray(tris, dtype='int32')
        quads = np.asarray(quads, dtype='int32')
        #print('tris =', tris)
        #print('quads =', quads)

        nnodes = int(su2_file.readline().split('=')[1])
        nodes = np.zeros((nelem, 2), dtype='int32')
        for inode in range(nnodes):
            sline = su2_file.readline().split()
            assert len(sline) == 3, sline
            x, y, unused_z = sline
            #print(x, y, z)
            nodes[inode, :] = [float(x), float(y)]

        # boundary conditions
        nmark = int(su2_file.readline().split('=')[1])
        for unused_imark in range(nmark):
            unused_marker = su2_file.readline().split('=')[1].strip()
            nelements_mark = int(su2_file.readline().split('=')[1])

            if ndim == 2:
                lines = []
                #np.zeros((nelements_mark, 2), dtype='int32')
                for unused_ne in range(nelements_mark):
                    # what are the 3 slots?
                    #Line          (2D)     3
                    #Triangle      (3D)     5
                    #Quadrilateral (3D)     9
                    sline = su2_file.readline().split()
                    etype = int(sline[0])
                    if etype == 3:
                        etype, n1, n2 = sline
                        lines.append([n1, n2])
                    else:
                        raise NotImplementedError(etype)

        assert len(tris) > 0 or len(quads) > 0
        elements = {
            5 : tris,
            9 : quads,
        }
        regions = {3 : lines}
        return nodes, elements, regions

    def read_3d(self, su2_file, unused_ndim):
        nelem = int(su2_file.readline().split('=')[1])
        tets = []
        hexs = []
        wedges = []
        #pents = []
        pyramids = []
        for unused_ne in range(nelem):
            data = su2_file.readline().split()[1:-1]
            etype = data[0]
            nodes = data[1:]
            if etype == '10':
                tets.append(nodes)
            elif etype == '12':
                hexs.append(nodes)
            elif etype == '13':
                wedges.append(nodes)
            elif etype == '14':
                pyramids.append(nodes)
            else:
                raise NotImplementedError(etype)
        tets = np.array(tets, dtype='int32')
        #pents = np.array(pents, dtype='int32')
        wedges = np.array(wedges, dtype='int32')
        pyramids = np.array(pyramids, dtype='int32')
        elements = {
            10 : tets,
            12 : hexs,
            13 : wedges,
            14 : pyramids,
        }

        nnodes = int(su2_file.readline().split('=')[1])
        nodes = np.zeros((nelem, 2), dtype='int32')
        for inode in range(nnodes):
            x, y, z = su2_file.readline().split()[:-1]
            nodes[inode, :] = [x, y, z]

        nmark = int(su2_file.readline().split('=')[1])
        for unused_imark in range(nmark):
            unused_marker = su2_file.readline().split('=')[1].strip()
            nelements_mark = int(su2_file.readline().split('=')[1])
            tris = []
            quads = []
            for unused_ne in range(nelements_mark):
                data = su2_file.readline().split()[1:-1]
                etype = data[0]
                nodes = data[1:]
                if etype == '9':
                    quads.append(nodes)
                elif etype == '5':
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
        with open(su2_filename, 'r') as su2_file:
            ndim = int(su2_file.readline().split('=')[1])

            if ndim == 2:
                nodes, elements, regions = self.read_2d(su2_file, ndim)
            elif ndim == 3:
                nodes, elements, regions = self.read_3d(su2_file, ndim)
            else:
                raise RuntimeError(ndim)
        return ndim, nodes, elements, regions

    def write_su2(self, su2_filename, nodes, elements, unused_regions):
        nnodes, ndim = nodes.shape
        nnodes, ndim = nodes.shape
        with open(su2_filename, 'wb') as su2_file:
            su2_file.write('NDIM = %i\n' % ndim)
            if ndim == 2:
                for etype, elementsi in sorted(elements.items()):
                    element_num = self.etype_nnodes_map[etype]
                    fmt = '%%s' + ' %%s' * (element_num-1) + '\n'
                    for element in elementsi:
                        su2_file.write(fmt % element)

                su2_file.write('NPOINTS = %i\n' % nnodes)
                for inode, node in enumerate(nodes):
                    su2_file.write('%i %i %i\n' % (node[0], node[1], inode))
            elif ndim == 3:
                su2_file.write('NPOINTS = %i\n' % nnodes)
                for inode, node in enumerate(nodes):
                    su2_file.write('%i %i %i %i\n' % (node[0], node[1], node[2], inode))
