from __future__ import print_function
from pyNastran.op2.fortran_format import FortranFormat
from numpy import array, zeros, hstack, vstack

class FGridReader(object):
    def __init__(self, log=None, debug=False):
        self.log = log
        self.debug = debug

    def read_fgrid(self, fgrid_filename, dimension_flag=3):
        self.f = open(fgrid_filename, 'r')
        nnodes, ntris, ntets = self.f.readline().split()
        nnodes = int(nnodes)
        ntris = int(ntris)
        ntets = int(ntets)

        #nodes = zeros((nnodes, 3), dtype='float32')
        nodesi = zeros((nnodes * 3), dtype='float32')
        #x = zeros(nnodes, dtype='float32')
        #y = zeros(nnodes, dtype='float32')
        #z = zeros(nnodes, dtype='float32')
        tris = zeros((ntris, 3), dtype='int32')

        i0 = 0
        # I think this goes xxx, yyy, zzz
        # instead of x, y, z
        #            x, y, z
        for i in range(nnodes):
            #nodes[i, :] = self.f.readline().split()
            #nodes[i0:i1] = self.f.readline().split()
            sline = self.f.readline().split()
            #x[i] = sline[0]
            #y[i] = sline[1]
            #z[i] = sline[2]
            nodesi[i0] = sline[0]
            nodesi[i0 + 1] = sline[1]
            nodesi[i0 + 2] = sline[2]
            i0 += 3

        #nodes = vstack([x, y, z]).T
        #nodes = hstack([x, y, z]).reshape((3, nnodes)).T
        #nodes = hstack([x, y, z]).reshape((3, nnodes)).T
        #nodes = nodesi.reshape((nnodes, 3))
        nodes = nodesi.reshape((3, nnodes)).T
        print('nodes.shape', nodes.shape)
        ntri_lines = ntris // 2
        j = 0
        for i in range(ntri_lines):
            sline = self.f.readline().split()
            tris[j, :] = sline[:3]
            tris[j + 1, :] = sline[3:]
            j += 2

        k = 0
        if ntris % 2: # if there's a leftover line
            tris[j] = sline[:3]
            #bcs[:3] = sline[3:]
            k = 3

        nquestion = ntris // 3 // 2
        for i in range(nquestion):
            sline = self.f.readline().split()

        #sline = self.f.readline().split()
        # nnodes=16391
        # ntris=4512
        # ntets = 90892
        # 60595 lines * 6/4 = 90892.5
        #752, 753 = 4512/3/2
        #assert sline == ['15770','15767','16242','15772','15765','15762'], sline
        ntet_lines = ntets * 6 // 44

        tets = zeros(ntets * 4, dtype='int32')
        tets[:] = self.f.read().split()
        tets = tets.reshape((ntets, 4))

        self.nodes = nodes
        self.tris = tris
        self.tets = tets

def main():
    f = FGridReader()
    f.read_fgrid()

if __name__ == '__main__':
    main()
