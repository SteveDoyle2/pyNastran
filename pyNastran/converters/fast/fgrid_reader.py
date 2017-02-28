from __future__ import print_function
import os
#from pyNastran.op2.fortran_format import FortranFormat
import numpy as np
from pyNastran.utils import print_bad_path

class FGridReader(object):
    def __init__(self, log=None, debug=False):
        self.log = log
        self.debug = debug
        self.nodes = None
        self.tris = None
        self.tets = None

    def read_fgrid(self, fgrid_filename, dimension_flag=3):
        assert os.path.exists(fgrid_filename), print_bad_path(fgrid_filename)
        with open(fgrid_filename, 'r') as fgrid:
            nnodes, ntris, ntets = fgrid.readline().split()
            nnodes = int(nnodes)
            ntris = int(ntris)
            ntets = int(ntets)

            nodesi = np.zeros((nnodes * 3), dtype='float32')
            tris = np.zeros((ntris, 3), dtype='int32')
            assert nnodes > 0, nnodes
            i0 = 0
            # I think this goes xxx, yyy, zzz
            # instead of x, y, z
            #            x, y, z
            for i in range(nnodes):
                #nodes[i, :] = fgrid.readline().split()
                #nodes[i0:i1] = fgrid.readline().split()
                sline = fgrid.readline().split()
                #x[i] = sline[0]
                #y[i] = sline[1]
                #z[i] = sline[2]
                nodesi[i0] = sline[0]
                nodesi[i0 + 1] = sline[1]
                nodesi[i0 + 2] = sline[2]
                i0 += 3

            # we want a contiguous array
            nodes = nodesi.reshape((3, nnodes)).T.ravel().reshape(nnodes, 3)

            ntri_lines = ntris // 2
            j = 0
            for i in range(ntri_lines):
                sline = fgrid.readline().split()
                tris[j, :] = sline[:3]
                tris[j + 1, :] = sline[3:]
                j += 2

            #k = 0
            if ntris % 2: # if there's a leftover line
                tris[j] = sline[:3]
                #bcs[:3] = sline[3:]
                #k = 3

            nquestion = ntris // 3 // 2
            for i in range(nquestion):
                sline = fgrid.readline().split()

            #sline = fgrid.readline().split()
            # nnodes=16391
            # ntris=4512
            # ntets = 90892
            # 60595 lines * 6/4 = 90892.5
            #752, 753 = 4512/3/2
            #assert sline == ['15770','15767','16242','15772','15765','15762'], sline

            #ntet_lines = ntets * 6 // 44

            tets = np.zeros(ntets * 4, dtype='int32')
            tets[:] = fgrid.read().split()
            tets = tets.reshape((ntets, 4))

        self.nodes = nodes
        self.tris = tris
        self.tets = tets
