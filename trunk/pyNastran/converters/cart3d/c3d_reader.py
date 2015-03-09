from struct import unpack, Struct

from numpy import zeros

from pyNastran.op2.fortran_format import FortranFormat

from pyNastran.utils import is_binary
from pyNastran.utils.log import get_logger


class C3D_Reader(FortranFormat):
    def __init__(self, log=None, debug=False):
        self.log = get_logger(log, 'debug' if debug else 'info')
        FortranFormat.__init__(self)

    def read_c3d(self, c3d_filename):
        self.f = open(c3d_filename, 'r')
        data = self.f.read(32)
        nVolHexes, nCutHexes, nSplitCells, nFlowFaces, facesX, facesY, facesZ, nCutFaces = unpack(b'8i', data)
        print('nVolHexes=%s nCutHexes=%s nSplitCells=%s nFlowFaces=%s' % (nVolHexes, nCutHexes, nSplitCells, nFlowFaces))
        print('facesX=%s facesY=%s facesZ=%s nCutFaces=%s' % (facesX, facesY, facesZ, nCutFaces))

        #facesXYZ
        #nCutFaces
        self.n += 32

        fmt = 'Iccc'
        s = Struct(fmt)
        n = 10
        vol_hexs = zeros((n, 3), dtype='int32')
        for j in range(nVolHexes):
            data = self.f.read(7)
            self.n += 7
            n, x, y, z = s.unpack(data)
            #print('n=%s x=%s y=%s z=%s' % (n, x, y, z))
        self.show(50, types='iflILq')

        self.f.close()

def main():
    c3d_filename = 'Mesh.R.c3d'
    c3d = C3D_Reader()
    c3d.read_c3d(c3d_filename)

if __name__ == '__main__':
    main()