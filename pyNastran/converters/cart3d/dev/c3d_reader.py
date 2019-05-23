from struct import unpack, Struct

from numpy import zeros
from cpylog import get_logger2

#from pyNastran.utils import is_binary_file


class C3D_Reader:
    def __init__(self, log=None, debug=False):
        self.log = get_logger2(log, debug=debug)

    def read_c3d(self, c3d_filename):
        self.f = open(c3d_filename, 'r')
        data = self.f.read(32)
        nvol_hexes, ncut_hexes, nsplit_cells, nflow_faces, faces_x, faces_y, faces_z, ncut_faces = unpack(b'8i', data)
        print('nvol_hexes=%s ncut_hexes=%s nsplit_cells=%s nflow_faces=%s' % (
            nvol_hexes, ncut_hexes, nsplit_cells, nflow_faces))
        print('faces_x=%s faces_y=%s faces_z=%s ncut_faces=%s' % (
            faces_x, faces_y, faces_z, ncut_faces))

        #facesXYZ
        #ncut_faces
        self.n += 32

        fmt = 'Iccc'
        s = Struct(fmt)
        n = 10
        vol_hexs = zeros((n, 3), dtype='int32')
        for j in range(nvol_hexes):
            data = self.f.read(7)
            self.n += 7
            n, x, y, z = s.unpack(data)
            #print('n=%s x=%s y=%s z=%s' % (n, x, y, z))
        self.show(50, types='iflILq')

        self.f.close()

def main():  # pragma: no cover
    c3d_filename = 'Mesh.R.c3d'
    c3d = C3D_Reader()
    c3d.read_c3d(c3d_filename)

if __name__ == '__main__':  # pragma: no cover
    main()
