from pyNastran.op2.fortran_format import FortranFormat
from struct import unpack
from numpy import degrees, array, allclose, zeros, hstack, vstack


class ADB_Reader(FortranFormat):
    def __init__(self, log=None, debug=None):
        self.log = log
        self.debug = debug
        self.n = 0
        #self.p_inf =
        self.rho_inf = 0.002377
        self.v_inf = 100.
        self.q_inf = 0.5 * self.rho_inf * self.v_inf ** 2

    @property
    def nnodes(self):
        return self.nodes.shape[0] #+ self.wake_xyz.shape[0]

    def read_adb(self, adb_filename):
        self.f = open(adb_filename, 'rb')
        n = 4
        data = self.f.read(4)
        version_code_small, = unpack('>i', data)
        version_code_large, = unpack('<i', data)
        if version_code_small == -123789456:
            endian = '>'
        elif version_code_large == -123789456:
            endian = '<'
        else:
            print(version_code_small, version_code_large)

        self.n = 4
        n = 40
        data = self.f.read(n)
        nnodes, ntris, nmach, nalpha, sref, cref, bref, xcg, ycg, zcg = unpack(endian + '4i6f', data)
        print('nnodes=%s ntris=%s nmach=%s nalpha=%s' % (nnodes, ntris, nmach, nalpha))
        print('sref=%s cref=%s bref=%s cg=(%s, %s, %s)' % (sref, cref, bref, xcg, ycg, zcg))
        self.n += 40
        assert self.n == self.f.tell(), 'n=%s tell=%s' % (self.n, self.f.tell())

        nbeta = 1

        fmt = '%s%if' % (endian, nmach)
        data = self.f.read(4 * nmach)
        self.machs = unpack(fmt, data)
        self.n += 4 * nmach

        fmt = '%s%if' % (endian, nalpha)
        data = self.f.read(4 * nalpha)
        self.alphas = degrees(unpack(fmt, data))
        self.n += 4 * nalpha
        assert self.n == self.f.tell(), 'n=%s tell=%s' % (self.n, self.f.tell())

        fmt = '%s%if' % (endian, nbeta)
        data = self.f.read(4 * nbeta)
        self.betas = degrees(unpack(fmt, data))
        self.n += 4 * nbeta
        print('machs = %s' % self.machs)
        print('alphas = %s' % self.alphas)
        print('betas = %s' % self.betas)

        data = self.f.read(4)
        self.n += 4
        nwings, = unpack(endian + 'i', data)

        assert self.n == self.f.tell(), 'n=%s tell=%s' % (self.n, self.f.tell())

        fmt = endian + 'i'
        assert nwings < 10, nwings
        print('nwings=%s' % nwings)
        for i in range(nwings):
            data = self.f.read(4)
            j, = unpack(endian + 'i', data)
            print('j = %r' % j)
            assert j < 1000, j
            self.n += 4

            #data = self.f.read(j * 4)
            #fmt = '%s%is' % (endian, j * 4)
            #wing_name = unpack(fmt, data)
            #print('wing_name = %r' % wing_name)

            data = self.f.read(100)
            fmt = '%s%is' % (endian, 100)
            self.n += 100
            wing_name = unpack(fmt, data)
            print('wing_name2 = %r' % wing_name)
            #self.n += j * 4
        assert self.n == self.f.tell(), 'n=%s tell=%s' % (self.n, self.f.tell())

        data = self.f.read(4)
        self.n += 4
        nbodies, = unpack(endian + 'i', data)
        print('nbodies=%s' % nbodies)
        assert nbodies < 10, nbodies
        for i in range(nbodies):
            data = self.f.read(4)
            j, = unpack(endian + 'i', data)
            print('j = %r' % j)
            self.n += 4
            assert j < 1000, j

            data = self.f.read(j * 4)
            fmt = '%s%is' % (endian, j * 4)
            body_name = unpack(fmt, data)
            self.n += j * 4


        self.tris = zeros((ntris, 3), dtype='int32')
        self.surf_id = zeros(ntris, dtype='int32')
        self.area = zeros(ntris, dtype='float32')
        for i in range(ntris):
            data = self.f.read(24)
            n1, n2, n3, surf_type, surf_id, area = unpack(endian + '5if', data)
            self.tris[i, :] = [n1, n2, n3]
            self.surf_id[i] = surf_id
            self.area[i] = area
            self.n += 24
        assert self.n == self.f.tell(), 'n=%s tell=%s' % (self.n, self.f.tell())

        nxyz = nnodes * 3
        data = self.f.read(nxyz * 4)
        fmt = '%s%if' % (endian, nxyz)
        self.nodes = array(unpack(fmt, data), dtype='float32')
        self.nodes = self.nodes.reshape((nnodes, 3))
        self.n += nxyz * 4

        #self.show_ndata(200, types='ifs')
        n = 8
        data = self.f.read(n)
        self.n += n
        #p0, q0 = unpack(endian + '2f', data)
        #print('p0=%s q0=%s qinf=%s' % (p0, q0, self.q_inf))

        fmt = '%s%if' % (endian, ntris)
        data = self.f.read(4 * ntris)
        self.Cp = array(unpack(fmt, data), dtype='float32')
        print('Cp =', self.Cp.max(), self.Cp.min())

        data = self.f.read(4)
        self.n += 4
        nwakes, = unpack(endian + 'i', data)
        X = []
        Y = []
        Z = []
        wakes = []
        iwake_node = nnodes
        for i in range(nwakes):
            data = self.f.read(4)
            j, = unpack(endian + 'i', data)
            print('j = %r' % j)
            self.n += 4
            assert 0 < j < 1000, j
            for ji in range(j):
                data = self.f.read(12)
                x, y, z = unpack(endian + '3f', data)
                self.n += 12
                X.append(x)
                Y.append(y)
                Z.append(z)
            wake = [iwake_node, iwake_node + j]
            iwake_node += j + 1
            wakes.append(wake)
        X = array(X, dtype='float32')
        Y = array(Y, dtype='float32')
        Z = array(Z, dtype='float32')
        self.wake_xyz = vstack([X, Y, Z]).T
        self.wake_elements = wakes
        print('X', X)
        print('Z', Z)
        assert self.wake_xyz.shape[1] == 3, self.wake_xyz.shape
        #print(X.max(), X.min())
        #print(Y.max(), Y.min())
        #print(Z.max(), Z.min())



        n = 64
        data = self.f.read(n)
        self.n += n
        self.show_data(data, types='ifs')
        print('--------')

        if 0:
            npath = 68
            data = self.f.read(npath)
            #self.show_data(data, types='s')
            p = unpack('68s', data)
            print('p = %r' % p)
            print('--------')

            n = 36
            data = self.f.read(n)
            self.show_data(data, types='ifs')

            print('--------')
            data = self.f.read(npath)
            p = unpack('68s', data)
            print('p = %r' % p)
            self.show_data(data, types='s')
            print('----------------')

            self.show_ndata(300, types='ifs')
        return self.nodes, self.tris

def main():
    adb_filename = 'model_DegenGeom.adb'
    a = ADB_Reader()
    a.read_adb(adb_filename)

if __name__ == '__main__':
    main()
#WingGeom\x00tran/pyNastran/OpenVSP-3.1.0-win32/model_DegenGeom.history