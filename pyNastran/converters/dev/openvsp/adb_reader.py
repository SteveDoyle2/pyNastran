from __future__ import print_function
from struct import unpack
from numpy import degrees, array, zeros, vstack, where
from pyNastran.op2.fortran_format import FortranFormat
from pyNastran.utils.log import get_logger2


class ADB_Reader(FortranFormat):
    def __init__(self, log=None, debug=None, batch=False):
        self.log = get_logger2(log, debug)
        self.debug = debug
        self.n = 0
        #self.p_inf =
        self.rho_inf = 0.002377
        self.v_inf = 100.
        self.q_inf = 0.5 * self.rho_inf * self.v_inf ** 2
        self.batch = batch

        self.cg = None
        self.machs = None
        self.alphas = None
        self.betas = None

        self.nodes = None
        self.Cp = None
        self.tris = None
        self.surf_id = None

        self.edges = None
        self.edge_surf_id = None
        self.area = None
        self.wake_xyz = None
        self.wake_elements = None

    @property
    def nnodes(self):
        return self.nodes.shape[0] #+ self.wake_xyz.shape[0]

    def read_adb(self, adb_filename):
        with open(adb_filename, 'rb') as adb_file:
            self.f = adb_file
            n = 4
            data = self.f.read(4)
            version_code_big, = unpack('>i', data)
            version_code_little, = unpack('<i', data)
            if version_code_big == -123789456:
                endian = '>'
            elif version_code_little == -123789456:
                endian = '<'
            else:
                print(version_code_little, version_code_big)
                raise RuntimeError('incorrect endian')

            self.n = 4
            n = 40
            data = self.f.read(n)
            header = unpack(endian + '4i6f', data)
            nnodes, ntris, nmach, nalpha, sref, cref, bref, xcg, ycg, zcg = header
            print('nnodes=%s ntris=%s nmach=%s nalpha=%s' % (nnodes, ntris, nmach, nalpha))
            print('sref=%s cref=%s bref=%s cg=(%s, %s, %s)' % (sref, cref, bref, xcg, ycg, zcg))
            self.cg = array([xcg, ycg, zcg], dtype='float32')
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
            p0, q0 = unpack(endian + '2f', data)
            print('p0=%s q0=%s qinf=%s' % (p0, q0, self.q_inf))
            assert p0 == 1000., p0 # this is hardcoded in the VSP solver
            assert q0 == 100000., q0 # this is hardcoded in the VSP solver

            fmt = '%s%if' % (endian, ntris)
            data = self.f.read(4 * ntris)
            self.Cp = array(unpack(fmt, data), dtype='float32')
            print('Cp =', self.Cp.max(), self.Cp.min())

            data = self.f.read(4)
            self.n += 4
            nwakes, = unpack(endian + 'i', data)
            #nwakes -= 2
            X = []
            Y = []
            Z = []
            wakes = []
            iwake_node = nnodes
            self.log.debug('nwakes = %s' % nwakes)
            for iwake in range(nwakes):
                data = self.f.read(4)
                nsub_vortex_nodes, = unpack(endian + 'i', data)
                print('nsub_vortex_nodes = %r' % nsub_vortex_nodes)
                self.n += 4
                assert 0 < nsub_vortex_nodes < 1000, nsub_vortex_nodes
                for ji in range(nsub_vortex_nodes - 1):
                    data = self.f.read(12)
                    x, y, z = unpack(endian + '3f', data)
                    self.n += 12
                    X.append(x)
                    Y.append(y)
                    Z.append(z)

                if 1:
                    # there is a trailing infinity node that I don't want to save
                    data = self.f.read(12)
                    x, y, z = unpack(endian + '3f', data)
                    self.n += 12
                    nsub_vortex_nodes -= 1
                #wake = [iwake_node - 1, iwake_node + nsub_vortex_nodes - 1]
                wake = [iwake_node, iwake_node + nsub_vortex_nodes]
                iwake_node += nsub_vortex_nodes #+ 1
                wakes.append(wake)
            assert iwake_node - nnodes == len(X), 'inode=%s nx=%s' % (iwake_node, len(X))
            X = array(X, dtype='float32')
            Y = array(Y, dtype='float32')
            Z = array(Z, dtype='float32')
            ix = where(X > 1000.)[0]
            iy = where(Y > 1000.)[0]
            iz = where(Z > 1000.)[0]
            if self.batch:
                print('ix=%s' % ix)
                print('iy=%s' % iy)
                print('iz=%s' % iz)
            else:
                X[ix] = 0.
                Y[iy] = 0.
                Z[iz] = 0.
            self.wake_xyz = vstack([X, Y, Z]).T
            print('X.max=%s X.min=%s' % (X.max(), X.min()))
            print('Y.max=%s Y.min=%s' % (Y.max(), Y.min()))
            print('Z.max=%s Z.min=%s' % (Z.max(), Z.min()))
            self.wake_elements = wakes
            #print('X', X)
            #print('Z', Z)
            assert self.wake_xyz.shape[1] == 3, self.wake_xyz.shape

            # shift to the reference point
            self.nodes -= self.cg
            self.wake_xyz -= self.cg
            #print(X.max(), X.min())
            #print(Y.max(), Y.min())
            #print(Z.max(), Z.min())

            n = 4
            data = self.f.read(n)
            self.n += n
            npropulsion_elements, = unpack(endian + 'i', data)
            if npropulsion_elements != 0:
                raise NotImplementedError(npropulsion_elements)

            n = 4
            data = self.f.read(n)
            self.n += n
            edges = []
            edge_surf_ids = []
            nmesh_levels, = unpack(endian + 'i', data)
            if nmesh_levels != 0:
                data = self.f.read(4)
                self.n += 4
                nedges, = unpack(endian + 'i', data)
                for iedge in range(nedges):
                    data = self.f.read(12)
                    self.n += 12
                    surf_idi, n1, n2 = unpack(endian + '3i', data)
                    edge = [n1, n2]
                    edges.append(edge)
                    edge_surf_ids.append(surf_idi)
            print(edge_surf_ids)
            print(edges)

            self.edges = array(edges, dtype='int32')
            self.edge_surf_id = array(edge_surf_ids, dtype='int32')

            n = 64
            data = self.f.read(n)
            self.n += n
            self.show_data(data, types='ifs')
            self.show_ndata(200, types='ifs')
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

def main():  # pragma: no cover
    adb_filename = 'model_DegenGeom.adb'
    a = ADB_Reader(batch=True)
    a.read_adb(adb_filename)

if __name__ == '__main__':  # pragma: no cover
    main()
#WingGeom\x00tran/pyNastran/OpenVSP-3.1.0-win32/model_DegenGeom.history
