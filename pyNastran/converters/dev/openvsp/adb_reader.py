import sys
from struct import unpack
from numpy import degrees, array, zeros, vstack, where
from cpylog import get_logger2


class ADB_Reader:
    def __init__(self, log=None, debug: str | bool | None=None, batch=False):
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
        self._uendian = '<'
        self.f = None

    @property
    def nnodes(self):
        return self.nodes.shape[0] #+ self.wake_xyz.shape[0]

    def read_adb(self, adb_filename):
        with open(adb_filename, 'rb') as adb_file:
            self.f = adb_file
            self._read_adb_file(adb_file)
        return self.nodes, self.tris

    def _read_adb_file(self, adb_file):
        n = 4
        data = adb_file.read(4)
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
        data = adb_file.read(n)
        header = unpack(endian + '4i6f', data)
        nnodes, ntris, nmach, nalpha, sref, cref, bref, xcg, ycg, zcg = header
        print('nnodes=%s ntris=%s nmach=%s nalpha=%s' % (nnodes, ntris, nmach, nalpha))
        print('sref=%s cref=%s bref=%s cg=(%s, %s, %s)' % (sref, cref, bref, xcg, ycg, zcg))
        self.cg = array([xcg, ycg, zcg], dtype='float32')
        self.n += 40
        assert self.n == adb_file.tell(), 'n=%s tell=%s' % (self.n, adb_file.tell())

        nbeta = 1

        fmt = '%s%if' % (endian, nmach)
        data = adb_file.read(4 * nmach)
        self.machs = unpack(fmt, data)
        self.n += 4 * nmach

        fmt = '%s%if' % (endian, nalpha)
        data = adb_file.read(4 * nalpha)
        self.alphas = degrees(unpack(fmt, data))
        self.n += 4 * nalpha
        assert self.n == adb_file.tell(), 'n=%s tell=%s' % (self.n, adb_file.tell())

        fmt = '%s%if' % (endian, nbeta)
        data = adb_file.read(4 * nbeta)
        self.betas = degrees(unpack(fmt, data))
        self.n += 4 * nbeta
        print('machs = %s' % self.machs)
        print('alphas = %s' % self.alphas)
        print('betas = %s' % self.betas)

        data = adb_file.read(4)
        self.n += 4
        nwings, = unpack(endian + 'i', data)

        assert self.n == adb_file.tell(), 'n=%s tell=%s' % (self.n, adb_file.tell())

        fmt = endian + 'i'
        assert nwings < 10, nwings
        print('nwings=%s' % nwings)
        for i in range(nwings):
            data = adb_file.read(4)
            j, = unpack(endian + 'i', data)
            print('j = %r' % j)
            assert j < 1000, j
            self.n += 4

            #data = adb_file.read(j * 4)
            #fmt = '%s%is' % (endian, j * 4)
            #wing_name = unpack(fmt, data)
            #print('wing_name = %r' % wing_name)

            data = adb_file.read(100)
            fmt = '%s%is' % (endian, 100)
            self.n += 100
            wing_name = unpack(fmt, data)
            print('wing_name2 = %r' % wing_name)
            #self.n += j * 4
        assert self.n == adb_file.tell(), 'n=%s tell=%s' % (self.n, adb_file.tell())

        data = adb_file.read(4)
        self.n += 4
        nbodies, = unpack(endian + 'i', data)
        print('nbodies=%s' % nbodies)
        assert nbodies < 10, nbodies
        for i in range(nbodies):
            data = adb_file.read(4)
            j, = unpack(endian + 'i', data)
            print('j = %r' % j)
            self.n += 4
            assert j < 1000, j

            data = adb_file.read(j * 4)
            fmt = '%s%is' % (endian, j * 4)
            unused_body_name = unpack(fmt, data)
            self.n += j * 4


        self.tris = zeros((ntris, 3), dtype='int32')
        self.surf_id = zeros(ntris, dtype='int32')
        self.area = zeros(ntris, dtype='float32')
        for i in range(ntris):
            data = adb_file.read(24)
            n1, n2, n3, unused_surf_type, surf_id, area = unpack(endian + '5if', data)
            self.tris[i, :] = [n1, n2, n3]
            self.surf_id[i] = surf_id
            self.area[i] = area
            self.n += 24
        assert self.n == adb_file.tell(), 'n=%s tell=%s' % (self.n, adb_file.tell())

        nxyz = nnodes * 3
        data = adb_file.read(nxyz * 4)
        fmt = '%s%if' % (endian, nxyz)
        self.nodes = array(unpack(fmt, data), dtype='float32')
        self.nodes = self.nodes.reshape((nnodes, 3))
        self.n += nxyz * 4

        #self.show_ndata(200, types='ifs')
        n = 8
        data = adb_file.read(n)
        self.n += n
        p0, q0 = unpack(endian + '2f', data)
        print('p0=%s q0=%s qinf=%s' % (p0, q0, self.q_inf))
        assert p0 == 1000., p0 # this is hardcoded in the VSP solver
        assert q0 == 100000., q0 # this is hardcoded in the VSP solver

        fmt = '%s%if' % (endian, ntris)
        data = adb_file.read(4 * ntris)
        self.Cp = array(unpack(fmt, data), dtype='float32')
        print('Cp =', self.Cp.max(), self.Cp.min())

        data = adb_file.read(4)
        self.n += 4
        nwakes, = unpack(endian + 'i', data)
        #nwakes -= 2
        X = []
        Y = []
        Z = []
        wakes = []
        iwake_node = nnodes
        self.log.debug('nwakes = %s' % nwakes)
        for unused_iwake in range(nwakes):
            data = adb_file.read(4)
            nsub_vortex_nodes, = unpack(endian + 'i', data)
            print('nsub_vortex_nodes = %r' % nsub_vortex_nodes)
            self.n += 4
            assert 0 < nsub_vortex_nodes < 1000, nsub_vortex_nodes
            for unused_ji in range(nsub_vortex_nodes - 1):
                data = adb_file.read(12)
                x, y, z = unpack(endian + '3f', data)
                self.n += 12
                X.append(x)
                Y.append(y)
                Z.append(z)

            if 1:
                # there is a trailing infinity node that I don't want to save
                data = adb_file.read(12)
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
        data = adb_file.read(n)
        self.n += n
        npropulsion_elements, = unpack(endian + 'i', data)
        if npropulsion_elements != 0:
            raise NotImplementedError(npropulsion_elements)

        n = 4
        data = adb_file.read(n)
        self.n += n
        edges = []
        edge_surf_ids = []
        nmesh_levels, = unpack(endian + 'i', data)
        if nmesh_levels != 0:
            data = adb_file.read(4)
            self.n += 4
            nedges, = unpack(endian + 'i', data)
            for unused_iedge in range(nedges):
                data = adb_file.read(12)
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
        data = adb_file.read(n)
        self.n += n
        self.show_data(data, types='ifs')
        self.show_ndata(200, types='ifs')
        print('--------')

        if 0:
            npath = 68
            data = adb_file.read(npath)
            #self.show_data(data, types='s')
            p = unpack('68s', data)
            print('p = %r' % p)
            print('--------')

            n = 36
            data = adb_file.read(n)
            self.show_data(data, types='ifs')

            print('--------')
            data = adb_file.read(npath)
            p = unpack('68s', data)
            print('p = %r' % p)
            self.show_data(data, types='s')
            print('----------------')

            self.show_ndata(300, types='ifs')
        return

    def show_data(self, data, types='ifs', endian=None):  # pragma: no cover
        """
        Shows a data block as various types

        Parameters
        ----------
        data : bytes
            the binary string bytes
        types : str; default='ifs'
            i - int
            f - float
            s - string
            d - double (float64; 8 bytes)
            q - long long (int64; 8 bytes)

            l - long (int; 4 bytes)
            I - unsigned int (int; 4 bytes)
            L - unsigned long (int; 4 bytes)
            Q - unsigned long long (int; 8 bytes)
        endian : str; default=None -> auto determined somewhere else in the code
            the big/little endian {>, <}

        .. warning:: 's' is apparently not Python 3 friendly

        """
        return self._write_data(sys.stdout, data, types=types, endian=endian)

    def _write_data(self, file_obj, data, types='ifs', endian=None):  # pragma: no cover
        """
        Useful function for seeing what's going on locally when debugging.

        Parameters
        ----------
        data : bytes
            the binary string bytes
        types : str; default='ifs'
            i - int
            f - float
            s - string
            d - double (float64; 8 bytes)
            q - long long (int64; 8 bytes)

            l - long (int; 4 bytes)
            I - unsigned int (int; 4 bytes)
            L - unsigned long (int; 4 bytes)
            Q - unsigned long long (int; 8 bytes)
        endian : str; default=None -> auto determined somewhere else in the code
            the big/little endian {>, <}

        """
        n = len(data)
        nints = n // 4
        ndoubles = n // 8
        strings = None
        ints = None
        floats = None
        longs = None

        if endian is None:
            endian = self._uendian
            assert endian is not None, endian

        file_obj.write('\nndata = %s:\n' % n)
        for typei in types:
            assert typei in 'sifdq lIL', 'type=%r is invalid' % typei

        if 's' in types:
            strings = unpack('%s%is' % (endian, n), data)
            file_obj.write("  strings = %s\n" % str(strings))
        if 'i' in types:
            ints = unpack('%s%ii' % (endian, nints), data)
            file_obj.write("  ints    = %s\n" % str(ints))
        if 'f' in types:
            floats = unpack('%s%if' % (endian, nints), data)
            file_obj.write("  floats  = %s\n" % str(floats))
        if 'd' in types:
            doubles = unpack('%s%id' % (endian, ndoubles), data[:ndoubles*8])
            file_obj.write("  doubles (float64) = %s\n" % str(doubles))

        if 'l' in types:
            longs = unpack('%s%il' % (endian, nints), data)
            file_obj.write("  long  = %s\n" % str(longs))
        if 'I' in types:
            ints2 = unpack('%s%iI' % (endian, nints), data)
            file_obj.write("  unsigned int = %s\n" % str(ints2))
        if 'L' in types:
            longs2 = unpack('%s%iL' % (endian, nints), data)
            file_obj.write("  unsigned long = %s\n" % str(longs2))
        if 'q' in types:
            longs = unpack('%s%iq' % (endian, ndoubles), data[:ndoubles*8])
            file_obj.write("  long long (int64) = %s\n" % str(longs))
        file_obj.write('\n')
        return strings, ints, floats

    def show_ndata(self, n, types='ifs'):  # pragma: no cover
        return self._write_ndata(sys.stdout, n, types=types)

    def _write_ndata(self, file_obj, n, types='ifs'):  # pragma: no cover
        """Useful function for seeing what's going on locally when debugging."""
        nold = self.n
        data = file_obj.read(n)
        self.n = nold
        file_obj.seek(self.n)
        return self._write_data(file_obj, data, types=types)

def main():  # pragma: no cover
    adb_filename = 'model_DegenGeom.adb'
    model = ADB_Reader(batch=True)
    model.read_adb(adb_filename)

if __name__ == '__main__':  # pragma: no cover
    main()
#WingGeom\x00tran/pyNastran/OpenVSP-3.1.0-win32/model_DegenGeom.history
