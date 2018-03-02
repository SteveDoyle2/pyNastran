from six import iteritems, integer_types
from numpy import zeros, empty

from pyNastran.op2.result_objects.op2_objects import ScalarObject
from pyNastran.f06.f06_formatting import (
    write_floats_13e, write_floats_10e, get_key0,
    _eigenvalue_header)


class GridPointStressesArray(ScalarObject):
    """
    '                                  S T R E S S E S   A T   G R I D   P O I N T S   - -     S U R F A C E       5\n',
    '0                       SURFACE X-AXIS X  NORMAL(Z-AXIS)  Z         REFERENCE COORDINATE SYSTEM FOR SURFACE DEFINITION CID        0\n',
    '     GRID      ELEMENT            STRESSES IN SURFACE SYSTEM           PRINCIPAL STRESSES            MAX             \n',
    '     ID          ID    FIBRE   NORMAL-X   NORMAL-Y   SHEAR-XY     ANGLE      MAJOR      MINOR      SHEAR     VON MISES\n']
    '0     13683          3736    TRIAX6         4.996584E+00   0.0            1.203093E+02   0.0            0.0            0.0'
    '      13683          3737    TRIAX6        -4.996584E+00   0.0           -1.203093E+02   0.0            0.0            0.0'
    '      13683                  *TOTALS*       6.366463E-12   0.0           -1.364242E-12   0.0            0.0            0.0'
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase, apply_data_code=True)
        self.ntotal = 0
        self.ntimes = 0
        self.nelements = 0
        self.itotal = 0
        self.ielement = 0
        self.data = None
        self.itime = None
        self.node_element = None

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def build(self):
        """sizes the vectorized attributes of the GridPointStressesArray"""
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        self.itime = 0
        self.ielement = 0
        self.itotal = 0

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        self.nelements //= self.ntimes

        self.node_element = zeros((self.ntotal, 2), dtype='int32')
        #oxx, oyy, txy, angle, major, minor, ovm
        self.data = zeros((self.ntimes, self.ntotal, 8), dtype='float32')
        self.location = empty(self.ntotal, dtype='U8')
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'

        self._times = zeros(self.ntimes, dtype=dtype)
        self.is_built = True

    #def build_dataframe(self):
        #headers = self.get_headers()
        #element_node = [self.element_node[:, 0], self.element_node[:, 1]]
        #if self.nonlinear_factor is not None:
            #column_names, column_values = self._build_dataframe_transient_header()
            #self.data_frame = pd.Panel(self.data, items=column_values, major_axis=element_node, minor_axis=headers).to_frame()
            #self.data_frame.columns.names = column_names
        #else:
            #self.data_frame = pd.Panel(self.data, major_axis=element_node, minor_axis=headers).to_frame()
            #self.data_frame.columns.names = ['Static']
        #self.data_frame.index.names = ['NodeID', 'ElementID', 'Item']

    def add_sort1(self, dt, nid, eid, fiber, nx, ny, txy, angle, majorP, minorP, tmax, ovm):
        """unvectorized method for adding SORT1 transient data"""
        self._times[self.itime] = dt
        self.node_element[self.itotal, :] = [nid, eid]
        self.location[self.itotal] = fiber
        self.data[self.itime, self.itotal, :] = [nx, ny, txy, angle, majorP, minorP, tmax, ovm]
        self.itotal += 1

    def get_stats(self, short=False):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        ntimes, nelements, _ = self.data.shape
        assert self.ntimes == ntimes, 'ntimes=%s expected=%s' % (self.ntimes, ntimes)
        assert self.nelements == nelements, 'nelements=%s expected=%s' % (self.nelements, nelements)

        msg = []
        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i\n'
                       % (self.__class__.__name__, ntimes, nelements))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i\n'
                       % (self.__class__.__name__, nelements))
            ntimes_word = '1'
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  node_element.shape = %s\n' % str(self.node_element.shape).replace('L', ''))
        msg.append('  location.shape = %s\n' % str(self.location.shape).replace('L', ''))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg += self.get_data_code()

        #msg = self.get_data_code()
        #if self.nonlinear_factor is not None:  # transient
            #ntimes = len(self.nx)
            #times0 = get_key0(self.nx)
            #nelements = len(self. nx[times0])
            #msg.append('  type=%s ntimes=%s nelements=%s\n'
                       #% (self.__class__.__name__, ntimes, nelements))
            #ntimes_word = 'ntimes'
        #else:
            #nelements = len(self.nx)
            #msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     #nelements))
            #ntimes_word = '1'

        #headers = self.get_headers()
        #n = len(headers)
        #msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (ntimes_word, n, n,
                                                                 #', '.join(headers)))
        #msg.append('  grid_element.shape=%s\n' % str(self.grid_element.shape).replace('L', ''))
        #msg.append('  data.shape=%s\n' % str(self.data.shape).replace('L', ''))
        return msg

    def get_headers(self):
        headers = ['nx', 'ny', 'txy', 'angle', 'majorP', 'minorP', 'tmax', 'ovm']
        return headers

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []

        i = self.ID
        cid = self.refid
        axis_int = self.oCoord
        axis_map = {0 : 'X', 1 : 'Y', 2 : 'Z'}
        axis = axis_map[axis_int]
        msg = [
            '                                  S T R E S S E S   A T   G R I D   P O I N T S   - -     S U R F A C E     %s\n' % i,
            '0                       SURFACE X-AXIS X  NORMAL(Z-AXIS)  %s         REFERENCE COORDINATE SYSTEM FOR SURFACE DEFINITION CID        %s\n' % (axis, cid),
            '     GRID      ELEMENT            STRESSES IN SURFACE SYSTEM           PRINCIPAL STRESSES            MAX             \n',
            '     ID          ID    FIBRE   NORMAL-X   NORMAL-Y   SHEAR-XY     ANGLE      MAJOR      MINOR      SHEAR     VON MISES\n']
           #'0     13683          3736    TRIAX6         4.996584E+00   0.0            1.203093E+02   0.0            0.0            0.0'
           #'      13683          3737    TRIAX6        -4.996584E+00   0.0           -1.203093E+02   0.0            0.0            0.0'
           #'      13683                  *TOTALS*       6.366463E-12   0.0           -1.364242E-12   0.0            0.0            0.0'

        ntimes = self.data.shape[0]

        nids = self.node_element[:, 0]
        eids = self.node_element[:, 1]
        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg))

            nx = self.data[itime, :, 0]
            ny = self.data[itime, :, 1]
            txy = self.data[itime, :, 2]
            angle = self.data[itime, :, 3]
            majorp = self.data[itime, :, 4]
            minorp = self.data[itime, :, 5]
            tmax = self.data[itime, :, 6]
            ovm = self.data[itime, :, 7]
            fibers = self.location
            nid_old = -1
            for (nid, eid, fiber, nxi, nyi, txyi, anglei, majorpi, minorpi, tmaxi, ovmi) in zip(
                nids, eids, fibers, nx, ny, txy, angle, majorp, minorp, tmax, ovm):
                [nxi, nyi, txyi, majorpi, minorpi, tmaxi, ovmi] = write_floats_10e([
                    nxi, nyi, txyi, majorpi, minorpi, tmaxi, ovmi])
                if nid > nid_old:
                    f06_file.write(
                        '0%8s  %8s   %4s    %s %s %s  %8.4f %10s %10s %10s  %s\n' % (
                            nid, eid, fiber, nxi, nyi, txyi, anglei, majorpi, minorpi,
                            tmaxi, ovmi))
                else:
                    f06_file.write(
                        ' %8s  %8s   %4s    %s %s %s  %8.4f %10s %10s %10s  %s\n' % (
                            '', '', fiber, nxi, nyi, txyi, anglei, majorpi, minorpi,
                            tmaxi, ovmi))
                nid_old = nid
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

        #for ekey, nxs in sorted(iteritems(self.nx)):
            #ekey2 = ekey
            #zero = '0'
            #for iLoad, nx in enumerate(nxs):
                #ny = self.ny[ekey][iLoad]
                #txy = self.txy[ekey][iLoad]
                #angle = self.angle[ekey][iLoad]
                #majorP = self.majorP[ekey][iLoad]
                #minorP = self.minorP[ekey][iLoad]
                #tmax = self.tmax[ekey][iLoad]
                #ovm = self.ovm[ekey][iLoad]

                #elem_name = self.elemName[ekey][iLoad]
                #eid = self.eids[ekey][iLoad]
                #vals = [nx, ny, txy, majorP, minorP, tmax, ovm]
                #vals2 = write_floats_10e(vals)
                #[nx, ny, txy, majorP, minorP, tmax, ovm] = vals2
                #if eid == 0:
                    #eid = zero
                #angle = write_floats_8p4f([angle])
                #anglei = angle[0]
                #msg.append('%s%8s  %8s   %4s    %s %s %s   %8s %10s %10s %10s  %s\n' % (
                    #zero, ekey2, eid, elem_name, nx, ny, txy, anglei, majorP, minorP, tmax, ovm))
                #zero = ' '
                #ekey2 = ' '
        #msg.append(page_stamp % page_num)
        #f.write(''.join(msg))
        #return page_num


class GridPointStressesVolume(ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
        self.nx = {}
        self.ny = {}
        self.nz = {}
        self.txy = {}
        self.tyz = {}
        self.txz = {}
        self.pressure = {}
        self.ovm = {}

        self.elemName = {}
        self.eids = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.add_sort2

    def get_stats(self, short=False):
        msg = self.get_data_code()
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.nx)
            times0 = get_key0(self.nx)
            nelements = len(self. nx[times0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self. nx)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  nx, ny, nz, txy, tyz, txz, pressure, ovm\n')
        return msg

    def add_new_transient(self, dt):  # ekey
        """initializes the transient variables"""
        self.nx[dt] = {}
        self.ny[dt] = {}
        self.nz[dt] = {}
        self.txy[dt] = {}
        self.tyz[dt] = {}
        self.txz[dt] = {}
        self.pressure[dt] = {}
        self.ovm[dt] = {}

        self.elemName = {}
        self.eids = {}

    def add(self, dt, ekey, nx, ny, nz, txy, tyz, txz, pressure, ovm):
        if ekey not in self.nx:
            #self.eids[ekey] = []
            #self.elemName[ekey] = []
            self.nx[ekey] = []
            self.ny[ekey] = []
            self.nz[ekey] = []
            self.txy[ekey] = []
            self.tyz[ekey] = []
            self.txz[ekey] = []
            self.pressure[ekey] = []
            self.ovm[ekey] = []
        self.nx[ekey].append(nx)
        self.ny[ekey].append(ny)
        self.nz[ekey].append(nz)
        self.txy[ekey].append(txy)
        self.tyz[ekey].append(tyz)
        self.txz[ekey].append(txz)
        self.pressure[ekey].append(pressure)
        self.ovm[ekey].append(ovm)

        #self.elemName[ekey].append(elemName)
        #self.eids[ekey].append(eid)

    def add_sort1(self, dt, ekey, nx, ny, nz, txy, tyz, txz, pressure, ovm):
        """unvectorized method for adding SORT1 transient data"""
        if dt not in self.nx:
            self.add_new_transient(dt)

        #print "%s=%s ekey=%s eid=%s elemName=%s f1=%s" %(self.data_code['name'],dt,ekey,eid,elemName,f1)
        if ekey not in self.nx[dt]:
            #self.eids[ekey] = []
            #self.elemName[ekey] = []
            self.nx[dt][ekey] = []
            self.ny[dt][ekey] = []
            self.nz[dt][ekey] = []
            self.txy[dt][ekey] = []
            self.tyz[dt][ekey] = []
            self.txz[dt][ekey] = []
            self.pressure[ekey] = []
            self.ovm[dt][ekey] = []
        self.eids[ekey].append(eid)
        #self.elemName[ekey].append(elemName)

        self.nx[dt][ekey].append(nx)
        self.ny[dt][ekey].append(ny)
        self.nz[dt][ekey].append(nz)
        self.txy[dt][ekey].append(txy)
        self.tyz[dt][ekey].append(tyz)
        self.txz[dt][ekey].append(txz)
        self.pressure[dt][ekey].append(pressure)
        self.ovm[dt][ekey].append(ovm)

    def delete_transient(self, dt):
        del self.nx[dt]
        del self.ny[dt]
        del self.nz[dt]
        del self.txy[dt]
        del self.tyz[dt]
        del self.txz[dt]
        del self.pressure[dt]
        del self.ovm[dt]

    def get_transients(self):
        k = self.nx.keys()
        k.sort()
        return k

    #def cleanupObj(self):
        #k = self.elemName.keys()
        #self.elemName = self.elemName[k[0]]
        #self.eids = self.eids[k[0]]

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        f06_file.write('GridPointStressesVolume write_f06 is not implemented...\n')
        return page_num

        #raise NotImplementedError()
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f06_file, is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        msg = header + ['                                  S T R E S S E S   A T   G R I D   P O I N T S   - -     S U R F A C E       5\n',
                        '0                       SURFACE X-AXIS X  NORMAL(Z-AXIS)  Z         REFERENCE COORDINATE SYSTEM FOR SURFACE DEFINITION CID        0\n',
                        '     GRID      ELEMENT            STRESSES IN SURFACE SYSTEM           PRINCIPAL STRESSES            MAX             \n',
                        '     ID          ID    FIBRE   NORMAL-X   NORMAL-Y   SHEAR-XY     ANGLE      MAJOR      MINOR      SHEAR     VON MISES\n']
              #'0     13683          3736    TRIAX6         4.996584E+00   0.0            1.203093E+02   0.0            0.0            0.0'
              #'      13683          3737    TRIAX6        -4.996584E+00   0.0           -1.203093E+02   0.0            0.0            0.0'
              #'      13683                  *TOTALS*       6.366463E-12   0.0           -1.364242E-12   0.0            0.0            0.0'
        for ekey, nxs in sorted(iteritems(self.nx)):
            #ekey2 = ekey
            zero = '0'
            for iLoad, nx in enumerate(nxs):
                ny = self.ny[ekey][iLoad]
                nz = self.nz[ekey][iLoad]
                txy = self.txy[ekey][iLoad]
                tyz = self.tyz[ekey][iLoad]
                txz = self.txz[ekey][iLoad]
                pressure = self.pressure[ekey][iLoad]
                ovm = self.ovm[ekey][iLoad]

                #(elemName) = self.elemName[ekey][iLoad]
                #eid = self.eids[ekey][iLoad]
                vals = [nx, ny, nz, txy, tyz, txz, pressure, ovm]
                vals2 = write_floats_10e(vals)
                [nx, ny, nz, txy, tyz, txz, pressure, ovm] = vals2
                f06_file.write('%s%8s  %s %s %s   %s %s %s %s  %-s\n' % (
                    zero, ekey, nx, ny, nz, txy, tyz, txz, pressure, ovm.rstrip()))
                zero = ' '
                #ekey2 = ' '

        f06_file.write(page_stamp % page_num)
        return page_num

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        f.write('GridPointStressesVolume _write_f06_transient is not implemented...\n')
        return page_num
        #raise NotImplementedError()
        msg = header + ['                                  S T R E S S E S   A T   G R I D   P O I N T S   - -     S U R F A C E       5\n',
                        '0                       SURFACE X-AXIS X  NORMAL(Z-AXIS)  Z         REFERENCE COORDINATE SYSTEM FOR SURFACE DEFINITION CID        0\n',
                        '     GRID      ELEMENT            STRESSES IN SURFACE SYSTEM           PRINCIPAL STRESSES            MAX             \n',
                        '     ID          ID    FIBRE   NORMAL-X   NORMAL-Y   SHEAR-XY     ANGLE      MAJOR      MINOR      SHEAR     VON MISES\n']
              #'0     13683          3736    TRIAX6         4.996584E+00   0.0            1.203093E+02   0.0            0.0            0.0'
              #'      13683          3737    TRIAX6        -4.996584E+00   0.0           -1.203093E+02   0.0            0.0            0.0'
              #'      13683                  *TOTALS*       6.366463E-12   0.0           -1.364242E-12   0.0            0.0            0.0'
        for dt, Forces in sorted(iteritems(self.forces)):
            for ekey, force in sorted(iteritems(Forces)):
                zero = '0'
                for iLoad, f in enumerate(force):
                    (f1, f2, f3) = f
                    (m1, m2, m3) = self.moments[dt][ekey][iLoad]
                    elem_name = self.elemName[ekey][iLoad]
                    eid = self.eids[ekey][iLoad]

                    vals = [f1, f2, f3, m1, m2, m3]
                    vals2 = write_floats_13e(vals)
                    [f1, f2, f3, m1, m2, m3] = vals2
                    if eid == 0:
                        eid = ''

                    msg.append('%s  %8s    %10s    %8s      %s  %s  %s  %s  %s  %-s\n' % (zero, ekey, eid, elem_name, f1, f2, f3, m1, m2, m3))
                    zero = ' '

            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            msg = ['']
            page_num += 1
        return page_num - 1
