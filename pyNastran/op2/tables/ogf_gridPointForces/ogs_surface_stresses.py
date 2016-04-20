from six import iteritems
from numpy import zeros

from pyNastran.op2.result_objects.op2_objects import ScalarObject
from pyNastran.f06.f06_formatting import write_floats_13e, writeFloats10E, writeFloats8p4F, get_key0


class GridPointStressesArray(ScalarObject):
    """
        msg = header + ['                                  S T R E S S E S   A T   G R I D   P O I N T S   - -     S U R F A C E       5\n',
                        '0                       SURFACE X-AXIS X  NORMAL(Z-AXIS)  Z         REFERENCE COORDINATE SYSTEM FOR SURFACE DEFINITION CID        0\n',
                        '     GRID      ELEMENT            STRESSES IN SURFACE SYSTEM           PRINCIPAL STRESSES            MAX             \n',
                        '     ID          ID    FIBRE   NORMAL-X   NORMAL-Y   SHEAR-XY     ANGLE      MAJOR      MINOR      SHEAR     VON MISES\n']
              #'0     13683          3736    TRIAX6         4.996584E+00   0.0            1.203093E+02   0.0            0.0            0.0'
              #'      13683          3737    TRIAX6        -4.996584E+00   0.0           -1.203093E+02   0.0            0.0            0.0'
              #'      13683                  *TOTALS*       6.366463E-12   0.0           -1.364242E-12   0.0            0.0            0.0'
    """
    def __init__(self, data_code, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
        self.ntotal = 0
        self.ntimes = 0

    def build(self):
        self.grid_element = zeros((self.ntotal, 2), dtype='int32')
        #oxx, oyy, txy, angle, major, minor, ovm
        self.data = zeros((self.ntimes, self.ntotal, 7), dtype='float32')

    def build_dataframe(self):
        headers = self.get_headers()
        element_node = [self.element_node[:, 0], self.element_node[:, 1]]
        if self.nonlinear_factor is not None:
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=element_node, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
        else:
            self.data_frame = pd.Panel(self.data, major_axis=element_node, minor_axis=headers).to_frame()
            self.data_frame.columns.names = ['Static']
        self.data_frame.index.names = ['NodeID', 'ElementID', 'Item']

    def add_sort1(self, dt, ekey, eid, elemName, nx, ny, txy, angle, majorP, minorP, tmax, ovm):
        self.times[self.itime] = dt
        self.grid_element[self.ntotal, :] = [ekey, eid]
        self.data[self.itime, self.ntotal, :] = [nx, ny, txy, angle, majorP, minorP, tmax, ovm]

    def get_stats(self):
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
        msg.append('  nx, ny, txy, angle, majorP, minorP, tmax, ovm\n')
        return msg


class GridPointStresses(ScalarObject):

    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
        self.nx = {}
        self.ny = {}
        self.txy = {}
        self.angle = {}
        self.majorP = {}
        self.minorP = {}
        self.tmax = {}
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

    def get_stats(self):
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
        msg.append('  nx, ny, txy, angle, majorP, minorP, tmax, ovm\n')
        return msg

    def add_new_transient(self, dt):  # ekey
        """initializes the transient variables"""
        self.nx[dt] = {}
        self.ny[dt] = {}
        self.txy[dt] = {}
        self.angle[dt] = {}
        self.majorP[dt] = {}
        self.minorP[dt] = {}
        self.tmax[dt] = {}
        self.ovm[dt] = {}

        self.elemName = {}
        self.eids = {}

    def add(self, dt, ekey, eid, elemName, nx, ny, txy, angle, majorP, minorP, tmax, ovm):
        if ekey not in self.nx:
            self.eids[ekey] = []
            self.elemName[ekey] = []
            self.nx[ekey] = []
            self.ny[ekey] = []
            self.txy[ekey] = []
            self.angle[ekey] = []
            self.majorP[ekey] = []
            self.minorP[ekey] = []
            self.tmax[ekey] = []
            self.ovm[ekey] = []
        self.nx[ekey].append(nx)
        self.ny[ekey].append(ny)
        self.txy[ekey].append(txy)
        self.angle[ekey].append(angle)
        self.majorP[ekey].append(majorP)
        self.minorP[ekey].append(minorP)
        self.tmax[ekey].append(tmax)
        self.ovm[ekey].append(ovm)

        self.elemName[ekey].append(elemName)
        self.eids[ekey].append(eid)

    def add_sort1(self, dt, ekey, eid, elemName, nx, ny, txy, angle, majorP, minorP, tmax, ovm):
        if dt not in self.nx:
            self.add_new_transient(dt)

        #print "%s=%s ekey=%s eid=%s elemName=%s f1=%s" %(self.data_code['name'],dt,ekey,eid,elemName,f1)
        if ekey not in self.nx[dt]:
            self.eids[ekey] = []
            self.elemName[ekey] = []
            self.nx[dt][ekey] = []
            self.ny[dt][ekey] = []
            self.txy[dt][ekey] = []
            self.angle[dt][ekey] = []
            self.majorP[dt][ekey] = []
            self.minorP[dt][ekey] = []
            self.tmax[dt][ekey] = []
            self.ovm[dt][ekey] = []
        self.eids[ekey].append(eid)
        self.elemName[ekey].append(elemName)

        self.nx[dt][ekey].append(nx)
        self.ny[dt][ekey].append(ny)
        self.txy[dt][ekey].append(txy)
        self.angle[dt][ekey].append(angle)
        self.majorP[dt][ekey].append(majorP)
        self.minorP[dt][ekey].append(minorP)
        self.tmax[dt][ekey].append(tmax)
        self.ovm[dt][ekey].append(ovm)

    def delete_transient(self, dt):
        del self.nx[dt]
        del self.ny[dt]
        del self.txy[dt]
        del self.angle[dt]
        del self.majorP[dt]
        del self.minorP[dt]
        del self.tmax[dt]
        del self.ovm[dt]

    def get_transients(self):
        k = self.nx.keys()
        k.sort()
        return k

    #def cleanupObj(self):
        #k = self.elemName.keys()
        #self.elemName = self.elemName[k[0]]
        #self.eids = self.eids[k[0]]

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        msg = header + ['                                  S T R E S S E S   A T   G R I D   P O I N T S   - -     S U R F A C E       5\n',
                        '0                       SURFACE X-AXIS X  NORMAL(Z-AXIS)  Z         REFERENCE COORDINATE SYSTEM FOR SURFACE DEFINITION CID        0\n',
                        '     GRID      ELEMENT            STRESSES IN SURFACE SYSTEM           PRINCIPAL STRESSES            MAX             \n',
                        '     ID          ID    FIBRE   NORMAL-X   NORMAL-Y   SHEAR-XY     ANGLE      MAJOR      MINOR      SHEAR     VON MISES\n']
              #'0     13683          3736    TRIAX6         4.996584E+00   0.0            1.203093E+02   0.0            0.0            0.0'
              #'      13683          3737    TRIAX6        -4.996584E+00   0.0           -1.203093E+02   0.0            0.0            0.0'
              #'      13683                  *TOTALS*       6.366463E-12   0.0           -1.364242E-12   0.0            0.0            0.0'
        for ekey, nxs in sorted(iteritems(self.nx)):
            ekey2 = ekey
            zero = '0'
            for iLoad, nx in enumerate(nxs):
                ny = self.ny[ekey][iLoad]
                txy = self.txy[ekey][iLoad]
                angle = self.angle[ekey][iLoad]
                majorP = self.majorP[ekey][iLoad]
                minorP = self.minorP[ekey][iLoad]
                tmax = self.tmax[ekey][iLoad]
                ovm = self.ovm[ekey][iLoad]

                (elemName) = self.elemName[ekey][iLoad]
                eid = self.eids[ekey][iLoad]
                vals = [nx, ny, txy, majorP, minorP, tmax, ovm]
                (vals2, is_all_zeros) = writeFloats10E(vals)
                [nx, ny, txy, majorP, minorP, tmax, ovm] = vals2
                if eid == 0:
                    eid = zero
                angle, isAllZero = writeFloats8p4F([angle])
                anglei = angle[0]
                msg.append('%s%8s  %8s   %4s    %s %s %s   %8s %10s %10s %10s  %s\n' % (zero, ekey2, eid, elemName, nx, ny, txy, anglei, majorP, minorP, tmax, ovm))
                zero = ' '
                ekey2 = ' '
        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        f.write('GridPointStressesObject write_f06 is not implemented...\n')
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
                    (elemName) = self.elemName[ekey][iLoad]
                    eid = self.eids[ekey][iLoad]

                    vals = [f1, f2, f3, m1, m2, m3]
                    vals2 = write_floats_13e(vals)
                    [f1, f2, f3, m1, m2, m3] = vals2
                    if eid == 0:
                        eid = ''
                    msg.append('%s  %8s    %10s    %8s      %10s  %10s  %10s  %10s  %10s  %s\n' % (zero, ekey, eid, elemName, f1, f2, f3, m1, m2, m3))
                    zero = ' '

            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            msg = ['']
            page_num += 1
        return page_num - 1


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

    def get_stats(self):
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

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        f.write('GridPointStressesVolumeObject write_f06 is not implemented...\n')
        return page_num

        #raise NotImplementedError()
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        msg = header + ['                                  S T R E S S E S   A T   G R I D   P O I N T S   - -     S U R F A C E       5\n',
                        '0                       SURFACE X-AXIS X  NORMAL(Z-AXIS)  Z         REFERENCE COORDINATE SYSTEM FOR SURFACE DEFINITION CID        0\n',
                        '     GRID      ELEMENT            STRESSES IN SURFACE SYSTEM           PRINCIPAL STRESSES            MAX             \n',
                        '     ID          ID    FIBRE   NORMAL-X   NORMAL-Y   SHEAR-XY     ANGLE      MAJOR      MINOR      SHEAR     VON MISES\n']
              #'0     13683          3736    TRIAX6         4.996584E+00   0.0            1.203093E+02   0.0            0.0            0.0'
              #'      13683          3737    TRIAX6        -4.996584E+00   0.0           -1.203093E+02   0.0            0.0            0.0'
              #'      13683                  *TOTALS*       6.366463E-12   0.0           -1.364242E-12   0.0            0.0            0.0'
        for ekey, nxs in sorted(iteritems(self.nx)):
            ekey2 = ekey
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
                (vals2, is_all_zeros) = writeFloats10E(vals)
                [nx, ny, nz, txy, tyz, txz, pressure, ovm] = vals2
                msg.append('%s%8s  %s %s %s   %s %s %s %s  %-s\n' % (zero, ekey, nx, ny, nz, txy, tyz, txz, pressure, ovm.rstrip()))
                zero = ' '
                ekey2 = ' '

        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
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
                    (elemName) = self.elemName[ekey][iLoad]
                    eid = self.eids[ekey][iLoad]

                    vals = [f1, f2, f3, m1, m2, m3]
                    vals2 = write_floats_13e(vals)
                    [f1, f2, f3, m1, m2, m3] = vals2
                    if eid == 0:
                        eid = ''

                    msg.append('%s  %8s    %10s    %8s      %s  %s  %s  %s  %s  %-s\n' % (zero, ekey, eid, elemName, f1, f2, f3, m1, m2, m3))
                    zero = ' '

            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            msg = ['']
            page_num += 1
        return page_num - 1
