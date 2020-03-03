from typing import List
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.result_objects.op2_objects import ScalarObject, get_times_dtype
from pyNastran.f06.f06_formatting import (
    write_floats_10e, _eigenvalue_header)


class GridPointSurfaceStressesArray(ScalarObject):
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
        self.location = None
        self._times = None

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    @property
    def is_real(self):
        return True
    @property
    def is_complex(self):
        return False

    def build(self):
        """sizes the vectorized attributes of the GridPointStressesArray"""
        if self.is_built:
            return
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        self.itime = 0
        self.ielement = 0
        self.itotal = 0

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        self.nelements //= self.ntimes

        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size)
        self.node_element = np.zeros((self.ntotal, 2), dtype=idtype)
        #oxx, oyy, txy, angle, major, minor, ovm
        self.data = np.zeros((self.ntimes, self.nelements, 8), dtype=fdtype)
        self.location = np.empty(self.ntotal, dtype='U8')

        self._times = np.zeros(self.ntimes, dtype=dtype)
        self.is_built = True

    #def build_dataframe(self):
        #"""creates a pandas dataframe"""
        #import pandas as pd
        #headers = self.get_headers()
        #element_node = [self.element_node[:, 0], self.element_node[:, 1]]
        #if self.nonlinear_factor not in (None, np.nan):
            #column_names, column_values = self._build_dataframe_transient_header()
            #self.data_frame = pd.Panel(self.data, items=column_values, major_axis=element_node, minor_axis=headers).to_frame()
            #self.data_frame.columns.names = column_names
        #else:
            #self.data_frame = pd.Panel(self.data, major_axis=element_node, minor_axis=headers).to_frame()
            #self.data_frame.columns.names = ['Static']
        #self.data_frame.index.names = ['NodeID', 'ElementID', 'Item']

    def add_sort1(self, dt, nid, eid, fiber, nx, ny, txy, angle, majorP, minorP, tmax, ovm):
        """unvectorized method for adding SORT1 transient data"""
        #assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.node_element[self.itotal, :] = [nid, eid]
        self.location[self.itotal] = fiber
        self.data[self.itime, self.itotal, :] = [nx, ny, txy, angle, majorP, minorP, tmax, ovm]
        self.itotal += 1

    def get_stats(self, short=False) -> List[str]:
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
        if self.nonlinear_factor not in (None, np.nan):  # transient
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
        return msg

    def get_headers(self) -> List[str]:
        headers = ['nx', 'ny', 'txy', 'angle', 'majorP', 'minorP', 'tmax', 'ovm']
        return headers

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []

        cid = self.refid
        axis_int = self.oCoord
        axis_map = {0 : 'X', 1 : 'Y', 2 : 'Z'}
        axis = axis_map[axis_int]
        msg = [
            '                                  S T R E S S E S   A T   G R I D   P O I N T S   - -     S U R F A C E    %s\n' % self.ogs_id,
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
                        '0%8s  %8s   %4s    %-10s %-10s %-10s  %8.4f %10s %10s %10s  %s\n' % (
                            nid, eid, fiber, nxi, nyi, txyi, anglei, majorpi, minorpi,
                            tmaxi, ovmi))
                else:
                    f06_file.write(
                        ' %8s  %8s   %4s    %-10s %-10s %-10s  %8.4f %10s %10s %10s  %s\n' % (
                            '', '', fiber, nxi, nyi, txyi, anglei, majorpi, minorpi,
                            tmaxi, ovmi))
                nid_old = nid
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            i = 0
            if self.is_sort1:
                for itime in range(ntimes):
                    for inid, (nid, eid) in enumerate(self.node_element):
                        t1 = self.data[itime, inid, :]
                        t2 = table.data[itime, inid, :]
                        (nx1, ny1, txy1, majorp1, minorp1, tmax1, ovm1) = t1
                        (nx2, ny2, txy2, majorp2, minorp2, tmax2, ovm2) = t2
                        if not np.allclose(t1, t2):
                        #if not np.array_equal(t1, t2):
                            msg += '%s %s\n  (%s, %s, %s, %s, %s, %s, %s)\n  (%s, %s, %s, %s, %s, %s, %s)\n' % (
                                nid, eid,
                                nx1, ny1, txy1, majorp1, minorp1, tmax1, ovm1,
                                nx2, ny2, txy2, majorp2, minorp2, tmax2, ovm2)
                            i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
            else:
                raise NotImplementedError(self.is_sort2)
            if i > 0:
                print(msg)
                raise ValueError(msg)
        return True


class GridPointStressesVolumePrincipalArray(ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase, apply_data_code=True)
        self.ntotal = 0
        self.ntimes = 0
        self.nelements = 0
        self.itotal = 0
        self.ielement = 0
        self.data = None
        self.itime = None
        self._times = None

    def get_headers(self) -> List[str]:
        headers = [
            'lxa', 'lxb', 'lxc',
            'lya', 'lyb', 'lyc',
            'lza', 'lzb', 'lzc',
            'sa', 'sb', 'sc',
            'epr', 'ovm']
        return headers

    def __eq__(self, table):  # pragma: no cover
        self._eq_header(table)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            i = 0
            if self.is_sort1:
                for itime in range(ntimes):
                    for inid, nid in enumerate(self.node):
                        t1 = self.data[itime, inid, :]
                        t2 = table.data[itime, inid, :]
                        (lxa1, lxb1, lxc1, lya1, lyb1, lyc1, lza1, lzb1, lzc1, sa1, sb1, sc1, epr1, ovm1) = t1
                        (lxa2, lxb2, lxc2, lya2, lyb2, lyc2, lza2, lzb2, lzc2, sa2, sb2, sc2, epr2, ovm2) = t2
                        if not np.allclose(t1, t2):
                        #if not np.array_equal(t1, t2):
                            msg += '%s\n  (%s, %s, %s, %s, %s, %s, %s)\n  (%s, %s, %s, %s, %s, %s, %s)\n' % (
                                nid,
                                lxa1, lxb1, lxc1, lya1, lyb1, lyc1, lza1,
                                lxa2, lxb2, lxc2, lya2, lyb2, lyc2, lza2)
                            i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
            else:
                raise NotImplementedError(self.is_sort2)
            if i > 0:
                print(msg)
                raise ValueError(msg)
        return True

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    @property
    def is_real(self):
        return True
    @property
    def is_complex(self):
        return False

    def build(self):
        """sizes the vectorized attributes of the GridPointStressesArray"""
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        #print('self.IDs', self.data)
        self.itime = 0
        self.ielement = 0
        self.itotal = 0

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        self.nelements //= self.ntimes

        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size)
        self.node = np.zeros(self.ntotal, dtype=idtype)
        #lxa, lxb, lxc, lya, lyb, lyc, lza, lzb, lzc, sa, sb, sc, epr, ovm
        self.data = np.zeros((self.ntimes, self.ntotal, 14), dtype=fdtype)
        self.location = np.empty(self.ntotal, dtype='U8')

        self._times = np.zeros(self.ntimes, dtype=dtype)
        self.is_built = True

    def get_stats(self, short=False) -> List[str]:
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
        if self.nonlinear_factor not in (None, np.nan):  # transient
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
        msg.append('  node.shape = %s\n' % str(self.node.shape).replace('L', ''))
        msg.append('  location.shape = %s\n' % str(self.location.shape).replace('L', ''))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg += self.get_data_code()
        return msg

    def add_sort1(self, dt, nid, lxa, lxb, lxc, lya, lyb, lyc, lza, lzb, lzc, sa, sb, sc, epr, ovm):
        assert isinstance(nid, int) and nid > 0, 'dt=%s nid=%s' % (dt, nid)
        self._times[self.itime] = dt
        self.node[self.itotal] = nid
        self.data[self.itime, self.itotal, :] = [lxa, lxb, lxc, lya, lyb, lyc, lza, lzb, lzc, sa, sb, sc, epr, ovm]
        self.itotal += 1

    #def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  #page_num=1, is_mag_phase=False, is_sort1=True):
        #pass


class GridPointStressesVolumeDirectArray(ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase, apply_data_code=True)
        self.ntotal = 0
        self.ntimes = 0
        self.nelements = 0
        self.itotal = 0
        self.ielement = 0
        self.data = None
        self.itime = None
        self._times = None

    def get_headers(self) -> List[str]:
        headers = ['ox', 'oy', 'oz', 'txy', 'tyz', 'txz', 'pressure', 'ovm']
        return headers

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    @property
    def is_real(self):
        return True
    @property
    def is_complex(self):
        return False

    def build(self):
        """sizes the vectorized attributes of the GridPointStressesArray"""
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        #print('self.IDs', self.data)
        self.itime = 0
        self.ielement = 0
        self.itotal = 0

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        self.nelements //= self.ntimes

        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size)
        self.node = np.zeros(self.ntotal, dtype=idtype)
        #oxx, oyy, txy, angle, major, minor, ovm
        self.data = np.zeros((self.ntimes, self.ntotal, 8), dtype=fdtype)
        self.location = np.empty(self.ntotal, dtype='U8')
        self._times = np.zeros(self.ntimes, dtype=dtype)
        self.is_built = True

    def get_stats(self, short=False) -> List[str]:
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
        if self.nonlinear_factor not in (None, np.nan):  # transient
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
        msg.append('  node.shape = %s\n' % str(self.node.shape).replace('L', ''))
        msg.append('  location.shape = %s\n' % str(self.location.shape).replace('L', ''))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg += self.get_data_code()
        return msg

    def add_sort1(self, dt, nid, nx, ny, nz, txy, tyz, txz, pressure, ovm):
        assert isinstance(nid, int) and nid > 0, 'dt=%s nid=%s' % (dt, nid)
        self._times[self.itime] = dt
        self.node[self.itotal] = nid
        self.data[self.itime, self.itotal, :] = [nx, ny, nz, txy, tyz, txz, pressure, ovm]
        self.itotal += 1

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        """
        '    D I R E C T   S T R E S S E S   A T   G R I D   P O I N T S   - -       V O L U M E      101'
        '        OUTPUT COORDINATE SYSTEM =       0  BASIC   '
        '    GRID            NORMAL-X    NORMAL-Y    NORMAL-Z      SHEAR-XY    SHEAR-YZ    SHEAR-ZX        MEAN      VON MISES'
        '    ID                                                                                           PRESSURE'
        '        1           1.455E+03  -1.548E+02  -2.927E+02    -1.573E+01   3.326E+01  -3.438E+03     -3.357E+02   6.188E+03'
        '        2           1.093E+03  -1.996E+02  -1.682E+02     1.542E+02   5.962E+01  -4.104E+03     -2.417E+02   7.227E+03'
        """
        if header is None:
            header = []


        cid = self.refid
        #axis_int = self.oCoord
        #axis_map = {0 : 'X', 1 : 'Y', 2 : 'Z'}
        #axis = axis_map[axis_int]
        msg = [
            '                    D I R E C T   S T R E S S E S   A T   G R I D   P O I N T S   - -       V O L U M E      %3i\n'
            '                              OUTPUT COORDINATE SYSTEM = %7i  ELEMENT \n'
            '     GRID            NORMAL-X    NORMAL-Y    NORMAL-Z      SHEAR-XY    SHEAR-YZ    SHEAR-ZX        MEAN      VON MISES\n'
            '     ID                                                                                           PRESSURE\n' % (
            #'     8086           6.136E-02   2.131E-01   8.353E-02    -2.268E+00  -2.274E-13   1.525E-13     -1.193E-01   3.930E+00'
            self.ogs_id, cid)
        ]

        ntimes = self.data.shape[0]

        nids = self.node
        zero = ' '
        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg))

            nx = self.data[itime, :, 0]
            ny = self.data[itime, :, 1]
            nz = self.data[itime, :, 2]
            txy = self.data[itime, :, 3]
            tyz = self.data[itime, :, 4]
            txz = self.data[itime, :, 5]
            pressure = self.data[itime, :, 6]
            ovm = self.data[itime, :, 7]
            for (nid, nxi, nyi, nzi, txyi, tyzi, txzi, pressurei, ovmi) in zip(
                nids, nx, ny, nz, txy, tyz, txz, pressure, ovm):
                [nxi, nyi, nzi, txyi, tyzi, txzi, pressurei, ovmi] = write_floats_10e([
                    nxi, nyi, nzi, txyi, tyzi, txzi, pressurei, ovmi])

                f06_file.write('%s%8s          %-10s  %-10s  %-10s    %-10s  %-10s  %-10s     %-10s  %-s\n' % (
                    zero, nid, nxi, nyi, nzi, txyi, tyzi, txzi, pressurei, ovmi.rstrip()))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            i = 0
            if self.is_sort1:
                for itime in range(ntimes):
                    for inid, nid in enumerate(self.node):
                        t1 = self.data[itime, inid, :]
                        t2 = table.data[itime, inid, :]
                        (nx1, ny1, nz1, txy1, tyz1, txz1, pressure1, ovm1) = t1
                        (nx2, ny2, nz2, txy2, tyz2, txz2, pressure2, ovm2) = t2
                        if not np.allclose(t1, t2):
                        #if not np.array_equal(t1, t2):
                            msg += '%s\n  (%s, %s, %s, %s, %s, %s, %s, %s)\n  (%s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                                nid,
                                nx1, ny1, nz1, txy1, tyz1, txz1, pressure1, ovm1,
                                nx2, ny2, nz2, txy2, tyz2, txz2, pressure2, ovm2)
                            i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
            else:
                raise NotImplementedError(self.is_sort2)
            if i > 0:
                print(msg)
                raise ValueError(msg)
        return True

#msg = [
    #' P R I N C I P A L   G R I D   P O I N T   S T R E S S   D I S C O N T I N U I T I E S  - -       V O L U M E       %s\n'
    #'                              OUTPUT COORDINATE SYSTEM = %7i  ELEMENT \n'
    #'                              GRID       PRINCIPAL STRESS DISCONTINUITY     MEAN      VON MISES      ERROR\n'
    #'                              ID           A          B          C          PRESSURE                 EST.\n' % (
        #ivolume, cid)
    #'                              8086         5.448E-09  9.886E-08  2.026E-15  2.484E-09  1.086E-07   5.716E-08'
#]
# not sure what result this is for
#zero = '                              '
#f06_file.write('%s%8s  %-10s %-10s %-10s   %-10s %-10s %-10s %-10s  %-s\n' % (
    #zero, nid, nxi, nyi, nzi, txyi, tyzi, txzi, pressurei, ovmi.rstrip()))

GridPointStressesVolumeDiscontinutiesArray = None # tCode=34

class GridPointStressesSurfaceDiscontinutiesArray(ScalarObject): # tCode=35
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase, apply_data_code=True)
        self.ntotal = 0
        self.ntimes = 0
        self.nelements = 0
        self.itotal = 0
        self.ielement = 0
        self.data = None
        self.itime = None
        #self.node_element = None
        self._times = None

    def get_headers(self) -> List[str]:
        headers = ['oxx', 'oyy', 'ozz', 'txy', 'pressure']
        return headers

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    @property
    def is_real(self):
        return True
    @property
    def is_complex(self):
        return False

    def build(self):
        """sizes the vectorized attributes of the GridPointStressesArray"""
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        #print('self.IDs', self.data)
        self.itime = 0
        self.ielement = 0
        self.itotal = 0

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        self.nelements //= self.ntimes

        self.node = np.zeros(self.ntotal, dtype='int32')
        #oxx, oyy, ozz, txy, pressure
        self.data = np.zeros((self.ntimes, self.ntotal, 5), dtype='float32')
        self.location = np.empty(self.ntotal, dtype='U8')
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'

        self._times = np.zeros(self.ntimes, dtype=dtype)
        self.is_built = True

    def get_stats(self, short=False) -> List[str]:
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
        if self.nonlinear_factor not in (None, np.nan):  # transient
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
        msg.append('  node.shape = %s\n' % str(self.node.shape).replace('L', ''))
        msg.append('  location.shape = %s\n' % str(self.location.shape).replace('L', ''))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg += self.get_data_code()
        return msg

    def add_sort1(self, dt, nid, oxx, oyy, ozz, txy, pressure):
        assert isinstance(nid, int) and nid > 0, 'dt=%s nid=%s' % (dt, nid)
        self._times[self.itime] = dt
        self.node[self.itotal] = nid
        self.data[self.itime, self.itotal, :] = [oxx, oyy, ozz, txy, pressure]
        self.itotal += 1
