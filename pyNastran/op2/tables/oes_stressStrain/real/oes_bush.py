from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
import numpy as np
from numpy import zeros, array_equal
from itertools import count

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import write_floats_13e, _eigenvalue_header
try:
    import pandas as pd
except ImportError:
    pass


class RealBushArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        self.eType = {}
        #self.code = [self.format_code, self.sort_code, self.s_code]
        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.ielement = 0
        self.nelements = 0  # result specific
        #print('RealBushArray.nonlinear_factor =', self.nonlinear_factor)

    def is_real(self):
        return True

    def is_complex(self):
        return False

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def _get_msgs(self):
        raise NotImplementedError('%s needs to implement _get_msgs' % self.__class__.__name__)

    def get_headers(self):
        raise NotImplementedError('%s needs to implement get_headers' % self.__class__.__name__)
        return headers

    def build(self):
        if self.is_built:
            return
        #print("self.ielement =", self.ielement)
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal

        if self.element_type == 102:
            nnodes_per_element = 1
        else:
            raise NotImplementedError(self.element_type)

        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("***name=%s type=%s nnodes_per_element=%s ntimes=%s nelements=%s ntotal=%s" % (
            #self.element_name, self.element_type, nnodes_per_element, self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, int):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.ntotal, dtype='int32')

        # [tx, ty, tz, rx, ry, rz]
        self.data = zeros((self.ntimes, self.ntotal, 6), dtype='float32')

    def build_dataframe(self):
        headers = self.get_headers()
        if self.nonlinear_factor is not None:
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
            self.data_frame.index.names = ['ElementID', 'Item']
        else:
            self.data_frame = pd.Panel(self.data, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = ['Static']
            self.data_frame.index.names = ['ElementID', 'Item']

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()
        self._eq_header(table)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            i = 0
            if self.is_sort1():
                for itime in range(ntimes):
                    for ieid, eid, in enumerate(self.element):
                        t1 = self.data[itime, ieid, :]
                        t2 = table.data[itime, ieid, :]
                        (fx1, fy1, fz1, mx1, my1, mz1) = t1
                        (fx2, fy2, fz2, mx2, my2, mz2) = t2
                        if not allclose(t1, t2):
                        #if not np.array_equal(t1, t2):
                            msg += '%s\n  (%s, %s, %s)\n  (%s, %s, %s)\n' % (
                                eid,
                                fx1, fy1, fz1, mx1, my1, mz1,
                                fx2, fy2, fz2, mx2, my2, mz2)
                            i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
            else:
                raise NotImplementedError(self.is_sort2())
            if i > 0:
                print(msg)
                raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, tx, ty, tz, rx, ry, rz):
        assert isinstance(eid, int)
        self._times[self.itime] = dt
        self.element[self.itotal] = eid
        self.data[self.itime, self.itotal, :] = [tx, ty, tz, rx, ry, rz]
        self.itotal += 1
        self.ielement += 1

    def get_stats(self):
        if not self.is_built:
            return ['<%s>\n' % self.__class__.__name__,
                    '  ntimes: %i\n' % self.ntimes,
                    '  ntotal: %i\n' % self.ntotal,
                    ]

        nelements = self.ntotal
        ntimes = self.ntimes
        ntotal = self.ntotal
        nelements = self.ntotal

        msg = []
        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i\n'
                       % (self.__class__.__name__, ntimes, nelements))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i\n'
                       % (self.__class__.__name__, nelements))
            ntimes_word = 1
        headers = self.get_headers()

        n = len(headers)
        assert n == self.data.shape[2], 'nheaders=%s shape=%s' % (n, str(self.data.shape))
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element)  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        ind = ravel([searchsorted(self.element == eid) for eid in eids])
        return ind

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg = self._get_msgs()
        (ntimes, ntotal) = self.data.shape[:2]
        eids = self.element

        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg))
            #[tx, ty, tz, rx, ry, rz]
            tx = self.data[itime, :, 0]
            ty = self.data[itime, :, 1]
            tz = self.data[itime, :, 2]
            rx = self.data[itime, :, 3]
            ry = self.data[itime, :, 4]
            rz = self.data[itime, :, 5]

            for eid, txi, tyi, tzi, rxi, ryi, rzi in zip(
                eids, tx, ty, tz, rx, ry, rz):

                vals = [txi, tyi, tzi, rxi, ryi, rzi]
                vals2 = write_floats_13e(vals)
                [txi, tyi, tzi, rxi, ryi, rzi] = vals2
                f.write('0                   %8i     %-13s %-13s %-13s %-13s %-13s %s\n' % (
                        eid, txi, tyi, tzi, rxi, ryi, rzi))
            f.write(page_stamp % page_num)
            page_num += 1
        if self.nonlinear_factor is None:
            page_num -= 1
        return page_num


class RealBushStressArray(RealBushArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealBushArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        headers = ['tx', 'ty', 'tz', 'rx', 'ry', 'rz']
        return headers

    def _get_msgs(self):
        if self.element_type == 102:
            pass
        else:
            raise NotImplementedError(self.element_type)

        msg = [
            '                                  S T R E S S E S   I N   B U S H   E L E M E N T S        ( C B U S H )\n \n',
            '                  ELEMENT-ID        STRESS-TX     STRESS-TY     STRESS-TZ    STRESS-RX     STRESS-RY     STRESS-RZ \n',
        ]
        return msg

class RealBushStrainArray(RealBushArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealBushArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        headers = ['tx', 'ty', 'tz', 'rx', 'ry', 'rz']
        return headers

    def _get_msgs(self):
        if self.element_type == 102:
            pass
        else:
            raise NotImplementedError(self.element_type)

        msg = [
            '                                    S T R A I N S   I N   B U S H   E L E M E N T S        ( C B U S H )\n'
            ' \n'
            '                  ELEMENT-ID        STRAIN-TX     STRAIN-TY     STRAIN-TZ    STRAIN-RX     STRAIN-RY     STRAIN-RZ \n'
        ]
        return msg
