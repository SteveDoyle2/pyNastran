from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import write_floats_13e, write_imag_floats_13e
import numpy as np
from numpy import concatenate, zeros
try:
    import pandas as pd
except ImportError:
    pass


class ComplexBarArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        self.result_flag = 0
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.itime = 0
        self.nelements = 0  # result specific
        #self.cid = {}  # gridGauss

        if is_sort1:
            #sort1
            pass
        else:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def is_real(self):
        return False

    def is_complex(self):
        return True

    def build(self):
        #print('ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (self.ntimes, self.nelements, self.ntotal, self.subtitle))
        if self.is_built:
            return
        nnodes = 1

        #self.names = []
        #self.nelements //= nnodes
        self.nelements //= self.ntimes
        #self.ntotal //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        self.is_built = True
        #print('ntotal=%s ntimes=%s nelements=%s' % (self.ntotal, self.ntimes, self.nelements))
        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        self._times = zeros(self.ntimes, 'float32')
        #self.element = array(self.nelements, dtype='|S8')

        #self.ntotal = self.nelements * nnodes
        self.element = zeros(self.ntotal, 'int32')

        # the number is messed up because of the offset for the element's properties
        if not self.nelements * nnodes == self.ntotal:
            msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (self.ntimes,
                                                                           self.nelements, nnodes,
                                                                           self.nelements * nnodes,
                                                                           self.ntotal)
            raise RuntimeError(msg)

        #[s1a, s2a, s3a, s4a, axial, s2a, s2b, s2c, s2d]
        self.data = zeros((self.ntimes, self.ntotal, 9), 'complex64')

    def build_dataframe(self):
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
        self.data_frame.columns.names = column_names
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
                    for ieid, eid in enumerate(self.element):
                        t1 = self.data[itime, ieid, :]
                        t2 = table.data[itime, ieid, :]
                        (s1a1, s2a1, s3a1, s4a1, axial1, s2a1, s2b1, s2c1, s2d1) = t1
                        (s1a2, s2a2, s3a2, s4a2, axial2, s2a2, s2b2, s2c2, s2d2) = t2
                        #d = t1 - t2
                        if not np.allclose(
                            [s1a1.real, s2a1.real, s3a1.real, s4a1.real, axial1.real, s2a1.real, s2b1.real, s2c1.real, s2d1.real],
                            [s1a2.real, s2a2.real, s3a2.real, s4a2.real, axial2.real, s2a2.real, s2b2.real, s2c2.real, s2d2.real], atol=0.0001):
                        #if not np.array_equal(t1, t2):
                            msg += '%-4s  (%s, %s, %s, %s, %s, %s, %s, %s, %s)\n      (%s, %s, %s, %s, %s, %s, %s, %s, %s)\m' % (
                                eid,
                                s1a1.real, s2a1.real, s3a1.real, s4a1.real, axial1.real, s2a1.real, s2b1.real, s2c1.real, s2d1.real,
                                s1a2.real, s2a2.real, s3a2.real, s4a2.real, axial2.real, s2a2.real, s2b2.real, s2c2.real, s2d2.real,
                                )
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

    def add_new_eid_sort1(self, dt, eid, e1a, e2a, e3a, e4a, axial,
                          e1b, e2b, e3b, e4b,):
        #self.e1[dt][eid] = [e1a, e1b]
        #self.e2[dt][eid] = [e2a, e2b]
        #self.e3[dt][eid] = [e3a, e3b]
        #self.e4[dt][eid] = [e4a, e4b]
        #self.axial[dt][eid] = axial

        self._times[self.itime] = dt
        #[sa1, sa2, sa3, sa4, axial, sb1, sb2, sb3, sb4]
        self.data[self.itime, self.itotal, :] = [e1a, e2a, e3a, e4a, axial,
                                                 e1b, e2b, e3b, e4b,]
        self.element[self.itotal] = eid
        self.itotal += 1

    def get_stats(self):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        #ntotal = self.ntotal
        nelements2 = self.element.shape[0]
        assert nelements, nelements2
        msg = []

        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%i\n' % (self.__class__.__name__, nelements))
        msg.append('  eType, cid\n')
        msg.append('  data: [ntimes, nnodes, 6] where 6=[%s]\n' % str(', '.join(self.get_headers())))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))

        msg.append('  CBAR\n  ')
        msg += self.get_data_code()
        return msg

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        #msg_temp, nnodes = get_f06_header(self, is_mag_phase, is_sort1)

        # write the f06

        if is_mag_phase:
            mag_phase = '                                                          (MAGNITUDE/PHASE)'
        else:
            mag_phase = '                                                          (REAL/IMAGINARY)'

        # force
        #msg = header + [
            #'                             C O M P L E X   F O R C E S   I N   B A R   E L E M E N T S   ( C B A R )',
            #mag_phase,
            #' ',
            #'                     BEND-MOMENT-END-A            BEND-MOMENT-END-B                  SHEAR',
            #'   FREQUENCY       PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE',
        #]
        name = self.data_code['name']
        if self.is_sort1():
            line1 = '            ELEMENT                    LOCATION       LOCATION       LOCATION       LOCATION             AVERAGE\n'
            line2 = '              ID.                          1              2              3              4             AXIAL STRESS\n'
        else:
            line1 = '                                       LOCATION       LOCATION       LOCATION       LOCATION             AVERAGE\n'
            if name == 'freq':
                line2 = '           FREQUENCY                       1              2              3              4             AXIAL STRESS\n'
            else:
                msg = 'name=%r\n\n%s' % (name, self.code_information())
                raise RuntimeError(msg)

        if self.is_stress():
            stress_strain = '                             C O M P L E X   S T R E S S E S   I N   B A R   E L E M E N T S   ( C B A R )'
        else:
            stress_strain = '                             C O M P L E X   S T R A I N S   I N   B A R   E L E M E N T S   ( C B A R )'

        msg_temp = [
            stress_strain,
            mag_phase,
            ' ',
            line1,
            line2,
        ]
        if self.is_sort1():
            self._write_sort1_as_sort1(f, name, header, page_stamp, msg_temp, page_num,
                                       is_mag_phase=is_mag_phase)
        else:
            raise NotImplementedError()
        return page_num - 1

    def _write_sort1_as_sort1(self, f, name, header, page_stamp, msg_temp, page_num,
                              is_mag_phase=False):
        ntimes = self.data.shape[0]
        for itime in range(ntimes):
            dt = self._times[itime]

            dt_line = ' %14s = %12.5E\n' % (name, dt)
            header[1] = dt_line
            msg = header + msg_temp
            f.write('\n'.join(msg))

            sa1 = self.data[itime, :, 0]
            sa2 = self.data[itime, :, 1]
            sa3 = self.data[itime, :, 2]
            sa4 = self.data[itime, :, 3]
            axial = self.data[itime, :, 4]
            sb1 = self.data[itime, :, 5]
            sb2 = self.data[itime, :, 6]
            sb3 = self.data[itime, :, 7]
            sb4 = self.data[itime, :, 8]
            #[sa1, sa2, sa3, sa4, axial, sb1, sb2, sb3, sb4]

            eids = self.element
            for eid, s1ai, s2ai, s3ai, s4ai, axiali, s2ai, s2bi, s2ci, s2di in zip(eids, sa1, sa2, sa3, sa4, axial, sb1, sb2, sb3, sb4):
                vals = (s1ai, s2ai, s3ai, s4ai, axiali,
                        s2ai, s2bi, s2ci, s2di)
                vals2 = write_imag_floats_13e(vals, is_mag_phase)
                (s1ar, s2ar, s3ar, s4ar, axialr,
                 s1br, s2br, s3br, s4br,
                 s1ai, s2ai, s3ai, s4ai, axiali,
                 s1bi, s2bi, s3bi, s4bi) = vals2

                msg.append('0%8i   %-13s  %-13s  %-13s  %-13s  %s\n' % (eid, s1ar, s2ar, s3ar, s4ar, axialr))
                msg.append(' %8s   %-13s  %-13s  %-13s  %-13s  %s\n' % ('', s1ai, s2ai, s3ai, s4ai, axiali))

                msg.append(' %8s   %-13s  %-13s  %-13s  %s\n' % ('', s1br, s2br, s3br, s4br))
                msg.append(' %8s   %-13s  %-13s  %-13s  %s\n' % ('', s1bi, s2bi, s3bi, s4bi))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num

class ComplexBarStressArray(ComplexBarArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexBarArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        headers = ['s1a', 's1b', 's1c', 's1d', 'axial',
                   's2a', 's2b', 's2c', 's2d', ]
        #headers = ['s1a', 's1b', 's1c', 's1d', 's1e', 'axial',
                   #'s2a', 's2b', 's2c', 's2d', 's2e', ]
        return headers

class ComplexBarStrainArray(ComplexBarArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexBarArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        headers = ['e1a', 'e1b', 'e1c', 'e1d', 'axial',
                   'e2a', 'e2b', 'e2c', 'e2d', ]
        #headers = ['e1a', 'e1b', 'e1c', 'e1d', 'e1e', 'axial',
                   #'e2a', 'e2b', 'e2c', 'e2d', 'e2e', ]
        return headers
