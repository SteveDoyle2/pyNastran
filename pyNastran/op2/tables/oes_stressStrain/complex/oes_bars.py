from typing import List

import numpy as np
from numpy import zeros

from pyNastran.op2.result_objects.op2_objects import get_complex_times_dtype
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import (
    StressObject, StrainObject, OES_Object)
from pyNastran.f06.f06_formatting import write_imag_floats_13e


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

    @property
    def is_real(self) -> bool:
        return False

    @property
    def is_complex(self) -> bool:
        return True

    @property
    def nnodes_per_element(self) -> int:
        return 1

    def build(self):
        """sizes the vectorized attributes of the ComplexCBarArray"""
        #print('ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (
            #self.ntimes, self.nelements, self.ntotal, self.subtitle))
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
        dtype, idtype, cfdtype = get_complex_times_dtype(self.nonlinear_factor, self.size)
        self._times = zeros(self.ntimes, dtype=dtype)
        #self.element = array(self.nelements, dtype='|S8')

        #self.ntotal = self.nelements * nnodes
        self.element = zeros(self.ntotal, dtype=idtype)

        # the number is messed up because of the offset for the element's properties
        #if self.nelements * nnodes != self.ntotal:
            #msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (self.ntimes,
                                                                           #self.nelements, nnodes,
                                                                           #self.nelements * nnodes,
                                                                           #self.ntotal)
            #raise RuntimeError(msg)

        #[s1a, s2a, s3a, s4a, axial, s2a, s2b, s2c, s2d]
        self.data = zeros((self.ntimes, self.ntotal, 9), cfdtype)

    def build_dataframe(self):
        """creates a pandas dataframe"""
        # Freq            0.00001  10.00000 20.00000 30.00000 40.00000 50.00000 60.00000
        # ElementID Item
        # 13        s1a         0j       0j       0j       0j       0j       0j       0j
        #           s1b         0j       0j       0j       0j       0j       0j       0j
        #           s1c         0j       0j       0j       0j       0j       0j       0j
        #           s1d         0j       0j       0j       0j       0j       0j       0j
        #           axial  (-0+0j)  (-0+0j)  (-0+0j)  (-0+0j)  (-0+0j)  (-0+0j)  (-0+0j)
        #           s2a    (-0+0j)  (-0+0j)  (-0+0j)  (-0+0j)  (-0+0j)  (-0+0j)  (-0+0j)
        #           s2b    (-0+0j)  (-0+0j)  (-0+0j)  (-0+0j)  (-0+0j)  (-0+0j)  (-0+0j)
        #           s2c    (-0+0j)  (-0+0j)  (-0+0j)  (-0+0j)  (-0+0j)  (-0+0j)  (-0+0j)
        #           s2d    (-0+0j)  (-0+0j)  (-0+0j)  (-0+0j)  (-0+0j)  (-0+0j)  (-0+0j)
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        data_frame = self._build_pandas_transient_elements(column_values, column_names,
                                                           headers, self.element, self.data)
        #data_frame = pd.Panel(self.data, items=column_values,
                              #major_axis=self.element, minor_axis=headers).to_frame()
        #data_frame.columns.names = column_names
        #data_frame.index.names = ['ElementID', 'Item']
        self.data_frame = data_frame

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
                            msg += '%-4s  (%s, %s, %s, %s, %s, %s, %s, %s, %s)\n      (%s, %s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                                eid,
                                s1a1.real, s2a1.real, s3a1.real, s4a1.real, axial1.real, s2a1.real, s2b1.real, s2c1.real, s2d1.real,
                                s1a2.real, s2a2.real, s3a2.real, s4a2.real, axial2.real, s2a2.real, s2b2.real, s2c2.real, s2d2.real,
                                )
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

    def get_stats(self, short=False) -> List[str]:
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

        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i; table_name=%r\n' % (
                self.__class__.__name__, ntimes, nelements, self.table_name))
        else:
            msg.append('  type=%s nelements=%i; table_name=%r\n' % (
                self.__class__.__name__, nelements, self.table_name))
        msg.append(
            '  eType, cid\n'
            '  data: [ntimes, nnodes, 6] where 6=[%s]\n'
            '  element.shape = %s\n'
            '  data.shape = %s\n' % (
                ', '.join(self.get_headers()),
                str(self.element.shape).replace('L', ''),
                str(self.data.shape).replace('L', ''),
        ))

        msg.append('  CBAR\n')
        msg += self.get_data_code()
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
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
        if self.is_sort1:
            line1 = '            ELEMENT                    LOCATION       LOCATION       LOCATION       LOCATION             AVERAGE\n'
            line2 = '              ID.                          1              2              3              4             AXIAL STRESS\n'
        else:
            line1 = '                                       LOCATION       LOCATION       LOCATION       LOCATION             AVERAGE\n'
            if name == 'freq':
                line2 = '           FREQUENCY                       1              2              3              4             AXIAL STRESS\n'
            else:
                msg = 'name=%r\n\n%s' % (name, self.code_information())
                raise RuntimeError(msg)

        if self.is_stress:
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
        if self.is_sort1:
            self._write_sort1_as_sort1(f06_file, name, header, page_stamp, msg_temp, page_num,
                                       is_mag_phase=is_mag_phase)
        else:
            raise NotImplementedError()
        return page_num - 1

    def _write_sort1_as_sort1(self, f06_file, name, header, page_stamp, msg_temp, page_num,
                              is_mag_phase=False):
        ntimes = self.data.shape[0]
        for itime in range(ntimes):
            dt = self._times[itime]

            dt_line = ' %14s = %12.5E\n' % (name, dt)
            header[1] = dt_line
            msg = header + msg_temp
            f06_file.write('\n'.join(msg))

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

                f06_file.write(
                    '0%8i   %-13s  %-13s  %-13s  %-13s  %s\n'
                    ' %8s   %-13s  %-13s  %-13s  %-13s  %s\n' % (
                        eid, s1ar, s2ar, s3ar, s4ar, axialr,
                        '', s1ai, s2ai, s3ai, s4ai, axiali))

                f06_file.write(
                    ' %8s   %-13s  %-13s  %-13s  %s\n'
                    ' %8s   %-13s  %-13s  %-13s  %s\n' % (
                        '', s1br, s2br, s3br, s4br,
                        '', s1bi, s2bi, s3bi, s4bi))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num

    def write_op2(self, op2, op2_ascii, itable, new_result,
                  date, is_mag_phase=False, endian='>'):
        """writes an OP2"""
        import inspect
        from struct import Struct, pack
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write('%s.write_op2: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        if itable == -1:
            self._write_table_header(op2, op2_ascii, date)
            itable = -3

        #if isinstance(self.nonlinear_factor, float):
            #op2_format = '%sif' % (7 * self.ntimes)
            #raise NotImplementedError()
        #else:
            #op2_format = 'i21f'
        #s = Struct(op2_format)

        eids = self.element
        eids_device = eids * 10 + self.device_code

        # table 4 info
        #ntimes = self.data.shape[0]
        #nnodes = self.data.shape[1]
        nelements = self.data.shape[1]

        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)

        ntotali = self.num_wide
        ntotal = ntotali * nelements

        #print('shape = %s' % str(self.data.shape))
        #assert self.ntimes == 1, self.ntimes

        op2_ascii.write('  ntimes = %s\n' % self.ntimes)

        #fmt = '%2i %6f'
        #print('ntotal=%s' % (ntotal))
        #assert ntotal == 193, ntotal

        if self.is_sort1:
            struct1 = Struct(endian + b'i 18f')
        else:
            raise NotImplementedError('SORT2')

        op2_ascii.write('%s-nelements=%i\n' % (self.element_name, nelements))
        for itime in range(self.ntimes):
            self._write_table_3(op2, op2_ascii, new_result, itable, itime)

            # record 4
            #print('stress itable = %s' % itable)
            itable -= 1
            header = [4, itable, 4,
                      4, 1, 4,
                      4, 0, 4,
                      4, ntotal, 4,
                      4 * ntotal]
            op2.write(pack('%ii' % len(header), *header))
            op2_ascii.write('r4 [4, 0, 4]\n')
            op2_ascii.write('r4 [4, %s, 4]\n' % (itable))
            op2_ascii.write('r4 [4, %i, 4]\n' % (4 * ntotal))

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
            for eid_device, s1ai, s2ai, s3ai, s4ai, axiali, s2ai, s2bi, s2ci, s2di in zip(eids_device, sa1, sa2, sa3, sa4, axial, sb1, sb2, sb3, sb4):

                data = [eid_device,
                        s1ai.real, s2ai.real, s3ai.real, s4ai.real, axiali.real, s2ai.real, s2bi.real, s2ci.real, s2di.real,
                        s1ai.imag, s2ai.imag, s3ai.imag, s4ai.imag, axiali.imag, s2ai.imag, s2bi.imag, s2ci.imag, s2di.imag]
                op2_ascii.write('  eid_device=%s data=%s\n' % (eid_device, str(data)))
                op2.write(struct1.pack(*data))

            itable -= 1
            header = [4 * ntotal,]
            op2.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable

class ComplexBarStressArray(ComplexBarArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexBarArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        headers = ['s1a', 's1b', 's1c', 's1d', 'axial',
                   's2a', 's2b', 's2c', 's2d', ]
        #headers = ['s1a', 's1b', 's1c', 's1d', 's1e', 'axial',
                   #'s2a', 's2b', 's2c', 's2d', 's2e', ]
        return headers

class ComplexBarStrainArray(ComplexBarArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexBarArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        headers = ['e1a', 'e1b', 'e1c', 'e1d', 'axial',
                   'e2a', 'e2b', 'e2c', 'e2d', ]
        #headers = ['e1a', 'e1b', 'e1c', 'e1d', 'e1e', 'axial',
                   #'e2a', 'e2b', 'e2c', 'e2d', 'e2e', ]
        return headers
