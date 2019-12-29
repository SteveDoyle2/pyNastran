from typing import List
import numpy as np
from numpy import zeros, searchsorted, allclose


from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import (
    StressObject, StrainObject, OES_Object)
from pyNastran.f06.f06_formatting import write_imag_floats_13e, _eigenvalue_header


class ComplexCBushArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        #self.code = [self.format_code, self.sort_code, self.s_code]
        self.nelements = 0  # result specific

    @property
    def is_real(self):
        return False

    @property
    def is_complex(self):
        return True

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def _get_msgs(self):
        raise NotImplementedError()

    def get_headers(self):
        raise NotImplementedError()

    def build(self):
        """sizes the vectorized attributes of the ComplexCBushArray"""
        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        self.nelements //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.nelements, dtype='int32')

        #[tx, ty, tz, rx, ry, rz]
        self.data = zeros((self.ntimes, self.nelements, 6), dtype='complex64')

    def build_dataframe(self):
        """creates a pandas dataframe"""
        # Freq                          0.0                 2.5
        # ElementID Item
        # 10210     tx    0.010066-0.000334j  0.010066-0.000334j
        #           ty    0.000000+0.000000j  0.000000+0.000000j
        #           tz    0.431447-0.009564j  0.431461-0.009564j
        #           rx    0.000000+0.000000j  0.000000+0.000000j
        #           ry    0.000000+0.000000j  0.000000+0.000000j
        #           rz    0.000000+0.000000j  0.000000+0.000000j
        # 10211     tx   -0.000002+0.000000j -0.000002+0.000000j
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        self.data_frame = self._build_pandas_transient_elements(
            column_values, column_names,
            headers, self.element, self.data)

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
                        (tx1, ty1, tz1, rx1, ry1, rz1) = t1
                        (tx2, ty2, tz2, rx2, ry2, rz2) = t2
                        d = t1 - t2
                        if not allclose([tx1.real, tx1.imag, ty1.real, ty1.imag],
                                        [tx2.real, tx2.imag, ty2.real, ty2.imag], atol=0.0001):
                        #if not np.array_equal(t1, t2):
                            msg += '%-4s  (%s, %sj, %s, %sj)\n      (%s, %sj, %s, %sj)\n  dt12=(%s, %sj, %s, %sj)\n' % (
                                eid,
                                tx1.real, tx1.imag, ty1.real, ty1.imag,
                                tx2.real, tx2.imag, ty2.real, ty2.imag,
                                d[0].real, d[0].imag, d[1].real, d[1].imag,)
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

    def add_sort1(self, dt, eid, tx, ty, tz, rx, ry, rz):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [tx, ty, tz, rx, ry, rz]
        self.ielement += 1

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
            msg.append('  type=%s ntimes=%i nelements=%i; table_name=%r\n' % (
                self.__class__.__name__, ntimes, nelements, self.table_name))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i; table_name=%r\n' % (
                       self.__class__.__name__, nelements, self.table_name))
            ntimes_word = '1'
        msg.append('  eType\n')
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  element.shape = %s\n' % str(self.element.shape).replace('L', ''))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_element_index(self, eids):
        itot = searchsorted(eids, self.element)
        return itot

    def eid_to_element_node_index(self, eids):
        ind = searchsorted(eids, self.element)
        return ind

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = self.get_f06_header(is_mag_phase)

        if self.is_sort1:
            page_num = self._write_sort1_as_sort1(
                header, page_stamp, page_num, f06_file, msg_temp, is_mag_phase)
        else:
            raise NotImplementedError()
        return page_num

    def _write_sort1_as_sort1(self, header, page_stamp, page_num, f06_file, msg_temp, is_mag_phase):
        r"""
        C:\\Users\\sdoyle\\Dropbox\\move_tpl\\ofprand1.op2

        '                         C O M P L E X   F O R C E S   I N   B U S H   E L E M E N T S   ( C B U S H ) '
        '                                                          (REAL/IMAGINARY)'
        ' '
        '                  FREQUENCY         FORCE-X       FORCE-Y       FORCE-Z      MOMENT-X      MOMENT-Y      MOMENT-Z  '
        '0               0.0               1.006599E-02  0.0           4.314467E-01  0.0           0.0           0.0'
        '                                 -3.338092E-04  0.0          -9.563536E-03  0.0           0.0           0.0'
        """
        ntimes = self.data.shape[0]

        eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            tx = self.data[itime, :, 0]
            ty = self.data[itime, :, 1]
            tz = self.data[itime, :, 2]
            rx = self.data[itime, :, 3]
            ry = self.data[itime, :, 4]
            rz = self.data[itime, :, 5]
            for eid, itx, ity, itz, irx, iry, irz in zip(eids, tx, ty, tz, rx, ry, rz):
                [txr, tyr, tzr, rxr, ryr, rzr,
                 txi, tyi, tzi, rxi, ryi, rzi] = write_imag_floats_13e([itx, ity, itz, irx, iry, irz], is_mag_phase)
                #'0               1.000000E-01      0.0           2.912573E+00  0.0           0.0           0.0           0.0'
                #'                                    0.0         179.9942        0.0           0.0           0.0           0.0'
                f06_file.write('                %8i    %-13s %-13s %-13s %-13s %-13s %s\n'
                               '                %8s    %-13s %-13s %-13s %-13s %-13s %s\n' % (
                                   eid, txr, tyr, tzr, rxr, ryr, rzr,
                                   '', txi, tyi, tzi, rxi, ryi, rzi))

            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def write_op2(self, op2, op2_ascii, itable, new_result, date,
                  is_mag_phase=False, endian='>'):
        """writes an OP2"""
        # see TestOP2.test_op2_other_01
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

        device_code = self.device_code
        op2_ascii.write('  ntimes = %s\n' % self.ntimes)

        eids_device = self.element * 10 + self.device_code

        #fmt = '%2i %6f'
        #print('ntotal=%s' % (ntotal))
        #assert ntotal == 193, ntotal

        if self.is_sort1:
            struct1 = Struct(endian + b'i12f')
        else:
            raise NotImplementedError('SORT2')

        op2_ascii.write('nelements=%i\n' % nelements)

        for itime in range(self.ntimes):
            #print('3, %s' % itable)
            self._write_table_3(op2, op2_ascii, new_result, itable, itime)

            # record 4
            #print('stress itable = %s' % itable)
            itable -= 1
            #print('4, %s' % itable)
            header = [4, itable, 4,
                      4, 1, 4,
                      4, 0, 4,
                      4, ntotal, 4,
                      4 * ntotal]
            op2.write(pack('%ii' % len(header), *header))
            op2_ascii.write('r4 [4, 0, 4]\n')
            op2_ascii.write('r4 [4, %s, 4]\n' % (itable))
            op2_ascii.write('r4 [4, %i, 4]\n' % (4 * ntotal))

            tx = self.data[itime, :, 0]
            ty = self.data[itime, :, 1]
            tz = self.data[itime, :, 2]
            rx = self.data[itime, :, 3]
            ry = self.data[itime, :, 4]
            rz = self.data[itime, :, 5]
            for eid, itx, ity, itz, irx, iry, irz in zip(eids, tx, ty, tz, rx, ry, rz):
                [txr, tyr, tzr, rxr, ryr, rzr,
                 txi, tyi, tzi, rxi, ryi, rzi] = write_imag_floats_13e([itx, ity, itz, irx, iry, irz], is_mag_phase)
                data = [
                    eid,
                    itx.real, ity.real, itz.real, irx.real, iry.real, irz.real,
                    itx.imag, ity.imag, itz.imag, irx.imag, iry.imag, irz.imag]
                op2_ascii.write('  eid=%s data=%s\n' % (eids_device, str(data)))
                op2.write(struct1.pack(*data))

            itable -= 1
            header = [4 * ntotal,]
            op2.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable


class ComplexCBushStressArray(ComplexCBushArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexCBushArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        headers = ['tx', 'ty', 'tz', 'rx', 'ry', 'rz']
        return headers

    def get_f06_header(self, is_mag_phase=True):
        #'                         C O M P L E X   S T R E S S E S   I N   B U S H   E L E M E N T S   ( C B U S H ) '
        #' '
        #'                  FREQUENCY         STRESS-TX     STRESS-TY     STRESS-TZ    STRESS-RX     STRESS-RY     STRESS-RZ '
        #'0               1.000000E-01      0.0           2.912573E+00  0.0           0.0           0.0           0.0'
        #'                                    0.0         179.9942        0.0           0.0           0.0           0.0'

        if self.element_type == 102:
            element_header = '                         C O M P L E X   S T R E S S E S   I N   B U S H   E L E M E N T S   ( C B U S H ) \n'
            #''
            #' '
            #
            #'0               1.000000E-01      0.0           2.912573E+00  0.0           0.0           0.0           0.0'
            #'                                    0.0         179.9942        0.0           0.0           0.0           0.0'
        else:
            raise NotImplementedError('element_name=%r element_type=%s' % (self.element_name, self.element_type))

        if is_mag_phase:
            mag_phase = '                                                         (MAGNITUDE/PHASE)\n \n'
        else:
            mag_phase = '                                                          (REAL/IMAGINARY)\n \n'  # not tested

        words = [
            element_header,
            mag_phase,
            '                  FREQUENCY         STRESS-TX     STRESS-TY     STRESS-TZ    STRESS-RX     STRESS-RY     STRESS-RZ \n',
            '                   ID.                               FORCE\n',]
        return words


class ComplexCBushStrainArray(ComplexCBushArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexCBushArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        headers = ['tx', 'ty', 'tz', 'rx', 'ry', 'rz'] # tx, ty, tz, rx, ry, rz
        return headers

    def get_f06_header(self, is_mag_phase=True):
        """C:\\Users\\sdoyle\\Dropbox\\move_tpl\\ofprand1.op2"""
        if self.element_type == 102:
            element_header = '                         C O M P L E X   F O R C E S   I N   B U S H   E L E M E N T S   ( C B U S H ) \n'
        else:
            raise NotImplementedError('element_name=%r element_type=%s' % (self.element_name, self.element_type))

        if is_mag_phase:
            mag_phase = '                                                          (MAG/PHASE)\n'  # not tested
        else:
            mag_phase = '                                                          (REAL/IMAGINARY)\n'

        words = [
            element_header,
            mag_phase,
            ' \n',
            ' ',
            '                  FREQUENCY         FORCE-X       FORCE-Y       FORCE-Z      MOMENT-X      MOMENT-Y      MOMENT-Z  \n',
        ]
        return words
