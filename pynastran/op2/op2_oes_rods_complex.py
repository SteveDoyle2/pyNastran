from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
from six.moves import range
import numpy as np
from numpy import zeros, array_equal, allclose

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import write_imag_floats_13e, get_key0, _eigenvalue_header
try:
    import pandas as pd
except ImportError:
    pass


class ComplexRodArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        self.eType = {}
        #self.code = [self.format_code, self.sort_code, self.s_code]

        self.nelements = 0  # result specific

    def is_real(self):
        return False

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
        if self.is_built:
            return

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
        if isinstance(self.nonlinear_factor, int):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.nelements, dtype='int32')

        #[axial, torsion]
        self.data = zeros((self.ntimes, self.nelements, 2), dtype='complex64')

    def build_dataframe(self):
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
        self.data_frame.columns.names = column_names
        self.data_frame.index.names = ['ElementID', 'Item']

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()
        self._eq_header(table)
        if not np.array_equal(self.element, table.element):
            assert self.element.shape == table.element.shape, 'shape=%s element.shape=%s' % (self.element.shape, table.element.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            for eid, eid2 in zip(self.element, table.element):
                msg += '%s, %s\n' % (eid, eid2)
            print(msg)
            raise ValueError(msg)
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
                        (axial1, torsion1) = t1
                        (axial2, torsion2) = t2
                        d = t1 - t2
                        if not allclose([axial1.real, axial1.imag, torsion1.real, torsion1.imag],
                                        [axial2.real, axial2.imag, torsion2.real, torsion2.imag], atol=0.0001):
                        #if not np.array_equal(t1, t2):
                            msg += '%-4s  (%s, %sj, %s, %sj)\n      (%s, %sj, %s, %sj)\n  dt12=(%s, %sj, %s, %sj)\n' % (
                                eid,
                                axial1.real, axial1.imag, torsion1.real, torsion1.imag,
                                axial2.real, axial2.imag, torsion2.real, torsion2.imag,
                                d[0].real, d[0].imag, d[1].real, d[1].imag,)
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

    def add_sort1(self, dt, eid, axial, torsion):
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [axial, torsion]
        self.ielement += 1

    def get_stats(self):
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
            ntimes_word = 1
        msg.append('  eType\n')
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_element_index(self, eids):
        itot = searchsorted(eids, self.element)
        return itot

    def eid_to_element_node_index(self, eids):
        ind = searchsorted(eids, self.element)
        return ind

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = self.get_f06_header(is_mag_phase)

        if self.is_sort1():
            page_num = self._write_sort1_as_sort1(header, page_stamp, page_num, f, msg_temp, is_mag_phase)
        else:
            raise NotImplementedError()
        return page_num

    def _write_sort1_as_sort1(self, header, page_stamp, page_num, f, msg_temp, is_mag_phase):
        ntimes = self.data.shape[0]

        eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            axial = self.data[itime, :, 0]
            torsion = self.data[itime, :, 1]
            for eid, iaxial, itorsion in zip(eids, axial, torsion):
                [axialr, torsionr, axiali, torsioni] = write_imag_floats_13e([iaxial, itorsion], is_mag_phase)
                f.write('                %8i                 %-13s / %-13s                 %-13s / %s\n' % (
                    eid, axialr, axiali, torsionr, torsioni))

            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class ComplexRodStressArray(ComplexRodArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexRodArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        headers = ['axial', 'torsion']
        return headers

    def get_f06_header(self, is_mag_phase=True):
        if self.element_type == 1:
            element_header = '                           C O M P L E X   S T R E S S E S   I N   R O D   E L E M E N T S   ( C R O D )\n'
        elif self.element_type == 3:
            element_header = '                          C O M P L E X   S T R E S S E S   I N   R O D   E L E M E N T S   ( C T U B E )\n'
        elif self.element_type == 10:
            element_header = '                         C O M P L E X   S T R E S S E S   I N   R O D   E L E M E N T S   ( C O N R O D )\n'
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
            '                 ELEMENT                             AXIAL                                         TORQUE\n',
            '                   ID.                               FORCE\n',]
           #'                       1                 -2.459512E+05 /  3.377728E+04                  0.0          /  0.0\n',]
        return words


class ComplexRodStrainArray(ComplexRodArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexRodArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        headers = ['axial', 'torsion']
        return headers

    def get_f06_header(self, is_mag_phase=True):
        if self.element_type == 1:
            element_header = '                           C O M P L E X    S T R A I N S    I N   R O D   E L E M E N T S   ( C R O D )\n'
        elif self.element_type == 3:
            element_header = '                          C O M P L E X    S T R A I N S    I N   R O D   E L E M E N T S   ( C T U B E )\n'
        elif self.element_type == 10:
            element_header = '                         C O M P L E X    S T R A I N S    I N   R O D   E L E M E N T S   ( C O N R O D )\n'
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
            '                 ELEMENT                             AXIAL                                         TORQUE\n',
            '                   ID.                               FORCE\n',
           #'                       1                 -2.459512E+05 /  3.377728E+04                  0.0          /  0.0\n',
        ]
        return words
