import numpy as np
from numpy import zeros, allclose

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.result_objects.op2_objects import get_complex_times_dtype
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import (
    StressObject, StrainObject, OES_Object,
    oes_complex_data_code, get_scode, set_element_case,
    set_freq_case, set_complex_modes_case)
from pyNastran.f06.f06_formatting import write_imag_floats_13e, _eigenvalue_header # get_key0,


ELEMENT_NAME_TO_ELEMENT_TYPE = {
    'CROD' : 1,
    'CONROD' : 10,
    'CTUBE' : 3,
}

class ComplexRodArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        #self.code = [self.format_code, self.sort_code, self.s_code]
        self.nelements = 0  # result specific

    @property
    def is_real(self) -> bool:
        return False

    @property
    def is_complex(self) -> bool:
        return True

    @property
    def nnodes_per_element(self) -> int:
        return 1

    def _reset_indices(self) -> None:
        self.itotal = 0
        self.ielement = 0

    def _get_msgs(self):
        raise NotImplementedError()

    def get_headers(self):
        raise NotImplementedError()

    def build(self):
        """sizes the vectorized attributes of the ComplexRodArray"""
        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        self.nelements //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        idtype, cfdtype = get_complex_times_dtype(self.size)

        self._times = zeros(self.ntimes, dtype=self.analysis_fmt)
        self.element = zeros(self.nelements, dtype=idtype)

        #[axial, torsion]
        self.data = zeros((self.ntimes, self.nelements, 2), dtype=cfdtype)

    def build_dataframe(self):
        """creates a pandas dataframe"""
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()

        #Freq              0.00001  10.00000 20.00000 30.00000 40.00000 50.00000 60.00000
        #ElementID Item
        #14        axial         0j       0j       0j       0j       0j       0j       0j
        #          torsion       0j       0j       0j       0j       0j       0j       0j
        #15        axial         0j       0j       0j       0j       0j       0j       0j
        #          torsion       0j       0j       0j       0j       0j       0j       0j
        self.data_frame = self._build_pandas_transient_elements(
            column_values, column_names,
            headers, self.element, self.data)

    @classmethod
    def _add_case(cls,
                  table_name, element_name, isubcase,
                  is_sort1, is_random, is_msc,
                  random_code, title, subtitle, label):
        is_strain = 'Strain' in cls.__name__
        assert isinstance(is_strain, bool), is_strain
        assert isinstance(element_name, str), element_name
        num_wide = 5
        data_code = oes_complex_data_code(
            table_name,
            element_name, num_wide,
            is_sort1=is_sort1, is_random=is_random,
            random_code=random_code, title=title, subtitle=subtitle, label=label,
            is_msc=is_msc)
        #data_code['modes'] = modes
        #data_code['eigns'] = eigenvalues
        #data_code['mode_cycles'] = mode_cycles
        #
        #stress_bits[1] = 1  # strain bit (vs. stress)
        #stress_bits[2] = 1  # curvature bit (vs. fiber)     ---> always 0
        #stress_bits[3] = 1  # strain bit (vs. stress)
        #stress_bits[4] = 1  # von mises bit (vs. max shear) ---> always 0
        if is_strain:
            # fiber   # 2  =0
            # strain  # 1,3=1
            #stress_bits[2] == 0
            stress_bits = [1, 1, 0, 1, 0]
            #data_code['s_code'] = 1 # strain?
        else:
            # fiber   # 2   =0
            # stress  # 1, 3=0
            stress_bits = [0, 0, 0, 0, 0]
            #data_code['s_code'] = 0

        s_code = get_scode(stress_bits)
        data_code['stress_bits'] = stress_bits
        data_code['s_code'] = s_code

        assert stress_bits[1] == stress_bits[3]  # strain

        element_type = ELEMENT_NAME_TO_ELEMENT_TYPE[element_name.upper()]
        #data_code['element_name'] = element_name.upper()
        data_code['element_type'] = element_type
        return data_code

    @classmethod
    def add_freq_case(cls, table_name, element, data, isubcase,
                      freqs,
                      element_name: str,
                      is_sort1=True, is_random=False, is_msc=True,
                      random_code=0, title='', subtitle='', label=''):
        data_code = cls._add_case(
            table_name, element_name, isubcase,
            is_sort1, is_random, is_msc,
            random_code, title, subtitle, label)

        obj = set_freq_case(cls, is_sort1, isubcase, data_code,
                            set_element_case, (element, data), freqs)
        #obj.stress_bits = stress_bits
        return obj

    @classmethod
    def add_complex_modes_case(cls, table_name, element, data, isubcase,
                               modes, eigrs, eigis,
                               element_name: str,
                               is_sort1=True, is_random=False, is_msc=True,
                               random_code=0, title='', subtitle='', label=''):
        data_code = cls._add_case(
            table_name, element_name, isubcase,
            is_sort1, is_random, is_msc,
            random_code, title, subtitle, label)

        obj = set_complex_modes_case(cls, is_sort1, isubcase, data_code,
                                     set_element_case, (element, data), modes, eigrs, eigis)
        #obj.stress_bits = stress_bits
        return obj

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
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
            if self.is_sort1:
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
                raise NotImplementedError(self.is_sort2)
            if i > 0:
                print(msg)
                raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, axial, torsion):
        """unvectorized method for adding SORT1 transient data"""
        assert self.sort_method == 1, self
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [axial, torsion]
        self.ielement += 1

    def get_stats(self, short: bool=False) -> list[str]:
        if not self.is_built:
            return [
                f'<{self.__class__.__name__}>; table_name={self.table_name!r}\n',
                f'  ntimes: {self.ntimes:d}\n',
                f'  ntotal: {self.ntotal:d}\n',
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
            msg.append('  type=%s nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, nelements, self.table_name))
            ntimes_word = '1'
        msg.append('  eType\n')
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append(f'  element.shape = {self.element.shape}\n')
        msg.append(f'  data.shape = {self.data.shape}\n')
        msg.append(f'  element type: {self.element_name}-{self.element_type}\n')
        msg += self.get_data_code()
        return msg

    def get_element_index(self, eids):
        itot = np.searchsorted(eids, self.element)
        return itot

    def eid_to_element_node_index(self, eids):
        ind = np.searchsorted(eids, self.element)
        return ind

    def get_f06_header(self, is_mag_phase=True):
        raise NotImplementedError('overwrite this')

    def write_f06(self, f06_file, header=None, page_stamp: str='PAGE %s', page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        if header is None:
            header = []
        msg_temp = self.get_f06_header(is_mag_phase)

        if self.is_sort1:
            page_num = self._write_sort1_as_sort1(header, page_stamp, page_num, f06_file, msg_temp, is_mag_phase)
        else:
            raise NotImplementedError()
        return page_num

    def _write_sort1_as_sort1(self, header, page_stamp, page_num, f06_file, msg_temp, is_mag_phase):
        ntimes = self.data.shape[0]

        eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            axial = self.data[itime, :, 0]
            torsion = self.data[itime, :, 1]
            for eid, iaxial, itorsion in zip(eids, axial, torsion):
                [axialr, torsionr, axiali, torsioni] = write_imag_floats_13e([iaxial, itorsion], is_mag_phase)
                f06_file.write('                %8i                 %-13s / %-13s                 %-13s / %s\n' % (
                    eid, axialr, axiali, torsionr, torsioni))

            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def write_op2(self, op2_file, op2_ascii, itable, new_result, date,
                  is_mag_phase=False, endian='>'):
        """writes an OP2"""
        import inspect
        from struct import Struct, pack
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write(f'{self.__class__.__name__}.write_op2: {call_frame[1][3]}\n')

        if itable == -1:
            self._write_table_header(op2_file, op2_ascii, date)
            itable = -3

        eids = self.element

        # table 4 info
        #ntimes = self.data.shape[0]
        #nnodes = self.data.shape[1]
        nelements = self.data.shape[1]

        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)

        ntotali = self.num_wide
        ntotal = ntotali * nelements

        device_code = self.device_code
        op2_ascii.write(f'  ntimes = {self.ntimes}\n')

        eids_device = self.element * 10 + self.device_code

        if self.is_sort1:
            struct1 = Struct(endian + b'i4f')
        else:
            raise NotImplementedError('SORT2')

        op2_ascii.write(f'nelements={nelements:d}\n')

        for itime in range(self.ntimes):
            self._write_table_3(op2_file, op2_ascii, new_result, itable, itime)

            # record 4
            itable -= 1
            header = [4, itable, 4,
                      4, 1, 4,
                      4, 0, 4,
                      4, ntotal, 4,
                      4 * ntotal]
            op2_file.write(pack('%ii' % len(header), *header))
            op2_ascii.write('r4 [4, 0, 4]\n')
            op2_ascii.write(f'r4 [4, {itable:d}, 4]\n')
            op2_ascii.write(f'r4 [4, {4 * ntotal:d}, 4]\n')

            axial = self.data[itime, :, 0]
            torsion = self.data[itime, :, 1]

            for eid_device, axiali, torsioni in zip(eids_device, axial, torsion):
                data = [eid_device, axiali.real, torsioni.real, axiali.imag, torsioni.imag]
                op2_ascii.write('  eid_device=%s data=%s\n' % (eid_device, tuple(data)))
                op2_file.write(struct1.pack(*data))

            itable -= 1
            header = [4 * ntotal,]
            op2_file.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable


class ComplexRodStressArray(ComplexRodArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexRodArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> list[str]:
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

    def get_headers(self) -> list[str]:
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
