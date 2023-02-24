import numpy as np
from numpy import zeros

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.result_objects.op2_objects import get_complex_times_dtype
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import (
    StressObject, StrainObject, OES_Object, get_scode,
    oes_complex_data_code, set_freq_case, set_complex_modes_case,
    set_element_case)
from pyNastran.f06.f06_formatting import write_imag_floats_13e, _eigenvalue_header

ELEMENT_NAME_TO_ELEMENT_TYPE = {
    'CELAS1': 11,
    'CELAS2': 12,
    'CELAS3': 13,
    'CELAS4': 14,
}

class ComplexSpringDamperArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        if not is_sort1:
            raise NotImplementedError('SORT2')

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

    #def get_headers(self):
        #headers = ['axial', 'torque']
        #return headers

    def build(self):
        """sizes the vectorized attributes of the ComplexSpringDamperArray"""
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        self.nelements //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        idtype, cfdtype = get_complex_times_dtype(self.size)
        self._times = zeros(self.ntimes, dtype=self.analysis_fmt)
        self.element = zeros(self.nelements, dtype=idtype)

        #[spring_stress]
        self.data = zeros((self.ntimes, self.ntotal, 1), dtype=cfdtype)

    def build_dataframe(self):
        """creates a pandas dataframe"""
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
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
        num_wide = 3
        data_code = oes_complex_data_code(
            table_name,
            element_name, num_wide,
            is_sort1=is_sort1, is_random=is_random,
            random_code=random_code, title=title, subtitle=subtitle, label=label,
            is_msc=is_msc)

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
        data_code['element_name'] = element_name.upper()
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
        return obj

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ielem, eid in enumerate(self.element):
                    t1 = self.data[itime, ielem, :]
                    t2 = table.data[itime, ielem, :]
                    if not np.array_equal(t1, t2):
                        msg += '%s    (%s, %s)  (%s, %s)\n' % (
                            eid,
                            t1.real, t1.imag,
                            t2.real, t2.imag)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, stress):
        """unvectorized method for adding SORT1 transient data"""
        assert self.sort_method == 1, self
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, 0] = stress
        self.ielement += 1

    def get_stats(self, short: bool=False) -> list[str]:
        if not self.is_built:
            return [
                f'<{self.__class__.__name__}>; table_name={self.table_name!r}\n',
                f'  ntimes: {self.ntimes:d}\n',
                f'  ntotal: {self.ntotal:d}\n',
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        #ntotal = self.ntotal

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
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (
            ntimes_word, n, n, str(', '.join(headers))))
        msg.append(f'  element.shape = {self.element.shape}\n')
        msg.append(f'  data.shape = {self.data.shape}\n')
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True, is_sort1=True):
        # 11-CELAS1, 12-CELAS2, 13-CELAS3, 14-CELAS4

        #'            FREQUENCY                   STRESS                        FREQUENCY                   STRESS'
        if self.element_type == 11:
            msg = ['                       C O M P L E X   S T R E S S E S   I N   S C A L A R   S P R I N G S   ( C E L A S 1 )\n']
        elif self.element_type == 12:
            msg = ['                       C O M P L E X   S T R E S S E S   I N   S C A L A R   S P R I N G S   ( C E L A S 2 )\n']
        elif self.element_type == 13:
            msg = ['                       C O M P L E X   S T R E S S E S   I N   S C A L A R   S P R I N G S   ( C E L A S 3 )\n']
        elif self.element_type == 14:
            msg = ['                       C O M P L E X   S T R E S S E S   I N   S C A L A R   S P R I N G S   ( C E L A S 4 )\n']
        #elif self.element_type == 20: # CDAMP1
            #msg = ['                         C O M P L E X   F O R C E S   I N   S C A L A R   D A M P E R S   ( C D A M P 1 )\n']
        #elif self.element_type == 21: # CDAMP2
            #msg = ['                         C O M P L E X   F O R C E S   I N   S C A L A R   D A M P E R S   ( C D A M P 2 )\n']
        else:
            raise NotImplementedError('element_name=%s element_type=%s' % (self.element_name, self.element_type))

        if is_mag_phase:
            msg += ['                                                          (MAGNITUDE/PHASE)\n \n']
        else:
            msg += ['                                                          (REAL/IMAGINARY)\n \n']

        if is_sort1:
            msg += [
                '                ELEMENT                                                   ELEMENT\n'
                '                  ID.                    STRESS                             ID.                    STRESS\n'
            ]
            #'                      14                  0.0          /  0.0                           0.0          /  0.0'
        else:
            msg += ['            FREQUENCY                    STRESS                       FREQUENCY                    STRESS\n']
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        if header is None:
            header = []
        msg_temp = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        ntimes = self.data.shape[0]

        eids = self.element
        #is_odd = False
        #nwrite = len(eids)
        #if len(eids) % 2 == 1:
            #nwrite -= 1
            #is_odd = True

        #print('len(eids)=%s nwrite=%s is_odd=%s' % (len(eids), nwrite, is_odd))
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            spring_force = self.data[itime, :, 0]

            for eid, spring_forcei in zip(eids, spring_force):
                [rspring, ispring] = write_imag_floats_13e([spring_forcei], is_mag_phase)
                #ELEMENT                             AXIAL                                       TORSIONAL
                    #ID.                              STRESS                                         STRESS
                    #14                  0.0          /  0.0                           0.0          /  0.0


                f06_file.write('      %8i   %-13s / %-13s\n' % (eid, rspring, ispring))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def write_op2(self, op2_file, op2_ascii, itable, new_result,
                  date, is_mag_phase=False, endian='>'):
        """writes an OP2"""
        import inspect
        from struct import Struct, pack
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write(f'{self.__class__.__name__}.write_op2: {call_frame[1][3]}\n')

        if itable == -1:
            self._write_table_header(op2_file, op2_ascii, date)
            itable = -3

        #eids = self.element

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

        #device_code = self.device_code
        op2_ascii.write(f'  ntimes = {self.ntimes}\n')

        eids_device = self.element * 10 + self.device_code

        #print('ntotal=%s' % (ntotal))
        #assert ntotal == 193, ntotal

        if self.is_sort1:
            struct1 = Struct(endian + b'i2f')
        else:
            raise NotImplementedError('SORT2')

        op2_ascii.write('%s-nelements=%i\n' % (self.element_name, nelements))
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

            from pyNastran.op2.op2_interface.utils import to_mag_phase
            stress = self.data[itime, :, 0]
            reals, imags = to_mag_phase(stress, is_mag_phase)

            for eid, stress_real, stress_imag in zip(eids_device, reals, imags):
                data = [eid, stress_real, stress_imag]
                op2_ascii.write(f'  eid={eid} stress={[stress_real, stress_imag]}\n')
                op2_file.write(struct1.pack(*data))

            itable -= 1
            header = [4 * ntotal,]
            op2_file.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable


class ComplexSpringStressArray(ComplexSpringDamperArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexSpringDamperArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> list[str]:
        headers = ['spring_stress']
        return headers


class ComplexSpringStrainArray(ComplexSpringDamperArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexSpringDamperArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> list[str]:
        headers = ['spring_strain']
        return headers
