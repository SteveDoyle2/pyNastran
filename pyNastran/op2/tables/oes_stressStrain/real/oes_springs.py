from itertools import count
from typing import List

import numpy as np
from numpy import zeros

from pyNastran.utils.numpy_utils import integer_types, float_types
from pyNastran.op2.result_objects.op2_objects import get_times_dtype
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import (
    StressObject, StrainObject, OES_Object, oes_data_code)
from pyNastran.f06.f06_formatting import write_float_13e, _eigenvalue_header

ELEMENT_NAME_TO_ELEMENT_TYPE = {
    'CELAS1': 11,
    'CELAS2': 12,
    'CELAS3': 13,
    'CELAS4': 14,
}

class RealSpringArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)

        self.nelements = 0  # result specific

        self.itime = 0
        self.itotal = 0
        self.ielement = 0
        self.element = None

    @classmethod
    def add_static_case(cls, table_name, element_name, element, data, isubcase,
                        is_sort1=True, is_random=False, is_stress=True, is_msc=True,
                        random_code=0, title='', subtitle='', label=''):

        analysis_code = 1 # static
        data_code = oes_data_code(table_name, analysis_code,
                                  is_sort1=is_sort1, is_random=is_random,
                                  random_code=random_code,
                                  title=title, subtitle=subtitle, label=label,
                                  is_msc=is_msc)
        data_code['lsdvmns'] = [0] # TODO: ???
        data_code['data_names'] = []

        # I'm only sure about the 1s in the strains and the
        # corresponding 0s in the stresses.
        if is_stress:
            data_code['stress_bits'] = [0, 0, 0, 0]
            data_code['s_code'] = 0
        else:
            data_code['stress_bits'] = [0, 1, 0, 1]
            data_code['s_code'] = 1 # strain?

        element_type = ELEMENT_NAME_TO_ELEMENT_TYPE[element_name]
        data_code['element_name'] = element_name
        data_code['element_type'] = element_type
        data_code['load_set'] = 1

        ntimes = data.shape[0]
        nnodes = data.shape[1]
        dt = None
        obj = cls(data_code, is_sort1, isubcase, dt)
        obj.element = element
        obj.data = data

        obj.ntimes = ntimes
        obj.ntotal = nnodes
        obj._times = [None]
        obj.is_built = True
        return obj

    @property
    def is_real(self) -> bool:
        return True

    @property
    def is_complex(self) -> bool:
        return False

    @property
    def nnodes_per_element(self) -> int:
        return 1

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        raise NotImplementedError()

    #def __mul__(self, factor):
        #"""in-place multiplication"""
        #assert isinstance(factor, float_types), 'factor=%s and must be a float' % (factor)
        #self.data *= factor
    #def __rmul__(self, factor):
        #assert isinstance(factor, float_types), 'factor=%s and must be a float' % (factor)
        #self.data *= factor

    #def __sub__(self, factor):
        #if isinstance(factor, float_types):
            #self.data -= factor
        #else:
            ## TODO: should support arrays
            #raise TypeError('factor=%s and must be a float' % (factor))
    #def __add__(self, factor):
        #"""[C] = [A] + b"""
        #if isinstance(factor, float_types):
            #self.data += factor
        #else:
            ## TODO: should support arrays
            #raise TypeError('factor=%s and must be a float' % (factor))

    #def __radd__(self, factor):
        #"""[C] = b + [A]"""
        #return self.__add__(factor)

    def update_data_components(self):
        pass

    def __iadd__(self, factor):
        """[A] += b"""
        if isinstance(factor, float_types):
            self.data += factor
        else:
            # TODO: should support arrays
            raise TypeError('factor=%s and must be a float' % (factor))
        self.update_data_components()

    def __isub__(self, factor):
        """[A] -= b"""
        if isinstance(factor, float_types):
            self.data -= factor
        else:
            # TODO: should support arrays
            raise TypeError('factor=%s and must be a float' % (factor))
        self.update_data_components()

    def __imul__(self, factor):
        """[A] *= b"""
        assert isinstance(factor, float_types), 'factor=%s and must be a float' % (factor)
        self.data *= factor
        self.update_data_components()

    def __idiv__(self, factor):
        """[A] *= b"""
        assert isinstance(factor, float_types), 'factor=%s and must be a float' % (factor)
        self.data *= 1. / factor
        self.update_data_components()

    #def linear_combination(a, coeffs):
        #import numexpr as ne
        #local_vars = locals()
        #letters = [
            #'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',
            #'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z']
        #expr = ''
        #for ai, coeff, letter in zip(a, coeffs, letters):
            #expr += '%s*%s' % (coeff, letter)
            #local_vars[letter] = ai.data
        #c = ne.evaluate(expr)
        #return c

    def build(self):
        """sizes the vectorized attributes of the RealSpringArray"""
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
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        float_fmt = 'float32' if self.size == 4 else 'float64'
        dtype = float_fmt
        if isinstance(self.nonlinear_factor, integer_types):
            int_fmt = 'int32' if self.size == 4 else 'int64'
            dtype = int_fmt
        self.build_data(self.ntimes, self.nelements, dtype)

    def build_data(self, ntimes, nelements, dtype):
        """actually performs the build step"""
        self.ntimes = ntimes
        self.nelements = nelements
        _times = zeros(ntimes, dtype=dtype)
        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size)
        element = zeros(nelements, dtype=idtype)

        #[stress]
        data = zeros((ntimes, nelements, 1), dtype=fdtype)

        if self.load_as_h5:
            #for key, value in sorted(self.data_code.items()):
                #print(key, value)
            group = self._get_result_group()
            self._times = group.create_dataset('_times', data=_times)
            self.element = group.create_dataset('element', data=element)
            self.data = group.create_dataset('data', data=data)
        else:
            self._times = _times
            self.element = element
            self.data = data

    def build_dataframe(self):
        """creates a pandas dataframe

        v 0.24
        Static                     0
        ElementID Item
        30        spring_stress  0.0
        31        spring_stress  0.0
        32        spring_stress  0.0
        33        spring_stress  0.0

        v 0.25 for test_bdf_op2_elements_01
        Static  ElementID  spring_stress
        0              30            0.0
        1              31            0.0
        2              32            0.0
        3              33            0.0
        ...
        """
        import pandas as pd

        headers = self.get_headers()
        if self.nonlinear_factor not in (None, np.nan):
            # Mode                                1             2             3
            # Freq                     1.482246e-10  3.353940e-09  1.482246e-10
            # Eigenvalue              -8.673617e-19  4.440892e-16  8.673617e-19
            # Radians                  9.313226e-10  2.107342e-08  9.313226e-10
            # ElementID Item
            # 30        spring_stress           0.0          -0.0          -0.0
            # 31        spring_stress           0.0          -0.0          -0.0
            # 32        spring_stress           0.0           0.0           0.0
            # 33        spring_stress           0.0           0.0           0.0
            column_names, column_values = self._build_dataframe_transient_header()
            data_frame = self._build_pandas_transient_elements(
                column_values, column_names,
                headers, self.element, self.data)
        else:
            #Static     spring_stress
            #ElementID
            #30                   0.0
            #31                   0.0
            #32                   0.0
            #33                   0.0
            data_frame = pd.DataFrame(self.data[0], columns=headers, index=self.element)
            data_frame.index.name = 'ElementID'
            data_frame.columns.names = ['Static']
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
                        (force1, stress1) = t1
                        (force2, stress2) = t2
                        if not np.allclose(t1, t2):
                        #if not np.array_equal(t1, t2):
                            msg += '%s\n  (%s, %s)\n  (%s, %s)\n' % (
                                eid,
                                force1, stress1,
                                force2, stress2)
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

    def add_sort1(self, dt, eid, stress):
        self._times[self.itime] = dt
        #if self.itime == 0:
        #print('itime=%s eid=%s' % (self.itime, eid))
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [stress]
        self.ielement += 1

    def get_stats(self, short=False) -> List[str]:
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        #print(self.data.shape[:1])
        #ntimes, nelements = self.data.shape[:1]
        ntimes = self.data.shape[0]
        nelements = self.data.shape[1]
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
        # elements are always sorted; nodes are not
        itot = np.searchsorted(eids, self.element)  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.element_node[:, 0] == eid) for eid in eids])
        ind = np.searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = self.get_f06_header(is_mag_phase)

        if self.is_sort1:
            page_num = self._write_sort1_as_sort1(header, page_stamp, page_num, f06_file, msg_temp)
        else:
            raise NotImplementedError(self.code_information())
            #page_num = self._write_sort2_as_sort2(header, page_stamp, page_num, f06_file, msg_temp)
        return page_num

    def _write_sort1_as_sort1(self, header, page_stamp, page_num, f06_file, msg_temp):
        ntimes = self.data.shape[0]

        eids = self.element
        nwrite = len(eids)
        nrows = nwrite // 4
        nleftover = nwrite - nrows * 4

        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))
            stress = self.data[itime, :, 0]

            out = []
            for eid, stressi in zip(eids, stress):
                out.append([eid, write_float_13e(stressi)])

            for i in range(0, nrows * 4, 4):
                f06_file.write('    %10i  %13s    %10i  %13s    %10i  %13s    %10i  %13s\n' % (
                    tuple(out[i] + out[i + 1] + out[i + 2] + out[i + 3])))

            i = nrows * 4
            if nleftover == 3:
                f06_file.write('    %10i  %13s    %10i  %13s    %10i  %13s\n' % (
                    tuple(out[i] + out[i + 1] + out[i + 2])))
            elif nleftover == 2:
                f06_file.write('    %10i  %13s    %10i  %13s\n' % (
                    tuple(out[i] + out[i + 1])))
            elif nleftover == 1:
                f06_file.write('    %10i  %13s\n' % tuple(out[i]))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

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
        op2_ascii.write('  ntimes = %s\n' % self.ntimes)

        eids_device = self.element * 10 + self.device_code

        #print('ntotal=%s' % (ntotal))
        #assert ntotal == 193, ntotal

        if self.is_sort1:
            struct1 = Struct(endian + b'if')
        else:
            raise NotImplementedError('SORT2')

        op2_ascii.write('%s-nelements=%i\n' % (self.element_name, nelements))
        for itime in range(self.ntimes):
            self._write_table_3(op2, op2_ascii, new_result, itable, itime)

            # record 4
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

            stress = self.data[itime, :, 0]

            for eid, stressi in zip(eids_device, stress):
                data = [eid, stressi]
                op2_ascii.write('  eid=%s force=%s\n' % tuple(data))
                op2.write(struct1.pack(*data))

            itable -= 1
            header = [4 * ntotal,]
            op2.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable


class RealSpringStressArray(RealSpringArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealSpringArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        headers = ['spring_stress']
        return headers

    def get_f06_header(self, is_mag_phase=True):
        if self.element_type == 11:  # CELAS1
            msg = ['                              S T R E S S E S   I N   S C A L A R   S P R I N G S        ( C E L A S 1 )\n']
        elif self.element_type == 12:  # CELAS2
            msg = ['                              S T R E S S E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )\n']
        elif self.element_type == 13:  # CELAS3
            msg = ['                              S T R E S S E S   I N   S C A L A R   S P R I N G S        ( C E L A S 3 )\n']
        elif self.element_type == 14:  # CELAS4
            msg = ['                              S T R E S S E S   I N   S C A L A R   S P R I N G S        ( C E L A S 4 )\n']
        else:
            raise NotImplementedError('element_name=%s element_type=%s' % (self.element_name, self.element_type))

        msg += [
            '      ELEMENT         STRESS           ELEMENT         STRESS           ELEMENT         STRESS           ELEMENT         STRESS\n'
            '        ID.                              ID.                              ID.                              ID.\n'
        ]
        return msg



class RealSpringStrainArray(RealSpringArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealSpringArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        headers = ['spring_strain']
        return headers

    def get_f06_header(self, is_mag_phase=True):
        if self.element_type == 11:  # CELAS1
            msg = ['                               S T R A I N S    I N   S C A L A R   S P R I N G S        ( C E L A S 1 )\n']
        elif self.element_type == 12:  # CELAS2
            msg = ['                               S T R A I N S    I N   S C A L A R   S P R I N G S        ( C E L A S 2 )\n']
        elif self.element_type == 13:  # CELAS3
            msg = ['                               S T R A I N S    I N   S C A L A R   S P R I N G S        ( C E L A S 3 )\n']
        elif self.element_type == 14:  # CELAS4
            msg = ['                               S T R A I N S    I N   S C A L A R   S P R I N G S        ( C E L A S 4 )\n']
        else:
            raise NotImplementedError('element_name=%s element_type=%s' % (self.element_name, self.element_type))

        msg += [
            '      ELEMENT         STRAIN           ELEMENT         STRAIN           ELEMENT         STRAIN           ELEMENT         STRAIN\n'
            '        ID.                              ID.                              ID.                              ID.\n'
        ]
        return msg


class RealNonlinearSpringStressArray(OES_Object):
    """
    ::

      #ELEMENT-ID =     102
                               #N O N L I N E A R   S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )
        #TIME          AXIAL STRESS         EQUIVALENT         TOTAL STRAIN       EFF. STRAIN          EFF. CREEP        LIN. TORSIONAL
                                             #STRESS                             PLASTIC/NLELAST          STRAIN              STRESS
      #2.000E-02        1.941367E+01        1.941367E+01        1.941367E-04        0.0                 0.0                 0.0
      #3.000E-02        1.941367E+01        1.941367E+01        1.941367E-04        0.0                 0.0                 0.0
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=True)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        self.nelements = 0  # result specific

        self.itime = 0
        self.itotal = 0
        self.ielement = 0
        self.element = None

        if is_sort1:
            pass
        else:
            raise NotImplementedError('SORT2')

    @property
    def is_real(self):
        return True

    @property
    def is_complex(self):
        return False

    @property
    def is_stress(self):
        return True

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def _get_msgs(self):
        raise NotImplementedError()

    def get_headers(self) -> List[str]:
        headers = ['force', 'stress']
        return headers

    def build(self):
        """sizes the vectorized attributes of the RealNonlinearSpringStressArray"""
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
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'
        _times = zeros(self.ntimes, dtype=dtype)
        element = zeros(self.nelements, dtype='int32')

        #[force, stress]
        data = zeros((self.ntimes, self.nelements, 2), dtype='float32')

        if self.load_as_h5:
            #for key, value in sorted(self.data_code.items()):
                #print(key, value)
            group = self._get_result_group()
            self._times = group.create_dataset('_times', data=_times)
            self.element = group.create_dataset('element', data=element)
            self.data = group.create_dataset('data', data=data)
        else:
            self._times = _times
            self.element = element
            self.data = data

    def __eq__(self, table):  # pragma: no cover
        self._eq_header(table)
        assert self.is_sort1 == table.is_sort1
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            i = 0
            if self.is_sort1:
                for itime in range(ntimes):
                    for ieid, eid, in enumerate(self.element):
                        t1 = self.data[itime, ieid, :]
                        t2 = table.data[itime, ieid, :]
                        (force1, stress1) = t1
                        (force2, stress2) = t2
                        if not np.allclose(t1, t2):
                        #if not np.array_equal(t1, t2):
                            msg += '%s\n  (%s, %s)\n  (%s, %s)\n' % (
                                eid,
                                force1, stress1,
                                force2, stress2)
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

    def add_sort1(self, dt, eid, force, stress):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [force, stress]
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
            msg.append('  type=%s ntimes=%i nelements=%i\n'
                       % (self.__class__.__name__, ntimes, nelements))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i\n'
                       % (self.__class__.__name__, nelements))
            ntimes_word = '1'
        msg.append('  eType\n')
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        if self.is_sort1:
            if self.element_type == 224:
                nspring = 1 # CELAS1
            elif self.element_type == 225:
                nspring = 3 # CELAS3
            else:
                raise NotImplementedError('type=%s name=%s' % (self.element_type, self.element_name))

            msg = [
                '          N O N L I N E A R   F O R C E S  A N D  S T R E S S E S  I N   S C A L A R   S P R I N G S    ( C E L A S %i )\n'
                ' \n'
                '         ELEMENT-ID          FORCE         STRESS                    ELEMENT-ID          FORCE         STRESS\n' % nspring
                #'         5.000000E-02     2.000000E+01   1.000000E+01                1.000000E-01     4.000000E+01   2.000000E+01'
            ]
            page_num = self._write_sort1_as_sort1(header, page_stamp, page_num, f06_file, msg)
        else:
            msg = [
                '          N O N L I N E A R   F O R C E S  A N D  S T R E S S E S  I N   S C A L A R   S P R I N G S    ( C E L A S %i )\n'
                ' \n'
                '             STEP            FORCE         STRESS                        STEP            FORCE         STRESS\n' % nspring
                #'         5.000000E-02     2.000000E+01   1.000000E+01                1.000000E-01     4.000000E+01   2.000000E+01'
            ]
            raise NotImplementedError('RealNonlinearSpringStressArray-sort2')
        return page_num

    def _write_sort1_as_sort1(self, header, page_stamp, page_num, f06_file, msg_temp):
        """
        ::
              ELEMENT-ID =      20
                  N O N L I N E A R   F O R C E S  A N D  S T R E S S E S  I N   S C A L A R   S P R I N G S    ( C E L A S 1 )

                     STEP            FORCE         STRESS                        STEP            FORCE         STRESS
                 5.000000E-02     2.000000E+01   1.000000E+01                1.000000E-01     4.000000E+01   2.000000E+01
                 1.500000E-01     6.000000E+01   3.000000E+01                2.000000E-01     8.000000E+01   4.000000E+01
        """
        ntimes = self.data.shape[0]

        eids = self.element
        neids = len(eids)
        is_odd = neids % 2 == 1
        if is_odd:
            neids -= 1

        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))
            force = self.data[itime, :, 0]
            stress = self.data[itime, :, 1]
            for i, eid, forcei, stressi, in zip(count(step=2), eids[:neids:2], force[:neids:2], stress[:neids:2]):
                f06_file.write(
                    '         %-13i   %-13s  %-13s                %-13s   %-13s  %s\n' % (
                        eid,
                        write_float_13e(forcei),
                        write_float_13e(stressi),
                        eids[i + 1],
                        write_float_13e(force[i + 1]),
                        write_float_13e(stress[i + 1])
                ))
            if is_odd:
                f06_file.write('         %-13i   %-13s  %s\n' % (
                    eids[neids],
                    write_float_13e(force[neids]),
                    write_float_13e(stress[neids])
                ))

            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1
