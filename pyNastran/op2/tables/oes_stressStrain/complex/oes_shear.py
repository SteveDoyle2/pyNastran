import copy
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.result_objects.op2_objects import get_complex_times_dtype
from pyNastran.op2.result_objects.utils_pandas import build_dataframe_transient_header, build_pandas_transient_elements
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import (
    StressObject, StrainObject, OES_Object,
    oes_complex_data_code, get_scode,
    set_element_case, set_freq_case, set_complex_modes_case,
    slice_eids_by_index)
from pyNastran.f06.f06_formatting import write_imag_floats_13e


class ComplexShearArray(OES_Object):
    """
    Common class for:
     - ComplexShearStressArray
     - ComplexShearStrainArray
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)   ## why???
        #self.element_node = None
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        #self.itime = 0
        self.nelements = 0  # result specific

        if is_sort1:
            pass
        else:
            raise NotImplementedError('SORT2')

    def slice_by_element_id(self, eids: np.ndarray,
                            assume_exists: bool=False, inplace: bool=False):
        eids, ieid, neid2 = slice_eids_by_index(self.element, eids, assume_exists)
        obj = self if inplace else copy.deepcopy(self)
        obj.element = obj.element[ieid]
        obj.data = obj.data[:, ieid, :]
        obj.nelements = neid2
        return obj

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

    #def get_nnodes(self):
        #return get_nnodes(self)

    def build(self):
        """sizes the vectorized attributes of the ComplexShearArray"""
        if not hasattr(self, 'subtitle'):
            self.subtitle = self.data_code['subtitle']
        #print('ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (
            #self.ntimes, self.nelements, self.ntotal, self.subtitle))
        nnodes = 1

        #self.names = []
        #self.nelements //= nnodes
        self.nelements //= self.ntimes
        self.ntotal = self.nelements * nnodes# * 2
        #self.ntotal
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #print('ntotal=%s ntimes=%s nelements=%s' % (self.ntotal, self.ntimes, self.nelements))

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        idtype, cfdtype = get_complex_times_dtype(self.size)
        self._times = np.zeros(self.ntimes, dtype=self.analysis_fmt)
        #self.ntotal = self.nelements * nnodes

        self.element = np.zeros(self.nelements, dtype=idtype)

        # the number is messed up because of the offset for the element's properties
        if self.nelements != self.ntotal:
            msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (
                self.ntimes, self.nelements, nnodes, self.nelements * nnodes,
                self.ntotal)
            raise RuntimeError(msg)

        # [max_shear, avg_shear]
        self.data = np.zeros((self.ntimes, self.ntotal, 2), dtype=cfdtype)

    def build_dataframe(self):
        """creates a pandas dataframe"""
        #Mode                                          1                   2
        #EigenvalueReal                             -0.0                -0.0
        #EigenvalueImag                             -0.0                -0.0
        #Damping                                     0.0                 0.0
        #ElementID Item
        #22        max_shear  5.855954e-09+0.000000e+00j  0.000000+0.000000j
        #          avg_shear  5.855954e-09+0.000000e+00j  0.000000+0.000000j
        #import pandas as pd
        column_names, column_values = build_dataframe_transient_header(self)
        self.data_frame = build_pandas_transient_elements(
            self, column_values, column_names,
            self.headers, self.element, self.data)


    @classmethod
    def _add_case(cls,
                  table_name, element_name, isubcase,
                  is_sort1, is_random, is_msc,
                  random_code, title, subtitle, label):
        is_strain = 'Strain' in cls.__name__
        assert isinstance(is_strain, bool), is_strain
        assert isinstance(element_name, str), element_name
        is_strain = 'Strain' in cls.__name__
        num_wide = 3
        data_code = oes_complex_data_code(
            table_name,
            element_name, num_wide,
            is_sort1=is_sort1, is_random=is_random,
            random_code=random_code, title=title, subtitle=subtitle, label=label,
            is_msc=is_msc)

        ELEMENT_NAME_TO_ELEMENT_TYPE = {
            'CSHEAR': 4,
        }
        #
        #stress_bits[0] = 1  # Henky von-mises bit (vs. octehedral)????
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
    def add_freq_case(cls, table_name, element_name: str,
                      element, data, isubcase,
                      freqs,
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
                                     set_element_case, (element, data),
                                     modes, eigrs, eigis)
        return obj

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
                        (tx1, ty1) = t1
                        (tx2, ty2) = t2
                        d = t1 - t2
                        if not np.allclose([tx1.real, tx1.imag, ty1.real, ty1.imag],
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

    def add_sort1(self, dt, eid, max_shear, avg_shear):
        """unvectorized method for adding SORT1 transient data"""
        assert self.sort_method == 1, self
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.data[self.itime, self.itotal] = [max_shear, avg_shear]
        self.element[self.itotal] = eid
        #self.ielement += 1
        self.itotal += 1

    def get_stats(self, short: bool=False) -> list[str]:
        if not self.is_built:
            return [
                f'<{self.__class__.__name__}>; table_name={self.table_name!r}\n',
                f'  ntimes: {self.ntimes:d}\n',
                f'  ntotal: {self.ntotal:d}\n',
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        nnodes = self.element.shape[0]
        #ntotal = self.ntotal
        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i nnodes=%i; table_name=%r\n' % (
                self.__class__.__name__, ntimes, nelements, nnodes, self.table_name))
        else:
            msg.append('  type=%s nelements=%i nnodes=%i; table_name=%r\n' % (
                self.__class__.__name__, nelements, nnodes, self.table_name))
        msg.append('  data: [ntimes, nnodes, 2] where 2=[%s]\n' % str(', '.join(self._get_headers())))
        msg.append(f'  element.shape = {self.element.shape}\n')
        msg.append(f'  data.shape = {self.data.shape}\n')
        msg.append('  %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        """
                C O M P L E X   F O R C E S   A C T I N G   O N   S H E A R   P A N E L   E L E M E N T S   (CSHEAR)
                                                          (REAL/IMAGINARY)

                  ====== POINT  1 ======      ====== POINT  2 ======      ====== POINT  3 ======      ====== POINT  4 ======
 ELEMENT          F-FROM-4      F-FROM-2      F-FROM-1      F-FROM-3      F-FROM-2      F-FROM-4      F-FROM-3      F-FROM-1
         ID               KICK-1       SHEAR-12       KICK-2       SHEAR-23       KICK-3       SHEAR-34       KICK-4       SHEAR-41
            28  0.0           0.0           0.0           0.0           0.0           0.0           0.0           0.0
                0.0           0.0           0.0           0.0           0.0           0.0           0.0           0.0
                       0.0           0.0           0.0           0.0           0.0           0.0           0.0           0.0
                       0.0           0.0           0.0           0.0           0.0           0.0           0.0           0.0
        """
        #if header is None:
            #header = []
        #f.write(self.code_information())
        #return page_num
        msg_temp = _get_cshear_msg(is_mag_phase, is_sort1)

        ntimes = self.data.shape[0]
        eids = self.element
        if self.is_sort1:
            if is_sort1:
                for itime in range(ntimes):
                    dt = self._times[itime]

                    dt_line = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
                    header[1] = dt_line
                    msg = header + msg_temp
                    f06_file.write('\n'.join(msg))

                    max_shear = self.data[itime, :, 0]
                    avg_shear = self.data[itime, :, 1]
                    assert len(eids) == len(max_shear)
                    assert len(max_shear) > 0, max_shear
                    for eid, max_sheari, avg_sheari in zip(eids, max_shear, avg_shear):
                        assert isinstance(eid, integer_types), 'eid=%s type=%s' % (eid, type(eid))
                        [rmax_shear, imax_shear, ravg_shear, iavg_shear
                         ,] = write_imag_floats_13e([max_sheari, avg_sheari], is_mag_phase)

                        #f.write('                      28                  0.0          /  0.0                           0.0          /  0.0\n')
                        f06_file.write(
                            '%24s                 %-13s / %-13s                 %-13s / %-13s\n' % (
                                eid, rmax_shear, imax_shear, ravg_shear, iavg_shear))
                    f06_file.write(page_stamp % page_num)
                    page_num += 1
            else:
                # TODO: write in SORT2
                times = self._times
                for ieid, eid in enumerate(eids):
                    max_shear = self.data[:, ieid, 0].ravel()
                    avg_shear = self.data[:, ieid, 1].ravel()
                    for itime, max_sheari, avg_sheari in zip(times, max_shear, avg_shear):
                        [rmax_shear, imax_shear, ravg_shear, iavg_shear
                         ] = write_imag_floats_13e([max_sheari, avg_sheari], is_mag_phase)

                        #f06_file.write(
                            #'   %6s   %-13s / %-13s     %-13s / %-13s\n' % (
                                #eid, rmax_shear, imax_shear, ravg_shear, iavg_shear))
                        f06_file.write(
                            '%24s                 %-13s / %-13s                 %-13s / %-13s\n' % (
                                eid, rmax_shear, imax_shear, ravg_shear, iavg_shear))
                    f06_file.write(page_stamp % page_num)
                    page_num += 1
        else:
            raise NotImplementedError('ComplexShearArray-sort2')
        return page_num - 1

    def write_op2(self, op2_file, op2_ascii, itable, new_result, date,
                  is_mag_phase=False, endian: bytes=b'>'):
        """
        writes an OP2

        record = [eid_device, max_shear_r, max_shear_i, avg_shear_r, avg_shear_i]
        (interleaved; see oes_cshear_complex_5)
        """
        import inspect
        from struct import Struct
        from pyNastran.op2.op2_interface.write_utils import view_idtype_as_fdtype
        from pyNastran.op2.writer.utils import fdtype_from_data
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write(f'{self.__class__.__name__}.write_op2: {call_frame[1][3]}\n')
        if not self.is_sort1:  # pragma: no cover
            raise NotImplementedError('SORT2')

        if itable == -1:
            self._write_table_header(op2_file, op2_ascii, date)
            itable = -3

        nelements = self.data.shape[1]
        ntotali = self.num_wide
        assert ntotali == 5, ntotali
        ntotal = ntotali * nelements
        op2_ascii.write(f'  ntimes = {self.ntimes}\n')
        op2_ascii.write(f'nelements={nelements:d}\n')

        eids_device = self.element * 10 + self.device_code
        fdtype = fdtype_from_data(self.data.real, self.size)
        eid_device_floats = view_idtype_as_fdtype(eids_device, fdtype)

        struct_13i = Struct(endian + b'13i')
        struct_i = Struct(endian + b'i')
        for itime in range(self.ntimes):
            self._write_table_3(op2_file, op2_ascii, new_result, itable, itime)

            # record 4
            itable -= 1
            header = [4, itable, 4,
                      4, 1, 4,
                      4, 0, 4,
                      4, ntotal, 4,
                      4 * ntotal]
            op2_file.write(struct_13i.pack(*header))
            op2_ascii.write('r4 [4, 0, 4]\n')
            op2_ascii.write(f'r4 [4, {itable:d}, 4]\n')
            op2_ascii.write(f'r4 [4, {4 * ntotal:d}, 4]\n')

            max_shear = self.data[itime, :, 0]
            avg_shear = self.data[itime, :, 1]
            datai = np.column_stack([
                eid_device_floats,
                max_shear.real, max_shear.imag,
                avg_shear.real, avg_shear.imag]).astype(fdtype, copy=False)
            assert ntotal == datai.size
            op2_file.write(np.ascontiguousarray(datai).tobytes())

            itable -= 1
            op2_file.write(struct_i.pack(4 * ntotal))
            op2_ascii.write('footer = %s\n' % [4 * ntotal])
            new_result = False
        return itable

    @property
    def headers(self) -> list[str]:
        return self._get_headers()

    @property
    def headers(self) -> list[str]:
        return self.headers

class ComplexShearStressArray(ComplexShearArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexShearArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    @property
    def headers(self) -> list[str]:
        return ['max_shear', 'avg_shear']


class ComplexShearStrainArray(ComplexShearArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexShearArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)
        assert self.is_strain, self.stress_bits

    @property
    def headers(self) -> list[str]:
        return ['max_shear', 'avg_shear']

def _get_cshear_msg(is_mag_phase, is_sort1):
    if is_mag_phase:
        raise NotImplementedError()
    else:
        out = [
            '                         C O M P L E X   S T R E S S E S   I N   S H E A R P A N E L S   ( C S H E A R )\n'
            '                                                          (REAL/IMAGINARY)\n'
            ' \n'
            '                 ELEMENT                            MAXIMUM                                       AVERAGE\n'
            '                   ID.                               SHEAR                                         SHEAR\n'
            #'                      28                  9.089907E+01 /  1.709346E+02                 -9.089907E+01 / -1.709346E+02\n'
        ]
    return out
        #out = [
            #'                C O M P L E X   F O R C E S   A C T I N G   O N   S H E A R   P A N E L   E L E M E N T S   (CSHEAR)\n',
            #'                                                          (REAL/IMAGINARY)\n',
            #' \n',
            #'                  ====== POINT  1 ======      ====== POINT  2 ======      ====== POINT  3 ======      ====== POINT  4 ======\n',
            #' ELEMENT          F-FROM-4      F-FROM-2      F-FROM-1      F-FROM-3      F-FROM-2      F-FROM-4      F-FROM-3      F-FROM-1\n',
            #'         ID               KICK-1       SHEAR-12       KICK-2       SHEAR-23       KICK-3       SHEAR-34       KICK-4       SHEAR-41\n',
        #]
    #else:
        #out = [
            #'                C O M P L E X   F O R C E S   A C T I N G   O N   S H E A R   P A N E L   E L E M E N T S   (CSHEAR)\n',
            #'                                                          (REAL/IMAGINARY)\n',
            #' \n',
            #'                  ====== POINT  1 ======      ====== POINT  2 ======      ====== POINT  3 ======      ====== POINT  4 ======\n',
            #' ELEMENT          F-FROM-4      F-FROM-2      F-FROM-1      F-FROM-3      F-FROM-2      F-FROM-4      F-FROM-3      F-FROM-1\n',
            #'         ID               KICK-1       SHEAR-12       KICK-2       SHEAR-23       KICK-3       SHEAR-34       KICK-4       SHEAR-41\n',
        #]
    #return out
