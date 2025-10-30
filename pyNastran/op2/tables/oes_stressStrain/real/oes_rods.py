from typing import TextIO
import numpy as np
from numpy import zeros, searchsorted, allclose

from pyNastran.op2.result_objects.op2_objects import get_times_dtype
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import (
    StressObject, StrainObject, OES_Object, oes_real_data_code,
    set_element_case, set_static_case, set_modal_case,
    set_transient_case, set_post_buckling_case)
from pyNastran.op2.op2_interface.write_utils import view_dtype, view_idtype_as_fdtype
from pyNastran.f06.f06_formatting import (
    write_floats_13e, write_floats_13e_long,
    _eigenvalue_header) #, get_key0

ELEMENT_NAME_TO_ELEMENT_TYPE = {
    'CROD' : 1,
    'CTUBE' : 3,
    'CONROD' : 10,
}

# there is a bug for mode_solid_shell_bar.op2 for multiple times
class RealRodArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        self.nelements = 0  # result specific

        self.itime = 0
        self.itotal = 0
        self.ielement = 0
        self.element = None

    @classmethod
    def _add_case(cls,
                  table_name, element_name, isubcase,
                  is_sort1, is_random, is_msc,
                  random_code, title, subtitle, label):
        num_wide = 3
        is_strain = 'Strain' in cls.__name__
        data_code = oes_real_data_code(table_name,
                                       element_name, num_wide,
                                       is_sort1=is_sort1, is_random=is_random,
                                       random_code=random_code,
                                       title=title, subtitle=subtitle, label=label,
                                       is_msc=is_msc)

        # I'm only sure about the 1s in the strains and the
        # corresponding 0s in the stresses.
        #stress / strain -> 1, 3
        if is_strain:
            stress_bits = [0, 1, 0, 1]
            s_code = 1
        else:
            stress_bits = [0, 0, 0, 0]
            s_code = 0
            assert stress_bits[1] == 0, 'stress_bits=%s' % (stress_bits)

        # stress
        assert stress_bits[1] == stress_bits[3], 'stress_bits=%s' % (stress_bits)
        data_code['stress_bits'] = stress_bits
        data_code['s_code'] = s_code

        element_type = ELEMENT_NAME_TO_ELEMENT_TYPE[element_name]
        data_code['element_name'] = element_name
        data_code['element_type'] = element_type
        return data_code

    @classmethod
    def add_static_case(cls, table_name, element_name, element, data, isubcase,
                        is_sort1=True, is_random=False, is_msc=True,
                        random_code=0, title='', subtitle='', label=''):
        data_code = cls._add_case(
            table_name, element_name,
            isubcase, is_sort1, is_random, is_msc,
            random_code, title, subtitle, label)
        obj = set_static_case(cls, is_sort1, isubcase, data_code,
                              set_element_case, (element, data))
        return obj

    @classmethod
    def add_modal_case(cls, table_name, element_name: str, element, data, isubcase,
                       modes, eigns, cycles,
                       is_sort1=True, is_random=False, is_msc=True,
                       random_code=0, title='', subtitle='', label=''):
        data_code = cls._add_case(
            table_name, element_name,
            isubcase, is_sort1, is_random, is_msc,
            random_code, title, subtitle, label)
        obj = set_modal_case(cls, is_sort1, isubcase, data_code,
                             set_element_case, (element, data),
                             modes, eigns, cycles)
        return obj

    @classmethod
    def add_transient_case(cls, table_name, element_name, element, data, isubcase,
                           times,
                           is_sort1=True, is_random=False, is_msc=True,
                           random_code=0, title='', subtitle='', label=''):
        data_code = cls._add_case(
            table_name, element_name,
            isubcase, is_sort1, is_random, is_msc,
            random_code, title, subtitle, label)
        obj = set_transient_case(cls, is_sort1, isubcase, data_code,
                                 set_element_case, (element, data),
                                 times)
        return obj

    @classmethod
    def add_post_buckling_case(cls, table_name, element_name, element, data, isubcase,
                               modes, eigrs, eigis,
                               is_sort1=True, is_random=False, is_msc=True,
                               random_code=0, title='', subtitle='', label=''):
        data_code = cls._add_case(
            table_name, element_name,
            isubcase, is_sort1, is_random, is_msc,
            random_code, title, subtitle, label)
        obj = set_post_buckling_case(cls, is_sort1, isubcase, data_code,
                                     set_element_case, (element, data),
                                     modes, eigrs, eigis)
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

    def _reset_indices(self) -> None:
        self.itotal = 0
        self.ielement = 0

    def _get_msgs(self):
        raise NotImplementedError()

    def get_headers(self):
        raise NotImplementedError()

    def build(self):
        """sizes the vectorized attributes of the RealRodArray"""
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
        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size, self.analysis_fmt)
        self.build_data(self.ntimes, self.nelements, dtype, idtype, fdtype)

    def build_data(self, ntimes, nelements, dtype, idtype, fdtype):
        """actually performs the build step"""
        self.ntimes = ntimes
        self.nelements = nelements
        _times = zeros(ntimes, dtype=self.analysis_fmt)
        element = zeros(nelements, dtype=idtype)

        #[axial, torsion, SMa, SMt]
        data = zeros((ntimes, nelements, 4), dtype=fdtype)

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
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()
        if self.nonlinear_factor not in (None, np.nan):
            column_names, column_values = self._build_dataframe_transient_header()
            data_frame = self._build_pandas_transient_elements(
                column_values, column_names,
                headers, self.element, self.data)
        else:
            #Static     axial           SMa  torsion           SMt
            #ElementID
            #14           0.0  1.401298e-45      0.0  1.401298e-45
            #15           0.0  1.401298e-45      0.0  1.401298e-45
            data_frame = pd.DataFrame(self.data[0], columns=headers, index=self.element)
            data_frame.index.name = 'ElementID'
            data_frame.columns.names = ['Static']
            #Static                        0
            #ElementID Item
            #14        axial    0.000000e+00
            #          SMa      1.401298e-45
            #          torsion  0.000000e+00
            #          SMt      1.401298e-45
            #15        axial    0.000000e+00
            #          SMa      1.401298e-45
            #          torsion  0.000000e+00
            #          SMt      1.401298e-45
            #data_frame = pd.Panel(self.data, major_axis=self.element, minor_axis=headers).to_frame()
            #data_frame.columns.names = ['Static']
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
                        (axial1, torsion1, sma1, smt1) = t1
                        (axial2, torsion2, sma2, smt2) = t2
                        if not allclose(t1, t2):
                        #if not np.array_equal(t1, t2):
                            msg += '%s\n  (%s, %s, %s, %s)\n  (%s, %s, %s, %s)\n' % (
                                eid,
                                axial1, torsion1, sma1, smt1,
                                axial2, torsion2, sma2, smt2)
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

    def add_sort1(self, dt, eid, axial, SMa, torsion, SMt):
        assert self.sort_method == 1, self
        self._times[self.itime] = dt
        #if self.itime == 0:
        #print('itime=%s eid=%s' % (self.itime, eid))
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [axial, SMa, torsion, SMt]
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
        msg.append(f'  element.shape = {self.element.shape}\n')
        msg.append(f'  data.shape = {self.data.shape}\n')
        msg.append(f'  element type: {self.element_name}-{self.element_type}\n')
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True):
        crod_msg, conrod_msg, ctube_msg = self._get_msgs()
        if 'CROD' in self.element_name:
            msg = crod_msg
        elif 'CONROD' in self.element_name:
            msg = conrod_msg
        elif 'CTUBE' in self.element_name:
            msg = ctube_msg
        else:
            raise NotImplementedError(self.element_name)
        return self.element_name, msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element)  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.element_node[:, 0] == eid) for eid in eids])
        ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_csv(self, csv_file: TextIO,
                  is_exponent_format: bool=False,
                  is_mag_phase: bool=False, is_sort1: bool=True,
                  write_header: bool=True):
        """
        Stress Table - CROD/CONROD/CTUBE
        --------------------------------
        Flag, SubcaseID, iTime, EID,  BLANK,  BLANK,  Sxx,       Syy,  Szz,      Sxy,  Syz,  Szx
        10,           1,     0, 312,      0,      0,  1642.503,    0,    0,  167.541,    0,    0
        10,           1,     0, 313,      0,      0,  3937.541,    0,    0,   66.171,    0,    0
        """
        name = str(self.__class__.__name__)
        if write_header:
            csv_file.write('# %s\n' % name)
            headers = ['Flag', 'SubcaseID', 'iTime', 'Eid', 'BLANK', 'BLANK', 'Sxx', 'Syy', 'Szz', 'Sxy', 'Syz', 'Sxz']
            csv_file.write('# ' + ','.join(headers) + '\n')

        # stress vs. strain
        flag = 10 if 'Stress' in name else 11

        isubcase = self.isubcase
        #times = self._times

        # write the f06
        ntimes = self.data.shape[0]
        eids = self.element
        eid_len = '%d' % len(str(eids.max()))

        zero = ' 0.000000E+00'
        for itime in range(ntimes):
            dt = self._times[itime]
            #header = _eigenvalue_header(self, header, itime, ntimes, dt)

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            axial = self.data[itime, :, 0]
            torsion = self.data[itime, :, 2]
            for eid, axiali, torsioni in zip(eids, axial, torsion):
                if is_exponent_format:
                    [axiali, torsioni] = write_floats_13e_long([axiali, torsioni])
                csv_file.write(f'{flag}, {isubcase}, {itime}, {eid:{eid_len}d}, 0, 0, '
                               f'{axiali}, {zero}, {zero}, {torsioni}, {zero}, {zero}\n')
        return

    def write_f06(self, f06_file: TextIO, header=None, page_stamp: str='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        if header is None:
            header = []
        (unused_elem_name, msg_temp) = self.get_f06_header(is_mag_phase)
        if self.is_sort1:
            page_num = self._write_sort1_as_sort1(header, page_stamp, page_num, f06_file, msg_temp)
        return page_num

    def _write_sort1_as_sort1(self, header, page_stamp, page_num, f06_file: TextIO, msg_temp):
        ntimes = self.data.shape[0]

        eids = self.element
        is_odd = False
        nwrite = len(eids)
        if len(eids) % 2 == 1:
            nwrite -= 1
            is_odd = True

        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            axial = self.data[itime, :, 0]
            SMa = self.data[itime, :, 1]
            torsion = self.data[itime, :, 2]
            SMt = self.data[itime, :, 3]

            out = []
            for eid, axiali, SMai, torsioni, SMti in zip(eids, axial, SMa, torsion, SMt):
                [axiali, torsioni, SMai, SMti] = write_floats_13e([axiali, torsioni, SMai, SMti])
                out.append([eid, axiali, SMai, torsioni, SMti])

            for i in range(0, nwrite, 2):
                f06_file.write(
                    '      %8i %-13s  %-13s %-13s  %-13s %-8i   %-13s  %-13s %-13s  %-s\n' % (
                        tuple(out[i] + out[i + 1])))
            if is_odd:
                f06_file.write('      %8i %-13s  %-13s %-13s  %13s\n' % (tuple(out[-1])))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def write_op2(self, op2_file, op2_ascii, itable, new_result, date,
                  is_mag_phase=False, endian='>'):
        """writes an OP2"""
        import inspect
        from struct import pack # Struct,
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write(f'{self.__class__.__name__}.write_op2: {call_frame[1][3]}\n')

        if itable == -1:
            self._write_table_header(op2_file, op2_ascii, date)
            itable = -3

        #if isinstance(self.nonlinear_factor, float):
            #op2_format = '%sif' % (7 * self.ntimes)
            #raise NotImplementedError()
        #else:
            #op2_format = 'i21f'
        #s = Struct(op2_format)

        unused_eids = self.element

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

        op2_ascii.write(f'  ntimes = {self.ntimes}\n')

        eids_device = self.element * 10 + self.device_code

        #fmt = '%2i %6f'
        #print('ntotal=%s' % (ntotal))
        #assert ntotal == 193, ntotal

        if not self.is_sort1:
            raise NotImplementedError('SORT2')
        #struct1 = Struct(endian + b'i4f')

        fdtype = self.data.dtype
        if self.size == fdtype.itemsize:
            pass
        else:
            # print(f'downcasting {self.class_name}...')
            #idtype = np.int32(1)
            fdtype = np.float32(1.0)

        # [eid, axial, SMa, torsion, SMt]
        data_out = np.empty((nelements, 5), dtype=fdtype)
        data_out[:, 0] = view_idtype_as_fdtype(eids_device, fdtype)

        op2_ascii.write(f'nelements={nelements:d}\n')

        for itime in range(self.ntimes):
            #print('3, %s' % itable)
            self._write_table_3(op2_file, op2_ascii, new_result, itable, itime)

            # record 4
            #print('stress itable = %s' % itable)
            itable -= 1
            #print('4, %s' % itable)
            header = [4, itable, 4,
                      4, 1, 4,
                      4, 0, 4,
                      4, ntotal, 4,
                      4 * ntotal]
            op2_file.write(pack('%ii' % len(header), *header))
            op2_ascii.write('r4 [4, 0, 4]\n')
            op2_ascii.write(f'r4 [4, {itable:d}, 4]\n')
            op2_ascii.write(f'r4 [4, {4 * ntotal:d}, 4]\n')

            # [eid, axial, SMa, torsion, SMt]
            data_out[:, 1:] = self.data[itime, :, :]
            op2_file.write(data_out)

            itable -= 1
            header = [4 * ntotal,]
            op2_file.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable


class RealBushStressArray(RealRodArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealRodArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> list[str]:
        headers = ['axial', 'SMa', 'torsion', 'SMt']
        return headers

    def _get_msgs(self):
        raise NotImplementedError()
        #base_msg = ['       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                    #'         ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN\n']
        #crod_msg   = ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )\n', ]
        #conrod_msg = ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C O N R O D )\n', ]
        #ctube_msg  = ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C T U B E )\n', ]
        ##cbush_msg  = ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C B U S H )\n', ]
        #crod_msg += base_msg
        #conrod_msg += base_msg
        #ctube_msg += base_msg
        ##cbush_msg += base_msg
        #return crod_msg, conrod_msg, ctube_msg


class RealRodStressArray(RealRodArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealRodArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> list[str]:
        headers = ['axial', 'SMa', 'torsion', 'SMt']
        return headers

    def _get_msgs(self):
        base_msg = ['       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                    '         ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN\n']
        crod_msg = ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )\n', ]
        conrod_msg = ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C O N R O D )\n', ]
        ctube_msg = ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C T U B E )\n', ]
        crod_msg += base_msg
        conrod_msg += base_msg
        ctube_msg += base_msg
        return crod_msg, conrod_msg, ctube_msg

class RealRodStrainArray(RealRodArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealRodArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> list[str]:
        headers = ['axial', 'SMa', 'torsion', 'SMt']
        return headers

    def _get_msgs(self):
        base_msg = ['       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                    '         ID.        STRAIN       MARGIN        STRAIN      MARGIN         ID.        STRAIN       MARGIN        STRAIN      MARGIN\n']
        crod_msg = ['                                       S T R A I N S   I N   R O D   E L E M E N T S      ( C R O D )\n', ]
        conrod_msg = ['                                       S T R A I N S   I N   R O D   E L E M E N T S      ( C O N R O D )\n', ]
        ctube_msg = ['                                       S T R A I N S   I N   R O D   E L E M E N T S      ( C T U B E )\n', ]
        crod_msg += base_msg
        conrod_msg += base_msg
        ctube_msg += base_msg
        return crod_msg, conrod_msg, ctube_msg
