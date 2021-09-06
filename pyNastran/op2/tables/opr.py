from __future__ import annotations
from collections import defaultdict
import numpy as np
from typing import TYPE_CHECKING

import numpy as np
from pyNastran.op2.result_objects.table_object import RealTableArray # , ComplexTableArray
from pyNastran.utils.numpy_utils import integer_types

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


class RealPressureArray(RealTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        if header is None:
            header = []
        words = ['                                                   P R E S S U R E\n', ]
        #' \n',
        #'      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #if self.table_name in ['OUGV1', 'OUGV2', 'BOUGV1', 'OUPV1', 'OUGV1PAT', 'OUG1', 'OUG2']:
            #pass
        #elif  self.table_name in ['OUXY1', 'OUXY2']:
            #words = ['                                       D I S P L A C E M E N T   V E C T O R   (SOLUTION SET)']
        #elif self.table_name in ['ROUGV1', 'ROUGV2']:
            #words += ['                                                (RELATIVE TO ENFORCED MOTION INPUT)']
        if self.table_name in ['OPRATO1', 'OPRATO2']:
            words += ['                                                 ( AUTO-CORRELATION FUNCTION )']
        elif self.table_name in ['OPRPSD1', 'OPRPSD2']:
            words += ['                                             ( POWER SPECTRAL DENSITY FUNCTION )']
        elif self.table_name in ['OPRRMS1', 'OPRRMS2']:
            words += ['                                                     ( ROOT MEAN SQUARE )']
        elif self.table_name in ['OPRCRM1', 'OPRCRM2']:
            words += ['                                               ( CUMULATIVE ROOT MEAN SQUARE )']
        elif self.table_name in ['OPRNO1', 'OPRNO2']:
            words += ['                                                 ( NUMBER OF ZERO CROSSINGS )']
        #elif self.table_name in ['OCRUG']:
            #words += ['                                                 ( OCRUG??? )']
        else:
            raise NotImplementedError('table_name=%r' % self.table_name)
        #words += self.get_table_marker()
        write_words = True
        if self.nonlinear_factor not in (None, np.nan):
            return self._write_f06_transient_block(words, header, page_stamp, page_num, f06_file, write_words,
                                                   is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        return self._write_f06_block(words, header, page_stamp, page_num, f06_file, write_words,
                                     is_mag_phase=is_mag_phase, is_sort1=is_sort1)


class OPR:
    def __init__(self, op2: OP2):
        self.op2 = op2
        self.obj = None

    def _read_opr1_3(self, data: bytes, ndata: int):
        """
        reads table 3 (the header table)

        Word Name Type Description
        1 ACODE(C)     I     Device code + 10*Approach code
        2 TCODE(C)     I     Table code 92
        3 UNDEF(3)     None
        6 DCYCLE       I     Design cycle number
        7 ROBJ         RS    Objective value
        8 RCON         RS    Critical constraint value
        9 UNDEF        None
        10 NUMWDE      I     Number of words per entry in DATA record (always 2)
        11 UNDEF(40)   None
        51 TITLE(32)   CHAR4 Title character string (TITLE)
        83 SUBTITL(32) CHAR4 Subtitle character string (SUBTITLE)
        115 LABEL(32)  CHAR4 LABEL character string (LABEL)
        """
        op2 = self.op2
        op2.to_nx('; found OPR (pressure) table')
        """reads table 3 (the header table)"""
        op2 = self.op2
        #self._set_times_dtype()
        op2.nonlinear_factor = np.nan
        op2.is_table_1 = True
        op2.is_table_2 = False
        unused_three = op2.parse_approach_code(data)
        op2.words = [
            'approach_code', 'table_code', '???', 'isubcase',
            '???', '???', '???', 'random_code',
            'format_code', 'num_wide', '???', '???',
            'acoustic_flag', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', 'thermal', '???',
            '???', 'Title', 'subtitle', 'label']

        ## random code
        op2.random_code = op2.add_data_parameter(data, 'random_code', b'i', 8, False)

        ## format code
        op2.format_code = op2.add_data_parameter(data, 'format_code', b'i', 9, False)

        ## number of words per entry in record
        op2.num_wide = op2.add_data_parameter(data, 'num_wide', b'i', 10, False)

        ## acoustic pressure flag
        op2.acoustic_flag = op2.add_data_parameter(data, 'acoustic_flag', b'i', 13, False)

        ## thermal flag; 1 for heat transfer, 0 otherwise
        op2.thermal = op2.add_data_parameter(data, 'thermal', b'i', 23, False)

        if op2.analysis_code == 1:   # statics / displacement / heat flux
            # load set number
            op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
            op2.setNullNonlinearFactor()
        elif op2.analysis_code == 2:  # real eigenvalues
            # mode number
            op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            # eigenvalue
            op2.eign = op2.add_data_parameter(data, 'eign', b'f', 6, False)
            # mode or cycle .. todo:: confused on the type - F1???
            # float - C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\mftank.op2
            #op2.mode_cycle = op2.add_data_parameter(data, 'mode_cycle', b'i', 7, False)  # nope...
            op2.mode_cycle = op2.add_data_parameter(data, 'mode_cycle', b'f', 7, False) # radians
            self.update_mode_cycle('mode_cycle')
            op2.data_names = op2.apply_data_code_value('data_names', ['mode', 'eign', 'mode_cycle'])
        #elif op2.analysis_code == 3: # differential stiffness
            #op2.lsdvmn = self.get_values(data, b'i', 5) ## load set number
            #op2.data_code['lsdvmn'] = op2.lsdvmn
        #elif op2.analysis_code == 4: # differential stiffness
            #op2.lsdvmn = self.get_values(data, b'i', 5) ## load set number
        elif op2.analysis_code == 5:   # frequency
            # frequency
            op2.freq = op2.add_data_parameter(data, 'freq', b'f', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['freq'])
        elif op2.analysis_code == 6:  # transient
            # time step
            op2.dt = op2.add_data_parameter(data, 'dt', b'f', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['dt'])
        elif op2.analysis_code == 7:  # pre-buckling
            # load set number
            op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
        elif op2.analysis_code == 8:  # post-buckling
            # load set number
            op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            # real eigenvalue
            op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn', 'eigr'])
        elif op2.analysis_code == 9:  # complex eigenvalues
            # mode number
            op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            # real eigenvalue
            op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            # imaginary eigenvalue
            op2.eigi = op2.add_data_parameter(data, 'eigi', b'f', 7, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['mode', 'eigr', 'eigi'])
        elif op2.analysis_code == 10:  # nonlinear statics
            # load step
            op2.lftsfq = op2.add_data_parameter(data, 'lftsfq', b'f', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['lftsfq'])
        elif op2.analysis_code == 11:  # old geometric nonlinear statics
            # load set number
            op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
        elif op2.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            # load set number
            op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
        else:  # pragma: no cover
            msg = f'invalid analysis_code...analysis_code={op2.analysis_code}\ndata={op2.data_code}'
            raise RuntimeError(msg)

        #print op2.code_information()
        op2._fix_oug_format_code()
        op2._parse_thermal_code()
        if op2.is_debug_file:
            op2.binary_debug.write('  approach_code  = %r\n' % op2.approach_code)
            op2.binary_debug.write('  tCode          = %r\n' % op2.tCode)
            op2.binary_debug.write('  isubcase       = %r\n' % op2.isubcase)
        op2._read_title(data)
        op2._write_debug_bits()
        #self._correct_eigenvalue()

    def _read_opr2_3(self, data: bytes, ndata: int) -> None:
        """reads the SORT2 version of table 4 (the data table)"""
        op2 = self.op2
        op2.to_nx('; found OPR (pressure) table')
        #self._set_times_dtype()
        #return self._read_oug1_3(data)
        op2.nonlinear_factor = np.nan

        op2.is_table_1 = False
        op2.is_table_2 = True
        unused_three = op2.parse_approach_code(data)
        op2.words = [
            'approach_code', 'table_code', '???', 'isubcase',
            '???', '???', '???', 'random_code',
            'format_code', 'num_wide', '???', '???',
            'acoustic_flag', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', 'thermal', '???',
            '???', 'Title', 'subtitle', 'label']

        ## random code
        op2.random_code = op2.add_data_parameter(data, 'random_code', b'i', 8, False)

        ## format code
        op2.format_code = op2.add_data_parameter(data, 'format_code', b'i', 9, False)

        ## number of words per entry in record
        op2.num_wide = op2.add_data_parameter(data, 'num_wide', b'i', 10, False)

        ## acoustic pressure flag
        op2.acoustic_flag = op2.add_data_parameter(data, 'acoustic_flag', b'i', 13, False)

        ## thermal flag; 1 for heat transfer, 0 otherwise
        op2.thermal = op2.add_data_parameter(data, 'thermal', b'i', 23, False)

        op2.node_id = op2.add_data_parameter(data, 'node_id', b'i', 5, fix_device_code=True)
        #if op2.analysis_code == 1:  # statics / displacement / heat flux
            ## load set number
            #op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5, False)
            #op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            #op2.setNullNonlinearFactor()

        if op2.analysis_code == 1:  # static...because reasons.
            op2._analysis_code_fmt = b'f'
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            op2.apply_data_code_value('analysis_method', 'N/A')
        elif op2.analysis_code == 2:  # real eigenvalues
            ## mode number
            #op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            op2._analysis_code_fmt = b'i'
            # real eigenvalue
            op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            ## mode or cycle .. todo:: confused on the type - F1???
            # float - C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\mftank.op2
             # mode or cycle .. todo:: confused on the type - F1???
            op2.mode_cycle = op2.add_data_parameter(data, 'mode_cycle', b'i', 7, False)
            #op2.mode_cycle = op2.add_data_parameter(data, 'mode_cycle', b'f', 7, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id', 'eigr', 'mode_cycle'])
            op2.apply_data_code_value('analysis_method', 'mode')
        #elif op2.analysis_code == 3: # differential stiffness
            #op2.lsdvmn = self.get_values(data, b'i', 5) ## load set number
            #op2.data_code['lsdvmn'] = op2.lsdvmn
        #elif op2.analysis_code == 4: # differential stiffness
            #op2.lsdvmn = self.get_values(data, b'i', 5) ## load set number
        elif op2.analysis_code == 5:   # frequency
            # frequency
            #op2.freq = op2.add_data_parameter(data, 'freq', b'f', 5)
            op2._analysis_code_fmt = b'f'
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            op2.apply_data_code_value('analysis_method', 'freq')
        elif op2.analysis_code == 6:  # transient
            ## time step
            #op2.dt = op2.add_data_parameter(data, 'dt', b'f', 5)
            op2._analysis_code_fmt = b'f'
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            op2.apply_data_code_value('analysis_method', 'dt')
        elif op2.analysis_code == 7:  # pre-buckling
            ## load set number
            #op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            op2._analysis_code_fmt = b'i'
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            op2.apply_data_code_value('analysis_method', 'lsdvmn')
        elif op2.analysis_code == 8:  # post-buckling
            ## load set number
            #op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            op2._analysis_code_fmt = b'f'
            ## real eigenvalue
            op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id', 'eigr'])
            op2.apply_data_code_value('analysis_method', 'eigr')
        elif op2.analysis_code == 9:  # complex eigenvalues
            ## mode number
            #op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            op2._analysis_code_fmt = b'i'
            ## real eigenvalue
            op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            ## imaginary eigenvalue
            op2.eigi = op2.add_data_parameter(data, 'eigi', b'f', 7, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id', 'eigr', 'eigi'])
            op2.apply_data_code_value('analysis_method', 'mode')
        elif op2.analysis_code == 10:  # nonlinear statics
            ## load step
            #op2.lftsfq = op2.add_data_parameter(data, 'lftsfq', b'f', 5)
            op2._analysis_code_fmt = b'f'
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            op2.apply_data_code_value('analysis_method', 'lftsfq')
        elif op2.analysis_code == 11:  # old geometric nonlinear statics
            ## load set number
            #op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
        elif op2.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## load set number
            #op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            op2.apply_data_code_value('analysis_method', 'lsdvmn')
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % op2.analysis_code
            raise RuntimeError(msg)

        op2._fix_oug_format_code()
        op2._parse_thermal_code()
        if op2.is_debug_file:
            op2.binary_debug.write('  %-14s = %r %s\n' % ('approach_code', op2.approach_code,
                                                           op2.approach_code_str(op2.approach_code)))
            op2.binary_debug.write('  %-14s = %r\n' % ('tCode', op2.tCode))
            op2.binary_debug.write('  %-14s = %r\n' % ('isubcase', op2.isubcase))
        op2._read_title(data)
        op2._write_debug_bits()
        assert isinstance(op2.nonlinear_factor, integer_types), op2.nonlinear_factor

    def _read_opr_ato(self, data: bytes, ndata: int) -> int:
        n = self._read_opr(data, ndata, 'ato.')
        return n

    def _read_opr_crm(self, data: bytes, ndata: int) -> int:
        n = self._read_opr(data, ndata, 'crm.')
        return n

    def _read_opr_psd(self, data: bytes, ndata: int) -> int:
        n = self._read_opr(data, ndata, 'psd.')
        return n

    def _read_opr_no(self, data: bytes, ndata: int) -> int:
        n = self._read_opr(data, ndata, 'no.')
        return n

    def _read_opr(self, data: bytes, ndata: int, prefix: str) -> int:
        """
      POINT ID.   TYPE           P              P              P              P              P              P
         50000      S      5.821456E-03
         50003      S      1.288847E-02
         50024      S      2.896047E-03
         50040      S      4.176410E-03
         50042      S      3.067493E-03
         50075      S      4.594804E-03
         50088      S      5.080107E-03

        ints    = (500001, 2, 1102221187, 0, 0, 0, 0, 0,
                   500031, 2, 1102162789, 0, 0, 0, 0, 0,
                   500241, 2, 1102331348, 0, 0, 0, 0, 0,
                   500401, 2, 1102175796, 0, 0, 0, 0, 0, 500421, 2, 1102263146, 0, 0, 0, 0, 0, 500751, 2, 1102573063, 0, 0, 0, 0, 0, 500881, 2, 1102274036, 0, 0, 0, 0, 0)
        floats  = (500001, 2, 22.3200740814209, 0.0, 0.0, 0.0, 0.0, 0.0, 7.006926724148026e-40, 2.802596928649634e-45, 22.208688735961914, 0.0, 0.0, 0.0, 0.0, 0.0, 7.009869450923108e-40, 2.802596928649634e-45, 22.530189514160156, 0.0, 0.0, 0.0, 0.0, 0.0, 7.012111528466028e-40, 2.802596928649634e-45, 22.233497619628906, 0.0, 0.0, 0.0, 0.0, 0.0, 7.012391788158893e-40, 2.802596928649634e-45, 22.400104522705078, 0.0, 0.0, 0.0, 0.0, 0.0, 7.017016073091165e-40, 2.802596928649634e-45, 22.99122428894043, 0.0, 0.0, 0.0, 0.0, 0.0, 7.018837761094787e-40, 2.802596928649634e-45, 22.420875549316406, 0.0, 0.0, 0.0, 0.0, 0.0)
        """
        op2 = self.op2
        if op2.thermal == 0:
            if op2.table_code == 501:
                # displacement
                assert op2.table_name in [b'OPRCRM1', b'OPRCRM2'], 'op2.table_name=%r' % op2.table_name
                result_name = prefix + 'pressures'
                obj = RealPressureArray
            elif op2.table_code == 601:
                # displacement
                assert op2.table_name in [b'OPRPSD1', b'OPRPSD2'], 'op2.table_name=%r' % op2.table_name
                result_name = prefix + 'pressures'
                obj = RealPressureArray
            elif op2.table_code == 701:
                # displacement
                assert op2.table_name in [b'OPRATO1', b'OPRATO2'], 'op2.table_name=%r' % op2.table_name
                result_name = prefix + 'pressures'
                obj = RealPressureArray
            elif op2.table_code == 901:
                # displacement
                assert op2.table_name in [b'OPRNO1', b'OPRNO2'], 'op2.table_name=%r' % op2.table_name
                result_name = prefix + 'pressures'
                obj = RealPressureArray
            else:
                raise RuntimeError(op2.code_information())
        else:
            raise RuntimeError(op2.code_information())

        #if data is None:
            #assert isinstance(ndata, int), ndata
            #return ndata

        op2 = self.op2
        #idata = np.frombuffer(data, dtype=op2.idtype8)
        #ndata = len(idata)
        #nrows = ndata // 8
        #idata = idata.reshape(nrows, 8)
        #fdata = np.frombuffer(data, dtype=op2.fdtype8).reshape(nrows, 8)[:, 2:]
        #[50000 50003 50024 50040 50042 50075 50088]

        if op2._results.is_not_saved(result_name):
            return ndata
        op2._results._found_result(result_name)
        storage_obj = op2.get_result(result_name)

        n = op2._read_table_vectorized(data, ndata, result_name, storage_obj,
                                       RealPressureArray, None,
                                       'node', random_code=op2.random_code)
        #n = op2._read_scalar_table_vectorized(data, ndata, result_name, storage_obj,
                                              #RealPressureArray, None,
                                              #'node', random_code=op2.random_code,
                                              #is_cid=False)

        #eid = idata[:, 0] // 10
        #gridtype = idata[:, 1]
        #print(eid)
        #print(idata)
        #print(fdata)

        #op2.show_data(data)
        return ndata
        #asdf
    def _read_oprs_4(self, data: bytes, ndata: int) -> int:
        """
        Word Name Type Description
        1 EKEY  I Device code + 10*Element ID
        2 VALUE RS Scalar value for element
        """
        if not data:
            return ndata
        op2 = self.op2

        fdata = np.frombuffer(data, dtype=op2.fdtype8)
        idata = np.frombuffer(data, dtype=op2.idtype8)
        ndata = len(idata)
        #op2.log.warning(f'ndata={ndata}')
        eids = idata[::2] // 10
        data = fdata[1::2]
        #print(f'eids = {eids}')
        #print(f'data = {data}')
        self.obj.eids = eids
        self.obj.data = data

        key = op2.isubcase
        responses = op2.op2_results.responses
        if responses.normalized_mass_density is None:
            responses.normalized_mass_density = defaultdict(list)
        responses.normalized_mass_density[key].append(self.obj)
        return ndata
