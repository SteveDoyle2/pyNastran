"""
This file defines the OUG Table, which contains:
 * Real/Complex SPC Forces
   - SPCFORCE = ALL
 * Real/Complex MPC Forces
   - MPCFORCE = ALL
 * Real Temperature Gradient & Flux
   - FLUX = ALL
"""
from __future__ import print_function
from six import integer_types
from pyNastran.op2.op2_interface.op2_common import OP2Common

from pyNastran.op2.tables.oqg_constraintForces.oqg_spc_forces import (
    RealSPCForcesArray, ComplexSPCForcesArray,
    #RealSPCForces, ComplexSPCForces
)
from pyNastran.op2.tables.oqg_constraintForces.oqg_mpc_forces import (
    RealMPCForcesArray, ComplexMPCForcesArray,
    #RealMPCForces, ComplexMPCForces
)
from pyNastran.op2.tables.oqg_constraintForces.oqg_thermal_gradient_and_flux import (
    #RealTemperatureGradientAndFlux,
    RealTemperatureGradientAndFluxArray)

from pyNastran.op2.tables.table_reader import TableReader


class OQG(TableReader):
    """
    Defines the OQS class that is used to read SPC/MPC forces
    """
    def __init__(self, op2_reader, op2):
        TableReader.__init__(self, op2_reader, op2)

    def _read_oqg1_3(self, data, ndata):
        op2_reader = self.op2_reader
        op2_reader._data_factor = 1
        op2_reader.nonlinear_factor = None
        op2_reader.is_table_1 = True
        op2_reader.is_table_2 = False
        unused_three = self.parse_approach_code(data)
        op2_reader.words = [
            'analysis_code', 'table_code', '???', 'isubcase',
            '???', '???', '???', 'random_code',
            'format_code', 'num_wide', '11', '???',
            'acoustic_flag', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', 'thermal', '???',
            '???', 'Title', 'subtitle', 'label']

        ## random code
        op2_reader.random_code = self.add_data_parameter(data, 'random_code', b'i', 8, False)

        ## format code
        op2_reader.format_code = self.add_data_parameter(data, 'format_code', b'i', 9, False)

        ## number of words per entry in record
        op2_reader.num_wide = self.add_data_parameter(data, 'num_wide', b'i', 10, False)

        ## acoustic pressure flag
        op2_reader.acoustic_flag = self.add_data_parameter(data, 'acoustic_flag', b'f', 13, False)

        ## thermal flag; 1 for heat transfer, 0 otherwise
        op2_reader.thermal = self.add_data_parameter(data, 'thermal', b'i', 23, False)

        if not self.is_sort1:
            raise NotImplementedError('SORT2; code_info=\n%s' % self.code_information())
        #assert self.isThermal()==False,self.thermal

        ## assuming tCode=1
        analysis_code = op2_reader.analysis_code
        if analysis_code == 1:   # statics / displacement / heat flux
            ## load set number
            op2_reader.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5, False)
            op2_reader.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
            op2_reader.set_null_nonlinear_factor()
        elif analysis_code == 2:  # real eigenvalues
            ## mode number
            op2_reader.mode = self.add_data_parameter(data, 'mode', b'i', 5)
            ## eigenvalue
            op2_reader.eign = self.add_data_parameter(data, 'eign', b'f', 6, False)
            ## mode or cycle .. todo:: confused on the type - F1???
            op2_reader.mode_cycle = self.add_data_parameter(data, 'mode_cycle', b'f', 7, False)
            self.update_mode_cycle('mode_cycle')
            op2_reader.data_names = self.apply_data_code_value('data_names',
                                                               ['mode', 'eign', 'mode_cycle'])
        #elif self.analysis_code == 3: # differential stiffness
            #self.lsdvmn = self.get_values(data, b'i', 5) ## load set number
            #self.data_names = self.data_code['lsdvmn'] = self.lsdvmn
        #elif self.analysis_code == 4: # differential stiffness
            #self.lsdvmn = self.get_values(data, b'i', 5) ## load set number
        elif analysis_code == 5:   # frequency
            ## frequency
            self.freq = self.add_data_parameter(data, 'freq', b'f', 5)
            self.data_names = self.apply_data_code_value('data_names', ['freq'])
        elif analysis_code == 6:  # transient
            ## time step
            self.dt = self.add_data_parameter(data, 'dt', b'f', 5)
            self.data_names = self.apply_data_code_value('data_names', ['dt'])
        elif self.analysis_code == 7:  # pre-buckling
            ## load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5)
            op2_reader.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
        elif analysis_code == 8:  # post-buckling
            ## load set number
            op2_reader.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5)
            ## real eigenvalue
            op2_reader.eigr = self.add_data_parameter(data, 'eigr', b'f', 6, False)
            op2_reader.data_names = self.apply_data_code_value('data_names', ['lsdvmn', 'eigr'])
        elif analysis_code == 9:  # complex eigenvalues
            ## mode number
            op2_reader.mode = self.add_data_parameter(data, 'mode', b'i', 5)
            ## real eigenvalue
            op2_reader.eigr = self.add_data_parameter(data, 'eigr', b'f', 6, False)
            ## imaginary eigenvalue
            op2_reader.eigi = self.add_data_parameter(data, 'eigi', b'f', 7, False)
            op2_reader.data_names = self.apply_data_code_value('data_names', ['mode', 'eigr', 'eigi'])
        elif analysis_code == 10:  # nonlinear statics
            ## load step
            op2_reader.lftsfq = self.add_data_parameter(data, 'lftsfq', b'f', 5)
            op2_reader.data_names = self.apply_data_code_value('data_names', ['lftsfq'])
        elif analysis_code == 11:  # old geometric nonlinear statics
            ## load set number
            op2_reader.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5)
            op2_reader.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
        elif analysis_code == 12:
            # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## load set number
            op2_reader.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5)
            op2_reader.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % analysis_code
            raise RuntimeError(msg)

        self.fix_format_code()
        if op2_reader.num_wide == 8:
            op2_reader.format_code = 1
            op2_reader.data_code['format_code'] = 1
        else:
            #self.fix_format_code()
            if op2_reader.format_code == 1:
                op2_reader.format_code = 2
                op2_reader.data_code['format_code'] = 2
            assert op2_reader.format_code in [2, 3], op2_reader.code_information()

        self._parse_thermal_code()
        if self.is_debug_file:
            self.binary_debug.write('  approach_code  = %r\n' % self.approach_code)
            self.binary_debug.write('  tCode          = %r\n' % self.tCode)
            self.binary_debug.write('  isubcase       = %r\n' % self.isubcase)
        op2_reader._read_title(data)
        op2_reader._write_debug_bits()

    def _read_oqg2_3(self, data, ndata):
        """reads the SORT2 version of table 4 (the data table)"""
        op2_reader = self.op2_reader
        op2_reader._data_factor = 1
        op2_reader.nonlinear_factor = None
        op2_reader.is_table_1 = False
        op2_reader.is_table_2 = True
        unused_three = self.parse_approach_code(data)
        op2_reader.words = [
            'analysis_code', 'table_code', '???', 'isubcase',
            '???', '???', '???', 'random_code',
            'format_code', 'num_wide', '11', '???',
            'acoustic_flag', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', 'thermal', '???',
            '???', 'Title', 'subtitle', 'label']

        ## random code
        op2_reader.random_code = self.add_data_parameter(data, 'random_code', b'i', 8, False)

        ## format code
        op2_reader.format_code = self.add_data_parameter(data, 'format_code', b'i', 9, False)

        ## number of words per entry in record
        op2_reader.num_wide = self.add_data_parameter(data, 'num_wide', b'i', 10, False)

        ## acoustic pressure flag
        op2_reader.acoustic_flag = self.add_data_parameter(data, 'acoustic_flag', b'f', 13, False)

        ## thermal flag; 1 for heat transfer, 0 otherwise
        op2_reader.thermal = self.add_data_parameter(data, 'thermal', b'i', 23, False)

        #assert self.isThermal() == False, self.thermal

        self.node_id = self.add_data_parameter(data, 'node_id', b'i', 5, fix_device_code=True)
        #if self.analysis_code == 1:  # statics / displacement / heat flux
            ## load set number
            #self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5, False)
            #self.data_names = self.apply_data_code_value('data_names', ['node_id'])
            #op2_reader.set_null_nonlinear_factor()

        analysis_code = op2_reader.analysis_code
        if analysis_code == 1:  # static...because reasons.
            self._analysis_code_fmt = b'i'
            op2_reader.data_names = self.apply_data_code_value('data_names', ['node_id'])
            self.apply_data_code_value('analysis_method', 'N/A')
        elif analysis_code == 2:  # real eigenvalues
            ## mode number
            op2_reader.mode = self.add_data_parameter(data, 'mode', b'i', 5)
            self._analysis_code_fmt = b'i'
            ## real eigenvalue
            op2_reader.eigr = self.add_data_parameter(data, 'eigr', b'f', 6, False)
            ## mode or cycle .. todo:: confused on the type - F1???
            op2_reader.mode_cycle = self.add_data_parameter(data, 'mode_cycle', b'f', 7, False)
            op2_reader.data_names = self.apply_data_code_value('data_names',
                                                         ['node_id', 'eigr', 'mode_cycle'])
            self.apply_data_code_value('analysis_method', 'mode')
        #elif self.analysis_code == 3: # differential stiffness
            #op2_reader.lsdvmn = self.get_values(data, b'i', 5) ## load set number
            #op2_reader.data_names = self.data_code['lsdvmn'] = self.lsdvmn
        #elif self.analysis_code == 4: # differential stiffness
            #self.lsdvmn = self.get_values(data, b'i', 5) ## load set number
        elif self.analysis_code == 5:   # frequency
            ## frequency
            #op2_reader.freq = self.add_data_parameter(data, 'freq', b'f', 5)
            self._analysis_code_fmt = b'f'
            op2_reader.data_names = self.apply_data_code_value('data_names', ['node_id'])
            self.apply_data_code_value('analysis_method', 'freq')
        elif self.analysis_code == 6:  # transient
            ## time step
            #op2_reader.dt = self.add_data_parameter(data, 'dt', b'f', 5)
            self._analysis_code_fmt = b'f'
            op2_reader.data_names = self.apply_data_code_value('data_names', ['node_id'])
            self.apply_data_code_value('analysis_method', 'dt')
        elif self.analysis_code == 7:  # pre-buckling
            ## load set number
            #self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5)
            self._analysis_code_fmt = b'i'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
            self.apply_data_code_value('analysis_method', 'lsdvmn')
        elif self.analysis_code == 8:  # post-buckling
            ## load set number
            #self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5)
            self._analysis_code_fmt = b'i'
            ## real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', b'f', 6, False)
            self.data_names = self.apply_data_code_value('data_names', ['node_id', 'eigr'])
            self.apply_data_code_value('analysis_method', 'eigr')
        elif self.analysis_code == 9:  # complex eigenvalues
            ## mode number
            self.mode = self.add_data_parameter(data, 'mode', b'i', 5)
            self._analysis_code_fmt = b'i'
            ## real eigenvalue
            #self.eigr = self.add_data_parameter(data, 'eigr', b'f', 6, False)
            ## imaginary eigenvalue
            self.eigi = self.add_data_parameter(data, 'eigi', b'f', 7, False)
            self.data_names = self.apply_data_code_value('data_names', ['node_id', 'eigr', 'eigi'])
            self.apply_data_code_value('analysis_method', 'mode')
        elif self.analysis_code == 10:  # nonlinear statics
            ## load step
            #self.lftsfq = self.add_data_parameter(data, 'lftsfq', b'f', 5)
            self._analysis_code_fmt = b'f'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
            self.apply_data_code_value('analysis_method', 'lftsfq')
        elif self.analysis_code == 11:  # old geometric nonlinear statics
            ## load set number
            #self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5)
            self._analysis_code_fmt = b'f'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
        elif self.analysis_code == 12:
            # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## load set number
            #self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5)
            self._analysis_code_fmt = b'i'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
            self.apply_data_code_value('analysis_method', 'lsdvmn')
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % self.analysis_code
            raise RuntimeError(msg)

        self.fix_format_code()
        if self.num_wide == 8:
            self.format_code = 1
            self.data_code['format_code'] = 1
        else:
            #self.fix_format_code()
            if self.format_code == 1:
                self.format_code = 2
                self.data_code['format_code'] = 2
            assert self.format_code in [2, 3], self.code_information()

        self._parse_thermal_code()
        if self.is_debug_file:
            self.binary_debug.write('  approach_code  = %r\n' % self.approach_code)
            self.binary_debug.write('  tCode          = %r\n' % self.tCode)
            self.binary_debug.write('  isubcase       = %r\n' % self.isubcase)
        self._read_title(data)
        self._write_debug_bits()

    def _read_oqg_4(self, data, ndata):
        """
        This function is called by:

        MSC : SPC/MPC forces
         - SPC forces table_code = 3  (OQGPSD1 ???)
         - MPC forces table_code = 39 (OQGPSD1 ???)

        NX  : SPC forces
         - SPC forces table_code = 3  (OQGPSD1 ???)
         - MPC forces table_code = 3  (OQGMPSD1 ???)

        elif table_code == 3:
            table = "OQG - SPC Force vector"
        elif table_code == 39:
            table = "OQG - MPC Forces"

        """
        result_name = 'constraint_forces'
        if self._results.is_not_saved(result_name):
            return ndata

        op2_reader = self.op2_reader
        table_name = op2_reader.table_name
        table_code = op2_reader.table_code
        if table_name in [b'OQGCRM1', b'OQGCRM2']:
            if table_code not in [3, 503]:
                msg = 'table_name=%s table_code=%s' % (table_name, table_code)
                raise AssertionError(msg)
            n = self._read_oqg_spc_crm(data, ndata)
        elif table_name in [b'OQGPSD1', b'OQGPSD2']:
            if table_code not in [3, 603]:  # was 3
                msg = 'table_name=%s table_code=%s' % (table_name, table_code)
                raise AssertionError(msg + '\n%s' % self.code_information())
            n = self._read_oqg_spc_psd(data, ndata)
        elif table_name in [b'OQGATO1', b'OQGATO2']:
            if table_code not in [3, 703]:
                msg = 'table_name=%s table_code=%s' % (table_name, table_code)
                raise AssertionError(msg)
            n = self._read_oqg_spc_ato(data, ndata)
        elif table_name in [b'OQGRMS1', b'OQGRMS2']:
            if table_code not in [3, 803]:
                msg = 'table_name=%s table_code=%s' % (table_name, table_code)
                raise AssertionError(msg)
            n = self._read_oqg_spc_rms(data, ndata)
        elif table_name in [b'OQGNO1', b'OQGNO2']:
            if table_code not in [3, 903]:
                msg = 'table_name=%s table_code=%s' % (table_name, table_code)
                raise AssertionError(msg)
            n = self._read_oqg_spc_no(data, ndata)

        elif table_code == 3:   # SPC Forces
            assert self.table_name in [b'OQG1', b'OQGV1', b'OQP1', b'OQG2'], op2_reader.code_information()
            n = self._read_spc_forces(data, ndata)
        elif table_code == 39:  # MPC Forces
            assert table_name in [b'OQMG1', b'OQMG2'], op2_reader.code_information() # , b'OQMPSD1', b'OQMPSD2'
            n = self._read_oqg_mpc_forces(data, ndata)
        elif table_name in [b'RAQEATC', b'RAQCONS']:
            # self.table_code == 5 and
            n = self._read_oqg_mpc_forces(data, ndata)
        else:
            raise RuntimeError(op2_reader.code_information())
        return n

    def _read_spc_forces(self, data, ndata):
        """
        table_code = 3
        """
        if self.table_name in [b'OQG1', b'OQG2', b'OQGV1', b'OQP1']:
            pass
        else:
            msg = 'spc_forces; table_name=%s' % self.table_name
            raise NotImplementedError(msg)

        op2_reader = self.op2_reader
        thermal = op2_reader.thermal
        if thermal == 0:
            op2_reader._setup_op2_subcase('SPCFORCES')

            result_name = 'spc_forces'
            storage_obj = self.get_result(result_name)
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            n = op2_reader._read_table_vectorized(data, ndata, result_name, storage_obj,
                                                  RealSPCForcesArray, ComplexSPCForcesArray,
                                                  'node', random_code=op2_reader.random_code)
        elif thermal == 1:
            #'finite element temperature gradients and fluxes'
            op2_reader._setup_op2_subcase('FLUX')

            result_name = 'thermal_gradient_and_flux'
            storage_obj = self.get_result(result_name)
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            n = op2_reader._read_table_vectorized(data, ndata, result_name, storage_obj,
                                                  RealTemperatureGradientAndFluxArray, None,
                                                  'node', random_code=op2_reader.random_code)
        elif self.thermal == 8:  # 4 ?
            result_name = 'spc_forces_scaled_response_spectra_NRL'
            storage_obj = self.get_result(result_name)
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            n = op2_reader._read_table_vectorized(data, ndata, result_name, storage_obj,
                                                  RealSPCForcesArray, ComplexSPCForcesArray,
                                                  'node', random_code=op2_reader.random_code)
        else:
            raise RuntimeError(op2_reader.code_information())
            #msg = 'thermal=%s' % self.thermal
            #return self._not_implemented_or_skip(data, ndata, msg)
        return n

    def _read_oqg_mpc_forces(self, data, ndata):
        """
        table_code = 39
        """
        self._setup_op2_subcase('MPCFORCES')

        if self.table_name in [b'OQMG1', b'OQMG2']:
            result_name = 'mpc_forces'
        elif self.table_name == b'RAQEATC':
            result_name = 'mpc_forces_RAQEATC'
        elif self.table_name == b'RAQCONS':
            result_name = 'mpc_forces_RAQCONS'
        elif self.table_name == b'RAQEATC':
            result_name = 'mpc_forces_RAQEATC'
        else:
            msg = 'mpc_forces; table_name=%s' % self.table_name
            raise NotImplementedError(msg)

        storage_obj = self.get_result(result_name)
        if self.thermal == 0:
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            n = self._read_table_vectorized(data, ndata, result_name, storage_obj,
                                            RealMPCForcesArray, ComplexMPCForcesArray,
                                            'node', random_code=self.random_code)
        else:
            raise RuntimeError(self.code_information())
            #msg = 'thermal=%s' % self.thermal
            #return self._not_implemented_or_skip(data, ndata, msg)
        assert isinstance(n, integer_types), 'table_name=%s n=%s' % (self.table_name, n)
        return n

    def _read_oqg_spc_psd(self, data, ndata):
        """
        table_code = 601/610/611
        """
        if self.thermal == 0:
            if self.table_code in [3, 603]:
                result_name = 'psd.spc_forces'
            else:
                raise RuntimeError(self.code_information())
            obj = RealSPCForcesArray
        else:
            raise RuntimeError(self.code_information())

        if self._results.is_not_saved(result_name):
            return ndata
        self._results._found_result(result_name)

        storage_obj = self.get_result(result_name)
        n = self._read_random_table(data, ndata, result_name, storage_obj,
                                    obj, 'node',
                                    random_code=self.random_code)
        return n

    def _read_oqg_spc_rms(self, data, ndata):
        """
        table_code = 3/???/?10/?11
        """
        if self.thermal == 0:
            if self.table_code in [3, 803]:
                result_name = 'rms.spc_forces'
                if self._results.is_not_saved(result_name):
                    return ndata
                storage_obj = self.get_result(result_name)
                self._results._found_result(result_name)
                n = self._read_random_table(data, ndata, result_name, storage_obj,
                                            RealSPCForcesArray, 'node',
                                            random_code=self.random_code)
            else:
                raise RuntimeError(self.code_information())
        else:
            raise RuntimeError(self.code_information())
        return n

    def _read_oqg_spc_ato(self, data, ndata):
        """
        table_code = ???
        """
        if self.thermal == 0:
            if self.table_code in [3]:
                assert self.table_name in [b'OQGATO2'], 'self.table_name=%r' % self.table_name
                result_name = 'ato.spc_forces'
            #elif self.table_code in [603]:
                #assert self.table_name in [b'OQGATO2'], 'self.table_name=%r' % self.table_name
                #result_name = 'psd.mpc_forces'
            else:
                print(self.table_code)
                raise RuntimeError(self.code_information())

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            storage_obj = self.get_result(result_name)
            n = self._read_random_table(data, ndata, result_name, storage_obj,
                                        RealSPCForcesArray, 'node',
                                        random_code=self.random_code)
        else:
            raise RuntimeError(self.code_information())
        assert n is not None, n
        return n

    def _read_oqg_spc_crm(self, data, ndata):
        """
        table_code = 3/???/?10/?11
        """
        if self.thermal == 0:
            if self.table_code in [3, 503]:
                result_name = 'crm.spc_forces'
                storage_obj = self.get_result(result_name)
                if self._results.is_not_saved(result_name):
                    return ndata
                self._results._found_result(result_name)
                n = self._read_random_table(data, ndata, result_name, storage_obj,
                                            RealSPCForcesArray, 'node',
                                            random_code=self.random_code)
            else:
                raise RuntimeError(self.code_information())
        else:
            raise RuntimeError(self.code_information())
        return n

    def _read_oqg_spc_no(self, data, ndata):
        """
        table_code = 3/???/?10/?11
        """
        if self.thermal == 0:
            if self.table_code in [3, 903]:
                assert self.table_name in [b'OQGNO1', b'OQGNO2'], 'self.table_name=%r' % self.table_name
                result_name = 'no.spc_forces'
                obj = RealSPCForcesArray
            else:
                raise RuntimeError(self.code_information())
        else:
            raise RuntimeError(self.code_information())

        if self._results.is_not_saved(result_name):
            return ndata
        self._results._found_result(result_name)

        storage_obj = self.get_result(result_name)
        n = self._read_random_table(data, ndata, result_name, storage_obj,
                                    obj, 'node',
                                    random_code=self.random_code)
        return n

    def _read_oqg_mpc_psd(self, data, ndata):
        """
        table_code = 603
        """
        if self.thermal == 0:
            if self.table_code in [39]:
                assert self.table_name in [b'OQMPSD1', b'OQMPSD2'], 'self.table_name=%r' % self.table_name
                result_name = 'psd.mpc_forces'
            elif self.table_code in [603]:
                assert self.table_name in [b'OQMPSD1', b'OQMPSD2'], 'self.table_name=%r' % self.table_name
                result_name = 'psd.mpc_forces'
            else:
                print(self.table_code)
                raise RuntimeError(self.code_information())
            obj = RealMPCForcesArray
        else:
            raise RuntimeError(self.code_information())

        if self._results.is_not_saved(result_name):
            return ndata
        self._results._found_result(result_name)

        storage_obj = self.get_result(result_name)
        n = self._read_random_table(data, ndata, result_name, storage_obj,
                                    obj, 'node',
                                    random_code=self.random_code)

        assert n is not None, n
        return n

    def _read_oqg_mpc_ato(self, data, ndata):
        """
        table_code = ???
        """
        if self.thermal == 0:
            if self.table_code in [39]:
                assert self.table_name in [b'OQMATO1', b'OQMATO2'], 'self.table_name=%r' % self.table_name
                result_name = 'ato.mpc_forces'
                obj = RealMPCForcesArray
            #elif self.table_code in [603]:
                #assert self.table_name in [b''], 'self.table_name=%r' % self.table_name
                #result_name = 'psd.mpc_forces'
            else:
                print(self.table_code)
                raise RuntimeError(self.code_information())
        else:
            raise RuntimeError(self.code_information())

        if self._results.is_not_saved(result_name):
            return ndata
        self._results._found_result(result_name)

        storage_obj = self.get_result(result_name)
        n = self._read_random_table(data, ndata, result_name, storage_obj,
                                    obj, 'node',
                                    random_code=self.random_code)

        assert n is not None, n
        return n

    def _read_oqg_mpc_crm(self, data, ndata):
        """
        table_code = ???
        """
        if self.thermal == 0:
            if self.table_code in [39]:
                assert self.table_name in [b'OQMCRM1', b'OQMCRM2'], 'self.table_name=%r' % self.table_name
                result_name = 'crm.mpc_forces'
                obj = RealMPCForcesArray
            #elif self.table_code in [603]:
                #assert self.table_name in [b''], 'self.table_name=%r' % self.table_name
                #result_name = 'psd.mpc_forces'
            else:
                print(self.table_code)
                raise RuntimeError(self.code_information())
        else:
            raise RuntimeError(self.code_information())

        if self._results.is_not_saved(result_name):
            return ndata
        self._results._found_result(result_name)

        storage_obj = self.get_result(result_name)
        n = self._read_random_table(data, ndata, result_name, storage_obj,
                                    obj, 'node',
                                    random_code=self.random_code)
        assert n is not None, n
        return n

    def _read_oqg_mpc_rms(self, data, ndata):
        """
        table_code = ???
        """
        if self.thermal == 0:
            if self.table_code in [39]:
                assert self.table_name in [b'OQMRMS1', b'OQMRMS2'], 'self.table_name=%r' % self.table_name
                result_name = 'rms.mpc_forces'
                obj = RealMPCForcesArray
            #elif self.table_code in [603]:
                #assert self.table_name in [b''], 'self.table_name=%r' % self.table_name
                #result_name = 'psd.mpc_forces'
            else:
                print(self.table_code)
                raise RuntimeError(self.code_information())
        else:
            raise RuntimeError(self.code_information())

        if self._results.is_not_saved(result_name):
            return ndata
        self._results._found_result(result_name)

        storage_obj = self.get_result(result_name)
        n = self._read_random_table(data, ndata, result_name, storage_obj,
                                    obj, 'node',
                                    random_code=self.random_code)

        assert n is not None, n
        return n

    def _read_oqg_mpc_no(self, data, ndata):
        """
        table_code = ???
        """
        if self.thermal == 0:
            if self.table_code in [39]:
                assert self.table_name in [b'OQMNO1', b'OQMNO2'], 'self.table_name=%r' % self.table_name
                result_name = 'no.mpc_forces'
                obj = RealMPCForcesArray
            #elif self.table_code in [603]:
                #assert self.table_name in [b''], 'self.table_name=%r' % self.table_name
                #result_name = 'psd.mpc_forces'
            else:
                print(self.table_code)
                raise RuntimeError(self.code_information())
        else:
            raise RuntimeError(self.code_information())

        if self._results.is_not_saved(result_name):
            return ndata
        self._results._found_result(result_name)

        storage_obj = self.get_result(result_name)
        n = self._read_random_table(data, ndata, result_name, storage_obj,
                                    obj, 'node',
                                    random_code=self.random_code)
        assert n is not None, n
        return n
