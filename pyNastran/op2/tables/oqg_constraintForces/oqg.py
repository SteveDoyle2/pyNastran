"""
This file defines the OUG Table, which contains:
 * Real/Complex SPC Forces
   - SPCFORCE = ALL
 * Real/Complex MPC Forces
   - MPCFORCE = ALL
 * Real Temperature Gradient & Flux
   - FLUX = ALL
"""
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


class OQG(OP2Common):
    def __init__(self):
        OP2Common.__init__(self)

    def _read_oqg1_3(self, data, ndata):
        self.nonlinear_factor = None
        self.is_table_1 = True
        self.is_table_2 = False
        three = self.parse_approach_code(data)
        self.words = [
            'analysis_code', 'table_code', '???', 'isubcase',
            '???', '???', '???', 'random_code',
            'format_code', 'num_wide', '11', '???',
            'acoustic_flag', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', 'thermal', '???',
            '???', 'Title', 'subtitle', 'label']

        ## random code
        self.random_code = self.add_data_parameter(data, 'random_code', 'i', 8, False)

        ## format code
        self.format_code = self.add_data_parameter(data, 'format_code', 'i', 9, False)

        ## number of words per entry in record
        self.num_wide = self.add_data_parameter(data, 'num_wide', 'i', 10, False)

        ## acoustic pressure flag
        self.acoustic_flag = self.add_data_parameter(data, 'acoustic_flag', 'f', 13, False)

        ## thermal flag; 1 for heat transfer, 0 otherwise
        self.thermal = self.add_data_parameter(data, 'thermal', 'i', 23, False)

        if not self.is_sort1():
            raise NotImplementedError('SORT2; code_info=\n%s' % self.code_information())
        #assert self.isThermal()==False,self.thermal

        ## assuming tCode=1
        if self.analysis_code == 1:   # statics / displacement / heat flux
            ## load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5, False)
            self.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysis_code == 2:  # real eigenvalues
            ## mode number
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)
            ## eigenvalue
            self.eign = self.add_data_parameter(data, 'eign', 'f', 6, False)
            ## mode or cycle .. todo:: confused on the type - F1???
            self.mode_cycle = self.add_data_parameter(data, 'mode_cycle', 'f', 7, False)
            self.update_mode_cycle('mode_cycle')
            self.data_names = self.apply_data_code_value('data_names',
                                                         ['mode', 'eign', 'mode_cycle'])
        #elif self.analysis_code == 3: # differential stiffness
            #self.lsdvmn = self.get_values(data,'i',5) ## load set number
            #self.data_names = self.data_code['lsdvmn'] = self.lsdvmn
        #elif self.analysis_code == 4: # differential stiffness
            #self.lsdvmn = self.get_values(data,'i',5) ## load set number
        elif self.analysis_code == 5:   # frequency
            ## frequency
            self.freq = self.add_data_parameter(data, 'freq', 'f', 5)
            self.data_names = self.apply_data_code_value('data_names', ['freq'])
        elif self.analysis_code == 6:  # transient
            ## time step
            self.dt = self.add_data_parameter(data, 'dt', 'f', 5)
            self.data_names = self.apply_data_code_value('data_names', ['dt'])
        elif self.analysis_code == 7:  # pre-buckling
            ## load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
        elif self.analysis_code == 8:  # post-buckling
            ## load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            ## real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            self.data_names = self.apply_data_code_value('data_names', ['lsdvmn', 'eigr'])
        elif self.analysis_code == 9:  # complex eigenvalues
            ## mode number
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            ## imaginary eigenvalue
            self.eigi = self.add_data_parameter(data, 'eigi', 'f', 7, False)
            self.data_names = self.apply_data_code_value('data_names', ['mode', 'eigr', 'eigi'])
        elif self.analysis_code == 10:  # nonlinear statics
            ## load step
            self.lftsfq = self.add_data_parameter(data, 'lftsfq', 'f', 5)
            self.data_names = self.apply_data_code_value('data_names', ['lftsfq'])
        elif self.analysis_code == 11:  # old geometric nonlinear statics
            ## load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
        elif self.analysis_code == 12:
            # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
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

    def _read_oqg2_3(self, data, ndata):
        self.nonlinear_factor = None
        self.is_table_1 = False
        self.is_table_2 = True
        three = self.parse_approach_code(data)
        self.words = [
            'analysis_code', 'table_code', '???', 'isubcase',
            '???', '???', '???', 'random_code',
            'format_code', 'num_wide', '11', '???',
            'acoustic_flag', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', 'thermal', '???',
            '???', 'Title', 'subtitle', 'label']

        ## random code
        self.random_code = self.add_data_parameter(data, 'random_code', 'i', 8, False)

        ## format code
        self.format_code = self.add_data_parameter(data, 'format_code', 'i', 9, False)

        ## number of words per entry in record
        self.num_wide = self.add_data_parameter(data, 'num_wide', 'i', 10, False)

        ## acoustic pressure flag
        self.acoustic_flag = self.add_data_parameter(data, 'acoustic_flag', 'f', 13, False)

        ## thermal flag; 1 for heat transfer, 0 otherwise
        self.thermal = self.add_data_parameter(data, 'thermal', 'i', 23, False)

        #assert self.isThermal() == False, self.thermal

        self.node_id = self.add_data_parameter(data, 'node_id', 'i', 5, fix_device_code=True)
        #if self.analysis_code == 1:  # statics / displacement / heat flux
            ## load set number
            #self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5, False)
            #self.data_names = self.apply_data_code_value('data_names', ['node_id'])
            #self.setNullNonlinearFactor()
        if self.analysis_code == 1:  # static...because reasons.
            self._analysis_code_fmt = 'i'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
        elif self.analysis_code == 2:  # real eigenvalues
            ## mode number
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)
            self._analysis_code_fmt = 'i'
            ## real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            ## mode or cycle .. todo:: confused on the type - F1???
            self.mode_cycle = self.add_data_parameter(data, 'mode_cycle', 'f', 7, False)
            self.data_names = self.apply_data_code_value('data_names',
                                                         ['node_id', 'eigr', 'mode_cycle'])
        #elif self.analysis_code == 3: # differential stiffness
            #self.lsdvmn = self.get_values(data, 'i', 5) ## load set number
            #self.data_names = self.data_code['lsdvmn'] = self.lsdvmn
        #elif self.analysis_code == 4: # differential stiffness
            #self.lsdvmn = self.get_values(data, 'i', 5) ## load set number
        elif self.analysis_code == 5:   # frequency
            ## frequency
            #self.freq = self.add_data_parameter(data, 'freq', 'f', 5)
            self._analysis_code_fmt = 'f'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
        elif self.analysis_code == 6:  # transient
            ## time step
            #self.dt = self.add_data_parameter(data, 'dt', 'f', 5)
            self._analysis_code_fmt = 'f'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
        elif self.analysis_code == 7:  # pre-buckling
            ## load set number
            #self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self._analysis_code_fmt = 'i'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
        elif self.analysis_code == 8:  # post-buckling
            ## load set number
            #self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self._analysis_code_fmt = 'i'
            ## real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            self.data_names = self.apply_data_code_value('data_names', ['node_id', 'eigr'])
        elif self.analysis_code == 9:  # complex eigenvalues
            ## mode number
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)
            self._analysis_code_fmt = 'i'
            ## real eigenvalue
            #self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            ## imaginary eigenvalue
            self.eigi = self.add_data_parameter(data, 'eigi', 'f', 7, False)
            self.data_names = self.apply_data_code_value('data_names', ['node_id', 'eigr', 'eigi'])
        elif self.analysis_code == 10:  # nonlinear statics
            ## load step
            #self.lftsfq = self.add_data_parameter(data, 'lftsfq', 'f', 5)
            self._analysis_code_fmt = 'f'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
        elif self.analysis_code == 11:  # old geometric nonlinear statics
            ## load set number
            #self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self._analysis_code_fmt = 'f'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
        elif self.analysis_code == 12:
            # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## load set number
            #self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self._analysis_code_fmt = 'i'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
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

        if self.table_name in [b'OQGCRM1', b'OQGCRM2']:
            if self.table_code not in [3, 503]:
                msg = 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
                raise AssertionError(msg)
            n = self._read_oqg_spc_crm(data, ndata)
        elif self.table_name in [b'OQGPSD1', b'OQGPSD2']:
            if self.table_code not in [3, 603]:  # was 3
                msg = 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
                raise AssertionError(msg + '\n%s' % self.code_information())
            n = self._read_oqg_spc_psd(data, ndata)
        elif self.table_name in [b'OQGATO1', b'OQGATO2']:
            if self.table_code not in [3, 703]:
                msg = 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
                raise AssertionError(msg)
            n = self._read_oqg_spc_ato(data, ndata)
        elif self.table_name in [b'OQGRMS1', b'OQGRMS2']:
            if self.table_code not in [3, 803]:
                msg = 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
                raise AssertionError(msg)
            n = self._read_oqg_spc_rms(data, ndata)
        elif self.table_name in [b'OQGNO1', b'OQGNO2']:
            if self.table_code not in [3, 903]:
                msg = 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
                raise AssertionError(msg)
            n = self._read_oqg_spc_no(data, ndata)

        elif self.table_code == 3:   # SPC Forces
            assert self.table_name in [b'OQG1', b'OQGV1', b'OQP1', b'OQG2'], self.code_information()
            n = self._read_spc_forces(data, ndata)
        elif self.table_code == 39:  # MPC Forces
            assert self.table_name in [b'OQMG1', b'OQMG2'], self.code_information() # , b'OQMPSD1', b'OQMPSD2'
            n = self._read_oqg_mpc_forces(data, ndata)
        elif self.table_name in [b'RAQEATC', b'RAQCONS']:
            # self.table_code == 5 and
            n = self._read_oqg_mpc_forces(data, ndata)
        else:
            raise RuntimeError(self.code_information())
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

        if self.thermal == 0:
            if self.isubcase not in self.case_control_deck.subcases:
                self.subcase = self.case_control_deck.create_new_subcase(self.isubcase)
            self.subcase.add_op2_data(self.data_code, 'SPCFORCES', self.log)

            result_name = 'spc_forces'
            storage_obj = self.spc_forces
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            n = self._read_table_vectorized(data, ndata, result_name, storage_obj,
                                            RealSPCForcesArray, ComplexSPCForcesArray,
                                            'node', random_code=self.random_code)
        elif self.thermal == 1:
            #'finite element temperature gradients and fluxes'
            if self.isubcase not in self.case_control_deck.subcases:
                self.subcase = self.case_control_deck.create_new_subcase(self.isubcase)
            self.subcase.add_op2_data(self.data_code, 'FLUX', self.log)

            result_name = 'thermal_gradient_and_flux'
            storage_obj = self.thermal_gradient_and_flux
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            n = self._read_table_vectorized(data, ndata, result_name, storage_obj,
                                            RealTemperatureGradientAndFluxArray, None,
                                            'node', random_code=self.random_code)
        elif self.thermal == 8:  # 4 ?
            result_name = 'spc_forces_scaled_response_spectra_NRL'
            storage_obj = self.spc_forces_scaled_response_spectra_NRL
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            n = self._read_table_vectorized(data, ndata, result_name, storage_obj,
                                            RealSPCForcesArray, ComplexSPCForcesArray,
                                            'node', random_code=self.random_code)
        else:
            raise RuntimeError(self.code_information())
            #msg = 'thermal=%s' % self.thermal
            #return self._not_implemented_or_skip(data, ndata, msg)
        return n

    def _read_oqg_mpc_forces(self, data, ndata):
        """
        table_code = 39
        """
        if self.isubcase not in self.case_control_deck.subcases:
            self.subcase = self.case_control_deck.create_new_subcase(self.isubcase)
        self.subcase.add_op2_data(self.data_code, 'MPCFORCES', self.log)

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

        storage_obj = getattr(self, result_name)
        if self.thermal == 0:
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            n = self._read_table_vectorized(data, ndata, result_name, storage_obj,
                                            RealMPCForcesArray, ComplexMPCForcesArray,
                                            'node', random_code=self.random_code)
        #elif self.thermal == 1:
            #raise NotImplementedError(self.thermal)
            #n = self._read_table(data, result_name, storage_obj,
                                 #None, None,
                                 #None, None, 'node')
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
                result_name = 'spc_forces_PSD'
                storage_obj = self.spc_forces_PSD
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

    def _read_oqg_spc_rms(self, data, ndata):
        """
        table_code = 3/???/?10/?11
        """
        if self.thermal == 0:
            if self.table_code in [3, 803]:
                result_name = 'spc_forces_RMS'
                storage_obj = self.spc_forces_RMS
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

    def _read_oqg_spc_ato(self, data, ndata):
        """
        table_code = ???
        """
        if self.thermal == 0:
            if self.table_code in [3]:
                assert self.table_name in [b'OQGATO2'], 'self.table_name=%r' % self.table_name
                result_name = 'spc_forces_ATO'
            #elif self.table_code in [603]:
                #assert self.table_name in [b'OQGATO2'], 'self.table_name=%r' % self.table_name
                #result_name = 'mpc_forces_PSD'
            else:
                print(self.table_code)
                raise RuntimeError(self.code_information())

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            storage_obj = getattr(self, result_name)
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
                result_name = 'spc_forces_CRM'
                storage_obj = self.spc_forces_CRM
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
                result_name = 'spc_forces_NO'
                storage_obj = self.spc_forces_NO
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

    def _read_oqg_mpc_psd(self, data, ndata):
        """
        table_code = 603
        """
        if self.thermal == 0:
            if self.table_code in [39]:
                assert self.table_name in [b'OQMPSD1', b'OQMPSD2'], 'self.table_name=%r' % self.table_name
                result_name = 'mpc_forces_PSD'
            elif self.table_code in [603]:
                assert self.table_name in [b'OQMPSD1', b'OQMPSD2'], 'self.table_name=%r' % self.table_name
                result_name = 'mpc_forces_PSD'
            else:
                print(self.table_code)
                raise RuntimeError(self.code_information())

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            storage_obj = getattr(self, result_name)
            n = self._read_random_table(data, ndata, result_name, storage_obj,
                                        RealMPCForcesArray, 'node',
                                        random_code=self.random_code)
        else:
            raise RuntimeError(self.code_information())
        assert n is not None, n
        return n

    def _read_oqg_mpc_ato(self, data, ndata):
        """
        table_code = ???
        """
        if self.thermal == 0:
            if self.table_code in [39]:
                assert self.table_name in [b'OQMATO1', b'OQMATO2'], 'self.table_name=%r' % self.table_name
                result_name = 'mpc_forces_ATO'
            #elif self.table_code in [603]:
                #assert self.table_name in [b''], 'self.table_name=%r' % self.table_name
                #result_name = 'mpc_forces_PSD'
            else:
                print(self.table_code)
                raise RuntimeError(self.code_information())

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            storage_obj = getattr(self, result_name)
            n = self._read_random_table(data, ndata, result_name, storage_obj,
                                        RealMPCForcesArray, 'node',
                                        random_code=self.random_code)
        else:
            raise RuntimeError(self.code_information())
        assert n is not None, n
        return n

    def _read_oqg_mpc_crm(self, data, ndata):
        """
        table_code = ???
        """
        if self.thermal == 0:
            if self.table_code in [39]:
                assert self.table_name in [b'OQMCRM1', b'OQMCRM2'], 'self.table_name=%r' % self.table_name
                result_name = 'mpc_forces_CRM'
            #elif self.table_code in [603]:
                #assert self.table_name in [b''], 'self.table_name=%r' % self.table_name
                #result_name = 'mpc_forces_PSD'
            else:
                print(self.table_code)
                raise RuntimeError(self.code_information())

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            storage_obj = getattr(self, result_name)
            n = self._read_random_table(data, ndata, result_name, storage_obj,
                                        RealMPCForcesArray, 'node',
                                        random_code=self.random_code)
        else:
            raise RuntimeError(self.code_information())
        assert n is not None, n
        return n

    def _read_oqg_mpc_rms(self, data, ndata):
        """
        table_code = ???
        """
        if self.thermal == 0:
            if self.table_code in [39]:
                assert self.table_name in [b'OQMRMS1', b'OQMRMS2'], 'self.table_name=%r' % self.table_name
                result_name = 'mpc_forces_RMS'
            #elif self.table_code in [603]:
                #assert self.table_name in [b''], 'self.table_name=%r' % self.table_name
                #result_name = 'mpc_forces_PSD'
            else:
                print(self.table_code)
                raise RuntimeError(self.code_information())

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            storage_obj = getattr(self, result_name)
            n = self._read_random_table(data, ndata, result_name, storage_obj,
                                        RealMPCForcesArray, 'node',
                                        random_code=self.random_code)
        else:
            raise RuntimeError(self.code_information())
        assert n is not None, n
        return n

    def _read_oqg_mpc_no(self, data, ndata):
        """
        table_code = ???
        """
        if self.thermal == 0:
            if self.table_code in [39]:
                assert self.table_name in [b'OQMNO1', b'OQMNO2'], 'self.table_name=%r' % self.table_name
                result_name = 'mpc_forces_CRM'
            #elif self.table_code in [603]:
                #assert self.table_name in [b''], 'self.table_name=%r' % self.table_name
                #result_name = 'mpc_forces_PSD'
            else:
                print(self.table_code)
                raise RuntimeError(self.code_information())

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            storage_obj = getattr(self, result_name)
            n = self._read_random_table(data, ndata, result_name, storage_obj,
                                        RealMPCForcesArray, 'node',
                                        random_code=self.random_code)
        else:
            raise RuntimeError(self.code_information())
        assert n is not None, n
        return n
