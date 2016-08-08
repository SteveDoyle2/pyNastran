"""
This file defines the OUG Table, which contains:
 * Real/Complex SPC Forces
   - SPCFORCE = ALL
 * Real/Complex MPC Forces
   - MPCFORCE = ALL
 * Real Temperature Gradient & Flux
   - FLUX = ALL
"""
from pyNastran.op2.op2_common import OP2Common

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
            raise NotImplementedError('sort2...')
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
            ## real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            ## mode or cycle .. todo:: confused on the type - F1???
            self.mode_cycle = self.add_data_parameter(data, 'mode_cycle', 'f', 7, False)
            self.update_mode_cycle('mode_cycle')
            self.data_names = self.apply_data_code_value('data_names', ['mode', 'eigr', 'mode_cycle'])
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
        elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
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
        if self.analysis_code == 2:  # real eigenvalues
            ## mode number
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)
            self._analysis_code_fmt = 'i'
            ## real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            ## mode or cycle .. todo:: confused on the type - F1???
            self.mode_cycle = self.add_data_parameter(data, 'mode_cycle', 'f', 7, False)
            self.data_names = self.apply_data_code_value('data_names', ['node_id', 'eigr', 'mode_cycle'])
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
        elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
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
        result_name = 'constraint_forces'
        if self._results.is_not_saved(result_name):
            return ndata

        if self.is_msc:
            if self.table_name == b'OQGPSD2':
                if self.table_code not in [3]:
                    msg = 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
                    raise AssertionError(msg)
                n = self._read_oqg_spc_psd(data, ndata)

            elif self.table_code == 3:   # SPC Forces
                assert self.table_name in [b'OQG1', b'OQGV1', b'OQP1'], self.code_information()
                n = self._read_spc_forces(data, ndata)
            elif self.table_code == 39:  # MPC Forces
                assert self.table_name in [b'OQMG1'], self.code_information()
                n = self._read_mpc_forces(data, ndata)
            else:
                raise RuntimeError(self.code_information())
                #msg = self.code_information()
                #return self._not_implemented_or_skip(data, ndata, msg)
        elif self.is_nx:
            if self.table_name == b'OQMPSD2':
                if self.table_code not in [603]:
                    msg = 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
                    raise AssertionError(msg)
                n = self._read_oqg_mpc_psd(data, ndata)

            elif self.table_code == 3:   # SPC Forces
                if self.table_name in [b'OQG1', b'OQG2', b'OQGV1', b'OQP1']:
                    n = self._read_spc_forces(data, ndata)
                elif self.table_name in [b'OQMG1']:
                    n = self._read_mpc_forces(data, ndata)
                else:
                    raise RuntimeError(self.code_information())
                    #msg = self.code_information()
                    #return self._not_implemented_or_skip(data, ndata, msg)
            elif self.table_code == 5 and self.table_name in [b'RAQEATC', b'RAQCONS']:
                n = self._read_mpc_forces(data, ndata)
            else:
                raise RuntimeError(self.code_information())
        else:
            raise RuntimeError(self.code_information())
            #msg = 'table_code=%s' % self.table_code
            #return self._not_implemented_or_skip(data, ndata, msg)
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
            self.subcase.add_op2_data(self.data_code, 'SPCFORCES, self.log')

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
            msg = 'thermal=%s' % self.thermal
            return self._not_implemented_or_skip(data, ndata, msg)
        return n

    def _read_mpc_forces(self, data, ndata):
        """
        table_code = 39
        """
        if self.isubcase not in self.case_control_deck.subcases:
            self.subcase = self.case_control_deck.create_new_subcase(self.isubcase)
        self.subcase.add_op2_data(self.data_code, 'MPCFORCES', self.log)

        if self.table_name == b'OQMG1':
            result_name = 'mpc_forces'
        elif self.table_name == b'RAQEATC':
            result_name = 'mpc_forces_RAQEATC'
        elif self.table_name == b'RAQCONS':
            result_name = 'mpc_forces_RAQCONS'
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
            msg = 'thermal=%s' % self.thermal
            return self._not_implemented_or_skip(data, ndata, msg)
        return n

    def _read_oqg_spc_psd(self, data, ndata):
        """
        table_code = 601/610/611
        """
        if self.thermal == 0:
            if self.table_code in [3]:
                result_name = 'spc_forcesPSD'
                storage_obj = self.displacementsPSD
                if self._results.is_not_saved(result_name):
                    return ndata
                self._results._found_result(result_name)
                n = self._read_table_vectorized(data, ndata, result_name, storage_obj,
                                                RealSPCForcesArray, ComplexSPCForcesArray,
                                                'node', random_code=self.random_code)
            #elif self.table_code == 610:
                #result_name = 'velocitiesPSD'
                #storage_obj = self.velocitiesPSD
                #if self._results.is_not_saved(result_name):
                    #return ndata
                #self._results._found_result(result_name)
                #n = self._read_table_vectorized(data, ndata, result_name, storage_obj,
                                                #RealVelocityArray, ComplexVelocityArray,
                                                #'node', random_code=self.random_code)
            #elif self.table_code == 611:
                #result_name = 'accelerationsPSD'
                #storage_obj = self.accelerationsPSD
                #if self._results.is_not_saved(result_name):
                    #return ndata
                #n = self._read_table_vectorized(data, ndata, result_name, storage_obj,
                                                #RealAccelerationArray, ComplexAccelerationArray,
                                                #'node', random_code=self.random_code)
            #elif self.table_code in [1]:
                #if self.format_code == 2:
                    #self.format_code = 1
                    #self.data['format_code'] = 1
                #result_name = 'displacements'
                #storage_obj = self.displacements
                #if self._results.is_not_saved(result_name):
                    #return ndata
                #self._results._found_result(result_name)
                #n = self._read_table_vectorized(data, ndata, result_name, storage_obj,
                                     #RealDisplacementArray, ComplexDisplacementArray,
                                     #'node', random_code=self.random_code)
            else:
                raise RuntimeError(self.code_information())
        else:
            raise RuntimeError(self.code_information())

    def _read_oqg_mpc_psd(self, data, ndata):
        """
        table_code = 603
        """
        if self.thermal == 0:
            if self.table_code in [603]:
                result_name = 'mpc_forcesPSD'
                storage_obj = self.mpc_forcesPSD
                if self._results.is_not_saved(result_name):
                    return ndata
                self._results._found_result(result_name)
                n = self._read_table_vectorized(data, ndata, result_name, storage_obj,
                                                RealMPCForcesArray, ComplexMPCForcesArray,
                                                'node', random_code=self.random_code)
            else:
                raise RuntimeError(self.code_information())
        else:
            raise RuntimeError(self.code_information())

        assert n is not None, n
        return n
