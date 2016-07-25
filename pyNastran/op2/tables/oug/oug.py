#pylint: disable=C0301,C0103
"""
This file defines the OUG Table, which contains:
 * Real/Complex Displacement
   - DISPLACEMENT = ALL
 * Real/Complex Acceleration
   - ACCELERATION = ALL
 * Real/Complex Velocity
   - VELOCITY = ALL
 * Real/Complex Eigenvectors
   - DISPLACEMENT = ALL
 * Real Temperature
   - DISPLACEMENT = ALL
"""
import numpy as np
from pyNastran import is_release
from pyNastran.op2.op2_common import OP2Common

from pyNastran.op2.tables.oug.oug_displacements import (
    RealDisplacementArray, ComplexDisplacementArray)

from pyNastran.op2.tables.oug.oug_velocities import (
    RealVelocityArray, ComplexVelocityArray)

from pyNastran.op2.tables.oug.oug_accelerations import (
    RealAccelerationArray, ComplexAccelerationArray)

from pyNastran.op2.tables.oug.oug_temperatures import (
    RealTemperatureArray)

from pyNastran.op2.tables.oug.oug_eigenvectors import (
    RealEigenvectorArray, ComplexEigenvectorArray,
)

from pyNastran.op2.tables.opg_appliedLoads.opg_load_vector import RealThermalVelocityVectorArray


class OUG(OP2Common):

    def __init__(self):
        OP2Common.__init__(self)

    def update_mode_cycle(self, name):
        value = getattr(self, name)
        if value == 0.0:
            #print('table_name=%r mode=%s eigr=%s' % (self.table_name, self.mode, self.eigr))
            value = np.sqrt(np.abs(self.eigr)) / (2. * np.pi)
            setattr(self, name, value)
            self.data_code[name] = value

    def _read_oug1_3(self, data, ndata):
        #self._set_times_dtype()
        self.nonlinear_factor = None
        self.is_table_1 = True
        self.is_table_2 = False
        three = self.parse_approach_code(data)
        self.words = [
            'approach_code', 'table_code', '???', 'isubcase',
            '???', '???', '???', 'random_code',
            'format_code', 'num_wide', '???', '???',
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
        self.acoustic_flag = self.add_data_parameter(data, 'acoustic_flag', 'i', 13, False)

        ## thermal flag; 1 for heat transfer, 0 otherwise
        self.thermal = self.add_data_parameter(data, 'thermal', 'i', 23, False)

        if self.analysis_code == 1:   # statics / displacement / heat flux
            # load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5, False)
            self.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysis_code == 2:  # real eigenvalues
            # mode number
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)
            # real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            # mode or cycle .. todo:: confused on the type - F1???
            self.mode_cycle = self.add_data_parameter(data, 'mode_cycle', 'i', 7, False)
            self.update_mode_cycle('mode_cycle')
            self.data_names = self.apply_data_code_value('data_names', ['mode', 'eigr', 'mode_cycle'])
        #elif self.analysis_code == 3: # differential stiffness
            #self.lsdvmn = self.get_values(data, 'i', 5) ## load set number
            #self.data_code['lsdvmn'] = self.lsdvmn
        #elif self.analysis_code == 4: # differential stiffness
            #self.lsdvmn = self.get_values(data, 'i', 5) ## load set number
        elif self.analysis_code == 5:   # frequency
            # frequency
            self.freq = self.add_data_parameter(data, 'freq', 'f', 5)
            self.data_names = self.apply_data_code_value('data_names', ['freq'])
        elif self.analysis_code == 6:  # transient
            # time step
            self.dt = self.add_data_parameter(data, 'dt', 'f', 5)
            self.data_names = self.apply_data_code_value('data_names', ['dt'])
        elif self.analysis_code == 7:  # pre-buckling
            # load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
        elif self.analysis_code == 8:  # post-buckling
            # load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            # real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            self.data_names = self.apply_data_code_value('data_names', ['lsdvmn', 'eigr'])
        elif self.analysis_code == 9:  # complex eigenvalues
            # mode number
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)
            # real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            # imaginary eigenvalue
            self.eigi = self.add_data_parameter(data, 'eigi', 'f', 7, False)
            self.data_names = self.apply_data_code_value('data_names', ['mode', 'eigr', 'eigi'])
        elif self.analysis_code == 10:  # nonlinear statics
            # load step
            self.lftsfq = self.add_data_parameter(data, 'lftsfq', 'f', 5)
            self.data_names = self.apply_data_code_value('data_names', ['lftsfq'])
        elif self.analysis_code == 11:  # old geometric nonlinear statics
            # load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
        elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            # load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % self.analysis_code
            raise RuntimeError(msg)

        #print self.code_information()
        #
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


    def _read_oug2_3(self, data, ndata):
        #self._set_times_dtype()
        #return self._read_oug1_3(data)
        self.nonlinear_factor = None

        self.is_table_1 = False
        self.is_table_2 = True
        three = self.parse_approach_code(data)
        self.words = [
            'approach_code', 'table_code', '???', 'isubcase',
            '???', '???', '???', 'random_code',
            'format_code', 'num_wide', '???', '???',
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
        self.acoustic_flag = self.add_data_parameter(data, 'acoustic_flag', 'i', 13, False)

        ## thermal flag; 1 for heat transfer, 0 otherwise
        self.thermal = self.add_data_parameter(data, 'thermal', 'i', 23, False)

        self.node_id = self.add_data_parameter(data, 'node_id', 'i', 5, fix_device_code=True)
        #if self.analysis_code == 1:  # statics / displacement / heat flux
            # load set number
            #self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5, False)
            #self.data_names = self.apply_data_code_value('data_names', ['node_id'])
            #self.setNullNonlinearFactor()
        if self.analysis_code == 1:  # static...because reasons.
            self._analysis_code_fmt = 'i'
        elif self.analysis_code == 2:  # real eigenvalues
            # mode number
            #self.mode = self.add_data_parameter(data, 'mode', 'i', 5)
            self._analysis_code_fmt = 'i'
            # real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            self.mode_cycle = self.add_data_parameter(data, 'mode_cycle', 'i', 7, False)  # mode or cycle .. todo:: confused on the type - F1???
            self.data_names = self.apply_data_code_value('data_names', ['node_id', 'eigr', 'mode_cycle'])
        #elif self.analysis_code == 3: # differential stiffness
            #self.lsdvmn = self.get_values(data, 'i', 5) ## load set number
            #self.data_code['lsdvmn'] = self.lsdvmn
        #elif self.analysis_code == 4: # differential stiffness
            #self.lsdvmn = self.get_values(data, 'i', 5) ## load set number
        elif self.analysis_code == 5:   # frequency
            # frequency
            #self.freq = self.add_data_parameter(data, 'freq', 'f', 5)
            self._analysis_code_fmt = 'f'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
        elif self.analysis_code == 6:  # transient
            # time step
            #self.dt = self.add_data_parameter(data, 'dt', 'f', 5)
            self._analysis_code_fmt = 'f'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
        elif self.analysis_code == 7:  # pre-buckling
            # load set number
            #self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self._analysis_code_fmt = 'i'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
        elif self.analysis_code == 8:  # post-buckling
            # load set number
            #self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self._analysis_code_fmt = 'f'
            # real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            self.data_names = self.apply_data_code_value('data_names', ['node_id', 'eigr'])
        elif self.analysis_code == 9:  # complex eigenvalues
            # mode number
            #self.mode = self.add_data_parameter(data, 'mode', 'i', 5)
            self._analysis_code_fmt = 'i'
            # real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            # imaginary eigenvalue
            self.eigi = self.add_data_parameter(data, 'eigi', 'f', 7, False)
            self.data_names = self.apply_data_code_value('data_names', ['node_id', 'eigr', 'eigi'])
        elif self.analysis_code == 10:  # nonlinear statics
            # load step
            #self.lftsfq = self.add_data_parameter(data, 'lftsfq', 'f', 5)
            self._analysis_code_fmt = 'f'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
        elif self.analysis_code == 11:  # old geometric nonlinear statics
            # load set number
            #self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
        elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            # load set number
            #self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
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
            self.binary_debug.write('  %-14s = %r %s\n' % ('approach_code', self.approach_code,
                                                           self.approach_code_str(self.approach_code)))
            self.binary_debug.write('  %-14s = %r\n' % ('tCode', self.tCode))
            self.binary_debug.write('  %-14s = %r\n' % ('isubcase', self.isubcase))
        self._read_title(data)
        self._write_debug_bits()
        assert isinstance(self.nonlinear_factor, int), self.nonlinear_factor

    def _read_oug_4(self, data, ndata):
        if self.table_name == b'OUGPSD2':
            if self.table_code not in [1, 601, 610, 611]:
                msg = 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
                raise AssertionError(msg)
            n = self._read_oug_psd(data, ndata)
        elif self.table_name in [b'OUGRMS1', b'OUGRMS2']:
            if self.table_code not in [1, 801]:
                msg = 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
                raise AssertionError(msg)
            n = self._read_oug_rms(data, ndata)
        elif self.table_name in [b'OUGNO1', b'OUGNO2']:
            if self.table_code not in [1, 901]:
                msg = 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
                raise AssertionError(msg)
            n = self._read_oug_no(data, ndata)
        elif self.table_name in [b'OUGATO1', b'OUGATO2']:
            if self.table_code not in [1]:
                msg = 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
                raise AssertionError(msg)
            n = self._read_oug_ato(data, ndata)
        elif self.table_name in [b'OUGCRM2']:
            if self.table_code not in [1]:
                msg = 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
                raise AssertionError(msg)
            n = self._not_implemented_or_skip(data, ndata)

        elif self.table_code == 1:   # Displacements
            if self.table_name not in [b'OUG1', b'BOUGV1', b'OUGV1', b'OUGV1PAT', b'TOUGV1',
                                       b'OUGV2',
                                       b'OUPV1', b'ROUGV1']:
                msg = 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
                raise AssertionError(msg)
            is_cid = False
            if self.table_name == b'OUGV1PAT':
                is_cid = True
            n = self._read_displacement(data, ndata, is_cid)
        elif self.table_code == 7:
            n = self._read_eigenvector(data, ndata)
        elif self.table_code == 10:
            n = self._read_velocity(data, ndata)
        elif self.table_code == 11:
            n = self._read_acceleration(data, ndata)
        else:
            raise NotImplementedError(self.code_information())
        #else:
            #self._not_implemented_or_skip(data, ndata, 'bad OUG table')
        return n

    def _read_eigenvector_displacement_solution_set(self, data, ndata):
        """
        table_code = 14
        """
        raise NotImplementedError()

    def _read_displacement_solution_set(self, data, ndata):
        """
        table_code = 15
        """
        raise NotImplementedError()

    def _read_velocity_solution_set(self, data, ndata):
        """
        table_code = 16
        """
        raise NotImplementedError()

    def _read_acceleration_solution_set(self, data, ndata):
        """
        table_code = 17
        """
        raise NotImplementedError()

    def _read_displacement(self, data, ndata, is_cid):
        if self.read_mode == 1:
            if self.isubcase not in self.case_control_deck.subcases:
                self.subcase = self.case_control_deck.create_new_subcase(self.isubcase)
            self.subcase.add_op2_data(self.data_code, 'displacement', self.log)

        if self.table_name in [b'OUG1', b'OUGV1', b'OUGV2', b'OUGV1PAT', b'BOUGV1', b'ROUGV1']:
            assert self.thermal in [0, 1], self.code_information()
            if self.thermal == 0:
                result_name = 'displacements'
            elif self.thermal == 1:
                result_name = 'temperatures'
            else:
                msg = 'displacements; table_name=%s' % self.table_name
                raise NotImplementedError(msg)

        elif self.table_name in [b'OUPV1']:
            result_name = 'temperatures'
            assert self.thermal in [2, 4, 8], self.code_information()
            if self.thermal == 2:
                result_name = 'displacement_scaled_response_spectra_ABS'
            elif self.thermal == 4:
                result_name = 'displacement_scaled_response_spectra_SRSS'
            elif self.thermal == 8:
                result_name = 'displacement_scaled_response_spectra_NRL'
            else:
                msg = 'displacements; table_name=%s' % self.table_name
                raise NotImplementedError(msg)

        elif self.table_name in [b'TOUGV1']:
            result_name = 'temperatures'
            assert self.thermal == 1, self.code_information()
        else:
            msg = 'displacements; table_name=%s' % self.table_name
            raise NotImplementedError(msg)

        if self.thermal == 0:
            result_name = 'displacements'
            storage_obj = self.displacements
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            n = self._read_table_vectorized(data, ndata, result_name, storage_obj,
                                            RealDisplacementArray, ComplexDisplacementArray,
                                            'node', random_code=self.random_code,
                                            is_cid=is_cid)
        elif self.thermal == 1:
            result_name = 'temperatures'
            storage_obj = self.temperatures
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            n = self._read_scalar_table_vectorized(data, ndata, result_name, storage_obj,
                                                   RealTemperatureArray, None,
                                                   'node', random_code=self.random_code,
                                                   is_cid=is_cid)
        elif self.thermal == 2:
            result_name = 'displacement_scaled_response_spectra_ABS'
            storage_obj = self.displacement_scaled_response_spectra_ABS
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            n = self._read_table_vectorized(data, ndata, result_name, storage_obj,
                                            RealDisplacementArray, ComplexDisplacementArray,
                                            'node', random_code=self.random_code)
        elif self.thermal == 4:
            result_name = 'displacement_scaled_response_spectra_SRSS'
            storage_obj = self.displacement_scaled_response_spectra_SRSS
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            n = self._read_table_vectorized(data, ndata, result_name, storage_obj,
                                            RealDisplacementArray, ComplexDisplacementArray,
                                            'node', random_code=self.random_code)
        elif self.thermal == 8:  # 4 ?
            result_name = 'displacement_scaled_response_spectra_NRL'
            storage_obj = self.displacement_scaled_response_spectra_NRL
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            n = self._read_table_vectorized(data, ndata, result_name, storage_obj,
                                            RealDisplacementArray, ComplexDisplacementArray,
                                            'node', random_code=self.random_code)
            #return self._not_implemented_or_skip(data, ndata, msg='thermal=4')
        else:
            raise RuntimeError(self.code_information())
            n = self._not_implemented_or_skip(data, ndata, 'bad thermal=%r table' % self.thermal)
        #else:
            #raise NotImplementedError(self.thermal)
        return n

    def _read_velocity(self, data, ndata):
        """
        table_code = 10
        """
        if self.read_mode == 1:
            if self.isubcase not in self.case_control_deck.subcases:
                self.subcase = self.case_control_deck.create_new_subcase(self.isubcase)
            self.subcase.add_op2_data(self.data_code, 'velocity', self.log)
        if self.table_name in [b'OUGV1', b'OUGV2', b'ROUGV1']:
            result_name = 'velocities'
        elif self.table_name == b'OUPV1':
            assert self.thermal in [2, 4], self.thermal
            if self.thermal == 2:
                result_name = 'velocity_scaled_response_spectra_ABS'
            elif self.thermal == 4:
                result_name = 'velocity_scaled_response_spectra_NRL'
            else:
                msg = 'velocities; table_name=%s' % self.table_name
                raise NotImplementedError(msg)
        else:
            msg = 'velocities; table_name=%s' % self.table_name
            raise NotImplementedError(msg)

        #result_name = 'velocities'
        #storage_obj = self.velocities
        storage_obj = getattr(self, result_name)
        if self.thermal == 0:
            result_name = 'velocities'
            storage_obj = self.velocities
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            n = self._read_table_vectorized(data, ndata, result_name, storage_obj,
                                            RealVelocityArray, ComplexVelocityArray,
                                            'node', random_code=self.random_code)
        elif self.thermal == 1:
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            #n = self._read_table(data, ndata, result_name, storage_obj,
                                 #RealThermalVelocityVector, None,
                                 #None, None,
                                 #'node', random_code=self.random_code)
            n = self._read_scalar_table_vectorized(data, ndata, result_name, storage_obj,
                                            RealThermalVelocityVectorArray, None,
                                            'node', random_code=self.random_code)

        elif self.thermal == 2:
            result_name = 'velocity_scaled_response_spectra_ABS'
            storage_obj = self.velocity_scaled_response_spectra_ABS
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            n = self._read_table_vectorized(data, ndata, result_name, storage_obj,
                                            RealVelocityArray, ComplexVelocityArray,
                                            'node', random_code=self.random_code)
            #n = self._not_implemented_or_skip(data, ndata, msg='thermal=2')
        else:
            raise NotImplementedError(self.thermal)
        return n

    def _read_acceleration(self, data, ndata):
        """
        table_code = 11
        """
        if self.read_mode == 1:
            if self.isubcase not in self.case_control_deck.subcases:
                self.subcase = self.case_control_deck.create_new_subcase(self.isubcase)
            self.subcase.add_op2_data(self.data_code, 'acceleration', self.log)

        if self.table_name in [b'OUGV1', b'OUGV2']:
            result_name = 'accelerations'
            assert self.thermal == 0, self.code_information()
        elif self.table_name == b'OUPV1':
            assert self.thermal in [2, 4], self.thermal
            if self.thermal == 2:
                result_name = 'acceleration_scaled_response_spectra_ABS'
            elif self.thermal == 4:
                result_name = 'acceleration_scaled_response_spectra_NRL'
            else:
                msg = 'accelerations; table_name=%s' % self.table_name
                raise NotImplementedError(msg)
        else:
            msg = 'accelerations; table_name=%s' % self.table_name
            raise NotImplementedError(msg)

        if self.thermal == 0:
            result_name = 'accelerations'
            storage_obj = self.accelerations
            if self._results.is_not_saved(result_name):
                return ndata
            n = self._read_table_vectorized(data, ndata, result_name, storage_obj,
                                            RealAccelerationArray,
                                            ComplexAccelerationArray,
                                            'node', random_code=self.random_code)
        elif self.thermal == 1:
            result_name = 'accelerations'
            storage_obj = self.accelerations
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            raise NotImplementedError(self.code_information())
            n = self._read_table(data, ndata, result_name, storage_obj,
                                 None, None,
                                 None, None, 'node', random_code=self.random_code)
        elif self.thermal == 2:
            result_name = 'acceleration_scaled_response_spectra_ABS'
            storage_obj = self.acceleration_scaled_response_spectra_ABS
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            n = self._read_table_vectorized(data, ndata, result_name, storage_obj,
                                            RealAccelerationArray, ComplexAccelerationArray,
                                            'node', random_code=self.random_code)
            #n = self._not_implemented_or_skip(data, ndata, msg='thermal=2')
        elif self.thermal == 4:
            result_name = 'acceleration_scaled_response_spectra_NRL'
            storage_obj = self.acceleration_scaled_response_spectra_NRL
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            n = self._read_table_vectorized(data, ndata, result_name, storage_obj,
                                            RealAccelerationArray, ComplexAccelerationArray,
                                            'node', random_code=self.random_code)
            #n = self._not_implemented_or_skip(data, ndata, msg='thermal=4')
        else:
            raise NotImplementedError(self.thermal)
        return n

    def _read_eigenvector(self, data, ndata):
        """
        table_code = 7
        """
        if self.isubcase not in self.case_control_deck.subcases:
            self.subcase = self.case_control_deck.create_new_subcase(self.isubcase)
        self.subcase.add_op2_data(self.data_code, 'VECTOR', self.log)

        if self.table_name in [b'OUGV1', b'OUGV2', b'BOUGV1', b'BOPHIG']:
            result_name = 'eigenvectors'
        elif self.table_name == b'RADCONS':
            result_name = 'eigenvectors_RADCONS'
        elif self.table_name == b'RADEFFM':
            result_name = 'eigenvectors_RADEFFM'
        elif self.table_name == b'RADEATC':
            result_name = 'eigenvectors_RADEATC'
        else:
            msg = 'eigenvectors; table_name=%s' % self.table_name
            raise NotImplementedError(msg)
        assert self.thermal == 0, self.code_information()

        storage_obj = getattr(self, result_name)
        if self.thermal == 0:
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            n = self._read_table_vectorized(data, ndata, result_name, storage_obj,
                                            RealEigenvectorArray, ComplexEigenvectorArray,
                                            'node', random_code=self.random_code)
        elif self.thermal == 1:
            n = self._not_implemented_or_skip(data, ndata, msg='thermal=1')
            #n = self._read_table(data, result_name, storage_obj,
            #                     None, None,
            #                     None, None, 'node', random_code=self.random_code)
        else:
            raise NotImplementedError(self.thermal)
        return n

    def _read_oug_psd(self, data, ndata):
        """
        table_code = 601/610/611
        """
        if self.thermal == 0:
            if self.table_code in [601]:
                result_name = 'displacementsPSD'
                storage_obj = self.displacementsPSD
                if self._results.is_not_saved(result_name):
                    return ndata
                self._results._found_result(result_name)
                n = self._read_table_vectorized(
                    data, ndata, result_name, storage_obj,
                    RealDisplacementArray, ComplexDisplacementArray,
                    'node', random_code=self.random_code)
            elif self.table_code == 610:
                result_name = 'velocitiesPSD'
                storage_obj = self.velocitiesPSD
                if self._results.is_not_saved(result_name):
                    return ndata
                self._results._found_result(result_name)
                n = self._read_table_vectorized(
                    data, ndata, result_name, storage_obj,
                    RealVelocityArray, ComplexVelocityArray,
                    'node', random_code=self.random_code)
            elif self.table_code == 611:
                result_name = 'accelerationsPSD'
                storage_obj = self.accelerationsPSD
                if self._results.is_not_saved(result_name):
                    return ndata
                n = self._read_table_vectorized(
                    data, ndata, result_name, storage_obj,
                    RealAccelerationArray, ComplexAccelerationArray,
                    'node', random_code=self.random_code)
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
                n = self._not_implemented_or_skip(data, ndata, self.code_information())
                #raise RuntimeError(self.code_information())
        #elif self.thermal == 1:
            #result_name = 'accelerations'
            #storage_obj = self.accelerations
            #if self._results.is_not_saved(result_name):
                #return ndata
            #self._results._found_result(result_name)
            #n = self._read_table(data, ndata, result_name, storage_obj,
                                 #None, None,
                                 #None, None, 'node', random_code=self.random_code)
        #elif self.thermal == 2:
            #result_name = 'acceleration_scaled_response_spectra_ABS'
            #storage_obj = self.acceleration_scaled_response_spectra_ABS
            #if self._results.is_not_saved(result_name):
                #return ndata
            #self._results._found_result(result_name)
            #n = self._read_table(data, ndata, result_name, storage_obj,
                                 #RealAcceleration, ComplexAcceleration,
                                 #RealAccelerationArray, ComplexAccelerationArray,
                                 #'node', random_code=self.random_code)
            ##n = self._not_implemented_or_skip(data, ndata, msg='thermal=2')
        #elif self.thermal == 4:
            #result_name = 'acceleration_scaled_response_spectra_NRL'
            #storage_obj = self.acceleration_scaled_response_spectra_NRL
            #if self._results.is_not_saved(result_name):
                #return ndata
            #self._results._found_result(result_name)
            #n = self._read_table(data, ndata, result_name, storage_obj,
                                 #RealAcceleration, ComplexAcceleration,
                                 #RealAccelerationArray, ComplexAccelerationArray,
                                 #'node', random_code=self.random_code)
            ##n = self._not_implemented_or_skip(data, ndata, msg='thermal=4')
        else:
            raise NotImplementedError(self.thermal)
        return n

    def _read_oug_rms(self, data, ndata):
        """
        table_code = 801  # /610/611
        """
        if self.thermal == 0:
            if self.table_code in [801]:
                result_name = 'displacementsRMS'
                storage_obj = self.displacementsRMS
                if self._results.is_not_saved(result_name):
                    return ndata
                self._results._found_result(result_name)
                n = self._read_table_vectorized(
                    data, ndata, result_name, storage_obj,
                    RealDisplacementArray, ComplexDisplacementArray,
                    'node', random_code=self.random_code)
            #elif self.table_code == 610:
                #result_name = 'velocitiesPSD'
                #storage_obj = self.velocitiesPSD
                #if self._results.is_not_saved(result_name):
                    #return ndata
                #self._results._found_result(result_name)
                #n = self._read_table_vectorized(
                    #data, ndata, result_name, storage_obj,
                    #RealVelocityArray, ComplexVelocityArray,
                    #'node', random_code=self.random_code)
            #elif self.table_code == 611:
                #result_name = 'accelerationsPSD'
                #storage_obj = self.accelerationsPSD
                #if self._results.is_not_saved(result_name):
                    #return ndata
                #n = self._read_table_vectorized(
                    #data, ndata, result_name, storage_obj,
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
                n = self._not_implemented_or_skip(data, ndata, self.code_information())
                #raise RuntimeError(self.code_information())
        #elif self.thermal == 1:
            #result_name = 'accelerations'
            #storage_obj = self.accelerations
            #if self._results.is_not_saved(result_name):
                #return ndata
            #self._results._found_result(result_name)
            #n = self._read_table(data, ndata, result_name, storage_obj,
                                 #None, None,
                                 #None, None, 'node', random_code=self.random_code)
        #elif self.thermal == 2:
            #result_name = 'acceleration_scaled_response_spectra_ABS'
            #storage_obj = self.acceleration_scaled_response_spectra_ABS
            #if self._results.is_not_saved(result_name):
                #return ndata
            #self._results._found_result(result_name)
            #n = self._read_table(data, ndata, result_name, storage_obj,
                                 #RealAcceleration, ComplexAcceleration,
                                 #RealAccelerationArray, ComplexAccelerationArray,
                                 #'node', random_code=self.random_code)
            ##n = self._not_implemented_or_skip(data, ndata, msg='thermal=2')
        #elif self.thermal == 4:
            #result_name = 'acceleration_scaled_response_spectra_NRL'
            #storage_obj = self.acceleration_scaled_response_spectra_NRL
            #if self._results.is_not_saved(result_name):
                #return ndata
            #self._results._found_result(result_name)
            #n = self._read_table(data, ndata, result_name, storage_obj,
                                 #RealAcceleration, ComplexAcceleration,
                                 #RealAccelerationArray, ComplexAccelerationArray,
                                 #'node', random_code=self.random_code)
            ##n = self._not_implemented_or_skip(data, ndata, msg='thermal=4')
        else:
            raise NotImplementedError(self.thermal)
        return n

    def _read_oug_no(self, data, ndata):
        """
        table_code = 901  # /610/611
        """
        if self.thermal == 0:
            if self.table_code in [901]:
                result_name = 'displacementsNO'
                storage_obj = self.displacementsNO
                if self._results.is_not_saved(result_name):
                    return ndata
                self._results._found_result(result_name)
                n = self._read_table_vectorized(
                    data, ndata, result_name, storage_obj,
                    RealDisplacementArray, ComplexDisplacementArray,
                    'node', random_code=self.random_code)
            #elif self.table_code == 610:
                #result_name = 'velocitiesPSD'
                #storage_obj = self.velocitiesPSD
                #if self._results.is_not_saved(result_name):
                    #return ndata
                #self._results._found_result(result_name)
                #n = self._read_table_vectorized(
                    #data, ndata, result_name, storage_obj,
                    #RealVelocityArray, ComplexVelocityArray,
                    #'node', random_code=self.random_code)
            #elif self.table_code == 611:
                #result_name = 'accelerationsPSD'
                #storage_obj = self.accelerationsPSD
                #if self._results.is_not_saved(result_name):
                    #return ndata
                #n = self._read_table_vectorized(
                    #data, ndata, result_name, storage_obj,
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
                n = self._not_implemented_or_skip(data, ndata, self.code_information())
                #raise RuntimeError(self.code_information())
        #elif self.thermal == 1:
            #result_name = 'accelerations'
            #storage_obj = self.accelerations
            #if self._results.is_not_saved(result_name):
                #return ndata
            #self._results._found_result(result_name)
            #n = self._read_table(data, ndata, result_name, storage_obj,
                                 #None, None,
                                 #None, None, 'node', random_code=self.random_code)
        #elif self.thermal == 2:
            #result_name = 'acceleration_scaled_response_spectra_ABS'
            #storage_obj = self.acceleration_scaled_response_spectra_ABS
            #if self._results.is_not_saved(result_name):
                #return ndata
            #self._results._found_result(result_name)
            #n = self._read_table(data, ndata, result_name, storage_obj,
                                 #RealAcceleration, ComplexAcceleration,
                                 #RealAccelerationArray, ComplexAccelerationArray,
                                 #'node', random_code=self.random_code)
            ##n = self._not_implemented_or_skip(data, ndata, msg='thermal=2')
        #elif self.thermal == 4:
            #result_name = 'acceleration_scaled_response_spectra_NRL'
            #storage_obj = self.acceleration_scaled_response_spectra_NRL
            #if self._results.is_not_saved(result_name):
                #return ndata
            #self._results._found_result(result_name)
            #n = self._read_table(data, ndata, result_name, storage_obj,
                                 #RealAcceleration, ComplexAcceleration,
                                 #RealAccelerationArray, ComplexAccelerationArray,
                                 #'node', random_code=self.random_code)
            ##n = self._not_implemented_or_skip(data, ndata, msg='thermal=4')
        else:
            raise NotImplementedError(self.thermal)
        return n

    def _read_oug_ato(self, data, ndata):
        """
        table_code = 901  # /610/611
        """
        if self.thermal == 0:
            if self.table_code in []:
                result_name = 'displacementsATO'
                storage_obj = self.displacementsATO
                if self._results.is_not_saved(result_name):
                    return ndata
                self._results._found_result(result_name)
                n = self._read_table_vectorized(
                    data, ndata, result_name, storage_obj,
                    RealDisplacementArray, ComplexDisplacementArray,
                    'node', random_code=self.random_code)
            #elif self.table_code == 610:
                #result_name = 'velocitiesPSD'
                #storage_obj = self.velocitiesPSD
                #if self._results.is_not_saved(result_name):
                    #return ndata
                #self._results._found_result(result_name)
                #n = self._read_table_vectorized(
                    #data, ndata, result_name, storage_obj,
                    #RealVelocityArray, ComplexVelocityArray,
                    #'node', random_code=self.random_code)
            #elif self.table_code == 611:
                #result_name = 'accelerationsPSD'
                #storage_obj = self.accelerationsPSD
                #if self._results.is_not_saved(result_name):
                    #return ndata
                #n = self._read_table_vectorized(
                    #data, ndata, result_name, storage_obj,
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
                n = self._not_implemented_or_skip(data, ndata, self.code_information())
                #raise RuntimeError(self.code_information())
        #elif self.thermal == 1:
            #result_name = 'accelerations'
            #storage_obj = self.accelerations
            #if self._results.is_not_saved(result_name):
                #return ndata
            #self._results._found_result(result_name)
            #n = self._read_table(data, ndata, result_name, storage_obj,
                                 #None, None,
                                 #None, None, 'node', random_code=self.random_code)
        #elif self.thermal == 2:
            #result_name = 'acceleration_scaled_response_spectra_ABS'
            #storage_obj = self.acceleration_scaled_response_spectra_ABS
            #if self._results.is_not_saved(result_name):
                #return ndata
            #self._results._found_result(result_name)
            #n = self._read_table(data, ndata, result_name, storage_obj,
                                 #RealAcceleration, ComplexAcceleration,
                                 #RealAccelerationArray, ComplexAccelerationArray,
                                 #'node', random_code=self.random_code)
            ##n = self._not_implemented_or_skip(data, ndata, msg='thermal=2')
        #elif self.thermal == 4:
            #result_name = 'acceleration_scaled_response_spectra_NRL'
            #storage_obj = self.acceleration_scaled_response_spectra_NRL
            #if self._results.is_not_saved(result_name):
                #return ndata
            #self._results._found_result(result_name)
            #n = self._read_table(data, ndata, result_name, storage_obj,
                                 #RealAcceleration, ComplexAcceleration,
                                 #RealAccelerationArray, ComplexAccelerationArray,
                                 #'node', random_code=self.random_code)
            ##n = self._not_implemented_or_skip(data, ndata, msg='thermal=4')
        else:
            raise NotImplementedError(self.thermal)
        return n
