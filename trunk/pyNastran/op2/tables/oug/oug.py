#pylint: disable=C0301,C0103
"""
This file defines the OUG Table, which contains:
 * Real/Complex Displacement
 * Real/Complex Acceleration
 * Real/Complex Velocity
 * Real/Complex Eigenvectors
 * Real Temperature
"""
from struct import unpack

from pyNastran.op2.op2_common import OP2Common

from pyNastran.op2.tables.oug.oug_displacements import (
    RealDisplacementArray,
    ComplexDisplacementArray,
    RealDisplacement,          # table_code=1     format_code=1 sort_code=0
    ComplexDisplacement)       # analysis_code=5  format_code=3 sort_code=1

# table_code=10 format_code=1 sort_code=0
from pyNastran.op2.tables.oug.oug_velocities import (
    RealVelocityArray, ComplexVelocityArray,
    RealVelocity, ComplexVelocity)

# table_code=11 format_code=1 sort_code=0
from pyNastran.op2.tables.oug.oug_accelerations import (
    RealAccelerationArray, ComplexAccelerationArray,
    RealAcceleration, ComplexAcceleration)

# table_code=1 format_code=1 sort_code=0
from pyNastran.op2.tables.oug.oug_temperatures import (
    RealTemperatureArray, RealTemperature)

from pyNastran.op2.tables.oug.oug_eigenvectors import (
    RealEigenvectorArray,
    ComplexEigenvectorArray,

    #ComplexEigenvector              # analysis_code=5, sort_code=1 format_code=1 table_code=7
     Eigenvector,                    # analysis_code=2, sort_code=0 format_code   table_code=7
     ComplexEigenvector,             # analysis_code=5, sort_code=1 format_code=1 table_code=7
    #RealEigenvector,                # analysis_code=9, sort_code=1 format_code=1 table_code=7
)

from pyNastran.op2.tables.opg_appliedLoads.opg_loadVector import RealThermalVelocityVector


class OUG(OP2Common):

    def __init__(self):
        OP2Common.__init__(self)

    def _read_oug1_3(self, data):
        three = self.parse_approach_code(data)
        self.words = [
            'approach_code','table_code','???',         'isubcase',
             '???',         '???',      '???',          'random_code'
             'format_code', 'num_wide', '???',          '???',
            'acoustic_flag','???',      '???',          '???',
             '???',         '???',      '???',          '???',
             '???',         '???',      'thermal',      '???',
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
            self.dataNames = self.apply_data_code_value('dataNames', ['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysis_code == 2:  # real eigenvalues
            # mode number
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)
            # real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            self.mode_cycle = self.add_data_parameter(data, 'mode_cycle', 'i', 7, False)  # mode or cycle .. todo:: confused on the type - F1???
            self.dataNames = self.apply_data_code_value('dataNames', ['mode', 'eigr', 'mode_cycle'])
        #elif self.analysis_code == 3: # differential stiffness
            #self.lsdvmn = self.get_values(data, 'i', 5) ## load set number
            #self.data_code['lsdvmn'] = self.lsdvmn
        #elif self.analysis_code == 4: # differential stiffness
            #self.lsdvmn = self.get_values(data, 'i', 5) ## load set number
        elif self.analysis_code == 5:   # frequency
            # frequency
            self.freq = self.add_data_parameter(data, 'freq', 'f', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['freq'])
        elif self.analysis_code == 6:  # transient
            # time step
            self.dt = self.add_data_parameter(data, 'dt', 'f', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['dt'])
        elif self.analysis_code == 7:  # pre-buckling
            # load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['lsdvmn'])
        elif self.analysis_code == 8:  # post-buckling
            # load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            # real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            self.dataNames = self.apply_data_code_value('dataNames', ['lsdvmn', 'eigr'])
        elif self.analysis_code == 9:  # complex eigenvalues
            # mode number
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)
            # real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            # imaginary eigenvalue
            self.eigi = self.add_data_parameter(data, 'eigi', 'f', 7, False)
            self.dataNames = self.apply_data_code_value('dataNames', ['mode', 'eigr', 'eigi'])
        elif self.analysis_code == 10:  # nonlinear statics
            # load step
            self.lftsfq = self.add_data_parameter(data, 'lftsfq', 'f', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['lftsfq'])
        elif self.analysis_code == 11:  # old geometric nonlinear statics
            # load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['lsdvmn'])
        elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            # load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['lsdvmn'])
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % self.analysis_code
            raise RuntimeError(msg)

        #print self.code_information()
        if self.debug:
            self.binary_debug.write('  approach_code = %r\n' % self.approach_code)
            self.binary_debug.write('  tCode    = %r\n' % self.tCode)
            self.binary_debug.write('  isubcase = %r\n' % self.isubcase)
        self._read_title(data)
        self._write_debug_bits()

    def _read_oug1_4(self, data):
        if self.is_sort2():
            raise NotImplementedError('SORT2!!!!!')

        if self.table_code == 1:   # Displacements
            assert self.table_name in ['OUG1', 'BOUGV1', 'OUGV1', 'OUPV1'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
            n = self._read_displacement(data)
        elif self.table_code == 7:
            n = self._read_eigenvector(data)
        elif self.table_code == 10:
            n = self._read_velocity(data)
        elif self.table_code == 11:
            n = self._read_acceleration(data)
        else:
            raise NotImplementedError(self.table_code)
        #else:
            #self._not_implemented_or_skip(data, 'bad OUG table')
        return n

    def _read_eigenvector_displacement_solution_set(self, data):
        """
        table_code = 14
        """
        raise NotImplementedError()

    def _read_displacement_solution_set(self, data):
        """
        table_code = 15
        """
        raise NotImplementedError()

    def _read_velocity_solution_set(self, data):
        """
        table_code = 16
        """
        raise NotImplementedError()

    def _read_acceleration_solution_set(self, data):
        """
        table_code = 17
        """
        raise NotImplementedError()

    def _read_displacement(self, data):
        if self.thermal == 0:
            result_name = 'displacements'
            storage_obj = self.displacements
            if result_name not in self._saved_results:
                return len(data)
            n = self._read_table(data, result_name, storage_obj,
                                 RealDisplacement, ComplexDisplacement,
                                 RealDisplacementArray, ComplexDisplacementArray, 'node', random_code=self.random_code)
        elif self.thermal == 1:
            result_name = 'temperatures'
            storage_obj = self.temperatures
            if result_name not in self._saved_results:
                return len(data)
            n = self._read_table(data, result_name, storage_obj,
                                 RealTemperature, None,
                                 RealTemperatureArray, None, 'node', random_code=self.random_code)
        elif self.thermal == 2:
            result_name = 'displacement_scaled_response_spectra_ABS'
            storage_obj = self.displacement_scaled_response_spectra_ABS
            if result_name not in self._saved_results:
                return len(data)
            n = self._read_table(data, result_name, storage_obj,
                                 RealDisplacement, ComplexDisplacement,
                                 RealDisplacementArray, ComplexDisplacementArray, 'node', random_code=self.random_code)
        elif self.thermal == 4:
            result_name = 'displacement_scaled_response_spectra_NRL'
            storage_obj = self.displacement_scaled_response_spectra_NRL
            if result_name not in self._saved_results:
                return len(data)
            n = self._read_table(data, result_name, storage_obj,
                                 RealDisplacement, ComplexDisplacement,
                                 RealDisplacementArray, ComplexDisplacementArray, 'node', random_code=self.random_code)
            #return self._not_implemented_or_skip(data, msg='thermal=4')
        else:
            n = self._not_implemented_or_skip(data, 'bad thermal=%r table' % self.thermal)
        #else:
            #raise NotImplementedError(self.thermal)
        return n

    def _read_velocity(self, data):
        """
        table_code = 10
        """
        result_name = 'velocities'
        storage_obj = self.velocities
        if self.thermal == 0:
            real_obj = RealVelocity
            complex_obj = ComplexVelocity
            if result_name not in self._saved_results:
                return len(data)
            n = self._read_table(data, result_name, storage_obj,
                                 RealVelocity, ComplexVelocity,
                                 RealVelocityArray, ComplexVelocityArray,
                                 'node', random_code=self.random_code)
        elif self.thermal == 1:
            real_obj = RealThermalVelocityVector
            complex_obj = None
            if result_name not in self._saved_results:
                return len(data)
            n = self._read_table(data, result_name, storage_obj,
                                 RealThermalVelocityVector, None,
                                 None, None, 'node', random_code=self.random_code)
        elif self.thermal == 2:
            result_name = 'velocity_scaled_response_spectra_ABS'
            storage_obj = self.velocity_scaled_response_spectra_ABS
            if result_name not in self._saved_results:
                return len(data)
            n = self._read_table(data, result_name, storage_obj,
                                 RealVelocity, ComplexVelocity,
                                 RealVelocityArray, ComplexVelocityArray,
                                 'node', random_code=self.random_code)
            #n = self._not_implemented_or_skip(data, msg='thermal=2')
        else:
            raise NotImplementedError(self.thermal)
        return n

    def _read_acceleration(self, data):
        """
        table_code = 11
        """
        if self.thermal == 0:
            result_name = 'accelerations'
            storage_obj = self.accelerations
            if result_name not in self._saved_results:
                return len(data)
            n = self._read_table(data, result_name, storage_obj,
                                 RealAcceleration, ComplexAcceleration,
                                 RealAccelerationArray, ComplexAccelerationArray,
                                 'node', random_code=self.random_code)
        elif self.thermal == 1:
            result_name = 'accelerations'
            storage_obj = self.accelerations
            if result_name not in self._saved_results:
                return len(data)
            n = self._read_table(data, result_name, storage_obj,
                                 None, None,
                                 None, None, 'node', random_code=self.random_code)
        elif self.thermal == 2:
            result_name = 'acceleration_scaled_response_spectra_ABS'
            storage_obj = self.acceleration_scaled_response_spectra_ABS
            if result_name not in self._saved_results:
                return len(data)
            n = self._read_table(data, result_name, storage_obj,
                                 RealAcceleration, ComplexAcceleration,
                                 RealAccelerationArray, ComplexAccelerationArray,
                                 'node', random_code=self.random_code)
            #n = self._not_implemented_or_skip(data, msg='thermal=2')
        elif self.thermal == 4:
            result_name = 'acceleration_scaled_response_spectra_NRL'
            storage_obj = self.acceleration_scaled_response_spectra_NRL
            if result_name not in self._saved_results:
                return len(data)
            n = self._read_table(data, result_name, storage_obj,
                                 RealAcceleration, ComplexAcceleration,
                                 RealAccelerationArray, ComplexAccelerationArray,
                                 'node', random_code=self.random_code)
            #n = self._not_implemented_or_skip(data, msg='thermal=4')
        else:
            raise NotImplementedError(self.thermal)
        return n

    def _read_eigenvector(self, data):
        """
        table_code = 7
        """
        result_name = 'eigenvectors'
        storage_obj = self.eigenvectors
        if self.thermal == 0:
            if result_name not in self._saved_results:
                return len(data)
            n = self._read_table(data, result_name, storage_obj,
                                 Eigenvector, ComplexEigenvector,
                                 RealEigenvectorArray, ComplexEigenvectorArray, 'node')
        elif self.thermal == 1:
            n = self._not_implemented_or_skip(data, msg='thermal=1')
            #n = self._read_table(data, result_name, storage_obj,
            #                     None, None,
            #                     None, None, 'node', random_code=self.random_code)
        else:
            raise NotImplementedError(self.thermal)
        return n

    def _read_eigenvector_backup(self, data):
        """
        table_code = 7
        """
        result_name = 'eigenvectors'
        storage_obj = self.eigenvectors
        if self.thermal == 0:
            if result_name not in self._saved_results:
                return len(data)
            #if self.isRandomResponse():
            if self.isRandomResponse():
                if self.format_code == 1 and self.num_wide == 8:  # real/random
                    real_obj = Eigenvector
                    #assert real_obj is not None
                    nnodes = len(data) // 32  # 8*4
                    auto_return = self._create_table_object(result_name, nnodes, storage_obj, real_obj, real_vector)
                    if auto_return:
                        return len(data)
                    n = self._read_real_table(data, result_name, node_elem)
                elif self.format_code in [1, 2, 3] and self.num_wide == 14:  # real (fsi.op2 odd...) or real/imaginary or mag/phase
                    complex_obj = ComplexEigenvector
                    #assert complex_obj is not None
                    nnodes = len(data) // 56  # 14*4
                    auto_return = self._create_table_object(result_name, nnodes, storage_obj, complex_obj, complex_vector)
                    if auto_return:
                        return len(data)
                    n = self._read_complex_table(data, result_name, node_elem)
                else:
                    msg = 'only num_wide=8 or 14 is allowed  num_wide=%s' % self.num_wide
                    n = self._not_implemented_or_skip(data, msg)
            else:
                msg = 'invalid random_code=%s num_wide=%s' % (self.random_code, self.num_wide)
                n = self._not_implemented_or_skip(data, msg)
            #n = self._read_table(data, result_name, storage_obj,
            #                     Eigenvector, ComplexEigenvector,
            #                     RealEigenvectorArray, ComplexEigenvectorArray, 'node')
        elif self.thermal == 1:
            n = self._not_implemented_or_skip(data, msg='thermal=1')
            #n = self._read_table(data, result_name, storage_obj,
            #                     None, None,
            #                     None, None, 'node', random_code=self.random_code)
        else:
            raise NotImplementedError(self.thermal)
        return n
