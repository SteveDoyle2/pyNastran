from struct import unpack

from pyNastran.op2.dev.op2_common import OP2Common

from pyNastran.op2.tables.oqg_constraintForces.oqg_spcForces import (
    SPCForcesObject, ComplexSPCForcesObject)
from pyNastran.op2.tables.oqg_constraintForces.oqg_mpcForces import (
    MPCForcesObject, ComplexMPCForcesObject)
from pyNastran.op2.tables.oqg_constraintForces.oqg_thermalGradientAndFlux import (
    TemperatureGradientAndFluxObject)


class OQG(OP2Common):
    def __init__(self):
        OP2Common.__init__(self)

    def _read_oqg1_3(self, data):
        three = self.parse_approach_code(data)
        self.words = [
            'aCode',       'tCode',    '???',           'isubcase',
             '???',         '???',      '???',          'random_code'
             'format_code', 'num_wide', '11',           '???',
            'acoustic_flag','???',      '???',          '???',
             '???',         '???',      '???',          '???',
             '???',         '???',      'thermal',      '???',
             '???', 'Title', 'subtitle', 'label']

        ## random code
        self.random_code = self.add_data_parameter( data, 'random_code', 'i', 8, False)

        ## format code
        self.format_code = self.add_data_parameter( data, 'format_code', 'i', 9, False)

        ## number of words per entry in record
        self.num_wide = self.add_data_parameter(data, 'num_wide', 'i', 10, False)

        ## acoustic pressure flag
        self.acoustic_flag = self.add_data_parameter(data, 'acoustic_flag', 'f', 13, False)

        ## thermal flag; 1 for heat ransfer, 0 otherwise
        self.thermal = self.add_data_parameter(data, 'thermal', 'i', 23, False)

        if not self.is_sort1():
            raise NotImplementedError('sort2...')
        #assert self.isThermal()==False,self.thermal

        #self.print_block(data) # on
        ## assuming tCode=1
        if self.analysis_code == 1:   # statics / displacement / heat flux
            ## load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5, False)
            self.dataNames = self.apply_data_code_value('dataNames', ['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysis_code == 2:  # real eigenvalues
            ## mode number
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            ## mode or cycle .. todo:: confused on the type - F1???
            self.mode_cycle = self.add_data_parameter(data, 'mode_cycle', 'f', 7, False)
            self.dataNames = self.apply_data_code_value('dataNames', ['mode', 'eigr', 'mode_cycle'])
        #elif self.analysis_code==3: # differential stiffness
            #self.lsdvmn = self.get_values(data,'i',5) ## load set number
            #self.dataNames = self.data_code['lsdvmn'] = self.lsdvmn
        #elif self.analysis_code==4: # differential stiffness
            #self.lsdvmn = self.get_values(data,'i',5) ## load set number
        elif self.analysis_code == 5:   # frequency
            ## frequency
            self.freq = self.add_data_parameter(data, 'freq', 'f', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['freq'])
        elif self.analysis_code == 6:  # transient
            ## time step
            self.dt = self.add_data_parameter(data, 'dt', 'f', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['dt'])
        elif self.analysis_code == 7:  # pre-buckling
            ## load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['lsdvmn'])
        elif self.analysis_code == 8:  # post-buckling
            ## load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            ## real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            self.dataNames = self.apply_data_code_value('dataNames', ['lsdvmn', 'eigr'])
        elif self.analysis_code == 9:  # complex eigenvalues
            ## mode number
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            ## imaginary eigenvalue
            self.eigi = self.add_data_parameter(data, 'eigi', 'f', 7, False)
            self.dataNames = self.apply_data_code_value('dataNames', ['mode', 'eigr', 'eigi'])
        elif self.analysis_code == 10:  # nonlinear statics
            ## load step
            self.lftsfq = self.add_data_parameter(data, 'lftsfq', 'f', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['lftsfq'])
        elif self.analysis_code == 11:  # old geometric nonlinear statics
            ## load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['lsdvmn'])
        elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['lsdvmn'])
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % (self.analysis_code)
            raise RuntimeError(msg)

        #print self.code_information()
        if self.debug:
            self.binary_debug.write('  aCode    = %r\n' % self.aCode)
            self.binary_debug.write('  tCode    = %r\n' % self.tCode)
            self.binary_debug.write('  isubcase = %r\n' % self.isubcase)
        self._read_title(data)
        self._write_debug_bits()

    def _read_oqg1_4(self, data):
        if self.table_code == 3:   # SPC Forces
            assert self.table_name in ['OQG1', 'OQGV1', 'OQP1'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
            self._read_spc_forces(data)
        elif self.table_code == 39:  # MPC Forces
            assert self.table_name in ['OQMG1'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
            self._read_mpc_forces(data)
        else:
            raise NotImplementedError(self.table_code)
        #else:
            #self.not_implemented_or_skip('bad OQG table')

    def _read_spc_forces(self, data):
        """
        table_code = 3
        """
        if self.thermal == 0:
            real_obj = SPCForcesObject
            complex_obj = ComplexSPCForcesObject
            storage_obj = self.spcForces
            self._read_table(data, storage_obj, real_obj, complex_obj, 'node')
        elif self.thermal == 1:
            result_name = 'thermalGradientAndFlux' #'finite element temperature gradients and fluxes'
            storage_obj =  self.thermalGradientAndFlux
            real_obj = TemperatureGradientAndFluxObject
            complex_obj = None
            self._read_table(data, storage_obj, real_obj, complex_obj, 'node')
        else:
            raise NotImplementedError(self.thermal)

    def _read_mpc_forces(self, data):
        """
        table_code = 39
        """
        result_name = 'mpcForces'
        storage_obj = self.mpcForces
        if self.thermal == 0:
            real_obj = MPCForcesObject
            complex_obj = ComplexMPCForcesObject
            self._read_table(data, storage_obj, real_obj, complex_obj, 'node')
        elif self.thermal == 1:
            real_obj = None
            complex_obj = None
            self._read_table(data, storage_obj, real_obj, complex_obj, 'node')
        else:
            raise NotImplementedError(self.thermal)
