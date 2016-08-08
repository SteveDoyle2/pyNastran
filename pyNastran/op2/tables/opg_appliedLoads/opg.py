#pylint: disable=C0301,C0103
"""
Defines the Real/Complex Forces created by:
    OLOAD = ALL
"""

#from pyNastran.op2.tables.opg_appliedLoads.opg_objects import (#RealAppliedLoads,  #ComplexAppliedLoads,
                                                               #RealAppliedLoadsVectorArray, ComplexAppliedLoadsVectorArray)
from pyNastran.op2.tables.opg_appliedLoads.opg_load_vector import (
    #RealLoadVector, ComplexLoadVector,
    RealLoadVectorArray, ComplexLoadVectorArray,
    #RealThermalLoadVector,
    RealTemperatureVectorArray,
    #RealThermalVelocityArray
)
from pyNastran.op2.tables.opg_appliedLoads.opnl_force_vector import RealForceVectorArray#, ComplexForceVectorArray

from pyNastran.op2.op2_common import OP2Common


class OPG(OP2Common):
    def __init__(self):
        pass

    def _read_opg1_3(self, data, ndata):
        self.words = [
            'aCode', 'tCode', '???', 'isubcase',
            '???', '???', '???', 'dLoadID',
            'format_code', 'num_wide', 'o_code', '???',
            'acoustic_flag', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', 'thermal', '???',
            '???', 'Title', 'subtitle', 'label']

        self.parse_approach_code(data)
        #isubcase = self.get_values(data,'i',4)

        ## dynamic load set ID/random code
        self.random_code = self.add_data_parameter(data, 'random_code', 'i', 8, False)

        ## format code
        self.format_code = self.add_data_parameter(data, 'format_code', 'i', 9, False)

        ## number of words per entry in record
        ## .. note:: is this needed for this table ???
        self.num_wide = self.add_data_parameter(data, 'num_wide', 'i', 10, False)

        ## undefined in DMAP...
        self.oCode = self.add_data_parameter(data, 'oCode', 'i', 11, False)

        ## thermal flag; 1 for heat transfer, 0 otherwise
        self.thermal = self.add_data_parameter(data, 'thermal', 'i', 23, False)

        #print "dLoadID(8)=%s format_code(9)=%s num_wide(10)=%s oCode(11)=%s thermal(23)=%s" %(self.dLoadID,self.format_code,self.num_wide,self.oCode,self.thermal)
        if not self.is_sort1():
            raise NotImplementedError('sort2...')
        #assert self.isThermal()==False,self.thermal

        ## assuming tCode=1
        if self.analysis_code == 1:   # statics
            ## load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5, False)
            self.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysis_code == 2:  # normal modes/buckling (real eigenvalues)
            ## mode number
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            ## mode or cycle .. todo:: confused on the type - F1???
            self.mode2 = self.add_data_parameter(data, 'mode2', 'i', 7, False)
            self.cycle = self.add_data_parameter(data, 'cycle', 'f', 7, False)
            self.update_mode_cycle('cycle')
            self.data_names = self.apply_data_code_value('data_names', ['mode', 'eigr', 'mode2', 'cycle', ])
        #elif self.analysis_code == 3: # differential stiffness
        #    ## load set number
        #    self.lsdvmn = self.get_values(data,'i',5)
        #elif self.analysis_code == 4: # differential stiffness
        #    ## load set number
        #    self.lsdvmn = self.get_values(data,'i',5)
        elif self.analysis_code == 5:   # frequency
            ## frequency
            self.freq = self.add_data_parameter(data, 'freq', 'f', 5)
            self.data_names = self.apply_data_code_value('data_names', ['freq'])
        elif self.analysis_code == 6:  # transient
            ## time step
            self.time = self.add_data_parameter(data, 'time', 'f', 5)
            self.data_names = self.apply_data_code_value('data_names', ['time'])
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
            raise RuntimeError('invalid analysis_code...analysis_code=%s' %
                               self.analysis_code)

        # tCode=2
        #if self.analysis_code==2: # sort2
        #    self.lsdvmn = self.get_values(data,'i',5) ## load set, Mode number

        #print "*isubcase=%s"%(self.isubcase)
        #print "analysis_code=%s table_code=%s thermal=%s" %(self.analysis_code,self.table_code,self.thermal)

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


    def _read_opg1_4(self, data, ndata):
        if self.table_code == 2:  # load vector
            assert self.table_name in [b'OPG1', b'OPGV1'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
            n = self._read_load_vector(data, ndata)
        elif self.table_code == 12:  # ???
            n = self._read_force_vector(data, ndata)
        #else:
            #n = self._not_implemented_or_skip('bad OPG table')
        else:
            raise NotImplementedError(self.table_code)
        return n

    def _read_load_vector(self, data, ndata):
        """
        table_code = 2
        """
        if self.thermal == 0:
            result_name = 'load_vectors'
            storage_obj = self.load_vectors
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            n = self._read_table_vectorized(data, ndata, result_name, storage_obj,
                                            RealLoadVectorArray, ComplexLoadVectorArray,
                                            'node', random_code=self.random_code)
        elif self.thermal == 1:
            result_name = 'thermal_load_vectors'
            storage_obj = self.thermal_load_vectors

            #RealThermalLoadVectorVector = None
            #ComplexThermalLoadVectorVector = None
            ComplexThermalLoadVectorArray = None
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            #n = self._read_table(data, ndata, result_name, storage_obj,
                                 #RealThermalLoadVector, ComplexThermalLoadVectorVector,
                                 #RealTemperatureVectorArray, ComplexThermalLoadVectorArray,
                                 #'node', random_code=self.random_code)
            n = self._read_scalar_table_vectorized(data, ndata, result_name, storage_obj,
                                                   RealTemperatureVectorArray, ComplexThermalLoadVectorArray,
                                                   'node', random_code=self.random_code)

        else:
            raise NotImplementedError(self.thermal)
        return n

    def _read_force_vector(self, data, ndata):
        """
        table_code = 12
        """
        if self.thermal == 0:
            result_name = 'force_vectors'
            storage_obj = self.force_vectors
            #ForceVectorVector = None
            ComplexForceVectorArray = None
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            #n = self._read_table(data, ndata, result_name, storage_obj,
                                 #RealForceVector, ComplexForceVector,
                                 #RealForceVectorArray, ComplexForceVectorArray,
                                 #'node', random_code=self.random_code)
            n = self._read_table_vectorized(data, ndata, result_name, storage_obj,
                                            RealForceVectorArray, ComplexForceVectorArray,
                                            'node', random_code=self.random_code)

        #elif self.thermal == 1:
            #result_name = 'thermal_force_vectors'
            #storage_obj = self.thermal_force_vectors
            #RealThermalForceVector = None
            #ComplexThermalForceVector = None
            #RealThermalForceVectorVector = None
            #ComplexThermalForceVectorArray = None
            #n = self._read_table(data, result_name, storage_obj,
                                 #RealThermalForceVector, ComplexThermalForceVector,
                                 #RealThermalForceVectorArray, ComplexThermalForceVectorArray,
                                 #'node', random_code=self.random_code)
        else:
            raise NotImplementedError(self.thermal)
        return n
