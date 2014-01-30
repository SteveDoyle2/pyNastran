#pylint: disable=C0301,C0103
from struct import unpack

from pyNastran.op2.tables.opg_appliedLoads.opg_Objects import AppliedLoadsObject  # ComplexAppliedLoadsObject
from pyNastran.op2.tables.opg_appliedLoads.opg_loadVector import LoadVectorObject, ComplexLoadVectorObject, ThermalLoadVectorObject
from pyNastran.op2.tables.opg_appliedLoads.opnl_forceVector import ForceVectorObject, ComplexForceVectorObject

class OPG(object):
    def __init__(self):
        pass


    def read_opg1_3(self, data):
        three = self.parse_approach_code(data)
        self.words = [
            'aCode',       'tCode',    '???',           'isubcase',
             '???',         '???',      '???',          'dLoadID'
             'format_code', 'num_wide', 'o_code',       '???',
            'acoustic_flag','???',      '???',          '???',
             '???',         '???',      '???',          '???',
             '???',         '???',      'thermal',      '???',
             '???', 'Title', 'subtitle', 'label']

        self.parse_approach_code(data)
        #isubcase = self.get_values(data,'i',4)

        ## dynamic load set ID/random code
        self.dLoadID = self.add_data_parameter(data, 'dLoadID', 'i', 8, False)

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
            self.add_data_parameter(data, 'lsdvmn', 'i', 5, False)
            self.apply_data_code_value('dataNames', ['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysis_code == 2:  # normal modes/buckling (real eigenvalues)
            ## mode number
            self.add_data_parameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            self.add_data_parameter(data, 'eign', 'f', 6, False)
            ## mode or cycle .. todo:: confused on the type - F1???
            self.add_data_parameter(data, 'mode_cycle', 'f', 7, False)
            self.apply_data_code_value('dataNames', ['mode', 'eign', 'mode_cycle'])
        #elif self.analysis_code == 3: # differential stiffness
        #    ## load set number
        #    self.lsdvmn = self.get_values(data,'i',5)
        #elif self.analysis_code == 4: # differential stiffness
        #    ## load set number
        #    self.lsdvmn = self.get_values(data,'i',5)
        elif self.analysis_code == 5:   # frequency
            ## frequency
            self.add_data_parameter(data, 'freq', 'f', 5)
            self.apply_data_code_value('dataNames', ['freq'])
        elif self.analysis_code == 6:  # transient
            ## time step
            self.add_data_parameter(data, 'time', 'f', 5)
            self.apply_data_code_value('dataNames', ['time'])
        elif self.analysis_code == 7:  # pre-buckling
            ## load set number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.apply_data_code_value('dataNames', ['lsdvmn'])
        elif self.analysis_code == 8:  # post-buckling
            ## load set number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            ## real eigenvalue
            self.add_data_parameter(data, 'eigr', 'f', 6, False)
            self.apply_data_code_value('dataNames', ['lsdvmn', 'eigr'])
        elif self.analysis_code == 9:  # complex eigenvalues
            ## mode number
            self.add_data_parameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            self.add_data_parameter(data, 'eigr', 'f', 6, False)
            ## imaginary eigenvalue
            self.add_data_parameter(data, 'eigi', 'f', 7, False)
            self.apply_data_code_value('dataNames', ['mode', 'eigr', 'eigi'])
        elif self.analysis_code == 10:  # nonlinear statics
            ## load step
            self.add_data_parameter(data, 'lftsfq', 'f', 5)
            self.apply_data_code_value('dataNames', ['lftsfq'])
        elif self.analysis_code == 11:  # old geometric nonlinear statics
            ## load set number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.apply_data_code_value('dataNames', ['lsdvmn'])
        elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## load set number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.apply_data_code_value('dataNames', ['lsdvmn'])
        else:
            raise RuntimeError('invalid analysis_code...analysis_code=%s' %
                               self.analysis_code)

        # tCode=2
        #if self.analysis_code==2: # sort2
        #    self.lsdvmn = self.get_values(data,'i',5) ## load set, Mode number

        #print "*isubcase=%s"%(self.isubcase)
        #print "analysis_code=%s table_code=%s thermal=%s" %(self.analysis_code,self.table_code,self.thermal)

        #print self.code_information()
        if self.debug:
            self.binary_debug.write('  aCode    = %r\n' % self.aCode)
            self.binary_debug.write('  tCode    = %r\n' % self.tCode)
            self.binary_debug.write('  isubcase = %r\n' % self.isubcase)
        self._read_title(data)
        self._write_debug_bits()


    def read_opg1_4(self, data):
        if self.table_code == 2:  # load vector
            assert self.table_name in ['OPG1', 'OPGV1'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
            self.read_load_vector(data)
        elif self.table_code == 12:  # ???
            self.read_force_vector(data)
            pass
        #else:
            #self.not_implemented_or_skip('bad OPG table')
        else:
            raise NotImplementedError(self.table_code)

    def read_load_vector(self, data):
        """
        table_code = 2
        """
        if self.thermal == 0:
            result_name = 'loadVectors'
            storage_obj = self.loadVectors
            real_obj = LoadVectorObject
            complex_obj = ComplexLoadVectorObject
            self._read_table(data, storage_obj, real_obj, complex_obj, 'node')
        elif self.thermal == 1:
            result_name = 'thermalLoadVectors'
            storage_obj = self.thermalLoadVectors
            real_obj = ThermalLoadVectorObject
            complex_obj = None
            self._read_table(data, storage_obj, real_obj, complex_obj, 'node')
        else:
            raise NotImplementedError(self.thermal)

    def read_force_vector(self, data):
        """
        table_code = 12
        """
        if self.thermal == 0:
            result_name = 'forceVectors'
            storage_obj = self.forceVectors
            real_obj = ForceVectorObject
            complex_obj = ComplexForceVectorObject
            self._read_table(data, storage_obj, real_obj, complex_obj, 'node')
        #elif self.thermal == 1:
            #result_name = 'thermalForceVectors'
            #storage_obj = self.thermalForceVectors
            #real_obj = ThermalForceVectorObject
            #complex_obj = None
            #self.read_table(data, storage_obj, real_obj, complex_obj, 'node')
        else:
            raise NotImplementedError(self.thermal)