#pylint: disable=C0301,C0103
"""
Defines the Real/Complex Forces created by:
    OLOAD = ALL
"""
import numpy as np
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

from pyNastran.op2.op2_interface.op2_common import OP2Common


class OPG(OP2Common):
    def __init__(self):
        pass

    def _read_opg1_3(self, data: bytes, ndata: int) -> None:
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
        self.random_code = self.add_data_parameter(data, 'random_code', b'i', 8, False)

        ## format code
        self.format_code = self.add_data_parameter(data, 'format_code', b'i', 9, False)

        ## number of words per entry in record
        ## .. note:: is this needed for this table ???
        self.num_wide = self.add_data_parameter(data, 'num_wide', b'i', 10, False)

        ## undefined in DMAP...
        self.oCode = self.add_data_parameter(data, 'oCode', b'i', 11, False)

        ## thermal flag; 1 for heat transfer, 0 otherwise
        self.thermal = self.add_data_parameter(data, 'thermal', b'i', 23, False)

        #print "dLoadID(8)=%s format_code(9)=%s num_wide(10)=%s oCode(11)=%s thermal(23)=%s" %(self.dLoadID,self.format_code,self.num_wide,self.oCode,self.thermal)
        if not self.is_sort1:
            raise NotImplementedError('sort2...')
        #assert self.isThermal()==False,self.thermal

        ## assuming tCode=1
        if self.analysis_code == 1:   # statics
            ## load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5, False)
            self.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysis_code == 2:  # normal modes/buckling (real eigenvalues)
            ## mode number
            self.mode = self.add_data_parameter(data, 'mode', b'i', 5)
            ## eigenvalue
            self.eign = self.add_data_parameter(data, 'eign', b'f', 6, False)
            ## mode or cycle .. todo:: confused on the type - F1???
            self.mode2 = self.add_data_parameter(data, 'mode2', b'i', 7, False)
            self.cycle = self.add_data_parameter(data, 'cycle', b'f', 7, False)
            self.update_mode_cycle('cycle')
            self.data_names = self.apply_data_code_value('data_names', ['mode', 'eign', 'mode2', 'cycle', ])
        #elif self.analysis_code == 3: # differential stiffness
        #    ## load set number
        #    self.lsdvmn = self.get_values(data,'i',5)
        #elif self.analysis_code == 4: # differential stiffness
        #    ## load set number
        #    self.lsdvmn = self.get_values(data,'i',5)
        elif self.analysis_code == 5:   # frequency
            ## frequency
            self.freq = self.add_data_parameter(data, 'freq', b'f', 5)
            self.data_names = self.apply_data_code_value('data_names', ['freq'])
        elif self.analysis_code == 6:  # transient
            ## time step
            self.time = self.add_data_parameter(data, 'time', b'f', 5)
            self.data_names = self.apply_data_code_value('data_names', ['time'])
        elif self.analysis_code == 7:  # pre-buckling
            ## load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5)
            self.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
        elif self.analysis_code == 8:  # post-buckling
            ## load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5)
            ## real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', b'f', 6, False)
            self.data_names = self.apply_data_code_value('data_names', ['lsdvmn', 'eigr'])
        elif self.analysis_code == 9:  # complex eigenvalues
            ## mode number
            self.mode = self.add_data_parameter(data, 'mode', b'i', 5)
            ## real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', b'f', 6, False)
            ## imaginary eigenvalue
            self.eigi = self.add_data_parameter(data, 'eigi', b'f', 7, False)
            self.data_names = self.apply_data_code_value('data_names', ['mode', 'eigr', 'eigi'])
        elif self.analysis_code == 10:  # nonlinear statics
            ## load step
            self.lftsfq = self.add_data_parameter(data, 'lftsfq', b'f', 5)
            self.data_names = self.apply_data_code_value('data_names', ['lftsfq'])
        elif self.analysis_code == 11:  # old geometric nonlinear statics
            ## load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5)
            self.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
        elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5)
            self.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
        else:
            raise RuntimeError('invalid analysis_code...analysis_code=%s' %
                               self.analysis_code)

        # tCode=2
        #if self.analysis_code==2: # sort2
        #    self.lsdvmn = self.get_values(data, b'i', 5) ## load set, Mode number

        #print "*isubcase=%s"%(self.isubcase)
        #print "analysis_code=%s table_code=%s thermal=%s" % (
            #self.analysis_code, self.table_code, self.thermal)

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

    def _read_opg2_3(self, data, ndata):
        """reads the SORT2 version of table 4 (the data table)"""
        self.nonlinear_factor = np.nan
        self.is_table_1 = False
        self.is_table_2 = True
        unused_three = self.parse_approach_code(data)
        self.words = [
            'analysis_code', 'table_code', '???', 'isubcase',
            '???', '???', '???', 'random_code',
            'format_code', 'num_wide', '11', '???',
            'acoustic_flag', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', 'thermal', '???',
            '???', 'Title', 'subtitle', 'label']

        ## random code
        self.random_code = self.add_data_parameter(data, 'random_code', b'i', 8, False)

        ## format code
        self.format_code = self.add_data_parameter(data, 'format_code', b'i', 9, False)

        ## number of words per entry in record
        self.num_wide = self.add_data_parameter(data, 'num_wide', b'i', 10, False)

        ## acoustic pressure flag
        self.acoustic_flag = self.add_data_parameter(data, 'acoustic_flag', b'f', 13, False)

        ## thermal flag; 1 for heat transfer, 0 otherwise
        self.thermal = self.add_data_parameter(data, 'thermal', b'i', 23, False)

        #assert self.isThermal() == False, self.thermal

        self.node_id = self.add_data_parameter(data, 'node_id', b'i', 5, fix_device_code=True)
        #if self.analysis_code == 1:  # statics / displacement / heat flux
            ## load set number
            #self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5, False)
            #self.data_names = self.apply_data_code_value('data_names', ['node_id'])
            #self.setNullNonlinearFactor()

        if self.analysis_code == 1:  # static...because reasons.
            self._analysis_code_fmt = b'i'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
            self.apply_data_code_value('analysis_method', 'N/A')
        elif self.analysis_code == 2:  # real eigenvalues
            ## mode number
            self.mode = self.add_data_parameter(data, 'mode', b'i', 5)
            self._analysis_code_fmt = b'i'
            ## real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', b'f', 6, False)
            ## mode or cycle .. todo:: confused on the type - F1???
            self.mode_cycle = self.add_data_parameter(data, 'mode_cycle', b'f', 7, False)
            self.data_names = self.apply_data_code_value('data_names',
                                                         ['node_id', 'eigr', 'mode_cycle'])
            self.apply_data_code_value('analysis_method', 'mode')
        #elif self.analysis_code == 3: # differential stiffness
            #self.lsdvmn = self.get_values(data, b'i', 5) ## load set number
            #self.data_names = self.data_code['lsdvmn'] = self.lsdvmn
        #elif self.analysis_code == 4: # differential stiffness
            #self.lsdvmn = self.get_values(data, b'i', 5) ## load set number
        elif self.analysis_code == 5:   # frequency
            ## frequency
            #self.freq = self.add_data_parameter(data, 'freq', b'f', 5)
            self._analysis_code_fmt = b'f'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
            self.apply_data_code_value('analysis_method', 'freq')
        elif self.analysis_code == 6:  # transient
            ## time step
            #self.dt = self.add_data_parameter(data, 'dt', b'f', 5)
            self._analysis_code_fmt = b'f'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
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

    def _read_opg1_4(self, data, ndata):
        if self.table_code == 2:  # load vector
            prefixs = {
                b'BOPG1' : '',
                b'OPG1' : '',
                b'OPG2' : '',
                b'OPGV1' : '',
                b'OCRPG' : '',
                b'OPGPSD1' : 'psd.',
                b'OPGPSD2' : 'psd.',
                b'OPGATO1' : 'ato.',
                b'OPGATO2' : 'ato.',
                b'OPGCRM1' : 'crm.',
                b'OPGCRM2' : 'crm.',
                b'OPGRMS1' : 'rms.',
                b'OPGRMS2' : 'rms.',
                b'OPGNO1' : 'no.',
                b'OPGNO2' : 'no.',
            }
            postfixs = {
                b'BOPG1' : '',
                b'OPG1' : '',
                b'OPG2' : '',
                #b'OPGV1' : '_v',
                b'OPGV1' : '',
                b'OCRPG' : '',
                b'OPGPSD1' : '',
                b'OPGPSD2' : '',
                b'OPGATO1' : '',
                b'OPGATO2' : '',
                b'OPGCRM1' : '',
                b'OPGCRM2' : '',
                b'OPGRMS1' : '',
                b'OPGRMS2' : '',
                b'OPGNO1' : '',
                b'OPGNO2' : '',
            }
            keys = list(prefixs.keys())
            prefix = prefixs[self.table_name]
            postfix = postfixs[self.table_name]
            assert self.table_name in keys, 'tables=%s table_name=%s table_code=%s' % (keys, self.table_name, self.table_code)
            n = self._read_load_vector(data, ndata, prefix=prefix, postfix=postfix)
        elif self.table_code == 12:  # ???
            n = self._read_force_vector(data, ndata)
        elif self.table_code == 502:  # load vector
            assert self.table_name in [b'OPGCRM2'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
            n = self._read_load_vector(data, ndata, prefix='crm.')
        elif self.table_code == 602:  # load vector
            assert self.table_name in [b'OPGPSD2'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
            n = self._read_load_vector(data, ndata, prefix='psd.')
        elif self.table_code == 802:  # load vector
            assert self.table_name in [b'OPGRMS1'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
            n = self._read_load_vector(data, ndata, prefix='rms.')
        elif self.table_code == 902:  # load vector
            assert self.table_name in [b'OPGNO1'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
            n = self._read_load_vector(data, ndata, prefix='no.')
        #else:
            #n = self._not_implemented_or_skip('bad OPG table')
        else:
            raise NotImplementedError('table_name=%r table_code=%r' % (self.table_name, self.table_code))
        return n

    def _read_load_vector(self, data, ndata, prefix='', postfix=''):
        """
        table_code = 2
        """
        if self.thermal == 0:
            result_name = prefix + 'load_vectors' + postfix
            storage_obj = self.get_result(result_name)
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            n = self._read_table_vectorized(data, ndata, result_name, storage_obj,
                                            RealLoadVectorArray, ComplexLoadVectorArray,
                                            'node', random_code=self.random_code)
        elif self.thermal == 1:
            result_name = prefix + 'thermal_load_vectors' + postfix
            storage_obj = self.get_result(result_name)

            #RealThermalLoadVectorVector = None
            #ComplexThermalLoadVectorVector = None
            ComplexThermalLoadVectorArray = None
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
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
            storage_obj = self.get_result(result_name)
            #ForceVectorVector = None
            ComplexForceVectorArray = None
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
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
