#pylint: disable=C0301,C0103
"""
Defines the Real/Complex Forces created by:
    OLOAD = ALL
"""
from __future__ import annotations
from typing import TYPE_CHECKING
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

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


class OPG:
    def __init__(self, OP2: OP2):
        self.op2 = OP2

    @property
    def size(self) -> int:
        return self.op2.size
    @property
    def factor(self) -> int:
        return self.op2.factor

    def _read_opg1_3(self, data: bytes, ndata: int) -> None:
        op2 = self.op2
        op2.words = [
            'aCode', 'tCode', '???', 'isubcase',
            '???', '???', '???', 'dLoadID',
            'format_code', 'num_wide', 'o_code', '???',
            'acoustic_flag', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', 'thermal', '???',
            '???', 'Title', 'subtitle', 'label']

        op2.parse_approach_code(data)
        #isubcase = self.get_values(data,'i',4)

        ## dynamic load set ID/random code
        op2.random_code = op2.add_data_parameter(data, 'random_code', b'i', 8, False)

        ## format code
        op2.format_code = op2.add_data_parameter(data, 'format_code', b'i', 9, False)

        ## number of words per entry in record
        ## .. note:: is this needed for this table ???
        op2.num_wide = op2.add_data_parameter(data, 'num_wide', b'i', 10, False)

        ## undefined in DMAP...
        op2.oCode = op2.add_data_parameter(data, 'oCode', b'i', 11, False)

        ## thermal flag; 1 for heat transfer, 0 otherwise
        op2.thermal = op2.add_data_parameter(data, 'thermal', b'i', 23, False)

        #print "dLoadID(8)=%s format_code(9)=%s num_wide(10)=%s oCode(11)=%s thermal(23)=%s" %(self.dLoadID,op2.format_code,op2.num_wide,op2.oCode,op2.thermal)
        if not op2.is_sort1:
            raise NotImplementedError('sort2...')
        #assert self.isThermal()==False,op2.thermal

        ## assuming tCode=1
        if op2.analysis_code == 1:   # statics
            ## load set number
            op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
            op2.setNullNonlinearFactor()
        elif op2.analysis_code == 2:  # normal modes/buckling (real eigenvalues)
            ## mode number
            op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            ## eigenvalue
            op2.eign = op2.add_data_parameter(data, 'eign', b'f', 6, False)
            ## mode or cycle .. todo:: confused on the type - F1???
            op2.mode2 = op2.add_data_parameter(data, 'mode2', b'i', 7, False)
            op2.cycle = op2.add_data_parameter(data, 'cycle', b'f', 7, False)
            op2._op2_readers.reader_oug.update_mode_cycle('cycle')
            op2.data_names = op2.apply_data_code_value('data_names', ['mode', 'eign', 'mode2', 'cycle', ])
        #elif op2.analysis_code == 3: # differential stiffness
        #    ## load set number
        #    op2.lsdvmn = self.get_values(data,'i',5)
        #elif op2.analysis_code == 4: # differential stiffness
        #    ## load set number
        #    op2.lsdvmn = self.get_values(data,'i',5)
        elif op2.analysis_code == 5:   # frequency
            ## frequency
            op2.freq = op2.add_data_parameter(data, 'freq', b'f', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['freq'])
        elif op2.analysis_code == 6:  # transient
            ## time step
            self.time = op2.add_data_parameter(data, 'time', b'f', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['time'])
        elif op2.analysis_code == 7:  # pre-buckling
            ## load set number
            op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
        elif op2.analysis_code == 8:  # post-buckling
            ## load set number
            op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            ## real eigenvalue
            op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn', 'eigr'])
        elif op2.analysis_code == 9:  # complex eigenvalues
            ## mode number
            op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            ## real eigenvalue
            op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            ## imaginary eigenvalue
            op2.eigi = op2.add_data_parameter(data, 'eigi', b'f', 7, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['mode', 'eigr', 'eigi'])
        elif op2.analysis_code == 10:  # nonlinear statics
            ## load step
            op2.lftsfq = op2.add_data_parameter(data, 'lftsfq', b'f', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['lftsfq'])
        elif op2.analysis_code == 11:  # old geometric nonlinear statics
            ## load set number
            op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
        elif op2.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## load set number
            op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
        else:
            raise RuntimeError('invalid analysis_code...analysis_code=%s' %
                               op2.analysis_code)

        # tCode=2
        #if op2.analysis_code==2: # sort2
        #    op2.lsdvmn = self.get_values(data, b'i', 5) ## load set, Mode number

        #print "*isubcase=%s"%(op2.isubcase)
        #print "analysis_code=%s table_code=%s thermal=%s" % (
            #op2.analysis_code, op2.table_code, op2.thermal)

        op2._fix_oug_format_code()
        op2._parse_thermal_code()
        if op2.is_debug_file:
            op2.binary_debug.write('  approach_code  = %r\n' % op2.approach_code)
            op2.binary_debug.write('  tCode          = %r\n' % op2.tCode)
            op2.binary_debug.write('  isubcase       = %r\n' % op2.isubcase)
        op2._read_title(data)
        op2._write_debug_bits()

    def _read_opg2_3(self, data: bytes, ndata: int):
        """reads the SORT2 version of table 4 (the data table)"""
        op2 = self.op2
        op2.nonlinear_factor = np.nan
        op2.is_table_1 = False
        op2.is_table_2 = True
        unused_three = op2.parse_approach_code(data)
        op2.words = [
            'analysis_code', 'table_code', '???', 'isubcase',
            '???', '???', '???', 'random_code',
            'format_code', 'num_wide', '11', '???',
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
        op2.acoustic_flag = op2.add_data_parameter(data, 'acoustic_flag', b'f', 13, False)

        ## thermal flag; 1 for heat transfer, 0 otherwise
        op2.thermal = op2.add_data_parameter(data, 'thermal', b'i', 23, False)

        #assert self.isThermal() == False, op2.thermal

        op2.node_id = op2.add_data_parameter(data, 'node_id', b'i', 5, fix_device_code=True)
        #if op2.analysis_code == 1:  # statics / displacement / heat flux
            ## load set number
            #op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5, False)
            #op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            #op2.setNullNonlinearFactor()

        if op2.analysis_code == 1:  # static...because reasons.
            op2._analysis_code_fmt = b'i'
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            op2.apply_data_code_value('analysis_method', 'N/A')
        elif op2.analysis_code == 2:  # real eigenvalues
            ## mode number
            op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            op2._analysis_code_fmt = b'i'
            ## real eigenvalue
            op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            ## mode or cycle .. todo:: confused on the type - F1???
            op2.mode_cycle = op2.add_data_parameter(data, 'mode_cycle', b'f', 7, False)
            op2.data_names = op2.apply_data_code_value('data_names',
                                                         ['node_id', 'eigr', 'mode_cycle'])
            op2.apply_data_code_value('analysis_method', 'mode')
        #elif op2.analysis_code == 3: # differential stiffness
            #op2.lsdvmn = self.get_values(data, b'i', 5) ## load set number
            #op2.data_names = self.data_code['lsdvmn'] = op2.lsdvmn
        #elif op2.analysis_code == 4: # differential stiffness
            #op2.lsdvmn = self.get_values(data, b'i', 5) ## load set number
        elif op2.analysis_code == 5:   # frequency
            ## frequency
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
            op2._analysis_code_fmt = b'i'
            ## real eigenvalue
            op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id', 'eigr'])
            op2.apply_data_code_value('analysis_method', 'eigr')
        elif op2.analysis_code == 9:  # complex eigenvalues
            ## mode number
            op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            op2._analysis_code_fmt = b'i'
            ## real eigenvalue
            #op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
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
            op2._analysis_code_fmt = b'f'
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
        elif op2.analysis_code == 12:
            # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## load set number
            #op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            op2._analysis_code_fmt = b'i'
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            op2.apply_data_code_value('analysis_method', 'lsdvmn')
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % op2.analysis_code
            raise RuntimeError(msg)

        op2._fix_oug_format_code()
        op2._parse_thermal_code()
        if op2.is_debug_file:
            op2.binary_debug.write('  approach_code  = %r\n' % op2.approach_code)
            op2.binary_debug.write('  tCode          = %r\n' % op2.tCode)
            op2.binary_debug.write('  isubcase       = %r\n' % op2.isubcase)
        op2._read_title(data)
        op2._write_debug_bits()

    def _read_opg1_4(self, data: bytes, ndata: int) -> int:
        op2 = self.op2
        if op2.table_code == 2:  # load vector
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
            prefix = prefixs[op2.table_name]
            postfix = postfixs[op2.table_name]
            assert op2.table_name in keys, 'tables=%s table_name=%s table_code=%s' % (keys, op2.table_name, op2.table_code)
            n = self._read_load_vector(data, ndata, prefix=prefix, postfix=postfix)
        elif op2.table_code == 12:  # ???
            n = self._read_force_vector(data, ndata)
        elif op2.table_code == 502:  # load vector
            assert op2.table_name in [b'OPGCRM2'], f'table_name={op2.table_name} table_code={op2.table_code}'
            n = self._read_load_vector(data, ndata, prefix='crm.')
        elif op2.table_code == 602:  # load vector
            assert op2.table_name in [b'OPGPSD2'], f'table_name={op2.table_name} table_code={op2.table_code}'
            n = self._read_load_vector(data, ndata, prefix='psd.')
        elif op2.table_code == 802:  # load vector
            assert op2.table_name in [b'OPGRMS1'], f'table_name={op2.table_name} table_code={op2.table_code}'
            n = self._read_load_vector(data, ndata, prefix='rms.')
        elif op2.table_code == 902:  # load vector
            assert op2.table_name in [b'OPGNO1'], f'table_name={op2.table_name} table_code={op2.table_code}'
            n = self._read_load_vector(data, ndata, prefix='no.')
        #else:
            #n = self._not_implemented_or_skip('bad OPG table')
        else:
            raise NotImplementedError('table_name=%r table_code=%r' % (op2.table_name, op2.table_code))
        return n

    def _read_load_vector(self, data, ndata, prefix='', postfix='') -> int:
        """
        table_code = 2
        """
        op2 = self.op2
        if op2.thermal == 0:
            result_name = prefix + 'load_vectors' + postfix
            storage_obj = op2.get_result(result_name)
            if op2._results.is_not_saved(result_name):
                return ndata
            op2._results._found_result(result_name)
            n = op2._read_table_vectorized(data, ndata, result_name, storage_obj,
                                           RealLoadVectorArray, ComplexLoadVectorArray,
                                           'node', random_code=op2.random_code)
        elif op2.thermal == 1:
            result_name = prefix + 'thermal_load_vectors' + postfix
            storage_obj = op2.get_result(result_name)

            #RealThermalLoadVectorVector = None
            #ComplexThermalLoadVectorVector = None
            ComplexThermalLoadVectorArray = None
            if op2._results.is_not_saved(result_name):
                return ndata
            op2._results._found_result(result_name)
            n = op2._read_scalar_table_vectorized(data, ndata, result_name, storage_obj,
                                                  RealTemperatureVectorArray, ComplexThermalLoadVectorArray,
                                                  'node', random_code=op2.random_code)

        else:
            raise NotImplementedError(op2.thermal)
        return n

    def _read_force_vector(self, data: bytes, ndata: int) -> int:
        """
        table_code = 12
        """
        op2 = self.op2
        if op2.thermal == 0:
            result_name = 'force_vectors'
            storage_obj = op2.get_result(result_name)
            #ForceVectorVector = None
            ComplexForceVectorArray = None
            if op2._results.is_not_saved(result_name):
                return ndata
            op2._results._found_result(result_name)
            n = op2._read_table_vectorized(data, ndata, result_name, storage_obj,
                                           RealForceVectorArray, ComplexForceVectorArray,
                                           'node', random_code=op2.random_code)

        #elif op2.thermal == 1:
            #result_name = 'thermal_force_vectors'
            #storage_obj = op2.thermal_force_vectors
            #RealThermalForceVector = None
            #ComplexThermalForceVector = None
            #RealThermalForceVectorVector = None
            #ComplexThermalForceVectorArray = None
            #n = self._read_table(data, result_name, storage_obj,
                                 #RealThermalForceVector, ComplexThermalForceVector,
                                 #RealThermalForceVectorArray, ComplexThermalForceVectorArray,
                                 #'node', random_code=op2.random_code)
        else:
            raise NotImplementedError(op2.thermal)
        return n
