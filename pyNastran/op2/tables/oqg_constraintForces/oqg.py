"""
This file defines the OUG Table, which contains:
 * Real/Complex SPC Forces
   - SPCFORCE = ALL
 * Real/Complex MPC Forces
   - MPCFORCE = ALL
 * Real Temperature Gradient & Flux
   - FLUX = ALL

"""
from __future__ import annotations
from struct import Struct
from typing import Dict, Any, TYPE_CHECKING
import numpy as np
from pyNastran.op2.op2_interface.op2_reader import mapfmt

from pyNastran.op2.tables.oug.oug import get_shock_prefix_postfix
from pyNastran.op2.tables.oqg_constraintForces.separation_distance import (
    SeparationDistanceArray)
from pyNastran.op2.tables.oqg_constraintForces.oqg_spc_forces import (
    RealSPCForcesArray, ComplexSPCForcesArray,)
from pyNastran.op2.tables.oqg_constraintForces.oqg_mpc_forces import (
    RealMPCForcesArray, ComplexMPCForcesArray,)
from pyNastran.op2.tables.oqg_constraintForces.oqg_contact_forces import RealContactForcesArray
from pyNastran.op2.tables.oqg_constraintForces.oqg_thermal_gradient_and_flux import (
    RealTemperatureGradientAndFluxArray)
from pyNastran.op2.result_objects.contact_traction_and_pressure import RealContactTractionAndPressureArray
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2

class OQG:
    def __init__(self, op2: OP2):
        self.op2 = op2

    @property
    def size(self) -> int:
        return self.op2.size
    @property
    def factor(self) -> int:
        return self.op2.factor

    def _read_opsdi1_3(self, data: bytes, ndata: int) -> None:
        """Initial separation distance"""
        op2 = self.op2
        op2.to_nx(f' because table_name={op2.table_name} was found')
        self._read_oqg1_3(data, ndata)

    def _read_opsdi1_4(self, data: bytes, ndata: int) -> int:
        """Initial separation distance"""
        op2 = self.op2
        return self._read_opsds1_intial_final_4(data, ndata,
                                                'separation_initial',
                                                op2.op2_results.separation_initial,
                                                SeparationDistanceArray)

    def _read_opsds1_3(self, data: bytes, ndata: int) -> None:
        """Final separation distance"""
        op2 = self.op2
        op2.to_nx(f' because table_name={op2.table_name} was found')
        self._read_oqg1_3(data, ndata)

    def _read_opsds1_4(self, data: bytes, ndata: int) -> int:
        """Final separation distance"""
        op2 = self.op2
        return self._read_opsds1_intial_final_4(data, ndata,
                                                'separation_final',
                                                op2.op2_results.separation_final,
                                                SeparationDistanceArray)

    def _read_opsds1_intial_final_4(self, data: bytes, ndata: int,
                                    result_name: str,
                                    slot: Dict[Any, Any],
                                    obj_vector_real: SeparationDistanceArray) -> int:
        """
        $                    D E F O R M E D  C O N T A C T  S E P A R A T I O N  D I S T A N C E
        $
        $      POINT ID.   TYPE       DISTANCE
        $             1      G      1.010815E-01
        $             6      G      9.858034E-02
        """
        op2 = self.op2
        assert op2.num_wide == 2, op2.code_information()
        #if op2._results.is_not_saved(result_name):
            #return ndata
        #op2._results._found_result(result_name)
        #slot = op2.get_result(result_name)
        ntotal = 2 * 4 * self.factor
        nnodes = ndata // ntotal

        auto_return, is_vectorized = op2._create_node_object4(
            nnodes, result_name, slot, obj_vector_real)
        if auto_return:
            return nnodes * ntotal

        obj = op2.obj
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nnodes * 4 * op2.num_wide
            itotal = obj.itotal
            itotal2 = obj.itotal + nnodes
            ints = np.frombuffer(data, op2.idtype8)[::2]
            floats = np.frombuffer(data, op2.fdtype8)[1::2]
            obj.node = ints // 10
            obj.data[obj.itime, itotal:itotal2, 0] = floats
        else:
            n = 0
            dt = op2.nonlinear_factor
            structi = Struct(op2._endian + mapfmt(b'if', self.size))
            for unused_i in range(nnodes):
                edata = data[n:n+ntotal]
                nid_code, distance = structi.unpack(edata)
                nid = nid_code // 10
                obj.add_sort1(dt, nid, distance)
                n += ntotal
        return ndata

    def _read_oqg1_3(self, data: bytes, ndata: int) -> None:
        op2 = self.op2
        op2.nonlinear_factor = np.nan #None
        op2.is_table_1 = True
        op2.is_table_2 = False
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

        if not op2.is_sort1:
            raise NotImplementedError('SORT2; code_info=\n%s' % op2.code_information())
        #assert self.isThermal()==False,op2.thermal

        ## assuming tCode=1
        if op2.analysis_code == 1:   # statics / displacement / heat flux
            ## load set number
            op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
            op2.setNullNonlinearFactor()
        elif op2.analysis_code == 2:  # real eigenvalues
            ## mode number
            op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            ## eigenvalue
            op2.eign = op2.add_data_parameter(data, 'eign', b'f', 6, False)
            ## mode or cycle .. todo:: confused on the type - F1???
            op2.mode_cycle = op2.add_data_parameter(data, 'mode_cycle', b'f', 7, False)
            op2.reader_oug.update_mode_cycle('mode_cycle')
            op2.data_names = op2.apply_data_code_value('data_names',
                                                         ['mode', 'eign', 'mode_cycle'])
        #elif op2.analysis_code == 3: # differential stiffness
            #op2.lsdvmn = self.get_values(data, b'i', 5) ## load set number
            #op2.data_names = op2.data_code['lsdvmn'] = op2.lsdvmn
        #elif op2.analysis_code == 4: # differential stiffness
            #op2.lsdvmn = self.get_values(data, b'i', 5) ## load set number
        elif op2.analysis_code == 5:   # frequency
            ## frequency
            op2.freq = op2.add_data_parameter(data, 'freq', b'f', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['freq'])
        elif op2.analysis_code == 6:  # transient
            ## time step
            op2.dt = op2.add_data_parameter(data, 'dt', b'f', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['dt'])
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
        elif op2.analysis_code == 12:
            # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## load set number
            op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
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

    def _read_oqg2_3(self, data: bytes, ndata: int) -> None:
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

        self.node_id = op2.add_data_parameter(data, 'node_id', b'i', 5, fix_device_code=True)
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
            #op2.data_names = op2.data_code['lsdvmn'] = op2.lsdvmn
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

    def _read_oqg_4(self, data: bytes, ndata: int):
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
        op2 = self.op2
        result_name = 'constraint_forces'
        if op2._results.is_not_saved(result_name):
            return ndata

        if op2.table_name in [b'OQGCRM1', b'OQGCRM2']:
            if op2.table_code not in [3, 503]:
                msg = f'table_name={op2.table_name} table_code={op2.table_code}'
                raise AssertionError(msg)
            n = self._read_oqg_spc_crm(data, ndata)
        elif op2.table_name in [b'OQGPSD1', b'OQGPSD2']:
            if op2.table_code not in [3, 603]:  # was 3
                msg = f'table_name={op2.table_name} table_code={op2.table_code}'
                raise AssertionError(msg + '\n%s' % op2.code_information())
            n = self._read_oqg_spc_psd(data, ndata)
        elif op2.table_name in [b'OQGATO1', b'OQGATO2']:
            if op2.table_code not in [3, 703]:
                msg = f'table_name={op2.table_name} table_code={op2.table_code}'
                raise AssertionError(msg)
            n = self._read_oqg_spc_ato(data, ndata)
        elif op2.table_name in [b'OQGRMS1', b'OQGRMS2']:
            if op2.table_code not in [3, 803]:
                msg = f'table_name={op2.table_name} table_code={op2.table_code}'
                raise AssertionError(msg)
            n = self._read_oqg_spc_rms(data, ndata)
        elif op2.table_name in [b'OQGNO1', b'OQGNO2']:
            if op2.table_code not in [3, 903]:
                msg = f'table_name={op2.table_name} table_code={op2.table_code}'
                raise AssertionError(msg)
            n = self._read_oqg_spc_no(data, ndata)

        elif op2.table_code == 3:   # SPC Forces
            assert op2.table_name in [b'OQG1', b'OQGV1', b'OQP1', b'OQG2'], op2.code_information()
            n = self._read_spc_forces(data, ndata)
        elif op2.table_code == 39:  # MPC Forces
            assert op2.table_name in [b'OQMG1', b'OQMG2'], op2.code_information() # , b'OQMPSD1', b'OQMPSD2'
            n = self._read_oqg_mpc_forces(data, ndata)
        elif op2.table_name in [b'RAQEATC', b'RAQCONS']:
            # op2.table_code == 5 and
            n = self._read_oqg_mpc_forces(data, ndata)
        elif op2.table_code == 63:   # Contact Forces
            op2.to_nx(f' because table_name={op2.table_name} was found')
            assert op2.table_name in [b'OQGCF1', b'OQGCF2'], op2.code_information()
            n = self._oqg_read_contact_forces(data, ndata)
        elif op2.table_code == 67:   # Glue Forces
            assert op2.table_name in [b'OQGGF1', b'OQGGF2'], op2.code_information()
            n = self._oqg_read_glue_forces(data, ndata)
        else:
            raise RuntimeError(op2.code_information())
        return n

    def _read_spc_forces(self, data: bytes, ndata: int):
        """
        table_code = 3
        """
        op2 = self.op2
        if op2.table_name in [b'OQG1', b'OQG2', b'OQGV1', b'OQP1']:
            pass
        else:  # pragma: no cover
            msg = 'spc_forces; table_name=%s' % op2.table_name
            raise NotImplementedError(msg)

        if op2.thermal == 0:
            op2._setup_op2_subcase('SPCFORCES')
            if op2.table_name in [b'OQG1', b'OQG2']:
                result_name = 'spc_forces'
            #elif op2.table_name in [b'OQGV1']:
                #result_name = 'spc_forces_v'
            elif op2.table_name in [b'OQGV1', b'OQP1']:
                result_name = 'spc_forces'
            else:
                raise NotImplementedError(op2.table_name_str)
            if op2._results.is_not_saved(result_name):
                return ndata
            op2._results._found_result(result_name)
            storage_obj = getattr(op2, result_name)
            n = op2._read_table_vectorized(data, ndata, result_name, storage_obj,
                                           RealSPCForcesArray, ComplexSPCForcesArray,
                                           'node', random_code=op2.random_code)
        elif op2.thermal == 1:
            #'finite element temperature gradients and fluxes'
            op2._setup_op2_subcase('FLUX')

            result_name = 'thermal_gradient_and_flux'
            storage_obj = op2.thermal_gradient_and_flux
            if op2._results.is_not_saved(result_name):
                return ndata
            op2._results._found_result(result_name)
            n = op2._read_table_vectorized(data, ndata, result_name, storage_obj,
                                            RealTemperatureGradientAndFluxArray, None,
                                            'node', random_code=op2.random_code)
        elif op2.thermal == 8:  # 4 ?
            result_name0 = 'spc_forces'
            prefix, postfix = get_shock_prefix_postfix(op2.thermal)
            result_name = prefix + result_name0 + postfix
            if op2._results.is_not_saved(result_name):
                return ndata
            op2._results._found_result(result_name)
            storage_obj = op2.get_result(result_name)
            n = op2._read_table_vectorized(data, ndata, result_name, storage_obj,
                                           RealSPCForcesArray, ComplexSPCForcesArray,
                                           'node', random_code=op2.random_code)
        else:
            raise RuntimeError(op2.code_information())
            #msg = 'thermal=%s' % op2.thermal
            #return self._not_implemented_or_skip(data, ndata, msg)
        return n

    def _read_oqg_mpc_forces(self, data: bytes, ndata: int):
        """
        table_code = 39
        """
        op2 = self.op2
        op2._setup_op2_subcase('MPCFORCES')

        if op2.table_name in [b'OQMG1', b'OQMG2']:
            result_name = 'mpc_forces'
        elif op2.table_name == b'RAQEATC':
            result_name = 'mpc_forces_RAQEATC'
        elif op2.table_name == b'RAQCONS':
            result_name = 'mpc_forces_RAQCONS'
        elif op2.table_name == b'RAQEATC':
            result_name = 'mpc_forces_RAQEATC'
        else:  # pragma: no cover
            msg = 'mpc_forces; table_name=%s' % op2.table_name
            raise NotImplementedError(msg)

        storage_obj = op2.get_result(result_name)
        if op2.thermal == 0:
            if op2._results.is_not_saved(result_name):
                return ndata
            op2._results._found_result(result_name)
            n = op2._read_table_vectorized(data, ndata, result_name, storage_obj,
                                           RealMPCForcesArray, ComplexMPCForcesArray,
                                           'node', random_code=op2.random_code)
        else:
            raise RuntimeError(op2.code_information())
            #msg = 'thermal=%s' % op2.thermal
            #return self._not_implemented_or_skip(data, ndata, msg)
        assert isinstance(n, int), 'table_name=%s n=%s' % (op2.table_name, n)
        return n

    def _oqg_read_contact_forces(self, data: bytes, ndata: int):
        """
        table_code = 63
        """
        op2 = self.op2
        op2.to_nx(f' because table_name={op2.table_name} was found')
        #op2._setup_op2_subcase('MPCFORCES')
        if op2.table_name in [b'OQGCF1', b'OQGCF2']:
            result_name = 'contact_forces'
        else:  # pragma: no cover
            msg = 'contact_forces; table_name=%s' % op2.table_name
            raise NotImplementedError(msg)

        storage_obj = op2.get_result(result_name)
        if op2.thermal == 0:
            if op2._results.is_not_saved(result_name):
                return ndata
            op2._results._found_result(result_name)
            n = op2._read_table_vectorized(data, ndata, result_name, storage_obj,
                                           RealContactForcesArray, None,
                                           'node', random_code=op2.random_code)
        else:
            raise RuntimeError(op2.code_information())
            #msg = 'thermal=%s' % op2.thermal
            #return self._not_implemented_or_skip(data, ndata, msg)
        assert isinstance(n, int), 'table_name=%s n=%s' % (op2.table_name, n)
        return n

    def _oqg_read_glue_forces(self, data: bytes, ndata: int):
        """
        table_code = 67
        """
        op2 = self.op2
        op2.to_nx(f' because table_name={op2.table_name} was found')
        #op2._setup_op2_subcase('MPCFORCES')
        if op2.table_name in [b'OQGGF1', b'OQGGF2']:
            result_name = 'glue_forces'
        else:  # pragma: no cover
            msg = 'glue_forces; table_name=%s' % op2.table_name
            raise NotImplementedError(msg)

        storage_obj = op2.get_result(result_name)
        if op2.thermal == 0:
            if op2._results.is_not_saved(result_name):
                return ndata
            op2._results._found_result(result_name)
            n = op2._read_table_vectorized(data, ndata, result_name, storage_obj,
                                           RealContactForcesArray, None,
                                           'node', random_code=op2.random_code)
        else:
            raise RuntimeError(op2.code_information())
            #msg = 'thermal=%s' % op2.thermal
            #return self._not_implemented_or_skip(data, ndata, msg)
        assert isinstance(n, int), 'table_name=%s n=%s' % (op2.table_name, n)
        return n

    def _read_oqg_spc_psd(self, data: bytes, ndata: int):
        """
        table_code = 601/610/611
        """
        op2 = self.op2
        if op2.thermal == 0:
            if op2.table_code in [3, 603]:
                result_name = 'psd.spc_forces'
            else:
                raise RuntimeError(op2.code_information())
            obj = RealSPCForcesArray
        else:
            raise RuntimeError(op2.code_information())

        if op2._results.is_not_saved(result_name):
            return ndata
        op2._results._found_result(result_name)

        storage_obj = op2.get_result(result_name)
        n = op2._read_random_table(data, ndata, result_name, storage_obj,
                                   obj, 'node',
                                   random_code=op2.random_code)
        return n

    def _read_oqg_spc_rms(self, data: bytes, ndata: int):
        """
        table_code = 3/???/?10/?11
        """
        op2 = self.op2
        if op2.thermal == 0:
            if op2.table_code in [3, 803]:
                result_name = 'rms.spc_forces'
                if op2._results.is_not_saved(result_name):
                    return ndata
                storage_obj = op2.get_result(result_name)
                op2._results._found_result(result_name)
                n = op2._read_random_table(data, ndata, result_name, storage_obj,
                                           RealSPCForcesArray, 'node',
                                           random_code=op2.random_code)
            else:
                raise RuntimeError(op2.code_information())
        else:
            raise RuntimeError(op2.code_information())
        return n

    def _read_oqg_spc_ato(self, data: bytes, ndata: int):
        """
        table_code = ???
        """
        op2 = self.op2
        if op2.thermal == 0:
            if op2.table_code in [3]:
                assert op2.table_name in [b'OQGATO1', b'OQGATO2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'ato.spc_forces'
            #elif op2.table_code in [603]:
                #assert op2.table_name in [b'OQGATO2'], 'op2.table_name=%r' % op2.table_name
                #result_name = 'psd.mpc_forces'
            else:
                print(op2.table_code)
                raise RuntimeError(op2.code_information())

            if op2._results.is_not_saved(result_name):
                return ndata
            op2._results._found_result(result_name)
            storage_obj = op2.get_result(result_name)
            n = op2._read_random_table(data, ndata, result_name, storage_obj,
                                       RealSPCForcesArray, 'node',
                                       random_code=op2.random_code)
        else:
            raise RuntimeError(op2.code_information())
        assert n is not None, n
        return n

    def _read_oqg_spc_crm(self, data: bytes, ndata: int) -> int:
        """
        table_code = 3/???/?10/?11
        """
        op2 = self.op2
        if op2.thermal == 0:
            if op2.table_code in [3, 503]:
                result_name = 'crm.spc_forces'
                storage_obj = op2.get_result(result_name)
                if op2._results.is_not_saved(result_name):
                    return ndata
                op2._results._found_result(result_name)
                n = op2._read_random_table(data, ndata, result_name, storage_obj,
                                           RealSPCForcesArray, 'node',
                                           random_code=op2.random_code)
            else:
                raise RuntimeError(op2.code_information())
        else:
            raise RuntimeError(op2.code_information())
        return n

    def _read_oqg_spc_no(self, data: bytes, ndata: int) -> int:
        """
        table_code = 3/???/?10/?11
        """
        op2 = self.op2
        if op2.thermal == 0:
            if op2.table_code in [3, 903]:
                assert op2.table_name in [b'OQGNO1', b'OQGNO2'], f'op2.table_name={op2.table_name!r}'
                result_name = 'no.spc_forces'
                obj = RealSPCForcesArray
            else:
                raise RuntimeError(op2.code_information())
        else:
            raise RuntimeError(op2.code_information())

        if op2._results.is_not_saved(result_name):
            return ndata
        op2._results._found_result(result_name)

        storage_obj = op2.get_result(result_name)
        n = op2._read_random_table(data, ndata, result_name, storage_obj,
                                   obj, 'node',
                                   random_code=op2.random_code)
        return n

    def _read_oqg_mpc_psd(self, data: bytes, ndata: int) -> int:
        """
        table_code = 603
        """
        op2 = self.op2
        if op2.thermal == 0:
            if op2.table_code in [39]:
                assert op2.table_name in [b'OQMPSD1', b'OQMPSD2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'psd.mpc_forces'
            elif op2.table_code in [603]:
                assert op2.table_name in [b'OQMPSD1', b'OQMPSD2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'psd.mpc_forces'
            else:
                print(op2.table_code)
                raise RuntimeError(op2.code_information())
            obj = RealMPCForcesArray
        else:
            raise RuntimeError(op2.code_information())

        if op2._results.is_not_saved(result_name):
            return ndata
        op2._results._found_result(result_name)

        storage_obj = op2.get_result(result_name)
        n = op2._read_random_table(data, ndata, result_name, storage_obj,
                                    obj, 'node',
                                    random_code=op2.random_code)

        assert n is not None, n
        return n

    def _read_oqg_mpc_ato(self, data: bytes, ndata: int) -> int:
        """
        table_code = ???
        """
        op2 = self.op2
        if op2.thermal == 0:
            if op2.table_code in [39]:
                assert op2.table_name in [b'OQMATO1', b'OQMATO2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'ato.mpc_forces'
                obj = RealMPCForcesArray
            #elif op2.table_code in [603]:
                #assert op2.table_name in [b''], 'op2.table_name=%r' % op2.table_name
                #result_name = 'psd.mpc_forces'
            else:
                print(op2.table_code)
                raise RuntimeError(op2.code_information())
        else:
            raise RuntimeError(op2.code_information())

        if op2._results.is_not_saved(result_name):
            return ndata
        op2._results._found_result(result_name)

        storage_obj = op2.get_result(result_name)
        n = op2._read_random_table(data, ndata, result_name, storage_obj,
                                   obj, 'node',
                                   random_code=op2.random_code)

        assert n is not None, n
        return n

    def _read_oqg_mpc_crm(self, data: bytes, ndata: int) -> int:
        """
        table_code = ???
        """
        op2 = self.op2
        if op2.thermal == 0:
            if op2.table_code in [39]:
                assert op2.table_name in [b'OQMCRM1', b'OQMCRM2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'crm.mpc_forces'
                obj = RealMPCForcesArray
            #elif op2.table_code in [603]:
                #assert op2.table_name in [b''], 'op2.table_name=%r' % op2.table_name
                #result_name = 'psd.mpc_forces'
            else:
                print(op2.table_code)
                raise RuntimeError(op2.code_information())
        else:
            raise RuntimeError(op2.code_information())

        if op2._results.is_not_saved(result_name):
            return ndata
        op2._results._found_result(result_name)

        storage_obj = op2.get_result(result_name)
        n = op2._read_random_table(data, ndata, result_name, storage_obj,
                                   obj, 'node',
                                   random_code=op2.random_code)
        assert n is not None, n
        return n

    def _read_oqg_mpc_rms(self, data: bytes, ndata: int) -> int:
        """
        table_code = ???
        """
        op2 = self.op2
        if op2.thermal == 0:
            if op2.table_code in [39]:
                assert op2.table_name in [b'OQMRMS1', b'OQMRMS2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'rms.mpc_forces'
                obj = RealMPCForcesArray
            #elif op2.table_code in [603]:
                #assert op2.table_name in [b''], 'op2.table_name=%r' % op2.table_name
                #result_name = 'psd.mpc_forces'
            else:
                raise RuntimeError(op2.code_information())
        else:
            raise RuntimeError(op2.code_information())

        if op2._results.is_not_saved(result_name):
            return ndata
        op2._results._found_result(result_name)

        storage_obj = op2.get_result(result_name)
        n = op2._read_random_table(data, ndata, result_name, storage_obj,
                                   obj, 'node',
                                   random_code=op2.random_code)

        assert n is not None, n
        return n

    def _read_oqg_mpc_no(self, data: bytes, ndata: int) -> int:
        """
        table_code = ???
        """
        op2 = self.op2
        if op2.thermal == 0:
            if op2.table_code in [39]:
                assert op2.table_name in [b'OQMNO1', b'OQMNO2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'no.mpc_forces'
                obj = RealMPCForcesArray
            #elif op2.table_code in [603]:
                #assert op2.table_name in [b''], 'op2.table_name=%r' % op2.table_name
                #result_name = 'psd.mpc_forces'
            else:
                raise RuntimeError(op2.code_information())
        else:
            raise RuntimeError(op2.code_information())

        if op2._results.is_not_saved(result_name):
            return ndata
        op2._results._found_result(result_name)

        storage_obj = op2.get_result(result_name)
        n = op2._read_random_table(data, ndata, result_name, storage_obj,
                                   obj, 'node',
                                   random_code=op2.random_code)
        assert n is not None, n
        return n

    def _read_obc1_3(self, data: bytes, ndata: int):
        op2 = self.op2
        op2.to_nx(f' because table_name={op2.table_name} was found')
        op2.nonlinear_factor = np.nan #None
        op2.is_table_1 = True
        op2.is_table_2 = False
        unused_three = op2.parse_approach_code(data)
        op2.words = [
            'analysis_code', 'table_code', '???', 'isubcase',
            '???', '???', '???', 'random_code',
            'format_code', 'num_wide', '11', '12',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', 'Title', 'subtitle', 'label']

        ## random code
        op2.random_code = op2.add_data_parameter(data, 'random_code', b'i', 8, False)

        ## format code
        op2.format_code = op2.add_data_parameter(data, 'format_code', b'i', 9, False)

        ## number of words per entry in record
        op2.num_wide = op2.add_data_parameter(data, 'num_wide', b'i', 10, False)

        #12 PID I Physical property
        #28 BLTSEQID I Bolt sequence
        op2.pid = op2.add_data_parameter(data, 'pid', b'i', 12, False)
        op2.bolt_seq_id = op2.add_data_parameter(data, 'bolt_seq_id', b'i', 28, False)

        if not op2.is_sort1:
            raise NotImplementedError('SORT2; code_info=\n%s' % op2.code_information())
        #assert self.isThermal()==False,op2.thermal

        ## assuming tCode=1
        if op2.analysis_code == 1:   # statics / displacement / heat flux
            ## load set number
            op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
            op2.setNullNonlinearFactor()
        #elif op2.analysis_code == 2:  # real eigenvalues
            ### mode number
            #op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            ### eigenvalue
            #op2.eign = op2.add_data_parameter(data, 'eign', b'f', 6, False)
            ### mode or cycle .. todo:: confused on the type - F1???
            #op2.mode_cycle = op2.add_data_parameter(data, 'mode_cycle', b'f', 7, False)
            #op2.reader_oug.update_mode_cycle('mode_cycle')
            #op2.data_names = op2.apply_data_code_value('data_names',
                                                         #['mode', 'eign', 'mode_cycle'])
        #elif op2.analysis_code == 3: # differential stiffness
            #op2.lsdvmn = self.get_values(data, b'i', 5) ## load set number
            #op2.data_names = op2.data_code['lsdvmn'] = op2.lsdvmn
        #elif op2.analysis_code == 4: # differential stiffness
            #op2.lsdvmn = self.get_values(data, b'i', 5) ## load set number
        #elif op2.analysis_code == 5:   # frequency
            ### frequency
            #op2.freq = op2.add_data_parameter(data, 'freq', b'f', 5)
            #op2.data_names = op2.apply_data_code_value('data_names', ['freq'])
        elif op2.analysis_code == 6:  # transient
            ## time step
            op2.dt = op2.add_data_parameter(data, 'dt', b'f', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['dt'])
        #elif op2.analysis_code == 7:  # pre-buckling
            ### load set number
            #op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            #op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
        #elif op2.analysis_code == 8:  # post-buckling
            ### load set number
            #op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            ### real eigenvalue
            #op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            #op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn', 'eigr'])
        #elif op2.analysis_code == 9:  # complex eigenvalues
            ### mode number
            #op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            ### real eigenvalue
            #op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            ### imaginary eigenvalue
            #op2.eigi = op2.add_data_parameter(data, 'eigi', b'f', 7, False)
            #op2.data_names = op2.apply_data_code_value('data_names', ['mode', 'eigr', 'eigi'])
        elif op2.analysis_code == 10:  # nonlinear statics
            ## load step
            op2.lftsfq = op2.add_data_parameter(data, 'lftsfq', b'f', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['lftsfq'])
        #elif op2.analysis_code == 11:  # old geometric nonlinear statics
            ### load set number
            #op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            #op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
        #elif op2.analysis_code == 12:
            ## contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ### load set number
            #op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            #op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % op2.analysis_code
            raise RuntimeError(msg)

        #self.fix_format_code()
        #if op2.num_wide == 8:
            #op2.format_code = 1
            #op2.data_code['format_code'] = 1
        #else:
            ##self.fix_format_code()
            #if op2.format_code == 1:
                #op2.format_code = 2
                #op2.data_code['format_code'] = 2
            #assert op2.format_code in [2, 3], op2.code_information()

        #op2._parse_thermal_code()
        if op2.is_debug_file:
            op2.binary_debug.write('  approach_code  = %r\n' % op2.approach_code)
            op2.binary_debug.write('  tCode          = %r\n' % op2.tCode)
            op2.binary_debug.write('  isubcase       = %r\n' % op2.isubcase)
        op2._read_title(data)
        op2._write_debug_bits()

    def _read_obc1_4(self, data: bytes, ndata: int) -> int:
        """C O N T A C T  P R E S S U R E  A N D  T R A C T I O N S"""
        op2 = self.op2
        result_name = 'contact_tractions_and_pressure'
        #if op2._results.is_not_saved(result_name):
            #print('return A')
            #return ndata
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        ntotal = 20 * self.factor
        nnodes = ndata // ntotal

        op2.data_code['_times_dtype'] = 'float32'
        obj_vector_real = RealContactTractionAndPressureArray
        #auto_return, is_vectorized = self._create_oes_object4(
            #nnodes, result_name, slot, obj_vector_real)
        auto_return = op2._create_table_vector(
            result_name, nnodes, slot, obj_vector_real, is_cid=False)

        if auto_return:
            return nnodes * ntotal
        obj = op2.obj
        #ints    = (327, 1233577877, 0, 0, -988060411,
                   #357, 1233066655, 0, 0, -965162024,
                   #367, 1233710164, 0, 0, -993503801,
                   #987, 1233710164, 0, 0, 1153979847,
                   #997, 1233445590, 0, 0, 1163042599,
                   #1007, 1233066655, 0, 0, 1182321624)
        #floats  = (4.582245978342152e-43, 1105138.625, 0.0, 0.0, -2486.313720703125,
                   #5.002635517639597e-43, 1044905.9375, 0.0, 0.0, -15922.9609375,
                   #5.142765364072079e-43, 1121674.5, 0.0, 0.0, -1602.6805419921875,
                   #1.3830815842885944e-42, 1121674.5, 0.0, 0.0, 1602.6805419921875,
                   #1.3970945689318426e-42, 1088602.75, 0.0, 0.0, 3369.947021484375,
                   #1.4111075535750908e-42, 1044905.9375, 0.0, 0.0, 15922.9609375)
        dt = op2.nonlinear_factor
        from struct import Struct
        struct1 = Struct(mapfmt(op2._endian + b'i fff f', self.size))

        is_vectorized = False
        if op2.use_vector and is_vectorized:
            # POINT ID. TYPE  PRESSURE  S1  S2  S3
            itime = obj.itime
            n = nnodes * ntotal
            itotal = obj.itotal
            itotal2 = itotal + nnodes

            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nnodes, 4)

                nids = ints[:, 0] // 10
                assert nids.min() > 0, nids.min()
                obj.node_gridtype[itotal:itotal2, 0] = nids
                #obj.node_gridtype[itotal:itotal2, 1] = ints[:, 1].copy()

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nnodes, 4)
            obj.data[obj.itime, obj.itotal:itotal2, :] = floats[:, 2:].copy()
            obj._times[itime] = dt
            obj.itotal = itotal2
        else:
            n = 0
            for unused_i in range(nnodes):
                edata = data[n:n+ntotal]
                nid_device, pressure, s1, s2, s3 = struct1.unpack(edata)
                unused_nid = nid_device // 10
                #out2 = [nid, pressure, s1, s2, s3]
                #obj.add_sort1(nid, pressure, s1, s2, s3)
                #print(out2)
                n += ntotal
            #self.show_data(data)

        #if self._table4_count == 0:
            #self.log.warning('no output for "Contact Pressure and Tractions"')
        #self._table4_count += 1
        return ndata
