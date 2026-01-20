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
from typing import TYPE_CHECKING
import numpy as np
from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.tables.utils import get_is_slot_saved

from pyNastran.op2.tables.oqg_constraintForces.separation_distance import (
    SeparationDistanceArray)
#from pyNastran.op2.tables.oqg_constraintForces.oqg_contact_forces import RealContactForcesArray
from pyNastran.op2.result_objects.contact_traction_and_pressure import RealContactTractionAndPressureArray
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2

class OBC:
    def __init__(self, op2: OP2):
        self.op2 = op2

    def read_sort1_3(self, data: bytes, ndata: int):
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

    def read_4(self, data: bytes, ndata: int) -> int:
        """C O N T A C T  P R E S S U R E  A N D  T R A C T I O N S"""
        op2 = self.op2
        result_name = 'contact_tractions_and_pressure'
        is_saved, slot = get_is_slot_saved(op2, result_name)
        if not is_saved:
            return ndata

        ntotal = 20 * op2.factor
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

        is_vectorized = True
        if op2.use_vector and is_vectorized:
            # POINT ID. TYPE  PRESSURE  S1  S2  S3
            itime = obj.itime
            n = nnodes * ntotal
            itotal = obj.itotal
            itotal2 = itotal + nnodes

            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nnodes, 5)

                nids = ints[:, 0] // 10
                assert nids.min() > 0, nids.min()
                obj.node_gridtype[itotal:itotal2, 0] = nids
                #obj.node_gridtype[itotal:itotal2, 1] = ints[:, 1].copy()

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nnodes, 5)
            obj.data[obj.itime, obj.itotal:itotal2, :] = floats[:, 1:].copy()
            obj._times[itime] = dt
            obj.itotal = itotal2
        else:
            n = obc_real_5(op2, obj, data, ntotal, nnodes, dt)

        #if self._table4_count == 0:
        #self._table4_count += 1
        assert n == ndata
        return ndata


def obc_real_5(op2: OP2,
               obj: RealContactTractionAndPressureArray,
               data: bytes, ntotal: int, nnodes: int, dt) -> int:
    n = 0
    struct1 = Struct(mapfmt(op2._endian + b'i fff f', op2.size))
    for unused_i in range(nnodes):
        edata = data[n:n+ntotal]
        nid_device, pressure, s1, s2, s3 = struct1.unpack(edata)
        nid = nid_device // 10
        #out2 = [nid, pressure, s1, s2, s3]
        obj.add_sort1(dt, nid, pressure, s1, s2, s3)
        #print(nid, pressure, s1, s2, s3)
        n += ntotal
    return n
