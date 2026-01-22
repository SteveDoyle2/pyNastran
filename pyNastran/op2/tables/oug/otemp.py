#pylint: disable=C0301,C0103
"""
This file defines the OTEMP1 Table, which contains:
 * Real Temperature
   - DISPLACEMENT = ALL
"""
from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.op2.tables.utils import get_is_slot_saved
from pyNastran.op2.tables.oug.oug_temperatures import RealTemperatureArray

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


class OTEMP:
    """
    OTEMP1: New version of NX tempererature output for SOL 401?
    """
    def __init__(self, op2: OP2):
        self.op2 = op2

    def _read_otemp1_3(self, data: bytes, ndata: int):
        """SOL 401 table"""
        op2 = self.op2
        op2.nonlinear_factor = np.nan
        op2.is_table_1 = True
        op2.is_table_2 = False
        op2.parse_approach_code(data)  # field 3
        op2.words = [
            'approach_code', 'table_code', '???', 'isubcase',
            '???', '???', '???', '???',
            'format_code', 'num_wide', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', 'thermal', '???',
            '???', 'Title', 'subtitle', 'label']

        ## format code
        op2.format_code = op2.add_data_parameter(data, 'format_code', b'i', 9, False)

        ## number of words per entry in record
        op2.num_wide = op2.add_data_parameter(data, 'num_wide', b'i', 10, False)

        ## nBolt sequence number for SOL 401 preloaded bolts
        self.bolt_seq_id = op2.add_data_parameter(data, 'bolt_seq_id', b'i', 28, False)

        if op2.analysis_code == 6:  # transient
            # time step
            op2.dt = op2.add_data_parameter(data, 'dt', b'f', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['dt'])
        elif op2.analysis_code == 10:  # nonlinear statics
            # load step
            op2.lftsfq = op2.add_data_parameter(data, 'lftsfq', b'f', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['lftsfq'])
        else:  # pragma: no cover
            msg = f'invalid analysis_code...analysis_code={op2.analysis_code}\ndata={op2.data_code}'
            raise RuntimeError(msg)

        if op2.is_debug_file:
            op2.binary_debug.write('  approach_code  = %r\n' % op2.approach_code)
            op2.binary_debug.write('  tCode          = %r\n' % op2.tCode)
            op2.binary_debug.write('  isubcase       = %r\n' % op2.isubcase)
        op2._read_title(data)
        op2._write_debug_bits()
        #print(op2.code_information())

    def _read_otemp1_4(self, data: bytes, ndata: int):
        """SOL 401 table"""
        op2 = self.op2
        nfields = ndata // 4
        nnodes = nfields // 2
        result_name = 'temperatures'
        real_vector = RealTemperatureArray
        is_cid = False
        op2.data_code['_times_dtype'] = 'float32'
        #self._times_dtype = 'float32'
        is_saved, slot = get_is_slot_saved(op2, result_name)
        if not is_saved:
            return ndata
        auto_return = op2._create_table_vector(
            result_name, nnodes, slot, real_vector, is_cid=is_cid)
        if auto_return:
            return ndata

        #print(op2.obj)
        #print(op2.code_information())
        floats = np.frombuffer(data, dtype=op2.fdtype).reshape(nnodes, 2).copy()
        ints = np.frombuffer(data, dtype=op2.idtype).reshape(nnodes, 2) // 10
        #print(op2.obj.get_stats())
        nids = ints[:, 0]
        temps = floats[:, 1]
        op2.obj.node_gridtype[:, 0] = nids
        op2.obj.data[op2.obj.itime, :, 0] = temps
        return ndata
