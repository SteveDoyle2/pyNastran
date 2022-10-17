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
from __future__ import annotations
#from struct import Struct
from typing import TYPE_CHECKING
import numpy as np
#from pyNastran import DEV
#from pyNastran.utils.numpy_utils import integer_types
#from pyNastran.op2.op2_interface.op2_reader import mapfmt

from pyNastran.op2.tables.oug.oug_displacements import RealDisplacementArray
from pyNastran.op2.tables.oug.oug_velocities import RealVelocityArray
from pyNastran.op2.tables.oug.oug_accelerations import RealAccelerationArray

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


class OUGPK:
    """
    OUGPK1 : Peak to Peak Displacement?

    Output U in the global frame

    U is:
     - Displacement
     - Velocity
     - Accelerations

    The global frame is:
     - the analysis coordinate frame, not the 0 coordinate frame
     """
    def __init__(self, op2: OP2):
        self.op2 = op2

    @property
    def size(self) -> int:
        return self.op2.size
    @property
    def factor(self) -> int:
        return self.op2.factor

    def _read_ougpk1_3(self, data: bytes, ndata: int):
        """reads table 3 (the header table)"""
        op2 = self.op2
        assert ndata == 146 * op2.size
        #self._set_times_dtype()
        op2.nonlinear_factor = np.nan
        op2.is_table_1 = True
        op2.is_table_2 = False
        unused_three = op2.parse_approach_code(data)
        op2.words = [
            'approach_code', 'table_code', '???', 'isubcase',
            '???', '???', '???', 'random_code',
            'format_code', 'num_wide', '???', '???',
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
        op2.acoustic_flag = op2.add_data_parameter(data, 'acoustic_flag', b'i', 13, False)

        ## thermal flag; 1 for heat transfer, 0 otherwise
        op2.thermal = op2.add_data_parameter(data, 'thermal', b'i', 23, False)

        if op2.analysis_code == 5:   # frequency
            # frequency
            op2.freq = op2.add_data_parameter(data, 'freq', b'f', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['freq'])
        else:  # pragma: no cover
            msg = f'invalid analysis_code...analysis_code={op2.analysis_code}\ndata={op2.data_code}'
            raise RuntimeError(msg)

        #print(op2.code_information())
        op2._fix_oug_format_code()
        op2._parse_thermal_code()
        if op2.is_debug_file:
            op2.binary_debug.write('  approach_code  = %r\n' % op2.approach_code)
            op2.binary_debug.write('  tCode          = %r\n' % op2.tCode)
            op2.binary_debug.write('  isubcase       = %r\n' % op2.isubcase)
        op2._read_title(data)
        op2._write_debug_bits()
        #op2.show_data(data[:200], types='if')
        #self._correct_eigenvalue()

    #def _correct_eigenvalue(self):
        #"""Nastran 95 gets the frequency wrong"""
        #op2 = self.op2
        #if op2._nastran_format == 'nasa95' and op2.analysis_code == 2:  # real eigenvalues
            ##print(op2.mode, op2.eign, op2.mode_cycle)
            ## sqrt(lambda) = omega = 2*pi*f
            #freq = (op2.eign) ** 0.5 / (2 * np.pi)
            #op2.mode_cycle = freq
            #op2.data_code['mode_cycle'] = freq

    def _read_ougpk1_4(self, data: bytes, ndata: int):
        """reads the SORT1 version of table 4 (the data table)"""
        op2 = self.op2
        #table_name_bytes = op2.table_name

        #is_vectorized = True
        if op2.table_code == 401:
            # Peak Displacement Vector
            result_name = 'rms.displacements'
            obj = RealDisplacementArray
        elif op2.table_code == 410:
            # Peak Velocity Vector
            result_name = 'rms.velocities'
            obj = RealVelocityArray
        elif op2.table_code == 411:
            # Peak Acceleration Vector
            obj = RealAccelerationArray
            result_name = 'rms.accelerations'
        else:
            raise RuntimeError(op2.table_code)

        if op2._results.is_not_saved(result_name):
            return ndata
        op2._results._found_result(result_name)
        storage_obj = op2.get_result(result_name)

        assert op2.format_code == 1, op2.format_code
        assert op2.num_wide == 8, op2.num_wide

        op2.random_code
        n = op2._read_random_table(data, ndata, result_name, storage_obj,
                                    obj, 'node',
                                    random_code=op2.random_code)
        return n
