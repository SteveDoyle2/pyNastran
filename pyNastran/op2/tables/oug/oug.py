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
from struct import Struct
from typing import TYPE_CHECKING
import numpy as np
from pyNastran import DEV
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.op2_interface.op2_reader import mapfmt

from pyNastran.op2.tables.oug.oug_displacements import (
    RealDisplacementArray, ComplexDisplacementArray)

from pyNastran.op2.tables.oug.oug_velocities import (
    RealVelocityArray, ComplexVelocityArray)

from pyNastran.op2.tables.oug.oug_accelerations import (
    RealAccelerationArray, ComplexAccelerationArray)

from pyNastran.op2.tables.oug.oug_temperatures import (
    RealTemperatureArray)

from pyNastran.op2.tables.oug.oug_eigenvectors import (
    RealEigenvectorArray, ComplexEigenvectorArray,
)

from pyNastran.op2.tables.opg_appliedLoads.opg_load_vector import RealThermalVelocityVectorArray
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


class OUG:
    """
    OUG : Output U in the global frame

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

    def update_mode_cycle(self, name):
        op2 = self.op2
        value = getattr(op2, name)
        if value == 0.0:
            #print('table_name=%r mode=%s eigr=%s' % (op2.table_name, op2.mode, op2.eigr))
            value = np.sqrt(np.abs(op2.eign)) / (2. * np.pi)
            setattr(op2, name, value)
            op2.data_code[name] = value

    def _read_otemp1_3(self, data: bytes, ndata: int):
        """SOL 401 table"""
        op2 = self.op2
        op2.nonlinear_factor = np.nan
        op2.is_table_1 = True
        op2.is_table_2 = False
        unused_three = op2.parse_approach_code(data)
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
        storage_obj = op2.temperatures
        real_vector = RealTemperatureArray
        is_cid = False
        op2.data_code['_times_dtype'] = 'float32'
        #self._times_dtype = 'float32'
        auto_return = op2._create_table_vector(
            result_name, nnodes, storage_obj, real_vector, is_cid=is_cid)
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

    def _read_oug1_3(self, data: bytes, ndata: int):
        """reads table 3 (the header table)"""
        op2 = self.op2
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

        if op2.analysis_code == 1:   # statics / displacement / heat flux
            # load set number
            op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
            op2.setNullNonlinearFactor()
        elif op2.analysis_code == 2:  # real eigenvalues
            # mode number
            op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            # eigenvalue
            op2.eign = op2.add_data_parameter(data, 'eign', b'f', 6, False)
            # mode or cycle .. todo:: confused on the type - F1???
            # float - C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\mftank.op2
            #op2.mode_cycle = op2.add_data_parameter(data, 'mode_cycle', b'i', 7, False)  # nope...
            op2.mode_cycle = op2.add_data_parameter(data, 'mode_cycle', b'f', 7, False) # radians
            self.update_mode_cycle('mode_cycle')
            op2.data_names = op2.apply_data_code_value('data_names', ['mode', 'eign', 'mode_cycle'])
        #elif op2.analysis_code == 3: # differential stiffness
            #op2.lsdvmn = self.get_values(data, b'i', 5) ## load set number
            #op2.data_code['lsdvmn'] = op2.lsdvmn
        #elif op2.analysis_code == 4: # differential stiffness
            #op2.lsdvmn = self.get_values(data, b'i', 5) ## load set number
        elif op2.analysis_code == 5:   # frequency
            # frequency
            op2.freq = op2.add_data_parameter(data, 'freq', b'f', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['freq'])
        elif op2.analysis_code == 6:  # transient
            # time step
            op2.dt = op2.add_data_parameter(data, 'dt', b'f', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['dt'])
        elif op2.analysis_code == 7:  # pre-buckling
            # load set number
            op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
        elif op2.analysis_code == 8:  # post-buckling
            # load set number
            op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            # real eigenvalue
            op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn', 'eigr'])
        elif op2.analysis_code == 9:  # complex eigenvalues
            # mode number
            op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            # real eigenvalue
            op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            # imaginary eigenvalue
            op2.eigi = op2.add_data_parameter(data, 'eigi', b'f', 7, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['mode', 'eigr', 'eigi'])
        elif op2.analysis_code == 10:  # nonlinear statics
            # load step
            op2.lftsfq = op2.add_data_parameter(data, 'lftsfq', b'f', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['lftsfq'])
        elif op2.analysis_code == 11:  # old geometric nonlinear statics
            # load set number
            op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
        elif op2.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            # load set number
            op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
        else:  # pragma: no cover
            msg = f'invalid analysis_code...analysis_code={op2.analysis_code}\ndata={op2.data_code}'
            raise RuntimeError(msg)

        #print op2.code_information()
        op2._fix_oug_format_code()
        op2._parse_thermal_code()
        if op2.is_debug_file:
            op2.binary_debug.write('  approach_code  = %r\n' % op2.approach_code)
            op2.binary_debug.write('  tCode          = %r\n' % op2.tCode)
            op2.binary_debug.write('  isubcase       = %r\n' % op2.isubcase)
        op2._read_title(data)
        op2._write_debug_bits()
        self._correct_eigenvalue()

    def _correct_eigenvalue(self):
        """Nastran 95 gets the frequency wrong"""
        op2 = self.op2
        if op2._nastran_format == 'nasa95' and op2.analysis_code == 2:  # real eigenvalues
            #print(op2.mode, op2.eign, op2.mode_cycle)
            # sqrt(lambda) = omega = 2*pi*f
            freq = (op2.eign) ** 0.5 / (2 * np.pi)
            op2.mode_cycle = freq
            op2.data_code['mode_cycle'] = freq

    def _read_oug2_3(self, data: bytes, ndata: int) -> None:
        """reads the SORT2 version of table 4 (the data table)"""
        op2 = self.op2
        #self._set_times_dtype()
        #return self._read_oug1_3(data)
        op2.nonlinear_factor = np.nan

        op2.is_table_1 = False
        op2.is_table_2 = True
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

        op2.node_id = op2.add_data_parameter(data, 'node_id', b'i', 5, fix_device_code=True)
        #if op2.analysis_code == 1:  # statics / displacement / heat flux
            ## load set number
            #op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5, False)
            #op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            #op2.setNullNonlinearFactor()

        if op2.analysis_code == 1:  # static...because reasons.
            op2._analysis_code_fmt = b'f'
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            op2.apply_data_code_value('analysis_method', 'N/A')
        elif op2.analysis_code == 2:  # real eigenvalues
            ## mode number
            #op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            op2._analysis_code_fmt = b'i'
            # real eigenvalue
            op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            ## mode or cycle .. todo:: confused on the type - F1???
            # float - C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\mftank.op2
             # mode or cycle .. todo:: confused on the type - F1???
            op2.mode_cycle = op2.add_data_parameter(data, 'mode_cycle', b'i', 7, False)
            #op2.mode_cycle = op2.add_data_parameter(data, 'mode_cycle', b'f', 7, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id', 'eigr', 'mode_cycle'])
            op2.apply_data_code_value('analysis_method', 'mode')
        #elif op2.analysis_code == 3: # differential stiffness
            #op2.lsdvmn = self.get_values(data, b'i', 5) ## load set number
            #op2.data_code['lsdvmn'] = op2.lsdvmn
        #elif op2.analysis_code == 4: # differential stiffness
            #op2.lsdvmn = self.get_values(data, b'i', 5) ## load set number
        elif op2.analysis_code == 5:   # frequency
            # frequency
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
            op2._analysis_code_fmt = b'f'
            ## real eigenvalue
            op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id', 'eigr'])
            op2.apply_data_code_value('analysis_method', 'eigr')
        elif op2.analysis_code == 9:  # complex eigenvalues
            ## mode number
            #op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            op2._analysis_code_fmt = b'i'
            ## real eigenvalue
            op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
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
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
        elif op2.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## load set number
            #op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            op2.apply_data_code_value('analysis_method', 'lsdvmn')
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % op2.analysis_code
            raise RuntimeError(msg)

        op2._fix_oug_format_code()
        op2._parse_thermal_code()
        if op2.is_debug_file:
            op2.binary_debug.write('  %-14s = %r %s\n' % ('approach_code', op2.approach_code,
                                                           op2.approach_code_str(op2.approach_code)))
            op2.binary_debug.write('  %-14s = %r\n' % ('tCode', op2.tCode))
            op2.binary_debug.write('  %-14s = %r\n' % ('isubcase', op2.isubcase))
        op2._read_title(data)
        op2._write_debug_bits()
        assert isinstance(op2.nonlinear_factor, integer_types), op2.nonlinear_factor

    def _read_ougpc1_3(self, data: bytes, ndata: int):
        """reads table 3 (the header table)"""
        op2 = self.op2
        op2.to_nx(f' because table_name={op2.table_name} was found')
        #self._set_times_dtype()
        op2.nonlinear_factor = np.nan
        op2.is_table_1 = True
        op2.is_table_2 = False
        unused_three = op2.parse_approach_code(data)
        op2.words = [
            'approach_code', 'table_code', '???', 'isubcase',
            '???', '???', '???', '???',
            'format_code', 'num_wide', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', 'thermal', '???',
            '???', 'Title', 'subtitle', 'label']

        ## pcode
        # Panel contribution code: +/-1=abs, +/-2=norm

        ## data_type
        #Acoustic dof code (10*grid ID + direction)
        #Direction has the following values:
        #=0, Pressure
        #=1, X-displacement
        #=2, Y-displacement
        #=3, Z-displacement
        #=4, RX-displacement
        #=5, RY-displacement
        #=6, RZ-displacement
        op2.dcode = op2.add_data_parameter(data, 'dcode', b'i', 5, False)

        ## Panel name (0 for TOTAL)
        panel_name1 = op2.add_data_parameter(data, 'data_type', b'4s', 6, False,
                                             add_to_dict=False)
        panel_name2 = op2.add_data_parameter(data, 'panel_name', b'4s', 7, False,
                                             add_to_dict=False)
        self.panel_name = panel_name1 + panel_name2
        op2.data_code['panel_name'] = self.panel_name

        ## data_type
        ## (1=pressure, 2=first derivative, 3=second derivative)
        op2.data_type = op2.add_data_parameter(data, 'data_type', b'i', 8, False)

        ## format code
        op2.format_code = op2.add_data_parameter(data, 'format_code', b'i', 9, False)

        ## number of words per entry in record
        op2.num_wide = op2.add_data_parameter(data, 'num_wide', b'i', 10, False)


        # 1 ACODE(C) I Device code + 10*Approach Code
        # 2 TCODE(C) I Table Code
        # 3 PCODE I Panel contribution code: +/-1=abs, +/-2=norm
        # 4 SUBCASE I Subcase number
        # 5 DCODE I Acoustic dof code (10*grid ID + direction)
        # Direction has the following values:
        # =0, Pressure
        # =1, X-displacement
        # =2, Y-displacement
        # =3, Z-displacement
        # =4, RX-displacement
        # =5, RY-displacement
        # =6, RZ-displacement
        # TCODE,1=01 Sort 1
        #   ACODE,4=05 Frequency
        #   6 FREQ RS Frequency (Hz)
        #   End ACODE,4
        # TCODE,1=02 Sort 2
        # 6 PNAME(2) CHAR4 Panel name (0 for TOTAL)
        # End TCODE,1
        # 8 DATTYP I Data Type (1=pressure, 2=first derivative, 3=second derivative)
        if op2.analysis_code == 1 and 0:   # statics / displacement / heat flux
            # load set number
            op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
            op2.setNullNonlinearFactor()
        #elif op2.analysis_code == 2:  # real eigenvalues
            ## mode number
            #op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            ## eigenvalue
            #op2.eign = op2.add_data_parameter(data, 'eign', b'f', 6, False)
            ## mode or cycle .. todo:: confused on the type - F1???
            ## float - C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\mftank.op2
            ##op2.mode_cycle = op2.add_data_parameter(data, 'mode_cycle', b'i', 7, False)  # nope...
            #op2.mode_cycle = op2.add_data_parameter(data, 'mode_cycle', b'f', 7, False) # radians
            #self.update_mode_cycle('mode_cycle')
            #op2.data_names = op2.apply_data_code_value('data_names', ['mode', 'eign', 'mode_cycle'])
        #elif op2.analysis_code == 3: # differential stiffness
            #op2.lsdvmn = self.get_values(data, b'i', 5) ## load set number
            #op2.data_code['lsdvmn'] = op2.lsdvmn
        #elif op2.analysis_code == 4: # differential stiffness
            #op2.lsdvmn = self.get_values(data, b'i', 5) ## load set number
        elif op2.analysis_code == 5:   # frequency
            # frequency
            op2.freq = op2.add_data_parameter(data, 'freq', b'i', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['freq'])
        #elif op2.analysis_code == 6:  # transient
            ## time step
            #op2.dt = op2.add_data_parameter(data, 'dt', b'f', 5)
            #op2.data_names = op2.apply_data_code_value('data_names', ['dt'])
        #elif op2.analysis_code == 7:  # pre-buckling
            ## load set number
            #op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            #op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
        #elif op2.analysis_code == 8:  # post-buckling
            ## load set number
            #op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            ## real eigenvalue
            #op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            #op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn', 'eigr'])
        #elif op2.analysis_code == 9:  # complex eigenvalues
            ## mode number
            #op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            ## real eigenvalue
            #op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            ## imaginary eigenvalue
            #op2.eigi = op2.add_data_parameter(data, 'eigi', b'f', 7, False)
            #op2.data_names = op2.apply_data_code_value('data_names', ['mode', 'eigr', 'eigi'])
        #elif op2.analysis_code == 10:  # nonlinear statics
            ## load step
            #op2.lftsfq = op2.add_data_parameter(data, 'lftsfq', b'f', 5)
            #op2.data_names = op2.apply_data_code_value('data_names', ['lftsfq'])
        #elif op2.analysis_code == 11:  # old geometric nonlinear statics
            ## load set number
            #op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            #op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
        #elif op2.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## load set number
            #op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            #op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
        else:  # pragma: no cover
            msg = f'invalid analysis_code...analysis_code={op2.analysis_code}\ndata={op2.data_code}'
            raise RuntimeError(msg)

        #print op2.code_information()
        op2._fix_oug_format_code()
        if op2.is_debug_file:
            op2.binary_debug.write('  approach_code  = %r\n' % op2.approach_code)
            op2.binary_debug.write('  tCode          = %r\n' % op2.tCode)
            op2.binary_debug.write('  isubcase       = %r\n' % op2.isubcase)
        op2._read_title(data)
        op2._write_debug_bits()
        self.warn_skip_table()

    def warn_skip_table(self):
        op2 = self.op2
        if op2._table4_count == 0:
            op2.log.warning(f'skipping {op2.table_name}')
            op2._table4_count += 1

    def _read_ougpc2_3(self, data: bytes, ndata: int):
        """reads the SORT2 version of table 4 (the data table)"""
        op2 = self.op2
        #self._set_times_dtype()
        op2.nonlinear_factor = np.nan

        op2.is_table_1 = False
        op2.is_table_2 = True
        unused_three = op2.parse_approach_code(data)
        op2.words = [
            'approach_code', 'table_code', '???', 'isubcase',
            '???', '???', '???', '???',
            'format_code', 'num_wide', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', 'Title', 'subtitle', 'label']

        ## data_type
        ## (1=pressure, 2=first derivative, 3=second derivative)
        op2.data_type = op2.add_data_parameter(data, 'data_type', b'i', 8, False)

        ## format code
        op2.format_code = op2.add_data_parameter(data, 'format_code', b'i', 9, False)

        ## number of words per entry in record
        op2.num_wide = op2.add_data_parameter(data, 'num_wide', b'i', 10, False)

        op2.node_id = op2.add_data_parameter(data, 'node_id', b'i', 5, fix_device_code=True)
        #if op2.analysis_code == 1:  # statics / displacement / heat flux
            ## load set number
            #op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5, False)
            #op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            #op2.setNullNonlinearFactor()

        #if op2.analysis_code == 1:  # static...because reasons.
            #op2._analysis_code_fmt = b'f'
            #op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            #op2.apply_data_code_value('analysis_method', 'N/A')
        #elif op2.analysis_code == 2:  # real eigenvalues
            ### mode number
            ##op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            #op2._analysis_code_fmt = b'i'
            ## real eigenvalue
            #op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            ### mode or cycle .. todo:: confused on the type - F1???
            ## float - C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\mftank.op2
             ## mode or cycle .. todo:: confused on the type - F1???
            #op2.mode_cycle = op2.add_data_parameter(data, 'mode_cycle', b'i', 7, False)
            ##op2.mode_cycle = op2.add_data_parameter(data, 'mode_cycle', b'f', 7, False)
            #op2.data_names = op2.apply_data_code_value('data_names', ['node_id', 'eigr', 'mode_cycle'])
            #op2.apply_data_code_value('analysis_method', 'mode')
        #elif op2.analysis_code == 3: # differential stiffness
            #op2.lsdvmn = self.get_values(data, b'i', 5) ## load set number
            #op2.data_code['lsdvmn'] = op2.lsdvmn
        #elif op2.analysis_code == 4: # differential stiffness
            #op2.lsdvmn = self.get_values(data, b'i', 5) ## load set number

        if op2.analysis_code == 5:   # frequency
            # frequency
            ## Panel name (0 for TOTAL)
            panel_name1 = op2.add_data_parameter(data, 'panel_name1', b'4s', 6, False,
                                                  add_to_dict=False)
            panel_name2 = op2.add_data_parameter(data, 'panel_name2', b'4s', 7, False,
                                                  add_to_dict=False)
            self.panel_name = panel_name1 + panel_name2
            op2.data_code['panel_name'] = self.panel_name
            #print(self.panel_name)


            #op2.freq = op2.add_data_parameter(data, 'freq', b'f', 5)
            op2._analysis_code_fmt = b'f'
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            op2.apply_data_code_value('analysis_method', 'freq')
        #elif op2.analysis_code == 6:  # transient
            ### time step
            ##op2.dt = op2.add_data_parameter(data, 'dt', b'f', 5)
            #op2._analysis_code_fmt = b'f'
            #op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            #op2.apply_data_code_value('analysis_method', 'dt')
        #elif op2.analysis_code == 7:  # pre-buckling
            ### load set number
            ##op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            #op2._analysis_code_fmt = b'i'
            #op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            #op2.apply_data_code_value('analysis_method', 'lsdvmn')
        #elif op2.analysis_code == 8:  # post-buckling
            ### load set number
            ##op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            #op2._analysis_code_fmt = b'f'
            ### real eigenvalue
            #op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            #op2.data_names = op2.apply_data_code_value('data_names', ['node_id', 'eigr'])
            #op2.apply_data_code_value('analysis_method', 'eigr')
        #elif op2.analysis_code == 9:  # complex eigenvalues
            ### mode number
            ##op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            #op2._analysis_code_fmt = b'i'
            ### real eigenvalue
            #op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            ### imaginary eigenvalue
            #op2.eigi = op2.add_data_parameter(data, 'eigi', b'f', 7, False)
            #op2.data_names = op2.apply_data_code_value('data_names', ['node_id', 'eigr', 'eigi'])
            #op2.apply_data_code_value('analysis_method', 'mode')
        #elif op2.analysis_code == 10:  # nonlinear statics
            ### load step
            ##op2.lftsfq = op2.add_data_parameter(data, 'lftsfq', b'f', 5)
            #op2._analysis_code_fmt = b'f'
            #op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            #op2.apply_data_code_value('analysis_method', 'lftsfq')
        #elif op2.analysis_code == 11:  # old geometric nonlinear statics
            ### load set number
            ##op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            #op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
        #elif op2.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ### load set number
            ##op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            #op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            #op2.apply_data_code_value('analysis_method', 'lsdvmn')
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % op2.analysis_code
            raise RuntimeError(msg)

        op2._fix_oug_format_code()
        if op2.is_debug_file:
            op2.binary_debug.write('  %-14s = %r %s\n' % ('approach_code', op2.approach_code,
                                                           op2.approach_code_str(op2.approach_code)))
            op2.binary_debug.write('  %-14s = %r\n' % ('tCode', op2.tCode))
            op2.binary_debug.write('  %-14s = %r\n' % ('isubcase', op2.isubcase))
        op2._read_title(data)
        op2._write_debug_bits()
        assert isinstance(op2.nonlinear_factor, integer_types), op2.nonlinear_factor
        self.warn_skip_table()

    def _read_ougpc_4(self, data: bytes, ndata: int):
        """reads table 4 (the results table)"""
        op2 = self.op2
        assert op2.table_code == 49, op2.code_information()
        if op2.read_mode == 1:
            return ndata
        #self.show_data(data)
        #print(f'data_type = {op2.data_type}')
        assert op2.data_type == 1, op2.data_type

        n = 0
        ntotal = 16 * self.factor
        npanels = ndata // ntotal
        assert ndata % ntotal == 0, f'ndata={ndata} ntotal={ntotal} data_type={op2.data_type}'
        if op2.sort_method == 1:
            struct1 = Struct(op2._endian + b'8s ff')
            for unused_i in range(npanels):
                # Panel name (0 for TOTAL)
                edata = data[n:n+ntotal]
                name, real, imag = struct1.unpack(edata)
                #print(name, real, imag)
                n += ntotal
        else:
            struct1 = Struct(op2._endian + b'ff ff')
            for unused_i in range(npanels):
                # Panel name (0 for TOTAL)
                edata = data[n:n+ntotal]
                #self.show_data(edata)
                freq, null, real, imag = struct1.unpack(edata)
                #print(freq, null, real, imag)
                n += ntotal
        #strings = (b'0       \xed\xf0{AE\x9d\\APANEL3  O\xa1\x82A\x81\xc5eAPANEL5  a\xeb"\xbe\xbf/\xa4\xbePANEL6  \xda\xe2\xa1\xbe\x0c\x1a!\xbcPANEL2  \xe3~\r\xbetE\xab\xbdPANEL1  ep\x1a\xbd\x05@\xb8\xbd',)
        #ints    = (538976304, 538976288, 1098641645, 1096588613, 1162756432, 538981196, 1099080015, 1097188737, 1162756432, 538981708, -1105007775, -1096536129, 1162756432, 538981964, -1096686886, -1138681332, 1162756432, 538980940, -1106411805, -1112849036, 1162756432, 538980684, -1122340763, -1111998459)
        #floats  = (1.3563177106455426e-19, 1.3563156426940112e-19, 15.746319770812988, 13.788395881652832, 3300.08203125, 1.3569499868262628e-19, 16.328763961791992, 14.360718727111816, 3300.08203125, 1.357016161275267e-19, -0.15910102427005768, -0.3206767737865448, 3300.08203125, 1.3570492484997691e-19, -0.316183865070343, -0.009832870215177536, 3300.08203125, 1.3569168996017607e-19, -0.13817934691905975, -0.0836285650730133, 3300.08203125, 1.3568838123772585e-19, -0.037704844027757645, -0.08996585756540298)
        return ndata

    def _read_ougmc_4(self, data: bytes, ndata: int):
        op2 = self.op2
        if op2.table_code == 44:   # Displacements
            if op2.table_name in [b'OUGMC1', b'OUGMC2']:
                assert op2.thermal == 0, op2.code_information()
                result_name = 'modal_contribution.displacements'
            else:
                raise NotImplementedError(op2.code_information())
        elif op2.table_code == 48:   # spc_forces
            if op2.table_name in [b'OQGMC1', b'OQGMC2']:
                assert op2.thermal == 0, op2.code_information()
                result_name = 'modal_contribution.spc_forces'
            else:
                raise NotImplementedError(op2.code_information())
        else:
            raise NotImplementedError(op2.code_information())

        n = 0
        if op2.table_name in [b'OUGMC1', b'OQGMC1']:
            if op2.read_mode == 1:
                return ndata
            from struct import Struct
            ntotal = 16 * self.factor  # 4*4
            nnodes = ndata // ntotal
            fmt = mapfmt(op2._endian + b'i 3f', self.size)
            struct1 = Struct(fmt)
            for unused_inode in range(nnodes):
                edata = data[n:n+ntotal]
                out = struct1.unpack(edata)
                #print(out)
                n += ntotal
        else:
            raise NotImplementedError(op2.code_information())

        self.warn_skip_table()
        return n

    def _read_oug_4(self, data: bytes, ndata: int):
        """reads the SORT1 version of table 4 (the data table)"""
        op2 = self.op2
        table_name_bytes = op2.table_name
        if op2.table_code == 1:   # Displacements
            if table_name_bytes in [b'OUGV1', b'OUGV2',
                                    b'OUG1',
                                    b'BOUGV1',
                                    b'OUPV1', b'OUG1F']:
                # OUG1F - acoustic displacements?
                #msg = f'table_name={op2.table_name} table_code={op2.table_code}'
                #raise AssertionError(msg)
                n = self._read_oug_displacement(data, ndata, is_cid=False)
            elif table_name_bytes in [b'ROUGV1', b'ROUGV2', b'TOUGV1',
                                      b'OUGF1', b'OUGF2',
                                      b'BOUGF1', ]:
                op2.to_nx(f' because table_name={op2.table_name} was found')
                n = self._read_oug_displacement(data, ndata, is_cid=False)
            elif table_name_bytes == b'OUGV1PAT':
                n = self._read_oug_displacement(data, ndata, is_cid=True)
            elif table_name_bytes == b'OAG1':
                n = self._read_oug_acceleration(data, ndata)
            elif table_name_bytes == b'OCRUG':
                n = self._read_oug_displacement(data, ndata, is_cid=False)
            else:
                raise NotImplementedError(op2.code_information())
        elif op2.table_code == 7:
            n = self._read_oug_eigenvector(data, ndata)
        elif op2.table_code == 10:
            n = self._read_oug_velocity(data, ndata)
        elif op2.table_code == 11:
            n = self._read_oug_acceleration(data, ndata)

        elif op2.table_code == 14:  # eigenvector (solution set)
            assert table_name_bytes in [b'OPHSA'], op2.table_name
            op2.to_nx(f' because table_name={op2.table_name} was found')
            n = self._read_oug_eigenvector(data, ndata)
        elif op2.table_code == 15:  # displacement (solution set)
            assert table_name_bytes in [b'OUXY1', b'OUXY2'], op2.table_name
            op2.to_nx(f' because table_name={op2.table_name} was found')
            n = self._read_oug_displacement(data, ndata, is_cid=False)
        elif op2.table_code == 16:  # velocity (solution set)
            assert table_name_bytes in [b'OUXY1', b'OUXY2'], op2.table_name
            op2.to_nx(f' because table_name={op2.table_name} was found')
            n = self._read_oug_velocity(data, ndata)
        elif op2.table_code == 17:  # acceleration (solution set)
            assert table_name_bytes in [b'OUXY1', b'OUXY2'], op2.table_name
            op2.to_nx(f' because table_name={op2.table_name} was found')
            n = self._read_oug_acceleration(data, ndata)
        elif op2.table_code == 44:   # Displacements
            assert table_name_bytes in [b'OUGMC1', b'OUGMC2'], op2.table_name
            op2.to_nx(f' because table_name={op2.table_name} was found')
            n = self._read_oug_displacement(data, ndata, is_cid=False)
        else:
            raise NotImplementedError(op2.code_information())
        return n

    #def _read_eigenvector_displacement_solution_set(self, data: bytes, ndata: int):
        #"""
        #table_code = 14
        #"""
        #raise NotImplementedError()

    #def _read_displacement_solution_set(self, data: bytes, ndata: int):
        #"""
        #table_code = 15
        #"""
        #raise NotImplementedError()

    #def _read_velocity_solution_set(self, data: bytes, ndata: int):
        #"""
        #table_code = 16
        #"""
        #raise NotImplementedError()

    #def _read_acceleration_solution_set(self, data: bytes, ndata: int):
        #"""
        #table_code = 17
        #"""
        #raise NotImplementedError()

    def _read_oug_displacement(self, data: bytes, ndata: int, is_cid: bool) -> int:
        """
        Table     Description
        -----     -----------
        OUG1      displacements in the global? frame
        OUGV1/2   displacements in the global frame
        OUGV1PAT  displacements in the global? frame
        BOUGV1    displacments in the basic frame
        ROUGV1/2  relative displacments in the global frame
        OUPV1     ???
        OUXY1/2   eigenvectors in the basic frame
        TOUGV1/2  temperature
        OCRUG     ???
        OUG1F     acoustic displacements
        OUGF1     acoustic displacements

        """
        op2 = self.op2
        op2._setup_op2_subcase('displacement')

        if op2.table_name in [b'ROUGV1', b'ROUGV2']:
            assert op2.thermal in [0], op2.code_information()
            result_name = 'displacements_ROUGV1'

        elif op2.table_name in [b'OUG1', b'OUGV1', b'OUGV2', b'OUGV1PAT', b'BOUGV1']:
            # OUG1F - acoustic displacements
            assert op2.thermal in [0, 1], op2.code_information()
            # NX THERMAL
            # 1: heat transfer
            # 2: axisymmetric Fourier
            # 3: for cyclic symmetric;
            # 0: otherwise
            if op2.thermal == 0:
                result_name = 'displacements'
            elif op2.thermal == 1:
                result_name = 'temperatures'
            else:  # pragma: no cover
                msg = 'displacements; table_name=%s' % op2.table_name
                raise NotImplementedError(msg)
        elif op2.table_name in [b'OUXY1', b'OUXY2']:
            assert op2.thermal == 0, op2.code_information()
            result_name = 'solution_set.displacements'
        elif op2.table_name == b'OUPV1':
            #result_name = 'temperatures'
            result_name0 = 'displacements' # is this right?
            prefix, postfix = _oug_get_prefix_postfix(op2.thermal)
            result_name = prefix + result_name0 + postfix

        elif op2.table_name in [b'TOUGV1', b'TOUGV2']:
            result_name = 'temperatures'
            assert op2.thermal == 1, op2.code_information()
        elif op2.table_name in [b'OCRUG']:
            result_name = 'displacements'
            assert op2.thermal == 0, op2.code_information()
        elif op2.table_name in [b'OUG1F', b'OUGF1', b'OUGF2', b'BOUGF1']:
            result_name = 'acoustic.displacements'  # acoustic displacements
            assert op2.thermal == 0, op2.code_information()
        else:  # pragma: no cover
            msg = 'displacements; table_name=%s' % op2.table_name
            raise NotImplementedError(msg)

        if op2._results.is_not_saved(result_name):
            return ndata
        op2._results._found_result(result_name)
        storage_obj = op2.get_result(result_name)
        if op2.thermal == 0:
            #result_name = 'displacements'
            #storage_obj = self.displacements
            assert op2.table_name in [b'BOUGV1', b'ROUGV1', b'ROUGV2', b'OUGV1', b'OUGV2',
                                       b'OUG1', b'OCRUG', b'OUGV1PAT', b'OUXY1', b'OUXY2',
                                       b'OUG1F',
                                       b'OUGF1', b'OUGF2',
                                       b'BOUGF1', ], op2.table_name
            n = op2._read_table_vectorized(data, ndata, result_name, storage_obj,
                                            RealDisplacementArray, ComplexDisplacementArray,
                                            'node', random_code=op2.random_code,
                                            is_cid=is_cid)
        elif op2.thermal == 1:
            assert op2.table_name in [b'OUGV1', b'OUGV2', b'TOUGV1', b'TOUGV2', b'OUG1'], op2.table_name
            n = op2._read_scalar_table_vectorized(data, ndata, result_name, storage_obj,
                                                  RealTemperatureArray, None,
                                                  'node', random_code=op2.random_code,
                                                  is_cid=is_cid)
        elif op2.thermal == 2:
            assert op2.table_name in [b'OUPV1'], op2.table_name
            n = op2._read_table_vectorized(data, ndata, result_name, storage_obj,
                                           RealDisplacementArray, ComplexDisplacementArray,
                                           'node', random_code=op2.random_code)
        elif op2.thermal == 4:
            # F:\work\pyNastran\examples\Dropbox\move_tpl\ms103.op2
            assert op2.table_name in [b'OUPV1'], op2.table_name
            n = op2._read_table_vectorized(data, ndata, result_name, storage_obj,
                                           RealDisplacementArray, ComplexDisplacementArray,
                                           'node', random_code=op2.random_code)
        elif op2.thermal == 8:  # 4 ?
            assert op2.table_name in [b'OUPV1'], op2.table_name
            n = op2._read_table_vectorized(data, ndata, result_name, storage_obj,
                                           RealDisplacementArray, ComplexDisplacementArray,
                                           'node', random_code=op2.random_code)
        else:
            raise RuntimeError(op2.code_information())
            #n = op2._not_implemented_or_skip(data, ndata, 'bad thermal=%r table' % op2.thermal)
        #else:
            #raise NotImplementedError(op2.thermal)
        return n

    def _read_oug_velocity(self, data: bytes, ndata: int):
        """
        table_code = 10
        """
        op2 = self.op2
        op2._setup_op2_subcase('velocity')
        if op2.table_name in [b'OUGV1', b'OUGV2', b'BOUGV1', b'OVG1']:
            assert op2.thermal in [0, 1], op2.code_information()
            result_name = 'velocities'
        elif op2.table_name in [b'OUXY1', b'OUXY2']:
            op2.to_nx(f' because table_name={op2.table_name} was found')
            assert op2.thermal == 0, op2.code_information()
            result_name = 'solution_set.velocities'
        elif op2.table_name in [b'ROUGV1', b'ROUGV2']:
            op2.to_nx(f' because table_name={op2.table_name} was found')
            result_name = 'velocities_ROUGV1'
            assert op2.thermal == 0, op2.code_information()
        elif op2.table_name == b'OUPV1':
            result_name0 = 'velocities'
            prefix, postfix = _oug_get_prefix_postfix(op2.thermal)
            result_name = prefix + result_name0 + postfix
            assert op2.thermal in [2, 4], op2.thermal
        else:  # pragma: no cover
            msg = 'velocities; table_name=%s' % op2.table_name
            raise NotImplementedError(msg)

        #result_name = 'velocities'
        #storage_obj = self.velocities
        if op2._results.is_not_saved(result_name):
            return ndata
        op2._results._found_result(result_name)
        storage_obj = op2.get_result(result_name)
        if op2.thermal == 0:
            n = op2._read_table_vectorized(data, ndata, result_name, storage_obj,
                                            RealVelocityArray, ComplexVelocityArray,
                                            'node', random_code=op2.random_code)
        elif op2.thermal == 1:
            n = op2._read_scalar_table_vectorized(data, ndata, result_name, storage_obj,
                                                   RealThermalVelocityVectorArray, None,
                                                   'node', random_code=op2.random_code)

        elif op2.thermal == 2:
            n = op2._read_table_vectorized(data, ndata, result_name, storage_obj,
                                           RealVelocityArray, ComplexVelocityArray,
                                           'node', random_code=op2.random_code)
        else:
            raise NotImplementedError(op2.thermal)
        return n

    def _read_oug_acceleration(self, data: bytes, ndata: int):
        """
        table_code = 11
        """
        op2 = self.op2
        op2._setup_op2_subcase('acceleration')

        result_name = None
        if op2.table_name in [b'OUGV1', b'OUGV2', b'OAG1', b'BOUGV1']:
            result_name = 'accelerations'
            assert op2.thermal == 0, op2.code_information()
        elif op2.table_name in [b'OUXY1', b'OUXY2']:
            op2.to_nx(f' because table_name={op2.table_name} was found')
            assert op2.thermal == 0, op2.code_information()
            result_name = 'solution_set.accelerations'
        elif op2.table_name in [b'ROUGV1', b'ROUGV2']:
            op2.to_nx(f' because table_name={op2.table_name} was found')
            result_name = 'accelerations_ROUGV1'
            assert op2.thermal == 0, op2.code_information()
        elif op2.table_name in [b'OAGPSD1', b'OAGPSD2',
                                b'OAGRMS1', b'OAGRMS2',
                                b'OACRM1', b'OAGCRM2',
                                b'OAGNO1', b'OAGNO2']:
            assert op2.thermal == 0, op2.code_information()
            pass
        elif op2.table_name == b'OUPV1':
            assert op2.thermal in [0, 2, 4], op2.thermal
            result_name0 = 'accelerations'
            prefix, postfix = _oug_get_prefix_postfix(op2.thermal)
            result_name = prefix + result_name0 + postfix
        else:  # pragma: no cover
            msg = 'accelerations; table_name=%s' % op2.table_name
            raise NotImplementedError(msg)

        if op2.thermal == 0:
            if op2.table_name in [b'OUGV1', b'OUGV2', b'ROUGV1', b'ROUGV2', b'OAG1', b'BOUGV1', b'OUXY1', b'OUXY2', b'OUPV1']:
                assert result_name is not None, op2.table_name
                if op2._results.is_not_saved(result_name):
                    return ndata
                storage_obj = op2.get_result(result_name)
                n = op2._read_table_vectorized(data, ndata, result_name, storage_obj,
                                                RealAccelerationArray,
                                                ComplexAccelerationArray,
                                                'node', random_code=op2.random_code)
            elif op2.table_name in [b'OAGPSD1', b'OAGPSD2']:
                n = self._read_oug_psd(data, ndata)
            elif op2.table_name in [b'OAGRMS1', b'OAGRMS2']:
                n = self._read_oug_rms(data, ndata)
            elif op2.table_name in [b'OACRM1', b'OAGCRM2']:
                n = self._read_oug_crm(data, ndata)
            elif op2.table_name in [b'OAGNO1', b'OAGNO2']:
                n = self._read_oug_no(data, ndata)
            else:
                raise NotImplementedError(op2.code_information())
        elif op2.thermal == 1:
            result_name = 'accelerations'
            storage_obj = op2.accelerations
            if op2._results.is_not_saved(result_name):
                return ndata
            op2._results._found_result(result_name)
            raise NotImplementedError(op2.code_information())
        elif op2.thermal == 2:
            result_name = 'abs.accelerations'
            storage_obj = op2.get_result(result_name)
            if op2._results.is_not_saved(result_name):
                return ndata
            op2._results._found_result(result_name)
            n = op2._read_table_vectorized(data, ndata, result_name, storage_obj,
                                            RealAccelerationArray, ComplexAccelerationArray,
                                            'node', random_code=op2.random_code)
        elif op2.thermal == 4:
            result_name = 'srss.accelerations'
            storage_obj = op2.get_result(result_name)
            if op2._results.is_not_saved(result_name):
                return ndata
            op2._results._found_result(result_name)
            n = op2._read_table_vectorized(data, ndata, result_name, storage_obj,
                                            RealAccelerationArray, ComplexAccelerationArray,
                                            'node', random_code=op2.random_code)
        else:
            raise NotImplementedError(op2.thermal)
        return n

    def _read_oug_eigenvector(self, data: bytes, ndata: int):
        """
        table_code = 7
        """
        op2 = self.op2
        # NX THERMAL
        # 1: heat transfer
        # 2: axisymmetric Fourier
        # 3: for cyclic symmetric;
        # 0: otherwise
        assert op2.thermal in [0, 2, 3], op2.code_information()
        if op2.table_name in [b'OUGV1', b'OUGV2', b'OUG1',
                               b'BOUGV1',
                               b'OPHIG', b'BOPHIG', ]:
            op2._setup_op2_subcase('VECTOR')
            result_name = 'eigenvectors'
        elif op2.table_name in [b'OUGF1', b'OUGF2',
                                 b'BOUGF1',
                                 b'BOPHIGF']:
            op2._setup_op2_subcase('VECTOR')
            result_name = 'eigenvectors_fluid'

        elif op2.table_name == b'OPHSA':
            op2.to_nx(f' because table_name={op2.table_name} was found')
            op2._setup_op2_subcase('SVECTOR')
            assert op2.thermal == 0, op2.code_information()
            result_name = 'solution_set.eigenvectors'

        elif op2.table_name == b'RADCONS':
            op2.to_nx(f' because table_name={op2.table_name} was found')
            op2._setup_op2_subcase('VECTOR')
            result_name = 'RADCONS.eigenvectors'
        elif op2.table_name == b'RADEFFM':
            op2.to_nx(f' because table_name={op2.table_name} was found')
            op2._setup_op2_subcase('VECTOR')
            result_name = 'RADEFFM.eigenvectors'
        elif op2.table_name == b'RADEATC':
            op2.to_nx(f' because table_name={op2.table_name} was found')
            op2._setup_op2_subcase('VECTOR')
            result_name = 'RADEATC.eigenvectors'
        elif op2.table_name in [b'ROUGV1', 'ROUGV2']:
            op2.to_nx(f' because table_name={op2.table_name} was found')
            op2._setup_op2_subcase('VECTOR')
            result_name = 'ROUGV1.eigenvectors'
        else:  # pragma: no cover
            msg = 'eigenvectors; table_name=%s' % op2.table_name
            raise NotImplementedError(msg)
        assert op2.thermal in [0, 2, 3], op2.code_information()

        if op2._results.is_not_saved(result_name):
            return ndata
        op2._results._found_result(result_name)
        storage_obj = op2.get_result(result_name)

        # NX THERMAL
        # 1: heat transfer
        # 2: axisymmetric Fourier
        # 3: for cyclic symmetric;
        # 0: otherwise
        if op2.thermal in [0, 2, 3]:
            n = op2._read_table_vectorized(data, ndata, result_name, storage_obj,
                                           RealEigenvectorArray, ComplexEigenvectorArray,
                                           'node', random_code=op2.random_code)
        elif op2.thermal == 1:
            n = op2._not_implemented_or_skip(data, ndata, msg='thermal=1')
        else:
            raise NotImplementedError(op2.thermal)
        return n

    def _read_oug_psd(self, data: bytes, ndata: int):
        """
        table_code = 601/610/611

        +-----+-------------+---------+
        | Bit |     0       |    1    |
        +-----+-------------+---------+
        |  0  | Not Random  | Random  |
        |  1  | SORT1       | SORT2   |
        |  2  | Real        | Complex |
        +-----+-------------+---------+

          sort_code = 0 -> sort_bits = [0,0,0]  #         sort1, real
          sort_code = 1 -> sort_bits = [0,0,1]  #         sort1, complex
          sort_code = 2 -> sort_bits = [0,1,0]  #         sort2, real
          sort_code = 3 -> sort_bits = [0,1,1]  #         sort2, complex
          sort_code = 4 -> sort_bits = [1,0,0]  # random, sort1, real
          sort_code = 5 -> sort_bits = [1,0,1]  # random, sort1, real
          sort_code = 6 -> sort_bits = [1,1,0]  # random, sort2, real
          sort_code = 7 -> sort_bits = [1,1,1]  # random, sort2, complex
          # random, sort2, complex <- [1, 1, 1]

          sort_bits[0] = 0 -> isSorted=True isRandom=False
          sort_bits[1] = 0 -> is_sort1=True is_sort2=False
          sort_bits[2] = 0 -> isReal=True   isReal/Imaginary=False
        """
        op2 = self.op2
        #self.sort_code = 6
        #self.sort_bits = [1, 1, 0]
        #if op2.table_code < 50:
        #    op2.table_code += 600

        if op2.thermal == 0:
            if op2.table_code == 1:
                # displacement
                assert op2.table_name in [b'OUGPSD1', b'OUGPSD2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'psd.displacements'
                obj = RealDisplacementArray
            elif op2.table_code == 10:
                # velocity
                assert op2.table_name in [b'OVGPSD1', b'OVGPSD2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'psd.velocities'
                obj = RealVelocityArray
            elif op2.table_code == 11:
                # acceleration
                assert op2.table_name in [b'OAGPSD1', b'OAGPSD2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'psd.accelerations'
                obj = RealAccelerationArray

            elif op2.table_code == 601:
                # displacement
                assert op2.table_name in [b'OUGPSD1', b'OUGPSD2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'psd.displacements'
                obj = RealDisplacementArray
            elif op2.table_code == 610:
                # velocity
                assert op2.table_name in [b'OUGPSD1', b'OUGPSD2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'psd.velocities'
                obj = RealVelocityArray
            elif op2.table_code == 611:
                # acceleration
                assert op2.table_name in [b'OUGPSD1', b'OUGPSD2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'psd.accelerations'
                obj = RealAccelerationArray
            else:
                n = op2._not_implemented_or_skip(data, ndata, op2.code_information())
                return n

            if op2._results.is_not_saved(result_name):
                return ndata
            op2._results._found_result(result_name)

            storage_obj = op2.get_result(result_name)
            n = op2._read_random_table(data, ndata, result_name, storage_obj,
                                        obj, 'node',
                                        random_code=op2.random_code)
        #elif op2.thermal == 1:
            #result_name = 'accelerations'
            #storage_obj = op2.accelerations
            #if op2._results.is_not_saved(result_name):
                #return ndata
            #op2._results._found_result(result_name)
            #n = self._read_table(data, ndata, result_name, storage_obj,
                                 #None, None,
                                 #None, None, 'node', random_code=op2.random_code)
        #elif op2.thermal == 2:
            #result_name = 'acceleration_scaled_response_spectra_abs'
            #storage_obj = op2.acceleration_scaled_response_spectra_abs
            #if op2._results.is_not_saved(result_name):
                #return ndata
            #op2._results._found_result(result_name)
            #n = self._read_table(data, ndata, result_name, storage_obj,
                                 #RealAcceleration, ComplexAcceleration,
                                 #RealAccelerationArray, ComplexAccelerationArray,
                                 #'node', random_code=op2.random_code)
            ##n = op2._not_implemented_or_skip(data, ndata, msg='thermal=2')
        #elif op2.thermal == 4:
            #result_name = 'acceleration_scaled_response_spectra_nrl'
            #storage_obj = op2.acceleration_scaled_response_spectra_nrl
            #if op2._results.is_not_saved(result_name):
                #return ndata
            #op2._results._found_result(result_name)
            #n = self._read_table(data, ndata, result_name, storage_obj,
                                 #RealAcceleration, ComplexAcceleration,
                                 #RealAccelerationArray, ComplexAccelerationArray,
                                 #'node', random_code=op2.random_code)
            ##n = op2._not_implemented_or_skip(data, ndata, msg='thermal=4')
        else:
            raise NotImplementedError(op2.thermal)
        return n

    def _read_oug_rms(self, data: bytes, ndata: int):
        """
        table_code = 801  # /610/611
        """
        op2 = self.op2
        #self.sort_code = 6
        #self.sort_bits = [1, 1, 0]
        #if op2.table_code < 50:
        #    op2.table_code += 800

        if op2.thermal == 0:
            if op2.table_code == 1:
                # displacement
                assert op2.table_name in [b'OUGRMS1', b'OUGRMS2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'rms.displacements'
                obj = RealDisplacementArray
            elif op2.table_code == 10:
                # velocity
                assert op2.table_name in [b'OVGRMS1', b'OVGRMS2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'rms.velocities'
                obj = RealVelocityArray
            elif op2.table_code == 11:
                # acceleration
                assert op2.table_name in [b'OAGRMS1', b'OAGRMS2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'rms.accelerations'
                obj = RealAccelerationArray
            elif op2.table_code == 801:
                result_name = 'rms.displacements'
                assert op2.table_name in [b'OUGRMS1', b'OUGRM2'], 'op2.table_name=%r' % op2.table_name
                obj = RealDisplacementArray
            elif op2.table_code == 810:
                assert op2.table_name in [b'OUGRMS1', b'OUGRM2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'rms.velocities'
                obj = RealVelocityArray
            elif op2.table_code == 811:
                assert op2.table_name in [b'OUGRMS1', b'OUGRMS2'], 'op2.table_name=%r' % op2.table_name # , b'OAGRMS1', b'OAGRMS2'
                result_name = 'rms.accelerations'
                obj = RealAccelerationArray
            else:
                if DEV:  # pragma: no cover
                    raise RuntimeError(op2.code_information())
                n = op2._not_implemented_or_skip(data, ndata, op2.code_information())
                return n

            if op2._results.is_not_saved(result_name):
                return ndata
            op2._results._found_result(result_name)

            storage_obj = op2.get_result(result_name)
            n = op2._read_random_table(data, ndata, result_name, storage_obj,
                                       obj, 'node',
                                       random_code=op2.random_code)
            #n = self._read_table_sort1_real(data, ndata, result_name, storage_obj,
                                            #RealDisplacementArray, 'node',
                                            #random_code=op2.random_code)
            #n = op2._read_table_vectorized(data, ndata, result_name, storage_obj,
                                            #RealDisplacementArray, ComplexDisplacementArray,
                                            #'node')

        #elif op2.thermal == 1:
            #result_name = 'accelerations'
            #storage_obj = op2.accelerations
            #if op2._results.is_not_saved(result_name):
                #return ndata
            #op2._results._found_result(result_name)
            #n = self._read_table(data, ndata, result_name, storage_obj,
                                 #None, None,
                                 #None, None, 'node', random_code=op2.random_code)
        #elif op2.thermal == 2:
            #result_name = 'acceleration_scaled_response_spectra_abs'
            #storage_obj = op2.acceleration_scaled_response_spectra_abs
            #if op2._results.is_not_saved(result_name):
                #return ndata
            #op2._results._found_result(result_name)
            #n = self._read_table(data, ndata, result_name, storage_obj,
                                 #RealAcceleration, ComplexAcceleration,
                                 #RealAccelerationArray, ComplexAccelerationArray,
                                 #'node', random_code=op2.random_code)
            ##n = op2._not_implemented_or_skip(data, ndata, msg='thermal=2')
        #elif op2.thermal == 4:
            #result_name = 'acceleration_scaled_response_spectra_nrl'
            #storage_obj = op2.acceleration_scaled_response_spectra_nrl
            #if op2._results.is_not_saved(result_name):
                #return ndata
            #op2._results._found_result(result_name)
            #n = self._read_table(data, ndata, result_name, storage_obj,
                                 #RealAcceleration, ComplexAcceleration,
                                 #RealAccelerationArray, ComplexAccelerationArray,
                                 #'node', random_code=op2.random_code)
            ##n = op2._not_implemented_or_skip(data, ndata, msg='thermal=4')
        else:
            raise NotImplementedError(op2.thermal)
        return n

    def _read_oug_no(self, data: bytes, ndata: int):
        """
        table_code = 901  # /610/611
        """
        op2 = self.op2
        if op2.thermal == 0:
            if op2.table_code == 1:
                # displacement
                assert op2.table_name in [b'OUGNO1', b'OUGNO2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'no.displacements'
                obj = RealDisplacementArray
            elif op2.table_code == 10:
                # velocity
                assert op2.table_name in [b'OVGNO1', b'OVGNO2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'no.velocities'
                obj = RealVelocityArray
            elif op2.table_code == 11:
                # acceleration
                assert op2.table_name in [b'OAGNO1', b'OAGNO2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'no.accelerations'
                obj = RealAccelerationArray

            elif op2.table_code == 901:
                assert op2.table_name in [b'OUGNO1', b'OUGNO2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'no.displacements'
                obj = RealDisplacementArray
            elif op2.table_code == 910:
                assert op2.table_name in [b'OUGNO1', b'OUGNO2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'no.velocities'
                obj = RealVelocityArray
            elif op2.table_code == 911:
                assert op2.table_name in [b'OUGNO1', b'OUGNO2', b'OAGNO1', b'OAGNO2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'no.accelerations'
                obj = RealAccelerationArray
            else:
                n = op2._not_implemented_or_skip(data, ndata, op2.code_information())
                return n
            if op2._results.is_not_saved(result_name):
                return ndata
            op2._results._found_result(result_name)

            storage_obj = op2.get_result(result_name)
            n = op2._read_random_table(data, ndata, result_name, storage_obj,
                                       obj, 'node',
                                       random_code=op2.random_code)

        #elif op2.thermal == 1:
            #result_name = 'accelerations'
            #storage_obj = op2.accelerations
            #if op2._results.is_not_saved(result_name):
                #return ndata
            #op2._results._found_result(result_name)
            #n = self._read_table(data, ndata, result_name, storage_obj,
                                 #None, None,
                                 #None, None, 'node', random_code=op2.random_code)
        #elif op2.thermal == 2:
            #result_name = 'acceleration_scaled_response_spectra_abs'
            #storage_obj = op2.acceleration_scaled_response_spectra_abs
            #if op2._results.is_not_saved(result_name):
                #return ndata
            #op2._results._found_result(result_name)
            #n = self._read_table(data, ndata, result_name, storage_obj,
                                 #RealAcceleration, ComplexAcceleration,
                                 #RealAccelerationArray, ComplexAccelerationArray,
                                 #'node', random_code=op2.random_code)
            ##n = op2._not_implemented_or_skip(data, ndata, msg='thermal=2')
        #elif op2.thermal == 4:
            #result_name = 'acceleration_scaled_response_spectra_nrl'
            #storage_obj = op2.acceleration_scaled_response_spectra_nrl
            #if op2._results.is_not_saved(result_name):
                #return ndata
            #op2._results._found_result(result_name)
            #n = self._read_table(data, ndata, result_name, storage_obj,
                                 #RealAcceleration, ComplexAcceleration,
                                 #RealAccelerationArray, ComplexAccelerationArray,
                                 #'node', random_code=op2.random_code)
            ##n = op2._not_implemented_or_skip(data, ndata, msg='thermal=4')
        else:
            raise NotImplementedError(op2.thermal)
        return n

    def _read_oug_ato(self, data: bytes, ndata: int):
        """
        table_code = 901  # /610/611
        """
        op2 = self.op2
        if op2.thermal == 0:
            if op2.table_code == 1:
                result_name = 'ato.displacements'
                obj = RealDisplacementArray
                assert op2.table_name in [b'OUGATO1', b'OUGATO2'], 'op2.table_name=%r' % op2.table_name
            elif op2.table_code == 10:
                result_name = 'ato.velocities'
                obj = RealVelocityArray
                assert op2.table_name in [b'OVGATO1', b'OVGATO2'], 'op2.table_name=%r' % op2.table_name
            elif op2.table_code == 11:
                result_name = 'ato.accelerations'
                obj = RealAccelerationArray
                assert op2.table_name in [b'OAGATO1', b'OAGATO2'], 'op2.table_name=%r' % op2.table_name
            else:
                n = op2._not_implemented_or_skip(data, ndata, op2.code_information())
                return n
        else:
            raise NotImplementedError(op2.thermal)

        if op2._results.is_not_saved(result_name):
            return ndata
        op2._results._found_result(result_name)
        storage_obj = op2.get_result(result_name)
        n = op2._read_random_table(data, ndata, result_name, storage_obj,
                                   obj, 'node',
                                   random_code=op2.random_code)
        return n

    def _read_oug_crm(self, data: bytes, ndata: int):
        """
        table_code = 501  # /510/511
        """
        op2 = self.op2
        #self.sort_code = 6
        #self.sort_bits = [1, 1, 0]
        #if op2.table_code < 50:
        #    op2.table_code += 800

        if op2.thermal == 0:
            if op2.table_code == 1:
                assert op2.table_name in [b'OUGCRM1', b'OUGCRM2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'crm.displacements'
                obj = RealDisplacementArray
            elif op2.table_code == 10:
                # velocity
                assert op2.table_name in [b'OVGCRM1', b'OVGCRM2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'crm.velocities'
                obj = RealVelocityArray
            elif op2.table_code == 11:
                # acceleration
                assert op2.table_name in [b'OAGCRM1', b'OAGCRM2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'crm.accelerations'
                obj = RealAccelerationArray
            elif op2.table_code == 501:
                assert op2.table_name in [b'OUGCRM1', b'OUGCRM2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'crm.displacements'
                obj = RealDisplacementArray
            elif op2.table_code == 510:
                # velocity
                assert op2.table_name in [b'OUGCRM1', b'OUGCRM2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'crm.velocities'
                obj = RealVelocityArray
            elif op2.table_code == 511:
                # acceleration
                assert op2.table_name in [b'OUGCRM1', b'OUGCRM2'], 'op2.table_name=%r' % op2.table_name
                result_name = 'crm.accelerations'
                obj = RealAccelerationArray
            else:
                n = op2._not_implemented_or_skip(data, ndata, op2.code_information())
                #raise RuntimeError(op2.code_information())
                return n

            if op2._results.is_not_saved(result_name):
                return ndata
            op2._results._found_result(result_name)

            storage_obj = op2.get_result(result_name)
            n = op2._read_random_table(data, ndata, result_name, storage_obj,
                                        obj, 'node',
                                        random_code=op2.random_code)

                #n = self._read_table_sort1_real(data, ndata, result_name, storage_obj,
                                                #RealDisplacementArray, 'node',
                                                #random_code=op2.random_code)
                #n = op2._read_table_vectorized(data, ndata, result_name, storage_obj,
        else:
            raise NotImplementedError(op2.thermal)
        return n


def _oug_get_prefix_postfix(thermal: int) -> Tuple[str, str]:
    prefix = ''
    postfix = ''
    if thermal == 0:
        pass
    if thermal == 2:
        prefix = 'abs.'
    elif thermal == 4:
        prefix = 'srss.'
    elif thermal == 8:
        prefix = 'nrl.'
    else:  # pragma: no cover
        msg = 'thermal=%s' % thermal
        raise NotImplementedError(msg)

    #assert op2.thermal in [0, 2, 4, 8], op2.code_information()
    #if op2.thermal == 0:
        #result_name = 'displacement' # is this right?
    #elif op2.thermal == 2:
        #result_name = 'abs.displacement' # displacement_scaled_response_spectra_abs
    #elif op2.thermal == 4:
        #result_name = 'srss.displacement' # displacement_scaled_response_spectra_srss
    #elif op2.thermal == 8:
        #result_name = 'nrl.displacement'  # displacement_scaled_response_spectra_nrl
    #else:  # pragma: no cover
        #msg = 'displacements; table_name=%s' % op2.table_name
        #raise NotImplementedError(msg)

    return prefix, postfix
