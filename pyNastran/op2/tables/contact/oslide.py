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
from pyNastran.op2.tables.utils import get_is_slot_saved

from pyNastran.op2.tables.contact.slide_objects import (
    RealSlideDistanceArray, ) # ComplexDisplacementArray

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


class OSLIDE:
    """
      TIME =  1.000000E+00
                               C O N T A C T     S L I D E     D I S T A N C E
                                         INCREMENTAL                                      TOTAL
      POINT ID.   TYPE          T1             T2             T3             T1             T2             T3
            57      G     -1.783432E-03   1.065814E-14   0.0           -1.783432E-03   1.065814E-14  -7.993606E-15

    #OUG : Output U in the global frame

    #U is:
     #- Displacement
     #- Velocity
     #- Accelerations

    #The global frame is:
     #- the analysis coordinate frame, not the 0 coordinate frame
     """
    def __init__(self, op2: OP2):
        self.op2 = op2

    @property
    def size(self) -> int:
        return self.op2.size
    @property
    def factor(self) -> int:
        return self.op2.factor

    #def update_mode_cycle(self, name):
        #op2 = self.op2
        #value = getattr(op2, name)
        #if value == 0.0:
            ##print('table_name=%r mode=%s eigr=%s' % (op2.table_name, op2.mode, op2.eigr))
            #value = np.sqrt(np.abs(op2.eign)) / (2. * np.pi)
            #setattr(op2, name, value)
            #op2.data_code[name] = value

    def read_sort1_3(self, data: bytes, ndata: int):
        """reads table 3 (the header table)"""
        op2 = self.op2
        assert ndata == 146 * op2.size
        #op2.show_data(data)
        #self._set_times_dtype()
        op2.nonlinear_factor = np.nan
        op2.is_table_1 = True
        op2.is_table_2 = False
        op2.parse_approach_code(data)  # field 3
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
            op2.show_data(data)
            #op2._write_data(op2.binary_debug, data, types='ifs')
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
        #self._correct_eigenvalue()

    #def _correct_eigenvalue(self):
        #"""Nastran 95 gets the frequency wrong"""
        #op2 = self.op2

    #def read_sort2_3(self, data: bytes, ndata: int) -> None:
        #"""reads the SORT2 version of table 4 (the data table)"""
        #op2 = self.op2
        ##self._set_times_dtype()
        ##return self._read_oug1_3(data)
        #op2.nonlinear_factor = np.nan

        #op2.is_table_1 = False
        #op2.is_table_2 = True
        #unused_three = op2.parse_approach_code(data)
        #op2.words = [
            #'approach_code', 'table_code', '???', 'isubcase',
            #'???', '???', '???', 'random_code',
            #'format_code', 'num_wide', '???', '???',
            #'acoustic_flag', '???', '???', '???',
            #'???', '???', '???', '???',
            #'???', '???', 'thermal', '???',
            #'???', 'Title', 'subtitle', 'label']

        ### random code
        #op2.random_code = op2.add_data_parameter(data, 'random_code', b'i', 8, False)

        ### format code
        #op2.format_code = op2.add_data_parameter(data, 'format_code', b'i', 9, False)

        ### number of words per entry in record
        #op2.num_wide = op2.add_data_parameter(data, 'num_wide', b'i', 10, False)

        ### acoustic pressure flag
        #op2.acoustic_flag = op2.add_data_parameter(data, 'acoustic_flag', b'i', 13, False)

        ### thermal flag; 1 for heat transfer, 0 otherwise
        #op2.thermal = op2.add_data_parameter(data, 'thermal', b'i', 23, False)

        #op2.node_id = op2.add_data_parameter(data, 'node_id', b'i', 5, fix_device_code=True)
        ##if op2.analysis_code == 1:  # statics / displacement / heat flux
            ### load set number
            ##op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5, False)
            ##op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            ##op2.setNullNonlinearFactor()

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
        ##elif op2.analysis_code == 3: # differential stiffness
            ##op2.lsdvmn = self.get_values(data, b'i', 5) ## load set number
            ##op2.data_code['lsdvmn'] = op2.lsdvmn
        ##elif op2.analysis_code == 4: # differential stiffness
            ##op2.lsdvmn = self.get_values(data, b'i', 5) ## load set number
        #elif op2.analysis_code == 5:   # frequency
            ## frequency
            ##op2.freq = op2.add_data_parameter(data, 'freq', b'f', 5)
            #op2._analysis_code_fmt = b'f'
            #op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            #op2.apply_data_code_value('analysis_method', 'freq')
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
        #else:
            #msg = 'invalid analysis_code...analysis_code=%s' % op2.analysis_code
            #raise RuntimeError(msg)

        #op2._fix_oug_format_code()
        #op2._parse_thermal_code()
        #if op2.is_debug_file:
            #op2.binary_debug.write('  %-14s = %r %s\n' % ('approach_code', op2.approach_code,
                                                           #op2.approach_code_str(op2.approach_code)))
            #op2.binary_debug.write('  %-14s = %r\n' % ('tCode', op2.tCode))
            #op2.binary_debug.write('  %-14s = %r\n' % ('isubcase', op2.isubcase))
        #op2._read_title(data)
        #op2._write_debug_bits()
        #assert isinstance(op2.nonlinear_factor, integer_types), op2.nonlinear_factor

    def read_4(self, data: bytes, ndata: int):
        """reads table 4 (the results table)"""
        op2 = self.op2
        assert op2.table_code == 74, op2.code_information()
        #if op2.read_mode == 1:
            #return ndata
        #self.show_data(data)
        #print(f'data_type = {op2.data_type}')
        #assert op2.data_type == 1, op2.data_type

        #print('op2.analysis_code =', op2.analysis_code)
        #op2.show_data(data, types='ifs', endian=None, force=False)

        if op2.table_name == b'OSLIDEG1':
            result_name = 'glue_contact_slide_distance'
        elif op2.table_name == b'OSLIDE1':
            result_name = 'contact_slide_distance'
        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())

        is_saved, slot = get_is_slot_saved(op2, result_name)
        if not is_saved:
            return ndata

        if op2.num_wide == 7:
            factor = op2.factor
            ntotal = 28 * factor
            nnodes = ndata // ntotal  # 8*4
            auto_return = op2._create_table_vector(
                result_name, nnodes, slot, RealSlideDistanceArray, is_cid=False)
            if auto_return:
                return ndata

            assert op2.format_code == 2, op2.code_information
            ints = np.frombuffer(data, dtype=op2.idtype8)
            nfields = len(ints)
            nrows = nfields // 7
            assert nfields % 7 == 0
            ints = ints.reshape(nrows, 7)
            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nrows, 7)

            if op2.is_sort1 or op2.analysis_code in {1, 2, 3, 4, 7, 8, 9, 11, 12}:
                # 1: static
                # 2: modes
                # 3: diff. stiff 0
                # 4: diff stiff 1
                # 7: pre-buckling
                # 8: post-buckling
                # 9: complex modes
                # 11: geometric nonlinear statics
                # 12 contran
                node_id = ints[:, 0] // 10
                datai = floats[:, 1:]
                obj = op2.obj
                if op2.analysis_code not in {1}:
                    obj._times[obj.itime] = obj.nonlinear_factor
                obj.node_gridtype[:, 0] = node_id
                obj.node_gridtype[:, 1] = 1
                obj.data[obj.itime, :, :] = datai
            else:  # pragma: no cover
                # 5: freq
                # 6: time step
                # 10 nonlinear statics
                raise RuntimeError(op2.analysis_code)

            #print('node_id', node_id)
            #print(datai)
            #reals = floats[:, [0, 1, 2]]
            #imags = floats[:, [3, 4, 5]]
        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
        #raise NotImplementedError(op2.code_information())
        return ndata
