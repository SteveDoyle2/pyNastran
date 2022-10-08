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

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


class OUGPK1:
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

    #def update_mode_cycle(self, name):
        #op2 = self.op2
        #value = getattr(op2, name)
        #if value == 0.0:
            ##print('table_name=%r mode=%s eigr=%s' % (op2.table_name, op2.mode, op2.eigr))
            #value = np.sqrt(np.abs(op2.eign)) / (2. * np.pi)
            #setattr(op2, name, value)
            #op2.data_code[name] = value

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
        table_name_bytes = op2.table_name

        if data is None:
            return ndata

        is_vectorized = True
        if op2.table_code == 401:
            restype = 'displacement?'
        elif op2.table_code == 411:
            restype = 'velocity?'
        else:
            raise RuntimeError(op2.table_code)

        #ints    = (11, 1, 1057254660, 0, 0, 0, 0, 0,
                   #21, 1, 0, 0, 0, 0, 0, 0,
                   #31, 1, 1057273077, 690366181, 688497915, 0, 709531316, 734276092,
                   #41, 1, 1057257192, 753975636, 756267685, 0, 731678033, 720573197,
                   #51, 1, 1057260305, 754263159, 774410412, 0, 726903076, 705517750,
                   #61, 1, 1057263994, 751504871, 773197589, 0, 727850162, 715943374,
                   #71, 1, 1057268254, 758442129, 746517448, 0, 721989502, 726473951)
        #floats  = (11, 1, 0.5172884464263916, 0.0, 0.0, 0.0, 0.0, 0.0, 2.942726775082116e-44, 1.401298464324817e-45, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.344025239406933e-44, 1.401298464324817e-45, 0.51838618516922, 3.689314004577844e-14, 3.056320862083843e-14, 0.0, 1.7993489727710643e-13, 1.3938290625847838e-12, 5.74532370373175e-44, 1.401298464324817e-45, 0.5174393653869629, 6.842673996865223e-12, 8.397425983741602e-12, 0.0, 1.112146941380232e-12, 4.3178690332414116e-13, 7.146622168056567e-44, 1.401298464324817e-45, 0.5176249146461487, 6.967367221361043e-12, 3.832727191177554e-11, 0.0, 7.519698839303368e-13, 1.2554093506943198e-13, 8.547920632381384e-44, 1.401298464324817e-45, 0.5178447961807251, 5.771150484584764e-12, 3.41194468511663e-11, 0.0, 8.033115188668671e-13, 3.0629529945355727e-13, 9.949219096706201e-44, 1.401298464324817e-45, 0.5180987119674683, 1.0283455510740058e-11, 3.623089675497404e-12, 0.0, 4.856045036569223e-13, 7.287069710669447e-13)

        print('----------------------------------------------')
        print(restype)
        #print(op2.code_information())
        #op2.show_data(data)
        #fmt = mapfmt(b'2i 6f', self.size)
        ntotal = 32 * self.factor  # 32=4*8
        nnodes = ndata // ntotal
        n = ndata
        assert ndata % ntotal == 0, 'ndata=%s ntotal=%s' % (ndata, ntotal)
        if op2.use_vector and is_vectorized:
            #itime = obj.itime
            #n = nnodes * ntotal
            #itotal = obj.itotal
            #itotal2 = itotal + nnodes

            if 1 or obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nnodes, 8)

                nids = ints[:, 0] // 10
                assert nids.min() > 0, nids.min()
                print('nids', nids)
                print('gridtype', ints[:, 1])
                #obj.node_gridtype[itotal:itotal2, 0] = nids
                #obj.node_gridtype[itotal:itotal2, 1] = ints[:, 1].copy()

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nnodes, 8)
            np.set_printoptions(threshold=100000)
            print('floats', floats[:, 2:])
            #obj.data[obj.itime, obj.itotal:itotal2, :] = floats[:, 2:].copy()
            #obj._times[itime] = dt
            #obj.itotal = itotal2
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
            prefix, postfix = get_shock_prefix_postfix(op2.thermal)
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
        elif op2.thermal == 2:  # ABS
            assert op2.table_name in [b'OUPV1'], op2.table_name
            n = op2._read_table_vectorized(data, ndata, result_name, storage_obj,
                                           RealDisplacementArray, ComplexDisplacementArray,
                                           'node', random_code=op2.random_code)
        elif op2.thermal == 4:  # SRSS
            # F:\work\pyNastran\examples\Dropbox\move_tpl\ms103.op2
            assert op2.table_name in [b'OUPV1'], op2.table_name
            n = op2._read_table_vectorized(data, ndata, result_name, storage_obj,
                                           RealDisplacementArray, ComplexDisplacementArray,
                                           'node', random_code=op2.random_code)
        elif op2.thermal == 8:  # NRL
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
            assert op2.thermal in {2, 4, 8}, op2.code_information()
            prefix, postfix = get_shock_prefix_postfix(op2.thermal)
            result_name = prefix + result_name0 + postfix
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
            assert op2.thermal in {2, 4, 8}, op2.code_information() # should 0 be here?
            result_name0 = 'accelerations'
            prefix, postfix = get_shock_prefix_postfix(op2.thermal)
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
                               b'OPHIG', b'BOPHIG', b'OUGV1PAT']:
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
