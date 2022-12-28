#pylint: disable=C0301,W0212,C0103,W0201
"""
Defines the Real/Complex Forces created by:
    FORCE = ALL

NX Case Control  Block         Description
===============  ==========    ===========
FORCE            OEF1          Element forces or heat flux (linear elements only)
FORCE            OEF1X         Element forces (nonlinear elements only)
????             HOEF1         ???
FORCE            DOEF1         Scaled Response Spectra
MODCON           OEFMC         Modal contributions
FORCE            OEF1X         Element forces with intermediate (CBAR and CBEAM)
                               station forces and forces on nonlinear elements
FLUX             HOEF1         Element heat flux

"""
from __future__ import annotations
from struct import Struct
from typing import Union, Any, TYPE_CHECKING

import numpy as np
from numpy import frombuffer, vstack, sin, cos, radians, array, hstack, zeros

from pyNastran.op2.op2_interface.function_codes import func1, func7
from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.tables.utils import get_eid_dt_from_eid_device
from pyNastran.op2.op2_helper import polar_to_real_imag
from pyNastran.op2.op2_interface.utils import apply_mag_phase, reshape_bytes_block_strip
from pyNastran.op2.op2_interface.msc_tables import MSC_OEF_REAL_MAPPER, MSC_OEF_IMAG_MAPPER
from pyNastran.op2.op2_interface.nx_tables import NX_OEF_REAL_MAPPER, NX_OEF_IMAG_MAPPER
from pyNastran.op2.op2_interface.op2_codes import SORT1_TABLES_BYTES, TABLES_BYTES

from pyNastran.op2.tables.oef_forces.oef_thermal_objects import (
    Real1DHeatFluxArray,
    RealHeatFlux_2D_3DArray,
    RealChbdyHeatFluxArray,
    RealConvHeatFluxArray,
)
from pyNastran.op2.tables.oef_forces.oef_force_objects import (
    FailureIndicesArray,
    RealRodForceArray, RealViscForceArray,
    RealCBarForceArray, RealCBar100ForceArray,
    RealCFastForceArrayNX, RealCWeldForceArray,
    RealCFastForceArrayMSC, RealCBushForceArray, RealCBearForceArray,
    RealCWeldForceArrayMSC,
    RealPlateForceArray,
    RealPlateBilinearForceArray,
    RealSpringForceArray, RealDamperForceArray,
    RealCShearForceArray,
    RealCGapForceArray,
    RealConeAxForceArray,
    RealSolidPressureForceArray,
    RealCBeamForceArray,
    RealBendForceArray,
)
from pyNastran.op2.tables.oef_forces.oef_complex_force_objects import (
    ComplexRodForceArray,
    ComplexCBarForceArray, ComplexCWeldForceArray, ComplexCWeldForceArrayMSC,
    ComplexCBeamForceArray,
    ComplexCBushForceArray, ComplexCFastForceArrayMSC,
    ComplexCBearForceArray,
    ComplexCShearForceArray,
    ComplexSpringForceArray,
    ComplexDamperForceArray,
    ComplexViscForceArray,
    ComplexPlateForceArray,
    ComplexPlate2ForceArray,
    ComplexSolidPressureForceArray,
    ComplexCBendForceArray,
)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


class OEF:
    """Defines OEFx table reading for element forces/heat flux"""
    def __init__(self, op2: OP2):
        self.op2 = op2

    @property
    def size(self) -> int:
        return self.op2.size
    @property
    def factor(self) -> int:
        return self.op2.factor

    def _create_oes_object4(self, *args, **kwargs):
        return self.op2._create_oes_object4(*args, **kwargs)

    def get_oef_prefix_postfix(self):
        """
        NX Case Control  Block         Description
        ===============  ==========    ===========
        NLSTRESS         OESNLXR       Nonlinear static stresses
        BOUTPUT          OESNLBR       Slideline stresses
        STRESS           OESNLXD       Nonlinear Transient Stresses
        STRESS           OES1C/OSTR1C  Ply stresses/strains
        STRESS           OES1X         Element stresses with intermediate (CBAR and CBEAM)
                                       station stresses and stresses on nonlinear elements
        STRESS           OES/OESVM     Element stresses (linear elements only)
        STRAIN           OSTR1         Element strains
        STRESS/STRAIN    DOES1/DOSTR1  Scaled Response Spectra
        MODCON           OSTRMC        Modal contributions
        """
        op2 = self.op2
        prefix = ''
        postfix = ''
        table_name_bytes = op2.table_name
        assert isinstance(table_name_bytes, bytes), table_name_bytes
        is_sort1 = table_name_bytes in SORT1_TABLES_BYTES
        assert table_name_bytes in TABLES_BYTES, table_name_bytes

        if table_name_bytes in [b'OEF1X', b'OEF1', b'OEF2']:
            if op2.thermal == 0:
                prefix = 'force.'
            elif op2.thermal == 1:
                prefix = 'thermal_load.'
            else:
                raise NotImplementedError(op2.code_information())
        elif table_name_bytes in [b'HOEF1']:
            postfix = '_flux'
        #elif op2.table_name in ['OESNLXR']:
            #prefix = 'sideline_'
        #elif op2.table_name in ['OESNLXD', 'OESNL1X', 'OESNLXR']:
            #prefix = 'nonlinear_'
        #elif op2.table_name == 'OESNLBR':
            #prefix = 'sideline_'
        #elif op2.table_name == 'OESRT':
            #prefix = 'strength_ratio.'
        #elif op2.table_name in ['OESCP', 'OESTRCP']:
            #pass # TODO: update
        elif table_name_bytes in [b'OEFCRM1', b'OEFCRM2']:
            assert op2.table_code in [4, 504], op2.code_information()
            prefix = 'crm.'
            op2._op2_readers.reader_oes._set_as_random()
        elif table_name_bytes in [b'OEFPSD1', b'OEFPSD2']:
            assert op2.table_code in [4, 604], op2.code_information()
            op2._op2_readers.reader_oes._set_as_random()
            prefix = 'psd.'
        elif table_name_bytes in [b'OEFRMS1', b'OEFRMS2', b'OEFPK1']:
            assert op2.table_code in [4, 404, 804], op2.code_information()
            op2._op2_readers.reader_oes._set_as_random()
            is_sort1 = True
            op2._analysis_code_fmt = b'i'
            prefix = 'rms.'
        elif table_name_bytes in [b'OEFNO1', b'OEFNO2']:
            assert op2.table_code in [4, 904], op2.code_information()
            op2._op2_readers.reader_oes._set_as_random()
            op2.sort_method = 1
            op2.data_code['nonlinear_factor'] = None
            op2._analysis_code_fmt = b'i'
            assert op2.sort_method == 1, op2.code_information()
            prefix = 'no.'
        elif table_name_bytes in [b'OEFATO1', b'OEFATO2']:
            assert op2.table_code in [4], op2.code_information()
            prefix = 'ato.'

        elif table_name_bytes in [b'RAFCONS']:
            prefix = 'RAFCONS.'
        elif table_name_bytes in [b'RAFEATC']:
            prefix = 'RAFEATC.'
        elif table_name_bytes in [b'DOEF1']:
            assert op2.table_code in [4], op2.code_information()
            prefix = shock_response_prefix(op2.thermal)
        elif table_name_bytes in [b'OEFIT']:
            assert op2.table_code in [25], op2.code_information()
            prefix = 'failure_indices.'
            #raise NotImplementedError(op2.code_information())
        elif table_name_bytes in [b'OEFITSTN']: # composite failure indicies
            assert op2.table_code in [25], op2.code_information()
            prefix = 'failure_indices.'
        else:
            raise NotImplementedError('%r' % op2.table_name)

        op2.sort_bits.is_sort1 = is_sort1 # sort2
        op2.sort_method = 1 if is_sort1 else 2
        return prefix, postfix

    def _oef_force_code(self):
        """
        Gets the numwide codes for the element to determine if
        the real or complex result should be found.
        The format and sort codes do not always give the right answer...
        """
        op2 = self.op2
        if op2.is_nx:
            real_mapper = NX_OEF_REAL_MAPPER
            imag_mapper = NX_OEF_IMAG_MAPPER
        else:
            real_mapper = MSC_OEF_REAL_MAPPER
            imag_mapper = MSC_OEF_IMAG_MAPPER


        try:
            real = real_mapper[op2.element_type]
        except KeyError:
            real = None

        try:
            imag = imag_mapper[op2.element_type]
        except KeyError:
            imag = None
        return (real, imag)

    def _read_oef1_3(self, data: bytes, ndata: int):
        """Table 3 parser for OEF1 table"""
        op2 = self.op2
        op2._analysis_code_fmt = b'i'
        op2._data_factor = 1
        op2.words = [
            'aCode', 'tCode', 'element_type', 'isubcase',
            '???', '???', '???', '???',
            'format_code', 'num_wide', 'o_code', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', 'Title', 'subtitle', 'label']

        op2.parse_approach_code(data)

        #: element type
        op2.element_type = op2.add_data_parameter(data, 'element_type', b'i', 3, False)

        # dynamic load set ID/random code
        #self.dLoadID = op2.add_data_parameter(data, 'dLoadID', b'i', 8, False)

        #: format code
        op2.format_code = op2.add_data_parameter(data, 'format_code', b'i', 9, False)

        #: number of words per entry in record
        #: .. note: is this needed for this table ???
        op2.num_wide = op2.add_data_parameter(data, 'num_wide', b'i', 10, False)

        #: undefined in DMAP...
        op2.o_code = op2.add_data_parameter(data, 'o_code', b'i', 11, False)

        #: thermal flag; 1 for heat ransfer, 0 otherwise
        op2.thermal = op2.add_data_parameter(data, 'thermal', b'i', 23, False)

        ## assuming tCode=1
        if op2.analysis_code == 1:   # statics
            op2.loadID = op2.add_data_parameter(data, 'loadID', b'i', 5, False)  # load set ID number
            op2.data_names = op2.apply_data_code_value('data_names', ['loadID'])
            op2.setNullNonlinearFactor()
        elif op2.analysis_code == 2:  # normal modes/buckling (real eigenvalues)
            #: mode number
            op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            #: eigenvalue
            op2.eign = op2.add_data_parameter(data, 'eign', b'f', 6, False)
            op2.cycle = 0.
            op2._op2_readers.reader_oug.update_mode_cycle('cycle')
            op2.data_names = op2.apply_data_code_value('data_names', ['mode', 'eign', 'cycle'])
            # TODO: mode_cycle is not defined?
            #op2.data_names = op2.apply_data_code_value('data_names', ['mode', 'eign', 'mode_cycle'])
        elif op2.analysis_code == 3:  # differential stiffness 0
            #: load set ID number
            op2.loadID = op2.add_data_parameter(data, 'loadID', b'i', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['loadID'])
        elif op2.analysis_code == 4:  # differential stiffness 1
            #: load set ID number
            op2.loadID = op2.add_data_parameter(data, 'loadID', b'i', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['loadID'])
        elif op2.analysis_code == 5:   # frequency
            self.freq = op2.add_data_parameter(data, 'freq', b'f', 5)  # frequency
            op2.data_names = op2.apply_data_code_value('data_names', ['freq'])
        elif op2.analysis_code == 6:  # transient
            self.time = op2.add_data_parameter(data, 'time', b'f', 5)  # time step
            op2.data_names = op2.apply_data_code_value('data_names', ['time'])
        elif op2.analysis_code == 7:  # pre-buckling
            #: load set ID number
            op2.loadID = op2.add_data_parameter(data, 'loadID', b'i', 5)
            #op2.apply_data_code_value('data_names',['lsdvmn'])
            op2.data_names = op2.apply_data_code_value('data_names', ['loadID'])
        elif op2.analysis_code == 8:  # post-buckling
            #: load set ID number
            op2.loadID = op2.add_data_parameter(data, 'loadID', b'i', 5)
            #: real eigenvalue
            op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['loadID', 'eigr'])
        elif op2.analysis_code == 9:  # complex eigenvalues
            #: mode number
            op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            #: real eigenvalue
            op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            #: imaginary eigenvalue
            op2.eigi = op2.add_data_parameter(data, 'eigi', b'f', 7, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['mode', 'eigr', 'eigi'])
        elif op2.analysis_code == 10:  # nonlinear statics
            #: load step
            self.load_step = op2.add_data_parameter(data, 'load_step', b'f', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['load_step'])
        elif op2.analysis_code == 11:  # geometric nonlinear statics
            #: load set ID number
            op2.loadID = op2.add_data_parameter(data, 'loadID', b'i', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['loadID'])
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % str(op2.analysis_code)
            raise RuntimeError(msg)

        op2.fix_format_code()
        op2._parse_thermal_code()
        self._set_force_stress_element_name()

        if op2.is_debug_file:
            op2.binary_debug.write('  %-14s = %r\n' % ('element_name', op2.element_name))
            op2.binary_debug.write('  %-14s = %r %s\n' % ('approach_code', op2.approach_code,
                                                           op2.approach_code_str(op2.approach_code)))
            op2.binary_debug.write('  %-14s = %r\n' % ('tCode', op2.tCode))
            op2.binary_debug.write('  %-14s = %r\n' % ('isubcase', op2.isubcase))


        op2._read_title(data)
        if op2.element_type not in op2.element_mapper:
            msg = 'element_type = %s' % op2.element_type
            return op2._not_implemented_or_skip(data, ndata, msg)
        op2._write_debug_bits()
        assert op2.num_wide != 146, op2.code_information()
        #print('OEF-%s' % op2.element_name)
        #self._check_result_type()

    def _set_force_stress_element_name(self):
        """
        Not all cards can have OES/OEF output, so if they do,
        we have in the wrong solver, specifically:
        - RBAR

        """
        op2 = self.op2
        try:
            op2.element_name = op2.element_mapper[op2.element_type]
        except KeyError:
            op2.log.error(op2.code_information())
            raise

        if op2.element_type == 227 and op2.element_name == 'RBAR' and op2.is_msc:
            op2.to_nx(' because element_type=227 was found')
            op2.element_name = op2.element_mapper[op2.element_type]
        assert op2.element_name != '', op2.code_information()

        op2.data_code['element_name'] = op2.element_name

    def _read_oef2_3(self, data, unused_ndata):
        """Table 3 parser for OEF2 table"""
        op2 = self.op2
        op2._data_factor = 1
        op2.words = [
            'aCode', 'tCode', 'element_type', 'isubcase',
            '???', '???', '???', '???',
            'format_code', 'num_wide', 'o_code', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', 'Title', 'subtitle', 'label']

        op2.parse_approach_code(data)  # 3
        op2.sort_method = 2

        #: element type
        op2.element_type = op2.add_data_parameter(data, 'element_type', b'i', 3, False)

        # dynamic load set ID/random code
        #self.dLoadID = op2.add_data_parameter(data, 'dLoadID', b'i', 8, False)

        #: format code
        op2.format_code = op2.add_data_parameter(data, 'format_code', b'i', 9, False)

        #: number of words per entry in record
        #: .. note: is this needed for this table ???
        op2.num_wide = op2.add_data_parameter(data, 'num_wide', b'i', 10, False)

        #: undefined in DMAP...
        op2.o_code = op2.add_data_parameter(data, 'o_code', b'i', 11, False)

        #: thermal flag; 1 for heat ransfer, 0 otherwise
        op2.thermal = op2.add_data_parameter(data, 'thermal', b'i', 23, False)

        op2.element_id = op2.add_data_parameter(data, 'element_id', b'i', 5, fix_device_code=True)
        op2._element_id = op2.add_data_parameter(data, '_element_id', b'f', 5, apply_nonlinear_factor=False, add_to_dict=True)

        if op2.analysis_code == 1:  # static...because reasons.
            op2._analysis_code_fmt = b'f'
            op2.data_names = op2.apply_data_code_value('data_names', ['element_id'])
        elif op2.analysis_code == 2:  # real eigenvalues
            op2._analysis_code_fmt = b'i'
            op2.eign = op2.add_data_parameter(data, 'eign', b'f', 6, False)
            op2.mode_cycle = op2.add_data_parameter(data, 'mode_cycle', b'i', 7, False)  # mode or cycle .. todo:: confused on the type - F1???
            op2.data_names = op2.apply_data_code_value('data_names', ['element_id', 'eign', 'mode_cycle'])
        elif op2.analysis_code == 5:   # frequency
            op2._analysis_code_fmt = b'f'
            op2.data_names = op2.apply_data_code_value('data_names', ['element_id'])
            op2.apply_data_code_value('analysis_method', 'freq')
        elif op2.analysis_code == 6:  # transient
            op2._analysis_code_fmt = b'f'
            op2.data_names = op2.apply_data_code_value('data_names', ['element_id'])
            op2.apply_data_code_value('analysis_method', 'time')
        elif op2.analysis_code == 7:  # pre-buckling
            op2._analysis_code_fmt = b'i'
            op2.data_names = op2.apply_data_code_value('data_names', ['element_id'])
            op2.apply_data_code_value('analysis_method', 'lsdvmn')
        elif op2.analysis_code == 8:  # post-buckling
            op2._analysis_code_fmt = b'f'
            op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['element_id', 'eigr'])
            op2.apply_data_code_value('analysis_method', 'eigr')
        elif op2.analysis_code == 9:  # complex eigenvalues
            # mode number
            op2._analysis_code_fmt = b'i'
            op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            op2.eigi = op2.add_data_parameter(data, 'eigi', b'f', 7, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['element_id', 'eigr', 'eigi'])
            op2.apply_data_code_value('analysis_method', 'mode')
        elif op2.analysis_code == 10:  # nonlinear statics
            # load step
            op2._analysis_code_fmt = b'f'
            op2.data_names = op2.apply_data_code_value('data_names', ['element_id'])
            op2.apply_data_code_value('analysis_method', 'lftsfq')
        elif op2.analysis_code == 11:  # old geometric nonlinear statics
            # load set number
            op2.data_names = op2.apply_data_code_value('data_names', ['element_id'])
        elif op2.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            op2.data_names = op2.apply_data_code_value('data_names', ['element_id'])
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % op2.analysis_code
            raise RuntimeError(msg)

        op2.fix_format_code()
        op2._op2_readers.reader_oes._fix_oes_sort2(data)
        self._set_force_stress_element_name()
        #assert isinstance(op2.nonlinear_factor, int), op2.nonlinear_factor
        #self._check_result_type()

    def _read_oef1_4(self, data: bytes, ndata: int):
        """Table 4 parser for OEF1 table"""
        op2 = self.op2
        if op2.thermal == 0:
            op2._setup_op2_subcase('FORCE')
            n = self._read_oef1_loads(data, ndata)
        elif op2.thermal == 1:
            n = self._read_oef1_thermal(data, ndata)
        elif op2.thermal in [2, 4, 8]: # 2=ABS, 4=SRSS, 8=NRL
            #C:\NASA\m4\formats\git\examples\move_tpl\ms103.op2 # SRSS
            n = self._read_oef1_loads(data, ndata)
        else:
            n = op2._not_implemented_or_skip(data, ndata, 'thermal=%s' % op2.thermal)
        return n

    def _read_oef2_4(self, data: bytes, ndata: int):
        op2 = self.op2
        if op2.thermal == 0: # and op2.element_type not in [77]:
            op2._setup_op2_subcase('FORCE')
            n = self._read_oef1_loads(data, ndata)
        else:
            n = op2._not_implemented_or_skip(data, ndata, 'thermal=%s' % op2.thermal)
        return n

    def _read_oef1_thermal(self, data: bytes, ndata: int):
        """Table 4 parser for OEF1 thermal table"""
        op2 = self.op2
        if op2._results.is_not_saved('element_forces'):
            return ndata
        prefix, postfix = self.get_oef_prefix_postfix()

        n = 0
        #thermal
        #is_magnitude_phase = op2.is_magnitude_phase()
        dt = op2.nonlinear_factor

        #flag = 'element_id'
        if op2.element_type in [1, 2, 3, 10, 34, 69]:  # ROD,BEAM,TUBE,CONROD,BAR,BEND
            n, nelements, ntotal = self._thermal_1d_heat_flux(data, ndata, dt, prefix, postfix)
            if nelements is None:
                return n

        elif op2.element_type in [33, 53, 64, 74, 75,  # CQUAD4, CTRIAX6, CQUAD8, CTRIA3, CTRIA6
                                   39, 67, 68]:  # TETRA, HEXA, PENTA
            n, nelements, ntotal = self._thermal_2d_3d_heat_flux(data, ndata, dt, prefix, postfix)

        elif op2.element_type in [107, 108, 109]:  # CHBDYE, CHBDYG, CHBDYP
            n, nelements, ntotal = self._thermal_chbdy(data, ndata, dt, prefix, postfix)

        elif op2.element_type == 110:  # CONV
            n, nelements, ntotal = self._thermal_conv(data, ndata, dt, prefix, postfix)

        elif op2.element_type in [145, 146, 147,  # VUHEXA,VUPENTA,VUTETRA
                                   189, 190,  # VUQUAD,VUTRIA
                                   191]:  # VUBEAM
            # removed by msc/nx
            msg = f'{op2.element_name}-{op2.element_type} has been removed'
            return op2._not_implemented_or_skip(data, ndata, msg)
        else:
            msg = 'OEF sort1 thermal Type=%s num=%s' % (op2.element_name, op2.element_type)
            return op2._not_implemented_or_skip(data, ndata, msg)
        if nelements is None:
            return n


        assert op2.thermal == 1, op2.thermal
        assert ndata > 0, ndata
        try:
            assert nelements > 0, 'nelements=%r element_type=%s element_name=%r\n%s' % (nelements, op2.element_type, op2.element_name, op2.code_information())
        except UnboundLocalError:
            raise UnboundLocalError('element_name=%r' % op2.element_name)
        #assert ndata % ntotal == 0, '%s n=%s nwide=%s len=%s ntotal=%s' % (op2.element_name, ndata % ntotal, ndata % op2.num_wide, ndata, ntotal)
        assert op2.num_wide * 4 == ntotal, 'numwide*4=%s ntotal=%s' % (op2.num_wide*4, ntotal)
        assert n > 0, 'n=%s element_type=%s element_name=%s numwide=%s' % (
            n, op2.element_type, op2.element_name, op2.num_wide)
        return n

    def _thermal_1d_heat_flux(self, data, ndata, dt, prefix, postfix):
        """
        1-CROD
        2-CBEAM
        3-CTUBE
        10-CONROD
        34-CBAR
        69-CBEND
        """
        op2 = self.op2
        n = 0
        obj_vector_real = Real1DHeatFluxArray
        #if op2.element_type == 1: # CROD
        if op2.element_type == 1:
            result_name = prefix + 'crod_thermal_load' + postfix
        elif op2.element_type == 2:
            result_name = prefix + 'cbeam_thermal_load' + postfix
        elif op2.element_type == 3:
            result_name = prefix + 'ctube_thermal_load' + postfix
        elif op2.element_type == 10:
            result_name = prefix + 'conrod_thermal_load' + postfix
        elif op2.element_type == 34:
            result_name = prefix + 'cbar_thermal_load' + postfix
        elif op2.element_type == 69:
            result_name = prefix + 'cbend_thermal_load' + postfix
        else:
            raise NotImplementedError('element_type=%s element_name=%s' % (
                op2.element_type, op2.element_name))

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        if op2.format_code == 1 and op2.num_wide == 9:  # real
            ntotal = 36
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return nelements * op2.num_wide * 4, None, None
            obj = op2.obj
            #if op2.is_debug_file:
                #op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                #op2.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
                #op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * 4 * op2.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 9)
                obj._times[obj.itime] = dt

                strings = frombuffer(data, dtype=op2._uendian + 'S4').reshape(nelements, 9)
                s = array([s1+s2 for s1, s2 in zip(strings[:, 1], strings[:, 2])])
                #print(s)
                #print('itime = ', obj.itime)
                #print('---------')
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, 9)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids
                    obj.element_data_type[itotal:itotal2] = s
                    #obj.element_type[obj.itime, itotal:itotal2, :] = strings[:, 3:]

                #[etype, xgrad, ygrad, zgrad, xflux, yflux, zflux]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 3:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                s = Struct(op2._endian + op2._analysis_code_fmt + b'8s6f')
                add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]
                    out = s.unpack(edata)
                    (eid_device, etype, xgrad, ygrad, zgrad, xflux, yflux, zflux) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, op2.nonlinear_factor, op2.sort_method)
                    add_sort_x(dt, eid, etype, xgrad, ygrad, zgrad, xflux, yflux, zflux)
                    n += ntotal
        else:  # pragma: no cover
            msg = op2.code_information()
            return op2._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _thermal_2d_3d_heat_flux(self, data, ndata, dt, prefix, postfix):
        """
        33-QUAD4-centroidal
        53-TRIAX6
        64-QUAD8
        74-TRIA3
        75-TRIA6

        39-TETRA
        67-HEXA
        68-PENTA
        """
        op2 = self.op2
        n = 0
        if op2.element_type == 33:
            result_name = prefix + 'cquad4_thermal_load' + postfix
        elif op2.element_type == 53:
            result_name = prefix + 'ctriax6_thermal_load' + postfix
        elif op2.element_type == 64:
            result_name = prefix + 'cquad8_thermal_load' + postfix
        elif op2.element_type == 74:
            result_name = prefix + 'ctria3_thermal_load' + postfix
        elif op2.element_type == 75:
            result_name = prefix + 'ctria6_thermal_load' + postfix
        elif op2.element_type == 39:
            result_name = prefix + 'ctetra_thermal_load' + postfix
        elif op2.element_type == 67:
            result_name = prefix + 'chexa_thermal_load' + postfix
        elif op2.element_type == 68:
            result_name = prefix + 'cpenta_thermal_load' + postfix
        else:
            raise NotImplementedError('element_type=%s element_name=%s' % (
                op2.element_type, op2.element_name))

        if op2._results.is_not_saved(result_name):
            return ndata, None, None

        obj_vector_real = RealHeatFlux_2D_3DArray
        #if op2.element_type == 1: # CROD
        #result_name = 'thermalLoad_2D_3D'

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)
        if op2.format_code == 1 and op2.num_wide == 9:  # real - 2D
            # [33, 53, 64, 74, 75]
            ntotal = 4 * op2.num_wide * self.factor
            ntotal = 36 * self.factor
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return nelements * ntotal, None, None
            obj = op2.obj
            #if op2.is_debug_file:
                #op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                #op2.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
                #op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 9)
                obj._times[obj.itime] = dt
                #if obj.itime == 0:
                ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, 9)
                eids = ints[:, 0] // 10
                assert eids.min() > 0, eids.min()
                obj.element[itotal:itotal2] = eids
                strings = frombuffer(data, dtype=op2._uendian + 'S4').reshape(nelements, 9)
                obj.element_data_type[itotal:itotal2] = array([s1+s2 for s1, s2 in zip(strings[:, 1], strings[:, 2])])

                #[etype, xgrad, ygrad, zgrad, xflux, yflux, zflux]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 3:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                # no zed on this element for some reason...
                if self.size == 4:
                    fmt = op2._endian + op2._analysis_code_fmt + b'8s 6f'
                else:
                    fmt = op2._endian + mapfmt(op2._analysis_code_fmt, self.size) + b'16s 6d'
                s = Struct(fmt)
                add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]
                    n += ntotal
                    out = s.unpack(edata)
                    (eid_device, etype, xgrad, ygrad, zgrad, xflux, yflux, zflux) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, op2.nonlinear_factor, op2.sort_method)
                    add_sort_x(dt, eid, etype, xgrad, ygrad, zgrad, xflux, yflux, zflux)

        elif op2.format_code == 1 and op2.num_wide == 10:  # real - 3D
            # [39, 67, 68]:  # HEXA,PENTA
            ntotal = 40 * self.factor
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return nelements * ntotal, None, None
            obj = op2.obj
            assert nelements > 0, 'ndata=%s ntotal=%s' % (ndata, ntotal)
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 10)
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, 10)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids
                    strings = frombuffer(data, dtype=op2._uendian + 'S4').reshape(nelements, 10)
                    obj.element_data_type[itotal:itotal2] = array([s1+s2 for s1, s2 in zip(strings[:, 1], strings[:, 2])])

                #[etype, xgrad, ygrad, zgrad, xflux, yflux, zflux, zed]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 3:-1].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                s = Struct(op2._endian + op2._analysis_code_fmt + b'8s6fi')
                add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]
                    n += ntotal
                    out = s.unpack(edata)
                    (eid_device, etype, xgrad, ygrad, zgrad, xflux, yflux, zflux, unused_zed) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, op2.nonlinear_factor, op2.sort_method)
                    add_sort_x(dt, eid, etype, xgrad, ygrad, zgrad, xflux, yflux, zflux)
        else:
            raise RuntimeError(op2.code_information())
        return n, nelements, ntotal

    def _thermal_chbdy(self, data, ndata, dt, prefix, postfix):
        """
        107-CHBDYE
        108-CHBDYG
        109-CHBDYP

        """
        op2 = self.op2
        n = 0
        #if op2.table_name in ['OEF1X']:
        if op2.element_type == 107:
            result_name = prefix + 'chbdye_thermal_load' + postfix
        elif op2.element_type == 108:
            result_name = prefix + 'chbdyg_thermal_load' + postfix
        elif op2.element_type == 109:
            result_name = prefix + 'chbdyp_thermal_load' + postfix
        else:
            raise NotImplementedError('element_type=%s element_name=%s' % (
                op2.element_type, op2.element_name))

        if op2._results.is_not_saved(result_name):
            return ndata, None, None

        #elif op2.table_name in ['HOEF1']:
            #if op2.element_type == 107:
                #result_name = 'chbdye_thermal_flux'
            #elif op2.element_type == 108:
                #result_name = 'chbdyg_thermal_flux'
            #elif op2.element_type == 109:
                #result_name = 'chbdyp_thermal_flux'
            #else:
                #raise NotImplementedError('element_type=%s element_name=%s' % (
                    #op2.element_type, op2.element_name))
        #else:
            #raise NotImplementedError(msg)

        if op2.format_code == 1 and op2.num_wide == 8:  # real
            #result_name = 'thermalLoad_CHBDY'
            if op2._results.is_not_saved(result_name):
                return ndata, None, None
            op2._results._found_result(result_name)
            slot = op2.get_result(result_name)

            if op2.format_code == 1 and op2.num_wide == 8:  # real
                obj_vector_real = RealChbdyHeatFluxArray
                ntotal = 32
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    return nelements * op2.num_wide * 4, None, None
                obj = op2.obj
                #if op2.is_debug_file:
                    #op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    #op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    #op2.binary_debug.write('  #elementi = [eid_device, etype, fapplied, free_conv, force_conv, frad, ftotal]\n')
                    #op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                if op2.use_vector and is_vectorized and op2.sort_method == 1:
                    n = nelements * 4 * op2.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 8)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, 8)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids
                        #obj.element_type[obj.itime, itotal:itotal2, :] = strings[:, 3:]

                    #[fapplied, free_conv, force_conv, frad, ftotal]
                    obj.data[obj.itime, itotal:itotal2, :] = floats[:, 3:].copy()
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s1 = Struct(op2._endian + op2._analysis_code_fmt + b'8s5f')
                    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
                    for unused_i in range(nelements):
                        edata = data[n:n+32]
                        n += ntotal
                        out = s1.unpack(edata)
                        (eid_device, etype, fapplied, free_conv, force_conv, frad, ftotal) = out
                        eid, dt = get_eid_dt_from_eid_device(
                            eid_device, op2.nonlinear_factor, op2.sort_method)

                        if op2.is_debug_file:
                            op2.binary_debug.write('  %s -> [%s, %s, %s, %s, %s, %s, %s]\n'
                                                    % (eid, eid_device, etype, fapplied, free_conv, force_conv, frad, ftotal))
                        add_sort_x(dt, eid, etype, fapplied, free_conv, force_conv, frad, ftotal)
        else:  # pragma: no cover
            msg = op2.code_information()
            return op2._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _thermal_conv(self, data, ndata, dt, prefix, postfix):
        op2 = self.op2
        n = 0
        # 110-CONV
        result_name = prefix + 'conv_thermal_load' + postfix
        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        if op2.format_code == 1 and op2.num_wide == 4:
            ntotal = 16
            nelements = ndata // ntotal

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, RealConvHeatFluxArray)
            if auto_return:
                return nelements * op2.num_wide * 4, None, None
            obj = op2.obj
            #if op2.is_debug_file:
                #op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                #op2.binary_debug.write('  #elementi = [eid_device, etype, fapplied, free_conv, force_conv, frad, ftotal]\n')
                #op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * 4 * op2.num_wide
                ielement = obj.ielement
                ielement2 = ielement + nelements

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 4).copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, 4).copy()
                    eids = ints[:, 0] // 10
                    nids = ints[:, 2]
                    assert eids.min() > 0, eids.min()
                    assert nids.min() >= 0, nids.min()
                    obj.element_node[ielement:ielement2, 0] = eids
                    obj.element_node[ielement:ielement2, 1] = nids

                #[eid, free_conv, cntl_node, free_conv_k]
                obj.data[obj.itime, ielement:ielement2, :] = floats[:, [1, 3]]
                obj.itotal = ielement2
                obj.ielement = ielement2
            else:
                s1 = Struct(op2._endian + op2._analysis_code_fmt + b'fif')
                add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
                for unused_i in range(nelements):
                    edata = data[n:n+16]
                    n += 16
                    out = s1.unpack(edata)
                    (eid_device, free_conv, cntl_node, free_conv_k) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, op2.nonlinear_factor, op2.sort_method)
                    assert cntl_node >= 0, cntl_node
                    add_sort_x(dt, eid, cntl_node, free_conv, free_conv_k)
        else:  # pragma: no cover
            msg = op2.code_information()
            return op2._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal


    def _print_obj_name_on_crash(func):
        """
        Decorator debugging function to print the object name and an needed parameters
        """
        def new_func(self, data):
            """
            The actual function exec'd by the decorated function.
            """
            try:
                n = func(self, data)
            except Exception:
                raise
                #print("----------")
                #try:
                    #print(op2.obj)
                #except Exception:
                    #print("error printing %r" % op2.obj.__class__.__name__)
                #print(op2.data_code)
                #if op2.obj is not None:
                    ##from pyNastran.utils import object_attributes
                    ##print object_attributes(op2.obj)
                    #print(op2.obj.data_code)
                #print("----------")
                #raise
            return n
        return new_func

    def _read_oef1_loads_nasa95(self, data: bytes, ndata: int):
        """Reads the OEF1 table for NASA 95 Nastran"""
        op2 = self.op2
        if op2._results.is_not_saved('element_forces'):
            return ndata

        prefix, postfix = self.get_oef_prefix_postfix()
        _sort_method = func1(op2.tCode)
        result_type = op2.result_type # func7(op2.tCode)

        #print('prefix=%r postfix=%s element_name=%s' % (prefix, postfix, op2.element_name))
        (num_wide_real, num_wide_imag) = self._oef_force_code()
        if op2.is_debug_file:
            op2.binary_debug.write('  num_wide_real = %r\n' % num_wide_real)
            op2.binary_debug.write('  num_wide_imag = %r\n' % num_wide_imag)

        n = 0
        is_magnitude_phase = op2.is_magnitude_phase()
        dt = op2.nonlinear_factor

        if op2.element_type in [1, 3, 10]:  # rods
            n, nelements, ntotal = self._oef_crod(data, ndata, dt, is_magnitude_phase,
                                                  result_type, prefix, postfix)

        elif op2.element_type == 2:  # cbeam
            #2-CBEAM
            n, nelements, ntotal = self._oef_cbeam(data, ndata, dt, is_magnitude_phase,
                                                   result_type, prefix, postfix)

        elif op2.element_type in [11, 12, 13, 14,   # springs
                                   20, 21, 22, 23]:  # dampers
            # 11-CELAS1
            # 12-CELAS2
            # 13-CELAS3
            # 14-CELAS4

            # 20-CDAMP1
            # 21-CDAMP2
            # 22-CDAMP3
            # 23-CDAMP4
            n, nelements, ntotal = self._oef_celas_cdamp(data, ndata, dt, is_magnitude_phase,
                                                         result_type, prefix, postfix)

        elif op2.element_type == 24:  # CVISC
            n, nelements, ntotal = self._oef_cvisc(data, ndata, dt, is_magnitude_phase,
                                                   result_type, prefix, postfix)

        elif op2.element_type == 34:  # cbar
            # 34-CBAR
            n, nelements, ntotal = self._oef_cbar_34(data, ndata, dt, is_magnitude_phase,
                                                     result_type, prefix, postfix)

        elif op2.element_type in [83]: # centroidal shells
            # 33-CQUAD4???
            # 83-CTRIA3
            n, nelements, ntotal = self._oef_shells_centroidal(data, ndata, dt, is_magnitude_phase,
                                                               result_type, prefix, postfix)
        elif op2.element_type == 4:  # cshear
            n, nelements, ntotal = self._oef_cshear(data, ndata, dt, is_magnitude_phase,
                                                    result_type, prefix, postfix)

        elif op2.element_type == 35:  # cconeax
            n, nelements, ntotal = self._oef_cconeax(data, ndata, dt, is_magnitude_phase,
                                                     result_type, prefix, postfix)
        else:
            return op2._not_implemented_or_skip(data, ndata, op2.code_information())

        if nelements is None:
            return n

        #assert op2.thermal == 0, op2.thermal
        assert ndata > 0, ndata
        assert nelements > 0, 'nelements=%r element_type=%s element_name=%r num_wide=%s' % (
            nelements, op2.element_type, op2.element_name, op2.num_wide)
        #assert ndata % ntotal == 0, '%s n=%s nwide=%s len=%s ntotal=%s' % (op2.element_name, ndata % ntotal, ndata % op2.num_wide, ndata, ntotal)
        assert op2.num_wide * 4 == ntotal, 'numwide*4=%s ntotal=%s' % (op2.num_wide*4, ntotal)
        assert n > 0, n
        return n

    # @_print_obj_name_on_crash
    def _read_oef1_loads(self, data: bytes, ndata: int):
        """Reads the OEF1 table; stores the element forces/heat flux."""
        op2 = self.op2
        #self._apply_oef_ato_crm_psd_rms_no('') # TODO: just testing
        if op2._results.is_not_saved('element_forces'):
            return ndata

        prefix, postfix = self.get_oef_prefix_postfix()
        _sort_method = func1(op2.tCode)
        result_type = op2.result_type # func7(op2.tCode)

        #print('prefix=%r postfix=%s element_name=%s' % (prefix, postfix, op2.element_name))
        unused_flag = 'element_id'
        (num_wide_real, num_wide_imag) = self._oef_force_code()
        if op2.is_debug_file:
            op2.binary_debug.write(f'  num_wide_real = {num_wide_real!r}\n')
            op2.binary_debug.write(f'  num_wide_imag = {num_wide_imag!r}\n')

        n = 0
        is_magnitude_phase = op2.is_magnitude_phase()
        dt = op2.nonlinear_factor

        element_type = op2.element_type
        if element_type in [1, 3, 10]:  # rods
            n, nelements, ntotal = self._oef_crod(data, ndata, dt, is_magnitude_phase,
                                                  result_type, prefix, postfix)

        elif element_type == 2:  # cbeam
            #2-CBEAM
            n, nelements, ntotal = self._oef_cbeam(data, ndata, dt, is_magnitude_phase,
                                                   result_type, prefix, postfix)

        elif element_type in [11, 12, 13, 14,   # springs
                                  20, 21, 22, 23]:  # dampers
            # 11-CELAS1
            # 12-CELAS2
            # 13-CELAS3
            # 14-CELAS4

            # 20-CDAMP1
            # 21-CDAMP2
            # 22-CDAMP3
            # 23-CDAMP4
            n, nelements, ntotal = self._oef_celas_cdamp(data, ndata, dt, is_magnitude_phase,
                                                         result_type, prefix, postfix)

        elif element_type == 24:  # CVISC
            n, nelements, ntotal = self._oef_cvisc(data, ndata, dt, is_magnitude_phase,
                                                   result_type, prefix, postfix)

        elif element_type == 34:  # cbar
            # 34-CBAR
            n, nelements, ntotal = self._oef_cbar_34(data, ndata, dt, is_magnitude_phase,
                                                     result_type, prefix, postfix)

        elif element_type == 100:  # cbar
            #100-BARS
            n, nelements, ntotal = self._oef_cbar_100(data, ndata, dt, is_magnitude_phase,
                                                      result_type, prefix, postfix)

        elif element_type in [33, 74]:  # centroidal shells
            # 33-CQUAD4
            # 74-CTRIA3
            n, nelements, ntotal = self._oef_shells_centroidal(data, ndata, dt, is_magnitude_phase,
                                                               result_type, prefix, postfix)
        elif op2.is_nx and element_type in [227, 228]:  # centroidal shells
            # 227-CTRIAR? (C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\cqrdbx102.op2)
            # 228-CQUADR
            n, nelements, ntotal = self._oef_shells_centroidal(data, ndata, dt, is_magnitude_phase,
                                                               result_type, prefix, postfix)

        elif element_type in [64, 70, 75, 82, 144]: # bilinear shells
            # 64-CQUAD8
            # 70-CTRIAR
            # 75-CTRIA6
            # 82-CQUADR
            # 144-CQUAD4-bilinear
            n, nelements, ntotal = self._oef_shells_nodal(data, ndata, dt, is_magnitude_phase,
                                                          result_type, prefix, postfix)

        elif element_type in [95, 96, 97, 98]: # composites
            # 95 - CQUAD4
            # 96 - CQUAD8
            # 97 - CTRIA3
            # 98 - CTRIA6 (composite)
            n, nelements, ntotal = self._oef_shells_composite(data, ndata, dt, is_magnitude_phase,
                                                              result_type, prefix, postfix)
        elif op2.is_nx and element_type in [232, 233]: # composites
            # 232 - CQUADR
            # 233 - CTRIAR
            n, nelements, ntotal = self._oef_shells_composite(data, ndata, dt, is_magnitude_phase,
                                                              result_type, prefix, postfix)

        elif element_type in [39, 67, 68]: # solids
            # 39-CTETRA
            # 67-CHEXA
            # 68-CPENTA
            if op2.read_mode == 1:
                return ndata
            #op2._results._found_result('solid_forces')
            raise RuntimeError(op2.code_information())
            #if op2.format_code == 1 and op2.num_wide == 0:  # real
                ##self.create_transient_object(result_name, self.solidForces, RealCSolidForce)
                #raise RuntimeError(op2.code_information())
            #else:
                #msg = op2.code_information()
                #return op2._not_implemented_or_skip(data, ndata, msg)

        elif element_type == 53:  # ctriax6
            # 53-CTRIAX6
            op2._results._found_result('ctriax_force')
            #if op2.format_code == 1 and op2.num_wide == 0:  # real
                #pass
                #self.create_transient_object(self.ctriax_force, RealCTriaxForce)  # undefined
            #else:  # pragma: no cover
            raise NotImplementedError(op2.code_information())
                #msg = op2.code_information()
                #return op2._not_implemented_or_skip(data, ndata, msg)
            #return ndata

        elif element_type == 4:  # cshear
            n, nelements, ntotal = self._oef_cshear(data, ndata, dt, is_magnitude_phase,
                                                    result_type, prefix, postfix)

        elif element_type == 35:  # coneax
            n, nelements, ntotal = self._oef_cconeax(data, ndata, dt, is_magnitude_phase,
                                                     result_type, prefix, postfix)

        elif element_type == 38:  # cgap
            n, nelements, ntotal = self._oef_cgap(data, ndata, dt, is_magnitude_phase,
                                                  result_type, prefix, postfix)

        elif element_type == 69:  # cbend
            n, nelements, ntotal = self._oef_cbend(data, ndata, dt, is_magnitude_phase,
                                                   result_type, prefix, postfix)

        elif element_type in [76, 77, 78, 79]:
            # 76-HEXPR
            # 77-PENPR
            # 78-TETPR
            # 79-CPYRAM
            n, nelements, ntotal = self._oef_csolid_pressure(data, ndata, dt, is_magnitude_phase,
                                                             result_type, prefix, postfix)

        elif element_type in [102, 280]:
            # 102: cbush
            # 280: cbear
            n, nelements, ntotal = self._oef_cbush(data, ndata, dt, is_magnitude_phase,
                                                   result_type, prefix, postfix)

        elif element_type in [145, 146, 147]:
            # 145-VUHEXA
            # 146-VUPENTA
            # 147-VUTETRA
            #if op2.read_mode == 1:
                #return ndata
            return ndata

        elif element_type == 126 and op2.is_msc:
            # 119-CFAST-MSC
            n, nelements, ntotal = self._oef_cbush(data, ndata, dt, is_magnitude_phase,
                                                   result_type, prefix, postfix)
        elif element_type == 119 and op2.is_nx:
            # 119-CFAST-NX
            n, nelements, ntotal = self._oef_cbar_34(data, ndata, dt, is_magnitude_phase,
                                                     result_type, prefix, postfix)
        elif element_type in [117, 200]:
            # 117-CWELDC
            # 200-CWELD
            n, nelements, ntotal = self._oef_cbar_34(data, ndata, dt, is_magnitude_phase,
                                                     result_type, prefix, postfix)
        #elif element_type == 119 and op2.is_msc:
            #raise NotImplementedError(op2.code_information())
        elif element_type == 235:
            # 235-CQUADR
            return op2._not_implemented_or_skip(data, ndata, op2.code_information())
        elif element_type in [189, 190,  # VUQUAD, VTRIA - order?
                                   191]: # VUBEAM
            # removed from msc/nx
            #n, nelements, ntotal = self._oef_vu_shell(data, ndata, dt, is_magnitude_phase,
                                                      #result_type, prefix, postfix)
            #n, nelements, ntotal = self._oef_vu_beam(data, ndata, dt, is_magnitude_phase,
                                                     #result_type, prefix, postfix)
            msg = f'{op2.element_name}-{element_type} has been removed'
            return op2._not_implemented_or_skip(data, ndata, msg)
        elif op2.is_nx:
            if element_type in [118, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352,
                                     356, 357, 363]:
                # 118 CWELDP
                # 343 CTRIA6 SOL 401
                #344 CQUAD8 SOL 401
                #345 CTRIAR SOL 401
                # 347 CBAR SOL 401
                # 348 CBEAM SOL 401
                #346 CQUADR SOL 401
                # 349 CBUSH1D SOL 401
                # 350 CELAS1 SOL 401
                # 351 CELAS2 SOL 401
                # 352 CBUSH SOL 401
                #356 Composite quadrilateral shell element (CQUAD8); SOL 402?
                #357 Composite triangular shell element (CTRIAR); SOL 402?
                # 363 CROD SOL 402
                return op2._not_implemented_or_skip(data, ndata, op2.code_information())
            #print(op2.code_information())
            #nx_missing
            return op2._not_implemented_or_skip(data, ndata, op2.code_information())
        else:
            if element_type == 118:
                # 118-WELDP
                n, nelements, ntotal = self._oef_cbar_34(
                    data, ndata, dt, is_magnitude_phase,
                    result_type, prefix, postfix)
            else:
                #print(op2.code_information())
                #msc_missing
                return op2._not_implemented_or_skip(data, ndata, op2.code_information())

        if nelements is None:
            return n
        if element_type != 2:
            # CBEAM-2: has a finalize step
            op2._op2_readers.reader_oes.check_element_ids()

        #assert op2.thermal == 0, op2.thermal
        assert ndata > 0, ndata
        assert nelements > 0, 'nelements=%r element_type=%s element_name=%r num_wide=%s' % (
            nelements, op2.element_type, op2.element_name, op2.num_wide)
        #assert ndata % ntotal == 0, '%s n=%s nwide=%s len=%s ntotal=%s' % (op2.element_name, ndata % ntotal, ndata % op2.num_wide, ndata, ntotal)
        assert op2.num_wide * 4 * self.factor == ntotal, f'numwide*4={op2.num_wide*4} ntotal={ntotal}'
        assert n is not None and n > 0, op2.code_information()
        return n

    def _oef_crod(self, data, ndata, dt, is_magnitude_phase,
                  result_type, prefix, postfix):
        """
        1-CROD
        3-CTUBE
        10-CONROD

        """
        op2 = self.op2
        n = 0
        obj_real = RealRodForceArray
        obj_complex = ComplexRodForceArray
        if op2.element_type == 1: # CROD
            result_name = prefix + 'crod_force' + postfix
        elif op2.element_type == 3:  # CTUBE
            result_name = prefix + 'ctube_force' + postfix
        elif op2.element_type == 10:  # CONROD
            result_name = prefix + 'conrod_force' + postfix
        else:
            raise NotImplementedError(op2.code_information())
            #msg = 'sort1 Type=%s num=%s' % (op2.element_name, op2.element_type)
            #return op2._not_implemented_or_skip(data, ndata, msg)

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)

        slot = op2.get_result(result_name)
        if op2.format_code == 1 and op2.num_wide == 3: # real
            ntotal = 3 * 4 * self.factor
            nelements = ndata // ntotal

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_real)
            if auto_return:
                return nelements * op2.num_wide * 4, None, None
            obj = op2.obj
            if op2.is_debug_file:
                op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                op2.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
                op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 3)
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype8).reshape(nelements, 3)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                #[axial, torsion]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                n = oef_crod_real_3(op2, data, obj,
                                    nelements, ntotal)

        elif op2.format_code in [2, 3] and op2.num_wide == 5: # imag
            ntotal = 20 * self.factor  # 5*4
            nelements = ndata // ntotal

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_complex)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.is_debug_file:
                op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                op2.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
                op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 5).copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype8).reshape(nelements, 5)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                #[axial_force, torque]
                #(eid_device, axial_real, torque_real, axial_imag, torque_imag) = out
                real_imag = apply_mag_phase(floats, is_magnitude_phase, [1, 2], [3, 4])
                obj.data[obj.itime, itotal:itotal2, :] = real_imag
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                n = oef_crod_imag_5(op2, data, obj,
                                    nelements, ntotal,
                                    is_magnitude_phase)

        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
            #msg = op2.code_information()
            #print(msg)
            #return op2._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oef_cbeam(self, data: bytes, ndata: int, dt, is_magnitude_phase: bool,
                   result_type: str, prefix: str, postfix: str) -> int:
        """2-CBEAM"""
        op2 = self.op2
        n = 0
        result_name = prefix + 'cbeam_force' + postfix
        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)

        #if op2.format_code == 1 and op2.num_wide == 9:  # real centroid ???
            #raise RuntimeError('is this used?')
            #auto_return, is_vectorized = self._create_oes_object4(
                #nelements, result_name, slot, RealCBeamForceArray)
            #if auto_return:
                #return nelements * op2.num_wide * 4

            #obj = op2.obj
            ##is_vectorized = False
            #if op2.use_vector and is_vectorized and op2.sort_method == 1:
                #n = nelements * 4 * op2.num_wide
                #itotal = obj.itotal
                #itotal2 = obj.itotal + nelements
                #ielement2 = obj.ielement + nelements

                #floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 9)[:, 1:]
                #obj._times[obj.itime] = dt
                #if obj.itime == 0:
                    #ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, 9)
                    #eids = ints[:, 0] // 10
                    #assert eids.min() > 0, eids.min()
                    #assert 0 not in eids, eids

                    #obj.element[itotal:itotal2] = eids
                    #obj.element_node[itotal:itotal2, 0] = eids
                    ##obj.element_node[itotal:itotal2, 1] = nids

                ##[sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq]
                #obj.data[obj.itime, itotal:itotal2, :] = floats.copy()
                #obj.itotal = itotal2
                #obj.ielement = ielement2
            #else:
                #s = Struct(op2._endian + b'i8f')  # 36
                #ntotal = 36
                #nelements = ndata // ntotal
                #obj = op2.obj
                #for i in range(nelements):
                    #edata = data[n:n+36]
                    #out = s.unpack(edata)
                    #if op2.is_debug_file:
                        #op2.binary_debug.write('OEF_Beam - %s\n' % (str(out)))
                    #(eid_device, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq) = out
                    #eid, dt = get_eid_dt_from_eid_device(
                        #eid_device, op2.nonlinear_factor, op2.sort_method)
                    #n += 36

        if result_type in [0, 2] and op2.num_wide == 100:  # real/random
            if op2.sort_method == 2:
                msg = op2.code_information()
                if op2.read_mode == 2:
                    return op2._not_implemented_or_skip(data, ndata, msg), None, None
                return ndata, None, None
            # real - format_code == 1
            # random - format_code == 2
            #result_name, is_random = self._apply_oef_ato_crm_psd_rms_no(result_name)
            slot = op2.get_result(result_name)
            ntotal = 400 * self.factor  # 1+(10-1)*11=100 ->100*4 = 400
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, RealCBeamForceArray)
            if auto_return:
                op2._data_factor = 11
                return nelements * op2.num_wide * 4, None, None
            obj = op2.obj

            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.itotal
                itotal2 = obj.itotal + nelements * 11
                #ielement = obj.ielement
                ielement2 = obj.ielement + nelements

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 100)[:, 1:]
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype8).reshape(nelements, 100)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    assert 0 not in eids, eids
                    eids2 = np.repeat(eids, 11)

                    ints2 = ints[:, 1:].reshape(nelements * 11, 9)
                    nids = ints2[:, 0]

                    obj.element[itotal:itotal2] = eids2
                    obj.element_node[itotal:itotal2, 0] = eids2
                    obj.element_node[itotal:itotal2, 1] = nids

                #[nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq]
                floats2 = floats.reshape(nelements * 11, 9)[:, 1:].copy()
                #sd = floats2[:, 0]
                #obj.data[obj.itime, itotal:itotal2, :] = sd
                obj.data[obj.itime, itotal:itotal2, :] = floats2

                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                n = oef_cbeam_real_100(op2, data, obj,
                                       nelements, ntotal, dt)

        elif result_type == 1 and op2.num_wide == 177: # imag
            slot = op2.get_result(result_name)
            ntotal = 708 * self.factor # 177*4
            nelements = ndata // ntotal

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, ComplexCBeamForceArray)
            if auto_return:
                op2._data_factor = 11
                return nelements * op2.num_wide * 4, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.itotal
                itotal2 = obj.itotal + nelements * 11
                #ielement = obj.ielement
                ielement2 = obj.ielement + nelements

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 177)[:, 1:]
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype8).reshape(nelements, 177)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    assert 0 not in eids, eids
                    eids2 = np.repeat(eids, 11)

                    ints2 = ints[:, 1:].reshape(nelements * 11, 16)
                    nids = ints2[:, 0].copy()

                    obj.element[itotal:itotal2] = eids2
                    obj.element_node[itotal:itotal2, 0] = eids2
                    obj.element_node[itotal:itotal2, 1] = nids

                #[nid, sd, bm1r, bm2r, ts1r, ts2r, afr, ttrqr, wtrqr,
                #          bm1i, bm2i, ts1i, ts2i, afi, ttrqi, wtrqi]
                floats2 = floats.reshape(nelements * 11, 16)[:, 1:].copy()
                isave1 = slice(1, 8)
                isave2 = slice(8, None)
                real_imag = apply_mag_phase(floats2, is_magnitude_phase, isave1, isave2)

                sd = floats2[:, 0]
                obj.data[obj.itime, itotal:itotal2, 0] = sd
                obj.data[obj.itime, itotal:itotal2, 1:] = real_imag

                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if op2.is_debug_file:
                    op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    #op2.binary_debug.write('  #elementi = [eid_device, force]\n')
                    #op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                ntotal = 708 * self.factor # (16*11+1)*4 = 177*4
                nelements = ndata // ntotal
                n = oef_cbeam_imag_177(op2, data, obj, nelements, ntotal, is_magnitude_phase)
        else:
            raise RuntimeError(op2.code_information())
            #msg = op2.code_information()
            #raise RuntimeError(msg)
            #print(msg)
            #return op2._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oef_celas_cdamp(self, data, ndata, dt, is_magnitude_phase,
                         result_type, prefix, postfix):
        """
        11-CELAS1
        12-CELAS2
        13-CELAS3
        14-CELAS4

        20-CDAMP1
        21-CDAMP2
        22-CDAMP3
        23-CDAMP4

        """
        op2 = self.op2
        n = 0
        if op2.element_type == 11:
            result_name = prefix + 'celas1_force' + postfix
            obj_real = RealSpringForceArray
            obj_complex = ComplexSpringForceArray
        elif op2.element_type == 12:
            result_name = prefix + 'celas2_force' + postfix
            obj_real = RealSpringForceArray
            obj_complex = ComplexSpringForceArray
        elif op2.element_type == 13:
            result_name = prefix + 'celas3_force' + postfix
            obj_real = RealSpringForceArray
            obj_complex = ComplexSpringForceArray
        elif op2.element_type == 14:
            result_name = prefix + 'celas4_force' + postfix
            obj_real = RealSpringForceArray
            obj_complex = ComplexSpringForceArray

        elif op2.element_type == 20:
            result_name = prefix + 'cdamp1_force' + postfix
            obj_real = RealDamperForceArray
            obj_complex = ComplexDamperForceArray
        elif op2.element_type == 21:
            result_name = prefix + 'cdamp2_force' + postfix
            obj_real = RealDamperForceArray
            obj_complex = ComplexDamperForceArray
        elif op2.element_type == 22:
            result_name = prefix + 'cdamp3_force' + postfix
            obj_real = RealDamperForceArray
            obj_complex = ComplexDamperForceArray
        elif op2.element_type == 23:
            result_name = prefix + 'cdamp4_force' + postfix
            obj_real = RealDamperForceArray
            obj_complex = ComplexDamperForceArray
        else:
            raise NotImplementedError(op2.code_information())

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)
        if op2.format_code == 1 and op2.num_wide == 2:  # real
            ntotal = 8 * self.factor # 2 * 4
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_real)
            if auto_return:
                return nelements * op2.num_wide * 4, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 2)
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, 2)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                #(eid_device, force)
                obj.data[obj.itime, itotal:itotal2, 0] = floats[:, 1].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                n = oef_celas_cdamp_real_2(op2, data, obj,
                                           nelements, ntotal, dt)

        elif op2.format_code in [2, 3] and op2.num_wide == 3:  # imag
            ntotal = 12 * self.factor  # 3*4
            nelements = ndata // ntotal

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_complex)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.is_debug_file:
                op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                op2.binary_debug.write('  #elementi = [eid_device, force]\n')
                op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 3).copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype8).reshape(nelements, 3)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                #[spring_force]
                real_imag = apply_mag_phase(floats, is_magnitude_phase, 1, 2)
                obj.data[obj.itime, itotal:itotal2, 0] = real_imag
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                n = oef_celas_cdamp_imag_3(op2, data, obj,
                                           nelements, ntotal,
                                           is_magnitude_phase)
        else:
            raise RuntimeError(op2.code_information())
            #msg = 'OEF: element_name=%s element_type=%s' % (op2.element_name, op2.element_type)
            #msg = op2.code_information()
            #print(msg)
            #return op2._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oef_cvisc(self, data, ndata, dt, is_magnitude_phase,
                   result_type, prefix, postfix):
        op2 = self.op2
        result_name = prefix + 'cvisc_force' + postfix
        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        n = 0
        op2._results._found_result(result_name)

        slot = op2.get_result(result_name)
        obj_real = RealViscForceArray

        if op2.format_code == 1 and op2.num_wide == 3: # real
            ntotal = 12 * self.factor # 3 * 4
            nelements = ndata // ntotal

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_real)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 3)
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, 3)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                #(eid_device, axial, torque)
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                n = oef_cvisc_real_3(op2, data, obj,
                                     nelements, ntotal, dt)

        elif op2.format_code in [2, 3] and op2.num_wide == 5: # complex
            ntotal = 20 * self.factor # 5*4
            nelements = ndata // ntotal

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, ComplexViscForceArray)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.is_debug_file:
                op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                op2.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
                op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 5).copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, 5)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                #[axial_force, torque]
                #(eid_device, axial_real, torque_real, axial_imag, torque_imag) = out
                real_imag = apply_mag_phase(floats, is_magnitude_phase, [1, 2], [3, 4])
                obj.data[obj.itime, itotal:itotal2, :] = real_imag
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                n = oef_cvisc_imag_5(op2, data, obj,
                                     nelements, ntotal,
                                     is_magnitude_phase)
        else:
            raise RuntimeError(op2.code_information())
            #msg = op2.code_information()
            #print(msg)
            #return op2._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oef_cbar_34(self, data: bytes, ndata: int, dt: Any,
                     is_magnitude_phase: bool,
                     result_type: str, prefix: str, postfix: str) -> tuple[int, int, int]:
        """
        34-CBAR
        117-CWELDC
        119-CFAST-NX  (126-CFAST-MSC is like the CBUSH)
        118-WELDP-MSC

        """
        op2 = self.op2
        n = 0
        if op2.element_type == 34:
            result_name = prefix + 'cbar_force' + postfix
            obj_real = RealCBarForceArray
            obj_complex = ComplexCBarForceArray
        elif op2.element_type in [117, 200]:
            result_name = prefix + 'cweld_force' + postfix
            obj_real = RealCWeldForceArray
            obj_complex = ComplexCWeldForceArray
            assert op2.num_wide in [9, 17], op2.code_information()
        elif op2.element_type == 118:  # WELDP
            result_name = prefix + 'cweld_force' + postfix
            obj_real = RealCWeldForceArrayMSC
            obj_complex = ComplexCWeldForceArrayMSC
        elif op2.element_type == 119:
            result_name = prefix + 'cfast_force' + postfix
            obj_real = RealCFastForceArrayNX
            assert op2.num_wide == 9, op2.code_information()
        else:
            raise NotImplementedError(op2.element_type)

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)
        #print(result_type in [0, 2], op2.num_wide == 9)
        #print(result_type == 1, op2.num_wide == 17)
        if result_type in [0, 2] and op2.num_wide == 9: # real/random
            # real - format_code == 1
            # random - format_code == 3
            #result_name, is_random = self._apply_oef_ato_crm_psd_rms_no(result_name)
            slot = op2.get_result(result_name)

            ntotal = 36 * self.factor  # 9*4
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_real)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            #return nelements * op2.num_wide * 4
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 9)
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype8).reshape(nelements, 9)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                #[bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            # elif op2.use_vector and is_vectorized and op2.sort_method == 1:
            else:
                n = oef_cbar_real_9(op2, data, obj, nelements, ntotal)
        elif result_type == 1 and op2.num_wide == 17: # imag
            # TODO: vectorize
            ntotal = 68 * self.factor  # 17*4
            nelements = ndata // ntotal
            assert ndata % ntotal == 0

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_complex)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            n = oef_cbar_imag_17(op2, data, obj, nelements, ntotal, is_magnitude_phase)
        else:
            raise RuntimeError(op2.code_information())
            #print(op2.table_name)
            #msg = op2.code_information()
            #print(msg)
            #print(result_type)
            #raise NotImplementedError(op2.code_information())
            #aaa
            #return op2._not_implemented_or_skip(data, ndata, msg), None, None
        #print self.barForces
        return n, nelements, ntotal

    def _oef_cbar_100(self, data, ndata, dt, unused_is_magnitude_phase,
                      result_type, prefix, postfix):
        op2 = self.op2
        n = 0
        #100-BARS
        result_name = prefix + 'cbar_force' + postfix  # _10nodes
        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        if op2.format_code == 1 and op2.num_wide == 8:  # real
            ntotal = 32 * self.factor # 8*4
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, RealCBar100ForceArray)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 8)
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, 8)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                #[axial, torsion, SMa, SMt]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                n = oef_cbar_100_real_8(op2, data, obj,
                                        nelements, ntotal)

        #elif op2.format_code in [2, 3] and op2.num_wide == 14:  # imag
        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
            #msg = op2.code_information()
            #print(msg)
            #return op2._not_implemented_or_skip(data, ndata, msg)
        return n, nelements, ntotal

    def _oef_shells_centroidal(self, data, ndata, dt, is_magnitude_phase,
                               result_type, prefix, postfix):
        """
        33-CQUAD4
        74-CTRIA3
        227-CTRIAR
        228-CQUADR

        """
        op2 = self.op2
        assert op2.element_name != 'RBAR', op2.code_information()
        n = 0
        if op2.element_type == 33:
            result_name = prefix + 'cquad4_force' + postfix
        elif op2.element_type in [74, 83]:
            result_name = prefix + 'ctria3_force' + postfix
        elif op2.element_type == 227:
            result_name = prefix + 'ctriar_force' + postfix
        elif op2.element_type == 228:
            result_name = prefix + 'cquadr_force' + postfix
        else:
            #msg = 'sort1 Type=%s num=%s' % (op2.element_name, op2.element_type)
            #return op2._not_implemented_or_skip(data, ndata, msg)
            raise NotImplementedError(op2.code_information())
        #result_name, is_random = self._apply_oef_ato_crm_psd_rms_no(result_name)

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        assert op2._data_factor == 1, op2._data_factor
        if op2.format_code in [1, 2] and op2.num_wide == 9:
            # real - format_code == 1
            # random - format_code == 2
            ntotal = 36 * self.factor # 9*4
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, RealPlateForceArray)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal
                ielement = obj.ielement
                ielement2 = ielement + nelements

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 9)
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype8).reshape(nelements, 9)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[ielement:ielement2] = eids

                #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
                obj.data[obj.itime, ielement:ielement2, :] = floats[:, 1:].copy()
                obj.itotal = ielement2
                obj.ielement = ielement2
            else:
                n = oef_cquad4_33_real_9(op2, data, obj,
                                         nelements, ntotal)

        elif op2.format_code in [2, 3] and op2.num_wide == 17:  # imag
            ntotal = 68 * self.factor
            nelements = ndata // ntotal

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, ComplexPlateForceArray)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal
                ielement = obj.ielement
                ielement2 = ielement + nelements
                itotal = obj.itotal
                itotal2 = itotal + nelements

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 17).copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype8).reshape(nelements, 17)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
                isave1 = [1, 2, 3, 4, 5, 6, 7, 8]
                isave2 = [9, 10, 11, 12, 13, 14, 15, 16]
                real_imag = apply_mag_phase(floats, is_magnitude_phase, isave1, isave2)
                obj.data[obj.itime, itotal:itotal2, :] = real_imag
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                n = oef_cquad4_33_imag_17(op2, data, ndata, obj,
                                          nelements, ntotal,
                                          is_magnitude_phase)

        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
            #msg = op2.code_information()
            #print(msg)
            #return op2._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oef_shells_nodal(self, data, ndata, dt, is_magnitude_phase,
                          result_type, prefix, postfix):
        """
        64-CQUAD8
        70-CTRIAR
        75-CTRIA6
        82-CQUADR
        144-CQUAD4-bilinear

        """
        op2 = self.op2
        n = 0
        if op2.element_type == 64:
            result_name = prefix + 'cquad8_force' + postfix
        elif op2.element_type == 70:
            result_name = prefix + 'ctriar_force' + postfix
        elif op2.element_type == 75:
            result_name = prefix + 'ctria6_force' + postfix
        elif op2.element_type == 82:
            result_name = prefix + 'cquadr_force' + postfix
        elif op2.element_type == 144:
            result_name = prefix + 'cquad4_force' + postfix
        else:
            raise NotImplementedError(op2.code_information())
            #msg = op2.code_information()
            #return op2._not_implemented_or_skip(data, ndata, msg)

        if op2.element_type in [70, 75]:  # CTRIAR,CTRIA6
            nnodes = 3
        elif op2.element_type in [64, 82, 144]:  # CQUAD8,CQUADR,CQUAD4-bilinear
            nnodes = 4
        else:
            raise NotImplementedError(op2.code_information())
            #msg = 'name=%r type=%r' % (op2.element_name, op2.element_type)
            #return op2._not_implemented_or_skip(data, ndata, msg), None, None

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)

        slot = op2.get_result(result_name)
        nnodes_all = nnodes + 1
        numwide_real = 2 + nnodes_all * 9 # centroidal node is the + 1
        numwide_imag = 2 + nnodes_all * 17

        if op2.format_code == 1 and op2.num_wide == numwide_real:  # real
            obj_real = RealPlateBilinearForceArray

            ntotal = (8 + nnodes_all * 36) * self.factor # centroidal node is the + 1
            assert ntotal == op2.num_wide * 4 * self.factor, 'ntotal=%s numwide=%s' % (ntotal, op2.num_wide * 4)

            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_real)
            if auto_return:
                op2._data_factor = nnodes_all
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                nlayers = nelements * nnodes_all
                n = nelements * op2.num_wide * 4

                istart = obj.itotal
                iend = istart + nlayers
                obj._times[obj.itime] = dt

                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype8).reshape(nelements, numwide_real).copy()
                    # Nastran makes this a 4 for CQUAD4s instead
                    # of 0 like the bilinear stress element...
                    ints[:, 2] = 0

                    nids = ints[:, 2:].reshape(nlayers, 9)[:, 0]
                    eids = ints[:, 0] // 10
                    eids2 = vstack([eids] * nnodes_all).T.ravel()
                    obj.element_node[istart:iend, 0] = eids2
                    obj.element_node[istart:iend, 1] = nids

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, numwide_real)
                results = floats[:, 2:].reshape(nlayers, 9)[:, 1:].copy()
                #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
                obj.data[obj.itime, istart:iend, :] = results
            else:
                n = oef_cquad4_144_real_9(op2, data, obj,
                                          nelements, nnodes)

        elif op2.format_code in [2, 3] and op2.num_wide == numwide_imag: # complex
            ntotal = numwide_imag * 4 * self.factor
            nelements = ndata // ntotal

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, ComplexPlate2ForceArray)
            if auto_return:
                op2._data_factor = nnodes_all
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.itotal
                ielement = obj.ielement
                ielement2 = obj.ielement + nelements
                itotal2 = obj.itotal + nelements * nnodes_all

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, numwide_imag)
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype8).reshape(nelements, numwide_imag).copy()
                    ints[:, 2] = 0
                    ints2 = ints[:, 2:].reshape(nelements * nnodes_all, 17)

                    eids = ints[:, 0] // 10
                    nids = ints2[:, 0]
                    assert eids.min() > 0, eids.min()
                    eids2 = vstack([eids] * nnodes_all).T.ravel()
                    obj.element[ielement:ielement2] = eids
                    obj.element_node[itotal:itotal2, 0] = eids2
                    obj.element_node[itotal:itotal2, 1] = nids

                #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
                floats2 = floats[:, 2:].reshape(nelements * nnodes_all, 17).copy()
                isave1 = [1, 2, 3, 4, 5, 6, 7, 8]
                isave2 = [9, 10, 11, 12, 13, 14, 15, 16]
                real_imag = apply_mag_phase(floats2, is_magnitude_phase, isave1, isave2)
                obj.data[obj.itime, itotal:itotal2, :] = real_imag
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                n = oef_cquad4_imag_17(op2, data, ndata, obj,
                                       nelements, nnodes,
                                       is_magnitude_phase)

        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
            #msg = op2.code_information()
            #print(msg)
            #return op2._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oef_shells_composite(self, data, ndata, dt, unused_is_magnitude_phase,
                              result_type, prefix, postfix):
        """
        95 - CQUAD4
        96 - CQUAD8
        97 - CTRIA3
        98 - CTRIA6 (composite)
        232 - CQUADR
        233 - CTRIAR

        """
        op2 = self.op2
        if op2.element_type == 95:
            result_name = prefix + 'cquad4_composite_force' + postfix
        elif op2.element_type == 96:
            result_name = prefix + 'cquad8_composite_force' + postfix
        elif op2.element_type == 97:
            result_name = prefix + 'ctria3_composite_force' + postfix
        elif op2.element_type == 98:
            result_name = prefix + 'ctria6_composite_force' + postfix
        elif op2.element_type == 232:
            result_name = prefix + 'cquadr_composite_force' + postfix
        elif op2.element_type == 233:
            result_name = prefix + 'ctriar_composite_force' + postfix
        else:  # pragma: no cover
            raise NotImplementedError(op2.code_information())
        if op2._results.is_not_saved(result_name):
            return ndata, None, None

        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        n = 0
        if op2.format_code == 1 and op2.num_wide == 9:  # real
            ntotal = 36 * self.factor # 9 * 4
            nelements = ndata // ntotal


            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, FailureIndicesArray)
            #print('read_mode ', op2.read_mode, auto_return, is_vectorized)
            if auto_return:
                #op2._data_factor = nnodes_all
                return nelements * ntotal, None, None

            obj = op2.obj
            nelements = ndata // ntotal

            ## TODO: add
            #if op2.is_debug_file:
                #op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                ##op2.binary_debug.write('  #centeri = [eid_device, j, grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
                ##op2.binary_debug.write('  #                                fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,)]\n')
                ##op2.binary_debug.write('  #nodeji = [eid, ilayer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)]\n')
                #op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)
            n = oef_shells_composite_real_9(op2, data, obj, nelements, ntotal, dt)
        elif op2.format_code in {2, 3} and op2.num_wide == 9:  # complex
            #  device_code   = 1   Print
            #  analysis_code = 5   Frequency
            #  table_code    = 25  OEFIT-OEF - Composite failure indices
            #  format_code   = 2   Real/Imaginary
            #  result_type   = 1   Complex
            #  sort_method   = 1
            #  sort_code     = 1
            #      sort_bits   = (1, 0, 1)
            #      data_format = 1   Real/Imaginary
            #      sort_type   = 0   Sort1
            #      is_random   = 1   Random Responses
            #  random_code   = 0
            #  element_type  = 95  QUAD4LC-composite
            #  s_code        = None ???
            #  thermal       = 0   isHeatTransfer = False
            #  thermal_bits  = [0, 0, 0, 0, 0]
            #  num_wide      = 9
            #  isubcase      = 1
            #  MSC Nastran
            msg = op2.code_information()
            msg = (f'etype={op2.element_name} ({op2.element_type}) '
                   f'{op2.table_name_str}-COMP-random-numwide={op2.num_wide} '
                   f'numwide_real=11 numwide_imag=9 result_type={result_type}')
            if data is None:
                return op2._not_implemented_or_skip(data, ndata, msg), None, None
            '      ELEMENT-ID =      11'
            '          F A I L U R E   I N D I C E S   F O R   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )'
            '    PLY     FAILURE              FP=FAILURE INDEX FOR PLY    FB=FAILURE INDEX FOR BONDING   FAILURE INDEX FOR ELEMENT      FLAG'
            '     ID      THEORY      FREQ  (DIRECT STRESSES/STRAINS)     (INTER-LAMINAR STRESSES)      MAX OF FP,FB FOR ALL PLIES  '
            '         1   HFABRIC  2.0000E+01        0.1217    1 '
            '                                                                           0.0000                                               '
            '                      4.0000E+01        0.5472   -1 '
            '                                                                           0.0001                                               '
            '                      6.0000E+01        0.2454    1 '
            '                                                                           0.0000                                               '
            '                      8.0000E+01        0.3217    1 '
            '                                                                           0.0000                                               '
            '                      2.0000E+02        0.7284    1 '
            '                                                                           0.0000                        0.7284                 '
            assert op2.sort_method == 1, op2.code_information()
            ntotal = 36 * self.factor
            nelements = ndata // ntotal
            sf = Struct(op2._endian + b'i8s if i ff 4s')  # if if i
            si = Struct(op2._endian + b'i8s if i fi 4s')  # if if f
            sf2 = Struct(op2._endian + b'f')
            for unused_i in range(nelements):
                edata = data[n:n + ntotal]
                out = si.unpack(edata)
                (ply_id, failure_theory_bytes,
                 c, d,
                 e, f,
                 g,
                #ply_fp, failure_index_for_ply,
                #ply_fb, failure_index_for_bonding,
                #failure_index_for_element,
                end) = out
                failure_theory = failure_theory_bytes.decode(op2._encoding).rstrip()
                #print(f'ply_id={ply_id} failure_theory={failure_theory!r} ply_fp={ply_fp} failure_index_for_ply={failure_index_for_ply:.3g} '
                      #f'ply_fb={ply_fb} fi_bonding={failure_index_for_bonding:.3g} fi_for_element={failure_index_for_element:g}')
                #if failure_index_for_bonding != -1:
                    #failure_index_for_bonding = -1000.
                    #failure_index_for_element = f.unpack(edata[-8:-4])

                    #*junk, failure_index_for_bonding, failure_index_for_element, end = sf.unpack(edata)
                    #print(f'ply_id={ply_id} failure_theory={failure_theory!r} ply_fp={ply_fp} fi_ply={failure_index_for_ply:.3g} '
                          #f'ply_fb={ply_fb} fi_bonding={failure_index_for_bonding:.3g} fi_for_element={failure_index_for_element:g}')
                #else:
                if g != -1:
                    g, = sf2.unpack(edata[-8:-4])
                    #op2.show_data(edata)
                    out = (ply_id, failure_theory_bytes, c, d, e, f, g)
                #print(out)
                assert end in {b'    ', b'*** '}, end
                n += ntotal
            #aaa
            return op2._not_implemented_or_skip(data, ndata, msg), None, None
        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
            #msg = op2.code_information()
            #print(msg)
            #return op2._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oef_cshear(self, data, ndata, dt, is_magnitude_phase,
                    result_type, prefix, postfix):
        """4-CSHEAR"""
        op2 = self.op2
        result_name = prefix + 'cshear_force' + postfix
        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        n = 0
        if op2.format_code == 1 and op2.num_wide == 17:  # real
            ntotal = 68  # 17*4
            nelements = ndata // ntotal

            obj_real = RealCShearForceArray
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_real)
            if auto_return:
                return nelements * op2.num_wide * 4, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * 4 * op2.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 17)
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, 17)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                # [f41, f21, f12, f32, f23, f43, f34, f14, kf1,
                #  s12, kf2, s23, kf3, s34, kf4, s41]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                n = oef_cshear_real_17(op2, data, obj,
                                       nelements, ntotal, dt)

        elif op2.format_code in [2, 3] and op2.num_wide == 33:  # imag
            ntotal = 132 * self.factor # 33*4
            nelements = ndata // ntotal

            obj_complex = ComplexCShearForceArray
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_complex)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            #if op2.is_debug_file:
                #op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                #op2.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
                #op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * 4 * op2.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 33).copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype8).reshape(nelements, 33)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                #[f41r, f21r, f12r, f32r, f23r, f43r, f34r, f14r
                # kf1r, s12r, kf2r, s23r, kf3r, s34r, kf4r, s41r
                # f41i, f21i, f12i, f32i, f23i, f43i, f34i, f14i
                # kf1i, s12i, kf2i, s23i, kf3i, s34i, kf4i, s41i]
                isave1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
                isave2 = [17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32]
                real_imag = apply_mag_phase(floats, is_magnitude_phase, isave1, isave2)
                obj.data[obj.itime, itotal:itotal2, :] = real_imag

                #obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:]
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                n = oef_cshear_imag_33(op2, data, obj,
                                       nelements, ntotal,
                                       is_magnitude_phase)
        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
            #msg = op2.code_information()
            #return op2._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oef_cconeax(self, data, ndata, dt, unused_is_magnitude_phase,
                     result_type, prefix, postfix):
        """35-CONEAX"""
        op2 = self.op2
        result_name = prefix + 'cconeax_force' + postfix
        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        n = 0
        if op2.format_code == 1 and op2.num_wide == 7:  # real
            ntotal = 28  # 7*4
            nelements = ndata // ntotal

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, RealConeAxForceArray)
            if auto_return:
                return nelements * op2.num_wide * 4, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * 4 * op2.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 7)
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, 7)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                # [hopa, bmu, bmv, tm, su, sv]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                n = oef_cconeax_real_7(op2, data, obj,
                                       nelements, ntotal, dt)

        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
            #msg = op2.code_information()
            #return op2._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oef_cgap(self, data, ndata, dt, unused_is_magnitude_phase,
                  result_type, prefix, postfix):
        """38-GAP"""
        op2 = self.op2
        result_name = prefix + 'cgap_force' + postfix
        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        n = 0
        if op2.format_code == 1 and op2.num_wide == 9:  # real
            ntotal = 36 *  self.factor # 9*4
            nelements = ndata // ntotal
            obj = op2.obj

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, RealCGapForceArray)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 9)
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, 9)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                # [fx, sfy, sfz, u, v, w, sv, sw]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                n = oef_cgap_real_9(op2, data, obj,
                                    nelements, ntotal)

        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
            #msg = op2.code_information()
            #print(msg)
            #return op2._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oef_cbend(self, data, ndata, dt, is_magnitude_phase,
                   result_type, prefix, postfix):
        """69-CBEND"""
        op2 = self.op2
        result_name = prefix + 'cbend_force' + postfix
        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        n = 0
        if op2.format_code == 1 and op2.num_wide == 15:  # real
            ntotal = 60  # 15*4
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, RealBendForceArray)
            if auto_return:
                return nelements * op2.num_wide * 4, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * 4 * op2.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 15).copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, 15).copy()
                    eids = ints[:, 0] // 10
                    nids_a = ints[:, 1]
                    nids_b = ints[:, 8]
                    assert eids.min() > 0, eids.min()
                    obj.element_node[itotal:itotal2, 0] = eids
                    obj.element_node[itotal:itotal2, 1] = nids_a
                    obj.element_node[itotal:itotal2, 2] = nids_b

                # [nid_a, bm1_a, bm2_a, ts1_a, ts2_a, af_a, trq_a,
                #  nid_b, bm1_b, bm2_b, ts1_b, ts2_b, af_b, trq_b]
                assert floats[:, 2:8].shape[1] == 6, floats[:, 2:8].shape
                assert floats[:, 9:].shape[1] == 6, floats[:, 9:].shape
                obj.data[obj.itime, itotal:itotal2, :6] = floats[:, 2:8]
                obj.data[obj.itime, itotal:itotal2, 6:] = floats[:, 9:]
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                s = Struct(op2._endian + op2._analysis_code_fmt + b' i6fi6f')
                add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]

                    out = s.unpack(edata)
                    if op2.is_debug_file:
                        op2.binary_debug.write('OEF_BEND-69 - %s\n' % (str(out)))
                    (eid_device,
                     nid_a, bm1_a, bm2_a, ts1_a, ts2_a, af_a, trq_a,
                     nid_b, bm1_b, bm2_b, ts1_b, ts2_b, af_b, trq_b) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, op2.nonlinear_factor, op2.sort_method)

                    add_sort_x(
                        dt, eid,
                        nid_a, bm1_a, bm2_a, ts1_a, ts2_a, af_a, trq_a,
                        nid_b, bm1_b, bm2_b, ts1_b, ts2_b, af_b, trq_b)
                    n += ntotal
        elif op2.format_code in [2, 3] and op2.num_wide == 27:  # imag
            # TODO: vectorize
            ntotal = 108  # 27*4
            nelements = ndata // ntotal

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, ComplexCBendForceArray)
            if auto_return:
                return nelements * op2.num_wide * 4, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * 4 * op2.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                ireal = [2, 3, 4, 5, 6, 7, 15, 16, 17, 18, 19, 20]
                iimag = [8, 9, 10, 11, 12, 13, 21, 22, 23, 24, 25, 26]
                # 0    1
                # eid, nidA,
                # 2      3      4       5     6     7
                # 8      9      10      11    12    13
                # bm1Ar, bm2Ar, ts1Ar, ts2Ar, afAr, trqAr,
                # bm1Ai, bm2Ai, ts1Ai, ts2Ai, afAi, trqAi,
                # 14
                # nidB
                # 15     16     17     18     19    20
                # 21     22     23     24     25    26
                # bm1Br, bm2Br, ts1Br, ts2Br, afBr, trqBr,
                # bm1Bi, bm2Bi, ts1Bi, ts2Bi, afBi, trqBi
                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 27).copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, 27).copy()
                    eids = ints[:, 0] // 10
                    nids_a = ints[:, 1]
                    nids_b = ints[:, 14]
                    assert nids_a.min() > 0, nids_a
                    assert nids_b.min() > 0, nids_b
                    assert eids.min() > 0, eids.min()
                    #print(nids_b)
                    obj.element_node[itotal:itotal2, 0] = eids
                    obj.element_node[itotal:itotal2, 1] = nids_a
                    obj.element_node[itotal:itotal2, 2] = nids_b

                real_imag = apply_mag_phase(floats, is_magnitude_phase, ireal, iimag)
                obj.data[obj.itime, itotal:itotal2, :] = real_imag
                #assert floats[:, 1:6].shape[1] == 12, floats[:, 1:6].shape
                #assert floats[:, 7:].shape[1] == 6, floats[:, 7:].shape
                #obj.data[obj.itime, itotal:itotal2, :6] = floats[:, 1:6]
                #obj.data[obj.itime, itotal:itotal2, 6:] = floats[:, 7:]
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                s = Struct(op2._endian + op2._analysis_code_fmt + b' i12f i12f')
                add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
                for unused_i in range(nelements):
                    edata = data[n:n+108]
                    n += ntotal
                    out = s.unpack(edata)
                    if op2.is_debug_file:
                        op2.binary_debug.write('OEF_BEND-69 - %s\n' % (str(out)))
                    (eid_device, nid_a,
                     bm1_ar, bm2_ar, ts1_ar, ts2_ar, af_ar, trq_ar,
                     bm1_ai, bm2_ai, ts1_ai, ts2_ai, af_ai, trq_ai,
                     nid_b,
                     bm1_br, bm2_br, ts1_br, ts2_br, af_br, trq_br,
                     bm1_bi, bm2_bi, ts1_bi, ts2_bi, af_bi, trq_bi) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, op2.nonlinear_factor, op2.sort_method)

                    if is_magnitude_phase:
                        bm1_a = polar_to_real_imag(bm1_ar, bm1_ai)
                        bm1_b = polar_to_real_imag(bm1_br, bm1_bi)
                        bm2_a = polar_to_real_imag(bm2_ar, bm2_ai)
                        bm2_b = polar_to_real_imag(bm2_br, bm2_bi)
                        ts1_a = polar_to_real_imag(ts1_ar, ts1_ai)
                        ts1_b = polar_to_real_imag(ts1_br, ts1_bi)
                        ts2_a = polar_to_real_imag(ts2_ar, ts2_ai)
                        ts2_b = polar_to_real_imag(ts2_br, ts2_bi)
                        af_a = polar_to_real_imag(af_ar, af_ai)
                        af_b = polar_to_real_imag(af_br, af_bi)
                        trq_a = polar_to_real_imag(trq_ar, trq_ai)
                        trq_b = polar_to_real_imag(trq_br, trq_bi)
                    else:
                        bm1_a = complex(bm1_ar, bm1_ai)
                        bm1_b = complex(bm1_br, bm1_bi)
                        bm2_a = complex(bm2_ar, bm2_ai)
                        bm2_b = complex(bm2_br, bm2_bi)
                        ts1_a = complex(ts1_ar, ts1_ai)
                        ts1_b = complex(ts1_br, ts1_bi)
                        ts2_a = complex(ts2_ar, ts2_ai)
                        ts2_b = complex(ts2_br, ts2_bi)
                        af_a = complex(af_ar, af_ai)
                        af_b = complex(af_br, af_bi)
                        trq_a = complex(trq_ar, trq_ai)
                        trq_b = complex(trq_br, trq_bi)

                    add_sort_x(dt, eid,
                               nid_a, bm1_a, bm2_a, ts1_a, ts2_a, af_a, trq_a,
                               nid_b, bm1_b, bm2_b, ts1_b, ts2_b, af_b, trq_b)
        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
            #msg = op2.code_information()
            #return op2._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oef_csolid_pressure(self, data, ndata, dt, is_magnitude_phase,
                             result_type, prefix, postfix):
        """
        76-HEXPR
        77-PENPR
        78-TETPR
        """
        op2 = self.op2
        n = 0
        if op2.element_type == 76:
            result_name = prefix + 'chexa_pressure_force' + postfix
        elif op2.element_type == 77:
            result_name = prefix + 'cpenta_pressure_force' + postfix
        elif op2.element_type == 78:
            result_name = prefix + 'ctetra_pressure_force' + postfix
        elif op2.element_type == 79:
            result_name = prefix + 'cpyram_pressure_force' + postfix
        else:
            msg = op2.code_information()
            return op2._not_implemented_or_skip(data, ndata, msg), None, None
        slot = op2.get_result(result_name)

        op2._results._found_result(result_name)
        if op2.format_code == 1 and op2.num_wide == 10:  # real
            ntotal = 40 * self.factor
            nelements = ndata // ntotal
            #nelements = ndata // ntotal

            obj_real = RealSolidPressureForceArray
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_real)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            #if op2.is_debug_file:
                #op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                #op2.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
                #op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                # self.itime = 0
                # self.ielement = 0
                # self.itotal = 0
                #self.ntimes = 0
                #self.nelements = 0
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 10)
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, 10)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                #[axial_force, torque]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 3:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                n = oef_csolid_pressure_10(op2, data, obj,
                                           nelements, ntotal, dt)

        elif op2.format_code in [2, 3] and op2.num_wide == 16:  # imag
            ntotal = 64 * self.factor
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, ComplexSolidPressureForceArray)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            #if op2.is_debug_file:
                #op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                #op2.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
                #op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 16).copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype8).reshape(nelements, 16)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                #[xaccr, yaccr, zaccr, xvelr, yvelr, zvelr, pressure,
                # xacci, yacci, zacci, xveli, yveli, zveli]
                if is_magnitude_phase:
                    mag = floats[:, [3, 4, 5, 6, 7, 8, 9]]
                    phase = hstack([
                        floats[:, [10, 11, 12, 13, 14, 15]],
                        zeros((len(floats), 1), dtype='float32')
                    ])
                    rtheta = radians(phase)
                    real_imag = mag * (cos(rtheta) + 1.j * sin(rtheta))
                else:
                    real = floats[:, [3, 4, 5, 6, 7, 8, 9]]
                    imag = hstack([
                        floats[:, [10, 11, 12, 13, 14, 15]],
                        zeros((len(floats), 1), dtype='float32')
                    ])
                    real_imag = real + 1.j * imag
                obj.data[obj.itime, itotal:itotal2, :] = real_imag
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                n = oef_csolid_imag_16(op2, data, obj,
                                       nelements, ntotal,
                                       is_magnitude_phase)
        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
            #msg = op2.code_information()
            #return op2._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oef_cbush(self, data, ndata, dt, is_magnitude_phase,
                   result_type, prefix, postfix):
        """
        102-CBUSH
        126-CFAST-MSC
        280-CBEAR

        """
        op2 = self.op2
        num_wide = op2.num_wide
        if op2.element_type == 102:
            result_name = prefix + 'cbush_force' + postfix
            real_obj = RealCBushForceArray
            complex_obj = ComplexCBushForceArray
            assert num_wide in [7, 13], op2.code_information()
        elif op2.element_type == 126:
            result_name = prefix + 'cfast_force' + postfix
            real_obj = RealCFastForceArrayMSC
            complex_obj = ComplexCFastForceArrayMSC
            assert num_wide in [7, 13], op2.code_information()
        elif op2.element_type == 280:
            result_name = prefix + 'cbear_force' + postfix
            assert num_wide in [7, 13], op2.code_information()
            real_obj = RealCBearForceArray
            complex_obj = ComplexCBearForceArray
        else:
            raise NotImplementedError(op2.code_information())

        #print(op2.code_information())
        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        #result_name, is_random = self._apply_oef_ato_crm_psd_rms_no(result_name)
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        n = 0
        #op2.log.warning('dt=%s num_wide=%s result_type=%s', dt, num_wide, result_type)
        if result_type in {0, 2} and num_wide == 7:  # real/random
            numwide_real = 7
            # real - format_code == 1
            # random - format_code == 3
            ntotal = 28 *  self.factor # 7*4
            nelements = ndata // ntotal

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, real_obj)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                # self.itime = 0
                # self.ielement = 0
                # self.itotal = 0
                #self.ntimes = 0
                #self.nelements = 0
                n = nelements * ntotal

                istart = obj.itotal
                iend = istart + nelements
                obj._times[obj.itime] = dt

                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype8).reshape(nelements, numwide_real).copy()
                    eids = ints[:, 0] // 10
                    obj.element[istart:iend] = eids
                results = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, numwide_real)

                #[fx, fy, fz, mx, my, mz]
                obj.data[obj.itime, istart:iend, :] = results[:, 1:].copy()
            else:
                n = oef_cbush_real_7(op2, data, obj,
                                     nelements, ntotal, dt)
        elif result_type == 1 and op2.num_wide == 13:  # imag
            # TCODE,7 =1 Real/imaginary or magnitude/phase
            # 2 FXR RS Force x - real/mag. part
            # 3 FYR RS Force y - real/mag. part
            # 4 FZR RS Force z - real/mag. part
            # 5 MXR RS Moment x - real/mag. part
            # 6 MYR RS Moment y - real/mag. part
            # 7 MZR RS Moment z - real/mag. part
            # 8 FXI RS Force x - imag./phase part
            # 9 FYI RS Force y - imag./phase part
            # 10 FZI RS Force z - imag./phase part
            # 11 MXI RS Moment x - imag./phase part
            # 12 MYI RS Moment y - imag./phase part
            # 13 MZI RS Moment z - imag./phase part

            # TODO: vectorize
            ntotal = 52 * self.factor  # 13*4
            nelements = ndata // ntotal
            #result_name = prefix + 'cbush_force' + postfix
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, complex_obj)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            n = oef_cbush_imag_13(op2, data, obj,
                                  nelements, ntotal,
                                  is_magnitude_phase)
        #elif op2.format_code == 2 and op2.num_wide == 7:
            #op2.log.warning(op2.code_information())
        else:
            raise NotImplementedError(op2.code_information())
        #else:  # pragma: no cover
            #msg = op2.code_information()
            #print(msg)
            #return op2._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

def oef_crod_real_3(op2: OP2, data: bytes,
                    obj: RealRodForceArray,
                    nelements: int, ntotal: int) -> int:
    n = 0
    fmt = mapfmt(op2._endian + op2._analysis_code_fmt + b'ff', op2.size)  # 3
    s = Struct(fmt)
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = s.unpack(edata)
        (eid_device, axial, torque) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_Rod - %s\n' % (str(out)))
        add_sort_x(dt, eid, axial, torque)
        n += ntotal
    return n

def oef_crod_imag_5(op2: OP2, data: bytes,
                    obj: ComplexRodForceArray,
                    nelements: int, ntotal: int,
                    is_magnitude_phase: bool) -> int:
    n = 0
    fmt = mapfmt(op2._endian + op2._analysis_code_fmt + b'4f', op2.size)
    s = Struct(fmt)
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]

        out = s.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_CRod - %s\n' % (str(out)))
        (eid_device, axial_real, torque_real, axial_imag, torque_imag) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if is_magnitude_phase:
            axial = polar_to_real_imag(axial_real, axial_imag)
            torque = polar_to_real_imag(torque_real, torque_imag)
        else:
            axial = complex(axial_real, axial_imag)
            torque = complex(torque_real, torque_imag)

        add_sort_x(dt, eid, axial, torque)
        n += ntotal
    return n

def oef_celas_cdamp_imag_3(self, data: bytes,
                           obj: Union[ComplexSpringForceArray, ComplexDamperForceArray],
                           nelements: int, ntotal: int,
                           is_magnitude_phase: bool) -> int:
    op2 = self
    n = 0
    fmt = mapfmt(op2._endian + op2._analysis_code_fmt + b'2f', self.size)
    structi = Struct(fmt)
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = structi.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_SpringDamper - %s\n' % str(out))
        (eid_device, force_real, force_imag) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if is_magnitude_phase:
            force = polar_to_real_imag(force_real, force_imag)
        else:
            force = complex(force_real, force_imag)
        add_sort_x(dt, eid, force)
        n += ntotal
    return n

def oef_cgap_real_9(self, data: bytes,
                    obj: RealCGapForceArray,
                    nelements: int, ntotal: int) -> int:
    op2 = self
    n = 0
    s = Struct(op2._endian + op2._analysis_code_fmt + b'8f')
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+36]

        out = s.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_CGAP-38 - %s\n' % (str(out)))
        (eid_device, fx, sfy, sfz, u, v, w, sv, sw) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #data_in = [eid, fx, sfy, sfz, u, v, w, sv, sw]
        #print "%s" %(self.get_element_type(op2.element_type)), data_in
        #eid = obj.add_new_eid_sort1(out)
        add_sort_x(dt, eid, fx, sfy, sfz, u, v, w, sv, sw)
        n += ntotal
    return n

def oef_cbar_real_9(self, data: bytes,
                    obj: RealCBarForceArray,
                    nelements: int, ntotal: int) -> int:
    op2 = self
    n = 0
    fmt = mapfmt(op2._endian + op2._analysis_code_fmt + b'8f', self.size)
    s = Struct(fmt)
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]

        out = s.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_CBar - %s\n' % (str(out)))
        (eid_device, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #data_in = [eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq]
        add_sort_x(dt, eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq)
        n += ntotal
    return n

def oef_cbar_imag_17(self, data: bytes,
                     obj: ComplexCBarForceArray,
                     nelements: int, ntotal: int,
                     is_magnitude_phase: bool) -> int:
    op2 = self
    n = 0
    fmt = mapfmt(op2._endian + op2._analysis_code_fmt + b'16f', self.size)
    s = Struct(fmt)
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]

        out = s.unpack(edata)
        (eid_device,
         bm1ar, bm2ar, bm1br, bm2br, ts1r, ts2r, afr, trqr,
         bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi) = out
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_CBar - %s\n' % (str(out)))
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        if is_magnitude_phase:
            bm1a = polar_to_real_imag(bm1ar, bm1ai)
            bm2a = polar_to_real_imag(bm2ar, bm2ai)
            bm1b = polar_to_real_imag(bm1br, bm1bi)
            bm2b = polar_to_real_imag(bm2br, bm2bi)
            ts1 = polar_to_real_imag(ts1r, ts1i)
            ts2 = polar_to_real_imag(ts2r, ts2i)
            af = polar_to_real_imag(afr, afi)
            trq = polar_to_real_imag(trqr, trqi)
        else:
            bm1a = complex(bm1ar, bm1ai)
            bm2a = complex(bm2ar, bm2ai)
            bm1b = complex(bm1br, bm1bi)
            bm2b = complex(bm2br, bm2bi)
            ts1 = complex(ts1r, ts1i)
            ts2 = complex(ts2r, ts2i)
            af = complex(afr, afi)
            trq = complex(trqr, trqi)

        #data_in = [bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq]
        #print("eid_device=%s eid=%s dt=%s %s" % (
            #eid_device, eid, dt, self.get_element_type(op2.element_type)), data_in)
        add_sort_x(dt, eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq)
        n += ntotal
    return n

def oef_cbar_100_real_8(self, data: bytes,
                        obj: RealCBar100ForceArray,
                        nelements: int, ntotal: int) -> int:
    op2 = self
    n = 0
    s = Struct(op2._endian + op2._analysis_code_fmt + b'7f')
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]

        out = s.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_CBar100 - %s\n' % (str(out)))
        (eid_device, sd, bm1, bm2, ts1, ts2, af, trq) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        add_sort_x(dt, eid, sd, bm1, bm2, ts1, ts2, af, trq)
        n += ntotal
    return n

def oef_cbeam_real_100(self, data: bytes, obj: RealCBeamForceArray,
                       nelements: int, ntotal: int, dt: Any) -> int:
    op2 = self
    n = 0
    ntotal1 = 4 * self.factor
    ntotal2 = 36 * self.factor
    if self.size == 4:
        s1 = self.struct_i
        s2 = Struct(op2._endian + b'i8f')  # 36
    else:
        s1 = self.struct_q
        s2 = Struct(op2._endian + b'q8d')  # 36

    sort_method = op2.sort_method
    add_sort_x = getattr(obj, 'add_sort' + str(sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal1]
        eid_device, = s1.unpack(edata)
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, sort_method)
        n += ntotal1

        for istation in range(11):
            edata = data[n:n+ntotal2]
            out = s2.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('OEF_Beam - %s\n' % (str(out)))
            (nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq) = out

            if istation == 0 or sd > 0:
                add_sort_x(dt, eid, nid, sd, bm1, bm2, ts1, ts2,
                           af, ttrq, wtrq)
            n += ntotal2
    return n

def oef_cbeam_imag_177(self, data: bytes,
                       obj: ComplexCBeamForceArray,
                       nelements: int, ntotal: int,
                       is_magnitude_phase: bool) -> int:
    op2 = self
    n = 0
    #s1 = self.struct_i
    ntotal1 = 4 * self.factor
    ntotal2 = 64 * self.factor
    s1 = Struct(mapfmt(op2._endian + b'i', op2.size)) # self.struct_i
    s2 = Struct(mapfmt(op2._endian + b'i15f', op2.size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + ntotal1]
        eid_device, = s1.unpack(edata)
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        n += ntotal1
        for unused_istation in range(11):
            edata = data[n:n+ntotal2]
            n += ntotal2
            out = s2.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('OEF_Beam - %s\n' % (str(out)))
            (nid, sd,
             bm1r, bm2r, ts1r, ts2r, afr, ttrqr, wtrqr,
             bm1i, bm2i, ts1i, ts2i, afi, ttrqi, wtrqi) = out

            if is_magnitude_phase:
                bm1 = polar_to_real_imag(bm1r, bm1i)
                bm2 = polar_to_real_imag(bm2r, bm2i)
                ts1 = polar_to_real_imag(ts1r, ts1i)
                ts2 = polar_to_real_imag(ts2r, ts2i)
                af = polar_to_real_imag(afr, afi)
                ttrq = polar_to_real_imag(ttrqr, ttrqi)
                wtrq = polar_to_real_imag(wtrqr, wtrqi)
            else:
                bm1 = complex(bm1r, bm1i)
                bm2 = complex(bm2r, bm2i)
                ts1 = complex(ts1r, ts1i)
                ts2 = complex(ts2r, ts2i)
                af = complex(afr, afi)
                ttrq = complex(ttrqr, ttrqi)
                wtrq = complex(wtrqr, wtrqi)

            #if i == 0:
                #obj.add_new_element_sort1(
                    #dt, eid, nid, sd, bm1, bm2, ts1, ts2,
                    #af, ttrq, wtrq)
            #elif sd > 0.:
            add_sort_x(
                dt, eid, nid, sd, bm1, bm2, ts1, ts2,
                af, ttrq, wtrq)
            #else:
                ## don't add this field
                #pass
                #raise RuntimeError('CBEAM error; i=%s sd=%s' % (i, sd))
    return n

def oef_shells_composite_real_9(self, data: bytes,
                                obj: FailureIndicesArray,
                                nelements: int, ntotal: int,
                                dt: Any) -> int:
    """
    2 THEORY(2) CHA/R4
    4 PLY       I
    5 DIRECT    RS
    6 INDEX     CHA/R4
    7 LAMIN     RS (I)
    8 MAX       RS (I)
    9 FLAG      CHA
    """
    op2 = self
    n = 0
    size = self.size
    if size == 4:
        #                                5 6  7 8-i/f 9
        s1 = Struct(op2._endian + b'i8sif  i f i     4s')
        s2 = Struct(op2._endian + b'i8sif 4s f f     4s')
    elif size == 8:
        s1 = Struct(op2._endian + b'q16sqd  q d q     8s')
        s2 = Struct(op2._endian + b'q16sqd 8s d d     8s')
    else:  # pragma: no cover
        raise RuntimeError(size)

    eid_old = None
    #print(op2.element_type)
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        #2 THEORY(2) CHAR4 Theory
        #4 LAMID     I Lamina number

        #5 FP       RS Failure index for direct stresses
        #6 FM       RS Failure mode for maximum strain theory
        #7 FB       RS Failure index for interlaminar shear stress or -1
        #8 FMAX     RS Maximum of FP and FB or -1.
        #9 FFLAG CHAR4 Failure flag
        edata = data[n:n+ntotal]  # 4*9
        #print(self.show_data(edata[4+8+4:-4-4], types='ifs'))
        out1 = s1.unpack(edata)
        out2 = s2.unpack(edata)

        # failure_stress_for_ply = failure_strain_for_ply = failure_index_for_ply???
        # i    8s               i      f
        (eid, failure_theoryb, ply_id, failure_stress_for_ply,
         # 4s   7-f/i                8-f/i      9-4s
         flagi, interlaminar_stress, max_value, failure_flagb,
         #failure_index_for_bonding,
         #failure_index_for_element,
         #flag,
         #direct_stress_or_strain,
         #interlaminar_stress,
         #max_of_fb_fp_for_all_plies
        ) = out1

        (_eid, _failure_theoryb, _ply_id, _failure_stress_for_ply,
         # 4s   7-f/i                8-f/i      9-4s
         flagb, _interlaminar_stress, max_value_float, _failure_flagb,
         #failure_index_for_bonding,
         #failure_index_for_element,
         #flag,
         #direct_stress_or_strain,
         #interlaminar_stress,
         #max_of_fb_fp_for_all_plies
        ) = out2
        #print(interlaminar_stress, _interlaminar_stress)

        #print('failure_flagb = %r' % failure_flagb)
        if flagi == 0:
            flag = ''
        else:
            #print('flagb = %r' % flagb)
            flag = reshape_bytes_block_strip(flagb, size=size)

        failure_theory = reshape_bytes_block_strip(failure_theoryb, size=size)
        failure_flag = reshape_bytes_block_strip(failure_flagb, size=size)
        #print('flag = %r' % flag)
        #print('failure_flag = %r' % failure_flag)

        if max_value == -1:
            max_value = np.nan
        else:
            max_value = max_value_float # out2[6]

        if eid == -1:
            #print(f'  ply_id={ply_id} failure_stress_for_ply={failure_stress_for_ply:g} '
                  #f'flag={flag!r} interlaminar_stress={interlaminar_stress} '
                  #f'max_value={max_value:g} failure_flag={failure_flag!r}')
            eid = eid_old
        else:
            #print(f"eid={eid} ft={failure_theory!r}\n"
                  #f'  ply_id={ply_id} failure_stress_for_ply={failure_stress_for_ply:g} '
                  #f'flag={flag!r} interlaminar_stress={interlaminar_stress} '
                  #f'max_value={max_value} failure_flag={failure_flag!r}')
            eid_old = eid
        assert flag in ['', '-1', '-2', '-12', 'IN'], f'flag={flag!r} flagb={flagb!r}'

        # 'HILL' for the Hill theory.
        # 'HOFF' for the Hoffman theory.
        # 'TSAI' for the Tsai-Wu theory.
        # 'STRN' for the Maximum Strain theory.
        # 'HFAIL' for the Hashin failure criterion
        # 'HTAPE' for the Hashin tape criterion
        # 'HFABR' for the Hashin fabric criterion
        assert failure_theory in ['TSAI-WU', 'STRAIN', 'HILL', 'HOFFMAN', 'HFAIL', 'HFABRIC', 'HTAPE', ''], f'failure_theory={failure_theory!r}'
        assert failure_flag in ['', '***'], 'failure_flag=%r' % failure_flag
        add_sort_x(dt, eid, failure_theory, ply_id, failure_stress_for_ply, flag,
                   interlaminar_stress, max_value, failure_flag)
        n += ntotal

    #add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    #s = Struct(op2._endian + b'i8si4f4s')
    #for i in range(nelements):
        #if i % 10000 == 0:
            #print 'i = ', i
        #edata = data[n:n+ntotal]  # 4*9
        #out = s.unpack(edata)
        #(eid_device, theory, lamid, failure_index_direct_stress, failure_mode_max_shear,
                 #failure_index_interlaminar_shear, fmax, failure_flag) = out
        #eid, dt = get_eid_dt_from_eid_device(
            #eid_device, op2.nonlinear_factor, op2.sort_method)
        #if op2.is_debug_file:
            #if eid > 0:
                #op2.binary_debug.write('  eid=%i; C=[%s]\n' % (', '.join(['%r' % di for di in out]) ))
            #else:
                #op2.binary_debug.write('      %s  C=[%s]\n' % (' ' * len(str(eid)), ', '.join(['%r' % di for di in out]) ))

        #if eid > 0:
            #obj.add_new_eid_sort1(eType, dt, eid, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)
        #else:
            #add_sort_x(dt, eid, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)
        #n += ntotal
    return n

def oef_cbush_real_7(self, data: bytes,
                     obj: Union[RealCBushForceArray, RealCFastForceArrayMSC, RealCBearForceArray],
                     nelements: int, ntotal: int, dt: Any) -> int:
    op2 = self
    n = 0
    struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'6f', op2.size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_CBUSH-102 - %s\n' % (str(out)))
        (eid_device, fx, fy, fz, mx, my, mz) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #op2.log.debug('eid=%s dt=%s', eid, dt)
        add_sort_x(dt, eid, fx, fy, fz, mx, my, mz)
        n += ntotal
    return n

def oef_cbush_imag_13(self, data: bytes,
                      obj: Union[ComplexCBushForceArray, ComplexCFastForceArrayMSC],
                      nelements: int, ntotal: int,
                      is_magnitude_phase: bool) -> int:
    """
    102-CBUSH
    126-CFAST-MSC
    280-CBEAR

    """
    op2 = self
    n = 0
    s = Struct(op2._endian + op2._analysis_code_fmt + b'12f')
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = s.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_CBUSH-102 - %s\n' % (str(out)))
        (eid_device,
         fxr, fyr, fzr, mxr, myr, mzr,
         fxi, fyi, fzi, mxi, myi, mzi) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        if is_magnitude_phase:
            fx = polar_to_real_imag(fxr, fxi)
            mx = polar_to_real_imag(mxr, mxi)
            fy = polar_to_real_imag(fyr, fyi)
            my = polar_to_real_imag(myr, myi)
            fz = polar_to_real_imag(fzr, fzi)
            mz = polar_to_real_imag(mzr, mzi)
        else:
            fx = complex(fxr, fxi)
            mx = complex(mxr, mxi)
            fy = complex(fyr, fyi)
            my = complex(myr, myi)
            fz = complex(fzr, fzi)
            mz = complex(mzr, mzi)

        add_sort_x(dt, eid, fx, fy, fz, mx, my, mz)
        n += ntotal
    return n

def oef_cvisc_real_3(self, data: bytes, obj: RealViscForceArray,
                     nelements: int, ntotal: int, dt: Any) -> int:
    op2 = self
    n = 0
    fmt = mapfmt(op2._endian + op2._analysis_code_fmt + b'ff', self.size)
    s = Struct(fmt)
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]

        out = s.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_CVisc - %s\n' % (str(out)))
        (eid_device, axial, torque) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        add_sort_x(dt, eid, axial, torque)
        n += ntotal
    return n

def oef_cvisc_imag_5(self, data: bytes,
                     obj: ComplexViscForceArray,
                     nelements: int, ntotal: int,
                     is_magnitude_phase: bool) -> int:
    op2 = self
    n = 0
    s = Struct(op2._endian + op2._analysis_code_fmt + b'4f')  # 5
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+20]

        out = s.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_CVisc - %s\n' % (str(out)))
        (eid_device, axial_real, torque_real, axial_imag, torque_imag) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if is_magnitude_phase:
            axial = polar_to_real_imag(axial_real, axial_imag)
            torque = polar_to_real_imag(torque_real, torque_imag)
        else:
            axial = complex(axial_real, axial_imag)
            torque = complex(torque_real, torque_imag)

        add_sort_x(dt, eid, axial, torque)
        n += ntotal
    return n

def oef_celas_cdamp_real_2(self, data: bytes,
                           obj: Union[RealSpringForceArray, RealDamperForceArray],
                           nelements: int, ntotal: int, dt: Any) -> int:
    op2 = self
    n = 0
    fmt = mapfmt(op2._endian + op2._analysis_code_fmt + b'f', self.size)
    s = Struct(fmt)  # 2
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = s.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_SpringDamper - %s\n' % str(out))
        (eid_device, force) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        add_sort_x(dt, eid, force)
        n += ntotal
    return n

def oef_cshear_real_17(self, data: bytes,
                       obj: RealCShearForceArray,
                       nelements: int, ntotal: int, dt: Any) -> int:
    op2 = self
    n = 0
    s = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'16f', op2.size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]

        out = s.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_Shear - %s\n' % (str(out)))
        (eid_device,
         f41, f21, f12, f32, f23, f43, f34, f14, kf1,
         s12, kf2, s23, kf3, s34, kf4, s41) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        #data_in = [eid,
                   #f41, f21, f12, f32, f23, f43, f34,
                   #f14, kf1, s12, kf2, s23, kf3, s34, kf4, s41]
        #print "%s" % (self.get_element_type(op2.element_type)), data_in
        add_sort_x(dt, eid,
                   f41, f21, f12, f32, f23, f43, f34,
                   f14, kf1, s12, kf2, s23, kf3, s34, kf4, s41)
        n += ntotal
    return n

def oef_cshear_imag_33(self, data: bytes,
                       obj: ComplexCShearForceArray,
                       nelements: int, ntotal: int,
                       is_magnitude_phase: bool) -> int:
    op2 = self
    #ntotal1 = 132 * self.factor
    n = 0
    s = Struct(mapfmt(op2._endian + op2._analysis_code_fmt + b'32f', op2.size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        n += ntotal
        out = s.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_Shear - %s\n' % (str(out)))
        (eid_device,
         f41r, f21r, f12r, f32r, f23r, f43r, f34r, f14r, # 8
         kf1r, s12r, kf2r, s23r, kf3r, s34r, kf4r, s41r, # 16
         f41i, f21i, f12i, f32i, f23i, f43i, f34i, f14i,
         kf1i, s12i, kf2i, s23i, kf3i, s34i, kf4i, s41i) = out
        if is_magnitude_phase:
            f41 = polar_to_real_imag(f41r, f41i)
            kf1 = polar_to_real_imag(kf1r, kf1i)
            f21 = polar_to_real_imag(f21r, f21i)
            kf2 = polar_to_real_imag(kf2r, kf2i)
            f12 = polar_to_real_imag(f12r, f12i)
            kf3 = polar_to_real_imag(kf3r, kf3i)
            f23 = polar_to_real_imag(f23r, f23i)
            kf4 = polar_to_real_imag(kf4r, kf4i)
            f32 = polar_to_real_imag(f32r, f32i)
            s12 = polar_to_real_imag(s12r, s12i)
            f43 = polar_to_real_imag(f43r, f43i)
            s23 = polar_to_real_imag(s23r, s23i)
            f34 = polar_to_real_imag(f34r, f34i)
            s34 = polar_to_real_imag(s34r, s34i)
            f14 = polar_to_real_imag(f14r, f14i)
            s41 = polar_to_real_imag(s41r, s41i)
        else:
            f41 = complex(f41r, f41i)
            kf1 = complex(kf1r, kf1i)
            f21 = complex(f21r, f21i)
            kf2 = complex(kf2r, kf2i)
            f12 = complex(f12r, f12i)
            kf3 = complex(kf3r, kf3i)
            f23 = complex(f23r, f23i)
            kf4 = complex(kf4r, kf4i)
            f32 = complex(f32r, f32i)
            s12 = complex(s12r, s12i)
            f43 = complex(f43r, f43i)
            s23 = complex(s23r, s23i)
            f34 = complex(f34r, f34i)
            s34 = complex(s34r, s34i)
            f14 = complex(f14r, f14i)
            s41 = complex(s41r, s41i)

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        add_sort_x(dt, eid,
                   f41, f21, f12, f32, f23, f43, f34, f14,
                   kf1, s12, kf2, s23, kf3, s34, kf4, s41)
    return n

def oef_cquad4_33_real_9(self, data: bytes,
                         obj: RealPlateForceArray,
                         nelements: int, ntotal: int) -> int:
    op2 = self
    n = 0
    s = Struct(mapfmt(op2._endian + op2._analysis_code_fmt + b'8f', op2.size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = s.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('real_OEF_Plate-%s - %s\n' % (op2.element_type, str(out)))
        (eid_device, mx, my, mxy, bmx, bmy, bmxy, tx, ty) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        add_sort_x(dt, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)
        n += ntotal
    return n

def oef_cquad4_33_imag_17(self, data: bytes, ndata: int,
                          obj: ComplexPlateForceArray,
                          nelements: int, ntotal: int,
                          is_magnitude_phase: bool) -> int:
    op2 = self
    n = 0
    s = Struct(mapfmt(op2._endian + op2._analysis_code_fmt + b'16f', op2.size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = s.unpack(edata)
        (eid_device,
         mxr, myr, mxyr, bmxr, bmyr, bmxyr, txr, tyr,
         mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('complex_OEF_Plate-%s - %s\n' % (op2.element_type, str(out)))

        if is_magnitude_phase:
            mx = polar_to_real_imag(mxr, mxi)
            my = polar_to_real_imag(myr, myi)
            mxy = polar_to_real_imag(mxyr, mxyi)
            bmx = polar_to_real_imag(bmxr, bmxi)
            bmy = polar_to_real_imag(bmyr, bmyi)
            bmxy = polar_to_real_imag(bmxyr, bmxyi)
            tx = polar_to_real_imag(txr, txi)
            ty = polar_to_real_imag(tyr, tyi)
        else:
            mx = complex(mxr, mxi)
            my = complex(myr, myi)
            mxy = complex(mxyr, mxyi)
            bmx = complex(bmxr, bmxi)
            bmy = complex(bmyr, bmyi)
            bmxy = complex(bmxyr, bmxyi)
            tx = complex(txr, txi)
            ty = complex(tyr, tyi)
        add_sort_x(dt, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)
        n += ntotal
    return n

def oef_cquad4_144_real_9(op2: OP2, data: bytes,
                          obj: RealPlateBilinearForceArray,
                          nelements: int, nnodes: int) -> int:
    n = 0
    n44 = 44 * op2.factor
    n36 = 36 * op2.factor
    if op2.size == 4:
        fmt1 = op2._endian + op2._analysis_code_fmt + b'4si8f'  # 8+36
        fmt2 = op2._endian + b'i8f' # 36
    else:
        fmt1 = op2._endian + mapfmt(op2._analysis_code_fmt, 8) + b'8sq8d'
        fmt2 = op2._endian + b'q8d'
    s1 = Struct(fmt1)
    s2 = Struct(fmt2)

    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + n44]

        out = s1.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_Plate2-%s - %s\n' % (op2.element_type, str(out)))
        (eid_device, term, _nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty) = out
        #term= 'CEN\'
        #_nid = 4
        # -> CEN/4
        nid = 0
        inode = 0
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        add_sort_x(dt, eid, term,
                   inode, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)
        n += n44
        for jnode in range(nnodes):
            edata = data[n : n + n36]
            out = s2.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('    %s\n' % (str(out)))
            (nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty) = out
            assert nid > 0, 'nid=%s' % nid
            add_sort_x(dt, eid, term,
                       jnode+1, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)
            n += n36
    return n

def oef_cquad4_imag_17(self, data: bytes, ndata: int,
                       obj: ComplexPlate2ForceArray,
                       nelements: int, nnodes: int,
                       is_magnitude_phase: bool) -> int:
    op2 = self
    n = 0
    if self.size == 4:
        s1 = Struct(op2._endian + b'i4s17f')  # 2+17=19 * 4 = 76
        s2 = Struct(op2._endian + b'i16f')  # 17 * 4 = 68
    else:
        s1 = Struct(op2._endian + b'q8s17d')  # 2+17=19 * 4 = 768
        s2 = Struct(op2._endian + b'q16d')  # 17 * 4 = 68
    ntotal = (8 + (nnodes + 1) * 68) * self.factor
    ntotal1 = 76 * self.factor
    ntotal2 = 68 * self.factor

    nelements = ndata // ntotal
    obj = op2.obj
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + ntotal1]
        n += ntotal1

        out = s1.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_Plate2-%s - %s\n' % (op2.element_type, str(out)))
        (eid_device, term, nid,
         mxr, myr, mxyr, bmxr, bmyr, bmxyr, txr, tyr,
         mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi) = out
        #term = 'CEN\'

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if is_magnitude_phase:
            mx = polar_to_real_imag(mxr, mxi)
            my = polar_to_real_imag(myr, myi)
            mxy = polar_to_real_imag(mxyr, mxyi)
            bmx = polar_to_real_imag(bmxr, bmxi)
            bmy = polar_to_real_imag(bmyr, bmyi)
            bmxy = polar_to_real_imag(bmxyr, bmxyi)
            tx = polar_to_real_imag(txr, txi)
            ty = polar_to_real_imag(tyr, tyi)
        else:
            mx = complex(mxr, mxi)
            my = complex(myr, myi)
            mxy = complex(mxyr, mxyi)
            bmx = complex(bmxr, bmxi)
            bmy = complex(bmyr, bmyi)
            bmxy = complex(bmxyr, bmxyi)
            tx = complex(txr, txi)
            ty = complex(tyr, tyi)
        obj.add_new_element_sort1(dt, eid, term, nid, mx, my, mxy,
                                  bmx, bmy, bmxy, tx, ty)

        for unused_j in range(nnodes):  # .. todo:: fix crash...
            edata = data[n:n+ntotal2]
            n += ntotal2
            out = s2.unpack(edata)
            (nid,
             mxr, myr, mxyr, bmxr, bmyr, bmxyr, txr, tyr,
             mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi) = out
            if is_magnitude_phase:
                mx = polar_to_real_imag(mxr, mxi)
                my = polar_to_real_imag(myr, myi)
                mxy = polar_to_real_imag(mxyr, mxyi)
                bmx = polar_to_real_imag(bmxr, bmxi)
                bmy = polar_to_real_imag(bmyr, bmyi)
                bmxy = polar_to_real_imag(bmxyr, bmxyi)
                tx = polar_to_real_imag(txr, txi)
                ty = polar_to_real_imag(tyr, tyi)
            else:
                mx = complex(mxr, mxi)
                my = complex(myr, myi)
                mxy = complex(mxyr, mxyi)
                bmx = complex(bmxr, bmxi)
                bmy = complex(bmyr, bmyi)
                bmxy = complex(bmxyr, bmxyi)
                tx = complex(txr, txi)
                ty = complex(tyr, tyi)
            if op2.is_debug_file:
                op2.binary_debug.write('OEF_Plate2 - eid=%i nid=%s out=%s\n' % (
                    eid, nid, str(out)))
            add_sort_x(dt, eid, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)
    return n

def oef_cconeax_real_7(self, data: bytes,
                       obj: Union[RealConeAxForceArray],
                       nelements: int, ntotal: int, dt: Any) -> int:
    op2 = self
    n = 0
    s = Struct(op2._endian + op2._analysis_code_fmt + b'6f')
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = s.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_CONEAX-35 - %s\n' % (str(out)))
        (eid_device, hopa, bmu, bmv, tm, su, sv) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        add_sort_x(dt, eid, hopa, bmu, bmv, tm, su, sv)
        n += ntotal
    return n

def oef_csolid_pressure_10(self, data: bytes,
                           obj: Union[RealSolidPressureForceArray],
                           nelements: int, ntotal: int, dt: Any) -> int:
    op2 = self
    n = 0
    if self.size == 4:
        fmt = op2._endian + op2._analysis_code_fmt + b'8s7f'
    else:
        fmt = op2._endian + mapfmt(op2._analysis_code_fmt, self.size) + b'16s7d'

    s = Struct(fmt)
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n : n + ntotal]
        n += ntotal
        out = s.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_PentaPressure-%s %s\n' % (op2.element_type, str(out)))
        (eid_device, ename, ax, ay, az, vx, vy, vz, pressure) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        add_sort_x(dt, eid, ename, ax, ay, az, vx, vy, vz, pressure)
    return n

def oef_csolid_imag_16(self, data: bytes,
                       obj: ComplexSolidPressureForceArray,
                       nelements: int, ntotal: int,
                       is_magnitude_phase: bool) -> int:
    op2 = self
    n = 0
    if self.size == 4:
        s = Struct(op2._endian + op2._analysis_code_fmt + b'8s 13f')
    else:
        s = Struct(mapfmt(op2._endian + op2._analysis_code_fmt, self.size) + b'16s 13d')

    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        n += ntotal

        #print(len(edata))
        out = s.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('_oef_csolid_pressure-%s %s\n' % (op2.element_type, str(out)))
        (eid_device, ename,
         axr, ayr, azr, vxr, vyr, vzr, pressure,
         axi, ayi, azi, vxi, vyi, vzi) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        ename = ename.decode('utf-8').strip()

        if is_magnitude_phase:
            ax = polar_to_real_imag(axr, axi)
            vx = polar_to_real_imag(vxr, vxi)
            ay = polar_to_real_imag(ayr, ayi)
            vy = polar_to_real_imag(vyr, vyi)
            az = polar_to_real_imag(azr, azi)
            vz = polar_to_real_imag(vzr, vzi)
        else:
            ax = complex(axr, axi)
            vx = complex(vxr, vxi)
            ay = complex(ayr, ayi)
            vy = complex(vyr, vyi)
            az = complex(azr, azi)
            vz = complex(vzr, vzi)
        cpressure = complex(pressure, 0.)
        add_sort_x(dt, eid, ename, ax, ay, az, vx, vy, vz, cpressure)
    return n

def shock_response_prefix(thermal: int) -> str:
    #if thermal == 0:
        #prefix = 'abs.'
    if thermal == 2:
        prefix = 'abs.'  # Scaled response spectra ABS
    elif thermal == 4:
        #D:\NASA\git\examples\move_tpl\ms103.op2
        #C:\NASA\m4\formats\git\examples\move_tpl\ms103.op2
        prefix = 'srss.' #   # Scaled response spectra SRSS
    elif thermal == 8:
        prefix = 'nrl.'  # Scaled response spectra NRL
    else:
        assert thermal in [2, 4, 8], thermal # , op2.code_information() # abs
    return prefix
