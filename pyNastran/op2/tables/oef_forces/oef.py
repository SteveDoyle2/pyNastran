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
from struct import Struct
from typing import Any

import numpy as np
from numpy import frombuffer, vstack, sin, cos, radians, array, hstack, zeros

from pyNastran.op2.op2_interface.op2_reader import mapfmt, reshape_bytes_block_strip
from pyNastran.op2.tables.utils import get_eid_dt_from_eid_device
from pyNastran.op2.op2_helper import polar_to_real_imag
from pyNastran.op2.op2_interface.op2_common import OP2Common
from pyNastran.op2.op2_interface.utils import apply_mag_phase
from pyNastran.op2.op2_interface.msc_tables import MSC_OEF_REAL_MAPPER, MSC_OEF_IMAG_MAPPER
from pyNastran.op2.op2_interface.nx_tables import NX_OEF_REAL_MAPPER, NX_OEF_IMAG_MAPPER


from pyNastran.op2.tables.oef_forces.oef_thermal_objects import (
    Real1DHeatFluxArray,
    RealHeatFlux_2D_3DArray,
    RealChbdyHeatFluxArray,
    RealConvHeatFluxArray,

    RealHeatFluxVUBeamArray,
    RealHeatFluxVU3DArray,
    RealHeatFluxVUShellArray,
)
from pyNastran.op2.tables.oef_forces.oef_force_objects import (
    FailureIndicesArray,
    RealRodForceArray, RealViscForceArray,
    RealCBarForceArray, RealCBar100ForceArray,
    RealCFastForceArrayNX, RealCWeldForceArray,
    RealCFastForceArrayMSC, RealCBushForceArray, RealCBearForceArray,
    RealPlateForceArray,
    RealPlateBilinearForceArray,
    RealSpringForceArray, RealDamperForceArray,
    RealCShearForceArray,
    RealCGapForceArray,
    RealConeAxForceArray,
    RealSolidPressureForceArray,
    RealCBeamForceArray,
    RealBendForceArray,
    RealForceVU2DArray,
    RealCBeamForceVUArray,
)
from pyNastran.op2.tables.oef_forces.oef_complex_force_objects import (
    ComplexRodForceArray,
    ComplexCBarForceArray, ComplexCWeldForceArray,
    ComplexCBeamForceArray,
    ComplexCBushForceArray,
    ComplexCBearForceArray,
    ComplexCShearForceArray,
    ComplexSpringForceArray,
    ComplexDamperForceArray,
    ComplexViscForceArray,
    ComplexPlateForceArray,
    ComplexPlate2ForceArray,
    ComplexSolidPressureForceArray,
    ComplexCBendForceArray,
    ComplexForceVU_2DArray,
    ComplexCBeamForceVUArray,
)


class OEF(OP2Common):
    """Defines OEFx table reading for element forces/heat flux"""
    def __init__(self):
        OP2Common.__init__(self)

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
        prefix = ''
        postfix = ''
        table_name_bytes = self.table_name
        if table_name_bytes in [b'OEF1X', b'OEF1', b'OEF2']:
            if self.thermal == 0:
                prefix = 'force.'
            elif self.thermal == 1:
                prefix = 'thermal_load.'
            else:
                raise NotImplementedError(self.code_information())
        elif table_name_bytes in [b'HOEF1']:
            postfix = '_flux'
        #elif self.table_name in ['OESNLXR']:
            #prefix = 'sideline_'
        #elif self.table_name in ['OESNLXD', 'OESNL1X', 'OESNLXR']:
            #prefix = 'nonlinear_'
        #elif self.table_name == 'OESNLBR':
            #prefix = 'sideline_'
        #elif self.table_name == 'OESRT':
            #prefix = 'strength_ratio.'
        #elif self.table_name in ['OESCP', 'OESTRCP']:
            #pass # TODO: update
        elif table_name_bytes in [b'OEFCRM1']:
            assert self.table_code in [4, 504], self.code_information()
            prefix = 'crm.'
        elif table_name_bytes in [b'OEFCRM2']:
            assert self.table_code in [4, 504], self.code_information()
            # sort2, random
            self.format_code = 1 # real
            self.sort_bits[0] = 0 # real
            self.sort_bits[1] = 1 # sort2
            self.sort_bits[2] = 1 # random
            self.sort_method = 2
            prefix = 'crm.'
        elif table_name_bytes in [b'OEFPSD1']:
            assert self.table_code in [4, 604], self.code_information()
            prefix = 'psd.'
        elif table_name_bytes in [b'OEFPSD2']:
            assert self.table_code in [4, 604], self.code_information()
            self.format_code = 1
            self.sort_bits[0] = 0 # real
            self.sort_bits[1] = 1 # sort2
            self.sort_bits[2] = 1 # random
            prefix = 'psd.'
        elif table_name_bytes in [b'OEFRMS1', b'OEFRMS2']:
            assert self.table_code in [4, 804], self.code_information()
            #self.format_code = 1
            self.sort_bits[0] = 0 # real
            self.sort_bits[1] = 0 # sort1
            self.sort_bits[2] = 1 # random
            self.sort_method = 1
            self._analysis_code_fmt = b'i'
            prefix = 'rms.'
        elif table_name_bytes in [b'OEFNO1', b'OEFNO2']:
            assert self.table_code in [4, 904], self.code_information()
            self.format_code = 1
            self.sort_bits[0] = 0 # real
            self.sort_bits[1] = 0 # sort1
            #self.sort_bits[0] = 1 # sort2
            self.sort_bits[2] = 1 # random
            self.sort_method = 1
            self.data_code['nonlinear_factor'] = None
            self._analysis_code_fmt = b'i'
            assert self.sort_method == 1, self.code_information()
            prefix = 'no.'
        elif table_name_bytes in [b'OEFATO1', b'OEFATO2']:
            assert self.table_code in [4], self.code_information()
            prefix = 'ato.'

        elif table_name_bytes in [b'RAFCONS']:
            prefix = 'RAFCONS.'
        elif table_name_bytes in [b'RAFEATC']:
            prefix = 'RAFEATC.'
        elif table_name_bytes in [b'DOEF1']:
            assert self.table_code in [4], self.code_information()
            if self.thermal == 0:
                postfix = '_abs'
            elif self.thermal == 2:
                postfix = '_abs'  # Scaled response spectra ABS
            elif self.thermal == 4:
                #C:\NASA\m4\formats\git\examples\move_tpl\ms103.op2
                postfix = '_srss' #   # Scaled response spectra SRSS
            elif self.thermal == 8:
                postfix = '_nrl'  # Scaled response spectra NRL
            else:
                assert self.thermal in [2, 4, 8], self.code_information() # abs
        elif table_name_bytes in [b'OEFIT']:
            assert self.table_code in [25], self.code_information()
            prefix = 'failure_indices.'
            #raise NotImplementedError(self.code_information())
        elif table_name_bytes in [b'OEFITSTN']: # composite failure indicies
            assert self.table_code in [25], self.code_information()
            prefix = 'failure_indices.'
        else:
            raise NotImplementedError('%r' % self.table_name)
        return prefix, postfix

    def _oef_force_code(self):
        """
        Gets the numwide codes for the element to determine if
        the real or complex result should be found.
        The format and sort codes do not always give the right answer...
        """
        if self.is_nx:
            real_mapper = NX_OEF_REAL_MAPPER
            imag_mapper = NX_OEF_IMAG_MAPPER
        else:
            real_mapper = MSC_OEF_REAL_MAPPER
            imag_mapper = MSC_OEF_IMAG_MAPPER


        try:
            real = real_mapper[self.element_type]
        except KeyError:
            real = None

        try:
            imag = imag_mapper[self.element_type]
        except KeyError:
            imag = None
        return (real, imag)

    def _read_oef1_3(self, data, ndata):
        """Table 3 parser for OEF1 table"""
        self._analysis_code_fmt = b'i'
        self._data_factor = 1
        self.words = [
            'aCode', 'tCode', 'element_type', 'isubcase',
            '???', '???', '???', '???',
            'format_code', 'num_wide', 'o_code', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', 'Title', 'subtitle', 'label']

        self.parse_approach_code(data)

        #: element type
        self.element_type = self.add_data_parameter(data, 'element_type', b'i', 3, False)

        # dynamic load set ID/random code
        #self.dLoadID = self.add_data_parameter(data, 'dLoadID', b'i', 8, False)

        #: format code
        self.format_code = self.add_data_parameter(data, 'format_code', b'i', 9, False)

        #: number of words per entry in record
        #: .. note: is this needed for this table ???
        self.num_wide = self.add_data_parameter(data, 'num_wide', b'i', 10, False)

        #: undefined in DMAP...
        self.o_code = self.add_data_parameter(data, 'o_code', b'i', 11, False)

        #: thermal flag; 1 for heat ransfer, 0 otherwise
        self.thermal = self.add_data_parameter(data, 'thermal', b'i', 23, False)

        ## assuming tCode=1
        if self.analysis_code == 1:   # statics
            self.loadID = self.add_data_parameter(data, 'loadID', b'i', 5, False)  # load set ID number
            self.data_names = self.apply_data_code_value('data_names', ['loadID'])
            self.setNullNonlinearFactor()
        elif self.analysis_code == 2:  # normal modes/buckling (real eigenvalues)
            #: mode number
            self.mode = self.add_data_parameter(data, 'mode', b'i', 5)
            #: eigenvalue
            self.eign = self.add_data_parameter(data, 'eign', b'f', 6, False)
            self.cycle = 0.
            self.update_mode_cycle('cycle')
            self.data_names = self.apply_data_code_value('data_names', ['mode', 'eign', 'cycle'])
            # TODO: mode_cycle is not defined?
            #self.data_names = self.apply_data_code_value('data_names', ['mode', 'eign', 'mode_cycle'])
        elif self.analysis_code == 3:  # differential stiffness 0
            #: load set ID number
            self.loadID = self.add_data_parameter(data, 'loadID', b'i', 5)
            self.data_names = self.apply_data_code_value('data_names', ['loadID'])
        elif self.analysis_code == 4:  # differential stiffness 1
            #: load set ID number
            self.loadID = self.add_data_parameter(data, 'loadID', b'i', 5)
            self.data_names = self.apply_data_code_value('data_names', ['loadID'])
        elif self.analysis_code == 5:   # frequency
            self.freq = self.add_data_parameter(data, 'freq', b'f', 5)  # frequency
            self.data_names = self.apply_data_code_value('data_names', ['freq'])
        elif self.analysis_code == 6:  # transient
            self.time = self.add_data_parameter(data, 'time', b'f', 5)  # time step
            self.data_names = self.apply_data_code_value('data_names', ['time'])
        elif self.analysis_code == 7:  # pre-buckling
            #: load set ID number
            self.loadID = self.add_data_parameter(data, 'loadID', b'i', 5)
            #self.apply_data_code_value('data_names',['lsdvmn'])
            self.data_names = self.apply_data_code_value('data_names', ['loadID'])
        elif self.analysis_code == 8:  # post-buckling
            #: load set ID number
            self.loadID = self.add_data_parameter(data, 'loadID', b'i', 5)
            #: real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', b'f', 6, False)
            self.data_names = self.apply_data_code_value('data_names', ['loadID', 'eigr'])
        elif self.analysis_code == 9:  # complex eigenvalues
            #: mode number
            self.mode = self.add_data_parameter(data, 'mode', b'i', 5)
            #: real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', b'f', 6, False)
            #: imaginary eigenvalue
            self.eigi = self.add_data_parameter(data, 'eigi', b'f', 7, False)
            self.data_names = self.apply_data_code_value('data_names', ['mode', 'eigr', 'eigi'])
        elif self.analysis_code == 10:  # nonlinear statics
            #: load step
            self.load_step = self.add_data_parameter(data, 'load_step', b'f', 5)
            self.data_names = self.apply_data_code_value('data_names', ['load_step'])
        elif self.analysis_code == 11:  # geometric nonlinear statics
            #: load set ID number
            self.loadID = self.add_data_parameter(data, 'loadID', b'i', 5)
            self.data_names = self.apply_data_code_value('data_names', ['loadID'])
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % str(self.analysis_code)
            raise RuntimeError(msg)

        self.fix_format_code()
        self._parse_thermal_code()
        try:
            self.element_name = self.element_mapper[self.element_type]
        except KeyError:
            self.log.error(self.code_information())
            raise
        assert self.element_name != '', self.code_information()

        #self.element_name = self.element_mapper[self.element_type]
        self.data_code['element_name'] = self.element_name

        if self.is_debug_file:
            self.binary_debug.write('  %-14s = %r\n' % ('element_name', self.element_name))
            self.binary_debug.write('  %-14s = %r %s\n' % ('approach_code', self.approach_code,
                                                           self.approach_code_str(self.approach_code)))
            self.binary_debug.write('  %-14s = %r\n' % ('tCode', self.tCode))
            self.binary_debug.write('  %-14s = %r\n' % ('isubcase', self.isubcase))


        self._read_title(data)
        if self.element_type not in self.element_mapper:
            msg = 'element_type = %s' % self.element_type
            return self._not_implemented_or_skip(data, ndata, msg)
        self._write_debug_bits()
        assert self.num_wide != 146, self.code_information()
        #print('OEF-%s' % self.element_name)

    def _read_oef2_3(self, data, unused_ndata):
        """Table 3 parser for OEF2 table"""
        self._data_factor = 1
        self.words = [
            'aCode', 'tCode', 'element_type', 'isubcase',
            '???', '???', '???', '???',
            'format_code', 'num_wide', 'o_code', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', 'Title', 'subtitle', 'label']

        self.parse_approach_code(data)  # 3
        self.sort_method = 2

        #: element type
        self.element_type = self.add_data_parameter(data, 'element_type', b'i', 3, False)

        # dynamic load set ID/random code
        #self.dLoadID = self.add_data_parameter(data, 'dLoadID', b'i', 8, False)

        #: format code
        self.format_code = self.add_data_parameter(data, 'format_code', b'i', 9, False)

        #: number of words per entry in record
        #: .. note: is this needed for this table ???
        self.num_wide = self.add_data_parameter(data, 'num_wide', b'i', 10, False)

        #: undefined in DMAP...
        self.o_code = self.add_data_parameter(data, 'o_code', b'i', 11, False)

        #: thermal flag; 1 for heat ransfer, 0 otherwise
        self.thermal = self.add_data_parameter(data, 'thermal', b'i', 23, False)

        self.element_id = self.add_data_parameter(data, 'element_id', b'i', 5, fix_device_code=True)
        self._element_id = self.add_data_parameter(data, '_element_id', b'f', 5, apply_nonlinear_factor=False, add_to_dict=True)

        if self.analysis_code == 1:  # static...because reasons.
            self._analysis_code_fmt = b'f'
            self.data_names = self.apply_data_code_value('data_names', ['element_id'])
        elif self.analysis_code == 2:  # real eigenvalues
            self._analysis_code_fmt = b'i'
            self.eign = self.add_data_parameter(data, 'eign', b'f', 6, False)
            self.mode_cycle = self.add_data_parameter(data, 'mode_cycle', b'i', 7, False)  # mode or cycle .. todo:: confused on the type - F1???
            self.data_names = self.apply_data_code_value('data_names', ['element_id', 'eign', 'mode_cycle'])
        elif self.analysis_code == 5:   # frequency
            self._analysis_code_fmt = b'f'
            self.data_names = self.apply_data_code_value('data_names', ['element_id'])
            self.apply_data_code_value('analysis_method', 'freq')
        elif self.analysis_code == 6:  # transient
            self._analysis_code_fmt = b'f'
            self.data_names = self.apply_data_code_value('data_names', ['element_id'])
            self.apply_data_code_value('analysis_method', 'time')
        elif self.analysis_code == 7:  # pre-buckling
            self._analysis_code_fmt = b'i'
            self.data_names = self.apply_data_code_value('data_names', ['element_id'])
            self.apply_data_code_value('analysis_method', 'lsdvmn')
        elif self.analysis_code == 8:  # post-buckling
            self._analysis_code_fmt = b'f'
            self.eigr = self.add_data_parameter(data, 'eigr', b'f', 6, False)
            self.data_names = self.apply_data_code_value('data_names', ['element_id', 'eigr'])
            self.apply_data_code_value('analysis_method', 'eigr')
        elif self.analysis_code == 9:  # complex eigenvalues
            # mode number
            self._analysis_code_fmt = b'i'
            self.eigr = self.add_data_parameter(data, 'eigr', b'f', 6, False)
            self.eigi = self.add_data_parameter(data, 'eigi', b'f', 7, False)
            self.data_names = self.apply_data_code_value('data_names', ['element_id', 'eigr', 'eigi'])
            self.apply_data_code_value('analysis_method', 'mode')
        elif self.analysis_code == 10:  # nonlinear statics
            # load step
            self._analysis_code_fmt = b'f'
            self.data_names = self.apply_data_code_value('data_names', ['element_id'])
            self.apply_data_code_value('analysis_method', 'lftsfq')
        elif self.analysis_code == 11:  # old geometric nonlinear statics
            # load set number
            self.data_names = self.apply_data_code_value('data_names', ['element_id'])
        elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            self.data_names = self.apply_data_code_value('data_names', ['element_id'])
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % self.analysis_code
            raise RuntimeError(msg)

        self._fix_oes_sort2(data)
        self._set_element_name()
        #assert isinstance(self.nonlinear_factor, int), self.nonlinear_factor

    def _read_oef1_4(self, data, ndata):
        """Table 4 parser for OEF1 table"""
        if self.thermal == 0:
            self._setup_op2_subcase('FORCE')
            n = self._read_oef1_loads(data, ndata)
        elif self.thermal == 1:
            n = self._read_oef1_thermal(data, ndata)
        elif self.thermal in [2, 4, 8]: # 2=ABS, 4=SRSS, 8=NRL
            #C:\NASA\m4\formats\git\examples\move_tpl\ms103.op2 # SRSS
            n = self._read_oef1_loads(data, ndata)
        else:
            n = self._not_implemented_or_skip(data, ndata, 'thermal=%s' % self.thermal)
        return n

    def _read_oef2_4(self, data, ndata):
        if self.thermal == 0: # and self.element_type not in [77]:
            self._setup_op2_subcase('FORCE')
            n = self._read_oef1_loads(data, ndata)
        else:
            n = self._not_implemented_or_skip(data, ndata, 'thermal=%s' % self.thermal)
        return n

    def _read_oef1_thermal(self, data, ndata):
        """Table 4 parser for OEF1 thermal table"""
        if self._results.is_not_saved('element_forces'):
            return ndata
        prefix, postfix = self.get_oef_prefix_postfix()

        n = 0
        #thermal
        #is_magnitude_phase = self.is_magnitude_phase()
        dt = self.nonlinear_factor

        #flag = 'element_id'
        if self.element_type in [1, 2, 3, 10, 34, 69]:  # ROD,BEAM,TUBE,CONROD,BAR,BEND
            n, nelements, ntotal = self._thermal_1d_heat_flux(data, ndata, dt, prefix, postfix)
            if nelements is None:
                return n

        elif self.element_type in [33, 53, 64, 74, 75,  # CQUAD4, CTRIAX6, CQUAD8, CTRIA3, CTRIA6
                                   39, 67, 68]:  # TETRA, HEXA, PENTA
            n, nelements, ntotal = self._thermal_2d_3d_heat_flux(data, ndata, dt, prefix, postfix)

        elif self.element_type in [107, 108, 109]:  # CHBDYE, CHBDYG, CHBDYP
            n, nelements, ntotal = self._thermal_chbdy(data, ndata, dt, prefix, postfix)

        elif self.element_type == 110:  # CONV
            n, nelements, ntotal = self._thermal_conv(data, ndata, dt, prefix, postfix)

        elif self.element_type in [145, 146, 147]:  # VUHEXA,VUPENTA,VUTETRA
            n, nelements, ntotal = self._thermal_vu_solid(data, ndata, dt, prefix, postfix)

        elif self.element_type in [189, 190]:  # VUQUAD,VUTRIA
            n, nelements, ntotal = self._thermal_vu_shell(data, ndata, dt, prefix, postfix)

        elif self.element_type == 191:  # VUBEAM
            n, nelements, ntotal = self._thermal_vu_beam(data, ndata, dt, prefix, postfix)
        else:
            msg = 'OEF sort1 thermal Type=%s num=%s' % (self.element_name, self.element_type)
            return self._not_implemented_or_skip(data, ndata, msg)
        if nelements is None:
            return n


        assert self.thermal == 1, self.thermal
        assert ndata > 0, ndata
        try:
            assert nelements > 0, 'nelements=%r element_type=%s element_name=%r\n%s' % (nelements, self.element_type, self.element_name, self.code_information())
        except UnboundLocalError:
            raise UnboundLocalError('element_name=%r' % self.element_name)
        #assert ndata % ntotal == 0, '%s n=%s nwide=%s len=%s ntotal=%s' % (self.element_name, ndata % ntotal, ndata % self.num_wide, ndata, ntotal)
        assert self.num_wide * 4 == ntotal, 'numwide*4=%s ntotal=%s' % (self.num_wide*4, ntotal)
        assert n > 0, 'n=%s element_type=%s element_name=%s numwide=%s' % (
            n, self.element_type, self.element_name, self.num_wide)
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
        n = 0
        obj_vector_real = Real1DHeatFluxArray
        #if self.element_type == 1: # CROD
        if self.element_type == 1:
            result_name = prefix + 'crod_thermal_load' + postfix
        elif self.element_type == 2:
            result_name = prefix + 'cbeam_thermal_load' + postfix
        elif self.element_type == 3:
            result_name = prefix + 'ctube_thermal_load' + postfix
        elif self.element_type == 10:
            result_name = prefix + 'conrod_thermal_load' + postfix
        elif self.element_type == 34:
            result_name = prefix + 'cbar_thermal_load' + postfix
        elif self.element_type == 69:
            result_name = prefix + 'cbend_thermal_load' + postfix
        else:
            raise NotImplementedError('element_type=%s element_name=%s' % (
                self.element_type, self.element_name))

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)

        if self.format_code == 1 and self.num_wide == 9:  # real
            ntotal = 36
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return nelements * self.num_wide * 4, None, None
            obj = self.obj
            #if self.is_debug_file:
                #self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                #self.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
                #self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * 4 * self.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 9)
                obj._times[obj.itime] = dt

                strings = frombuffer(data, dtype=self._uendian + 'S4').reshape(nelements, 9)
                s = array([s1+s2 for s1, s2 in zip(strings[:, 1], strings[:, 2])])
                #print(s)
                #print('itime = ', obj.itime)
                #print('---------')
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 9)
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
                s = Struct(self._endian + self._analysis_code_fmt + b'8s6f')
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]
                    out = s.unpack(edata)
                    (eid_device, eType, xgrad, ygrad, zgrad, xflux, yflux, zflux) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    obj.add_sort1(dt, eid, eType, xgrad, ygrad, zgrad, xflux, yflux, zflux)
                    n += ntotal
        else:  # pragma: no cover
            msg = self.code_information()
            return self._not_implemented_or_skip(data, ndata, msg), None, None
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
        n = 0
        if self.element_type == 33:
            result_name = prefix + 'cquad4_thermal_load' + postfix
        elif self.element_type == 53:
            result_name = prefix + 'ctriax6_thermal_load' + postfix
        elif self.element_type == 64:
            result_name = prefix + 'cquad8_thermal_load' + postfix
        elif self.element_type == 74:
            result_name = prefix + 'ctria3_thermal_load' + postfix
        elif self.element_type == 75:
            result_name = prefix + 'ctria6_thermal_load' + postfix
        elif self.element_type == 39:
            result_name = prefix + 'ctetra_thermal_load' + postfix
        elif self.element_type == 67:
            result_name = prefix + 'chexa_thermal_load' + postfix
        elif self.element_type == 68:
            result_name = prefix + 'cpenta_thermal_load' + postfix
        else:
            raise NotImplementedError('element_type=%s element_name=%s' % (
                self.element_type, self.element_name))

        if self._results.is_not_saved(result_name):
            return ndata, None, None

        obj_vector_real = RealHeatFlux_2D_3DArray
        #if self.element_type == 1: # CROD
        #result_name = 'thermalLoad_2D_3D'

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)
        if self.format_code == 1 and self.num_wide == 9:  # real - 2D
            # [33, 53, 64, 74, 75]
            ntotal = 4 * self.num_wide * self.factor
            ntotal = 36 * self.factor
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return nelements * ntotal, None, None
            obj = self.obj
            #if self.is_debug_file:
                #self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                #self.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
                #self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 9)
                obj._times[obj.itime] = dt
                #if obj.itime == 0:
                ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 9)
                eids = ints[:, 0] // 10
                assert eids.min() > 0, eids.min()
                obj.element[itotal:itotal2] = eids
                strings = frombuffer(data, dtype=self._uendian + 'S4').reshape(nelements, 9)
                obj.element_data_type[itotal:itotal2] = array([s1+s2 for s1, s2 in zip(strings[:, 1], strings[:, 2])])

                #[etype, xgrad, ygrad, zgrad, xflux, yflux, zflux]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 3:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                # no zed on this element for some reason...
                if self.size == 4:
                    fmt = self._endian + self._analysis_code_fmt + b'8s 6f'
                else:
                    fmt = self._endian + mapfmt(self._analysis_code_fmt, self.size) + b'16s 6d'
                s = Struct(fmt)
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]
                    n += ntotal
                    out = s.unpack(edata)
                    (eid_device, etype, xgrad, ygrad, zgrad, xflux, yflux, zflux) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    obj.add_sort1(dt, eid, etype, xgrad, ygrad, zgrad, xflux, yflux, zflux)

        elif self.format_code == 1 and self.num_wide == 10:  # real - 3D
            # [39, 67, 68]:  # HEXA,PENTA
            ntotal = 40 * self.factor
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return nelements * ntotal, None, None
            obj = self.obj
            assert nelements > 0, 'ndata=%s ntotal=%s' % (ndata, ntotal)
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 10)
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 10)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids
                    strings = frombuffer(data, dtype=self._uendian + 'S4').reshape(nelements, 10)
                    obj.element_data_type[itotal:itotal2] = array([s1+s2 for s1, s2 in zip(strings[:, 1], strings[:, 2])])

                #[etype, xgrad, ygrad, zgrad, xflux, yflux, zflux, zed]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 3:-1].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                s = Struct(self._endian + self._analysis_code_fmt + b'8s6fi')
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]
                    n += ntotal
                    out = s.unpack(edata)
                    (eid_device, etype, xgrad, ygrad, zgrad, xflux, yflux, zflux, unused_zed) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    obj.add_sort1(dt, eid, etype, xgrad, ygrad, zgrad, xflux, yflux, zflux)
        else:
            raise RuntimeError(self.code_information())
        return n, nelements, ntotal

    def _thermal_chbdy(self, data, ndata, dt, prefix, postfix):
        """
        107-CHBDYE
        108-CHBDYG
        109-CHBDYP

        """
        n = 0
        #if self.table_name in ['OEF1X']:
        if self.element_type == 107:
            result_name = prefix + 'chbdye_thermal_load' + postfix
        elif self.element_type == 108:
            result_name = prefix + 'chbdyg_thermal_load' + postfix
        elif self.element_type == 109:
            result_name = prefix + 'chbdyp_thermal_load' + postfix
        else:
            raise NotImplementedError('element_type=%s element_name=%s' % (
                self.element_type, self.element_name))

        if self._results.is_not_saved(result_name):
            return ndata, None, None

        #elif self.table_name in ['HOEF1']:
            #if self.element_type == 107:
                #result_name = 'chbdye_thermal_flux'
            #elif self.element_type == 108:
                #result_name = 'chbdyg_thermal_flux'
            #elif self.element_type == 109:
                #result_name = 'chbdyp_thermal_flux'
            #else:
                #raise NotImplementedError('element_type=%s element_name=%s' % (
                    #self.element_type, self.element_name))
        #else:
            #raise NotImplementedError(msg)

        if self.format_code == 1 and self.num_wide == 8:  # real
            #result_name = 'thermalLoad_CHBDY'
            if self._results.is_not_saved(result_name):
                return ndata, None, None
            self._results._found_result(result_name)
            slot = self.get_result(result_name)

            if self.format_code == 1 and self.num_wide == 8:  # real
                obj_vector_real = RealChbdyHeatFluxArray
                ntotal = 32
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    return nelements * self.num_wide * 4, None, None
                obj = self.obj
                #if self.is_debug_file:
                    #self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    #self.binary_debug.write('  #elementi = [eid_device, etype, fapplied, free_conv, force_conv, frad, ftotal]\n')
                    #self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                if self.use_vector and is_vectorized and self.sort_method == 1:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 8)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 8)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids
                        #obj.element_type[obj.itime, itotal:itotal2, :] = strings[:, 3:]

                    #[fapplied, free_conv, force_conv, frad, ftotal]
                    obj.data[obj.itime, itotal:itotal2, :] = floats[:, 3:].copy()
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s1 = Struct(self._endian + self._analysis_code_fmt + b'8s5f')
                    for unused_i in range(nelements):
                        edata = data[n:n+32]
                        n += ntotal
                        out = s1.unpack(edata)
                        (eid_device, etype, fapplied, free_conv, force_conv, frad, ftotal) = out
                        eid, dt = get_eid_dt_from_eid_device(
                            eid_device, self.nonlinear_factor, self.sort_method)

                        if self.is_debug_file:
                            self.binary_debug.write('  %s -> [%s, %s, %s, %s, %s, %s, %s]\n'
                                                    % (eid, eid_device, etype, fapplied, free_conv, force_conv, frad, ftotal))
                        obj.add_sort1(dt, eid, etype, fapplied, free_conv, force_conv, frad, ftotal)
        else:  # pragma: no cover
            msg = self.code_information()
            return self._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _thermal_conv(self, data, ndata, dt, prefix, postfix):
        n = 0
        # 110-CONV
        result_name = prefix + 'conv_thermal_load' + postfix
        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)

        if self.format_code == 1 and self.num_wide == 4:
            ntotal = 16
            nelements = ndata // ntotal

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, RealConvHeatFluxArray)
            if auto_return:
                return nelements * self.num_wide * 4, None, None
            obj = self.obj
            #if self.is_debug_file:
                #self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                #self.binary_debug.write('  #elementi = [eid_device, etype, fapplied, free_conv, force_conv, frad, ftotal]\n')
                #self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * 4 * self.num_wide
                ielement = obj.ielement
                ielement2 = ielement + nelements

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 4).copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 4).copy()
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
                s1 = Struct(self._endian + self._analysis_code_fmt + b'fif')
                for unused_i in range(nelements):
                    edata = data[n:n+16]
                    n += 16
                    out = s1.unpack(edata)
                    (eid_device, free_conv, cntl_node, free_conv_k) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    assert cntl_node >= 0, cntl_node
                    obj.add_sort1(dt, eid, cntl_node, free_conv, free_conv_k)
        else:  # pragma: no cover
            msg = self.code_information()
            return self._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _thermal_vu_shell(self, data, ndata, dt, prefix, postfix):
        """
        189-VUQUAD
        190-VUTRIA

        TODO: vectorize
        """
        n = 0
        if self.format_code == 1 and self.num_wide == 27:  # real
            result_name = prefix + 'thermalLoad_VU' + postfix
            if self._results.is_not_saved(result_name):
                return ndata, None, None
            self._results._found_result(result_name)
            slot = self.get_result(result_name)

            if self.element_type == 189:
                nnodes = 4
            elif self.element_type == 190:
                nnodes = 3
            elif self.element_type == 191:
                nnodes = 2
            else:
                raise NotImplementedError(self.code_information())

            numwide = 6 + 7 * nnodes
            ntotal = (24 + 28 * nnodes) * self.factor
            assert ntotal == numwide * 4 * self.factor
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, RealHeatFluxVUShellArray)
            if auto_return:
                self._data_factor = nnodes
                return nelements * ntotal, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.itotal
                itotal2 = itotal + nelements * nnodes
                ielement = obj.ielement
                ielement2 = obj.ielement + nelements

                ints = frombuffer(data, dtype=self.idtype).reshape(nelements, numwide).copy()
                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, numwide)
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    eids = ints[:, 0] // 10
                    parent = ints[:, 1]
                    coord = ints[:, 2]
                    # icord - 4s
                    theta = ints[:, 4]
                    assert eids.min() > 0, eids.min()
                    obj.element[ielement:ielement2] = eids
                    obj.element_parent_coord_icord[ielement:ielement2, 0] = eids
                    obj.element_parent_coord_icord[ielement:ielement2, 1] = parent
                    obj.element_parent_coord_icord[ielement:ielement2, 2] = coord
                    obj.element_parent_coord_icord[ielement:ielement2, 3] = theta
                #obj.data[itotal:itotal2, 0] = floats[:, 4]

                #ints2 = ints[:, 6:].reshape(nelements * nnodes, 7)
                floats2 = floats[:, 6:].reshape(nelements * nnodes, 7)

                #[vugrid]
                #[xgrad, ygrad, zgrad, xflux, yflux, zflux]
                #obj.int_data[obj.itime, itotal:itotal2, 0] = ints2[:, 0]
                obj.data[obj.itime, itotal:itotal2, :] = floats2[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                s1 = Struct(self._endian + b'3i4s2i')
                s2 = Struct(self._endian + b'i6f')
                ntotal1 = 24 * self.factor # 6*4
                ntotal2 = 28 * self.factor # 7*4
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal1]
                    n += ntotal1

                    out = s1.unpack(edata)
                    (eid_device, parent, coord, icord, theta, unused_null) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    #self.log.debug('RealHeatFluxVUArray = %s' % str(out))

                    for unused_j in range(nnodes):
                        edata = data[n:n+ntotal2]
                        n += ntotal2
                        out = s2.unpack(edata)
                        #print(out)
                        vugrid, xgrad, ygrad, zgrad, xflux, yflux, zflux = out
                        obj.add_sort1(dt, eid, parent, coord, icord, theta,
                                      xgrad, ygrad, zgrad, xflux, yflux, zflux)
        else:  # pragma: no cover
            msg = self.code_information()
            return self._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _thermal_vu_solid(self, data, ndata, dt, prefix, postfix):
        """
        145-VUHEXA
        146-VUPENTA
        147-VUTETRA

        TODO: vectorize
        """
        n = 0
        if self.element_type == 147:  # VUTETRA
            nnodes = 4
        elif self.element_type == 146:  # VUPENTA
            nnodes = 6
        elif self.element_type == 145:  # VUHEXA
            nnodes = 8
        else:
            raise NotImplementedError(msg)
            #msg = self.code_information()
            #return self._not_implemented_or_skip(data, ndata, msg), None, None

        result_name = prefix + 'thermalLoad_VU_3D' + postfix
        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)

        numwide_real = 2 + 7 * nnodes
        if self.format_code == 1 and self.num_wide == numwide_real:  # real
            ntotal = 8 + 28 * nnodes
            nelements = ndata // ntotal
            #self.create_transient_object(result_name, self.thermalLoad_VU_3D, HeatFluxVU_3DArray)
            #is_vectorized = False

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, RealHeatFluxVU3DArray)
            if auto_return:
                self._data_factor = nnodes
                return nelements * self.num_wide * 4, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * 4 * self.num_wide
                ielement = obj.ielement
                ielement2 = ielement + nelements
                itotal = obj.itotal
                itotal2 = itotal + nelements * nnodes

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, numwide_real)
                floats2 = floats[:, 2:].reshape(nelements * nnodes, 7)
                obj._times[obj.itime] = dt
                #if obj.itime == 0:
                ints = frombuffer(data, dtype=self.idtype).reshape(nelements, numwide_real).copy()
                ints2 = ints[:, 2:].reshape(nelements * nnodes, 7)
                eids = ints[:, 0] // 10
                parent = ints[:, 1]
                assert eids.min() > 0, eids.min()
                obj.element_parent[ielement:ielement2, 0] = eids
                obj.element_parent[ielement:ielement2, 1] = parent

                #[vugrid, xgrad, ygrad, zgrad, xflux, yflux, zflux]
                obj.vugrid[obj.itime, itotal:itotal2] = ints2[:, 0]
                obj.data[obj.itime, itotal:itotal2, :] = floats2[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                s1 = self.struct_2i
                s2 = Struct(self._endian + self._analysis_code_fmt + b'6f')
                for unused_i in range(nelements):
                    out = s1.unpack(data[n:n+8])
                    n += 8
                    (eid_device, parent) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)

                    grad_fluxes = []
                    for unused_j in range(nnodes):
                        out = s2.unpack(data[n:n+28])
                        grad_fluxes.append(out)
                        n += 28
                    obj.add_sort1(dt, eid, parent, grad_fluxes)
        else:  # pragma: no cover
            msg = self.code_information()
            return self._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _thermal_vu_beam(self, data, ndata, dt, prefix, postfix):
        """191-VUBEAM"""
        n = 0
        # TODO: vectorize
        nnodes = 2
        numwide_real = 4 + 7 * nnodes

        result_name = prefix + 'vu_beam_thermal_load' + postfix
        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)
        if self.format_code == 1 and self.num_wide == numwide_real:  # real
            ntotal = 16 + 28 * nnodes
            nelements = ndata // ntotal

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, RealHeatFluxVUBeamArray)
            if auto_return:
                self._data_factor = nnodes
                return nelements * self.num_wide * 4, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * 4 * self.num_wide
                itotal = obj.itotal
                itotal2 = itotal + nelements * nnodes
                ielement = obj.ielement
                ielement2 = obj.ielement + nelements

                ints = frombuffer(data, dtype=self.idtype).reshape(nelements, numwide_real).copy()
                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, numwide_real)
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    eids = ints[:, 0] // 10
                    parent = ints[:, 1]
                    coord = ints[:, 2]
                    # icord - 4s
                    assert eids.min() > 0, eids.min()
                    obj.element_parent_coord[ielement:ielement2, 0] = eids
                    obj.element_parent_coord[ielement:ielement2, 1] = parent
                    obj.element_parent_coord[ielement:ielement2, 2] = coord
                #obj.data[itotal:itotal2, 0] = floats[:, 4]

                ints2 = ints[:, 4:].reshape(nelements * nnodes, 7)
                floats2 = floats[:, 4:].reshape(nelements * nnodes, 7)

                #[vugrid]
                #[xgrad, ygrad, zgrad, xflux, yflux, zflux]
                obj.vugrid[obj.itime, itotal:itotal2, 0] = ints2[:, 0]
                obj.data[obj.itime, itotal:itotal2, :] = floats2[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                s1 = Struct(self._endian + b'iii4s')
                s2 = Struct(self._endian + b'i6f')
                for unused_i in range(nelements):
                    edata = data[n:n+16]  # 4*4
                    n += 16

                    out = s1.unpack(edata)
                    (eid_device, parent, coord, icord) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    data_in = [eid, parent, coord, icord]
                    self.log.debug('VUBeam %s' % data_in)

                    grad_fluxes = []
                    for unused_j in range(nnodes):
                        edata = data[n:n+28]  # 7*4
                        n += 28
                        out = s2.unpack(edata)
                        grad_fluxes.append(out)
                    data_in.append(grad_fluxes)
                    obj.add_sort1(dt, eid, parent, coord, icord, grad_fluxes)

        elif self.element_type == 118:  # CWELDP
            msg = 'OEF sort1 thermal Type=%s num=%s' % (self.element_name, self.element_type)
            return self._not_implemented_or_skip(data, ndata, msg)
        else:
            msg = 'OEF sort1 thermal Type=%s num=%s' % (self.element_name, self.element_type)
            return self._not_implemented_or_skip(data, ndata, msg)
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
            except:
                raise
                #print("----------")
                #try:
                    #print(self.obj)
                #except:
                    #print("error printing %r" % self.obj.__class__.__name__)
                #print(self.data_code)
                #if self.obj is not None:
                    ##from pyNastran.utils import object_attributes
                    ##print object_attributes(self.obj)
                    #print(self.obj.data_code)
                #print("----------")
                #raise
            return n
        return new_func

    def _read_oef1_loads_nasa95(self, data, ndata):
        """Reads the OEF1 table for NASA 95 Nastran"""
        if self._results.is_not_saved('element_forces'):
            return ndata
        prefix, postfix = self.get_oef_prefix_postfix()
        #print('prefix=%r postfix=%s element_name=%s' % (prefix, postfix, self.element_name))
        (num_wide_real, num_wide_imag) = self._oef_force_code()
        if self.is_debug_file:
            self.binary_debug.write('  num_wide_real = %r\n' % num_wide_real)
            self.binary_debug.write('  num_wide_imag = %r\n' % num_wide_imag)

        n = 0
        is_magnitude_phase = self.is_magnitude_phase()
        dt = self.nonlinear_factor

        if self.element_type in [1, 3, 10]:  # rods
            n, nelements, ntotal = self._oef_crod(data, ndata, dt, is_magnitude_phase, prefix, postfix)

        elif self.element_type == 2:  # cbeam
            #2-CBEAM
            n, nelements, ntotal = self._oef_cbeam(data, ndata, dt, is_magnitude_phase, prefix, postfix)

        elif self.element_type in [11, 12, 13, 14,   # springs
                                   20, 21, 22, 23]:  # dampers
            # 11-CELAS1
            # 12-CELAS2
            # 13-CELAS3
            # 14-CELAS4

            # 20-CDAMP1
            # 21-CDAMP2
            # 22-CDAMP3
            # 23-CDAMP4
            n, nelements, ntotal = self._oef_celas_cdamp(data, ndata, dt, is_magnitude_phase, prefix, postfix)

        elif self.element_type == 24:  # CVISC
            n, nelements, ntotal = self._oef_cvisc(data, ndata, dt, is_magnitude_phase, prefix, postfix)

        elif self.element_type == 34:  # cbar
            # 34-CBAR
            n, nelements, ntotal = self._oef_cbar_34(data, ndata, dt, is_magnitude_phase, prefix, postfix)

        elif self.element_type in [83]: # centroidal shells
            # 33-CQUAD4???
            # 83-CTRIA3
            n, nelements, ntotal = self._oef_shells_centroidal(data, ndata, dt, is_magnitude_phase, prefix, postfix)
        elif self.element_type == 4:  # cshear
            n, nelements, ntotal = self._oef_cshear(data, ndata, dt, is_magnitude_phase, prefix, postfix)

        elif self.element_type == 35:  # cconeax
            n, nelements, ntotal = self._oef_cconeax(data, ndata, dt, is_magnitude_phase, prefix, postfix)
        else:
            return self._not_implemented_or_skip(data, ndata, self.code_information())

        if nelements is None:
            return n

        #assert self.thermal == 0, self.thermal
        assert ndata > 0, ndata
        assert nelements > 0, 'nelements=%r element_type=%s element_name=%r num_wide=%s' % (
            nelements, self.element_type, self.element_name, self.num_wide)
        #assert ndata % ntotal == 0, '%s n=%s nwide=%s len=%s ntotal=%s' % (self.element_name, ndata % ntotal, ndata % self.num_wide, ndata, ntotal)
        assert self.num_wide * 4 == ntotal, 'numwide*4=%s ntotal=%s' % (self.num_wide*4, ntotal)
        assert n > 0, n
        return n

    # @_print_obj_name_on_crash
    def _read_oef1_loads(self, data, ndata):
        """Reads the OEF1 table; stores the element forces/heat flux."""
        #self._apply_oef_ato_crm_psd_rms_no('') # TODO: just testing
        if self._results.is_not_saved('element_forces'):
            return ndata
        prefix, postfix = self.get_oef_prefix_postfix()
        #print('prefix=%r postfix=%s element_name=%s' % (prefix, postfix, self.element_name))
        unused_flag = 'element_id'
        (num_wide_real, num_wide_imag) = self._oef_force_code()
        if self.is_debug_file:
            self.binary_debug.write('  num_wide_real = %r\n' % num_wide_real)
            self.binary_debug.write('  num_wide_imag = %r\n' % num_wide_imag)

        n = 0
        is_magnitude_phase = self.is_magnitude_phase()
        dt = self.nonlinear_factor

        if self.element_type in [1, 3, 10]:  # rods
            n, nelements, ntotal = self._oef_crod(data, ndata, dt, is_magnitude_phase, prefix, postfix)

        elif self.element_type == 2:  # cbeam
            #2-CBEAM
            n, nelements, ntotal = self._oef_cbeam(data, ndata, dt, is_magnitude_phase, prefix, postfix)

        elif self.element_type in [11, 12, 13, 14,   # springs
                                   20, 21, 22, 23]:  # dampers
            # 11-CELAS1
            # 12-CELAS2
            # 13-CELAS3
            # 14-CELAS4

            # 20-CDAMP1
            # 21-CDAMP2
            # 22-CDAMP3
            # 23-CDAMP4
            n, nelements, ntotal = self._oef_celas_cdamp(data, ndata, dt, is_magnitude_phase, prefix, postfix)

        elif self.element_type == 24:  # CVISC
            n, nelements, ntotal = self._oef_cvisc(data, ndata, dt, is_magnitude_phase, prefix, postfix)

        elif self.element_type == 34:  # cbar
            # 34-CBAR
            n, nelements, ntotal = self._oef_cbar_34(data, ndata, dt, is_magnitude_phase, prefix, postfix)

        elif self.element_type == 100:  # cbar
            #100-BARS
            n, nelements, ntotal = self._oef_cbar_100(data, ndata, dt, is_magnitude_phase, prefix, postfix)

        elif self.element_type in [33, 74, 227, 228]: # centroidal shells
            # 33-CQUAD4
            # 74-CTRIA3
            # 227-CTRIAR? (C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\cqrdbx102.op2)
            # 228-CQUADR
            n, nelements, ntotal = self._oef_shells_centroidal(data, ndata, dt, is_magnitude_phase, prefix, postfix)

        elif self.element_type in [64, 70, 75, 82, 144]: # bilinear shells
            # 64-CQUAD8
            # 70-CTRIAR
            # 75-CTRIA6
            # 82-CQUADR
            # 144-CQUAD4-bilinear
            n, nelements, ntotal = self._oef_shells_nodal(data, ndata, dt, is_magnitude_phase, prefix, postfix)

        elif self.element_type in [95, 96, 97, 98, 232, 233]: # composites
            # 95 - CQUAD4
            # 96 - CQUAD8
            # 97 - CTRIA3
            # 98 - CTRIA6 (composite)
            # 232 - CQUADR
            # 233 - CTRIAR
            n, nelements, ntotal = self._oef_shells_composite(data, ndata, dt, is_magnitude_phase,
                                                              prefix, postfix)

        elif self.element_type in [39, 67, 68]: # solids
            # 39-CTETRA
            # 67-CHEXA
            # 68-CPENTA
            if self.read_mode == 1:
                return ndata
            #self._results._found_result('solid_forces')
            raise RuntimeError(self.code_information())
            #if self.format_code == 1 and self.num_wide == 0:  # real
                ##self.create_transient_object(result_name, self.solidForces, RealCSolidForce)
                #raise RuntimeError(self.code_information())
            #else:
                #msg = self.code_information()
                #return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type == 53:  # ctriax6
            # 53-CTRIAX6
            self._results._found_result('ctriax_force')
            #if self.format_code == 1 and self.num_wide == 0:  # real
                #pass
                #self.create_transient_object(self.ctriax_force, RealCTriaxForce)  # undefined
            #else:  # pragma: no cover
            raise NotImplementedError(self.code_information())
                #msg = self.code_information()
                #return self._not_implemented_or_skip(data, ndata, msg)
            #return ndata

        elif self.element_type == 4:  # cshear
            n, nelements, ntotal = self._oef_cshear(data, ndata, dt, is_magnitude_phase, prefix, postfix)

        elif self.element_type == 35:  # coneax
            n, nelements, ntotal = self._oef_cconeax(data, ndata, dt, is_magnitude_phase, prefix, postfix)

        elif self.element_type == 38:  # cgap
            n, nelements, ntotal = self._oef_cgap(data, ndata, dt, is_magnitude_phase, prefix, postfix)

        elif self.element_type == 69:  # cbend
            n, nelements, ntotal = self._oef_cbend(data, ndata, dt, is_magnitude_phase, prefix, postfix)

        elif self.element_type in [76, 77, 78, 79]:
            # 76-HEXPR
            # 77-PENPR
            # 78-TETPR
            # 79-CPYRAM
            n, nelements, ntotal = self._oef_csolid_pressure(data, ndata, dt, is_magnitude_phase, prefix, postfix)

        elif self.element_type in [102, 280]:
            # 102: cbush
            # 280: cbear
            n, nelements, ntotal = self._oef_cbush(data, ndata, dt, is_magnitude_phase, prefix, postfix)

        elif self.element_type in [145, 146, 147]:
            # 145-VUHEXA
            # 146-VUPENTA
            # 147-VUTETRA
            #if self.read_mode == 1:
                #return ndata
            return ndata

        elif self.element_type in [189, 190]:
            n, nelements, ntotal = self._oef_vu_shell(data, ndata, dt, is_magnitude_phase, prefix, postfix)

        elif self.element_type == 191:
            n, nelements, ntotal = self._oef_vu_beam(data, ndata, dt, is_magnitude_phase, prefix, postfix)

        elif self.element_type == 126 and self.is_msc:
            # 119-CFAST-MSC
            n, nelements, ntotal = self._oef_cbush(data, ndata, dt, is_magnitude_phase, prefix, postfix)
        elif self.element_type == 119 and self.is_nx:
            # 119-CFAST-NX
            n, nelements, ntotal = self._oef_cbar_34(data, ndata, dt, is_magnitude_phase,
                                                     prefix, postfix)
        elif self.element_type in [117, 200]:
            # 117-CWELDC
            # 200-CWELD
            n, nelements, ntotal = self._oef_cbar_34(data, ndata, dt, is_magnitude_phase,
                                                     prefix, postfix)
        #elif self.element_type == 119 and self.is_msc:
            #raise NotImplementedError(self.code_information())
        elif self.element_type == 235:
            # 235-CQUADR
            return self._not_implemented_or_skip(data, ndata, self.code_information())
        elif self.is_nx:
            if self.element_type in [118, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352,
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
                #print(self.code_information())
                return self._not_implemented_or_skip(data, ndata, self.code_information())
            #print(self.code_information())
            #nx_missing
            return self._not_implemented_or_skip(data, ndata, self.code_information())
        else:
            #print(self.code_information())
            #msc_missing
            return self._not_implemented_or_skip(data, ndata, self.code_information())

        if nelements is None:
            return n

        #assert self.thermal == 0, self.thermal
        assert ndata > 0, ndata
        assert nelements > 0, 'nelements=%r element_type=%s element_name=%r num_wide=%s' % (
            nelements, self.element_type, self.element_name, self.num_wide)
        #assert ndata % ntotal == 0, '%s n=%s nwide=%s len=%s ntotal=%s' % (self.element_name, ndata % ntotal, ndata % self.num_wide, ndata, ntotal)
        assert self.num_wide * 4 * self.factor == ntotal, 'numwide*4=%s ntotal=%s' % (self.num_wide*4, ntotal)
        assert n > 0, n
        return n

    def _oef_crod(self, data, ndata, dt, is_magnitude_phase, prefix, postfix):
        """
        1-CROD
        3-CTUBE
        10-CONROD
        """
        n = 0
        obj_real = RealRodForceArray
        obj_complex = ComplexRodForceArray
        if self.element_type == 1: # CROD
            result_name = prefix + 'crod_force' + postfix
        elif self.element_type == 3:  # CTUBE
            result_name = prefix + 'ctube_force' + postfix
        elif self.element_type == 10:  # CONROD
            result_name = prefix + 'conrod_force' + postfix
        else:
            raise NotImplementedError(self.code_information())
            #msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
            #return self._not_implemented_or_skip(data, ndata, msg)

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)

        slot = self.get_result(result_name)
        if self.format_code == 1 and self.num_wide == 3: # real
            ntotal = 3 * 4 * self.factor
            nelements = ndata // ntotal

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_real)
            if auto_return:
                return nelements * self.num_wide * 4, None, None
            obj = self.obj
            if self.is_debug_file:
                self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                self.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
                self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 3)
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 3)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                #[axial, torsion]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                fmt = mapfmt(self._endian + self._analysis_code_fmt + b'ff', self.size)  # 3
                s = Struct(fmt)
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]
                    out = s.unpack(edata)
                    (eid_device, axial, torque) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    if self.is_debug_file:
                        self.binary_debug.write('OEF_Rod - %s\n' % (str(out)))
                    obj.add_sort1(dt, eid, axial, torque)
                    n += ntotal

        elif self.format_code in [2, 3] and self.num_wide == 5: # imag
            ntotal = 20 * self.factor  # 5*4
            nelements = ndata // ntotal

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_complex)
            if auto_return:
                return nelements * ntotal, None, None

            obj = self.obj
            if self.is_debug_file:
                self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                self.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
                self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, 5).copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype8).reshape(nelements, 5)
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
                fmt = mapfmt(self._endian + self._analysis_code_fmt + b'4f', self.size)
                s = Struct(fmt)
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]

                    out = s.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('OEF_CRod - %s\n' % (str(out)))
                    (eid_device, axial_real, torque_real, axial_imag, torque_imag) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    if is_magnitude_phase:
                        axial = polar_to_real_imag(axial_real, axial_imag)
                        torque = polar_to_real_imag(torque_real, torque_imag)
                    else:
                        axial = complex(axial_real, axial_imag)
                        torque = complex(torque_real, torque_imag)

                    obj.add_sort1(dt, eid, axial, torque)
                    n += ntotal
        else:  # pragma: no cover
            msg = self.code_information()
            print(msg)
            return self._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oef_cbeam(self, data, ndata, dt, is_magnitude_phase, prefix, postfix):
        """2-CBEAM"""
        n = 0
        result_name = prefix + 'cbeam_force' + postfix
        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)

        if self.format_code == 1 and self.num_wide == 9:  # real centroid ???
            raise RuntimeError('is this used?')
            #auto_return, is_vectorized = self._create_oes_object4(
                #nelements, result_name, slot, RealCBeamForceArray)
            #if auto_return:
                #return nelements * self.num_wide * 4

            #obj = self.obj
            ##is_vectorized = False
            #if self.use_vector and is_vectorized and self.sort_method == 1:
                #n = nelements * 4 * self.num_wide
                #itotal = obj.itotal
                #itotal2 = obj.itotal + nelements
                #ielement2 = obj.ielement + nelements

                #floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 9)[:, 1:]
                #obj._times[obj.itime] = dt
                #if obj.itime == 0:
                    #ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 9)
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
                #s = Struct(self._endian + b'i8f')  # 36
                #ntotal = 36
                #nelements = ndata // ntotal
                #obj = self.obj
                #for i in range(nelements):
                    #edata = data[n:n+36]
                    #out = s.unpack(edata)
                    #if self.is_debug_file:
                        #self.binary_debug.write('OEF_Beam - %s\n' % (str(out)))
                    #(eid_device, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq) = out
                    #eid, dt = get_eid_dt_from_eid_device(
                        #eid_device, self.nonlinear_factor, self.sort_method)
                    #n += 36

        elif self.format_code in [1, 2] and self.num_wide == 100:  # real/random
            # real - format_code == 1
            # random - format_code == 2
            #result_name, is_random = self._apply_oef_ato_crm_psd_rms_no(result_name)
            slot = self.get_result(result_name)
            ntotal = 400  # 1+(10-1)*11=100 ->100*4 = 400
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, RealCBeamForceArray)
            if auto_return:
                self._data_factor = 11
                return nelements * self.num_wide * 4, None, None
            obj = self.obj

            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * 4 * self.num_wide
                itotal = obj.itotal
                itotal2 = obj.itotal + nelements * 11
                #ielement = obj.ielement
                ielement2 = obj.ielement + nelements

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 100)[:, 1:]
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 100)
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
                n = oef_cbeam_real_100(self, data, obj,
                                       nelements, ntotal, dt)
        elif self.format_code in [2, 3] and self.num_wide == 177: # imag
            slot = self.get_result(result_name)
            ntotal = 708 * self.factor # 177*4
            nelements = ndata // ntotal

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, ComplexCBeamForceArray)
            if auto_return:
                self._data_factor = 11
                return nelements * self.num_wide * 4, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.itotal
                itotal2 = obj.itotal + nelements * 11
                #ielement = obj.ielement
                ielement2 = obj.ielement + nelements

                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, 177)[:, 1:]
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype8).reshape(nelements, 177)
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
                if self.is_debug_file:
                    self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    #self.binary_debug.write('  #elementi = [eid_device, force]\n')
                    #self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                ntotal = 708 * self.factor # (16*11+1)*4 = 177*4
                nelements = ndata // ntotal
                n = oef_cbeam_imag_177(self, data, obj, nelements, ntotal, is_magnitude_phase)
        else:
            msg = self.code_information()
            print(msg)
            return self._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oef_celas_cdamp(self, data, ndata, dt, is_magnitude_phase, prefix, postfix):
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
        n = 0
        if self.element_type == 11:
            result_name = prefix + 'celas1_force' + postfix
            obj_real = RealSpringForceArray
            obj_complex = ComplexSpringForceArray
        elif self.element_type == 12:
            result_name = prefix + 'celas2_force' + postfix
            obj_real = RealSpringForceArray
            obj_complex = ComplexSpringForceArray
        elif self.element_type == 13:
            result_name = prefix + 'celas3_force' + postfix
            obj_real = RealSpringForceArray
            obj_complex = ComplexSpringForceArray
        elif self.element_type == 14:
            result_name = prefix + 'celas4_force' + postfix
            obj_real = RealSpringForceArray
            obj_complex = ComplexSpringForceArray

        elif self.element_type == 20:
            result_name = prefix + 'cdamp1_force' + postfix
            obj_real = RealDamperForceArray
            obj_complex = ComplexDamperForceArray
        elif self.element_type == 21:
            result_name = prefix + 'cdamp2_force' + postfix
            obj_real = RealDamperForceArray
            obj_complex = ComplexDamperForceArray
        elif self.element_type == 22:
            result_name = prefix + 'cdamp3_force' + postfix
            obj_real = RealDamperForceArray
            obj_complex = ComplexDamperForceArray
        elif self.element_type == 23:
            result_name = prefix + 'cdamp4_force' + postfix
            obj_real = RealDamperForceArray
            obj_complex = ComplexDamperForceArray
        else:
            raise NotImplementedError(self.code_information())

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)
        if self.format_code == 1 and self.num_wide == 2:  # real
            ntotal = 8 * self.factor # 2 * 4
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_real)
            if auto_return:
                return nelements * self.num_wide * 4, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 2)
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 2)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                #(eid_device, force)
                obj.data[obj.itime, itotal:itotal2, 0] = floats[:, 1].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                fmt = mapfmt(self._endian + self._analysis_code_fmt + b'f', self.size)
                s = Struct(fmt)  # 2
                for unused_i in range(nelements):
                    edata = data[n:n + ntotal]
                    out = s.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('OEF_SpringDamper - %s\n' % str(out))
                    (eid_device, force) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    obj.add_sort1(dt, eid, force)
                    n += ntotal
        elif self.format_code in [2, 3] and self.num_wide == 3:  # imag
            ntotal = 12 * self.factor  # 3*4
            nelements = ndata // ntotal

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_complex)
            if auto_return:
                return nelements * ntotal, None, None

            obj = self.obj
            if self.is_debug_file:
                self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                self.binary_debug.write('  #elementi = [eid_device, force]\n')
                self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, 3).copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype8).reshape(nelements, 3)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                #[spring_force]
                real_imag = apply_mag_phase(floats, is_magnitude_phase, 1, 2)
                obj.data[obj.itime, itotal:itotal2, 0] = real_imag
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                fmt = mapfmt(self._endian + self._analysis_code_fmt + b'2f', self.size)
                s = Struct(fmt)
                for unused_i in range(nelements):
                    edata = data[n:n + ntotal]
                    out = s.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('OEF_SpringDamper - %s\n' % str(out))
                    (eid_device, force_real, force_imag) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    if is_magnitude_phase:
                        force = polar_to_real_imag(force_real, force_imag)
                    else:
                        force = complex(force_real, force_imag)
                    obj.add_sort1(dt, eid, force)
                    n += ntotal
        else:
            #msg = 'OEF: element_name=%s element_type=%s' % (self.element_name, self.element_type)
            msg = self.code_information()
            print(msg)
            return self._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oef_cvisc(self, data, ndata, dt, is_magnitude_phase, prefix, postfix):
        result_name = prefix + 'cvisc_force' + postfix
        if self._results.is_not_saved(result_name):
            return ndata, None, None
        n = 0
        self._results._found_result(result_name)

        slot = self.get_result(result_name)
        obj_real = RealViscForceArray

        if self.format_code == 1 and self.num_wide == 3: # real
            ntotal = 12 * self.factor # 3 * 4
            nelements = ndata // ntotal

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_real)
            if auto_return:
                return nelements * ntotal, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 3)
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 3)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                #(eid_device, axial, torque)
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                fmt = mapfmt(self._endian + self._analysis_code_fmt + b'ff', self.size)
                s = Struct(fmt)
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]

                    out = s.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('OEF_CVisc - %s\n' % (str(out)))
                    (eid_device, axial, torque) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    obj.add_sort1(dt, eid, axial, torque)
                    n += ntotal
        elif self.format_code in [2, 3] and self.num_wide == 5: # complex
            ntotal = 20 * self.factor # 5*4
            nelements = ndata // ntotal

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, ComplexViscForceArray)
            if auto_return:
                return nelements * ntotal, None, None

            obj = self.obj
            if self.is_debug_file:
                self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                self.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
                self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 5).copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 5)
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
                s = Struct(self._endian + self._analysis_code_fmt + b'4f')  # 5
                for unused_i in range(nelements):
                    edata = data[n:n+20]

                    out = s.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('OEF_CVisc - %s\n' % (str(out)))
                    (eid_device, axial_real, torque_real, axial_imag, torque_imag) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    if is_magnitude_phase:
                        axial = polar_to_real_imag(axial_real, axial_imag)
                        torque = polar_to_real_imag(torque_real, torque_imag)
                    else:
                        axial = complex(axial_real, axial_imag)
                        torque = complex(torque_real, torque_imag)

                    obj.add_sort1(dt, eid, axial, torque)
                    n += ntotal
        else:
            msg = self.code_information()
            print(msg)
            return self._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oef_cbar_34(self, data, ndata, dt, is_magnitude_phase, prefix, postfix):
        """
        34-CBAR
        117-CWELDC
        119-CFAST-NX  (126-CFAST-MSC is like the CBUSH)

        """
        n = 0
        if self.element_type == 34:
            result_name = prefix + 'cbar_force' + postfix
            obj_real = RealCBarForceArray
            obj_complex = ComplexCBarForceArray
        elif self.element_type in [117, 200]:
            result_name = prefix + 'cweld_force' + postfix
            obj_real = RealCWeldForceArray
            obj_complex = ComplexCWeldForceArray
            assert self.num_wide in [9, 17], self.code_information()
        elif self.element_type == 119:
            result_name = prefix + 'cfast_force' + postfix
            obj_real = RealCFastForceArrayNX
            assert self.num_wide == 9, self.code_information()
        else:
            raise NotImplementedError(self.element_type)

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)

        if self.format_code in [1, 2] and self.num_wide == 9:
            # real - format_code == 1
            # random - format_code == 2
            #result_name, is_random = self._apply_oef_ato_crm_psd_rms_no(result_name)
            slot = self.get_result(result_name)

            ntotal = 36 * self.factor  # 9*4
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_real)
            if auto_return:
                return nelements * ntotal, None, None

            obj = self.obj
            #return nelements * self.num_wide * 4
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, 9)
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype8).reshape(nelements, 9)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                #[bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            # elif self.use_vector and is_vectorized and self.sort_method == 1:
            else:
                n = oef_cbar_real(self, data, obj, nelements, ntotal)
        elif self.format_code in [2, 3] and self.num_wide == 17: # imag
            # TODO: vectorize
            ntotal = 68 * self.factor  # 17*4
            nelements = ndata // ntotal
            assert ndata % ntotal == 0

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_complex)
            if auto_return:
                return nelements * ntotal, None, None

            obj = self.obj
            n = oef_cbar_imag(self, data, obj, nelements, ntotal, is_magnitude_phase)
        else:
            msg = self.code_information()
            print(msg)
            return self._not_implemented_or_skip(data, ndata, msg)
        #print self.barForces
        return n, nelements, ntotal

    def _oef_cbar_100(self, data, ndata, dt, unused_is_magnitude_phase,
                      prefix, postfix):
        n = 0
        #100-BARS
        result_name = prefix + 'cbar_force' + postfix  # _10nodes
        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)

        if self.format_code == 1 and self.num_wide == 8:  # real
            ntotal = 32  # 8*4
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, RealCBar100ForceArray)
            if auto_return:
                return nelements * self.num_wide * 4, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * 4 * self.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 8)
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 8)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                #[axial, torsion, SMa, SMt]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                s = Struct(self._endian + self._analysis_code_fmt + b'7f')
                for unused_i in range(nelements):
                    edata = data[n:n+32]

                    out = s.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('OEF_CBar100 - %s\n' % (str(out)))
                    (eid_device, sd, bm1, bm2, ts1, ts2, af, trq) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    obj.add_sort1(dt, eid, sd, bm1, bm2, ts1, ts2, af, trq)
                    n += 32
        #elif self.format_code in [2, 3] and self.num_wide == 14:  # imag
        else:  # pragma: no cover
            msg = self.code_information()
            print(msg)
            return self._not_implemented_or_skip(data, ndata, msg)
        return n, nelements, ntotal

    def _oef_shells_centroidal(self, data, ndata, dt, is_magnitude_phase, prefix, postfix):
        """
        33-CQUAD4
        74-CTRIA3
        227-CTRIAR
        228-CQUADR

        """
        n = 0
        if self.element_type == 33:
            result_name = prefix + 'cquad4_force' + postfix
        elif self.element_type in [74, 83]:
            result_name = prefix + 'ctria3_force' + postfix
        elif self.element_type == 227:
            result_name = prefix + 'ctriar_force' + postfix
        elif self.element_type == 228:
            result_name = prefix + 'cquadr_force' + postfix
        else:
            #msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
            #return self._not_implemented_or_skip(data, ndata, msg)
            raise NotImplementedError(self.code_information())
        #result_name, is_random = self._apply_oef_ato_crm_psd_rms_no(result_name)

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)

        assert self._data_factor == 1, self._data_factor
        if self.format_code in [1, 2] and self.num_wide == 9:
            # real - format_code == 1
            # random - format_code == 2
            ntotal = 36 * self.factor # 9*4
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, RealPlateForceArray)
            if auto_return:
                return nelements * ntotal, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * ntotal
                ielement = obj.ielement
                ielement2 = ielement + nelements

                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, 9)
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype8).reshape(nelements, 9)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[ielement:ielement2] = eids

                #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
                obj.data[obj.itime, ielement:ielement2, :] = floats[:, 1:].copy()
                obj.itotal = ielement2
                obj.ielement = ielement2
            else:
                s = Struct(mapfmt(self._endian + self._analysis_code_fmt + b'8f', self.size))
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]
                    out = s.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('real_OEF_Plate-%s - %s\n' % (self.element_type, str(out)))
                    (eid_device, mx, my, mxy, bmx, bmy, bmxy, tx, ty) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    obj.add_sort1(dt, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)
                    n += ntotal
        elif self.format_code in [2, 3] and self.num_wide == 17:  # imag
            ntotal = 68 * self.factor
            nelements = ndata // ntotal

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, ComplexPlateForceArray)
            if auto_return:
                return nelements * ntotal, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * ntotal
                ielement = obj.ielement
                ielement2 = ielement + nelements
                itotal = obj.itotal
                itotal2 = itotal + nelements

                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, 17).copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype8).reshape(nelements, 17)
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
                s = Struct(mapfmt(self._endian + self._analysis_code_fmt + b'16f', self.size))
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]
                    out = s.unpack(edata)
                    (eid_device,
                     mxr, myr, mxyr, bmxr, bmyr, bmxyr, txr, tyr,
                     mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    if self.is_debug_file:
                        self.binary_debug.write('complex_OEF_Plate-%s - %s\n' % (self.element_type, str(out)))

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
                    obj.add_sort1(dt, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)
                    n += ntotal
        else:  # pragma: no cover
            msg = self.code_information()
            print(msg)
            return self._not_implemented_or_skip(data, ndata, msg)
        return n, nelements, ntotal

    def _oef_shells_nodal(self, data, ndata, dt, is_magnitude_phase, prefix, postfix):
        """
        64-CQUAD8
        70-CTRIAR
        75-CTRIA6
        82-CQUADR
        144-CQUAD4-bilinear

        """
        n = 0
        if self.element_type == 64:
            result_name = prefix + 'cquad8_force' + postfix
        elif self.element_type == 70:
            result_name = prefix + 'ctriar_force' + postfix
        elif self.element_type == 75:
            result_name = prefix + 'ctria6_force' + postfix
        elif self.element_type == 82:
            result_name = prefix + 'cquadr_force' + postfix
        elif self.element_type == 144:
            result_name = prefix + 'cquad4_force' + postfix
        else:
            raise NotImplementedError(self.code_information())
            #msg = self.code_information()
            #return self._not_implemented_or_skip(data, ndata, msg)

        if self.element_type in [70, 75]:  # CTRIAR,CTRIA6
            nnodes = 3
        elif self.element_type in [64, 82, 144]:  # CQUAD8,CQUADR,CQUAD4-bilinear
            nnodes = 4
        else:
            raise NotImplementedError(self.code_information())
            #msg = 'name=%r type=%r' % (self.element_name, self.element_type)
            #return self._not_implemented_or_skip(data, ndata, msg), None, None

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)

        slot = self.get_result(result_name)
        nnodes_all = nnodes + 1
        numwide_real = 2 + nnodes_all * 9 # centroidal node is the + 1
        numwide_imag = 2 + nnodes_all * 17

        if self.format_code == 1 and self.num_wide == numwide_real:  # real
            obj_real = RealPlateBilinearForceArray

            ntotal = (8 + nnodes_all * 36) * self.factor # centroidal node is the + 1
            assert ntotal == self.num_wide * 4 * self.factor, 'ntotal=%s numwide=%s' % (ntotal, self.num_wide * 4)

            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_real)
            if auto_return:
                self._data_factor = nnodes_all
                return nelements * ntotal, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                nlayers = nelements * nnodes_all
                n = nelements * self.num_wide * 4

                istart = obj.itotal
                iend = istart + nlayers
                obj._times[obj.itime] = dt

                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype8).reshape(nelements, numwide_real).copy()
                    # Nastran makes this a 4 for CQUAD4s instead
                    # of 0 like the bilinear stress element...
                    ints[:, 2] = 0

                    nids = ints[:, 2:].reshape(nlayers, 9)[:, 0]
                    eids = ints[:, 0] // 10
                    eids2 = vstack([eids] * nnodes_all).T.ravel()
                    obj.element_node[istart:iend, 0] = eids2
                    obj.element_node[istart:iend, 1] = nids

                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, numwide_real)
                results = floats[:, 2:].reshape(nlayers, 9)[:, 1:].copy()
                #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
                obj.data[obj.itime, istart:iend, :] = results
            else:

                n44 = 44 * self.factor
                n36 = 36 * self.factor
                if self.size == 4:
                    fmt1 = self._endian + b'i4si8f'  # 8+36
                    fmt2 = self._endian + b'i8f' # 36
                else:
                    fmt1 = self._endian + b'q8sq8d'
                    fmt2 = self._endian + b'q8d'
                s1 = Struct(fmt1)
                s2 = Struct(fmt2)

                for unused_i in range(nelements):
                    edata = data[n:n + n44]

                    out = s1.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('OEF_Plate2-%s - %s\n' % (self.element_type, str(out)))
                    (eid_device, term, _nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty) = out
                    #term= 'CEN\'
                    #_nid = 4
                    # -> CEN/4
                    nid = 0
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    obj.add_sort1(dt, eid, term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)
                    n += n44
                    for unused_j in range(nnodes):
                        edata = data[n : n + n36]
                        out = s2.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('    %s\n' % (str(out)))
                        (nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty) = out
                        assert nid > 0, 'nid=%s' % nid
                        obj.add_sort1(dt, eid, term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)
                        n += n36
        elif self.format_code in [2, 3] and self.num_wide == numwide_imag: # complex
            ntotal = numwide_imag * 4 * self.factor
            nelements = ndata // ntotal

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, ComplexPlate2ForceArray)
            if auto_return:
                self._data_factor = nnodes_all
                return nelements * ntotal, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.itotal
                ielement = obj.ielement
                ielement2 = obj.ielement + nelements
                itotal2 = obj.itotal + nelements * nnodes_all

                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, numwide_imag)
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype8).reshape(nelements, numwide_imag).copy()
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
                if self.size == 4:
                    s1 = Struct(self._endian + b'i4s17f')  # 2+17=19 * 4 = 76
                    s2 = Struct(self._endian + b'i16f')  # 17 * 4 = 68
                else:
                    s1 = Struct(self._endian + b'q8s17d')  # 2+17=19 * 4 = 768
                    s2 = Struct(self._endian + b'q16d')  # 17 * 4 = 68
                ntotal = (8 + (nnodes + 1) * 68) * self.factor
                ntotal1 = 76 * self.factor
                ntotal2 = 68 * self.factor

                nelements = ndata // ntotal
                obj = self.obj
                for unused_i in range(nelements):
                    edata = data[n:n + ntotal1]
                    n += ntotal1

                    out = s1.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('OEF_Plate2-%s - %s\n' % (self.element_type, str(out)))
                    (eid_device, term, nid,
                     mxr, myr, mxyr, bmxr, bmyr, bmxyr, txr, tyr,
                     mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi) = out
                    #term = 'CEN\'

                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
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
                        if self.is_debug_file:
                            self.binary_debug.write('OEF_Plate2 - eid=%i nid=%s out=%s\n' % (
                                eid, nid, str(out)))
                        obj.add_sort1(dt, eid, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)
        else:  # pragma: no cover
            msg = self.code_information()
            print(msg)
            return self._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oef_shells_composite(self, data, ndata, dt, unused_is_magnitude_phase,
                              prefix, postfix):
        """
        95 - CQUAD4
        96 - CQUAD8
        97 - CTRIA3
        98 - CTRIA6 (composite)
        232 - CQUADR
        233 - CTRIAR

        """
        if self.element_type == 95:
            result_name = prefix + 'cquad4_composite_force' + postfix
        elif self.element_type == 96:
            result_name = prefix + 'cquad8_composite_force' + postfix
        elif self.element_type == 97:
            result_name = prefix + 'ctria3_composite_force' + postfix
        elif self.element_type == 98:
            result_name = prefix + 'ctria6_composite_force' + postfix
        elif self.element_type == 232:
            result_name = prefix + 'cquadr_composite_force' + postfix
        elif self.element_type == 233:
            result_name = prefix + 'ctriar_composite_force' + postfix
        else:  # pragma: no cover
            raise NotImplementedError(self.code_information())
        if self._results.is_not_saved(result_name):
            return ndata, None, None

        self._results._found_result(result_name)
        slot = self.get_result(result_name)

        n = 0
        if self.format_code == 1 and self.num_wide == 9:  # real
            ntotal = 36 # 9 * 4
            nelements = ndata // ntotal


            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, FailureIndicesArray)
            #print('read_mode ', self.read_mode, auto_return, is_vectorized)
            if auto_return:
                #self._data_factor = nnodes_all
                return nelements * self.num_wide * 4, None, None

            obj = self.obj
            nelements = ndata // ntotal


            ## TODO: add
            #if self.is_debug_file:
                #self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                ##self.binary_debug.write('  #centeri = [eid_device, j, grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
                ##self.binary_debug.write('  #                                fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,)]\n')
                ##self.binary_debug.write('  #nodeji = [eid, ilayer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)]\n')
                #self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            n = oef_shells_composite_real_9(self, data, obj, nelements, ntotal, dt)
        else:  # pragma: no cover
            msg = self.code_information()
            print(msg)
            return self._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oef_cshear(self, data, ndata, dt, is_magnitude_phase, prefix, postfix):
        """4-CSHEAR"""
        result_name = prefix + 'cshear_force' + postfix
        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)

        n = 0
        if self.format_code == 1 and self.num_wide == 17:  # real
            ntotal = 68  # 17*4
            nelements = ndata // ntotal

            obj_real = RealCShearForceArray
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_real)
            if auto_return:
                return nelements * self.num_wide * 4, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * 4 * self.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 17)
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 17)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                # [f41, f21, f12, f32, f23, f43, f34, f14, kf1,
                #  s12, kf2, s23, kf3, s34, kf4, s41]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                s = Struct(self._endian + self._analysis_code_fmt + b'16f')
                for unused_i in range(nelements):
                    edata = data[n:n+68]

                    out = s.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('OEF_Shear - %s\n' % (str(out)))
                    (eid_device,
                     f41, f21, f12, f32, f23, f43, f34, f14, kf1,
                     s12, kf2, s23, kf3, s34, kf4, s41) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)

                    #data_in = [eid,
                               #f41, f21, f12, f32, f23, f43, f34,
                               #f14, kf1, s12, kf2, s23, kf3, s34, kf4, s41]
                    #print "%s" % (self.get_element_type(self.element_type)), data_in
                    obj.add_sort1(dt, eid,
                                  f41, f21, f12, f32, f23, f43, f34,
                                  f14, kf1, s12, kf2, s23, kf3, s34, kf4, s41)
                    n += ntotal

        elif self.format_code in [2, 3] and self.num_wide == 33:  # imag
            ntotal = 132 * self.factor # 33*4
            nelements = ndata // ntotal

            obj_complex = ComplexCShearForceArray
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_complex)
            if auto_return:
                return nelements * ntotal, None, None

            obj = self.obj
            #if self.is_debug_file:
                #self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                #self.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
                #self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * 4 * self.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, 33).copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype8).reshape(nelements, 33)
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
                #self.create_transient_object(self.cshear_force, ComplexCShearForce)
                s = Struct(mapfmt(self._endian + self._analysis_code_fmt + b'32f', self.size))
                #ntotal1 = 132 * self.factor
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]
                    n += ntotal
                    out = s.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('OEF_Shear - %s\n' % (str(out)))
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
                        eid_device, self.nonlinear_factor, self.sort_method)
                    obj.add_sort1(dt, eid,
                                  f41, f21, f12, f32, f23, f43, f34, f14,
                                  kf1, s12, kf2, s23, kf3, s34, kf4, s41)
        else:  # pragma: no cover
            msg = self.code_information()
            return self._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oef_cconeax(self, data, ndata, dt, unused_is_magnitude_phase,
                     prefix, postfix):
        """35-CONEAX"""
        result_name = prefix + 'cconeax_force' + postfix
        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)

        n = 0
        if self.format_code == 1 and self.num_wide == 7:  # real
            ntotal = 28  # 7*4
            nelements = ndata // ntotal

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, RealConeAxForceArray)
            if auto_return:
                return nelements * self.num_wide * 4, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * 4 * self.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 7)
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 7)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                # [hopa, bmu, bmv, tm, su, sv]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                s = Struct(self._endian + self._analysis_code_fmt + b'6f')
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]
                    out = s.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('OEF_CONEAX-35 - %s\n' % (str(out)))
                    (eid_device, hopa, bmu, bmv, tm, su, sv) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    obj.add_sort1(dt, eid, hopa, bmu, bmv, tm, su, sv)
                    n += ntotal
        else:  # pragma: no cover
            msg = self.code_information()
            return self._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oef_cgap(self, data, ndata, dt, unused_is_magnitude_phase,
                  prefix, postfix):
        """38-GAP"""
        result_name = prefix + 'cgap_force' + postfix
        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)

        n = 0
        if self.format_code == 1 and self.num_wide == 9:  # real
            ntotal = 36 # 9*4
            nelements = ndata // ntotal
            obj = self.obj

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, RealCGapForceArray)
            if auto_return:
                return nelements * self.num_wide * 4, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * 4 * self.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 9)
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 9)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                # [fx, sfy, sfz, u, v, w, sv, sw]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                s = Struct(self._endian + self._analysis_code_fmt + b'8f')
                for unused_i in range(nelements):
                    edata = data[n:n+36]

                    out = s.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('OEF_CGAP-38 - %s\n' % (str(out)))
                    (eid_device, fx, sfy, sfz, u, v, w, sv, sw) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    #data_in = [eid, fx, sfy, sfz, u, v, w, sv, sw]
                    #print "%s" %(self.get_element_type(self.element_type)), data_in
                    #eid = obj.add_new_eid_sort1(out)
                    obj.add_sort1(dt, eid, fx, sfy, sfz, u, v, w, sv, sw)
                    n += ntotal
        else:  # pragma: no cover
            msg = self.code_information()
            print(msg)
            return self._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oef_cbend(self, data, ndata, dt, is_magnitude_phase, prefix, postfix):
        """69-CBEND"""
        result_name = prefix + 'cbend_force' + postfix
        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)

        n = 0
        if self.format_code == 1 and self.num_wide == 15:  # real
            ntotal = 60  # 15*4
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, RealBendForceArray)
            if auto_return:
                return nelements * self.num_wide * 4, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * 4 * self.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 15).copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 15).copy()
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
                s = Struct(self._endian + self._analysis_code_fmt + b' i6fi6f')
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]

                    out = s.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('OEF_BEND-69 - %s\n' % (str(out)))
                    (eid_device,
                     nid_a, bm1_a, bm2_a, ts1_a, ts2_a, af_a, trq_a,
                     nid_b, bm1_b, bm2_b, ts1_b, ts2_b, af_b, trq_b) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)

                    obj.add_sort1(
                        dt, eid,
                        nid_a, bm1_a, bm2_a, ts1_a, ts2_a, af_a, trq_a,
                        nid_b, bm1_b, bm2_b, ts1_b, ts2_b, af_b, trq_b)
                    n += ntotal
        elif self.format_code in [2, 3] and self.num_wide == 27:  # imag
            # TODO: vectorize
            ntotal = 108  # 27*4
            nelements = ndata // ntotal

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, ComplexCBendForceArray)
            if auto_return:
                return nelements * self.num_wide * 4, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * 4 * self.num_wide
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
                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 27).copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 27).copy()
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
                s = Struct(self._endian + self._analysis_code_fmt + b' i12f i12f')
                for unused_i in range(nelements):
                    edata = data[n:n+108]
                    n += ntotal
                    out = s.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('OEF_BEND-69 - %s\n' % (str(out)))
                    (eid_device, nid_a,
                     bm1_ar, bm2_ar, ts1_ar, ts2_ar, af_ar, trq_ar,
                     bm1_ai, bm2_ai, ts1_ai, ts2_ai, af_ai, trq_ai,
                     nid_b,
                     bm1_br, bm2_br, ts1_br, ts2_br, af_br, trq_br,
                     bm1_bi, bm2_bi, ts1_bi, ts2_bi, af_bi, trq_bi) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)

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

                    obj.add_sort1(dt, eid,
                                  nid_a, bm1_a, bm2_a, ts1_a, ts2_a, af_a, trq_a,
                                  nid_b, bm1_b, bm2_b, ts1_b, ts2_b, af_b, trq_b)
        else:  # pragma: no cover
            msg = self.code_information()
            return self._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oef_csolid_pressure(self, data, ndata, dt, is_magnitude_phase, prefix, postfix):
        """
        76-HEXPR
        77-PENPR
        78-TETPR
        """
        n = 0
        if self.element_type == 76:
            result_name = prefix + 'chexa_pressure_force' + postfix
        elif self.element_type == 77:
            result_name = prefix + 'cpenta_pressure_force' + postfix
        elif self.element_type == 78:
            result_name = prefix + 'ctetra_pressure_force' + postfix
        elif self.element_type == 79:
            result_name = prefix + 'cpyram_pressure_force' + postfix
        else:
            msg = self.code_information()
            return self._not_implemented_or_skip(data, ndata, msg), None, None
        slot = self.get_result(result_name)

        self._results._found_result(result_name)
        if self.format_code == 1 and self.num_wide == 10:  # real
            ntotal = 40 * self.factor
            nelements = ndata // ntotal
            #nelements = ndata // ntotal

            obj_real = RealSolidPressureForceArray
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_real)
            if auto_return:
                return nelements * ntotal, None, None

            obj = self.obj
            #if self.is_debug_file:
                #self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                #self.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
                #self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)
            if self.use_vector and is_vectorized and self.sort_method == 1:
                # self.itime = 0
                # self.ielement = 0
                # self.itotal = 0
                #self.ntimes = 0
                #self.nelements = 0
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 10)
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 10)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                #[axial_force, torque]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 3:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if self.size == 4:
                    fmt = self._endian + self._analysis_code_fmt + b'8s7f'
                else:
                    fmt = self._endian + mapfmt(self._analysis_code_fmt, self.size) + b'16s7d'
                s = Struct(fmt)
                for unused_i in range(nelements):
                    edata = data[n : n + ntotal]
                    n += ntotal
                    out = s.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('OEF_PentaPressure-%s %s\n' % (self.element_type, str(out)))
                    (eid_device, ename, ax, ay, az, vx, vy, vz, pressure) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    obj.add_sort1(dt, eid, ename, ax, ay, az, vx, vy, vz, pressure)

        elif self.format_code in [2, 3] and self.num_wide == 16:  # imag
            ntotal = 64 * self.factor
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, ComplexSolidPressureForceArray)
            if auto_return:
                return nelements * ntotal, None, None

            obj = self.obj
            #if self.is_debug_file:
                #self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                #self.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
                #self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, 16).copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype8).reshape(nelements, 16)
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
                if self.size == 4:
                    s = Struct(self._endian + self._analysis_code_fmt + b'8s 13f')
                else:
                    s = Struct(mapfmt(self._endian + self._analysis_code_fmt, self.size) + b'16s 13d')
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]
                    n += ntotal

                    #print(len(edata))
                    out = s.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('_oef_csolid_pressure-%s %s\n' % (self.element_type, str(out)))
                    (eid_device, ename,
                     axr, ayr, azr, vxr, vyr, vzr, pressure,
                     axi, ayi, azi, vxi, vyi, vzi) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
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
                    obj.add_sort1(dt, eid, ename, ax, ay, az, vx, vy, vz, cpressure)
        else:  # pragma: no cover
            msg = self.code_information()
            return self._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oef_cbush(self, data, ndata, dt, is_magnitude_phase, prefix, postfix):
        """
        102-CBUSH
        126-CFAST-MSC
        280-CBEAR
        """
        if self.element_type == 102:
            result_name = prefix + 'cbush_force' + postfix
            real_obj = RealCBushForceArray
            complex_obj = ComplexCBushForceArray
        elif self.element_type == 126:
            result_name = prefix + 'cfast_force' + postfix
            real_obj = RealCFastForceArrayMSC
            complex_obj = None
        elif self.element_type == 280:
            result_name = prefix + 'cbear_force' + postfix
            assert self.num_wide in [7, 13], self.code_information()
            real_obj = RealCBearForceArray
            complex_obj = ComplexCBearForceArray
        else:
            raise NotImplementedError(self.code_information())

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        #result_name, is_random = self._apply_oef_ato_crm_psd_rms_no(result_name)
        self._results._found_result(result_name)
        slot = self.get_result(result_name)

        n = 0
        numwide_real = 7
        if self.format_code in [1, 2] and self.num_wide == 7:  # real/random
            # real - format_code == 1
            # random - format_code == 2
            ntotal = 28 # 7*4
            nelements = ndata // ntotal

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, real_obj)
            if auto_return:
                return nelements * self.num_wide * 4, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                # self.itime = 0
                # self.ielement = 0
                # self.itotal = 0
                #self.ntimes = 0
                #self.nelements = 0
                n = nelements * self.num_wide * 4

                istart = obj.itotal
                iend = istart + nelements
                obj._times[obj.itime] = dt

                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, numwide_real).copy()
                    eids = ints[:, 0] // 10
                    obj.element[istart:iend] = eids
                results = frombuffer(data, dtype=self.fdtype).reshape(nelements, numwide_real)

                #[fx, fy, fz, mx, my, mz]
                obj.data[obj.itime, istart:iend, :] = results[:, 1:].copy()
            else:
                s = Struct(self._endian + self._analysis_code_fmt + b'6f')
                for unused_i in range(nelements):
                    edata = data[n:n+28]
                    out = s.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('OEF_CBUSH-102 - %s\n' % (str(out)))
                    (eid_device, fx, fy, fz, mx, my, mz) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    obj.add_sort1(dt, eid, fx, fy, fz, mx, my, mz)
                    n += ntotal
        elif self.format_code in [2, 3] and self.num_wide == 13:  # imag
            # TODO: vectorize
            ntotal = 52  # 13*4
            nelements = ndata // ntotal
            #result_name = prefix + 'cbush_force' + postfix
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, complex_obj)
            if auto_return:
                return nelements * self.num_wide * 4, None, None

            s = Struct(self._endian + self._analysis_code_fmt + b'12f')

            obj = self.obj
            for unused_i in range(nelements):
                edata = data[n:n + 52]

                out = s.unpack(edata)
                if self.is_debug_file:
                    self.binary_debug.write('OEF_CBUSH-102 - %s\n' % (str(out)))
                (eid_device,
                 fxr, fyr, fzr, mxr, myr, mzr,
                 fxi, fyi, fzi, mxi, myi, mzi) = out
                eid, dt = get_eid_dt_from_eid_device(
                    eid_device, self.nonlinear_factor, self.sort_method)

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

                obj.add_sort1(dt, eid, fx, fy, fz, mx, my, mz)
                n += ntotal
        #elif self.format_code == 2 and self.num_wide == 7:
            #self.log.warning(self.code_information())
        else:  # pragma: no cover
            msg = self.code_information()
            print(msg)
            return self._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oef_vu_shell(self, data, ndata, dt, is_magnitude_phase, prefix, postfix):
        """
        189-VUQUAD
        190-VUTRIA

        """
        n = 0
        if self.element_type == 189:  # VUQUAD
            nnodes = 4
            etype = 'VUQUAD4'
            result_name = prefix + 'vu_quad_force' + postfix
        elif self.element_type == 190:  # VUTRIA
            nnodes = 3
            etype = 'VUTRIA3'
            result_name = prefix + 'vu_tria_force' + postfix
        else:
            raise NotImplementedError(self.code_information())

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)

        slot = self.get_result(result_name)
        numwide_real = 6 + 13 * nnodes
        numwide_imag = 6 + 25 * nnodes

        if self.format_code == 1 and self.num_wide == numwide_real:  # real
            # real - format_code == 1
            # random - format_code == 2

            #ntotal = (6 + nnodes * 13) * 4 # 6+n*13
            ntotal = (24 + 52 * nnodes) * self.factor
            nelements = ndata // ntotal
            nlayers = nelements * nnodes
            #result_name = 'force_VU_2D'
            self._results._found_result(result_name)

            auto_return, is_vectorized = self._create_oes_object4(
                nlayers, result_name, slot, RealForceVU2DArray)
            if auto_return:
                self._data_factor = nnodes  # number of "layers" for an element
                return nelements * ntotal, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and 0: ## TODO: vectorize
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
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, numwide_real)
                    eids = ints[:, 0] // 10
                    obj.element_node[istart:iend] = eids_nids
                results = frombuffer(data, dtype=self.fdtype).reshape(nelements, numwide_real)

                #[fx, fy, fz, mx, my, mz]
                obj.data[obj.itime, istart:iend, :] = results[:, 1:].copy()
            else:
                # 6+n*13
                if self.size == 4:
                    s1 = Struct(self._endian + b'3i4s2i') # 6
                else:
                    s1 = Struct(self._endian + b'3q8s2q') # 6
                s2 = Struct(mapfmt(self._endian + b'i3f3i5fi', self.size)) # 13
                ntotal1 = 24 * self.factor # 6*4
                ntotal2 = 52 * self.factor # 13*4
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal1]
                    n += ntotal1

                    out = s1.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('OEF_Force_%s-%s - %s\n' % (
                            etype, self.element_type, str(out)))
                    (eid_device, parent, coord, icord, theta, _) = out

                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)

                    #vugrids = []
                    #forces = []
                    for unused_j in range(nnodes):
                        edata = data[n:n+ntotal2]
                        n += ntotal2
                        out = s2.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('%s\n' % (str(out)))
                        (vugrid, mfx, mfy, mfxy, unused_ai, unused_bi, unused_ci, bmx, bmy,
                         bmxy, syz, szx, unused_di) = out
                        #out2 = (vugrid, mfx, mfy, mfxy, bmx, bmy, bmxy, syz, szx)
                        obj.add_sort1(dt, eid, parent, coord, icord, theta,
                                      vugrid, mfx, mfy, mfxy, bmx, bmy, bmxy, syz, szx)
                        #vugrids.append(vugrid)
                        #forces.append(out2)
                    #data_in = [vugrid,mfx,mfy,mfxy,a,b,c,bmx,bmy,bmxy,syz,szx,d]
                    #obj.add_sort1(dt, eid, parent, coord, icord, theta, vugrids, forces)

        elif self.format_code in [2, 3] and self.num_wide == numwide_imag:  # imag
            if self._results.is_not_saved(result_name):
                return ndata, None, None
            self._results._found_result(result_name)

            ntotal = numwide_imag #6 + 25 * nnodes; nnodes=[3, 4]
            nelements = ndata // ntotal
            assert ndata % ntotal == 0

            nlayers = nelements * nnodes
            auto_return, is_vectorized = self._create_oes_object4(
                nlayers, result_name, slot, ComplexForceVU_2DArray)
            if auto_return:
                self._data_factor = nnodes  # number of "layers" for an element
                return nelements * self.num_wide * 4, None, None

            obj = self.obj
            #print('dt=%s, itime=%s' % (obj.itime, dt))
            is_vectorized = False
            if self.use_vector and is_vectorized and self.sort_method == 1:
                raise NotImplementedError()
            else:
                ntotal = (24 + 100 * nnodes) * self.factor
                if self.size == 4:
                    s1 = Struct(self._endian + b'iii4sii')
                else:
                    s1 = Struct(self._endian + b'qqq8sqq')
                s2 = Struct(mapfmt(self._endian + b'i3f3i5fi3f3i5fi', self.size))
                nelements = ndata // ntotal
                ntotal1 = 24 * self.factor # 6*4
                ntotal2 = 100 * self.factor # 13*4
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal1]
                    n += ntotal1

                    out = s1.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('OEF_Force_%s-%s - %s\n' % (
                            etype, self.element_type, str(out)))
                    (eid_device, parent, coord, icord, theta, _) = out

                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)

                    vugrids = []
                    forces = []
                    for unused_j in range(nnodes):
                        edata = data[n:n+ntotal2]
                        n += ntotal2
                        out = s2.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('OEF_Force_%s-%s - %s\n' % (
                                etype, self.element_type, str(out)))
                        [vugrid, mfxr, mfyr, mfxyr, unused_ai1, unused_bi1, unused_ci1,
                         bmxr, bmyr, bmxyr, syzr, szxr, unused_di1,
                         mfxi, mfyi, mfxyi, unused_ai2, unused_bi2, unused_ci2,
                         bmxi, bmyi, bmxyi, syzi, szxi, unused_di2] = out

                        if is_magnitude_phase:
                            mfx = polar_to_real_imag(mfxr, mfxi)
                            mfy = polar_to_real_imag(mfyr, mfyi)
                            mfxy = polar_to_real_imag(mfxyr, mfxyi)
                            bmx = polar_to_real_imag(bmxr, bmxi)
                            bmy = polar_to_real_imag(bmyr, bmyi)
                            bmxy = polar_to_real_imag(bmxyr, bmxyi)
                            syz = polar_to_real_imag(syzr, syzi)
                            szx = polar_to_real_imag(szxr, szxi)
                        else:
                            mfx = complex(mfxr, mfxi)
                            mfy = complex(mfyr, mfyi)
                            mfxy = complex(mfxyr, mfxyi)
                            bmx = complex(bmxr, bmxi)
                            bmy = complex(bmyr, bmyi)
                            bmxy = complex(bmxyr, bmxyi)
                            syz = complex(syzr, syzi)
                            szx = complex(szxr, szxi)
                        vugrids.append(vugrid)
                        forcei = [mfx, mfy, mfxy, bmx, bmy, bmxy, syz, szx]
                        forces.append(forcei)

                    #data_in = [vugrid,mfxr,mfyr,mfxyr,bmxr,bmyr,bmxyr,syzr,szxr,
                                     #mfxi,mfyi,mfxyi,bmxi,bmyi,bmxyi,syzi,szxi]
                    obj.add_sort1(nnodes, dt, eid, parent, coord, icord, theta, vugrids, forces)
        else:  # pragma: no cover
            msg = self.code_information()
            return self._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oef_vu_beam(self, data, ndata, dt, is_magnitude_phase, prefix, postfix):
        """191-VUBEAM"""
        n = 0
        result_name = prefix + 'cbeam_force_vu' + postfix
        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)

        if self.format_code == 1 and self.num_wide == 20:  # real
            #ELTYPE = 191 Beam view element (VUBEAM)
            #---------------------------------------
            #2 PARENT I     Parent p-element identification number
            #3 COORD  I     Coordinate system identification number
            #4 ICORD  CHAR4 Flat/curved and so on

            #TCODE,7 = 0 Real
            #5 VUGRID   I  VU grid ID for output grid
            #6 POSIT    RS x/L position of VU grid identification number
            #7 FORCEX   RS Force x
            #8 SHEARY   RS Shear force y
            #9 SHEARZ   RS Shear force z
            #10 TORSION RS Torsional moment x
            #11 BENDY   RS Bending moment y
            #12 BENDZ   RS Bending moment z

            ntotal = self.num_wide * 4 # 80
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, RealCBeamForceVUArray)
            if auto_return:
                self._data_factor = 2
                return nelements * self.num_wide * 4, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                # self.itime = 0
                # self.ielement = 0
                # self.itotal = 0
                #self.ntimes = 0
                #self.nelements = 0
                n = nelements * self.num_wide * 4
                itotal = obj.itotal
                ielement = obj.ielement
                ielement2 = ielement + nelements
                itotal2 = itotal + nelements * 2

                # 20 values
                ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 20).copy()
                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 20)[:, 4:]
                ints2 = ints[:, 4:].reshape(nelements * 2, 8)
                assert floats.shape[1] == 16, floats.shape

                obj._times[obj.itime] = dt
                #if obj.itime == 0:
                eids = ints[:, 0] // 10
                parent = ints[:, 1]
                coord = ints[:, 2]
                #icord = ints[:, 3]

                ints2 = ints[:, 4:].reshape(nelements*2, 8)
                nids = ints2[:, 0]
                eids2 = np.repeat(eids, 2)

                #icord = ints[:, 3]  is this a string?
                assert eids.min() > 0, eids.min()
                obj.element_node[itotal:itotal2, 0] = eids2
                obj.element_node[itotal:itotal2, 1] = nids
                obj.parent_coord[ielement:ielement2, 1] = parent
                obj.parent_coord[ielement:ielement2, 1] = coord

                floats2 = floats.reshape(nelements*2, 8)
                #[xxb, fx. fy, fz, mx, my, mz]
                obj.data[obj.itime, itotal:itotal2, :] = floats2[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                s1 = Struct(self._endian + b'3i 4s')
                s2 = Struct(self._endian + b'i7f')
                nnodes = 2
                for unused_i in range(nelements):
                    edata = data[n:n+16]  # 8*4
                    n += 16

                    out = s1.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('OEF_Force_VU-191 - %s\n' % (str(out)))
                    (eid_device, parent, coord, icord) = out

                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    #data_in = [eid, parent, coord, icord]

                    #forces = []
                    for unused_j in range(nnodes):
                        edata = data[n:n+32]  # 8*4
                        n += 32
                        out = s2.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('%s\n' % str(out))
                        nid, xxb, fx, fy, fz, mx, my, mz = out
                        obj._add_sort1(dt, eid, parent, coord, icord, nid, xxb, fx, fy, fz, mx, my, mz)
                    #data_in = [vugrid, posit, forceX, shearY, shearZ, torsion, bendY, bendZ]
        elif self.format_code == 1 and self.num_wide == 32:  # random
            # TODO: vectorize
            return ndata, None, None
        elif self.format_code in [2, 3] and self.num_wide == 32:  # imag
            #TCODE,7 = 1 Real/imaginary or magnitude/phase
            #5 VUGRID   I  VU grid identification number for output grid
            #6 POSIT    RS x/L position of VU grid identification number

            #7 FORCEXR  RS Force x real/mag.
            #8 SHEARYR  RS Shear force y real/mag.
            #9 SHEARZR  RS Shear force z real/mag.
            #10 TORSINR RS Torsional moment x real/mag.
            #11 BENDYR  RS Bending moment y real/mag.
            #12 BENDZR  RS Bending moment z real/mag.

            #13 FORCEXI RS Force x imag./phase
            #14 SHEARYI RS Shear force y imag./phase
            #15 SHEARZI RS Shear force z imag./phase
            #16 TORSINI RS Torsional moment x imag./phase
            #17 BENDYI  RS Bending moment y imag./phase
            #18 BENDZI  RS Bending moment z imag./phase
            #Words 5 through max repeat 2 times

            # 32 = 4 + 56/4 * 2 = 4 + 14 * 2 = 4 + 28
            #nnodes = 2
            #ntotal = 16 + 56 * nnodes
            ntotal = self.num_wide * 4 * self.factor
            nelements = ndata // ntotal

            # TODO: vectorize
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, ComplexCBeamForceVUArray)
            if auto_return:
                self._data_factor = 2
                return ndata, None, None
            obj = self.obj

            #is_vectorized = False
            if self.use_vector and is_vectorized and self.sort_method == 1:
                #ntotal = nelements * 2
                #raise NotImplementedError('ComplexBeamForceVUArray')
                n = nelements * self.num_wide * 4 * self.factor

                ielement = obj.ielement
                ielement2 = ielement + nelements

                itotal = obj.itotal
                itotal2 = itotal + nelements*2
                obj._times[obj.itime] = dt

                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 32).copy()

                    eids = ints[:, 0] // 10
                    eids2 = np.repeat(eids, 2)
                    parent = ints[:, 1]
                    coord = ints[:, 2]
                    #parent = ints[:, 3]

                    ints2 = ints[:, 4:].reshape(nelements*2, 14)
                    nids = ints2[:, 0]
                    obj.element_node[itotal:itotal2, 0] = eids2
                    obj.element_node[itotal:itotal2, 1] = nids
                    obj.parent_coord[ielement:ielement2, 1] = parent
                    obj.parent_coord[ielement:ielement2, 1] = coord

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 32)[:, 4:]
                assert floats.shape[1] == 28, floats.shape
                # skipping [form1, form2]
                floats2 = floats.reshape(nelements*2, 14)
                isave1 = [2, 3, 4, 5, 6, 7]
                isave2 = [8, 9, 10, 11, 12, 13]
                #[xxb, force_x, shear_y, shear_z, torsion, bending_y, bending_z]
                real_imag = apply_mag_phase(floats2, is_magnitude_phase, isave1, isave2)
                obj.data[obj.itime, itotal:itotal2, 0] = floats2[:, 1].copy()
                obj.data[obj.itime, itotal:itotal2, 1:] = real_imag
            else:
                nnodes = 2
                if self.size == 4:
                    s1 = Struct(self._endian + b'i2i4s')
                else:
                    s1 = Struct(self._endian + b'q2q8s')
                s2 = Struct(mapfmt(self._endian + b'i13f', self.size))
                n = 0
                obj._times[obj.itime] = dt
                ntotal1 = 16 * self.factor  # 4*4
                ntotal2 = 56 * self.factor  # 14*4
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal1]
                    n += ntotal1

                    out = s1.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('OEF_Force_191-%s - %s\n' % (
                            self.element_type, str(out)))
                    (eid_device, parent, coord, icord) = out

                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)

                    unused_forces = []
                    for unused_i in range(nnodes):
                        edata = data[n:n+ntotal2]
                        n += ntotal2
                        out = s2.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('%s\n' % str(out))
                        [vugrid, posit,
                         force_xr, shear_yr, shear_zr, torsionr, bending_yr, bending_zr,
                         force_xi, shear_yi, shear_zi, torsioni, bending_yi, bending_zi] = out

                        if is_magnitude_phase:
                            force_x = polar_to_real_imag(force_xr, force_xi)
                            shear_y = polar_to_real_imag(shear_yr, shear_yi)
                            shear_z = polar_to_real_imag(shear_zr, shear_zi)
                            torsion = polar_to_real_imag(torsionr, torsioni)
                            bending_y = polar_to_real_imag(bending_yr, bending_yi)
                            bending_z = polar_to_real_imag(bending_zr, bending_zi)
                        else:
                            force_x = complex(force_xr, force_xi)
                            shear_y = complex(shear_yr, shear_yi)
                            shear_z = complex(shear_zr, shear_zi)
                            torsion = complex(torsionr, torsioni)
                            bending_y = complex(bending_yr, bending_yi)
                            bending_z = complex(bending_zr, bending_zi)

                        obj._add_sort1(dt, eid, parent, coord, icord,
                                       vugrid, posit, force_x, shear_y, shear_z, torsion, bending_y, bending_z)
        else:  # pragma: no cover
            msg = self.code_information()
            return self._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

def oef_cbar_real(self, data, obj: RealCBarForceArray, nelements, ntotal) -> int:
    n = 0
    fmt = mapfmt(self._endian + self._analysis_code_fmt + b'8f', self.size)
    s = Struct(fmt)
    if self.is_sort1:
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]

            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('OEF_CBar - %s\n' % (str(out)))
            (eid_device, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq) = out
            eid, dt = get_eid_dt_from_eid_device(
                eid_device, self.nonlinear_factor, self.sort_method)
            #data_in = [eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq]
            obj.add_sort1(dt, eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq)
            n += ntotal
    else:
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]

            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('OEF_CBar - %s\n' % (str(out)))
            (eid_device, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq) = out
            eid, dt = get_eid_dt_from_eid_device(
                eid_device, self.nonlinear_factor, self.sort_method)
            #data_in = [eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq]
            obj.add_sort2(dt, eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq)
            n += ntotal
    return n

def oef_cbar_imag(self, data, obj: ComplexCBarForceArray, nelements, ntotal, is_magnitude_phase) -> int:
    n = 0

    fmt = mapfmt(self._endian + self._analysis_code_fmt + b'16f', self.size)
    s = Struct(fmt)
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]

        out = s.unpack(edata)
        (eid_device,
         bm1ar, bm2ar, bm1br, bm2br, ts1r, ts2r, afr, trqr,
         bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi) = out
        if self.is_debug_file:
            self.binary_debug.write('OEF_CBar - %s\n' % (str(out)))
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)

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
            #eid_device, eid, dt, self.get_element_type(self.element_type)), data_in)
        obj.add_sort1(dt, eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq)
        n += ntotal
    return n

def oef_cbeam_real_100(self, data: bytes, obj: RealCBeamForceArray,
                       nelements: int, ntotal: int, dt: Any) -> int:
    n = 0
    ntotal1 = 4
    ntotal2 = 36
    s1 = self.struct_i
    s2 = Struct(self._endian + b'i8f')  # 36
    sort_method = self.sort_method
    for unused_i in range(nelements):
        edata = data[n:n+ntotal1]
        eid_device, = s1.unpack(edata)
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, sort_method)
        n += ntotal1

        for istation in range(11):
            edata = data[n:n+ntotal2]
            out = s2.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('OEF_Beam - %s\n' % (str(out)))
            (nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq) = out

            if istation == 0:  # isNewElement
                obj.add_new_element_sort1(
                    dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq)
            elif sd > 0.:
                obj.add_sort1(
                    dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq)
            n += ntotal2
    return n

def oef_cbeam_imag_177(self, data, obj: ComplexCBeamForceArray,
                       nelements, ntotal, is_magnitude_phase) -> int:
    n = 0
    #s1 = self.struct_i
    ntotal1 = 4 * self.factor
    ntotal2 = 64 * self.factor
    s1 = Struct(mapfmt(self._endian + b'i', self.size)) # self.struct_i
    s2 = Struct(mapfmt(self._endian + b'i15f', self.size))
    for unused_i in range(nelements):
        edata = data[n:n + ntotal1]
        eid_device, = s1.unpack(edata)
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)

        n += ntotal1
        for istation in range(11):
            edata = data[n:n+ntotal2]
            n += ntotal2
            out = s2.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('OEF_Beam - %s\n' % (str(out)))
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
            obj.add_sort1(
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
    n = 0
    size = self.size
    #                                5 6  7 8-i/f 9
    s1 = Struct(self._endian + b'i8sif 4s f i     4s')
    s2 = Struct(self._endian + b'i8sif 4s f f     4s')
    eid_old = None
    for unused_i in range(nelements):
        #2 THEORY(2) CHAR4 Theory
        #4 LAMID     I Lamina number
        #5 FP       RS Failure index for direct stresses
        #6 FM       RS Failure mode for maximum strain theory
        #7 FB       RS Failure index for interlaminar shear stress or -1
        #8 FMAX     RS Maximum of FP and FB or -1.
        #9 FFLAG CHAR4 Failure flag
        edata = data[n:n+ntotal]  # 4*9
        out = s1.unpack(edata)

        # failure_stress_for_ply = failure_strain_for_ply = failure_index_for_ply???
        # i    8s               i      f
        (eid, failure_theoryb, ply_id, failure_stress_for_ply,
         # 4s   7-f/i                8-f/i      9-4s
         flagb, interlaminar_stress, max_value, failure_flagb,
         #failure_index_for_bonding,
         #failure_index_for_element,
         #flag,
         #direct_stress_or_strain,
         #interlaminar_stress,
         #max_of_fb_fp_for_all_plies
        ) = out

        failure_theory = reshape_bytes_block_strip(failure_theoryb, size=size)
        flag = reshape_bytes_block_strip(flagb, size=size)
        failure_flag = reshape_bytes_block_strip(failure_flagb, size=size)

        if max_value == -1:
            max_value = np.nan
        else:
            max_value = s2.unpack(edata)[6]

        if eid == -1:
            #print(f'  ply_id={ply_id} failure_stress_for_ply={failure_stress_for_ply} '
                  #f'flag={flag} interlaminar_stress={interlaminar_stress} max_value={max_value} failure_flag={failure_flag}')
            eid = eid_old
        else:
            #print(f"eid={eid} ft='{failure_theory}'\n"
                  #f'  ply_id={ply_id} failure_stress_for_ply={failure_stress_for_ply} '
                  #f'flag={flag} interlaminar_stress={interlaminar_stress} max_value={max_value} failure_flag={failure_flag}')
            eid_old = eid
        assert flag in ['', '-1', '-2', '-12', 'IN'], f'flag={flag!r} flagb={flagb!r}'
        assert failure_theory in ['TSAI-WU', 'STRAIN', 'HILL', 'HOFFMAN', ''], 'failure_theory={failure_theory!r}'
        assert  failure_flag in ['', '***'], 'failure_flag=%r' % failure_flag
        obj.add_sort1(dt, eid, failure_theory, ply_id, failure_stress_for_ply, flag,
                      interlaminar_stress, max_value, failure_flag)
        n += ntotal

    #s = Struct(self._endian + b'i8si4f4s')
    #for i in range(nelements):
        #if i % 10000 == 0:
            #print 'i = ', i
        #edata = data[n:n+ntotal]  # 4*9
        #out = s.unpack(edata)
        #(eid_device, theory, lamid, failure_index_direct_stress, failure_mode_max_shear,
                 #failure_index_interlaminar_shear, fmax, failure_flag) = out
        #eid, dt = get_eid_dt_from_eid_device(
            #eid_device, self.nonlinear_factor, self.sort_method)
        #if self.is_debug_file:
            #if eid > 0:
                #self.binary_debug.write('  eid=%i; C=[%s]\n' % (', '.join(['%r' % di for di in out]) ))
            #else:
                #self.binary_debug.write('      %s  C=[%s]\n' % (' ' * len(str(eid)), ', '.join(['%r' % di for di in out]) ))

        #if eid > 0:
            #obj.add_new_eid_sort1(eType, dt, eid, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)
        #else:
            #obj.add_sort1(dt, eid, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)
        #n += ntotal
    return n
