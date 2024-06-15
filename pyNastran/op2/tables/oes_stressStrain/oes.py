#pylint: disable=C0301,W0201,R0911,R0915,R0914,C0103,W0212
"""
Defines the Real/Complex Stresses/Strains created by:
    STRESS = ALL
    STRAIN = ALL

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
from __future__ import annotations
from struct import Struct
from typing import Union, Any, TYPE_CHECKING
from numpy import fromstring, frombuffer, radians, sin, cos, vstack, repeat, array
import numpy as np

from pyNastran.op2.op2_interface.op2_codes import SORT1_TABLES_BYTES, TABLES_BYTES
from pyNastran.op2.op2_interface.utils import (
    mapfmt, mapfmt8,
    apply_mag_phase, build_obj, reshape_bytes_block_strip
)
from pyNastran.op2.op2_helper import polar_to_real_imag
from pyNastran.op2.op2_interface.function_codes import func1, func7

from pyNastran.op2.tables.utils import get_eid_dt_from_eid_device
from pyNastran.op2.tables.oug.oug import get_shock_prefix_postfix
from pyNastran.op2.tables.oes_stressStrain.real.oes_bars import RealBarStressArray, RealBarStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_bars100 import RealBar10NodesStressArray, RealBar10NodesStrainArray

from pyNastran.op2.tables.oes_stressStrain.real.oes_beams import (RealBeamStressArray, RealBeamStrainArray,
                                                                  RealNonlinearBeamStressArray)
from pyNastran.op2.tables.oes_stressStrain.real.oes_bush import RealBushStressArray, RealBushStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_bush1d import RealBush1DStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_composite_plates import RealCompositePlateStressArray, RealCompositePlateStrainArray
#RealCompositePlateStressStrengthRatioArray, RealCompositePlateStrainStrengthRatioArray = None, None
#RealCompositePlateStrainStrengthRatioArray = None
from pyNastran.op2.tables.oes_stressStrain.real.oes_composite_plates_strength_ratio import RealCompositePlateStressStrengthRatioArray # , RealCompositePlateStrainStrengthRatioArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_gap import NonlinearGapStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_plates import RealPlateStressArray, RealPlateStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_plate_strain import RealCPLSTRNPlateStressArray, RealCPLSTRNPlateStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_plate_strain_nx import RealCPLSTRNPlateStressNXArray, RealCPLSTRNPlateStrainNXArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_rods import RealRodStressArray, RealRodStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_shear import RealShearStrainArray, RealShearStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_solids import RealSolidStrainArray, RealSolidStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_solids_nx import RealSolidStressArrayNx, RealSolidStrainArrayNx
from pyNastran.op2.tables.oes_stressStrain.real.oes_solids_composite_nx import RealSolidCompositeStressArray, RealSolidCompositeStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_springs import (RealSpringStressArray, RealSpringStrainArray,
                                                                    RealNonlinearSpringStressArray)
from pyNastran.op2.tables.oes_stressStrain.real.oes_triax import RealTriaxStressArray, RealTriaxStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_bend import RealBendStressArray, RealBendStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_weld import RealWeldStressArray, RealWeldStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_fast import RealFastStressArray, RealFastStrainArray


from pyNastran.op2.tables.oes_stressStrain.complex.oes_bars import ComplexBarStressArray, ComplexBarStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_beams import ComplexBeamStressArray, ComplexBeamStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_bush import (ComplexCBushStressArray, ComplexCBushStrainArray)
from pyNastran.op2.tables.oes_stressStrain.complex.oes_fast import ComplexFastStressArray, ComplexFastStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_weld import ComplexWeldStressArray, ComplexWeldStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_bush1d import ComplexCBush1DStressArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_plates import ComplexPlateStressArray, ComplexPlateStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_composite_plates import ComplexLayeredCompositeStressArray, ComplexLayeredCompositeStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_plates_vm import (
    ComplexPlateVMStressArray, ComplexPlateVMStrainArray)
from pyNastran.op2.tables.oes_stressStrain.complex.oes_triax import ComplexTriaxStressArray, ComplexTriaxStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_rods import ComplexRodStressArray, ComplexRodStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_shear import ComplexShearStressArray, ComplexShearStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_solids import ComplexSolidStressArray, ComplexSolidStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_solids_vm import ComplexSolidStressVMArray, ComplexSolidStrainVMArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_springs import ComplexSpringStressArray, ComplexSpringStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_bend import ComplexBendStressArray, ComplexBendStrainArray

from pyNastran.op2.tables.oes_stressStrain.random.oes_rods import RandomRodStressArray, RandomRodStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_bars import RandomBarStressArray, RandomBarStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_beams import RandomBeamStressArray, RandomBeamStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_bend import RandomBendStressArray, RandomBendStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_plates import RandomPlateStressArray, RandomPlateStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_plates_vm import RandomPlateVMStressArray, RandomPlateVMStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_solids import RandomSolidStressArray, RandomSolidStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_solids_vm import RandomSolidVMStressArray, RandomSolidVMStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_shear import RandomShearStressArray, RandomShearStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_composite_plates import (
    RandomCompositePlateStressArray, RandomCompositePlateStrainArray,
    RandomCompositePlateStressVMArray, RandomCompositePlateStrainVMArray)

from pyNastran.op2.tables.oes_stressStrain.oes_nonlinear_rod import RealNonlinearRodArray
from pyNastran.op2.tables.oes_stressStrain.oes_nonlinear_bush import RealNonlinearBushArray
from pyNastran.op2.tables.oes_stressStrain.oes_hyperelastic import (
    HyperelasticQuadArray)
from pyNastran.op2.tables.oes_stressStrain.oes_nonlinear import RealNonlinearPlateArray, RealNonlinearSolidArray
from pyNastran.op2.tables.oes_stressStrain.cplstn import (
    oes_cplstn3_real_6, oes_cplstn4_real_32, oes_cplstn6_real_26,
)
from pyNastran.op2.tables.oes_stressStrain.utils import (
    oes_cbar_complex_19,
    oes_cbar100_real_10,
    oesrt_comp_shell_real_9, oesrt_cquad4_95,
    oes_cbar_real_16, oes_cbar_random_10,
    oes_cbeam_complex_111, oes_cbeam_random_67, oes_cbeam_real_111,
    oes_cbend_real_21,
    oes_cbush_complex_13, oes_cbush_real_7,
    oes_celas_complex_3, oes_celas_real_2,
    oes_cshear_real_4, oes_cshear_complex_5, oes_cshear_random_3,
    oes_comp_shell_real_11,
    oes_weldp_msc_real_8, oes_weldp_msc_complex_15,
    oes_fastp_msc_real_7, oes_fastp_msc_complex_13,
    oes_cquad4_33_complex_15, oes_cquad4_33_random_9,
    oes_cquad4_33_complex_vm_17, oes_cquad4_33_random_vm_11, oes_cquad4_random_vm_57,
    oes_cquad4_144_complex_77, oes_cquad4_144_random, oes_cquad4_144_real,
    oes_cquad4_complex_vm_87, oes_crod_real_5, oes_crod_complex_5, oes_crod_random_3,
    oes_csolid_real, oes_csolid_complex, oes_csolid_random, _oes_csolid2_real,
    oes_csolid_composite_real,
    oes_csolid_linear_hyperelastic_cosine_real, oes_csolid_linear_hyperelastic_real,
    oes_ctria3_real_17, oes_ctria3_complex_vm_17, oes_ctria3_random_9,
    oes_ctriax6_real_33, oes_ctriax_complex_37,
    oes_ctria3_random_vm_11,
    oes_quad4_33_real_17,
    oes_shell_composite_complex_11, oes_shell_composite_complex_13,
)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2
NX_TABLES_BYTES = [b'OESVM1', b'OESVM2']
NASA_TABLES_BYTES = [b'OESC1']

class OP2Common2:
    def __init__(self, op2: OP2):
        self.op2 = op2

    @property
    def size(self) -> int:
        return self.op2.size
    @property
    def factor(self) -> int:
        return self.op2.factor

class OES(OP2Common2):
    """
    Defines  the OES class that is used to read stress/strain data
    """
    def __init__(self, op2: OP2):
        super().__init__(op2)
        self.ntotal = 0

    def _read_oes1_3(self, data, unused_ndata):
        """
        reads OES1 subtable 3
        """
        op2 = self.op2
        op2._analysis_code_fmt = b'i'
        op2._data_factor = 1
        op2.words = [
            'aCode', 'tCode', 'element_type', 'isubcase',
            '???', '???', '???', 'load_set'
            'format_code', 'num_wide', 's_code', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', 'Title', 'subtitle', 'label', '???']

        op2.parse_approach_code(data)  # 3

        ## element type
        op2.element_type = op2.add_data_parameter(data, 'element_type', b'i', 3, False)

        ## load set ID
        op2.load_set = op2.add_data_parameter(data, 'load_set', b'i', 8, False)

        ## format code
        op2.format_code = op2.add_data_parameter(data, 'format_code', b'i', 9, False)

        ## number of words per entry in record
        ## .. note:: is this needed for this table ???
        op2.num_wide = op2.add_data_parameter(data, 'num_wide', b'i', 10, False)

        ## stress/strain codes
        op2.s_code = op2.add_data_parameter(data, 's_code', b'i', 11, False)

        ## thermal flag; 1 for heat ransfer, 0 otherwise
        op2.thermal = op2.add_data_parameter(data, 'thermal', b'i', 23, False)

        ## assuming tCode=1
        analysis_code = op2.analysis_code
        if analysis_code == 1:   # statics / displacement / heat flux
            op2._analysis_code_fmt = b'i'
            ## load set number
            op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
            op2.setNullNonlinearFactor()
        elif analysis_code == 2:  # real eigenvalues
            op2._analysis_code_fmt = b'i'
            #: mode number
            op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            #: eigenvalue
            op2.eign = op2.add_data_parameter(data, 'eign', b'f', 6, False)
            #: mode or cycle TODO confused on the type - F1 means float/int???
            op2.mode2 = op2.add_data_parameter(data, 'mode2', b'i', 7, False)
            op2.cycle = op2.add_data_parameter(data, 'cycle', b'f', 7, False)
            op2._op2_readers.reader_oug.update_mode_cycle('cycle')
            op2.data_names = op2.apply_data_code_value('data_names', ['mode', 'eign', 'mode2', 'cycle'])
        #elif analysis_code == 3: # differential stiffness
            #self.lsdvmn = self.get_values(data,'i',5) ## load set number
            #op2.data_code['lsdvmn'] = self.lsdvmn
        #elif analysis_code == 4: # differential stiffness
        #    self.lsdvmn = self.get_values(data,'i',5) ## load set number
        elif analysis_code == 5:   # frequency
            ## frequency
            op2.freq = op2.add_data_parameter(data, 'freq', b'f', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['freq'])
        elif analysis_code == 6:  # transient
            ## time step
            op2.dt = op2.add_data_parameter(data, 'dt', b'f', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['dt'])
        elif analysis_code == 7:  # pre-buckling
            ## load set
            op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
        elif analysis_code == 8:  # post-buckling
            ## mode number
            op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)  # real eigenvalue
            op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn', 'eigr'])
        elif analysis_code == 9:  # complex eigenvalues
            ## mode number
            op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            ## real eigenvalue
            op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            ## imaginary eigenvalue
            op2.eigi = op2.add_data_parameter(data, 'eigi', b'f', 7, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['mode', 'eigr', 'eigi'])
        elif analysis_code == 10:  # nonlinear statics
            ## load step
            self.lftsfq = op2.add_data_parameter(data, 'lftsfq', b'f', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['lftsfq'])
        elif analysis_code == 11:  # old geometric nonlinear statics
            ## load set number
            self.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
        elif analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## Time step ??? --> straight from DMAP
            op2.dt = op2.add_data_parameter(data, 'dt', b'f', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['dt'])
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % op2.analysis_code
            raise RuntimeError(msg)
        # tcode=2
        #if op2.analysis_code==2: # sort2
        #    self.lsdvmn = self.get_values(data,'i',5)

        #print(f'tCode={op2.tCode} -> result_type={result_type}')
        #print(op2.code_information())

        op2.fix_format_code()
        self.get_oes_prefix_postfix()
        op2._parse_thermal_code()
        op2._op2_readers.reader_oef._set_force_stress_element_name()
        if op2.is_debug_file:
            op2.binary_debug.write('  element_name   = %r\n' % op2.element_name)
            op2.binary_debug.write('  approach_code  = %r\n' % op2.approach_code)
            op2.binary_debug.write('  tCode          = %r\n' % op2.tCode)
            op2.binary_debug.write('  isubcase       = %r\n' % op2.isubcase)

        op2._read_title(data)

        try:
            op2.element_name = op2.element_mapper[op2.element_type]
        except KeyError:
            op2.log.error(op2.code_information())
            raise
        assert op2.element_name != '', op2.code_information()
        #if op2.element_type not in self.element_mapper:
            #return op2._not_implemented_or_skip(data, ndata, op2.code_information())

        self._parse_stress_code_to_stress_bits()
        op2._write_debug_bits()
        assert isinstance(op2.format_code, int), op2.format_code
        #print('op2.nonlinear_factor =', op2.nonlinear_factor)
        #assert op2.num_wide != 146, op2.code_information()
        #self._check_result_type()

    def _check_result_type(self):
        op2 = self.op2
        sort_method = func1(op2.tCode)
        if sort_method == 1:  # SORT1
            assert op2.sort_bits.is_sort1 == 1, op2.code_information()
        elif sort_method == 2:  # SORT2
            assert op2.sort_bits.is_sort1 == 0, op2.code_information()
        else:
            raise NotImplementedError(sort_method)


        result_type = op2.result_type #  func7(op2.tCode)
        if result_type == 0:
            assert op2.sort_bits.is_real == 1, op2.code_information()
        elif result_type == 1:
            assert op2.sort_bits.is_complex == 1, op2.code_information()
        elif result_type == 2:
            assert op2.sort_bits.is_random == 1, op2.code_information()
        else:
            raise NotImplementedError(result_type)

    def _parse_stress_code_to_stress_bits(self):
        """
        s_code =  0 -> stress_bits = [0,0,0,0,0]
        s_code =  1 -> stress_bits = [0,0,0,0,1]
        s_code =  2 -> stress_bits = [0,0,0,1,0]
        s_code =  3 -> stress_bits = [0,0,0,1,1]
        etc.
        s_code = 32 -> stress_bits = [1,1,1,1,1]

        stress_bits[0] = 0 -> is_max_shear=True       isVonMises=False
        stress_bits[0] = 1 -> is_max_shear=False      isVonMises=True

        stress_bits[1] = 0 -> is_stress=True        is_strain=False
        stress_bits[2] = 0 -> isFiberCurvature=True isFiberDistance=False
        stress_bits[3] = 0 -> duplicate of Bit[1] (stress/strain)
        stress_bits[4] = 0 -> material coordinate system flag

        """
        op2 = self.op2
        bits = [0, 0, 0, 0, 0]

        s_code = op2.s_code
        i = 4
        while s_code > 0:
            value = s_code % 2
            s_code = (s_code - value) // 2
            bits[i] = value
            i -= 1
        op2.stress_bits = bits
        op2.data_code['stress_bits'] = op2.stress_bits

    def _read_oes2_4(self, data: bytes, ndata: int):
        """
        Reads the Stress Table 4
        """
        op2 = self.op2
        if op2.table_name in NX_TABLES_BYTES:
            op2.to_nx(f' because table_name={op2.table_name}')
        elif op2.table_name in NASA_TABLES_BYTES:
            op2.to_nasa(f' because table_name={op2.table_name}')

        #assert self.isubtable == -4, self.isubtable
        #if op2.is_debug_file:
            #op2.binary_debug.write('  element_name = %r\n' % op2.element_name)
        #print("element_name =", op2.element_name)
        assert isinstance(op2.format_code, int), op2.format_code
        assert op2.is_stress is True, op2.code_information()
        op2.data_code['is_stress_flag'] = True
        op2.data_code['is_strain_flag'] = False

        op2._setup_op2_subcase('STRESS/STRAIN')
        elements_to_read = [
            1, 3, 10,  # CROD, CTUBE, CTUBE
            11, 12, 13,  # CELAS1, CELAS2, CELAS3,
            2, 4, 34, 33, 74,  # CBEAM, CSHEAR, CBAR, CQUAD4, CTRIA3,
            75, 64, 70, 82, 144,  # CTRIA6, CQUAD8, CTRIAR, CQUADR, CQUAD4
            69,  # CBEND
            67, 68, 95, 102,  # #CHEXA, CPENTA, QUAD4-comp, CBUSH
            39,  #CTETRA
            86,  # GAPNL
            88,  # TRIA3-nonlinear
            89,  # ROD-nonlinear
            90,  # QUAD4-nonlinear
            91,  # PENTANL
            93,  # HEXANL
            97,  # CTRIA3-C
            96,  # QUAD8-nonlinear
            98,  # TRIA6-nonlinear
            100,  # CBAR-100
            228,  # CQUADR
            232,  # CQUADR-composite
            243,  # CQUADX4
            189,  # VUQUAD
            190,  # VUTRIA
            191,  # VUBEAM
            256,  # CPYRAM
            227,  # CTRIAR
            275,  # CPLSTS3
        ]
        if op2.element_type in elements_to_read:
            n = self._read_oes_4_sort(data, ndata)
        else:
            msg = op2.code_information()
            n = op2._not_implemented_or_skip(data, ndata, msg)
        return n

    def _read_oes1_4(self, data: bytes, ndata: int):
        """
        Reads the Stress Table 4
        """
        op2 = self.op2
        if op2.table_name in NX_TABLES_BYTES:
            op2.to_nx(f' because table_name={op2.table_name}')
        #assert self.isubtable == -4, self.isubtable
        #if op2.is_debug_file:
            #op2.binary_debug.write('  element_name = %r\n' % op2.element_name)
        #print "element_name =", op2.element_name
        assert isinstance(op2.format_code, int), op2.format_code
        #assert op2.is_stress is True, op2.code_information()
        op2.data_code['is_stress_flag'] = True
        op2.data_code['is_strain_flag'] = False

        op2._setup_op2_subcase('STRESS/STRAIN')
        n = self._read_oes_4_sort(data, ndata)
        return n

    def _read_oes2_3(self, data, unused_ndata):
        """
        reads OES1 subtable 3
        """
        op2 = self.op2
        op2._data_factor = 1
        op2.words = [
            'aCode', 'tCode', 'element_type', 'isubcase',
            '???', '???', '???', 'load_set'
            'format_code', 'num_wide', 's_code', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', 'Title', 'subtitle', 'label', '???']

        op2.parse_approach_code(data)  # 3
        op2.sort_method = 2

        ## element type
        op2.element_type = op2.add_data_parameter(data, 'element_type', b'i', 3, False)

        ## load set ID
        op2.load_set = op2.add_data_parameter(data, 'load_set', b'i', 8, False)

        ## format code
        op2.format_code = op2.add_data_parameter(data, 'format_code', b'i', 9, False)

        ## number of words per entry in record
        ## .. note:: is this needed for this table ???
        op2.num_wide = op2.add_data_parameter(data, 'num_wide', b'i', 10, False)

        ## stress/strain codes
        op2.s_code = op2.add_data_parameter(data, 's_code', b'i', 11, False)

        ## thermal flag; 1 for heat ransfer, 0 otherwise
        op2.thermal = op2.add_data_parameter(data, 'thermal', b'i', 23, False)
        op2.element_id = op2.add_data_parameter(data, 'element_id', b'i', 5, fix_device_code=True)
        op2._element_id = op2.add_data_parameter(data, '_element_id', b'f', 5, apply_nonlinear_factor=False, add_to_dict=True)

        if op2.analysis_code == 1:  # static...because reasons.
            op2._analysis_code_fmt = b'i'
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
            op2.apply_data_code_value('analysis_method', 'dt')
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
        self._parse_stress_code_to_stress_bits()
        self._fix_oes_sort2(data)
        op2._op2_readers.reader_oef._set_force_stress_element_name()
        #assert isinstance(op2.nonlinear_factor, int), op2.nonlinear_factor

        #def get_format_code(is_sort2, is_complex, is_random):
            #format_code = 0
            #if is_complex:
                #format_code += 1
            #if is_sort2:
                #format_code += 2
            #if is_random:
                #format_code += 4
            #if format_code > 5:
                #format_code = 5
            ##print('is_complex =', is_complex)
            ##print('is_sort2 =', is_sort2)
            ##print('is_random =', is_random)
            ##print('format_code =', format_code)
            #return format_code

        #is_sort2 = True
        #is_complex = False
        #is_random = True
        #format_code = get_format_code(is_sort2, is_complex, is_random)
        #op2.format_code = format_code
        #op2.data_code['format_code'] = format_code
        #self._check_result_type()


    def _fix_oes_sort2(self, data):
        op2 = self.op2
        op2.fix_format_code()
        #print('op2.format_code_original =', op2.format_code_original)
        #print('op2.format_code =', op2.format_code)
            #op2.fix_format_code()
            #if op2.format_code == 1:
                #op2.format_code = 2
                #op2.data_code['format_code'] = 2
            #assert op2.format_code in [2, 3], op2.code_information()

        if 1:
            op2.fix_format_code()
            #if op2.num_wide == 8:
                #op2.format_code = 1
                #op2.data_code['format_code'] = 1
            #else:
                ##op2.fix_format_code()
                #if op2.format_code == 1:
                    #op2.format_code = 2
                    #op2.data_code['format_code'] = 2
                #assert op2.format_code in [2, 3], op2.code_information()

        op2._parse_thermal_code()
        if op2.is_debug_file:
            op2.binary_debug.write('  %-14s = %r %s\n' % ('approach_code', op2.approach_code,
                                                          op2.approach_code_str(op2.approach_code)))
            op2.binary_debug.write('  %-14s = %r\n' % ('tCode', op2.tCode))
            op2.binary_debug.write('  %-14s = %r\n' % ('isubcase', op2.isubcase))
        op2._read_title(data)
        op2._write_debug_bits()
        #assert isinstance(op2.nonlinear_factor, int), op2.nonlinear_factor

    def _read_ostr1_4(self, data: bytes, ndata: int):
        """
        Reads the Strain Table 4
        """
        op2 = self.op2
        if op2.table_name in NX_TABLES_BYTES:
            op2.to_nx(f' because table_name={op2.table_name} was found')
        #assert self.isubtable == -4, self.isubtable
        #if op2.is_debug_file:
            #op2.binary_debug.write('  element_name = %r\n' % op2.element_name)
        #print "element_name =", op2.element_name
        assert op2.is_strain is True, op2.code_information()
        op2.data_code['is_stress_flag'] = False
        op2.data_code['is_strain_flag'] = True

        op2._setup_op2_subcase('STRESS/STRAIN')
        n = self._read_ostr_4_sort(data, ndata)
        return n

    def _read_ostr2_4(self, data: bytes, ndata: int):
        """
        Reads the Strain Table 4
        """
        op2 = self.op2
        if op2.table_name in NX_TABLES_BYTES:
            op2.to_nx(f' because table_name={op2.table_name} was found')
        #assert self.isubtable == -4, self.isubtable
        #if op2.is_debug_file:
            #op2.binary_debug.write('  element_name = %r\n' % op2.element_name)
        #print("element_name =", op2.element_name)
        assert op2.is_strain is True, op2.code_information()
        op2.data_code['is_stress_flag'] = False
        op2.data_code['is_strain_flag'] = True

        op2._setup_op2_subcase('STRESS/STRAIN')
        if op2.element_type in [1, 3, 10,  # CROD, CTUBE, CTUBE
                                11, 12, 13,  # CELAS1, CELAS2, CELAS3,
                                2, 4, 34, 33, 74,  # CBEAM, CSHEAR, CBAR, CQUAD4, CTRIA3,
                                75, 64, 70, 82, 144,  # CTRIA6, CQUAD8, CTRIAR, CQUADR, CQUAD4
                                69,  # CBEND
                                67, 68, 95, 102,  #CHEXA, CPENTA, QUAD4-comp, CBUSH
                                96,  # QUAD8-nonlinear
                                98,  # TRIA6-nonlinear
                                39,  #CTETRA
                                228,  #CQUADR
                                232,  #CQUADR-composite
                                233,  #CTRIAR-composite
                                97]:  # CTRIA3-C
            n = self._read_ostr_4_sort(data, ndata)
        else:
            msg = op2.code_information()
            n = op2._not_implemented_or_skip(data, ndata, msg)
        return n

    #def _autojit3(func):
        #"""
        #Debugging function to print the object name and an needed parameters
        #"""
        #def new_func(self, data):
            #"""
            #The actual function exec'd by the decorated function.
            #"""
            #n = func(self, data)
            #return n
        #return new_func

    def _print_obj_name_on_crash(func):
        """Debugging function to print the object name and an needed parameters"""
        def new_func(self, data):
            """The actual function exec'd by the decorated function."""
            try:
                n = func(self, data)
            except NameError:
                raise
            except AttributeError:
                raise
            #except Exception:
                #raise
                #print("----------")
                #print(op2.obj)
                #print(op2.data_code)
                #if op2.obj is not None:
                    ##from pyNastran.utils import object_attributes
                    ##print object_attributes(op2.obj)
                    #print(op2.obj.data_code)
                #print("----------")
                #raise
            return n
        return new_func

    #@_print_obj_name_on_crash
    def _read_oes_4_sort(self, data: bytes, ndata: int):
        """Reads OES1 subtable 4 for NX/MSC/Autodesk/Optistruct"""
        op2 = self.op2
        #if op2.num_wide == 146:
            #assert op2.num_wide != 146
            #assert ndata != 146, op2.code_information()
        assert isinstance(op2.format_code, int), op2.format_code
        if op2.thermal == 0:
            n = self._read_oes1_loads(data, ndata)
        elif op2.thermal == 1:
            n = self._read_oes1_thermal(data, ndata)
        elif op2.thermal == 4:  # abs/nrl/srss - shock response spectra
            n = self._read_oes1_loads(data, ndata)
        else:
            msg = 'thermal=%s' % op2.thermal
            n = op2._not_implemented_or_skip(data, ndata, msg)
        return n

    #@_print_obj_name_on_crash
    def _read_ostr_4_sort(self, data: bytes, ndata: int):
        """
        Reads OSTR1 subtable 4
        """
        op2 = self.op2
        #if op2.num_wide == 146:
            #assert op2.num_wide != 146
            #assert ndata != 146, op2.code_information()
        if op2.thermal == 0:
            n = self._read_oes1_loads(data, ndata)
        elif op2.thermal == 1:
            n = self._read_oes1_thermal(data, ndata)
        else:
            msg = 'thermal=%s' % op2.thermal
            n = op2._not_implemented_or_skip(data, ndata, msg)
        return n

    def _read_oes1_thermal(self, unused_data: bytes, ndata: int) -> int:
        """
        Reads OES op2.thermal=1 tables; uses a hackish method to just skip the table
        """
        return ndata

    def _read_ostr1_thermal(self, unused_data: bytes, ndata: int) -> int:
        """
        Reads OSTR op2.thermal=1 tables; uses a hackish method to just skip the table
        """
        return ndata

    def get_stress_mapper(self):
        stress_mapper = {
            # element_type, format_code, num_wide

            # rods
            (1, 1, 5, b'OES1') : ('crod_stress', RealRodStressArray), # real
            (1, 1, 5, b'OES1X') : ('crod_stress', RealRodStressArray), # real
            (1, 1, 5, b'OES1X1') : ('crod_stress', RealRodStressArray), # real
            (1, 2, 5, b'OES1X') : ('crod_stress', ComplexRodStressArray), # real/imag
            (1, 3, 5, b'OES1X') : ('crod_stress', ComplexRodStressArray), # mag/phase
            (1, 2, 5, b'OESVM1') : ('crod_stress', ComplexRodStressArray), # real/imag
            (1, 3, 5, b'OESVM1') : ('crod_stress', ComplexRodStressArray), # mag/phase

            (3, 1, 5, b'OES1X1') : ('ctube_stress', RealRodStressArray),
            (3, 1, 5, b'OES1X') : ('ctube_stress', RealRodStressArray),
            (3, 2, 5, b'OES1X') : ('ctube_stress', ComplexRodStressArray),
            (3, 2, 5, b'OESVM1') : ('ctube_stress', ComplexRodStressArray),  # freq nx
            (3, 3, 5, b'OESVM1') : ('ctube_stress', ComplexRodStressArray),  # freq nx

            #(3, 3, 5) : ('ctube_stress', ComplexRodStressArray),

            (10, 1, 5, b'OES1') : ('conrod_stress', RealRodStressArray),
            (10, 1, 5, b'OES1X') : ('conrod_stress', RealRodStressArray),
            (10, 2, 5, b'OES1X') : ('conrod_stress', ComplexRodStressArray),
            (10, 1, 5, b'OES1X1') : ('conrod_stress', RealRodStressArray),
            (10, 2, 5, b'OESVM1') : ('conrod_stress', ComplexRodStressArray),
            (10, 3, 5, b'OESVM1') : ('conrod_stress', ComplexRodStressArray),
            #(10, 2, 5) : ('conrod_stress', ComplexRodStressArray),
            #(10, 3, 5) : ('conrod_stress', ComplexRodStressArray),

            # beams
            (2, 1, 111, b'OES1X1') : ('cbeam_stress', RealBeamStressArray),
            (2, 1, 111, b'OES1X') : ('cbeam_stress', RealBeamStressArray),
            (2, 1, 111, b'OES1') : ('cbeam_stress', RealBeamStressArray),
            (2, 2, 111, b'OES1X') : ('cbeam_stress', ComplexBeamStressArray),
            (2, 3, 111, b'OES1X') : ('cbeam_stress', ComplexBeamStressArray),
            (2, 3, 111, b'OESVM1') : ('cbeam_stress', ComplexBeamStressArray),

            (4, 1, 4, b'OES1X1') : ('cshear_stress', RealShearStressArray),
            #(4, 2, 5) : ('cshear_stress', ComplexShearStressArray),
            #(4, 3, 5) : ('cshear_stress', ComplexShearStressArray),
            (4, 2, 5, b'OES1X') : ('cshear_stress', ComplexShearStressArray),
            (4, 2, 5, b'OESVM1') : ('cshear_stress', ComplexShearStressArray),
            (4, 3, 5, b'OESVM1') : ('cshear_stress', ComplexShearStressArray),
            #(4, 3, 3) : ('cshear_stress', RandomShearStressArray),

            (11, 1, 2, b'OES1X1'): ('celas1_stress', RealSpringStressArray), # real
            (11, 2, 3, b'OES1X'): ('celas1_stress', ComplexSpringStressArray), # real/imag
            (11, 3, 3, b'OES1X'): ('celas1_stress', ComplexSpringStressArray), # mag/phase
            (11, 2, 3, b'OESVM1'): ('celas1_stress', ComplexSpringStressArray), # mag/phase
            (11, 3, 3, b'OESVM1'): ('celas1_stress', ComplexSpringStressArray), # mag/phase

            (12, 1, 2, b'OES1X1'): ('celas2_stress', RealSpringStressArray),
            (12, 1, 2, b'OES1X'): ('celas2_stress', RealSpringStressArray),
            (12, 1, 2, b'OES1'): ('celas2_stress', RealSpringStressArray),
            (12, 2, 3, b'OES1X'): ('celas2_stress', ComplexSpringStressArray),
            (12, 3, 3, b'OES1X'): ('celas2_stress', ComplexSpringStressArray),
            (12, 2, 3, b'OESVM1'): ('celas2_stress', ComplexSpringStressArray),
            (12, 3, 3, b'OESVM1'): ('celas2_stress', ComplexSpringStressArray),

            (13, 1, 2, b'OES1X1'): ('celas3_stress', RealSpringStressArray),
            #(13, 2, 3) : ('celas3_stress', ComplexSpringStressArray),
            #(13, 3, 3) : ('celas3_stress', ComplexSpringStressArray),
            (13, 2, 3, b'OES1X'): ('celas3_stress', ComplexSpringStressArray),
            (13, 2, 3, b'OESVM1'): ('celas3_stress', ComplexSpringStressArray),
            (13, 3, 3, b'OESVM1'): ('celas3_stress', ComplexSpringStressArray),

            (14, 1, 2) : ('celas4_stress', RealSpringStressArray),
            (14, 2, 3) : ('celas4_stress', ComplexSpringStressArray),
            (14, 3, 3) : ('celas4_stress', ComplexSpringStressArray),

            (34, 1, 16, b'OES1X1'): ('cbar_stress', RealBarStressArray),
            (34, 1, 16, b'OES1X'): ('cbar_stress', RealBarStressArray),
            (34, 1, 16, b'OES1'): ('cbar_stress', RealBarStressArray),
            (34, 2, 19, b'OES1X'): ('cbar_stress', ComplexBarStressArray),
            (34, 1, 10, b'OESNO1'): ('cbar_stress', ComplexBarStressArray),
            (34, 2, 10, b'OESXRMS1'): ('cbar_stress', ComplexBarStressArray),

            (34, 1, 10, b'OESRMS2'): ('cbar_stress', RandomBarStressArray),

            (34, 2, 10, b'OESPSD2'): ('cbar_stress', RandomBarStressArray),
            (34, 2, 10, b'OESRMS2'): ('cbar_stress', RandomBarStressArray),
            (34, 2, 10, b'OESNO2'): ('cbar_stress', RandomBarStressArray),
            (34, 2, 10, b'OESATO2'): ('cbar_stress', RandomBarStressArray),
            (34, 2, 10, b'OESCRM2'): ('cbar_stress', RandomBarStressArray),

            # Missing stress_mapper key for OES1 table #501
            # see cbarao_random_x_mini.op2 for an example with OES1 and OES1X...
            # it looks to be an error in MSC [2008-2012)
            (34, 2, 19, b'OES1'): ('cbar_stress', ComplexBarStressArray),
            (34, 3, 19, b'OES1X'): ('cbar_stress', ComplexBarStressArray),
            (34, 3, 19, b'OESVM1'): ('cbar_stress', ComplexBarStressArray),
            #(34, 1, 19) : ('cbar_stress', RandomBarStressArray),
            (100, 1, 10, b'OES1X1'): ('cbar_stress_10nodes', RealBar10NodesStressArray),
            (100, 1, 10, b'OES1X'): ('cbar_stress_10nodes', RealBar10NodesStressArray),

            # solids
            (39, 1, 109, b'OES1X1'): ('ctetra_stress', RealSolidStressArray), # real
            (39, 1, 109, b'OES1X'): ('ctetra_stress', RealSolidStressArray), # real
            (39, 1, 109, b'OES1'): ('ctetra_stress', RealSolidStressArray), # real
            (39, 3, 74, b'OESVM1'): ('ctetra_stress', ComplexSolidStressArray), # mag/phase

            (67, 1, 193, b'OES1X1'): ('chexa_stress', RealSolidStressArray),
            (67, 1, 193, b'OES1X'): ('chexa_stress', RealSolidStressArray),
            (67, 1, 193, b'OES1'): ('chexa_stress', RealSolidStressArray),
            (67, 1, 193, b'RASCONS'): ('chexa_stress', RealSolidStressArray),

            (68, 1, 151, b'OES1X1'): ('cpenta_stress', RealSolidStressArray),
            (68, 1, 151, b'OES1X'): ('cpenta_stress', RealSolidStressArray),
            (68, 1, 151, b'OES1'): ('cpenta_stress', RealSolidStressArray),
            (68, 3, 102, b'OESVM1'): ('cpenta_stress', ComplexSolidStressArray),

            (39, 2, 69, b'OES1X'): ('ctetra_stress', ComplexSolidStressArray), # real/imag
            (39, 2, 69, b'OES1'): ('ctetra_stress', ComplexSolidStressArray),
            (39, 2, 74, b'OESVM1'): ('ctetra_stress', 'NA'), # real/imag
            #(39, 3, 69) : ('ctetra_stress', ComplexSolidStressArray), # mag/phase

            (67, 2, 121, b'OES1X'): ('chexa_stress', ComplexSolidStressArray),
            (67, 3, 121, b'OES1X'): ('chexa_stress', ComplexSolidStressArray),
            (67, 3, 130, b'OESVM1'): ('chexa_stress', ComplexSolidStressArray),
            (67, 2, 121, b'OES1'): ('chexa_stress', ComplexSolidStressArray),
            (67, 3, 121, b'OES1'): ('chexa_stress', ComplexSolidStressArray),

            (68, 2, 95, b'OES1X'): ('cpenta_stress', ComplexSolidStressArray),
            (68, 3, 95, b'OES1X'): ('cpenta_stress', ComplexSolidStressArray),
            (68, 2, 95, b'OES1'): ('cpenta_stress', ComplexSolidStressArray),

            (33, 1, 17, b'OES1X1'):  ('cquad4_stress', RealPlateStressArray),
            (33, 1, 17, b'OES1X'):  ('cquad4_stress', RealPlateStressArray),
            (33, 1, 17, b'OES1'):  ('cquad4_stress', RealPlateStressArray),
            (33, 2, 15, b'OES1X'):  ('cquad4_stress', ComplexPlateStressArray),
            (33, 3, 15, b'OES1X'):  ('cquad4_stress', ComplexPlateStressArray),
            #(33, 3, 0) :  ('cquad4_stress', RandomPlateStressArray),
            (33, 1, 9, b'OESNO1'): ('cquad4_stress', ComplexPlateStressArray),
            (33, 2, 11, b'OESXRMS1'): ('cquad4_stress', ComplexPlateStressArray),

            (33, 2, 9, b'OESATO2'): ('cquad4_stress', 'NA'),
            (33, 2, 9, b'OESCRM2'): ('cquad4_stress', 'NA'),
            (33, 2, 9, b'OESPSD2'): ('cquad4_stress', 'NA'),
            (33, 2, 9, b'OESNO2'): ('cquad4_stress', 'NA'),

            (33, 1, 9, b'OESRMS2'): ('cquad4_stress', 'NA'),
            (33, 2, 9, b'OESRMS2'): ('cquad4_stress', 'NA'),


            (74, 1, 17, b'OES1X1'): ('ctria3_stress', RealPlateStrainArray),
            (74, 1, 17, b'OES1X'): ('ctria3_stress', RealPlateStrainArray),
            (74, 1, 17, b'OES1'): ('ctria3_stress', RealPlateStrainArray),
            (74, 2, 15, b'OES1X'): ('ctria3_stress', ComplexPlateStrainArray),
            (74, 3, 15, b'OES1X'): ('ctria3_stress', ComplexPlateStrainArray),
            (74, 2, 11, b'OESXRMS1'): ('ctria3_stress', ComplexPlateStrainArray),
            (74, 1, 9, b'OESNO1'): ('ctria3_stress', ComplexPlateStrainArray),
            (74, 2, 17, b'OESVM1'): ('ctria3_stress', 'NA'),
            (74, 3, 17, b'OESVM1'): ('ctria3_stress', 'NA'),

            (74, 1, 9, b'OESRMS2'): ('ctria3_stress', 'NA'),

            #(74, 1, 9) : ('ctria3_stress', RandomPlateStressArray),

            (82, 1, 87, b'OES1X1'): ('cquadr_stress', RealPlateStressArray),
            (82, 1, 87, b'OES1X'): ('cquadr_stress', RealPlateStressArray),
            (82, 2, 77, b'OES1X'): ('cquadr_stress', ComplexPlateStressArray),
            (82, 3, 77, b'OES1X') : ('cquadr_stress', ComplexPlateStressArray),

            (64, 1, 87, b'OES1X1') : ('cquad8_stress', RealPlateStressArray), # real
            (64, 1, 87, b'OES1X'): ('cquad8_stress', RealPlateStressArray),
            (64, 1, 87, b'OES1'): ('cquad8_stress', RealPlateStressArray),
            (64, 2, 77, b'OES1'): ('cquad8_stress', ComplexPlateStressArray), # real/imag
            (64, 3, 77, b'OES1'): ('cquad8_stress', ComplexPlateStressArray), # mag/phase
            (64, 2, 77, b'OES1X'): ('cquad8_stress', ComplexPlateStressArray), # real/imag
            (64, 3, 77, b'OES1X'): ('cquad8_stress', ComplexPlateStressArray), # mag/phase
            (64, 2, 87, b'OESVM1'): ('cquad8_stress', ComplexPlateStressArray), # real/imag
            (64, 3, 87, b'OESVM1'): ('cquad8_stress', ComplexPlateStressArray), # mag/phase

            (70, 1, 70, b'OES1X1'): ('ctriar_stress', RealPlateStressArray),
            (70, 1, 70, b'OES1X'): ('ctriar_stress', RealPlateStressArray),
            (70, 2, 62, b'OES1X'): ('ctriar_stress', ComplexPlateStressArray),
            (70, 3, 62, b'OES1X'): ('ctriar_stress', ComplexPlateStressArray),

            (75, 1, 70, b'OES1X1'): ('ctria6_stress', RealPlateStressArray),
            (75, 2, 62, b'OES1X'): ('ctria6_stress', ComplexPlateStressArray),
            (75, 3, 62, b'OES1X'): ('ctria6_stress', ComplexPlateStressArray),
            (75, 2, 70, b'OESVM1'): ('ctria6_stress', ComplexPlateStressArray),
            (75, 3, 70, b'OESVM1'): ('ctria6_stress', ComplexPlateStressArray),

            (144, 1, 87, b'OES1X1'): ('cquad4_stress', RealPlateStressArray),
            (144, 1, 87, b'OES1'): ('cquad4_stress', RealPlateStressArray),
            (144, 1, 87, b'RASCONS'): ('cquad4_stress', RealPlateStressArray),

            (144, 2, 77, b'OES1X'): ('cquad4_stress', ComplexPlateStressArray),
            (144, 3, 77, b'OES1X'): ('cquad4_stress', ComplexPlateStressArray),
            (144, 3, 87, b'OESVM1'): ('cquad4_stress', ComplexPlateStressArray),
            #(144, 3, 77) : ('cquad4_stress', ComplexPlateStressArray),
            #(64, 1, 47) : ('cquad8_stress', RandomPlateStressArray), # random
            #(70, 1, 39) : ('ctriar_stress', RandomPlateStressArray),
            #(75, 1, 39) : ('ctria6_stress', RandomPlateStressArray),
            #(82, 1, 47) : ('cquadr_stress', RandomPlateStressArray),
            #(144, 1, 47) : ('cquad4_stress', RandomPlateStressArray),

            (88, 1, 13, b'OESNLXR'): ('nonlinear_ctria3_stress', RealNonlinearPlateArray), # real
            (88, 1, 25, b'OESNL1X'): ('nonlinear_ctria3_stress', RealNonlinearPlateArray), # real?
            (88, 1, 25, b'OESNLXR'): ('nonlinear_ctria3_stress', RealNonlinearPlateArray), # real?

            (90, 1, 13, b'OESNLXR'): ('nonlinear_cquad4_stress', RealNonlinearPlateArray),
            (90, 1, 25, b'OESNL1X'): ('nonlinear_cquad4_stress', RealNonlinearPlateArray),
            (90, 1, 25, b'OESNLXR'): ('nonlinear_cquad4_stress', RealNonlinearPlateArray),
            (90, 1, 25, b'OESNLXD'): ('nonlinear_cquad4_stress', RealNonlinearPlateArray),

            (95, 1, 11, b'OES1C'): ('cquad4_composite_stress', RealCompositePlateStressArray), # real
            (95, 1, 11, b'OESCP'): ('cquad4_composite_stress', RealCompositePlateStressArray), # real
            (95, 1, 9, b'OESRT'): ('cquad4_composite_stress', 'RandomCompositePlateStressArray'), # real
            (95, 2, 11, b'OESCP'): ('cquad4_composite_stress', RealCompositePlateStressArray), # real?
            (95, 2, 11, b'OESRT'): ('cquad4_composite_stress', RealCompositePlateStressArray), # real?
            #(95, 2, 9) : ('cquad4_composite_stress', ComplexCompositePlateStressArray), # real/imag
            #(95, 3, 9) : ('cquad4_composite_stress', ComplexCompositePlateStressArray), # mag/phase

            #(96, 1, 9) : ('cquad8_composite_stress', 'RandomCompositePlateStressArray'),
            (96, 1, 11, b'OES1C'): ('cquad8_composite_stress', RealCompositePlateStressArray),
            #(96, 1, 11) : ('cquad8_composite_stress', RealCompositePlateStressArray),
            #(96, 2, 9) : ('cquad8_composite_stress', ComplexCompositePlateStressArray),
            #(96, 3, 9) : ('cquad8_composite_stress', ComplexCompositePlateStressArray),

            (97, 1, 9, b'OESRT'): ('ctria3_composite_stress', 'RandomCompositePlateStressArray'),
            (97, 1, 11, b'OES1C'): ('ctria3_composite_stress', RealCompositePlateStressArray),
            (97, 1, 11, b'OESCP'): ('ctria3_composite_stress', RealCompositePlateStressArray),
            (97, 2, 11, b'OESCP'): ('ctria3_composite_stress', RealCompositePlateStressArray),
            #(97, 2, 9) : ('ctria3_composite_stress', ComplexCompositePlateStressArray),
            #(97, 3, 9) : ('ctria3_composite_stress', ComplexCompositePlateStressArray),

            (98, 1, 9, b'OESRT'): ('ctria6_composite_stress', 'RandomCompositePlateStressArray'),
            (98, 1, 11, b'OES1C'): ('ctria6_composite_stress', RealCompositePlateStressArray),
            #(98, 1, 11) : ('ctria6_composite_stress', RealCompositePlateStressArray),
            #(98, 2, 9) : ('ctria6_composite_stress', ComplexCompositePlateStressArray),
            #(98, 3, 9) : ('ctria6_composite_stress', ComplexCompositePlateStressArray),

            (53, 1, 33, b'OES1X1'): ('ctriax_stress', RealTriaxStressArray),
            (53, 1, 33, b'OES1X'): ('ctriax_stress', RealTriaxStressArray),
            (53, 2, 37, b'OES1X'): ('ctriax_stress', ComplexTriaxStressArray),
            #(53, 3, 37) : ('ctriax_stress', ComplexTriaxStressArray),

            (102, 1, 7, b'OES1X1'): ('cbush_stress', RealBushStressArray),
            (102, 1, 7, b'OES1X'): ('cbush_stress', RealBushStressArray),
            (102, 1, 7, b'OES1'): ('cbush_stress', RealBushStressArray),
            (102, 2, 13, b'OES1X'): ('cbush_stress', ComplexCBushStressArray),
            (102, 3, 13, b'OES1X'): ('cbush_stress', ComplexCBushStressArray),
            (102, 2, 13, b'OESVM1'): ('cbush_stress', 'NA'),
            (102, 2, 13, b'OES1'): ('cbush_stress', ComplexCBushStressArray),

            (40, 1, 8, b'OES1X1'): ('cbush1d_stress_strain', RealBushStressArray),
            (40, 1, 8, b'OESNLXD'): ('cbush1d_stress_strain', RealBushStressArray),
            #(40, 2, 9) : ('cbush1d_stress_strain', ComplexCBushStressArray),
            #(40, 3, 9) : ('cbush1d_stress_strain', ComplexCBushStressArray),

            (87, 1, 7, b'OESNL1X'): ('nonlinear_ctube_stress', RealNonlinearRodArray),
            (87, 1, 7, b'OESNLXR'): ('nonlinear_ctube_stress', RealNonlinearRodArray),
            (89, 1, 7, b'OESNL1X'): ('nonlinear_crod_stress', RealNonlinearRodArray),
            (89, 1, 7, b'OESNLXD'): ('nonlinear_crod_stress', RealNonlinearRodArray),
            (89, 1, 7, b'OESNLXR'): ('nonlinear_crod_stress', RealNonlinearRodArray),
            (92, 1, 7, b'OESNL1X'): ('nonlinear_conrod_stress', RealNonlinearRodArray),
            (92, 1, 7, b'OESNLXD'): ('nonlinear_conrod_stress', RealNonlinearRodArray),
            (92, 1, 7, b'OESNLXR'): ('nonlinear_conrod_stress', RealNonlinearRodArray),

            (224, 1, 3, b'OESNLXD'): ('nonlinear_celas1_stress', RealNonlinearSpringStressArray),
            (224, 1, 3, b'OESNLXR'): ('nonlinear_celas1_stress', RealNonlinearSpringStressArray),
            (225, 1, 3, b'OESNLXR'): ('nonlinear_celas3_stress', RealNonlinearSpringStressArray),

            (35, 1, 18, b'OES1X1'): ('NA', 'NA'), # CCONEAX
            (35, 1, 18, b'OES1'): ('NA', 'NA'), # CCONEAX

            (60, 1, 10, b'OES1X'): ('NA', 'NA'), # DUM8/CCRAC2D
            (61, 1, 10, b'OES1X'): ('NA', 'NA'), # DUM8/CCRAC3D

            (69, 1, 21, b'OES1X1'): ('NA', 'NA'), # CBEND
            (69, 2, 21, b'OES1X'): ('NA', 'NA'), # CBEND
            (69, 3, 21, b'OES1X'): ('NA', 'NA'), # CBEND

            (86, 1, 11, b'OESNL1X'): ('nonlinear_cgap_stress', NonlinearGapStressArray),
            (86, 1, 11, b'OESNLXR'): ('nonlinear_cgap_stress', NonlinearGapStressArray),
            (86, 1, 11, b'OESNLXD'): ('nonlinear_cgap_stress', NonlinearGapStressArray),
            (94, 1, 51, b'OESNL1X'): ('nonlinear_cbeam_stress', RealNonlinearBeamStressArray),
            (94, 1, 51, b'OESNLXR'): ('nonlinear_cbeam_stress', RealNonlinearBeamStressArray),

            (85, 1, 82, b'OESNLXR'): ('NA', 'NA'),  # TETRANL
            (91, 1, 114, b'OESNLXD'): ('NA', 'NA'),  # PENTANL
            (91, 1, 114, b'OESNLXR'): ('NA', 'NA'),  # PENTANL
            (93, 1, 146, b'OESNL1X'): ('NA', 'NA'),  # HEXANL
            (93, 1, 146, b'OESNLXD'): ('NA', 'NA'),  # HEXANL
            (93, 1, 146, b'OESNLXR'): ('NA', 'NA'),  # HEXANL

            # 101-AABSF
            (101, 2, 4, b'OES1X'): ('NA', 'NA'),

            # 140-HEXA8FD
            (140, 1, 162, b'OES1X1'): ('NA', 'NA'),
            #201-QUAD4FD
            (201, 1, 46, b'OESNLXD'): ('NA', 'NA'),
            (201, 1, 46, b'OESNLXR'): ('NA', 'NA'),

            # 145-VUHEXA  (8 nodes)
            (145, 1, 98, b'OES1X1'): ('NA', 'NA'),
            (145, 2, 106, b'OES1X'): ('NA', 'NA'),
            (145, 3, 106, b'OES1X'): ('NA', 'NA'),
            # 146-VUPENTA (6 nodes)
            (146, 1, 74, b'OES1X1') : ('NA', 'NA'),
            (146, 2, 80, b'OES1X'): ('NA', 'NA'),
            (146, 3, 80, b'OES1X'): ('NA', 'NA'),
            # 147-VUTETRA (4 nodes)
            (147, 1, 50, b'OES1X1'): ('NA', 'NA'),
            (147, 2, 54, b'OES1X'): ('NA', 'NA'),
            (147, 3, 54, b'OES1X'): ('NA', 'NA'),

            # 139-QUAD4FD
            # self.hyperelastic_cquad4_strain, HyperelasticQuad
            (139, 1, 30, b'OES1X1'): ('NA', 'NA'),

            # 189-VUQUAD
            (189, 1, 74, b'OES1X1'): ('NA', 'NA'),
            (189, 2, 114, b'OES1X'): ('NA', 'NA'),

            # 47-AXIF2
            (47, 2, 9, b'OES1X'): ('axif2', 'NA'),
            # 48-AXIF3
            (48, 2, 19, b'OES1X'): ('axif3', 'NA'),
            # 190-VUTRIA
            (190, 1, 57, b'OES1X1'): ('NA', 'NA'),
            (190, 2, 87, b'OES1X'): ('NA', 'NA'),
            (190, 3, 87, b'OES1X'): ('NA', 'NA'),

            # 191-VUBEAM
            #(191, 1, 60, b'OES1X1') : ('vubeam', 'NA'),
            #(191, 2, 80, b'OES1X') : ('vubeam', 'NA'),
            #(191, 3, 80, b'OES1X') : ('vubeam', 'NA'),

            # 203-SLIF1D?
            (203, 1, 14, b'OESNLBR'): ('slif1d', 'NA'),
            # 50-SLOT3
            (50, 2, 11, b'OES1X'): ('slot3', 'NA'),
            # 51-SLOT4
            (51, 2, 13, b'OES1X'): ('slot4', 'NA'),

            # 160-PENTA6FD
            (160, 1, 122, b'OES1X1'): ('cpenta', 'NA'),
            # 161-TETRA4FD
            (161, 1, 22, b'OES1X1'): ('ctetra', 'NA'),
            # 162-TRIA3FD
            (162, 1, 9, b'OES1X1'): ('ctria', 'NA'),
            # 163-HEXAFD
            (163, 1, 542, b'OES1X1'): ('chexa', 'NA'),
            # 164-QUADFD
            (164, 1, 65, b'OES1X1'): ('cquad', 'NA'),
            # 165-PENTAFD
            (165, 1, 422, b'OES1X1'): ('cpenta', 'NA'),
            # 166-TETRAFD
            (166, 1, 102, b'OES1X1'): ('ctetra', 'NA'),
            # 167-TRIAFD
            (167, 1, 23, b'OES1X1'): ('NA', 'NA'),
            # 168-TRIAX3FD
            (168, 1, 9, b'OES1X1'): ('ctriax3', 'NA'),
            # 169-TRIAXFD
            (169, 1, 23, b'OES1X1'): ('ctriax', 'NA'),
            # 170-QUADX4FD
            (170, 1, 30, b'OES1X1'): ('cquadx4fd', 'NA'),
            # 171-QUADXFD
            (171, 1, 65, b'OES1X1'): ('cquadx', 'NA'),
            # 172-QUADRNL
            (172, 1, 25, b'OESNLXR'): ('cquadrnl', 'NA'),
            # 202-HEXA8FD
            (202, 1, 122, b'OESNLXD'): ('chexa', 'NA'),
            (202, 1, 122, b'OESNLXR'): ('chexa', 'NA'),
            # 204-PENTA6FD
            (204, 1, 92, b'OESNLXR'): ('cpenta', 'NA'),
            # 211-TRIAFD
            (211, 1, 35, b'OESNLXR'): ('ctria3', 'NA'),
            # 213-TRIAXFD
            (213, 1, 35, b'OESNLXR'): ('ctriax', 'NA'),
            # 214-QUADX4FD
            (214, 1, 46, b'OESNLXR'): ('cquadx4', 'NA'),
            # 216-TETRA4FD
            (216, 1, 62, b'OESNLXD'): ('NA', 'NA'),
            (216, 1, 62, b'OESNLXR'): ('NA', 'NA'),
            # 217-TRIA3FD
            (217, 1, 35, b'OESNLXR') : ('ctria3', 'NA'),
            # 218-HEXAFD
            (218, 1, 122, b'OESNLXR'): ('chexa', 'NA'),
            # 219-QUADFD
            (219, 1, 46, b'OESNLXR'): ('cquad', 'NA'),
            # 220-PENTAFD
            (220, 1, 92, b'OESNLXR'): ('cpenta', 'NA'),
            # 221-TETRAFD
            (221, 1, 62, b'OESNLXR'): ('tetrafd', 'NA'),
            # 222-TRIAX3FD
            (222, 1, 35, b'OESNLXR'): ('ctriax3fd', 'NA'),
            # 223-CQUADXFD
            (223, 1, 46, b'OESNLXR'): ('cquadx', 'NA'),
            # 226-BUSH
            (226, 1, 19, b'OESNLXD'): ('cbush', 'NA'),
            (226, 1, 19, b'OESNLXR'): ('cbush', 'NA'),
            # 227-CTRIAR
            (227, 1, 17, b'OES1X1'): ('ctriar', 'NA'),
            (227, 1, 17, b'OES1X'): ('ctriar', 'NA'),
            # 228-CQUADR
            (228, 1, 17, b'OES1X1'): ('cquadr', 'NA'),
            (228, 1, 17, b'OES1X'): ('cquadr', 'NA'),
            # 232-QUADRLC
            (232, 1, 11, b'OES1C'): ('cquadr', 'NA'),
            (232, 1, 11, b'OESCP'): ('cquadr', 'NA'),
            (232, 2, 13, b'OESVM1C'): ('cquadr', 'NA'),  # freq nx
            (232, 3, 13, b'OESVM1C'): ('cquadr', 'NA'),  # freq nx
            #(234, 1, 11) : ('cquadr', 'NA'), # bad?
            # 233-TRIARLC
            (233, 1, 11, b'OES1C'): ('ctriar', 'NA'),
            (233, 2, 13, b'OESVM1C'): ('ctriar', 'NA'),  # freq nx
            (233, 3, 13, b'OESVM1C'): ('ctriar', 'NA'),  # freq nx
            # 235-CQUADR
            (235, 1, 17, b'OES1X1'): ('NA', 'NA'),
            (235, 2, 15, b'OES1X'): ('NA', 'NA'),

            # 242-CTRAX
            # 244-CTRAX6
            (242, 1, 34, b'OES1X1'): ('ctrax', 'NA'),
            (244, 1, 34, b'OES1X1'): ('ctrax6', 'NA'),

            # 243-CQUADX4
            # 245-CQUADX8
            (243, 1, 42, b'OES1X1'): ('cquadx4', 'NA'),
            (245, 1, 42, b'OES1X1'): ('cquadx8', 'NA'),

            #256-CPYRAM
            (255, 1, 130, b'OES1X1'): ('cpyram', 'NA'),
            (255, 2, 82, b'OES1X'): ('cpyram', 'NA'),
            (256, 1, 98, b'OESNLXD'): ('cpyram', 'NA'),

            # 271-CPLSTN3
            # 272-CPLSTN4
            (271, 1, 6, b'OES1X1'): ('cplstn3', 'NA'),
            (271, 1, 6, b'OES1X'): ('cplstn3', 'NA'),
            (272, 1, 32, b'OES1X1') : ('cplstn4', 'NA'),
            (272, 1, 32, b'OES1X'): ('cplstn4', 'NA'),
            (273, 1, 26, b'OES1X1'): ('cplstn6', 'NA'),
            (273, 1, 26, b'OES1X'): ('cplstn6', 'NA'),
            (274, 1, 32, b'OES1X1'): ('cplstn3', 'NA'),
            (274, 1, 32, b'OES1X'): ('cplstn3', 'NA'),

            # 275-CPLSTS3
            # 277-CPLSTS6
            (275, 1, 6, b'OES1X1'): ('cplsts3', 'NA'),
            (276, 1, 32, b'OES1X1'): ('cplsts4', 'NA'),
            (277, 1, 26, b'OES1X1'): ('cplsts6', 'NA'),
            (278, 1, 32, b'OES1X1'): ('cplsts8', 'NA'),

            (1, 2, 5, 'OESVM1'): ('crod', 'NA'),
            (10, 2, 5, 'OESVM1'): ('conrod', 'NA'),
            (10, 2, 5, 'OES1X'): ('conrod', 'NA'),

            (11, 2, 3, 'OESVM1'): ('celas1', 'NA'),
            (12, 2, 3, 'OESVM1'): ('celas2', 'NA'),

            (2, 2, 111, b'OESVM1'): ('cbeam', 'NA'),
            (34, 2, 19, b'OESVM1'): ('cbar', 'NA'),

            (4, 2, 5, 'OESVM1'): ('cshear', 'NA'),
            (4, 2, 5, 'OES1X'): ('cshear', 'NA'),
            (74, 2, 17, 'OESVM1'): ('ctria3', 'NA'),
            (144, 2, 87, b'OESVM1'): ('cquad4', 'NA'),

            (95, 2, 13, b'OESVM1C'): ('cquad4', 'NA'),
            (95, 3, 13, b'OESVM1C'): ('cquad4', 'NA'),
            (97, 2, 13, b'OESVM1C'): ('ctria3', 'NA'),
            (97, 3, 13, b'OESVM1C'): ('ctria3', 'NA'),

            (39, 2, 74, 'OESVM1'): ('ctetra', 'NA'),
            (67, 2, 130, b'OESVM1'): ('chexa', 'NA'),
            (68, 2, 102, b'OESVM1'): ('cpenta', 'NA'),
        }
        op2 = self.op2
        key = (op2.element_type, op2.format_code, op2.num_wide, op2.table_name)
        try:
            return stress_mapper[key]
        except KeyError:  # pragma: no cover
            op2.log.error(op2.code_information())
            msg = ('stress_mapper (~line 850 of oes.py) does not contain the '
                   'following key and must be added\n'
                   'key=(element_type=%r, format_code=%r, num_wide=%r, table_name=%r) ' % key)
            op2.log.error(msg)
            #raise KeyError(msg)
            raise
            #return None, None

    def get_oes_prefix_postfix(self) -> tuple[str, str]:
        """
        Creates the prefix/postfix that splits off ATO, CRM, PSD, nonlinear,
        etc. results.  We also fix some of the sort bits as typing:

            STRESS(PLOT,SORT1,RALL) = ALL

        will actually create the OESRMS2 table (depending on what else
        is in your case control).  However, it's in an OESATO2 table, so
        we know it's really SORT2.

        Also, if you're validating the sort_bit flags, *RMS2 and *NO2 are
        actually SORT1 tables.

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

        if table_name_bytes in [b'OES1X1', b'OES1X',
                                b'OES1C', b'OES1']:
            prefix = 'stress.'
            self._set_as_sort1()
        elif table_name_bytes == b'OES1A':
            prefix = 'stressa.'
        elif table_name_bytes in [b'OSTR1X', b'OSTR1', b'OSTR1C']:
            prefix = 'strain.'
            self._set_as_sort1()
        elif table_name_bytes in [b'OES2', b'OSTR2', b'OES2C', b'OSTR2C']:
            self._set_as_sort2()
        #elif table_name_bytes in ['OESNLXR']:
            #prefix = 'sideline_'
        elif table_name_bytes in [b'DOES1', b'DOSTR1']:
            self._set_as_sort1()
            assert op2.thermal in {0, 2, 4, 8}, op2.code_information()
            prefix, postfix = get_shock_prefix_postfix(op2.thermal)
        elif table_name_bytes in [b'OESNLXD', b'OESNL1X', b'OESNLXR', b'OESNL2']:
            prefix = 'nonlinear_'
        elif table_name_bytes in [b'OESNLXR2']:
            prefix = 'nonlinear_'
        elif table_name_bytes == b'OESNLBR':
            prefix = 'sideline_'
        elif table_name_bytes == b'OESRT':
            #OESRT: Table of composite element strength ratios
            prefix = 'strength_ratio.'
        elif table_name_bytes == b'OESRTN':
            op2.stress_bits[1] = 0
            #OESRTN: Table of composite element strength ratios nonlinear???
            prefix = 'strength_ratio.'
        elif table_name_bytes in [b'OESCP', b'OESTRCP']:
            # guessing
            pass
            #op2.sort_bits[0] = 0 # real; ???
            #op2.sort_bits[1] = 0 # sort1
            #op2.sort_bits[2] = 1 # random; ???
        elif table_name_bytes in [b'OESVM1C', b'OSTRVM1C', b'OESVM1', b'OSTRVM1',
                                  #b'OESVM1C', b'OSTRVM1C',
                                  b'OESVM2', b'OSTRVM2',]:
            prefix = 'modal_contribution.'
            op2.to_nx(f' because table_name={table_name_bytes} was found')

        #----------------------------------------------------------------
        elif table_name_bytes in [b'OSTRMS1C']:  #, b'OSTRMS1C']:
            op2.format_code = 1
            op2.sort_bits[0] = 0  # real
            prefix = 'rms.'

        elif table_name_bytes in [b'OESXRMS1']:
            op2._analysis_code_fmt = b'i'
            self._set_as_random()
            self._set_as_sort1()
            prefix = 'rms.'
        elif table_name_bytes in [b'OESXRMS2']:  # wrong?
            self._set_as_random()
            self._set_as_sort2()
            prefix = 'rms.'

        elif table_name_bytes in [b'OESXNO1']:
            self._set_as_random()
            self._set_as_sort1()
            prefix = 'no.'
        elif table_name_bytes in [b'OESXNO1C']:
            # - ply-by-ply Stresses including:
            #    - von Mises Stress for PSDF (OESPSD1C),
            #    - Cumulative Root Mean Square output (OESXNO1C)
            #    - Positive Crossing (OESCRM1C) output sets
            # - ply-by-ply Strains for:
            #    - PSDF (OSTPSD1C)
            #    - Cumulative Root Mean Square (OSTCRM1C) output sets
            self._set_as_random()
            self._set_as_sort1()
            prefix = 'crm.'
        elif table_name_bytes in [b'OESXRM1C']:
            self._set_as_random()
            self._set_as_sort1()
            prefix = 'rms.'
            #print(op2.code_information())

        elif table_name_bytes in [b'OESRMS1', b'OSTRRMS1']:
            op2._analysis_code_fmt = b'i'
            self._set_as_random()
            self._set_as_sort1()
            prefix = 'rms.'
        elif table_name_bytes in [b'OESRMS2', b'OSTRRMS2']:
            op2._analysis_code_fmt = b'i'
            self._set_as_random()
            #print(op2.code_information())
            self._set_as_sort1()  # it's not really SORT2...
            op2.sort_method = 1
            if table_name_bytes == b'OESRMS2':
                op2.table_name = b'OESRMS1'
            elif table_name_bytes == b'OSTRRMS2':
                op2.table_name = b'OSTRRMS1'
            else:
                raise NotImplementedError(table_name_bytes)
            #assert op2.sort_method == 2, op2.code_information()
            prefix = 'rms.'

        elif table_name_bytes in [b'OESNO1', b'OSTRNO1', b'OSTNO1C']:
            assert op2.sort_method == 1, op2.code_information()
            self._set_as_random()
            prefix = 'no.'
        elif table_name_bytes in [b'OESNO2', b'OSTRNO2']:
            self._set_as_random()
            self._set_as_sort1()
            op2.data_code['nonlinear_factor'] = None
            op2._analysis_code_fmt = b'i'
            prefix = 'no.'
        #----------------------------------------------------------------

        elif table_name_bytes in [b'OESPSD1', b'OSTRPSD1']:
            #op2.format_code = 1
            op2.sort_bits[0] = 0  # real
            op2.sort_bits[1] = 0  # sort1
            op2.sort_bits[2] = 1  # random
            prefix = 'psd.'
        elif table_name_bytes in [b'OESPSD2', b'OSTRPSD2',
                                  b'OESPSD2C', b'OSTPSD2C']:
            if 0:
                # TODO: the sort bits might not be right...isat_random
                #print(op2.code_information())
                #print(op2.sort_bits)
                op2.format_code = 1
                #op2.sort_bits[0] = 0 # real
                #op2.sort_bits[1] = 1 # sort2
                #op2.sort_bits[2] = 1 # random
                op2.sort_bits.is_real = 1
                op2.sort_bits.is_sort2 = 1
                op2.sort_bits.is_random = 1
                #print(op2.code_information())
                #print(op2.sort_bits)
            else:
                op2.format_code = 1 # real
                op2.result_type = 2 # random
                op2.sort_bits[0] = 0 # real
                op2.sort_bits[2] = 1 # random
            prefix = 'psd.'

        elif table_name_bytes in [b'OESATO1', b'OSTRATO1']:
            prefix = 'ato.'
        elif table_name_bytes in [b'OESATO2', b'OSTRATO2']:
            prefix = 'ato.'

        elif table_name_bytes in [b'OESCRM1', b'OSTRCRM1']:
            prefix = 'crm.'
            op2.result_type = 2 # random
            op2.sort_bits[2] = 1 # random
        elif table_name_bytes in [b'OESCRM2', b'OSTRCRM2']:
            # sort2, random
            op2.format_code = 1 # real
            op2.result_type = 2 # random
            op2.sort_bits[0] = 0 # real
            op2.sort_bits[1] = 1 # sort2
            op2.sort_bits[2] = 1 # random
            op2.sort_method = 2
            prefix = 'crm.'
        #elif op2.table_name in ['DOES1', 'DOSTR1']:
            #prefix = 'scaled_response_spectra_'
        #elif op2.table_name in ['OESCP']:

        elif table_name_bytes in [b'RASCONS']: #, b'OSTRMS1C']:
            op2.format_code = 1
            op2.sort_bits[0] = 0 # real
            prefix = 'RASCONS.'
        elif table_name_bytes in [b'RAECONS']: #, b'OSTRMS1C']:
            op2.format_code = 1
            op2.sort_bits[0] = 0 # real
            prefix = 'RAECONS.'
        elif table_name_bytes in [b'RAPCONS']: #, b'OSTRMS1C']:
            op2.format_code = 1
            op2.sort_bits[0] = 0 # real
            prefix = 'RAPCONS.'

        elif table_name_bytes in [b'RASEATC']: #, b'OSTRMS1C']:
            self._set_as_real()
            prefix = 'RASEATC.'
        elif table_name_bytes in [b'RAEEATC']: #, b'OSTRMS1C']:
            self._set_as_real()
            prefix = 'RAEEATC.'
        elif table_name_bytes in [b'RAPEATC']: #, b'OSTRMS1C']:
            self._set_as_real()
            prefix = 'RAPEATC.'
        elif table_name_bytes in [b'OESMC1', b'OSTRMC1']:
            prefix = 'modal_contribution.'
        elif table_name_bytes in [b'OESC1']:
            # NASA95
            prefix = 'stress.'
        else:
            raise NotImplementedError(op2.table_name)

        #if op2.analysis_code == 1:
            #op2.sort_bits[1] = 0 # sort1
            #op2.sort_method = 1

        op2.data_code['sort_bits'] = op2.sort_bits
        op2.data_code['nonlinear_factor'] = op2.nonlinear_factor
        return prefix, postfix

    def _set_as_real(self):
        op2 = self.op2
        op2.format_code = 1
        op2.result_type = 0
        op2.sort_bits[0] = 0 # real
        op2.sort_bits.is_real = True
        op2.sort_bits.is_random = False

    def _set_as_random(self):
        op2 = self.op2
        op2.format_code = 1  # real
        op2.result_type = 2  # random
        op2.sort_bits.is_real = True
        op2.sort_bits.is_random = True

    def _set_as_sort1(self):
        op2 = self.op2
        op2.sort_bits[1] = 0 # sort1
        op2.sort_method = 1

    def _set_as_sort2(self):
        op2 = self.op2
        op2.sort_bits[1] = 1 # sort2
        op2.sort_method = 2

    def _read_oesmc_4(self, data: bytes, ndata: int) -> int:
        op2 = self.op2
        n = 0
        log = op2.log
        if op2.element_type == 1:
            assert op2.num_wide == 4, op2.code_information()
            if op2.read_mode == 1:
                return ndata
            ntotal = 16 * self.factor # 4*4
            nelements = ndata // ntotal
            fmt = mapfmt(op2._endian + b'i3f', self.size)
            struct1 = Struct(fmt)
            for ielem in range(nelements):
                edata = data[n:n+ntotal]
                out = struct1.unpack(edata)
                #print(out)
                n += ntotal
            log.warning(f'skipping {op2.table_name} with {op2.element_name}-{op2.element_type}')
        else:
            raise NotImplementedError(op2.code_information())
        return n

    def _read_oes1_loads_nasa95(self, data, ndata: int) -> tuple[int, Any, Any]:
        """Reads OES1 subtable 4 for NASA 95"""
        op2 = self.op2
        prefix, postfix = self.get_oes_prefix_postfix()
        result_type = op2.result_type
        #self._apply_oes_ato_crm_psd_rms_no('') # TODO: just testing
        n = 0
        is_magnitude_phase = op2.is_magnitude_phase()
        dt = op2.nonlinear_factor

        if op2.is_stress:
            result_name = 'stress'
            #stress_name = 'STRESS'
        else:
            result_name = 'strain'
            #stress_name = 'STRAIN'

        if op2._results.is_not_saved(result_name):
            return ndata
        if op2.element_type in [1, 3, 10]:  # rods
            # 1-CROD
            # 3-CTUBE
            # 10-CONROD
            n, nelements, ntotal = self._oes_crod(data, ndata, dt, is_magnitude_phase,
                                                  result_type, prefix, postfix)

        #elif op2.element_type == 2: # CBEAM
            #n, nelements, ntotal = self._oes_cbeam(data, ndata, dt, is_magnitude_phase,
                                                   #result_type, prefix, postfix)

        elif op2.element_type == 4: # CSHEAR
            n, nelements, ntotal = self._oes_cshear(data, ndata, dt, is_magnitude_phase,
                                                    result_type, prefix, postfix)

        elif op2.element_type in {11, 12, 13, 14}:  # springs
            # 11-CELAS1
            # 12-CELAS2
            # 13-CELAS3
            # 14-CELAS4
            n, nelements, ntotal = self._oes_celas(data, ndata, dt, is_magnitude_phase,
                                                   result_type, prefix, postfix)

        elif op2.element_type == 19:
            # 19-CQUAD1
            n, nelements, ntotal = self._oes_cquad4_33(data, ndata, dt, is_magnitude_phase,
                                                       result_type, prefix, postfix)
        elif op2.element_type == 34: # CBAR
            n, nelements, ntotal = self._oes_cbar_34(data, ndata, dt, is_magnitude_phase,
                                                     result_type, prefix, postfix)

        elif op2.element_type == 64: # CQUAD4
            op2.element_type = 33  # faking...
            n, nelements, ntotal = self._oes_cquad4_33(data, ndata, dt, is_magnitude_phase,
                                                       result_type, prefix, postfix)

        elif op2.element_type == 83:
            # 83: TRIA3
            n, nelements, ntotal = self._oes_ctria3(data, ndata, dt, is_magnitude_phase,
                                                    result_type, prefix, postfix)

        #elif op2.element_type in [64, 70, 75, 82, 144]:  # bilinear plates
            # 64-CQUAD8
            # 70-CTRIAR
            # 75-CTRIA6
            # 82-CQUADR
            # 144-CQUAD4-bilinear
            #n, nelements, ntotal = self._oes_cquad4_144(data, ndata, dt, is_magnitude_phase,
                                                        #result_type, prefix, postfix)
        else:
            #msg = 'sort1 Type=%s num=%s' % (op2.element_name, op2.element_type)
            msg = op2.code_information()
            print(msg)
            return op2._not_implemented_or_skip(data, ndata, msg)

        if nelements is None:
            return n
        assert ndata > 0, ndata
        assert nelements > 0, 'nelements=%r element_type=%s element_name=%r' % (nelements, op2.element_type, op2.element_name)
        #assert ndata % ntotal == 0, '%s n=%s nwide=%s len=%s ntotal=%s' % (op2.element_name, ndata % ntotal, ndata % op2.num_wide, ndata, ntotal)
        assert op2.num_wide * 4 == ntotal, 'numwide*4=%s ntotal=%s' % (op2.num_wide * 4, ntotal)
        assert op2.thermal == 0, "thermal = %%s" % op2.thermal
        assert n > 0, f'n = {n} result_name={result_name}'
        return n

    def _read_oes1_loads(self, data, ndata: int):
        """Reads OES op2.thermal=0 stress/strain"""
        op2 = self.op2
        log = op2.log
        prefix, postfix = self.get_oes_prefix_postfix()
        result_type = op2.result_type

        #self._apply_oes_ato_crm_psd_rms_no('') # TODO: just testing
        n = 0
        is_magnitude_phase = op2.is_magnitude_phase()
        dt = op2.nonlinear_factor

        #flag = 'element_id'
        if op2.is_stress:
            result_name = 'stress'
            stress_name = 'STRESS'
        else:
            result_name = 'strain'
            stress_name = 'STRAIN'

        #if op2.is_stress:
            #_result_name, _class_obj = self.get_stress_mapper()
        if op2.table_name_str == 'OESXRMS1':
            assert op2.sort_method == 1, op2.code_information()

        if op2._results.is_not_saved(result_name):
            return ndata
        if op2.element_type in [1, 3, 10]:  # rods
            # 1-CROD
            # 3-CTUBE
            # 10-CONROD
            n, nelements, ntotal = self._oes_crod(data, ndata, dt, is_magnitude_phase,
                                                  result_type, prefix, postfix)

        elif op2.element_type == 2: # CBEAM
            n, nelements, ntotal = self._oes_cbeam(data, ndata, dt, is_magnitude_phase,
                                                   result_type, prefix, postfix)

        elif op2.element_type == 4: # CSHEAR
            n, nelements, ntotal = self._oes_cshear(data, ndata, dt, is_magnitude_phase,
                                                    result_type, prefix, postfix)

        elif op2.element_type in [11, 12, 13, 14]:  # springs
            # 11-CELAS1
            # 12-CELAS2
            # 13-CELAS3
            # 14-CELAS4
            n, nelements, ntotal = self._oes_celas(data, ndata, dt, is_magnitude_phase,
                                                   result_type, prefix, postfix)

        elif op2.element_type == 34: # CBAR
            n, nelements, ntotal = self._oes_cbar_34(data, ndata, dt, is_magnitude_phase,
                                                     result_type, prefix, postfix)

        elif op2.element_type in [39, 67, 68, 255]: # solid stress
            # 39-CTETRA
            # 67-CHEXA
            # 68-CPENTA
            # 255-CPYRAM
            n, nelements, ntotal = self._oes_csolid(data, ndata, dt, is_magnitude_phase,
                                                    result_type, prefix, postfix)
        elif op2.element_type in [140]:
            # 144-CHEXAFD
            #TestOP2.test_bdf_op2_other_23
            '             S T R E S S E S   I N   H Y P E R E L A S T I C   H E X A H E D R O N   E L E M E N T S  ( H E X A F D ) '
            '              GRID/    POINT               --------CAUCHY STRESSES---------              DIR.  COSINES       MEAN'
            '  ELEMENT-ID  GAUSS      ID         NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE'
            '0      211    GAUS'
            '                           1  X   2.627481E+02  XY   8.335709E+01   A   3.040629E+02  LX 0.94-0.34-0.00  -8.630411E+01'
            '                              Y  -2.599833E+00  YZ  -1.138583E+01   B  -5.463896E+01  LY 0.25 0.68 0.69'
            '                              Z  -1.235891E+00  ZX   7.851368E+01   C   9.488372E+00  LZ 0.23 0.64-0.73'
            '                           2  X   1.541617E+02  XY   4.154493E+01   A   1.964021E+02  LX 0.88-0.47-0.00  -8.630411E+01'
            '                              Y   5.412691E+01  YZ  -3.499344E+00   B   6.221669E+00  LY 0.25 0.46-0.85'
            '                              Z   5.062376E+01  ZX   6.725376E+01   C   5.628857E+01  LZ 0.40 0.75 0.53'
            n, nelements, ntotal = self._oes_csolid_linear_hyperelastic_cosine(data, ndata, dt, is_magnitude_phase,
                                                                               result_type, prefix, postfix)

        #elif op2.element_type in [160, 163, 166,
                                   #161, # centroid
                                   #165,]:
            # nonlinear hyperelastic solids
            # 160-CPENTAFD
            # 163-CHEXAFD
            # 166-CTETRAFD

            # centroid??
            # 161-CTETRAFD

            # many nodes?
            # 165-CPENTAFD
            #n, nelements, ntotal = self._oes_csolid_linear_hyperelastic(data, ndata, dt, is_magnitude_phase,
                                                                        #result_type, prefix, postfix)
        elif op2.element_type in [202, 204,
                                   216, 218, 220, 221]:
            # nonlinear hyperelastic solids
            # 202-CHEXAFD
            # 204-CPENTAFD

            # also nonlinear hyperelastic solid, but somewhat different
            '   N O N L I N E A R   S T R E S S E S   I N   H Y P E R E L A S T I C   H E X A H E D R O N   E L E M E N T S  ( HEXA8FD )'
            ' '
            '  ELEMENT GRID/   POINT                       CAUCHY STRESSES/ LOG STRAINS                        PRESSURE    VOL. STRAIN'
            '     ID   GAUSS     ID       X           Y           Z           XY          YZ          ZX'
            '0     401 GRID      401  1.9128E+03  6.2729E+02 -3.4828E+02 -7.5176E+01  7.8259E+00 -2.5001E+02  7.3060E+02  7.3060E-03'
            '                         6.8270E-01 -6.5437E-04 -1.2874E+00 -3.9645E-02 -2.9882E-03 -5.9975E-02'
            # 216-TETRAFD
            # 218-HEXAFD
            # 220-PENTAFD
            # 221-TETRAFD
            n, nelements, ntotal = self._oes_csolid_nonlinear_hyperelastic(data, ndata, dt, is_magnitude_phase,
                                                                           result_type, prefix, postfix)

        elif op2.element_type in [300, 301, 302, 303]: # solid stress
            # solids without stress eigenvectors
            # 300-CHEXA
            # 301-CPENTA
            # 302-CTETRA
            # 303-CPYRAM
            n, nelements, ntotal = self._oes_csolid2(data, ndata, dt, is_magnitude_phase,
                                                     result_type, prefix, postfix)
        elif op2.element_type in [306, 307]:
            # 306-CHEXALN
            # 307-CPENTALN
            n, nelements, ntotal = self._oes_csolid_composite(data, ndata, dt, is_magnitude_phase,
                                                              result_type, prefix, postfix)
        #=========================
        # plates
        elif op2.element_type in [33, 228]:
            # 33: CQUAD4-centroidal
            # 228: CQUADR-centroidal
            n, nelements, ntotal = self._oes_cquad4_33(data, ndata, dt, is_magnitude_phase,
                                                       result_type, prefix, postfix)

        elif op2.element_type in [74, 227]: # 229???
            # 74: TRIA3
            # 227: TRIAR
            n, nelements, ntotal = self._oes_ctria3(data, ndata, dt, is_magnitude_phase,
                                                    result_type, prefix, postfix)

        elif op2.element_type in [64, 70, 75, 82, 144]:  # bilinear plates
            # 64-CQUAD8
            # 70-CTRIAR
            # 75-CTRIA6
            # 82-CQUADR
            # 144-CQUAD4-bilinear
            n, nelements, ntotal = self._oes_cquad4_144(data, ndata, dt, is_magnitude_phase,
                                                        result_type, prefix, postfix)

        elif op2.element_type in [88, 90]: # nonlinear shells
            # 88-CTRIA3NL
            # 90-CQUAD4NL
            n, nelements, ntotal = self._oes_shells_nonlinear(data, ndata, dt, is_magnitude_phase,
                                                              result_type, prefix, postfix)

        elif op2.element_type in [95, 96, 97, 98, 232, 233]: # composite shell
            # 95 - CQUAD4
            # 96 - CQUAD8
            # 97 - CTRIA3
            # 98 - CTRIA6 (composite)
            # 232 - QUADRLC (CQUADR-composite)
            # 233 - TRIARLC (CTRIAR-composite)
            n, nelements, ntotal = self._oes_shells_composite(data, ndata, dt, is_magnitude_phase,
                                                              result_type, prefix, postfix)

        elif op2.element_type == 53: # axial plates - ctriax6
            n, nelements, ntotal = self._oes_ctriax6(data, ndata, dt, is_magnitude_phase,
                                                     result_type, prefix, postfix)

        elif op2.element_type == 102: # cbush
            n, nelements, ntotal = self._oes_cbush(data, ndata, dt, is_magnitude_phase,
                                                   result_type, prefix, postfix)

        elif op2.element_type == 40:  # cbush1d
            n, nelements, ntotal = self._oes_cbush1d(data, ndata, dt, is_magnitude_phase,
                                                     result_type, prefix, postfix)

        elif op2.element_type in [87, 89, 92]:  # nonlinear rods
            # 87-CTUBENL
            # 89-RODNL
            # 92-CONRODNL
            n, nelements, ntotal = self._oes_crod_nonlinear(data, ndata, dt, is_magnitude_phase,
                                                            result_type, prefix, postfix)

        elif op2.element_type in [224, 225]: # nonlinear spring
            # 224-CELAS1
            # 225-CELAS3
            # NonlinearSpringStress
            n, nelements, ntotal = self._oes_celas_nonlinear(data, ndata, dt, is_magnitude_phase,
                                                             result_type, prefix, postfix)

        elif op2.element_type == 69:  # cbend
            # 69-CBEND
            n, nelements, ntotal = self._oes_cbend(data, ndata, dt, is_magnitude_phase,
                                                   result_type, prefix, postfix)

        elif op2.element_type == 86:  # cgap
            # 86-GAPNL
            n, nelements, ntotal = self._oes_cgap_nonlinear(data, ndata, dt, is_magnitude_phase,
                                                            result_type, prefix, postfix)

        elif op2.element_type == 94:
            # 94-BEAMNL
            n, nelements, ntotal = self._oes_cbeam_nonlinear(data, ndata, dt, is_magnitude_phase,
                                                             result_type, prefix, postfix)

        elif op2.element_type in [85, 91, 93, 256]:
            # 256-PYRAM
            n, nelements, ntotal = self._oes_csolid_nonlinear(data, ndata, dt, is_magnitude_phase,
                                                              result_type, prefix, postfix)

        elif op2.element_type == 100:  # bars
            # 100-BARS
            n, nelements, ntotal = self._oes_cbar_100(data, ndata, dt, is_magnitude_phase,
                                                      result_type, prefix, postfix)

        #-----------------------------------------------------------------------

        elif op2.element_type == 139:
            n, nelements, ntotal = self._oes_hyperelastic_quad(data, ndata, dt, is_magnitude_phase,
                                                               result_type, prefix, postfix)

        elif op2.element_type == 226:
            # 226-BUSHNL
            n, nelements, ntotal = self._oes_cbush_nonlinear(data, ndata, dt, is_magnitude_phase,
                                                             result_type, prefix, postfix)
        elif op2.element_type in [275]: # 271,
            #271 CPLSTN3 Triangle plane strain linear format (Center Only)
            #272 CPLSTN4 Quadrilateral plane strain linear format (Center and Corners)
            #273 CPLSTN6 Triangle plane strain linear format (Center and Corners)
            #274 CPLSTN8 Quadrilateral plane strain linear format (Center and Corners)

            #275 CPLSTS3 Triangle plane stress linear Format (Center Only)
            #276 CPLSTS4 Quadrilateral plane stress linear format (Center and Corners)
            #277 CPLSTS6 Triangle plane stress linear format (Center and Corners)
            #278 CPLSTS8 Quadrilateral plane stress linear format (Center and Corners)

            # 271-CPLSTN3
            # 275-CPLSTS3
            n, nelements, ntotal = self._oes_plate_stress_34(data, ndata, dt, is_magnitude_phase,
                                                             stress_name, prefix, postfix)

        elif op2.element_type in [276, 277, 278]:
            # 276-CPLSTS4
            # 277-CPLSTS6
            # 278-CPLSTS8
            n, nelements, ntotal = self._oes_plate_stress_68(data, ndata, dt, is_magnitude_phase,
                                                             stress_name, prefix, postfix)


        elif op2.element_type == 35: # CON
            return ndata

        elif op2.element_type in [60, 61]:
            # 60-DUM8
            # 61-DUM9
            return ndata

        elif op2.element_type == 101: # AABSF
            return ndata

        elif op2.element_type in [47, 48, 189, 190]:
            # 47-AXIF2
            # 48-AXIF3
            # 189-???
            # 190-VUTRIA
            return ndata
        elif op2.element_type in [50, 51, 203]:
            # 203-SLIF1D?
            # 50-SLOT3
            # 51-SLOT4
            return ndata
        elif op2.element_type in [162, 164, 167, 168,
                                   169, 170, 171, 172,
                                   218, 211, 213, 214,
                                   217, 219, 222, 223,
                                   232, 235]:
            # 162-TRIA3FD
            # 164-QUADFD
            # 167-TRIAFD
            # 168-TRIAX3FD
            # 169-TRIAXFD
            # 170-QUADX4FD
            # 171-QUADXFD
            # 172-QUADRNL
            # 211-TRIAFD
            # 213-TRIAXFD
            # 214-QUADX4FD
            # 217-TRIA3FD
            # 219-QUADFD
            # 223-QUADXFD
            # 222-TRIAX3FD
            # 232-QUADRLC
            # 235-CQUADR
            return op2._not_implemented_or_skip(data, ndata, op2.code_information())
        elif op2.element_type in [145, 146, 147, # VU-solid
                                   189,  # VUQUAD
                                   191]: # VUBEAM
            msg = f'{op2.element_name}-{op2.element_type} has been removed'
            return op2._not_implemented_or_skip(data, ndata, msg)
        elif op2.element_type == 118:  # WELDP-MSC
            #'                                S T R E S S E S   I N   W E L D   E L E M E N T S   ( C W E L D P ) '
            #' '
            #'    ELEMENT          AXIAL         MAX  STRESS      MIN  STRESS      MAX  STRESS      MIN  STRESS        MAXIMUM          BEARING '
            #'      ID             STRESS           END-A            END-A            END-B            END-B        SHEAR  STRESS       STRESS'
            #'        179      -3.153108E+00     8.089753E+02    -8.152815E+02     7.946552E+02    -8.009614E+02     2.852777E+01     1.179798E+01'

            # ELEMENT-ID =     100
            #     S T R A I N S   I N   W E L D   E L E M E N T S   ( C W E L D P )
            #
            #     AXIAL         MAX  STRAIN      MIN  STRAIN      MAX  STRAIN      MIN  STRAIN        MAXIMUM
            #     TIME          STRAIN           END-A            END-A            END-B            END-B        SHEAR  STRAIN
            # 0.0           0.0              0.0              0.0              0.0              0.0              0.0              0.0
            # 1.000000E-01  0.0              0.0              0.0              0.0              0.0              0.0              0.0
            # 2.000000E-01  1.652614E-02     2.381662E+02    -2.381332E+02     2.381623E+02    -2.381293E+02     5.678050E+01     0.0
            # 3.000000E-01  6.468190E-03     4.706443E+01    -4.705150E+01     4.703462E+01    -4.702168E+01     1.121626E+01     0.0

            #ints    = (1001, -0.0007072892040014267, 0.6948937773704529, -0.6963083744049072, 0.6948915123939514, -0.6963061094284058, 6.161498617984762e-07, 0)
            #floats  = (1001, -0.0007072892040014267, 0.6948937773704529, -0.6963083744049072, 0.6948915123939514, -0.6963061094284058, 6.161498617984762e-07, 0.0)
            #if data:
                #self.show_data(data)
            n, nelements, ntotal = self._oes_weld_118(
                data, ndata, dt, is_magnitude_phase,
                result_type, prefix, postfix)
        elif op2.element_type == 126:  # FASTP
            #C:\MSC.Software\msc_nastran_runs\cf103e.op2
            # S T R E S S E S   I N   F A S T E N E R   E L E M E N T S   ( C F A S T )
            #
            # ELEMENT-ID         FORCE-X          FORCE-Y          FORCE-Z         MOMENT-X         MOMENT-Y         MOMENT-Z
            # data  = (301, -4.547473508864641e-09, 1.8571810755929619e-09, -7.94031507211912e-10, -0.0, -0.0, 0.0,
            #         401, -4.547473508864641e-09, -2.0263790645458357e-09, 1.1617373729677638e-09, -0.0, 0.0, 0.0)
            n, nelements, ntotal = self._oes_fast_126(
                data, ndata, dt, is_magnitude_phase,
                result_type, prefix, postfix)
        elif op2.is_nx and op2.element_type in [269, 270]:
            # 269-CHEXAL
            # 270-PENTAL
            n, nelements, ntotal = self._oes_composite_solid_nx(
                data, ndata, dt, is_magnitude_phase,
                result_type, prefix, postfix)


        elif op2.element_type in [271, 272, 273, 274]:
            # 271 CPLSTN3
            # 272 CPLSTN4
            # 273 CPLSTN6
            # 274 CPLSTN3 -> 8???
            n, nelements, ntotal = self._oes_cplstn_nx(
                data, ndata, dt, is_magnitude_phase,
                result_type, prefix, postfix)
        elif op2.element_type in [159, 184,
                                  200, 201, 236, 237, 242, 243, 244, 245]:
            # 159-SEAMP
            # 184-CBEAM3
            #
            # 200-WELD
            # 201 CQUAD4FD
            # 236 CTRIAR-corner
            # 237 CTRIAR-center
            # 242-CHEXA?
            # 243 CQUADX4
            # 244 CTRAX6
            # 245 CQUADX8
            log.warning(f'skipping {op2.element_name}-{op2.element_type}')
            return ndata
        elif op2.element_type in [312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323,
                                   343, 344, 345, 346, 347, 348, 349,
                                   350, 351, 352, 355, 356, 357, 358, 363]:
            #
            # 312 TRAX3
            # 313 QUADX4
            # 314 TRAX6
            # 315 QUADX8
            # 316 PLSTN3
            # 317 PLSTN4
            # 318 PLSTN6
            # 319 PLSTN8
            # 320 PLSTS3
            # 321 PLSTS4
            # 322 PLSTS6
            # 323 PLSTS8
            #
            # 343 CTRIA6 SOL 401
            # 344 CQUAD8 SOL 401
            # 345 CTRIAR SOL 401
            # 346 CQUADR SOL 401
            # 347 CBAR SOL 401
            # 348 CBEAM SOL 401
            # 349 CBUSH1D SOL 401
            #
            # 350 CELAS1 SOL 401
            # 351 CELAS2 SOL 401
            # 352 CBUSH SOL 401
            # 355 Composite triangular shell element (CTRIA6); SOL 402?
            # 356 Composite quadrilateral shell element (CQUAD8); SOL 402?
            # 357 Composite triangular shell element (CTRIAR); SOL 402?
            # 358 Composite quadrilateral shell element (CQUADR); SOL 402?
            # 363 CROD SOL 402
            log.warning(f'skipping {op2.element_name}-{op2.element_type}')
            return ndata
        else:
            #msg = 'sort1 Type=%s num=%s' % (op2.element_name, op2.element_type)
            msg = op2.code_information()
            #raise NotImplementedError(msg)
            return op2._not_implemented_or_skip(data, ndata, msg)

        try:
            nelements
        except NameError:
            raise RuntimeError(op2.code_information())

        if nelements is None:
            return n

        #self.check_element_ids()
        assert ndata > 0, ndata
        assert nelements > 0, f'nelements={nelements} element_type={op2.element_type} element_name={op2.element_name!r}'
        #assert ndata % ntotal == 0, '%s n=%s nwide=%s len=%s ntotal=%s' % (op2.element_name, ndata % ntotal, ndata % op2.num_wide, ndata, ntotal)
        assert op2.num_wide * 4 * self.factor == ntotal, f'numwide*4={op2.num_wide*4} ntotal={ntotal} element_name={op2.element_name!r}\n{op2.code_information()}'
        #assert op2.thermal == 0, "thermal = %%s" % op2.thermal
        assert n is not None and n > 0, f'n={n} result_name={result_name}\n{op2.code_information()}'

        #if self.is_sort2:
            #assert len(np.unique(op2.obj._times)) == len(op2.obj._times), f'{op2.obj._times.tolist()}\n{op2.code_information()}'
        return n

    def check_element_ids(self):
        op2 = self.op2
        if op2.read_mode == 1:
            return
        if op2.is_sort1:
            obj = op2.obj
            if obj is None:
                raise RuntimeError('obj is None...\n' + op2.code_information())
            if hasattr(obj, 'element_node'):
                eids = obj.element_node[:, 0]
            elif hasattr(obj, 'element_layer'):
                eids = obj.element_layer[:, 0]
            elif hasattr(obj, 'element'):
                eids = obj.element
            else:
                print(op2.code_information())
                raise RuntimeError(''.join(obj.get_stats()))
            if eids.min() <= 0:
                #print(obj.code_information())
                print(''.join(obj.get_stats()))
                raise RuntimeError(f'{op2.element_name}-{op2.element_type}: {eids}')
        #else:
            #assert._times

    def _create_nodes_object(self, nnodes, result_name, slot, obj_vector):
        """same as _create_oes_object4 except it adds to the nnodes parameter"""
        op2 = self.op2
        auto_return = False
        #is_vectorized = True
        is_vectorized = op2._is_vectorized(obj_vector)
        #print("vectorized...read_mode=%s...%s; %s" % (op2.read_mode, result_name, is_vectorized))

        if is_vectorized:
            if op2.read_mode == 1:
                #print('oes-op2.nonlinear_factor =', op2.nonlinear_factor)
                #print(op2.data_code)
                op2.create_transient_object(result_name, slot, obj_vector)
                #print("read_mode 1; ntimes=%s" % op2.obj.ntimes)
                op2.result_names.add(result_name)
                #print('op2.obj =', op2.obj)
                op2.obj.nnodes += nnodes
                auto_return = True
            elif op2.read_mode == 2:
                self.code = op2._get_code()
                #op2.log.info("code = %s" % str(self.code))
                #print("code = %s" % str(self.code))

                # if this is failing, you probably set obj_vector to None...
                try:
                    op2.obj = slot[self.code]
                except KeyError:
                    msg = 'Could not find key=%s in result=%r\n' % (self.code, result_name)
                    msg += "There's probably an extra check for read_mode=1..."
                    self.op2.log.error(msg)
                    raise
                #op2.obj.update_data_code(op2.data_code)
                build_obj(op2.obj)

            else:  # not vectorized
                auto_return = True
        else:
            auto_return = True
        return auto_return, is_vectorized

    def _create_ntotal_object(self, ntotal, result_name, slot, obj_vector):
        """same as _create_oes_object4 except it adds to the ntotal parameter"""
        op2 = self.op2
        auto_return = False
        #is_vectorized = True
        is_vectorized = op2._is_vectorized(obj_vector)
        #print("vectorized...read_mode=%s...%s; %s" % (op2.read_mode, result_name, is_vectorized))

        if is_vectorized:
            if op2.read_mode == 1:
                #print('oes-op2.nonlinear_factor =', op2.nonlinear_factor)
                #print(op2.data_code)
                op2.create_transient_object(result_name, slot, obj_vector)
                #print("read_mode 1; ntimes=%s" % op2.obj.ntimes)
                op2.result_names.add(result_name)
                #print('op2.obj =', op2.obj)
                op2.obj.ntotal += ntotal
                auto_return = True
            elif op2.read_mode == 2:
                self.code = op2._get_code()
                #op2.log.info("code = %s" % str(self.code))
                #print("code = %s" % str(self.code))

                # if this is failing, you probably set obj_vector to None...
                try:
                    op2.obj = slot[self.code]
                except KeyError:
                    msg = 'Could not find key=%s in result=%r\n' % (self.code, result_name)
                    msg += "There's probably an extra check for read_mode=1..."
                    op2.log.error(msg)
                    raise
                #op2.obj.update_data_code(op2.data_code)
                build_obj(op2.obj)

            else:  # not vectorized
                auto_return = True
        else:
            auto_return = True
        return auto_return, is_vectorized

    def _oes_celas(self, data, ndata, dt, is_magnitude_phase,
                   result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 11 : CELAS1
         - 12 : CELAS2
         - 13 : CELAS3
         - 14 : CELAS4
        """
        op2 = self.op2
        n = 0
        if op2.is_stress:
            if prefix == '' and postfix == '':
                prefix = 'stress.'
            obj_real = RealSpringStressArray
            obj_complex = ComplexSpringStressArray
            if op2.element_type == 11:
                result_name = prefix + 'celas1_stress' + postfix
            elif op2.element_type == 12:
                result_name = prefix + 'celas2_stress' + postfix
            elif op2.element_type == 13:
                result_name = prefix + 'celas3_stress' + postfix
            elif op2.element_type == 14:
                result_name = prefix + 'celas4_stress' + postfix
            else:
                raise RuntimeError(op2.element_type)
        else:
            if prefix == '' and postfix == '':
                prefix = 'strain.'
            obj_real = RealSpringStrainArray
            obj_complex = ComplexSpringStrainArray
            if op2.element_type == 11:
                result_name = prefix + 'celas1_strain' + postfix
            elif op2.element_type == 12:
                result_name = prefix + 'celas2_strain' + postfix
            elif op2.element_type == 13:
                result_name = prefix + 'celas3_strain' + postfix
            elif op2.element_type == 14:
                result_name = prefix + 'celas4_strain' + postfix
            else:
                raise RuntimeError(op2.element_type)

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        log = op2.log
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        if op2.format_code == 1 and op2.num_wide == 2:  # real
            ntotal = 8 * self.factor # 2 * 4
            nelements = ndata // ntotal
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_real)
            if auto_return:
                #assert ntotal == op2.num_wide * 4
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * 4 * op2.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 2)
                obj._times[obj.itime] = dt
                self.obj_set_element(obj, itotal, itotal2, data, nelements)

                #(eid_device, stress)
                obj.data[obj.itime, itotal:itotal2, 0] = floats[:, 1].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    log.debug('vectorize CELASx real SORT%s' % op2.sort_method)
                n = oes_celas_real_2(op2, data, obj, nelements, ntotal, dt)

        elif op2.format_code in [2, 3] and op2.num_wide == 3:  # imag
            ntotal = 12 * self.factor
            nelements = ndata // ntotal

            nelements = ndata // ntotal
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_complex)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            assert obj is not None, op2.code_information()
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 3).copy()
                obj._times[obj.itime] = dt
                self.obj_set_element(obj, itotal, itotal2, data, nelements)

                if is_magnitude_phase:
                    mag = floats[:, 1]
                    phase = floats[:, 2]
                    rtheta = radians(phase)
                    real_imag = mag * (cos(rtheta) + 1.j * sin(rtheta))
                else:
                    real = floats[:, 1]
                    imag = floats[:, 2]
                    real_imag = real + 1.j * imag
                obj.data[obj.itime, itotal:itotal2, 0] = real_imag

                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    log.debug('vectorize CELASx imag SORT%s' % op2.sort_method)

                n = oes_celas_complex_3(op2, data, obj,
                                        nelements, ntotal,
                                        dt, is_magnitude_phase)
        elif op2.format_code == 1 and op2.num_wide == 3: # random
            raise RuntimeError(op2.code_information())
            #msg = op2.code_information()
            #return op2._not_implemented_or_skip(data, ndata, msg)
        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
        return n, nelements, ntotal

    def _oes_crod(self, data, ndata, dt, is_magnitude_phase,
                  result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 1 : CROD
         - 3 : CTUBE
         - 10 : CONROD

        """
        op2 = self.op2
        n = 0
        if op2.is_stress:
            obj_vector_real = RealRodStressArray
            obj_vector_complex = ComplexRodStressArray
            obj_vector_random = RandomRodStressArray
            if op2.element_type == 1: # CROD
                result_name = prefix + 'crod_stress' + postfix
            elif op2.element_type == 3:  # CTUBE
                result_name = prefix + 'ctube_stress' + postfix
            elif op2.element_type == 10:  # CONROD
                result_name = prefix + 'conrod_stress' + postfix
            else:  # pragma: no cover
                msg = op2.code_information()
                return op2._not_implemented_or_skip(data, ndata, msg)
        else:
            obj_vector_real = RealRodStrainArray
            obj_vector_complex = ComplexRodStrainArray
            obj_vector_random = RandomRodStrainArray
            if op2.element_type == 1: # CROD
                result_name = prefix + 'crod_strain' + postfix
            elif op2.element_type == 3:  # CTUBE
                result_name = prefix + 'ctube_strain' + postfix
            elif op2.element_type == 10:  # CONROD
                result_name = prefix + 'conrod_strain' + postfix
            else:  # pragma: no cover
                msg = op2.code_information()
                return op2._not_implemented_or_skip(data, ndata, msg)

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)

        #result_name, unused_is_random = self._apply_oes_ato_crm_psd_rms_no(result_name)
        slot = op2.get_result(result_name)
        if result_type == 0 and op2.num_wide == 5:  # real
            ntotal = 5 * 4 * self.factor
            nelements = ndata // ntotal

            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                #assert ntotal == op2.num_wide * 4
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * 4 * op2.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 5)
                obj._times[obj.itime] = dt
                self.obj_set_element(obj, itotal, itotal2, data, nelements)

                #[axial, torsion, SMa, SMt]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    op2.log.debug('vectorize CROD real SORT%s' % op2.sort_method)
                if op2.is_debug_file:
                    op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    op2.binary_debug.write('  #elementi = [eid_device, axial, axial_margin, torsion, torsion_margin]\n')
                    op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)
                n = oes_crod_real_5(op2, data, obj, nelements, ntotal, dt)

        elif result_type == 1 and op2.num_wide == 5: # imag
            ntotal = 20 * self.factor
            nelements = ndata // ntotal
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_complex)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 5)
                obj._times[obj.itime] = dt
                self.obj_set_element(obj, itotal, itotal2, data, nelements)

                real_imag = apply_mag_phase(floats, is_magnitude_phase, [1, 3], [2, 4])
                obj.data[obj.itime, itotal:itotal2, :] = real_imag

                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    op2.log.debug('vectorize CROD imag SORT%s' % op2.sort_method)
                n = oes_crod_complex_5(op2, data, obj, nelements, ntotal, dt, is_magnitude_phase)

        #elif op2.format_code in [2, 3] and op2.num_wide == 8:  # is this imag ???
            #ntotal = 32
            #s = self.self.struct_i
            #nelements = ndata // ntotal
            #for i in range(nelements):
                #edata = data[n:n + 4]
                #eid_device, = s.unpack(edata)
                #assert eid > 0, eid
                #n += ntotal
        elif result_type == 2 and op2.num_wide == 3: # random
            ntotal = 3 * 4 * self.factor
            nelements = ndata // ntotal

            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_random)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    op2.log.debug('vectorize CROD random SORT%s' % op2.sort_method)
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 3)
                obj._times[obj.itime] = dt
                self.obj_set_element(obj, itotal, itotal2, data, nelements)

                #[axial, torsion, SMa, SMt]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                n = oes_crod_random_3(op2, data, ndata, obj, nelements, ntotal)

        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
        assert op2.num_wide * 4 * self.factor == ntotal, f'numwide*4={op2.num_wide*4} ntotal={ntotal} element_name={op2.element_name!r}\n{op2.code_information()}'
        return n, nelements, ntotal

    def _oes_cbeam(self, data, ndata, dt, is_magnitude_phase,
                   result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 2 : CBEAM

        """
        op2 = self.op2
        n = 0
        ## TODO: fix method to follow correct pattern...regarding???

        if op2.is_stress:
            result_name = prefix + 'cbeam_stress' + postfix
        else:
            result_name = prefix + 'cbeam_strain' + postfix
        table_name_bytes = op2.table_name
        assert isinstance(table_name_bytes, bytes), table_name_bytes
        assert table_name_bytes in TABLES_BYTES, table_name_bytes

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)

        slot = op2.get_result(result_name)
        if result_type == 0 and op2.num_wide == 111:  # real
            # TODO: vectorize
            ntotal = 444 * self.factor # 44 + 10*40  (11 nodes)

            if op2.is_stress:
                obj_vector_real = RealBeamStressArray
            else:
                obj_vector_real = RealBeamStrainArray

            nelements = ndata // ntotal
            nlayers = nelements * 11
            auto_return, is_vectorized = op2._create_oes_object4(
                nlayers, result_name, slot, obj_vector_real)
            if auto_return:
                op2._data_factor = 11
                #assert ntotal == op2.num_wide * 4
                return nelements * ntotal, None, None
            obj = op2.obj

            ntotal = op2.num_wide * 4 * self.factor
            nelements = ndata // ntotal
            if op2.use_vector and is_vectorized and 0:
                raise NotImplementedError('CBEAM-2-real not vectorized')
            else:
                if op2.use_vector and is_vectorized and op2.sort_method == 1:
                    obj._times[obj.itime] = dt

                    n = nelements * ntotal
                    itotal = obj.itotal
                    itotal2 = itotal + nelements * 11

                    ints = frombuffer(data, dtype=op2.idtype8).reshape(nelements, 111)
                    floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 111)
                    #print(ints[:2, :].tolist())
                    #CBEAM    6       2       6       8       0.      1.      0.
                    #CBEAM    7       2       8       9       0.      1.      0.
                    #CBEAM    8       2       9       10      0.      1.      0.
                    #CBEAM    9       2       10      11      0.      1.      0.
                    #CBEAM    10      2       11      12      0.      1.      0.
                    #[[61,
                    #      6, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    #      8, 1065353216, 0, 0, 0, 0, 0, 0, 1, 1],
                    # [71,
                    #      8, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    #      9, 1065353216, 0, 0, 0, 0, 0, 0, 1, 1]]
                    eids = ints[:, 0] // 10
                    ints2 = ints[:, 1:].reshape(nelements * 11, 10)
                    #print('floats[:, 1:].shape', floats[:, 1:].shape)  # (5,110)
                    floats2 = floats[:, 1:].reshape(nelements * 11, 10)

                    xxb = floats2[:, 1]
                    # ints2 = ints[:, :2]
                    #print(ints2[0, :])
                    nids = ints2[:, 0]
                    #print("eids =", eids)
                    #print("nids =", nids.tolist())
                    #print("xxb =", xxb)

                    eids2 = array([eids] * 11, dtype=op2.idtype8).T.ravel()
                    assert len(eids2) == len(nids)
                    obj.element_node[itotal:itotal2, 0] = eids2
                    obj.element_node[itotal:itotal2, 1] = nids
                    obj.xxb[itotal:itotal2] = xxb
                    obj.data[obj.itime, itotal:itotal2, :] = floats2[:, 2:]
                    #self.data[self.itime, self.itotal, :] = [sxc, sxd, sxe, sxf,
                                                             #smax, smin, mst, msc]
                else:
                    if op2.use_vector:
                        op2.log.debug('vectorize CBEAM real SORT%s' % op2.sort_method)
                    n = oes_cbeam_real_111(op2, data,
                                           obj,
                                           nelements, dt)

        elif result_type == 1 and op2.num_wide == 111:  # imag and random?
            # definitely complex results for MSC Nastran 2016.1

            ntotal = 444 * self.factor # 44 + 10*40  (11 nodes)
            nelements = ndata // ntotal

            if op2.is_stress:
                obj_vector_complex = ComplexBeamStressArray
            else:
                obj_vector_complex = ComplexBeamStrainArray

            nlayers = nelements * 11
            auto_return, is_vectorized = op2._create_oes_object4(
                nlayers, result_name, slot, obj_vector_complex)
            if auto_return:
                op2._data_factor = 11
                return nelements * ntotal, None, None

            obj = op2.obj

            nnodes = 10  # 11-1
            #ntotal = op2.num_wide * 4
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.itotal
                itotal2 = itotal + nelements * 11

                # chop off eid
                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 111)[:, 1:]
                floats2 = floats.reshape(nelements * 11, 10).copy()

                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype8).reshape(nelements, 111)
                    eids = ints[:, 0] // 10
                    eids2 = array([eids] * 11, dtype=op2.idtype8).T.ravel()

                    ints2 = ints[:, 1:].reshape(nelements * 11, 10)

                    nids = ints2[:, 0]
                    assert eids.min() > 0, eids.min()
                    assert nids.min() >= 0, nids.min()
                    obj.element_node[itotal:itotal2, 0] = eids2
                    obj.element_node[itotal:itotal2, 1] = nids

                #  0    1   2  3  4  5  6   7   8   9
                # grid, sd, c, d, e, f, c2, d2, e2, f2
                real_imag = apply_mag_phase(floats2, is_magnitude_phase, [2, 3, 4, 5], [6, 7, 8, 9])
                obj.data[obj.itime, itotal:itotal2, :] = real_imag
                obj.sd[itotal:itotal2] = floats2[:, 1]

                obj.itotal = itotal2
                #obj.ielement = ielement2
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    op2.log.debug('vectorize CBEAM imag SORT%s' % op2.sort_method)
                n = oes_cbeam_complex_111(op2, data, obj,
                                          nelements, nnodes,
                                          is_magnitude_phase)

        elif result_type == 2 and op2.num_wide == 67: # random
            # TODO: vectorize
            ntotal = 268 # 1 + 11*6  (11 nodes)

            if op2.is_stress:
                obj_vector_random = RandomBeamStressArray
            else:
                obj_vector_random = RandomBeamStrainArray

            nelements = ndata // ntotal
            nlayers = nelements * 11
            auto_return, is_vectorized = op2._create_oes_object4(
                nlayers, result_name, slot, obj_vector_random)
            if auto_return:
                op2._data_factor = 11
                assert ntotal == op2.num_wide * 4
                return nelements * ntotal, None, None
            obj = op2.obj

            nnodes = 10  # 11-1
            ntotal = op2.num_wide * 4
            nelements = ndata // ntotal
            if op2.use_vector and is_vectorized and 0:
                raise NotImplementedError('CBEAM-2-random not vectorized')
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    op2.log.debug('vectorize CBEAM random SORT%s' % op2.sort_method)
                n = oes_cbeam_random_67(op2, data, obj,
                                        nelements, nnodes, dt)

        elif result_type == 1 and op2.num_wide in [67] and table_name_bytes in [b'OESXNO1']:  # CBEAM
            # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1081x.op2
            msg = 'skipping freq CBEAM; numwide=67'
            n = op2._not_implemented_or_skip(data, ndata, msg)
            nelements = None
            ntotal = None
        elif result_type == 2 and op2.num_wide in [67] and table_name_bytes in [b'OESXNO1']:  # CBEAM
            #C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1081x.op2
            msg = 'skipping random CBEAM; numwide=67'
            n = op2._not_implemented_or_skip(data, ndata, msg)
            nelements = None
            ntotal = None
        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
        return n, nelements, ntotal

    def _oes_cshear(self, data, ndata, dt, is_magnitude_phase,
                    result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 4 : CSHEAR
        """
        op2: OP2 = self.op2
        n = 0
        # 4-CSHEAR
        if op2.is_stress:
            obj_vector_real = RealShearStressArray
            obj_vector_complex = ComplexShearStressArray
            obj_vector_random = RandomShearStressArray
            result_name = prefix + 'cshear_stress' + postfix
        else:
            obj_vector_real = RealShearStrainArray
            obj_vector_complex = ComplexShearStrainArray
            obj_vector_random = RandomShearStrainArray
            result_name = prefix + 'cshear_strain' + postfix

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)

        slot = op2.get_result(result_name)
        if result_type == 0 and op2.num_wide == 4:  # real
            ntotal = 16 * self.factor # 4*4
            nelements = ndata // ntotal
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return nelements * ntotal, None, None

            obj: Union[RealShearStressArray, RealShearStrainArray] = op2.obj
            assert obj is not None
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 4)
                itime = obj.itime
                obj._times[itime] = dt
                self.obj_set_element(obj, itotal, itotal2, data, nelements)

                #[max_strain, avg_strain, margin]
                obj.data[itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    op2.log.debug('vectorize CSHEAR real SORT%s' % op2.sort_method)
                n = oes_cshear_real_4(op2, data, obj,
                                      ntotal, nelements, dt)

        elif result_type == 1 and op2.num_wide == 5:  # imag
            ntotal = 20 * self.factor # 4*5
            nelements = ndata // ntotal
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_complex)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 5).copy()
                obj._times[obj.itime] = dt
                self.obj_set_element(obj, itotal, itotal2, data, nelements)

                #(eid_device, etmaxr, etmaxi, etavgr, etavgi)
                real_imag = apply_mag_phase(floats, is_magnitude_phase, [1, 3], [2, 4])
                obj.data[obj.itime, itotal:itotal2, :] = real_imag
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    op2.log.debug('vectorize CSHEAR imag SORT%s' % op2.sort_method)
                n = oes_cshear_complex_5(op2, data, obj,
                                         nelements, ntotal, is_magnitude_phase)
        elif result_type == 2 and op2.num_wide == 3: # random
            ntotal = 12  # 3*4
            nelements = ndata // ntotal
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_random)
            if auto_return:
                assert ntotal == op2.num_wide * 4
                return nelements * ntotal, None, None

            obj = op2.obj
            assert obj is not None
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * 4 * op2.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 3)
                itime = obj.itime
                obj._times[itime] = dt
                self.obj_set_element(obj, itotal, itotal2, data, nelements)

                #[max_strain, avg_strain, margin]
                obj.data[itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    op2.log.debug('vectorize CSHEAR random SORT%s' % op2.sort_method)
                n = oes_cshear_random_3(op2, data, obj,
                                        nelements, ntotal)

        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
        return n, nelements, ntotal

    def _oes_cbar_34(self, data: bytes, ndata: int, dt: Any,
                     is_magnitude_phase: bool,
                     result_type: str, prefix: str, postfix: str) -> tuple[int, int, int]:
        """
        reads stress/strain for element type:
         - 34 : CBAR

        """
        op2 = self.op2
        #if isinstance(op2.nonlinear_factor, float):
            #op2.sort_bits[0] = 1 # sort2
            #op2.sort_method = 2

        n = 0
        if op2.is_stress:
            result_name = prefix + 'cbar_stress' + postfix
        else:
            result_name = prefix + 'cbar_strain' + postfix

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)
        if result_type == 0 and op2.num_wide == 16:  # real
            if op2.is_stress:
                obj_vector_real = RealBarStressArray
            else:
                obj_vector_real = RealBarStrainArray

            ntotal = 64 * self.factor  # 16*4
            nelements = ndata // ntotal
            #print('CBAR nelements =', nelements)

            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return ndata, None, None

            if op2.is_debug_file:
                op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                op2.binary_debug.write('  #elementi = [eid_device, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,\n')
                op2.binary_debug.write('                           s1b, s2b, s3b, s4b, smaxb, sminb,        MSc]\n')
                op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                # self.itime = 0
                # self.ielement = 0
                # self.itotal = 0
                #self.ntimes = 0
                #self.nelements = 0
                n = nelements * op2.num_wide * 4

                ielement = obj.ielement
                ielement2 = ielement + nelements
                obj._times[obj.itime] = dt
                self.obj_set_element(obj, ielement, ielement2, data, nelements)

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 16)

                #[s1a, s2a, s3a, s4a, axial, smaxa, smina, margin_tension,
                # s1b, s2b, s3b, s4b,        smaxb, sminb, margin_compression]
                obj.data[obj.itime, ielement:ielement2, :] = floats[:, 1:].copy()
                obj.itotal = ielement2
                obj.ielement = ielement2
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    op2.log.debug('vectorize CBAR real SORT%s' % op2.sort_method)
                n = oes_cbar_real_16(op2, data, obj, nelements, ntotal, dt)
        elif result_type == 1 and op2.num_wide == 19:  # imag
            if op2.is_stress:
                obj_vector_complex = ComplexBarStressArray
            else:
                obj_vector_complex = ComplexBarStrainArray

            ntotal = 76 * self.factor
            nelements = ndata // ntotal
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_complex)
            if auto_return:
                return ndata, None, None

            if op2.is_debug_file:
                op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                op2.binary_debug.write('  #elementi = [eid_device, s1a, s2a, s3a, s4a, axial,\n')
                op2.binary_debug.write('                           s1b, s2b, s3b, s4b]\n')
                op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.itotal
                itotal2 = itotal + nelements
                ielement2 = itotal2

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 19).copy()
                obj._times[obj.itime] = dt
                self.obj_set_element(obj, itotal, itotal2, data, nelements)

                isave1 = [1, 2, 3, 4, 5, 11, 12, 13, 14]
                isave2 = [6, 7, 8, 9, 10, 15, 16, 17, 18]
                real_imag = apply_mag_phase(floats, is_magnitude_phase, isave1, isave2)
                obj.data[obj.itime, itotal:itotal2, :] = real_imag

                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    op2.log.debug('vectorize CBAR imag SORT%s' % op2.sort_method)
                n = oes_cbar_complex_19(op2, data, obj, nelements, ntotal, is_magnitude_phase)
        elif result_type == 2 and op2.num_wide == 19: # random strain?
            raise RuntimeError(op2.code_information())

        elif result_type in [1, 2] and op2.num_wide == 10:  # random
            # random stress/strain per example
            #
            # DMAP says random stress has num_wide=10 and
            # random strain has numwide=19, but it's wrong...maybe???
            #
            # format_code = 1 - NO/RMS (SORT1 regardless of whether this is a SORT2 table or not)
            # format_code = 2 - ATO/PSD/CRM (actually SORT2)
            #
            element_id = op2.nonlinear_factor
            if op2.is_stress:
                obj_vector_random = RandomBarStressArray
            else:
                obj_vector_random = RandomBarStrainArray
            op2.data_code['nonlinear_factor'] = element_id

            ntotal = 10 * self.size
            nelements = ndata // ntotal
            #print(f'CBAR* nelements={nelements}')
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_random)
            if auto_return:
                return ndata, None, None

            if op2.is_debug_file:
                op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                op2.binary_debug.write('  #elementi = [eid_device, s1a, s2a, s3a, s4a, axial,\n')
                op2.binary_debug.write('                           s1b, s2b, s3b, s4b]\n')
                op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            obj = op2.obj
            if op2.use_vector and is_vectorized and 0:  # pragma: no cover
                # self.itime = 0
                # self.ielement = 0
                # self.itotal = 0
                #self.ntimes = 0
                #self.nelements = 0
                n = nelements * ntotal

                ielement = obj.ielement
                ielement2 = ielement + nelements
                obj._times[obj.itime] = dt
                self.obj_set_element(obj, itotal, itotal2, data, nelements)

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 10)

                #[s1a, s2a, s3a, s4a, axial,
                # s1b, s2b, s3b, s4b]
                obj.data[obj.itime, ielement:ielement2, :] = floats[:, 1:].copy()
                obj.itotal = ielement2
                obj.ielement = ielement2
            else:
                if is_vectorized and op2.use_vector and obj.itime == 0:  # pragma: no cover
                    op2.log.debug('vectorize CBAR random SORT%s' % op2.sort_method)
                n = oes_cbar_random_10(op2, data, obj,
                                       nelements, ntotal)
        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
        return n, nelements, ntotal

    def _oes_csolid(self, data, ndata: int, dt, is_magnitude_phase: bool,
                    result_type: str, prefix: str, postfix: str):
        """
        reads stress/strain for element type:
         - 39 : CTETRA
         - 67 : CHEXA
         - 68 : CPENTA
         - 255: CPYRAM

        """
        op2 = self.op2
        n = 0
        etype_map = {
            #element_type : (element_base, nnodes_expected, element_name)
            39 : ('ctetra', 5, 'CTETRA4'),
            67 : ('chexa', 9, 'CHEXA8'),
            68 : ('cpenta', 7, 'CPENTA6'),
            255 : ('cpyram', 6, 'CPYRAM5'),
        }
        element_base, nnodes_expected, element_name = etype_map[op2.element_type]

        if op2.is_stress:
            stress_strain = 'stress'
            obj_vector_real = RealSolidStressArray
            obj_vector_complex = ComplexSolidStressArray
            obj_vector_complex_vm = ComplexSolidStressVMArray
            obj_vector_random = RandomSolidStressArray
            obj_vector_random_vm = RandomSolidVMStressArray
        else:
            stress_strain = 'strain'
            obj_vector_real = RealSolidStrainArray
            obj_vector_complex = ComplexSolidStrainArray
            obj_vector_complex_vm = ComplexSolidStrainVMArray
            obj_vector_random = RandomSolidStrainArray
            obj_vector_random_vm = RandomSolidVMStrainArray

        if prefix == '' and postfix == '':
            prefix = stress_strain + '.'

        # stress.chexa_stress
        result_name = prefix + f'{element_base}_{stress_strain}' + postfix

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        numwide_real = 4 + 21 * nnodes_expected
        numwide_imag = 4 + (17 - 4) * nnodes_expected
        numwide_imag2 = 4 + 14 * nnodes_expected
        numwide_random = 4 + (11 - 4) * nnodes_expected
        numwide_random2_vm = 4 + 14 * nnodes_expected
        numwide_random3 = 4 + 8 * nnodes_expected
        preline1 = '%s-%s' % (op2.element_name, op2.element_type)
        preline2 = ' ' * len(preline1)

        #print('numwide real=%s imag=%s random=%s' % (numwide_real, numwide_imag, numwide_random2))
        op2._data_factor = nnodes_expected
        log = op2.log
        # print(op2.format_code, op2.table_name, op2.num_wide, op2.element_name)
        # print(op2.code_information())
        if op2.format_code == 1 and op2.num_wide == numwide_real:  # real
            op2.log.debug(f'numwide_real={numwide_real}')
            ntotal = (16 + 84 * nnodes_expected) * self.factor
            nelements = ndata // ntotal
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            #sort_method = op2.sort_method2()
            #sort_method = self.sort_method2
            sort_method = op2.sort_method
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                itotali = obj.itotal + nelements
                itotal2 = obj.itotal + nelements * nnodes_expected
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    # (eid_device, cid, abcd, nnodes)
                    ints = frombuffer(data, dtype=op2.idtype8).copy()
                    try:
                        ints1 = ints.reshape(nelements, numwide_real)
                    except ValueError:
                        msg = 'ints.shape=%s; size=%s ' % (str(ints.shape), ints.size)
                        msg += 'nelements=%s numwide_real=%s nelements*numwide=%s' % (
                            nelements, numwide_real, nelements * numwide_real)
                        raise ValueError(msg)
                    eids = ints1[:, 0] // 10
                    cids = ints1[:, 1]
                    #nids = ints1[:, 4]
                    assert eids.min() > 0, eids.min()
                    obj.element_node[itotal:itotal2, 0] = repeat(eids, nnodes_expected)
                    ints2 = ints1[:, 4:].reshape(nelements * nnodes_expected, 21)
                    grid_device = ints2[:, 0]#.reshape(nelements, nnodes_expected)

                    #print('%s-grid_device=%s' % (op2.element_name, grid_device))
                    unused_grid_device2 = repeat(grid_device, nnodes_expected)
                    try:
                        obj.element_node[itotal:itotal2, 1] = grid_device
                    except ValueError:
                        msg = '%s; nnodes=%s\n' % (op2.element_name, nnodes_expected)
                        msg += 'itotal=%s itotal2=%s\n' % (itotal, itotal2)
                        msg += 'grid_device.shape=%s; size=%s\n' % (str(grid_device.shape), grid_device.size)
                        #msg += 'nids=%s' % nids
                        raise ValueError(msg)
                    #log.debug(f'cids = {np.unique(cids)}')
                    obj.element_cid[itotal:itotali, 0] = eids
                    obj.element_cid[itotal:itotali, 1] = cids

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, numwide_real)[:, 4:]
                # 1     9    15   2    10   16  3   11  17   8
                #[oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, ovm]
                #isave = [1, 9, 15, 2, 10, 16, 3, 11, 17, 8]
                #(grid_device,
                    #sxx, sxy, s1, a1, a2, a3, pressure, svm,
                    #syy, syz, s2, b1, b2, b3,
                    #szz, sxz, s3, c1, c2, c3)
                floats1 = floats.reshape(nelements * nnodes_expected, 21)#[:, 1:] # drop grid_device

                # o1/o2/o3 is not max/mid/min.  They are not consistently ordered, so we force it.
                max_mid_min = np.vstack([
                    floats1[:, 3],
                    floats1[:, 11],
                    floats1[:, 17],
                ]).T
                max_mid_min.sort(axis=1)
                assert max_mid_min.shape == (nelements * nnodes_expected, 3), max_mid_min.shape
                obj.data[obj.itime, itotal:itotal2, 6:9] = max_mid_min[:, [2, 1, 0]]

                #obj.data[obj.itime, itotal:itotal2, :] = floats1[:, isave]
                obj.data[obj.itime, itotal:itotal2, :6] = floats1[:, [1, 9, 15, 2, 10, 16]]
                obj.data[obj.itime, itotal:itotal2, 9] = floats1[:, 8]
                obj.itotal = itotal2
                obj.ielement = itotali
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    log.debug(f'vectorize CSolid real SORT{sort_method} from {op2.table_name}')
                n = oes_csolid_real(op2, data, obj,
                                    nelements, dt,
                                    element_name, nnodes_expected,
                                    preline1, preline2)

        elif op2.format_code in [2, 3] and op2.num_wide == numwide_imag:  # complex
            op2.log.debug(f'numwide_imag={numwide_imag}')
            ntotal = numwide_imag * 4 * self.factor
            nelements = ndata // ntotal
            self.ntotal += nelements * nnodes_expected
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_complex)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj

            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal
                ielement = obj.ielement
                ielement2 = ielement + nelements
                itotal = obj.itotal
                itotal2 = itotal + nelements * nnodes_expected

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, numwide_imag)
                floats1 = floats[:, 4:].reshape(nelements * nnodes_expected, 13).copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype8).reshape(nelements, numwide_imag)
                    ints1 = ints[:, 4:].reshape(nelements * nnodes_expected, 13)
                    eids = ints[:, 0] // 10
                    cids = ints[:, 1]
                    nids = ints1[:, 0]
                    # TODO: ctype, nodef not considered
                    assert eids.min() > 0, eids.min()
                    assert nids.min() >= 0, nids.min()
                    eids2 = np.vstack([eids] * nnodes_expected).T.ravel()
                    #nids2 = np.vstack([nids] * nnodes_expected).T.ravel()
                    #print(nids2)
                    obj.element_node[itotal:itotal2, 0] = eids2
                    obj.element_node[itotal:itotal2, 1] = nids

                    obj.element_cid[ielement:ielement2, 0] = eids
                    obj.element_cid[ielement:ielement2, 1] = cids

                # 0 is nid
                isave1 = [1, 2, 3, 4, 5, 6]
                isave2 = [7, 8, 9, 10, 11, 12]
                real_imag = apply_mag_phase(floats1, is_magnitude_phase, isave1, isave2)
                obj.data[obj.itime, itotal:itotal2, :] = real_imag

                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    log.debug(f'vectorize CSolid imag SORT{op2.sort_method}')
                n = oes_csolid_complex(op2, data, obj,
                                       nelements, # nnodes,
                                       element_name, nnodes_expected,
                                       is_magnitude_phase)

        elif op2.format_code == 1 and op2.num_wide == numwide_random: # random
            op2.log.debug(f'numwide_random={numwide_random}')
            if not op2.is_sort1:
                log.debug(f'support CSolid random SORT{op2.sort_method}')
                return ndata, None, None

            ntotal = numwide_random * 4
            nelements = ndata // ntotal
            assert ndata % ntotal == 0, ndata
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_random)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and 0:  # pragma: no cover
                n = nelements * 4 * op2.num_wide
                itotal = obj.ielement
                itotali = obj.itotal + nelements
                itotal2 = obj.itotal + nelements * nnodes_expected
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    # (eid_device, cid, abcd, nnodes)
                    ints = frombuffer(data, dtype=op2.idtype).copy()
                    try:
                        ints1 = ints.reshape(nelements, numwide_real)
                    except ValueError:
                        msg = 'ints.shape=%s; size=%s ' % (str(ints.shape), ints.size)
                        msg += 'nelements=%s numwide_real=%s nelements*numwide=%s' % (
                            nelements, numwide_real, nelements * numwide_real)
                        raise ValueError(msg)
                    eids = ints1[:, 0] // 10
                    cids = ints1[:, 1]
                    #nids = ints1[:, 4]
                    assert eids.min() > 0, eids.min()
                    obj.element_node[itotal:itotal2, 0] = repeat(eids, nnodes_expected)
                    ints2 = ints1[:, 4:].reshape(nelements * nnodes_expected, 21)
                    grid_device = ints2[:, 0]#.reshape(nelements, nnodes_expected)

                    #print('%s-grid_device=%s' % (op2.element_name, grid_device))
                    unused_grid_device2 = repeat(grid_device, nnodes_expected)
                    try:
                        obj.element_node[itotal:itotal2, 1] = grid_device
                    except ValueError:
                        msg = '%s; nnodes=%s\n' % (op2.element_name, nnodes_expected)
                        msg += 'itotal=%s itotal2=%s\n' % (itotal, itotal2)
                        msg += 'grid_device.shape=%s; size=%s\n' % (str(grid_device.shape), grid_device.size)
                        #msg += 'nids=%s' % nids
                        raise ValueError(msg)
                    obj.element_cid[itotal:itotali, 0] = eids
                    obj.element_cid[itotal:itotali, 1] = cids

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, numwide_real)[:, 4:]
                # 1     9    15   2    10   16  3   11  17   8
                #[oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, ovm]
                #isave = [1, 9, 15, 2, 10, 16, 3, 11, 17, 8]
                #(grid_device,
                    #sxx, sxy, s1, a1, a2, a3, pressure, svm,
                    #syy, syz, s2, b1, b2, b3,
                    #szz, sxz, s3, c1, c2, c3)
                floats1 = floats.reshape(nelements * nnodes_expected, 21)#[:, 1:] # drop grid_device

                # o1/o2/o3 is not max/mid/min.  They are not consistently ordered, so we force it.
                max_mid_min = np.vstack([
                    floats1[:, 3],
                    floats1[:, 11],
                    floats1[:, 17],
                ]).T
                max_mid_min.sort(axis=1)
                assert max_mid_min.shape == (nelements * nnodes_expected, 3), max_mid_min.shape
                obj.data[obj.itime, itotal:itotal2, 6:9] = max_mid_min[:, [2, 1, 0]]

                #obj.data[obj.itime, itotal:itotal2, :] = floats1[:, isave]
                obj.data[obj.itime, itotal:itotal2, :6] = floats1[:, [1, 9, 15, 2, 10, 16]]
                obj.data[obj.itime, itotal:itotal2, 9] = floats1[:, 8]
                obj.itotal = itotal2
                obj.ielement = itotali
            else:
                if is_vectorized and op2.use_vector and obj.itime == 0:  # pragma: no cover
                    log.debug(f'vectorize CSolid random SORT{op2.sort_method}')
                n = oes_csolid_random(op2, data, obj, nelements,
                                      element_name, nnodes_expected,
                                      preline1, preline2)

        elif op2.format_code in [2, 3] and op2.num_wide == numwide_random2_vm and \
                op2.table_name in {b'OESVM1', b'OSTRVM1', b'OESVM2', b'OSTRVM2'}:
            #op2.log.debug(f'numwide_random2={numwide_random2}')
            #raise RuntimeError(op2.code_information())
            ## a = 18
            ## b = 14
            ## a + b * nnodes = numwide_random3
            ## a + b * 4 = 74  # CTETRA
            ## a + b * 6 = 102 # CPENTA
            ## a + b * 8 = 130 # CHEXA-67
            #msg = 'OES-CHEXA-random-numwide=%s numwide_real=%s numwide_imag=%s numwide_random=%s' % (
                #op2.num_wide, numwide_real, numwide_imag, numwide_random)
            #return op2._not_implemented_or_skip(data, ndata, msg)

            #print('numwide real=%s imag=%s random=%s' % (numwide_real, numwide_imag, numwide_random))
            unused_num_wide_random = 4 + nnodes_expected * (17 - 4)

            #print('random2=%s' % num_wide_random)
            #print(op2.code_information())

            #print('numwide=%s numwide_random=%s attempt2=%s subcase=%s' % (
                #op2.num_wide, numwide_random, num_wide_random, op2.isubcase))
            ##msg = op2.code_information()
            #ntotal = 130
            ntotal = op2.num_wide * op2.size
            nelements = ndata // ntotal

            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_complex_vm)
            if auto_return:
                return nelements * ntotal, None, None
            ## cid, coord_type, nactive_pnts,
            ##      nid, oxx, oyy, ozz, txy, tyz, txz
            ntotal1 = 4 * self.size # 3*4
            ntotal2 = 14 * self.size
            ntotali = ntotal1 + nnodes_expected * ntotal2
            assert ntotali == ntotal, (ntotali, ntotal)

            # 2 CID I Stress Coordinate System
            # 3 CTYPE CHAR4 Coordinate System Type (BCD)
            # 4 NODEF I Number of Active Points
            #
            # 5 GRID I External grid identification number (0=center)
            # 6 SXR RS Normal in x
            # 7 SYR RS Normal in y
            # 8 SZR RS Normal in z
            # 9 TXYR RS Shear in xy
            # 10 TYZR RS Shear in yz
            # 11 TZXR RS Shear in zx
            # 12 SXI RS Normal in x
            # 13 SYI RS Normal in y
            # 14 SZI RS Normal in z
            # 15 TXYI RS Shear in xy
            # 16 TYZI RS Shear in yz
            # 17 TZXI RS Shear in zx
            # 18 VM RS von Mises
            # Words 5 through 18 repeat 005 times
            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                itotali = obj.itotal + nelements
                itotal2 = obj.itotal + nelements * nnodes_expected
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    # (eid_device, cid, abcd, nnodes)
                    ints = frombuffer(data, dtype=op2.idtype8).copy()
                    try:
                        ints1 = ints.reshape(nelements, numwide_random2_vm)
                    except ValueError:
                        msg = 'ints.shape=%s; size=%s ' % (str(ints.shape), ints.size)
                        msg += 'nelements=%s numwide_real=%s nelements*numwide=%s' % (
                            nelements, numwide_real, nelements * numwide_real)
                        raise ValueError(msg)

                    eids = ints1[:, 0] // 10
                    cids = ints1[:, 1]
                    #print(f'eids = {eids}')
                    assert eids.min() > 0, eids.min()
                    obj.element_node[itotal:itotal2, 0] = repeat(eids, nnodes_expected)
                    ints2 = ints1[:, 4:].reshape(nelements * nnodes_expected, 14)
                    grid_device = ints2[:, 0]#.reshape(nelements, nnodes_expected)
                    #print(ints1)
                    #print(ints2)
                    #print(nelements, nnodes_expected, nelements*nnodes_expected)
                    #print(f'nids = {grid_device}')
                    #print('ints2.shape =', ints2.shape)

                    #print('%s-grid_device=%s' % (op2.element_name, grid_device))
                    unused_grid_device2 = repeat(grid_device, nnodes_expected)
                    try:
                        obj.element_node[itotal:itotal2, 1] = grid_device
                    except ValueError:
                        msg = '%s; nnodes=%s\n' % (op2.element_name, nnodes_expected)
                        msg += 'itotal=%s itotal2=%s\n' % (itotal, itotal2)
                        msg += 'grid_device.shape=%s; size=%s\n' % (str(grid_device.shape), grid_device.size)
                        #msg += 'nids=%s' % nids
                        raise ValueError(msg)
                    obj.element_cid[itotal:itotali, 0] = eids
                    obj.element_cid[itotal:itotali, 1] = cids

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, numwide_random2_vm)[:, 4:]
                # 1    2    3    4    5    6    7 - verify...
                #[oxx, oyy, ozz, txy, tyz, txz, ovm]
                # sxr, syr, szr, txyr, tyzr, txzr, \
                #     sxi, syi, szi, txyi, tyzi, txzi, ovm = out
                isave1 = [1, 2, 3,  4,  5,  6]
                isave2 = [7, 8, 9, 10, 11, 12]
                floats1 = floats.reshape(nelements * nnodes_expected, 14)
                real_imag = apply_mag_phase(floats1, is_magnitude_phase, isave1, isave2)
                #floats2 = floats1[:, 1:] # drop grid_device

                # o1/o2/o3 is not max/mid/min.  They are not consistently ordered, so we force it.
                obj.data[obj.itime, itotal:itotal2, -1] = floats1[:, -1]
                obj.data[obj.itime, itotal:itotal2, :-1] = real_imag
                obj.itotal = itotal2
                obj.ielement = itotali
            else:
                #op2.log.warning(f'vectorize complex solid vm SORT1; numwide={op2.num_wide}')
                struct1 = Struct(op2._endian + b'2i 4s i')
                struct2 = Struct(op2._endian + b'i 13f')
                element_num = op2.element_type
                element_type = op2.element_name
                for i in range(nelements):
                    edata = data[n:n+ntotal1]
                    out = struct1.unpack(edata)
                    (eid_device, cid, ctype, nodef) = out
                    #print(eid_device, cid, ctype, nodef)
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, op2.nonlinear_factor, op2.sort_method)
                    #obj.add_sort1(eid, dt, cid)
                    obj.add_eid_sort1(element_num, element_type, dt, eid, cid, ctype, nodef)

                    if op2.is_debug_file:
                        op2.binary_debug.write('%s - eid=%i; %s\n' % (preline1, eid, str(out)))
                    n += ntotal1
                    for inode in range(nnodes_expected):  # nodes pts, +1 for centroid (???)
                        grid, *out = struct2.unpack(data[n:n + ntotal2])
                        #print(' ', grid, out)
                        #if op2.is_debug_file:
                            #op2.binary_debug.write('%s - %s\n' % (preline2, str(out)))
                        #(grid_device, sxx, syy, sz, txy, tyz, txz) = out
                        #     ELEMENT-ID    GRID-ID    NORMAL-X       NORMAL-Y       NORMAL-Z       SHEAR-XY       SHEAR-YZ       SHEAR-ZX      VON MISES
                        # 0                     1048   2.823925E+03  -8.389856E+03  -1.373890E+03  -3.886024E+03  -9.440013E+02  -8.069132E+02
                        #                             -1.142355E+02   3.635427E+02   6.020642E+01   1.693545E+02   4.121148E+01   3.410083E+01   1.210361E+04

                        sxr, syr, szr, txyr, tyzr, txzr, \
                        sxi, syi, szi, txyi, tyzi, txzi, ovm = out
                        if is_magnitude_phase:
                            sx = polar_to_real_imag(sxr, sxi)
                            sy = polar_to_real_imag(syr, syi)
                            sz = polar_to_real_imag(szr, szi)
                            txy = polar_to_real_imag(txyr, txyi)
                            tyz = polar_to_real_imag(tyzr, tyzi)
                            txz = polar_to_real_imag(txzr, txzi)
                        else:
                            sx = complex(sxr, sxi)
                            sy = complex(syr, syi)
                            sz = complex(szr, szi)
                            txy = complex(txyr, txyi)
                            tyz = complex(tyzr, tyzi)
                            txz = complex(txzr, txzi)
                        #del grid, sx, sy, sz, txy, tyz, txz
                        #if eid == 1048:
                            #print(f'grid={grid}', sx, sy, sz, txy, tyz, txz)
                            #out2 = [round(val, 2) for val in out]
                            #print(grid, out2)
                        obj.add_node_sort1(dt, eid, grid, inode, sx, sy, sz, txy, tyz, txz, ovm)
                        n += ntotal2
            #etype_num = f'{op2.element_name}-{op2.element_type:d};'
            #msg = f'OES-{etype_num:<11s} complex2; table_name={op2.table_name}; numwide={op2.num_wide}'
            #return op2._not_implemented_or_skip(data, ndata, msg), None, None
        elif op2.format_code in [1, 2] and op2.num_wide == 67 and op2.element_name == 'CHEXA':  # CHEXA
            op2.log.debug(f'numwide complex HEXA?')
            # CTETRA:
            # CPYRAM:
            # CPENTA:
            # CHEXA:  67
            #44 = 5 * 8 + 4  (TETRA)
            #52 = 6 * 8 + 4  (PYRAM)
            #60 = 7 * 8 + 4  (PENTA)
            #76 = 9 * 8 + 4  (HEXA)
            msg = 'skipping random CHEXA; numwide=67'
            #print(op2.code_information())
            #asdf
            n = op2._not_implemented_or_skip(data, ndata, msg)
            nelements = None
            ntotal = None
        elif op2.format_code in [1, 2, 3] and op2.num_wide in [60] and op2.table_name in [b'OESXRMS1', b'OESXNO1']:  # CPENTA
            op2.log.debug(f'numwide complex PENTA?')
            # bad
            #if op2.read_mode == 2:
                #ints    = (68011, 0, 805.28, 6,
                           #0,   0, 0, 0, 0, 0, 0, 0,
                           #120000, 0, 0, 0, 0, 0, 0, 0,
                           #120001, 0, 0, 0, 0, 0, 0, 0,
                           #120010, 0, 0, 0, 0, 0, 0, 0,
                           #120100, 0, 0, 0, 0, 0, 0, 0,
                           #120101, 0, 0, 0, 0, 0, 0, 0,
                           #120110, 0, 0, 0, 0, 0, 0, 0,

                           #68111, 0, 1145655879, 15,
                           #0,      1080284864, 1080296481, 1080284990, 1080279426, 1080285072, 1080285570, 1080285750,
                           #130000, 1080285656, 1080285656, 1080285656, 1080285162, 1080287537, 1080285308, 1080285794,
                           #130002, 1080285656, 1080285656, 1080285656, 1080284551, 1080287537, 1080285308, 1080285794,
                           #130020, 1080285656, 1080285656, 1080285656, 1080289401, 1080287537, 1080285308, 1080285794,
                           #130200, 1080285656, 1080285656, 1080285656, 1080285269, 1080287537, 1080285308, 1080285794,
                           #130202, 1080288409, 1080287759, 1080323139, 1080285308, 1080285512, 1080285308, 1080285874,
                           #130220, 1080285333, 1080285373, 1080285450, 1080287537, 1080287537, 1080285625, 1080285771)
                #floats  = (68011, 0.0, 805.28, 6,
                           #0.0,                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           #1.681558157189705e-40, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           #1.681572170174437e-40, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           #1.681698287036213e-40, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           #1.682959455654153e-40, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           #1.682973468638786e-40, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           #1.683099585500578e-40, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,

                           #9.544383970362762e-41, 0.0, 805.28, 2.1019476964872256e-44,
                           #0.0,                    3.559982299804687, 3.562752008438110, 3.5600123405456, 3.558685779571533, 3.560031890869140, 3.560150623321533, 3.56019353866577,
                           #1.8216880036222622e-40, 3.560171127319336, 3.560171127319336, 3.5601711273193, 3.560053348541259, 3.560619592666626, 3.560088157653808, 3.56020402908325,
                           #1.8217160295915487e-40, 3.560171127319336, 3.560171127319336, 3.5601711273193, 3.559907674789428, 3.560619592666626, 3.560088157653808, 3.56020402908325,
                           #1.8219682633151272e-40, 3.560171127319336, 3.560171127319336, 3.5601711273193, 3.561064004898071, 3.560619592666626, 3.560088157653808, 3.56020402908325,
                           #1.8244906005509118e-40, 3.560171127319336, 3.560171127319336, 3.5601711273193, 3.560078859329223, 3.560619592666626, 3.560088157653808, 3.56020402908325,
                           #1.8245186265201983e-40, 3.560827493667602, 3.560672521591186, 3.5691077709198, 3.560088157653808, 3.560136795043945, 3.560088157653808, 3.56022310256958,
                           #1.8247708602437768e-40, 3.560094118118286, 3.56010365486145,  3.5601220130920, 3.560619592666626, 3.560619592666626, 3.560163736343384, 3.56019854545593)
                #self.show_data(data, types='ifs')

            #C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1081x.op2
            msg = 'skipping random CPENTA; numwide=60'
            n = op2._not_implemented_or_skip(data, ndata, msg)
            nelements = None
            ntotal = None
        elif op2.format_code in [1, 2, 3] and op2.num_wide == numwide_random3 and op2.table_name in [b'OESXRMS1', b'OESXNO1']:
            op2.log.debug(f'numwide_random3={numwide_random3}')
            # CTETRA: 44
            # CPYRAM: 52
            # CHEXA:  76
            # '      FREQUENCY =  0.000000E+00'
            # '                      S T R E S S E S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )'
            # '                                    ( ROOT MEAN SQUARE; RMSSF SCALE FACTOR =  1.00E+00 )'
            # '0                   CORNER      --------------------------CENTER AND CORNER POINT STRESSES---------------------------'
            # '      ELEMENT-ID   GRID-ID      NORMAL-X    NORMAL-Y    NORMAL-Z      SHEAR-XY    SHEAR-YZ    SHEAR-ZX   VON MISES'
            # '0            6701        0GRID CS  8 GP'
            # '0                   CENTER     6.413E+01   2.752E-01   6.413E+01     6.496E+00   6.496E+00   3.842E+01   9.359E+01'
            # '0                    30000     5.367E+00   1.069E+01   5.367E+00     1.134E-15   1.134E-15   4.473E-15   5.321E+00'
            # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1081x.op2

            # 2 CID       I Stress Coordinate System
            # 3 CTYPE CHAR4 Coordinate System Type (BCD)
            # 4 NODEF     I Number of Active Points
            # 5 GRID      I External grid identification number (0=center)
            msg = f'skipping random {element_name}; numwide={numwide_random3}'

            ntotal = op2.num_wide * 4
            nelements = ndata // ntotal
            assert nelements > 0, nelements
            assert ndata % ntotal == 0, ndata

            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_random_vm)
            if auto_return:
                return nelements * ntotal, None, None
            if op2.read_mode == 1:
                return ndata, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized:  # pragma: no cover
                is_finished = True
                n = nelements * 4 * op2.num_wide
                itotal = obj.ielement
                itotali = obj.itotal + nelements
                itotal2 = obj.itotal + nelements * nnodes_expected
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    # (eid_device, cid, abcd, nnodes)
                    ints = frombuffer(data, dtype=op2.idtype).copy()
                    try:
                        ints1 = ints.reshape(nelements, numwide_random3)
                    except ValueError:
                        msg = 'ints.shape=%s; size=%s ' % (str(ints.shape), ints.size)
                        msg += 'nelements=%s numwide_random3=%s nelements*numwide=%s' % (
                            nelements, numwide_random3, nelements * numwide_random3)
                        raise ValueError(msg)
                    eids = ints1[:, 0] // 10
                    cids = ints1[:, 1]
                    #nids = ints1[:, 4]
                    assert eids.min() > 0, eids.min()
                    obj.element_node[itotal:itotal2, 0] = repeat(eids, nnodes_expected)
                    ints2 = ints1[:, 4:].reshape(nelements * nnodes_expected, 8)
                    grid_device = ints2[:, 0]#.reshape(nelements, nnodes_expected)

                    #print('%s-grid_device=%s' % (op2.element_name, grid_device))
                    unused_grid_device2 = repeat(grid_device, nnodes_expected)
                    print(f'grid_device = {grid_device.tolist()}')
                    print(f'eids = {eids.tolist()}')
                    print(f'cids = {cids.tolist()}')
                    try:
                        obj.element_node[itotal:itotal2, 1] = grid_device
                    except ValueError:
                        msg = '%s; nnodes=%s\n' % (op2.element_name, nnodes_expected)
                        msg += 'itotal=%s itotal2=%s\n' % (itotal, itotal2)
                        msg += 'grid_device.shape=%s; size=%s\n' % (str(grid_device.shape), grid_device.size)
                        #msg += 'nids=%s' % nids
                        raise ValueError(msg)
                    obj.element_cid[itotal:itotali, 0] = eids
                    obj.element_cid[itotal:itotali, 1] = cids

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, numwide_random3)[:, 4:]

                #(grid_device, oxx, oyy, ozz, txy, tyz, txz, von_mises)
                floats1 = floats.reshape(nelements * nnodes_expected, 8)[:, 1:] # drop grid_device
                assert floats1.shape == (nelements * nnodes_expected, 7), floats1.shape

                obj.data[obj.itime, itotal:itotal2, :] = floats1
                obj.itotal = itotal2
                obj.ielement = itotali
            else:
                if is_vectorized and op2.use_vector and obj.itime == 0:  # pragma: no cover
                    log.debug(f'vectorize CSolid random SORT{op2.sort_method}')
                n = oes_csolid_random3(op2, data, obj, nelements,
                                       element_name, nnodes_expected,
                                       ntotal)
        #elif op2.format_code in [2, 3] and op2.num_wide == 76:  # imag
            # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1081x.op2
            # analysis_code = 5   Frequency
            # table_code    = 905 OESXNO1-OESXNO1C - Cumulative Root Mean Square output
            # format_code   = 3   Magnitude/Phase
            # result_type   = 3   Magnitude/Phase
            # sort_method   = 1
            # random_code   = 0
            # element_type  = 67  CHEXA
            #msg = op2.code_information()
            #return op2._not_implemented_or_skip(data, ndata, msg), None, None
        else:  # pragma: no cover
            raise RuntimeError(op2.code_information() +
                               f'\nnumwide real={numwide_real} imag={numwide_imag} '
                               f'random2={numwide_random2} random3={numwide_random3}')
        return n, nelements, ntotal

    def _oes_csolid2(self, data, ndata, dt, is_magnitude_phase,
                     result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 300 : CHEXA
         - 301 : CPENTA
         - 302 : CTETRA
         - 303 : CPYRAM
        """
        op2 = self.op2
        n = 0
        if op2.is_stress:
            stress_strain = 'stress'
            obj_real = RealSolidStressArrayNx
            #obj_vector_complex = ComplexSolidStressArray
            #obj_vector_random = RandomSolidStressArray
        else:
            obj_real = RealSolidStrainArrayNx
            #obj_vector_complex = ComplexSolidStrainArray
            #obj_vector_random = RandomSolidStrainArray
            stress_strain = 'strain'

        if prefix == '' and postfix == '':
            prefix = stress_strain + '.'

        etype_map = {
            300 : ('chexa', 8, 'CHEXA8'),
            301 : ('cpenta', 6, 'CPENTA6'),
            302 : ('ctetra', 4, 'CTETRA4'),
            303 : ('cpyram', 5, 'CPYRAM5'),
        }
        element_base, nnodes_expected, element_name = etype_map[op2.element_type]
        # chexa_stress
        result_name = prefix + f'{element_base}_{stress_strain}' + postfix

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        numwide_real = 3 + 8 * nnodes_expected
        #numwide_imag = 4 + (17 - 4) * nnodes_expected
        #numwide_random = 4 + (11 - 4) * nnodes_expected
        #numwide_random2 = 18 + 14 * (nnodes_expected - 1)

        #print('nnodes_expected =', nnodes_expected)
        #print('numwide real=%s imag=%s random=%s' % (numwide_real, numwide_imag, numwide_random2))
        op2._data_factor = nnodes_expected
        if op2.format_code == 1 and op2.num_wide == numwide_real:  # real
            ntotal = (12 + 32 * nnodes_expected) * self.factor
            nelements = ndata // ntotal
            #auto_return, is_vectorized = op2._create_oes_object4(
                #nelements, result_name, slot, obj_vector_real)
            #auto_return = op2.read_mode == 1
            #is_vectorized = False
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_real)

            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:  # pragma: no cover
                n = nelements * ntotal
                itotal = obj.ielement
                itotali = obj.itotal + nelements
                itotal2 = obj.itotal + nelements * nnodes_expected
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    # (eid_device, cid, abcd, nnodes)
                    ints = frombuffer(data, dtype=op2.idtype8).copy()
                    try:
                        ints1 = ints.reshape(nelements, numwide_real)
                    except ValueError:
                        msg = 'ints.shape=%s; size=%s ' % (str(ints.shape), ints.size)
                        msg += 'nelements=%s numwide_real=%s nelements*numwide=%s' % (
                            nelements, numwide_real, nelements * numwide_real)
                        raise ValueError(msg)

                    eids = ints1[:, 0] // 10
                    cids = ints1[:, 1]
                    #nids = ints1[:, 4]
                    assert eids.min() > 0, eids.min()
                    obj.element_node[itotal:itotal2, 0] = repeat(eids, nnodes_expected)
                    ints2 = ints1[:, 3:].reshape(nelements * nnodes_expected, 8)
                    grid_device = ints2[:, 0]#.reshape(nelements, nnodes_expected)

                    #print('%s-grid_device=%s' % (op2.element_name, grid_device))
                    unused_grid_device2 = repeat(grid_device, nnodes_expected)
                    try:
                        obj.element_node[itotal:itotal2, 1] = grid_device
                    except ValueError:
                        msg = '%s; nnodes=%s\n' % (op2.element_name, nnodes_expected)
                        msg += 'itotal=%s itotal2=%s\n' % (itotal, itotal2)
                        msg += 'grid_device.shape=%s; size=%s\n' % (str(grid_device.shape), grid_device.size)
                        #msg += 'nids=%s' % nids
                        raise ValueError(msg)
                    obj.element_cid[itotal:itotali, 0] = eids
                    obj.element_cid[itotal:itotali, 1] = cids

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, numwide_real)[:, 3:]
                # 1    2    3    4    5    6    7 - verify...
                #[oxx, oyy, ozz, txy, tyz, txz, ovm]
                #isave = [1, 9, 15, 2, 10, 16, 3, 11, 17, 8]
                #(grid_device,
                    #sxx, sxy, s1, a1, a2, a3, pressure, svm,
                    #syy, syz, s2, b1, b2, b3,
                    #szz, sxz, s3, c1, c2, c3)
                floats1 = floats.reshape(nelements * nnodes_expected, 8)[:, 1:] # drop grid_device

                # o1/o2/o3 is not max/mid/min.  They are not consistently ordered, so we force it.

                obj.data[obj.itime, itotal:itotal2, :] = floats1
                obj.itotal = itotal2
                obj.ielement = itotali
            else:
                n = _oes_csolid2_real(op2, data, n,
                                      obj,
                                      nnodes_expected,
                                      nelements,
                                      element_name,
                                      stress_strain=stress_strain)
        else:  # pragma: no cover
            raise NotImplementedError(op2.code_information())
        assert isinstance(n, int), n
        assert isinstance(ntotal, int), ntotal
        assert isinstance(nelements, int), nelements

        #assert n == ntotal * nelements, f'n={n} ntotal={ntotal*nelements}'
        return n, nelements, ntotal

    def _oes_csolid_composite(self, data, ndata, dt, is_magnitude_phase: bool,
                              result_type: int, prefix: str, postfix: str) -> int:
        """
        306: Nonlinear composite HEXA element (CHEXALN)
        307: Nonlinear composite PENTA element (CPENTALN)

        reads stress/strain for element type:
         - 306 : CHEXALN
         - 307 : CPENTA
         #- 302 : CTETRA
         #- 303 : CPYRAM

        """
        op2 = self.op2
        n = 0
        if op2.is_stress:
            #obj_vector_real = RealSolidStressArray
            #obj_vector_complex = ComplexSolidStressArray
            #obj_vector_random = RandomSolidStressArray
            word = 'stress'
            if op2.element_type == 306:  # CHEXALN
                nedges = 4 # quad
                nnodes_expected = 8
                result_name = prefix + 'chexa_stress' + postfix
                element_name = 'CHEXA8'
                # real=67
            elif op2.element_type == 307:  # CPENTALN
                nedges = 3 # tri
                nnodes_expected = 6
                result_name = prefix + 'cpenta_stress' + postfix
                element_name = 'CPENTA6'
            #elif op2.element_type == 302:  # CTETRA
                #nnodes_expected = 4
                #result_name = prefix + 'ctetra_stress' + postfix
                #element_name = 'CTETRA4'
            #elif op2.element_type == 303:  # CPYRAM
                #nnodes_expected = 5
                #result_name = prefix + 'cpyram_stress' + postfix
                #element_name = 'CPYRAM5'
            else:  # pragma: no cover
                raise RuntimeError(op2.code_information())
        else:
            #obj_vector_real = RealSolidStrainArray
            #obj_vector_complex = ComplexSolidStrainArray
            #obj_vector_random = RandomSolidStrainArray
            word = 'strain'
            if op2.element_type == 306:  # CHEXALN
                nedges = 4 # quad
                nnodes_expected = 8
                result_name = prefix + 'chexa_strain' + postfix
                element_name = 'CHEXA8'
            elif op2.element_type == 307:  # CPENTA
                nedges = 3 # tri
                nnodes_expected = 6
                result_name = prefix + 'cpenta_strain' + postfix
                element_name = 'CPENTA6'
            #elif op2.element_type == 302:  # CTETRA
                #nnodes_expected = 4
                #result_name = prefix + 'ctetra_strain' + postfix
                #element_name = 'CTETRA4'
            #elif op2.element_type == 303:  # CPYRAM
                #nnodes_expected = 5
                #result_name = prefix + 'cpyram_strain' + postfix
                #element_name = 'CPYRAM5'
            else:  # pragma: no cover
                raise NotImplementedError(op2.code_information())
                #msg = 'sort1 Type=%s num=%s' % (op2.element_name, op2.element_type)
                #return op2._not_implemented_or_skip(data, ndata, msg)

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        numwide_real = 3 + 8 * nedges # 3 + 8*4 = 35
        numwide_imag = 3 + 7 * 14 # 3 + 7 * 14 = 101
        #print(op2.code_information())
        #print(f'{op2.element_name} numwide_real={numwide_real} numwide_imag={numwide_imag} -> {op2.num_wide}')
        #numwide_real = 3 + 8 * nnodes_expected
        #numwide_imag = 4 + (17 - 4) * nnodes_expected
        #numwide_random = 4 + (11 - 4) * nnodes_expected
        #numwide_random2 = 18 + 14 * (nnodes_expected - 1)
        preline1 = '%s-%s' % (op2.element_name, op2.element_type)
        preline2 = ' ' * len(preline1)

        #print('nnodes_expected =', nnodes_expected)
        #print('numwide real=%s imag=%s random=%s' % (numwide_real, numwide_imag, numwide_random2))
        op2._data_factor = nedges
        if result_type == 0 and op2.num_wide == numwide_real:  # real
            op2.log.warning(f'skipping {op2.table_name_str}: {op2.element_name}-{op2.element_type} {word} csolid composite')
            ntotal = 12 + 32 * nedges
            nelements = ndata // ntotal
            #auto_return, is_vectorized = op2._create_oes_object4(
                #nelements, result_name, slot, obj_vector_real)
            auto_return = op2.read_mode == 1
            is_vectorized = False
            if auto_return:
                assert ntotal == op2.num_wide * 4
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1 and 0:  # pragma: no cover
                n = nelements * 4 * op2.num_wide
                itotal = obj.ielement
                itotali = obj.itotal + nelements
                itotal2 = obj.itotal + nelements * nnodes_expected
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    # (eid_device, cid, abcd, nnodes)
                    ints = frombuffer(data, dtype=op2.idtype).copy()
                    try:
                        ints1 = ints.reshape(nelements, numwide_real)
                    except ValueError:
                        msg = 'ints.shape=%s; size=%s ' % (str(ints.shape), ints.size)
                        msg += 'nelements=%s numwide_real=%s nelements*numwide=%s' % (
                            nelements, numwide_real, nelements * numwide_real)
                        raise ValueError(msg)
                    eids = ints1[:, 0] // 10
                    cids = ints1[:, 1]
                    #nids = ints1[:, 4]
                    assert eids.min() > 0, eids.min()
                    obj.element_node[itotal:itotal2, 0] = repeat(eids, nnodes_expected)
                    ints2 = ints1[:, 4:].reshape(nelements * nnodes_expected, 21)
                    grid_device = ints2[:, 0]#.reshape(nelements, nnodes_expected)

                    #print('%s-grid_device=%s' % (op2.element_name, grid_device))
                    unused_grid_device2 = repeat(grid_device, nnodes_expected)
                    try:
                        obj.element_node[itotal:itotal2, 1] = grid_device
                    except ValueError:
                        msg = '%s; nnodes=%s\n' % (op2.element_name, nnodes_expected)
                        msg += 'itotal=%s itotal2=%s\n' % (itotal, itotal2)
                        msg += 'grid_device.shape=%s; size=%s\n' % (str(grid_device.shape), grid_device.size)
                        #msg += 'nids=%s' % nids
                        raise ValueError(msg)
                    obj.element_cid[itotal:itotali, 0] = eids
                    obj.element_cid[itotal:itotali, 1] = cids

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, numwide_real)[:, 4:]
                # 1     9    15   2    10   16  3   11  17   8
                #[oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, ovm]
                #isave = [1, 9, 15, 2, 10, 16, 3, 11, 17, 8]
                #(grid_device,
                    #sxx, sxy, s1, a1, a2, a3, pressure, svm,
                    #syy, syz, s2, b1, b2, b3,
                    #szz, sxz, s3, c1, c2, c3)
                floats1 = floats.reshape(nelements * nnodes_expected, 21)#[:, 1:] # drop grid_device

                # o1/o2/o3 is not max/mid/min.  They are not consistently ordered, so we force it.
                max_mid_min = np.vstack([
                    floats1[:, 3],
                    floats1[:, 11],
                    floats1[:, 17],
                ]).T
                max_mid_min.sort(axis=1)
                assert max_mid_min.shape == (nelements * nnodes_expected, 3), max_mid_min.shape
                obj.data[obj.itime, itotal:itotal2, 6:9] = max_mid_min[:, [2, 1, 0]]

                #obj.data[obj.itime, itotal:itotal2, :] = floats1[:, isave]
                obj.data[obj.itime, itotal:itotal2, :6] = floats1[:, [1, 9, 15, 2, 10, 16]]
                obj.data[obj.itime, itotal:itotal2, 9] = floats1[:, 8]
                obj.itotal = itotal2
                obj.ielement = itotali
            else:
                #if is_vectorized and op2.use_vector:  # pragma: no cover
                    #op2.log.debug('vectorize CSolid real SORT%s' % op2.sort_method)
                n = oes_csolid_composite_real(op2, data, obj,
                                              nelements, nedges,
                                              element_name, preline1, preline2, dt)

        elif result_type == 1 and op2.num_wide == numwide_imag:  # complex
            # 1 PLY I Lamina number
            # 2 FLOC I Fiber location (BOT, MID, TOP)
            #
            # 3 GRID I Edge grid ID
            # 4 EX1R RS Normal strain in the 1-direction
            # 5 EY1R RS Normal strain in the 2-direction
            # 6 EZ1R RS Normal strain in the 3-direction
            # 7 ET1R RS Shear strain in the 12-plane
            # 8 EL2R RS Shear strain in the 23-plane
            # 9 EL1R RS Shear strain in the 13-plane
            # 10 EX1I RS Normal strain in the 1-direction
            # 11 EY1I RS Normal strain in the 2-direction
            # 12 EZ1I RS Normal strain in the 3-direction
            # 13 ET1I RS Shear strain in the 12-plane
            # 14 EL2I RS Shear strain in the 23-plane
            # 15 EL1I RS Shear strain in the 13-plane
            # 16 ETMAX1 RS von Mises strain
            # For each fiber location requested (PLSLOC), words 3 through 16 repeat 4 times.
            if op2.read_mode == 1:
                return ndata, None, None
            self.show_data(data[n:n+4*op2.num_wide])
            aaa
        elif result_type == 0 and op2.num_wide == 68:  # real
            msg = (f'etype={op2.element_name} ({op2.element_type}) '
                   f'{op2.table_name_str}-OES-CSOLID-random-numwide={op2.num_wide} '
                   f'numwide_real=11 numwide_imag=9 result_type={result_type}')
            return op2._not_implemented_or_skip(data, ndata, msg), None, None
        elif op2.element_type == 306 and op2.num_wide in {52, 101}:
             #HEXALN
            msg = (f'etype={op2.element_name} ({op2.element_type}) '
                   f'{op2.table_name_str}-OES-CSOLID-real-numwide={op2.num_wide} '
                   f'numwide_real={numwide_real} numwide_imag={numwide_imag} result_type={result_type}')
            return op2._not_implemented_or_skip(data, ndata, msg), None, None
        elif op2.element_type == 307 and op2.num_wide == 52:
            # CPENTA
            msg = (f'etype={op2.element_name} ({op2.element_type}) '
                   f'{op2.table_name_str}-OES-CSOLID-real-numwide={op2.num_wide} '
                   f'numwide_real={numwide_real} numwide_imag={numwide_imag} result_type={result_type}')
            return op2._not_implemented_or_skip(data, ndata, msg), None, None
        else:  # pragma: no cover
            raise NotImplementedError(op2.code_information())
        assert n == ntotal * nelements, f'n={n} ntotal={ntotal*nelements}'
        return n, nelements, ntotal

    def _oes_csolid_linear_hyperelastic_cosine(self, data, ndata, dt, unused_is_magnitude_phase,
                                               result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 140 :CHEXAFD

        """
        op2 = self.op2
        n = 0
        log = op2.log
        if op2.is_stress:
            #obj_vector_real = RealSolidStressArray
            #obj_vector_complex = ComplexSolidStressArray
            #obj_vector_random = RandomSolidStressArray
            word = 'stress'
            prefix = word + '.'
            if op2.element_type == 140:  # CHEXA
                nnodes_expected = 8
                result_name = prefix + 'chexa_stress' + postfix
                #element_name = 'CHEXA8'
                # real=122
            #elif op2.element_type == 160:  # CPENTA
                #nnodes_expected = 6
                #result_name = prefix + 'cpenta_stress' + postfix
                #element_name = 'CPENTA6'
            #elif op2.element_type == 165:  # CPENTA
                #nnodes_expected = 21
                #result_name = prefix + 'cpenta_stress' + postfix
                #element_name = 'CPENTA6'

            #elif op2.element_type == 161:  # CTETRA
                #nnodes_expected = 1
                #result_name = prefix + 'ctetra_stress' + postfix
                #element_name = 'CTETRA4'
            #elif op2.element_type == 166:  # CTETRA
                #nnodes_expected = 5
                #result_name = prefix + 'ctetra_stress' + postfix
                #element_name = 'CTETRA4'
            #elif op2.element_type == 303:  # CPYRAM
                #nnodes_expected = 5
                #result_name = prefix + 'cpyram_stress' + postfix
                #element_name = 'CPYRAM5'
            else:  # pragma: no cover
                raise RuntimeError(op2.code_information())
        else:
            #obj_vector_real = RealSolidStrainArray
            #obj_vector_complex = ComplexSolidStrainArray
            #obj_vector_random = RandomSolidStrainArray
            word = 'strain'
            #if op2.element_type == 202:  # CHEXA
                #nnodes_expected = 8
                #result_name = prefix + 'chexa_strain' + postfix
                #element_name = 'CHEXA8'
            #elif op2.element_type == 301:  # CPENTA
                #nnodes_expected = 6
                #result_name = prefix + 'cpenta_strain' + postfix
                #element_name = 'CPENTA6'
            #elif op2.element_type == 302:  # CTETRA
                #nnodes_expected = 4
                #result_name = prefix + 'ctetra_strain' + postfix
                #element_name = 'CTETRA4'
            #elif op2.element_type == 303:  # CPYRAM
                #nnodes_expected = 5
                #result_name = prefix + 'cpyram_strain' + postfix
                #element_name = 'CPYRAM5'
            #else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
                #msg = 'sort1 Type=%s num=%s' % (op2.element_name, op2.element_type)
                #return op2._not_implemented_or_skip(data, ndata, msg)

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        numwide_real = 2 + 20 * nnodes_expected
        #numwide_real = 162  # CHEXA
        #print(op2.num_wide, numwide_real)
        assert numwide_real == op2.num_wide, numwide_real
        #numwide_imag = 4 + (17 - 4) * nnodes_expected
        #numwide_random = 4 + (11 - 4) * nnodes_expected
        #numwide_random2 = 18 + 14 * (nnodes_expected - 1)
        preline1 = '%s-%s' % (op2.element_name, op2.element_type)
        preline2 = ' ' * len(preline1)

        #print('nnodes_expected =', nnodes_expected)
        #print('numwide real=%s imag=%s random=%s' % (numwide_real, numwide_imag, numwide_random2))
        op2._data_factor = nnodes_expected

        if op2.format_code == 1 and op2.num_wide == numwide_real:  # real
            ntotal = 8 + 80 * nnodes_expected
            #ntotal = numwide_real * 4
            nelements = ndata // ntotal
            assert ndata % ntotal == 0
            #auto_return, is_vectorized = op2._create_oes_object4(
                #nelements, result_name, slot, obj_vector_real)
            auto_return = op2.read_mode == 1
            is_vectorized = False
            if auto_return:
                assert ntotal == op2.num_wide * 4
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1 and 0:  # pragma: no cover
                n = nelements * 4 * op2.num_wide
                itotal = obj.ielement
                itotali = obj.itotal + nelements
                itotal2 = obj.itotal + nelements * nnodes_expected
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    # (eid_device, cid, abcd, nnodes)
                    ints = frombuffer(data, dtype=op2.idtype).copy()
                    try:
                        ints1 = ints.reshape(nelements, numwide_real)
                    except ValueError:
                        msg = 'ints.shape=%s; size=%s ' % (str(ints.shape), ints.size)
                        msg += 'nelements=%s numwide_real=%s nelements*numwide=%s' % (
                            nelements, numwide_real, nelements * numwide_real)
                        raise ValueError(msg)
                    eids = ints1[:, 0] // 10
                    cids = ints1[:, 1]
                    #nids = ints1[:, 4]
                    assert eids.min() > 0, eids.min()
                    obj.element_node[itotal:itotal2, 0] = repeat(eids, nnodes_expected)
                    ints2 = ints1[:, 4:].reshape(nelements * nnodes_expected, 21)
                    grid_device = ints2[:, 0]#.reshape(nelements, nnodes_expected)

                    #print('%s-grid_device=%s' % (op2.element_name, grid_device))
                    unused_grid_device2 = repeat(grid_device, nnodes_expected)
                    try:
                        obj.element_node[itotal:itotal2, 1] = grid_device
                    except ValueError:
                        msg = '%s; nnodes=%s\n' % (op2.element_name, nnodes_expected)
                        msg += 'itotal=%s itotal2=%s\n' % (itotal, itotal2)
                        msg += 'grid_device.shape=%s; size=%s\n' % (str(grid_device.shape), grid_device.size)
                        #msg += 'nids=%s' % nids
                        raise ValueError(msg)
                    obj.element_cid[itotal:itotali, 0] = eids
                    obj.element_cid[itotal:itotali, 1] = cids

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, numwide_real)[:, 4:]
                # 1     9    15   2    10   16  3   11  17   8
                #[oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, ovm]
                #isave = [1, 9, 15, 2, 10, 16, 3, 11, 17, 8]
                #(grid_device,
                    #sxx, sxy, s1, a1, a2, a3, pressure, svm,
                    #syy, syz, s2, b1, b2, b3,
                    #szz, sxz, s3, c1, c2, c3)
                floats1 = floats.reshape(nelements * nnodes_expected, 21)#[:, 1:] # drop grid_device

                # o1/o2/o3 is not max/mid/min.  They are not consistently ordered, so we force it.
                max_mid_min = np.vstack([
                    floats1[:, 3],
                    floats1[:, 11],
                    floats1[:, 17],
                ]).T
                max_mid_min.sort(axis=1)
                assert max_mid_min.shape == (nelements * nnodes_expected, 3), max_mid_min.shape
                obj.data[obj.itime, itotal:itotal2, 6:9] = max_mid_min[:, [2, 1, 0]]

                #obj.data[obj.itime, itotal:itotal2, :] = floats1[:, isave]
                obj.data[obj.itime, itotal:itotal2, :6] = floats1[:, [1, 9, 15, 2, 10, 16]]
                obj.data[obj.itime, itotal:itotal2, 9] = floats1[:, 8]
                obj.itotal = itotal2
                obj.ielement = itotali
            else:
                #if is_vectorized and op2.use_vector:  # pragma: no cover
                    #op2.log.debug('vectorize CSolid real SORT%s' % op2.sort_method)
                n = oes_csolid_linear_hyperelastic_cosine_real(
                    op2, data,
                    nelements, nnodes_expected,
                    preline1, preline2)
            log.warning(f'skipping {op2.table_name_str}: {op2.element_name}-{op2.element_type} linear hyperelastic cosine {word}')
            return n, None, None
        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
        assert n == ntotal * nelements, f'n={n} ntotal={ntotal*nelements}'
        return n, nelements, ntotal

    def _oes_csolid_linear_hyperelastic(self, data, ndata, dt, is_magnitude_phase,
                                        result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 160 :CPENTAFD
         - 163 : CHEXAFD
         - 166 : CTETRAFD

        # centroid??
         - 161 : CTETRAFD

        # more nodes???
         - 165 : CPENTAFD
        """
        op2 = self.op2
        n = 0
        if op2.is_stress:
            #obj_vector_real = RealSolidStressArray
            #obj_vector_complex = ComplexSolidStressArray
            #obj_vector_random = RandomSolidStressArray
            word = 'stress'
            if op2.element_type == 163:  # CHEXA
                nnodes_expected = 27
                result_name = prefix + 'chexa_stress' + postfix
                #element_name = 'CHEXA8'
                # real=122
            elif op2.element_type == 160:  # CPENTA
                nnodes_expected = 6
                result_name = prefix + 'cpenta_stress' + postfix
                #element_name = 'CPENTA6'
            elif op2.element_type == 165:  # CPENTA
                nnodes_expected = 21
                result_name = prefix + 'cpenta_stress' + postfix
                #element_name = 'CPENTA6'

            elif op2.element_type == 161:  # CTETRA
                nnodes_expected = 1
                result_name = prefix + 'ctetra_stress' + postfix
                #element_name = 'CTETRA4'
            elif op2.element_type == 166:  # CTETRA
                nnodes_expected = 5
                result_name = prefix + 'ctetra_stress' + postfix
                #element_name = 'CTETRA4'
            #elif op2.element_type == 303:  # CPYRAM
                #nnodes_expected = 5
                #result_name = prefix + 'cpyram_stress' + postfix
                #element_name = 'CPYRAM5'
            else:  # pragma: no cover
                raise RuntimeError(op2.code_information())
        else:
            #obj_vector_real = RealSolidStrainArray
            #obj_vector_complex = ComplexSolidStrainArray
            #obj_vector_random = RandomSolidStrainArray
            word = 'strain'
            #if op2.element_type == 202:  # CHEXA
                #nnodes_expected = 8
                #result_name = prefix + 'chexa_strain' + postfix
                #element_name = 'CHEXA8'
            #elif op2.element_type == 301:  # CPENTA
                #nnodes_expected = 6
                #result_name = prefix + 'cpenta_strain' + postfix
                #element_name = 'CPENTA6'
            #elif op2.element_type == 302:  # CTETRA
                #nnodes_expected = 4
                #result_name = prefix + 'ctetra_strain' + postfix
                #element_name = 'CTETRA4'
            #elif op2.element_type == 303:  # CPYRAM
                #nnodes_expected = 5
                #result_name = prefix + 'cpyram_strain' + postfix
                #element_name = 'CPYRAM5'
            #else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
                #msg = 'sort1 Type=%s num=%s' % (op2.element_name, op2.element_type)
                #return op2._not_implemented_or_skip(data, ndata, msg)

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        numwide_real = 2 + 20 * nnodes_expected
        #numwide_real = 122  # CHEXA
        #print(op2.num_wide, numwide_real)
        assert numwide_real == op2.num_wide, numwide_real
        numwide_imag = 4 + (17 - 4) * nnodes_expected
        #numwide_random = 4 + (11 - 4) * nnodes_expected
        numwide_random2 = 18 + 14 * (nnodes_expected - 1)
        preline1 = '%s-%s' % (op2.element_name, op2.element_type)
        preline2 = ' ' * len(preline1)

        #print('nnodes_expected =', nnodes_expected)
        #print('numwide real=%s imag=%s random=%s' % (numwide_real, numwide_imag, numwide_random2))
        op2._data_factor = nnodes_expected

        if op2.format_code == 1 and op2.num_wide == numwide_real:  # real
            ntotal = 8 + 80 * nnodes_expected
            #ntotal = numwide_real * 4
            nelements = ndata // ntotal
            assert ndata % ntotal == 0
            #auto_return, is_vectorized = op2._create_oes_object4(
                #nelements, result_name, slot, obj_vector_real)
            auto_return = op2.read_mode == 1
            is_vectorized = False
            if auto_return:
                assert ntotal == op2.num_wide * 4
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1 and 0:  # pragma: no cover
                n = nelements * 4 * op2.num_wide
                itotal = obj.ielement
                itotali = obj.itotal + nelements
                itotal2 = obj.itotal + nelements * nnodes_expected
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    # (eid_device, cid, abcd, nnodes)
                    ints = frombuffer(data, dtype=op2.idtype).copy()
                    try:
                        ints1 = ints.reshape(nelements, numwide_real)
                    except ValueError:
                        msg = 'ints.shape=%s; size=%s ' % (str(ints.shape), ints.size)
                        msg += 'nelements=%s numwide_real=%s nelements*numwide=%s' % (
                            nelements, numwide_real, nelements * numwide_real)
                        raise ValueError(msg)
                    eids = ints1[:, 0] // 10
                    cids = ints1[:, 1]
                    #nids = ints1[:, 4]
                    assert eids.min() > 0, eids.min()
                    obj.element_node[itotal:itotal2, 0] = repeat(eids, nnodes_expected)
                    ints2 = ints1[:, 4:].reshape(nelements * nnodes_expected, 21)
                    grid_device = ints2[:, 0]#.reshape(nelements, nnodes_expected)

                    #print('%s-grid_device=%s' % (op2.element_name, grid_device))
                    unused_grid_device2 = repeat(grid_device, nnodes_expected)
                    try:
                        obj.element_node[itotal:itotal2, 1] = grid_device
                    except ValueError:
                        msg = '%s; nnodes=%s\n' % (op2.element_name, nnodes_expected)
                        msg += 'itotal=%s itotal2=%s\n' % (itotal, itotal2)
                        msg += 'grid_device.shape=%s; size=%s\n' % (str(grid_device.shape), grid_device.size)
                        #msg += 'nids=%s' % nids
                        raise ValueError(msg)
                    obj.element_cid[itotal:itotali, 0] = eids
                    obj.element_cid[itotal:itotali, 1] = cids

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, numwide_real)[:, 4:]
                # 1     9    15   2    10   16  3   11  17   8
                #[oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, ovm]
                #isave = [1, 9, 15, 2, 10, 16, 3, 11, 17, 8]
                #(grid_device,
                    #sxx, sxy, s1, a1, a2, a3, pressure, svm,
                    #syy, syz, s2, b1, b2, b3,
                    #szz, sxz, s3, c1, c2, c3)
                floats1 = floats.reshape(nelements * nnodes_expected, 21)#[:, 1:] # drop grid_device

                # o1/o2/o3 is not max/mid/min.  They are not consistently ordered, so we force it.
                max_mid_min = np.vstack([
                    floats1[:, 3],
                    floats1[:, 11],
                    floats1[:, 17],
                ]).T
                max_mid_min.sort(axis=1)
                assert max_mid_min.shape == (nelements * nnodes_expected, 3), max_mid_min.shape
                obj.data[obj.itime, itotal:itotal2, 6:9] = max_mid_min[:, [2, 1, 0]]

                #obj.data[obj.itime, itotal:itotal2, :] = floats1[:, isave]
                obj.data[obj.itime, itotal:itotal2, :6] = floats1[:, [1, 9, 15, 2, 10, 16]]
                obj.data[obj.itime, itotal:itotal2, 9] = floats1[:, 8]
                obj.itotal = itotal2
                obj.ielement = itotali
            else:
                #if is_vectorized and op2.use_vector:  # pragma: no cover
                    #op2.log.debug('vectorize CSolid real SORT%s' % op2.sort_method)
                n = oes_csolid_linear_hyperelastic_real(op2, data, obj, nelements,
                                                        nnodes_expected,
                                                        preline1, preline2)

            op2.log.warning(f'skipping {op2.table_name_str}: {op2.element_name}-{op2.element_type} linear hyperelastic {word}')
        else:  # pragma: no cover
            raise RuntimeError(op2.code_information() +
                               '\nnumwide real=%s imag=%s random=%s' % (
                                   numwide_real, numwide_imag, numwide_random2))
        assert n == ntotal * nelements, f'n={n} ntotal={ntotal*nelements}'
        return n, nelements, ntotal

    def _oes_csolid_nonlinear_hyperelastic(self, data, ndata, dt, is_magnitude_phase,
                                           result_type: str, prefix: str, postfix: str):
        """
        reads stress/strain for element type:
         - 202 : CHEXAFD
         - 204 : PENTA6FD

            '   N O N L I N E A R   S T R E S S E S   I N   H Y P E R E L A S T I C   H E X A H E D R O N   E L E M E N T S  ( HEXA8FD )'
            ' '
            '  ELEMENT GRID/   POINT                       CAUCHY STRESSES/ LOG STRAINS                        PRESSURE    VOL. STRAIN'
            '     ID   GAUSS     ID       X           Y           Z           XY          YZ          ZX'
            '0     401 GRID      401  1.9128E+03  6.2729E+02 -3.4828E+02 -7.5176E+01  7.8259E+00 -2.5001E+02  7.3060E+02  7.3060E-03'
            '                         6.8270E-01 -6.5437E-04 -1.2874E+00 -3.9645E-02 -2.9882E-03 -5.9975E-02'
            '                    402  1.1024E+03  1.0686E+03  2.0832E+01 -1.7936E+00 -2.3656E-01 -1.1467E+02  7.3060E+02  7.3060E-03'
            '                         6.8201E-01  6.4335E-01 -1.2964E+00 -2.7195E-03 -1.0809E-02  6.3608E-02'
            '                   1402  1.1024E+03  1.0686E+03  2.0832E+01  1.7936E+00  2.3656E-01 -1.1467E+02  7.3060E+02  7.3060E-03'
            '                         6.8201E-01  6.4335E-01 -1.2964E+00 -2.7195E-03 -1.0809E-02  6.3608E-02'
            '                   1401  1.9128E+03  6.2729E+02 -3.4828E+02  7.5176E+01 -7.8259E+00 -2.5001E+02  7.3060E+02  7.3060E-03'
            '                         6.8201E-01  6.4335E-01 -1.2964E+00 -2.7195E-03 -1.0809E-02  6.3608E-02'
            '                    501  1.9159E+03  6.2332E+02 -3.4744E+02 -7.5730E+01  7.9009E+00 -2.5075E+02  7.3060E+02  7.3060E-03'
            '                         6.8201E-01  6.4335E-01 -1.2964E+00 -2.7195E-03 -1.0809E-02  6.3608E-02'
            '                    502  1.1004E+03  1.0667E+03  2.4631E+01 -1.7898E+00 -2.2971E-01 -1.1434E+02  7.3060E+02  7.3060E-03'
            '                         6.8201E-01  6.4335E-01 -1.2964E+00 -2.7195E-03 -1.0809E-02  6.3608E-02'
            '                   1502  1.1004E+03  1.0667E+03  2.4631E+01  1.7898E+00  2.2971E-01 -1.1434E+02  7.3060E+02  7.3060E-03'
            '                         6.8201E-01  6.4335E-01 -1.2964E+00 -2.7195E-03 -1.0809E-02  6.3608E-02'
            '                   1501  1.9159E+03  6.2332E+02 -3.4744E+02  7.5730E+01 -7.9009E+00 -2.5075E+02  7.3060E+02  7.3060E-03'
            '                         6.8201E-01  6.4335E-01 -1.2964E+00 -2.7195E-03 -1.0809E-02  6.3608E-02'

          # 216 TETRAFD
          # 218 HEXAFD
          # 220 (CPENTA)
          # 221 (CTETRA)
          ELEMENT GRID/   POINT                         STRESSES/ TOTAL STRAINS                          EQUIVALENT EFF. STRAIN  EFF. CREEP
     ID   GAUSS     ID       X           Y           Z           XY          YZ          ZX        STRESS   PLAS/NLELAS   STRAIN

        """
        op2 = self.op2
        n = 0
        if op2.is_stress:
            #obj_vector_real = RealSolidStressArray
            #obj_vector_complex = ComplexSolidStressArray
            #obj_vector_random = RandomSolidStressArray
            word = 'stress'
            if op2.element_type in [202, 218]:  # CHEXA
                nnodes_expected = 8
                result_name = prefix + 'chexa_stress_strain' + postfix
                element_name = 'CHEXA8'
                # real=122
            elif op2.element_type in [204, 220]:  # CPENTA
                nnodes_expected = 6
                result_name = prefix + 'cpenta_stress_strain' + postfix
                element_name = 'CPENTA6'
            elif op2.element_type in [216, 221]:  # CTETRA
                nnodes_expected = 4
                result_name = prefix + 'ctetra_stress_strain' + postfix
                element_name = 'CTETRA4'
            #elif op2.element_type == 303:  # CPYRAM
                #nnodes_expected = 5
                #result_name = prefix + 'cpyram_stress' + postfix
                #element_name = 'CPYRAM5'
            else:  # pragma: no cover
                raise RuntimeError(op2.code_information())
        else:
            #obj_vector_real = RealSolidStrainArray
            #obj_vector_complex = ComplexSolidStrainArray
            #obj_vector_random = RandomSolidStrainArray
            word = 'strain'
            #if op2.element_type == 202:  # CHEXA
                #nnodes_expected = 8
                #result_name = prefix + 'chexa_strain' + postfix
                #element_name = 'CHEXA8'
            #elif op2.element_type == 301:  # CPENTA
                #nnodes_expected = 6
                #result_name = prefix + 'cpenta_strain' + postfix
                #element_name = 'CPENTA6'
            #elif op2.element_type == 302:  # CTETRA
                #nnodes_expected = 4
                #result_name = prefix + 'ctetra_strain' + postfix
                #element_name = 'CTETRA4'
            #elif op2.element_type == 303:  # CPYRAM
                #nnodes_expected = 5
                #result_name = prefix + 'cpyram_strain' + postfix
                #element_name = 'CPYRAM5'
            #else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
                #msg = 'sort1 Type=%s num=%s' % (op2.element_name, op2.element_type)
                #return op2._not_implemented_or_skip(data, ndata, msg)

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        numwide_real = 2 + 15 * nnodes_expected
        #numwide_real = 122  # CHEXA
        assert numwide_real == op2.num_wide, numwide_real
        numwide_imag = 4 + (17 - 4) * nnodes_expected
        #numwide_random = 4 + (11 - 4) * nnodes_expected
        numwide_random2 = 18 + 14 * (nnodes_expected - 1)
        preline1 = '%s-%s' % (op2.element_name, op2.element_type)
        preline2 = ' ' * len(preline1)

        #print('nnodes_expected =', nnodes_expected)
        #print('numwide real=%s imag=%s random=%s' % (numwide_real, numwide_imag, numwide_random2))
        op2._data_factor = nnodes_expected

        if result_type == 0 and op2.num_wide == numwide_real:  # real
            ntotal = 8 + 60 * nnodes_expected
            #ntotal = numwide_real * 4
            nelements = ndata // ntotal
            assert ndata % ntotal == 0
            #auto_return, is_vectorized = op2._create_oes_object4(
                #nelements, result_name, slot, obj_vector_real)
            auto_return = op2.read_mode == 1
            is_vectorized = False
            if auto_return:
                assert ntotal == op2.num_wide * 4
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1 and 0:  # pragma: no cover
                n = nelements * 4 * op2.num_wide
                itotal = obj.ielement
                itotali = obj.itotal + nelements
                itotal2 = obj.itotal + nelements * nnodes_expected
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    # (eid_device, cid, abcd, nnodes)
                    ints = frombuffer(data, dtype=op2.idtype).copy()
                    try:
                        ints1 = ints.reshape(nelements, numwide_real)
                    except ValueError:
                        msg = 'ints.shape=%s; size=%s ' % (str(ints.shape), ints.size)
                        msg += 'nelements=%s numwide_real=%s nelements*numwide=%s' % (
                            nelements, numwide_real, nelements * numwide_real)
                        raise ValueError(msg)
                    eids = ints1[:, 0] // 10
                    cids = ints1[:, 1]
                    #nids = ints1[:, 4]
                    assert eids.min() > 0, eids.min()
                    obj.element_node[itotal:itotal2, 0] = repeat(eids, nnodes_expected)
                    ints2 = ints1[:, 4:].reshape(nelements * nnodes_expected, 21)
                    grid_device = ints2[:, 0]#.reshape(nelements, nnodes_expected)

                    #print('%s-grid_device=%s' % (op2.element_name, grid_device))
                    unused_grid_device2 = repeat(grid_device, nnodes_expected)
                    try:
                        obj.element_node[itotal:itotal2, 1] = grid_device
                    except ValueError:
                        msg = '%s; nnodes=%s\n' % (op2.element_name, nnodes_expected)
                        msg += 'itotal=%s itotal2=%s\n' % (itotal, itotal2)
                        msg += 'grid_device.shape=%s; size=%s\n' % (str(grid_device.shape), grid_device.size)
                        #msg += 'nids=%s' % nids
                        raise ValueError(msg)
                    obj.element_cid[itotal:itotali, 0] = eids
                    obj.element_cid[itotal:itotali, 1] = cids

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, numwide_real)[:, 4:]
                # 1     9    15   2    10   16  3   11  17   8
                #[oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, ovm]
                #isave = [1, 9, 15, 2, 10, 16, 3, 11, 17, 8]
                #(grid_device,
                    #sxx, sxy, s1, a1, a2, a3, pressure, svm,
                    #syy, syz, s2, b1, b2, b3,
                    #szz, sxz, s3, c1, c2, c3)
                floats1 = floats.reshape(nelements * nnodes_expected, 21)#[:, 1:] # drop grid_device

                # o1/o2/o3 is not max/mid/min.  They are not consistently ordered, so we force it.
                max_mid_min = np.vstack([
                    floats1[:, 3],
                    floats1[:, 11],
                    floats1[:, 17],
                ]).T
                max_mid_min.sort(axis=1)
                assert max_mid_min.shape == (nelements * nnodes_expected, 3), max_mid_min.shape
                obj.data[obj.itime, itotal:itotal2, 6:9] = max_mid_min[:, [2, 1, 0]]

                #obj.data[obj.itime, itotal:itotal2, :] = floats1[:, isave]
                obj.data[obj.itime, itotal:itotal2, :6] = floats1[:, [1, 9, 15, 2, 10, 16]]
                obj.data[obj.itime, itotal:itotal2, 9] = floats1[:, 8]
                obj.itotal = itotal2
                obj.ielement = itotali
            else:
                #if is_vectorized and op2.use_vector:  # pragma: no cover
                    #op2.log.debug('vectorize CSolid real SORT%s' % op2.sort_method)

                n = oes_csolid_nonlinear_hyperelastic_real(
                    op2, data, obj, nnodes_expected, nelements, ntotal,
                    element_name, preline1, preline2)
            op2.log.warning(f'skipping {op2.table_name_str}: {op2.element_name}-{op2.element_type} nonlinear hyperelastic {word}')
            return n, None, None
        else:  # pragma: no cover
            raise RuntimeError(op2.code_information() +
                               '\nnumwide real=%s imag=%s random=%s' % (
                                   numwide_real, numwide_imag, numwide_random2))
        assert n == ntotal * nelements, f'n={n} ntotal={ntotal*nelements}'
        return n, nelements, ntotal

    def _oes_csolid_nonlinear(self, data, ndata, dt, unused_is_magnitude_phase,
                              result_type, prefix, postfix):
        """
        reads stress/strain for element type:
        - 85-TETRANL
        - 91-PENTANL
        - 93-HEXANL
        - 256-PYRAMNL

        2  CTYPE CHAR4
        3  NODEF 1 Number of active GRID points

        4  GRID  I Grid / Gauss
        5  SX    RS Stress in x
        6  SY    RS Stress in y
        7  SZ    RS Stress in z
        8  SXY   RS Stress in xy
        9  SYZ   RS Stress in yz
        10 SZX   RS Stress in zx
        11 SE    RS Equivalent stress
        12 EPS   RS Effective plastic strain
        13 ECS   RS Effective creep strain
        14 EX    RS Strain in x
        15 EY    RS Strain in y
        16 EZ    RS Strain in z
        17 EXY   RS Strain in xy
        18 EYZ   RS Strain in yz
        19 EZX   RS Strain in zx
        Words 3 through 19 repeat 005 times
        """
        op2 = self.op2
        #real
        #85:  2 + (18 - 2) * 5,  # Nonlinear CTETRA
        #256: 4 + (18 - 2) * 6,  # Nonlinear CHEXA -> ???

        # random
        #91:  4 + (25 - 4) * 7,  # Nonlinear CPENTA
        #93:  4 + (25 - 4) * 9,  # Nonlinear CHEXA -> 584 (can cause a crash)
        #256: 4 + (25 - 4) * 6,  # Nonlinear CHEXA -> ???

        # the nodes are nnodes + 1
        if op2.element_type == 85:
            etype = 'CTETRANL'
            nnodes = 5
            result_name = prefix + 'ctetra_stress_strain' + postfix
        elif op2.element_type == 91:
            etype = 'CPENTANL'
            nnodes = 7
            result_name = prefix + 'cpenta_stress_strain' + postfix
        elif op2.element_type == 93:
            etype = 'CHEXANL'
            nnodes = 9
            result_name = prefix + 'chexa_stress_strain' + postfix
        elif op2.element_type == 256:
            etype = 'CPYRAMNL'
            nnodes = 6
            result_name = prefix + 'chexa_stress_strain' + postfix

        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())

        numwide_real = 4 + (25 - 4) * nnodes  # real???
        numwide_random = 2 + (18 - 2) * nnodes  # imag???
        #op2.log.debug("format_code=%s numwide=%s numwide_real=%s numwide_random=%s" % (
            #op2.format_code, op2.num_wide, numwide_real, numwide_random))

        #numwide_real = 0
        #numwide_imag = 2 + 16 * nnodes
        #ntotal = 8 + 64 * nnodes

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)

        slot = op2.get_result(result_name)

        if op2.format_code == 1 and op2.num_wide == numwide_real:
            #if op2.read_mode == 1:
                #return ndata, None, None

            ntotal = numwide_real * 4 * self.factor
            #if op2.is_stress:
                #op2.create_transient_object(self.nonlinearPlateStress, NonlinearSolid)
            #else:
                #op2.create_transient_object(self.nonlinearPlateStrain, NonlinearSolid)
            #self.handle_results_buffer(self.OES_CQUAD4NL_90, resultName, name)
            raise RuntimeError('OES_CSOLIDNL_90')
        elif op2.format_code == 1 and op2.num_wide == numwide_random:  # random
            # 82 : CTETRA_NL (etype=85)
            # 146 : CHEXA_NL (etype=93)
            #raise RuntimeError(op2.code_information())
        #elif op2.format_code in [2, 3] and op2.num_wide == numwide_imag:  # imag

            ntotal = numwide_random * 4 * self.factor

            nelements = ndata // ntotal
            self.ntotal += nelements * nnodes
            #print(op2.read_mode, RealNonlinearSolidArray)
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, RealNonlinearSolidArray)
            if auto_return:
                op2._data_factor = nnodes
                return nelements * ntotal, None, None

            nelements = ndata // ntotal
            obj = op2.obj
            n = oes_csolidnl_real(op2, data, obj, etype, nnodes, nelements, ntotal)

        else:  # pragma: no cover
            #msg = op2.code_information()
            msg = "format_code=%s numwide=%s numwide_real=%s numwide_random=%s\n" % (
                op2.format_code, op2.num_wide, numwide_real, numwide_random)
            #return op2._not_implemented_or_skip(data, ndata, msg)
            raise RuntimeError(msg + op2.code_information())
        return n, nelements, ntotal


    def _oes_weld_118(self, data, ndata, dt, is_magnitude_phase,
                      result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 118 : WELDP

        """
        op2 = self.op2
        #if isinstance(op2.nonlinear_factor, float):
            #op2.sort_bits[0] = 1 # sort2
            #op2.sort_method = 2

        n = 0
        if op2.is_stress:
            result_name = prefix + 'cweld_stress' + postfix
        else:
            result_name = prefix + 'cweld_strain' + postfix

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)
        if result_type == 0 and op2.num_wide == 8:  # real
            if op2.is_stress:
                obj_vector_real = RealWeldStressArray
            else:
                obj_vector_real = RealWeldStrainArray

            ntotal = 32 * self.factor  # 8*4
            nelements = ndata // ntotal
            #print('WELDP nelements =', nelements)

            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return ndata, None, None

            if op2.is_debug_file:
                op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                op2.binary_debug.write('  #elementi = [eid_device, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,\n')
                op2.binary_debug.write('                           s1b, s2b, s3b, s4b, smaxb, sminb,        MSc]\n')
                op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                # self.itime = 0
                # self.ielement = 0
                # self.itotal = 0
                #self.ntimes = 0
                #self.nelements = 0
                n = nelements * op2.num_wide * 4

                ielement = obj.ielement
                ielement2 = ielement + nelements
                obj._times[obj.itime] = dt
                self.obj_set_element(obj, ielement, ielement2, data, nelements)

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 8)

                #[axial, maxa, mina, maxb, minb, max_shear, bearing]
                obj.data[obj.itime, ielement:ielement2, :] = floats[:, 1:].copy()
                obj.itotal = ielement2
                obj.ielement = ielement2
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    op2.log.debug('vectorize WELDP real SORT%s' % op2.sort_method)
                n = oes_weldp_msc_real_8(op2, data, obj, nelements, ntotal, dt)
        elif result_type == 1 and op2.num_wide == 15:  # complex
            if op2.is_stress:
                obj_vector_complex = ComplexWeldStressArray
            else:
                obj_vector_complex = ComplexWeldStrainArray

            ntotal = 60 * self.factor  # 15*4
            nelements = ndata // ntotal
            #print('WELDP nelements =', nelements)

            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_complex)
            if auto_return:
                return ndata, None, None

            if op2.is_debug_file:
                op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                op2.binary_debug.write('  #elementi = [eid_device, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,\n')
                op2.binary_debug.write('                           s1b, s2b, s3b, s4b, smaxb, sminb,        MSc]\n')
                op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1 and 0:
                # self.itime = 0
                # self.ielement = 0
                # self.itotal = 0
                #self.ntimes = 0
                #self.nelements = 0
                n = nelements * op2.num_wide * 4

                ielement = obj.ielement
                ielement2 = ielement + nelements
                obj._times[obj.itime] = dt
                self.obj_set_element(obj, ielement, ielement2, data, nelements)

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 8)

                #[axial, maxa, mina, maxb, minb, max_shear, bearing]
                obj.data[obj.itime, ielement:ielement2, :] = floats[:, 1:].copy()
                obj.itotal = ielement2
                obj.ielement = ielement2
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    op2.log.debug('vectorize WELDP real SORT%s' % op2.sort_method)
                n = oes_weldp_msc_complex_15(op2, data, obj, nelements, ntotal, is_magnitude_phase, dt)
            return n, None, None
        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
        assert op2.obj.element_name == op2.element_name, op2.obj
        assert n > 0
        return n, nelements, ntotal

    def _oes_fast_126(self, data, ndata, dt, is_magnitude_phase,
                      result_type, prefix, postfix):
        r"""
        reads stress/strain for element type:
         - 126 : FASTP

        C:\MSC.Software\msc_nastran_runs\cf103e.op2
         S T R E S S E S   I N   F A S T E N E R   E L E M E N T S   ( C F A S T )

         ELEMENT-ID         FORCE-X          FORCE-Y          FORCE-Z         MOMENT-X         MOMENT-Y         MOMENT-Z
         data  = (301, -4.547473508864641e-09, 1.8571810755929619e-09, -7.94031507211912e-10, -0.0, -0.0, 0.0,
                  401, -4.547473508864641e-09, -2.0263790645458357e-09, 1.1617373729677638e-09, -0.0, 0.0, 0.0)

        """
        op2 = self.op2
        #if isinstance(op2.nonlinear_factor, float):
            #op2.sort_bits[0] = 1 # sort2
            #op2.sort_method = 2

        n = 0
        if op2.is_stress:
            result_name = prefix + 'cfast_stress' + postfix
        else:
            result_name = prefix + 'cfast_strain' + postfix

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)
        if result_type == 0 and op2.num_wide == 7:  # real
            if op2.is_stress:
                obj_vector_real = RealFastStressArray
            else:
                obj_vector_real = RealFastStrainArray

            ntotal = 28 * self.factor  # 7*4
            nelements = ndata // ntotal
            assert ndata % ntotal == 0
            #print('WELDP nelements =', nelements)

            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return ndata, None, None

            if op2.is_debug_file:
                op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                op2.binary_debug.write('  #elementi = [eid_device, force_x, force_y, force_z, moment_x, moment_y, moment_z]\n')
                op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                # self.itime = 0
                # self.ielement = 0
                # self.itotal = 0
                #self.ntimes = 0
                #self.nelements = 0
                n = nelements * op2.num_wide * 4

                ielement = obj.ielement
                ielement2 = ielement + nelements
                obj._times[obj.itime] = dt
                self.obj_set_element(obj, ielement, ielement2, data, nelements)

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 7)

                #[force_x, force_y, force_z, moment_x, moment_y, moment_z]
                obj.data[obj.itime, ielement:ielement2, :] = floats[:, 1:].copy()
                obj.itotal = ielement2
                obj.ielement = ielement2
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    op2.log.debug('vectorize real-FASTP real SORT%s' % op2.sort_method)
                n = oes_fastp_msc_real_7(op2, data, obj, nelements, ntotal, dt)
        elif result_type == 1 and op2.num_wide == 13:  # complex
            if op2.is_stress:
                obj_vector_complex = ComplexFastStressArray
            else:
                obj_vector_complex = ComplexFastStrainArray
            ntotal = 52 * self.factor  # 13*4
            nelements = ndata // ntotal
            assert ndata % ntotal == 0
            #print('FASTP nelements =', nelements)

            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_complex)
            if auto_return:
                return ndata, None, None

            if op2.is_debug_file:
                op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                op2.binary_debug.write('  #elementi = [eid_device, disp_x, disp_y, disp_z, rotation_x, rotation_y, rotation_z]\n')
                op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                # self.itime = 0
                # self.ielement = 0
                # self.itotal = 0
                #self.ntimes = 0
                #self.nelements = 0
                n = nelements * op2.num_wide * 4

                ielement = obj.ielement
                ielement2 = ielement + nelements
                obj._times[obj.itime] = dt
                self.obj_set_element(obj, ielement, ielement2, data, nelements)

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 13)# [:, 1:]
                isave_real = [1, 2, 3, 4, 5, 6]
                isave_imag = [7, 8, 9, 10, 11, 12]
                real_imag = apply_mag_phase(floats, is_magnitude_phase, isave_real, isave_imag)


                #[force_x, force_y, force_z, moment_x, moment_y, moment_z]
                obj.data[obj.itime, ielement:ielement2, :] = real_imag # .copy()
                obj.itotal = ielement2
                obj.ielement = ielement2
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    op2.log.debug('vectorize complex-FASTP real SORT%s' % op2.sort_method)
                n = oes_fastp_msc_complex_13(op2, data, obj, nelements, ntotal, is_magnitude_phase, dt)
        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
        assert op2.obj.element_name == op2.element_name, op2.obj
        assert n > 0
        return n, nelements, ntotal

    def _oes_cquad4_33(self, data, ndata: int, dt, is_magnitude_phase: bool,
                       result_type: int, prefix: str, postfix: str) -> tuple[int, Any, Any]:
        """
        reads stress/strain for element type:
         - 33 : CQUAD4-centroidal
         - 64 : CQUAD4-centroidal (nastran95)
         - 228 : CQUADR-centroidal

        """
        op2 = self.op2
        n = 0
        factor = self.factor
        size = self.size
        #print('_oes_cquad4_33')
        if op2.is_stress:
            obj_vector_real = RealPlateStressArray
            obj_vector_complex = ComplexPlateStressArray
            if op2.element_type == 33:
                result_name = prefix + 'cquad4_stress' + postfix
            elif op2.element_type == 228:
                result_name = prefix + 'cquadr_stress' + postfix
                assert op2.num_wide in [17, 15], op2.code_information()

            elif op2._nastran_format == 'nasa95':
                if op2.element_type == 19:
                    result_name = prefix + 'cquad1_stress' + postfix
                elif op2.element_type == 64:
                    result_name = prefix + 'cquad4_stress' + postfix
                else:  # pragma: no cover
                    raise NotImplementedError(op2.code_information())
            else:  # pragma: no cover
                raise NotImplementedError(op2.code_information())
        else:
            obj_vector_real = RealPlateStrainArray
            obj_vector_complex = ComplexPlateStrainArray
            if op2.element_type == 33:
                result_name = prefix + 'cquad4_strain' + postfix
            elif op2.element_type == 228:
                result_name = prefix + 'cquadr_strain' + postfix
                assert op2.num_wide in [17, 15], op2.code_information()
            elif op2._nastran_format == 'nasa95':
                if op2.element_type == 19:
                    result_name = prefix + 'cquad1_strain' + postfix
                elif op2.element_type == 64:
                    result_name = prefix + 'cquad4_strain' + postfix
                else:  # pragma: no cover
                    raise NotImplementedError(op2.code_information())
            else:  # pragma: no cover
                raise NotImplementedError(op2.code_information())

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        numwide_real = 17
        sort_method = op2.sort_method
        if result_type == 0 and op2.num_wide == 17:  # real
            ntotal = 68 * factor  # 4*17
            nelements = ndata // ntotal
            nlayers = nelements * 2  # 2 layers per node
            #op2.log.info(f'CQUAD4-33: len(data)={ndata} numwide={op2.num_wide} nelements={nelements} nlayers={nlayers}')

            auto_return, is_vectorized = op2._create_oes_object4(
                nlayers, result_name, slot, obj_vector_real)
            if auto_return:
                op2._data_factor = 2  # number of "layers" for an element
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and sort_method == 1:
                nnodes_expected = 2
                n = nelements * ntotal
                ielement = obj.ielement
                ielement2 = ielement + nelements
                itotal = obj.itotal
                itotal2 = itotal + nelements * nnodes_expected
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype8)
                    ints1 = ints.reshape(nelements, numwide_real)
                    eids = ints1[:, 0] // 10
                    eids = np.vstack([eids, eids]).T.ravel()
                    assert eids.min() > 0, eids.min()
                    obj.element_node[itotal:itotal2, 0] = eids

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, numwide_real)[:, 1:]

                #fd, sx, sy, txy, angle, major, minor, max_shear
                floats1 = floats.reshape(nelements * nnodes_expected, 8)
                obj.data[obj.itime, itotal:itotal2, :] = floats1.copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    op2.log.debug(f'vectorize centroidal quad: {op2.element_name}-{op2.element_type} real '
                                   f'SORT{sort_method}')
                n = oes_quad4_33_real_17(op2, data, obj, ntotal, nelements, dt)
            if op2.is_sort1:
                assert obj.element_node[:, 0].min() > 0, obj.element_node[:, 0].shape

        elif result_type == 1 and op2.num_wide == 15:  # imag
            #op2.to_nx(f' because CQUAD4-33 (numwide=15) was found')
            #nnodes = 0  # centroid + 4 corner points
            ntotal = op2.num_wide * size
            #op2.log.info(f'CQUAD4-33: len(data)={ndata} numwide={op2.num_wide} nelements={nelements} nlayers={nlayers}')

            nelements = ndata // ntotal
            nlayers = nelements * 2
            auto_return, is_vectorized = op2._create_oes_object4(
                nlayers, result_name, slot, obj_vector_complex)
            if auto_return:
                op2._data_factor = 2
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and sort_method == 1:
                n = nelements * ntotal
                nnodes_all = 1
                itotal = obj.itotal
                itotal2 = itotal + 2 * nelements * nnodes_all
                ielement = obj.ielement
                ielement2 = ielement + nelements

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 15 * nnodes_all)
                floats1 = floats[:, 1:].reshape(nelements * nnodes_all * 2, 7).copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, 15 * nnodes_all).copy()
                    eids = ints[:, 0] // 10
                    ints[:, 0] = 0
                    ints1 = ints.reshape(nelements * nnodes_all, 15)
                    nids = ints[:, 0]
                    assert eids.min() > 0, eids.min()
                    eids2 = np.vstack([eids, eids]).T.ravel()
                    nids2 = np.vstack([nids, nids]).T.ravel()
                    obj.element_node[itotal:itotal2, 0] = eids2
                    obj.element_node[itotal:itotal2, 1] = nids2

                #[fd, sxr, sxi, syr, syi, txyr, txyi]
                isave1 = [1, 3, 5]
                isave2 = [2, 4, 6]
                real_imag = apply_mag_phase(floats1, is_magnitude_phase, isave1, isave2)

                obj.fiber_distance[itotal:itotal2] = floats1[:, 0]
                obj.data[obj.itime, itotal:itotal2, :] = real_imag
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    op2.log.debug(f'vectorize CQUAD4-33 imag SORT{sort_method}')

                n = oes_cquad4_33_complex_15(
                    op2, data, obj,
                    nelements, ntotal,
                    is_magnitude_phase)

        elif result_type in [1, 2] and op2.num_wide == 9: # random msc
            # _oes_cquad4 is the same as _oes_ctria3
            element_id = op2.nonlinear_factor
            if op2.is_stress:
                obj_vector_random = RandomPlateStressArray
            else:
                obj_vector_random = RandomPlateStrainArray
            op2.data_code['nonlinear_factor'] = element_id

            if op2._results.is_not_saved(result_name):
                op2._data_factor = 2
                return ndata, None, None
            op2._results._found_result(result_name)
            slot = op2.get_result(result_name)

            ntotal = 36 * self.factor  # 4*9
            nelements = ndata // ntotal
            nlayers = nelements * 2
            nnodes_expected = 1
            #if op2.table_name_str.startswith('OSTRRMS'):
                #print(f'{op2.table_name_str} {result_name}: {nelements} {ntotal}; sort_method={sort_method}')

            auto_return, is_vectorized = op2._create_oes_object4(
                nlayers, result_name, slot, obj_vector_random)
            if auto_return:
                op2._data_factor = 2
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized:
                n = nelements * ntotal
                itotal = obj.itotal
                itotal2 = itotal + nelements * 2

                if sort_method == 1:
                    ielement = obj.ielement
                    ielement2 = ielement + nelements

                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = frombuffer(data, dtype=op2.idtype)
                        ints1 = ints.reshape(nelements, 9)
                        eids = ints1[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        #print(eids)
                        eids2 = np.vstack([eids, eids]).T.ravel()
                        #print(eids2)

                        # TODO: what about layer 1/2?
                        #print(op2.code_information())
                        #print(f'eids.shape={eids.shape} obj.element_node.shape={obj.element_node.shape}')
                        obj.element_node[itotal:itotal2, 0] = eids2

                    floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 9)[:, 1:]
                    #fd1, sx1, sy1, txy1, fd2, fx2, fy2, txy2
                    floats2 = floats.reshape(nelements * nnodes_expected, 8)

                    #[eid_device, fd1, sx1, sy1, txy1,
                    #             fd2, sx2, sy2, txy2,]
                    nf2 = floats2.shape[0]
                    floats3 = floats2.reshape(nf2*2, 4)

                    #print(obj)
                    #print(op2.code_information())
                    obj.fiber_distance[itotal:itotal2] = floats3[:, 0].copy()
                    obj.data[obj.itime, itotal:itotal2, :] = floats3[:, 1:].copy()
                    obj.itotal = itotal2
                    obj.ielement = ielement2

                elif sort_method == 2 and op2._analysis_code_fmt == b'f':
                    ielement = obj.itime
                    ie_upper = 2 * ielement
                    ie_lower = 2 * ielement + 1

                    obj.element_node[ie_upper, 0] = dt
                    obj.element_node[ie_lower, 0] = dt
                    #obj._times[obj.itime] = dt

                    floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 9)# [:, 1:]
                    # itime is actually ielement
                    # we grab the element id from the ints for all times
                    if op2._analysis_code_fmt == b'i' and obj.itime == 0:
                        #print('analysis_code ', op2.analysis_code, op2._analysis_code_fmt)
                        ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, 9)
                        eids = ints[:, 0] // 10
                        #nids = np.zeros(len(eids), dtype='int32')
                        #print(eids)
                        #eids = np.vstack([eids, nids]).T.ravel()
                        #print(eids.shape)
                        #print(eids)
                        #print(obj.element)
                        #assert eids.min() > 0, eids.min()
                        #obj.element[itotal:itotal2, 0] = eids
                        obj._times[itotal:itotal2] = eids
                        aaa
                    elif op2._analysis_code_fmt == b'f' and obj.itime == 0:
                        #print(floats[:, 0])
                        #print(floats[:, 0].shape, obj._times.shape)
                        obj._times[itotal:itotal2] = floats[:, 0]

                    floats1 = floats[:, 1:]
                    #print(floats1)
                    #print(floats1.shape)
                    #fd, sx, sy, txy,
                    floats2 = floats1.reshape(nelements * nnodes_expected, 8)
                    nf2 = floats2.shape[0]
                    # reshape it into 2 layers
                    floats3 = floats2.reshape(nf2*2, 4)
                    # we only need to grab the first two fiber/curvature values
                    # as they're duplicated many times for the same element
                    obj.fiber_distance[2*obj.itime:2*obj.itime+2] = floats3[:2, 0].copy()
                    # we apply the data across 2 rows because we have 2 layers
                    obj.data[:, ie_upper, :] = floats3[::2, 1:].copy()
                    obj.data[:, ie_lower, :] = floats3[1::2, 1:].copy()
                else:
                    raise NotImplementedError(op2.code_information())
                obj.itotal = itotal2
                #obj.ielement = ielement2
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    op2.log.debug(f'vectorize CQUAD4-33 random numwide=9 SORT{sort_method}')
                n = oes_cquad4_33_random_9(op2, data, obj, nelements, ntotal)

        elif result_type in [1, 2] and op2.num_wide == 11: # random
            #2 FD1 RS Z1 = Fibre Distance
            #3 SX1 RS Normal in x at Z1
            #4 SY1 RS Normal in y at Z1
            #5 TXY1 RS Shear in xy at Z1
            #6 RMSVM1 RS RMS von Mises at Z1

            #7 FD2 RS Z2 = Fibre Distance
            #8 SX2 RS Normal in x at Z2
            #9 SY2 RS Normal in y at Z2
            #10 TXY2 RS Shear in xy at Z2
            #11 RMSVM2 RS RMS von Mises at Z2

            element_id = op2.nonlinear_factor
            if op2.is_stress:
                obj_vector_random = RandomPlateVMStressArray
            else:
                obj_vector_random = RandomPlateVMStrainArray
            op2.data_code['nonlinear_factor'] = element_id

            if op2._results.is_not_saved(result_name):
                return ndata, None, None
            op2._results._found_result(result_name)
            slot = op2.get_result(result_name)

            ntotal = 44 * factor # 4*11
            nelements = ndata // ntotal
            nlayers = nelements * 2
            nnodes_expected = 1

            auto_return, is_vectorized = op2._create_oes_object4(
                nlayers, result_name, slot, obj_vector_random)
            if auto_return:
                op2._data_factor = 2
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and 0:  # pragma: no cover
                n = nelements * 4 * op2.num_wide
                ielement = obj.ielement
                ielement2 = ielement + nelements
                itotal = obj.itotal
                itotal2 = itotal + nelements * nnodes_expected
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype)
                    ints1 = ints.reshape(nelements, 9)
                    eids = ints1[:, 0] // 10
                    print(eids)
                    eids = np.vstack([eids, eids]).T.ravel()
                    print(eids.shape)
                    print(eids)
                    print(obj.element)
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2, 0] = eids

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 11)[:, 1:]
                print(floats.shape)
                #fd, sx, sy, txy,
                floats1 = floats.reshape(nelements * nnodes_expected, 10)
                obj.data[obj.itime, itotal:itotal2, :] = floats1.copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                n = oes_cquad4_33_random_vm_11(op2, data, obj, nelements, ntotal)

        elif result_type == 1 and op2.num_wide == 17 and op2.table_name in [b'OESVM1', b'OESVM2', b'OSTRVM1', b'OSTRVM2']: # freq
            # Table of element stresses for frequency response analysis that includes
            # von Mises stress output in SORT1 format.
            element_id = op2.nonlinear_factor
            if op2.is_stress:  # TODO: add new complex type
                obj_vector_complex = ComplexPlateVMStressArray
            else:
                obj_vector_complex = ComplexPlateVMStrainArray
            op2.data_code['nonlinear_factor'] = element_id

            if op2._results.is_not_saved(result_name):
                return ndata, None, None
            op2._results._found_result(result_name)
            slot = op2.get_result(result_name)

            ntotal = 68 * self.factor  # 4*17
            nelements = ndata // ntotal
            nlayers = nelements * 2
            nnodes_expected = 1
            #op2.log.info(f'CQUAD4-33: len(data)={ndata} numwide={op2.num_wide} nelements={nelements} nlayers={nlayers}')

            auto_return, is_vectorized = op2._create_oes_object4(
                nlayers, result_name, slot, obj_vector_complex)
            if auto_return:
                op2._data_factor = 2
                return nelements * ntotal, None, None

            #self.show_data(data)
            #   ELEMENT      FIBER                                     - STRESSES IN ELEMENT  COORDINATE SYSTEM -
            #      ID.       DISTANCE              NORMAL-X                       NORMAL-Y                      SHEAR-XY               VON MISES
            #0     101  -5.000000E-01  -8.152692E-01 /  0.0           -1.321875E+00 /  0.0           -3.158517E+00 /  0.0            5.591334E+00
            #            5.000000E-01   1.728573E+00 /  0.0           -7.103837E+00 /  0.0            2.856040E+00 /  0.0            9.497519E+00
                #floats  = (1011,
                           #-0.5, -0.8152692317962646, 0.0, -1.321874737739563, 0.0, -3.1585168838500977, 0.0, 5.591334342956543,
                           #0.5,   1.7285730838775635, 0.0, -7.103837490081787, 0.0,  2.8560397624969482, 0.0, 9.497518539428711)
            obj = op2.obj
            if is_vectorized and op2.use_vector and op2.sort_method == 1 and 0:  # pragma: no cover
                raise NotImplementedError(op2.table_name_str)
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    op2.log.debug('vectorize CQUAD4-33 complex '
                                   f'{op2.table_name_str} SORT{op2.sort_method}')
                n = oes_cquad4_33_complex_vm_17(op2, data, obj, nelements, ntotal,
                                                is_magnitude_phase)

        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
        assert op2.obj.element_name == op2.element_name, op2.obj
        assert n > 0
        return n, nelements, ntotal

    def _oes_ctria3(self, data, ndata: int, dt, is_magnitude_phase: bool,
                    result_type: int, prefix: str, postfix: str) -> tuple[int, Any, Any]:
        """
        reads stress/strain for element type:
         - 74 : CTRIA3-centroidal
         - 83 : CTRIA3-centroidal  (NASA 95)
         - 227: TRIAR-centroidal

        """
        op2 = self.op2
        #print('_oes_ctria3')
        n = 0

        etype_map = {
            #element_type : (element_base, element_name)
            74 : ('ctria3', 'CTRIA3'),
            83 : ('ctria3', 'CTRIA3'),  # NASA-95
            227 : ('ctriar', 'CTRIAR'),
        }
        if op2.is_stress:
            stress_strain = 'stress'
            obj_vector_real = RealPlateStressArray
            obj_vector_complex = ComplexPlateStressArray
        else:
            stress_strain = 'strain'
            obj_vector_real = RealPlateStrainArray
            obj_vector_complex = ComplexPlateStrainArray

        #if prefix == '' and postfix == '':
            #prefix = stress_strain + '.'

        element_base, element_name = etype_map[op2.element_type]
        # stress.ctria3_stress
        result_name = prefix + f'{element_base}_{stress_strain}' + postfix

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)

        slot = op2.get_result(result_name)
        #print(op2.element_name, result_name, op2.format_code, op2.num_wide)
        #table_names = [
        #    b'OES1', b'OES1X', b'OES1X1', b'OSTR1X',
        #    b'OES2', b'OSTR2',
        #    b'OESVM1', b'OESVM2', b'OSTRVM1', b'OSTRVM2',
        #    b'OESPSD1', b'OESRMS1',
        #    b'OESPSD2',  b'OESATO2',  b'OESCRM2',  b'OESNO1',  b'OESXRMS1', b'OESXNO1',
        #    b'OSTRPSD2', b'OSTRATO2', b'OSTRCRM2', b'OSTRNO1', b'OSTRRMS1',
        #]
        #assert op2.table_name in table_names, op2.table_name
        sort_method = op2.sort_method
        element_name_type = f'{op2.element_name}-{op2.element_type}'

        if op2.format_code in [1, 3] and op2.num_wide == 17:  # real
            ntotal = 68 * self.factor # 4*17
            nelements = ndata // ntotal
            nlayers = nelements * 2  # 2 layers per node
            #if self.code[0] == 100:
                #print(f'\nnelements={nelements} nlayers={nlayers} {self.code in slot}')
            auto_return, is_vectorized = op2._create_oes_object4(
                nlayers, result_name, slot, obj_vector_real)
            if auto_return:
                op2._data_factor = 2
                return nelements * ntotal, None, None

            if op2.is_debug_file:
                op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                op2.binary_debug.write('  #elementi = [eid_device, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
                op2.binary_debug.write('  #                        fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,]\n')
                op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            obj = op2.obj
            if op2.use_vector and is_vectorized and sort_method == 1:
                nfields = 17 * nelements
                nbytes = nfields * 4
                itotal = obj.itotal
                iend = obj.itotal + nlayers

                itime = obj.itime
                if itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype8).reshape(nelements, 17)
                    eids = ints[:, 0] // 10
                    #ilayers = ints[:, 1]
                    #ints2 = ints[:, 1:].reshape(nlayers, 8)
                    assert eids.min() > 0, eids
                    obj._times[obj.itime] = dt
                    obj.element_node[itotal:iend:2, 0] = eids
                    obj.element_node[itotal+1:iend+1:2, 0] = eids
                    #obj.element_node[itotal:iend, 1] = 0

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 17)
                floats1 = floats[:, 1:].reshape(nlayers, 8).copy()
                obj.data[obj.itime, itotal:iend, :] = floats1
                obj._times[obj.itime] = dt
                obj.itotal += nlayers
                n = nbytes
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    op2.log.debug(f'vectorize centroidal tri: {element_name_type} real SORT{sort_method}')
                n = oes_ctria3_real_17(op2, data, obj,
                                       ntotal, nelements, dt)
            if op2.is_sort1:
                assert obj.element_node[:, 0].min() > 0, obj.element_node[:, 0]

        elif op2.format_code in [2, 3] and op2.num_wide == 15:  # imag
            ntotal = 60 * self.factor # 4*15
            nelements = ndata // ntotal
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_complex)
            if auto_return:
                op2._data_factor = 2
                #assert ntotal == op2.num_wide * 4
                return nelements * ntotal, None, None
            obj = op2.obj

            if op2.use_vector and is_vectorized and sort_method == 1:
                n = nelements * ntotal
                itotal = obj.itotal
                itotal2 = itotal + nelements * 2
                ielement = obj.ielement
                ielement2 = ielement + nelements

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 15)
                floats1 = floats[:, 1:].reshape(nelements * 2, 7).copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype8).reshape(nelements, 15).copy()
                    eids = ints[:, 0] // 10
                    ints[:, 0] = 0
                    unused_ints1 = ints.reshape(nelements, 15)
                    nids = ints[:, 0]

                    assert eids.min() > 0, eids.min()
                    eids2 = np.vstack([eids, eids]).T.ravel()
                    nids2 = np.vstack([nids, nids]).T.ravel()
                    obj.element_node[itotal:itotal2, 0] = eids2
                    obj.element_node[itotal:itotal2, 1] = nids2

                #[fd, sxr, sxi, syr, syi, txyr, txyi]
                isave1 = [1, 3, 5]
                isave2 = [2, 4, 6]
                real_imag = apply_mag_phase(floats1, is_magnitude_phase, isave1, isave2)

                obj.fiber_distance[itotal:itotal2] = floats1[:, 0]
                obj.data[obj.itime, itotal:itotal2, :] = real_imag
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    op2.log.debug(f'vectorize CTRIA3 imag SORT{sort_method}')
                struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'14f', self.size))
                cen = 0 # CEN/3
                for unused_i in range(nelements):
                    edata = data[n:n + ntotal]
                    out = struct1.unpack(edata)
                    (eid_device,
                     fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
                     fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i,) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, op2.nonlinear_factor, op2.sort_method)

                    if op2.is_debug_file:
                        op2.binary_debug.write('  OESC %s - eid=%i; C=[%s]\n' % (
                            element_name_type, eid,
                            ', '.join(['%r' % di for di in out])))

                    if is_magnitude_phase:
                        sx1 = polar_to_real_imag(sx1r, sx1i)
                        sy1 = polar_to_real_imag(sy1r, sy1i)
                        sx2 = polar_to_real_imag(sx2r, sx2i)
                        sy2 = polar_to_real_imag(sy2r, sy2i)
                        txy1 = polar_to_real_imag(txy1r, txy1i)
                        txy2 = polar_to_real_imag(txy2r, txy2i)
                    else:
                        sx1 = complex(sx1r, sx1i)
                        sy1 = complex(sy1r, sy1i)
                        sx2 = complex(sx2r, sx2i)
                        sy2 = complex(sy2r, sy2i)
                        txy1 = complex(txy1r, txy1i)
                        txy2 = complex(txy2r, txy2i)
                    obj.add_sort1(dt, eid, cen,
                                  fd1, sx1, sy1, txy1,
                                  fd2, sx2, sy2, txy2)
                    n += ntotal
        #elif op2.format_code == 1 and op2.num_wide == 9: # random?
            #msg = op2.code_information()
            #return op2._not_implemented_or_skip(data, ndata, msg), None, None
        elif op2.format_code in [2, 3] and op2.num_wide == 17 and op2.table_name in [b'OESVM1', b'OESVM2', b'OSTRVM1', b'OSTRVM2']:
            # freq:
            # # random; CTRIA3
            assert op2.table_name in [b'OESVM1', b'OESVM2', b'OSTRVM1', b'OSTRVM2'], op2.code_information()

            element_id = op2.nonlinear_factor
            if op2.is_stress:  # TODO: add new complex type
                obj_vector_complex = ComplexPlateVMStressArray
            else:
                obj_vector_complex = ComplexPlateVMStrainArray
            op2.data_code['nonlinear_factor'] = element_id

            if op2._results.is_not_saved(result_name):
                return ndata, None, None
            op2._results._found_result(result_name)
            slot = op2.get_result(result_name)

            ntotal = 68 * self.factor  # 17*4
            nelements = ndata // ntotal
            nlayers = nelements * 2
            nnodes_expected = 1

            auto_return, is_vectorized = op2._create_oes_object4(
                nlayers, result_name, slot, obj_vector_complex)
            if auto_return:
                op2._data_factor = 2
                return nelements * ntotal, None, None

            obj = op2.obj
            # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\shlthk14.op2
            nelements = ndata // ntotal
            assert ndata % ntotal == 0

            #                   C O M P L E X   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )
            #                                                          (REAL/IMAGINARY)
            #
            #  ELEMENT       FIBER                                     - STRESSES IN ELEMENT  COORDINATE SYSTEM -
            #    ID.        DISTANCE              NORMAL-X                       NORMAL-Y                      SHEAR-XY               VON MISES
            #0       1  -4.359080E+00  -1.391918E+00 /  2.474756E-03  -1.423926E+00 /  2.530494E-03   2.655153E-02 / -5.158625E-05   1.408948E+00
            #            4.359080E+00   1.391918E+00 / -2.474756E-03   1.423926E+00 / -2.530494E-03  -2.655153E-02 /  5.158625E-05   1.408948E+00
            n = oes_ctria3_complex_vm_17(op2, data, obj, nelements, ntotal, dt,
                                         is_magnitude_phase)
            assert n is not None, n

        elif op2.format_code in [1, 2, 3] and op2.num_wide == 11: # random; CTRIA3
            #2 FD1 RS Z1 = Fibre Distance
            #3 SX1 RS Normal in x at Z1
            #4 SY1 RS Normal in y at Z1
            #5 TXY1 RS Shear in xy at Z1
            #6 RMSVM1 RS RMS von Mises at Z1

            #7 FD2 RS Z2 = Fibre Distance
            #8 SX2 RS Normal in x at Z2
            #9 SY2 RS Normal in y at Z2
            #10 TXY2 RS Shear in xy at Z2
            #11 RMSVM2 RS RMS von Mises at Z2

            element_id = op2.nonlinear_factor
            if op2.is_stress:
                obj_vector_random = RandomPlateVMStressArray
            else:
                obj_vector_random = RandomPlateVMStrainArray
            op2.data_code['nonlinear_factor'] = element_id

            if op2._results.is_not_saved(result_name):
                return ndata, None, None
            op2._results._found_result(result_name)
            slot = op2.get_result(result_name)

            ntotal = 44 * self.factor  # 4*11
            nelements = ndata // ntotal
            nlayers = nelements * 2
            nnodes_expected = 1

            auto_return, is_vectorized = op2._create_oes_object4(
                nlayers, result_name, slot, obj_vector_random)
            if auto_return:
                op2._data_factor = 2
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and 0:  # pragma: no cover
                n = nelements * 4 * op2.num_wide
                ielement = obj.ielement
                ielement2 = ielement + nelements
                itotal = obj.itotal
                itotal2 = itotal + nelements * nnodes_expected
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype)
                    ints1 = ints.reshape(nelements, 9)
                    eids = ints1[:, 0] // 10
                    print(eids)
                    eids = np.vstack([eids, eids]).T.ravel()
                    print(eids.shape)
                    print(eids)
                    print(obj.element)
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2, 0] = eids

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 11)[:, 1:]
                print(floats.shape)
                #fd, sx, sy, txy,
                floats1 = floats.reshape(nelements * nnodes_expected, 10)
                obj.data[obj.itime, itotal:itotal2, :] = floats1.copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if is_vectorized and op2.use_vector and obj.itime == 0:  # pragma: no cover
                    op2.log.debug(f'vectorize {element_name_type} random numwide=11 SORT{sort_method}')
                n = oes_ctria3_random_vm_11(op2, data, obj, nelements, ntotal)

        elif op2.format_code in [1, 2] and op2.num_wide == 9: # random MSC stress/strain; CTRIA3
            # _oes_cquad4 is the same as _oes_ctria3
            element_id = op2.nonlinear_factor
            if op2.is_stress:
                obj_vector_random = RandomPlateStressArray
                #result_name = prefix + 'ctria3_stress' + postfix
            else:
                obj_vector_random = RandomPlateStrainArray
                #result_name = prefix + 'ctria3_strain' + postfix

            op2.data_code['nonlinear_factor'] = element_id

            if op2._results.is_not_saved(result_name):
                return ndata, None, None
            op2._results._found_result(result_name)
            slot = op2.get_result(result_name)

            ntotal = 36 * self.factor # 4*9
            nelements = ndata // ntotal
            nlayers = nelements * 2
            nnodes_expected = 1

            auto_return, is_vectorized = op2._create_oes_object4(
                nlayers, result_name, slot, obj_vector_random)
            if auto_return:
                op2._data_factor = 2
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and 0:  # pragma: no cover
                n = nelements * 4 * op2.num_wide
                ielement = obj.ielement
                ielement2 = ielement + nelements
                itotal = obj.itotal
                itotal2 = itotal + nelements * nnodes_expected
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype)
                    ints1 = ints.reshape(nelements, 9)
                    eids = ints1[:, 0] // 10
                    #print(eids)
                    eids = np.vstack([eids, eids]).T.ravel()
                    #print(eids.shape)
                    #print(eids)
                    #print(obj.element)
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2, 0] = eids

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 9)[:, 1:]
                #print(floats.shape)
                #fd, sx, sy, txy,
                floats1 = floats.reshape(nelements * nnodes_expected, 8)
                obj.data[obj.itime, itotal:itotal2, :] = floats1.copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    op2.log.debug(f'vectorize {element_name_type} random2 SORT{sort_method}')
                n = oes_ctria3_random_9(op2, data, obj, nelements, ntotal)

        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
        assert n is not None, op2.code_information()
        return n, nelements, ntotal

    def _oes_cquad4_144(self, data: bytes, ndata: int, dt, is_magnitude_phase: bool,
                        result_type: int, prefix: str, postfix: str) -> tuple[int, Any, Any]:
        """
        reads stress/strain for element type:
         - 64 : CQUAD8
         - 70 : CTRIAR
         - 75 : CTRIA6
         - 82 : CQUADR
         - 144 : CQUAD4-bilinear

        """
        op2 = self.op2
        n = 0
        size = self.size

        etype_map = {
            #element_type : (element_base, nnodes_expected, element_name)
            64 : ('cquad8', 4, 'CQUAD8'),
            70 : ('ctriar', 3, 'CTRIAR'),
            75 : ('ctria6', 3, 'CTRIA6'),
            82 : ('cquadr', 4, 'CQUADR'),
            144 : ('cquad4', 4, 'CQUAD4-bilinear'),
        }
        if op2.is_stress:
            stress_strain = 'stress'
            obj_vector_real = RealPlateStressArray
            obj_vector_complex = ComplexPlateStressArray
            obj_vector_random = RandomPlateStressArray
        else:
            stress_strain = 'strain'
            obj_vector_real = RealPlateStrainArray
            obj_vector_complex = ComplexPlateStrainArray
            obj_vector_random = RandomPlateStrainArray

        # centroid not included in nnodes
        element_base, nnodes, element_name = etype_map[op2.element_type]
        #if prefix == '' and postfix == '':
            #prefix = stress_strain + '.'

        # stress.cquad4_stress
        result_name = prefix + f'{element_base}_{stress_strain}' + postfix

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        log = op2.log

        nnodes_all = nnodes + 1 # adding the centroid

        slot = op2.get_result(result_name)
        numwide_real = 2 + 17 * nnodes_all
        numwide_imag = 2 + 15 * nnodes_all
        numwide_random = 2 + 9 * nnodes_all

        #numwide_imag2 = 2 + 16 * nnodes_all
        #print('%s real=%s imag=%s imag2=%s random=%s' % (
            #op2.element_name, numwide_real, numwide_imag, numwide_imag2, numwide_random
        #))
        #etype = op2.element_name
        #grid_center = 'CEN/%i' % nnodes

        # OESVM1/2 (complex)
        #   87 - CQUAD8

        # OSTRNO1
        #   47 - CQUAD4

        # CQUAD8 real=87 imag=77 imag2=82 random=47
        # CQUAD4 ???=
        sort_method = op2.sort_method
        element_name_type = f'{op2.element_name}-{op2.element_type}'
        #print(op2.code_information())

        #if result_type == 0 and op2.num_wide == numwide_real and op2._nastran_format == 'optistruct' and op2.element_type == 75:  # real
            #could be a quad or tri
            #x = 1

        if result_type == 0 and op2.num_wide == numwide_real:  # real
            ntotal = 4 * (2 + 17 * nnodes_all) * self.factor
            nelements = ndata // ntotal
            assert ndata % ntotal == 0
            nlayers = 2 * nelements * nnodes_all  # 2 layers per node

            auto_return, is_vectorized = op2._create_oes_object4(
                nlayers, result_name, slot, obj_vector_real)
            if auto_return:
                op2._data_factor = 2 * nnodes_all  # number of "layers" for an element
                return nelements * ntotal, None, None

            obj = op2.obj
            #print('dt=%s, itime=%s' % (obj.itime, dt))
            if op2.use_vector and is_vectorized and sort_method == 1:
                # self.itime = 0
                # self.ielement = 0
                # self.itotal = 0
                #self.ntimes = 0
                #self.nelements = 0
                n = nelements * ntotal

                istart = obj.itotal
                iend = istart + nlayers
                obj._times[obj.itime] = dt

                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype8).reshape(nelements, numwide_real)
                    ints1 = ints[:, 2:].reshape(nlayers//2, 17)[:, 0].reshape(nelements, nnodes_all).copy()
                    ints1[:, 0] = 0.
                    nids = ints1.ravel()

                    eids = ints[:, 0] // 10
                    eids2 = array([eids] * (nnodes_all * 2), dtype=op2.idtype8).T.ravel()
                    nids2 = vstack([nids, nids]).T.ravel()
                    obj.element_node[istart:iend, 0] = eids2
                    obj.element_node[istart:iend, 1] = nids2
                    #assert obj.element_node[:iend, 0].min() > 0, eids2
                    if obj.nonlinear_factor is not None:
                        float_mask = np.arange(nelements * numwide_real, dtype=np.int32).reshape(nelements, numwide_real)
                        float_mask1 = float_mask[:, 2:].reshape(nlayers // 2, 17)[:, 1:].reshape(nlayers, 8)
                        obj.float_mask = float_mask1

                if obj.nonlinear_factor is not None:
                    results = frombuffer(data, dtype=op2.fdtype8)[obj.float_mask].copy()
                else:
                    floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, numwide_real)
                    floats1 = floats[:, 2:].reshape(nlayers // 2, 17)
                    results = floats1[:, 1:].reshape(nlayers, 8).copy()

                #[fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm]
                obj.data[obj.itime, istart:iend, :] = results
                assert obj.element_node[:, 0].min() > 0, obj.element_node[:, 0]
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    log.debug(f'vectorize nodal shell: {element_name_type}... real SORT{sort_method}')

                n = oes_cquad4_144_real(op2, data, ndata, obj,
                                        nelements, nnodes, dt)
            if op2.is_sort1:
                assert obj.element_node[:, 0].min() > 0, obj.element_node[:, 0]

        elif result_type == 1 and op2.num_wide == numwide_imag:  # complex
            ntotal = numwide_imag * 4 * self.factor
            #assert op2.num_wide * 4 == ntotal, 'numwide*4=%s ntotal=%s' % (op2.num_wide*4, ntotal)
            nelements = ndata // ntotal
            nlayers = nelements * 2 * nnodes_all
            #print(element_name_type)
            #print('ndata', ndata)
            #print('ntotal', ntotal)
            #print('nelements', nelements)
            #print('nlayers', nlayers)

            auto_return, is_vectorized = op2._create_oes_object4(
                nlayers, result_name, slot, obj_vector_complex)
            if auto_return:
                op2._data_factor = 2 * nnodes_all
                #if op2.num_wide == 77:
                    #print('ntotal =', ndata, ntotal, op2.num_wide)
                    #print('nelements * ntotal =', nelements * ntotal)
                return nelements * ntotal, None, None
            obj = op2.obj

            if op2.use_vector and is_vectorized and sort_method == 1:
                n = nelements * ntotal
                itotal = obj.itotal
                itotal2 = itotal + nelements * (nnodes_all * 2)
                ielement = obj.ielement
                ielement2 = ielement + nelements

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, numwide_imag)
                floats1 = floats[:, 2:].reshape(nelements * nnodes_all, 15)
                floats2 = floats1[:, 1:].reshape(nelements * nnodes_all * 2, 7).copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype8).reshape(nelements, numwide_imag).copy()
                    ints[:, 2] = 0  # set center node to 0
                    ints1 = ints[:, 2:].reshape(nelements * nnodes_all, 15)
                    eids = ints[:, 0] // 10
                    nids = ints1[:, 0]
                    eids2 = np.vstack([eids] * (nnodes_all * 2)).T.ravel()
                    nids2 = np.vstack([nids, nids]).T.ravel()
                    assert eids.min() > 0, eids.min()
                    obj.element_node[itotal:itotal2, 0] = eids2
                    obj.element_node[itotal:itotal2, 1] = nids2

                #[fd, sxr, sxi, syr, syi, txyr, txyi]
                isave1 = [1, 3, 5]
                isave2 = [2, 4, 6]
                real_imag = apply_mag_phase(floats2, is_magnitude_phase, isave1, isave2)

                obj.fiber_distance[itotal:itotal2] = floats2[:, 0]
                obj.data[obj.itime, itotal:itotal2, :] = real_imag
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    log.debug(f'vectorize CQUAD4-144/{element_name_type}... imag SORT{sort_method}')
                # nnodes_cquad4 = 5
                # nelements = 3
                # nlayers = nelements * nodes_cquad4 * 2 = 3*5*2 = 30
                # ntotal = nlayers
                n = oes_cquad4_144_complex_77(op2, data, obj,
                                              nelements, nnodes,
                                              dt, is_magnitude_phase)
        #elif op2.format_code == 1 and op2.num_wide == numwide_random: # random
            #msg = op2.code_information()
            #msg += '  numwide=%s numwide_real=%s numwide_imag=%s numwide_random=%s' % (
                #op2.num_wide, numwide_real, numwide_imag, numwide_random)
            #return op2._not_implemented_or_skip(data, ndata, msg), None, None
        elif result_type == 2 and op2.num_wide == numwide_random: # random
            # 47 - CQUAD8-64
            # 38 - CTRIAR-70
            ntotal = op2.num_wide * 4 * self.factor
            nelements = ndata // ntotal
            assert ndata % ntotal == 0
            nlayers = 2 * nelements * nnodes_all  # 2 layers per node
            #op2.log.info(f'random quad-144 ntotal={ntotal} ndata={ndata} ntotal={ntotal} numwide={op2.num_wide} -> nelements={nelements}')

            #if op2.read_mode == 1:
                #msg = ''
                #return op2._not_implemented_or_skip(data, ndata, msg), None, None

            auto_return, is_vectorized = op2._create_oes_object4(
                nlayers, result_name, slot, obj_vector_random)
            if auto_return:
                op2._data_factor = 2 * nnodes_all  # number of "layers" for an element
                return nelements * ntotal, None, None

            obj = op2.obj
            #print('dt=%s, itime=%s' % (obj.itime, dt))
            if op2.use_vector and is_vectorized and 0:
                # self.itime = 0
                # self.ielement = 0
                # self.itotal = 0
                #self.ntimes = 0
                #self.nelements = 0
                n = nelements * ntotal

                istart = obj.itotal
                iend = istart + nlayers
                obj._times[obj.itime] = dt

                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, numwide_real)
                    ints1 = ints[:, 2:].reshape(nlayers//2, 17)[:, 0].reshape(nelements, nnodes_all).copy()
                    ints1[:, 0] = 0.
                    nids = ints1.ravel()

                    eids = ints[:, 0] // 10
                    eids2 = array([eids] * (nnodes_all * 2), dtype='int32').T.ravel()
                    nids2 = vstack([nids, nids]).T.ravel()
                    obj.element_node[istart:iend, 0] = eids2
                    obj.element_node[istart:iend, 1] = nids2
                    #assert obj.element_node[:iend, 0].min() > 0, eids2
                    if obj.nonlinear_factor is not None:
                        float_mask = np.arange(nelements * numwide_real, dtype=np.int32).reshape(nelements, numwide_real)
                        float_mask1 = float_mask[:, 2:].reshape(nlayers // 2, 17)[:, 1:].reshape(nlayers, 8)
                        obj.float_mask = float_mask1

                if obj.nonlinear_factor is not None:
                    results = frombuffer(data, dtype=op2.fdtype)[obj.float_mask].copy()
                else:
                    floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, numwide_real)
                    floats1 = floats[:, 2:].reshape(nlayers // 2, 17)
                    results = floats1[:, 1:].reshape(nlayers, 8).copy()

                #[fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm]
                obj.data[obj.itime, istart:iend, :] = results
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    log.debug(f'vectorize CQUAD4-144/{element_name_type}... random SORT{sort_method}')
                #numwide_random = 2 + 9 * nnodes_all
                n = oes_cquad4_144_random(op2, data, obj, nelements, nnodes, ndata)

            #if op2.read_mode == 1:
                #msg = ''
                #return op2._not_implemented_or_skip(data, ndata, msg), None, None
            ##self.show_data(data[:44])
            #ntotal = 44
            #struct1 = Struct(op2._endian + b'i4si 8f')
            #for i in ra
            #for i in range(20):
                #edata = data[n:n+ntotal]
                #out = struct1.unpack(edata)
                #self.show_data(edata)
                #print(out)
                #n += ntotal

            ## 47 - CQUAD8-64
            ##msg = op2.code_information()
            #msg = '%s-CQUAD4-numwide=%s format_code=%s;\n numwide_real=%s numwide_imag=%s numwide_random=%s' % (
                #op2.table_name_str, op2.num_wide, op2.format_code,
                #numwide_real, numwide_imag, numwide_random)
            #return op2._not_implemented_or_skip(data, ndata, msg), None, None
        elif op2.table_name in [b'OESVM1', b'OESVM2', b'OSTRVM1', b'OSTRVM2']:
            # 82  CQUADR -> 87
            # CQUAD8 sort_method=2 ntotal=348 nelements=3
            # CTRIA6 sort_method=2 ntotal=348 nelements=2

            msg = op2.code_information()
            if result_type == 1: # complex
                # ndata = 3828
                # ???
                #
                # ndata=1044
                assert op2.num_wide in [70, 87], op2.code_information()
                ntotal = op2.num_wide * self.size # 87*4
                nelements = ndata // ntotal
                nlayers = nelements * 2
                assert ndata % ntotal == 0

                if op2.is_stress:
                    obj_vector_complex = ComplexPlateVMStressArray
                else:
                    obj_vector_complex = ComplexPlateVMStrainArray

                auto_return, is_vectorized = op2._create_oes_object4(
                    nlayers, result_name, slot, obj_vector_complex)

                if auto_return:
                    op2._data_factor = 2 * nnodes_all  # number of "layers" for an element
                    #op2._data_factor = 2
                    return nelements * ntotal, None, None
                #if op2.read_mode == 1:
                    #return op2._not_implemented_or_skip(data, ndata, msg), None, None

                #self.show_data(data)
                #print(ndata, ntotal)
                obj = op2.obj
                n = oes_cquad4_complex_vm_87(op2, data, obj, nelements, nnodes_all,
                                             is_magnitude_phase)

            #if result_type == 1 and op2.num_wide in [70, 87]:
                # 70 - CTRIA6-75
                # 87 - CQUAD4-144
                #pass
            else:
                #msg = (f'skipping {op2.table_name_str}-{op2.element_name}: numwide={op2.num_wide} '
                       #f'result_type={op2.result_type} (complex);\n numwide_real={numwide_real} '
                       #f'numwide_imag={numwide_imag} numwide_random={numwide_random}')
                #print(msg)
                raise NotImplementedError(msg)
            #return op2._not_implemented_or_skip(data, ndata, msg), None, None

        #elif op2.format_code in [2, 3] and op2.num_wide == 70:
            ## 87 - CQUAD4-144
            ##msg = op2.code_information()
            #msg = '%s-CTRIA6-numwide=%s numwide_real=%s numwide_imag=%s numwide_random=%s' % (
                #op2.table_name_str, op2.num_wide, numwide_real, numwide_imag, numwide_random)
            #return op2._not_implemented_or_skip(data, ndata, msg), None, None

        elif result_type == 2 and op2.table_name in [b'OESXRMS1', b'OESXNO1']:
            # CTRIA6-75  numwide=46  C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1081x.op2
            # CQUAD8-64  numwide=64  C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1081x.op2

            # corners + centroid
            if op2.element_type == 75:  # CTRIA6  (numwide=46)
                # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1081x.op2
                nnodes = 4
            elif op2.element_type in [64, 144]:
                # CQUAD8     (numwide=57) C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1081x.op2
                # CQUAD4-144 (numwide=57) C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\plate_111o.op2
                nnodes = 5
            else:
                raise RuntimeError(op2.code_information())

            assert 2 + 11 * nnodes == op2.num_wide, op2.code_information()
            if op2.is_stress:
                obj_vector_random = RandomPlateVMStressArray
            else:
                obj_vector_random = RandomPlateVMStrainArray

            ntotal = op2.num_wide * size
            nelements = ndata // ntotal
            nlayers = nelements * nnodes_all * 2
            assert ndata % ntotal == 0
            #print(f'selement_name={op2.element_name} sort_method={op2.sort_method} ntotal={ntotal} num_wide={op2.num_wide} nelements={nelements} ndata={ndata} nlayers={nlayers}')
            #print('nnodes_all =', nnodes_all)
            # 57 = 2 + 11 * 5
            #ints    = (64011, 'CEN/',
            #           5, -1127428915, 1168048724, 1174717287, 1170166698, 1179261635, 0.025, 1167390785, 1175199013, 1169871798, 1179273875,
            #           80000, -1127428915, 1144450062, 1161227278, 1155248064, 1165738662, 0.025, 1149212546, 1165989762, 1154256422, 1167150968,
            #           80002, -1127428915, 1182003448, 1198780664, 1170939895, 1197489549, 0.025, 1182601858, 1199379074, 1169253917, 1197953979,
            #           80202, -1127428915, 1171034189, 1187612615, 1180541053, 1191448609, 0.025, 1169771716, 1186411799, 1181780048, 1191501390,
            #           80200, -1127428915, 1178323054, 1132616020, 1157628875, 1178722187, 0.025, 1176075105, 1152142782, 1150736003, 1177059240,
            #           64021, 'CEN/',
            #           5,
            # -1127428915, 1180051065, 1181683390, 1180538582, 1189027141,
            #  0.025, 1176941913, 1182077843, 1179498554, 1188072152,
            #           80200, -1127428915, 1179821552, 1170149790, 1168886008, 1181456870, 0.025, 1178168847, 1171555762, 1167658184, 1179691422,
            #           80202, -1127428915, 1180244047, 1197120002, 1184242951, 1198270294, 0.025, 1180604356, 1197449345, 1183498771, 1198177490,
            #           80402, -1127428915, 1192586977, 1184669070, 1184219971, 1194803525, 0.025, 1193976276, 1184750091, 1182770117, 1194901195,
            #           80400, -1127428915, 1199592178, 1184219800, 1179170385, 1198821390, 0.025, 1198216710, 1184079548, 1179988640, 1197735836,
            #           64111, 'CEN/',
            #           8,     -1127428915, 1176450101, 1187579168, 1179109877, 1190548891, 0.025, 1176466945, 1187779876, 1179033188, 1190633996,
            #           90000, -1127428915, 1163199692, 1169764570, 1159841955, 1171614247, 0.025, 1163848412, 1172543886, 1158998260, 1173065416,
            #           90002, -1127428915, 1183839433, 1201603141, 1152572706, 1200652206, 0.025, 1184231034, 1201902140, 1154931831, 1200916849,
            #           90202, -1127428915, 1156439069, 1187515149, 1193045713, 1200691352, 0.025, 1130954657, 1188796022, 1193252160, 1200929208,
            #           90200, -1127428915, 1188552094, 1155874499, 1172127542, 1189476821, 0.025, 1187333567, 1120409103, 1169694799, 1188368246,
            #           64121, 'CEN/',
            #           8,     -1127428915, 1188544334, 1178584960, 1188275834, 1196282134, 0.025, 1186800673, 1178575754, 1187794850, 1195571664,
            #           90200, -1127428915, 1189840387, 1170594118, 1173995469, 1190183343, 0.025, 1188877086, 1165971753, 1173556714, 1189607807,
            #           90202, -1127428915, 1140886249, 1189485297, 1193242194, 1200973710, 0.025, 1161620989, 1191030157, 1193534822, 1201299068,
            #           90402, -1127428915, 1176972259, 1174260793, 1195749260, 1202470249, 0.025, 1185547735, 1170629093, 1194910619, 1201965268,
            #           90400, -1127428915, 1202771906, 1185377334, 1175183109, 1201858966, 0.025, 1202359853, 1184055805, 1173072162, 1201506363)
            #floats  = (8.969851599989587e-41, 'CEN/',
            #           5,                      -0.025, 5088.291015625, 8496.8505859375, 6122.4580078125, 12934.6904296875,
            #                                    0.025, 4767.03173828125, 8967.2861328125, 5978.4638671875, 12946.6435546875,
            #           1.1210387714598537e-40, -0.025, 731.6883544921875, 2926.75341796875, 1757.4921875, 4028.16552734375,
            #                                    0.025, 1022.3673095703125, 4089.46923828125, 1636.442138671875, 4649.93359375,
            #           1.1210667974291402e-40, -0.025, 15612.2421875, 62448.96875, 6499.99560546875, 57405.55078125,
            #                                    0.025, 16196.626953125, 64786.5078125, 5676.76416015625, 59219.73046875,
            #           1.1238693943577898e-40, -0.025, 6546.03759765625, 25795.888671875, 14184.1220703125, 33808.12890625,
            #                                    0.025, 5929.595703125, 23450.544921875, 15394.078125, 34014.3046875,
            #           1.1238413683885033e-40, -0.025, 12018.107421875, 260.6978759765625, 2048.237060546875, 12407.8857421875,
            #                                    0.025, 9822.8447265625, 1378.429443359375, 1206.7034912109375, 10783.9140625,
            #           8.971252898453911e-41, 'CEN/',
            #           5,                      -0.025, 13705.6181640625, 15299.685546875, 14181.708984375, 28558.634765625,
            #                                    0.025, 10669.3369140625, 15684.8935546875, 13166.056640625, 26693.421875,
            #           1.1238413683885033e-40, -0.025, 13481.484375, 6114.2021484375, 5497.12109375, 15078.474609375,
            #                                    0.025, 11867.5146484375, 6800.7119140625, 4897.59765625, 13354.404296875,
            #           1.1238693943577898e-40, -0.025, 13894.0771484375, 55962.0078125, 19214.513671875, 60455.3359375,
            #                                    0.025, 14245.94140625, 57248.50390625, 17761.037109375, 60092.8203125,
            #           1.1266719912864394e-40, -0.025, 38254.87890625, 20046.77734375, 19169.630859375, 46913.26953125,
            #                                    0.025, 43681.828125, 20205.021484375, 16360.9423828125, 47294.79296875,
            #           1.126643965317153e-40,  -0.025, 65701.890625, 19169.296875, 12845.5791015625, 62608.0546875,
            #                                    0.025, 60246.0234375, 18895.3671875, 13644.65625, 58367.609375,
            #           8.983864584632835e-41, 'CEN/',
            #           8,                      -0.025, 10189.0517578125, 25730.5625, 12786.4892578125, 31530.802734375,
            #                                    0.025, 10205.5009765625, 26122.5703125, 12711.59765625, 31697.0234375,
            #           1.2611686178923354e-40, -0.025, 3408.2998046875, 5926.1064453125, 2588.539794921875, 6829.26904296875,
            #                                    0.025, 3566.6787109375, 7283.1943359375, 2382.5595703125, 7537.84765625,
            #           1.2611966438616219e-40, -0.025, 18426.392578125, 81412.5390625, 1430.910400390625, 73983.359375,
            #                                    0.025, 19191.23828125, 83748.46875, 1718.8895263671875, 76050.8828125,
            #           1.2639992407902715e-40, -0.025, 1902.8785400390625, 25605.525390625, 40046.81640625, 74289.1875,
            #                                    0.025, 232.99855041503906, 28107.23046875, 40853.25, 76147.4375,
            #           1.263971214820985e-40,  -0.025, 27630.80859375, 1833.9613037109375, 7079.9013671875, 29436.916015625,
            #                                    0.025, 25250.873046875, 100.04308319091797, 5892.03857421875, 27271.73046875,
            #           8.98526588309716e-41, 'CEN/',
            #           8,                     -0.025, 27615.65234375, 12273.875, 27091.23828125, 52689.0859375,
            #                                   0.025, 24210.064453125, 12264.884765625, 26151.81640625, 49913.8125,
            #           1.263971214820985e-40, -0.025, 30147.005859375, 6331.1591796875, 7991.97509765625, 30816.841796875,
            #                                   0.025, 28265.55859375, 4085.072509765625, 7777.7392578125, 29692.748046875,
            #           1.2639992407902715e-40,-0.025, 514.1704711914062, 29453.470703125, 40814.3203125, 76495.109375,
            #                                   0.025, 3022.874267578125, 32470.775390625, 41957.3984375, 79036.96875,
            #           1.2668018377189211e-40,-0.025, 10698.9716796875, 8121.52783203125, 50607.546875, 88186.8203125,
            #                                   0.025, 21762.919921875, 6348.23681640625, 47331.60546875, 84241.65625,
            #           1.2667738117496346e-40,-0.025, 90543.515625, 21430.10546875, 8951.7548828125, 83411.171875,
            #                                   0.025, 87324.3515625, 18848.994140625, 7541.1416015625, 80656.4609375)
            #if op2.read_mode == 2:
                #self.show_data(data)
                #ddd
            auto_return, is_vectorized = op2._create_oes_object4(
                nlayers, result_name, slot, obj_vector_random)
            if auto_return:
                op2._data_factor = 2 * nnodes_all  # number of "layers" for an element
                return nelements * ntotal, None, None

            obj = op2.obj
            n = oes_cquad4_random_vm_57(op2, data, op2.obj, nelements, ntotal, nnodes,
                                        dt)
        elif result_type == 1 and op2.table_name in [b'OESPSD1'] and op2.element_type == 82 and op2.num_wide == 47: # QUADR
            msg = op2.code_information()
            return op2._not_implemented_or_skip(data, ndata, msg), None, None
        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
        return n, nelements, ntotal

    def _oes_shells_nonlinear(self, data, ndata, dt, is_magnitude_phase,
                              result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 88 : CTRIA3NL
         - 90 : CQUAD4NL

        """
        op2 = self.op2
        n = 0
        if op2.is_stress:
            if op2.element_type == 88:
                result_name = prefix + 'ctria3_stress' # + postfix  nonlinear_
            elif op2.element_type == 90:
                result_name = prefix + 'cquad4_stress'  + postfix # nonlinear_
            else:
                raise RuntimeError(op2.element_type)
        else:
            if op2.element_type == 88:
                result_name = prefix + 'ctria3_strain' + postfix # nonlinear_
            elif op2.element_type == 90:
                result_name = prefix + 'cquad4_strain' + postfix # nonlinear_
            else:
                raise RuntimeError(op2.element_type)

        slot = op2.get_result(result_name)
        op2._results._found_result(result_name)
        #print(op2.code_information())

        log = op2.log
        if op2.format_code == 1 and op2.num_wide == 13 and op2.element_type in [88, 90]:  # real
            # single layered hyperelastic (???) ctria3, cquad4
            ntotal = 52 * self.factor  # 4*13
            nelements = ndata // ntotal

            obj_vector_real = RealNonlinearPlateArray
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * op2.num_wide * 4

                ielement = obj.ielement
                ielement2 = ielement + nelements
                obj._times[obj.itime] = dt

                self.obj_set_element(obj, ielement, ielement2, data, nelements)

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 13).copy()

                #[fiber_distance, oxx, oyy, ozz, txy, exx, eyy, ezz, exy, es, eps, ecs]
                floats[:, 1] = 0
                obj.data[obj.itime, ielement:ielement2, :] = floats[:, 1:]
                obj.ielement = ielement2
                obj.itotal = ielement2
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    log.debug('vectorize CTRIA3/CQUAD4_NL real SORT%s' % op2.sort_method)
                n = oes_cshellnl_real_13(op2, data, obj, etype, nnodes, nelements, ntotal)

        elif op2.format_code == 1 and op2.num_wide == 25 and op2.element_type in [88, 90]:
            # TODO: vectorize
            #     ELEMENT      FIBER                        STRESSES/ TOTAL STRAINS                     EQUIVALENT    EFF. STRAIN     EFF. CREEP
            #        ID      DISTANCE           X              Y             Z               XY           STRESS    PLASTIC/NLELAST     STRAIN
            # 0       721  -7.500000E+00   5.262707E+02   2.589492E+02   0.000000E+00  -2.014457E-14   4.557830E+02   5.240113E-02   0.0
            #                              4.775555E-02  -2.775558E-17  -4.625990E-02  -7.197441E-18
            #               7.500000E+00   5.262707E+02   2.589492E+02   0.000000E+00   1.308169E-14   4.557830E+02   5.240113E-02   0.0
            #                              4.775555E-02  -1.387779E-17  -4.625990E-02   4.673947E-18
            # 0       722  -7.500000E+00   5.262707E+02   2.589492E+02   0.000000E+00   2.402297E-13   4.557830E+02   5.240113E-02   0.0
            #                              4.775555E-02  -2.081668E-17  -4.625990E-02   8.583152E-17
            #               7.500000E+00   5.262707E+02   2.589492E+02   0.000000E+00   2.665485E-14   4.557830E+02   5.240113E-02   0.0
            #                              4.775555E-02  -2.081668E-17  -4.625990E-02   9.523495E-18
            #
            ntotal = 100 * self.factor  # 4*25
            nelements = ndata // ntotal
            obj_vector_real = RealNonlinearPlateArray

            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                op2._data_factor = 2
                return nelements * ntotal, None, None

            #return nelements * op2.num_wide * 4
            obj = op2.obj
            is_vectorized = False
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * op2.num_wide * 4

                ielement = obj.ielement
                ielement2 = ielement + nelements
                itotal = obj.itotal
                itotal2 = itotal + nelements * 2
                obj._times[obj.itime] = dt
                #print('ielement=%s:%s' % (ielement, ielement2))
                #print('itotal=%s:%s' % (itotal, itotal2))

                if obj.itime == 0:
                    try:
                        ints = fromstring(data, dtype=op2.idtype).reshape(nelements, 25)
                    except ValueError:
                        unused_values = fromstring(data, dtype=op2.idtype)

                    eids = ints[:, 0] // 10
                    #eids2 = vstack([eids, eids]).T.ravel()
                    #print(eids.tolist())
                    obj.element[ielement:ielement2] = eids  # 150
                     #print(obj.element_node[:10, :])

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 25)[:, 1:]

                #[fiber_distance, oxx, oyy, ozz, txy, exx, eyy, ezz, exy, es, eps, ecs]
                #floats[:, 1] = 0
                obj.data[obj.itime, itotal:itotal2, :] = floats.reshape(nelements * 2, 12).copy()
                #obj.data[obj.itime, ielement:ielement2, :] = floats[:, 1:]
                obj.ielement = ielement2
                obj.itotal = itotal2
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    log.debug('vectorize CTRIA3/CQUAD4_NL imag SORT%s' % op2.sort_method)

                etype = op2.element_type
                n = oes_cshellnl_real_25(op2, data, obj, etype, nelements, ntotal)

        elif op2.format_code == 1 and op2.num_wide == 0: # random
            msg = op2.code_information()
            return op2._not_implemented_or_skip(data, ndata, msg)
        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
        return n, nelements, ntotal

    def _oes_shells_composite(self, data, ndata, dt, is_magnitude_phase,
                              result_type: int, prefix: str, postfix: str) -> tuple[int, Any, Any]:
        """
        reads stress/strain for element type:
         - 95 : CQUAD4
         - 96 : CQUAD8
         - 97 : CTRIA3
         - 98 : CTRIA6 (composite)
         - 232 : QUADRLC (CQUADR-composite)
         - 233 : TRIARLC (CTRIAR-composite)

        """
        op2 = self.op2  # type: OP2
        table_name_bytes = op2.table_name
        assert isinstance(table_name_bytes, bytes), table_name_bytes
        n = 0

        composite_element_name_map = {
            95: 'cquad4',
            96: 'cquad8',
            97: 'ctria3',
            98: 'ctria6',
            232: 'cquadr',
            233: 'ctriar',
        }
        try:
            element_name = composite_element_name_map[op2.element_type]
        except KeyError:  # pragma: no cover
            raise KeyError(op2.code_information())

        if op2.is_stress:
            stress_strain = 'stress'
            obj_vector_real = RealCompositePlateStressArray
            obj_vector_strength = RealCompositePlateStressStrengthRatioArray
            #obj_vector_complex = ComplexCompositePlateStressArray
            layered_cls = ComplexLayeredCompositeStressArray
        else:
            stress_strain = 'strain'
            obj_vector_real = RealCompositePlateStrainArray
            obj_vector_strength = None # RealCompositePlateStrainStrengthRatioArray
            #obj_vector_complex = ComplexCompositePlateStrainArray
            #obj_vector_random = RandomCompositePlateStrainArray
            layered_cls = ComplexLayeredCompositeStrainArray

        result_name = prefix + f'{element_name}_composite_{stress_strain}' + postfix
        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        etype = op2.element_name
        sort_method = op2.sort_method

        num_wide = op2.num_wide
        #print('table_name =', table_name_bytes)
        if result_type == 0 and num_wide == 11:  # real
            #                    S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A R )
            #   ELEMENT      PLY STRESSES IN FIBER AND MATRIX DIRECTIONS   INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      MAX
            #     ID          ID   NORMAL-1     NORMAL-2     SHEAR-12    SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        SHEAR
            #      7070        1  7.50000E-01  3.00000E+00  9.86167E-08  -6.58903E-08  3.00000E+00   90.00  3.00000E+00  7.50000E-01  1.12500E+00
            #      7070        2 -7.50000E-01 -3.00000E+00 -9.86167E-08   0.0          0.0           -0.00 -7.50000E-01 -3.00000E+00  1.12500E+00
            ntotal = 44 * self.factor
            nelements = ndata // ntotal
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.is_debug_file:
                op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                op2.binary_debug.write('  element1 = [eid_device, layer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)]\n')
                op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if op2.use_vector and is_vectorized and sort_method == 1:
                n = nelements * op2.num_wide * 4

                istart = obj.itotal
                iend = istart + nelements
                obj._times[obj.itime] = dt

                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype8).reshape(nelements, 11).copy()
                    eids = ints[:, 0] // 10
                    nids = ints[:, 1]
                    obj.element_layer[istart:iend, 0] = eids
                    obj.element_layer[istart:iend, 1] = nids

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 11)
                #[o1, o2, t12, t1z, t2z, angle, major, minor, ovm]
                obj.data[obj.itime, istart:iend, :] = floats[:, 2:].copy()
            else:
                if is_vectorized and op2.use_vector:  # pragma: no cover
                    op2.log.debug(f'vectorize COMP_SHELL real SORT{sort_method}')

                n = oes_comp_shell_real_11(op2, data, ndata, obj,
                                           ntotal, nelements, etype, dt)

        #elif result_type == 1 and op2.num_wide == 9:  # TODO: imag? - not done...
            # TODO: vectorize
            #raise NotImplementedError('imaginary composite stress?')
            #msg = op2.code_information()
            #nelements = ndata // ntotal
            #obj_vector_complex = None
            #auto_return, is_vectorized = op2._create_oes_object4(
                #nelements, result_name, slot, obj_vector_complex)
            #if auto_return:
                #assert ntotal == op2.num_wide * 4
                #return nelements * ntotal, None, None

            ## TODO: this is an OEF result???
            ##    furthermore the actual table is calle dout as
            ##    'i8si4f4s', not 'i8si3fi4s'
            #ntotal = 36
            #nelements = ndata // ntotal
            #s = self.struct_i
            #s2 = Struct(op2._endian + b'8si3fi4s')
            #s3 = Struct(op2._endian + b'8si4f4s')
            #for i in range(nelements):
                ##out = s.unpack(data[n:n + ntotal])
                #eid_device, = s.unpack(data[n:n+4])
                ##t, = s.unpack(data[n:n+4])

                #if eid_device > 0:
                    #out = s2.unpack(data[n+4:n+ntotal])
                #else:
                    #unused_out1 = s2.unpack(data[n+4:n+ntotal])
                    #out = s3.unpack(data[n+4:n+ntotal])
                #(theory, lamid, fp, fm, fb, fmax, fflag) = out

                #if op2.is_debug_file:
                    #op2.binary_debug.write('%s-%s - (%s) + %s\n' % (op2.element_name, op2.element_type, eid_device, str(out)))
                #obj.add_new_eid_sort1(dt, eid, theory, lamid, fp, fm, fb, fmax, fflag)
                #n += ntotal
            #raise NotImplementedError('this is a really weird case...')
        elif result_type == 1 and num_wide == 11 and table_name_bytes in [b'OESCP', b'OESTRCP']:  # complex
            # OESCP - STRAINS IN LAYERED COMPOSITE ELEMENTS (QUAD4)
            ntotal = 44 * self.factor
            nelements = ndata // ntotal

            op2.log.warning(f'skipping complex {op2.table_name_str}-PCOMP')
            complex_obj = ComplexLayeredCompositeStressArray if op2.is_stress else ComplexLayeredCompositeStrainArray
            return nelements * ntotal, None, None


            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, complex_obj)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            n = oes_shell_composite_complex_11(op2, data, obj,
                                               ntotal, nelements, sort_method,
                                               dt, is_magnitude_phase)
            return nelements * ntotal, None, None

        elif table_name_bytes in {b'OESRT', b'OESRTN'}:
            n, nelements, ntotal = self._oes_shells_composite_oesrt(
                result_name, slot,
                result_type, sort_method, obj_vector_strength,
                data, ndata, dt)

        elif result_type == 1 and num_wide == 13 and table_name_bytes in [b'OESVM1C', b'OSTRVM1C']: # complex
            #op2.log.warning(f'skipping complex {op2.table_name_str}-PCOMP (numwide=13)')
            ntotal = 52 * self.factor
            nelements = ndata // ntotal
            #return nelements * ntotal, None, None
            op2.table_name = table_name_bytes
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, layered_cls)
            if auto_return:
                return nelements * ntotal, None, None

            if is_vectorized and op2.use_vector:  # pragma: no cover
                op2.log.warning(f'vectorize COMP_SHELL complex SORT{sort_method} (numwide=13)')

            obj = op2.obj
            n = oes_shell_composite_complex_13(op2, data, obj,
                                               ntotal, nelements, sort_method,
                                               dt, is_magnitude_phase)
            #return nelements * ntotal, None, None

        elif num_wide == 7 and (result_type == 2 or table_name_bytes == b'OSTRMS1C'): # random (no VM)
            # TCODE,7 =0 Real
            # 2 PLY I Lamina Number
            # 3 EX1 RS Normal-1
            # 4 EY1 RS Normal-2
            # 5 ET1 RS Shear-12
            # 6 EL1 RS Shear-1Z
            # 7 EL2 RS Shear-2Z
            # 8 A1 RS Shear angle
            # 9 EMJRP1 RS Major Principal
            # 10 EMNRP1 RS Minor Principal
            # 11 ETMAX1 RS von Mises or Maximum shear
            #
            # TCODE,7 =1 Real/imaginary
            # 2 PLY I Lamina Number
            # 3 EX1 RS Normal-1
            # 4 EY1 RS Normal-2
            # 5 ET1 RS Shear-12
            # 6 EL1 RS Shear-1Z
            # 7 EL2 RS Shear-2Z
            # 8 EX1I RS Normal-1
            # 9 EY1I RS Normal-2
            # 10 ET1I RS Shear-12
            # 11 EL1I RS Shear-1Z
            # 12 EL2I RS Shear-2Z
            #
            # TCODE,7 =2 Random Response
            # 2 PLY I Lamina Number
            # 3 EX1 RS Normal-1
            # 4 EY1 RS Normal-2
            # 5 ET1 RS Shear-12
            # 6 EL1 RS Shear-1Z
            # 7 EL2 RS Shear-2Z
            ntotal = 28 * self.factor
            nelements = ndata // ntotal

            if op2.is_stress:
                obj_vector_random = RandomCompositePlateStressArray
            else:
                obj_vector_random = RandomCompositePlateStrainArray

            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_random)
            if auto_return:
                assert ntotal == op2.num_wide * 4
                return nelements * ntotal, None, None

            obj = op2.obj
            n = oes_composite_shells_nx_random_7(op2, data, obj, nelements, ntotal)

        elif result_type == 2 and num_wide == 8:  # random VM
            # analysis_code = 5   Frequency
            # table_code    = 805 OESXRM1C-OESXRMS1 - element RMS stresses for random analysis that includes von Mises stress output.
            # format_code   = 1   Real
            # result_type   = 2   Random
            # sort_method   = 1
            # sort_code     = 4
            #     sort_bits   = (0, 0, 1)
            #     data_format = 0   Real
            #     sort_type   = 0   Sort1
            #     is_random   = 1   Random Responses
            # random_code   = 0
            # element_type  = 95  QUAD4-nonlinear
            # num_wide      = 8
            # freq          = 0.0
            # NX Nastran
            if op2.is_stress:
                obj_vector_random = RandomCompositePlateStressVMArray
            else:
                obj_vector_random = RandomCompositePlateStrainVMArray
            msg = op2.code_information()
            msg = (f'etype={op2.element_name} ({op2.element_type}) '
                   f'{op2.table_name_str}-COMP-random-numwide={op2.num_wide} '
                   f'numwide_real=11 numwide_imag=9 result_type={result_type}')
            return op2._not_implemented_or_skip(data, ndata, msg), None, None

        elif result_type == 1 and num_wide == 12 and op2.is_msc:
            # analysis_code = 5   Frequency
            # table_code    = 5   OES1C-OES - Element Stress
            # format_code   = 2   Real/Imaginary
            # result_type   = 1   Complex
            # sort_method   = 1
            # sort_code     = 1
            #     sort_bits   = (1, 0, 1)
            #     data_format = 1   Real/Imaginary
            #     sort_type   = 0   Sort1
            #     is_random   = 1   Random Responses
            # random_code   = 0
            # element_type  = 95  QUAD4LC-composite
            # s_code        = 0   Coordinate Element - Stress Max Shear (Octahedral)
            # thermal       = 0   isHeatTransfer = False
            # thermal_bits  = [0, 0, 0, 0, 0]
            # num_wide      = 12
            # freq          = 990.0
            # MSC Nastran
            msg = (f'etype={op2.element_name} ({op2.element_type}) '
                   f'{op2.table_name_str}-COMP-complex-numwide={num_wide} '
                   #f'numwide_real=11 numwide_imag=9 result_type={result_type}'
                   )
            #op2.log.warning(f'skipping complex {op2.table_name_str}-PCOMP-12')
            return op2._not_implemented_or_skip(data, ndata, msg), None, None

        elif result_type == 1 and num_wide == 11:
            # analysis_code = 9   Complex eigenvalues
            # table_code    = 5   OESCP-OES - Element Stress
            # format_code   = 2   Real/Imaginary
            # result_type   = 1   Complex
            # sort_method   = 1
            # sort_code     = 0
            #     sort_bits   = (1, 0, 0)
            #     data_format = 1   Real/Imaginary
            #     sort_type   = 0   Sort1
            #     is_random   = 0   Sorted Responses
            # random_code   = 0
            # element_type  = 95  QUAD4-nonlinear
            # num_wide      = 11
            # mode          = 0
            # eigr          = 0.0
            # eigi          = 0.0
            # NX Nastran
            msg = op2.code_information()
            msg = (f'etype={op2.element_name} ({op2.element_type}) '
                   f'{op2.table_name_str}-COMP-random-numwide={num_wide} '
                   f'numwide_real=11 numwide_imag=9 result_type={result_type}')
            return op2._not_implemented_or_skip(data, ndata, msg), None, None
        else:
            raise RuntimeError(op2.code_information())
            #msg = op2.code_information()
            #msg = (f'etype={op2.element_name} ({op2.element_type}) '
                   #f'{op2.table_name_str}-COMP-random-numwide={op2.num_wide} '
                   #f'numwide_real=11 numwide_imag=9 result_type={result_type}')
            #return op2._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal


    def _oes_shells_composite_oesrt(self, result_name: str, slot: dict[Any, Any],
                                    result_type: int, sort_method: int,
                                    obj_vector_strength,
                                    data: bytes, ndata: int, dt: Any):
        """Table of composite element ply strength ratio. Output by SDRCOMP"""
        op2 = self.op2
        n = 0
        if result_type == 0 and op2.num_wide == 9: # real
            # strength_ratio.cquad4_composite_stress
            ntotal = 36 * self.factor
            nelements = ndata // ntotal

            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_strength)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.is_debug_file:
                op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                op2.binary_debug.write('  element1 = [eid_device, failure_theory, ply_id, strength_ratio_ply, failure_index_bonding, strength_ratio_bonding, flag, flag2)]\n')
                op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            element_type = op2.element_type
            if op2.use_vector and is_vectorized and sort_method == 1 and 0:
                n = nelements * op2.num_wide * 4
                asdf
            else:
                op2.log.warning(f'need to vectorize oes_shell_composite; {op2.element_name}-{op2.element_type} '
                                 f'(numwide={op2.num_wide}) {op2.table_name_str}')
                n = oesrt_comp_shell_real_9(op2, data, ndata, obj,
                                            ntotal, nelements, element_type, dt)
        elif result_type == 1 and op2.num_wide == 9: # complex
            # '          S T R E N G T H   R A T I O S   F O R   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )'
            # '   ELEMENT  FAILURE        PLY  SRP-STRENGTH RATIO FOR PLY  SRB-STRENGTH RATIO FOR BONDING  STRENGTH RATIO FOR ELEMENT     FLAG'
            # '     ID      THEORY         ID  (DIRECT STRESSES/STRAINS)     (INTER-LAMINAR STRESSES)      MIN OF SRP,SRB FOR ALL PLIES'
            # '         3   HILL            1      4.510528E-02      '
            # '                                                                     4.160837E+01                                               '
            # '                             2      1.414114E-01      '
            # '                                                                     4.807957E+01                                               '
            # '                             3      2.013647E-01      '
            # '                                                                                                   4.510528E-02             *** '
            ntotal = 36 * self.factor
            nelements = ndata // ntotal
            if data is None:
                return ndata, None, None #   eid hill ply_id ? i f i blank

                auto_return, is_vectorized = op2._create_oes_object4(
                    nelements, result_name, slot, obj_vector_strength)
                if auto_return:
                    return nelements * ntotal, None, None

            is_vectorized = False
            op2.log.warning(f'OESRT: faking result; {op2.element_name}-{op2.element_type}')
            if op2.use_vector and is_vectorized and sort_method == 1 and 0:
                n = nelements * op2.num_wide * 4
                asdf
            else:
                op2.log.warning(f'need to vectorize oes_shell_composite; {op2.element_name}-{op2.element_type} '
                                 f'(numwide={op2.num_wide}) {op2.table_name_str}')

                structi = Struct(op2._endian + b'i   8s   i      f i f i 4s')
                structf = Struct(op2._endian + b'i   8s   i      f i f f 4s')
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]  # 4*9
                    outs = structi.unpack(edata)
                    eid, hill_bytes, ply_id, f1, i1, f2, minus_1, blank_bytes = outs
                    hill = hill_bytes.decode('latin1').strip()
                    blank = blank_bytes.decode('latin1').strip()
                    if minus_1 != -1:
                        outs = structf.unpack(edata)
                        eid, hill_bytes, ply_id, f1, i1, f2, minus_1, blank_bytes = outs
                        #op2.log.error(f'minus1={minus_1}')
                    assert hill in ['HILL', ''], hill
                    assert blank in ['', '***'], blank
                    #op2.show_data(edata)
                    if eid != -1:
                        print(f'eid={eid} hill={hill!r} ply={ply_id} i1={i1} f2={f2:.4e} minus1={minus_1} blank={blank!r}')
                    else:
                        print(f'    hill={hill!r} ply={ply_id} i1={i1} f2={f2:.4e} minus1={minus_1:.4e} blank={blank!r}')
                    n += ntotal
        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
        return n, nelements, ntotal

    def _oes_ctriax6(self, data, ndata, dt, is_magnitude_phase,
                     result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 53 : CTRIAX6
        """
        op2 = self.op2
        n = 0
        if op2.is_stress:
            result_name = prefix + 'ctriax_stress' + postfix
        else:
            result_name = prefix + 'ctriax_strain' + postfix

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)

        slot = op2.get_result(result_name)
        if result_type == 0 and op2.num_wide == 33: # real
            if op2.is_stress:
                obj_vector_real = RealTriaxStressArray
            else:
                obj_vector_real = RealTriaxStrainArray

            ntotal = 132 * self.factor  # (1+8*4)*4 = 33*4 = 132
            nelements = ndata // ntotal
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                op2._data_factor = 4
                assert ntotal == op2.num_wide * 4
                return nelements * ntotal, None, None

            obj = op2.obj
            nnodes_all = 4
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * op2.num_wide * 4

                itotal = obj.itotal
                itotal2 = itotal + nelements * nnodes_all
                ielement = obj.ielement
                ielement2 = ielement + nelements

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 33)
                floats1 = floats[:, 1:].reshape(nelements * nnodes_all, 8).copy()

                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, 33).copy()
                    ints1 = ints[:, 1:].reshape(nelements * nnodes_all, 8)
                    eids = ints[:, 0] // 10
                    ints[:, 0] = 0
                    nids = ints1[:, 0]
                    eids2 = np.vstack([eids] * nnodes_all).T.ravel()
                    assert eids.min() > 0, eids.min()
                    obj.element_node[itotal:itotal2, 0] = eids2
                    obj.element_node[itotal:itotal2, 1] = nids

                #[loc, rs, azs, As, ss, maxp, tmax, octs]
                obj.data[obj.itime, itotal:itotal2, :] = floats1[:, 1:]
                obj.ielement = ielement2
                obj.itotal = itotal2
            else:
                n = oes_ctriax6_real_33(op2, data, obj,
                                        nelements, ntotal, dt)

        elif result_type == 1 and op2.num_wide == 37: # imag
            # TODO: vectorize object
            if op2.is_stress:
                #print('op2.element_type', op2.element_type)
                #print('op2.element_name', op2.element_name)
                obj_vector_complex = ComplexTriaxStressArray
            else:
                obj_vector_complex = ComplexTriaxStrainArray

            num_wide = 1 + 4 * 9
            ntotal = num_wide * 4 * self.factor
            assert num_wide == op2.num_wide, num_wide
            nelements = ndata // ntotal  # (1+8*4)*4 = 33*4 = 132
            leftover = ndata % ntotal
            assert leftover == 0, 'ntotal=%s nelements=%s leftover=%s' % (ntotal, nelements, leftover)

            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_complex)

            if data is None:
                return ndata, None, None
            auto_return = False
            is_vectorized = False
            if auto_return:
                op2._data_factor = 4
                return nelements * ntotal, None, None

            obj = op2.obj
            nnodes_all = 4
            if op2.use_vector and is_vectorized and 0:
                n = nelements * ntotal
                itotal = obj.itotal
                itotal2 = itotal + nelements * nnodes_all
                ielement = obj.ielement
                ielement2 = ielement + nelements

                numwide_imag = 37
                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, numwide_imag)
                floats1 = floats[:, 1:].reshape(nelements * nnodes_all, 9).copy()

                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, numwide_imag)
                    ints1 = ints[:, 1:].reshape(nelements * nnodes_all, 9).copy()
                    eids = ints[:, 0] // 10
                    ints[:, 0] = 0
                    nids = ints1[:, 0]
                    eids2 = np.vstack([eids] * nnodes_all).T.ravel()
                    assert eids.min() > 0, eids.min()
                    obj.element_node[itotal:itotal2, 0] = eids2
                    obj.element_node[itotal:itotal2, 1] = nids

                # [loc, rsr, rsi, azsr, azsi, Asr, Asi, ssr, ssi]
                isave1 = [1, 3, 5, 7]
                isave2 = [2, 4, 6, 9]
                real_imag = apply_mag_phase(floats1, is_magnitude_phase, isave1, isave2)

                obj.data[obj.itime, itotal:itotal2, :] = real_imag
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                n = oes_ctriax_complex_37(op2, data, obj,
                                          nelements,
                                          is_magnitude_phase)
        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
            #msg = op2.code_information()
            #raise NotImplementedError(msg)
            #return op2._not_implemented_or_skip(data, ndata, msg)
        return n, nelements, ntotal

    def _oes_cbush(self, data, ndata, dt, is_magnitude_phase,
                   result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 102 : CBUSH

        """
        op2 = self.op2
        n = 0
        if op2.is_stress:
            result_name = prefix + 'cbush_stress' + postfix
        else:
            result_name = prefix + 'cbush_strain' + postfix

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)

        slot = op2.get_result(result_name)
        table_name = op2.table_name_str
        if result_type in [0, 2] and op2.num_wide == 7:  # real, random
            if op2.is_stress:
                obj_vector_real = RealBushStressArray
            else:
                obj_vector_real = RealBushStrainArray

            assert op2.num_wide == 7, "num_wide=%s not 7" % op2.num_wide
            ntotal = 28 * self.factor # 4*7

            nelements = ndata // ntotal
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return nelements * ntotal, None, None
            obj = op2.obj

            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal

                istart = obj.ielement
                iend = istart + nelements
                obj._times[obj.itime] = dt

                self.obj_set_element(obj, istart, iend, data, nelements)

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 7)
                #[tx, ty, tz, rx, ry, rz]
                obj.data[obj.itime, istart:iend, :] = floats[:, 1:].copy()
            else:
                n = oes_cbush_real_7(op2, data, obj,
                                     nelements, ntotal, dt)
        elif result_type == 1 and op2.num_wide == 13:  # imag
            if op2.is_stress:
                obj_complex = ComplexCBushStressArray
            else:
                obj_complex = ComplexCBushStrainArray

            ntotal = 52 * self.factor  # 4*13
            nelements = ndata // ntotal
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_complex)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 13).copy()
                obj._times[obj.itime] = dt
                self.obj_set_element(obj, itotal, itotal2, data, nelements)

                isave1 = [1, 2, 3, 4, 5, 6]
                isave2 = [7, 8, 9, 10, 11, 12]
                real_imag = apply_mag_phase(floats, is_magnitude_phase, isave1, isave2)
                obj.data[obj.itime, itotal:itotal2, :] = real_imag

                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                n = oes_cbush_complex_13(op2, data, obj,
                                         nelements, ntotal,
                                         is_magnitude_phase)
        elif result_type == 1 and op2.num_wide == 7 and table_name == 'OESPSD1':
            result_type = 2
            n, nelements, ntotal = self._oes_cbush(data, ndata, dt, is_magnitude_phase,
                                                   result_type, prefix, postfix)
        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
            #msg = op2.code_information()
            #raise NotImplementedError(msg)
            #return op2._not_implemented_or_skip(data, ndata, msg)
        return n, nelements, ntotal

    def _oes_cbush1d(self, data, ndata, dt, is_magnitude_phase,
                     result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 40 : CBUSH1D
        """
        op2 = self.op2
        n = 0
        if op2.is_stress:
            result_name = prefix + 'cbush1d_stress_strain' + postfix
        else:
            result_name = prefix + 'cbush1d_stress_strain' + postfix

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)
        if result_type == 0 and op2.num_wide == 8:  # real
            if op2.is_stress:
                obj_vector_real = RealBush1DStressArray
            else:
                #op2.create_transient_object(self.cbush1d_stress_strain, Bush1DStrain)  # undefined
                raise NotImplementedError('cbush1d_stress_strain; numwide=8')

            ntotal = 32 * self.factor # 4*8
            nelements = ndata // ntotal
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal

                itotal = obj.itotal
                itotal2 = itotal + nelements
                itime = obj.itime
                obj._times[itime] = dt

                if 1: #obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, 8).copy()
                    eids = ints[:, 0] // 10
                    fail = ints[:, 7]
                    obj.element[itotal:itotal2] = eids
                    obj.is_failed[itime, itotal:itotal2, 0] = fail

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 8)
                #[xxx, fe, ue, ve, ao, ae, ep, xxx]
                obj.data[itime, itotal:itotal2, :] = floats[:, 1:7].copy()

                obj.ielement = itotal2
                obj.itotal = itotal2
            else:
                struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'6fi')
                for unused_i in range(nelements):
                    edata = data[n:n + ntotal]
                    out = struct1.unpack(edata)  # num_wide=25
                    if op2.is_debug_file:
                        op2.binary_debug.write('CBUSH1D-40 - %s\n' % (str(out)))
                    (eid_device, fe, ue, ve, ao, ae, ep, fail) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, op2.nonlinear_factor, op2.sort_method)

                    # axial_force, axial_displacement, axial_velocity, axial_stress,
                    # axial_strain, plastic_strain, is_failed
                    obj.add_sort1(dt, eid, fe, ue, ve, ao, ae, ep, fail)
                    n += ntotal
        elif result_type == 1 and op2.num_wide == 9:  # imag
            # TODO: vectorize object
            ntotal = 36 * self.factor  # 4*9
            nelements = ndata // ntotal

            if op2.is_stress:
                auto_return, is_vectorized = op2._create_oes_object4(
                    nelements, result_name, slot, ComplexCBush1DStressArray)
            else:
                raise NotImplementedError('self.cbush1d_stress_strain; complex strain')

            if auto_return:
                assert ntotal == op2.num_wide * 4
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * op2.num_wide * 4

                itotal = obj.itotal
                itotal2 = itotal + nelements
                itime = obj.itime
                obj._times[itime] = dt
                self.obj_set_element(obj, itotal, itotal2, data, nelements)

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 9).copy()
                #[fer, uer, aor, aer,
                # fei, uei, aoi, aei]
                isave1 = [1, 3, 5, 7]
                isave2 = [2, 4, 6, 8]
                real_imag = apply_mag_phase(floats, is_magnitude_phase, isave1, isave2)
                obj.data[obj.itime, itotal:itotal2, :] = real_imag

                obj.ielement = itotal2
                obj.itotal = itotal2
            else:
                struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'8f')
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]

                    out = struct1.unpack(edata)  # num_wide=25
                    (eid_device,
                     fer, uer, aor, aer,
                     fei, uei, aoi, aei) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, op2.nonlinear_factor, op2.sort_method)

                    if is_magnitude_phase:
                        fe = polar_to_real_imag(fer, fei)
                        ue = polar_to_real_imag(uer, uei)
                        ao = polar_to_real_imag(aor, aoi)
                        ae = polar_to_real_imag(aer, aei)
                    else:
                        fe = complex(fer, fei)
                        ue = complex(uer, uei)
                        ao = complex(aor, aoi)
                        ae = complex(aer, aei)
                    obj.add_new_eid(op2.element_type, dt, eid, fe, ue, ao, ae)
        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
            #msg = op2.code_information()
            #raise NotImplementedError(msg)
            #return op2._not_implemented_or_skip(data, ndata, msg)
        return n, nelements, ntotal

    def _oes_crod_nonlinear(self, data, ndata, dt, unused_is_magnitude_phase,
                            result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 87 : CTUBENL
         - 89 : RODNL
         - 92 : CONRODNL
        """
        op2 = self.op2
        n = 0
        #prefix = 'nonlinear_'
        if op2.is_stress:
            if op2.element_type == 87:
                result_name = prefix + 'ctube_stress' + postfix
                name = 'CTUBENL-87'
            elif op2.element_type == 89:
                result_name = prefix + 'crod_stress' + postfix
                name = 'RODNL-89'
            elif op2.element_type == 92:
                result_name = prefix + 'conrod_stress' + postfix
                name = 'CONRODNL-92'
            else:  # pragma: no cover
                raise RuntimeError(op2.code_information())
        else:
            if op2.element_type == 87:
                result_name = prefix + 'ctube_strain' + postfix
                name = 'CTUBENL-87'
            elif op2.element_type == 89:
                result_name = prefix + 'crod_strain' + postfix
                name = 'RODNL-89'
            elif op2.element_type == 92:
                result_name = prefix + 'conrod_strain' + postfix
                name = 'CONRODNL-92'
            else:  # pragma: no cover
                raise RuntimeError(op2.code_information())

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)
        if result_type == 0 and op2.num_wide == 7:  # real
            ntotal = 28 * self.factor #  7*4 = 28
            nelements = ndata // ntotal
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, RealNonlinearRodArray)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            #if op2.is_debug_file:
                #op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                #op2.binary_debug.write('  element1 = [eid_device, layer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)]\n')
                #op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * op2.num_wide * 4
                istart = obj.itotal
                iend = istart + nelements
                obj._times[obj.itime] = dt

                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype8).reshape(nelements, 7).copy()
                    eids = ints[:, 0] // 10
                    obj.element[istart:iend] = eids
                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 7)
                #[axial_stress, equiv_stress, total_strain,
                # eff_plastic_creep_strain, eff_creep_strain, linear_torsional_stresss]
                obj.data[obj.itime, istart:iend, :] = floats[:, 1:].copy()
            else:
                struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'6f', self.size))  # 1+6=7
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]
                    out = struct1.unpack(edata)

                    (eid_device, axial_stress, equiv_stress, total_strain,
                     eff_plastic_creep_strain, eff_creep_strain, linear_torsional_stresss) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, op2.nonlinear_factor, op2.sort_method)
                    if op2.is_debug_file:
                        op2.binary_debug.write('%s - %s\n' % (name, str(out)))
                    obj.add_sort1(dt, eid, axial_stress, equiv_stress, total_strain,
                                  eff_plastic_creep_strain, eff_creep_strain, linear_torsional_stresss)
                    n += ntotal
        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
        return n, nelements, ntotal

    def _oes_celas_nonlinear(self, data, ndata, dt, unused_is_magnitude_phase,
                             result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 224 : CELAS1
         - 226 : CELAS3

        """
        op2 = self.op2
        # 224-CELAS1
        # 225-CELAS3
        # NonlinearSpringStress
        n = 0
        numwide_real = 3
        if op2.is_stress:
            if op2.element_type == 224:
                result_name = prefix + 'celas1_stress' + postfix # nonlinear_
            elif op2.element_type == 225:
                result_name = prefix + 'celas3_stress' + postfix # nonlinear_
        else:
            raise NotImplementedError('NonlinearSpringStrain')

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)

        slot = op2.get_result(result_name)
        if result_type == 0 and op2.num_wide == numwide_real:
            assert op2.num_wide == 3, "num_wide=%s not 3" % op2.num_wide
            ntotal = 12 * self.factor  # 4*3
            nelements = ndata // ntotal

            if op2.is_stress:
                auto_return, is_vectorized = op2._create_oes_object4(
                    nelements, result_name, slot, RealNonlinearSpringStressArray)
            else:
                raise NotImplementedError('NonlinearSpringStrainArray') # undefined

            if auto_return:
                assert ntotal == op2.num_wide * 4
                return nelements * ntotal, None, None
            obj = op2.obj

            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * 4 * op2.num_wide
                unused_itotal = obj.ielement
                ielement = obj.ielement
                ielement2 = obj.ielement + nelements
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, numwide_real).copy()
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[ielement:ielement2] = eids

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, numwide_real)

                #[force, stress]
                obj.data[obj.itime, ielement:ielement2, :] = floats[:, 1:].copy()
                obj.itotal = ielement2
                obj.ielement = ielement2
            else:
                struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'2f')
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]
                    out = struct1.unpack(edata)  # num_wide=3
                    (eid_device, force, stress) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, op2.nonlinear_factor, op2.sort_method)
                    if op2.is_debug_file:
                        op2.binary_debug.write('%s-%s - %s\n' % (op2.element_name, op2.element_type, str(out)))
                    obj.add_sort1(dt, eid, force, stress)
                    n += ntotal
        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
        return n, nelements, ntotal

    def _oes_cbush_nonlinear(self, data, ndata, dt, unused_is_magnitude_phase,
                             result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 226 : CBUSHNL
        """
        op2 = self.op2
        n = 0
        if op2.is_stress:
            if op2.element_type == 226:
                result_name = prefix + 'cbush_force_stress_strain' + postfix
                name = 'CBUSHNL-226'
            else:  # pragma: no cover
                raise RuntimeError(op2.code_information())
        else:
            if op2.element_type == 226:
                result_name = prefix + 'nonlinear_cbush_strain' + postfix
                name = 'CBUSHNL-226'
            else:  # pragma: no cover
                raise RuntimeError(op2.code_information())

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)
        if result_type == 0 and op2.num_wide == 19:  # real
            ntotal = 76 * self.factor  #  19*4 = 76
            nelements = ndata // ntotal
            assert ndata % ntotal == 0
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, RealNonlinearBushArray)
            if auto_return:
                assert ntotal == op2.num_wide * 4
                return nelements * ntotal, None, None

            obj = op2.obj
            #if op2.is_debug_file:
                #op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                #op2.binary_debug.write('  element1 = [eid_device, layer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)]\n')
                #op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * op2.num_wide * 4
                istart = obj.itotal
                iend = istart + nelements
                obj._times[obj.itime] = dt

                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, 19).copy()
                    eids = ints[:, 0] // 10
                    obj.element[istart:iend] = eids
                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 19)
                #[fx, fy, fz, otx, oty, otz, etx, ety, etz,
                # mx, my, mz, orx, ory, orz, erx, ery, erz]
                obj.data[obj.itime, istart:iend, :] = floats[:, 1:].copy()
            else:
                #             N O N L I N E A R   F O R C E S  A N D  S T R E S S E S  I N   B U S H   E L E M E N T S    ( C B U S H )
                #
                #                           F O R,C E                               S T R E S S                             S T R A I N
                # ELEMENT ID.   FORCE-X      FORCE-Y      FORCE-Z       STRESS-TX    STRESS-TY    STRESS-TZ     STRAIN-TX    STRAIN-TY    STRAIN-TZ
                #               MOMENT-X     MOMENT-Y     MOMENT-Z      STRESS-RX    STRESS-RY    STRESS-RZ     STRAIN-RX    STRAIN-RY    STRAIN-RZ
                #        6      0.0          0.0          0.0           0.0          0.0          0.0           0.0          0.0          0.0
                #               0.0          0.0          0.0           0.0          0.0          0.0           0.0          0.0          0.0
                struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'18f')
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]
                    out = struct1.unpack(edata)

                    (eid_device, fx, fy, fz, otx, oty, otz, etx, ety, etz,
                                 mx, my, mz, orx, ory, orz, erx, ery, erz) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, op2.nonlinear_factor, op2.sort_method)
                    if op2.is_debug_file:
                        op2.binary_debug.write('%s - %s\n' % (name, str(out)))
                    obj.add_sort1(dt, eid, fx, fy, fz, otx, oty, otz, etx, ety, etz,
                                  mx, my, mz, orx, ory, orz, erx, ery, erz)
                    n += ntotal
        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
        return n, nelements, ntotal

    def _oes_cplstn_nx(self, data, ndata: int, dt, is_magnitude_phase: bool,
                       result_type: str, prefix: str, postfix: str):
        """
        reads stress/strain for element type:
        - 271 CPLSTN3
        - 272 CPLSTN4
        - 273 CPLSTN6
        - 274 CPLSTN8

        """
        op2 = self.op2
        n = 0
        if op2.is_stress:
            stress_strain = 'stress'
            obj_vector_real = RealCPLSTRNPlateStressNXArray
        else:
            stress_strain = 'strain'
            obj_vector_real = RealCPLSTRNPlateStrainNXArray

        if op2.element_type == 271:    # CPLSTN3
            nnodes = 3
            nnodes_cen = 1
        elif op2.element_type == 272:    # CPLSTN4
            nnodes = 4
            nnodes_cen = 5
        elif op2.element_type == 273:  # CPLSTN6
            nnodes = 6
            nnodes_cen = 4
        elif op2.element_type == 274:  # CPLSTN8
            nnodes = 8
            nnodes_cen = 5
        else:  # pragma: no cover
            #raise RuntimeError(op2.element_type)
            raise RuntimeError(op2.code_information())
        result_name = f'{stress_strain}.cplstn{nnodes}_{stress_strain}'
        name = f'CPLSTN{nnodes}'

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        num_wide = op2.num_wide
        ntotal = num_wide * self.size
        num_wide_real = 32  # CPLSTN4, CPLSTN3
        num_wide_real = 26  # CPLSTN6
        nelements = ndata // ntotal
        assert ndata % ntotal == 0
        #print(f'result_name = {result_name}')

        if result_type == 0 and num_wide == 6 and nnodes == 3:  # real; CPLSTN3
            nlayers = nelements
            auto_return, is_vectorized = op2._create_oes_object4(
                nlayers, result_name, slot, obj_vector_real)
            if auto_return:
                assert ntotal == op2.num_wide * 4
                return nelements * ntotal, None, None

            obj: RealCPLSTRNPlateStressNXArray = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * op2.num_wide * 4
                istart = obj.itotal
                iend = istart + nlayers
                obj._times[obj.itime] = dt

                if obj.itime == 0:
                    ints1 = frombuffer(data, dtype=op2.idtype).reshape(nelements, 6).copy()
                    eids = ints1[:, 0] // 10
                    obj.element[istart:iend] = eids
                    obj.element_node[istart:iend, 0] = eids
                    obj.element_node[istart:iend, 1] = 0
                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 6)[:, 1:]
                #[oxx, oyy, ozz, txy, von_mises]
                obj.data[obj.itime, istart:iend, :] = floats.copy()
            else:
                n = oes_cplstn3_real_6(
                    op2, data, obj,
                    nelements, ntotal, dt)

        elif result_type == 0 and num_wide == 32 and nnodes in {4, 8}:  # real; CPLSTN4
            nlayers = nelements * nnodes_cen
            auto_return, is_vectorized = op2._create_oes_object4(
                nlayers, result_name, slot, obj_vector_real)
            if auto_return:
                op2._data_factor = nnodes_cen
                assert ntotal == op2.num_wide * 4
                return nelements * ntotal, None, None

            obj: RealCPLSTRNPlateStressNXArray = op2.obj
            #if op2.is_debug_file:
                #op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                #op2.binary_debug.write('  element1 = [eid_device, layer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)]\n')
                #op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * op2.num_wide * 4
                iestart = obj.ielement
                ieend = iestart + nelements
                istart = obj.itotal
                iend = istart + nlayers
                obj._times[obj.itime] = dt

                if obj.itime == 0:
                    ints1 = frombuffer(data, dtype=op2.idtype).reshape(nelements, 32).copy()
                    ints2 = ints1[:, 2:].reshape(nlayers, 6)
                    eids = ints1[:, 0] // 10
                    nids = ints2[:, 0]
                    obj.element[iestart:ieend] = eids
                    eids2 = np.vstack([eids] * nnodes_cen).T.ravel()
                    obj.element_node[istart:iend, 0] = eids2
                    obj.element_node[istart:iend, 1] = nids
                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 32)[:, 2:]
                floats2 = floats.reshape(nlayers, 6)
                #[oxx, oyy, ozz, txy, von_mises]
                obj.data[obj.itime, istart:iend, :] = floats2[:, 1:].copy()
            else:
                n = oes_cplstn4_real_32(
                    op2, data, obj,
                    nelements, ntotal, dt)

        elif result_type == 0 and num_wide == 26 and nnodes == 6:  # real; CPLSTN6
            nlayers = nelements * nnodes_cen
            auto_return, is_vectorized = op2._create_oes_object4(
                nlayers, result_name, slot, obj_vector_real)
            if auto_return:
                op2._data_factor = nnodes_cen
                assert ntotal == op2.num_wide * 4
                return nelements * ntotal, None, None

            obj: RealCPLSTRNPlateStressNXArray = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * op2.num_wide * 4
                iestart = obj.ielement
                ieend = iestart + nelements
                istart = obj.itotal
                iend = istart + nlayers
                obj._times[obj.itime] = dt

                if obj.itime == 0:
                    ints1 = frombuffer(data, dtype=op2.idtype).reshape(nelements, 26).copy()
                    ints2 = ints1[:, 2:].reshape(nlayers, 6)
                    eids = ints1[:, 0] // 10
                    nids = ints2[:, 0]
                    obj.element[iestart:ieend] = eids
                    eids2 = np.vstack([eids] * nnodes_cen).T.ravel()
                    obj.element_node[istart:iend, 0] = eids2
                    obj.element_node[istart:iend, 1] = nids
                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 26)[:, 2:]
                floats2 = floats.reshape(nlayers, 6)
                #[oxx, oyy, ozz, txy, von_mises]
                obj.data[obj.itime, istart:iend, :] = floats2[:, 1:].copy()
            else:
                n = oes_cplstn6_real_26(
                    op2, data, obj,
                    nelements, ntotal, dt)

        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
        return n, nelements, ntotal

    def _oes_composite_solid_nx(self, data, ndata: int, dt, is_magnitude_phase: bool,
                                result_type: str, prefix: str, postfix: str):
        """
        269: Composite HEXA element (CHEXAL)
        270: Composite PENTA element (CPENTAL)

        NX PCOMPS (CHEXAL, CPENTAL)
        Linear

        Etype        Corner Center
        ======       ====== ======
        CHEXAL-269   11     43
        CPENTAL-???  11???  3+6*5???

        TCODE,7 =0 Real
        Q4CSTR=0 Center option
           2 PLY I Lamina number (11)
           3 FLOC CHAR4 Fiber location (BOT, MID, TOP)
           4 GRID I Edge grid ID (center=0)
           5 E11 RS Normal strain in the 1-direction
           6 E22 RS Normal strain in the 2-direction
           7 E33 RS Normal strain in the 3-direction
           8 E12 RS Shear strain in the 12-plane
           9 E23 RS Shear strain in the 23-plane
           10 E13 RS Shear strain in the 13-plane
           11 ETMAX1 RS von Mises strain
        Q4CSTR=1 Center and Corner option (3+8*5=43)
           2 PLY I Lamina number
           3 FLOC CHAR4 Fiber location (BOT, MID, TOP)
           4 GRID I Edge grid ID (center=0)
           5 E11 RS Normal strain in the 1-direction
           6 E22 RS Normal strain in the 2-direction
           7 E33 RS Normal strain in the 3-direction
           8 E12 RS Shear strain in the 12-plane
           9 E23 RS Shear strain in the 23-plane
           10 E13 RS Shear strain in the 13-plane
           11 ETMAX1 RS Von Mises strain
           For each fiber location requested (PLSLOC), words 4 through 11 repeat 5 times.

        Complex
        TCODE,7 =1 Real/imaginary
        Q4CSTR=0 Center option
          2 PLY I Lamina number
          3 FLOC CHAR4 Fiber location (BOT, MID, TOP)
          4 GRID I Edge grid ID (Center = 0)
          5 E11r RS Normal strain in the 1-direction, real part
          6 E22r RS Normal strain in the 2-direction, real part
          7 E33r RS Normal strain in the 3-direction, real part
          8 E12r RS Shear strain in the 12-plane, real part
          9 E23r RS Shear strain in the 23-plane, real part
          10 E13r RS Shear strain in the 13-plane, real part
          11 E11i RS Normal strain in the 1-direction, imaginary part
          12 E22i RS Normal strain in the 2-direction, imaginary part
          13 E33i RS Normal strain in the 3-direction, imaginary part
          14 E12i RS Shear strain in the 12-plane, imaginary part
          15 EL23i RS Shear strain in the 23-plane, imaginary part
          16 EL13i RS Shear strain in the 13-plane, imaginary part
        Q4CSTR=1 Center and Corner option
          2 PLY I Lamina number
          3 FLOC CHAR4 Fiber location (BOT, MID, TOP)
          4 GRID I Edge grid ID (Center = 0)
          5 E11r RS Normal strain in the 1-direction, real part
          6 E22r RS Normal strain in the 2-direction, real part
          7 E33r RS Normal strain in the 3-direction, real part
          8 E12r RS Shear strain in the 12-plane, real part
          9 E23r RS Shear strain in the 23-plane, real part
          10 E13r RS Shear strain in the 13-plane, real part
          11 E11i RS Normal strain in the 1-direction, imaginary part
          12 E22i RS Normal strain in the 2-direction, imaginary part
          13 E33i RS Normal strain in the 3-direction, imaginary part
          14 E12i RS Shear strain in the 12-plane, imaginary part
          15 E23i RS Shear strain in the 23-plane, imaginary part
          16 E13i RS Shear strain in the 13-plane, imaginary part
          For each fiber location requested (PLSLOC), words 4 through 16 repeat 5 times.
        """
        op2 = self.op2
        n = 0
        #assert op2.is_stress is True, op2.code_information()

        stress_strain = 'stress' if op2.is_stress else 'strain'
        if op2.element_type == 269:
            result_name = f'{stress_strain}.chexa_composite_{stress_strain}'
        elif op2.element_type == 270:
            result_name = f'{stress_strain}.cpenta_composite_{stress_strain}'
        else:
            raise NotImplementedError(op2.code_information())

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        #if result_type == 0 and op2.num_wide == 43:  # real
            #op2.log.warning(f'skipping corner option for composite solid-{op2.element_name}-{op2.element_type}')
            #struct9 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'i 4s i 5f', self.size)) # 9
        obj_vector_real = RealSolidCompositeStressArray if op2.is_stress else RealSolidCompositeStrainArray
        #obj_vector_real = RealSolidCompositeStressArray

        if result_type == 0 and op2.num_wide == 11:  # real; center
            #op2.log.warning(f'skipping center option for composite solid-{op2.element_name}-{op2.element_type}')
            ntotal = 44 * self.factor # 11 * 4
            nelements = ndata // ntotal

            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                assert ntotal == op2.num_wide * 4
                return nelements * ntotal, None, None

            obj = op2.obj  # type: RealCompositeSolidStressArray
            n = oes_composite_solid_nx_real_center(op2, data, obj, nelements, ntotal)

        elif result_type == 0 and op2.num_wide == 43:  # real; center
            #op2.log.warning(f'skipping center option for composite solid-{op2.element_name}-{op2.element_type}')
            ntotal = 172 * self.factor  # 43*4
            nelements = ndata // ntotal

            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                assert ntotal == op2.num_wide * 4
                return nelements * ntotal, None, None

            obj = op2.obj  # type: RealCompositeSolidStressArray
            n = oes_composite_solid_nx_real_172(op2, data, obj, nelements, ntotal)
        else:
            raise NotImplementedError(op2.code_information())
        return n, nelements, ntotal

    def _oes_cbend(self, data, ndata, dt, is_magnitude_phase,
                   result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 69 : CBEND

        """
        op2 = self.op2
        if op2.is_stress:
            result_name = prefix + 'cbend_stress' + postfix
            obj_vector_real = RealBendStressArray
            obj_vector_complex = ComplexBendStressArray
            obj_vector_random = RandomBendStressArray
        else:
            result_name = prefix + 'cbend_strain' + postfix
            obj_vector_real = RealBendStrainArray
            obj_vector_complex = ComplexBendStrainArray
            obj_vector_random = RandomBendStrainArray

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        #print(op2.code_information())
        if result_type == 0 and op2.num_wide == 21:  # real
            #TCODE,7 =0 Real
            #2 GRID I External Grid Point identification number
            #3 CA RS Circumferential Angle
            #4 EC RS Long. strain at Point C
            #5 ED RS Long. strain at Point D
            #6 EE RS Long. strain at Point E
            #7 EF RS Long. strain at Point F
            #8 EMAX RS Maximum strain
            #9 EMIN RS Minimum strain
            #10 MST RS Margin of Safety in Tension
            #11 MSC RS Margin of Safety in Compression
            #Words 2 through 11 repeat 002 times
            n = 0
            ntotal = 84 * self.factor  # 4*21
            nelements = ndata // ntotal
            assert ndata % ntotal == 0, 'ndata=%s ntotal=%s nelements=%s error=%s' % (ndata, ntotal, nelements, ndata % ntotal)

            #nlayers = nelements * 2
            if op2.sort_method == 2:
                op2.log.warning('real cbend stress/strain for SORT2 is not supported')
                print(op2.code_information())
                return nelements * ntotal, None, None

            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            assert obj is not None
            if op2.use_vector and is_vectorized and op2.sort_method == 1 and 0:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 4)
                itime = obj.itime
                obj._times[itime] = dt
                if itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, 4)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                #[max_strain, avg_strain, margin]
                obj.data[itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                n = oes_cbend_real_21(op2, data, obj,
                                      nelements, ntotal, dt)

            #msg = ''
            #if op2.read_mode == 2:
                #msg = op2.code_information()
            #n = op2._not_implemented_or_skip(data, ndata, msg)
            #return n, None, None
        elif result_type == 1 and op2.num_wide == 21:  # complex
            n = 0
            ntotal = 84 * self.factor  # 4*21
            nelements = ndata // ntotal
            assert ndata % ntotal == 0, 'ndata=%s ntotal=%s nelements=%s error=%s' % (ndata, ntotal, nelements, ndata % ntotal)
            #TCODE,7 =1 Real / Imaginary
            #2 GRID I External Grid Point identification number
            #3 CA RS Circumferential Angle
            #4 SCR RS Long. Stress at Point C
            #5 SDR RS Long. Stress at Point D
            #6 SER RS Long. Stress at Point E
            #7 SFR RS Long. Stress at Point F
            #8 SCI RS Long. Stress at Point C
            #9 SDI RS Long. Stress at Point D
            #10 SEI RS Long. Stress at Point E
            #11 SFI RS Long. Stress at Point F
            #Words 2 through 11 repeat 002 times

            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_complex)
            if auto_return:
                assert ntotal == op2.num_wide * 4
                return nelements * ntotal, None, None

            obj = op2.obj
            assert obj is not None
            if op2.use_vector and is_vectorized and op2.sort_method == 1 and 0:
                n = nelements * 4 * op2.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 4)
                itime = obj.itime
                obj._times[itime] = dt

                if itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, 4)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                #[max_strain, avg_strain, margin]
                obj.data[itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                n = oes_cbend_complex_21(op2, data, obj, nelements, ntotal,
                                         is_magnitude_phase)

        elif result_type == 2 and op2.num_wide == 13:
            n = 0
            ntotal = 52 * self.factor  # 4*13
            nelements = ndata // ntotal
            #TCODE,7 =2 Real
            #2 GRID I External Grid Point identification number
            #3 CA RS Circumferential Angle
            #4 SC RS Long. Stress at Point C
            #5 SD RS Long. Stress at Point D
            #6 SE RS Long. Stress at Point E
            #7 SF RS Long. Stress at Point F
            #Words 2 through 7 repeat 002 times
            #if op2.table_name != "OESPSD2":
                #msg = ''
                #if op2.read_mode == 2:
                    #msg = op2.code_information()
                #n = op2._not_implemented_or_skip(data, ndata, msg)
                #return n, None, None

            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_random)
            if auto_return:
                assert ntotal == op2.num_wide * 4
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1 and 0:
                n = nelements * 4 * op2.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 4)
                itime = obj.itime
                obj._times[itime] = dt
                if itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, 4)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                #[max_strain, avg_strain, margin]
                obj.data[itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                ntotali = 24
                struct1 = Struct(op2._endian + op2._analysis_code_fmt)
                struct2 = Struct(op2._endian + b'i5f')

                for unused_i in range(nelements):
                    edata = data[n:n + 4]
                    #self.show_data(edata)
                    eid_device, = struct1.unpack(edata)
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, op2.nonlinear_factor, op2.sort_method)

                    n += 4
                    for unused_i in range(2):
                        edata = data[n:n + ntotali]
                        out = struct2.unpack(edata)
                        if op2.is_debug_file:
                            op2.binary_debug.write('BEND-69 - eid=%s dt=%s %s\n' % (eid, dt, str(out)))
                        #print('BEND-69 - eid=%s dt=%s %s\n' % (eid, dt, str(out)))

                        (grid, angle, sc, sd, se, sf) = out
                        obj.add_sort1(dt, eid, grid, angle, sc, sd, se, sf)
                        n += ntotali

        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
        return n, nelements, ntotal

    def _oes_cgap_nonlinear(self, data, ndata, dt, is_magnitude_phase,
                            result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 86 : GAPNL
        """
        op2 = self.op2
        n = 0
        if op2.is_stress:
            result_name = prefix + 'cgap_stress' + postfix # nonlinear_
        else:
            result_name = prefix + 'cgap_strain' + postfix # nonlinear_

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        if result_type == 0 and op2.num_wide == 11:  # real?
            if op2.is_stress:
                obj_vector_real = NonlinearGapStressArray
            else:
                raise NotImplementedError('NonlinearGapStrain')

            ntotal = 44 * self.factor # 4*11
            nelements = ndata // ntotal
            assert ndata % ntotal == 0
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return nelements * ntotal, None, None

            obj = op2.obj
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal

                ielement = obj.ielement
                ielement2 = ielement + nelements
                obj._times[obj.itime] = dt

                self.obj_set_element(obj, ielement, ielement2, data, nelements)

                #if obj.itime == 0:
                    #ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, 11).copy()
                    #eids = ints[:, 0] // 10
                    #obj.element[ielement:ielement2] = eids

                #print(data, len(data))
                if self.size == 4:
                    strings = frombuffer(data, dtype='|S4').reshape(nelements, 11)[:, [9, 10]]
                else:
                    strings = frombuffer(data, dtype='|S8').reshape(nelements, 11)[:, [9, 10]]
                form = [reshape_bytes_block_strip(string[0] + string[1], size=op2.size) for string in strings]

                #form = np.concatenate(strings) # strings[:, 0] + strings[:, 1] # np.hstack(strings)
                #form = np.char.join(strings[:, 0], strings[:, 1]) # np.hstack(strings)
                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 11)
                obj.form[ielement:ielement2] = form
                # skipping [form1, form2]
                #[cpx, shy, shz, au, shv, shw, slv, slp]
                obj.data[obj.itime, ielement:ielement2, :] = floats[:, 1:9].copy()
            else:
                n = oes_cgapnl_real_11(op2, data, obj, nelements, ntotal)
        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
        return n, nelements, ntotal

    def _oes_cbeam_nonlinear(self, data, ndata, dt, is_magnitude_phase,
                             result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 94 : BEAMNL

        """
        op2 = self.op2
        n = 0
        numwide_real = 51
        numwide_random = 0

        if op2.is_stress:
            result_name = prefix + 'cbeam_stress' + postfix
        else:
            result_name = prefix + 'cbeam_strain' + postfix

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        if result_type == 0 and op2.num_wide == numwide_real:
            msg = result_name
            if op2.is_stress:
                obj_vector_real = RealNonlinearBeamStressArray
            else:
                raise NotImplementedError('Nonlinear CBEAM Strain...this should never happen')

            ntotal = numwide_real * 4 * self.factor  # 204
            nelements = ndata // ntotal

            nlayers = nelements * 8
            auto_return, is_vectorized = op2._create_oes_object4(
                nlayers, result_name, slot, obj_vector_real)
            if auto_return:
                op2._data_factor = 8
                return ndata, None, None
            obj = op2.obj
            if op2.is_debug_file:
                op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                #op2.binary_debug.write('  #elementi = [eid_device, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,\n')
                #op2.binary_debug.write('                           s1b, s2b, s3b, s4b, smaxb, sminb,        MSc]\n')
                #op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if self.size == 4:
                struct1 = Struct(op2._endian + b'2i 4s5f 4s5f 4s5f 4s5f i 4s5f 4s5f 4s5f 4s5f')  # 2 + 6*8 + 1 = 51
            else:
                assert self.size == 8, self.size
                struct1 = Struct(op2._endian + b'2q 8s5d 8s5d 8s5d 8s5d q 8s5d 8s5d 8s5d 8s5d')  # 2 + 6*8 + 1 = 51

            for unused_i in range(nelements):  # num_wide=51
                edata = data[n:n + ntotal]
                out = struct1.unpack(edata)

                if op2.is_debug_file:
                    op2.binary_debug.write('BEAMNL-94 - %s\n' % str(out))

                #gridA, CA, long_CA, eqS_CA, tE_CA, eps_CA, ecs_CA,
                #       DA, long_DA, eqS_DA, tE_DA, eps_DA, ecs_DA,
                #       EA, long_EA, eqS_EA, tE_EA, eps_EA, ecs_EA,
                #       FA, long_FA, eqS_FA, tE_FA, eps_FA, ecs_FA,
                #gridB, CB, long_CB, eqS_CB, tE_CB, eps_CB, ecs_CB,
                #       DB, long_DB, eqS_DB, tE_DB, eps_DB, ecs_DB,
                #       EB, long_EB, eqS_EB, tE_EB, eps_EB, ecs_EB,
                #       FB, long_FB, eqS_FB, tE_FB, eps_FB, ecs_FB,
                # A
                assert out[3-1].rstrip() == b'   C', out[3-1]
                assert out[9-1].rstrip() == b'   D', out[9-1]
                assert out[15-1].rstrip() == b'   E', out[15-1]
                assert out[21-1].rstrip() == b'   F', out[21-1]

                # B
                assert out[28-1].rstrip() == b'   C', out[28-1]
                assert out[34-1].rstrip() == b'   D', out[34-1]
                assert out[40-1].rstrip() == b'   E', out[40-1]
                assert out[46-1].rstrip() == b'   F', out[46-1]

                eid_device = out[0]
                eid, dt = get_eid_dt_from_eid_device(
                    eid_device, op2.nonlinear_factor, op2.sort_method)
                obj.add_new_eid_sort1(dt, eid, *out[1:])
                n += ntotal

        elif result_type == 2 and op2.num_wide == numwide_random:  # random
            msg = op2.code_information()
            raise NotImplementedError(msg)
            #return op2._not_implemented_or_skip(data, ndata, msg)
        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
        return n, nelements, ntotal

    def _oes_cbar_100(self, data: bytes, ndata: int, dt, is_magnitude_phase: bool,
                      result_type: str, prefix: str, postfix: str):
        """
        reads stress/strain for element type:
         - 100 : BARS
        """
        op2 = self.op2
        n = 0
        if op2.is_stress:
            result_name = prefix + 'cbar_stress_10nodes' + postfix
        else:
            result_name = prefix + 'cbar_strain_10nodes' + postfix

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        if result_type == 0 and op2.num_wide == 10:  # real
            if op2.is_stress:
                obj_vector_real = RealBar10NodesStressArray
            else:
                obj_vector_real = RealBar10NodesStrainArray

            ntotal = 10 * 4 * self.factor
            nelements = ndata // ntotal

            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return ndata, None, None

            if op2.is_debug_file:
                op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                op2.binary_debug.write('  #elementi = [eid_device, sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS]\n')
                op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)
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

                self.obj_set_element(obj, istart, iend, data, nelements)

                floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 10)
                #[sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS]
                obj.data[obj.itime, istart:iend, :] = floats[:, 1:].copy()
            else:
                n = oes_cbar100_real_10(op2, data, obj, nelements, ntotal, dt)

        elif result_type == 1 and op2.num_wide == 16:  # complex
            msg = op2.code_information()
            return op2._not_implemented_or_skip(data, ndata, msg), None, None
        else:  # pragma: no cover
            raise RuntimeError(op2.code_information())
        return n, nelements, ntotal

    def _oes_hyperelastic_quad(self, data, ndata, dt, unused_is_magnitude_phase,
                               result_type, prefix, postfix):
        """
        139-QUAD4FD
        """
        op2 = self.op2
        #if op2.is_stress:
        result_name_base = 'hyperelastic_cquad4_'
        if op2.is_stress:
            flag = 'stress'
        else:
            flag = 'strain'
        result_name = prefix + result_name_base + flag + postfix

        if op2._results.is_not_saved(result_name):
            return ndata, None, None

        if result_type == 0 and op2.num_wide == 30:
            obj_vector_real = HyperelasticQuadArray

            op2._results._found_result(result_name)
            slot = op2.get_result(result_name)
            #op2.create_transient_object(result_name, slot, obj_vector_real)

            ntotal = 120 * self.factor # 36+28*3
            nelements = ndata // ntotal

            #print(op2.code_information())
            #print(op2.table_name_str)
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)

            if auto_return:
                op2._data_factor = 4  # number of "layers" for an element
                return nelements * ntotal, None, None
                #return ndata, None, None

            #if op2.is_debug_file:
                #op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                ##op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                #op2.binary_debug.write('  #elementi = [eid_device, sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS]\n')
                #op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)
            obj = op2.obj

            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                # self.itime = 0
                # self.ielement = 0
                # self.itotal = 0
                #self.ntimes = 0
                #self.nelements = 0
                n = nelements * ntotal

                istart = obj.itotal
                iend = istart + nelements * 4
                obj._times[obj.itime] = dt

                #if obj.itime == 0:
                # 30 = 2 + 28 = 2 + 7*4
                ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, 30).copy()
                #strs = frombuffer(data, dtype=self.sdtype)
                ints2 = ints[:, 2:].reshape(nelements * 4, 7)

                #strings = frombuffer(data, dtype=???)
                eids = ints[:, 0] // 10
                nids = ints2[:, 0]

                eids2 = np.vstack([eids, eids, eids, eids]).T.ravel()
                obj.element_node[istart:iend, 0] = eids2
                obj.element_node[istart:iend, 1] = nids
                #obj.element[istart:iend] = eids

                # dropping off eid and the string word (some kind of Type)
                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 30)[:, 2:].copy()
                floats2 = floats.reshape(nelements * 4, 7)
                #[oxx, oyy, txy, angle, majorp, minorp]
                obj.data[obj.itime, istart:iend, :] = floats2[:, 1:]
            else:
                n = 0
                # (2 + 7*4)*4 = 30*4 = 120
                ntotal1 = 36 * self.factor  # 4*9
                ntotal2 = 28 * self.factor  # 4*7
                s1 = Struct(op2._endian + op2._analysis_code_fmt + b'4s i6f')  # 1 + 4+1+6 = 12
                s2 = Struct(op2._endian + b'i6f')
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal1]
                    out = s1.unpack(edata)
                    if op2.is_debug_file:
                        op2.binary_debug.write('CQUAD4FD-139A- %s\n' % (str(out)))

                    (eid_device, etype, nid, sx, sy, sxy, angle, smj, smi) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, op2.nonlinear_factor, op2.sort_method)

                    obj._add_new_eid_sort1(dt, eid, etype, nid, sx, sy, sxy, angle, smj, smi)
                    n += ntotal1

                    for unused_i in range(3):  # TODO: why is this not 4?
                        edata = data[n:n + ntotal2]
                        out = s2.unpack(edata)
                        if op2.is_debug_file:
                            op2.binary_debug.write('               %s\n' % (str(out)))
                        (nid, sx, sy, sxy, angle, smj, smi) = out
                        obj._add_sort1(dt, eid, etype, nid, sx, sy, sxy, angle, smj, smi)
                        n += ntotal2
        else:
            raise RuntimeError(op2.code_information())
            #msg = 'numwide=%s element_num=%s etype=%s' % (
                #op2.num_wide, op2.element_type, op2.element_name)
            #return op2._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oes_plate_stress_34(self, data, ndata, unused_dt, unused_is_magnitude_phase,
                             unused_stress_name, unused_prefix, unused_postfix):
        """
        271-CPLSTN3
        275-CPLSTS3
        """
        op2 = self.op2
        msg = op2.code_information()
        return op2._not_implemented_or_skip(data, ndata, msg), None, None
        #if op2.element_type == 271:
            #result_name = 'cplstn3'
            #unused_nnodes = 1
            #ntotal = 4 * 6
        #elif op2.element_type == 275:
            #result_name = 'cplsts3'
            #unused_nnodes = 1
            #ntotal = 4 * 6
        #else:  # pragma: no cover
            #raise RuntimeError(op2.code_information())
        #if op2.is_stress:
            #obj_vector_real = RealCPLSTRNPlateStressArray
            #result_name += '_stress'
        #else:
            #obj_vector_real = RealCPLSTRNPlateStrainArray
            #result_name += '_strain'

        #numwide_real = ntotal // 4
        #if op2.format_code == 1 and op2.num_wide == numwide_real:
            ##ntotal = 4 * (1 + 6 * (nnodes))
            #nelements = ndata // ntotal

            ##op2._data_factor = 10  # TODO: why is this 10?
            #if op2.is_stress:
                #obj_vector_real = RealCPLSTRNPlateStressArray
                ##result_name = 'cplstn3_stress'
            #else:
                #obj_vector_real = RealCPLSTRNPlateStressArray
                ##result_name = 'cplstn3_strain'
            #slot = op2.get_result(result_name)

            #auto_return, is_vectorized = op2._create_oes_object4(
                #nelements, result_name, slot, obj_vector_real)
            #if auto_return:
                #assert ntotal == op2.num_wide * 4
                #return nelements * ntotal, None, None

            #obj = op2.obj
            ##if op2.use_vector and is_vectorized and op2.sort_method == 1:
            ##n = nelements * op2.num_wide * 4

            #istart = obj.itotal
            #iend = istart + nelements
            #obj._times[obj.itime] = dt

            #self.obj_set_element(obj, istart, iend, data, nelements)
            #floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, numwide_real)
            #results = floats[:, 1:].copy()
            ##print('results.shape', results.shape)

            ##[oxx, oyy, ozz, txy, ovm]
            #obj.data[obj.itime, istart:iend, :] = results
        #else:
            #msg = 'sort1 Type=%s num=%s' % (op2.element_name, op2.element_type)
            #return op2._not_implemented_or_skip(data, ndata, msg), None, None

    def _oes_plate_stress_68(self, data, ndata, unused_dt, unused_is_magnitude_phase,
                             unused_stress_name, unused_prefix, unused_postfix):
        # 276-CPLSTS4
        # 277-CPLSTS6
        # 278-CPLSTS8
        op2 = self.op2
        msg = op2.code_information()
        return op2._not_implemented_or_skip(data, ndata, msg), None, None
        #if op2.element_type == 276:
            #result_name = 'cplsts4'
            #nnodes = 5  # 4 + 1
            #ntotal = 4 * 32
        #elif op2.element_type == 277:
            #result_name = 'cplsts6'
            #nnodes = 4
            #ntotal = 4 * 26
        #elif op2.element_type == 278:
            #result_name = 'cplsts8'
            #nnodes = 5
            #ntotal = 4 * 32
        #else:
            #raise RuntimeError(op2.code_information())

        #if op2.is_stress:
            #obj_vector_real = RealCPLSTRNPlateStressArray
            #result_name += '_stress'
        #else:
            #obj_vector_real = RealCPLSTRNPlateStrainArray
            #result_name += '_strain'

        #numwide_real = 2 + 6 * (nnodes)
        #assert ntotal // 4 == numwide_real, 'notal/4=%s numwide_real=%s\n%s' % (
            #ntotal // 4, numwide_real, op2.code_information())

        #ntotal = numwide_real * 4
        #if op2.format_code == 1 and op2.num_wide == numwide_real:
            #nelements = ndata // ntotal

            ##op2._data_factor = 10  # TODO: why is this 10?
            #if op2.is_stress:
                #obj_vector_real = RealCPLSTRNPlateStressArray
                ##result_name = 'cplstn3_stress'
            #else:
                #obj_vector_real = RealCPLSTRNPlateStressArray
                ##result_name = 'cplstn3_strain'
            #slot = getattr(op2, result_name)

            #nlayers = nelements * nnodes
            #auto_return, is_vectorized = op2._create_oes_object4(
                #nlayers, result_name, slot, obj_vector_real)
            #if auto_return:
                #op2._data_factor = nnodes
                #assert ntotal == op2.num_wide * 4
                #return nelements * ntotal, None, None

            #obj = op2.obj
            ##if op2.use_vector and is_vectorized and op2.sort_method == 1:
            #n = nlayers * op2.num_wide * 4

            #istart = obj.itotal
            #iend = istart + nlayers
            #obj._times[obj.itime] = dt

            #if obj.itime == 0:
                #print(frombuffer(data, dtype=op2.idtype).size)
                #print('nelements=%s numwide=%s' % (nelements, numwide_real))
                #ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, numwide_real)
                #eids = ints[:, 0] // 10
                ##obj.element[istart:iend] = eids

            #floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, numwide_real).copy()
            #print('floats[:, 2:].shape', floats[:, 2:].shape)
            #print('nnelements=%s nnodes=%s numwide//nodes=%s' % (nelements, nnodes, (numwide_real-2) / nnodes))
            #results = floats[:, 2:].reshape(nelements, nnodes * 6)

            ##[oxx, oyy, ozz, txy, ovm]
            #obj.data[obj.itime, istart:iend, :] = results
        #else:
            #msg = 'sort1 Type=%s num=%s' % (op2.element_name, op2.element_type)
            #return op2._not_implemented_or_skip(data, ndata, msg)

    def obj_set_element(self, obj, ielement, ielement2, data, nelements):
        op2 = self.op2
        if obj.itime == 0:
            ints = frombuffer(data, dtype=op2.idtype8).reshape(nelements, op2.num_wide).copy()
            eids = ints[:, 0] // 10
            assert eids.min() > 0, eids.min()
            obj.element[ielement:ielement2] = eids

def oes_cgapnl_real_11(op2: OP2, data: bytes,
                       obj: NonlinearGapStressArray,
                       nelements: int, ntotal: int) -> int:
    n = 0
    size = op2.size
    if size == 4:
        struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'8f8s')
    else:
        struct1 = Struct(op2._endian + mapfmt8(op2._analysis_code_fmt) + b'8d16s')

    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        #op2.show_data(edata, 'ifsdq')
        out = struct1.unpack(edata)  # num_wide=25
        (eid_device, cpx, shy, shz, au, shv, shw, slv, slp, form_bytes) = out
        form = reshape_bytes_block_strip(form_bytes, size=size)
        assert form in ['OPEN', 'SLIDE', 'STICK', 'SLIP'], form
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('CGAPNL-86 - %s\n' % str(out))
        #form = (form1 + form2).rstrip(b' ').decode('utf8')
        #print((dt, eid, cpx, shy, shz, au, shv, shw, slv, slp, form))
        add_sort_x(dt, eid, cpx, shy, shz, au, shv, shw, slv, slp, form)
        n += ntotal
    return n


def oes_csolid_random3(op2: OP2, data: bytes, obj: RandomSolidVMStressArray, nelements: int,
                       element_name: str, nnodes_expected: int,
                       ntotal: int) -> int:
    # print(msg)
    n = 0
    assert op2.size == 4, op2.size
    fmt = b'ii4si' + b'i7f' * nnodes_expected
    structi = Struct(fmt)

    for ielement in range(nelements):
        datai = data[n:n + ntotal]
        # op2.show_data(datai, types='if')
        out = structi.unpack(datai)
        eid_device, cid, grid, nnodes, *other = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #print(eid, cid, grid, nnodes)
        # [oxx, oyy, ozz, txy, tyz, txz, von_mises]
        # 6 SX   RS Normal in x
        # 7 SY   RS Normal in y
        # 8 SZ   RS Normal in z
        # 9 TXY  RS Shear in xy
        # 10 TYZ RS Shear in yz
        # 11 TZX RS Shear in zx
        # 12 RMSVM RS RMS von Mises
        for inode in range(nnodes_expected):
            #print('   ', other[8 * inode], other[8 * inode + 1:8 * (inode + 1)])
            node_id = other[8 * inode]
            oxx, oyy, ozz, txy, tyz, txz, vm = other[8 * inode + 1:8 * (inode + 1)]
            obj.add_node_sort1(dt, eid, inode, node_id, oxx, oyy, ozz, txy, tyz, txz, vm)
        n += ntotal
    return n

def oes_csolidnl_real(op2: OP2, data: bytes,
                      obj: RealNonlinearSolidArray,
                      etype: str, nnodes: int,
                      nelements: int, ntotal: int) -> int:
    n = 0
    size = op2.size
    if size == 4:
        s1 = Struct(op2._endian + op2._analysis_code_fmt + b'4s')
        s2 = Struct(op2._endian + b'i15f')
    else:
        s1 = Struct(op2._endian + mapfmt8(op2._analysis_code_fmt) + b'8s')
        s2 = Struct(op2._endian + b'q15d')
    n1 = 8 * op2.factor
    n2 = 64 * op2.factor

    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):  # 2+16*9 = 146 -> 146*4 = 584
        edata = data[n:n+n1]
        n += n1

        out = s1.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('%s-%s - %s\n' % (etype, op2.element_type, str(out)))

        (eid_device, unused_ctype) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #print('%s-%s -eid=%s dt=%s %s\n' % (etype, op2.element_type, eid, dt, str(out)))

        for unused_j in range(nnodes):
            edata = data[n:n+n2]
            n += n2
            out = s2.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('%s-%sB - %s\n' % (etype, op2.element_type, str(out)))
            #print('%s-%sB - %s\n' % (etype, op2.element_type, str(out)))
            assert len(out) == 16

            (grid,
             sx, sy, sz, sxy, syz, sxz, se, eps, ecs,
             ex, ey, ez, exy, eyz, exz) = out
            add_sort_x(dt, eid, grid,
                       sx, sy, sz, sxy, syz, sxz, se, eps, ecs,
                       ex, ey, ez, exy, eyz, exz)
    return n


def oes_cshellnl_real_13(op2: OP2, data: bytes,
                         obj: RealNonlinearPlateArray,
                         etype,
                         nelements: int, ntotal: int) -> int:
    n = 0
    #size = op2.size
    struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'12f')  # 1+12=13
    add_new_eid_sort_x = getattr(obj, 'add_new_eid_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('CQUADNL-90 - %s\n' % str(out))

        (eid_device, fd1,
         sx1, sy1, sz1, txy1, es1, eps1, ecs1,
         ex1, ey1, ez1, exy1) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        add_new_eid_sort_x(
            dt, eid, op2.element_type, fd1,
            sx1, sy1, sz1, txy1, es1, eps1, ecs1,
            ex1, ey1, ez1, exy1)
        n += ntotal
    return n

def oes_cshellnl_real_25(op2: OP2, data: bytes,
                         obj: RealNonlinearPlateArray,
                         etype,
                         nelements: int, ntotal: int) -> int:
    n = 0
    #size = op2.size
    struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'24f', op2.size)) # 1+24=25

    add_new_eid_sort_x = getattr(obj, 'add_new_eid_sort' + str(op2.sort_method))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)
        if op2.is_debug_file:
            eid = out[0] // 10
            op2.binary_debug.write('CQUADNL-90 - %s : %s\n' % (eid, str(out)))

        (eid_device,
         fd1, sx1, sy1, undef1, txy1, es1, eps1, ecs1, ex1, ey1, undef2, etxy1,
         fd2, sx2, sy2, undef3, txy2, es2, eps2, ecs2, ex2, ey2, undef4, etxy2) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        add_new_eid_sort_x(
            dt, eid, etype,
            fd1, sx1, sy1, undef1, txy1, es1, eps1, ecs1, ex1, ey1, undef2, etxy1)
        add_sort_x(
            dt, eid, etype,
            fd2, sx2, sy2, undef3, txy2, es2, eps2, ecs2, ex2, ey2, undef4, etxy2)
        n += ntotal
    return n


def oes_composite_solid_nx_real_center(op2: OP2, data: bytes,
                                       obj: Union[RandomCompositePlateStressArray, RandomCompositePlateStrainArray],
                                       nelements: int, ntotal: int) -> int:
    n = 0
    size = op2.size

    struct11 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'i 4s i 7f', size)) # 11
    #sort_method = op2.sort_method
    #add_eid_sort_x = getattr(obj, 'add_eid_sort' + str(op2.sort_method))
    #add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))

    for unused_i in range(nelements):
        edata = data[n:n+ntotal]  # 4*11
        out = struct11.unpack(edata)
        #print(out)
        (eid_device, layer, location_bytes, grid, o11, o22, o33, t12, t23, t13, ovm) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        location = location_bytes.strip().decode('latin1')
        assert location == 'MID', out
        #print(f'eid,layer=({eid},{layer}) location={location!r} grid={grid} o11={o11:g} o22={o22:g} o33={o33:g} t12={t12:g} t1z={t13:g} t2z={t23:g} ovm={ovm:g}')
        obj.add_sort1(dt, eid, layer, location, grid, o11, o22, o33, t12, t23, t13, ovm)
        n += ntotal
    return n

def oes_composite_solid_nx_real_172(op2: OP2, data: bytes,
                                    obj: Union[RealSolidCompositeStressArray, RealSolidCompositeStrainArray],
                                    nelements: int, ntotal: int) -> int:
    n = 0
    #size = op2.size

    structa = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'i 4s', op2.size)) # 3
    structb = Struct(op2._endian + mapfmt(b'i 7f', op2.size)) # 8
    #sort_method = op2.sort_method
    #add_eid_sort_x = getattr(obj, 'add_eid_sort' + str(op2.sort_method))
    #add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    ntotal1 = 12 * op2.factor # 4*3
    ntotal2 = 32 * op2.factor # 4*8
    for unused_i in range(nelements):
        edata = data[n:n+ntotal1]  # 4*3, 4*40 = 4*43
        #op2.show_data(edata)
        out = structa.unpack(edata)
        #(13, 1, b' MID')
        eid_device, layer, location_bytes = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        location = location_bytes.strip().decode('latin1')
        assert location == 'MID', out
        #print(out)
        n += ntotal1
        for unused_j in range(5):
            edata = data[n:n+ntotal2]
            out = structb.unpack(edata)

            #print('  %s' % str(out))
            (grid, o11, o22, o33, t12, t23, t13, ovm) = out
            #print(f'eid,layer=({eid},{layer}) location={location!r} grid={grid} o11={o11:g} o22={o22:g} o33={o33:g} t12={t12:g} t1z={t13:g} t2z={t23:g} ovm={ovm:g}')
            obj.add_sort1(dt, eid, layer, location, grid, o11, o22, o33, t12, t23, t13, ovm)
            n += ntotal2
            #print(out)
        #(eid_device, layer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm) = out
        #print(out)
        #n += ntotal
    return n

def oes_csolid_nonlinear_hyperelastic_real(op2: OP2, data: bytes,
                                           obj: RealNonlinearPlateArray,
                                           nnodes_expected: int, nelements: int, ntotal: int,
                                           element_name: str,
                                           preline1: str, preline2: str) -> int:
    if obj is not None:
        assert isinstance(obj, RealNonlinearPlateArray), type(obj)
    n = 0
    #size = op2.size

    # 2 TYPE CHAR4 Grid or Gaus
    #
    # 3 ID I
    # 4 SX RS
    # 5 SY RS
    # 6 SZ RS
    # 7 SXY RS
    # 8 SYZ RS
    # 9 SZX RS
    # 10 PRESSURE RS
    # 11 VOLSTR RS
    # 12 EX RS
    # 13 EY RS
    # 14 EZ RS
    # 15 EXY RS
    # 16 EYZ RS
    # 17 EZX RS
    # Words 3 through 17 repeat 008 times
    struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'4s')
    struct2 = Struct(op2._endian + b'i14f')
    if op2.is_debug_file:
        msg = (
            f'{op2.element_name}-{op2.element_type} nelements={nelements} '
            f'nnodes={nnodes_expected}; '
            'C=[sxx, syy, szz, txy, tyz, txz, pressure, '
            'evol, exx, eyy, ezz, exy, eyz, exz]\n')
        op2.binary_debug.write(msg)

    for unused_i in range(nelements):
        edata = data[n:n+8]
        out = struct1.unpack(edata)
        (eid_device, unused_abcd, ) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        if op2.is_debug_file:
            op2.binary_debug.write('%s - eid=%i; %s\n' % (preline1, eid, str(out)))
        #assert nnodes < 21, 'print_block(data[n:n+16])'  #self.print_block(data[n:n+16])

        n += 8
        for inode in range(nnodes_expected):  # nodes pts, no centroid
            out = struct2.unpack(data[n:n + 60]) # 4*15 = 60
            if op2.is_debug_file:
                op2.binary_debug.write('%s - %s\n' % (preline2, str(out)))
            (grid_device, sxx, syy, szz, txy, tyz, txz, pressure,
             evol, exx, eyy, ezz, exy, eyz, exz) = out
            #print(out)

            if op2.is_debug_file:
                op2.binary_debug.write('  eid=%s inode=%i; C=[%s]\n' % (
                    eid, grid_device, ', '.join(['%r' % di for di in out])))

            grid = grid_device
            if 0:  # pragma: no cover
                if inode == 0:
                    #  this is correct, but fails
                    #element_name = op2.element_name + str(nnodes)
                    obj.add_eid_sort1(element_name, cid, dt, eid, grid,
                                      sxx, syy, szz, txy, tyz, txz, ovm)
                else:
                    obj.add_node_sort1(dt, eid, inode, grid,
                                       sxx, syy, szz, txy, tyz, txz, ovm)
            n += 60
    return n

def oes_composite_shells_nx_random_7(op2: OP2, data: bytes,
                                     obj: Union[RandomCompositePlateStressArray, RandomCompositePlateStrainArray],
                                     nelements: int, ntotal: int) -> int:
    n = 0
    size = op2.size
    assert size == 4, size
    struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'i5f')
    add_sort_x_7words = getattr(obj, f'add_sort{op2.sort_method}_7words')
    #obj.add_sort1_7words
    #print(f'random nelements={nelements}')
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device, ply_id, oxx, oyy, txy, txz, tyz) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #print(eid, out)

        #if op2.is_debug_file:
            #op2.binary_debug.write('%s-%s - (%s) + %s\n' % (op2.element_name, op2.element_type, eid_device, str(out)))
        #print(obj)
        add_sort_x_7words(dt, eid, ply_id, oxx, oyy, txy, txz, tyz)
        n += ntotal
    #print(f'all_times n={len(obj.all_times)}')
    #print(f'all_elements n={len(obj.all_elements)}')
    return n

def oes_cbend_complex_21(op2: OP2, data: bytes,
                         obj: Union[ComplexBendStressArray, ComplexBendStrainArray],
                         nelements: int, ntotal: int, is_magnitude_phase: bool) -> int:
    n = 0
    #size = op2.size
    ntotali = 40
    struct1 = Struct(op2._endian + op2._analysis_code_fmt)
    struct2 = Struct(op2._endian + b'i9f')
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))

    for unused_i in range(nelements):
        edata = data[n:n + 4]
        eid_device, = struct1.unpack(edata)
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        n += 4
        for unused_j in range(2):
            edata = data[n:n + ntotali]
            out = struct2.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('BEND-69 - eid=%s %s\n' % (eid, str(out)))
            #print('BEND-69 - eid=%s %s\n' % (eid, str(out)))

            (grid, angle, scr, sdr, ser, sfr,
             sci, sdi, sei, sfi) = out

            if is_magnitude_phase:
                sc = polar_to_real_imag(scr, sci)
                sd = polar_to_real_imag(sdr, sdi)
                se = polar_to_real_imag(ser, sei)
                sf = polar_to_real_imag(sfr, sfi)
            else:
                sc = complex(scr, sci)
                sd = complex(sdr, sdi)
                se = complex(ser, sei)
                sf = complex(sfr, sfi)
            add_sort_x(dt, eid, grid, angle, sc, sd, se, sf)
            n += ntotali
    return n
