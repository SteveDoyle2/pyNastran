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
from struct import Struct
from typing import Tuple, Union, Any
from numpy import fromstring, frombuffer, radians, sin, cos, vstack, repeat, array
import numpy as np

from pyNastran.op2.op2_interface.op2_codes import SORT1_TABLES_BYTES, TABLES_BYTES
from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.op2_interface.op2_common import OP2Common
from pyNastran.op2.op2_interface.utils import apply_mag_phase, build_obj
from pyNastran.op2.op2_helper import polar_to_real_imag
from pyNastran.op2.op2_interface.function_codes import func1, func7

from pyNastran.op2.tables.utils import get_eid_dt_from_eid_device
from pyNastran.op2.tables.oes_stressStrain.real.oes_bars import RealBarStressArray, RealBarStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_bars100 import RealBar10NodesStressArray, RealBar10NodesStrainArray

from pyNastran.op2.tables.oes_stressStrain.real.oes_beams import (RealBeamStressArray, RealBeamStrainArray,
                                                                  RealNonlinearBeamStressArray)
from pyNastran.op2.tables.oes_stressStrain.real.oes_bush import RealBushStressArray, RealBushStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_bush1d import RealBush1DStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_composite_plates import RealCompositePlateStressArray, RealCompositePlateStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_gap import NonlinearGapStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_plates import RealPlateStressArray, RealPlateStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_plate_strain import RealCPLSTRNPlateStressArray, RealCPLSTRNPlateStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_rods import RealRodStressArray, RealRodStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_shear import RealShearStrainArray, RealShearStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_solids import RealSolidStrainArray, RealSolidStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_springs import (RealSpringStressArray, RealSpringStrainArray,
                                                                    RealNonlinearSpringStressArray)
from pyNastran.op2.tables.oes_stressStrain.real.oes_triax import RealTriaxStressArray, RealTriaxStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_bend import RealBendStressArray, RealBendStrainArray


from pyNastran.op2.tables.oes_stressStrain.complex.oes_bars import ComplexBarStressArray, ComplexBarStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_beams import ComplexBeamStressArray, ComplexBeamStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_bush import (ComplexCBushStressArray, ComplexCBushStrainArray)
from pyNastran.op2.tables.oes_stressStrain.complex.oes_bush1d import ComplexCBush1DStressArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_plates import (
    ComplexPlateStressArray, ComplexPlateStrainArray, ComplexLayeredCompositesArray)
from pyNastran.op2.tables.oes_stressStrain.complex.oes_plates_vm import (
    ComplexPlateVMStressArray, ComplexPlateVMStrainArray)
from pyNastran.op2.tables.oes_stressStrain.complex.oes_triax import ComplexTriaxStressArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_rods import ComplexRodStressArray, ComplexRodStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_shear import ComplexShearStressArray, ComplexShearStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_solids import ComplexSolidStressArray, ComplexSolidStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_springs import ComplexSpringStressArray, ComplexSpringStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_bend import ComplexBendStressArray, ComplexBendStrainArray

from pyNastran.op2.tables.oes_stressStrain.random.oes_rods import RandomRodStressArray, RandomRodStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_bars import RandomBarStressArray, RandomBarStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_beams import RandomBeamStressArray, RandomBeamStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_bend import RandomBendStressArray, RandomBendStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_plates import RandomPlateStressArray, RandomPlateStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_plates_vm import RandomPlateVMStressArray, RandomPlateVMStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_solids import RandomSolidStressArray, RandomSolidStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_shear import RandomShearStressArray, RandomShearStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_composite_plates import RandomCompositePlateStressArray, RandomCompositePlateStrainArray

from pyNastran.op2.tables.oes_stressStrain.oes_nonlinear_rod import RealNonlinearRodArray
from pyNastran.op2.tables.oes_stressStrain.oes_nonlinear_bush import RealNonlinearBushArray
from pyNastran.op2.tables.oes_stressStrain.oes_hyperelastic import (
    HyperelasticQuadArray)
from pyNastran.op2.tables.oes_stressStrain.oes_nonlinear import RealNonlinearPlateArray, RealNonlinearSolidArray

NX_TABLES_BYTES = [b'OESVM1', b'OESVM2']

class OES(OP2Common):
    """
    Defines  the OES class that is used to read stress/strain data
    """
    def __init__(self):
        OP2Common.__init__(self)
        self.ntotal = 0

    def _read_oes1_3(self, data, unused_ndata):
        """
        reads OES1 subtable 3
        """
        self._analysis_code_fmt = b'i'
        self._data_factor = 1
        self.words = [
            'aCode', 'tCode', 'element_type', 'isubcase',
            '???', '???', '???', 'load_set'
            'format_code', 'num_wide', 's_code', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', 'Title', 'subtitle', 'label', '???']

        self.parse_approach_code(data)  # 3

        ## element type
        self.element_type = self.add_data_parameter(data, 'element_type', b'i', 3, False)

        ## load set ID
        self.load_set = self.add_data_parameter(data, 'load_set', b'i', 8, False)

        ## format code
        self.format_code = self.add_data_parameter(data, 'format_code', b'i', 9, False)

        ## number of words per entry in record
        ## .. note:: is this needed for this table ???
        self.num_wide = self.add_data_parameter(data, 'num_wide', b'i', 10, False)

        ## stress/strain codes
        self.s_code = self.add_data_parameter(data, 's_code', b'i', 11, False)

        ## thermal flag; 1 for heat ransfer, 0 otherwise
        self.thermal = self.add_data_parameter(data, 'thermal', b'i', 23, False)

        ## assuming tCode=1
        analysis_code = self.analysis_code
        if analysis_code == 1:   # statics / displacement / heat flux
            ## load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5, False)
            self.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
            self.setNullNonlinearFactor()
        elif analysis_code == 2:  # real eigenvalues
            #: mode number
            self.mode = self.add_data_parameter(data, 'mode', b'i', 5)
            #: eigenvalue
            self.eign = self.add_data_parameter(data, 'eign', b'f', 6, False)
            #: mode or cycle TODO confused on the type - F1 means float/int???
            self.mode2 = self.add_data_parameter(data, 'mode2', b'i', 7, False)
            self.cycle = self.add_data_parameter(data, 'cycle', b'f', 7, False)
            self.update_mode_cycle('cycle')
            self.data_names = self.apply_data_code_value('data_names', ['mode', 'eign', 'mode2', 'cycle'])
        #elif analysis_code == 3: # differential stiffness
            #self.lsdvmn = self.get_values(data,'i',5) ## load set number
            #self.data_code['lsdvmn'] = self.lsdvmn
        #elif analysis_code == 4: # differential stiffness
        #    self.lsdvmn = self.get_values(data,'i',5) ## load set number
        elif analysis_code == 5:   # frequency
            ## frequency
            self.freq = self.add_data_parameter(data, 'freq', b'f', 5)
            self.data_names = self.apply_data_code_value('data_names', ['freq'])
        elif analysis_code == 6:  # transient
            ## time step
            self.dt = self.add_data_parameter(data, 'dt', b'f', 5)
            self.data_names = self.apply_data_code_value('data_names', ['dt'])
        elif analysis_code == 7:  # pre-buckling
            ## load set
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5)
            self.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
        elif analysis_code == 8:  # post-buckling
            ## mode number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5)
            self.eigr = self.add_data_parameter(data, 'eigr', b'f', 6, False)  # real eigenvalue
            self.data_names = self.apply_data_code_value('data_names', ['lsdvmn', 'eigr'])
        elif analysis_code == 9:  # complex eigenvalues
            ## mode number
            self.mode = self.add_data_parameter(data, 'mode', b'i', 5)
            ## real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', b'f', 6, False)
            ## imaginary eigenvalue
            self.eigi = self.add_data_parameter(data, 'eigi', b'f', 7, False)
            self.data_names = self.apply_data_code_value('data_names', ['mode', 'eigr', 'eigi'])
        elif analysis_code == 10:  # nonlinear statics
            ## load step
            self.lftsfq = self.add_data_parameter(data, 'lftsfq', b'f', 5)
            self.data_names = self.apply_data_code_value('data_names', ['lftsfq'])
        elif analysis_code == 11:  # old geometric nonlinear statics
            ## load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5)
            self.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
        elif analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## Time step ??? --> straight from DMAP
            self.dt = self.add_data_parameter(data, 'dt', b'f', 5)
            self.data_names = self.apply_data_code_value('data_names', ['dt'])
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % self.analysis_code
            raise RuntimeError(msg)
        # tcode=2
        #if self.analysis_code==2: # sort2
        #    self.lsdvmn = self.get_values(data,'i',5)

        #print(f'tCode={self.tCode} -> result_type={result_type}')
        #print(self.code_information())

        self.fix_format_code()
        self.get_oes_prefix_postfix()
        self._parse_thermal_code()
        self._set_force_stress_element_name()
        if self.is_debug_file:
            self.binary_debug.write('  element_name   = %r\n' % self.element_name)
            self.binary_debug.write('  approach_code  = %r\n' % self.approach_code)
            self.binary_debug.write('  tCode          = %r\n' % self.tCode)
            self.binary_debug.write('  isubcase       = %r\n' % self.isubcase)

        self._read_title(data)


        try:
            self.element_name = self.element_mapper[self.element_type]
        except KeyError:
            self.log.error(self.code_information())
            raise
        assert self.element_name != '', self.code_information()
        #if self.element_type not in self.element_mapper:
            #return self._not_implemented_or_skip(data, ndata, self.code_information())

        self._parse_stress_code_to_stress_bits()
        self._write_debug_bits()
        assert isinstance(self.format_code, int), self.format_code
        #print('self.nonlinear_factor =', self.nonlinear_factor)
        #assert self.num_wide != 146, self.code_information()
        #self._check_result_type()

    def _check_result_type(self):
        sort_method = func1(self.tCode)
        if sort_method == 1:  # SORT1
            assert self.sort_bits.is_sort1 == 1, self.code_information()
        elif sort_method == 2:  # SORT2
            assert self.sort_bits.is_sort1 == 0, self.code_information()
        else:
            raise NotImplementedError(sort_method)


        result_type = self.result_type #  func7(self.tCode)
        if result_type == 0:
            assert self.sort_bits.is_real == 1, self.code_information()
        elif result_type == 1:
            assert self.sort_bits.is_complex == 1, self.code_information()
        elif result_type == 2:
            assert self.sort_bits.is_random == 1, self.code_information()
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
        bits = [0, 0, 0, 0, 0]

        s_code = self.s_code
        i = 4
        while s_code > 0:
            value = s_code % 2
            s_code = (s_code - value) // 2
            bits[i] = value
            i -= 1
        self.stress_bits = bits
        self.data_code['stress_bits'] = self.stress_bits

    def _read_oes2_4(self, data, ndata):
        """
        Reads the Stress Table 4
        """
        if self.table_name in NX_TABLES_BYTES:
            self.to_nx(f'because table_name={self.table_name}')

        #assert self.isubtable == -4, self.isubtable
        #if self.is_debug_file:
            #self.binary_debug.write('  element_name = %r\n' % self.element_name)
        #print("element_name =", self.element_name)
        assert isinstance(self.format_code, int), self.format_code
        assert self.is_stress is True, self.code_information()
        self.data_code['is_stress_flag'] = True
        self.data_code['is_strain_flag'] = False

        self._setup_op2_subcase('STRESS/STRAIN')
        elements_to_read = [
            1, 3, 10, # CROD, CTUBE, CTUBE
            11, 12, 13, # CELAS1, CELAS2, CELAS3,
            2, 4, 34, 33, 74, # CBEAM, CSHEAR, CBAR, CQUAD4, CTRIA3,
            75, 64, 70, 82, 144, # CTRIA6, CQUAD8, CTRIAR, CQUADR, CQUAD4
            69, # CBEND
            67, 68, 95, 102, #  #CHEXA, CPENTA, QUAD4-comp, CBUSH
            39, #CTETRA
            86, # GAPNL
            88, # TRIA3-nonlinear
            89, # ROD-nonlinear
            90, # QUAD4-nonlinear
            91, # PENTANL
            93, # HEXANL
            97, # CTRIA3-C
            96, # QUAD8-nonlinear
            98, # TRIA6-nonlinear
            100, # CBAR-100
            228, # CQUADR
            232, # CQUADR-composite
            243, # CQUADX4
            189, # VUQUAD
            190, # VUTRIA
            191, # VUBEAM
            256, # CPYRAM
            227, # CTRIAR
            275, # CPLSTS3
        ]
        if self.element_type in elements_to_read:
            n = self._read_oes_4_sort(data, ndata)
        else:
            msg = self.code_information()
            n = self._not_implemented_or_skip(data, ndata, msg)
        return n

    def _read_oes1_4(self, data, ndata):
        """
        Reads the Stress Table 4
        """
        if self.table_name in NX_TABLES_BYTES:
            self.to_nx(f'because table_name={self.table_name}')
        #assert self.isubtable == -4, self.isubtable
        #if self.is_debug_file:
            #self.binary_debug.write('  element_name = %r\n' % self.element_name)
        #print "element_name =", self.element_name
        assert isinstance(self.format_code, int), self.format_code
        assert self.is_stress is True, self.code_information()
        self.data_code['is_stress_flag'] = True
        self.data_code['is_strain_flag'] = False

        self._setup_op2_subcase('STRESS/STRAIN')
        n = self._read_oes_4_sort(data, ndata)
        return n

    def _read_oes2_3(self, data, unused_ndata):
        """
        reads OES1 subtable 3
        """
        self._data_factor = 1
        self.words = [
            'aCode', 'tCode', 'element_type', 'isubcase',
            '???', '???', '???', 'load_set'
            'format_code', 'num_wide', 's_code', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', 'Title', 'subtitle', 'label', '???']

        self.parse_approach_code(data)  # 3
        self.sort_method = 2

        ## element type
        self.element_type = self.add_data_parameter(data, 'element_type', b'i', 3, False)

        ## load set ID
        self.load_set = self.add_data_parameter(data, 'load_set', b'i', 8, False)

        ## format code
        self.format_code = self.add_data_parameter(data, 'format_code', b'i', 9, False)

        ## number of words per entry in record
        ## .. note:: is this needed for this table ???
        self.num_wide = self.add_data_parameter(data, 'num_wide', b'i', 10, False)

        ## stress/strain codes
        self.s_code = self.add_data_parameter(data, 's_code', b'i', 11, False)

        ## thermal flag; 1 for heat ransfer, 0 otherwise
        self.thermal = self.add_data_parameter(data, 'thermal', b'i', 23, False)
        self.element_id = self.add_data_parameter(data, 'element_id', b'i', 5, fix_device_code=True)
        self._element_id = self.add_data_parameter(data, '_element_id', b'f', 5, apply_nonlinear_factor=False, add_to_dict=True)

        if self.analysis_code == 1:  # static...because reasons.
            self._analysis_code_fmt = b'i'
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
            self.apply_data_code_value('analysis_method', 'dt')
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

        self.fix_format_code()
        self._parse_stress_code_to_stress_bits()
        self._fix_oes_sort2(data)
        self._set_force_stress_element_name()
        #assert isinstance(self.nonlinear_factor, int), self.nonlinear_factor

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
        #self.format_code = format_code
        #self.data_code['format_code'] = format_code
        #self._check_result_type()


    def _fix_oes_sort2(self, data):
        self.fix_format_code()
        #print('self.format_code_original =', self.format_code_original)
        #print('self.format_code =', self.format_code)
            #self.fix_format_code()
            #if self.format_code == 1:
                #self.format_code = 2
                #self.data_code['format_code'] = 2
            #assert self.format_code in [2, 3], self.code_information()

        if 1:
            self.fix_format_code()
            #if self.num_wide == 8:
                #self.format_code = 1
                #self.data_code['format_code'] = 1
            #else:
                ##self.fix_format_code()
                #if self.format_code == 1:
                    #self.format_code = 2
                    #self.data_code['format_code'] = 2
                #assert self.format_code in [2, 3], self.code_information()

        self._parse_thermal_code()
        if self.is_debug_file:
            self.binary_debug.write('  %-14s = %r %s\n' % ('approach_code', self.approach_code,
                                                           self.approach_code_str(self.approach_code)))
            self.binary_debug.write('  %-14s = %r\n' % ('tCode', self.tCode))
            self.binary_debug.write('  %-14s = %r\n' % ('isubcase', self.isubcase))
        self._read_title(data)
        self._write_debug_bits()
        #assert isinstance(self.nonlinear_factor, int), self.nonlinear_factor

    def _read_ostr1_4(self, data, ndata):
        """
        Reads the Strain Table 4
        """
        if self.table_name in NX_TABLES_BYTES:
            self.to_nx(f' because table_name={self.table_name} was found')
        #assert self.isubtable == -4, self.isubtable
        #if self.is_debug_file:
            #self.binary_debug.write('  element_name = %r\n' % self.element_name)
        #print "element_name =", self.element_name
        assert self.is_strain is True, self.code_information()
        self.data_code['is_stress_flag'] = False
        self.data_code['is_strain_flag'] = True

        self._setup_op2_subcase('STRESS/STRAIN')
        n = self._read_ostr_4_sort(data, ndata)
        return n

    def _read_ostr2_4(self, data, ndata):
        """
        Reads the Strain Table 4
        """
        if self.table_name in NX_TABLES_BYTES:
            self.to_nx(f' because table_name={self.table_name} was found')
        #assert self.isubtable == -4, self.isubtable
        #if self.is_debug_file:
            #self.binary_debug.write('  element_name = %r\n' % self.element_name)
        #print("element_name =", self.element_name)
        assert self.is_strain is True, self.code_information()
        self.data_code['is_stress_flag'] = False
        self.data_code['is_strain_flag'] = True

        self._setup_op2_subcase('STRESS/STRAIN')
        if self.element_type in [1, 3, 10, # CROD, CTUBE, CTUBE
                                 11, 12, 13, # CELAS1, CELAS2, CELAS3,
                                 2, 4, 34, 33, 74, # CBEAM, CSHEAR, CBAR, CQUAD4, CTRIA3,
                                 75, 64, 70, 82, 144, # CTRIA6, CQUAD8, CTRIAR, CQUADR, CQUAD4
                                 69, # CBEND
                                 67, 68, 95, 102,#CHEXA, CPENTA, QUAD4-comp, CBUSH
                                 96, # QUAD8-nonlinear
                                 98, # TRIA6-nonlinear
                                 39, #CTETRA
                                 228, #CQUADR
                                 232, #CQUADR-composite
                                 233, #CTRIAR-composite
                                 97]:  # CTRIA3-C
            n = self._read_ostr_4_sort(data, ndata)
        else:
            msg = self.code_information()
            n = self._not_implemented_or_skip(data, ndata, msg)
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
            #except:
                #raise
                #print("----------")
                #print(self.obj)
                #print(self.data_code)
                #if self.obj is not None:
                    ##from pyNastran.utils import object_attributes
                    ##print object_attributes(self.obj)
                    #print(self.obj.data_code)
                #print("----------")
                #raise
            return n
        return new_func

    #@_print_obj_name_on_crash
    def _read_oes_4_sort(self, data, ndata):
        """Reads OES1 subtable 4 for NX/MSC/Autodesk/Optistruct"""
        #if self.num_wide == 146:
            #assert self.num_wide != 146
            #assert ndata != 146, self.code_information()
        assert isinstance(self.format_code, int), self.format_code
        if self.thermal == 0:
            n = self._read_oes1_loads(data, ndata)
        elif self.thermal == 1:
            n = self._read_oes1_thermal(data, ndata)
        else:
            msg = 'thermal=%s' % self.thermal
            n = self._not_implemented_or_skip(data, ndata, msg)
        return n

    #@_print_obj_name_on_crash
    def _read_ostr_4_sort(self, data, ndata):
        """
        Reads OSTR1 subtable 4
        """
        #if self.num_wide == 146:
            #assert self.num_wide != 146
            #assert ndata != 146, self.code_information()
        if self.thermal == 0:
            n = self._read_oes1_loads(data, ndata)
        elif self.thermal == 1:
            n = self._read_oes1_thermal(data, ndata)
        else:
            msg = 'thermal=%s' % self.thermal
            n = self._not_implemented_or_skip(data, ndata, msg)
        return n

    def _read_oes1_thermal(self, unused_data, ndata):
        """
        Reads OES self.thermal=1 tables; uses a hackish method to just skip the table
        """
        return ndata

    def _read_ostr1_thermal(self, unused_data, ndata):
        """
        Reads OSTR self.thermal=1 tables; uses a hackish method to just skip the table
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

            (11, 1, 2, b'OES1X1') : ('celas1_stress', RealSpringStressArray), # real
            (11, 2, 3, b'OES1X') : ('celas1_stress', ComplexSpringStressArray), # real/imag
            (11, 3, 3, b'OES1X') : ('celas1_stress', ComplexSpringStressArray), # mag/phase
            (11, 2, 3, b'OESVM1') : ('celas1_stress', ComplexSpringStressArray), # mag/phase
            (11, 3, 3, b'OESVM1') : ('celas1_stress', ComplexSpringStressArray), # mag/phase

            (12, 1, 2, b'OES1X1') : ('celas2_stress', RealSpringStressArray),
            (12, 1, 2, b'OES1X') : ('celas2_stress', RealSpringStressArray),
            (12, 1, 2, b'OES1') : ('celas2_stress', RealSpringStressArray),
            (12, 2, 3, b'OES1X') : ('celas2_stress', ComplexSpringStressArray),
            (12, 3, 3, b'OES1X') : ('celas2_stress', ComplexSpringStressArray),
            (12, 2, 3, b'OESVM1') : ('celas2_stress', ComplexSpringStressArray),
            (12, 3, 3, b'OESVM1') : ('celas2_stress', ComplexSpringStressArray),

            (13, 1, 2, b'OES1X1') : ('celas3_stress', RealSpringStressArray),
            #(13, 2, 3) : ('celas3_stress', ComplexSpringStressArray),
            #(13, 3, 3) : ('celas3_stress', ComplexSpringStressArray),
            (13, 2, 3, b'OES1X') : ('celas3_stress', ComplexSpringStressArray),
            (13, 2, 3, b'OESVM1') : ('celas3_stress', ComplexSpringStressArray),
            (13, 3, 3, b'OESVM1') : ('celas3_stress', ComplexSpringStressArray),

            (14, 1, 2) : ('celas4_stress', RealSpringStressArray),
            (14, 2, 3) : ('celas4_stress', ComplexSpringStressArray),
            (14, 3, 3) : ('celas4_stress', ComplexSpringStressArray),

            (34, 1, 16, b'OES1X1') : ('cbar_stress', RealBarStressArray),
            (34, 1, 16, b'OES1X') : ('cbar_stress', RealBarStressArray),
            (34, 1, 16, b'OES1') : ('cbar_stress', RealBarStressArray),
            (34, 2, 19, b'OES1X') : ('cbar_stress', ComplexBarStressArray),
            (34, 1, 10, b'OESNO1') : ('cbar_stress', ComplexBarStressArray),
            (34, 2, 10, b'OESXRMS1') : ('cbar_stress', ComplexBarStressArray),

            (34, 1, 10, b'OESRMS2') : ('cbar_stress', RandomBarStressArray),

            (34, 2, 10, b'OESPSD2') : ('cbar_stress', RandomBarStressArray),
            (34, 2, 10, b'OESRMS2') : ('cbar_stress', RandomBarStressArray),
            (34, 2, 10, b'OESNO2') : ('cbar_stress', RandomBarStressArray),
            (34, 2, 10, b'OESATO2') : ('cbar_stress', RandomBarStressArray),
            (34, 2, 10, b'OESCRM2') : ('cbar_stress', RandomBarStressArray),

            # Missing stress_mapper key for OES1 table #501
            # see cbarao_random_x_mini.op2 for an example with OES1 and OES1X...
            # it looks to be an error in MSC [2008-2012)
            (34, 2, 19, b'OES1') : ('cbar_stress', ComplexBarStressArray),
            (34, 3, 19, b'OES1X') : ('cbar_stress', ComplexBarStressArray),
            (34, 3, 19, b'OESVM1') : ('cbar_stress', ComplexBarStressArray),
            #(34, 1, 19) : ('cbar_stress', RandomBarStressArray),
            (100, 1, 10, b'OES1X1') : ('cbar_stress_10nodes', RealBar10NodesStressArray),
            (100, 1, 10, b'OES1X') : ('cbar_stress_10nodes', RealBar10NodesStressArray),

            # solids
            (39, 1, 109, b'OES1X1') : ('ctetra_stress', RealSolidStressArray), # real
            (39, 1, 109, b'OES1X') : ('ctetra_stress', RealSolidStressArray), # real
            (39, 1, 109, b'OES1') : ('ctetra_stress', RealSolidStressArray), # real
            (39, 3, 74, b'OESVM1') : ('ctetra_stress', ComplexSolidStressArray), # mag/phase

            (67, 1, 193, b'OES1X1') : ('chexa_stress', RealSolidStressArray),
            (67, 1, 193, b'OES1X') : ('chexa_stress', RealSolidStressArray),
            (67, 1, 193, b'OES1') : ('chexa_stress', RealSolidStressArray),
            (67, 1, 193, b'RASCONS') : ('chexa_stress', RealSolidStressArray),

            (68, 1, 151, b'OES1X1') : ('cpenta_stress', RealSolidStressArray),
            (68, 1, 151, b'OES1X') : ('cpenta_stress', RealSolidStressArray),
            (68, 1, 151, b'OES1') : ('cpenta_stress', RealSolidStressArray),
            (68, 3, 102, b'OESVM1') : ('cpenta_stress', ComplexSolidStressArray),

            (39, 2, 69, b'OES1X') : ('ctetra_stress', ComplexSolidStressArray), # real/imag
            (39, 2, 69, b'OES1') : ('ctetra_stress', ComplexSolidStressArray),
            (39, 2, 74, b'OESVM1') : ('ctetra_stress', 'NA'), # real/imag
            #(39, 3, 69) : ('ctetra_stress', ComplexSolidStressArray), # mag/phase

            (67, 2, 121, b'OES1X') : ('chexa_stress', ComplexSolidStressArray),
            (67, 3, 121, b'OES1X') : ('chexa_stress', ComplexSolidStressArray),
            (67, 3, 130, b'OESVM1') : ('chexa_stress', ComplexSolidStressArray),
            (67, 2, 121, b'OES1') : ('chexa_stress', ComplexSolidStressArray),
            (67, 3, 121, b'OES1') : ('chexa_stress', ComplexSolidStressArray),

            (68, 2, 95, b'OES1X') : ('cpenta_stress', ComplexSolidStressArray),
            (68, 3, 95, b'OES1X') : ('cpenta_stress', ComplexSolidStressArray),
            (68, 2, 95, b'OES1') : ('cpenta_stress', ComplexSolidStressArray),

            (33, 1, 17, b'OES1X1') :  ('cquad4_stress', RealPlateStressArray),
            (33, 1, 17, b'OES1X') :  ('cquad4_stress', RealPlateStressArray),
            (33, 1, 17, b'OES1') :  ('cquad4_stress', RealPlateStressArray),
            (33, 2, 15, b'OES1X') :  ('cquad4_stress', ComplexPlateStressArray),
            (33, 3, 15, b'OES1X') :  ('cquad4_stress', ComplexPlateStressArray),
            #(33, 3, 0) :  ('cquad4_stress', RandomPlateStressArray),
            (33, 1, 9, b'OESNO1') : ('cquad4_stress', ComplexPlateStressArray),
            (33, 2, 11, b'OESXRMS1') : ('cquad4_stress', ComplexPlateStressArray),

            (33, 2, 9, b'OESATO2') : ('cquad4_stress', 'NA'),
            (33, 2, 9, b'OESCRM2') : ('cquad4_stress', 'NA'),
            (33, 2, 9, b'OESPSD2') : ('cquad4_stress', 'NA'),
            (33, 2, 9, b'OESNO2') : ('cquad4_stress', 'NA'),

            (33, 1, 9, b'OESRMS2') : ('cquad4_stress', 'NA'),
            (33, 2, 9, b'OESRMS2') : ('cquad4_stress', 'NA'),


            (74, 1, 17, b'OES1X1') : ('ctria3_stress', RealPlateStrainArray),
            (74, 1, 17, b'OES1X') : ('ctria3_stress', RealPlateStrainArray),
            (74, 1, 17, b'OES1') : ('ctria3_stress', RealPlateStrainArray),
            (74, 2, 15, b'OES1X') : ('ctria3_stress', ComplexPlateStrainArray),
            (74, 3, 15, b'OES1X') : ('ctria3_stress', ComplexPlateStrainArray),
            (74, 2, 11, b'OESXRMS1') : ('ctria3_stress', ComplexPlateStrainArray),
            (74, 1, 9, b'OESNO1') : ('ctria3_stress', ComplexPlateStrainArray),
            (74, 2, 17, b'OESVM1') : ('ctria3_stress', 'NA'),
            (74, 3, 17, b'OESVM1') : ('ctria3_stress', 'NA'),

            (74, 1, 9, b'OESRMS2') : ('ctria3_stress', 'NA'),

            #(74, 1, 9) : ('ctria3_stress', RandomPlateStressArray),

            (82, 1, 87, b'OES1X1') : ('cquadr_stress', RealPlateStressArray),
            (82, 1, 87, b'OES1X') : ('cquadr_stress', RealPlateStressArray),
            (82, 2, 77, b'OES1X') : ('cquadr_stress', ComplexPlateStressArray),
            (82, 3, 77, b'OES1X') : ('cquadr_stress', ComplexPlateStressArray),

            (64, 1, 87, b'OES1X1') : ('cquad8_stress', RealPlateStressArray), # real
            (64, 1, 87, b'OES1X') : ('cquad8_stress', RealPlateStressArray),
            (64, 1, 87, b'OES1') : ('cquad8_stress', RealPlateStressArray),
            (64, 2, 77, b'OES1') : ('cquad8_stress', ComplexPlateStressArray), # real/imag
            (64, 3, 77, b'OES1') : ('cquad8_stress', ComplexPlateStressArray), # mag/phase
            (64, 2, 77, b'OES1X') : ('cquad8_stress', ComplexPlateStressArray), # real/imag
            (64, 3, 77, b'OES1X') : ('cquad8_stress', ComplexPlateStressArray), # mag/phase
            (64, 2, 87, b'OESVM1') : ('cquad8_stress', ComplexPlateStressArray), # real/imag
            (64, 3, 87, b'OESVM1') : ('cquad8_stress', ComplexPlateStressArray), # mag/phase

            (70, 1, 70, b'OES1X1') : ('ctriar_stress', RealPlateStressArray),
            (70, 1, 70, b'OES1X') : ('ctriar_stress', RealPlateStressArray),
            (70, 2, 62, b'OES1X') : ('ctriar_stress', ComplexPlateStressArray),
            (70, 3, 62, b'OES1X') : ('ctriar_stress', ComplexPlateStressArray),

            (75, 1, 70, b'OES1X1') : ('ctria6_stress', RealPlateStressArray),
            (75, 2, 62, b'OES1X') : ('ctria6_stress', ComplexPlateStressArray),
            (75, 3, 62, b'OES1X') : ('ctria6_stress', ComplexPlateStressArray),
            (75, 2, 70, b'OESVM1') : ('ctria6_stress', ComplexPlateStressArray),
            (75, 3, 70, b'OESVM1') : ('ctria6_stress', ComplexPlateStressArray),

            (144, 1, 87, b'OES1X1') : ('cquad4_stress', RealPlateStressArray),
            (144, 1, 87, b'OES1') : ('cquad4_stress', RealPlateStressArray),
            (144, 1, 87, b'RASCONS') : ('cquad4_stress', RealPlateStressArray),

            (144, 2, 77, b'OES1X') : ('cquad4_stress', ComplexPlateStressArray),
            (144, 3, 77, b'OES1X') : ('cquad4_stress', ComplexPlateStressArray),
            (144, 3, 87, b'OESVM1') : ('cquad4_stress', ComplexPlateStressArray),
            #(144, 3, 77) : ('cquad4_stress', ComplexPlateStressArray),
            #(64, 1, 47) : ('cquad8_stress', RandomPlateStressArray), # random
            #(70, 1, 39) : ('ctriar_stress', RandomPlateStressArray),
            #(75, 1, 39) : ('ctria6_stress', RandomPlateStressArray),
            #(82, 1, 47) : ('cquadr_stress', RandomPlateStressArray),
            #(144, 1, 47) : ('cquad4_stress', RandomPlateStressArray),

            (88, 1, 13, b'OESNLXR') : ('nonlinear_ctria3_stress', RealNonlinearPlateArray), # real
            (88, 1, 25, b'OESNL1X') : ('nonlinear_ctria3_stress', RealNonlinearPlateArray), # real?
            (88, 1, 25, b'OESNLXR') : ('nonlinear_ctria3_stress', RealNonlinearPlateArray), # real?

            (90, 1, 13, b'OESNLXR') : ('nonlinear_cquad4_stress', RealNonlinearPlateArray),
            (90, 1, 25, b'OESNL1X') : ('nonlinear_cquad4_stress', RealNonlinearPlateArray),
            (90, 1, 25, b'OESNLXR') : ('nonlinear_cquad4_stress', RealNonlinearPlateArray),
            (90, 1, 25, b'OESNLXD') : ('nonlinear_cquad4_stress', RealNonlinearPlateArray),

            (95, 1, 11, b'OES1C') : ('cquad4_composite_stress', RealCompositePlateStressArray), # real
            (95, 1, 11, b'OESCP') : ('cquad4_composite_stress', RealCompositePlateStressArray), # real
            (95, 1, 9, b'OESRT') : ('cquad4_composite_stress', 'RandomCompositePlateStressArray'), # real
            (95, 2, 11, b'OESCP') : ('cquad4_composite_stress', RealCompositePlateStressArray), # real?
            (95, 2, 11, b'OESRT') : ('cquad4_composite_stress', RealCompositePlateStressArray), # real?
            #(95, 2, 9) : ('cquad4_composite_stress', ComplexCompositePlateStressArray), # real/imag
            #(95, 3, 9) : ('cquad4_composite_stress', ComplexCompositePlateStressArray), # mag/phase

            #(96, 1, 9) : ('cquad8_composite_stress', 'RandomCompositePlateStressArray'),
            (96, 1, 11, b'OES1C') : ('cquad8_composite_stress', RealCompositePlateStressArray),
            #(96, 1, 11) : ('cquad8_composite_stress', RealCompositePlateStressArray),
            #(96, 2, 9) : ('cquad8_composite_stress', ComplexCompositePlateStressArray),
            #(96, 3, 9) : ('cquad8_composite_stress', ComplexCompositePlateStressArray),

            (97, 1, 9, b'OESRT') : ('ctria3_composite_stress', 'RandomCompositePlateStressArray'),
            (97, 1, 11, b'OES1C') : ('ctria3_composite_stress', RealCompositePlateStressArray),
            (97, 1, 11, b'OESCP') : ('ctria3_composite_stress', RealCompositePlateStressArray),
            (97, 2, 11, b'OESCP') : ('ctria3_composite_stress', RealCompositePlateStressArray),
            #(97, 2, 9) : ('ctria3_composite_stress', ComplexCompositePlateStressArray),
            #(97, 3, 9) : ('ctria3_composite_stress', ComplexCompositePlateStressArray),

            (98, 1, 9, b'OESRT') : ('ctria6_composite_stress', 'RandomCompositePlateStressArray'),
            (98, 1, 11, b'OES1C') : ('ctria6_composite_stress', RealCompositePlateStressArray),
            #(98, 1, 11) : ('ctria6_composite_stress', RealCompositePlateStressArray),
            #(98, 2, 9) : ('ctria6_composite_stress', ComplexCompositePlateStressArray),
            #(98, 3, 9) : ('ctria6_composite_stress', ComplexCompositePlateStressArray),

            (53, 1, 33, b'OES1X1') : ('ctriax_stress', RealTriaxStressArray),
            (53, 1, 33, b'OES1X') : ('ctriax_stress', RealTriaxStressArray),
            (53, 2, 37, b'OES1X') : ('ctriax_stress', ComplexTriaxStressArray),
            #(53, 3, 37) : ('ctriax_stress', ComplexTriaxStressArray),

            (102, 1, 7, b'OES1X1') : ('cbush_stress', RealBushStressArray),
            (102, 1, 7, b'OES1X') : ('cbush_stress', RealBushStressArray),
            (102, 1, 7, b'OES1') : ('cbush_stress', RealBushStressArray),
            (102, 2, 13, b'OES1X') : ('cbush_stress', ComplexCBushStressArray),
            (102, 3, 13, b'OES1X') : ('cbush_stress', ComplexCBushStressArray),
            (102, 2, 13, b'OESVM1') : ('cbush_stress', 'NA'),
            (102, 2, 13, b'OES1'): ('cbush_stress', ComplexCBushStressArray),

            (40, 1, 8, b'OES1X1') : ('cbush1d_stress_strain', RealBushStressArray),
            (40, 1, 8, b'OESNLXD') : ('cbush1d_stress_strain', RealBushStressArray),
            #(40, 2, 9) : ('cbush1d_stress_strain', ComplexCBushStressArray),
            #(40, 3, 9) : ('cbush1d_stress_strain', ComplexCBushStressArray),

            (87, 1, 7, b'OESNL1X') : ('nonlinear_ctube_stress', RealNonlinearRodArray),
            (87, 1, 7, b'OESNLXR') : ('nonlinear_ctube_stress', RealNonlinearRodArray),
            (89, 1, 7, b'OESNL1X') : ('nonlinear_crod_stress', RealNonlinearRodArray),
            (89, 1, 7, b'OESNLXD') : ('nonlinear_crod_stress', RealNonlinearRodArray),
            (89, 1, 7, b'OESNLXR') : ('nonlinear_crod_stress', RealNonlinearRodArray),
            (92, 1, 7, b'OESNL1X') : ('nonlinear_conrod_stress', RealNonlinearRodArray),
            (92, 1, 7, b'OESNLXD') : ('nonlinear_conrod_stress', RealNonlinearRodArray),
            (92, 1, 7, b'OESNLXR') : ('nonlinear_conrod_stress', RealNonlinearRodArray),

            (224, 1, 3, b'OESNLXD') : ('nonlinear_celas1_stress', RealNonlinearSpringStressArray),
            (224, 1, 3, b'OESNLXR') : ('nonlinear_celas1_stress', RealNonlinearSpringStressArray),
            (225, 1, 3, b'OESNLXR') : ('nonlinear_celas3_stress', RealNonlinearSpringStressArray),

            (35, 1, 18, b'OES1X1') : ('NA', 'NA'), # CCONEAX
            (35, 1, 18, b'OES1') : ('NA', 'NA'), # CCONEAX

            (60, 1, 10, b'OES1X') : ('NA', 'NA'), # DUM8/CCRAC2D
            (61, 1, 10, b'OES1X') : ('NA', 'NA'), # DUM8/CCRAC3D

            (69, 1, 21, b'OES1X1') : ('NA', 'NA'), # CBEND
            (69, 2, 21, b'OES1X') : ('NA', 'NA'), # CBEND
            (69, 3, 21, b'OES1X') : ('NA', 'NA'), # CBEND

            (86, 1, 11, b'OESNL1X') : ('nonlinear_cgap_stress', NonlinearGapStressArray),
            (86, 1, 11, b'OESNLXR') : ('nonlinear_cgap_stress', NonlinearGapStressArray),
            (86, 1, 11, b'OESNLXD') : ('nonlinear_cgap_stress', NonlinearGapStressArray),
            (94, 1, 51, b'OESNL1X') : ('nonlinear_cbeam_stress', RealNonlinearBeamStressArray),
            (94, 1, 51, b'OESNLXR') : ('nonlinear_cbeam_stress', RealNonlinearBeamStressArray),

            (85, 1, 82, b'OESNLXR') : ('NA', 'NA'),  # TETRANL
            (91, 1, 114, b'OESNLXD') : ('NA', 'NA'),  # PENTANL
            (91, 1, 114, b'OESNLXR') : ('NA', 'NA'),  # PENTANL
            (93, 1, 146, b'OESNL1X') : ('NA', 'NA'),  # HEXANL
            (93, 1, 146, b'OESNLXD') : ('NA', 'NA'),  # HEXANL
            (93, 1, 146, b'OESNLXR') : ('NA', 'NA'),  # HEXANL

            # 101-AABSF
            (101, 2, 4, b'OES1X') : ('NA', 'NA'),

            # 140-HEXA8FD
            (140, 1, 162, b'OES1X1') : ('NA', 'NA'),
            #201-QUAD4FD
            (201, 1, 46, b'OESNLXD') : ('NA', 'NA'),
            (201, 1, 46, b'OESNLXR') : ('NA', 'NA'),

            # 145-VUHEXA  (8 nodes)
            (145, 1, 98, b'OES1X1') : ('NA', 'NA'),
            (145, 2, 106, b'OES1X') : ('NA', 'NA'),
            (145, 3, 106, b'OES1X') : ('NA', 'NA'),
            # 146-VUPENTA (6 nodes)
            (146, 1, 74, b'OES1X1') : ('NA', 'NA'),
            (146, 2, 80, b'OES1X') : ('NA', 'NA'),
            (146, 3, 80, b'OES1X') : ('NA', 'NA'),
            # 147-VUTETRA (4 nodes)
            (147, 1, 50, b'OES1X1') : ('NA', 'NA'),
            (147, 2, 54, b'OES1X') : ('NA', 'NA'),
            (147, 3, 54, b'OES1X') : ('NA', 'NA'),

            # 139-QUAD4FD
            # self.hyperelastic_cquad4_strain, HyperelasticQuad
            (139, 1, 30, b'OES1X1') : ('NA', 'NA'),

            # 189-VUQUAD
            (189, 1, 74, b'OES1X1') : ('NA', 'NA'),
            (189, 2, 114, b'OES1X') : ('NA', 'NA'),

            # 47-AXIF2
            (47, 2, 9, b'OES1X') : ('axif2', 'NA'),
            # 48-AXIF3
            (48, 2, 19, b'OES1X') : ('axif3', 'NA'),
            # 190-VUTRIA
            (190, 1, 57, b'OES1X1') : ('NA', 'NA'),
            (190, 2, 87, b'OES1X') : ('NA', 'NA'),
            (190, 3, 87, b'OES1X') : ('NA', 'NA'),

            # 191-VUBEAM
            (191, 1, 60, b'OES1X1') : ('vubeam', 'NA'),
            (191, 2, 80, b'OES1X') : ('vubeam', 'NA'),
            (191, 3, 80, b'OES1X') : ('vubeam', 'NA'),

            # 203-SLIF1D?
            (203, 1, 14, b'OESNLBR') : ('slif1d', 'NA'),
            # 50-SLOT3
            (50, 2, 11, b'OES1X') : ('slot3', 'NA'),
            # 51-SLOT4
            (51, 2, 13, b'OES1X') : ('slot4', 'NA'),

            # 160-PENTA6FD
            (160, 1, 122, b'OES1X1') : ('cpenta', 'NA'),
            # 161-TETRA4FD
            (161, 1, 22, b'OES1X1') : ('ctetra', 'NA'),
            # 162-TRIA3FD
            (162, 1, 9, b'OES1X1') : ('ctria', 'NA'),
            # 163-HEXAFD
            (163, 1, 542, b'OES1X1') : ('chexa', 'NA'),
            # 164-QUADFD
            (164, 1, 65, b'OES1X1') : ('cquad', 'NA'),
            # 165-PENTAFD
            (165, 1, 422, b'OES1X1') : ('cpenta', 'NA'),
            # 166-TETRAFD
            (166, 1, 102, b'OES1X1') : ('ctetra', 'NA'),
            # 167-TRIAFD
            (167, 1, 23, b'OES1X1') : ('NA', 'NA'),
            # 168-TRIAX3FD
            (168, 1, 9, b'OES1X1') : ('ctriax3', 'NA'),
            # 169-TRIAXFD
            (169, 1, 23, b'OES1X1') : ('ctriax', 'NA'),
            # 170-QUADX4FD
            (170, 1, 30, b'OES1X1') : ('cquadx4fd', 'NA'),
            # 171-QUADXFD
            (171, 1, 65, b'OES1X1') : ('cquadx', 'NA'),
            # 172-QUADRNL
            (172, 1, 25, b'OESNLXR') : ('cquadrnl', 'NA'),
            # 202-HEXA8FD
            (202, 1, 122, b'OESNLXD') : ('chexa', 'NA'),
            (202, 1, 122, b'OESNLXR') : ('chexa', 'NA'),
            # 204-PENTA6FD
            (204, 1, 92, b'OESNLXR') : ('cpenta', 'NA'),
            # 211-TRIAFD
            (211, 1, 35, b'OESNLXR') : ('ctria3', 'NA'),
            # 213-TRIAXFD
            (213, 1, 35, b'OESNLXR') : ('ctriax', 'NA'),
            # 214-QUADX4FD
            (214, 1, 46, b'OESNLXR') : ('cquadx4', 'NA'),
            # 216-TETRA4FD
            (216, 1, 62, b'OESNLXD') : ('NA', 'NA'),
            (216, 1, 62, b'OESNLXR') : ('NA', 'NA'),
            # 217-TRIA3FD
            (217, 1, 35, b'OESNLXR') : ('ctria3', 'NA'),
            # 218-HEXAFD
            (218, 1, 122, b'OESNLXR') : ('chexa', 'NA'),
            # 219-QUADFD
            (219, 1, 46, b'OESNLXR') : ('cquad', 'NA'),
            # 220-PENTAFD
            (220, 1, 92, b'OESNLXR') : ('cpenta', 'NA'),
            # 221-TETRAFD
            (221, 1, 62, b'OESNLXR') : ('tetrafd', 'NA'),
            # 222-TRIAX3FD
            (222, 1, 35, b'OESNLXR') : ('ctriax3fd', 'NA'),
            # 223-CQUADXFD
            (223, 1, 46, b'OESNLXR') : ('cquadx', 'NA'),
            # 226-BUSH
            (226, 1, 19, b'OESNLXD') : ('cbush', 'NA'),
            (226, 1, 19, b'OESNLXR') : ('cbush', 'NA'),
            # 227-CTRIAR
            (227, 1, 17, b'OES1X1') : ('ctriar', 'NA'),
            (227, 1, 17, b'OES1X') : ('ctriar', 'NA'),
            # 228-CQUADR
            (228, 1, 17, b'OES1X1') : ('cquadr', 'NA'),
            (228, 1, 17, b'OES1X') : ('cquadr', 'NA'),
            # 232-QUADRLC
            (232, 1, 11, b'OES1C') : ('cquadr', 'NA'),
            (232, 1, 11, b'OESCP') : ('cquadr', 'NA'),
            (232, 2, 13, b'OESVM1C') : ('cquadr', 'NA'),  # freq nx
            (232, 3, 13, b'OESVM1C') : ('cquadr', 'NA'),  # freq nx
            #(234, 1, 11) : ('cquadr', 'NA'), # bad?
            # 233-TRIARLC
            (233, 1, 11, b'OES1C') : ('ctriar', 'NA'),
            (233, 2, 13, b'OESVM1C') : ('ctriar', 'NA'),  # freq nx
            (233, 3, 13, b'OESVM1C') : ('ctriar', 'NA'),  # freq nx
            # 235-CQUADR
            (235, 1, 17, b'OES1X1') : ('NA', 'NA'),
            (235, 2, 15, b'OES1X') : ('NA', 'NA'),

            # 242-CTRAX
            # 244-CTRAX6
            (242, 1, 34, b'OES1X1') : ('ctrax', 'NA'),
            (244, 1, 34, b'OES1X1') : ('ctrax6', 'NA'),

            # 243-CQUADX4
            # 245-CQUADX8
            (243, 1, 42, b'OES1X1') : ('cquadx4', 'NA'),
            (245, 1, 42, b'OES1X1') : ('cquadx8', 'NA'),

            #256-CPYRAM
            (255, 1, 130, b'OES1X1') : ('cpyram', 'NA'),
            (255, 2, 82, b'OES1X') : ('cpyram', 'NA'),
            (256, 1, 98, b'OESNLXD') : ('cpyram', 'NA'),

            # 271-CPLSTN3
            # 272-CPLSTN4
            (271, 1, 6, b'OES1X1') : ('cplstn3', 'NA'),
            (271, 1, 6, b'OES1X') : ('cplstn3', 'NA'),
            (272, 1, 32, b'OES1X1') : ('cplstn4', 'NA'),
            (272, 1, 32, b'OES1X') : ('cplstn4', 'NA'),
            (273, 1, 26, b'OES1X1') : ('cplstn6', 'NA'),
            (273, 1, 26, b'OES1X') : ('cplstn6', 'NA'),
            (274, 1, 32, b'OES1X1') : ('cplstn3', 'NA'),
            (274, 1, 32, b'OES1X') : ('cplstn3', 'NA'),

            # 275-CPLSTS3
            # 277-CPLSTS6
            (275, 1, 6, b'OES1X1') : ('cplsts3', 'NA'),
            (276, 1, 32, b'OES1X1') : ('cplsts4', 'NA'),
            (277, 1, 26, b'OES1X1') : ('cplsts6', 'NA'),
            (278, 1, 32, b'OES1X1') : ('cplsts8', 'NA'),

            (1, 2, 5, 'OESVM1') : ('crod', 'NA'),
            (10, 2, 5, 'OESVM1') : ('conrod', 'NA'),
            (10, 2, 5, 'OES1X') : ('conrod', 'NA'),

            (11, 2, 3, 'OESVM1') : ('celas1', 'NA'),
            (12, 2, 3, 'OESVM1') : ('celas2', 'NA'),

            (2, 2, 111, b'OESVM1') : ('cbeam', 'NA'),
            (34, 2, 19, b'OESVM1') : ('cbar', 'NA'),

            (4, 2, 5, 'OESVM1') : ('cshear', 'NA'),
            (4, 2, 5, 'OES1X') : ('cshear', 'NA'),
            (74, 2, 17, 'OESVM1') : ('ctria3', 'NA'),
            (144, 2, 87, b'OESVM1') : ('cquad4', 'NA'),

            (95, 2, 13, b'OESVM1C') : ('cquad4', 'NA'),
            (95, 3, 13, b'OESVM1C') : ('cquad4', 'NA'),
            (97, 2, 13, b'OESVM1C') : ('ctria3', 'NA'),
            (97, 3, 13, b'OESVM1C') : ('ctria3', 'NA'),

            (39, 2, 74, 'OESVM1') : ('ctetra', 'NA'),
            (67, 2, 130, b'OESVM1') : ('chexa', 'NA'),
            (68, 2, 102, b'OESVM1') : ('cpenta', 'NA'),
        }

        key = (self.element_type, self.format_code, self.num_wide, self.table_name)
        try:
            return stress_mapper[key]
        except KeyError:  # pragma: no cover
            self.log.error(self.code_information())
            msg = ('stress_mapper (~line 850 of oes.py) does not contain the '
                   'following key and must be added\n'
                   'key=(element_type=%r, format_code=%r, num_wide=%r, table_name=%r) ' % key)
            self.log.error(msg)
            #raise KeyError(msg)
            raise
            #return None, None

    def get_oes_prefix_postfix(self):
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
        prefix = ''
        postfix = ''

        table_name_bytes = self.table_name
        assert isinstance(table_name_bytes, bytes), table_name_bytes
        is_sort1 = table_name_bytes in SORT1_TABLES_BYTES

        if table_name_bytes in [b'OES1X1', b'OES1X', b'OSTR1X', b'OSTR1',
                               b'OES1C', b'OSTR1C', b'OES1', ]:
            self._set_as_sort1()
        elif table_name_bytes in [b'OES2', b'OSTR2', b'OES2C', b'OSTR2C']:
            self._set_as_sort2()
        #elif table_name_bytes in ['OESNLXR']:
            #prefix = 'sideline_'
        elif table_name_bytes in [b'OESNLXD', b'OESNL1X', b'OESNLXR', b'OESNL2']:
            prefix = 'nonlinear_'
        elif table_name_bytes in [b'OESNLXR2']:
            prefix = 'nonlinear_'
        elif table_name_bytes == b'OESNLBR':
            prefix = 'sideline_'
        elif table_name_bytes == b'OESRT':
            prefix = 'strength_ratio.'
        elif table_name_bytes in [b'OESCP', b'OESTRCP']:
            # guessing
            pass
            #self.sort_bits[0] = 0 # real; ???
            #self.sort_bits[1] = 0 # sort1
            #self.sort_bits[2] = 1 # random; ???
        elif table_name_bytes in [b'OESVM1C', b'OSTRVM1C', b'OESVM1', b'OSTRVM1',
                                  #b'OESVM1C', b'OSTRVM1C',
                                  b'OESVM2', b'OSTRVM2',]:
            prefix = 'modal_contribution.'
            self.to_nx(f' because table_name={table_name_bytes} was found')

        #----------------------------------------------------------------
        elif table_name_bytes in [b'OSTRMS1C']: #, b'OSTRMS1C']:
            self.format_code = 1
            self.sort_bits[0] = 0 # real
            prefix = 'rms.'

        elif table_name_bytes in [b'OESXRMS1']:
            self._analysis_code_fmt = b'i'
            self._set_as_random()
            self._set_as_sort1()
            prefix = 'rms.'
        elif table_name_bytes in [b'OESXRMS2']: # wrong?
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
            #print(self.code_information())

        elif table_name_bytes in [b'OESRMS1', b'OSTRRMS1']:
            self._analysis_code_fmt = b'i'
            self._set_as_random()
            self._set_as_sort1()
            prefix = 'rms.'
        elif table_name_bytes in [b'OESRMS2', b'OSTRRMS2']:
            self._analysis_code_fmt = b'i'
            self._set_as_random()
            self._set_as_sort1()  # it's not really SORT2...
            self.sort_method = 1
            if table_name_bytes == b'OESRMS2':
                self.table_name = b'OESRMS1'
            elif table_name_bytes == b'OSTRRMS2':
                self.table_name = b'OSTRRMS1'
            else:
                raise NotImplementedError(table_name_bytes)
            #assert self.sort_method == 2, self.code_information()
            prefix = 'rms.'

        elif table_name_bytes in [b'OESNO1', b'OSTRNO1', b'OSTNO1C']:
            assert self.sort_method == 1, self.code_information()
            self._set_as_random()
            prefix = 'no.'
        elif table_name_bytes in [b'OESNO2', b'OSTRNO2']:
            self._set_as_random()
            self._set_as_sort1()
            self.data_code['nonlinear_factor'] = None
            self._analysis_code_fmt = b'i'
            prefix = 'no.'
        #----------------------------------------------------------------

        elif table_name_bytes in [b'OESPSD1', b'OSTRPSD1']:
            #self.format_code = 1
            self.sort_bits[0] = 0 # real
            self.sort_bits[1] = 0 # sort1
            self.sort_bits[2] = 1 # random
            prefix = 'psd.'
        elif table_name_bytes in [b'OESPSD2', b'OSTRPSD2',
                                  b'OESPSD2C', b'OSTPSD2C']:
            if 0:
                # TODO: the sort bits might not be right...isat_random
                #print(self.code_information())
                #print(self.sort_bits)
                self.format_code = 1
                #self.sort_bits[0] = 0 # real
                #self.sort_bits[1] = 1 # sort2
                #self.sort_bits[2] = 1 # random
                self.sort_bits.is_real = 1
                self.sort_bits.is_sort2 = 1
                self.sort_bits.is_random = 1
                #print(self.code_information())
                #print(self.sort_bits)
            else:
                self.format_code = 1 # real
                self.result_type = 2 # random
                self.sort_bits[0] = 0 # real
                self.sort_bits[2] = 1 # random
            prefix = 'psd.'

        elif table_name_bytes in [b'OESATO1', b'OSTRATO1']:
            prefix = 'ato.'
        elif table_name_bytes in [b'OESATO2', b'OSTRATO2']:
            prefix = 'ato.'

        elif table_name_bytes in [b'OESCRM1', b'OSTRCRM1']:
            prefix = 'crm.'
            self.result_type = 2 # random
            self.sort_bits[2] = 1 # random
        elif table_name_bytes in [b'OESCRM2', b'OSTRCRM2']:
            # sort2, random
            self.format_code = 1 # real
            self.result_type = 2 # random
            self.sort_bits[0] = 0 # real
            self.sort_bits[1] = 1 # sort2
            self.sort_bits[2] = 1 # random
            self.sort_method = 2
            prefix = 'crm.'
        #elif self.table_name in ['DOES1', 'DOSTR1']:
            #prefix = 'scaled_response_spectra_'
        #elif self.table_name in ['OESCP']:

        elif table_name_bytes in [b'RASCONS']: #, b'OSTRMS1C']:
            self.format_code = 1
            self.sort_bits[0] = 0 # real
            prefix = 'RASCONS.'
        elif table_name_bytes in [b'RAECONS']: #, b'OSTRMS1C']:
            self.format_code = 1
            self.sort_bits[0] = 0 # real
            prefix = 'RAECONS.'
        elif table_name_bytes in [b'RAPCONS']: #, b'OSTRMS1C']:
            self.format_code = 1
            self.sort_bits[0] = 0 # real
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
        else:
            raise NotImplementedError(self.table_name)

        #if self.analysis_code == 1:
            #self.sort_bits[1] = 0 # sort1
            #self.sort_method = 1

        self.data_code['sort_bits'] = self.sort_bits
        self.data_code['nonlinear_factor'] = self.nonlinear_factor
        return prefix, postfix

    def _set_as_real(self):
        self.format_code = 1
        self.result_type = 0
        self.sort_bits[0] = 0 # real
        self.sort_bits.is_real = True
        self.sort_bits.is_random = False

    def _set_as_random(self):
        self.format_code = 1  # real
        self.result_type = 2  # random
        self.sort_bits.is_real = True
        self.sort_bits.is_random = True

    def _set_as_sort1(self):
        self.sort_bits[1] = 0 # sort1
        self.sort_method = 1

    def _set_as_sort2(self):
        self.sort_bits[1] = 1 # sort2
        self.sort_method = 2

    def _read_oesmc_4(self, data, ndata):
        n = 0
        if self.element_type == 1:
            assert self.num_wide == 4, self.code_information()
            if self.read_mode == 1:
                return ndata
            ntotal = 16 * self.factor # 4*4
            nelements = ndata // ntotal
            fmt = mapfmt(self._endian + b'i3f', self.size)
            struct1 = Struct(fmt)
            for ielem in range(nelements):
                edata = data[n:n+ntotal]
                out = struct1.unpack(edata)
                #print(out)
                n += ntotal
            self.log.warning(f'skipping {self.table_name} with {self.element_name}-{self.element_type}')
        else:
            raise NotImplementedError(self.code_information())
        return n

    def _read_oes1_loads_nasa95(self, data, ndata: int) -> Tuple[int, Any, Any]:
        """Reads OES1 subtable 4 for NASA 95"""
        prefix, postfix = self.get_oes_prefix_postfix()
        result_type = self.result_type
        #self._apply_oes_ato_crm_psd_rms_no('') # TODO: just testing
        n = 0
        is_magnitude_phase = self.is_magnitude_phase()
        dt = self.nonlinear_factor

        if self.is_stress:
            result_name = 'stress'
            #stress_name = 'STRESS'
        else:
            result_name = 'strain'
            #stress_name = 'STRAIN'

        if self._results.is_not_saved(result_name):
            return ndata
        if self.element_type in [1, 3, 10]:  # rods
            # 1-CROD
            # 3-CTUBE
            # 10-CONROD
            n, nelements, ntotal = self._oes_crod(data, ndata, dt, is_magnitude_phase,
                                                  result_type, prefix, postfix)

        #elif self.element_type == 2: # CBEAM
            #n, nelements, ntotal = self._oes_cbeam(data, ndata, dt, is_magnitude_phase,
                                                   #result_type, prefix, postfix)

        elif self.element_type == 4: # CSHEAR
            n, nelements, ntotal = self._oes_cshear(data, ndata, dt, is_magnitude_phase,
                                                    result_type, prefix, postfix)

        elif self.element_type in [11, 12, 13, 14]:  # springs
            # 11-CELAS1
            # 12-CELAS2
            # 13-CELAS3
            # 14-CELAS4
            n, nelements, ntotal = self._oes_celas(data, ndata, dt, is_magnitude_phase,
                                                   result_type, prefix, postfix)

        elif self.element_type == 34: # CBAR
            n, nelements, ntotal = self._oes_cbar_34(data, ndata, dt, is_magnitude_phase,
                                                     result_type, prefix, postfix)

        elif self.element_type == 83:
            # 83: TRIA3
            n, nelements, ntotal = self._oes_ctria3(data, ndata, dt, is_magnitude_phase,
                                                    result_type, prefix, postfix)

        #elif self.element_type in [64, 70, 75, 82, 144]:  # bilinear plates
            # 64-CQUAD8
            # 70-CTRIAR
            # 75-CTRIA6
            # 82-CQUADR
            # 144-CQUAD4-bilinear
            #n, nelements, ntotal = self._oes_cquad4_144(data, ndata, dt, is_magnitude_phase,
                                                        #result_type, prefix, postfix)
        else:
            #msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
            msg = self.code_information()
            print(msg)
            return self._not_implemented_or_skip(data, ndata, msg)

        if nelements is None:
            return n
        assert ndata > 0, ndata
        assert nelements > 0, 'nelements=%r element_type=%s element_name=%r' % (nelements, self.element_type, self.element_name)
        #assert ndata % ntotal == 0, '%s n=%s nwide=%s len=%s ntotal=%s' % (self.element_name, ndata % ntotal, ndata % self.num_wide, ndata, ntotal)
        assert self.num_wide * 4 == ntotal, 'numwide*4=%s ntotal=%s' % (self.num_wide * 4, ntotal)
        assert self.thermal == 0, "thermal = %%s" % self.thermal
        assert n > 0, f'n = {n} result_name={result_name}'
        return n

    def _read_oes1_loads(self, data, ndata: int):
        """Reads OES self.thermal=0 stress/strain"""
        prefix, postfix = self.get_oes_prefix_postfix()
        result_type = self.result_type

        #self._apply_oes_ato_crm_psd_rms_no('') # TODO: just testing
        n = 0
        is_magnitude_phase = self.is_magnitude_phase()
        dt = self.nonlinear_factor

        #flag = 'element_id'
        if self.is_stress:
            result_name = 'stress'
            stress_name = 'STRESS'
        else:
            result_name = 'strain'
            stress_name = 'STRAIN'

        #if self.is_stress:
            #_result_name, _class_obj = self.get_stress_mapper()
        if self.table_name_str == 'OESXRMS1':
            assert self.sort_method == 1, self.code_information()

        if self._results.is_not_saved(result_name):
            return ndata
        if self.element_type in [1, 3, 10]:  # rods
            # 1-CROD
            # 3-CTUBE
            # 10-CONROD
            n, nelements, ntotal = self._oes_crod(data, ndata, dt, is_magnitude_phase,
                                                  result_type, prefix, postfix)

        elif self.element_type == 2: # CBEAM
            n, nelements, ntotal = self._oes_cbeam(data, ndata, dt, is_magnitude_phase,
                                                   result_type, prefix, postfix)

        elif self.element_type == 4: # CSHEAR
            n, nelements, ntotal = self._oes_cshear(data, ndata, dt, is_magnitude_phase,
                                                    result_type, prefix, postfix)

        elif self.element_type in [11, 12, 13, 14]:  # springs
            # 11-CELAS1
            # 12-CELAS2
            # 13-CELAS3
            # 14-CELAS4
            n, nelements, ntotal = self._oes_celas(data, ndata, dt, is_magnitude_phase,
                                                   result_type, prefix, postfix)

        elif self.element_type == 34: # CBAR
            n, nelements, ntotal = self._oes_cbar_34(data, ndata, dt, is_magnitude_phase,
                                                     result_type, prefix, postfix)

        elif self.element_type in [39, 67, 68, 255]: # solid stress
            # 39-CTETRA
            # 67-CHEXA
            # 68-CPENTA
            # 255-CPYRAM
            n, nelements, ntotal = self._oes_csolid(data, ndata, dt, is_magnitude_phase,
                                                    result_type, prefix, postfix)
        elif self.element_type in [140]:
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

        elif self.element_type in [160, 163, 166,
                                   161, # centroid
                                   165,]:
            # nonlinear hyperelastic solids
            # 160-CPENTAFD
            # 163-CHEXAFD
            # 166-CTETRAFD

            # centroid??
            # 161-CTETRAFD

            # many nodes?
            # 165-CPENTAFD
            n, nelements, ntotal = self._oes_csolid_linear_hyperelastic(data, ndata, dt, is_magnitude_phase,
                                                                        result_type, prefix, postfix)
        elif self.element_type in [202, 204,
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

        elif self.element_type in [300, 301, 302, 303]: # solid stress
            # solids without stress eigenvectors
            # 300-CHEXA
            # 301-CPENTA
            # 302-CTETRA
            # 303-CPYRAM
            n, nelements, ntotal = self._oes_csolid2(data, ndata, dt, is_magnitude_phase,
                                                     result_type, prefix, postfix)
        elif self.element_type in [306, 307]:
            # 306-CHEXALN
            # 307-CPENTALN
            n, nelements, ntotal = self._oes_csolid_composite(data, ndata, dt, is_magnitude_phase,
                                                              result_type, prefix, postfix)
        #=========================
        # plates
        elif self.element_type in [33, 228]:
            # 33: CQUAD4-centroidal
            # 228: CQUADR-centroidal
            n, nelements, ntotal = self._oes_cquad4_33(data, ndata, dt, is_magnitude_phase,
                                                       result_type, prefix, postfix)

        elif self.element_type in [74, 227, ]: # 229???
            # 74: TRIA3
            # 227: TRIAR
            n, nelements, ntotal = self._oes_ctria3(data, ndata, dt, is_magnitude_phase,
                                                    result_type, prefix, postfix)

        elif self.element_type in [64, 70, 75, 82, 144]:  # bilinear plates
            # 64-CQUAD8
            # 70-CTRIAR
            # 75-CTRIA6
            # 82-CQUADR
            # 144-CQUAD4-bilinear
            n, nelements, ntotal = self._oes_cquad4_144(data, ndata, dt, is_magnitude_phase,
                                                        result_type, prefix, postfix)

        elif self.element_type in [88, 90]: # nonlinear shells
            # 88-CTRIA3NL
            # 90-CQUAD4NL
            n, nelements, ntotal = self._oes_shells_nonlinear(data, ndata, dt, is_magnitude_phase,
                                                              result_type, prefix, postfix)

        elif self.element_type in [95, 96, 97, 98, 232, 233]: # composite shell
            # 95 - CQUAD4
            # 96 - CQUAD8
            # 97 - CTRIA3
            # 98 - CTRIA6 (composite)
            # 232 - QUADRLC (CQUADR-composite)
            # 233 - TRIARLC (CTRIAR-composite)
            n, nelements, ntotal = self._oes_shells_composite(data, ndata, dt, is_magnitude_phase,
                                                              result_type, prefix, postfix)

        elif self.element_type == 53: # axial plates - ctriax6
            n, nelements, ntotal = self._oes_ctriax6(data, ndata, dt, is_magnitude_phase,
                                                     result_type, prefix, postfix)

        elif self.element_type == 102: # cbush
            n, nelements, ntotal = self._oes_cbush(data, ndata, dt, is_magnitude_phase,
                                                   result_type, prefix, postfix)

        elif self.element_type == 40:  # cbush1d
            n, nelements, ntotal = self._oes_cbush1d(data, ndata, dt, is_magnitude_phase,
                                                     result_type, prefix, postfix)

        elif self.element_type in [87, 89, 92]:  # nonlinear rods
            # 87-CTUBENL
            # 89-RODNL
            # 92-CONRODNL
            n, nelements, ntotal = self._oes_crod_nonlinear(data, ndata, dt, is_magnitude_phase,
                                                            result_type, prefix, postfix)

        elif self.element_type in [224, 225]: # nonlinear spring
            # 224-CELAS1
            # 225-CELAS3
            # NonlinearSpringStress
            n, nelements, ntotal = self._oes_celas_nonlinear(data, ndata, dt, is_magnitude_phase,
                                                             result_type, prefix, postfix)

        elif self.element_type == 69:  # cbend
            # 69-CBEND
            n, nelements, ntotal = self._oes_cbend(data, ndata, dt, is_magnitude_phase,
                                                   result_type, prefix, postfix)

        elif self.element_type == 86:  # cgap
            # 86-GAPNL
            n, nelements, ntotal = self._oes_cgap_nonlinear(data, ndata, dt, is_magnitude_phase,
                                                            result_type, prefix, postfix)

        elif self.element_type == 94:
            # 94-BEAMNL
            n, nelements, ntotal = self._oes_cbeam_nonlinear(data, ndata, dt, is_magnitude_phase,
                                                             result_type, prefix, postfix)

        elif self.element_type in [85, 91, 93, 256]:
            # 256-PYRAM
            n, nelements, ntotal = self._oes_csolid_nonlinear(data, ndata, dt, is_magnitude_phase,
                                                              result_type, prefix, postfix)

        elif self.element_type == 100:  # bars
            # 100-BARS
            n, nelements, ntotal = self._oes_cbar_100(data, ndata, dt, is_magnitude_phase,
                                                      result_type, prefix, postfix)

        elif self.element_type in [145, 146, 147]:
            n, nelements, ntotal = self._oes_vu_solid(data, ndata, dt, is_magnitude_phase, stress_name,
                                                      result_type, prefix, postfix)

        #-----------------------------------------------------------------------

        elif self.element_type == 139:
            n, nelements, ntotal = self._oes_hyperelastic_quad(data, ndata, dt, is_magnitude_phase,
                                                               result_type, prefix, postfix)

        elif self.element_type == 189: # VUQUAD
            n, nelements, ntotal = self._oes_vu_quad(data, ndata, dt, is_magnitude_phase,
                                                     result_type, prefix, postfix)

        elif self.element_type == 226:
            # 226-BUSHNL
            n, nelements, ntotal = self._oes_cbush_nonlinear(data, ndata, dt, is_magnitude_phase,
                                                             result_type, prefix, postfix)
        elif self.element_type in [271, 275]:
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

        elif self.element_type in [276, 277, 278]:
            # 276-CPLSTS4
            # 277-CPLSTS6
            # 278-CPLSTS8
            n, nelements, ntotal = self._oes_plate_stress_68(data, ndata, dt, is_magnitude_phase,
                                                             stress_name, prefix, postfix)


        elif self.element_type == 35: # CON
            return ndata

        elif self.element_type in [60, 61]:
            # 60-DUM8
            # 61-DUM9
            return ndata

        elif self.element_type == 101: # AABSF
            return ndata

        #elif self.element_type in [140, 201]:
            ## 140-HEXA8FD, 201-QUAD4FD
            #return ndata

        elif self.element_type in [47, 48, 189, 190]:
            # 47-AXIF2
            # 48-AXIF3
            # 189-???
            # 190-VUTRIA
            return ndata
        elif self.element_type == 191:
            # 191-VUBEAM
            return ndata
        elif self.element_type in [50, 51, 203]:
            # 203-SLIF1D?
            # 50-SLOT3
            # 51-SLOT4
            return ndata
        elif self.element_type in [162, 164, 167, 168,
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
            return self._not_implemented_or_skip(data, ndata, self.code_information())
        else:
            #msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
            msg = self.code_information()
            return self._not_implemented_or_skip(data, ndata, msg)
        if nelements is None:
            return n

        #self.check_element_ids()
        assert ndata > 0, ndata
        assert nelements > 0, f'nelements={nelements} element_type={self.element_type} element_name={self.element_name!r}'
        #assert ndata % ntotal == 0, '%s n=%s nwide=%s len=%s ntotal=%s' % (self.element_name, ndata % ntotal, ndata % self.num_wide, ndata, ntotal)
        assert self.num_wide * 4 * self.factor == ntotal, f'numwide*4={self.num_wide*4} ntotal={ntotal} element_name={self.element_name!r}\n{self.code_information()}'
        assert self.thermal == 0, "thermal = %%s" % self.thermal
        assert n is not None and n > 0, f'n={n} result_name={result_name}\n{self.code_information()}'

        #if self.is_sort2:
            #assert len(np.unique(self.obj._times)) == len(self.obj._times), f'{self.obj._times.tolist()}\n{self.code_information()}'
        return n

    def check_element_ids(self):
        if self.read_mode == 1:
            return
        if self.is_sort1:
            obj = self.obj
            if obj is None:
                raise RuntimeError('obj is None...\n' + self.code_information())
            if hasattr(obj, 'element_node'):
                eids = obj.element_node[:, 0]
            elif hasattr(obj, 'element_layer'):
                eids = obj.element_layer[:, 0]
            elif hasattr(obj, 'element'):
                eids = obj.element
            else:
                print(self.code_information())
                raise RuntimeError(''.join(obj.get_stats()))
            if eids.min() <= 0:
                #print(obj.code_information())
                print(''.join(obj.get_stats()))
                raise RuntimeError(f'{self.element_name}-{self.element_type}: {eids}')
        #else:
            #assert._times

    def _create_nodes_object(self, nnodes, result_name, slot, obj_vector):
        """same as _create_oes_object4 except it adds to the nnodes parameter"""
        auto_return = False
        #is_vectorized = True
        is_vectorized = self._is_vectorized(obj_vector)
        #print("vectorized...read_mode=%s...%s; %s" % (self.read_mode, result_name, is_vectorized))

        if is_vectorized:
            if self.read_mode == 1:
                #print('oes-self.nonlinear_factor =', self.nonlinear_factor)
                #print(self.data_code)
                self.create_transient_object(result_name, slot, obj_vector)
                #print("read_mode 1; ntimes=%s" % self.obj.ntimes)
                self.result_names.add(result_name)
                #print('self.obj =', self.obj)
                self.obj.nnodes += nnodes
                auto_return = True
            elif self.read_mode == 2:
                self.code = self._get_code()
                #self.log.info("code = %s" % str(self.code))
                #print("code = %s" % str(self.code))

                # if this is failing, you probably set obj_vector to None...
                try:
                    self.obj = slot[self.code]
                except KeyError:
                    msg = 'Could not find key=%s in result=%r\n' % (self.code, result_name)
                    msg += "There's probably an extra check for read_mode=1..."
                    self.log.error(msg)
                    raise
                #self.obj.update_data_code(self.data_code)
                build_obj(self.obj)

            else:  # not vectorized
                auto_return = True
        else:
            auto_return = True
        return auto_return, is_vectorized

    def _create_ntotal_object(self, ntotal, result_name, slot, obj_vector):
        """same as _create_oes_object4 except it adds to the ntotal parameter"""
        auto_return = False
        #is_vectorized = True
        is_vectorized = self._is_vectorized(obj_vector)
        #print("vectorized...read_mode=%s...%s; %s" % (self.read_mode, result_name, is_vectorized))

        if is_vectorized:
            if self.read_mode == 1:
                #print('oes-self.nonlinear_factor =', self.nonlinear_factor)
                #print(self.data_code)
                self.create_transient_object(result_name, slot, obj_vector)
                #print("read_mode 1; ntimes=%s" % self.obj.ntimes)
                self.result_names.add(result_name)
                #print('self.obj =', self.obj)
                self.obj.ntotal += ntotal
                auto_return = True
            elif self.read_mode == 2:
                self.code = self._get_code()
                #self.log.info("code = %s" % str(self.code))
                #print("code = %s" % str(self.code))

                # if this is failing, you probably set obj_vector to None...
                try:
                    self.obj = slot[self.code]
                except KeyError:
                    msg = 'Could not find key=%s in result=%r\n' % (self.code, result_name)
                    msg += "There's probably an extra check for read_mode=1..."
                    self.log.error(msg)
                    raise
                #self.obj.update_data_code(self.data_code)
                build_obj(self.obj)

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
        n = 0
        if self.is_stress:
            if prefix == '' and postfix == '':
                prefix = 'stress.'
            obj_real = RealSpringStressArray
            obj_complex = ComplexSpringStressArray
            if self.element_type == 11:
                result_name = prefix + 'celas1_stress' + postfix
            elif self.element_type == 12:
                result_name = prefix + 'celas2_stress' + postfix
            elif self.element_type == 13:
                result_name = prefix + 'celas3_stress' + postfix
            elif self.element_type == 14:
                result_name = prefix + 'celas4_stress' + postfix
            else:
                raise RuntimeError(self.element_type)
        else:
            if prefix == '' and postfix == '':
                prefix = 'strain.'
            obj_real = RealSpringStrainArray
            obj_complex = ComplexSpringStrainArray
            if self.element_type == 11:
                result_name = prefix + 'celas1_strain' + postfix
            elif self.element_type == 12:
                result_name = prefix + 'celas2_strain' + postfix
            elif self.element_type == 13:
                result_name = prefix + 'celas3_strain' + postfix
            elif self.element_type == 14:
                result_name = prefix + 'celas4_strain' + postfix
            else:
                raise RuntimeError(self.element_type)

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
                n = nelements * 4 * self.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, 2)
                obj._times[obj.itime] = dt
                self.obj_set_element(obj, itotal, itotal2, data, nelements)

                #(eid_device, stress)
                obj.data[obj.itime, itotal:itotal2, 0] = floats[:, 1].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if is_vectorized and self.use_vector:  # pragma: no cover
                    self.log.debug('vectorize CELASx real SORT%s' % self.sort_method)
                fmt1 = mapfmt(self._endian + self._analysis_code_fmt + b'f', self.size)
                struct1 = Struct(fmt1)
                for i in range(nelements):
                    edata = data[n:n+ntotal]
                    out = struct1.unpack(edata)
                    (eid_device, ox) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    if eid <= 0: # pragma: no cover
                        msg = 'table_name=%s sort_method=%s eid_device=%s nonlinear_factor=%s'  % (
                            self.table_name_str, self.sort_method,
                            eid_device, self.nonlinear_factor)
                        raise RuntimeError(msg)
                    if self.is_debug_file:
                        self.binary_debug.write('  eid=%i result%i=[%i, %f]\n' % (
                            eid, i, eid_device, ox))
                    obj.add_sort1(dt, eid, ox)
                    n += ntotal
        elif self.format_code in [2, 3] and self.num_wide == 3:  # imag
            ntotal = 12 * self.factor
            nelements = ndata // ntotal

            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_complex)
            if auto_return:
                return nelements * ntotal, None, None

            obj = self.obj
            assert obj is not None, self.code_information()
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, 3).copy()
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
                if is_vectorized and self.use_vector:  # pragma: no cover
                    self.log.debug('vectorize CELASx imag SORT%s' % self.sort_method)
                struct1 = Struct(self._endian + mapfmt(self._analysis_code_fmt + b'2f', self.size))
                for i in range(nelements):
                    edata = data[n:n + ntotal]
                    out = struct1.unpack(edata)
                    (eid_device, axial_real, axial_imag) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)

                    if is_magnitude_phase:
                        axial = polar_to_real_imag(axial_real, axial_imag)
                    else:
                        axial = complex(axial_real, axial_imag)

                    if self.is_debug_file:
                        self.binary_debug.write('  eid=%i result%i=[%i, %f, %f]\n' % (
                            eid, i, eid_device, axial_real, axial_imag))
                    obj.add_sort1(dt, eid, axial)
                    n += ntotal
        elif self.format_code == 1 and self.num_wide == 3: # random
            raise RuntimeError(self.code_information())
            #msg = self.code_information()
            #return self._not_implemented_or_skip(data, ndata, msg)
        else:  # pragma: no cover
            raise RuntimeError(self.code_information())
        return n, nelements, ntotal

    def _oes_crod(self, data, ndata, dt, is_magnitude_phase,
                  result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 1 : CROD
         - 3 : CTUBE
         - 10 : CONROD

        """
        n = 0
        if self.is_stress:
            obj_vector_real = RealRodStressArray
            obj_vector_complex = ComplexRodStressArray
            obj_vector_random = RandomRodStressArray
            if self.element_type == 1: # CROD
                result_name = prefix + 'crod_stress' + postfix
            elif self.element_type == 3:  # CTUBE
                result_name = prefix + 'ctube_stress' + postfix
            elif self.element_type == 10:  # CONROD
                result_name = prefix + 'conrod_stress' + postfix
            else:  # pragma: no cover
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
        else:
            obj_vector_real = RealRodStrainArray
            obj_vector_complex = ComplexRodStrainArray
            obj_vector_random = RandomRodStrainArray
            if self.element_type == 1: # CROD
                result_name = prefix + 'crod_strain' + postfix
            elif self.element_type == 3:  # CTUBE
                result_name = prefix + 'ctube_strain' + postfix
            elif self.element_type == 10:  # CONROD
                result_name = prefix + 'conrod_strain' + postfix
            else:  # pragma: no cover
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)

        #result_name, unused_is_random = self._apply_oes_ato_crm_psd_rms_no(result_name)
        slot = self.get_result(result_name)
        if result_type == 0 and self.num_wide == 5:  # real
            ntotal = 5 * 4 * self.factor
            nelements = ndata // ntotal

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return nelements * self.num_wide * 4, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * 4 * self.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, 5)
                obj._times[obj.itime] = dt
                self.obj_set_element(obj, itotal, itotal2, data, nelements)

                #[axial, torsion, SMa, SMt]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if is_vectorized and self.use_vector:  # pragma: no cover
                    self.log.debug('vectorize CROD real SORT%s' % self.sort_method)
                if self.is_debug_file:
                    self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    self.binary_debug.write('  #elementi = [eid_device, axial, axial_margin, torsion, torsion_margin]\n')
                    self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                struct1 = Struct(self._endian + mapfmt(self._analysis_code_fmt + b'4f', self.size))
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]
                    out = struct1.unpack(edata)
                    (eid_device, axial, axial_margin, torsion, torsion_margin) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)

                    if self.is_debug_file:
                        self.binary_debug.write('  eid=%i; C=[%s]\n' % (
                            eid, ', '.join(['%r' % di for di in out])))
                    obj.add_sort1(dt, eid, axial, axial_margin, torsion, torsion_margin)
                    n += ntotal
        elif result_type == 1 and self.num_wide == 5: # imag
            ntotal = 20 * self.factor
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_complex)
            if auto_return:
                return nelements * ntotal, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, 5)
                obj._times[obj.itime] = dt
                self.obj_set_element(obj, itotal, itotal2, data, nelements)

                real_imag = apply_mag_phase(floats, is_magnitude_phase, [1, 3], [2, 4])
                obj.data[obj.itime, itotal:itotal2, :] = real_imag

                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if is_vectorized and self.use_vector:  # pragma: no cover
                    self.log.debug('vectorize CROD imag SORT%s' % self.sort_method)
                fmt = mapfmt(self._endian + self._analysis_code_fmt + b'4f', self.size)
                struct1 = Struct(fmt)
                for unused_i in range(nelements):
                    edata = data[n:n + ntotal]
                    out = struct1.unpack(edata)
                    (eid_device, axial_real, axial_imag, torsion_real, torsion_imag) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)

                    if is_magnitude_phase:
                        axial = polar_to_real_imag(axial_real, axial_imag)
                        torsion = polar_to_real_imag(torsion_real, torsion_imag)
                    else:
                        axial = complex(axial_real, axial_imag)
                        torsion = complex(torsion_real, torsion_imag)

                    obj.add_sort1(dt, eid, axial, torsion)
                    n += ntotal
        #elif self.format_code in [2, 3] and self.num_wide == 8:  # is this imag ???
            #ntotal = 32
            #s = self.self.struct_i
            #nelements = ndata // ntotal
            #for i in range(nelements):
                #edata = data[n:n + 4]
                #eid_device, = s.unpack(edata)
                #assert eid > 0, eid
                #n += ntotal
        elif result_type == 2 and self.num_wide == 3: # random
            ntotal = 3 * 4 * self.factor
            nelements = ndata // ntotal

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_random)
            if auto_return:
                return nelements * ntotal, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                if is_vectorized and self.use_vector:  # pragma: no cover
                    self.log.debug('vectorize CROD random SORT%s' % self.sort_method)
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 3)
                obj._times[obj.itime] = dt
                self.obj_set_element(obj, itotal, itotal2, data, nelements)

                #[axial, torsion, SMa, SMt]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                n = oes_crod_random_3(self, data, ndata, obj, nelements, ntotal)

        else:  # pragma: no cover
            raise RuntimeError(self.code_information())
        assert self.num_wide * 4 * self.factor == ntotal, f'numwide*4={self.num_wide*4} ntotal={ntotal} element_name={self.element_name!r}\n{self.code_information()}'
        return n, nelements, ntotal

    def _oes_cbeam(self, data, ndata, dt, is_magnitude_phase,
                   result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 2 : CBEAM

        """
        n = 0
        ## TODO: fix method to follow correct pattern...regarding???

        if self.is_stress:
            result_name = prefix + 'cbeam_stress' + postfix
        else:
            result_name = prefix + 'cbeam_strain' + postfix
        table_name_bytes = self.table_name
        assert isinstance(table_name_bytes, bytes), table_name_bytes
        assert table_name_bytes in TABLES_BYTES, table_name_bytes

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)

        slot = self.get_result(result_name)
        if result_type == 0 and self.num_wide == 111:  # real
            # TODO: vectorize
            ntotal = 444 * self.factor # 44 + 10*40  (11 nodes)

            if self.is_stress:
                obj_vector_real = RealBeamStressArray
            else:
                obj_vector_real = RealBeamStrainArray

            nelements = ndata // ntotal
            nlayers = nelements * 11
            auto_return, is_vectorized = self._create_oes_object4(
                nlayers, result_name, slot, obj_vector_real)
            if auto_return:
                self._data_factor = 11
                return nelements * self.num_wide * 4, None, None
            obj = self.obj

            ntotal = self.num_wide * 4 * self.factor
            nelements = ndata // ntotal
            if self.use_vector and is_vectorized and 0:
                raise NotImplementedError('CBEAM-2-real not vectorized')
            else:
                if is_vectorized and self.use_vector:  # pragma: no cover
                    self.log.debug('vectorize CBEAM real SORT%s' % self.sort_method)
                n = oes_cbeam_real_111(self, data,
                                       obj,
                                       nelements, dt)

        elif result_type == 1 and self.num_wide == 111:  # imag and random?
            # definitely complex results for MSC Nastran 2016.1

            ntotal = 444 * self.factor # 44 + 10*40  (11 nodes)
            nelements = ndata // ntotal

            if self.is_stress:
                obj_vector_complex = ComplexBeamStressArray
            else:
                obj_vector_complex = ComplexBeamStrainArray

            nlayers = nelements * 11
            auto_return, is_vectorized = self._create_oes_object4(
                nlayers, result_name, slot, obj_vector_complex)
            if auto_return:
                self._data_factor = 11
                return nelements * ntotal, None, None

            obj = self.obj

            nnodes = 10  # 11-1
            #ntotal = self.num_wide * 4
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.itotal
                itotal2 = itotal + nelements * 11

                # chop off eid
                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, 111)[:, 1:]
                floats2 = floats.reshape(nelements * 11, 10).copy()

                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype8).reshape(nelements, 111)
                    eids = ints[:, 0] // 10
                    eids2 = array([eids] * 11, dtype=self.idtype8).T.ravel()

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
                if is_vectorized and self.use_vector:  # pragma: no cover
                    self.log.debug('vectorize CBEAM imag SORT%s' % self.sort_method)
                n = oes_cbeam_complex_111(self, data, obj,
                                          nelements, nnodes,
                                          is_magnitude_phase)

        elif result_type == 2 and self.num_wide == 67: # random
            # TODO: vectorize
            ntotal = 268 # 1 + 11*6  (11 nodes)

            if self.is_stress:
                obj_vector_random = RandomBeamStressArray
            else:
                obj_vector_random = RandomBeamStrainArray

            nelements = ndata // ntotal
            nlayers = nelements * 11
            auto_return, is_vectorized = self._create_oes_object4(
                nlayers, result_name, slot, obj_vector_random)
            if auto_return:
                self._data_factor = 11
                return nelements * self.num_wide * 4, None, None
            obj = self.obj

            nnodes = 10  # 11-1
            ntotal = self.num_wide * 4
            nelements = ndata // ntotal
            if self.use_vector and is_vectorized and 0:
                raise NotImplementedError('CBEAM-2-random not vectorized')
            else:
                if is_vectorized and self.use_vector:  # pragma: no cover
                    self.log.debug('vectorize CBEAM random SORT%s' % self.sort_method)
                n = oes_cbeam_random_67(self, data, obj,
                                        nelements, nnodes, dt)

        elif result_type == 1 and self.num_wide in [67] and table_name_bytes in [b'OESXNO1']:  # CBEAM
            # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1081x.op2
            msg = 'skipping freq CBEAM; numwide=67'
            n = self._not_implemented_or_skip(data, ndata, msg)
            nelements = None
            ntotal = None
        elif result_type == 2 and self.num_wide in [67] and table_name_bytes in [b'OESXNO1']:  # CBEAM
            #C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1081x.op2
            msg = 'skipping random CBEAM; numwide=67'
            n = self._not_implemented_or_skip(data, ndata, msg)
            nelements = None
            ntotal = None
        else:  # pragma: no cover
            raise RuntimeError(self.code_information())
        return n, nelements, ntotal

    def _oes_cshear(self, data, ndata, dt, is_magnitude_phase,
                    result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 4 : CSHEAR
        """
        n = 0
        # 4-CSHEAR
        if self.is_stress:
            obj_vector_real = RealShearStressArray
            obj_vector_complex = ComplexShearStressArray
            obj_vector_random = RandomShearStressArray
            result_name = prefix + 'cshear_stress' + postfix
        else:
            obj_vector_real = RealShearStrainArray
            obj_vector_complex = ComplexShearStrainArray
            obj_vector_random = RandomShearStrainArray
            result_name = prefix + 'cshear_strain' + postfix

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)

        slot = self.get_result(result_name)
        if result_type == 0 and self.num_wide == 4:  # real
            ntotal = 16  # 4*4
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return nelements * self.num_wide * 4, None, None

            obj = self.obj
            assert obj is not None
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 4)
                itime = obj.itime
                obj._times[itime] = dt
                self.obj_set_element(obj, itotal, itotal2, data, nelements)

                #[max_strain, avg_strain, margin]
                obj.data[itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if is_vectorized and self.use_vector:  # pragma: no cover
                    self.log.debug('vectorize CSHEAR real SORT%s' % self.sort_method)
                struct1 = Struct(self._endian + self._analysis_code_fmt + b'3f')
                for unused_i in range(nelements):
                    edata = data[n:n + ntotal]
                    out = struct1.unpack(edata)  # num_wide=5
                    if self.is_debug_file:
                        self.binary_debug.write('CSHEAR-4 - %s\n' % str(out))

                    (eid_device, max_strain, avg_strain, margin) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    obj.add_sort1(dt, eid, max_strain, avg_strain, margin)
                    n += ntotal

        elif result_type == 1 and self.num_wide == 5:  # imag
            ntotal = 20 * self.factor # 4*5
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_complex)
            if auto_return:
                return nelements * ntotal, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, 5).copy()
                obj._times[obj.itime] = dt
                self.obj_set_element(obj, itotal, itotal2, data, nelements)

                #(eid_device, etmaxr, etmaxi, etavgr, etavgi)
                real_imag = apply_mag_phase(floats, is_magnitude_phase, [1, 3], [2, 4])
                obj.data[obj.itime, itotal:itotal2, :] = real_imag
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if is_vectorized and self.use_vector:  # pragma: no cover
                    self.log.debug('vectorize CSHEAR imag SORT%s' % self.sort_method)
                struct1 = Struct(self._endian + mapfmt(self._analysis_code_fmt + b'4f', self.size))
                for unused_i in range(nelements):
                    edata = data[n:n + ntotal]
                    out = struct1.unpack(edata)  # num_wide=5
                    if self.is_debug_file:
                        self.binary_debug.write('CSHEAR-4 - %s\n' % str(out))
                    (eid_device, etmaxr, etmaxi, etavgr, etavgi) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)

                    if is_magnitude_phase:
                        etmax = polar_to_real_imag(etmaxr, etmaxi)
                        etavg = polar_to_real_imag(etavgr, etavgi)
                    else:
                        etmax = complex(etmaxr, etmaxi)
                        etavg = complex(etavgr, etavgi)
                    obj.add_sort1(dt, eid, etmax, etavg)
                    n += ntotal
        elif result_type == 2 and self.num_wide == 3: # random
            ntotal = 12  # 3*4
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_random)
            if auto_return:
                return nelements * self.num_wide * 4, None, None

            obj = self.obj
            assert obj is not None
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * 4 * self.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 3)
                itime = obj.itime
                obj._times[itime] = dt
                self.obj_set_element(obj, itotal, itotal2, data, nelements)

                #[max_strain, avg_strain, margin]
                obj.data[itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if is_vectorized and self.use_vector:  # pragma: no cover
                    self.log.debug('vectorize CSHEAR random SORT%s' % self.sort_method)
                n = oes_cshear_random_3(self, data, obj,
                                        nelements, ntotal)

        else:  # pragma: no cover
            raise RuntimeError(self.code_information())
        return n, nelements, ntotal

    def _oes_cbar_34(self, data, ndata, dt, is_magnitude_phase,
                     result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 34 : CBAR

        """
        #if isinstance(self.nonlinear_factor, float):
            #self.sort_bits[0] = 1 # sort2
            #self.sort_method = 2

        n = 0
        if self.is_stress:
            result_name = prefix + 'cbar_stress' + postfix
        else:
            result_name = prefix + 'cbar_strain' + postfix

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)
        if result_type == 0 and self.num_wide == 16:  # real
            if self.is_stress:
                obj_vector_real = RealBarStressArray
            else:
                obj_vector_real = RealBarStrainArray

            ntotal = 64 * self.factor  # 16*4
            nelements = ndata // ntotal
            #print('CBAR nelements =', nelements)

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return ndata, None, None

            if self.is_debug_file:
                self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                self.binary_debug.write('  #elementi = [eid_device, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,\n')
                self.binary_debug.write('                           s1b, s2b, s3b, s4b, smaxb, sminb,        MSc]\n')
                self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                # self.itime = 0
                # self.ielement = 0
                # self.itotal = 0
                #self.ntimes = 0
                #self.nelements = 0
                n = nelements * self.num_wide * 4

                ielement = obj.ielement
                ielement2 = ielement + nelements
                obj._times[obj.itime] = dt
                self.obj_set_element(obj, ielement, ielement2, data, nelements)

                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, 16)

                #[s1a, s2a, s3a, s4a, axial, smaxa, smina, margin_tension,
                # s1b, s2b, s3b, s4b,        smaxb, sminb, margin_compression]
                obj.data[obj.itime, ielement:ielement2, :] = floats[:, 1:].copy()
                obj.itotal = ielement2
                obj.ielement = ielement2
            else:
                if is_vectorized and self.use_vector:  # pragma: no cover
                    self.log.debug('vectorize CBAR real SORT%s' % self.sort_method)
                n = oes_cbar_real_16(self, data, obj, nelements, ntotal, dt)
        elif result_type == 1 and self.num_wide == 19:  # imag
            if self.is_stress:
                obj_vector_complex = ComplexBarStressArray
            else:
                obj_vector_complex = ComplexBarStrainArray

            ntotal = 76 * self.factor
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_complex)
            if auto_return:
                return ndata, None, None

            if self.is_debug_file:
                self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                self.binary_debug.write('  #elementi = [eid_device, s1a, s2a, s3a, s4a, axial,\n')
                self.binary_debug.write('                           s1b, s2b, s3b, s4b]\n')
                self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.itotal
                itotal2 = itotal + nelements
                ielement2 = itotal2

                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, 19).copy()
                obj._times[obj.itime] = dt
                self.obj_set_element(obj, itotal, itotal2, data, nelements)

                isave1 = [1, 2, 3, 4, 5, 11, 12, 13, 14]
                isave2 = [6, 7, 8, 9, 10, 15, 16, 17, 18]
                real_imag = apply_mag_phase(floats, is_magnitude_phase, isave1, isave2)
                obj.data[obj.itime, itotal:itotal2, :] = real_imag

                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if is_vectorized and self.use_vector:  # pragma: no cover
                    self.log.debug('vectorize CBAR imag SORT%s' % self.sort_method)
                n = oes_cbar_complex_19(self, data, obj, nelements, ntotal, is_magnitude_phase)
        elif result_type == 2 and self.num_wide == 19: # random strain?
            raise RuntimeError(self.code_information())

        elif result_type in [1, 2] and self.num_wide == 10:  # random
            # random stress/strain per example
            #
            # DMAP says random stress has num_wide=10 and
            # random strain has numwide=19, but it's wrong...maybe???
            #
            # format_code = 1 - NO/RMS (SORT1 regardless of whether this is a SORT2 table or not)
            # format_code = 2 - ATO/PSD/CRM (actually SORT2)
            #
            element_id = self.nonlinear_factor
            if self.is_stress:
                obj_vector_random = RandomBarStressArray
            else:
                obj_vector_random = RandomBarStrainArray
            self.data_code['nonlinear_factor'] = element_id

            ntotal = 10 * self.size
            nelements = ndata // ntotal
            #print(f'CBAR* nelements={nelements}')
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_random)
            if auto_return:
                return ndata, None, None

            if self.is_debug_file:
                self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                self.binary_debug.write('  #elementi = [eid_device, s1a, s2a, s3a, s4a, axial,\n')
                self.binary_debug.write('                           s1b, s2b, s3b, s4b]\n')
                self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            obj = self.obj
            if self.use_vector and is_vectorized and 0:  # pragma: no cover
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

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 10)

                #[s1a, s2a, s3a, s4a, axial,
                # s1b, s2b, s3b, s4b]
                obj.data[obj.itime, ielement:ielement2, :] = floats[:, 1:].copy()
                obj.itotal = ielement2
                obj.ielement = ielement2
            else:
                if is_vectorized and self.use_vector and obj.itime == 0:  # pragma: no cover
                    self.log.debug('vectorize CBAR random SORT%s' % self.sort_method)
                n = oes_cbar_random_10(self, data, obj,
                                       nelements, ntotal)
        else:  # pragma: no cover
            raise RuntimeError(self.code_information())
        return n, nelements, ntotal

    def _oes_csolid(self, data, ndata, dt, is_magnitude_phase,
                    result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 39 : CTETRA
         - 67 : CHEXA
         - 68 : CPENTA

        """
        n = 0
        if self.is_stress:
            if prefix == '' and postfix == '':
                prefix = 'stress.'
            obj_vector_real = RealSolidStressArray
            obj_vector_complex = ComplexSolidStressArray
            obj_vector_random = RandomSolidStressArray
            if self.element_type == 39: # CTETRA
                nnodes_expected = 5  # 1 centroid + 4 corner points
                result_name = prefix + 'ctetra_stress' + postfix
                element_name = 'CTETRA4'
            elif self.element_type == 67:  # CHEXA
                nnodes_expected = 9
                result_name = prefix + 'chexa_stress' + postfix
                element_name = 'CHEXA8'
            elif self.element_type == 68:  # CPENTA
                nnodes_expected = 7
                result_name = prefix + 'cpenta_stress' + postfix
                element_name = 'CPENTA6'
            elif self.element_type == 255:  # CPYRAM
                nnodes_expected = 6
                result_name = prefix + 'cpyram_stress' + postfix
                element_name = 'CPYRAM5'
            else:  # pragma: no cover
                raise RuntimeError(self.code_information())
        else:
            if prefix == '' and postfix == '':
                prefix = 'strain.'
            obj_vector_real = RealSolidStrainArray
            obj_vector_complex = ComplexSolidStrainArray
            obj_vector_random = RandomSolidStrainArray

            if self.element_type == 39: # CTETRA
                nnodes_expected = 5  # 1 centroid + 4 corner points
                result_name = prefix + 'ctetra_strain' + postfix
                element_name = 'CTETRA4'
            elif self.element_type == 67:  # CHEXA
                nnodes_expected = 9
                result_name = prefix + 'chexa_strain' + postfix
                element_name = 'CHEXA8'
            elif self.element_type == 68:  # CPENTA
                nnodes_expected = 7
                result_name = prefix + 'cpenta_strain' + postfix
                element_name = 'CPENTA6'
            elif self.element_type == 255:  # CPYRAM
                nnodes_expected = 6
                result_name = prefix + 'cpyram_strain' + postfix
                element_name = 'CPYRAM5'
            else:  # pragma: no cover
                raise RuntimeError(self.code_information())
                #msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
                #return self._not_implemented_or_skip(data, ndata, msg)

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)

        numwide_real = 4 + 21 * nnodes_expected
        numwide_imag = 4 + (17 - 4) * nnodes_expected
        numwide_random = 4 + (11 - 4) * nnodes_expected
        numwide_random2 = 18 + 14 * (nnodes_expected - 1)
        preline1 = '%s-%s' % (self.element_name, self.element_type)
        preline2 = ' ' * len(preline1)

        #print('numwide real=%s imag=%s random=%s' % (numwide_real, numwide_imag, numwide_random2))
        self._data_factor = nnodes_expected
        if self.format_code == 1 and self.num_wide == numwide_real:  # real
            ntotal = (16 + 84 * nnodes_expected) * self.factor
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return nelements * ntotal, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                itotali = obj.itotal + nelements
                itotal2 = obj.itotal + nelements * nnodes_expected
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    # (eid_device, cid, abcd, nnodes)
                    ints = frombuffer(data, dtype=self.idtype8).copy()
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

                    #print('%s-grid_device=%s' % (self.element_name, grid_device))
                    unused_grid_device2 = repeat(grid_device, nnodes_expected)
                    try:
                        obj.element_node[itotal:itotal2, 1] = grid_device
                    except ValueError:
                        msg = '%s; nnodes=%s\n' % (self.element_name, nnodes_expected)
                        msg += 'itotal=%s itotal2=%s\n' % (itotal, itotal2)
                        msg += 'grid_device.shape=%s; size=%s\n' % (str(grid_device.shape), grid_device.size)
                        #msg += 'nids=%s' % nids
                        raise ValueError(msg)
                    #self.log.debug(f'cids = {np.unique(cids)}')
                    obj.element_cid[itotal:itotali, 0] = eids
                    obj.element_cid[itotal:itotali, 1] = cids

                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, numwide_real)[:, 4:]
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
                if is_vectorized and self.use_vector:  # pragma: no cover
                    self.log.debug(f'vectorize CSolid real SORT{self.sort_method}')
                n = oes_csolid_real(self, data, obj,
                                    nelements, dt,
                                    element_name, nnodes_expected,
                                    preline1, preline2)

        elif self.format_code in [2, 3] and self.num_wide == numwide_imag:  # complex
            ntotal = numwide_imag * 4 * self.factor
            nelements = ndata // ntotal
            self.ntotal += nelements * nnodes_expected
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_complex)
            if auto_return:
                return nelements * ntotal, None, None

            obj = self.obj

            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * ntotal
                ielement = obj.ielement
                ielement2 = ielement + nelements
                itotal = obj.itotal
                itotal2 = itotal + nelements * nnodes_expected

                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, numwide_imag)
                floats1 = floats[:, 4:].reshape(nelements * nnodes_expected, 13).copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype8).reshape(nelements, numwide_imag)
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
                if is_vectorized and self.use_vector:  # pragma: no cover
                    self.log.debug(f'vectorize CSolid imag SORT{self.sort_method}')
                n = oes_csolid_complex(self, data, obj,
                                       nelements, # nnodes,
                                       element_name, nnodes_expected,
                                       is_magnitude_phase)

        elif self.format_code == 1 and self.num_wide == numwide_random: # random
            if not self.is_sort1:
                self.log.debug(f'support CSolid random SORT{self.sort_method}')
                return ndata, None, None

            ntotal = numwide_random * 4
            nelements = ndata // ntotal
            assert ndata % ntotal == 0, ndata
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_random)
            if auto_return:
                return nelements * ntotal, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and 0:  # pragma: no cover
                n = nelements * 4 * self.num_wide
                itotal = obj.ielement
                itotali = obj.itotal + nelements
                itotal2 = obj.itotal + nelements * nnodes_expected
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    # (eid_device, cid, abcd, nnodes)
                    ints = frombuffer(data, dtype=self.idtype).copy()
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

                    #print('%s-grid_device=%s' % (self.element_name, grid_device))
                    unused_grid_device2 = repeat(grid_device, nnodes_expected)
                    try:
                        obj.element_node[itotal:itotal2, 1] = grid_device
                    except ValueError:
                        msg = '%s; nnodes=%s\n' % (self.element_name, nnodes_expected)
                        msg += 'itotal=%s itotal2=%s\n' % (itotal, itotal2)
                        msg += 'grid_device.shape=%s; size=%s\n' % (str(grid_device.shape), grid_device.size)
                        #msg += 'nids=%s' % nids
                        raise ValueError(msg)
                    obj.element_cid[itotal:itotali, 0] = eids
                    obj.element_cid[itotal:itotali, 1] = cids

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, numwide_real)[:, 4:]
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
                if is_vectorized and self.use_vector and obj.itime == 0:  # pragma: no cover
                    self.log.debug(f'vectorize CSolid random SORT{self.sort_method}')
                n = oes_csolid_random(self, data, obj, nelements,
                                      element_name, nnodes_expected,
                                      preline1, preline2)

        elif self.format_code in [2, 3] and self.num_wide == numwide_random2:
            #raise RuntimeError(self.code_information())
            ## a = 18
            ## b = 14
            ## a + b * nnodes = numwide_random3
            ## a + b * 4 = 74  # CTETRA
            ## a + b * 6 = 102 # CPENTA
            ## a + b * 8 = 130 # CHEXA-67
            #msg = 'OES-CHEXA-random-numwide=%s numwide_real=%s numwide_imag=%s numwide_random=%s' % (
                #self.num_wide, numwide_real, numwide_imag, numwide_random)
            #return self._not_implemented_or_skip(data, ndata, msg)

            #print('numwide real=%s imag=%s random=%s' % (numwide_real, numwide_imag, numwide_random))
            unused_num_wide_random = 4 + nnodes_expected * (17 - 4)

            #print('random2=%s' % num_wide_random)
            #print(self.code_information())

            #if self.num_wide ==
            if self.read_mode == 1:
                return ndata, None, None
            return ndata, None, None
            #print('numwide=%s numwide_random=%s attempt2=%s subcase=%s' % (
                #self.num_wide, numwide_random, num_wide_random, self.isubcase))
            ##msg = self.code_information()
            #ntotal = 130
            #nelements = ndata // ntotal

            ## cid, coord_type, nactive_pnts,
            ##      nid, oxx, oyy, ozz, txy, tyz, txz
            #struct1 = Struct(self._endian + b'2i 4s')
            #struct2 = Struct(self._endian + b'i6f')
            #for i in range(nelements):
                #edata = data[n:n+12]
                #out = struct1.unpack(edata)
                #(eid_device, cid, abcd) = out
                #eid, dt = get_eid_dt_from_eid_device(
                    #eid_device, self.nonlinear_factor, self.sort_method)
                #if self.is_debug_file:
                    #self.binary_debug.write('%s - eid=%i; %s\n' % (preline1, eid, str(out)))
                #n += 12
                #for inode in range(nnodes_expected):  # nodes pts, +1 for centroid (???)
                    #out = struct2.unpack(data[n:n + 28]) # 4*7 = 28
                    #if self.is_debug_file:
                        #self.binary_debug.write('%s - %s\n' % (preline2, str(out)))
                    #(grid_device, sxx, syy, sz, txy, tyz, txz) = out
            #msg = 'OES-CHEXA-random-numwide=%s numwide_real=%s numwide_imag=%s numwide_random=%s' % (
                #self.num_wide, numwide_real, numwide_imag, numwide_random)
            #return self._not_implemented_or_skip(data, ndata, msg)
        elif self.format_code in [1, 2] and self.num_wide == 67:  # CHEXA
            #44 = 5 * 8 + 4  (TETRA)
            #52 = 6 * 8 + 4  (PYRAM)
            #60 = 7 * 8 + 4  (PENTA)
            #76 = 9 * 8 + 4  (HEXA)
            msg = 'skipping random CHEXA; numwide=67'
            #print(self.code_information())
            #asdf
            n = self._not_implemented_or_skip(data, ndata, msg)
            nelements = None
            ntotal = None
        elif self.format_code in [1, 2, 3] and self.num_wide in [60] and self.table_name in [b'OESXRMS1', b'OESXNO1']:  # CPENTA
            # bad
            #if self.read_mode == 2:
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
            n = self._not_implemented_or_skip(data, ndata, msg)
            nelements = None
            ntotal = None
        elif self.format_code in [1, 2, 3] and self.num_wide in [76] and self.table_name in [b'OESXRMS1', b'OESXNO1']:  # CHEXA
            # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1081x.op2
            msg = 'skipping random CHEXA; numwide=76'
            n = self._not_implemented_or_skip(data, ndata, msg)
            nelements = None
            ntotal = None
        #elif self.format_code in [2, 3] and self.num_wide == 76:  # imag
            # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1081x.op2
            # analysis_code = 5   Frequency
            # table_code    = 905 OESXNO1-OESXNO1C - Cumulative Root Mean Square output
            # format_code   = 3   Magnitude/Phase
            # result_type   = 3   Magnitude/Phase
            # sort_method   = 1
            # random_code   = 0
            # element_type  = 67  CHEXA
            #msg = self.code_information()
            #return self._not_implemented_or_skip(data, ndata, msg), None, None
        else:  # pragma: no cover
            raise RuntimeError(self.code_information() +
                               '\nnumwide real=%s imag=%s random=%s' % (
                                   numwide_real, numwide_imag, numwide_random2))
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
        n = 0
        if self.is_stress:
            #obj_vector_real = RealSolidStressArray
            #obj_vector_complex = ComplexSolidStressArray
            #obj_vector_random = RandomSolidStressArray
            word = 'stress'
            if self.element_type == 300:  # CHEXA
                nnodes_expected = 8
                result_name = prefix + 'chexa_stress' + postfix
                element_name = 'CHEXA8'
                # real=67
            elif self.element_type == 301:  # CPENTA
                nnodes_expected = 6
                result_name = prefix + 'cpenta_stress' + postfix
                element_name = 'CPENTA6'
            elif self.element_type == 302:  # CTETRA
                nnodes_expected = 4
                result_name = prefix + 'ctetra_stress' + postfix
                element_name = 'CTETRA4'
            elif self.element_type == 303:  # CPYRAM
                nnodes_expected = 5
                result_name = prefix + 'cpyram_stress' + postfix
                element_name = 'CPYRAM5'
            else:  # pragma: no cover
                raise RuntimeError(self.code_information())
        else:
            #obj_vector_real = RealSolidStrainArray
            #obj_vector_complex = ComplexSolidStrainArray
            #obj_vector_random = RandomSolidStrainArray
            word = 'strain'
            if self.element_type == 300:  # CHEXA
                nnodes_expected = 8
                result_name = prefix + 'chexa_strain' + postfix
                element_name = 'CHEXA8'
            elif self.element_type == 301:  # CPENTA
                nnodes_expected = 6
                result_name = prefix + 'cpenta_strain' + postfix
                element_name = 'CPENTA6'
            elif self.element_type == 302:  # CTETRA
                nnodes_expected = 4
                result_name = prefix + 'ctetra_strain' + postfix
                element_name = 'CTETRA4'
            elif self.element_type == 303:  # CPYRAM
                nnodes_expected = 5
                result_name = prefix + 'cpyram_strain' + postfix
                element_name = 'CPYRAM5'
            else:  # pragma: no cover
                raise RuntimeError(self.code_information())
                #msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
                #return self._not_implemented_or_skip(data, ndata, msg)

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)

        numwide_real = 3 + 8 * nnodes_expected
        #numwide_imag = 4 + (17 - 4) * nnodes_expected
        #numwide_random = 4 + (11 - 4) * nnodes_expected
        #numwide_random2 = 18 + 14 * (nnodes_expected - 1)
        preline1 = '%s-%s' % (self.element_name, self.element_type)
        preline2 = ' ' * len(preline1)

        #print('nnodes_expected =', nnodes_expected)
        #print('numwide real=%s imag=%s random=%s' % (numwide_real, numwide_imag, numwide_random2))
        self._data_factor = nnodes_expected
        if self.format_code == 1 and self.num_wide == numwide_real:  # real
            ntotal = 12 + 32 * nnodes_expected
            nelements = ndata // ntotal
            #auto_return, is_vectorized = self._create_oes_object4(
                #nelements, result_name, slot, obj_vector_real)
            auto_return = self.read_mode == 1
            is_vectorized = False
            if auto_return:
                return nelements * self.num_wide * 4, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1 and 0:  # pragma: no cover
                n = nelements * 4 * self.num_wide
                itotal = obj.ielement
                itotali = obj.itotal + nelements
                itotal2 = obj.itotal + nelements * nnodes_expected
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    # (eid_device, cid, abcd, nnodes)
                    ints = frombuffer(data, dtype=self.idtype).copy()
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

                    #print('%s-grid_device=%s' % (self.element_name, grid_device))
                    unused_grid_device2 = repeat(grid_device, nnodes_expected)
                    try:
                        obj.element_node[itotal:itotal2, 1] = grid_device
                    except ValueError:
                        msg = '%s; nnodes=%s\n' % (self.element_name, nnodes_expected)
                        msg += 'itotal=%s itotal2=%s\n' % (itotal, itotal2)
                        msg += 'grid_device.shape=%s; size=%s\n' % (str(grid_device.shape), grid_device.size)
                        #msg += 'nids=%s' % nids
                        raise ValueError(msg)
                    obj.element_cid[itotal:itotali, 0] = eids
                    obj.element_cid[itotal:itotali, 1] = cids

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, numwide_real)[:, 4:]
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
                #if is_vectorized and self.use_vector:  # pragma: no cover
                    #self.log.debug('vectorize CSolid real SORT%s' % self.sort_method)

            # 2 CID I Coordinate System
            # 3 CTYPE CHAR4 Grid or Gauss
            #
            # 4 GRID I Corner grid ID
            # 5 EX RS Strain in X
            # 6 EY RS Strain in Y
            # 7 EZ RS Strain in Z
            # 8 EXY RS Strain in XY
            # 9 EYZ RS Strain in YZ
            # 10 EZX RS Strain in ZX
            # 11 EVM RS Von Mises strain
            # Words 4 through 11 repeat nnodes times.
                struct1 = Struct(self._endian + self._analysis_code_fmt + b'i4s')
                struct2 = Struct(self._endian + b'i7f')
                if self.is_debug_file:
                    msg = '%s-%s nelements=%s nnodes=%s; C=[sxx, syy, szz, txy, tyz, txz, ovm,\n' % (
                        self.element_name, self.element_type, nelements, nnodes_expected)
                    self.binary_debug.write(msg)

                for unused_i in range(nelements):
                    edata = data[n:n+12]
                    out = struct1.unpack(edata)
                    (eid_device, cid, unused_abcd) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)

                    if self.is_debug_file:
                        self.binary_debug.write('%s - eid=%i; %s\n' % (preline1, eid, str(out)))

                    #assert nnodes < 21, 'print_block(data[n:n+16])'  #self.print_block(data[n:n+16])

                    n += 12
                    for inode in range(nnodes_expected):  # nodes pts, no centroid
                        out = struct2.unpack(data[n:n + 32]) # 4*8 = 32
                        if self.is_debug_file:
                            self.binary_debug.write('%s - %s\n' % (preline2, str(out)))
                        (grid_device, sxx, syy, szz, txy, tyz, txz, ovm) = out

                        if self.is_debug_file:
                            self.binary_debug.write('  eid=%s inode=%i; C=[%s]\n' % (
                                eid, grid_device, ', '.join(['%r' % di for di in out])))

                        #if grid_device == 0:
                            #grid = 'CENTER'
                        #else:
                            ##grid = (grid_device - device_code) // 10
                            #grid = grid_device

                        grid = grid_device
                        #a_cos = [a1, a2, a3]
                        #b_cos = [b1, b2, b3]
                        #c_cos = [c1, c2, c3]
                        if 0:
                            if inode == 0:
                                #  this is correct, but fails
                                #element_name = self.element_name + str(nnodes)
                                obj.add_eid_sort1(element_name, cid, dt, eid, grid,
                                                  sxx, syy, szz, txy, tyz, txz, ovm)
                            else:
                                obj.add_node_sort1(dt, eid, inode, grid,
                                                   sxx, syy, szz, txy, tyz, txz, ovm)
                        n += 32
            self.log.warning(f'skipping {self.table_name_str}: {self.element_name}-{self.element_type} {word} csolid2')
        else:  # pragma: no cover
            raise NotImplementedError(self.code_information())
        assert n == ntotal * nelements, f'n={n} ntotal={ntotal*nelements}'
        return n, nelements, ntotal

    def _oes_csolid_composite(self, data, ndata, dt, is_magnitude_phase: bool,
                              result_type: int, prefix: str, postfix: str) -> int:
        """
        reads stress/strain for element type:
         - 306 : CHEXALN
         - 307 : CPENTA
         #- 302 : CTETRA
         #- 303 : CPYRAM

        """
        n = 0
        if self.is_stress:
            #obj_vector_real = RealSolidStressArray
            #obj_vector_complex = ComplexSolidStressArray
            #obj_vector_random = RandomSolidStressArray
            word = 'stress'
            if self.element_type == 306:  # CHEXALN
                nedges = 4 # quad
                #nnodes_expected = 8
                result_name = prefix + 'chexa_stress' + postfix
                element_name = 'CHEXA8'
                # real=67
            elif self.element_type == 307:  # CPENTALN
                nedges = 3 # tri
                #nnodes_expected = 6
                result_name = prefix + 'cpenta_stress' + postfix
                element_name = 'CPENTA6'
            #elif self.element_type == 302:  # CTETRA
                #nnodes_expected = 4
                #result_name = prefix + 'ctetra_stress' + postfix
                #element_name = 'CTETRA4'
            #elif self.element_type == 303:  # CPYRAM
                #nnodes_expected = 5
                #result_name = prefix + 'cpyram_stress' + postfix
                #element_name = 'CPYRAM5'
            else:  # pragma: no cover
                raise RuntimeError(self.code_information())
        else:
            #obj_vector_real = RealSolidStrainArray
            #obj_vector_complex = ComplexSolidStrainArray
            #obj_vector_random = RandomSolidStrainArray
            word = 'strain'
            if self.element_type == 306:  # CHEXALN
                nedges = 4 # quad
                #nnodes_expected = 8
                result_name = prefix + 'chexa_strain' + postfix
                element_name = 'CHEXA8'
            elif self.element_type == 307:  # CPENTA
                nedges = 3 # tri
                #nnodes_expected = 6
                result_name = prefix + 'cpenta_strain' + postfix
                element_name = 'CPENTA6'
            #elif self.element_type == 302:  # CTETRA
                #nnodes_expected = 4
                #result_name = prefix + 'ctetra_strain' + postfix
                #element_name = 'CTETRA4'
            #elif self.element_type == 303:  # CPYRAM
                #nnodes_expected = 5
                #result_name = prefix + 'cpyram_strain' + postfix
                #element_name = 'CPYRAM5'
            else:  # pragma: no cover
                raise NotImplementedError(self.code_information())
                #msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
                #return self._not_implemented_or_skip(data, ndata, msg)

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)

        numwide_real = 3 + 8 * nedges # 3 + 8*4 = 35
        numwide_imag = 3 + 7 * 14 # 3 + 7 * 14 = 101
        #print(self.code_information())
        #print(f'{self.element_name} numwide_real={numwide_real} numwide_imag={numwide_imag} -> {self.num_wide}')
        #numwide_real = 3 + 8 * nnodes_expected
        #numwide_imag = 4 + (17 - 4) * nnodes_expected
        #numwide_random = 4 + (11 - 4) * nnodes_expected
        #numwide_random2 = 18 + 14 * (nnodes_expected - 1)
        preline1 = '%s-%s' % (self.element_name, self.element_type)
        preline2 = ' ' * len(preline1)

        #print('nnodes_expected =', nnodes_expected)
        #print('numwide real=%s imag=%s random=%s' % (numwide_real, numwide_imag, numwide_random2))
        self._data_factor = nedges
        if result_type == 0 and self.num_wide == numwide_real:  # real
            ntotal = 12 + 32 * nedges
            nelements = ndata // ntotal
            #auto_return, is_vectorized = self._create_oes_object4(
                #nelements, result_name, slot, obj_vector_real)
            auto_return = self.read_mode == 1
            is_vectorized = False
            if auto_return:
                return nelements * self.num_wide * 4, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1 and 0:  # pragma: no cover
                n = nelements * 4 * self.num_wide
                itotal = obj.ielement
                itotali = obj.itotal + nelements
                itotal2 = obj.itotal + nelements * nnodes_expected
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    # (eid_device, cid, abcd, nnodes)
                    ints = frombuffer(data, dtype=self.idtype).copy()
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

                    #print('%s-grid_device=%s' % (self.element_name, grid_device))
                    unused_grid_device2 = repeat(grid_device, nnodes_expected)
                    try:
                        obj.element_node[itotal:itotal2, 1] = grid_device
                    except ValueError:
                        msg = '%s; nnodes=%s\n' % (self.element_name, nnodes_expected)
                        msg += 'itotal=%s itotal2=%s\n' % (itotal, itotal2)
                        msg += 'grid_device.shape=%s; size=%s\n' % (str(grid_device.shape), grid_device.size)
                        #msg += 'nids=%s' % nids
                        raise ValueError(msg)
                    obj.element_cid[itotal:itotali, 0] = eids
                    obj.element_cid[itotal:itotali, 1] = cids

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, numwide_real)[:, 4:]
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
                #if is_vectorized and self.use_vector:  # pragma: no cover
                    #self.log.debug('vectorize CSolid real SORT%s' % self.sort_method)
                n = oes_csolid_composite_real(self, data, obj,
                                              nelements, nedges,
                                              element_name, preline1, preline2, dt)
            self.log.warning(f'skipping {self.table_name_str}: {self.element_name}-{self.element_type} {word} csolid composite')

        elif result_type == 1 and self.num_wide == numwide_imag:  # complex
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
            if self.read_mode == 1:
                return ndata, None, None
            self.show_data(data[n:n+4*self.num_wide])
            aaa
        else:  # pragma: no cover
            raise NotImplementedError(self.code_information())
        assert n == ntotal * nelements, f'n={n} ntotal={ntotal*nelements}'
        return n, nelements, ntotal

    def _oes_csolid_linear_hyperelastic_cosine(self, data, ndata, dt, unused_is_magnitude_phase,
                                               result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 140 :CHEXAFD

        """
        n = 0
        if self.is_stress:
            #obj_vector_real = RealSolidStressArray
            #obj_vector_complex = ComplexSolidStressArray
            #obj_vector_random = RandomSolidStressArray
            word = 'stress'
            prefix = word + '.'
            if self.element_type == 140:  # CHEXA
                nnodes_expected = 8
                result_name = prefix + 'chexa_stress' + postfix
                #element_name = 'CHEXA8'
                # real=122
            #elif self.element_type == 160:  # CPENTA
                #nnodes_expected = 6
                #result_name = prefix + 'cpenta_stress' + postfix
                #element_name = 'CPENTA6'
            #elif self.element_type == 165:  # CPENTA
                #nnodes_expected = 21
                #result_name = prefix + 'cpenta_stress' + postfix
                #element_name = 'CPENTA6'

            #elif self.element_type == 161:  # CTETRA
                #nnodes_expected = 1
                #result_name = prefix + 'ctetra_stress' + postfix
                #element_name = 'CTETRA4'
            #elif self.element_type == 166:  # CTETRA
                #nnodes_expected = 5
                #result_name = prefix + 'ctetra_stress' + postfix
                #element_name = 'CTETRA4'
            #elif self.element_type == 303:  # CPYRAM
                #nnodes_expected = 5
                #result_name = prefix + 'cpyram_stress' + postfix
                #element_name = 'CPYRAM5'
            else:  # pragma: no cover
                raise RuntimeError(self.code_information())
        else:
            #obj_vector_real = RealSolidStrainArray
            #obj_vector_complex = ComplexSolidStrainArray
            #obj_vector_random = RandomSolidStrainArray
            word = 'strain'
            #if self.element_type == 202:  # CHEXA
                #nnodes_expected = 8
                #result_name = prefix + 'chexa_strain' + postfix
                #element_name = 'CHEXA8'
            #elif self.element_type == 301:  # CPENTA
                #nnodes_expected = 6
                #result_name = prefix + 'cpenta_strain' + postfix
                #element_name = 'CPENTA6'
            #elif self.element_type == 302:  # CTETRA
                #nnodes_expected = 4
                #result_name = prefix + 'ctetra_strain' + postfix
                #element_name = 'CTETRA4'
            #elif self.element_type == 303:  # CPYRAM
                #nnodes_expected = 5
                #result_name = prefix + 'cpyram_strain' + postfix
                #element_name = 'CPYRAM5'
            #else:  # pragma: no cover
            raise RuntimeError(self.code_information())
                #msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
                #return self._not_implemented_or_skip(data, ndata, msg)

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)

        numwide_real = 2 + 20 * nnodes_expected
        #numwide_real = 162  # CHEXA
        #print(self.num_wide, numwide_real)
        assert numwide_real == self.num_wide, numwide_real
        #numwide_imag = 4 + (17 - 4) * nnodes_expected
        #numwide_random = 4 + (11 - 4) * nnodes_expected
        #numwide_random2 = 18 + 14 * (nnodes_expected - 1)
        preline1 = '%s-%s' % (self.element_name, self.element_type)
        preline2 = ' ' * len(preline1)

        #print('nnodes_expected =', nnodes_expected)
        #print('numwide real=%s imag=%s random=%s' % (numwide_real, numwide_imag, numwide_random2))
        self._data_factor = nnodes_expected

        if self.format_code == 1 and self.num_wide == numwide_real:  # real
            ntotal = 8 + 80 * nnodes_expected
            #ntotal = numwide_real * 4
            nelements = ndata // ntotal
            assert ndata % ntotal == 0
            #auto_return, is_vectorized = self._create_oes_object4(
                #nelements, result_name, slot, obj_vector_real)
            auto_return = self.read_mode == 1
            is_vectorized = False
            if auto_return:
                return nelements * self.num_wide * 4, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1 and 0:  # pragma: no cover
                n = nelements * 4 * self.num_wide
                itotal = obj.ielement
                itotali = obj.itotal + nelements
                itotal2 = obj.itotal + nelements * nnodes_expected
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    # (eid_device, cid, abcd, nnodes)
                    ints = frombuffer(data, dtype=self.idtype).copy()
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

                    #print('%s-grid_device=%s' % (self.element_name, grid_device))
                    unused_grid_device2 = repeat(grid_device, nnodes_expected)
                    try:
                        obj.element_node[itotal:itotal2, 1] = grid_device
                    except ValueError:
                        msg = '%s; nnodes=%s\n' % (self.element_name, nnodes_expected)
                        msg += 'itotal=%s itotal2=%s\n' % (itotal, itotal2)
                        msg += 'grid_device.shape=%s; size=%s\n' % (str(grid_device.shape), grid_device.size)
                        #msg += 'nids=%s' % nids
                        raise ValueError(msg)
                    obj.element_cid[itotal:itotali, 0] = eids
                    obj.element_cid[itotal:itotali, 1] = cids

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, numwide_real)[:, 4:]
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
                #if is_vectorized and self.use_vector:  # pragma: no cover
                    #self.log.debug('vectorize CSolid real SORT%s' % self.sort_method)

                # ELTYPE =140 Hyperelastic 8-noded hexahedron element linear format
                # (HEXAFD)
                # 2 TYPE CHAR4 Gaus
                #
                # 3 ID I
                # 4 SX RS
                # 5 SXY RS
                # 6 PA RS
                # 7 AX RS
                # 8 AY RS
                # 9 AZ RS
                # 10 PRESSURE RS
                # 11 SY RS
                # 12 SYZ RS
                # 13 PB RS
                # 14 BX RS
                # 15 BY RS
                # 16 BZ RS
                # 17 SZ RS
                # 18 SZX RS
                # 19 PC RS
                # 20 CX RS
                # 21 CY RS
                # 22 CZ RS
                # Words 3 through 22 repeat 008 times
                struct1 = Struct(self._endian + self._analysis_code_fmt + b'4s')
                struct2 = Struct(self._endian + b'i19f')
                if self.is_debug_file:
                    msg = (
                        f'{self.element_name}-{self.element_type} nelements={nelements} '
                        f'nnodes={nnodes_expected}; '
                        'C=[oxx, oxy, pa, ax, ay, az, pressure, '
                        'oyy, oyz, pb, bx, by, bz, '
                        'ozz, oxz, pc, cx, cy, cz]\n')
                    self.binary_debug.write(msg)

                for unused_i in range(nelements):
                    edata = data[n:n+8]
                    out = struct1.unpack(edata)
                    (eid_device, grid_gauss, ) = out
                    #print(out)
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)

                    if self.is_debug_file:
                        self.binary_debug.write('%s - eid=%i; %s\n' % (preline1, eid, str(out)))
                    #assert nnodes < 21, 'print_block(data[n:n+16])'  #self.print_block(data[n:n+16])

                    n += 8
                    for inode in range(nnodes_expected):  # nodes pts, no centroid
                        out = struct2.unpack(data[n:n + 80]) # 4*20 = 80
                        if self.is_debug_file:
                            self.binary_debug.write('%s - %s\n' % (preline2, str(out)))
                        # nid, oxx, oxy, pa, ax, ay, az, pressure,
                        #      oyy, oyz, pb, bx, by, bz,
                        #      ozz, oxz, pc, cx, cy, cz
                        (grid_device,
                         oxx, oxy, pa, ax, ay, az, pressure,
                         oyy, oyz, pb, bx, by, bz,
                         ozz, oxz, pc, cx, cy, cz) = out
                        #print(out)

                        if self.is_debug_file:
                            self.binary_debug.write('  eid=%s inode=%i; C=[%s]\n' % (
                                eid, grid_device, ', '.join(['%r' % di for di in out])))

                        #if 0:
                            #if inode == 0:
                                #  this is correct, but fails
                                #element_name = self.element_name + str(nnodes)
                                #obj.add_eid_sort1(element_name, cid, dt, eid, grid,
                                                  #sxx, syy, szz, txy, tyz, txz, ovm)
                            #else:
                                #obj.add_node_sort1(dt, eid, inode, grid,
                                                   #sxx, syy, szz, txy, tyz, txz, ovm)
                        n += 80
            self.log.warning(f'skipping {self.table_name_str}: {self.element_name}-{self.element_type} linear hyperelastic cosine {word}')
            return n, None, None
        else:  # pragma: no cover
            raise RuntimeError(self.code_information())
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
        n = 0
        if self.is_stress:
            #obj_vector_real = RealSolidStressArray
            #obj_vector_complex = ComplexSolidStressArray
            #obj_vector_random = RandomSolidStressArray
            word = 'stress'
            if self.element_type == 163:  # CHEXA
                nnodes_expected = 27
                result_name = prefix + 'chexa_stress' + postfix
                #element_name = 'CHEXA8'
                # real=122
            elif self.element_type == 160:  # CPENTA
                nnodes_expected = 6
                result_name = prefix + 'cpenta_stress' + postfix
                #element_name = 'CPENTA6'
            elif self.element_type == 165:  # CPENTA
                nnodes_expected = 21
                result_name = prefix + 'cpenta_stress' + postfix
                #element_name = 'CPENTA6'

            elif self.element_type == 161:  # CTETRA
                nnodes_expected = 1
                result_name = prefix + 'ctetra_stress' + postfix
                #element_name = 'CTETRA4'
            elif self.element_type == 166:  # CTETRA
                nnodes_expected = 5
                result_name = prefix + 'ctetra_stress' + postfix
                #element_name = 'CTETRA4'
            #elif self.element_type == 303:  # CPYRAM
                #nnodes_expected = 5
                #result_name = prefix + 'cpyram_stress' + postfix
                #element_name = 'CPYRAM5'
            else:  # pragma: no cover
                raise RuntimeError(self.code_information())
        else:
            #obj_vector_real = RealSolidStrainArray
            #obj_vector_complex = ComplexSolidStrainArray
            #obj_vector_random = RandomSolidStrainArray
            word = 'strain'
            #if self.element_type == 202:  # CHEXA
                #nnodes_expected = 8
                #result_name = prefix + 'chexa_strain' + postfix
                #element_name = 'CHEXA8'
            #elif self.element_type == 301:  # CPENTA
                #nnodes_expected = 6
                #result_name = prefix + 'cpenta_strain' + postfix
                #element_name = 'CPENTA6'
            #elif self.element_type == 302:  # CTETRA
                #nnodes_expected = 4
                #result_name = prefix + 'ctetra_strain' + postfix
                #element_name = 'CTETRA4'
            #elif self.element_type == 303:  # CPYRAM
                #nnodes_expected = 5
                #result_name = prefix + 'cpyram_strain' + postfix
                #element_name = 'CPYRAM5'
            #else:  # pragma: no cover
            raise RuntimeError(self.code_information())
                #msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
                #return self._not_implemented_or_skip(data, ndata, msg)

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)

        numwide_real = 2 + 20 * nnodes_expected
        #numwide_real = 122  # CHEXA
        #print(self.num_wide, numwide_real)
        assert numwide_real == self.num_wide, numwide_real
        numwide_imag = 4 + (17 - 4) * nnodes_expected
        #numwide_random = 4 + (11 - 4) * nnodes_expected
        numwide_random2 = 18 + 14 * (nnodes_expected - 1)
        preline1 = '%s-%s' % (self.element_name, self.element_type)
        preline2 = ' ' * len(preline1)

        #print('nnodes_expected =', nnodes_expected)
        #print('numwide real=%s imag=%s random=%s' % (numwide_real, numwide_imag, numwide_random2))
        self._data_factor = nnodes_expected

        if self.format_code == 1 and self.num_wide == numwide_real:  # real
            ntotal = 8 + 80 * nnodes_expected
            #ntotal = numwide_real * 4
            nelements = ndata // ntotal
            assert ndata % ntotal == 0
            #auto_return, is_vectorized = self._create_oes_object4(
                #nelements, result_name, slot, obj_vector_real)
            auto_return = self.read_mode == 1
            is_vectorized = False
            if auto_return:
                return nelements * self.num_wide * 4, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1 and 0:  # pragma: no cover
                n = nelements * 4 * self.num_wide
                itotal = obj.ielement
                itotali = obj.itotal + nelements
                itotal2 = obj.itotal + nelements * nnodes_expected
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    # (eid_device, cid, abcd, nnodes)
                    ints = frombuffer(data, dtype=self.idtype).copy()
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

                    #print('%s-grid_device=%s' % (self.element_name, grid_device))
                    unused_grid_device2 = repeat(grid_device, nnodes_expected)
                    try:
                        obj.element_node[itotal:itotal2, 1] = grid_device
                    except ValueError:
                        msg = '%s; nnodes=%s\n' % (self.element_name, nnodes_expected)
                        msg += 'itotal=%s itotal2=%s\n' % (itotal, itotal2)
                        msg += 'grid_device.shape=%s; size=%s\n' % (str(grid_device.shape), grid_device.size)
                        #msg += 'nids=%s' % nids
                        raise ValueError(msg)
                    obj.element_cid[itotal:itotali, 0] = eids
                    obj.element_cid[itotal:itotali, 1] = cids

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, numwide_real)[:, 4:]
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
                #if is_vectorized and self.use_vector:  # pragma: no cover
                    #self.log.debug('vectorize CSolid real SORT%s' % self.sort_method)

                #ELTYPE =163 Linear form for hyperelastic 20 node HEXAFD
                #2 TYPE CHAR4 Gaus
                #
                #3 ID I
                #4 SX RS
                #5 SXY RS
                #6 PA RS
                #7 AX RS
                #8 AY RS
                #9 AZ RS
                #10 PRESSURE RS
                #11 SY RS
                #12 SYZ RS
                #13 PB RS
                #14 BX RS
                #15 BY RS
                #16 BZ RS
                #17 SZ RS
                #18 SZX RS
                #19 PC RS
                #20 CX RS
                #21 CY RS
                #22 CZ RS
                #Words 3 through 22 repeat 027 times
                struct1 = Struct(self._endian + self._analysis_code_fmt + b'4s')
                struct2 = Struct(self._endian + b'i19f')
                if self.is_debug_file:
                    #msg = (
                        #f'{self.element_name}-{self.element_type} nelements={nelements} '
                        #f'nnodes={nnodes_expected}; '
                        #'C=[sxx, syy, szz, txy, tyz, txz, pressure, '
                        #'evol, exx, eyy, ezz, exy, eyz, exz]\n')
                    self.binary_debug.write(msg)

                for unused_i in range(nelements):
                    edata = data[n:n+8]
                    out = struct1.unpack(edata)
                    (eid_device, unused_abcd, ) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)

                    if self.is_debug_file:
                        self.binary_debug.write('%s - eid=%i; %s\n' % (preline1, eid, str(out)))
                    #assert nnodes < 21, 'print_block(data[n:n+16])'  #self.print_block(data[n:n+16])

                    n += 8
                    for unused_inode in range(nnodes_expected):  # nodes pts, no centroid
                        out = struct2.unpack(data[n:n + 80]) # 4*20 = 80
                        if self.is_debug_file:
                            self.binary_debug.write('%s - %s\n' % (preline2, str(out)))
                        #(grid_device, sxx, syy, szz, txy, tyz, txz, pressure,
                         #evol, exx, eyy, ezz, exy, eyz, exz) = out
                        #print(out)

                        if self.is_debug_file:
                            self.binary_debug.write('  eid=%s inode=%i; C=[%s]\n' % (
                                eid, grid_device, ', '.join(['%r' % di for di in out])))

                        #grid = grid_device
                        #if 0:
                            #if inode == 0:
                                ##  this is correct, but fails
                                ##element_name = self.element_name + str(nnodes)
                                #obj.add_eid_sort1(element_name, cid, dt, eid, grid,
                                                  #sxx, syy, szz, txy, tyz, txz, ovm)
                            #else:
                                #obj.add_node_sort1(dt, eid, inode, grid,
                                                   #sxx, syy, szz, txy, tyz, txz, ovm)
                        n += 80
            self.log.warning(f'skipping {self.table_name_str}: {self.element_name}-{self.element_type} linear hyperelastic {word}')
        else:  # pragma: no cover
            raise RuntimeError(self.code_information() +
                               '\nnumwide real=%s imag=%s random=%s' % (
                                   numwide_real, numwide_imag, numwide_random2))
        assert n == ntotal * nelements, f'n={n} ntotal={ntotal*nelements}'
        return n, nelements, ntotal

    def _oes_csolid_nonlinear_hyperelastic(self, data, ndata, dt, is_magnitude_phase,
                                           result_type, prefix, postfix):
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
        n = 0
        if self.is_stress:
            #obj_vector_real = RealSolidStressArray
            #obj_vector_complex = ComplexSolidStressArray
            #obj_vector_random = RandomSolidStressArray
            word = 'stress'
            if self.element_type in [202, 218]:  # CHEXA
                nnodes_expected = 8
                result_name = prefix + 'chexa_stress_strain' + postfix
                element_name = 'CHEXA8'
                # real=122
            elif self.element_type in [204, 220]:  # CPENTA
                nnodes_expected = 6
                result_name = prefix + 'cpenta_stress_strain' + postfix
                element_name = 'CPENTA6'
            elif self.element_type in [216, 221]:  # CTETRA
                nnodes_expected = 4
                result_name = prefix + 'ctetra_stress_strain' + postfix
                element_name = 'CTETRA4'
            #elif self.element_type == 303:  # CPYRAM
                #nnodes_expected = 5
                #result_name = prefix + 'cpyram_stress' + postfix
                #element_name = 'CPYRAM5'
            else:  # pragma: no cover
                raise RuntimeError(self.code_information())
        else:
            #obj_vector_real = RealSolidStrainArray
            #obj_vector_complex = ComplexSolidStrainArray
            #obj_vector_random = RandomSolidStrainArray
            word = 'strain'
            #if self.element_type == 202:  # CHEXA
                #nnodes_expected = 8
                #result_name = prefix + 'chexa_strain' + postfix
                #element_name = 'CHEXA8'
            #elif self.element_type == 301:  # CPENTA
                #nnodes_expected = 6
                #result_name = prefix + 'cpenta_strain' + postfix
                #element_name = 'CPENTA6'
            #elif self.element_type == 302:  # CTETRA
                #nnodes_expected = 4
                #result_name = prefix + 'ctetra_strain' + postfix
                #element_name = 'CTETRA4'
            #elif self.element_type == 303:  # CPYRAM
                #nnodes_expected = 5
                #result_name = prefix + 'cpyram_strain' + postfix
                #element_name = 'CPYRAM5'
            #else:  # pragma: no cover
            raise RuntimeError(self.code_information())
                #msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
                #return self._not_implemented_or_skip(data, ndata, msg)

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)

        numwide_real = 2 + 15 * nnodes_expected
        #numwide_real = 122  # CHEXA
        assert numwide_real == self.num_wide, numwide_real
        numwide_imag = 4 + (17 - 4) * nnodes_expected
        #numwide_random = 4 + (11 - 4) * nnodes_expected
        numwide_random2 = 18 + 14 * (nnodes_expected - 1)
        preline1 = '%s-%s' % (self.element_name, self.element_type)
        preline2 = ' ' * len(preline1)

        #print('nnodes_expected =', nnodes_expected)
        #print('numwide real=%s imag=%s random=%s' % (numwide_real, numwide_imag, numwide_random2))
        self._data_factor = nnodes_expected

        if result_type == 0 and self.num_wide == numwide_real:  # real
            ntotal = 8 + 60 * nnodes_expected
            #ntotal = numwide_real * 4
            nelements = ndata // ntotal
            assert ndata % ntotal == 0
            #auto_return, is_vectorized = self._create_oes_object4(
                #nelements, result_name, slot, obj_vector_real)
            auto_return = self.read_mode == 1
            is_vectorized = False
            if auto_return:
                return nelements * self.num_wide * 4, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1 and 0:  # pragma: no cover
                n = nelements * 4 * self.num_wide
                itotal = obj.ielement
                itotali = obj.itotal + nelements
                itotal2 = obj.itotal + nelements * nnodes_expected
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    # (eid_device, cid, abcd, nnodes)
                    ints = frombuffer(data, dtype=self.idtype).copy()
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

                    #print('%s-grid_device=%s' % (self.element_name, grid_device))
                    unused_grid_device2 = repeat(grid_device, nnodes_expected)
                    try:
                        obj.element_node[itotal:itotal2, 1] = grid_device
                    except ValueError:
                        msg = '%s; nnodes=%s\n' % (self.element_name, nnodes_expected)
                        msg += 'itotal=%s itotal2=%s\n' % (itotal, itotal2)
                        msg += 'grid_device.shape=%s; size=%s\n' % (str(grid_device.shape), grid_device.size)
                        #msg += 'nids=%s' % nids
                        raise ValueError(msg)
                    obj.element_cid[itotal:itotali, 0] = eids
                    obj.element_cid[itotal:itotali, 1] = cids

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, numwide_real)[:, 4:]
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
                #if is_vectorized and self.use_vector:  # pragma: no cover
                    #self.log.debug('vectorize CSolid real SORT%s' % self.sort_method)

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
                struct1 = Struct(self._endian + self._analysis_code_fmt + b'4s')
                struct2 = Struct(self._endian + b'i14f')
                if self.is_debug_file:
                    msg = (
                        f'{self.element_name}-{self.element_type} nelements={nelements} '
                        f'nnodes={nnodes_expected}; '
                        'C=[sxx, syy, szz, txy, tyz, txz, pressure, '
                        'evol, exx, eyy, ezz, exy, eyz, exz]\n')
                    self.binary_debug.write(msg)

                for unused_i in range(nelements):
                    edata = data[n:n+8]
                    out = struct1.unpack(edata)
                    (eid_device, unused_abcd, ) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)

                    if self.is_debug_file:
                        self.binary_debug.write('%s - eid=%i; %s\n' % (preline1, eid, str(out)))
                    #assert nnodes < 21, 'print_block(data[n:n+16])'  #self.print_block(data[n:n+16])

                    n += 8
                    for inode in range(nnodes_expected):  # nodes pts, no centroid
                        out = struct2.unpack(data[n:n + 60]) # 4*15 = 60
                        if self.is_debug_file:
                            self.binary_debug.write('%s - %s\n' % (preline2, str(out)))
                        (grid_device, sxx, syy, szz, txy, tyz, txz, pressure,
                         evol, exx, eyy, ezz, exy, eyz, exz) = out
                        #print(out)

                        if self.is_debug_file:
                            self.binary_debug.write('  eid=%s inode=%i; C=[%s]\n' % (
                                eid, grid_device, ', '.join(['%r' % di for di in out])))

                        grid = grid_device
                        if 0:
                            if inode == 0:
                                #  this is correct, but fails
                                #element_name = self.element_name + str(nnodes)
                                obj.add_eid_sort1(element_name, cid, dt, eid, grid,
                                                  sxx, syy, szz, txy, tyz, txz, ovm)
                            else:
                                obj.add_node_sort1(dt, eid, inode, grid,
                                                   sxx, syy, szz, txy, tyz, txz, ovm)
                        n += 60
            self.log.warning(f'skipping {self.table_name_str}: {self.element_name}-{self.element_type} nonlinear hyperelastic {word}')
            return n, None, None
        else:  # pragma: no cover
            raise RuntimeError(self.code_information() +
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
        #real
        #85:  2 + (18 - 2) * 5,  # Nonlinear CTETRA
        #256: 4 + (18 - 2) * 6,  # Nonlinear CHEXA -> ???

        # random
        #91:  4 + (25 - 4) * 7,  # Nonlinear CPENTA
        #93:  4 + (25 - 4) * 9,  # Nonlinear CHEXA -> 584 (can cause a crash)
        #256: 4 + (25 - 4) * 6,  # Nonlinear CHEXA -> ???

        # the nodes are nnodes + 1
        if self.element_type == 85:
            etype = 'CTETRANL'
            nnodes = 5
            result_name = prefix + 'ctetra_stress_strain' + postfix
        elif self.element_type == 91:
            etype = 'CPENTANL'
            nnodes = 7
            result_name = prefix + 'cpenta_stress_strain' + postfix
        elif self.element_type == 93:
            etype = 'CHEXANL'
            nnodes = 9
            result_name = prefix + 'chexa_stress_strain' + postfix
        elif self.element_type == 256:
            etype = 'CPYRAMNL'
            nnodes = 6
            result_name = prefix + 'chexa_stress_strain' + postfix

        else:  # pragma: no cover
            raise RuntimeError(self.code_information())

        numwide_real = 4 + (25 - 4) * nnodes  # real???
        numwide_random = 2 + (18 - 2) * nnodes  # imag???
        #self.log.debug("format_code=%s numwide=%s numwide_real=%s numwide_random=%s" % (
            #self.format_code, self.num_wide, numwide_real, numwide_random))

        #numwide_real = 0
        #numwide_imag = 2 + 16 * nnodes
        #ntotal = 8 + 64 * nnodes

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)

        slot = self.get_result(result_name)

        if self.format_code == 1 and self.num_wide == numwide_real:
            #if self.read_mode == 1:
                #return ndata, None, None

            ntotal = numwide_real * 4 * self.factor
            #if self.is_stress:
                #self.create_transient_object(self.nonlinearPlateStress, NonlinearSolid)
            #else:
                #self.create_transient_object(self.nonlinearPlateStrain, NonlinearSolid)
            #self.handle_results_buffer(self.OES_CQUAD4NL_90, resultName, name)
            raise RuntimeError('OES_CSOLIDNL_90')
        elif self.format_code == 1 and self.num_wide == numwide_random:  # random
            # 82 : CTETRA_NL (etype=85)
            # 146 : CHEXA_NL (etype=93)
            #raise RuntimeError(self.code_information())
        #elif self.format_code in [2, 3] and self.num_wide == numwide_imag:  # imag

            ntotal = numwide_random * 4 * self.factor

            nelements = ndata // ntotal
            self.ntotal += nelements * nnodes
            #print(self.read_mode, RealNonlinearSolidArray)
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, RealNonlinearSolidArray)
            if auto_return:
                self._data_factor = nnodes
                return nelements * ntotal, None, None

            n = 0
            s1 = Struct(self._endian + self._analysis_code_fmt + b'4s')
            s2 = Struct(self._endian + b'i15f')
            nelements = ndata // ntotal
            obj = self.obj
            for unused_i in range(nelements):  # 2+16*9 = 146 -> 146*4 = 584
                edata = data[n:n+8]
                n += 8

                out = s1.unpack(edata)
                if self.is_debug_file:
                    self.binary_debug.write('%s-%s - %s\n' % (etype, self.element_type, str(out)))

                (eid_device, unused_ctype) = out
                eid, dt = get_eid_dt_from_eid_device(
                    eid_device, self.nonlinear_factor, self.sort_method)
                #print('%s-%s -eid=%s dt=%s %s\n' % (etype, self.element_type, eid, dt, str(out)))

                for unused_j in range(nnodes):
                    edata = data[n:n+64]
                    n += 64
                    out = s2.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('%s-%sB - %s\n' % (etype, self.element_type, str(out)))
                    #print('%s-%sB - %s\n' % (etype, self.element_type, str(out)))
                    assert len(out) == 16

                    (grid,
                     sx, sy, sz, sxy, syz, sxz, se, eps, ecs,
                     ex, ey, ez, exy, eyz, exz) = out
                    obj.add_sort1(dt, eid, grid,
                                  sx, sy, sz, sxy, syz, sxz, se, eps, ecs,
                                  ex, ey, ez, exy, eyz, exz)

        else:  # pragma: no cover
            #msg = self.code_information()
            msg = "format_code=%s numwide=%s numwide_real=%s numwide_random=%s\n" % (
                self.format_code, self.num_wide, numwide_real, numwide_random)
            #return self._not_implemented_or_skip(data, ndata, msg)
            raise RuntimeError(msg + self.code_information())
        return n, nelements, ntotal

    def _oes_cquad4_33(self, data, ndata: int, dt, is_magnitude_phase: bool,
                       result_type: int, prefix: str, postfix: str) -> Tuple[int, Any, Any]:
        """
        reads stress/strain for element type:
         - 33 : CQUAD4-centroidal
         - 228 : CQUADR-centroidal

        """
        n = 0
        factor = self.factor
        size = self.size
        #print('_oes_cquad4_33')
        if self.is_stress:
            obj_vector_real = RealPlateStressArray
            obj_vector_complex = ComplexPlateStressArray
            if self.element_type == 33:
                result_name = prefix + 'cquad4_stress' + postfix
            elif self.element_type == 228:
                result_name = prefix + 'cquadr_stress' + postfix
                assert self.num_wide in [17, 15], self.code_information()
            else:
                raise NotImplementedError(self.code_information())
        else:
            obj_vector_real = RealPlateStrainArray
            obj_vector_complex = ComplexPlateStrainArray
            if self.element_type == 33:
                result_name = prefix + 'cquad4_strain' + postfix
            elif self.element_type == 228:
                result_name = prefix + 'cquadr_strain' + postfix
                assert self.num_wide in [17, 15], self.code_information()
            else:
                raise NotImplementedError(self.code_information())

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)

        numwide_real = 17
        sort_method = self.sort_method
        if result_type == 0 and self.num_wide == 17:  # real
            ntotal = 68 * factor  # 4*17
            nelements = ndata // ntotal
            nlayers = nelements * 2  # 2 layers per node
            #self.log.info(f'CQUAD4-33: len(data)={ndata} numwide={self.num_wide} nelements={nelements} nlayers={nlayers}')

            auto_return, is_vectorized = self._create_oes_object4(
                nlayers, result_name, slot, obj_vector_real)
            if auto_return:
                self._data_factor = 2  # number of "layers" for an element
                return nelements * ntotal, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and sort_method == 1:
                nnodes_expected = 2
                n = nelements * ntotal
                ielement = obj.ielement
                ielement2 = ielement + nelements
                itotal = obj.itotal
                itotal2 = itotal + nelements * nnodes_expected
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype8)
                    ints1 = ints.reshape(nelements, numwide_real)
                    eids = ints1[:, 0] // 10
                    eids = np.vstack([eids, eids]).T.ravel()
                    assert eids.min() > 0, eids.min()
                    obj.element_node[itotal:itotal2, 0] = eids

                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, numwide_real)[:, 1:]

                #fd, sx, sy, txy, angle, major, minor, max_shear
                floats1 = floats.reshape(nelements * nnodes_expected, 8)
                obj.data[obj.itime, itotal:itotal2, :] = floats1.copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if is_vectorized and self.use_vector:  # pragma: no cover
                    self.log.debug(f'vectorize centroidal quad: {self.element_name}-{self.element_type} real '
                                   f'SORT{sort_method}')
                n = oes_quad4_33_real_17(self, data, obj, ntotal, nelements, dt)
            if self.is_sort1:
                assert obj.element_node[:, 0].min() > 0, obj.element_node[:, 0].shape

        elif result_type == 1 and self.num_wide == 15:  # imag
            #self.to_nx(f' because CQUAD4-33 (numwide=15) was found')
            #nnodes = 0  # centroid + 4 corner points
            ntotal = self.num_wide * size
            #self.log.info(f'CQUAD4-33: len(data)={ndata} numwide={self.num_wide} nelements={nelements} nlayers={nlayers}')

            nelements = ndata // ntotal
            nlayers = nelements * 2
            auto_return, is_vectorized = self._create_oes_object4(
                nlayers, result_name, slot, obj_vector_complex)
            if auto_return:
                self._data_factor = 2
                return nelements * ntotal, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and sort_method == 1:
                n = nelements * ntotal
                nnodes_all = 1
                itotal = obj.itotal
                itotal2 = itotal + 2 * nelements * nnodes_all
                ielement = obj.ielement
                ielement2 = ielement + nelements

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 15 * nnodes_all)
                floats1 = floats[:, 1:].reshape(nelements * nnodes_all * 2, 7).copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 15 * nnodes_all).copy()
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

                obj.fiber_curvature[itotal:itotal2] = floats1[:, 0]
                obj.data[obj.itime, itotal:itotal2, :] = real_imag
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if is_vectorized and self.use_vector:  # pragma: no cover
                    self.log.debug(f'vectorize CQUAD4-33 imag SORT{sort_method}')

                n = oes_cquad4_33_complex_15(
                    self, data, obj,
                    nelements, ntotal,
                    is_magnitude_phase)

        elif result_type in [1, 2] and self.num_wide == 9: # random msc
            # _oes_cquad4 is the same as _oes_ctria3
            element_id = self.nonlinear_factor
            if self.is_stress:
                obj_vector_random = RandomPlateStressArray
            else:
                obj_vector_random = RandomPlateStrainArray
            self.data_code['nonlinear_factor'] = element_id

            if self._results.is_not_saved(result_name):
                self._data_factor = 2
                return ndata, None, None
            self._results._found_result(result_name)
            slot = self.get_result(result_name)

            ntotal = 36 * self.factor  # 4*9
            nelements = ndata // ntotal
            nlayers = nelements * 2
            nnodes_expected = 1
            #if self.table_name_str.startswith('OSTRRMS'):
                #print(f'{self.table_name_str} {result_name}: {nelements} {ntotal}')

            auto_return, is_vectorized = self._create_oes_object4(
                nlayers, result_name, slot, obj_vector_random)
            if auto_return:
                self._data_factor = 2
                return nelements * ntotal, None, None

            obj = self.obj
            if self.use_vector and is_vectorized:
                n = nelements * ntotal
                itotal = obj.itotal
                itotal2 = itotal + nelements * 2

                if sort_method == 1:
                    ielement = obj.ielement
                    ielement2 = ielement + nelements

                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = frombuffer(data, dtype=self.idtype)
                        ints1 = ints.reshape(nelements, 9)
                        eids = ints1[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        #print(eids)
                        eids2 = np.vstack([eids, eids]).T.ravel()
                        #print(eids2)

                        # TODO: what about layer 1/2?
                        #print(self.code_information())
                        #print(f'eids.shape={eids.shape} obj.element_node.shape={obj.element_node.shape}')
                        obj.element_node[itotal:itotal2, 0] = eids2

                    floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 9)[:, 1:]
                    #fd1, sx1, sy1, txy1, fd2, fx2, fy2, txy2
                    floats2 = floats.reshape(nelements * nnodes_expected, 8)

                    #[eid_device, fd1, sx1, sy1, txy1,
                    #             fd2, sx2, sy2, txy2,]
                    nf2 = floats2.shape[0]
                    floats3 = floats2.reshape(nf2*2, 4)

                    obj.fiber_curvature[itotal:itotal2] = floats3[:, 0].copy()
                    obj.data[obj.itime, itotal:itotal2, :] = floats3[:, 1:].copy()
                    obj.itotal = itotal2
                    obj.ielement = ielement2

                elif sort_method == 2 and self._analysis_code_fmt == b'f':
                    ielement = obj.itime
                    ie_upper = 2 * ielement
                    ie_lower = 2 * ielement + 1

                    obj.element_node[ie_upper, 0] = dt
                    obj.element_node[ie_lower, 0] = dt
                    #obj._times[obj.itime] = dt

                    floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 9)# [:, 1:]
                    # itime is actually ielement
                    # we grab the element id from the ints for all times
                    if self._analysis_code_fmt == b'i' and obj.itime == 0:
                        #print('analysis_code ', self.analysis_code, self._analysis_code_fmt)
                        ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 9)
                        eids = ints[:, 0] // 10
                        #nids = np.zeros(len(eids), dtype='int32')
                        print(eids)
                        #eids = np.vstack([eids, nids]).T.ravel()
                        #print(eids.shape)
                        #print(eids)
                        #print(obj.element)
                        #assert eids.min() > 0, eids.min()
                        #obj.element[itotal:itotal2, 0] = eids
                        obj._times[itotal:itotal2] = eids
                        aaa
                    elif self._analysis_code_fmt == b'f' and obj.itime == 0:
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
                    obj.fiber_curvature[2*obj.itime:2*obj.itime+2] = floats3[:2, 0].copy()
                    # we apply the data across 2 rows because we have 2 layers
                    obj.data[:, ie_upper, :] = floats3[::2, 1:].copy()
                    obj.data[:, ie_lower, :] = floats3[1::2, 1:].copy()
                else:
                    raise NotImplementedError(self.code_information())
                obj.itotal = itotal2
                #obj.ielement = ielement2
            else:
                if is_vectorized and self.use_vector:  # pragma: no cover
                    self.log.debug(f'vectorize CQUAD4-33 random numwide=9 SORT{sort_method}')
                n = oes_cquad4_33_random_9(self, data, obj, nelements, ntotal)

        elif result_type in [1, 2] and self.num_wide == 11: # random
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

            element_id = self.nonlinear_factor
            if self.is_stress:
                obj_vector_random = RandomPlateVMStressArray
            else:
                obj_vector_random = RandomPlateVMStrainArray
            self.data_code['nonlinear_factor'] = element_id

            if self._results.is_not_saved(result_name):
                return ndata, None, None
            self._results._found_result(result_name)
            slot = self.get_result(result_name)

            ntotal = 44 * factor # 4*11
            nelements = ndata // ntotal
            nlayers = nelements * 2
            nnodes_expected = 1

            auto_return, is_vectorized = self._create_oes_object4(
                nlayers, result_name, slot, obj_vector_random)
            if auto_return:
                self._data_factor = 2
                return nelements * ntotal, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and 0:  # pragma: no cover
                n = nelements * 4 * self.num_wide
                ielement = obj.ielement
                ielement2 = ielement + nelements
                itotal = obj.itotal
                itotal2 = itotal + nelements * nnodes_expected
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype)
                    ints1 = ints.reshape(nelements, 9)
                    eids = ints1[:, 0] // 10
                    print(eids)
                    eids = np.vstack([eids, eids]).T.ravel()
                    print(eids.shape)
                    print(eids)
                    print(obj.element)
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2, 0] = eids

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 11)[:, 1:]
                print(floats.shape)
                #fd, sx, sy, txy,
                floats1 = floats.reshape(nelements * nnodes_expected, 10)
                obj.data[obj.itime, itotal:itotal2, :] = floats1.copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                n = oes_cquad4_33_random_vm_11(self, data, obj, nelements, ntotal)

        elif result_type == 1 and self.num_wide == 17 and self.table_name in [b'OESVM1', b'OESVM2', b'OSTRVM1', b'OSTRVM2']: # freq
            # Table of element stresses for frequency response analysis that includes
            # von Mises stress output in SORT1 format.
            element_id = self.nonlinear_factor
            if self.is_stress:  # TODO: add new complex type
                obj_vector_complex = ComplexPlateVMStressArray
            else:
                obj_vector_complex = ComplexPlateVMStrainArray
            self.data_code['nonlinear_factor'] = element_id

            if self._results.is_not_saved(result_name):
                return ndata, None, None
            self._results._found_result(result_name)
            slot = self.get_result(result_name)

            ntotal = 68 * self.factor  # 4*17
            nelements = ndata // ntotal
            nlayers = nelements * 2
            nnodes_expected = 1
            #self.log.info(f'CQUAD4-33: len(data)={ndata} numwide={self.num_wide} nelements={nelements} nlayers={nlayers}')

            auto_return, is_vectorized = self._create_oes_object4(
                nlayers, result_name, slot, obj_vector_complex)
            if auto_return:
                self._data_factor = 2
                return nelements * ntotal, None, None

            #self.show_data(data)
            #   ELEMENT      FIBER                                     - STRESSES IN ELEMENT  COORDINATE SYSTEM -
            #      ID.       DISTANCE              NORMAL-X                       NORMAL-Y                      SHEAR-XY               VON MISES
            #0     101  -5.000000E-01  -8.152692E-01 /  0.0           -1.321875E+00 /  0.0           -3.158517E+00 /  0.0            5.591334E+00
            #            5.000000E-01   1.728573E+00 /  0.0           -7.103837E+00 /  0.0            2.856040E+00 /  0.0            9.497519E+00
                #floats  = (1011,
                           #-0.5, -0.8152692317962646, 0.0, -1.321874737739563, 0.0, -3.1585168838500977, 0.0, 5.591334342956543,
                           #0.5,   1.7285730838775635, 0.0, -7.103837490081787, 0.0,  2.8560397624969482, 0.0, 9.497518539428711)
            obj = self.obj
            if is_vectorized and self.use_vector and self.sort_method == 1 and 0:  # pragma: no cover
                raise NotImplementedError(self.table_name_str)
            else:
                if is_vectorized and self.use_vector:  # pragma: no cover
                    self.log.debug('vectorize CQUAD4-33 complex '
                                   f'{self.table_name_str} SORT{self.sort_method}')
                n = oes_cquad4_33_complex_vm_17(self, data, obj, nelements, ntotal,
                                                is_magnitude_phase)

        else:  # pragma: no cover
            raise RuntimeError(self.code_information())
        assert self.obj.element_name == self.element_name, self.obj
        assert n > 0
        return n, nelements, ntotal

    def _oes_ctria3(self, data, ndata: int, dt, is_magnitude_phase: bool,
                    result_type: int, prefix: str, postfix: str) -> Tuple[int, Any, Any]:
        """
        reads stress/strain for element type:
         - 74 : CTRIA3-centroidal
         - 83 : CTRIA3-centroidal  (NASA 95)
         - 227: TRIAR-centroidal

        """
        #print('_oes_ctria3')
        n = 0
        if self.is_stress:
            obj_vector_real = RealPlateStressArray
            obj_vector_complex = ComplexPlateStressArray
            if self.element_type in [74, 83]:
                result_name = prefix + 'ctria3_stress' + postfix
            elif self.element_type in [227]:
                result_name = prefix + 'ctriar_stress' + postfix
            else:
                raise NotImplementedError(self.code_information())
        else:
            obj_vector_real = RealPlateStrainArray
            obj_vector_complex = ComplexPlateStrainArray
            if self.element_type in [74, 83]:
                result_name = prefix + 'ctria3_strain' + postfix
            elif self.element_type in [227]:
                result_name = prefix + 'ctriar_strain' + postfix
            else:
                raise NotImplementedError(self.code_information())

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)

        slot = self.get_result(result_name)
        #print(self.element_name, result_name, self.format_code, self.num_wide)
        #table_names = [
        #    b'OES1', b'OES1X', b'OES1X1', b'OSTR1X',
        #    b'OES2', b'OSTR2',
        #    b'OESVM1', b'OESVM2', b'OSTRVM1', b'OSTRVM2',
        #    b'OESPSD1', b'OESRMS1',
        #    b'OESPSD2',  b'OESATO2',  b'OESCRM2',  b'OESNO1',  b'OESXRMS1', b'OESXNO1',
        #    b'OSTRPSD2', b'OSTRATO2', b'OSTRCRM2', b'OSTRNO1', b'OSTRRMS1',
        #]
        #assert self.table_name in table_names, self.table_name
        sort_method = self.sort_method
        element_name_type = f'{self.element_name}-{self.element_type}'

        if self.format_code in [1, 3] and self.num_wide == 17:  # real
            ntotal = 68 * self.factor # 4*17
            nelements = ndata // ntotal
            nlayers = nelements * 2  # 2 layers per node
            #if self.code[0] == 100:
                #print(f'\nnelements={nelements} nlayers={nlayers} {self.code in slot}')
            auto_return, is_vectorized = self._create_oes_object4(
                nlayers, result_name, slot, obj_vector_real)
            if auto_return:
                self._data_factor = 2
                return nelements * ntotal, None, None

            if self.is_debug_file:
                self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                self.binary_debug.write('  #elementi = [eid_device, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
                self.binary_debug.write('  #                        fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,]\n')
                self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            obj = self.obj
            if self.use_vector and is_vectorized and sort_method == 1:
                nfields = 17 * nelements
                nbytes = nfields * 4
                itotal = obj.itotal
                iend = obj.itotal + nlayers

                itime = obj.itime
                if itime == 0:
                    ints = frombuffer(data, dtype=self.idtype8).reshape(nelements, 17)
                    eids = ints[:, 0] // 10
                    #ilayers = ints[:, 1]
                    #ints2 = ints[:, 1:].reshape(nlayers, 8)
                    assert eids.min() > 0, eids
                    obj._times[obj.itime] = dt
                    obj.element_node[itotal:iend:2, 0] = eids
                    obj.element_node[itotal+1:iend+1:2, 0] = eids
                    #obj.element_node[itotal:iend, 1] = 0

                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, 17)
                floats1 = floats[:, 1:].reshape(nlayers, 8).copy()
                obj.data[obj.itime, itotal:iend, :] = floats1
                obj._times[obj.itime] = dt
                obj.itotal += nlayers
                n = nbytes
            else:
                if is_vectorized and self.use_vector:  # pragma: no cover
                    self.log.debug(f'vectorize centroidal tri: {element_name_type} real SORT{sort_method}')
                n = oes_ctria3_real_17(self, data, obj,
                                       ntotal, nelements, dt)
            if self.is_sort1:
                assert obj.element_node[:, 0].min() > 0, obj.element_node[:, 0]

        elif self.format_code in [2, 3] and self.num_wide == 15:  # imag
            ntotal = 60 * self.factor # 4*15
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_complex)
            if auto_return:
                self._data_factor = 2
                return nelements * self.num_wide * 4, None, None
            obj = self.obj

            if self.use_vector and is_vectorized and sort_method == 1:
                n = nelements * ntotal
                itotal = obj.itotal
                itotal2 = itotal + nelements * 2
                ielement = obj.ielement
                ielement2 = ielement + nelements

                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, 15)
                floats1 = floats[:, 1:].reshape(nelements * 2, 7).copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype8).reshape(nelements, 15).copy()
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

                obj.fiber_curvature[itotal:itotal2] = floats1[:, 0]
                obj.data[obj.itime, itotal:itotal2, :] = real_imag
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if is_vectorized and self.use_vector:  # pragma: no cover
                    self.log.debug(f'vectorize CTRIA3 imag SORT{sort_method}')
                struct1 = Struct(self._endian + mapfmt(self._analysis_code_fmt + b'14f', self.size))
                cen = 0 # CEN/3
                for unused_i in range(nelements):
                    edata = data[n:n + ntotal]
                    out = struct1.unpack(edata)
                    (eid_device,
                     fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
                     fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i,) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)

                    if self.is_debug_file:
                        self.binary_debug.write('  OESC %s - eid=%i; C=[%s]\n' % (
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
        #elif self.format_code == 1 and self.num_wide == 9: # random?
            #msg = self.code_information()
            #return self._not_implemented_or_skip(data, ndata, msg), None, None
        elif self.format_code in [2, 3] and self.num_wide == 17 and self.table_name in [b'OESVM1', b'OESVM2', b'OSTRVM1', b'OSTRVM2']:
            # freq:
            # # random; CTRIA3
            assert self.table_name in [b'OESVM1', b'OESVM2', b'OSTRVM1', b'OSTRVM2'], self.code_information()

            element_id = self.nonlinear_factor
            if self.is_stress:  # TODO: add new complex type
                obj_vector_complex = ComplexPlateVMStressArray
            else:
                obj_vector_complex = ComplexPlateVMStrainArray
            self.data_code['nonlinear_factor'] = element_id

            if self._results.is_not_saved(result_name):
                return ndata, None, None
            self._results._found_result(result_name)
            slot = self.get_result(result_name)

            ntotal = 68 * self.factor  # 17*4
            nelements = ndata // ntotal
            nlayers = nelements * 2
            nnodes_expected = 1

            auto_return, is_vectorized = self._create_oes_object4(
                nlayers, result_name, slot, obj_vector_complex)
            if auto_return:
                self._data_factor = 2
                return nelements * ntotal, None, None

            obj = self.obj
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
            n = oes_ctria3_complex_vm_17(self, data, obj, nelements, ntotal, dt,
                                         is_magnitude_phase)
            assert n is not None, n

        elif self.format_code in [1, 2, 3] and self.num_wide == 11: # random; CTRIA3
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

            element_id = self.nonlinear_factor
            if self.is_stress:
                obj_vector_random = RandomPlateVMStressArray
            else:
                obj_vector_random = RandomPlateVMStrainArray
            self.data_code['nonlinear_factor'] = element_id

            if self._results.is_not_saved(result_name):
                return ndata, None, None
            self._results._found_result(result_name)
            slot = self.get_result(result_name)

            ntotal = 44 * self.factor  # 4*11
            nelements = ndata // ntotal
            nlayers = nelements * 2
            nnodes_expected = 1

            auto_return, is_vectorized = self._create_oes_object4(
                nlayers, result_name, slot, obj_vector_random)
            if auto_return:
                self._data_factor = 2
                return nelements * ntotal, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and 0:  # pragma: no cover
                n = nelements * 4 * self.num_wide
                ielement = obj.ielement
                ielement2 = ielement + nelements
                itotal = obj.itotal
                itotal2 = itotal + nelements * nnodes_expected
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype)
                    ints1 = ints.reshape(nelements, 9)
                    eids = ints1[:, 0] // 10
                    print(eids)
                    eids = np.vstack([eids, eids]).T.ravel()
                    print(eids.shape)
                    print(eids)
                    print(obj.element)
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2, 0] = eids

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 11)[:, 1:]
                print(floats.shape)
                #fd, sx, sy, txy,
                floats1 = floats.reshape(nelements * nnodes_expected, 10)
                obj.data[obj.itime, itotal:itotal2, :] = floats1.copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if is_vectorized and self.use_vector and obj.itime == 0:  # pragma: no cover
                    self.log.debug(f'vectorize {element_name_type} random numwide=11 SORT{sort_method}')
                n = oes_ctria3_random_vm_11(self, data, obj, nelements, ntotal)

        elif self.format_code in [1, 2] and self.num_wide == 9: # random MSC stress/strain; CTRIA3
            # _oes_cquad4 is the same as _oes_ctria3
            element_id = self.nonlinear_factor
            if self.is_stress:
                obj_vector_random = RandomPlateStressArray
                #result_name = prefix + 'ctria3_stress' + postfix
            else:
                obj_vector_random = RandomPlateStrainArray
                #result_name = prefix + 'ctria3_strain' + postfix

            self.data_code['nonlinear_factor'] = element_id

            if self._results.is_not_saved(result_name):
                return ndata, None, None
            self._results._found_result(result_name)
            slot = self.get_result(result_name)

            ntotal = 36 * self.factor # 4*9
            nelements = ndata // ntotal
            nlayers = nelements * 2
            nnodes_expected = 1

            auto_return, is_vectorized = self._create_oes_object4(
                nlayers, result_name, slot, obj_vector_random)
            if auto_return:
                self._data_factor = 2
                return nelements * ntotal, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and 0:  # pragma: no cover
                n = nelements * 4 * self.num_wide
                ielement = obj.ielement
                ielement2 = ielement + nelements
                itotal = obj.itotal
                itotal2 = itotal + nelements * nnodes_expected
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype)
                    ints1 = ints.reshape(nelements, 9)
                    eids = ints1[:, 0] // 10
                    print(eids)
                    eids = np.vstack([eids, eids]).T.ravel()
                    print(eids.shape)
                    print(eids)
                    print(obj.element)
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2, 0] = eids

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 9)[:, 1:]
                print(floats.shape)
                #fd, sx, sy, txy,
                floats1 = floats.reshape(nelements * nnodes_expected, 8)
                obj.data[obj.itime, itotal:itotal2, :] = floats1.copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if is_vectorized and self.use_vector:  # pragma: no cover
                    self.log.debug(f'vectorize {element_name_type} random2 SORT{sort_method}')
                n = oes_ctria3_random_9(self, data, obj, nelements, ntotal)

        else:  # pragma: no cover
            raise RuntimeError(self.code_information())
        assert n is not None, self.code_information()
        return n, nelements, ntotal

    def _oes_cquad4_144(self, data: bytes, ndata: int, dt, is_magnitude_phase: bool,
                        result_type: int, prefix: str, postfix: str) -> Tuple[int, Any, Any]:
        """
        reads stress/strain for element type:
         - 64 : CQUAD8
         - 70 : CTRIAR
         - 75 : CTRIA6
         - 82 : CQUADR
         - 144 : CQUAD4-bilinear

        """
        n = 0
        size = self.size
        if self.is_stress:
            obj_vector_real = RealPlateStressArray
            obj_vector_complex = ComplexPlateStressArray
            obj_vector_random = RandomPlateStressArray
            if self.element_type == 64: # CQUAD8
                result_name = prefix + 'cquad8_stress' + postfix
                #gridC = 'CEN/8'
            elif self.element_type == 70:  # CTRIAR
                result_name = prefix + 'ctriar_stress' + postfix
                #gridC = 'CEN/3'
            elif self.element_type == 75:  # CTRIA6
                result_name = prefix + 'ctria6_stress' + postfix
                #gridC = 'CEN/6'
            elif self.element_type == 82:  # CQUADR
                result_name = prefix + 'cquadr_stress' + postfix
                #gridC = 'CEN/4'
            elif self.element_type == 144:  # CQUAD4-bilinear
                # there's no nead to separate this with centroidal strain
                # because you can only have one in a given OP2
                result_name = prefix + 'cquad4_stress' + postfix
                #gridC = 'CEN/4'
            else:  # pragma: no cover
                raise RuntimeError(self.code_information())
                #msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
                #return self._not_implemented_or_skip(data, ndata, msg)
        else:
            obj_vector_real = RealPlateStrainArray
            obj_vector_complex = ComplexPlateStrainArray
            obj_vector_random = RandomPlateStrainArray

            if self.element_type == 64: # CQUAD8
                result_name = prefix + 'cquad8_strain' + postfix
                #gridC = 'CEN/8'
            elif self.element_type == 70:  # CTRIAR
                result_name = prefix + 'ctriar_strain' + postfix
                #gridC = 'CEN/3'
            elif self.element_type == 75:  # CTRIA6
                result_name = prefix + 'ctria6_strain' + postfix
                #gridC = 'CEN/6'
            elif self.element_type == 82: # CQUADR
                result_name = prefix + 'cquadr_strain' + postfix
                #gridC = 'CEN/4'
            elif self.element_type == 144: # CQUAD4-bilinear
                # there's no nead to separate this with centroidal strain
                # because you can only have one in a given OP2
                result_name = prefix + 'cquad4_strain' + postfix
                #gridC = 'CEN/4'
            else:  # pragma: no cover
                raise RuntimeError(self.code_information())

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)

        if self.element_type in [64, 82, 144]:
            nnodes = 4 # + 1 centroid
        elif self.element_type in [70, 75]:
            nnodes = 3 # + 1 centroid
        else:  # pragma: no cover
            raise RuntimeError(self.code_information())
        nnodes_all = nnodes + 1 # adding the centroid

        slot = self.get_result(result_name)
        numwide_real = 2 + 17 * nnodes_all
        numwide_imag = 2 + 15 * nnodes_all
        numwide_random = 2 + 9 * nnodes_all

        #numwide_imag2 = 2 + 16 * nnodes_all
        #print('%s real=%s imag=%s imag2=%s random=%s' % (
            #self.element_name, numwide_real, numwide_imag, numwide_imag2, numwide_random
        #))
        #etype = self.element_name
        #grid_center = 'CEN/%i' % nnodes

        # OESVM1/2 (complex)
        #   87 - CQUAD8

        # OSTRNO1
        #   47 - CQUAD4

        # CQUAD8 real=87 imag=77 imag2=82 random=47
        # CQUAD4 ???=
        sort_method = self.sort_method
        element_name_type = f'{self.element_name}-{self.element_type}'
        #print(self.code_information())
        if result_type == 0 and self.num_wide == numwide_real:  # real
            ntotal = 4 * (2 + 17 * nnodes_all) * self.factor
            nelements = ndata // ntotal
            assert ndata % ntotal == 0
            nlayers = 2 * nelements * nnodes_all  # 2 layers per node

            auto_return, is_vectorized = self._create_oes_object4(
                nlayers, result_name, slot, obj_vector_real)
            if auto_return:
                self._data_factor = 2 * nnodes_all  # number of "layers" for an element
                return nelements * ntotal, None, None

            obj = self.obj
            #print('dt=%s, itime=%s' % (obj.itime, dt))
            if self.use_vector and is_vectorized and sort_method == 1:
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
                    ints = frombuffer(data, dtype=self.idtype8).reshape(nelements, numwide_real)
                    ints1 = ints[:, 2:].reshape(nlayers//2, 17)[:, 0].reshape(nelements, nnodes_all).copy()
                    ints1[:, 0] = 0.
                    nids = ints1.ravel()

                    eids = ints[:, 0] // 10
                    eids2 = array([eids] * (nnodes_all * 2), dtype=self.idtype8).T.ravel()
                    nids2 = vstack([nids, nids]).T.ravel()
                    obj.element_node[istart:iend, 0] = eids2
                    obj.element_node[istart:iend, 1] = nids2
                    #assert obj.element_node[:iend, 0].min() > 0, eids2
                    if obj.nonlinear_factor is not None:
                        float_mask = np.arange(nelements * numwide_real, dtype=np.int32).reshape(nelements, numwide_real)
                        float_mask1 = float_mask[:, 2:].reshape(nlayers // 2, 17)[:, 1:].reshape(nlayers, 8)
                        obj.float_mask = float_mask1

                if obj.nonlinear_factor is not None:
                    results = frombuffer(data, dtype=self.fdtype8)[obj.float_mask].copy()
                else:
                    floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, numwide_real)
                    floats1 = floats[:, 2:].reshape(nlayers // 2, 17)
                    results = floats1[:, 1:].reshape(nlayers, 8).copy()

                #[fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm]
                obj.data[obj.itime, istart:iend, :] = results
                assert obj.element_node[:, 0].min() > 0, obj.element_node[:, 0]
            else:
                if is_vectorized and self.use_vector:  # pragma: no cover
                    self.log.debug(f'vectorize nodal shell: {element_name_type}... real SORT{sort_method}')

                n = oes_cquad4_144_real(self, data, ndata, obj,
                                        nelements, nnodes, dt)
            if self.is_sort1:
                assert obj.element_node[:, 0].min() > 0, obj.element_node[:, 0]

        elif result_type == 1 and self.num_wide == numwide_imag:  # complex
            ntotal = numwide_imag * 4 * self.factor
            #assert self.num_wide * 4 == ntotal, 'numwide*4=%s ntotal=%s' % (self.num_wide*4, ntotal)
            nelements = ndata // ntotal
            nlayers = nelements * 2 * nnodes_all
            #print(element_name_type)
            #print('ndata', ndata)
            #print('ntotal', ntotal)
            #print('nelements', nelements)
            #print('nlayers', nlayers)

            auto_return, is_vectorized = self._create_oes_object4(
                nlayers, result_name, slot, obj_vector_complex)
            if auto_return:
                self._data_factor = 2 * nnodes_all
                #if self.num_wide == 77:
                    #print('ntotal =', ndata, ntotal, self.num_wide)
                    #print('nelements * ntotal =', nelements * ntotal)
                return nelements * ntotal, None, None
            obj = self.obj

            if self.use_vector and is_vectorized and sort_method == 1:
                n = nelements * ntotal
                itotal = obj.itotal
                itotal2 = itotal + nelements * (nnodes_all * 2)
                ielement = obj.ielement
                ielement2 = ielement + nelements

                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, numwide_imag)
                floats1 = floats[:, 2:].reshape(nelements * nnodes_all, 15)
                floats2 = floats1[:, 1:].reshape(nelements * nnodes_all * 2, 7).copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype8).reshape(nelements, numwide_imag).copy()
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

                obj.fiber_curvature[itotal:itotal2] = floats2[:, 0]
                obj.data[obj.itime, itotal:itotal2, :] = real_imag
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if is_vectorized and self.use_vector:  # pragma: no cover
                    self.log.debug(f'vectorize CQUAD4-144/{element_name_type}... imag SORT{sort_method}')
                # nnodes_cquad4 = 5
                # nelements = 3
                # nlayers = nelements * nodes_cquad4 * 2 = 3*5*2 = 30
                # ntotal = nlayers
                n = oes_cquad4_144_complex_77(self, data, obj,
                                              nelements, nnodes,
                                              dt, is_magnitude_phase)
        #elif self.format_code == 1 and self.num_wide == numwide_random: # random
            #msg = self.code_information()
            #msg += '  numwide=%s numwide_real=%s numwide_imag=%s numwide_random=%s' % (
                #self.num_wide, numwide_real, numwide_imag, numwide_random)
            #return self._not_implemented_or_skip(data, ndata, msg), None, None
        elif result_type == 2 and self.num_wide == numwide_random: # random
            # 47 - CQUAD8-64
            # 38 - CTRIAR-70
            ntotal = self.num_wide * 4 * self.factor
            nelements = ndata // ntotal
            assert ndata % ntotal == 0
            nlayers = 2 * nelements * nnodes_all  # 2 layers per node
            #self.log.info(f'random quad-144 ntotal={ntotal} ndata={ndata} ntotal={ntotal} numwide={self.num_wide} -> nelements={nelements}')

            #if self.read_mode == 1:
                #msg = ''
                #return self._not_implemented_or_skip(data, ndata, msg), None, None

            auto_return, is_vectorized = self._create_oes_object4(
                nlayers, result_name, slot, obj_vector_random)
            if auto_return:
                self._data_factor = 2 * nnodes_all  # number of "layers" for an element
                return nelements * ntotal, None, None

            obj = self.obj
            #print('dt=%s, itime=%s' % (obj.itime, dt))
            if self.use_vector and is_vectorized and 0:
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
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, numwide_real)
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
                    results = frombuffer(data, dtype=self.fdtype)[obj.float_mask].copy()
                else:
                    floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, numwide_real)
                    floats1 = floats[:, 2:].reshape(nlayers // 2, 17)
                    results = floats1[:, 1:].reshape(nlayers, 8).copy()

                #[fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm]
                obj.data[obj.itime, istart:iend, :] = results
            else:
                if is_vectorized and self.use_vector:  # pragma: no cover
                    self.log.debug(f'vectorize CQUAD4-144/{element_name_type}... random SORT{sort_method}')
                #numwide_random = 2 + 9 * nnodes_all
                n = oes_cquad4_144_random(self, data, obj, nelements, nnodes, ndata)

            #if self.read_mode == 1:
                #msg = ''
                #return self._not_implemented_or_skip(data, ndata, msg), None, None
            ##self.show_data(data[:44])
            #ntotal = 44
            #struct1 = Struct(self._endian + b'i4si 8f')
            #for i in ra
            #for i in range(20):
                #edata = data[n:n+ntotal]
                #out = struct1.unpack(edata)
                #self.show_data(edata)
                #print(out)
                #n += ntotal

            ## 47 - CQUAD8-64
            ##msg = self.code_information()
            #msg = '%s-CQUAD4-numwide=%s format_code=%s;\n numwide_real=%s numwide_imag=%s numwide_random=%s' % (
                #self.table_name_str, self.num_wide, self.format_code,
                #numwide_real, numwide_imag, numwide_random)
            #return self._not_implemented_or_skip(data, ndata, msg), None, None
        elif self.table_name in [b'OESVM1', b'OESVM2', b'OSTRVM1', b'OSTRVM2']:
            # 82  CQUADR -> 87
            # CQUAD8 sort_method=2 ntotal=348 nelements=3
            # CTRIA6 sort_method=2 ntotal=348 nelements=2

            msg = self.code_information()
            if result_type == 1: # complex
                # ndata = 3828
                # ???
                #
                # ndata=1044
                assert self.num_wide in [70, 87], self.code_information()
                ntotal = self.num_wide * self.size # 87*4
                nelements = ndata // ntotal
                nlayers = nelements * 2
                assert ndata % ntotal == 0

                if self.is_stress:
                    obj_vector_complex = ComplexPlateVMStressArray
                else:
                    obj_vector_complex = ComplexPlateVMStrainArray

                auto_return, is_vectorized = self._create_oes_object4(
                    nlayers, result_name, slot, obj_vector_complex)

                if auto_return:
                    self._data_factor = 2 * nnodes_all  # number of "layers" for an element
                    #self._data_factor = 2
                    return nelements * ntotal, None, None
                #if self.read_mode == 1:
                    #return self._not_implemented_or_skip(data, ndata, msg), None, None

                #self.show_data(data)
                #print(ndata, ntotal)
                obj = self.obj
                n = oes_cquad4_complex_vm_87(self, data, obj, nelements, nnodes_all,
                                             is_magnitude_phase)

            #if result_type == 1 and self.num_wide in [70, 87]:
                # 70 - CTRIA6-75
                # 87 - CQUAD4-144
                #pass
            else:
                #msg = (f'skipping {self.table_name_str}-{self.element_name}: numwide={self.num_wide} '
                       #f'result_type={self.result_type} (complex);\n numwide_real={numwide_real} '
                       #f'numwide_imag={numwide_imag} numwide_random={numwide_random}')
                #print(msg)
                raise NotImplementedError(msg)
            #return self._not_implemented_or_skip(data, ndata, msg), None, None

        #elif self.format_code in [2, 3] and self.num_wide == 70:
            ## 87 - CQUAD4-144
            ##msg = self.code_information()
            #msg = '%s-CTRIA6-numwide=%s numwide_real=%s numwide_imag=%s numwide_random=%s' % (
                #self.table_name_str, self.num_wide, numwide_real, numwide_imag, numwide_random)
            #return self._not_implemented_or_skip(data, ndata, msg), None, None

        elif result_type == 2 and self.table_name in [b'OESXRMS1', b'OESXNO1']:
            # CTRIA6-75  numwide=46  C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1081x.op2
            # CQUAD8-64  numwide=64  C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1081x.op2

            # corners + centroid
            if self.element_type == 75:  # CTRIA6  (numwide=46)
                # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1081x.op2
                nnodes = 4
            elif self.element_type in [64, 144]:
                # CQUAD8     (numwide=57) C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1081x.op2
                # CQUAD4-144 (numwide=57) C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\plate_111o.op2
                nnodes = 5
            else:
                raise RuntimeError(self.code_information())

            assert 2 + 11 * nnodes == self.num_wide, self.code_information()
            if self.is_stress:
                obj_vector_random = RandomPlateVMStressArray
            else:
                obj_vector_random = RandomPlateVMStrainArray

            ntotal = self.num_wide * size
            nelements = ndata // ntotal
            nlayers = nelements * nnodes_all * 2
            assert ndata % ntotal == 0
            #print(f'selement_name={self.element_name} sort_method={self.sort_method} ntotal={ntotal} num_wide={self.num_wide} nelements={nelements} ndata={ndata} nlayers={nlayers}')
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
            #if self.read_mode == 2:
                #self.show_data(data)
                #ddd
            auto_return, is_vectorized = self._create_oes_object4(
                nlayers, result_name, slot, obj_vector_random)
            if auto_return:
                self._data_factor = 2 * nnodes_all  # number of "layers" for an element
                return nelements * ntotal, None, None

            obj = self.obj
            n = oes_cquad4_random_vm_57(self, data, self.obj, nelements, ntotal, nnodes,
                                        dt)
        else:  # pragma: no cover
            raise RuntimeError(self.code_information())
        return n, nelements, ntotal

    def _oes_shells_nonlinear(self, data, ndata, dt, is_magnitude_phase,
                              result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 88 : CTRIA3NL
         - 90 : CQUAD4NL

        """
        n = 0
        if self.is_stress:
            if self.element_type == 88:
                result_name = prefix + 'ctria3_stress' # + postfix  nonlinear_
            elif self.element_type == 90:
                result_name = prefix + 'cquad4_stress'  + postfix # nonlinear_
            else:
                raise RuntimeError(self.element_type)
        else:
            if self.element_type == 88:
                result_name = prefix + 'ctria3_strain' + postfix # nonlinear_
            elif self.element_type == 90:
                result_name = prefix + 'cquad4_strain' + postfix # nonlinear_
            else:
                raise RuntimeError(self.element_type)

        slot = self.get_result(result_name)
        self._results._found_result(result_name)
        #print(self.code_information())

        if self.format_code == 1 and self.num_wide == 13 and self.element_type in [88, 90]:  # real
            # single layered hyperelastic (???) ctria3, cquad4
            ntotal = 52 * self.factor  # 4*13
            nelements = ndata // ntotal

            obj_vector_real = RealNonlinearPlateArray
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return nelements * ntotal, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * self.num_wide * 4

                ielement = obj.ielement
                ielement2 = ielement + nelements
                obj._times[obj.itime] = dt

                self.obj_set_element(obj, ielement, ielement2, data, nelements)

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 13).copy()

                #[fiber_distance, oxx, oyy, ozz, txy, exx, eyy, ezz, exy, es, eps, ecs]
                floats[:, 1] = 0
                obj.data[obj.itime, ielement:ielement2, :] = floats[:, 1:]
                obj.ielement = ielement2
                obj.itotal = ielement2
            else:
                if is_vectorized and self.use_vector:  # pragma: no cover
                    self.log.debug('vectorize CTRIA3/CQUAD4_NL real SORT%s' % self.sort_method)
                struct1 = Struct(self._endian + self._analysis_code_fmt + b'12f')  # 1+12=13
                for unused_i in range(nelements):
                    edata = data[n:n + ntotal]
                    out = struct1.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('CQUADNL-90 - %s\n' % str(out))

                    (eid_device, fd1,
                     sx1, sy1, sz1, txy1, es1, eps1, ecs1,
                     ex1, ey1, ez1, exy1) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    obj.add_new_eid_sort1(
                        dt, eid, self.element_type, fd1,
                        sx1, sy1, sz1, txy1, es1, eps1, ecs1,
                        ex1, ey1, ez1, exy1)
                    n += ntotal
        elif self.format_code == 1 and self.num_wide == 25 and self.element_type in [88, 90]:
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

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                self._data_factor = 2
                return nelements * ntotal, None, None

            #return nelements * self.num_wide * 4
            obj = self.obj
            is_vectorized = False
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * self.num_wide * 4

                ielement = obj.ielement
                ielement2 = ielement + nelements
                itotal = obj.itotal
                itotal2 = itotal + nelements * 2
                obj._times[obj.itime] = dt
                #print('ielement=%s:%s' % (ielement, ielement2))
                #print('itotal=%s:%s' % (itotal, itotal2))

                if obj.itime == 0:
                    try:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 25)
                    except ValueError:
                        unused_values = fromstring(data, dtype=self.idtype)

                    eids = ints[:, 0] // 10
                    #eids2 = vstack([eids, eids]).T.ravel()
                    #print(eids.tolist())
                    obj.element[ielement:ielement2] = eids  # 150
                     #print(obj.element_node[:10, :])

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 25)[:, 1:]

                #[fiber_distance, oxx, oyy, ozz, txy, exx, eyy, ezz, exy, es, eps, ecs]
                #floats[:, 1] = 0
                obj.data[obj.itime, itotal:itotal2, :] = floats.reshape(nelements * 2, 12).copy()
                #obj.data[obj.itime, ielement:ielement2, :] = floats[:, 1:]
                obj.ielement = ielement2
                obj.itotal = itotal2
            else:
                if is_vectorized and self.use_vector:  # pragma: no cover
                    self.log.debug('vectorize CTRIA3/CQUAD4_NL imag SORT%s' % self.sort_method)
                etype = self.element_type
                struct1 = Struct(self._endian + mapfmt(self._analysis_code_fmt + b'24f', self.size)) # 1+24=25
                for unused_i in range(nelements):
                    edata = data[n:n + ntotal]
                    out = struct1.unpack(edata)
                    if self.is_debug_file:
                        eid = out[0] // 10
                        self.binary_debug.write('CQUADNL-90 - %s : %s\n' % (eid, str(out)))

                    (eid_device,
                     fd1, sx1, sy1, undef1, txy1, es1, eps1, ecs1, ex1, ey1, undef2, etxy1,
                     fd2, sx2, sy2, undef3, txy2, es2, eps2, ecs2, ex2, ey2, undef4, etxy2) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    obj.add_new_eid_sort1(
                        dt, eid, etype,
                        fd1, sx1, sy1, undef1, txy1, es1, eps1, ecs1, ex1, ey1, undef2, etxy1)
                    obj.add_sort1(
                        dt, eid, etype,
                        fd2, sx2, sy2, undef3, txy2, es2, eps2, ecs2, ex2, ey2, undef4, etxy2)
                    n += ntotal
        elif self.format_code == 1 and self.num_wide == 0: # random
            msg = self.code_information()
            return self._not_implemented_or_skip(data, ndata, msg)
        else:  # pragma: no cover
            raise RuntimeError(self.code_information())
        return n, nelements, ntotal

    def _oes_shells_composite(self, data, ndata, dt, is_magnitude_phase,
                              result_type: int, prefix: str, postfix: str) -> Tuple[int, Any, Any]:
        """
        reads stress/strain for element type:
         - 95 : CQUAD4
         - 96 : CQUAD8
         - 97 : CTRIA3
         - 98 : CTRIA6 (composite)
         - 232 : QUADRLC (CQUADR-composite)
         - 233 : TRIARLC (CTRIAR-composite)

        """
        table_name = self.table_name
        assert isinstance(table_name, bytes), table_name
        n = 0
        if self.is_stress:
            obj_vector_real = RealCompositePlateStressArray
            #obj_vector_complex = ComplexCompositePlateStressArray
            obj_vector_random = RandomCompositePlateStressArray
            if self.element_type == 95: # CQUAD4
                result_name = prefix + 'cquad4_composite_stress' + postfix
            elif self.element_type == 96:  # CQUAD8
                result_name = prefix + 'cquad8_composite_stress' + postfix
            elif self.element_type == 97:  # CTRIA3
                result_name = prefix + 'ctria3_composite_stress' + postfix
            elif self.element_type == 98:  # CTRIA6
                result_name = prefix + 'ctria6_composite_stress' + postfix
            elif self.element_type == 232:  # CQUADR
                result_name = prefix + 'cquadr_composite_stress' + postfix
            elif self.element_type == 233:  # CTRIAR
                result_name = prefix + 'ctriar_composite_stress' + postfix
            else:  # pragma: no cover
                raise RuntimeError(self.code_information())
        else:
            obj_vector_real = RealCompositePlateStrainArray
            #obj_vector_complex = ComplexCompositePlateStrainArray
            obj_vector_random = RandomCompositePlateStrainArray
            if self.element_type == 95: # CQUAD4
                result_name = prefix + 'cquad4_composite_strain' + postfix
            elif self.element_type == 96:  # CQUAD8
                result_name = prefix + 'cquad8_composite_strain' + postfix
            elif self.element_type == 97:  # CTRIA3
                result_name = prefix + 'ctria3_composite_strain' + postfix
            elif self.element_type == 98:  # CTRIA6
                result_name = prefix + 'ctria6_composite_strain' + postfix
            elif self.element_type == 232:  # CQUADR
                result_name = prefix + 'cquadr_composite_strain' + postfix
            elif self.element_type == 233:  # CTRIAR
                result_name = prefix + 'ctriar_composite_strain' + postfix
            else:  # pragma: no cover
                raise RuntimeError(self.code_information())

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)

        etype = self.element_name
        sort_method = self.sort_method
        if result_type == 0 and self.num_wide == 11:  # real
            #                    S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A R )
            #   ELEMENT      PLY STRESSES IN FIBER AND MATRIX DIRECTIONS   INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      MAX
            #     ID          ID   NORMAL-1     NORMAL-2     SHEAR-12    SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        SHEAR
            #      7070        1  7.50000E-01  3.00000E+00  9.86167E-08  -6.58903E-08  3.00000E+00   90.00  3.00000E+00  7.50000E-01  1.12500E+00
            #      7070        2 -7.50000E-01 -3.00000E+00 -9.86167E-08   0.0          0.0           -0.00 -7.50000E-01 -3.00000E+00  1.12500E+00
            ntotal = 44 * self.factor
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return nelements * ntotal, None, None

            obj = self.obj
            if self.is_debug_file:
                self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                self.binary_debug.write('  element1 = [eid_device, layer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)]\n')
                self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if self.use_vector and is_vectorized and sort_method == 1:
                n = nelements * self.num_wide * 4

                istart = obj.itotal
                iend = istart + nelements
                obj._times[obj.itime] = dt

                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype8).reshape(nelements, 11).copy()
                    eids = ints[:, 0] // 10
                    nids = ints[:, 1]
                    obj.element_layer[istart:iend, 0] = eids
                    obj.element_layer[istart:iend, 1] = nids

                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, 11)
                #[o1, o2, t12, t1z, t2z, angle, major, minor, ovm]
                obj.data[obj.itime, istart:iend, :] = floats[:, 2:].copy()
            else:
                if is_vectorized and self.use_vector:  # pragma: no cover
                    self.log.debug(f'vectorize COMP_SHELL real SORT{sort_method}')

                n = oes_comp_shell_real_11(self, data, ndata, obj,
                                           ntotal, nelements, etype, dt)

        #elif result_type == 1 and self.num_wide == 9:  # TODO: imag? - not done...
            # TODO: vectorize
            #raise NotImplementedError('imaginary composite stress?')
            #msg = self.code_information()
            #nelements = ndata // ntotal
            #obj_vector_complex = None
            #auto_return, is_vectorized = self._create_oes_object4(
                #nelements, result_name, slot, obj_vector_complex)
            #if auto_return:
                #return nelements * self.num_wide * 4, None, None

            ## TODO: this is an OEF result???
            ##    furthermore the actual table is calle dout as
            ##    'i8si4f4s', not 'i8si3fi4s'
            #ntotal = 36
            #nelements = ndata // ntotal
            #s = self.struct_i
            #s2 = Struct(self._endian + b'8si3fi4s')
            #s3 = Struct(self._endian + b'8si4f4s')
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

                #if self.is_debug_file:
                    #self.binary_debug.write('%s-%s - (%s) + %s\n' % (self.element_name, self.element_type, eid_device, str(out)))
                #obj.add_new_eid_sort1(dt, eid, theory, lamid, fp, fm, fb, fmax, fflag)
                #n += ntotal
            #raise NotImplementedError('this is a really weird case...')
        elif result_type == 1 and self.num_wide == 11 and table_name in [b'OESCP', b'OESTRCP']:  # complex
            # OESCP - STRAINS IN LAYERED COMPOSITE ELEMENTS (QUAD4)
            ntotal = 44 * self.factor
            nelements = ndata // ntotal

            self.log.warning(f'skipping complex {self.table_name_str}-PCOMP')
            return nelements * ntotal, None, None

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, ComplexLayeredCompositesArray)
            if auto_return:
                return nelements * ntotal, None, None

            obj = self.obj
            n = oes_shell_composite_complex_11(self, data, obj,
                                               ntotal, nelements, sort_method,
                                               dt, is_magnitude_phase)
            return nelements * ntotal, None, None

        elif result_type == 0 and self.num_wide == 9 and table_name == b'OESRT': # real
            # strength_ratio.cquad4_composite_stress
            ntotal = 36 * self.factor
            nelements = ndata // ntotal
            if self.read_mode == 1:
                return nelements * ntotal, None, None

            # not 100%
            self.log.warning(f'skipping oes_shell_composite; {self.element_name}-{self.element_type} '
                             f'(numwide={self.num_wide}) {self.table_name_str}')
            struct1 = Struct(self._endian + self._analysis_code_fmt + b' 8s i 3f if')
            for unused_i in range(nelements):
                edata = data[n:n+ntotal]
                #self.show_data(edata)
                out = struct1.unpack(edata)
                #print(out)

                #(eid_device, failure_theory, ply_id, strength_ratio_ply, failure_index_bonding, strength_ratio_bonding, flag, flag2) = out
                #(eid_device, failure_theory, ply_id, ratio_ply,          ratio_bonding,         ratio_element,          seven, eight, nine) = out
                #eid, dt = get_eid_dt_from_eid_device(
                    #eid_device, self.nonlinear_factor, sort_method)
                #print(eid, out)

                #if self.is_debug_file:
                    #self.binary_debug.write('%s-%s - (%s) + %s\n' % (self.element_name, self.element_type, eid_device, str(out)))
                #obj.add_new_eid_sort1(dt, eid, theory, lamid, fp, fm, fb, fmax, fflag)
                n += ntotal
            return nelements * ntotal, None, None
        elif result_type == 1 and self.num_wide == 13 and table_name in [b'OESVM1C', b'OSTRVM1C']: # complex
            self.log.warning(f'skipping complex {self.table_name_str}-PCOMP (numwide=13)')
            ntotal = 52 * self.factor
            nelements = ndata // ntotal
            return nelements * ntotal, None, None
            self.table_name = table_name
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, ComplexLayeredCompositesArray)
            if auto_return:
                return nelements * ntotal, None, None

            if is_vectorized and self.use_vector:  # pragma: no cover
                self.log.debug(f'vectorize COMP_SHELL random SORT{sort_method} (numwide=13)')

            obj = self.obj
            n = oes_shell_composite_complex_13(self, data, obj,
                                               ntotal, nelements, sort_method,
                                               dt, is_magnitude_phase)
            #return nelements * ntotal, None, None

        elif result_type == 2 and self.num_wide == 7:
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
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_random)
            if auto_return:
                return nelements * self.num_wide * 4, None, None

            obj = self.obj
            struct1 = Struct(self._endian + self._analysis_code_fmt + b'i5f')
            for unused_i in range(nelements):
                edata = data[n:n+ntotal]
                out = struct1.unpack(edata)

                (eid_device, ply_id, oxx, oyy, txy, txz, tyz) = out
                eid, dt = get_eid_dt_from_eid_device(
                    eid_device, self.nonlinear_factor, sort_method)
                #print(eid, out)

                #if self.is_debug_file:
                    #self.binary_debug.write('%s-%s - (%s) + %s\n' % (self.element_name, self.element_type, eid_device, str(out)))
                #print(obj)
                obj.add_sort1_7words(dt, eid, ply_id, oxx, oyy, txy, txz, tyz)
                n += ntotal

        elif result_type == 2 and self.num_wide == 8:
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
            msg = self.code_information()
            msg = (f'etype={self.element_name} ({self.element_type}) '
                   f'{self.table_name_str}-COMP-random-numwide={self.num_wide} '
                   f'numwide_real=11 numwide_imag=9 result_type={result_type}')
            return self._not_implemented_or_skip(data, ndata, msg), None, None

        else:
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
            raise RuntimeError(self.code_information())
            #msg = self.code_information()
            #msg = (f'etype={self.element_name} ({self.element_type}) '
                   #f'{self.table_name_str}-COMP-random-numwide={self.num_wide} '
                   #f'numwide_real=11 numwide_imag=9 result_type={result_type}')
            #return self._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oes_ctriax6(self, data, ndata, dt, is_magnitude_phase,
                     result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 53 : CTRIAX6
        """
        n = 0
        if self.is_stress:
            result_name = prefix + 'ctriax_stress' + postfix
        else:
            result_name = prefix + 'ctriax_strain' + postfix

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)

        slot = self.get_result(result_name)
        if result_type == 0 and self.num_wide == 33: # real
            if self.is_stress:
                obj_vector_real = RealTriaxStressArray
            else:
                obj_vector_real = RealTriaxStrainArray

            ntotal = 132 * self.factor  # (1+8*4)*4 = 33*4 = 132
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                self._data_factor = 4
                return nelements * self.num_wide * 4, None, None

            obj = self.obj
            nnodes_all = 4
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * self.num_wide * 4

                itotal = obj.itotal
                itotal2 = itotal + nelements * nnodes_all
                ielement = obj.ielement
                ielement2 = ielement + nelements

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 33)
                floats1 = floats[:, 1:].reshape(nelements * nnodes_all, 8).copy()

                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 33).copy()
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
                n = oes_ctriax6_real_33(self, data, obj,
                                        nelements, ntotal, dt)

        elif result_type == 1 and self.num_wide == 37: # imag
            # TODO: vectorize object
            if self.is_stress:
                #print('self.element_type', self.element_type)
                #print('self.element_name', self.element_name)
                #raise NotImplementedError('ComplexTriaxStressArray')
                obj_vector_complex = ComplexTriaxStressArray
            else:
                raise NotImplementedError('ComplexTriaxStrainArray')
                #obj_vector_complex = ComplexTriaxStrainArray

            num_wide = 1 + 4 * 9
            ntotal = num_wide * 4 * self.factor
            assert num_wide == self.num_wide, num_wide
            nelements = ndata // ntotal  # (1+8*4)*4 = 33*4 = 132
            leftover = ndata % ntotal
            assert leftover == 0, 'ntotal=%s nelements=%s leftover=%s' % (ntotal, nelements, leftover)

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_complex)

            if data is None:
                return ndata, None, None
            auto_return = False
            is_vectorized = False
            if auto_return:
                self._data_factor = 4
                return nelements * ntotal, None, None

            obj = self.obj
            nnodes_all = 4
            if self.use_vector and is_vectorized and 0:
                n = nelements * ntotal
                itotal = obj.itotal
                itotal2 = itotal + nelements * nnodes_all
                ielement = obj.ielement
                ielement2 = ielement + nelements

                numwide_imag = 37
                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, numwide_imag)
                floats1 = floats[:, 1:].reshape(nelements * nnodes_all, 9).copy()

                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, numwide_imag)
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
                n = oes_ctriax_complex_37(self, data, obj,
                                          nelements,
                                          is_magnitude_phase)
        else:  # pragma: no cover
            raise RuntimeError(self.code_information())
            #msg = self.code_information()
            #raise NotImplementedError(msg)
            #return self._not_implemented_or_skip(data, ndata, msg)
        return n, nelements, ntotal

    def _oes_cbush(self, data, ndata, dt, is_magnitude_phase,
                   result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 102 : CBUSH

        """
        n = 0
        if self.is_stress:
            result_name = prefix + 'cbush_stress' + postfix
        else:
            result_name = prefix + 'cbush_strain' + postfix

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)

        slot = self.get_result(result_name)
        if result_type in [0, 2] and self.num_wide == 7:  # real, random
            if self.is_stress:
                obj_vector_real = RealBushStressArray
            else:
                obj_vector_real = RealBushStrainArray

            assert self.num_wide == 7, "num_wide=%s not 7" % self.num_wide
            ntotal = 28 * self.factor # 4*7

            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return nelements * ntotal, None, None
            obj = self.obj

            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * ntotal

                istart = obj.ielement
                iend = istart + nelements
                obj._times[obj.itime] = dt

                self.obj_set_element(obj, istart, iend, data, nelements)

                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, 7)
                #[tx, ty, tz, rx, ry, rz]
                obj.data[obj.itime, istart:iend, :] = floats[:, 1:].copy()
            else:
                n = oes_cbush_real_7(self, data, obj,
                                     nelements, ntotal, dt)
        elif result_type == 1 and self.num_wide == 13:  # imag
            if self.is_stress:
                obj_complex = ComplexCBushStressArray
            else:
                obj_complex = ComplexCBushStrainArray

            ntotal = 52 * self.factor  # 4*13
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_complex)
            if auto_return:
                return nelements * ntotal, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, 13).copy()
                obj._times[obj.itime] = dt
                self.obj_set_element(obj, itotal, itotal2, data, nelements)

                isave1 = [1, 2, 3, 4, 5, 6]
                isave2 = [7, 8, 9, 10, 11, 12]
                real_imag = apply_mag_phase(floats, is_magnitude_phase, isave1, isave2)
                obj.data[obj.itime, itotal:itotal2, :] = real_imag

                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                n = oes_cbush_complex_13(self, data, obj,
                                         nelements, ntotal,
                                         is_magnitude_phase)
        else:  # pragma: no cover
            raise RuntimeError(self.code_information())
            #msg = self.code_information()
            #raise NotImplementedError(msg)
            #return self._not_implemented_or_skip(data, ndata, msg)
        return n, nelements, ntotal

    def _oes_cbush1d(self, data, ndata, dt, is_magnitude_phase,
                     result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 40 : CBUSH1D
        """
        n = 0
        if self.is_stress:
            result_name = prefix + 'cbush1d_stress_strain' + postfix
        else:
            result_name = prefix + 'cbush1d_stress_strain' + postfix

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)
        if result_type == 0 and self.num_wide == 8:  # real
            if self.is_stress:
                obj_vector_real = RealBush1DStressArray
            else:
                #self.create_transient_object(self.cbush1d_stress_strain, Bush1DStrain)  # undefined
                raise NotImplementedError('cbush1d_stress_strain; numwide=8')

            ntotal = 32 * self.factor # 4*8
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return nelements * ntotal, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * ntotal

                itotal = obj.itotal
                itotal2 = itotal + nelements
                itime = obj.itime
                obj._times[itime] = dt

                if 1: #obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 8).copy()
                    eids = ints[:, 0] // 10
                    fail = ints[:, 7]
                    obj.element[itotal:itotal2] = eids
                    obj.is_failed[itime, itotal:itotal2, 0] = fail

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 8)
                #[xxx, fe, ue, ve, ao, ae, ep, xxx]
                obj.data[itime, itotal:itotal2, :] = floats[:, 1:7].copy()

                obj.ielement = itotal2
                obj.itotal = itotal2
            else:
                struct1 = Struct(self._endian + self._analysis_code_fmt + b'6fi')
                for unused_i in range(nelements):
                    edata = data[n:n + ntotal]
                    out = struct1.unpack(edata)  # num_wide=25
                    if self.is_debug_file:
                        self.binary_debug.write('CBUSH1D-40 - %s\n' % (str(out)))
                    (eid_device, fe, ue, ve, ao, ae, ep, fail) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)

                    # axial_force, axial_displacement, axial_velocity, axial_stress,
                    # axial_strain, plastic_strain, is_failed
                    obj.add_sort1(dt, eid, fe, ue, ve, ao, ae, ep, fail)
                    n += ntotal
        elif result_type == 1 and self.num_wide == 9:  # imag
            # TODO: vectorize object
            ntotal = 36 * self.factor  # 4*9
            nelements = ndata // ntotal

            if self.is_stress:
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, ComplexCBush1DStressArray)
            else:
                raise NotImplementedError('self.cbush1d_stress_strain; complex strain')

            if auto_return:
                return nelements * self.num_wide * 4, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * self.num_wide * 4

                itotal = obj.itotal
                itotal2 = itotal + nelements
                itime = obj.itime
                obj._times[itime] = dt
                self.obj_set_element(obj, itotal, itotal2, data, nelements)

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 9).copy()
                #[fer, uer, aor, aer,
                # fei, uei, aoi, aei]
                isave1 = [1, 3, 5, 7]
                isave2 = [2, 4, 6, 8]
                real_imag = apply_mag_phase(floats, is_magnitude_phase, isave1, isave2)
                obj.data[obj.itime, itotal:itotal2, :] = real_imag

                obj.ielement = itotal2
                obj.itotal = itotal2
            else:
                struct1 = Struct(self._endian + self._analysis_code_fmt + b'8f')
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]

                    out = struct1.unpack(edata)  # num_wide=25
                    (eid_device,
                     fer, uer, aor, aer,
                     fei, uei, aoi, aei) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)

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
                    obj.add_new_eid(self.element_type, dt, eid, fe, ue, ao, ae)
        else:  # pragma: no cover
            raise RuntimeError(self.code_information())
            #msg = self.code_information()
            #raise NotImplementedError(msg)
            #return self._not_implemented_or_skip(data, ndata, msg)
        return n, nelements, ntotal

    def _oes_crod_nonlinear(self, data, ndata, dt, unused_is_magnitude_phase,
                            result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 87 : CTUBENL
         - 89 : RODNL
         - 92 : CONRODNL
        """
        n = 0
        #prefix = 'nonlinear_'
        if self.is_stress:
            if self.element_type == 87:
                result_name = prefix + 'ctube_stress' + postfix
                name = 'CTUBENL-87'
            elif self.element_type == 89:
                result_name = prefix + 'crod_stress' + postfix
                name = 'RODNL-89'
            elif self.element_type == 92:
                result_name = prefix + 'conrod_stress' + postfix
                name = 'CONRODNL-92'
            else:  # pragma: no cover
                raise RuntimeError(self.code_information())
        else:
            if self.element_type == 87:
                result_name = prefix + 'ctube_strain' + postfix
                name = 'CTUBENL-87'
            elif self.element_type == 89:
                result_name = prefix + 'crod_strain' + postfix
                name = 'RODNL-89'
            elif self.element_type == 92:
                result_name = prefix + 'conrod_strain' + postfix
                name = 'CONRODNL-92'
            else:  # pragma: no cover
                raise RuntimeError(self.code_information())

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)
        if result_type == 0 and self.num_wide == 7:  # real
            ntotal = 28 * self.factor #  7*4 = 28
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, RealNonlinearRodArray)
            if auto_return:
                return nelements * ntotal, None, None

            obj = self.obj
            #if self.is_debug_file:
                #self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                #self.binary_debug.write('  element1 = [eid_device, layer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)]\n')
                #self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * self.num_wide * 4
                istart = obj.itotal
                iend = istart + nelements
                obj._times[obj.itime] = dt

                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype8).reshape(nelements, 7).copy()
                    eids = ints[:, 0] // 10
                    obj.element[istart:iend] = eids
                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, 7)
                #[axial_stress, equiv_stress, total_strain,
                # eff_plastic_creep_strain, eff_creep_strain, linear_torsional_stresss]
                obj.data[obj.itime, istart:iend, :] = floats[:, 1:].copy()
            else:
                struct1 = Struct(self._endian + mapfmt(self._analysis_code_fmt + b'6f', self.size))  # 1+6=7
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]
                    out = struct1.unpack(edata)

                    (eid_device, axial_stress, equiv_stress, total_strain,
                     eff_plastic_creep_strain, eff_creep_strain, linear_torsional_stresss) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    if self.is_debug_file:
                        self.binary_debug.write('%s - %s\n' % (name, str(out)))
                    obj.add_sort1(dt, eid, axial_stress, equiv_stress, total_strain,
                                  eff_plastic_creep_strain, eff_creep_strain, linear_torsional_stresss)
                    n += ntotal
        else:  # pragma: no cover
            raise RuntimeError(self.code_information())
        return n, nelements, ntotal

    def _oes_celas_nonlinear(self, data, ndata, dt, unused_is_magnitude_phase,
                             result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 224 : CELAS1
         - 226 : CELAS3

        """
        # 224-CELAS1
        # 225-CELAS3
        # NonlinearSpringStress
        n = 0
        numwide_real = 3
        if self.is_stress:
            if self.element_type == 224:
                result_name = prefix + 'celas1_stress' + postfix # nonlinear_
            elif self.element_type == 225:
                result_name = prefix + 'celas3_stress' + postfix # nonlinear_
        else:
            raise NotImplementedError('NonlinearSpringStrain')

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)

        slot = self.get_result(result_name)
        if result_type == 0 and self.num_wide == numwide_real:
            assert self.num_wide == 3, "num_wide=%s not 3" % self.num_wide
            ntotal = 12 * self.factor  # 4*3
            nelements = ndata // ntotal

            if self.is_stress:
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, RealNonlinearSpringStressArray)
            else:
                raise NotImplementedError('NonlinearSpringStrainArray') # undefined

            if auto_return:
                return nelements * self.num_wide * 4, None, None
            obj = self.obj

            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * 4 * self.num_wide
                unused_itotal = obj.ielement
                ielement = obj.ielement
                ielement2 = obj.ielement + nelements
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, numwide_real).copy()
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[ielement:ielement2] = eids

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, numwide_real)

                #[force, stress]
                obj.data[obj.itime, ielement:ielement2, :] = floats[:, 1:].copy()
                obj.itotal = ielement2
                obj.ielement = ielement2
            else:
                struct1 = Struct(self._endian + self._analysis_code_fmt + b'2f')
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]
                    out = struct1.unpack(edata)  # num_wide=3
                    (eid_device, force, stress) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    if self.is_debug_file:
                        self.binary_debug.write('%s-%s - %s\n' % (self.element_name, self.element_type, str(out)))
                    obj.add_sort1(dt, eid, force, stress)
                    n += ntotal
        else:  # pragma: no cover
            raise RuntimeError(self.code_information())
        return n, nelements, ntotal

    def _oes_cbush_nonlinear(self, data, ndata, dt, unused_is_magnitude_phase,
                             result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 226 : CBUSHNL
        """
        n = 0
        if self.is_stress:
            if self.element_type == 226:
                result_name = prefix + 'cbush_force_stress_strain' + postfix
                name = 'CBUSHNL-226'
            else:  # pragma: no cover
                raise RuntimeError(self.code_information())
        else:
            if self.element_type == 226:
                result_name = prefix + 'nonlinear_cbush_strain' + postfix
                name = 'CBUSHNL-226'
            else:  # pragma: no cover
                raise RuntimeError(self.code_information())

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)
        if result_type == 0 and self.num_wide == 19:  # real
            ntotal = 76 * self.factor  #  19*4 = 76
            nelements = ndata // ntotal
            assert ndata % ntotal == 0
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, RealNonlinearBushArray)
            if auto_return:
                return nelements * self.num_wide * 4, None, None

            obj = self.obj
            #if self.is_debug_file:
                #self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                #self.binary_debug.write('  element1 = [eid_device, layer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)]\n')
                #self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * self.num_wide * 4
                istart = obj.itotal
                iend = istart + nelements
                obj._times[obj.itime] = dt

                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 19).copy()
                    eids = ints[:, 0] // 10
                    obj.element[istart:iend] = eids
                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 19)
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
                struct1 = Struct(self._endian + self._analysis_code_fmt + b'18f')
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]
                    out = struct1.unpack(edata)

                    (eid_device, fx, fy, fz, otx, oty, otz, etx, ety, etz,
                                 mx, my, mz, orx, ory, orz, erx, ery, erz) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    if self.is_debug_file:
                        self.binary_debug.write('%s - %s\n' % (name, str(out)))
                    obj.add_sort1(dt, eid, fx, fy, fz, otx, oty, otz, etx, ety, etz,
                                  mx, my, mz, orx, ory, orz, erx, ery, erz)
                    n += ntotal
        else:  # pragma: no cover
            raise RuntimeError(self.code_information())
        return n, nelements, ntotal

    def _oes_cbend(self, data, ndata, dt, is_magnitude_phase,
                   result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 69 : CBEND

        """
        if self.is_stress:
            result_name = prefix + 'cbend_stress' + postfix
            obj_vector_real = RealBendStressArray
            obj_vector_complex = ComplexBendStressArray
            obj_vector_random = RandomBendStressArray
        else:
            result_name = prefix + 'cbend_strain' + postfix
            obj_vector_real = RealBendStrainArray
            obj_vector_complex = ComplexBendStrainArray
            obj_vector_random = RandomBendStrainArray

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)

        #print(self.code_information())
        if result_type == 0 and self.num_wide == 21:  # real
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
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return nelements * ntotal, None, None

            obj = self.obj
            assert obj is not None
            if self.use_vector and is_vectorized and self.sort_method == 1 and 0:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 4)
                itime = obj.itime
                obj._times[itime] = dt
                if itime == 0:
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 4)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                #[max_strain, avg_strain, margin]
                obj.data[itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                n = oes_cbend_real_21(self, data, obj,
                                      nelements, ntotal, dt)

            #msg = ''
            #if self.read_mode == 2:
                #msg = self.code_information()
            #n = self._not_implemented_or_skip(data, ndata, msg)
            #return n, None, None
        elif result_type == 1 and self.num_wide == 21:  # complex
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

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_complex)
            if auto_return:
                return nelements * self.num_wide * 4, None, None

            obj = self.obj
            assert obj is not None
            if self.use_vector and is_vectorized and self.sort_method == 1 and 0:
                n = nelements * 4 * self.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 4)
                itime = obj.itime
                obj._times[itime] = dt

                if itime == 0:
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 4)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                #[max_strain, avg_strain, margin]
                obj.data[itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                ntotali = 40
                struct1 = Struct(self._endian + self._analysis_code_fmt)
                struct2 = Struct(self._endian + b'i9f')

                for unused_i in range(nelements):
                    edata = data[n:n + 4]
                    eid_device, = struct1.unpack(edata)
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)

                    n += 4
                    for unused_j in range(2):
                        edata = data[n:n + ntotali]
                        out = struct2.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('BEND-69 - eid=%s %s\n' % (eid, str(out)))
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
                        obj.add_sort1(dt, eid, grid, angle, sc, sd, se, sf)
                        n += ntotali

        elif result_type == 2 and self.num_wide == 13:
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
            #if self.table_name != "OESPSD2":
                #msg = ''
                #if self.read_mode == 2:
                    #msg = self.code_information()
                #n = self._not_implemented_or_skip(data, ndata, msg)
                #return n, None, None

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_random)
            if auto_return:
                return nelements * self.num_wide * 4, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1 and 0:
                n = nelements * 4 * self.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 4)
                itime = obj.itime
                obj._times[itime] = dt
                if itime == 0:
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 4)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                #[max_strain, avg_strain, margin]
                obj.data[itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                ntotali = 24
                struct1 = Struct(self._endian + self._analysis_code_fmt)
                struct2 = Struct(self._endian + b'i5f')

                for unused_i in range(nelements):
                    edata = data[n:n + 4]
                    #self.show_data(edata)
                    eid_device, = struct1.unpack(edata)
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)

                    n += 4
                    for unused_i in range(2):
                        edata = data[n:n + ntotali]
                        out = struct2.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('BEND-69 - eid=%s dt=%s %s\n' % (eid, dt, str(out)))
                        #print('BEND-69 - eid=%s dt=%s %s\n' % (eid, dt, str(out)))

                        (grid, angle, sc, sd, se, sf) = out
                        obj.add_sort1(dt, eid, grid, angle, sc, sd, se, sf)
                        n += ntotali

        else:  # pragma: no cover
            raise RuntimeError(self.code_information())
        return n, nelements, ntotal

    def _oes_cgap_nonlinear(self, data, ndata, dt, is_magnitude_phase,
                            result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 86 : GAPNL
        """
        n = 0
        if self.is_stress:
            result_name = prefix + 'cgap_stress' + postfix # nonlinear_
        else:
            result_name = prefix + 'cgap_strain' + postfix # nonlinear_

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)

        if result_type == 0 and self.num_wide == 11:  # real?
            if self.is_stress:
                obj_vector_real = NonlinearGapStressArray
            else:
                raise NotImplementedError('NonlinearGapStrain')

            ntotal = 44 * self.factor # 4*11
            nelements = ndata // ntotal
            assert ndata % ntotal == 0
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return nelements * ntotal, None, None

            obj = self.obj
            if self.use_vector and is_vectorized and self.sort_method == 1:
                n = nelements * ntotal

                ielement = obj.ielement
                ielement2 = ielement + nelements
                obj._times[obj.itime] = dt

                self.obj_set_element(obj, ielement, ielement2, data, nelements)

                #if obj.itime == 0:
                    #ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 11).copy()
                    #eids = ints[:, 0] // 10
                    #obj.element[ielement:ielement2] = eids

                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, 11)
                # skipping [form1, form2]
                #[cpx, shy, shz, au, shv, shw, slv, slp]
                obj.data[obj.itime, ielement:ielement2, :] = floats[:, 1:9].copy()
            else:
                if self.size == 4:
                    struct1 = Struct(self._endian + self._analysis_code_fmt + b'8f4s4s')
                else:
                    struct1 = Struct(self._endian + mapfmt(self._analysis_code_fmt, self.size) + b'8d8s8s')
                for unused_i in range(nelements):
                    edata = data[n:n + ntotal]
                    out = struct1.unpack(edata)  # num_wide=25
                    (eid_device, cpx, shy, shz, au, shv, shw, slv, slp, form1, form2) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)
                    if self.is_debug_file:
                        self.binary_debug.write('CGAPNL-86 - %s\n' % str(out))
                    obj.add_sort1(dt, eid, cpx, shy, shz, au, shv, shw, slv, slp, form1, form2)
                    n += ntotal
        else:  # pragma: no cover
            raise RuntimeError(self.code_information())
        return n, nelements, ntotal

    def _oes_cbeam_nonlinear(self, data, ndata, dt, is_magnitude_phase,
                             result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 94 : BEAMNL

        """
        n = 0
        numwide_real = 51
        numwide_random = 0

        if self.is_stress:
            result_name = prefix + 'cbeam_stress' + postfix
        else:
            result_name = prefix + 'cbeam_strain' + postfix

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)

        if result_type == 0 and self.num_wide == numwide_real:
            msg = result_name
            if self.is_stress:
                obj_vector_real = RealNonlinearBeamStressArray
            else:
                raise NotImplementedError('Nonlinear CBEAM Strain...this should never happen')

            ntotal = numwide_real * 4 * self.factor  # 204
            nelements = ndata // ntotal

            nlayers = nelements * 8
            auto_return, is_vectorized = self._create_oes_object4(
                nlayers, result_name, slot, obj_vector_real)
            if auto_return:
                self._data_factor = 8
                return ndata, None, None
            obj = self.obj
            if self.is_debug_file:
                self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                #self.binary_debug.write('  #elementi = [eid_device, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,\n')
                #self.binary_debug.write('                           s1b, s2b, s3b, s4b, smaxb, sminb,        MSc]\n')
                #self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if self.size == 4:
                struct1 = Struct(self._endian + b'2i 4s5f 4s5f 4s5f 4s5f i 4s5f 4s5f 4s5f 4s5f')  # 2 + 6*8 + 1 = 51
            else:
                assert self.size == 8, self.size
                struct1 = Struct(self._endian + b'2q 8s5d 8s5d 8s5d 8s5d q 8s5d 8s5d 8s5d 8s5d')  # 2 + 6*8 + 1 = 51

            for unused_i in range(nelements):  # num_wide=51
                edata = data[n:n + ntotal]
                out = struct1.unpack(edata)

                if self.is_debug_file:
                    self.binary_debug.write('BEAMNL-94 - %s\n' % str(out))

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
                    eid_device, self.nonlinear_factor, self.sort_method)
                obj.add_new_eid_sort1(dt, eid, *out[1:])
                n += ntotal

        elif result_type == 2 and self.num_wide == numwide_random:  # random
            msg = self.code_information()
            raise NotImplementedError(msg)
            #return self._not_implemented_or_skip(data, ndata, msg)
        else:  # pragma: no cover
            raise RuntimeError(self.code_information())
        return n, nelements, ntotal

    def _oes_cbar_100(self, data, ndata, dt, is_magnitude_phase,
                      result_type, prefix, postfix):
        """
        reads stress/strain for element type:
         - 100 : BARS
        """
        n = 0
        if self.is_stress:
            result_name = prefix + 'cbar_stress_10nodes' + postfix
        else:
            result_name = prefix + 'cbar_strain_10nodes' + postfix

        if self._results.is_not_saved(result_name):
            return ndata, None, None
        self._results._found_result(result_name)
        slot = self.get_result(result_name)

        if result_type == 0 and self.num_wide == 10:  # real
            if self.is_stress:
                obj_vector_real = RealBar10NodesStressArray
            else:
                obj_vector_real = RealBar10NodesStrainArray

            ntotal = 10 * 4 * self.factor
            nelements = ndata // ntotal

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return ndata, None, None

            if self.is_debug_file:
                self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                self.binary_debug.write('  #elementi = [eid_device, sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS]\n')
                self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)
            obj = self.obj

            if self.use_vector and is_vectorized and self.sort_method == 1:
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

                floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, 10)
                #[sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS]
                obj.data[obj.itime, istart:iend, :] = floats[:, 1:].copy()
            else:
                n = oes_cbar100_real_10(self, data, obj, nelements, ntotal, dt)
        else:  # pragma: no cover
            raise RuntimeError(self.code_information())
        return n, nelements, ntotal

    def _oes_vu_solid(self, data, ndata, dt, unused_is_magnitude_phase,
                      unused_stress_name, result_type, unused_prefix, unused_postfix):
        # TODO: vectorize
        if self.read_mode == 1:
            return ndata, None, None

        # 145-VUHEXA  (8 nodes)
        # 146-VUPENTA (6 nodes)
        # 147-VUTETRA (4 nodes)
        self.log.warning(f'skipping oes_vu_solid; {self.element_name}-{self.element_type} '
                         f'(numwide={self.num_wide})')
        if self.element_type == 147:
            etype = 'VUTETRA'
            nnodes = 4
            #numwide_a = 2 + (14 - 2) * nnodes  # 50
            #numwide_b = 2 + (9 - 2) * nnodes  # 30
            #numwide_c = 2 + 13 * nnodes  # 54
        elif self.element_type == 146:
            etype = 'VUPENTA'
            nnodes = 6
        elif self.element_type == 145:
            etype = 'VUHEXA'
            nnodes = 8
            # numwide=145
        else:  # pragma: no cover
            raise RuntimeError(self.code_information())

        #num_wideA = 2 + 12 * nnodes
        #ntotal = 8 + 48 * nnodes

        n = 0
        if self.format_code == 1: # real
            # assuming TETRA...
            # TODO: vectorize
            numwide_a = 2 + (14 - 2) * nnodes  # 50
            numwide_b = 2 + (9 - 2) * nnodes  # 30
            numwide_c = 2 + 13 * nnodes  # 54
            if self.num_wide == numwide_a:
                ntotal = numwide_a * 4
                s1 = self.struct_2i
                s2 = Struct(self._endian + b'i11f')
                nelements = ndata // ntotal  # 2+16*9 = 146 -> 146*4 = 584

                ntotal1 = 8 * self.factor
                ntotal2 = 48 * self.factor
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal1]
                    out = s1.unpack(edata)
                    (eid_device, parent_id) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)

                    for unused_j in range(nnodes):
                        edata = data[n:n+ntotal2]
                        out = s2.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('%s-%s - %s\n' % (etype, self.element_type, str(out)))
                        assert len(out) == 12
                        (grid, xnorm, ynorm, znorm, txy, tyz, txz,
                         prin1, prin2, prin3, smean, vono_roct) = out
                        del grid, xnorm, ynorm, znorm, txy, tyz, txz
                        del prin1, prin2, prin3, smean, vono_roct
                n = ndata
            elif self.num_wide == numwide_b:
                ntotal = numwide_b * 4
                nelements = ndata // ntotal
                n = nelements * ntotal
            elif self.num_wide == numwide_c:
                ntotal = numwide_c * 4
                nelements = ndata // ntotal
                n = nelements * ntotal
            else:
                msg = f'numwide={self.num_wide} A={numwide_a} B={numwide_b} C={numwide_c}'
                raise RuntimeError(self.code_information() + msg)
        elif result_type == 1:  # complex
            if self.num_wide == 54:  # VU-TETRA
                pass
            elif self.num_wide == 80:  # VU-PENTA
                pass
            elif self.num_wide == 106:  # VU-HEXA
                pass
            else:
                raise RuntimeError(self.code_information())
            # C:\NASA\m4\formats\git\examples\move_tpl\pe108p04.op2
            # analysis_code = 5   Frequency
            # element_type  = 146 VUPENTA
            # num_wide      = 80
        else:
            raise RuntimeError(self.code_information())
            #raise RuntimeError(self.code_information())
            #msg = self.code_information()
            #raise RuntimeError(msg)
            #return self._not_implemented_or_skip(data, ndata, msg), None, None
        return ndata, None, None
        return n, nelements, ntotal

    def _oes_hyperelastic_quad(self, data, ndata, dt, unused_is_magnitude_phase,
                               result_type, prefix, postfix):
        """
        139-QUAD4FD
        """
        #if self.is_stress:
        result_name = prefix + 'hyperelastic_cquad4_strain' + postfix
        if self._results.is_not_saved(result_name):
            return ndata, None, None

        if result_type == 0 and self.num_wide == 30:
            obj_vector_real = HyperelasticQuadArray

            self._results._found_result(result_name)
            slot = self.get_result(result_name)
            #self.create_transient_object(result_name, slot, obj_vector_real)

            ntotal = 120 * self.factor # 36+28*3
            nelements = ndata // ntotal

            #print(self.code_information())
            #print(self.table_name_str)
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)

            if auto_return:
                self._data_factor = 4  # number of "layers" for an element
                return nelements * ntotal, None, None
                #return ndata, None, None

            #if self.is_debug_file:
                #self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                ##self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                #self.binary_debug.write('  #elementi = [eid_device, sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS]\n')
                #self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)
            obj = self.obj

            if self.use_vector and is_vectorized and self.sort_method == 1:
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
                ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 30).copy()
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
                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 30)[:, 2:].copy()
                floats2 = floats.reshape(nelements * 4, 7)
                #[oxx, oyy, txy, angle, majorp, minorp]
                obj.data[obj.itime, istart:iend, :] = floats2[:, 1:]
            else:
                n = 0
                # (2 + 7*4)*4 = 30*4 = 120
                ntotal1 = 36 * self.factor  # 4*9
                ntotal2 = 28 * self.factor  # 4*7
                s1 = Struct(self._endian + self._analysis_code_fmt + b'4s i6f')  # 1 + 4+1+6 = 12
                s2 = Struct(self._endian + b'i6f')
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal1]
                    out = s1.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('CQUAD4FD-139A- %s\n' % (str(out)))

                    (eid_device, etype, nid, sx, sy, sxy, angle, smj, smi) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, self.nonlinear_factor, self.sort_method)

                    obj._add_new_eid_sort1(dt, eid, etype, nid, sx, sy, sxy, angle, smj, smi)
                    n += ntotal1

                    for unused_i in range(3):  # TODO: why is this not 4?
                        edata = data[n:n + ntotal2]
                        out = s2.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('               %s\n' % (str(out)))
                        (nid, sx, sy, sxy, angle, smj, smi) = out
                        obj._add_sort1(dt, eid, etype, nid, sx, sy, sxy, angle, smj, smi)
                        n += ntotal2
        else:
            raise RuntimeError(self.code_information())
            #msg = 'numwide=%s element_num=%s etype=%s' % (
                #self.num_wide, self.element_type, self.element_name)
            #return self._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oes_vu_quad(self, data, ndata, unused_dt, unused_is_magnitude_phase,
                     result_type, unused_prefix, unused_postfix):
        """Adds an adaptive VUQUAD"""
        n = 0
        if self.element_type == 189:  # VQUAD
            if self.read_mode == 1:
                return ndata, None, None
            #ntotal = 440  # 6+(33-7)*4 =  -> 110*4 = 440
            nnodes = 4    # 4 corner points + 1 centroid
            etype = 'VUQUAD4'
        #elif self.element_type == 144:  # CQUAD4
            #ntotal = 440  # 6+(33-7)*4 =  -> 110*4 = 440
            #nnodes = 4    # 4 corner points
            #etype = 'CQUAD4'
        #elif self.element_type == 64:  # CQUAD8
            #ntotal = 348  # 2+17*5 = 87 -> 87*4 = 348
            #nnodes = 4    # centroid + 4 corner points
            #etype = 'CQUAD8'
        #elif self.element_type == 82:  # CQUADR
            #ntotal = 348  # 2+17*5 = 87 -> 87*4 = 348
            #nnodes = 4    # centroid + 4 corner points
            #etype = 'CQUADR'
        #elif self.element_type == 75:  # CTRIA6
            #ntotal = 280  # 2+17*3 = 70 -> 70*4 = 280
            #nnodes = 3    # centroid + 3 corner points
            #etype = 'CTRIA6'
        #elif self.element_type == 70:  # CTRIAR
            #ntotal = 280  # 2+17*3 = 70 -> 70*4 = 280
            #nnodes = 3    # centroid + 3 corner points
            #etype = 'CTRIAR'
        else:
            return self._not_implemented_or_skip(data, ndata, self.code_information()), None, None

        numwide_real = 6 + (23 - 6) * nnodes  # real???
        numwide_imag = 6 + (33 - 6) * nnodes  # imag???

        if result_type == 0 and self.num_wide == numwide_real:  # real???
            ntotal = numwide_real * 4 * self.factor
            if self.size == 4:
                s2 = Struct(self._endian + b'3i4s2i')
            else:
                s2 = Struct(self._endian + b'3q8s2q')
            s3 = Struct(mapfmt(self._endian + b'i16f', self.size))
            nelements = ndata // ntotal
            ntotal1 = 24 * self.factor
            ntotal2 = 68 * self.factor
            for unused_i in range(nelements):
                out = s2.unpack(data[n:n + ntotal1])
                (eid_device, unused_parent, coord, icord, unused_theta, unused_itype) = out
                n += ntotal1
                eid, dt = get_eid_dt_from_eid_device(
                    eid_device, self.nonlinear_factor, self.sort_method)
                edata = data[n:n + ntotal2]
                out = s3.unpack(edata)  # len=17*4
                n += ntotal2

                if self.is_debug_file:
                    self.binary_debug.write('%s-%s - %s\n' % (etype, self.element_type, str(out)))

                #obj.add_new_node_sort1(dt, eid, parent, coord, icord, theta, itype)
                #obj.add_new_eid_sort1(eType, dt, eid, parent, coord, icord, theta, itype)
                for unused_node_id in range(nnodes - 1):  # nodes pts
                    edata = data[n:n + ntotal2]
                    n += ntotal2
                    out = s3.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('              %s\n' % (str(out)))

                    (unused_vuid, unused_dummy, unused_dummy2,
                     unused_msx, unused_msy, unused_mxy,
                     unused_dummy3, unused_dummy4, unused_dummy5,
                     unused_bcx, unused_bcy, unused_bcxy, unused_tyz, unused_tzx,
                     unused_dummy6, unused_dummy7, unused_dummy8) = out
                    #obj.add_sort1(vuid, dummy, dummy2, msx, msy, mxy,
                                   #dummy3, dummy4, dummy5,
                                   #bcx, bcy, bcxy, tyz, tzx,
                                   #dummy6, dummy7, dummy8)
            nelements = None
            ntotal = None
            self.log.warning(f'skipping oes_vu_quad; {self.element_name}-{self.element_type} '
                             f'(numwide={self.num_wide})')
        elif result_type == 1 and self.num_wide == numwide_imag:
            ntotal = (numwide_imag * 4) * self.factor
            nelements = ndata // ntotal
            n = nelements * ntotal
            nelements = None
            ntotal = None
            self.log.warning(f'skipping oes_vu_quad; {self.element_name}-{self.element_type} '
                             f'(numwide={self.num_wide})')
        else:
            raise RuntimeError(self.code_information())
            #msg = 'numwide=%s' % self.num_wide
            #return self._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _oes_plate_stress_34(self, data, ndata, unused_dt, unused_is_magnitude_phase,
                             unused_stress_name, unused_prefix, unused_postfix):
        """
        271-CPLSTN3
        275-CPLSTS3
        """
        msg = self.code_information()
        return self._not_implemented_or_skip(data, ndata, msg), None, None
        #if self.element_type == 271:
            #result_name = 'cplstn3'
            #unused_nnodes = 1
            #ntotal = 4 * 6
        #elif self.element_type == 275:
            #result_name = 'cplsts3'
            #unused_nnodes = 1
            #ntotal = 4 * 6
        #else:  # pragma: no cover
            #raise RuntimeError(self.code_information())
        #if self.is_stress:
            #obj_vector_real = RealCPLSTRNPlateStressArray
            #result_name += '_stress'
        #else:
            #obj_vector_real = RealCPLSTRNPlateStrainArray
            #result_name += '_strain'

        #numwide_real = ntotal // 4
        #if self.format_code == 1 and self.num_wide == numwide_real:
            ##ntotal = 4 * (1 + 6 * (nnodes))
            #nelements = ndata // ntotal

            ##self._data_factor = 10  # TODO: why is this 10?
            #if self.is_stress:
                #obj_vector_real = RealCPLSTRNPlateStressArray
                ##result_name = 'cplstn3_stress'
            #else:
                #obj_vector_real = RealCPLSTRNPlateStressArray
                ##result_name = 'cplstn3_strain'
            #slot = self.get_result(result_name)

            #auto_return, is_vectorized = self._create_oes_object4(
                #nelements, result_name, slot, obj_vector_real)
            #if auto_return:
                #return nelements * self.num_wide * 4, None, None

            #obj = self.obj
            ##if self.use_vector and is_vectorized and self.sort_method == 1:
            ##n = nelements * self.num_wide * 4

            #istart = obj.itotal
            #iend = istart + nelements
            #obj._times[obj.itime] = dt

            #self.obj_set_element(obj, istart, iend, data, nelements)
            #floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, numwide_real)
            #results = floats[:, 1:].copy()
            ##print('results.shape', results.shape)

            ##[oxx, oyy, ozz, txy, ovm]
            #obj.data[obj.itime, istart:iend, :] = results
        #else:
            #msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
            #return self._not_implemented_or_skip(data, ndata, msg), None, None

    def _oes_plate_stress_68(self, data, ndata, unused_dt, unused_is_magnitude_phase,
                             unused_stress_name, unused_prefix, unused_postfix):
        # 276-CPLSTS4
        # 277-CPLSTS6
        # 278-CPLSTS8
        msg = self.code_information()
        return self._not_implemented_or_skip(data, ndata, msg), None, None
        #if self.element_type == 276:
            #result_name = 'cplsts4'
            #nnodes = 5  # 4 + 1
            #ntotal = 4 * 32
        #elif self.element_type == 277:
            #result_name = 'cplsts6'
            #nnodes = 4
            #ntotal = 4 * 26
        #elif self.element_type == 278:
            #result_name = 'cplsts8'
            #nnodes = 5
            #ntotal = 4 * 32
        #else:
            #raise RuntimeError(self.code_information())

        #if self.is_stress:
            #obj_vector_real = RealCPLSTRNPlateStressArray
            #result_name += '_stress'
        #else:
            #obj_vector_real = RealCPLSTRNPlateStrainArray
            #result_name += '_strain'

        #numwide_real = 2 + 6 * (nnodes)
        #assert ntotal // 4 == numwide_real, 'notal/4=%s numwide_real=%s\n%s' % (
            #ntotal // 4, numwide_real, self.code_information())

        #ntotal = numwide_real * 4
        #if self.format_code == 1 and self.num_wide == numwide_real:
            #nelements = ndata // ntotal

            ##self._data_factor = 10  # TODO: why is this 10?
            #if self.is_stress:
                #obj_vector_real = RealCPLSTRNPlateStressArray
                ##result_name = 'cplstn3_stress'
            #else:
                #obj_vector_real = RealCPLSTRNPlateStressArray
                ##result_name = 'cplstn3_strain'
            #slot = getattr(self, result_name)

            #nlayers = nelements * nnodes
            #auto_return, is_vectorized = self._create_oes_object4(
                #nlayers, result_name, slot, obj_vector_real)
            #if auto_return:
                #self._data_factor = nnodes
                #return nelements * self.num_wide * 4

            #obj = self.obj
            ##if self.use_vector and is_vectorized and self.sort_method == 1:
            #n = nlayers * self.num_wide * 4

            #istart = obj.itotal
            #iend = istart + nlayers
            #obj._times[obj.itime] = dt

            #if obj.itime == 0:
                #print(frombuffer(data, dtype=self.idtype).size)
                #print('nelements=%s numwide=%s' % (nelements, numwide_real))
                #ints = frombuffer(data, dtype=self.idtype).reshape(nelements, numwide_real)
                #eids = ints[:, 0] // 10
                ##obj.element[istart:iend] = eids

            #floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, numwide_real).copy()
            #print('floats[:, 2:].shape', floats[:, 2:].shape)
            #print('nnelements=%s nnodes=%s numwide//nodes=%s' % (nelements, nnodes, (numwide_real-2) / nnodes))
            #results = floats[:, 2:].reshape(nelements, nnodes * 6)

            ##[oxx, oyy, ozz, txy, ovm]
            #obj.data[obj.itime, istart:iend, :] = results
        #else:
            #msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
            #return self._not_implemented_or_skip(data, ndata, msg)

    def obj_set_element(self, obj, ielement, ielement2, data, nelements):
        if obj.itime == 0:
            ints = frombuffer(data, dtype=self.idtype8).reshape(nelements, self.num_wide).copy()
            eids = ints[:, 0] // 10
            assert eids.min() > 0, eids.min()
            obj.element[ielement:ielement2] = eids


def oes_cquad4_33_complex_vm_17(self, data: bytes,
                                obj: Union[ComplexPlateVMStressArray, ComplexPlateVMStrainArray],
                                nelements: int, ntotal: int,
                                is_magnitude_phase: bool) -> int:
    """
    OESVM1/2 - Table of element stresses or strains with von Mises
    OSTRVM1/2  for frequency response results.

    """
    #self.log.info(f'nelements = {nelements}')
    n = 0
    struct1 = Struct(mapfmt(self._endian + self._analysis_code_fmt + b'16f', self.size))
    cen = 0 # CEN/4
    add_sort_x = getattr(obj, 'add_sort' + str(self.sort_method))

    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device,
         fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i, von_mises1,
         fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i, von_mises2) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)
        if self.is_debug_file:
            self.binary_debug.write('  eid=%i C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))

        if is_magnitude_phase:
            sx1 = polar_to_real_imag(sx1r, sx1i)
            sx2 = polar_to_real_imag(sx2r, sx2i)
            sy1 = polar_to_real_imag(sy1r, sy1i)
            sy2 = polar_to_real_imag(sy2r, sy2i)
            txy1 = polar_to_real_imag(txy1r, txy1i)
            txy2 = polar_to_real_imag(txy2r, txy2i)
        else:
            sx1 = complex(sx1r, sx1i)
            sx2 = complex(sx2r, sx2i)
            sy1 = complex(sy1r, sy1i)
            sy2 = complex(sy2r, sy2i)
            txy1 = complex(txy1r, txy1i)
            txy2 = complex(txy2r, txy2i)
        #print(dt, eid, cen, sx1, sy1, txy1, max_shear1)
        #print(dt, eid, cen, sx2, sy2, txy2, max_shear2)
        add_sort_x(dt, eid, cen,
                   fd1, sx1, sy1, txy1, von_mises1,
                   fd2, sx2, sy2, txy2, von_mises2)
        n += ntotal
    return n

def oes_cbeam_real_111(self, data: bytes,
                       obj: Union[RealBeamStressArray, RealBeamStrainArray],
                       nelements: int, dt: Any) -> int:
    n = 0
    nnodes = 10  # 11-1
    n1 = 44 * self.factor
    n2 = 40 * self.factor
    fmt1 = mapfmt(self._endian + self._analysis_code_fmt + b'i9f', self.size)
    fmt2 = mapfmt(self._endian + b'i9f', self.size)
    s1 = Struct(fmt1)
    s2 = Struct(fmt2)
    for unused_i in range(nelements):
        edata = data[n:n+n1]
        n += n1

        out = s1.unpack(edata)
        eid_device = out[0]
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)
        if self.is_debug_file:
            self.binary_debug.write('CBEAM-2 - eid=%i out=%s\n' % (eid, str(out)))

        #(grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = out
        obj.add_new_eid_sort1(dt, eid, *out[1:])

        for unused_inode in range(nnodes):
            edata = data[n:n+n2]
            n += n2
            out = s2.unpack(edata)
            # (grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = out
            obj.add_sort1(dt, eid, *out)
    return n

def oes_crod_random_3(self, data: bytes, ndata: int,
             obj: Union[RandomRodStressArray, RandomRodStrainArray],
             nelements: int, ntotal: int) -> int:
    n = 0
    if self.is_debug_file:
        self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
        self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
        self.binary_debug.write('  #elementi = [eid_device, axial, torsion]\n')
        self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

    struct1 = Struct(self._endian + self._analysis_code_fmt + b'2f')
    add_sort_x = getattr(obj, 'add_sort' + str(self.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)
        (eid_device, axial, torsion) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)
        if self.is_debug_file:
            self.binary_debug.write('  eid=%i; C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))
        add_sort_x(dt, eid, axial, torsion)
        n += ntotal
    return n

def oes_cbar100_real_10(self, data: bytes, obj, nelements: int, ntotal: int, dt: Any) -> int:
    n = 0
    struct1 = Struct(self._endian + mapfmt(self._analysis_code_fmt + b'9f', self.size))
    for i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)
        (eid_device, sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)
        if self.is_debug_file:
            self.binary_debug.write('  eid=%i; C%i=[%s]\n' % (
                eid, i, ', '.join(['%r' % di for di in out])))
        n += ntotal
        obj.add_new_eid_sort1(self.element_name, dt, eid,
                              sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS)
    return n

def oes_cbeam_complex_111(self, data: bytes,
                          obj: Union[ComplexBeamStressArray, ComplexBeamStrainArray],
                          nelements: int, nnodes: int,
                          is_magnitude_phase: bool) -> int:
    n = 0
    #itotal = obj.itotal
    n1 = 44 * self.factor
    n2 = 40 * self.factor

    s1 = Struct(self._endian + mapfmt(self._analysis_code_fmt + b'i9f', self.size))
    if self.size == 4:
        s2 = Struct(self._endian + b'i9f')
    else:
        s2 = Struct(self._endian + b'q9d')

    add_sort_x = getattr(obj, 'add_sort' + str(self.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+n1]
        n += n1
        out1 = s1.unpack(edata)
        eid_device = out1[0]
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)
        if self.is_debug_file:
            self.binary_debug.write('CBEAM-2 - eid=%i out1=%s\n' % (eid, str(out1)))

        (grid, sd,
         excr, exdr, exer, exfr,
         exci, exdi, exei, exfi) = out1[1:]

        if is_magnitude_phase:
            exc = polar_to_real_imag(excr, exci)
            exd = polar_to_real_imag(exdr, exdi)
            exe = polar_to_real_imag(exer, exei)
            exf = polar_to_real_imag(exfr, exfi)
        else:
            exc = complex(excr, exci)
            exd = complex(exdr, exdi)
            exe = complex(exer, exei)
            exf = complex(exfr, exfi)

        add_sort_x(dt, eid, grid, sd,
                   exc, exd, exe, exf)

        for unused_inode in range(nnodes):
            edata = data[n:n+n2]
            n += n2
            out2 = s2.unpack(edata)
            (grid, sd,
             excr, exdr, exer, exfr,
             exci, exdi, exei, exfi) = out2

            if is_magnitude_phase:
                exc = polar_to_real_imag(excr, exci)
                exd = polar_to_real_imag(exdr, exdi)
                exe = polar_to_real_imag(exer, exei)
                exf = polar_to_real_imag(exfr, exfi)
            else:
                exc = complex(excr, exci)
                exd = complex(exdr, exdi)
                exe = complex(exer, exei)
                exf = complex(exfr, exfi)

            add_sort_x(dt, eid, grid, sd,
                       exc, exd, exe, exf)
            if self.is_debug_file:
                self.binary_debug.write('CBEAM-2 - eid=%i out2=%s\n' % (eid, str(out2)))
    return n

def oes_cbeam_random_67(self, data: bytes,
                        obj: Union[RandomBeamStressArray, RandomBeamStrainArray],
                        nelements: int, nnodes: int, dt: Any) -> int:
    n = 0
    n1 = 28 * self.factor
    n2 = 24 * self.factor # 6*4
    s1 = Struct(self._endian + self._analysis_code_fmt + b'i5f')
    s2 = Struct(self._endian + b'i5f')
    add_sort_x = getattr(obj, 'add_sort' + str(self.sort_method))
    add_eid_sort_x = getattr(obj, 'add_eid_sort' + str(self.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+n1]
        n += n1

        out = s1.unpack(edata)
        eid_device = out[0]

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)
        if self.is_debug_file:
            self.binary_debug.write('CBEAM-2 - eid=%i out=%s\n' % (eid, str(out)))

        #(grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = out
        add_eid_sort_x(dt, eid, *out[1:])

        for unused_inode in range(nnodes):
            edata = data[n:n+n2]
            n += n2
            out = s2.unpack(edata)
            # (grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = out
            add_sort_x(dt, eid, *out)
    return n

def oes_cquad4_33_complex_15(self,
                             data: bytes,
                             obj: Union[ComplexPlateStressArray, ComplexPlateStrainArray],
                             nelements: int, ntotal: int,
                             is_magnitude_phase: bool) -> int:
    """
    NX/MSC 2018.2
    2 FD1    RS Z1 = Fibre distance
    3 EX1R   RS Normal in x at Z1
    4 EX1I   RS Normal in x at Z1
    5 EY1R   RS Normal in y at Z1
    6 EY1I   RS Normal in y at Z1
    7 EXY1R  RS Shear in xy at Z1
    8 EXY1I  RS Shear in xy at Z1
    9 FD2    RS Z2 = Fibre distance
    10 EX2R  RS Normal in x at Z2
    11 EX2I  RS Normal in x at Z2
    12 EY2R  RS Normal in y at Z2
    13 EY2I  RS Normal in y at Z2
    14 EXY2R RS Shear in xy at Z2
    15 EXY2I RS Shear in xy at Z2
    """
    n = 0
    s1 = Struct(self._endian + self._analysis_code_fmt + b'14f')
    cen = 0 # 'CEN/4'
    add_sort_x = getattr(obj, 'add_sort' + str(self.sort_method))

    for unused_i in range(nelements):
        edata = data[n:n+ntotal]  # 4*15=60
        n += ntotal
        out = s1.unpack(edata)  # 15
        (eid_device,
         fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
         fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)
        if self.is_debug_file:
            self.binary_debug.write('  eid=%i C=%s\n' % (eid, str(out)))

        if is_magnitude_phase:
            sx1 = polar_to_real_imag(sx1r, sx1i)
            sx2 = polar_to_real_imag(sx2r, sx2i)
            sy1 = polar_to_real_imag(sy1r, sy1i)
            sy2 = polar_to_real_imag(sy2r, sy2i)
            txy1 = polar_to_real_imag(txy1r, txy1i)
            txy2 = polar_to_real_imag(txy2r, txy2i)
        else:
            sx1 = complex(sx1r, sx1i)
            sx2 = complex(sx2r, sx2i)
            sy1 = complex(sy1r, sy1i)
            sy2 = complex(sy2r, sy2i)
            txy1 = complex(txy1r, txy1i)
            txy2 = complex(txy2r, txy2i)

        add_sort_x(dt, eid, cen,
                   fd1, sx1, sy1, txy1,
                   fd2, sx2, sy2, txy2)
    return n

def oes_cquad4_random_vm_57(self,
                            data: bytes,
                            obj: Union[ComplexPlateVMStressArray, ComplexPlateVMStrainArray],
                            nelements: int, ntotal: int, nnodes: int,
                            dt: Any) -> int:
    """
    name    ntotal numwide
    CQUAD8  228    57
    CTRIA6  184
    """

    #ntotal = 228 # 57*4
    n = 0
    cen = 0
    if self.size == 4:
        struct1 = Struct(self._endian + b'i4si 10f')
        struct2 = Struct(self._endian + b'i 10f')
    else:
        print('ntotal =', ntotal, 228*2)
        struct1 = Struct(self._endian + b'i4si 10f')
        struct2 = Struct(self._endian + b'i 10f')
        raise NotImplementedError('64-bit')

    nnodes_corners = nnodes - 1
    #print('random-57: nelements =', nelements)
    factor = self.factor

    ntotal1 = 52 * factor
    ntotal2 = 44 * factor
    for unused_i in range(nelements):
        ni = 0
        edata = data[n:n+ntotal]
        #self.show_data(edata)
        #out = struct1.unpack(edata)  # len=17*4
        (eid_device, grid1, four,
         fd1, ox1, oy1, txy1, vm1,
         fd2, ox2, oy2, txy2, vm2) = struct1.unpack(edata[ni:ni+ntotal1])  # 13*4=52
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)
        #print(eid, cen)
        obj.add_sort1(dt, eid, cen,
                      fd1, ox2, oy2, txy2, vm1,
                      fd2, ox2, oy2, txy2, vm2)
        ni += ntotal1
        #print(eid, grid1,
              #fd1, ox1, oy1, txy1, vm1,
              #fd2, ox2, oy2, txy2, vm2)

        for unused_jnode in range(nnodes_corners): # 4*44=176
            (grid,
             fd1, ox1, oy1, txy1, vm1,
             fd2, ox2, oy2, txy2, vm2) = struct2.unpack(edata[ni:ni+ntotal2])  # 11*4=44
            #print(eid, grid)
            ni += ntotal2
            #print(grid,
                  #fd1, ox1, oy1, txy1, vm1,
                  #fd2, ox2, oy2, txy2, vm2)
            obj.add_sort1(dt, eid, grid,
                          fd1, ox2, oy2, txy2, vm1,
                          fd2, ox2, oy2, txy2, vm2)
        assert ni == ntotal, f'ni={ni}'
        n += ntotal
        #print()
    #assert n == ndata
    #n = self._not_implemented_or_skip(data, ndata, msg)
    #nelements = None
    #ntotal = None
    return n

def oes_cquad4_144_complex_77(self,
                              data: bytes,
                              obj: Union[ComplexPlateStressArray, ComplexPlateStrainArray],
                              nelements: int, nnodes: int,
                              dt: Any,
                              is_magnitude_phase: bool) -> int:
    #print('nelements =', nelements)
    n = 0
    grid_center = 0
    if self.size == 4:
        s1 = self.struct_2i  # 2
    else:
        s1 = self.struct_2q  # 2
    s2 = Struct(self._endian + mapfmt(b'i14f', self.size)) # 15

    ntotal1 = 8 * self.factor
    ntotal2 = 60 * self.factor  # 4*15

    add_sort_x = getattr(obj, 'add_sort' + str(self.sort_method))
    for unused_i in range(nelements):
        (eid_device, _) = s1.unpack(data[n:n+ntotal1])
        n += ntotal1

        edata = data[n:n+ntotal2]
        n += ntotal2
        out = s2.unpack(edata)  # len=15*4

        eid = self._check_id(eid_device, 'element', out)
        if self.is_debug_file:
            self.binary_debug.write('%s\n' % (str(out)))
        (grid,
         fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
         fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i) = out
        #grid_center = 'CEN/%i' % grid   # this is correct, but fails

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

        add_sort_x(dt, eid, grid_center,
                   fd1, sx1, sy1, txy1,
                   fd2, sx2, sy2, txy2)

        for unused_node_id in range(nnodes):  # nodes pts
            edata = data[n:n + ntotal2]
            n += ntotal2
            out = s2.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('%s\n' % (str(out)))
            (grid,
             fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
             fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i) = out

            if is_magnitude_phase:
                sx1 = polar_to_real_imag(sx1r, sx1i)
                sx2 = polar_to_real_imag(sx2r, sx2i)
                sy1 = polar_to_real_imag(sy1r, sy1i)
                sy2 = polar_to_real_imag(sy2r, sy2i)
                txy1 = polar_to_real_imag(txy1r, txy1i)
                txy2 = polar_to_real_imag(txy2r, txy2i)
            else:
                sx1 = complex(sx1r, sx1i)
                sx2 = complex(sx2r, sx2i)
                sy1 = complex(sy1r, sy1i)
                sy2 = complex(sy2r, sy2i)
                txy1 = complex(txy1r, txy1i)
                txy2 = complex(txy2r, txy2i)

            add_sort_x(dt, eid, grid,
                       fd1, sx1, sy1, txy1,
                       fd2, sx2, sy2, txy2)
    return n

def oes_cquad4_complex_vm_87(self, data: bytes,
                             obj: Union[ComplexPlateVMStressArray, ComplexPlateVMStrainArray],
                             nelements: int, nnodes: int,
                             is_magnitude_phase: bool) -> int:
    """
    etype name      numwide nodes
    33  CQUAD4-33 -> 17      (not in this function)
    82  CQUADR    -> 87      (1+4 nodes)
    144 CQUAD144  -> 87      (1+4 nodes)
    ??? CTRIA6    -> 70?     (1+3 nodes)

    2 TERM    CHAR4
    3 GRID    I
    4 FD1     RS Fiber distance at z1
    5 EX1R    RS Normal in x at z1
    6 EX1I    RS Normal in x at z1
    7 EY1R    RS Normal in y at z1
    8 EY1I    RS Normal in y at z1
    9 ETXY1R  RS Shear in xy at z1
    10 ETXY1I RS Shear in xy at z1
    11 VM1    RS von Mises at z1
    12 FD2    RS Fiber distance at z2
    13 EX2R   RS Normal in x at z2
    14 EX2I   RS Normal in x at z2
    15 EY2R   RS Normal in y at z2
    16 EY2I   RS Normal in y at z2
    17 ETXY2R RS Shear in xy at z2
    18 ETXY2I RS Shear in xy at z2
    19 VM2    RS von Mises at z2

    """
    #strings = (b'20.0 CEN/ 4 \xcd\xccL\xbd\x91@\x1aD\xebX\xac\xc1b\x80FB\x96\xc9\xdd\xbf\xe7\xbe\n,\xaah\xc5\xa9b\x87\x14D\xcd\xccL=\x91@\x1a\xc4\xebX\xacAb\x80F\xc2\x96\xc9\xdd?\xe7\xbe\n\xac\xaah\xc5)b\x87\x14D\x01\x00\x00\x00\xcd\xccL\xbd\xc7\xd5\x9fD\xc9\x952\xc2b\x80FC\x96\xc9\xdd\xc0\xe7\xbe\n,\xaah\xc5\xa9\xca\x13\x95D\xcd\xccL=\xc7\xd5\x9f\xc4\xc9\x952Bb\x80F\xc3\x96\xc9\xdd@\xe7\xbe\n\xac\xaah\xc5)\xca\x13\x95D\x02\x00\x00\x00\xcd\xccL\xbd\xbf\xa62\xc2\xd4\x9b\xc7?b\x80\xc6\xc2\x96\xc9]@\xe7\xbe\n,\xaah\xc5\xa9nL\xacB\xcd\xccL=\xbf\xa62B\xd4\x9b\xc7\xbfb\x80\xc6B\x96\xc9]\xc0\xe7\xbe\n\xac\xaah\xc5)nL\xacB\x03\x00\x00\x00\xcd\xccL\xbd\xbf\xa62\xc2\xd4\x9b\xc7?b\x80\xc6\xc2\x96\xc9]@\xe7\xbe\n,\xaah\xc5\xa9nL\xacB\xcd\xccL=\xbf\xa62B\xd4\x9b\xc7\xbfb\x80\xc6B\x96\xc9]\xc0\xe7\xbe\n\xac\xaah\xc5)nL\xacB\x04\x00\x00\x00\xcd\xccL\xbd\xc7\xd5\x9fD\xc9\x952\xc2b\x80FC\x96\xc9\xdd\xc0\xe7\xbe\n,\xaah\xc5\xa9\xca\x13\x95D\xcd\xccL=\xc7\xd5\x9f\xc4\xc9\x952Bb\x80F\xc3\x96\xc9\xdd@\xe7\xbe\n\xac\xaah\xc5)\xca\x13\x95D',)
    # ints    = (20.0, CEN/, 4,
    #           -0.05, 617.0, -21.5, 49.6, -1.7, 738901735,              -1446680406,            1142196066,
    #           0.05, -617.0,  21.5, -49.6, 1.7, -1408581913, 700803242, 1142196066, 1, -1119040307, 1151325639, -1036872247, 1128693858, -1059206762, 738901735, -1446680406, 1150620618, 1028443341, -996158009, 1110611401, -1018789790, 1088276886, -1408581913, 700803242, 1150620618, 2, -1119040307, -1036867905, 1070046164, -1027178398, 1079888278, 738901735, -1446680406, 1118588014, 1028443341, 1110615743, -1077437484, 1120305250, -1067595370, -1408581913, 700803242, 1118588014, 3, -1119040307, -1036867905, 1070046164, -1027178398, 1079888278, 738901735, -1446680406, 1118588014, 1028443341, 1110615743, -1077437484, 1120305250, -1067595370, -1408581913, 700803242, 1118588014, 4, -1119040307, 1151325639, -1036872247, 1128693858, -1059206762, 738901735, -1446680406, 1150620618, 1028443341, -996158009, 1110611401, -1018789790, 1088276886, -1408581913, 700803242, 1150620618)
    #floats  = (20.0, CEN/, 4,
    #           -0.05,  617.0,   -21.5,   49.6,  -1.7,  1.9716951595721843e-12, -8.766713754677219e-14,  594.1,
    #            0.05, -617.0,    21.5,  -49.6,   1.7, -1.9716951595721843e-12,  8.766713754677219e-14,  594.1,
    #        1.4012984e-45,
    #           -0.05,  1278.68, -44.6,  198.5,  -6.9,  1.9716951595721843e-12, -8.766713754677219e-14, 1192.6,
    #            0.05, -1278.68,  44.6, -198.5,   6.9, -1.9716951595721843e-12,  8.766713754677219e-14, 1192.6,
    #        2.802596928649634e-45,
    #           -0.05, -44.6,      1.5,  -99.2,   3.4,  1.9716951595721843e-12, -8.766713754677219e-14, 86.149,
    #            0.05,  44.6,     -1.5,   99.2,  -3.4, -1.9716951595721843e-12,  8.766713754677219e-14, 86.149,
    #        4.203895392974451e-45,
    #           -0.05, -44.6,      1.5,  -99.2,   3.4,  1.9716951595721843e-12, -8.766713754677219e-14, 86.149,
    #            0.05,  44.6,     -1.5,   99.2,  -3.4, -1.9716951595721843e-12,  8.766713754677219e-14, 86.149,
    #        5.605193857299268e-45,
    #           -0.05, 1278.6,    -44.6, 198.5,  -6.9,  1.9716951595721843e-12, -8.766713754677219e-14, 1192.6,
    #           0.05, -1278.6,     44.6,-198.5,   6.9, -1.9716951595721843e-12,  8.766713754677219e-14, 1192.6)
    # -------------------------
    #ints    = (0,    'CEN/',
    #           5,     -0.025, 1131762990, -1056779073, 1137493831, -1051039598, 1133740074, -1055043205, 1142680183,
    #                   0.025, -1016740689, 1089802185, -1009242998, 1097240578, -1013972178, 1092196624, 1142689762,
    #           80000, -0.025, -1039337422, 1066818106, -1022560206, 1083595322, 1118381010, -1070307877, 1128400256,
    #                   0.025, 1111837962, -1076728084, 1128615178, -1059950868, -1029871408, 1076355832, 1130378350,
    #           80002, -0.025, 1144778127, -1043828805, 1161555343, -1027051589, 1134339017, -1054404396, 1160567412,
    #                   0.025, -1002241635, 1104149606, -985464419, 1120926822, -1014451576, 1091685317, 1160927605,
    #           80202, -0.025, -1013123432, 1093101862, -996473336, 1109743497, 1143660503, -1045020817, 1154190671,
    #                   0.025, 1133381240, -1055425923, 1150079339, -1038733090, -1002862704, 1103487198, 1154274046,
    #           80200, -0.025, 1141944942, -1046850567, 1095277553, -1093274063, 1120227570, -1068338410, 1142258188,
    #                   0.025, -1007929795, 1098641187, 1116060214, -1072783145, -1033402643, 1073147904, 1140969037,
    #           2.5,  'CEN/',
    #           5,     -0.025, 1131832225, -1056704978, 1137555116, -1050908425, 1133787322, -1054942075, 1142729036,
    #                   0.025, -1016671483, 1089950313, -1009181671, 1097371840, -1013924936, 1092297741, 1142738537,
    #           80000, -0.025, -1039289241, 1066921230, -1022512025, 1083698446, 1118448079, -1070164312, 1128472712,
    #                   0.025, 1111886248, -1076624736, 1128663464, -1059847520, -1029804359, 1076499355, 1130448349,
    #           80002, -0.025, 1144846119, -1043683269, 1161623335, -1026906053, 1134386659, -1054302422, 1160629636,
    #                   0.025, -1002173624, 1104295182, -985396408, 1121072398, -1014403979, 1091787195, 1160989873,
    #           80202, -0.025, -1013031023, 1093299683, -996403997, 1109891922, 1143714425, -1044905405, 1154266576,
    #                   0.025, 1133473626, -1055228152, 1150148641, -1038584742, -1002808755, 1103602670, 1154348580,
    #           80200, -0.025, 1141995401, -1046742558, 1095362185, -1093092785, 1120294524, -1068195089, 1142310151,
    #                   0.025, -1007828983, 1098856980, 1116049937, -1072805159, -1033268970, 1073434045, 1141020308,
    #           10.0, 'CEN/',
    #           5,     -0.025, 1132706983, -1055459734, 1138546383, -1048703317, 1134552107, -1053240134, 1143519784,
    #                   0.025, -1015552254, 1091479459, -1008189778, 1099242966, -1013160238, 1093999495, 1143528113,
    #           80000, -0.025, -1038510246, 1068653763, -1021733030, 1085430979, 1119535604, -1067741982, 1129646715,
    #                   0.025, 1112666820, -1074888834, 1129444036, -1058111618, -1028717126, 1078921059, 1131585244,
    #           80002, -0.025,  1145947464, -1041231474, 1162724680, -1024454258, 1135157762, -1052586493, 1161637710,
    #                   0.025, -1001071999, 1106747574, -984294783, 1123524790, -1013633545, 1093501695, 1161998610,
    #           80202, -0.025, -1011529318,  1096648297, -995279005, 1112398467, 1144586718, -1042964811, 1155498288,
    #                   0.025, 1134974980, -1051880286, 1151273087, -1036079365, -1001936046, 1105544152, 1155559868,
    #           80200, -0.025,  1142813256, -1044921274, 1096762940, -1090229644, 1121380799, -1065774838, 1143152359,
    #                   0.025, -1006413905, 1100701913, 1115879378, -1073189541, -1031449359, 1076004431, 1141852280)
    #floats  = (0.0,  'CEN/',
    #           5,       -0.025,  245.33273315429688, -8.176939964294434,   409.5568542480469, -13.650529861450195,  295.00128173828125,  -9.832392692565918,
    #           623.8668, 0.025, -229.76829528808594,  7.658176898956299, -432.34796142578125,  14.410158157348633, -288.02484130859375,   9.599868774414062, 624.2481689453125,
    #           1.12-40, -0.025, -35.24237823486328,   1.174628496170044, -140.96951293945312,   4.698513984680176,   84.56996154785156, -2.8187167644500732,
    #           194.4375, 0.025,  49.325233459472656, -1.644010066986084,  197.30093383789062,  -6.576040267944336,   -78.7047119140625,   2.623228073120117, 224.20480346679688,
    #           1.12106, -0.025,  751.7118530273438, -25.054555892944336,   3006.847412109375, -100.21822357177734,   313.2795715332031, -10.441608428955078,
    #           2765.65,  0.025, -780.0252075195312, 25.998241424560547, -3120.100830078125, 103.99296569824219, -273.394775390625, 9.112248420715332, 2853.591064453125,
    #           1.1238,  -0.025, -313.926513671875, 10.463171005249023, -1240.1884765625, 41.33548355102539, 683.4974975585938, -22.780973434448242,
    #           1628.4,   0.025, 284.050537109375, -9.46740436553955, 1126.5443115234375, -37.54772186279297, -742.1181640625, 24.734798431396484, 1638.593505859375,
    #           1.12,    -0.025, 578.7879638671875, -19.29100227355957, 12.538071632385254, -0.4178939163684845, 98.65809631347656, -3.2882742881774902,
    #           597.90,   0.025, -472.4237365722656, 15.745882987976074, 66.86369323730469, -2.228566884994507, -57.88176345825195, 1.92919921875, 519.2234497070312,
    #           2.5,  'CEN/',
    #           5,          -0.025, 246.38917541503906, -8.247602462768555, 411.4271240234375, -13.775626182556152, 296.44317626953125, -9.928837776184082,
    #           626.671875,  0.025, -230.8242950439453, 7.728809833526611, -434.2195129394531, 14.53533935546875, -289.466552734375, 9.696301460266113, 627.2251586914062,
    #           1.1217e-40, -0.025, -35.42617416381836, 1.1869218349456787, -141.70469665527344, 4.747687339782715, 85.08165740966797, -2.852945327758789,
    #           195.153125,  0.025, 49.509429931640625, -1.6563301086425781, 198.0377197265625, -6.6253204345703125, -79.21625518798828, 2.6574466228485107, 225.2729034423828,
    #           1.1212e-40, -0.025, 755.8617553710938, -25.332143783569336, 3023.447021484375, -101.32857513427734, 314.7334899902344, -10.538858413696289,
    #           2780.65625,  0.025, -784.17626953125, 26.27590560913086, -3136.705078125, 105.10362243652344, -274.8473205566406, 9.209406852722168, 2868.793212890625,
    #           1.1358e-40, -0.025, -316.7466125488281, 10.651827812194824, -1248.6527099609375, 41.90167999267578, 686.7886352539062, -23.0011043548584,
    #           1637.40625,  0.025, 286.86993408203125, -9.656013488769531, 1135.0040283203125, -38.113624572753906, -745.4109497070312, 24.95504379272461, 1647.69189453125,
    #           1.1383e-40, -0.025, 581.8677368164062, -19.497013092041016, 12.618782997131348, -0.42329642176628113, 99.16891479492188, -3.3224446773529053,
    #           601.471875,  0.025, -475.5002746582031, 15.951679229736328, 66.78528594970703, -2.223318338394165, -58.391685485839844, 1.96330988407135, 522.352783203125,
    #           10.0, 'CEN/',
    #           5,          -0.025, 263.4738464355469, -9.435159683227539, 441.6781921386719, -15.878581047058105, 319.7825622558594, -11.551935195922852,
    #           674.140625,  0.025, -247.90237426757812, 8.915926933288574, -464.48968505859375, 16.639568328857422, -312.80328369140625, 11.319220542907715, 675.4170532226562,
    #           1.1337e-40, -0.025, -38.397804260253906, 1.3934558629989624, -153.9121704101562, 5.57382345199585, 93.37881469726562, -3.430473804473877,
    #           213.310547,  0.025, 52.48707580566406, -1.8632657527923584, 209.94830322265625, -7.453063011169434, -87.51118469238281, 3.234825849533081, 242.62054443359375,
    #           1.1202e-40, -0.025, 823.08251953125, -30.00857162475586, 3292.330078125, -120.03428649902344, 338.26568603515625, -12.17529582977295,
    #           3026296875,  0.025, -851.4141235351562, 30.953472137451172, -3405.656494140625, 123.81388854980469, -298.3591003417969, 10.844481468200684, 3115.06689453125,
    #           1.1298e-40, -0.025, -362.57501220703125, 13.845314979553223, -1385.9808349609375, 51.4633903503418, 740.0291748046875, -26.70249366760254,
    #           1788109375,  0.025, 332.6876220703125, -12.848787307739258, 1272.2655029296875, -47.67087936401367, -798.6768798828125, 28.658126831054688, 1795.55419921875,
    #           1.1293e-40, -0.025, 631.78564453125, -22.970836639404297, 13.954647064208984, -0.517249345779419, 107.45653533935547, -3.899477481842041,
    #           652.410938,  0.025, -525.3700561523438, 19.42228889465332, 65.48402404785156, -2.1316745281219482, -66.66590118408203, 2.539447546005249, 573.13232421875)
    n = 0
    ntotal1 = 8 * self.factor  # 2*4
    ntotal2 = 68 * self.factor # 17*4
    add_sort_x = getattr(obj, 'add_sort' + str(self.sort_method))
    #print('ndata', len(data))
    #print('ntotal', ntotal1+nnodes*ntotal2)
    #print('nelements =', nelements)
    for i in range(nelements):
        #print(f'{i}/{nelements}')
        edata1 = data[n:n+ntotal1]
        #self.show_data(edata1)
        n += ntotal1

        fmt1 = mapfmt(self._endian + self._analysis_code_fmt + b'4s', self.size)
        #fmt2 = mapfmt(self._endian + b'i 5f2if 5f2if', self.size)
        fmt2 = mapfmt(self._endian + b'i 16f', self.size)
        s1 = Struct(fmt1)
        s2 = Struct(fmt2)

        eid_device, cen = s1.unpack(edata1)
        #self.show_data(data, types='ifs')
        assert cen == b'CEN/', (eid_device, cen)

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)

        for i in range(nnodes):
            edata2 = data[n:n+ntotal2]

            #DISTANCE             NORMAL-X                     NORMAL-Y                    SHEAR-XY             VON MISES
            out = s2.unpack(edata2)
            (nid,
             fd1, oxx1r, oxx1i, oyy1r, oyy1i, txy1r, txy1i, ovm1,
             fd2, oxx2r, oxx2i, oyy2r, oyy2i, txy2r, txy2i, ovm2) = out
            fd1, oxx1r, oxx1i, oyy1r, oyy1i, txy1r, txy1i, ovm1,
            if i == 0:
                nid = 0

            if is_magnitude_phase:
                oxx1 = polar_to_real_imag(oxx1r, oxx1i)
                oxx2 = polar_to_real_imag(oxx2r, oxx2i)
                oyy1 = polar_to_real_imag(oyy1r, oyy1i)
                oyy2 = polar_to_real_imag(oyy2r, oyy2i)
                txy1 = polar_to_real_imag(txy1r, txy1i)
                txy2 = polar_to_real_imag(txy2r, txy2i)
            else:
                oxx1 = complex(oxx1r, oxx1i)
                oxx2 = complex(oxx2r, oxx2i)
                oyy1 = complex(oyy1r, oyy1i)
                oyy2 = complex(oyy2r, oyy2i)
                txy1 = complex(txy1r, txy1i)
                txy2 = complex(txy2r, txy2i)

            #print((eid, dt), nid, f'{fd1:g}, {fd2:g}')
            #print('real', f'{oxx1r:g}, {oxx1i:g}, {oyy1r:g}, {oyy1i:g}, {txy1r:g}, {txy1i:g}, {ovm1:g}')
            #print('imag', oxx2i, oxx2i, oyy2r, oyy2i, txy2r, txy2i, ovm2)
            add_sort_x(dt, eid, nid,
                       fd1, oxx1, oyy1, txy1, ovm1,
                       fd2, oxx2, oyy2, txy2, ovm2)
            #print(out)
            #assert nid in (1, 2, 3, 4, 17), ((eid, dt), nid, fd1, fd2)
            n += ntotal2
        #print('-----------')
    return n

def oes_cquad4_33_random_9(self, data: bytes,
                           obj: Union[RandomPlateStressArray, RandomPlateStrainArray],
                           nelements: int, ntotal: int) -> int:
    n = 0
    add_sort_x = getattr(obj, 'add_sort' + str(self.sort_method))
    #print('cquad33_9 - SORT1')
    struct1 = Struct(self._endian + self._analysis_code_fmt + b'8f')
    #cen = 0 # CEN/4
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device,
         fd1, sx1, sy1, txy1,
         fd2, sx2, sy2, txy2,) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)
        #self.log.info(f'SORT1 no VM: eid={eid} dt={dt}')

        if self.is_debug_file:
            self.binary_debug.write('  eid=%i C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))
        add_sort_x(dt, eid, 0,
                   fd1, sx1, sy1, txy1,
                   fd2, sx2, sy2, txy2)
        n += ntotal
    return n

def oes_ctria3_real_17(self, data: bytes,
                       obj: Union[RealPlateStressArray, RealPlateStrainArray],
                       ntotal: int, nelements: int, dt: Any) -> int:
    n = 0
    cen = 0 # 'CEN/3'
    #assert self.sort_method == 1, self.code_information()

    add_new_eid_sort_x = getattr(obj, 'add_new_eid_sort' + str(self.sort_method))
    #print(add_new_eid_sort_x)
    struct1 = Struct(self._endian + mapfmt(self._analysis_code_fmt + b'16f', self.size))
    #print('nelements =', nelements, nelements*2, obj.ntotal)
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)
        (eid_device,
         fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
         fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)
        if self.is_debug_file:
            self.binary_debug.write('  OES CTRIA3-74 - eid=%i; C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))
        #print('  OES CTRIA3-74 - eid=%i; C=[%s]\n' % (
            #eid, ', '.join(['%r' % di for di in out])))

        add_new_eid_sort_x(dt, eid, cen,
                           fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                           fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2)
        n += ntotal
    return n

def oes_quad4_33_real_17(self, data: bytes,
                         obj: Union[RealPlateStressArray, RealPlateStrainArray],
                         ntotal: int, nelements: int, dt: Any) -> int:
    n = 0
    cen = 0 # CEN/4
    #assert self.sort_method == 1, self.code_information()
    add_new_eid_sort_x = getattr(obj, 'add_new_eid_sort' + str(self.sort_method))
    fmt = mapfmt(self._endian + self._analysis_code_fmt + b'16f', self.size)
    struct1 = Struct(fmt)
    #if self.sort_method == 1:
        #print('**nelements =', nelements, nelements*2, obj.ntotal)
    #else:
        #print('**ntimes =', nelements)
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device,
         fd1, sx1, sy1, txy1, angle1, major1, minor1, max_shear1,
         fd2, sx2, sy2, txy2, angle2, major2, minor2, max_shear2) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)
        if self.is_debug_file:
            self.binary_debug.write('  eid=%i C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))

        add_new_eid_sort_x(
            dt, eid, cen,
            fd1, sx1, sy1, txy1, angle1, major1, minor1, max_shear1,
            fd2, sx2, sy2, txy2, angle2, major2, minor2, max_shear2)
        n += ntotal
    return n

def oes_ctria3_complex_vm_17(self,
                             data: bytes,
                             obj: Union[ComplexPlateVMStressArray, ComplexPlateVMStrainArray],
                             nelements: int, ntotal: int,
                             dt,
                             is_magnitude_phase: bool) -> int:
    """
    #                   C O M P L E X   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )
    #                                                          (REAL/IMAGINARY)
    #
    #  ELEMENT       FIBER                                     - STRESSES IN ELEMENT  COORDINATE SYSTEM -
    #    ID.        DISTANCE              NORMAL-X                       NORMAL-Y                      SHEAR-XY               VON MISES
    #0       1  -4.359080E+00  -1.391918E+00 /  2.474756E-03  -1.423926E+00 /  2.530494E-03   2.655153E-02 / -5.158625E-05   1.408948E+00
    #            4.359080E+00   1.391918E+00 / -2.474756E-03   1.423926E+00 / -2.530494E-03  -2.655153E-02 /  5.158625E-05   1.408948E+00

    """
    n = 0
    cen = 0
    struct1 = Struct(mapfmt(self._endian + self._analysis_code_fmt + b'16f', self.size))
    add_sort_x = getattr(obj, 'add_sort' + str(self.sort_method))

    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)
        (eid_device, fd1, oxx1r, oxx1i, oyy1r, oyy1i, txy1r, txy1i, ovm1,
                     fd2, oxx2r, oxx2i, oyy2r, oyy2i, txy2r, txy2i, ovm2, ) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)

        if is_magnitude_phase:
            oxx1 = polar_to_real_imag(oxx1r, oxx1i)
            oxx2 = polar_to_real_imag(oxx2r, oxx2i)
            oyy1 = polar_to_real_imag(oyy1r, oyy1i)
            oyy2 = polar_to_real_imag(oyy2r, oyy2i)
            txy1 = polar_to_real_imag(txy1r, txy1i)
            txy2 = polar_to_real_imag(txy2r, txy2i)
        else:
            oxx1 = complex(oxx1r, oxx1i)
            oxx2 = complex(oxx2r, oxx2i)
            oyy1 = complex(oyy1r, oyy1i)
            oyy2 = complex(oyy2r, oyy2i)
            txy1 = complex(txy1r, txy1i)
            txy2 = complex(txy2r, txy2i)
        #print(dt, eid, cen, sx1, sy1, txy1, max_shear1)
        #print(dt, eid, cen, sx2, sy2, txy2, max_shear2)
        add_sort_x(dt, eid, cen,
                   fd1, oxx1, oyy1, txy1, ovm1,
                   fd2, oxx2, oyy2, txy2, ovm2)
        n += ntotal

    #msg = self.code_information()
    #self.log.warning(f'  skipping {self.table_name_str} {self.element_name}')
    #msg = '%s-%s' % (self.table_name_str, self.element_name)
    return n

def oes_ctria3_random_9(self, data: bytes,
                        obj: Union[RandomPlateStressArray, RandomPlateStrainArray],
                        nelements: int, ntotal: int) -> int:
    n = 0
    struct1 = Struct(self._endian + self._analysis_code_fmt + b'8f')
    cen = 0 # CEN/4
    add_sort_x = getattr(obj, 'add_sort' + str(self.sort_method))
    #print('ctria3_9 - SORT1')
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device,
         fd1, sx1, sy1, txy1,
         fd2, sx2, sy2, txy2,) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)

        if self.is_debug_file:
            self.binary_debug.write('  eid=%i C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))

        add_sort_x(dt, eid, cen,
                   fd1, sx1, sy1, txy1,
                   fd2, sx2, sy2, txy2)
        n += ntotal
    return n

def oes_cquad4_33_random_vm_11(self, data: bytes,
                               obj: Union[RandomPlateVMStressArray, RandomPlateVMStrainArray],
                               nelements: int, ntotal: int) -> int:
    struct1 = Struct(self._endian + self._analysis_code_fmt + b'10f')
    cen = 0 # CEN/4
    n = 0
    add_sort_x = getattr(obj, 'add_sort' + str(self.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device,
         fd1, sx1, sy1, txy1, ovm1,
         fd2, sx2, sy2, txy2, ovm2,) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)
        #self.log.info(f'SORT1 VM: eid={eid} dt={dt}')

        if self.is_debug_file:
            self.binary_debug.write('  eid=%i C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))

        add_sort_x(dt, eid, cen,
                   fd1, sx1, sy1, txy1, ovm1,
                   fd2, sx2, sy2, txy2, ovm2)
        n += ntotal
    return n

def oes_ctria3_random_vm_11(self, data: bytes,
                            obj: Union[RandomPlateStressArray, RandomPlateStrainArray],
                            nelements: int, ntotal: int) -> int:
    n = 0
    struct1 = Struct(self._endian + self._analysis_code_fmt + b'10f')

    add_sort_x = getattr(obj, 'add_sort' + str(self.sort_method))
    cen = 0 # CEN/4
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device,
         fd1, sx1, sy1, txy1, ovm1,
         fd2, sx2, sy2, txy2, ovm2,) = out
        #print('CTRIA3', self.nonlinear_factor, out)

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)
        #self.log.info(f'SORT1 VM: eid={eid} dt={dt}')

        if self.is_debug_file:
            self.binary_debug.write('  eid=%i C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))

        add_sort_x(dt, eid, cen,
                  fd1, sx1, sy1, txy1, ovm1,
                  fd2, sx2, sy2, txy2, ovm2)
        n += ntotal
    return n

def oes_cquad4_144_real(self, data: bytes, ndata: int,
                        obj: RealPlateStrainArray,
                        nelements: int, nnodes: int, dt: Any) -> int:
    n = 0
    if self.size == 4:
        center_format = self._endian + self._analysis_code_fmt + b'4si16f'
        node_format = self._endian + b'i16f'
    else:
        center_format = self._endian + mapfmt(self._analysis_code_fmt, self.size) + b'8sq16d'
        node_format = self._endian + b'q16d'
    cs = Struct(center_format)
    ns = Struct(node_format)

    add_new_eid_sort_x = getattr(obj, 'add_new_eid_sort' + str(self.sort_method))
    add_sort_x = getattr(obj, 'add_sort' + str(self.sort_method))
    if self.is_debug_file:
        self.binary_debug.write(
            '  [cap, element1, element2, ..., cap]\n'
            '  cap = %i  # assume 1 cap when there could have been multiple\n'
            '  #elementi = [centeri, node1i, node2i, node3i, node4i]\n'
            '  #centeri = [eid_device, j, grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n'
            '  #                                fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,)]\n'
            '  #nodeji = [grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n'
            '  #                fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,)]\n'
            '  nelements=%i; nnodes=%i # +1 centroid\n' % (ndata, nelements, nnodes))

    grid_center = 0
    n76 = 76 * self.factor
    n68 = 68 * self.factor
    for unused_i in range(nelements):
        edata = data[n:n+n76]

        out = cs.unpack(edata)  # len=17*4
        # j is the number of nodes, so CQUAD4 -> 4, but we don't need to save it...
        (eid_device, unused_j,
         grid,
         fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
         fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)
        #print(f'tri eid={eid} dt={dt} nid={grid_center}')
        #print(out[:3])
        if self.is_debug_file:
            self.binary_debug.write('  eid=%i; C=[%s]\n' % (eid, ', '.join(['%r' % di for di in out])))
        add_new_eid_sort_x(dt, eid, grid_center,
                           fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                           fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2)
        n += n76
        for inode in range(nnodes):
            out = ns.unpack(data[n:n + n68])
            (grid,
             fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
             fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out

            if self.is_debug_file:
                d = tuple([grid,
                           fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                           fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2])
                self.binary_debug.write('  node%i = [%s]\n' % (inode+1, ', '.join(['%r' % di for di in d])))
            assert isinstance(grid, int), grid
            assert grid > 0, grid

            #print(f' {inode}: eid={eid} dt={dt} nid={grid}')
            add_sort_x(dt, eid, grid,
                       fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                       fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2)
            n += n68
    return n

def oes_cquad4_144_random(self, data: bytes,
                          obj: Union[RandomPlateStressArray, RandomPlateStrainArray],
                          nelements: int, nnodes: int, ndata: int) -> int:
    n = 0
    center_format = self._endian + self._analysis_code_fmt + b'4s i8f'
    node_format = self._endian + b'i8f'
    cs = Struct(center_format)
    ns = Struct(node_format)

    if self.is_debug_file:
        self.binary_debug.write(
            '  [cap, element1, element2, ..., cap]\n'
            '  cap = %i  # assume 1 cap when there could have been multiple\n'
            '  #elementi = [centeri, node1i, node2i, node3i, node4i]\n'
            '  #centeri = [eid_device, j, grid, fd1, sx1, sy1, txy1,\n'
            '  #                                fd2, sx2, sy2, txy2,)]\n'
            '  #nodeji = [grid, fd1, sx1, sy1, txy1,\n'
            '  #                fd2, sx2, sy2, txy2,)]\n'
            '  nelements=%i; nnodes=%i # +1 centroid\n' % (ndata, nelements, nnodes))

    grid_center = 0
    add_sort_x = getattr(obj, 'add_sort' + str(self.sort_method))

    ntotal1 = 44 * self.factor
    ntotal2 = 36 * self.factor
    for unused_i in range(nelements):
        edata = data[n:n+ntotal1]
        #self.show_data(edata)
        out = cs.unpack(edata)  # len=17*4
        #print(out)
        # j is the number of nodes, so CQUAD4 -> 4, but we don't need to save it...
        eid_device = out[0]
        (eid_device, unused_j,
         grid,
         fd1, sx1, sy1, txy1,
         fd2, sx2, sy2, txy2,) = out
        #print(out)

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)

        if self.is_debug_file:
            self.binary_debug.write('  eid=%i; C=[%s]\n' % (eid, ', '.join(['%r' % di for di in out])))

        add_sort_x(dt, eid, grid_center,
                   fd1, sx1, sy1, txy1, fd2, sx2, sy2, txy2)
        n += ntotal1
        for inode in range(nnodes):
            out = ns.unpack(data[n:n + ntotal2])
            #print(out)
            (grid,
             fd1, sx1, sy1, txy1,
             fd2, sx2, sy2, txy2,) = out
            #print(out)
            if self.is_debug_file:
                d = tuple([grid,
                           fd1, sx1, sy1, txy1,
                           fd2, sx2, sy2, txy2])
                self.binary_debug.write('  node%i = [%s]\n' % (inode+1, ', '.join(['%r' % di for di in d])))
            assert isinstance(grid, int), grid
            assert grid > 0, grid

            # leaving off grid
            add_sort_x(dt, eid, grid,
                       fd1, sx1, sy1, txy1,
                       fd2, sx2, sy2, txy2)
            n += ntotal2
    return n

def oes_cbar_real_16(self, data: bytes,
                     obj: Union[RealBarStressArray, RealBarStrainArray],
                     nelements: int, ntotal: int, dt: Any) -> int:
    n = 0
    fmt = mapfmt(self._endian + self._analysis_code_fmt + b'15f', self.size)
    struct1 = Struct(fmt)
    add_sort_x = getattr(obj, 'add_new_eid_sort' + str(self.sort_method))
    for i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)
        (eid_device,
         s1a, s2a, s3a, s4a, axial, smaxa, smina, margin_tension,
         s1b, s2b, s3b, s4b, smaxb, sminb, margin_compression) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)
        if self.is_debug_file:
            self.binary_debug.write('  eid=%i; C%i=[%s]\n' % (
                eid, i, ', '.join(['%r' % di for di in out])))
        n += ntotal
        add_sort_x(
            dt, eid,
            s1a, s2a, s3a, s4a, axial, smaxa, smina, margin_tension,
            s1b, s2b, s3b, s4b, smaxb, sminb, margin_compression)

    if self.is_sort2:
        #print(add_sort_x)
        #print(''.join(obj.get_stats()))
        #print(f'{self.table_name} sort_method={self.sort_method}', obj._times)
        assert len(np.unique(obj._times)) == len(obj._times), obj._times.tolist()
    return n

def oes_cbar_complex_19(self,
                        data: bytes,
                        obj: Union[ComplexBarStressArray, ComplexBarStrainArray],
                        nelements: int, ntotal: int,
                        is_magnitude_phase: bool) -> int:
    n = 0
    struct1 = Struct(mapfmt(self._endian + self._analysis_code_fmt + b'18f', self.size))
    for i in range(nelements):
        edata = data[n:n+ntotal]
        n += ntotal
        out = struct1.unpack(edata)
        (eid_device,
         s1ar, s2ar, s3ar, s4ar, axialr,
         s1ai, s2ai, s3ai, s4ai, axiali,
         s1br, s2br, s3br, s4br,
         s1bi, s2bi, s3bi, s4bi) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)
        if self.is_debug_file:
            self.binary_debug.write('  eid=%i; C%i=[%s]\n' % (
                eid, i, ', '.join(['%r' % di for di in out])))
        if is_magnitude_phase:
            s1a = polar_to_real_imag(s1ar, s1ai)
            s1b = polar_to_real_imag(s1br, s1bi)
            s2a = polar_to_real_imag(s2ar, s2ai)
            s2b = polar_to_real_imag(s2br, s2bi)
            s3a = polar_to_real_imag(s3ar, s3ai)
            s3b = polar_to_real_imag(s3br, s3bi)
            s4a = polar_to_real_imag(s4ar, s4ai)
            s4b = polar_to_real_imag(s4br, s4bi)
            axial = polar_to_real_imag(axialr, axiali)
        else:
            s1a = complex(s1ar, s1ai)
            s1b = complex(s1br, s1bi)
            s2a = complex(s2ar, s2ai)
            s2b = complex(s2br, s2bi)
            s3a = complex(s3ar, s3ai)
            s3b = complex(s3br, s3bi)
            s4a = complex(s4ar, s4ai)
            s4b = complex(s4br, s4bi)
            axial = complex(axialr, axiali)

        obj.add_new_eid_sort1(dt, eid,
                              s1a, s2a, s3a, s4a, axial,
                              s1b, s2b, s3b, s4b)
    return n

def oes_cbar_random_10(self, data: bytes,
                       obj: Union[RandomBarStressArray, RandomBarStrainArray],
                       nelements: int, ntotal: int) -> int:
    n = 0
    #print(self.code_information())
    #print('self._analysis_code_fmt =', self._analysis_code_fmt)
    struct1 = Struct(self._endian + self._analysis_code_fmt + b'9f')
    add_sort_x = getattr(obj, 'add_sort' + str(self.sort_method))

    #self.log.info('self.nonlinear_factor = %s' % self.nonlinear_factor)
    #assert self.sort_method == 2, self.code_information()
    #if sort_method == 2:
        #obj.node_id = 42
    nonlinear_factor = self.nonlinear_factor
    #print(f'CBAR: nelements={nelements}')
    for i in range(nelements):
        edata = data[n:n+ntotal]
        n += ntotal

        out = struct1.unpack(edata)
        (eid_device,
         s1a, s2a, s3a, s4a, axial,
         s1b, s2b, s3b, s4b) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, nonlinear_factor, self.sort_method)
        #print(f'  eid_device={eid_device} eid={eid} dt={nonlinear_factor} nf={nonlinear_factor} -> {obj.data.shape}')
        #continue
        #print('  eid=%i; C%i=[%s]\n' % (eid, i, ', '.join(['%r' % di for di in out])))
        if self.table_name_str == 'OESXRMS1':
            #assert sort_method == 2
            assert self.sort_method == 1, self.code_information()

        if self.is_debug_file:
            self.binary_debug.write('  eid=%i; C%i=[%s]\n' % (
                eid, i, ', '.join(['%r' % di for di in out])))

        assert eid > 0, "dt=%s eid=%s" % (dt, eid)
        add_sort_x(dt, eid,
                   s1a, s2a, s3a, s4a, axial,
                   s1b, s2b, s3b, s4b)
    return n

def oes_cbush_real_7(self, data: bytes,
                     obj: Union[RealBushStressArray, RealBushStrainArray],
                     nelements: int, ntotal: int, dt) -> int:
    n = 0
    struct1 = Struct(self._endian + mapfmt(self._analysis_code_fmt + b'6f', self.size))
    add_sort_x = getattr(obj, 'add_sort' + str(self.sort_method))
    nonlinear_factor = self.nonlinear_factor
    #print(add_sort_x)
    #print('obj.is_sort1 =', obj.is_sort1, obj.table_name)
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        n += ntotal

        out = struct1.unpack(edata)  # num_wide=7
        if self.is_debug_file:
            self.binary_debug.write('CBUSH-102 - %s\n' % str(out))

        (eid_device, tx, ty, tz, rx, ry, rz) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, nonlinear_factor, self.sort_method)
        #print(f'CBUSH: eid_device={eid_device} eid={eid} dt={nonlinear_factor} nf={nonlinear_factor} -> {obj.data.shape}')

        add_sort_x(dt, eid, tx, ty, tz, rx, ry, rz)
    return n

def oes_cbush_complex_13(self,
                         data: bytes,
                         obj: Union[ComplexCBushStressArray, ComplexCBushStrainArray],
                         nelements: int, ntotal: int,
                         is_magnitude_phase: bool) -> int:
    n = 0
    struct1 = Struct(self._endian + mapfmt(self._analysis_code_fmt + b'12f', self.size))
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)  # num_wide=7
        if self.is_debug_file:
            self.binary_debug.write('CBUSH-102 - %s\n' % str(out))

        (eid_device,
         txr, tyr, tzr, rxr, ryr, rzr,
         txi, tyi, tzi, rxi, ryi, rzi) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)

        if is_magnitude_phase:
            tx = polar_to_real_imag(txr, txi)
            ty = polar_to_real_imag(tyr, tyi)
            tz = polar_to_real_imag(tzr, tzi)
            rx = polar_to_real_imag(rxr, rxi)
            ry = polar_to_real_imag(ryr, ryi)
            rz = polar_to_real_imag(rzr, rzi)
        else:
            tx = complex(txr, txi)
            ty = complex(tyr, tyi)
            tz = complex(tzr, tzi)
            rx = complex(rxr, rxi)
            ry = complex(ryr, ryi)
            rz = complex(rzr, rzi)
        obj.add_sort1(dt, eid, tx, ty, tz, rx, ry, rz)
        n += ntotal
    return ntotal

def oes_ctriax6_real_33(self, data: bytes,
                        obj: Union[RealTriaxStressArray, RealTriaxStrainArray],
                        nelements: int, ntotal: int, dt: Any) -> int:
    n = 0
    ntotal1 = 36 * self.factor
    ntotal2 = 32 * self.factor
    s1 = Struct(self._endian + b'2i7f')  # 36
    s2 = Struct(self._endian + b'i7f')
    for unused_i in range(nelements):
        out = s1.unpack(data[n:n + ntotal1])
        (eid_device, loc, rs, azs, As, ss, maxp, tmax, octs) = out
        if self.is_debug_file:
            self.binary_debug.write('CTRIAX6-53A - %s\n' % (str(out)))
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)

        obj.add_sort1(dt, eid, loc, rs, azs, As, ss, maxp, tmax, octs)
        n += ntotal1
        for unused_j in range(3):
            out = s2.unpack(data[n:n + ntotal2])
            (loc, rs, azs, As, ss, maxp, tmax, octs) = out
            if self.is_debug_file:
                self.binary_debug.write('CTRIAX6-53B - %s\n' % (str(out)))
            obj.add_sort1(dt, eid, loc, rs, azs, As, ss, maxp, tmax, octs)
            n += ntotal2
    return n

def oes_ctriax_complex_37(self,
                          data: bytes,
                          obj: ComplexTriaxStressArray,
                          nelements: int,
                          is_magnitude_phase: bool) -> int:
    n = 0
    ntotal1 = 40 * self.factor
    ntotal2 = 36 * self.factor
    s1 = Struct(self._endian + mapfmt(self._analysis_code_fmt + b'i8f', self.size)) # 10*4 = 40
    s2 = Struct(self._endian + mapfmt(b'i8f', self.size))  #  9*4 = 36

    add_sort_x = getattr(obj, 'add_sort' + str(self.sort_method))
    for unused_i in range(nelements):
        out = s1.unpack(data[n:n + ntotal1])
        (eid_device, loc, rsr, rsi, azsr, azsi, Asr, Asi, ssr, ssi) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)
        if self.is_debug_file:
            self.binary_debug.write('CTRIAX6-53 eid=%i\n    %s\n' % (eid, str(out)))
        #print('CTRIAX6-53 eid=%i\n    %s\n' % (eid, str(out)))

        if is_magnitude_phase:
            rs = polar_to_real_imag(rsr, rsi)
            azs = polar_to_real_imag(azsr, azsi)
            As = polar_to_real_imag(Asr, Asi)
            ss = polar_to_real_imag(ssr, ssi)
        else:
            rs = complex(rsr, rsi)
            azs = complex(azsr, azsi)
            As = complex(Asr, Asi)
            ss = complex(ssr, ssi)
        obj.add_element_sort1(dt, eid)
        obj.add_sort1(dt, eid, loc, rs, azs, As, ss)

        n += ntotal1
        for unused_j in range(3):
            out = s2.unpack(data[n:n + ntotal2])
            (loc, rsr, rsi, azsr, azsi, Asr, Asi, ssr, ssi) = out
            if self.is_debug_file:
                self.binary_debug.write('    %s\n' % (str(out)))
            #print("eid=%s loc=%s rs=%s azs=%s as=%s ss=%s" % (
                #eid, loc, rs, azs, As, ss))

            if is_magnitude_phase:
                rs = polar_to_real_imag(rsr, rsi)
                azs = polar_to_real_imag(azsr, azsi)
                As = polar_to_real_imag(Asr, Asi)
                ss = polar_to_real_imag(ssr, ssi)
            else:
                rs = complex(rsr, rsi)
                azs = complex(azsr, azsi)
                As = complex(Asr, Asi)
                ss = complex(ssr, ssi)
            add_sort_x(dt, eid, loc, rs, azs, As, ss)
            n += ntotal2  # 4*8
    return n

def oes_csolid_real(self, data: bytes,
                    obj: Union[RealSolidStressArray, RealSolidStrainArray],
                    nelements: int, dt: Any,
                    element_name: str, nnodes_expected: int,
                    preline1: str, preline2: str) -> int:
    """
    reads stress/strain for element type:
     - 39 : CTETRA
     - 67 : CHEXA
     - 68 : CPENTA

    """
    add_eid_sort_x = getattr(obj, 'add_eid_sort' + str(self.sort_method))
    add_node_sort_x = getattr(obj, 'add_node_sort' + str(self.sort_method))
    n = 0
    if self.size == 4:
        fmt1 = self._endian + self._analysis_code_fmt + b'i4si'
        fmt2 = self._endian + b'i20f'
    else:
        fmt1 = self._endian + mapfmt(self._analysis_code_fmt, self.size) + b'q8sq'
        fmt2 = self._endian + b'q20d'
    struct1 = Struct(fmt1)
    struct2 = Struct(fmt2)
    if self.is_debug_file:
        msg = '%s-%s nelements=%s nnodes=%s; C=[sxx, sxy, s1, a1, a2, a3, pressure, svm,\n' % (
            self.element_name, self.element_type, nelements, nnodes_expected)
        msg += '                                 syy, syz, s2, b1, b2, b3,\n'
        msg += '                                 szz, sxz, s3, c1, c2, c3]\n'
        self.binary_debug.write(msg)

    n16 = 16 * self.factor
    n84 = 84 * self.factor
    #if self.is_sort1:
        #print('nelements =', nelements)
    #else:
        #print('ntimes =', nelements)
    for unused_i in range(nelements):
        edata = data[n:n+n16]
        out = struct1.unpack(edata)
        (eid_device, cid, unused_abcd, nnodes) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)
        #print(f'eid_device={eid_device} eid={eid} dt={dt}')
        if self.is_debug_file:
            self.binary_debug.write('%s - eid=%i; %s\n' % (preline1, eid, str(out)))

        assert nnodes < 21, 'print_block(data[n:n+16])'  #self.print_block(data[n:n+16])

        n += n16
        for inode in range(nnodes_expected):  # nodes pts, +1 for centroid (???)
            out = struct2.unpack(data[n:n + n84]) # 4*21 = 84
            if self.is_debug_file:
                self.binary_debug.write('%s - %s\n' % (preline2, str(out)))
            (grid_device,
             sxx, sxy, s1, a1, a2, a3, pressure, svm,
             syy, syz, s2, b1, b2, b3,
             szz, sxz, s3, c1, c2, c3) = out

            if self.is_debug_file:
                self.binary_debug.write('  eid=%s inode=%i; C=[%s]\n' % (
                    eid, grid_device, ', '.join(['%r' % di for di in out])))

            #if grid_device == 0:
                #grid = 'CENTER'
            #else:
                ##grid = (grid_device - device_code) // 10
                #grid = grid_device

            grid = grid_device
            a_cos = [a1, a2, a3]
            b_cos = [b1, b2, b3]
            c_cos = [c1, c2, c3]

            if inode == 0:
                #  this is correct, but fails
                #element_name = self.element_name + str(nnodes)
                add_eid_sort_x(element_name, cid, dt, eid, grid,
                               sxx, syy, szz, sxy, syz, sxz, s1, s2, s3,
                               a_cos, b_cos, c_cos, pressure, svm)
            else:
                add_node_sort_x(dt, eid, inode, grid,
                                sxx, syy, szz, sxy, syz, sxz, s1, s2, s3,
                                a_cos, b_cos, c_cos, pressure, svm)
            n += n84
    return n

def oes_csolid_complex(self, data: bytes,
                       obj: Union[ComplexSolidStressArray, ComplexSolidStrainArray],
                       nelements: int, # nnodes: int,
                       element_name: str, nnodes_expected: int,
                       is_magnitude_phase: bool) -> int:
    n = 0
    if self.size == 4:
        s1 = Struct(self._endian + b'2i4si')
    else:
        s1 = Struct(self._endian + b'2q8sq')
    s2 = Struct(self._endian + mapfmt(b'i12f', self.size))
    ntotal1 = 16 * self.factor
    ntotal2 = 52 * self.factor
    for unused_i in range(nelements):
        edata = data[n:n+ntotal1]
        n += ntotal1
        out = s1.unpack(edata)
        (eid_device, cid, ctype, nodef) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)
        if self.is_debug_file:
            self.binary_debug.write('  eid=%i C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))

        #element_name = self.element_name + str(nodef)  # this is correct, but has problems...
        obj.add_eid_sort1(self.element_type, element_name, dt, eid, cid, ctype, nodef)
        for inode in range(nnodes_expected):
            edata = data[n:n+ntotal2]
            n += ntotal2
            out = s2.unpack(edata)
            (grid,
             exr, eyr, ezr, etxyr, etyzr, etzxr,
             exi, eyi, ezi, etxyi, etyzi, etzxi) = out
            #if grid == 0:
                #grid = 'CENTER'

            if is_magnitude_phase:
                ex = polar_to_real_imag(exr, exi)
                ey = polar_to_real_imag(eyr, eyi)
                ez = polar_to_real_imag(ezr, ezi)
                etxy = polar_to_real_imag(etxyr, etxyi)
                etyz = polar_to_real_imag(etyzr, etyzi)
                etzx = polar_to_real_imag(etzxr, etzxi)
            else:
                ex = complex(exr, exi)
                ey = complex(eyr, eyi)
                ez = complex(ezr, ezi)
                etxy = complex(etxyr, etxyi)
                etyz = complex(etyzr, etyzi)
                etzx = complex(etzxr, etzxi)

            if self.is_debug_file:
                self.binary_debug.write('       node%s=[%s]\n' % (
                    grid, ', '.join(['%r' % di for di in out])))
            obj.add_node_sort1(dt, eid, grid, inode,
                               ex, ey, ez, etxy, etyz, etzx)
    return n

def oes_csolid_random(self, data: bytes,
                      obj: Union[RandomSolidStressArray, RandomSolidStrainArray],
                      nelements: int,
                      element_name: str, nnodes_expected: int,
                      preline1: str, preline2: str) -> int:
    n = 0
    struct1 = Struct(self._endian + self._analysis_code_fmt + b'i4si')
    struct2 = Struct(self._endian + b'i6f')
    if self.is_debug_file and 0:
        msg = '%s-%s nelements=%s nnodes=%s; C=[sxx, sxy, s1, a1, a2, a3, pressure, svm,\n' % (
            self.element_name, self.element_type, nelements, nnodes_expected)
        msg += '                                 syy, syz, s2, b1, b2, b3,\n'
        msg += '                                 szz, sxz, s3, c1, c2, c3]\n'
        self.binary_debug.write(msg)

    ntotal1 = 16 * self.factor
    ntotal2 = 28 * self.factor
    for unused_i in range(nelements):
        edata = data[n:n+ntotal1]
        out = struct1.unpack(edata)
        (eid_device, cid, unused_abcd, grid) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)
        assert eid > 0, eid

        if self.is_debug_file and 0:
            self.binary_debug.write('%s - eid=%i; %s\n' % (preline1, eid, str(out)))

        #assert nnodes < 21, 'print_block(data[n:n+16])'  #self.print_block(data[n:n+16])

        n += ntotal1
        for inode in range(nnodes_expected):  # nodes pts, +1 for centroid (???)
            #self.show_data(data[n:n+48])
            out = struct2.unpack(data[n:n + ntotal2]) # 4*7 = 28
            if self.is_debug_file:
                self.binary_debug.write('%s - %s\n' % (preline2, str(out)))
            (grid_device, sxx, syy, szz, txy, tyz, txz) = out

            if self.is_debug_file:
                self.binary_debug.write('  eid=%s inode=%i; C=[%s]\n' % (
                    eid, grid_device, ', '.join(['%r' % di for di in out])))

            #if grid_device == 0:
                #grid = 'CENTER'
            #else:
                ##grid = (grid_device - device_code) // 10
                #grid = grid_device

            grid = grid_device
            if inode == 0:
                #  this is correct, but fails
                #element_name = self.element_name + str(nnodes)
                obj.add_eid_sort1(element_name, cid, dt, eid, grid,
                                  sxx, syy, szz, txy, tyz, txz)
            else:
                obj.add_node_sort1(dt, eid, inode, grid,
                                   sxx, syy, szz, txy, tyz, txz)
            n += ntotal2
    return n

def oes_comp_shell_real_11(self, data: bytes, ndata: int,
                           obj: Union[RealCompositePlateStressArray, RealCompositePlateStrainArray],
                           ntotal: int, nelements: int, etype: str, dt: Any) -> int:
    n = 0
    eid_old = 0
    struct1 = Struct(self._endian + mapfmt(self._analysis_code_fmt + b'i9f', self.size)) # 11
    sort_method = self.sort_method
    add_eid_sort_x = getattr(obj, 'add_eid_sort' + str(self.sort_method))
    add_sort_x = getattr(obj, 'add_sort' + str(self.sort_method))

    for unused_i in range(nelements):
        edata = data[n:n+ntotal]  # 4*11
        out = struct1.unpack(edata)
        (eid_device, layer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, sort_method)

        if self.is_debug_file:
            self.binary_debug.write('  eid=%i; layer=%i; C=[%s]\n' % (eid, layer, ', '.join(['%r' % di for di in out])))

        if eid != eid_old:
            # originally initialized to None, the buffer doesnt reset it, so it is the old value
            add_eid_sort_x(etype, dt, eid, layer, o1, o2,
                           t12, t1z, t2z, angle, major, minor, ovm)
        else:
            add_sort_x(dt, eid, layer, o1, o2,
                       t12, t1z, t2z, angle, major, minor, ovm)
        eid_old = eid
        n += ntotal
    return n

def oes_cshear_random_3(self, data: bytes,
                        obj: Union[RandomShearStressArray, RandomShearStrainArray],
                        nelements: int, ntotal: int) -> int:
    n = 0
    add_sort_x = getattr(obj, 'add_sort' + str(self.sort_method))
    struct1 = Struct(self._endian + self._analysis_code_fmt + b'2f')
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)  # num_wide=5
        if self.is_debug_file:
            self.binary_debug.write('CSHEAR-4 - %s\n' % str(out))

        (eid_device, max_strain, avg_strain) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)
        #print(eid, dt)
        add_sort_x(dt, eid, max_strain, avg_strain)
        n += ntotal
    return n

def oes_cbend_real_21(self, data: bytes,
                      obj: Union[RealBendStressArray, RealBendStrainArray],
                      nelements: int, ntotal: int, dt) -> int:
    n = 0
    ntotali = 40 * self.factor
    struct1 = Struct(self._endian + self._analysis_code_fmt)
    struct2 = Struct(self._endian + b'i9f')
    add_sort_x = getattr(obj, 'add_sort' + str(self.sort_method))

    #print('ntimes =', nelements)
    for unused_i in range(nelements):
        edata = data[n:n + 4]
        eid_device, = struct1.unpack(edata)
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)

        n += 4
        for unused_i in range(2):
            edata = data[n:n + ntotali]
            out = struct2.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('BEND-69 - eid=%s %s\n' % (eid, str(out)))
            #print('BEND-69 - eid=%s %s\n' % (eid, str(out)))

            (grid, angle, sc, sd, se, sf, omax, omin, mst, msc) = out

            add_sort_x(dt, eid, grid, angle, sc, sd, se, sf, omax, omin, mst, msc)
            n += ntotali
    return n


def oesrt_cquad4_95(self, data: bytes, ndata: int) -> int:
    """unsupported element"""
    assert self.num_wide == 9, f'num_wide={self.num_wide} not 9'
    ntotal = 36 * self.factor # 4*9
    #oesrt_cquad4_95

    n = 0
    struct1 = Struct(self._endian + self._analysis_code_fmt + b'8si3fi4s')
    nelements = ndata // ntotal
    #obj = self.obj
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)  # num_wide=9
        if self.is_debug_file:
            self.binary_debug.write('CQUAD4-95 - %s\n' % str(out))
        #eid, failure, ply, failureIndexPly, failureIndexBonding, failureIndexMax, flag
        # 3,TSAIWU,1,8.5640,0.0,None

        (eid, failure, ply, strength_ratio_ply, failure_index_bonding, strength_ratio_bonding, flag, flag2) = out
        #strength_ratio_ply
        #print("eid=%s failure=%r ply=%s failureIndexPly=%s  failure_index_bonding=%s strength_ratio_bonding=%s flag=%s flag2=%s" % (
        #    eid, failure.strip(), ply, failureIndexPly, failure_index_bonding, strength_ratio_bonding, flag, flag2))
        #print("eid=%s strength_ratio_ply=%g failure_index_bonding=%s strength_ratio_bonding=%s" % (
            #eid, strength_ratio_ply, failure_index_bonding, strength_ratio_bonding))
        #obj.add_new_eid(element_name, dt, eid, force, stress)
        n += ntotal
    return n

def oes_csolid_composite_real(self, data: bytes,
                              obj,
                              nelements: int, nedges: int,
                              element_name: str, preline1: str, preline2: str,
                              dt: Any) -> int:
    """
    1 PLY I Lamina number
    2 FLOC CHAR4 Fiber location (BOT, MID, TOP)
    3 GRID I Edge grid ID (center=0)

    4 EX1 RS Normal strain in the 1-direction
    5 EY1 RS Normal strain in the 2-direction
    6 EZ1 RS Normal strain in the 3-direction
    7 ET1 RS Shear strain in the 12-plane
    8 EL2 RS Shear strain in the 23-plane
    9 EL1 RS Shear strain in the 13-plane
    10 ETMAX1 RS von Mises strain
    For each fiber location requested (PLSLOC), words 3 through 10 repeat 4 times.
    """
    n = 0
    struct1 = Struct(self._endian + self._analysis_code_fmt + b'i4s')
    struct2 = Struct(self._endian + b'i7f')
    if self.is_debug_file:
        msg = '%s-%s nelements=%s nnodes=%s; C=[sxx, syy, szz, txy, tyz, txz, ovm,\n' % (
            self.element_name, self.element_type, nelements, nedges)
        self.binary_debug.write(msg)

    ntotal1 = 12 * self.factor
    ntotal2 = 32 * self.factor
    for unused_i in range(nelements):
        edata = data[n:n+ntotal1]
        out = struct1.unpack(edata)
        (eid_device, ply, fiber_location) = out
        #print(out)
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, self.sort_method)

        if self.is_debug_file:
            self.binary_debug.write('%s - eid=%i; %s\n' % (preline1, eid, str(out)))

        n += ntotal1
        for inode in range(nedges):  # nodes pts, no centroid
            out = struct2.unpack(data[n:n + ntotal2]) # 4*8 = 32
            if self.is_debug_file:
                self.binary_debug.write('%s - %s\n' % (preline2, str(out)))
            (grid_device, sxx, syy, szz, txy, tyz, txz, ovm) = out

            if self.is_debug_file:
                self.binary_debug.write('  eid=%s inode=%i; C=[%s]\n' % (
                    eid, grid_device, ', '.join(['%r' % di for di in out])))

            #if grid_device == 0:
                #grid = 'CENTER'
            #else:
                ##grid = (grid_device - device_code) // 10
                #grid = grid_device

            grid = grid_device
            #a_cos = [a1, a2, a3]
            #b_cos = [b1, b2, b3]
            #c_cos = [c1, c2, c3]
            if 0:
                if inode == 0:
                    #  this is correct, but fails
                    #element_name = self.element_name + str(nnodes)
                    obj.add_eid_sort1(element_name, cid, dt, eid, grid,
                                      sxx, syy, szz, txy, tyz, txz, ovm)
                else:
                    obj.add_node_sort1(dt, eid, inode, grid,
                                       sxx, syy, szz, txy, tyz, txz, ovm)
            n += ntotal2
    return n

def oes_shell_composite_complex_11(self,
                                   data: bytes,
                                   obj: ComplexLayeredCompositesArray,
                                   ntotal: int, nelements: int, sort_method: int,
                                   dt: Any, is_magnitude_phase: bool) -> int:
    """OESCP, OESTRCP"""
    n = 0
    struct1 = Struct(self._endian + self._analysis_code_fmt + b'i9f')
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device, ply_id, oxx, oyy, txy, txz, tyz, angle, omax, omin, max_shear) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, sort_method)
        #print(eid, out)

        #if self.is_debug_file:
            #self.binary_debug.write('%s-%s - (%s) + %s\n' % (self.element_name, self.element_type, eid_device, str(out)))
        #obj.add_sort1(dt, eid, theory, lamid, fp, fm, fb, fmax, fflag)
        obj.add_sort1(dt, eid, ply_id, oxx, oyy, txy, txz, tyz, angle, omax, omin, max_shear)
        n += ntotal
    return n

def oes_shell_composite_complex_13(self,
                                   data: bytes,
                                   obj: ComplexLayeredCompositesArray,
                                   ntotal: int, nelements: int, sort_method: int,
                                   dt: Any, is_magnitude_phase: bool) -> int:
    """OESVM1C, OSTRVM1C"""
    # OESCP - STRAINS IN LAYERED COMPOSITE ELEMENTS (QUAD4)
    n = 0
    struct1 = Struct(self._endian + self._analysis_code_fmt + b'i9f ff')
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device, ply_id,
         o1a, o2a, t12a, o1za, o2za,
         o1b, o2b, t12b, o1zb, e2zb, ovm,) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, self.nonlinear_factor, sort_method)
        #print(eid, out)

        #print('%s-%s - (%s) + %s\n' % (self.element_name, self.element_type, eid_device, str(out)))
        #if self.is_debug_file:
            #self.binary_debug.write('%s-%s - (%s) + %s\n' % (self.element_name, self.element_type, eid_device, str(out)))
        obj.add_sort1(dt, eid, ply_id,
                      o1a, o2a, t12a, o1za, o2za,
                      o1b, o2b, t12b, o1zb, e2zb, ovm)
        n += ntotal
    return n
