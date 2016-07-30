#pylint: disable=C0301,W0201,R0911,R0915,R0914
"""
Defines the Real/Complex Stresses/Strains created by:
    STRESS = ALL
    STRAIN = ALL
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from struct import Struct
from six import b
from six.moves import range
from numpy import fromstring, radians, sin, cos, vstack, repeat, array
import numpy as np

from pyNastran.op2.op2_common import OP2Common
from pyNastran.op2.op2_helper import polar_to_real_imag

from pyNastran.op2.tables.oes_stressStrain.real.oes_bars import RealBarStressArray, RealBarStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_bars100 import RealBar10NodesStressArray, RealBar10NodesStrainArray

from pyNastran.op2.tables.oes_stressStrain.real.oes_beams import (RealBeamStressArray, RealBeamStrainArray,
                                                                  RealNonlinearBeamStressArray)
from pyNastran.op2.tables.oes_stressStrain.real.oes_bush import RealBushStressArray, RealBushStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_bush1d import RealBush1DStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_composite_plates import RealCompositePlateStressArray, RealCompositePlateStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_gap import NonlinearGapStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_plates import (RealPlateStressArray, RealPlateStrainArray,
                                                                   RealCPLSTRNPlateStressArray, RealCPLSTRNPlateStrainArray)
from pyNastran.op2.tables.oes_stressStrain.real.oes_rods import RealRodStressArray, RealRodStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_shear import RealShearStrainArray, RealShearStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_solids import RealSolidStrainArray, RealSolidStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_springs import (RealSpringStressArray, RealSpringStrainArray,
                                                                    RealNonlinearSpringStressArray)
from pyNastran.op2.tables.oes_stressStrain.real.oes_triax import RealTriaxStressArray, RealTriaxStrainArray


from pyNastran.op2.tables.oes_stressStrain.complex.oes_bars import ComplexBarStressArray, ComplexBarStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_bush import (ComplexCBushStressArray, ComplexCBushStrainArray)
from pyNastran.op2.tables.oes_stressStrain.complex.oes_bush1d import ComplexCBush1DStressArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_plates import ComplexPlateStressArray, ComplexPlateStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_rods import ComplexRodStressArray, ComplexRodStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_shear import ComplexShearStressArray, ComplexShearStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_solids import ComplexSolidStressArray, ComplexSolidStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_springs import ComplexSpringStressArray, ComplexSpringStrainArray

# TODO: vectorize 1
from pyNastran.op2.tables.oes_stressStrain.oes_nonlinear import (RealNonlinearRodArray,
                                                                 HyperelasticQuad, # vectorize
                                                                 RealNonlinearPlateArray)


class OES(OP2Common):
    """
    Defines  the OES class that is used to read stress/strain data
    """
    def __init__(self):
        OP2Common.__init__(self)
        self.ntotal = 0

    #def _oes_cleanup():
    def _read_oes1_3(self, data, ndata):
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
            '???', 'Title', 'subtitle', 'label']

        self.parse_approach_code(data)  # 3

        ## element type
        self.element_type = self.add_data_parameter(data, 'element_type', 'i', 3, False)

        ## load set ID
        self.load_set = self.add_data_parameter(data, 'load_set', 'i', 8, False)

        ## format code
        self.format_code = self.add_data_parameter(data, 'format_code', 'i', 9, False)

        ## number of words per entry in record
        ## .. note:: is this needed for this table ???
        self.num_wide = self.add_data_parameter(data, 'num_wide', 'i', 10, False)

        ## stress/strain codes
        self.s_code = self.add_data_parameter(data, 's_code', 'i', 11, False)

        ## thermal flag; 1 for heat ransfer, 0 otherwise
        self.thermal = self.add_data_parameter(data, 'thermal', 'i', 23, False)

        ## assuming tCode=1
        if self.analysis_code == 1:   # statics / displacement / heat flux
            ## load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5, False)
            self.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysis_code == 2:  # real eigenvalues
            #: mode number
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)
            #: real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            #: mode or cycle TODO confused on the type - F1 means float/int???
            self.mode2 = self.add_data_parameter(data, 'mode2', 'i', 7, False)
            self.cycle = self.add_data_parameter(data, 'cycle', 'f', 7, False)
            self.update_mode_cycle('cycle')
            self.data_names = self.apply_data_code_value('data_names', ['mode', 'eigr', 'mode2', 'cycle'])
        #elif self.analysis_code==3: # differential stiffness
            #self.lsdvmn = self.get_values(data,'i',5) ## load set number
            #self.data_code['lsdvmn'] = self.lsdvmn
        #elif self.analysis_code==4: # differential stiffness
        #    self.lsdvmn = self.get_values(data,'i',5) ## load set number
        elif self.analysis_code == 5:   # frequency
            ## frequency
            self.freq = self.add_data_parameter(data, 'freq', 'f', 5)
            self.data_names = self.apply_data_code_value('data_names', ['freq'])
        elif self.analysis_code == 6:  # transient
            ## time step
            self.dt = self.add_data_parameter(data, 'dt', 'f', 5)
            self.data_names = self.apply_data_code_value('data_names', ['dt'])
        elif self.analysis_code == 7:  # pre-buckling
            ## load set
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
        elif self.analysis_code == 8:  # post-buckling
            ## mode number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)  # real eigenvalue
            self.data_names = self.apply_data_code_value('data_names', ['lsdvmn', 'eigr'])
        elif self.analysis_code == 9:  # complex eigenvalues
            ## mode number
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            ## imaginary eigenvalue
            self.eigi = self.add_data_parameter(data, 'eigi', 'f', 7, False)
            self.data_names = self.apply_data_code_value('data_names', ['mode', 'eigr', 'eigi'])
        elif self.analysis_code == 10:  # nonlinear statics
            ## load step
            self.lftsfq = self.add_data_parameter(data, 'lftsfq', 'f', 5)
            self.data_names = self.apply_data_code_value('data_names', ['lftsfq'])
        elif self.analysis_code == 11:  # old geometric nonlinear statics
            ## load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
        elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## Time step ??? --> straight from DMAP
            self.dt = self.add_data_parameter(data, 'dt', 'f', 5)
            self.data_names = self.apply_data_code_value('data_names', ['dt'])
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % self.analysis_code
            raise RuntimeError(msg)
        # tCode=2
        #if self.analysis_code==2: # sort2
        #    self.lsdvmn = self.get_values(data,'i',5)

        self.fix_format_code()
        self._parse_thermal_code()
        try:
            self.element_name = self.element_mapper[self.element_type]
        except KeyError:
            self.log.error(self.code_information())
            raise
        self.data_code['element_name'] = self.element_name
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

        self._parse_stress_code()
        self._write_debug_bits()
        if self.is_debug_file:
            self.binary_debug.write('catttt!')
        assert isinstance(self.format_code, int), self.format_code
        #print('self.nonlinear_factor =', self.nonlinear_factor)
        #assert self.num_wide != 146, self.code_information()

    def _parse_stress_code(self):
        """
        s_code =  0 -> stress_bits = [0,0,0,0,0]
        s_code =  1 -> stress_bits = [0,0,0,0,1]
        s_code =  2 -> stress_bits = [0,0,0,1,0]
        s_code =  3 -> stress_bits = [0,0,0,1,1]
        etc.
        s_code = 32 -> stress_bits = [1,1,1,1,1]

        stress_bits[0] = 0 -> isMaxShear=True       isVonMises=False
        stress_bits[0] = 1 -> isMaxShear=False      isVonMises=True

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

    def _read_oes1_4(self, data, ndata):
        """
        Reads the Stress Table 4
        """
        #assert self.isubtable == -4, self.isubtable
        #if self.is_debug_file:
            #self.binary_debug.write('  element_name = %r\n' % self.element_name)
        #print "element_name =", self.element_name
        assert isinstance(self.format_code, int), self.format_code
        assert self.is_stress() == True, self.code_information()
        self.data_code['is_stress_flag'] = True
        self.data_code['is_strain_flag'] = False

        if self.isubcase not in self.case_control_deck.subcases:
            self.subcase = self.case_control_deck.create_new_subcase(self.isubcase)
        self.subcase.add_op2_data(self.data_code, 'STRESS/STRAIN', self.log)

        if self.is_sort1():
            n = self._read_oes1_4_sort1(data, ndata)
        else:
            msg = self.code_information()
            n = self._not_implemented_or_skip(data, ndata, msg)
        return n

    def _read_oes2_3(self, data, ndata):
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
            '???', 'Title', 'subtitle', 'label']

        self.parse_approach_code(data)  # 3

        ## element type
        self.element_type = self.add_data_parameter(data, 'element_type', 'i', 3, False)

        ## load set ID
        self.load_set = self.add_data_parameter(data, 'load_set', 'i', 8, False)

        ## format code
        self.format_code = self.add_data_parameter(data, 'format_code', 'i', 9, False)

        ## number of words per entry in record
        ## .. note:: is this needed for this table ???
        self.num_wide = self.add_data_parameter(data, 'num_wide', 'i', 10, False)

        ## stress/strain codes
        self.s_code = self.add_data_parameter(data, 's_code', 'i', 11, False)

        ## thermal flag; 1 for heat ransfer, 0 otherwise
        self.thermal = self.add_data_parameter(data, 'thermal', 'i', 23, False)

        if self.analysis_code == 2:  # real eigenvalues
            self._analysis_code_fmt = 'i'
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            self.mode_cycle = self.add_data_parameter(data, 'mode_cycle', 'i', 7, False)  # mode or cycle .. todo:: confused on the type - F1???
            self.data_names = self.apply_data_code_value('data_names', ['node_id', 'eigr', 'mode_cycle'])
        elif self.analysis_code == 5:   # frequency
            self._analysis_code_fmt = 'f'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
        elif self.analysis_code == 6:  # transient
            self._analysis_code_fmt = 'f'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
        elif self.analysis_code == 7:  # pre-buckling
            self._analysis_code_fmt = 'i'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
        elif self.analysis_code == 8:  # post-buckling
            self._analysis_code_fmt = 'f'
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            self.data_names = self.apply_data_code_value('data_names', ['node_id', 'eigr'])
        elif self.analysis_code == 9:  # complex eigenvalues
            # mode number
            self._analysis_code_fmt = 'i'
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            self.eigi = self.add_data_parameter(data, 'eigi', 'f', 7, False)
            self.data_names = self.apply_data_code_value('data_names', ['node_id', 'eigr', 'eigi'])
        elif self.analysis_code == 10:  # nonlinear statics
            # load step
            self._analysis_code_fmt = 'f'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
        elif self.analysis_code == 11:  # old geometric nonlinear statics
            # load set number
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
        elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % self.analysis_code
            raise RuntimeError(msg)

        self._fix_sort2(data)

    def _fix_sort2(self, data):
        self.fix_format_code()
        if self.num_wide == 8:
            self.format_code = 1
            self.data_code['format_code'] = 1
        else:
            #self.fix_format_code()
            if self.format_code == 1:
                self.format_code = 2
                self.data_code['format_code'] = 2
            assert self.format_code in [2, 3], self.code_information()

        self._parse_thermal_code()
        if self.is_debug_file:
            self.binary_debug.write('  %-14s = %r %s\n' % ('approach_code', self.approach_code,
                                                           self.approach_code_str(self.approach_code)))
            self.binary_debug.write('  %-14s = %r\n' % ('tCode', self.tCode))
            self.binary_debug.write('  %-14s = %r\n' % ('isubcase', self.isubcase))
        self._read_title(data)
        self._write_debug_bits()
        #assert isinstance(self.nonlinear_factor, int), self.nonlinear_factor

    def _read_oes2_4(self, data, ndata):
        return self._table_passer(data, ndata)

    def _read_ostr1_4(self, data, ndata):
        """
        Reads the Strain Table 4
        """
        #assert self.isubtable == -4, self.isubtable
        #if self.is_debug_file:
            #self.binary_debug.write('  element_name = %r\n' % self.element_name)
        #print "element_name =", self.element_name
        assert self.is_strain() == True, self.code_information()
        self.data_code['is_stress_flag'] = False
        self.data_code['is_strain_flag'] = True

        if self.isubcase not in self.case_control_deck.subcases:
            self.subcase = self.case_control_deck.create_new_subcase(self.isubcase)
        self.subcase.add_op2_data(self.data_code, 'STRESS/STRAIN', self.log)

        if self.is_sort1():
            n = self._read_ostr1_4_sort1(data, ndata)
        else:
            msg = self.code_information()
            n = self._not_implemented_or_skip(data, ndata, msg)
        return n

    def _autojit3(func):
        """
        Debugging function to print the object name and an needed parameters
        """
        def new_func(self, data):
            """
            The actual function exec'd by the decorated function.
            """
            n = func(self, data)
            return n
        return new_func

    def _print_obj_name_on_crash(func):
        """
        Debugging function to print the object name and an needed parameters
        """
        def new_func(self, data):
            """
            The actual function exec'd by the decorated function.
            """
            try:
                n = func(self, data)
            except NameError:
                raise
            except AttributeError:
                raise
            except:
                raise
                print("----------")
                print(self.obj)
                print(self.data_code)
                if self.obj is not None:
                    #from pyNastran.utils import object_attributes
                    #print object_attributes(self.obj)
                    print(self.obj.data_code)
                print("----------")
                raise
            return n
        return new_func

    #@_print_obj_name_on_crash
    def _read_oes1_4_sort1(self, data, ndata):
        """
        Reads OES1 subtable 4
        """
        #if self.num_wide == 146:
            #assert self.num_wide != 146
            #assert ndata != 146, self.code_information()
        assert isinstance(self.format_code, int), self.format_code
        assert self.is_sort1() == True
        if self.thermal == 0:
            n = self._read_oes1_loads(data, ndata)
        elif self.thermal == 1:
            n = self._read_oes1_thermal(data, ndata)
        else:
            msg = 'thermal=%s' % self.thermal
            n = self._not_implemented_or_skip(data, ndata, msg)
        return n

    #@_print_obj_name_on_crash
    def _read_ostr1_4_sort1(self, data, ndata):
        """
        Reads OSTR1 subtable 4
        """
        #if self.num_wide == 146:
            #assert self.num_wide != 146
            #assert ndata != 146, self.code_information()
        assert self.is_sort1() == True
        if self.thermal == 0:
            n = self._read_oes1_loads(data, ndata)
        elif self.thermal == 1:
            n = self._read_oes1_thermal(data, ndata)
        else:
            msg = 'thermal=%s' % self.thermal
            n = self._not_implemented_or_skip(data, ndata, msg)
        return n

    def _read_oes1_thermal(self, data, ndata):
        """
        Reads OES self.thermal=1 tables; uses a hackish method to just skip the table.
        """
        return ndata

    def _read_ostr1_thermal(self, data, ndata):
        """
        Reads OSTR self.thermal=1 tables; uses a hackish method to just skip the table.
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

            (3, 1, 5, b'OES1X1') : ('ctube_stress', RealRodStressArray),
            (3, 1, 5, b'OES1X') : ('ctube_stress', RealRodStressArray),
            (3, 2, 5) : ('ctube_stress', ComplexRodStressArray),
            (3, 3, 5) : ('ctube_stress', ComplexRodStressArray),

            (10, 1, 5, b'OES1') : ('conrod_stress', RealRodStressArray),
            (10, 1, 5, b'OES1X') : ('conrod_stress', RealRodStressArray),
            (10, 1, 5, b'OES1X1') : ('conrod_stress', RealRodStressArray),
            (10, 2, 5) : ('conrod_stress', ComplexRodStressArray),
            (10, 3, 5) : ('conrod_stress', ComplexRodStressArray),

            # beams
            (2, 1, 111, b'OES1X1') : ('cbeam_stress', RealBeamStressArray),
            (2, 1, 111, b'OES1X') : ('cbeam_stress', RealBeamStressArray),
            (2, 1, 111, b'OES1') : ('cbeam_stress', RealBeamStressArray),
            (2, 2, 111, b'OES1X') : ('cbeam_stress', 'ComplexBeamStressArray'),
            (2, 3, 111, b'OES1X') : ('cbeam_stress', 'ComplexBeamStressArray'),

            (4, 1, 4, b'OES1X1') : ('cshear_stress', RealShearStressArray),
            (4, 2, 5) : ('cshear_stress', ComplexShearStressArray),
            (4, 3, 5) : ('cshear_stress', ComplexShearStressArray),
            #(4, 3, 3) : ('cshear_stress', RandomShearStressArray),

            (11, 1, 2, b'OES1X1') : ('celas1_stress', RealSpringStressArray), # real
            (11, 2, 3, b'OES1X') : ('celas1_stress', ComplexSpringStressArray), # real/imag
            (11, 3, 3, b'OES1X') : ('celas1_stress', ComplexSpringStressArray), # mag/phase

            (12, 1, 2, b'OES1X1') : ('celas2_stress', RealSpringStressArray),
            (12, 1, 2, b'OES1X') : ('celas2_stress', RealSpringStressArray),
            (12, 2, 3, b'OES1X') : ('celas2_stress', ComplexSpringStressArray),
            (12, 3, 3, b'OES1X') : ('celas2_stress', ComplexSpringStressArray),

            (13, 1, 2, b'OES1X1') : ('celas3_stress', RealSpringStressArray),
            (13, 2, 3) : ('celas3_stress', ComplexSpringStressArray),
            (13, 3, 3) : ('celas3_stress', ComplexSpringStressArray),

            (14, 1, 2) : ('celas4_stress', RealSpringStressArray),
            (14, 2, 3) : ('celas4_stress', ComplexSpringStressArray),
            (14, 3, 3) : ('celas4_stress', ComplexSpringStressArray),

            (34, 1, 16, b'OES1X1') : ('cbar_stress', RealBarStressArray),
            (34, 1, 16, b'OES1X') : ('cbar_stress', RealBarStressArray),
            (34, 1, 16, b'OES1') : ('cbar_stress', RealBarStressArray),
            (34, 2, 19, b'OES1X') : ('cbar_stress', ComplexBarStressArray),
            (34, 3, 19, b'OES1X') : ('cbar_stress', ComplexBarStressArray),
            #(34, 1, 19) : ('cbar_stress', RandomBarStressArray),
            (100, 1, 10, b'OES1X1') : ('cbar_stress_10nodes', RealBar10NodesStressArray),
            (100, 1, 10, b'OES1X') : ('cbar_stress_10nodes', RealBar10NodesStressArray),

            # solids
            (39, 1, 109, b'OES1X1') : ('ctetra_stress', RealSolidStressArray), # real
            (39, 1, 109, b'OES1X') : ('ctetra_stress', RealSolidStressArray), # real
            (39, 1, 109, b'OES1') : ('ctetra_stress', RealSolidStressArray), # real

            (67, 1, 193, b'OES1X1') : ('chexa_stress', RealSolidStressArray),
            (67, 1, 193, b'OES1X') : ('chexa_stress', RealSolidStressArray),
            (67, 1, 193, b'OES1') : ('chexa_stress', RealSolidStressArray),

            (68, 1, 151, b'OES1X1') : ('cpenta_stress', RealSolidStressArray),
            (68, 1, 151, b'OES1X') : ('cpenta_stress', RealSolidStressArray),
            (68, 1, 151, b'OES1') : ('cpenta_stress', RealSolidStressArray),

            (39, 2, 69, b'OES1X') : ('ctetra_stress', ComplexSolidStressArray), # real/imag
            #(39, 3, 69) : ('ctetra_stress', ComplexSolidStressArray), # mag/phase

            (67, 2, 121, b'OES1X') : ('chexa_stress', ComplexSolidStressArray),
            (67, 3, 121, b'OES1X') : ('chexa_stress', ComplexSolidStressArray),

            (68, 2, 95, b'OES1X') : ('cpenta_stress', ComplexSolidStressArray),
            (68, 3, 95, b'OES1X') : ('cpenta_stress', ComplexSolidStressArray),

            (33, 1, 17, b'OES1X1') :  ('cquad4_stress', RealPlateStressArray),
            (33, 1, 17, b'OES1X') :  ('cquad4_stress', RealPlateStressArray),
            (33, 1, 17, b'OES1') :  ('cquad4_stress', RealPlateStressArray),
            (33, 2, 15, b'OES1X') :  ('cquad4_stress', ComplexPlateStressArray),
            (33, 3, 15, b'OES1X') :  ('cquad4_stress', ComplexPlateStressArray),
            #(33, 3, 0) :  ('cquad4_stress', RandomPlateStressArray),

            (74, 1, 17, b'OES1X1') : ('ctria3_stress', RealPlateStrainArray),
            (74, 1, 17, b'OES1X') : ('ctria3_stress', RealPlateStrainArray),
            (74, 1, 17, b'OES1') : ('ctria3_stress', RealPlateStrainArray),
            (74, 2, 15, b'OES1X') : ('ctria3_stress', ComplexPlateStrainArray),
            (74, 3, 15, b'OES1X') : ('ctria3_stress', ComplexPlateStrainArray),
            #(74, 1, 9) : ('ctria3_stress', RandomPlateStressArray),

            (82, 1, 87, b'OES1X1') : ('cquadr_stress', RealPlateStressArray),
            (82, 1, 87, b'OES1X') : ('cquadr_stress', RealPlateStressArray),
            (82, 2, 77, b'OES1X') : ('cquadr_stress', ComplexPlateStressArray),
            (82, 3, 77, b'OES1X') : ('cquadr_stress', ComplexPlateStressArray),

            (64, 1, 87, b'OES1X1') : ('cquad8_stress', RealPlateStressArray), # real
            (64, 1, 87, b'OES1X') : ('cquad8_stress', RealPlateStressArray),
            (64, 1, 87, b'OES1') : ('cquad8_stress', RealPlateStressArray),
            (64, 2, 77, b'OES1X') : ('cquad8_stress', ComplexPlateStressArray), # real/imag
            (64, 3, 77, b'OES1X') : ('cquad8_stress', ComplexPlateStressArray), # mag/phase

            (70, 1, 70, b'OES1X1') : ('ctriar_stress', RealPlateStressArray),
            (70, 1, 70, b'OES1X') : ('ctriar_stress', RealPlateStressArray),
            (70, 2, 62, b'OES1X') : ('ctriar_stress', ComplexPlateStressArray),
            (70, 3, 62, b'OES1X') : ('ctriar_stress', ComplexPlateStressArray),

            (75, 1, 70, b'OES1X1') : ('ctria6_stress', RealPlateStressArray),
            (75, 2, 62, b'OES1X') : ('ctria6_stress', ComplexPlateStressArray),
            (75, 3, 62, b'OES1X') : ('ctria6_stress', ComplexPlateStressArray),

            (144, 1, 87, b'OES1X1') : ('cquad4_stress', RealPlateStressArray),
            (144, 1, 87, b'OES1') : ('cquad4_stress', RealPlateStressArray),
            (144, 2, 77, b'OES1X') : ('cquad4_stress', ComplexPlateStressArray),
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
            (53, 2, 37, b'OES1X') : ('ctriax_stress', 'ComplexTriaxStress'),
            (53, 3, 37) : ('ctriax_stress', 'ComplexTriaxStress'),

            (102, 1, 7, b'OES1X1') : ('cbush_stress', RealBushStressArray),
            (102, 1, 7, b'OES1X') : ('cbush_stress', RealBushStressArray),
            (102, 1, 7, b'OES1') : ('cbush_stress', RealBushStressArray),
            (102, 2, 13, b'OES1X') : ('cbush_stress', ComplexCBushStressArray),
            (102, 3, 13, b'OES1X') : ('cbush_stress', ComplexCBushStressArray),

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
            #(234, 1, 11) : ('cquadr', 'NA'), # bad?
            # 233-TRIARLC
            (233, 1, 11, b'OES1C') : ('ctriar', 'NA'),
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
        }

        try:
            return stress_mapper[self.element_type, self.format_code, self.num_wide, self.table_name]
        except KeyError:
            print(self.code_information())
            raise
            return None, None

    def _read_oes1_loads(self, data, ndata):
        """
        Reads OES self.thermal=0 stress/strain
        """
        n = 0
        is_magnitude_phase = self.is_magnitude_phase()
        dt = self.nonlinear_factor

        flag = 'element_id'
        if self.is_stress():
            result_name = 'stress'
            stress_name = 'STRESS'
        else:
            result_name = 'strain'
            stress_name = 'STRAIN'

        if self.is_stress():
            _result_name, _class_obj = self.get_stress_mapper()

        if self._results.is_not_saved(result_name):
            return ndata
        if self.element_type in [1, 3, 10]:  # rods
            # 1-CROD
            # 3-CTUBE
            # 10-CONROD
            if self.is_stress():
                obj_vector_real = RealRodStressArray
                obj_vector_complex = ComplexRodStressArray
                if self.element_type == 1: # CROD
                    result_name = 'crod_stress'
                elif self.element_type == 3:  # CTUBE
                    result_name = 'ctube_stress'
                elif self.element_type == 10:  # CONROD
                    result_name = 'conrod_stress'
                else:
                    msg = self.code_information()
                    return self._not_implemented_or_skip(data, ndata, msg)
            else:
                obj_vector_real = RealRodStrainArray
                obj_vector_complex = ComplexRodStrainArray
                if self.element_type == 1: # CROD
                    result_name = 'crod_strain'
                elif self.element_type == 3:  # CTUBE
                    result_name = 'ctube_strain'
                elif self.element_type == 10:  # CONROD
                    result_name = 'conrod_strain'
                else:
                    msg = self.code_information()
                    return self._not_implemented_or_skip(data, ndata, msg)

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)

            slot = getattr(self, result_name)
            if self.format_code == 1 and self.num_wide == 5:  # real
                ntotal = 5 * 4
                nelements = ndata // ntotal

                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.is_debug_file:
                    self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    self.binary_debug.write('  #elementi = [eid_device, axial, axial_margin, torsion, torsion_margin]\n')
                    self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 5)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 5)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids

                    #[axial, torsion, SMa, SMt]
                    obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:]
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    struct1 = Struct(b(self._endian + 'i4f'))
                    for i in range(nelements):
                        edata = data[n:n+ntotal]
                        out = struct1.unpack(edata)
                        (eid_device, axial, axial_margin, torsion, torsion_margin) = out
                        eid = eid_device // 10
                        if self.is_debug_file:
                            self.binary_debug.write('  eid=%i; C=[%s]\n' % (eid, ', '.join(['%r' % di for di in out])))
                        obj.add_new_eid(dt, eid, axial, axial_margin, torsion, torsion_margin)
                        n += ntotal
            elif self.format_code in [2, 3] and self.num_wide == 5: # imag
                ntotal = 20
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_complex)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 5)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 5)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids

                    if is_magnitude_phase:
                        mag = floats[:, [1, 3]]
                        phase = floats[:, [2, 4]]
                        rtheta = radians(phase)
                        real_imag = mag * (cos(rtheta) + 1.j * sin(rtheta))
                    else:
                        real = floats[:, [1, 3]]
                        imag = floats[:, [2, 4]]
                        real_imag = real + 1.j * imag
                    obj.data[obj.itime, itotal:itotal2, :] = real_imag

                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    struct1 = Struct(b(self._endian + 'i4f'))
                    for i in range(nelements):
                        edata = data[n:n + ntotal]
                        out = struct1.unpack(edata)
                        (eid_device, axial_real, axial_imag, torsion_real, torsion_imag) = out
                        eid = eid_device // 10

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
            elif self.format_code == 1 and self.num_wide == 3: # random
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
            else:
                raise RuntimeError(self.code_information())

        elif self.element_type == 2: # CBEAM
            # 2-CBEAM
            ## TODO: fix method to follow correct pattern...regarding???

            if self.is_stress():
                result_name = 'cbeam_stress'
            else:
                result_name = 'cbeam_strain'

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            slot = getattr(self, result_name)

            if self.format_code == 1 and self.num_wide == 111:  # real
                # TODO: vectorize
                ntotal = 444 # 44 + 10*40  (11 nodes)

                if self.is_stress():
                    obj_vector_real = RealBeamStressArray
                else:
                    obj_vector_real = RealBeamStrainArray

                nelements = ndata // ntotal
                nlayers = nelements * 11
                auto_return, is_vectorized = self._create_oes_object4(
                    nlayers, result_name, slot, obj_vector_real)
                if auto_return:
                    self._data_factor = 11
                    return nelements * self.num_wide * 4
                obj = self.obj
                #s = self.struct_i
                nnodes = 10  # 11-1
                ntotal = self.num_wide * 4
                n1 = 44
                n2 = 40
                s1 = Struct(b(self._endian + 'ii9f'))
                s2 = Struct(b(self._endian + 'i9f'))
                nelements = ndata // ntotal
                for i in range(nelements):
                    edata = data[n:n+n1]
                    n += n1

                    out = s1.unpack(edata)
                    eid_device = out[0]
                    eid = eid_device // 10
                    if self.is_debug_file:
                        self.binary_debug.write('CBEAM-2 - eid=%i out=%s\n' % (eid, str(out)))

                    #(grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = out
                    obj.add_new_eid(dt, eid, out[1:])

                    for inode in range(nnodes):
                        edata = data[n:n+n2]
                        n += n2
                        out = s2.unpack(edata)
                        # (grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = out
                        obj.add(dt, eid, out)
            elif self.format_code in [2, 3] and self.num_wide == 111:  # imag and random?
                # TODO: vectorize
                if self.read_mode == 1:
                    return ndata
                if self.is_debug_file:
                    self.binary_debug.write('skipping imag/random OES-CBEAM\n')

                ntotal = 444 # 44 + 10*40  (11 nodes)
                #if self.is_stress():
                    #self.create_transient_object(self.beamStress, RealBeamStress)
                #else:
                    #self.create_transient_object(self.beamStrain, RealBeamStrain)

                nelements = ndata // ntotal
                #s = self.struct_i

                nnodes = 10  # 11-1
                ntotal = self.num_wide * 4
                n1 = 44
                n2 = 40
                s1 = Struct(b(self._endian + 'ii9f'))
                s2 = Struct(b(self._endian + 'i9f'))

                nelements = ndata // ntotal
                for i in range(nelements):
                    edata = data[n:n+n1]
                    n += n1

                    out = s1.unpack(edata)
                    eid_device = out[0]
                    eid = eid_device // 10
                    if self.is_debug_file:
                        self.binary_debug.write('CBEAM-2 - eid=%i out=%s\n' % (eid, str(out)))

                    #(grid, sd, ercr, exdr, exer, exfr,
                    #           exci, exdi, exei, exfi) = out
                    #obj.add_new_eid(dt, eid, out[1:])

                    for inode in range(nnodes):
                        edata = data[n:n+n2]
                        n += n2
                        out = s2.unpack(edata)
                        #obj.add(dt, eid, out)
                return ndata
            elif self.format_code == 1 and self.num_wide == 67: # random
                msg = self.code_information()
                if self.is_debug_file:
                    self.binary_debug.write('skipping OES-CBEAM\n')
                return self._not_implemented_or_skip(data, ndata, msg)
            else:
                raise RuntimeError(self.code_information())

        elif self.element_type == 4: # CSHEAR
            # 4-CSHEAR
            if self.is_stress():
                obj_vector_real = RealShearStressArray
                obj_vector_complex = ComplexShearStressArray
                result_name = 'cshear_stress'
            else:
                obj_vector_real = RealShearStrainArray
                obj_vector_complex = ComplexShearStrainArray
                result_name = 'cshear_strain'

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)

            slot = getattr(self, result_name)
            if self.format_code == 1 and self.num_wide == 4:  # real
                ntotal = 16  # 4*4
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                assert obj is not None
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 4)
                    itime = obj.itime
                    obj._times[itime] = dt
                    if itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 4)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids

                    #[max_strain, avg_strain, margin]
                    obj.data[itime, itotal:itotal2, :] = floats[:, 1:]
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    struct1 = Struct(b(self._endian + 'i3f'))
                    for i in range(nelements):
                        edata = data[n:n + ntotal]
                        out = struct1.unpack(edata)  # num_wide=5
                        if self.is_debug_file:
                            self.binary_debug.write('CSHEAR-4 - %s\n' % str(out))

                        (eid_device, max_strain, avg_strain, margin) = out
                        eid = eid_device // 10
                        obj.add_new_eid(dt, eid, max_strain, avg_strain, margin)
                        n += ntotal

            elif self.format_code in [2, 3] and self.num_wide == 5:  # imag
                ntotal = 20  # 4*5
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_complex)
                if auto_return:
                    return nelements * self.num_wide * 4

                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 5)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 5)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids

                    #(eid_device, etmaxr, etmaxi, etavgr, etavgi)
                    if is_magnitude_phase:
                        real = floats[:, [1, 3]]
                        imag = floats[:, [2, 4]]
                        real_imag = real + 1.j * imag
                    else:
                        mag = floats[:, [1, 3]]
                        phase = floats[:, [2, 4]]
                        rtheta = radians(phase)
                        real_imag = mag * (cos(rtheta) + 1.j * sin(rtheta))
                    obj.data[obj.itime, itotal:itotal2, :] = real_imag
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    struct1 = Struct(b(self._endian + 'i4f'))
                    for i in range(nelements):
                        edata = data[n:n + ntotal]
                        out = struct1.unpack(edata)  # num_wide=5
                        if self.is_debug_file:
                            self.binary_debug.write('CSHEAR-4 - %s\n' % str(out))
                        (eid_device, etmaxr, etmaxi, etavgr, etavgi) = out
                        eid = eid_device // 10

                        if is_magnitude_phase:
                            etmax = polar_to_real_imag(etmaxr, etmaxi)
                            etavg = polar_to_real_imag(etavgr, etavgi)
                        else:
                            etmax = complex(etmaxr, etmaxi)
                            etavg = complex(etavgr, etavgi)
                        obj.add_new_eid_sort1(dt, eid, (etmax, etavg))
                        n += ntotal
            elif self.format_code == 1 and self.num_wide == 3: # random
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
            else:
                raise RuntimeError(self.code_information())

        elif self.element_type in [11, 12, 13, 14]:  # springs
            # 11-CELAS1
            # 12-CELAS2
            # 13-CELAS3
            # 14-CELAS4
            if self.is_stress():
                obj_real = RealSpringStressArray
                obj_complex = ComplexSpringStressArray
                if self.element_type == 11:
                    result_name = 'celas1_stress'
                elif self.element_type == 12:
                    result_name = 'celas2_stress'
                elif self.element_type == 13:
                    result_name = 'celas3_stress'
                elif self.element_type == 14:
                    result_name = 'celas4_stress'
                else:
                    raise RuntimeError(self.element_type)
            else:
                obj_real = RealSpringStrainArray
                obj_complex = ComplexSpringStrainArray
                if self.element_type == 11:
                    result_name = 'celas1_strain'
                elif self.element_type == 12:
                    result_name = 'celas2_strain'
                elif self.element_type == 13:
                    result_name = 'celas3_strain'
                elif self.element_type == 14:
                    result_name = 'celas4_strain'
                else:
                    raise RuntimeError(self.element_type)

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            slot = getattr(self, result_name)

            if self.format_code == 1 and self.num_wide == 2:  # real
                ntotal = 8 # 2 * 4
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_real)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 2)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 2)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids

                    #(eid_device, stress)
                    obj.data[obj.itime, itotal:itotal2, 0] = floats[:, 1]
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    struct1 = Struct(b(self._endian + 'if'))
                    for i in range(nelements):
                        edata = data[n:n+ntotal]
                        out = struct1.unpack(edata)
                        (eid_device, ox) = out
                        eid = eid_device // 10
                        if self.is_debug_file:
                            self.binary_debug.write('  eid=%i result%i=[%i, %f]\n' % (eid, i, eid_device, ox))
                        obj.add_new_eid(dt, eid, ox)
                        n += ntotal
            elif self.format_code in [2, 3] and self.num_wide == 3:  # imag
                ntotal = 12
                nelements = ndata // ntotal

                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_complex)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                assert obj is not None, self.code_information()
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 3)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 3)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids

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
                    struct1 = Struct(b(self._endian + 'i2f'))
                    for i in range(nelements):
                        edata = data[n:n + ntotal]
                        out = struct1.unpack(edata)
                        (eid_device, axial_real, axial_imag) = out
                        eid = eid_device // 10

                        if is_magnitude_phase:
                            axial = polar_to_real_imag(axial_real, axial_imag)
                        else:
                            axial = complex(axial_real, axial_imag)

                        if self.is_debug_file:
                            self.binary_debug.write('  eid=%i result%i=[%i, %f, %f]\n' % (eid, i, eid_device, axial_real, axial_imag))
                        obj.add_sort1(dt, eid, axial)
                        n += ntotal
            elif self.format_code == 1 and self.num_wide == 3: # random
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
            else:
                raise RuntimeError(self.code_information())

        elif self.element_type == 34: # CBAR
            if self.is_stress():
                result_name = 'cbar_stress'
            else:
                result_name = 'cbar_strain'

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            slot = getattr(self, result_name)
            slot_vector = slot

            if self.format_code == 1 and self.num_wide == 16:  # real
                if self.is_stress():
                    obj_vector_real = RealBarStressArray
                else:
                    obj_vector_real = RealBarStrainArray

                ntotal = 16 * 4
                nelements = ndata // ntotal

                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    return ndata

                if self.is_debug_file:
                    self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    self.binary_debug.write('  #elementi = [eid_device, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,\n')
                    self.binary_debug.write('                           s1b, s2b, s3b, s4b, smaxb, sminb,        MSc]\n')
                    self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                obj = self.obj
                #is_vectorized = False
                if self.use_vector and is_vectorized:
                    # self.itime = 0
                    # self.ielement = 0
                    # self.itotal = 0
                    #self.ntimes = 0
                    #self.nelements = 0
                    n = nelements * self.num_wide * 4

                    ielement = obj.ielement
                    ielement2 = ielement + nelements
                    obj._times[obj.itime] = dt

                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 16)
                        eids = ints[:, 0] // 10
                        obj.element[ielement:ielement2] = eids

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 16)

                    #[s1a, s2a, s3a, s4a, axial, smaxa, smina, margin_tension,
                    # s1b, s2b, s3b, s4b,        smaxb, sminb, margin_compression]
                    obj.data[obj.itime, ielement:ielement2, :] = floats[:, 1:]
                    obj.itotal = ielement2
                    obj.ielement = ielement2
                else:
                    struct1 = Struct(b(self._endian + 'i15f'))
                    for i in range(nelements):
                        edata = data[n:n+ntotal]
                        out = struct1.unpack(edata)
                        (eid_device,
                         s1a, s2a, s3a, s4a, axial, smaxa, smina, margin_tension,
                         s1b, s2b, s3b, s4b, smaxb, sminb, margin_compression) = out
                        eid = eid_device // 10
                        if self.is_debug_file:
                            self.binary_debug.write('  eid=%i; C%i=[%s]\n' % (eid, i, ', '.join(['%r' % di for di in out])))
                        n += ntotal
                        obj.add_new_eid(dt, eid,
                                        s1a, s2a, s3a, s4a, axial, smaxa, smina, margin_tension,
                                        s1b, s2b, s3b, s4b, smaxb, sminb, margin_compression)
            elif self.format_code in [2, 3] and self.num_wide == 19:  # imag
                if self.is_stress():
                    obj_vector_complex = ComplexBarStressArray
                else:
                    obj_vector_complex = ComplexBarStrainArray

                ntotal = 76
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_complex)
                if auto_return:
                    return ndata

                if self.is_debug_file:
                    self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    self.binary_debug.write('  #elementi = [eid_device, s1a, s2a, s3a, s4a, axial,\n')
                    self.binary_debug.write('                           s1b, s2b, s3b, s4b]\n')
                    self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                obj = self.obj
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.itotal
                    itotal2 = itotal + nelements
                    ielement2 = itotal2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 19)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 19)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids

                    isave1 = [1, 2, 3, 4, 5, 11, 12, 13, 14]
                    isave2 = [6, 7, 8, 9, 10, 15, 16, 17, 18]
                    if is_magnitude_phase:
                        mag = floats[:, isave1]
                        phase = floats[:, isave2]
                        rtheta = radians(phase)
                        real_imag = mag * (cos(rtheta) + 1.j * sin(rtheta))
                    else:
                        real = floats[:, isave1]
                        imag = floats[:, isave2]
                        real_imag = real + 1.j * imag
                    obj.data[obj.itime, itotal:itotal2, :] = real_imag

                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    struct1 = Struct(b(self._endian + 'i18f'))
                    for i in range(nelements):
                        edata = data[n:n+ntotal]
                        n += ntotal
                        out = struct1.unpack(edata)
                        (eid_device,
                         s1ar, s2ar, s3ar, s4ar, axialr,
                         s1ai, s2ai, s3ai, s4ai, axiali,
                         s1br, s2br, s3br, s4br,
                         s1bi, s2bi, s3bi, s4bi) = out

                        eid = eid_device // 10
                        if self.is_debug_file:
                            self.binary_debug.write('  eid=%i; C%i=[%s]\n' % (eid, i, ', '.join(['%r' % di for di in out])))
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
            elif self.format_code == 1 and self.num_wide == 19: # random
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
            else:
                raise RuntimeError(self.code_information())

        elif self.element_type in [39, 67, 68]: # solid stress
            # 39-CTETRA
            # 67-CHEXA
            # 68-CPENTA
            if self.is_stress():
                obj_vector_real = RealSolidStressArray
                obj_vector_complex = ComplexSolidStressArray
                if self.element_type == 39: # CTETRA
                    nnodes_expected = 5  # 1 centroid + 4 corner points
                    result_name = 'ctetra_stress'
                    element_name = 'CTETRA4'
                elif self.element_type == 67:  # CHEXA
                    nnodes_expected = 9
                    result_name = 'chexa_stress'
                    element_name = 'CHEXA8'
                elif self.element_type == 68:  # CPENTA
                    nnodes_expected = 7
                    result_name = 'cpenta_stress'
                    element_name = 'CPENTA6'
                else:
                    raise RuntimeError(self.code_information())
                    #msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
                    #return self._not_implemented_or_skip(data, ndata, msg)
            else:
                obj_vector_real = RealSolidStrainArray
                obj_vector_complex = ComplexSolidStrainArray

                if self.element_type == 39: # CTETRA
                    nnodes_expected = 5  # 1 centroid + 4 corner points
                    result_name = 'ctetra_strain'
                    element_name = 'CTETRA4'
                elif self.element_type == 67:  # CHEXA
                    nnodes_expected = 9
                    result_name = 'chexa_strain'
                    element_name = 'CHEXA8'
                elif self.element_type == 68:  # CPENTA
                    nnodes_expected = 7
                    result_name = 'cpenta_strain'
                    element_name = 'CPENTA6'
                else:
                    raise RuntimeError(self.code_information())
                    #msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
                    return self._not_implemented_or_skip(data, ndata, msg)

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            slot = getattr(self, result_name)


            numwide_real = 4 + 21 * nnodes_expected
            numwide_imag = 4 + (17 - 4) * nnodes_expected
            numwide_random = 4 + (11 - 4) * nnodes_expected
            preline1 = '%s-%s' % (self.element_name, self.element_type)
            preline2 = ' ' * len(preline1)

            self._data_factor = nnodes_expected
            if self.format_code == 1 and self.num_wide == numwide_real:  # real
                ntotal = 16 + 84 * nnodes_expected
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    itotali = obj.itotal + nelements
                    itotal2 = obj.itotal + nelements * nnodes_expected
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        # (eid_device, cid, abcd, nnodes)
                        ints = fromstring(data, dtype=self.idtype)
                        try:
                            ints1 = ints.reshape(nelements, numwide_real)
                        except ValueError:
                            msg = 'ints.shape=%s; size=%s ' % (str(ints.shape), ints.size)
                            msg += 'nelements=%s numwide_real=%s nelements*numwide=%s' % (nelements, numwide_real,
                                                                                          nelements * numwide_real)
                            raise ValueError(msg)
                        eids = ints1[:, 0] // 10
                        cids = ints1[:, 1]
                        #nids = ints1[:, 4]
                        assert eids.min() > 0, eids.min()
                        obj.element_node[itotal:itotal2, 0] = repeat(eids, nnodes_expected)
                        ints2 = ints1[:, 4:].reshape(nelements * nnodes_expected, 21)
                        grid_device = ints2[:, 0]#.reshape(nelements, nnodes_expected)

                        #print('%s-grid_device=%s' % (self.element_name, grid_device))
                        grid_device2 = repeat(grid_device, nnodes_expected)
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

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, numwide_real)[:, 4:]
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
                    struct1 = Struct(b(self._endian + 'ii4si'))
                    struct2 = Struct(b(self._endian + 'i20f'))
                    if self.is_debug_file:
                        msg = '%s-%s nelements=%s nnodes=%s; C=[sxx, sxy, s1, a1, a2, a3, pressure, svm,\n' % (
                            self.element_name, self.element_type, nelements, nnodes_expected)
                        msg += '                                 syy, syz, s2, b1, b2, b3,\n'
                        msg += '                                 szz, sxz, s3, c1, c2, c3]\n'
                        self.binary_debug.write(msg)

                    for i in range(nelements):
                        edata = data[n:n+16]
                        out = struct1.unpack(edata)
                        (eid_device, cid, abcd, nnodes) = out
                        eid = eid_device // 10

                        if self.is_debug_file:
                            self.binary_debug.write('%s - eid=%i; %s\n' % (preline1, eid, str(out)))

                        assert nnodes < 21, 'print_block(data[n:n+16])'  #self.print_block(data[n:n+16])

                        n += 16
                        for inode in range(nnodes_expected):  # nodes pts, +1 for centroid (???)
                            out = struct2.unpack(data[n:n + 84]) # 4*21 = 84
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
                                #element_name = self.element_name + str(nnodes) #  this is correct, but fails
                                obj.add_eid(element_name, cid, dt, eid, grid,
                                            sxx, syy, szz, sxy, syz, sxz, s1, s2, s3,
                                            a_cos, b_cos, c_cos, pressure, svm)
                            else:
                                obj.add_node(dt, eid, inode, grid,
                                             sxx, syy, szz, sxy, syz, sxz, s1, s2, s3,
                                             a_cos, b_cos, c_cos, pressure, svm)
                            n += 84

            elif self.format_code in [2, 3] and self.num_wide == numwide_imag:  # complex
                # TODO: vectorize
                ntotal = numwide_imag * 4
                nelements = ndata // ntotal
                self.ntotal += nelements * nnodes_expected
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_complex)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj

                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    ielement = obj.ielement
                    ielement2 = ielement + nelements
                    itotal = obj.itotal
                    itotal2 = itotal + nelements * nnodes_expected

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, numwide_imag)
                    floats1 = floats[:, 4:].reshape(nelements * nnodes_expected, 13)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, numwide_imag)
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
                    if is_magnitude_phase:
                        mag = floats1[:, isave1]
                        phase = floats1[:, isave2]
                        rtheta = radians(phase)
                        real_imag = mag * (cos(rtheta) + 1.j * sin(rtheta))
                    else:
                        real = floats1[:, isave1]
                        imag = floats1[:, isave2]
                        real_imag = real + 1.j * imag
                    obj.data[obj.itime, itotal:itotal2, :] = real_imag

                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s1 = Struct(b(self._endian + '2i4si'))
                    s2 = Struct(b(self._endian + 'i12f'))
                    for i in range(nelements):
                        edata = data[n:n+16]
                        n += 16
                        out = s1.unpack(edata)
                        (eid_device, cid, ctype, nodef) = out
                        eid = eid_device // 10
                        if self.is_debug_file:
                            self.binary_debug.write('  eid=%i C=[%s]\n' % (eid, ', '.join(['%r' % di for di in out])))

                        #element_name = self.element_name + str(nodef)  # this is correct, but has problems...
                        obj.add_eid_sort1(self.element_type, element_name, dt, eid, cid, ctype, nodef)
                        for inode in range(nnodes_expected):
                            edata = data[n:n+52]
                            n += 52
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
                                self.binary_debug.write('       node%s=[%s]\n' % (grid, ', '.join(['%r' % di for di in out])))
                            obj.add_node_sort1(dt, eid, grid, inode,
                                               ex, ey, ez, etxy, etyz, etzx)
            elif self.format_code == 1 and self.num_wide == numwide_random: # random
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
            elif self.format_code == 2 and self.num_wide == 130:
                # 130 - CHEXA-67
                #msg = 'OES-CHEXA-random-numwide=%s numwide_real=%s numwide_imag=%s numwide_random=%s' % (
                    #self.num_wide, numwide_real, numwide_imag, numwide_random)
                #return self._not_implemented_or_skip(data, ndata, msg)

                num_wide_random = 4 + nnodes_expected * (17 - 4)
                #if self.num_wide ==
                if self.read_mode == 1:
                    return ndata
                return ndata
                print('numwide=%s numwide_random=%s attempt2=%s subcase=%s' % (self.num_wide, numwide_random, num_wide_random, self.isubcase))
                #msg = self.code_information()
                ntotal = 130
                nelements = ndata // ntotal

                # cid, coord_type, nactive_pnts,
                #      nid, oxx, oyy, ozz, txy, tyz, txz
                struct1 = Struct(b('2i 4s'))
                struct2 = Struct(b('i6f'))
                for i in range(nelements):
                    edata = data[n:n+12]
                    out = struct1.unpack(edata)
                    (eid_device, cid, abcd) = out
                    eid = eid_device // 10
                    if self.is_debug_file:
                        self.binary_debug.write('%s - eid=%i; %s\n' % (preline1, eid, str(out)))
                    n += 12
                    for inode in range(nnodes_expected):  # nodes pts, +1 for centroid (???)
                        out = struct2.unpack(data[n:n + 28]) # 4*7 = 28
                        if self.is_debug_file:
                            self.binary_debug.write('%s - %s\n' % (preline2, str(out)))
                        (grid_device, sxx, syy, sz, txy, tyz, txz) = out
                #msg = 'OES-CHEXA-random-numwide=%s numwide_real=%s numwide_imag=%s numwide_random=%s' % (
                    #self.num_wide, numwide_real, numwide_imag, numwide_random)
                #return self._not_implemented_or_skip(data, ndata, msg)
            else:
                raise RuntimeError(self.code_information())

        #=========================
        # plates
        elif self.element_type == 33:  # CQUAD4-centroidal
            # 33-QUAD4-centroidal
            if self.is_stress():
                obj_vector_real = RealPlateStressArray
                obj_vector_complex = ComplexPlateStressArray
                result_name = 'cquad4_stress'
            else:
                obj_vector_real = RealPlateStrainArray
                obj_vector_complex = ComplexPlateStrainArray
                result_name = 'cquad4_strain'

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            slot = getattr(self, result_name)

            numwide_real = 17
            if self.format_code == 1 and self.num_wide == 17:  # real
                ntotal = 68  # 4*17
                nelements = ndata // ntotal
                nlayers = nelements * 2
                nnodes_expected = 2

                auto_return, is_vectorized = self._create_oes_object4(
                    nlayers, result_name, slot, obj_vector_real)
                if auto_return:
                    self._data_factor = 2
                    return nelements * ntotal

                obj = self.obj
                #print('dt=%s, itime=%s' % (obj.itime, dt))
                #is_vectorized = False
                assert obj.is_built == True, obj.is_built
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    ielement = obj.ielement
                    ielement2 = ielement + nelements
                    itotal = obj.itotal
                    itotal2 = itotal + nelements * nnodes_expected
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype)
                        ints1 = ints.reshape(nelements, numwide_real)
                        eids = ints1[:, 0] // 10
                        eids = np.vstack([eids, eids]).T.ravel()
                        assert eids.min() > 0, eids.min()
                        obj.element_node[itotal:itotal2, 0] = eids

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, numwide_real)[:, 1:]

                    #fd, sx, sy, txy, angle, major, minor, max_shear
                    floats1 = floats.reshape(nelements * nnodes_expected, 8)
                    obj.data[obj.itime, itotal:itotal2, :] = floats1
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    struct1 = Struct(b(self._endian + 'i16f'))
                    cen = 0 # CEN/4
                    for i in range(nelements):
                        edata = data[n:n+ntotal]
                        out = struct1.unpack(edata)

                        (eid_device,
                         fd1, sx1, sy1, txy1, angle1, major1, minor1, max_shear1,
                         fd2, sx2, sy2, txy2, angle2, major2, minor2, max_shear2) = out

                        eid = eid_device // 10
                        if self.is_debug_file:
                            self.binary_debug.write('  eid=%i C=[%s]\n' % (eid, ', '.join(['%r' % di for di in out])))

                        obj._add_new_eid(dt, eid, cen, fd1, sx1, sy1,
                                         txy1, angle1, major1, minor1, max_shear1)
                        obj._add(dt, eid, cen, fd2, sx2, sy2, txy2,
                                 angle2, major2, minor2, max_shear2)
                        n += ntotal
            elif self.format_code in [2, 3] and self.num_wide == 15:  # imag
                # TODO: vectorize
                nnodes = 0  # centroid + 4 corner points
                ntotal = 4 * (15 * (nnodes + 1))
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_complex)
                if auto_return:
                    self._data_factor = 2
                    return nelements * ntotal

                obj = self.obj
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    nnodes_all = (nnodes + 1)
                    itotal = obj.itotal
                    itotal2 = itotal + 2 * nelements * nnodes_all
                    ielement = obj.ielement
                    ielement2 = ielement + nelements

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 15 * nnodes_all)
                    floats1 = floats[:, 1:].reshape(nelements * nnodes_all * 2, 7)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 15 * nnodes_all)
                        eids = ints[:, 0] // 10
                        ints[:, 0] = 0
                        ints1 = ints.reshape(nelements * nnodes_all, 15)
                        nids = ints[:, 0]
                        #print(eids)
                        assert eids.min() > 0, eids.min()
                        eids2 = np.vstack([eids, eids]).T.ravel()
                        nids2 = np.vstack([nids, nids]).T.ravel()
                        #print(eids2, itotal, itotal2)
                        obj.element_node[itotal:itotal2, 0] = eids2
                        obj.element_node[itotal:itotal2, 1] = nids2

                    #[fd, sxr, sxi, syr, syi, txyr, txyi]
                    isave1 = [1, 3, 5]
                    isave2 = [2, 4, 6]
                    if is_magnitude_phase:
                        mag = floats1[:, isave1]
                        phase = floats1[:, isave2]
                        rtheta = radians(phase)
                        real_imag = mag * (cos(rtheta) + 1.j * sin(rtheta))
                    else:
                        real = floats1[:, isave1]
                        imag = floats1[:, isave2]
                        real_imag = real + 1.j * imag

                    obj.fiber_curvature[itotal:itotal2] = floats1[:, 0]
                    obj.data[obj.itime, itotal:itotal2, :] = real_imag
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s1 = Struct(b(self._endian + 'i14f'))
                    s2 = Struct(b(self._endian + 'i14f'))

                    cen = 0 # 'CEN/4'
                    for i in range(nelements):
                        edata = data[n:n+60]  # 4*15=60
                        n += 60
                        out = s1.unpack(edata)  # 15
                        (eid_device,
                         fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
                         fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i) = out

                        eid = eid_device // 10
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

                        obj.add_new_eid_sort1(dt, eid, cen, fd1, sx1, sy1, txy1)
                        obj.add_sort1(dt, eid, cen, fd2, sx2, sy2, txy2)

                        for node_id in range(nnodes):  # nodes pts
                            edata = data[n:n+60]  # 4*15=60
                            n += 60
                            out = s2.unpack(edata)
                            if self.is_debug_file:
                                self.binary_debug.write('  %s\n' % str(out))
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
                            obj.add_new_node_sort1(dt, eid, grid, fd1, sx1, sy1, txy1)
                            obj.add_sort1(dt, eid, grid, fd2, sx2, sy2, txy2)
            elif self.format_code == 1 and self.num_wide == 0: # random
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
            else:
                raise RuntimeError(self.code_information())

        elif self.element_type == 74:  # TRIA3
            if self.is_stress():
                obj_vector_real = RealPlateStressArray
                obj_vector_complex = ComplexPlateStressArray
                result_name = 'ctria3_stress'
            else:
                obj_vector_real = RealPlateStrainArray
                obj_vector_complex = ComplexPlateStrainArray
                result_name = 'ctria3_strain'

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)

            slot = getattr(self, result_name)
            if self.format_code == 1 and self.num_wide == 17:  # real
                ntotal = 68  # 4*17
                nelements = ndata // ntotal
                nlayers = nelements * 2  # 2 layers per node
                auto_return, is_vectorized = self._create_oes_object4(
                    nlayers, result_name, slot, obj_vector_real)
                if auto_return:
                    self._data_factor = 2
                    return nelements * ntotal

                if self.is_debug_file:
                    self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    self.binary_debug.write('  #elementi = [eid_device, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
                    self.binary_debug.write('  #                        fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,]\n')
                    self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                obj = self.obj
                if self.use_vector and is_vectorized:
                    nfields = 17 * nelements
                    nbytes = nfields * 4
                    itotal = obj.itotal
                    iend = obj.itotal + nlayers

                    itime = obj.itime
                    if itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 17)
                        eids = ints[:, 0] // 10
                        ilayers = ints[:, 1]
                        ints2 = ints[:, 1:].reshape(nlayers, 8)
                        assert eids.min() > 0, eids
                        obj._times[obj.itime] = dt
                        obj.element_node[itotal:iend:2, 0] = eids
                        obj.element_node[itotal+1:iend+1:2, 0] = eids
                        #obj.element_node[itotal:iend, 1] = 0
                        #print('obj.element_node\n', obj.element_node)

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 17)
                    floats1 = floats[:, 1:].reshape(nlayers, 8)
                    obj.data[obj.itime, itotal:iend, :] = floats1
                    obj._times[obj.itime] = dt
                    obj.itotal += nlayers
                    n = nbytes
                else:
                    cen = 0 # 'CEN/3'
                    struct1 = Struct(b(self._endian + 'i16f'))
                    for i in range(nelements):
                        edata = data[n:n + ntotal]
                        out = struct1.unpack(edata)
                        (eid_device,
                         fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                         fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out
                        eid = eid_device // 10
                        if self.is_debug_file:
                            self.binary_debug.write('  OES CTRIA3-74 - eid=%i; C=[%s]\n' % (eid, ', '.join(['%r' % di for di in out])))

                        obj._add_new_eid(dt, eid, cen, fd1, sx1, sy1,
                                         txy1, angle1, major1, minor1, vm1)
                        obj._add(dt, eid, cen, fd2, sx2, sy2, txy2,
                                 angle2, major2, minor2, vm2)
                        n += ntotal
                #ass
            elif self.format_code in [2, 3] and self.num_wide == 15:  # imag
                ntotal = 60  # 4*15
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_complex)
                if auto_return:
                    self._data_factor = 2
                    return nelements * self.num_wide * 4
                obj = self.obj

                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.itotal
                    itotal2 = itotal + nelements * 2
                    ielement = obj.ielement
                    ielement2 = ielement + nelements

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 15)
                    floats1 = floats[:, 1:].reshape(nelements * 2, 7)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 15)
                        eids = ints[:, 0] // 10
                        ints[:, 0] = 0
                        ints1 = ints.reshape(nelements, 15)
                        nids = ints[:, 0]

                        assert eids.min() > 0, eids.min()
                        eids2 = np.vstack([eids, eids]).T.ravel()
                        nids2 = np.vstack([nids, nids]).T.ravel()
                        obj.element_node[itotal:itotal2, 0] = eids2
                        obj.element_node[itotal:itotal2, 1] = nids2

                    #[fd, sxr, sxi, syr, syi, txyr, txyi]
                    isave1 = [1, 3, 5]
                    isave2 = [2, 4, 6]
                    if is_magnitude_phase:
                        mag = floats1[:, isave1]
                        phase = floats1[:, isave2]
                        rtheta = radians(phase)
                        real_imag = mag * (cos(rtheta) + 1.j * sin(rtheta))
                    else:
                        real = floats1[:, isave1]
                        imag = floats1[:, isave2]
                        real_imag = real + 1.j * imag

                    obj.fiber_curvature[itotal:itotal2] = floats1[:, 0]
                    obj.data[obj.itime, itotal:itotal2, :] = real_imag
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    struct1 = Struct(b(self._endian + 'i14f'))
                    cen = 0 # CEN/3
                    for i in range(nelements):
                        edata = data[n:n + ntotal]
                        out = struct1.unpack(edata)
                        (eid_device,
                         fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
                         fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i,) = out
                        eid = eid_device // 10

                        if self.is_debug_file:
                            self.binary_debug.write('  OESC CTRIA3-74 - eid=%i; C=[%s]\n' % (eid, ', '.join(['%r' % di for di in out])))

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
                        obj.add_new_eid_sort1(dt, eid, cen, fd1, sx1, sy1, txy1)
                        obj.add_sort1(dt, eid, cen, fd2, sx2, sy2, txy2)
                        n += ntotal
            elif self.format_code == 1 and self.num_wide == 9: # random
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
            else:
                raise RuntimeError(self.code_information())

        elif self.element_type in [64, 70, 75, 82, 144]:  # bilinear plates
            # 64-CQUAD8
            # 70-CTRIAR
            # 75-CTRIA6
            # 82-CQUADR
            # 144-CQUAD4-bilinear
            if self.is_stress():
                obj_vector_real = RealPlateStressArray
                obj_vector_complex = ComplexPlateStressArray
                if self.element_type == 64: # CQUAD8
                    result_name = 'cquad8_stress'
                    #gridC = 'CEN/8'
                elif self.element_type == 70:  # CTRIAR
                    result_name = 'ctriar_stress'
                    #gridC = 'CEN/3'
                elif self.element_type == 75:  # CTRIA6
                    result_name = 'ctria6_stress'
                    #gridC = 'CEN/6'
                elif self.element_type == 82:  # CQUADR
                    result_name = 'cquadr_stress'
                    #gridC = 'CEN/4'
                elif self.element_type == 144:  # CQUAD4-bilinear
                    # there's no nead to separate this with centroidal strain
                    # because you can only have one in a given OP2
                    result_name = 'cquad4_stress'
                    #gridC = 'CEN/4'
                else:
                    raise RuntimeError(self.code_information())
                    #msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
                    #return self._not_implemented_or_skip(data, ndata, msg)
            else:
                obj_vector_real = RealPlateStrainArray
                obj_vector_complex = ComplexPlateStrainArray

                if self.element_type == 64: # CQUAD8
                    result_name = 'cquad8_strain'
                    #gridC = 'CEN/8'
                elif self.element_type == 70:  # CTRIAR
                    result_name = 'ctriar_strain'
                    #gridC = 'CEN/3'
                elif self.element_type == 75:  # CTRIA6
                    result_name = 'ctria6_strain'
                    #gridC = 'CEN/6'
                elif self.element_type == 82: # CQUADR
                    result_name = 'cquadr_strain'
                    #gridC = 'CEN/4'
                elif self.element_type == 144: # CQUAD4-bilinear
                    # there's no nead to separate this with centroidal strain
                    # because you can only have one in a given OP2
                    result_name = 'cquad4_strain'
                    #gridC = 'CEN/4'
                else:
                    raise RuntimeError(self.code_information())

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)

            if self.element_type in [64, 82, 144]:
                nnodes = 4 # + 1 centroid
            elif self.element_type in [70, 75]:
                nnodes = 3 # + 1 centroid
            else:
                raise RuntimeError(self.code_information())
            nnodes_all = nnodes + 1 # adding the centroid

            slot = getattr(self, result_name)
            numwide_real = 2 + 17 * nnodes_all
            numwide_imag = 2 + 15 * nnodes_all
            numwide_random = 2 + 9 * nnodes_all

            etype = self.element_name
            #grid_center = 'CEN/%i' % nnodes
            if self.format_code == 1 and self.num_wide == numwide_real:  # real
                ntotal = 4 * (2 + 17 * nnodes_all)
                nelements = ndata // ntotal
                nlayers = 2 * nelements * nnodes_all  # 2 layers per node

                auto_return, is_vectorized = self._create_oes_object4(
                    nlayers, result_name, slot, obj_vector_real)
                if auto_return:
                    self._data_factor = 10  # TODO: why is this 10?
                    return nelements * self.num_wide * 4

                obj = self.obj
                #print('dt=%s, itime=%s' % (obj.itime, dt))
                if self.use_vector and is_vectorized:
                    # self.itime = 0
                    # self.ielement = 0
                    # self.itotal = 0
                    #self.ntimes = 0
                    #self.nelements = 0
                    n = nelements * self.num_wide * 4

                    istart = obj.itotal
                    iend = istart + nlayers
                    obj._times[obj.itime] = dt

                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, numwide_real)
                        ints1 = ints[:, 2:].reshape(nlayers//2, 17)[:, 0].reshape(nelements, nnodes_all)
                        ints1[:, 0] = 0.
                        nids = ints1.ravel()

                        eids = ints[:, 0] // 10
                        eids2 = array([eids] * (nnodes_all * 2), dtype='int32').T.ravel()
                        #nids2 = array([nids, nids], dtype='int32').T.ravel()
                        #eids2 = vstack([eids] * (nnodes_all * 2)).T.ravel()
                        nids2 = vstack([nids, nids]).T.ravel()
                        #eids2 = repeat(eids, 2 * nnodes_all)
                        # nids2 = repeat(nids, 2)
                        obj.element_node[istart:iend, 0] = eids2
                        obj.element_node[istart:iend, 1] = nids2
                        #assert obj.element_node[:iend, 0].min() > 0, eids2
                        if obj.nonlinear_factor is not None:
                            float_mask = np.arange(nelements * numwide_real, dtype=np.int32).reshape(nelements, numwide_real)
                            float_mask1 = float_mask[:, 2:].reshape(nlayers // 2, 17)[:, 1:].reshape(nlayers, 8)
                            obj.float_mask = float_mask1

                    if obj.nonlinear_factor is not None:
                        results = fromstring(data, dtype=self.fdtype)[obj.float_mask]
                    else:
                        floats = fromstring(data, dtype=self.fdtype).reshape(nelements, numwide_real)
                        floats1 = floats[:, 2:].reshape(nlayers // 2, 17)
                        results = floats1[:, 1:].reshape(nlayers, 8)

                    #[fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm]
                    obj.data[obj.itime, istart:iend, :] = results
                else:
                    n = 0
                    center_format = b'i4si16f'
                    node_format = b'i16f'
                    cs = Struct(center_format)
                    ns = Struct(node_format)

                    if self.is_debug_file:
                        self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                        self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                        self.binary_debug.write('  #elementi = [centeri, node1i, node2i, node3i, node4i]\n')
                        self.binary_debug.write('  #centeri = [eid_device, j, grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
                        self.binary_debug.write('  #                                fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,)]\n')
                        self.binary_debug.write('  #nodeji = [grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
                        self.binary_debug.write('  #                fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,)]\n')
                        self.binary_debug.write('  nelements=%i; nnodes=%i # +1 centroid\n' % (nelements, nnodes))

                    grid_center = 0
                    for i in range(nelements):
                        edata = data[n:n+76]

                        out = cs.unpack(edata)  # len=17*4
                        # j is the number of nodes, so CQUAD4 -> 4, but we don't need to save it...
                        (eid_device, j,
                         grid,
                         fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                         fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out
                        eid = eid_device // 10
                        #print(out[:3])
                        if self.is_debug_file:
                            self.binary_debug.write('  eid=%i; C=[%s]\n' % (eid, ', '.join(['%r' % di for di in out])))

                        obj._add_new_eid(dt, eid, grid_center, fd1, sx1, sy1,
                                         txy1, angle1, major1, minor1, vm1)
                        obj._add(dt, eid, grid_center, fd2, sx2, sy2, txy2,
                                 angle2, major2, minor2, vm2)
                        n += 76
                        for inode in range(nnodes):
                            out = ns.unpack(data[n:n + 68])
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

                            obj._add_new_node(dt, eid, grid, fd1, sx1, sy1,
                                              txy1, angle1, major1, minor1, vm1)
                            obj._add(dt, eid, grid, fd2, sx2, sy2,
                                     txy2, angle2, major2, minor2, vm2)
                            n += 68
            elif self.format_code in [2, 3] and self.num_wide == numwide_imag:  # imag
                ntotal = numwide_imag * 4
                assert self.num_wide * 4 == ntotal, 'numwide*4=%s ntotal=%s' % (self.num_wide*4, ntotal)
                nelements = ndata // ntotal

                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_complex)
                if auto_return:
                    self._data_factor = 2
                    return nelements * ntotal
                obj = self.obj

                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.itotal
                    itotal2 = itotal + nelements * (nnodes_all * 2)
                    ielement = obj.ielement
                    ielement2 = ielement + nelements

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, numwide_imag)
                    floats1 = floats[:, 2:].reshape(nelements * nnodes_all, 15)
                    floats2 = floats1[:, 1:].reshape(nelements * nnodes_all * 2, 7)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, numwide_imag)
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
                    if is_magnitude_phase:
                        mag = floats2[:, isave1]
                        phase = floats2[:, isave2]
                        rtheta = radians(phase)
                        real_imag = mag * (cos(rtheta) + 1.j * sin(rtheta))
                    else:
                        real = floats2[:, isave1]
                        imag = floats2[:, isave2]
                        real_imag = real + 1.j * imag

                    obj.fiber_curvature[itotal:itotal2] = floats2[:, 0]
                    obj.data[obj.itime, itotal:itotal2, :] = real_imag
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    grid_center = 0
                    s1 = self.struct_2i  # 2
                    s2 = Struct(b(self._endian + 'i14f')) # 15
                    for i in range(nelements):
                        (eid_device, _) = s1.unpack(data[n:n+8])
                        n += 8

                        edata = data[n:n+60]  # 4*15
                        n += 60
                        out = s2.unpack(edata)  # len=15*4

                        eid = self._check_id(eid_device, 'element', stress_name, out)
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

                        obj.add_new_eid_sort1(dt, eid, grid_center, fd1, sx1, sy1, txy1)
                        obj.add_sort1(dt, eid, grid_center, fd2, sx2, sy2, txy2)

                        for node_id in range(nnodes):  # nodes pts
                            edata = data[n:n + 60]  # 4*15=60
                            n += 60
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

                            obj.add_new_node_sort1(dt, eid, grid, fd1, sx1, sy1, txy1)
                            obj.add_sort1(dt, eid, grid, fd2, sx2, sy2, txy2)
            elif self.format_code == 1 and self.num_wide == numwide_random: # random
                msg = self.code_information()
                msg += '  numwide=%s numwide_real=%s numwide_imag=%s numwide_random=%s' % (
                    self.num_wide, numwide_real, numwide_imag, numwide_random)
                return self._not_implemented_or_skip(data, ndata, msg)
            elif self.format_code == 2 and self.num_wide == 87:
                # 87 - CQUAD4-144
                #msg = self.code_information()
                msg = 'OES-CQUAD4-random-numwide=%s numwide_real=%s numwide_imag=%s numwide_random=%s' % (
                    self.num_wide, numwide_real, numwide_imag, numwide_random)
                return self._not_implemented_or_skip(data, ndata, msg)
            else:
                raise RuntimeError(self.code_information())

        elif self.element_type in [88, 90]: # nonlinear shells
            # 88-CTRIA3NL
            # 90-CQUAD4NL
            if self.is_stress():
                if self.element_type == 88:
                    result_name = 'nonlinear_ctria3_stress'
                elif self.element_type == 90:
                    result_name = 'nonlinear_cquad4_stress'
                else:
                    raise RuntimeError(self.element_type)
            else:
                if self.element_type == 88:
                    result_name = 'nonlinear_ctria3_strain'
                elif self.element_type == 90:
                    result_name = 'nonlinear_cquad4_strain'
                else:
                    raise RuntimeError(self.element_type)

            slot = getattr(self, result_name)
            self._results._found_result(result_name)

            if self.format_code == 1 and self.num_wide == 13 and self.element_type in [88, 90]:  # real
                # TODO: vectorize
                # single layered hyperelastic (???) ctria3, cquad4
                ntotal = 52  # 4*13
                nelements = ndata // ntotal

                obj_vector_real = RealNonlinearPlateArray
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.use_vector and is_vectorized:
                    n = nelements * self.num_wide * 4

                    ielement = obj.ielement
                    ielement2 = ielement + nelements
                    obj._times[obj.itime] = dt

                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 13)
                        eids = ints[:, 0] // 10
                        obj.element_node[ielement:ielement2, 0] = eids

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 13)

                    #[fiber_distance, oxx, oyy, ozz, txy, exx, eyy, ezz, exy, es, eps, ecs]
                    #print(ints)
                    floats[:, 1] = 0
                    obj.data[obj.itime, ielement:ielement2, :] = floats[:, 1:]
                    obj.ielement = ielement2
                    obj.itotal = ielement2
                else:
                    struct1 = Struct(b(self._endian + 'i12f'))  # 1+12=13
                    for i in range(nelements):
                        edata = data[n:n + ntotal]
                        out = struct1.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('CQUADNL-90 - %s\n' % str(out))

                        (eid_device, fd1,
                         sx1, sy1, sz1, txy1, es1, eps1, ecs1,
                         ex1, ey1, ez1, exy1) = out
                        eid = eid_device // 10
                        obj.add_new_eid(
                            dt, eid, self.element_type, fd1,
                            sx1, sy1, sz1, txy1, es1, eps1, ecs1,
                            ex1, ey1, ez1, exy1)
                    n += ntotal
            elif self.format_code == 1 and self.num_wide == 25 and self.element_type in [88, 90]:
                # TODO: vectorize
                ntotal = 100  # 4*25
                nelements = ndata // ntotal
                obj_vector_real = RealNonlinearPlateArray
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    self._data_factor = 2
                    return nelements * self.num_wide * 4

                obj = self.obj
                is_vectorized = False
                if self.use_vector and is_vectorized:
                    n = nelements * self.num_wide * 4

                    ielement = obj.ielement
                    ielement2 = ielement + nelements
                    obj._times[obj.itime] = dt

                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 13)
                        eids = ints[:, 0] // 10
                        obj.element_node[ielement:ielement2, 0] = eids

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 13)

                    #[fiber_distance, oxx, oyy, ozz, txy, exx, eyy, ezz, exy, es, eps, ecs]
                    #print(ints)
                    floats[:, 1] = 0
                    obj.data[obj.itime, ielement:ielement2, :] = floats[:, 1:]
                    obj.ielement = ielement2
                    obj.itotal = ielement2
                else:
                    etype = self.element_type
                    struct1 = Struct(b(self._endian + 'i24f')) # 1+24=25
                    for i in range(nelements):
                        edata = data[n:n + ntotal]
                        out = struct1.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('CQUADNL-90 - %s\n' % str(out))
                        (eid_device,
                         fd1, sx1, sy1, undef1, txy1, es1, eps1, ecs1, ex1, ey1, undef2, etxy1,
                         fd2, sx2, sy2, undef3, txy2, es2, eps2, ecs2, ex2, ey2, undef4, etxy2) = out
                        eid = eid_device // 10
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
            else:
                raise RuntimeError(self.code_information())

        elif self.element_type in [95, 96, 97, 98]: # composite shell
            # 95 - CQUAD4
            # 96 - CQUAD8
            # 97 - CTRIA3
            # 98 - CTRIA6 (composite)
            if self.is_stress():
                obj_vector_real = RealCompositePlateStressArray
                #obj_vector_complex = ComplexCompositePlateStressArray
                if self.element_type == 95: # CQUAD4
                    result_name = 'cquad4_composite_stress'
                elif self.element_type == 96:  # CQUAD8
                    result_name = 'cquad8_composite_stress'
                elif self.element_type == 97:  # CTRIA3
                    result_name = 'ctria3_composite_stress'
                elif self.element_type == 98:  # CTRIA6
                    result_name = 'ctria6_composite_stress'
                #elif self.element_type == ???:  # CTRIA6
                    #result_name = 'ctriar_composite_stress'
                #elif self.element_type == 10:  # CTRIA6
                else:
                    raise RuntimeError(self.code_information())
            else:
                obj_vector_real = RealCompositePlateStrainArray
                #obj_vector_complex = ComplexCompositePlateStrainArray
                if self.element_type == 95: # CQUAD4
                    result_name = 'cquad4_composite_strain'
                elif self.element_type == 96:  # CQUAD8
                    result_name = 'cquad8_composite_strain'
                elif self.element_type == 97:  # CTRIA3
                    result_name = 'ctria3_composite_strain'
                elif self.element_type == 98:  # CTRIA6
                    result_name = 'ctria6_composite_strain'
                #elif self.element_type == ???:  # CTRIA6
                    #result_name = 'ctriar_composite_strain'
                else:
                    raise RuntimeError(self.code_information())

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            slot = getattr(self, result_name)

            etype = self.element_name
            if self.format_code == 1 and self.num_wide == 11:  # real
                ntotal = 44
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.is_debug_file:
                    self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    self.binary_debug.write('  element1 = [eid_device, layer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)]\n')
                    self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                if self.use_vector and is_vectorized:
                    n = nelements * self.num_wide * 4

                    istart = obj.itotal
                    iend = istart + nelements
                    obj._times[obj.itime] = dt

                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 11)
                        eids = ints[:, 0] // 10
                        nids = ints[:, 1]
                        obj.element_layer[istart:iend, 0] = eids
                        obj.element_layer[istart:iend, 1] = nids

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 11)
                    #[o1, o2, t12, t1z, t2z, angle, major, minor, ovm]
                    obj.data[obj.itime, istart:iend, :] = floats[:, 2:]
                else:
                    struct1 = Struct(b(self._endian + 'ii9f')) # 11
                    eid_old = 0
                    if hasattr(self, 'eid_old'):
                        eid_old = self.eid_old

                    for i in range(nelements):
                        edata = data[n:n+44]  # 4*11
                        out = struct1.unpack(edata)
                        (eid_device, layer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm) = out
                        eid = eid_device // 10

                        if self.is_debug_file:
                            self.binary_debug.write('  eid=%i; layer=%i; C=[%s]\n' % (eid, layer, ', '.join(['%r' % di for di in out])))

                        if eid != eid_old:
                            # originally initialized to None, the buffer doesnt reset it, so it is the old value
                            obj.add_new_eid_sort1(etype, dt, eid, layer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)
                        else:
                            obj.add_sort1(dt, eid, layer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)
                        eid_old = eid
                        n += 44
                    self.eid_old = eid_old
            elif self.format_code in [2, 3] and self.num_wide == 9:  # TODO: imag? - not done...
                # TODO: vectorize
                raise NotImplementedError('imaginary composite stress?')
                msg = self.code_information()
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_complex)
                if auto_return:
                    return nelements * self.num_wide * 4

                # TODO: this is an OEF result???
                #    furthermore the actual table is calle dout as
                #    'i8si4f4s', not 'i8si3fi4s'
                ntotal = 36
                nelements = ndata // ntotal
                s = self.struct_i
                s2 = Struct(b(self._endian + '8si3fi4s'))
                s3 = Struct(b(self._endian + '8si4f4s'))
                for i in range(nelements):
                    #out = s.unpack(data[n:n + ntotal])
                    eid_device, = s.unpack(data[n:n+4])
                    #t, = s.unpack(data[n:n+4])

                    if eid_device > 0:
                        out = s2.unpack(data[n+4:n+ntotal])
                    else:
                        out1 = s2.unpack(data[n+4:n+ntotal])
                        out = s3.unpack(data[n+4:n+ntotal])
                    (theory, lamid, fp, fm, fb, fmax, fflag) = out

                    if self.is_debug_file:
                        self.binary_debug.write('%s-%s - (%s) + %s\n' % (self.element_name, self.element_type, eid_device, str(out)))
                    obj.add_new_eid(dt, eid, theory, lamid, fp, fm, fb, fmax, fflag)
                    n += ntotal
                raise NotImplementedError('this is a really weird case...')
            else:
                #msg = self.code_information()
                msg = 'OES-COMP-random-numwide=%s numwide_real=11 numwide_imag=9' % (self.num_wide)
                return self._not_implemented_or_skip(data, ndata, msg)

        #=========================
        elif self.element_type == 53: # axial plates - ctriax6
            # 53 - CTRIAX6
            if self.is_stress():
                result_name = 'ctriax_stress'
            else:
                result_name = 'ctriax_strain'
            self._results._found_result(result_name)
            slot = getattr(self, result_name)

            if self.format_code == 1 and self.num_wide == 33: # real
                if self.is_stress():
                    obj_vector_real = RealTriaxStressArray
                else:
                    obj_vector_real = RealTriaxStrainArray

                ntotal = 132  # (1+8*4)*4 = 33*4 = 132
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    self._data_factor = 4
                    return nelements * self.num_wide * 4

                obj = self.obj
                nnodes_all = 4
                if self.use_vector and is_vectorized:
                    n = nelements * self.num_wide * 4

                    itotal = obj.itotal
                    itotal2 = itotal + nelements * nnodes_all
                    ielement = obj.ielement
                    ielement2 = ielement + nelements

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 33)
                    floats1 = floats[:, 1:].reshape(nelements * nnodes_all, 8)

                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 33)
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
                    s1 = Struct(b(self._endian + '2i7f'))  # 36
                    s2 = Struct(b(self._endian + 'i7f'))
                    for i in range(nelements):
                        out = s1.unpack(data[n:n + 36])
                        (eid_device, loc, rs, azs, As, ss, maxp, tmax, octs) = out
                        if self.is_debug_file:
                            self.binary_debug.write('CTRIAX6-53A - %s\n' % (str(out)))
                        eid = eid_device // 10

                        obj.add_sort1(dt, eid, loc, rs, azs, As, ss, maxp, tmax, octs)
                        n += 36
                        for i in range(3):
                            out = s2.unpack(data[n:n + 32])
                            (loc, rs, azs, As, ss, maxp, tmax, octs) = out
                            if self.is_debug_file:
                                self.binary_debug.write('CTRIAX6-53B - %s\n' % (str(out)))
                            obj.add_sort1(dt, eid, loc, rs, azs, As, ss, maxp, tmax, octs)
                            n += 32
            elif self.format_code in [2, 3] and self.num_wide == 37: # imag
                # TODO: vectorize object
                return ndata

                if self.is_stress():
                    raise NotImplementedError('ComplexTriaxStressArray')
                    obj_vector_complex = ComplexTriaxStressArray
                else:
                    raise NotImplementedError('ComplexTriaxStrainArray')
                    obj_vector_complex = ComplexTriaxStrainArray

                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_complex)
                if auto_return:
                    self._data_factor = 4
                    return nelements * self.num_wide * 4

                obj = self.obj
                nnodes_all = 4
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.itotal
                    itotal2 = itotal + nelements * nnodes_all
                    ielement = obj.ielement
                    ielement2 = ielement + nelements

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, numwide_imag)
                    floats1 = floats[:, 1:].reshape(nelements * nnodes_all, 9)

                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, numwide_imag)
                        ints1 = ints[:, 1:].reshape(nelements * nnodes_all, 9)
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
                    if is_magnitude_phase:
                        mag = floats1[:, isave1]
                        phase = floats1[:, isave2]
                        rtheta = radians(phase)
                        real_imag = mag * (cos(rtheta) + 1.j * sin(rtheta))
                    else:
                        real = floats1[:, isave1]
                        imag = floats1[:, isave2]
                        real_imag = real + 1.j * imag

                    obj.data[obj.itime, itotal:itotal2, :] = real_imag
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s1 = Struct(b(self._endian + 'ii8f'))
                    s2 = Struct(b(self._endian + 'i8f'))

                    num_wide = 1 + 4 * 9
                    ntotal = num_wide * 4
                    assert num_wide == self.num_wide, num_wide
                    nelements = ndata // ntotal  # (1+8*4)*4 = 33*4 = 132

                    for i in range(nelements):
                        out = s1.unpack(data[n:n + 40])
                        (eid_device, loc, rsr, rsi, azsr, azsi, Asr, Asi, ssr, ssi) = out
                        eid = eid_device // 10
                        if self.is_debug_file:
                            self.binary_debug.write('CTRIAX6-53 eid=%i\n    %s\n' % (eid, str(out)))

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
                        #obj.add_new_eid(dt, eid, loc, rs, azs, As, ss)

                        n += 40
                        for i in range(3):
                            out = s2.unpack(data[n:n + 36])
                            (loc, rsr, rsi, azsr, azsi, Asr, Asi, ssr, ssi) = out
                            if self.is_debug_file:
                                self.binary_debug.write('    %s\n' % (str(out)))
                            #print("eid=%s loc=%s rs=%s azs=%s as=%s ss=%s" % (
                                #id, loc, rs, azs, As, ss))

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
                            #obj.add(dt, eid, loc, rs, azs, As, ss)
                            n += 36  # 4*8
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
        elif self.element_type == 102: # bush
            # 102-CBUSH
            if self.is_stress():
                result_name = 'cbush_stress'
            else:
                result_name = 'cbush_strain'
            self._results._found_result(result_name)
            slot = getattr(self, result_name)

            if self.format_code == 1 and self.num_wide == 7:  # real
                if self.is_stress():
                    obj_vector_real = RealBushStressArray
                else:
                    obj_vector_real = RealBushStrainArray

                assert self.num_wide == 7, "num_wide=%s not 7" % self.num_wide
                ntotal = 28  # 4*7

                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    return nelements * self.num_wide * 4
                obj = self.obj

                if self.use_vector and is_vectorized:
                    n = nelements * self.num_wide * 4

                    istart = obj.ielement
                    iend = istart + nelements
                    obj._times[obj.itime] = dt

                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 7)
                        eids = ints[:, 0] // 10
                        obj.element[istart:iend] = eids

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 7)
                    #[tx, ty, tz, rx, ry, rz]
                    obj.data[obj.itime, istart:iend, :] = floats[:, 1:]
                else:
                    struct1 = Struct(b(self._endian + 'i6f'))
                    for i in range(nelements):
                        edata = data[n:n + ntotal]
                        out = struct1.unpack(edata)  # num_wide=7
                        if self.is_debug_file:
                            self.binary_debug.write('CBUSH-102 - %s\n' % str(out))

                        (eid_device, tx, ty, tz, rx, ry, rz) = out
                        eid = eid_device // 10

                        obj.add_sort1(dt, eid, tx, ty, tz, rx, ry, rz)
                        n += ntotal
            elif self.format_code in [2, 3] and self.num_wide == 13:  # imag
                if self.is_stress():
                    obj_complex = ComplexCBushStressArray
                else:
                    obj_complex = ComplexCBushStrainArray

                ntotal = 52 # 4*13
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_complex)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 13)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 13)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids

                    if is_magnitude_phase:
                        mag = floats[:, [1, 2, 3, 4, 5, 6]]
                        phase = floats[:, [7, 8, 9, 10, 11, 12]]
                        rtheta = radians(phase)
                        real_imag = mag * (cos(rtheta) + 1.j * sin(rtheta))
                    else:
                        real = floats[:, [1, 2, 3, 4, 5, 6]]
                        imag = floats[:, [7, 8, 9, 10, 11, 12]]
                        real_imag = real + 1.j * imag
                    obj.data[obj.itime, itotal:itotal2, :] = real_imag

                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    struct1 = Struct(b(self._endian + 'i12f'))
                    for i in range(nelements):
                        edata = data[n:n + ntotal]
                        out = struct1.unpack(edata)  # num_wide=7
                        if self.is_debug_file:
                            self.binary_debug.write('CBUSH-102 - %s\n' % str(out))

                        (eid_device,
                         txr, tyr, tzr, rxr, ryr, rzr,
                         txi, tyi, tzi, rxi, ryi, rzi) = out
                        eid = eid_device // 10

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
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type == 40:  # bush
            # 40-CBUSH1D
            if self.is_stress():
                result_name = 'cbush1d_stress_strain'
            else:
                result_name = 'cbush1d_stress_strain'
            self._results._found_result(result_name)
            slot = self.cbush1d_stress_strain

            if self.format_code == 1 and self.num_wide == 8:  # real
                if self.is_stress():
                    obj_vector_real = RealBush1DStressArray
                else:
                    #self.create_transient_object(self.cbush1d_stress_strain, Bush1DStrain)  # undefined
                    raise NotImplementedError('cbush1d_stress_strain; numwide=8')

                ntotal = 32  # 4*8
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.use_vector and is_vectorized:
                    n = nelements * self.num_wide * 4

                    itotal = obj.itotal
                    itotal2 = itotal + nelements
                    itime = obj.itime
                    obj._times[itime] = dt

                    if 1: #obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 8)
                        eids = ints[:, 0] // 10
                        fail = ints[:, 7]
                        obj.element[itotal:itotal2] = eids
                        obj.is_failed[itime, itotal:itotal2, 0] = fail

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 8)
                    #[xxx, fe, ue, ve, ao, ae, ep, xxx]
                    obj.data[itime, itotal:itotal2, :] = floats[:, 1:7]

                    obj.ielement = itotal2
                    obj.itotal = itotal2
                else:
                    struct1 = Struct(b(self._endian + 'i6fi'))
                    for i in range(nelements):
                        edata = data[n:n + 32]
                        out = struct1.unpack(edata)  # num_wide=25
                        if self.is_debug_file:
                            self.binary_debug.write('CBUSH1D-40 - %s\n' % (str(out)))
                        (eid_device, fe, ue, ve, ao, ae, ep, fail) = out
                        eid = eid_device // 10

                        # axial_force, axial_displacement, axial_velocity, axial_stress,
                        # axial_strain, plastic_strain, is_failed
                        obj.add_sort1(dt, eid, fe, ue, ve, ao, ae, ep, fail)
                        n += ntotal
            elif self.format_code in [2, 3] and self.num_wide == 9:  # imag
                # TODO: vectorize object
                ntotal = 36  # 4*9
                nelements = ndata // ntotal

                if self.is_stress():
                    auto_return, is_vectorized = self._create_oes_object4(
                        nelements, result_name, slot, ComplexCBush1DStressArray)
                else:
                    raise NotImplementedError('self.cbush1d_stress_strain; complex strain')

                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.use_vector and is_vectorized:
                    n = nelements * self.num_wide * 4

                    itotal = obj.itotal
                    itotal2 = itotal + nelements
                    itime = obj.itime
                    obj._times[itime] = dt

                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 9)
                        eids = ints[:, 0] // 10
                        obj.element[itotal:itotal2] = eids

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 9)
                    #[fer, uer, aor, aer,
                    # fei, uei, aoi, aei]
                    isave1 = [1, 3, 5, 7]
                    isave2 = [2, 4, 6, 8]
                    if is_magnitude_phase:
                        mag = floats[:, isave1]
                        phase = floats[:, isave2]
                        rtheta = radians(phase)
                        real_imag = mag * (cos(rtheta) + 1.j * sin(rtheta))
                    else:
                        real = floats[:, isave1]
                        imag = floats[:, isave2]
                        real_imag = real + 1.j * imag
                    obj.data[obj.itime, itotal:itotal2, :] = real_imag

                    obj.ielement = itotal2
                    obj.itotal = itotal2
                else:
                    struct1 = Struct(b(self._endian + 'i8f'))
                    for i in range(nelements):
                        edata = data[n:n+ntotal]

                        out = struct1.unpack(edata)  # num_wide=25
                        (eid_device,
                         fer, uer, aor, aer,
                         fei, uei, aoi, aei) = out
                        eid = eid_device // 10

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
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type in [87, 89, 92]:  # nonlinear rods
            # 87-CTUBENL
            # 89-RODNL
            # 92-CONRODNL
            if self.is_stress():
                if self.element_type == 87:
                    result_name = 'nonlinear_ctube_stress'
                    name = 'CTUBENL-87'
                elif self.element_type == 89:
                    result_name = 'nonlinear_crod_stress'
                    name = 'RODNL-89'
                elif self.element_type == 92:
                    result_name = 'nonlinear_conrod_stress'
                    name = 'CONRODNL-92'
                else:
                    raise RuntimeError(self.code_information())
            else:
                if self.element_type == 87:
                    result_name = 'nonlinear_ctube_strain'
                    name = 'CTUBENL-87'
                elif self.element_type == 89:
                    result_name = 'nonlinear_crod_strain'
                    name = 'RODNL-89'
                elif self.element_type == 92:
                    result_name = 'nonlinear_conrod_strain'
                    name = 'CONRODNL-92'
                else:
                    raise RuntimeError(self.code_information())
            self._results._found_result(result_name)
            slot = getattr(self, result_name)

            if self.format_code == 1 and self.num_wide == 7:  # real
                ntotal = 28  #  7*4 = 28
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, RealNonlinearRodArray)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                #if self.is_debug_file:
                    #self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    #self.binary_debug.write('  element1 = [eid_device, layer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)]\n')
                    #self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                if self.use_vector and is_vectorized:
                    n = nelements * self.num_wide * 4
                    istart = obj.itotal
                    iend = istart + nelements
                    obj._times[obj.itime] = dt

                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 7)
                        eids = ints[:, 0] // 10
                        obj.element[istart:iend] = eids
                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 7)
                    #[axial_stress, equiv_stress, total_strain,
                    # eff_plastic_creep_strain, eff_creep_strain, linear_torsional_stresss]
                    obj.data[obj.itime, istart:iend, :] = floats[:, 1:]
                else:
                    struct1 = Struct(b(self._endian + 'i6f'))  # 1+6=7
                    for i in range(nelements):
                        edata = data[n:n+ntotal]
                        out = struct1.unpack(edata)

                        (eid_device, axial_stress, equiv_stress, total_strain,
                         eff_plastic_creep_strain, eff_creep_strain, linear_torsional_stresss) = out
                        eid = eid_device // 10
                        if self.is_debug_file:
                            self.binary_debug.write('%s - %s\n' % (name, str(out)))
                        obj.add_sort1(dt, eid, axial_stress, equiv_stress, total_strain,
                                      eff_plastic_creep_strain, eff_creep_strain, linear_torsional_stresss)
                        n += ntotal
            else:
                raise RuntimeError(self.code_information())

        elif self.element_type in [224, 225]: # nonlinear spring
            # 224-CELAS1
            # 225-CELAS3
            # NonlinearSpringStress
            numwide_real = 3
            if self.is_stress():
                if self.element_type == 224:
                    result_name = 'nonlinear_celas1_stress'
                    slot = self.nonlinear_celas1_stress
                elif self.element_type == 225:
                    result_name = 'nonlinear_celas3_stress'
                    slot = self.nonlinear_celas3_stress
            else:
                raise NotImplementedError('NonlinearSpringStrain')

            self._results._found_result(result_name)
            if self.num_wide == numwide_real:
                assert self.num_wide == 3, "num_wide=%s not 3" % self.num_wide
                ntotal = 12  # 4*3
                nelements = ndata // ntotal

                if self.is_stress():
                    auto_return, is_vectorized = self._create_oes_object4(
                        nelements, result_name, slot, RealNonlinearSpringStressArray)
                else:
                    raise NotImplementedError('NonlinearSpringStrainArray') # undefined

                if auto_return:
                    return nelements * self.num_wide * 4
                obj = self.obj

                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement = obj.ielement
                    ielement2 = obj.ielement + nelements
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, numwide_real)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[ielement:ielement2] = eids

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, numwide_real)

                    #[force, stress]
                    obj.data[obj.itime, ielement:ielement2, :] = floats[:, 1:]
                    obj.itotal = ielement2
                    obj.ielement = ielement2
                else:
                    struct1 = Struct(b(self._endian + 'i2f'))
                    for i in range(nelements):
                        edata = data[n:n+ntotal]
                        out = struct1.unpack(edata)  # num_wide=3
                        (eid_device, force, stress) = out
                        eid = eid_device // 10
                        if self.is_debug_file:
                            self.binary_debug.write('%s-%s - %s\n' % (self.element_name, self.element_type, str(out)))
                        obj.add_sort1(dt, eid, force, stress)
                        n += ntotal
            else:
                raise RuntimeError(self.code_information())

        elif self.element_type == 35:
            # 35-CON
            return ndata
        elif self.element_type in [60, 61]:
            # 60-DUM8
            # 61-DUM9
            return ndata
        elif self.element_type == 69:  # cbend
            # 69-CBEND
            return ndata
        elif self.element_type == 86:  # cgap
            # 86-GAPNL
            if self.is_stress():
                result_name = 'nonlinear_cgap_stress'
            else:
                result_name = 'nonlinear_cgap_strain'
            self._results._found_result(result_name)
            slot = getattr(self, result_name)
            #print('self.nonlinear_factor = ', self.nonlinear_factor)
            if self.format_code == 1 and self.num_wide == 11:  # real?
                if self.is_stress():
                    obj_vector_real = NonlinearGapStressArray
                else:
                    raise NotImplementedError('NonlinearGapStrain')

                ntotal = 44  # 4*11
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.use_vector and is_vectorized:
                    n = nelements * self.num_wide * 4

                    ielement = obj.ielement
                    ielement2 = ielement + nelements
                    obj._times[obj.itime] = dt

                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 11)
                        eids = ints[:, 0] // 10
                        obj.element[ielement:ielement2] = eids

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 11)
                    # skipping [form1, form2]
                    #[cpx, shy, shz, au, shv, shw, slv, slp]
                    obj.data[obj.itime, ielement:ielement2, :] = floats[:, 1:9]
                else:
                    struct1 = Struct(b(self._endian + 'i8f4s4s'))
                    for i in range(nelements):
                        edata = data[n:n + ntotal]

                        out = struct1.unpack(edata)  # num_wide=25
                        (eid_device, cpx, shy, shz, au, shv, shw, slv, slp, form1, form2) = out
                        eid = eid_device // 10
                        if self.is_debug_file:
                            self.binary_debug.write('CGAPNL-86 - %s\n' % str(out))
                        obj.add_sort1(dt, eid, cpx, shy, shz, au, shv, shw, slv, slp, form1, form2)
                        n += ntotal
            else:
                raise RuntimeError(self.code_information())

        elif self.element_type == 94:
            if self.read_mode == 0:
                return ndata
            # 94-BEAMNL
            numwide_real = 51
            numwide_random = 0

            if self.is_stress():
                result_name = 'nonlinear_cbeam_stress'
                slot = self.nonlinear_cbeam_stress
            else:
                result_name = 'nonlinear_cbeam_strain'
                slot = self.nonlinear_cbeam_strain
            self._results._found_result(result_name)

            if self.format_code == 1 and self.num_wide == numwide_real:
                msg = result_name
                if self.is_stress():
                    obj_vector_real = RealNonlinearBeamStressArray
                else:
                    raise NotImplementedError('Nonlinear CBEAM Strain...this should never happen')

                ntotal = numwide_real * 4
                nelements = ndata // ntotal

                nlayers = nelements * 8
                auto_return, is_vectorized = self._create_oes_object4(
                    nlayers, result_name, slot, obj_vector_real)
                if auto_return:
                    self._data_factor = 8
                    return ndata
                obj = self.obj
                if self.is_debug_file:
                    self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    #self.binary_debug.write('  #elementi = [eid_device, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,\n')
                    #self.binary_debug.write('                           s1b, s2b, s3b, s4b, smaxb, sminb,        MSc]\n')
                    #self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)


                struct1 = Struct(b(self._endian + '2i 4s5f 4s5f 4s5f 4s5f i 4s5f 4s5f 4s5f 4s5f'))  # 2 + 6*8 + 1 = 51
                for i in range(nelements):  # num_wide=51
                    edata = data[n:n + 204]
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
                    assert out[3-1] == b'   C', out[3-1]
                    assert out[9-1] == b'   D', out[9-1]
                    assert out[15-1] == b'   E', out[15-1]
                    assert out[21-1] == b'   F', out[21-1]

                    # B
                    assert out[28-1] == b'   C', out[28-1]
                    assert out[34-1] == b'   D', out[34-1]
                    assert out[40-1] == b'   E', out[40-1]
                    assert out[46-1] == b'   F', out[46-1]

                    eid_device = out[0]
                    eid = eid_device // 10
                    obj.add_new_eid_sort1(dt, eid, out)
                    n += 204

            elif self.format_code == 1 and self.num_wide == numwide_random:  # random
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
            else:
                raise RuntimeError(self.code_information())

        elif self.element_type in [85, 91, 93]:
            if self.read_mode == 1:
                return ndata
            # 85-TETRANL
            # 91-PENTANL
            # 93-HEXANL

            #85: 2 + (18 - 2) * 5,  # Nonlinear CTETRA
            #91: 4 + (25 - 4) * 7,  # Nonlinear CPENTA
            #93: 4 + (25 - 4) * 9,  # Nonlinear CHEXA
            if self.element_type == 85:
                etype = 'CTETRANL'
                nnodes = 5
            elif self.element_type == 91:
                etype = 'CPENTANL'
                nnodes = 7
            elif self.element_type == 93:
                etype = 'CHEXANL'
                nnodes = 9
            else:
                raise RuntimeError(self.code_information())

            numwide_real = 4 + (25 - 4) * nnodes  # real???
            numwide_random = 2 + (18 - 2) * nnodes  # imag???
            #print("format_code=%s numwide=%s numwide_real=%s numwide_random=%s" % (self.format_code, self.num_wide, numwide_real, numwide_random))

            #numwide_real = 0
            #numwide_imag = 2 + 16 * nnodes
            #ntotal = 8 + 64 * nnodes

            if self.format_code == 1 and self.num_wide == numwide_real:
                ntotal = numwide_real * 4
                #if self.is_stress():
                    #self.create_transient_object(self.nonlinearPlateStress, NonlinearSolid)
                #else:
                    #self.create_transient_object(self.nonlinearPlateStrain, NonlinearSolid)
                #self.handle_results_buffer(self.OES_CQUAD4NL_90, resultName, name)
                raise RuntimeError('OES_CQUAD4NL_90')
            elif self.format_code == 1 and self.num_wide == numwide_random:  # random
                #raise RuntimeError(self.code_information())
            #elif self.format_code in [2, 3] and self.num_wide == numwide_imag:  # imag
                ntotal = numwide_random * 4
                #if self.is_stress():
                    #self.create_transient_object(self.nonlinearPlateStress, NonlinearSolid)
                #else:
                    #self.create_transient_object(self.nonlinearPlateStrain, NonlinearSolid)

                n = 0
                s1 = Struct(b(self._endian + 'i4s'))
                s2 = Struct(b(self._endian + 'i15f'))
                nelements = ndata // ntotal
                for i in range(nelements):  # 2+16*9 = 146 -> 146*4 = 584
                    edata = data[n:n+8]
                    n += 8

                    out = s1.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('%s-%s - %s\n' % (etype, self.element_type, str(out)))
                    (eid_device, ctype) = out
                    eid = eid_device // 10

                    for i in range(nnodes):
                        edata = data[n:n+64]
                        n += 64
                        out = s2.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('%s-%sB - %s\n' % (etype, self.element_type, str(out)))

                        assert len(out) == 16
                        (grid,
                         sx, sy, sz, sxy, syz, sxz, se, eps, ecs,
                         ex, ey, ez, exy, eyz, exz) = out
            else:
                #msg = self.code_information()
                msg = "format_code=%s numwide=%s numwide_real=%s numwide_random=%s" % (
                    self.format_code, self.num_wide, numwide_real, numwide_random)
                #return self._not_implemented_or_skip(data, ndata, msg)
                raise RuntimeError(self.code_information())

        elif self.element_type == 100:  # bars
            # 100-BARS
            if self.is_stress():
                result_name = 'cbar_stress_10nodes'
            else:
                result_name = 'cbar_strain_10nodes'

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            slot = getattr(self, result_name)

            if self.format_code == 1 and self.num_wide == 10:  # real
                if self.is_stress():
                    obj_vector_real = RealBar10NodesStressArray
                else:
                    obj_vector_real = RealBar10NodesStrainArray

                ntotal = 10 * 4
                nelements = ndata // ntotal

                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    return ndata

                if self.is_debug_file:
                    self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    self.binary_debug.write('  #elementi = [eid_device, sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS]\n')
                    self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)
                obj = self.obj

                if self.use_vector and is_vectorized:
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
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 10)
                        eids = ints[:, 0] // 10
                        obj.element[istart:iend] = eids

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 10)
                    #[sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS]
                    obj.data[obj.itime, istart:iend, :] = floats[:, 1:]
                else:
                    struct1 = Struct(b(self._endian + 'i9f'))
                    for i in range(nelements):
                        edata = data[n:n+ntotal]
                        out = struct1.unpack(edata)
                        (eid_device, sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS) = out
                        eid = eid_device // 10
                        if self.is_debug_file:
                            self.binary_debug.write('  eid=%i; C%i=[%s]\n' % (eid, i, ', '.join(['%r' % di for di in out])))
                        n += ntotal
                        obj.add_new_eid(self.element_name, dt, eid,
                                        sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS)
            else:
                raise RuntimeError(self.code_information())
        elif self.element_type == 101:
            # 101-AABSF
            return ndata
        elif self.element_type in [140, 201]:
            # 140-HEXA8FD, 201-QUAD4FD
            return ndata
        elif self.element_type in [145, 146, 147]:
            # TODO: vectorize
            if self.read_mode == 1:
                return ndata
            # 145-VUHEXA  (8 nodes)
            # 146-VUPENTA (6 nodes)
            # 147-VUTETRA (4 nodes)
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
            else:
                raise RuntimeError(self.code_information())

            #num_wideA = 2 + 12 * nnodes
            #ntotal = 8 + 48 * nnodes

            if self.format_code == -1 and self.num_wide == 33: # real
                # assuming TETRA...
                numwide_a = 2 + (14 - 2) * nnodes  # 50
                numwide_b = 2 + (9 - 2) * nnodes  # 30
                numwide_c = 2 + 13 * nnodes  # 54
                if self.num_wide == numwide_a:
                    ntotal = numwide_a * 4
                    s1 = self.struct_2i
                    s2 = Struct(b(self._endian + 'i11f'))
                    nelements = ndata // ntotal  # 2+16*9 = 146 -> 146*4 = 584
                    for i in range(nelements):
                        edata = data[n:n+8]
                        out = s1.unpack(edata)
                        (eid_device, parent_id) = out
                        eid = eid_device // 10

                        for i in range(nnodes):
                            edata = data[n:n+48]
                            out = s2.unpack(edata)
                            if self.is_debug_file:
                                self.binary_debug.write('%s-%s - %s\n' % (etype, self.element_type, str(out)))
                            assert len(out) == 12
                            (grid, xnorm, ynorm, znorm, txy, tyz, txz,
                             prin1, prin2, prin3, smean, vono_roct) = out
                    return ndata
                elif self.num_wide == numwide_b:
                    ntotal = numwide_b * 4
                    nelements = ndata // ntotal
                    n = nelements * ntotal
                elif self.num_wide == numwide_c:
                    ntotal = numwide_c * 4
                    nelements = ndata // ntotal
                    n = nelements * ntotal
                else:
                    msg = 'numwide=%s A=%s B=%s C=%s' % (self.num_wide, numwide_a, numwide_b, numwide_c)
                    raise RuntimeError(self.code_information())
            else:
                #raise RuntimeError(self.code_information())
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type == 139:
            # 139-QUAD4FD
            if self.read_mode == 1:
                return ndata
            if self.format_code == -4 and self.num_wide == 30:
                # TODO: vectorize
                if self.is_stress():
                    self.create_transient_object(self.hyperelastic_cquad4_strain, HyperelasticQuad)
                    result_name = 'hyperelastic_cquad4_strain'
                else:
                    msg = 'HyperelasticQuad???'
                    return self._not_implemented_or_skip(data, ndata, msg)

                self._results._found_result(result_name)
                n = 0
                ntotal = 120  # 36+28*3
                s1 = Struct(b(self._endian + 'i4si6f'))  # 1 + 4+1+6 = 12
                s2 = Struct(b(self._endian + 'i6f'))
                nelements = ndata // ntotal
                obj = self.obj
                for i in range(nelements):
                    edata = data[n:n+36]  # 4*9
                    out = s1.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('CQUAD4FD-139A- %s\n' % (str(out)))

                    (eid_device, Type, ID, sx, sy, sxy, angle, smj, smi) = out
                    eid = eid_device // 10
                    obj.add_new_eid(dt, [eid, Type, sx, sy, sxy, angle, smj, smi])
                    n += 36

                    for i in range(3):
                        edata = data[n:n + 28]  # 4*7
                        out = s2.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('               %s\n' % (str(out)))
                        (ID, sx, sy, sxy, angle, smj, smi) = out
                        obj.add(dt, eid, out)
                        n += 28
            else:
                msg = 'numwide=%s' % self.num_wide
                return self._not_implemented_or_skip(data, ndata, msg)
        elif self.element_type == 189:
            # 189-VUQUAD
            if self.element_type == 189:  # VQUAD
                if self.read_mode == 1:
                    return ndata
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
                #eType = 'CQUAD8'
            #elif self.element_type == 82:  # CQUADR
                #ntotal = 348  # 2+17*5 = 87 -> 87*4 = 348
                #nnodes = 4    # centroid + 4 corner points
                #eType = 'CQUADR'
            #elif self.element_type == 75:  # CTRIA6
                #ntotal = 280  # 2+17*3 = 70 -> 70*4 = 280
                #nnodes = 3    # centroid + 3 corner points
                #eType = 'CTRIA6'
            #elif self.element_type == 70:  # CTRIAR
                #ntotal = 280  # 2+17*3 = 70 -> 70*4 = 280
                #nnodes = 3    # centroid + 3 corner points
                #eType = 'CTRIAR'
            else:
                return self._not_implemented_or_skip(data, ndata, self.code_information())

            numwide_real = 6 + (23 - 6) * nnodes  # real???
            numwide_imag = 6 + (33 - 6) * nnodes  # imag???

            if self.format_code == 1 and self.num_wide == numwide_real:  # real???
                ntotal = numwide_real * 4
                s2 = Struct(b(self._endian + '3i4s2i'))
                s3 = Struct(b(self._endian + 'i16f'))
                nelements = ndata // ntotal
                for i in range(nelements):
                    out = s2.unpack(data[n:n + 24])
                    (eid_device, parent, coord, icord, theta, itype) = out
                    n += 24
                    eid = eid_device // 10
                    edata = data[n:n + 68]
                    out = s3.unpack(edata)  # len=17*4
                    n += 68

                    if self.is_debug_file:
                        self.binary_debug.write('%s-%s - %s\n' % (etype, self.element_type, str(out)))

                    #obj.add_new_node(dt, eid, parent, coord, icord, theta, itype)
                    #obj.add_new_eid(eType, dt, eid, parent, coord, icord, theta, itype)
                    for node_id in range(nnodes - 1):  # nodes pts
                        edata = data[n:n + 68]
                        n += 68
                        out = s3.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('              %s\n' % (str(out)))

                        (vuid, dummy, dummy2, msx, msy, mxy, dummy3, dummy4, dummy5,
                         bcx, bcy, bcxy, tyz, tzx, dummy6, dummy7, dummy8) = out
                        #obj.add(vuid, dummy, dummy2, msx, msy, mxy,
                                     #dummy3, dummy4, dummy5,
                                     #bcx, bcy, bcxy, tyz, tzx,
                                     #dummy6, dummy7, dummy8)
            elif self.num_wide == numwide_imag:
                ntotal = numwide_imag * 4
                nelements = ndata // ntotal
                n = nelements * ntotal
            else:
                msg = 'numwide=%s' % self.num_wide
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type in [47, 48, 189, 190]:
            # 47-AXIF2
            # 48-AXIF3
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
        elif self.element_type in [160, 161, 162, 163, 164, 165, 166, 167, 168,
                                   169, 170, 171, 172, 202,
                                   204, 218, 211, 213, 214,
                                   216, 217, 219, 220, 221, 222, 223,
                                   226, 232, 233, 235]:
            # 160-PENTA6FD
            # 161-TETRA4FD
            # 162-TRIA3FD
            # 163-HEXAFD
            # 164-QUADFD
            # 165-PENTAFD
            # 166-TETRAFD
            # 167-TRIAFD
            # 168-TRIAX3FD
            # 169-TRIAXFD
            # 170-QUADX4FD
            # 171-QUADXFD
            # 172-QUADRNL
            # 202-HEXA8FD
            # 204-PENTA6FD
            # 211-TRIAFD
            # 213-TRIAXFD
            # 214-QUADX4FD
            # 216-TETRA4FD
            # 217-TRIA3FD
            # 218-HEXAFD
            # 219-QUADFD
            # 220-PENTAFD
            # 221-TETRAFD
            # 223-QUADXFD
            # 222-TRIAX3FD
            # 226-BUSH
            # 232-QUADRLC
            # 233-TRIARLC
            # 235-CQUADR
            return ndata
        #elif self.element_type in [255]:
            #return ndata
        elif self.element_type in [271, 275]:
            # 271-CPLSTN3
            # 275-CPLSTS3
            msg = self.code_information()
            return self._not_implemented_or_skip(data, ndata, msg)
            if self.element_type == 271:
                result_name = 'cplstn3'
                nnodes = 1
                ntotal = 4 * 6
            elif self.element_type == 275:
                result_name = 'cplsts3'
                nnodes = 1
                ntotal = 4 * 6
            else:
                raise RuntimeError(self.code_information())
            if self.is_stress():
                obj_vector_real = RealCPLSTRNPlateStressArray
                result_name += '_stress'
            else:
                obj_vector_real = RealCPLSTRNPlateStrainArray
                result_name += '_strain'

            numwide_real = ntotal // 4
            if self.format_code == 1 and self.num_wide == numwide_real:
                #ntotal = 4 * (1 + 6 * (nnodes))
                nelements = ndata // ntotal

                #self._data_factor = 10  # TODO: why is this 10?
                if self.is_stress():
                    obj_vector_real = RealCPLSTRNPlateStressArray
                    #result_name = 'cplstn3_stress'
                else:
                    obj_vector_real = RealCPLSTRNPlateStressArray
                    #result_name = 'cplstn3_strain'
                slot = getattr(self, result_name)

                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                #if self.use_vector and is_vectorized:
                n = nelements * self.num_wide * 4

                istart = obj.itotal
                iend = istart + nelements
                obj._times[obj.itime] = dt

                if obj.itime == 0:
                    print(fromstring(data, dtype=self.idtype).size)
                    print('nelements=%s numwide=%s' % (nelements, numwide_real))
                    print('ndata=', ndata)
                    print('self.element_name=%s' % self.element_name)
                    ints = fromstring(data, dtype=self.idtype).reshape(nelements, numwide_real)
                    eids = ints[:, 0] // 10
                    obj.element[istart:iend] = eids

                floats = fromstring(data, dtype=self.fdtype).reshape(nelements, numwide_real)
                results = floats[:, 1:]
                print('results.shape', results.shape)

                #[oxx, oyy, ozz, txy, ovm]
                obj.data[obj.itime, istart:iend, :] = results
            else:
                msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type in [276, 277, 278]:
            # 276-CPLSTS4
            # 277-CPLSTS6
            # 278-CPLSTS8
            msg = self.code_information()
            return self._not_implemented_or_skip(data, ndata, msg)
            if self.element_type == 276:
                result_name = 'cplsts4'
                nnodes = 5  # 4 + 1
                ntotal = 4 * 32
            elif self.element_type == 277:
                result_name = 'cplsts6'
                nnodes = 4
                ntotal = 4 * 26
            elif self.element_type == 278:
                result_name = 'cplsts8'
                nnodes = 5
                ntotal = 4 * 32
            else:
                raise RuntimeError(self.code_information())

            if self.is_stress():
                obj_vector_real = RealCPLSTRNPlateStressArray
                result_name += '_stress'
            else:
                obj_vector_real = RealCPLSTRNPlateStrainArray
                result_name += '_strain'

            numwide_real = 2 + 6 * (nnodes)
            assert ntotal // 4 == numwide_real, 'notal/4=%s numwide_real=%s\n%s' % (
                ntotal // 4, numwide_real, self.code_information())

            ntotal = numwide_real * 4
            if self.format_code == 1 and self.num_wide == numwide_real:
                nelements = ndata // ntotal

                #self._data_factor = 10  # TODO: why is this 10?
                if self.is_stress():
                    obj_vector_real = RealCPLSTRNPlateStressArray
                    #result_name = 'cplstn3_stress'
                else:
                    obj_vector_real = RealCPLSTRNPlateStressArray
                    #result_name = 'cplstn3_strain'
                slot = getattr(self, result_name)

                nlayers = nelements * nnodes
                auto_return, is_vectorized = self._create_oes_object4(
                    nlayers, result_name, slot, obj_vector_real)
                if auto_return:
                    self._data_factor = nnodes
                    return nelements * self.num_wide * 4

                obj = self.obj
                #if self.use_vector and is_vectorized:
                n = nlayers * self.num_wide * 4

                istart = obj.itotal
                iend = istart + nlayers
                obj._times[obj.itime] = dt

                if obj.itime == 0:
                    print(fromstring(data, dtype=self.idtype).size)
                    print('nelements=%s numwide=%s' % (nelements, numwide_real))
                    ints = fromstring(data, dtype=self.idtype).reshape(nelements, numwide_real)
                    eids = ints[:, 0] // 10
                    #obj.element[istart:iend] = eids

                floats = fromstring(data, dtype=self.fdtype).reshape(nelements, numwide_real)
                print('floats[:, 2:].shape', floats[:, 2:].shape)
                print('nnelements=%s nnodes=%s numwide//nodes=%s' % (nelements, nnodes, (numwide_real-2) / nnodes))
                results = floats[:, 2:].reshape(nelements, nnodes * 6)

                #[oxx, oyy, ozz, txy, ovm]
                obj.data[obj.itime, istart:iend, :] = results
            else:
                msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
                return self._not_implemented_or_skip(data, ndata, msg)
        else:
            msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
            return self._not_implemented_or_skip(data, ndata, msg)

        assert ndata > 0, ndata
        assert nelements > 0, 'nelements=%r element_type=%s element_name=%r' % (nelements, self.element_type, self.element_name)
        #assert ndata % ntotal == 0, '%s n=%s nwide=%s len=%s ntotal=%s' % (self.element_name, ndata % ntotal, ndata % self.num_wide, ndata, ntotal)
        assert self.num_wide * 4 == ntotal, 'numwide*4=%s ntotal=%s' % (self.num_wide * 4, ntotal)
        assert self.thermal == 0, "thermal = %%s" % self.thermal
        assert n > 0, "n = %s result_name=%s" % (n, result_name)
        return n

    def OESRT_CQUAD4_95(self, data, ndata):
        """unsupported element"""
        assert self.num_wide == 9, "num_wide=%s not 9" % self.num_wide
        ntotal = 36  # 4*9

        n = 0
        struct1 = Struct(b(self._endian + 'i8si3fi4s'))
        nelements = ndata // ntotal
        obj = self.obj
        for i in range(nelements):
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
            print("eid=%s strength_ratio_ply=%g failure_index_bonding=%s strength_ratio_bonding=%s" % (
                eid, strength_ratio_ply, failure_index_bonding, strength_ratio_bonding))
            #obj.add_new_eid(element_name, dt, eid, force, stress)
            n += ntotal

    def _create_nodes_object(self, nnodes, result_name, slot, obj_vector):
        """same as _create_oes_object4 except it adds to the nnodes parameter"""
        auto_return = False
        #is_vectorized = True
        is_vectorized = self._is_vectorized(obj_vector, slot)
        #print("vectorized...read_mode=%s...%s; %s" % (self.read_mode, result_name, is_vectorized))

        if is_vectorized:
            if self.read_mode == 1:
                #print('oes-self.nonlinear_factor =', self.nonlinear_factor)
                #print(self.data_code)
                self.create_transient_object(slot, obj_vector)
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
                self.obj.build()

            else:  # not vectorized
                auto_return = True
        else:
            auto_return = True
        return auto_return, is_vectorized

    def _create_ntotal_object(self, ntotal, result_name, slot, obj_vector):
        """same as _create_oes_object4 except it adds to the ntotal parameter"""
        auto_return = False
        #is_vectorized = True
        is_vectorized = self._is_vectorized(obj_vector, slot)
        #print("vectorized...read_mode=%s...%s; %s" % (self.read_mode, result_name, is_vectorized))

        if is_vectorized:
            if self.read_mode == 1:
                #print('oes-self.nonlinear_factor =', self.nonlinear_factor)
                #print(self.data_code)
                self.create_transient_object(slot, obj_vector)
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
                self.obj.build()

            else:  # not vectorized
                auto_return = True
        else:
            auto_return = True
        return auto_return, is_vectorized

    def _is_vectorized(self, obj_vector, slot_vector):
        """
        Checks to see if the data array has been vectorized

        Parameters
        ----------
        obj_vector:  the object to check
            (obj or None; None happens when vectorization hasn't been implemented)
        slot_vector: the dictionary to put the object in
            (dict or None; None happens when obj hasn't been implemented)

        Returns
        -------
        is_vectorized : bool
            should the data object be vectorized

        .. note :: the Vectorized column refers to the setting given by the user
        """
        is_vectorized = False
        if self.is_vectorized:
            if obj_vector is not None:
                is_vectorized = True
        return is_vectorized
