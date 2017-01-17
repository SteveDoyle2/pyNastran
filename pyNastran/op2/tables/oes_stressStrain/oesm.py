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

from pyNastran.op2.op2_common import OP2Common, apply_mag_phase
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


class OESM(OP2Common):
    """
    Defines  the OES class that is used to read stress/strain data
    """
    def __init__(self):
        OP2Common.__init__(self)
        self.ntotal = 0

    def _read_oesm1_3(self, data, ndata):
        """
        reads OESM1 subtable 3
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

        self._parse_stress_code_to_stress_bits()
        self._write_debug_bits()
        if self.is_debug_file:
            self.binary_debug.write('catttt!')
        assert isinstance(self.format_code, int), self.format_code
        #print('self.nonlinear_factor =', self.nonlinear_factor)
        #assert self.num_wide != 146, self.code_information()

    def _read_oesm1_4(self, data, ndata):
        """
        Reads the Stress Table 4
        """
        #assert self.isubtable == -4, self.isubtable
        #if self.is_debug_file:
            #self.binary_debug.write('  element_name = %r\n' % self.element_name)
        #print "element_name =", self.element_name
        assert isinstance(self.format_code, int), self.format_code
        assert self.is_stress() is True, self.code_information()
        self.data_code['is_stress_flag'] = True
        self.data_code['is_strain_flag'] = False

        if self.isubcase not in self.case_control_deck.subcases:
            self.subcase = self.case_control_deck.create_new_subcase(self.isubcase)
        self.subcase.add_op2_data(self.data_code, 'STRESS/STRAIN', self.log)

        if self.is_sort1():
            n = self._read_oesm1_4_sort1(data, ndata)
        else:
            msg = self.code_information()
            n = self._not_implemented_or_skip(data, ndata, msg)
        return n

    def _read_oesm2_3(self, data, ndata):
        """
        reads OESM1 subtable 3
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

        self._parse_stress_code_to_stress_bits()
        #self._fix_sort2(data)

    def _read_oesm2_4(self, data, ndata):
        return self._table_passer(data, ndata)

    def _read_ostrm1_4(self, data, ndata):
        """
        Reads the Strain Table 4
        """
        #assert self.isubtable == -4, self.isubtable
        #if self.is_debug_file:
            #self.binary_debug.write('  element_name = %r\n' % self.element_name)
        #print "element_name =", self.element_name
        assert self.is_strain() is True, self.code_information()
        self.data_code['is_stress_flag'] = False
        self.data_code['is_strain_flag'] = True

        if self.isubcase not in self.case_control_deck.subcases:
            self.subcase = self.case_control_deck.create_new_subcase(self.isubcase)
        self.subcase.add_op2_data(self.data_code, 'STRESS/STRAIN', self.log)

        if self.is_sort1():
            n = self._read_ostrm1_4_sort1(data, ndata)
        else:
            msg = self.code_information()
            n = self._not_implemented_or_skip(data, ndata, msg)
        return n

    def _read_oesm1_4_sort1(self, data, ndata):
        """
        Reads OES1 subtable 4
        """
        #if self.num_wide == 146:
            #assert self.num_wide != 146
            #assert ndata != 146, self.code_information()
        assert isinstance(self.format_code, int), self.format_code
        assert self.is_sort1() is True
        if self.thermal == 0:
            n = self._read_oesm1_loads(data, ndata)
        elif self.thermal == 1:
            n = self._read_oesm1_thermal(data, ndata)
        else:
            msg = 'thermal=%s' % self.thermal
            n = self._not_implemented_or_skip(data, ndata, msg)
        return n

    #@_print_obj_name_on_crash
    def _read_ostrm1_4_sort1(self, data, ndata):
        """
        Reads OSTRM1 subtable 4
        """
        #if self.num_wide == 146:
            #assert self.num_wide != 146
            #assert ndata != 146, self.code_information()
        assert self.is_sort1() is True
        if self.thermal == 0:
            n = self._read_oesm1_loads(data, ndata)
        #elif self.thermal == 1:
            #n = self._read_oesm1_thermal(data, ndata)
        else:
            msg = 'thermal=%s' % self.thermal
            n = self._not_implemented_or_skip(data, ndata, msg)
        return n

    #def _read_oesm1_thermal(self, data, ndata):
        #"""
        #Reads OES self.thermal=1 tables; uses a hackish method to just skip the table.
        #"""
        #return ndata

    #def _read_ostrm1_thermal(self, data, ndata):
        #"""
        #Reads OSTR self.thermal=1 tables; uses a hackish method to just skip the table.
        #"""
        #return ndata

    def get_stress_mapper(self):
        stress_mapper = {
            # element_type, format_code, num_wide

            # rods
            (1, 2, 5, 'OESVM1') : ('crod', 'NA'),
            (10, 2, 5, 'OESVM1') : ('conrod', 'NA'),

            # springs
            (12, 2, 3, 'OESVM1') : ('celas2', 'NA'),

            # beam
            (2, 2, 111, b'OESVM1') : ('cbeam', 'NA'),

            # bar
            (34, 2, 19, b'OESVM1') : ('cbar', 'NA'),

            # shells
            (74, 2, 17, 'OESVM1') : ('ctria3', 'NA'),
            (144, 2, 87, b'OESVM1') : ('cquad4', 'NA'),

            # composites
            (95, 2, 13, b'OESVM1C') : ('cquad4', 'NA'),
            (97, 2, 13, b'OESVM1C') : ('ctria3', 'NA'),

            # bush
            (102, 2, 13, b'OESVM1') : ('cbush', 'NA'),

            # solids
            (39, 2, 74, 'OESVM1') : ('ctetra', 'NA'),
            (67, 2, 130, b'OESVM1') : ('chexa', 'NA'),
            (68, 2, 102, b'OESVM1') : ('cpenta', 'NA'),
        }

        try:
            return stress_mapper[self.element_type, self.format_code, self.num_wide, self.table_name]
        except KeyError:
            print(self.code_information())
            raise
            return None, None

    def _apply_oesm_ato_crm_psd_rms_no(self, result_name):
        """
        Appends a keyword onto the result_name in order to handle random results
        without 100 if loops.
        keywords = {_ATO, _CRM, _PSD, _RMS, _NO}

        Do this:
            result_name = 'cbar_stress'
            table_name = 'OEFPSD1'
            result_name, is_random = _apply_oef_ato_crm_psd_rms_no(self, result_name)
            slot = getattr(self, result_name)

        Or this:
            result_name = 'cbar_stress_PSD'
            slot = self.cbar_stress_PSD
        """
        is_random = True
        if self.table_name in [b'OES1', b'OES1X1', b'OES1X', b'OES1C', b'OESCP', b'OESRT',
                               b'OSTR1X', b'OSTR1C',
                               b'OESTRCP', b'OESNLXR', b'OESNLXD',
                               b'OESVM1', b'OSTRVM1', b'OESVM1C']:
            is_random = False
        elif self.table_name in [b'OESCRM1', b'OESCRM2']:
            assert self.table_code in [504], self.code_information()
            result_name += '_CRM'
        elif self.table_name in [b'OESPSD1', b'OESPSD2']:
            assert self.table_code in [604], self.code_information()
            result_name += '_PSD'
        elif self.table_name in [b'OESRMS1', b'OESRMS2']:
            assert self.table_code in [804], self.code_information()
            result_name += '_RMS'
        elif self.table_name in [b'OESNO1', b'OESNO2']:
            assert self.table_code in [904], self.code_information()
            result_name += '_NO'
        else:
            raise NotImplementedError(self.code_information())
        #print(result_name, self.table_name)
        return result_name, is_random

    def _read_oesm1_loads(self, data, ndata):
        """
        Reads OES self.thermal=0 stress/strain
        """
        #self._apply_oes_ato_crm_psd_rms_no('') # TODO: just testing
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

        #if self.is_stress():
            #_result_name, _class_obj = self.get_stress_mapper()

        if self._results.is_not_saved(result_name):
            return ndata
        if self.element_type in [1, 3, 10]:  # rods
            # 1-CROD
            # 3-CTUBE
            # 10-CONROD
            if self.is_stress():
                obj_vector_real = RealRodStressArray
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

            result_name, is_random = self._apply_oes_ato_crm_psd_rms_no(result_name)
            slot = getattr(self, result_name)
            #if self.format_code in [2, 3] and self.num_wide == 5: # imag
            if self.format_code == 2 and self.num_wide == 5:
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

            #if self.format_code == 2 and self.num_wide == 67: # random
                #msg = self.code_information()
                #if self.is_debug_file:
                    #self.binary_debug.write('skipping OES-CBEAM\n')
                #return self._not_implemented_or_skip(data, ndata, msg)
            if self.format_code == 2 and self.num_wide == 111:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
            else:
                raise RuntimeError(self.code_information())

        elif self.element_type == 4: # CSHEAR
            # 4-CSHEAR
            if self.is_stress():
                obj_vector_real = RealShearStressArray
                result_name = 'cshear_stress'
            else:
                obj_vector_real = RealShearStrainArray
                result_name = 'cshear_strain'

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)

            slot = getattr(self, result_name)
            if self.format_code == 2 and self.num_wide == 5: # random
                # checked
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

            if self.format_code == 2 and self.num_wide == 3: # random
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

            if self.format_code == 2 and self.num_wide == 19: # random
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


            #numwide_real = 4 + 21 * nnodes_expected
            #numwide_imag = 4 + (17 - 4) * nnodes_expected
            numwide_random = 4 + (11 - 4) * nnodes_expected
            numwide_random2 = 18 + 14 * (nnodes_expected - 1)
            preline1 = '%s-%s' % (self.element_name, self.element_type)
            preline2 = ' ' * len(preline1)

            #print('numwide real=%s imag=%s random=%s' % (numwide_real, numwide_imag, numwide_random2))
            self._data_factor = nnodes_expected

            if self.format_code == 1 and self.num_wide == numwide_random: # random
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
            elif self.format_code == 2 and self.num_wide == numwide_random2:
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
                num_wide_random = 4 + nnodes_expected * (17 - 4)

                #print('random2=%s' % num_wide_random)
                #print(self.code_information())

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
            if self.format_code == 1 and self.num_wide == 0: # random
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
            else:
                raise RuntimeError(self.code_information())

        elif self.element_type == 74:  # TRIA3
            if self.is_stress():
                obj_vector_real = RealPlateStressArray
                result_name = 'ctria3_stress'
            else:
                obj_vector_real = RealPlateStrainArray
                result_name = 'ctria3_strain'

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)

            slot = getattr(self, result_name)
            if self.format_code == 2 and self.num_wide == 17: # random; CTRIA3
                assert self.table_name in [b'OESVM1', b'OSTRVM1'], self.code_information()
                #OESVM
                #msg = self.code_information()
                msg = '%s-%s' % (self.table_name_str, self.element_name)
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
            if self.format_code == 1 and self.num_wide == numwide_random: # random
                msg = self.code_information()
                msg += '  numwide=%s numwide_real=%s numwide_imag=%s numwide_random=%s' % (
                    self.num_wide, numwide_real, numwide_imag, numwide_random)
                return self._not_implemented_or_skip(data, ndata, msg)
            elif self.format_code == 2 and self.num_wide == 87:
                # 87 - CQUAD4-144
                #msg = self.code_information()
                msg = '%s-CQUAD4-numwide=%s numwide_real=%s numwide_imag=%s numwide_random=%s' % (
                    self.table_name_str, self.num_wide, numwide_real, numwide_imag, numwide_random)
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

            if self.format_code == 1 and self.num_wide == 0: # random
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
            msg = '%s-COMP-random-numwide=%s numwide_real=11 numwide_imag=9' % (
                self.table_name_str, self.num_wide)
            msg = self.code_information()
            if self.format_code == 2 and self.num_wide == 13:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
            else:
                raise NotImplementedError(msg)

        #=========================
        elif self.element_type == 53: # axial plates - ctriax6
            # 53 - CTRIAX6
            if self.is_stress():
                result_name = 'ctriax_stress'
            else:
                result_name = 'ctriax_strain'
            self._results._found_result(result_name)
            slot = getattr(self, result_name)

            if self.format_code == 1 and self.num_wide == 0:
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

            if self.format_code == 1 and self.num_wide == 0:  # real
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

            if self.format_code == 1 and self.num_wide == 0:  # real
                pass
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
            if self.format_code == 1 and self.num_wide == 0:  # real?
                pass
            else:
                raise RuntimeError(self.code_information())

        elif self.element_type == 94:
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

            if self.format_code == 1 and self.num_wide == numwide_random:  # random
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

            numwide_random = 0
            if self.format_code == 1 and self.num_wide == numwide_random:  # random
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
                pass
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

            if self.format_code == -1 and self.num_wide == 0: # real
                pass
            else:
                #raise RuntimeError(self.code_information())
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type == 139:
            # 139-QUAD4FD
            raise NotImplementedError(self.code_information())
        elif self.element_type == 189:
            # 189-VUQUAD
            raise NotImplementedError(self.code_information())

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
            raise NotImplementedError(self.code_information())

        elif self.element_type in [276, 277, 278]:
            # 276-CPLSTS4
            # 277-CPLSTS6
            # 278-CPLSTS8
            msg = self.code_information()
            return self._not_implemented_or_skip(data, ndata, msg)

        assert ndata > 0, ndata
        assert nelements > 0, 'nelements=%r element_type=%s element_name=%r' % (nelements, self.element_type, self.element_name)
        #assert ndata % ntotal == 0, '%s n=%s nwide=%s len=%s ntotal=%s' % (self.element_name, ndata % ntotal, ndata % self.num_wide, ndata, ntotal)
        assert self.num_wide * 4 == ntotal, 'numwide*4=%s ntotal=%s' % (self.num_wide * 4, ntotal)
        assert self.thermal == 0, "thermal = %%s" % self.thermal
        assert n > 0, "n = %s result_name=%s" % (n, result_name)
        return n

    def oesrt_cquad4_95(self, data, ndata):
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
