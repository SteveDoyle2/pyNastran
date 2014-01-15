#pylint: disable=C0301,C0111,C0324,R0912,R0915,W0223,E1101
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys
from struct import unpack

from pyNastran import isRelease
from .real.elementsStressStrain import RealElementsStressStrain
from .real.oes_bars import BarStressObject, BarStrainObject
from .real.oes_beams import BeamStressObject, BeamStrainObject
from .real.oes_bush import BushStressObject #, BushStrainObject
from .real.oes_bush1d import Bush1DStressObject
from .real.oes_compositePlates import CompositePlateStressObject, CompositePlateStrainObject
from .real.oes_gap import NonlinearGapStressObject
from .real.oes_plates import PlateStressObject, PlateStrainObject
from .real.oes_rods import RodStressObject, RodStrainObject
from .real.oes_shear import ShearStressObject, ShearStrainObject
from .real.oes_solids import SolidStressObject, SolidStrainObject
from .real.oes_springs import CelasStressObject, CelasStrainObject, NonlinearSpringStressObject
from .real.oes_triax import TriaxStressObject, TriaxStrainObject


from .complex.elementsStressStrain import ComplexElementsStressStrain
from .complex.oes_bars import ComplexBarStressObject, ComplexBarStrainObject
from .complex.oes_bush import ComplexBushStressObject, ComplexBushStrainObject
from .complex.oes_bush1d import ComplexBush1DStressObject
from .complex.oes_plates import ComplexPlateStressObject, ComplexPlateStrainObject
from .complex.oes_rods import ComplexRodStressObject, ComplexRodStrainObject
from .complex.oes_springs import ComplexCelasStressObject, ComplexCelasStrainObject


from .oes_nonlinear import NonlinearRodObject, NonlinearQuadObject, HyperelasticQuadObject


class OES(RealElementsStressStrain, ComplexElementsStressStrain):
    """Table of stresses/strains"""

    def __init__(self):
        self._oes_complex_stress_map = {
            1 : ['CROD', self.rodStress, self.OES_Rod1_alt, ComplexRodStressObject, 'rodStress'],
            3 : ['CTUBE', self.rodStress, self.OES_Rod1_alt, ComplexRodStressObject, 'rodStress'],
            10 : ['CONROD', self.rodStress, self.OES_Rod1_alt, ComplexRodStressObject, 'rodStress'],
            # CTUBE
            # CONROD
        }
        self._oes_complex_strain_map = {
            1 : ['CROD', self.rodStrain, self.OES_Rod1_alt, ComplexRodStrainObject, 'rodStrain'],
            3 : ['CTUBE', self.rodStrain, self.OES_Rod1_alt, ComplexRodStrainObject, 'rodStrain'],
            10 : ['CONROD', self.rodStrain, self.OES_Rod1_alt, ComplexRodStrainObject, 'rodStrain'],
        }

        self._oes_real_stress_map = {
            # element_name, place to save data, extraction method, object_to_create, resultName (for debug)

            # rod
            1 : ['CROD', self.rodStress, self.OES_basicElement, RodStressObject, 'rodStress'],
            3 : ['CTUBE', self.rodStress, self.OES_basicElement, RodStressObject, 'rodStress'],
            10 : ['CONROD', self.rodStress, self.OES_basicElement, RodStressObject, 'rodStress'],
            #3 : ['CTUBE', self.ctubeStress, self.OES_basicElement, CtubeStrainObject, 'ctubeStress'],
            #10 : ['CONROD', self.conrodStress, self.OES_basicElement, ConrodStressObject, 'conrodStress'],

            # beam
            2 : ['CBEAM', self.beamStress, self.OES_CBEAM_2, BeamStressObject, 'beamStress'],

            # shear
            4 : ['CSHEAR', self.shearStress, self.OES_basicElement, ShearStressObject, 'shearStress'],

            # springs
            11 : ['CELAS1', self.celasStress, self.OES_basicElement, CelasStressObject, 'celasStress'],
            12 : ['CELAS2', self.celasStress, self.OES_basicElement, CelasStressObject, 'celasStress'],
            13 : ['CELAS3', self.celasStress, self.OES_basicElement, CelasStressObject, 'celasStress'],

            # dampers

            # bar
            34 : ['CBAR', self.barStress, self.OES_CBAR_34, BarStressObject, 'barStress'],

            # triax
            53 : ['CTRIAX6', self.ctriaxStress, self.OES_CTRIAX6_53, TriaxStressObject, 'ctriaxStress',],

            # plate
            33 : ['CQUAD4', self.plateStress, self.OES_CQUAD4_33, PlateStressObject, 'plateStress'],
            74 : ['CTRIA3', self.plateStress, self.OES_CTRIA3_74, PlateStressObject, 'plateStress'],

            # plate bilinear
            64  : ['CQUAD8', self.plateStress, self.OES_CQUAD4_144, PlateStressObject, 'plateStress'],
            144 : ['CQUAD4', self.plateStress, self.OES_CQUAD4_144, PlateStressObject, 'plateStress'],
            70  : ['CTRIAR', self.plateStress, self.OES_CQUAD4_144, PlateStressObject, 'plateStress'],
            75  : ['CTRIA6', self.plateStress, self.OES_CQUAD4_144, PlateStressObject, 'plateStress'],
            82  : ['CQUADR', self.plateStress, self.OES_CQUAD4_144, PlateStressObject, 'plateStress'],

            # bush
            40 : ['CBUSH1D', self.bush1dStressStrain, self.OES_CBUSH1D_40, Bush1DStressObject, 'bush1dStressStrain'],

            # solid
            39 : ['CTETRA', self.solidStress, self.OES_CSOLID_39_67_68, SolidStressObject, 'solidStress',],
            67 : ['CHEXA',  self.solidStress, self.OES_CSOLID_39_67_68, SolidStressObject, 'solidStress',],
            68 : ['CPENTA', self.solidStress, self.OES_CSOLID_39_67_68, SolidStressObject, 'solidStress',],

            # gap
            86 : ['CGAPNL', self.nonlinearGapStress, self.OES_CGAPNL_86, NonlinearGapStressObject, 'nonlinearGapStress',],

            88 : ['CTRIA3NL', self.nonlinearPlateStress, self.OES_CQUAD4NL_90, NonlinearQuadObject, 'nonlinearPlateStress',],
            90 : ['CQUAD4NL', self.nonlinearPlateStress, self.OES_CQUAD4NL_90, NonlinearQuadObject, 'nonlinearPlateStress',],

            # solid nonlinear
            85 : ['CTETRANL', self.solidStress, self.OES_TETRANL_85_PENTANL_91_CHEXANL_93, SolidStressObject, 'solidStress',],
            91 : ['CPENTANL', self.solidStress, self.OES_TETRANL_85_PENTANL_91_CHEXANL_93, SolidStressObject, 'solidStress',],
            93 : ['CHEXANL', self.solidStress, self.OES_TETRANL_85_PENTANL_91_CHEXANL_93, SolidStressObject, 'solidStress',],
        }

        self._oes_real_strain_map = {
            # rod
            1 : ['CROD', self.rodStrain, self.OES_basicElement, RodStrainObject, 'rodStrain'],
            3 : ['CTUBE', self.rodStrain, self.OES_basicElement, RodStrainObject, 'rodStrain'],
            10 : ['CONROD', self.rodStrain, self.OES_basicElement, RodStrainObject, 'rodStrain'],
            #3 : ['CTUBE', self.ctubeStrain, self.OES_basicElement, CtubeStrainObject, 'ctubeStrain'],
            #10 : ['CONROD', self.conrodStrain, self.OES_basicElement, ConrodStrainObject, 'conrodStrain'],

            # beam
            2 : ['CBEAM', self.beamStrain, self.OES_CBEAM_2, BeamStrainObject, 'beamStrain'],

            # shear
            4 : ['CSHEAR', self.shearStrain, self.OES_basicElement, ShearStrainObject, 'shearStrain'],

            # springs
            11 : ['CELAS1', self.celasStrain, self.OES_basicElement, CelasStrainObject, 'celasStrain'],
            12 : ['CELAS2', self.celasStrain, self.OES_basicElement, CelasStrainObject, 'celasStrain'],
            13 : ['CELAS3', self.celasStrain, self.OES_basicElement, CelasStrainObject, 'celasStrain'],
            # dampers

            # bar
            34 : ['CBAR', self.barStress, self.OES_CBAR_34, BarStressObject, 'barStress'],

            # triax
            53 : ['CTRIAX6', self.ctriaxStrain, self.OES_CTRIAX6_53, TriaxStrainObject, 'ctriaxStrain',],

            # plate
            33 : ['CQUAD4', self.plateStress, self.OES_CQUAD4_33, PlateStressObject, 'plateStrain'],
            74 : ['CTRIA3', self.plateStrain, self.OES_CTRIA3_74, PlateStrainObject, 'plateStrain'],

            # plate bilinear
            64  : ['CQUAD8', self.plateStrain, self.OES_CQUAD4_144, PlateStrainObject, 'plateStrain'],
            144 : ['CQUAD4', self.plateStrain, self.OES_CQUAD4_144, PlateStrainObject, 'plateStrain'],
            70  : ['CTRIAR', self.plateStrain, self.OES_CQUAD4_144, PlateStrainObject, 'plateStrain'],
            75  : ['CTRIA6', self.plateStrain, self.OES_CQUAD4_144, PlateStrainObject, 'plateStrain'],
            82  : ['CQUADR', self.plateStrain, self.OES_CQUAD4_144, PlateStrainObject, 'plateStrain'],

            # composite plate
            95 :  ['CQUAD4', self.compositePlateStrain, self.OES_CQUAD4_95, CompositePlateStrainObject, 'compositePlateStrain'],
            96 :  ['CQUAD8', self.compositePlateStrain, self.OES_CQUAD4_95, CompositePlateStrainObject, 'compositePlateStrain'],
            97 :  ['CTRIA3', self.compositePlateStrain, self.OES_CQUAD4_95, CompositePlateStrainObject, 'compositePlateStrain'],
            98 :  ['CTRIA6', self.compositePlateStrain, self.OES_CQUAD4_95, CompositePlateStrainObject, 'compositePlateStrain'],

            # bush
            #40 : ['CBUSH1D', self.bush1dStressStrain, self.OES_CBUSH1D_40, Bush1DStressObject, 'bush1dStressStrain'],

            # solid
            39 : ['CTETRA', self.solidStrain, self.OES_CSOLID_39_67_68, SolidStrainObject, 'solidStrain',],
            67 : ['CHEXA',  self.solidStrain, self.OES_CSOLID_39_67_68, SolidStrainObject, 'solidStrain',],
            68 : ['CPENTA', self.solidStrain, self.OES_CSOLID_39_67_68, SolidStrainObject, 'solidStrain',],

            # gap
            86 : ['CGAPNL', self.nonlinearGapStress, self.OES_CGAPNL_86, NonlinearGapStressObject, 'nonlinearGapStress',],

            88 : ['CTRIA3NL', self.nonlinearPlateStrain, self.OES_CQUAD4NL_90, NonlinearQuadObject, 'nonlinearPlateStrain',],
            90 : ['CQUAD4NL', self.nonlinearPlateStrain, self.OES_CQUAD4NL_90, NonlinearQuadObject, 'nonlinearPlateStrain',],

            # solid nonlinear
            85 : ['CTETRANL', self.solidStrain, self.OES_TETRANL_85_PENTANL_91_CHEXANL_93, SolidStrainObject, 'solidStrain',],
            91 : ['CPENTANL', self.solidStrain, self.OES_TETRANL_85_PENTANL_91_CHEXANL_93, SolidStrainObject, 'solidStrain',],
            93 : ['CHEXANL', self.solidStrain, self.OES_TETRANL_85_PENTANL_91_CHEXANL_93, SolidStrainObject, 'solidStrain',],
        }


    def readTable_OES(self):
        table3 = self.readTable_OES_3
        table4Data = self.readTable_OES_4_Data
        self.read_results_table(
            table3, table4Data, flag=1)  # flag=1 defines old style
        self._delete_attributes_OES()

    def _delete_attributes_OES(self):
        params = ['s_code', 'element_type', 'obj', 'markerStart', 'load_set', 'format_code', 's_code', 'thermal',
                  'lsdvmn', 'mode', 'eign', 'mode_cycle', 'freq', 'mode', 'eigr', 'eigi', 'dt']
        self._delete_attributes(params)

    def readTable_OES_3(self, iTable):
        #print "*iTable3 = ", iTable
        #if 0:
            #markers = self.read_markers([0,2])
            #print "markers=%s" % (markers)
            #block = self.read_block()
            #print "block = ",block
            #markers = self.read_markers([-1,7])
            #print "markers=%s" % (markers)
            #print self.print_section(200)

        buffer_words = self.get_buffer_words()

        data = self.get_data(4)
        buffer_size, = unpack(b'i', data)
        if self.make_op2_debug:
            self.op2_debug.write('buffer_size=%r\n' % buffer_size)

        data = self.get_data(4 * 50)
        #print(self.print_block(data))
        if self.make_op2_debug:
            self.op2_debug.write('readTable_OES_3\n')

        self.parse_approach_code(data)  # 3
        ## element type
        self.add_data_parameter(data, 'element_type', 'i', 3, False)
        ## load set ID
        self.add_data_parameter(data, 'load_set', 'i', 8, False)
        ## format code
        self.add_data_parameter(data, 'format_code', 'i', 9, False)
        ## number of words per entry in record
        ## .. note:: is this needed for this table ???
        self.add_data_parameter(data, 'num_wide', 'i', 10, False)
        ## stress/strain codes
        self.add_data_parameter(data, 's_code', 'i', 11, False)
        ## thermal flag; 1 for heat ransfer, 0 otherwise
        self.add_data_parameter(data, 'thermal', 'i', 23, False)

        #print "loadset=%s format_code=%s numWordsEntry=%s s_code=%s" % (self.load_set,self.format_code,self.num_wide,self.s_code)
        #print "thermal(23)=%s element_type(3)=%s" % (self.thermal,self.element_type)

        ## assuming tCode=1
        if self.analysis_code == 1:   # statics / displacement / heat flux
            ## load set number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5, False)
            self.apply_data_code_value('dataNames', ['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysis_code == 2:  # real eigenvalues
            ## mode number
            self.add_data_parameter(data, 'mode', 'i', 5)
            ## real eigenvalue            self.add_data_parameter(data, 'eign', 'f', 6, False)
            ## mode or cycle TODO confused on the type - F1???
            self.add_data_parameter(data, 'mode_cycle', 'f', 7, False)
            self.apply_data_code_value('dataNames', ['mode', 'eigr', 'mode_cycle'])
        #elif self.analysis_code==3: # differential stiffness
            #self.lsdvmn = self.get_values(data,'i',5) ## load set number
            #self.data_code['lsdvmn'] = self.lsdvmn
        #elif self.analysis_code==4: # differential stiffness
        #    self.lsdvmn = self.get_values(data,'i',5) ## load set number
        elif self.analysis_code == 5:   # frequency
            ## frequency
            self.add_data_parameter(data, 'freq', 'f', 5)
            self.apply_data_code_value('dataNames', ['freq'])
        elif self.analysis_code == 6:  # transient
            ## time step
            self.add_data_parameter(data, 'dt', 'f', 5)
            self.apply_data_code_value('dataNames', ['dt'])
        elif self.analysis_code == 7:  # pre-buckling
            ## load set
            self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.apply_data_code_value('dataNames', ['lsdvmn'])
        elif self.analysis_code == 8:  # post-buckling
            ## mode number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.add_data_parameter(data, 'eigr', 'f', 6, False)  # real eigenvalue
            self.apply_data_code_value('dataNames', ['lsdvmn', 'eigr'])
        elif self.analysis_code == 9:  # complex eigenvalues
            ## mode number
            self.add_data_parameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            self.add_data_parameter(data, 'eigr', 'f', 6, False)
            ## imaginary eigenvalue            self.add_data_parameter(data, 'eigi', 'f', 7, False)
            self.apply_data_code_value('dataNames', ['mode', 'eigr', 'eigi'])
        elif self.analysis_code == 10:  # nonlinear statics
            ## load step
            self.add_data_parameter(data, 'lftsfq', 'f', 5)
            self.apply_data_code_value('dataNames', ['lftsfq'])
        elif self.analysis_code == 11:  # old geometric nonlinear statics
            ## load set number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.apply_data_code_value('dataNames', ['lsdvmn'])
        elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## Time step ??? --> straight from DMAP
            self.add_data_parameter(data, 'dt', 'f', 5)
            self.apply_data_code_value('dataNames', ['dt'])
        else:
            raise RuntimeError('invalid analysis_code...analysis_code=%s' %
                               self.analysis_code)
        # tCode=2
        #if self.analysis_code==2: # sort2
        #    self.lsdvmn = self.get_values(data,'i',5)

        if not self.is_sort1() and not isRelease:
            raise NotImplementedError('sort2...')

        self.read_title()
        #print "n4 = ",self.n

    def parse_stress_code(self):
        """
        s_code =  0 -> stress_bits = [0,0,0,0,0]
        s_code =  1 -> stress_bits = [0,0,0,0,1]
        s_code =  2 -> stress_bits = [0,0,0,1,0]
        s_code =  3 -> stress_bits = [0,0,0,1,1]
        etc.
        s_code = 32 -> stress_bits = [1,1,1,1,1]

        stress_bits[0] = 0 -> isMaxShear=True       isVonMises=False
        stress_bits[0] = 1 -> isMaxShear=False      isVonMises=True

        stress_bits[1] = 0 -> isStress=True         isStrain=False
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

    def readTable_OES_4_Data(self, iTable):
        isTable4Done = False
        isBlockDone = False
        #print self.print_section(100)

        data = self.get_data(16)
        buffer_words, = unpack(b'i', data[4:8])
        if self.make_op2_debug:
            self.op2_debug.write('buffer_words=|%s|\n' % str(buffer_words))

        isBlockDone = not(buffer_words)

        #: table -4 is done, restarting table -3
        if self.isBufferDone:  # table is done when the buffer is done
            isTable4Done = True
            return (isTable4Done, isBlockDone)
        if buffer_words == 0:
            isTable4Done = True
            self.print_section(40)
            return (isTable4Done, isBlockDone)

        self._read_element_table()
        return isTable4Done, isBlockDone

    def _read_element_table(self):
        #print "**self.readElementTable"
        #print "*element_type = ",self.element_type
        #print "op2.tell=%s n=%s" % (self.op2.tell(),self.n)

        self.rewind(4)
        self.data = self.read_block()  # 348
        if self.make_op2_debug:
            self.op2_debug.write('reading big data block\n')
        #print self.print_block(self.data)

        #msg = 'element_type=%s -> %s' % (self.element_type, self.get_element_type(self.element_type))
        self.parse_stress_code()

        if not self.is_valid_subcase():  # lets the user skip a certain subcase
            self.log.debug("***skipping table=%s isubcase=%s" % (self.table_name, self.isubcase))
            print("***skipping table=%s isubcase=%s" % (self.table_name, self.isubcase))
            self.skipOES_Element()
            return

        if not self.is_sort1():
            if isRelease:
                self.not_implemented_or_skip('skipping OES SORT2')
            else:
                raise NotImplementedError('SORT2')


        if self.table_name in ['OES1', 'OES1X', 'OES1X1', 'OES1C',
                              'OESNLXR', 'OESNLXD', 'OESNL1X', 'OESCP', 'OESTRCP',
                              'OSTR1X', 'OSTR1C']:
            if self.thermal == 0:
                # Stress / Strain
                self.data_code['element_name'] = self.get_element_type(self.element_type)
                if self.table_code == 5:
                    assert self.table_name in ['OES1', 'OES1X', 'OES1X1', 'OES1C', 'OESNLXR', 'OESNLXD', 'OESNL1X', 'OESCP', 'OESTRCP',
                                              'OSTR1X', 'OSTR1C'], '%s is not supported' % (self.table_name)
                    self.readOES_Data()
                else:
                    self.not_implemented_or_skip('bad OES table')
            elif self.thermal == 1:
                self.OES_Thermal()
            else:
                self.not_implemented_or_skip('bad OES-thermal table')
        elif self.table_name in ['OESRT']:
            if isRelease:
                self.not_implemented_or_skip('bad OESRT table')
            else:
                if self.table_code == 25:  # failure indices for layered composite elements
                    self.readOESRT_Data()
                else:
                    self.not_implemented_or_skip('bad OESRT table')
        else:
            self.not_implemented_or_skip('bad OES table')

    def OES_StressStrainCode(self):
        """
        Gets the numwide codes for the element to determine if
        the real / complex / random result should be found.
        The format and sort codes do not always give the right answer...
        """
        realMapper = {
            0: 17,         # GRID - OES1G
            1: 5,          # CROD
            2: 1 + (11 - 1) * 11,  # CBEAM
            3: 5,          # CTUBE
            4: 4,          # CSHEAR
            10: 5,          # CONROD
            11: 2,          # CELAS1
            12: 2,          # CELAS2
            13: 2,          # CELAS3
            #14:  2,          # CELAS4
            24: None,       # CVISC
            33: 17,         # CQUAD4
            34: 16,         # CBAR
            35: 18,         # CCONEAX
            38: None,       # CGAP
            39: 4 + (25 - 4) * 5,  # CTETRA
            40: 8,          # CBUSH1D
            47: 5,          # CAXIF2
            48: 10,         # CAXIF3
            50: 6,          # CSLOT3
            51: 7,          # CSLOT4
            53: 1 + (9 - 1) * 4,  # CTRIAX6
            60: 10,         # CRAC2D, CDUM8
            61: 10,         # CRAC3D, CDUM9
            64: 2 + (19 - 2) * 5,  # CQUAD8
            67: 4 + (25 - 4) * 9,  # CHEXA
            68: 4 + (25 - 4) * 7,  # CPENTA
            69: 1 + (11 - 1) * 2,  # CBEND
            70: 2 + (19 - 2) * 4,  # CTRIAR
            74: 17,         # CTRIA3
            75: 2 + (19 - 2) * 4,  # CTRIA6
            82: 2 + (19 - 2) * 5,  # CQUADR
            85: 2 + (18 - 2) * 5,  # Nonlinear CTETRA
            86: 11,         # Nonlinear CGAP
            87: 7,          # Nonlinear CTUBE
            88: 13,         # Nonlinear CTRIA
            89: 7,          # Nonlinear CROD
            90: 13,         # Nonlinear QUAD
            91: 4 + (25 - 4) * 7,  # Nonlinear CPENTA
            92: 7,          # Nonlinear CONROD
            93: 4 + (25 - 4) * 9,  # Nonlinear CHEXA
            94: 51,         # Nonlinear CBEAM
            95: 11,         # Composite QUAD4
            96: 11,         # Composite QUAD8
            97: 11,         # Composite CTRIA3
            98: 11,         # Composite CTRIA6
            145: 98,         # VUHEXA
            146: 74,         # VUPENTA
            100: 10,         # CBAR 100
            101: 3,          # CAABSF
            102: 7,          # CBUSH
            139: 2 + (9 - 2) * 4,  # hyperelastic QUAD4FD
            140: 2 + (22 - 2) * 8,  # hyperelastic HEXAFD
            144: 2 + (19 - 2) * 5,  # bilinear CQUAD4
            147: 50,         # VUTETRA
            160: 2 + (22 - 2) * 6,  # linear hyperelastic PENTAFD - 6 nodes
            161: 2 + (22 - 2) * 1,  # linear hyperelastic TETRAFD - 4 nodes
            162: 2 + (9 - 2) * 1,  # linear hyperelastic TRIAFD
            163: 2 + (22 - 2) * 27,  # linear hyperelastic HEXAFD - 20 nodes
            164: 2 + (9 - 2) * 9,  # linear hyperelastic QUADFD
            165: 2 + (22 - 2) * 21,  # linear hyperelastic PENTAFD - 15 nodes
            166: 2 + (22 - 2) * 5,  # linear hyperelastic TETRAFD - 10 nodes
            167: 2 + (9 - 2) * 3,  # linear hyperelastic TRIA6FD
            168: 2 + (9 - 2) * 1,  # linear hyperelastic TRIAXFD
            169: 2 + (9 - 2) * 3,  # linear hyperelastic TRIAX6FD
            170: 2 + (9 - 2) * 4,  # linear hyperelastic QUADXFD
            171: 2 + (9 - 2) * 9,  # linear hyperelastic QUAD9XFD
            172: 13,         # Nonlinear CQUADR
            173: 13,         # Nonlinear CTRIAR
            189: 6 + (23 - 6) * 4,  # VUQUAD
            190: 6 + (23 - 6) * 3,  # VUTRIA
            191: 6 + ((12 - 6) * 2 + 4) * 2,  # VUBEAM
            200: 8,          # WELD
            201: 2 + (13 - 2) * 4,  # nonlinear hyperelastic QUAD
            202: 2 + (17 - 2) * 8,  # nonlinear hyperelastic HEXA8FD
            204: 2 + (17 - 2) * 8,  # hyperelastic HEXAFD - 20 nodes
            205: 2 + (17 - 2) * 5,  # hyperelastic PENTAFD - 4 nodes
            211: 2 + (13 - 2) * 3,  # hyperelastic TRIAFD
            212: 2 + (13 - 2) * 1,  # hyperelastic TRIAX3FD
            213: 2 + (13 - 2) * 3,  # hyperelastic TRIAXFD
            214: 2 + (13 - 2) * 4,  # hyperelastic QUAD4XFD
            215: 2 + (13 - 2) * 9,  # hyperelastic QUADXFD
            216: 2 + (17 - 2) * 4,  # hyperelastic TETRAFD -  4 nodes
            217: 2 + (13 - 2) * 3,  # hyperelastic TRIA3FD -  3 nodes
            218: 2 + (17 - 2) * 8,  # hyperelastic HEXAFD  - 20 nodes
            219: 2 + (13 - 2) * 4,  # hyperelastic QUAD4FD -  4 nodes
            220: 2 + (17 - 2) * 8,  # hyperelastic PENTAFD - 15 nodes
            221: 2 + (17 - 2) * 4,  # hyperelastic TETRAFD -  4 nodes
            222: 2 + (13 - 2) * 3,  # hyperelastic TRIAXFD -  3 nodes
            223: 2 + (13 - 2) * 4,  # hyperelastic QUADXFD -  4 nodes
            224: 3,          # nonlinear CELAS1
            225: 3,          # nonlinear CELAS3
            226: 19,         # nonlinear CBUSH
            232: 11,         # composite CQUADR (same as CQUAD4)
            233: 11,         # composite CTRIAR (same as CQUAD4)
            235: 16,         # punch plate CQUADR
            236: 16,         # punch plate CTRIAR
        }

        imagMapper = {
            0: None,       # GRID - OES1G
            1: 5,          # CROD
            2: 1 + (11 - 1) * 11,  # CBEAM
            3: 5,          # CTUBE
            4: 5,          # CSHEAR
            10: 5,          # CONROD
            11: 3,          # CELAS1
            12: 3,          # CELAS2
            13: 3,          # CELAS3
            14: 3,          # CELAS4
            24: 5,          # CVISC
            33: 15,         # CQUAD4
            34: 19,         # CBAR
            35: None,       # CCONEAX
            38: None,       # CGAP
            39: 4 + (17 - 4) * 5,  # CTETRA
            40: 9,          # CBUSH1D
            47: 9,          # CAXIF2
            48: 19,         # CAXIF3
            50: 11,         # CSLOT3
            51: 13,         # CSLOT4
            53: 1 + (10 - 1) * 4,  # CTRIAX6
            60: None,       # CRAC2D, CDUM8
            61: None,       # CRAC3D, CDUM9
            64: 2 + (17 - 2) * 5,  # CQUAD8
            67: 4 + (17 - 4) * 9,  # CHEXA
            68: 4 + (17 - 4) * 7,  # CPENTA
            69: 1 + (11 - 1) * 2,  # CBEND
            70: 2 + (17 - 2) * 4,  # CTRIAR
            74: 15,         # CTRIA3
            75: 2 + (17 - 2) * 4,  # CTRIA6
            82: 2 + (17 - 2) * 5,  # CQUADR
            85: None,       # Nonlinear CTETRA
            86: None,       # Nonlinear CGAP
            87: None,       # Nonlinear CTUBE
            88: 25,         # Nonlinear CTRIA
            89: None,       # Nonlinear CROD
            90: 25,         # Nonlinear QUAD
            91: None,       # Nonlinear CPENTA
            92: None,       # Nonlinear CONROD
            93: None,       # Nonlinear CHEXA
            94: None,       # Nonlinear CBEAM
            95: None,       # Composite QUAD4
            96: None,       # Composite QUAD8
            97: None,       # Composite CTRIA3
            98: None,       # Composite CTRIA6
            100: 16,         # CBAR 100
            101: 4,          # CAABSF
            102: 13,         # CBUSH
            139: None,       # hyperelastic QUADFD
            140: None,       # hyperelastic HEXAFD
            144: 2 + (17 - 2) * 5,  # bilinear CQUAD4
            145: 58,         # VUHEXA
            146: 44,         # VUPENTA
            147: 30,         # VUTETRA
            160: None,       # linear hyperelastic PENTAFD - 6 nodes
            161: None,       # linear hyperelastic TETRAFD - 4 nodes
            162: None,       # linear hyperelastic TRIAFD
            163: None,       # linear hyperelastic HEXAFD - 20 nodes
            164: None,       # linear hyperelastic QUADFD
            165: None,       # linear hyperelastic PENTAFD - 15 nodes
            166: None,       # linear hyperelastic TETRAFD - 10 nodes
            167: None,       # linear hyperelastic TRIA6FD
            168: None,       # linear hyperelastic TRIAXFD
            169: None,       # linear hyperelastic TRIAX6FD
            170: None,       # linear hyperelastic QUADXFD
            171: None,       # linear hyperelastic QUAD9XFD
            172: 25,         # Nonlinear CQUADR
            173: 25,         # Nonlinear CTRIAR
            189: 6 + (33 - 6) * 4,  # VUQUAD
            190: 6 + (33 - 6) * 3,  # VUTRIA
            191: 6 + ((15 - 6) * 2 + 4) * 2,  # VUBEAM
            200: None,       # WELD
            201: None,       # nonlinear hyperelastic QUAD
            202: None,       # nonlinear hyperelastic HEXA8FD
            204: None,       # hyperelastic HEXAFD - 20 nodes
            205: None,       # hyperelastic PENTAFD - 4 nodes
            211: None,       # hyperelastic TRIAFD
            212: None,       # hyperelastic TRIAX3FD
            213: None,       # hyperelastic TRIAXFD
            214: None,       # hyperelastic QUAD4XFD
            215: None,       # hyperelastic QUADXFD
            216: None,       # hyperelastic TETRAFD -  4 nodes
            217: None,       # hyperelastic TRIA3FD -  3 nodes
            218: None,       # hyperelastic HEXAFD  - 20 nodes
            219: None,       # hyperelastic QUAD4FD -  4 nodes
            220: None,       # hyperelastic PENTAFD - 15 nodes
            221: None,       # hyperelastic TETRAFD -  4 nodes
            222: None,       # hyperelastic TRIAXFD -  3 nodes
            223: None,       # hyperelastic QUADXFD -  4 nodes
                       224: None,       # nonlinear CELAS1
                       225: None,       # nonlinear CELAS3
                       226: None,       # nonlinear CBUSH
                       232: None,       # composite CQUADR (same as CQUAD4)
                       233: None,       # composite CTRIAR (same as CQUAD4)
                       235: 14,         # punch plate CQUADR
                       236: 14,         # punch plate CTRIAR
                 }

        randomMapper = {
                       0: None,       # GRID - OES1G
                       1: 3,          # CROD
                       2: 1 + (7 - 1) * 11,  # CBEAM
                       3: 3,          # CTUBE
                       4: 3,          # CSHEAR
                       10: 3,          # CONROD
                       11: 2,          # CELAS1
                       12: 2,          # CELAS2
                       13: 2,          # CELAS3
                      #14:  2,          # CELAS4
                       24: 3,          # CVISC   ## TODO:  ?????
                       33: 9,          # CQUAD4
                       34: None,       # CBAR   ## TODO: 19 stress; 10 for strain???
                       35: None,       # CCONEAX
                       38: None,       # CGAP
                       39: 4 + (11 - 4) * 5,  # CTETRA
                       40: 2 + (19 - 2) * 5,  # CBUSH1D
                       47: 5,          # CAXIF2
                       48: 10,         # CAXIF3
                       50: 6,          # CSLOT3
                       51: 7,          # CSLOT4
                       53: 1 + (9 - 1) * 4,  # CTRIAX6
                       60: None,       # CRAC2D, CDUM8
                       61: None,       # CRAC3D, CDUM9
                       64: 2 + (11 - 2) * 5,  # CQUAD8
                       67: 4 + (11 - 4) * 9,  # CHEXA
                       68: 4 + (11 - 4) * 7,  # CPENTA
                       69: 1 + (7 - 1) * 2,  # CBEND
                       70: 2 + (11 - 2) * 4,  # CTRIAR
                       74: 9,          # CTRIA3
                       75: 2 + (11 - 2) * 4,  # CTRIA6
                       82: 2 + (11 - 2) * 5,  # CQUADR
                       85: None,       # Nonlinear CTETRA
                       86: None,       # Nonlinear CGAP
                       87: None,       # Nonlinear CTUBE
                       88: None,       # Nonlinear CTRIA
                       89: None,       # Nonlinear CROD
                       90: None,       # Nonlinear QUAD
                       91: None,       # Nonlinear CPENTA
                       92: None,       # Nonlinear CONROD
                       93: None,       # Nonlinear CHEXA
                       94: None,       # Nonlinear CBEAM
                       95: None,       # Composite QUAD4
                       96: None,       # Composite QUAD8
                       97: None,       # Composite CTRIA3
                       98: None,       # Composite CTRIA6
                       100: 9,          # CBAR 100
                       101: 3,          # CAABSF
                       102: 7,          # CBUSH
                       139: None,       # hyperelastic QUADFD
                       140: None,       # hyperelastic HEXAFD
                       144: 2 + (11 - 2) * 5,  # bilinear CQUAD4
                       145: 106,        # VUHEXA
                       146: 80,         # VUPENTA
                       147: 54,         # VUTETRA
                       160: None,       # linear hyperelastic PENTAFD - 6 nodes
                       161: None,       # linear hyperelastic TETRAFD - 4 nodes
                       162: None,       # linear hyperelastic TRIAFD
                       163: None,       # linear hyperelastic HEXAFD - 20 nodes
                       164: None,       # linear hyperelastic QUADFD
                       165: None,       # linear hyperelastic PENTAFD - 15 nodes
                       166: None,       # linear hyperelastic TETRAFD - 10 nodes
                       167: None,       # linear hyperelastic TRIA6FD
                       168: None,       # linear hyperelastic TRIAXFD
                       169: None,       # linear hyperelastic TRIAX6FD
                       170: None,       # linear hyperelastic QUADXFD
                       171: None,       # linear hyperelastic QUAD9XFD
                       172: None,       # Nonlinear CQUADR
                       173: None,       # Nonlinear CTRIAR
                       189: None,       # VUQUAD
                       190: None,       # VUTRIA
                       191: None,       # VUBEAM
                       200: None,       # WELD
                       201: None,       # nonlinear hyperelastic QUAD
                       202: None,       # nonlinear hyperelastic HEXA8FD
                       204: None,       # hyperelastic HEXAFD - 20 nodes
                       205: None,       # hyperelastic PENTAFD - 4 nodes
                       211: None,       # hyperelastic TRIAFD
                       212: None,       # hyperelastic TRIAX3FD
                       213: None,       # hyperelastic TRIAXFD
                       214: None,       # hyperelastic QUAD4XFD
                       215: None,       # hyperelastic QUADXFD
                       216: None,       # hyperelastic TETRAFD -  4 nodes
                       217: None,       # hyperelastic TRIA3FD -  3 nodes
                       218: None,       # hyperelastic HEXAFD  - 20 nodes
                       219: None,       # hyperelastic QUAD4FD -  4 nodes
                       220: None,       # hyperelastic PENTAFD - 15 nodes
                       221: None,       # hyperelastic TETRAFD -  4 nodes
                       222: None,       # hyperelastic TRIAXFD -  3 nodes
                       223: None,       # hyperelastic QUADXFD -  4 nodes
                       224: None,       # nonlinear CELAS1
                       225: None,       # nonlinear CELAS3
                       226: None,       # nonlinear CBUSH
                       232: None,       # composite CQUADR (same as CQUAD4)
                       233: None,       # composite CTRIAR (same as CQUAD4)
                       235: None,       # punch plate CQUADR
                       236: None,       # punch plate CTRIAR
                 }

        Real = realMapper[self.element_type]
        Imag = imagMapper[self.element_type]
        Random = randomMapper[self.element_type]
        return (Real, Imag, Random)

    def readOESRT_Data(self):
        #msg = '%s-OES element_type=%-3s -> %-6s\n' % (self.table_name,self.element_type,self.get_element_type(self.element_type))
        msg = ''
        #if self.analysis_code not in [1,6,10]:
            #raise InvalidATFSCodeError('self.atfsCode=%s' % (self.atfsCode))

        readCase = True
        if self.isubcase in self.expected_times and len(self.expected_times[self.isubcase]) > 0:
            readCase = self.update_dt_map()

        if readCase == False:
            self.skipOES_Element()
            return

        #print 'self.element_type  = ',self.element_type
        (numWideReal, numWideImag, numWideRandom) = self.OES_StressStrainCode()
        print("numWideReal=%s numWideImag=%s numWideRandom=%s" %(numWideReal, numWideImag, numWideRandom))
        print('element_type=%s' % self.element_type)
        print('self.num_wide = %s' % self.num_wide)

        if self.element_type == 95:
            self.OESRT_CQUAD4_95()
        asdf

    def readOES_Data(self):
        #msg = '%s-OES element_type=%-3s -> %-6s\n' % (self.table_name, self.element_type, self.get_element_type(self.element_type))
        msg = ''

        readCase = True
        if self.isubcase in self.expected_times and len(self.expected_times[self.isubcase]) > 0:
            readCase = self.update_dt_map()

        if readCase == False:
            self.skipOES_Element()
            return

        (numWideReal, numWideImag, numWideRandom) = self.OES_StressStrainCode()

        if self.num_wide == numWideReal:
            if self.isStress():
                if self.element_type in self._oes_real_stress_map:
                    element_name, slot, extract_method, class_obj, resultName = self._oes_real_stress_map[self.element_type]

                    self.create_transient_object(slot, class_obj)
                    self.log.debug('stress %s' % resultName)
                    name = resultName + ': Subcase %s' % self.isubcase
                    self.handle_results_buffer(extract_method, resultName, name)
                    return
                else:
                    msg = 'need to add element_type=%s-%s to the stress_mapper' % (
                        self.data_code['element_name'], self.element_type)
                    self.log.debug(msg)
                    raise NotImplementedError(msg)
            else: # strain
                if self.element_type in self._oes_real_strain_map:
                    element_name, slot, extract_method, class_obj, resultName = self._oes_real_strain_map[self.element_type]
                    self.log.debug('strain %s' % resultName)
                    self.create_transient_object(slot, class_obj)
                    name = resultName + ': Subcase %s' % self.isubcase
                    self.handle_results_buffer(extract_method, resultName, name)
                    return
                else:
                    msg = 'need to add element_type=%s-%s to the real_strain_mapper' % (
                        self.data_code['element_name'], self.element_type)
                    self.log.debug(msg)
                    raise NotImplementedError(msg)
        elif self.num_wide == numWideImag:
            if self.isStress():
                if self.element_type in self._oes_complex_stress_map:
                    element_name, slot, extract_method, class_obj, resultName = self._oes_complex_stress_map[self.element_type]
                    self.log.debug('strain %s' % resultName)
                    self.create_transient_object(slot, class_obj)
                    name = resultName + ': Subcase %s' % self.isubcase
                    self.handle_results_buffer(extract_method, resultName, name)
                else:
                    msg = 'need to add element_type=%s-%s to the complex_stress_mapper' % (
                        self.data_code['element_name'], self.element_type)
                    self.log.debug(msg)
                    raise NotImplementedError(msg)
            else: # strain
                if self.element_type in self._oes_complex_strain_map:
                    element_name, slot, extract_method, class_obj, resultName = self._oes_complex_strain_map[self.element_type]
                    self.log.debug('strain %s' % resultName)
                    self.create_transient_object(slot, class_obj)
                    name = resultName + ': Subcase %s' % self.isubcase
                    self.handle_results_buffer(extract_method, resultName, name)
                    return
                else:
                    msg = 'need to add element_type=%s-%s to the complex_stress_mapper' % (
                        self.data_code['element_name'], self.element_type)
                    self.log.debug(msg)
                    raise NotImplementedError(msg)
        else:
            #self.not_implemented_or_skip()
            pass


        if self.element_type in [1, 3, 10]:  # crod/ctube/conrod
            #if self.element_type==3:    self.data_code['element_name'] = 'CTUBE'
            #if self.element_type==10:   self.data_code['element_name'] = 'CONROD'
            if self.num_wide == numWideReal:
                resultName = self.makeOES_Object(self.rodStress, RodStressObject, 'rodStress',
                                                 self.rodStrain, RodStrainObject, 'rodStrain')
                name = resultName + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_basicElement, resultName, name)
            elif self.num_wide == numWideImag:
                resultName = self.makeOES_Object(self.rodStress, ComplexRodStressObject, 'rodStress',
                                                 self.rodStrain, ComplexRodStrainObject, 'rodStrain')
                name = resultName + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_Rod1_alt, resultName, name)
            else:
                self.not_implemented_or_skip()

        elif self.element_type == 2:   # cbeam - these have the same num_wide
            self.data_code['element_name'] = 'CBEAM'
            #print self.code_information()
            #print "numWideReal=%s numWideImag=%s" % (numWideReal,numWideImag)
            if self.format_code == 1:  # Real
            #if self.num_wide==numWideReal:
                resultName = self.makeOES_Object(self.beamStress, BeamStressObject, 'beamStress',
                                                 self.beamStrain, BeamStrainObject, 'beamStrain')
                name = resultName + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_CBEAM_2, resultName, name)
            #elif self.format_code in [2,3]: # Imag
            ##elif self.num_wide==numWideImag:
            #    self.makeOES_Object(self.beamStress, ComplexBeamStressObject, 'beamStress',
            #                        self.beamStrain, ComplexBeamStrainObject, 'beamStrain')
            #    self.handle_results_buffer(self.OES_CBEAM_2_alt,resultName)
            else:
                self.not_implemented_or_skip()

        elif self.element_type in [4]:  # cshear
            self.data_code['element_name'] = 'CSHEAR'
            if self.num_wide == numWideReal:
                resultName = self.makeOES_Object(self.shearStress, ShearStressObject, 'shearStres',
                                                 self.shearStrain, ShearStrainObject, 'shearStrain')
                name = resultName + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_basicElement, resultName, name)
            #elif self.num_wide==numWideImag:
            #    resultName = self.makeOES_Object(self.shearStress, ComplexShearStressObject, 'shearStress',
            #                                     self.shearStrain, ComplexShearStrainObject, 'shearStrain')
            #    self.handle_results_buffer(self.OES_shear4_alt,resultName)
            else:
                self.not_implemented_or_skip()
        elif self.element_type in [11, 12, 13]:   # celas1/celas2/celas3
            #print "    found celas2_12"
            #if   self.element_type==11: self.data_code['element_name'] = 'CELAS1'
            #elif self.element_type==12: self.data_code['element_name'] = 'CELAS2'
            #elif self.element_type==13: self.data_code['element_name'] = 'CELAS3'
            #else:  raise Exception('not implemented error')

            if self.num_wide == numWideReal:
                resultName = self.makeOES_Object(self.celasStress, CelasStressObject, 'celasStress',
                                                 self.celasStrain, CelasStrainObject, 'celasStrain')
                name = resultName + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_basicElement, resultName, name)
            elif self.num_wide == numWideImag:
                resultName = self.makeOES_Object(self.celasStress, ComplexCelasStressObject, 'celasStress',
                                                 self.celasStrain, ComplexCelasStrainObject, 'celasStrain')
                name = resultName + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_Elas1_alt, resultName, name)
                #print self.obj
                #raise NotImplementedError('add printing to stress CELAS complex...')
            else:
                self.not_implemented_or_skip()
        elif self.element_type == 34:   # cbar
            self.data_code['element_name'] = 'CBAR'
            #print "numWideReal=%s numWideImag=%s" % (numWideReal,numWideImag)
            if self.num_wide == numWideReal:
                resultName = self.makeOES_Object(self.barStress, BarStressObject, 'barStress',
                                                 self.barStrain, BarStrainObject, 'barStrain')
                name = resultName + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_CBAR_34, resultName, name)
            elif self.num_wide == numWideImag:
                resultName = self.makeOES_Object(self.barStress, ComplexBarStressObject, 'barStress',
                                                 self.barStrain, ComplexBarStrainObject, 'barStrain')
                name = resultName + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_CBAR_34_alt, resultName, name)
            else:
                self.not_implemented_or_skip()

        elif self.element_type == 33:  # cquad4_33
            if self.num_wide == numWideReal:
                self.data_code['element_name'] = 'CQUAD4'
                resultName = self.makeOES_Object(self.plateStress, PlateStressObject, 'plateStress',
                                                 self.plateStrain, PlateStrainObject, 'plateStrain')
                name = resultName + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_CQUAD4_33, resultName, name)
            elif self.num_wide == numWideImag:
                resultName = self.makeOES_Object(self.plateStress, ComplexPlateStressObject, 'plateStress',
                                                 self.plateStrain, ComplexPlateStrainObject, 'plateStrain')
                name = resultName + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_CQUAD4_33_alt, resultName, name)
            else:
                self.not_implemented_or_skip()

        elif self.element_type == 53:  # ctriax6
            self.data_code['element_name'] = 'CTRIAX6'
            #print self.code_information()
            if self.num_wide == numWideReal:
                resultName = self.makeOES_Object(self.ctriaxStress, TriaxStressObject, 'ctriaxStress',
                                                 self.ctriaxStrain, TriaxStrainObject, 'ctriaxStrain')
                name = resultName + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_CTRIAX6_53, resultName, name)
            #elif self.num_wide==numWideImag:
            #    resultName = self.makeOES_Object(self.ctriaxStress,ComplexTriaxStressObject,'ctriaxStress',
            #                                     self.ctriaxStrain,ComplexTriaxStrainObject,'ctriaxStrain')
            #    self.handle_results_buffer(self.OES_CTRIAX6_53_alt,resultName)
            else:
                self.not_implemented_or_skip()

        elif self.element_type == 74:  # ctria
            self.data_code['element_name'] = 'CTRIA3'
            if self.num_wide == numWideReal:
                resultName = self.makeOES_Object(self.plateStress, PlateStressObject, 'plateStress',
                                                 self.plateStrain, PlateStrainObject, 'plateStrain')
                name = resultName + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_CTRIA3_74, resultName, name)
            elif self.num_wide == numWideImag:
                resultName = self.makeOES_Object(self.plateStress, ComplexPlateStressObject, 'plateStress',
                                                 self.plateStrain, ComplexPlateStrainObject, 'plateStrain')
                name = resultName + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_CTRIA3_74_alt, resultName, name)
            else:
                self.not_implemented_or_skip()
        elif self.element_type in [64, 144, 70, 75, 82]:  # 64-cquad8/cquad4/70-ctriar/ctria6/cquadr
            if     self.element_type == 64:
                self.data_code['element_name'] = 'CQUAD8'
            elif   self.element_type == 144:
                self.data_code['element_name'] = 'CQUAD4'
            elif   self.element_type == 70:
                self.data_code['element_name'] = 'CTRIAR'
            elif   self.element_type == 75:
                self.data_code['element_name'] = 'CTRIA6'
            elif   self.element_type == 82:
                self.data_code['element_name'] = 'CQUADR'
            else:
                msg = 'card not implemented element_type=%s' % (
                    self.element_type)
                raise NotImplementedError(msg)
            if self.num_wide == numWideReal:
                resultName = self.makeOES_Object(self.plateStress, PlateStressObject, 'plateStress',
                                                 self.plateStrain, PlateStrainObject, 'plateStrain')
                name = resultName + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_CQUAD4_144, resultName, name)
            elif self.num_wide == numWideImag:
                resultName = self.makeOES_Object(self.plateStress, ComplexPlateStressObject, 'plateStress',
                                                 self.plateStrain, ComplexPlateStrainObject, 'plateStrain')
                name = resultName + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_CQUAD4_144_alt, resultName, name)
            else:
                self.not_implemented_or_skip()

        elif self.element_type in [40]:   # CBUSH1D
            if self.num_wide == numWideReal:
                resultName = self.makeOES_Object(self.bush1dStressStrain, Bush1DStressObject, 'bush1dStressStrain',
                                                 self.bush1dStressStrain, Bush1DStressObject, 'bush1dStressStrain')
                name = resultName + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_CBUSH1D_40, resultName, name)
            elif self.num_wide==numWideImag:
                resultName = self.makeOES_Object(self.bush1dStressStrain, ComplexBush1DStressObject, 'bush1dStressStrain',
                                                 self.bush1dStressStrain, ComplexBush1DStressObject, 'bush1dStressStrain')
                name = resultName + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_CBUSH1D_40_alt, resultName, name)
            else:
                self.not_implemented_or_skip()

        elif self.element_type in [39, 67, 68]:   # ctetra/chexa/cpenta (linear)
            if self.num_wide == numWideReal:
                resultName = self.makeOES_Object(self.solidStress, SolidStressObject, 'solidStress',
                                                 self.solidStrain, SolidStrainObject, 'solidStrain')
                name = resultName + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_CSOLID_39_67_68, resultName, name)
            #elif self.num_wide==numWideImag:
            #    resultName = self.makeOES_Object(self.solidStress,ComplexSolidStressObject,'solidStress',
            #                                     self.solidStrain,ComplexSolidStrainObject,'solidStrain')
            #    self.handle_results_buffer(self.OES_CSOLID_67_alt,resultName)
            else:
                self.not_implemented_or_skip()

        elif self.element_type == 86:   # CGAPNL 86
            resultName = self.makeOES_Object(self.nonlinearGapStress, NonlinearGapStressObject, 'nonlinearGapStress',
                                             self.nonlinearGapStress, NonlinearGapStressObject, 'nonlinearGapStress')
            name = resultName + ': Subcase %s' % self.isubcase
            self.handle_results_buffer(self.OES_CGAPNL_86, resultName, name)

        elif self.element_type in [87, 89, 92]:   # CTUBENL, RODNL, CONRODNL
            resultName = self.makeOES_Object(self.nonlinearRodStress, NonlinearRodObject, 'nonlinearRodStress',
                                             self.nonlinearRodStrain, NonlinearRodObject, 'nonlinearRodStrain')
            name = resultName + ': Subcase %s' % self.isubcase
            self.handle_results_buffer(self.OES_RODNL_89_92, resultName, name)

        elif self.element_type in [88, 90]:  # CTRIA3NL, CQUAD4NL
            if self.num_wide == numWideReal:
                resultName = self.makeOES_Object(self.nonlinearPlateStress, NonlinearQuadObject, 'nonlinearPlateStress',
                                                 self.nonlinearPlateStrain, NonlinearQuadObject, 'nonlinearPlateStrain')
                name = resultName + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_CQUAD4NL_90, resultName, name)
            elif self.num_wide == numWideImag:  # TODO switch to ComplexNonlinearPlateObject
                resultName = self.makeOES_Object(self.nonlinearPlateStress, NonlinearQuadObject, 'nonlinearPlateStress',
                                                 self.nonlinearPlateStrain, NonlinearQuadObject, 'nonlinearPlateStrain')
                name = resultName + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_CQUAD4NL_90_alt, resultName, name)
            else:
                self.not_implemented_or_skip()

        elif self.element_type in [85, 91, 93]: # CTETRANL 85 / CPENTANL 91 / CHEXANL 93
            #print('not done...')
            resultName = self.makeOES_Object(self.solidStress, SolidStressObject, 'solidStress',
                                             self.solidStrain, SolidStrainObject, 'solidStrain')
            name = resultName + ': Subcase %s' % self.isubcase
            self.handle_results_buffer(self.OES_TETRANL_85_PENTANL_91_CHEXANL_93, resultName, name)
        elif self.element_type in [145, 146, 147]: # VUHEXA 145 / VUPENTA 146 / VUTETRA 147
            if isRelease:
                self.not_implemented_or_skip()
                return
            else:
                #print('only the first type read, not parsed...')
                resultName = self.makeOES_Object(self.solidStress, SolidStressObject, 'solidStress',
                                                 self.solidStrain, SolidStrainObject, 'solidStrain')
                name = resultName + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_VUHEXA_145_VUPENTA_146_VUTETRA_147, resultName, name)

        elif self.element_type in [95, 96, 97, 98]:  # CQUAD4, CQUAD8, CTRIA3, CTRIA6 (composite)
            self.eid2 = None  # stores the previous elementID
            if self.num_wide == numWideReal:
                resultName = self.makeOES_Object(self.compositePlateStress, CompositePlateStressObject, 'compositePlateStress',
                                                 self.compositePlateStrain, CompositePlateStrainObject, 'compositePlateStrain')
                name = resultName + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_CQUAD4_95, resultName, name)
            #elif self.num_wide==numWideImag:
            #    resultName = self.makeOES_Object(self.compositePlateStress,ComplexCompositePlateStressObject,'compositePlateStress',
            #                                     self.compositePlateStrain,ComplexCompositePlateStrainObject,'compositePlateStrain')
            #    self.handle_results_buffer(self.OES_CQUAD4_95_alt,resultName)
            else:
                self.not_implemented_or_skip()
            del self.eid2

        #elif self.element_type in [94]: # CBEAM (nonlinear)
            #print "    found a 94!"
            #self.eid2 = None # stores the previous elementID
            #self.makeOES_Object(self.beamStress,beamStressObject,
            #                    self.beamStrain,beamStrainObject)
            #self.handle_results_buffer(self.OES_CBEAM_94)
            #raise NotImplementedError('stoping at end of CBEAM_94')
            #del self.eid2

        elif self.element_type in [102]:   # CBUSH
            if self.num_wide == numWideReal:
                resultName = self.makeOES_Object(self.bushStress, BushStressObject, 'bushStress',
                                                 self.bushStrain, BushStressObject, 'bushStrain')
                name = resultName + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_CBUSH_102, resultName, name)
            elif self.num_wide==numWideImag:
                resultName = self.makeOES_Object(self.bushStress, ComplexBushStressObject, 'bushStress',
                                                 self.bushStrain, ComplexBushStrainObject, 'bushStrain')
                name = resultName + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_CBUSH_102_alt, resultName, name)
            else:
                self.not_implemented_or_skip()

        elif self.element_type in [139]:   # QUAD4FD (hyperelastic)
            if self.num_wide == numWideReal:
                resultName = self.makeOES_Object(self.hyperelasticPlateStress, HyperelasticQuadObject, 'hyperelasticPlateStress',
                                                 self.hyperelasticPlateStrain, HyperelasticQuadObject, 'hyperelasticPlateStrain')
                name = resultName + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_QUAD4FD_139, resultName, name)
            #elif self.num_wide==numWideImag:
            #    resultName = self.makeOES_Object(self.hyperelasticPlateStress,ComplexHyperelasticQuadObject,'hyperelasticPlateStress',
            #                                     self.hyperelasticPlateStrain,ComplexHyperelasticQuadObject,'hyperelasticPlateStrain')
            #    self.handle_results_buffer(self.OES_QUAD4FD_139,resultName)
            else:
                self.not_implemented_or_skip()

        elif self.element_type in [2189]:   # VUQUAD 189 ()
            #print "    found QUAD4FD_139"
            if self.num_wide == numWideReal:
                resultName = self.makeOES_Object(self.hyperelasticPlateStress, VUQuadObject, 'VUStress',
                                                 self.hyperelasticPlateStrain, VUQuadObject, 'VUStrain')
                name = resultName + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_QUAD4FD_139, resultName, name)
            #elif self.num_wide==numWideImag:
            #    resultName = self.makeOES_Object(self.hyperelasticPlateStress,ComplexHyperelasticQuadObject,'hyperelasticPlateStress',
            #                                     self.hyperelasticPlateStrain,ComplexHyperelasticQuadObject,'hyperelasticPlateStrain')
            #    self.handle_results_buffer(self.OES_QUAD4FD_139,resultName)
            else:
                self.not_implemented_or_skip()

        elif self.element_type in [224, 225]:   # CELAS1
            if self.element_type == 224:
                self.data_code['element_name'] = 'CELAS1'
            elif self.element_type == 225:
                self.data_code['element_name'] = 'CELAS3'
            else:
                raise NotImplementedError(self.element_type)
            resultName = self.makeOES_Object(self.nonlinearSpringStress, NonlinearSpringStressObject, 'nonlinearSpringStress',
                                             self.nonlinearSpringStress, NonlinearSpringStressObject, 'nonlinearSpringStress')
            name = resultName + ': Subcase %s' % self.isubcase
            self.handle_results_buffer(self.OES_CELAS_224_225, resultName, name)

        #elif self.element_type in [2,53,61,70,86,88,90,94,102,189,232,]:
            #element_type=53  -> TRIAX6  is not supported
            #element_type=61  -> DUM9    is not supported
            #element_type=70  -> TRIAR   is not supported
            #element_type=86  -> GAPNL   is not supported
            #element_type=88  -> TRIA3NL is not supported
            #element_type=90  -> QUAD4NL is not supported
            #element_type=94  -> BEAMNL  is not supported
            #element_type=189 -> VUQUAD  is not supported
            #element_type=232 -> QUADRLC is not supported
        #elif self.element_type in [100]:   # BARS
        #    self.makeOES_Object(self.barStress,barStressObject,
        #                        self.barStrain,barStrainObject)
        #    self.handle_results_buffer(self.OES_basicElement)
        #elif self.element_type in [75,89,90,92,93]:
        #    msg = '%s-OES format1_sort0 element_type=%-3s -> %s is not supported - fname=%s\n' % (self.table_name,self.element_type,self.get_element_type(self.element_type),self.op2FileName)
        #    raise AddNewElementError(msg)
        else:
            msg = '%s-OES format%s element_type=%-3s -> %s is not supported - fname=%s\n' % (self.table_name, self.format_code, self.element_type, self.get_element_type(self.element_type), self.op2FileName.strip())
            self.log.debug(msg)
            self.skippedCardsFile.write(msg)
            self.not_implemented_or_skip()

        #elif self.element_type == 1:    # crod     (done)
        #elif self.element_type == 2:    # cbeam    (done)
        #elif self.element_type == 3:    # ctube    (done)
        #elif self.element_type == 4:    # cshear   (done)
        #elif self.element_type == 10:   # conrod   (done)

        #elif self.element_type == 33:   # cquad4_33 (done)
        #elif self.element_type == 34:   # cbar      (done)
        #elif self.element_type == 35:   # cconeax
        #elif self.element_type == 38:   # cgap
        #elif self.element_type == 39:   # ctetra    (done)
        #elif self.element_type == 40:   # cbush1d

        #elif self.element_type == 47:   # caxif2
        #elif self.element_type == 48:   # caxif3
        #elif self.element_type == 49:   # caxif4
        #elif self.element_type == 50:   # cslot3
        #elif self.element_type == 51:   # cslot4
        #elif self.element_type == 53:   # ctriax6
        #elif self.element_type == 55:   # cdum3
        #elif self.element_type == 56:   # cdum4
        #elif self.element_type == 57:   # cdum5
        #elif self.element_type == 58:   # cdum6
        #elif self.element_type == 59:   # cdum7
        #elif self.element_type == 60:   # cdum8/crac2d
        #elif self.element_type == 64:   # cquad8   (done)
        #elif self.element_type == 67:   # chexa    (done)
        #elif self.element_type == 68:   # cpenta   (done)
        #elif self.element_type == 69:   # cbend
        #elif self.element_type == 70:   # ctriar
        #elif self.element_type == 74:   # ctria3  (done)
        #elif self.element_type == 75:   # ctria6  (done)
        #elif self.element_type == 82:   # cquadr  (done)
        #elif self.element_type == 95:   # composite CQUAD4    (done)
        #elif self.element_type == 96:   # composite CQUAD8    (done)
        #elif self.element_type == 97:   # composite CTRIA3    (done)
        #elif self.element_type == 98:   # composite CTRIA6    (done)
        #elif self.element_type == 100:  # cbar w/ cbarao or pload1
        #elif self.element_type == 102:  # cbush (done)
        #elif self.element_type == 144:  # cquad_144 - corner stresses (done)

        # rods/bars/beams/bush
        #elif self.element_type == 1:    # crod   (done)
        #elif self.element_type == 2:    # cbeam  (done)
        #elif self.element_type == 3:    # ctube  (done)
        #elif self.element_type == 10:   # conrod (done)
        #elif self.element_type == 34:   # cbar   (done)
        #elif self.element_type == 102:  # cbush  (done)
        #springs
        #elif self.element_type == 11:   # celas1 (done)
        #elif self.element_type == 12:   # celas2 (done)
        #elif self.element_type == 13:   # celas3 (done)
        #plate
        #elif self.element_type == 33:   # cquad_33 (done)
        #elif self.element_type == 74:   # ctria3_74  - (done)
        #elif self.element_type == 75:   # ctria6_75  - (in progress)
        #elif self.element_type == 64:   # cquad8_64  - corner stresses (done)
        #elif self.element_type == 144:  # cquad4_144 - corner stresses (done)
        #solid (???)
        #elif self.element_type == 39:  # ctetra (done)
        #elif self.element_type == 67:  # chexa  (done)
        #elif self.element_type == 68:  # cpenta (done)
        # composite plate
        #elif self.element_type == 95: # CQUAD4 (done)
        #elif self.element_type == 96: # CQUAD8 (done)
        #elif self.element_type == 97: # CTRIA3 (done)
        #elif self.element_type == 98: # CTRIA6 (done)
        # nonlinear
        #elif self.element_type == 85:   # tetra  (nonlinear,not integrated)
        #elif self.element_type == 86:   # gap    (nonlinear)
        #elif self.element_type == 87:   # ctube  (nonlinear, done)
        #elif self.element_type == 88:   # tria3  (nonlinear) - same as quad4
        #elif self.element_type == 89:   # crod   (nonlinear, done)
        #elif self.element_type == 90:   # quad4  (nonlinear, done)
        #elif self.element_type == 91:   # cpenta (nonlinear)
        #elif self.element_type == 92:   # conrod (nonlinear, done)
        #elif self.element_type == 93:   # chexa  (nonlinear,not integrated)
        #elif self.element_type == 94:   # cbeam  (nonlinear)
        # acoustic
        #elif self.element_type == 76:   # chexa  (acoustic)
        #elif self.element_type == 77:   # cpenta (acoustic)
        #elif self.element_type == 78:   # ctetra (acoustic)
        #elif self.element_type == 101:  # caabsf (acoustic)

    def makeOES_Object(self, stress, stressObject, stressName,
                             strain, strainObject, strainName):
        """
        Creates a stress/strain object if necessary
        """
        if self.isStress():
            self.create_transient_object(stress, stressObject)
            return stressName
        else:
            self.create_transient_object(strain, strainObject)
            return strainName
        #print "loading",self.obj.__class__.__name__
        #return self.obj
