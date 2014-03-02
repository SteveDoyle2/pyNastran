#pylint: disable=C0301,C0103,C0111,E1101
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys
from struct import unpack

from pyNastran import isRelease
from .real.elementsStressStrain import RealElementsStressStrain
from .real.oes_bush1d import Bush1DStressObject
from .real.oes_compositePlates import RealCompositePlateStress, RealCompositePlateStrain
from .real.oes_solids import RealSolidStress, RealSolidStrain


from .complex.elementsStressStrain import ComplexElementsStressStrain
from .complex.oes_bush import ComplexBushStress, ComplexBushStrain
from .complex.oes_bush1d import ComplexBush1DStressObject


from .oes_nonlinear import NonlinearQuadObject, HyperelasticQuad


class OES(RealElementsStressStrain, ComplexElementsStressStrain):
    """Table of stresses/strains"""

    def _delete_attributes_OES(self):
        params = ['s_code', 'element_type', 'obj', 'markerStart', 'load_set', 'format_code', 's_code', 'thermal',
                  'lsdvmn', 'mode', 'eign', 'mode_cycle', 'freq', 'mode', 'eigr', 'eigi', 'dt']
        self._delete_attributes(params)

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
        #print("numWideReal=%s numWideImag=%s numWideRandom=%s" % (numWideReal, numWideImag, numWideRandom))
        #print('element_type=%s' % self.element_type)
        #print('self.num_wide = %s' % self.num_wide)

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

        if 0:
            pass

        elif self.element_type in [40]:   # CBUSH1D
            if self.num_wide == numWideReal:
                result_name = self.makeOES_Object(self.bush1dStressStrain, Bush1DStressObject, 'bush1dStressStrain',
                                                 self.bush1dStressStrain, Bush1DStressObject, 'bush1dStressStrain')
                name = result_name + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_CBUSH1D_40, result_name, name)
            elif self.num_wide == numWideImag:
                result_name = self.makeOES_Object(self.bush1dStressStrain, ComplexBush1DStressObject, 'bush1dStressStrain',
                                                 self.bush1dStressStrain, ComplexBush1DStressObject, 'bush1dStressStrain')
                name = result_name + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_CBUSH1D_40_alt, result_name, name)
            else:
                self.not_implemented_or_skip()

        elif self.element_type in [88, 90]:  # CTRIA3NL, CQUAD4NL
            if self.num_wide == numWideReal:
                result_name = self.makeOES_Object(self.nonlinearPlateStress, NonlinearQuadObject, 'nonlinearPlateStress',
                                                 self.nonlinearPlateStrain, NonlinearQuadObject, 'nonlinearPlateStrain')
                name = result_name + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_CQUAD4NL_90, result_name, name)
            elif self.num_wide == numWideImag:  # TODO switch to ComplexNonlinearPlateObject
                result_name = self.makeOES_Object(self.nonlinearPlateStress, NonlinearQuadObject, 'nonlinearPlateStress',
                                                 self.nonlinearPlateStrain, NonlinearQuadObject, 'nonlinearPlateStrain')
                name = result_name + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_CQUAD4NL_90_alt, result_name, name)
            else:
                self.not_implemented_or_skip()

        elif self.element_type in [85, 91, 93]: # CTETRANL 85 / CPENTANL 91 / CHEXANL 93
            #print('not done...')
            result_name = self.makeOES_Object(self.solidStress, RealSolidStress, 'solidStress',
                                             self.solidStrain, RealSolidStrain, 'solidStrain')
            name = result_name + ': Subcase %s' % self.isubcase
            self.handle_results_buffer(self.OES_TETRANL_85_PENTANL_91_CHEXANL_93, result_name, name)

        elif self.element_type in [145, 146, 147]: # VUHEXA 145 / VUPENTA 146 / VUTETRA 147
            if isRelease:
                self.not_implemented_or_skip()
                return
            else:
                #print('only the first type read, not parsed...')
                result_name = self.makeOES_Object(self.solidStress, RealSolidStress, 'solidStress',
                                                 self.solidStrain, RealSolidStrain, 'solidStrain')
                name = result_name + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_VUHEXA_145_VUPENTA_146_VUTETRA_147, result_name, name)

        elif self.element_type in [95, 96, 97, 98]:  # CQUAD4, CQUAD8, CTRIA3, CTRIA6 (composite)
            self.eid2 = None  # stores the previous elementID
            if self.num_wide == numWideReal:
                result_name = self.makeOES_Object(self.compositePlateStress, RealCompositePlateStress, 'compositePlateStress',
                                                 self.compositePlateStrain, RealCompositePlateStrain, 'compositePlateStrain')
                name = result_name + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_CQUAD4_95, result_name, name)
            #elif self.num_wide == numWideImag:
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

        elif self.element_type in [139]:   # QUAD4FD (hyperelastic)
            if self.num_wide == numWideReal:
                result_name = self.makeOES_Object(self.hyperelasticPlateStress, HyperelasticQuad, 'hyperelasticPlateStress',
                                                 self.hyperelasticPlateStrain, HyperelasticQuad, 'hyperelasticPlateStrain')
                name = result_name + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_QUAD4FD_139, result_name, name)
            #elif self.num_wide == numWideImag:
            #    resultName = self.makeOES_Object(self.hyperelasticPlateStress,ComplexHyperelasticQuadObject,'hyperelasticPlateStress',
            #                                     self.hyperelasticPlateStrain,ComplexHyperelasticQuadObject,'hyperelasticPlateStrain')
            #    self.handle_results_buffer(self.OES_QUAD4FD_139,resultName)
            else:
                self.not_implemented_or_skip()

        elif self.element_type in [2189]:   # VUQUAD 189 ()
            #print "    found QUAD4FD_139"
            if self.num_wide == numWideReal:
                result_name = self.makeOES_Object(self.hyperelasticPlateStress, VUQuadObject, 'VUStress',
                                                 self.hyperelasticPlateStrain, VUQuadObject, 'VUStrain')
                name = result_name + ': Subcase %s' % self.isubcase
                self.handle_results_buffer(self.OES_QUAD4FD_139, result_name, name)
            #elif self.num_wide == numWideImag:
            #    resultName = self.makeOES_Object(self.hyperelasticPlateStress,ComplexHyperelasticQuadObject,'hyperelasticPlateStress',
            #                                     self.hyperelasticPlateStrain,ComplexHyperelasticQuadObject,'hyperelasticPlateStrain')
            #    self.handle_results_buffer(self.OES_QUAD4FD_139,resultName)
            else:
                self.not_implemented_or_skip()

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
        #    msg = '%s-OES format1_sort0 element_type=%-3s -> %s is not supported - fname=%s\n' % (self.table_name,self.element_type,self.get_element_type(self.element_type),self.op2_filename)
        #    raise AddNewElementError(msg)
        else:
            msg = '%s-OES format%s element_type=%-3s -> %s is not supported - fname=%s\n' % (self.table_name, self.format_code, self.element_type, self.get_element_type(self.element_type), self.op2_filename.strip())
            self.log.debug(msg)
            self.skippedCardsFile.write(msg)
            self.not_implemented_or_skip()
