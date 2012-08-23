#pylint: disable=C0301,C0111,C0324,R0912,R0915,W0223,E1101
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys
from struct import unpack

from .real.elementsStressStrain import RealElementsStressStrain
from .real.oes_rods import RodStressObject, RodStrainObject
from .real.oes_shear import ShearStressObject, ShearStrainObject
from .real.oes_bars import BarStressObject, BarStrainObject
from .real.oes_beams import BeamStressObject, BeamStrainObject
from .real.oes_solids import SolidStressObject, SolidStrainObject
from .real.oes_plates import PlateStressObject, PlateStrainObject
from .real.oes_springs import CelasStressObject, CelasStrainObject
from .real.oes_triax import TriaxStressObject, TriaxStrainObject
from .real.oes_compositePlates import CompositePlateStressObject, CompositePlateStrainObject


from .complex.elementsStressStrain import ComplexElementsStressStrain
from .complex.oes_rods import ComplexRodStressObject, ComplexRodStrainObject
from .complex.oes_springs import ComplexCelasStressObject, ComplexCelasStrainObject
from .complex.oes_bars import ComplexBarStressObject, ComplexBarStrainObject
from .complex.oes_plates import ComplexPlateStressObject, ComplexPlateStrainObject


from .oes_nonlinear import NonlinearRodObject, NonlinearQuadObject, HyperelasticQuadObject


class OES(RealElementsStressStrain, ComplexElementsStressStrain):
    """Table of stresses/strains"""

    def readTable_OES(self):
        table3 = self.readTable_OES_3
        table4Data = self.readTable_OES_4_Data
        self.readResultsTable(
            table3, table4Data, flag=1)  # flag=1 defines old style
        self.deleteAttributes_OES()

    def deleteAttributes_OES(self):
        params = ['sCode', 'elementType', 'obj', 'markerStart', 'loadSet', 'formatCode', 'sCode', 'thermal',
                  'lsdvmn', 'mode', 'eign', 'modeCycle', 'freq', 'mode', 'eigr', 'eigi', 'dt']
        self.deleteAttributes(params)

    def readTable_OES_3(self, iTable):
        #print "*iTable3 = ", iTable
        #if 0:
            #markers = self.readMarkers([0,2])
            #print "markers=%s" % (markers)
            #block = self.readBlock()
            #print "block = ",block
            #markers = self.readMarkers([-1,7])
            #print "markers=%s" % (markers)
            #print self.printSection(200)

        bufferWords = self.getBufferWords()

        data = self.getData(4)
        bufferSize, = unpack(b'i', data)
        if self.makeOp2Debug:
            self.op2Debug.write('bufferSize=|%s|\n' % (str(bufferSize)))

        data = self.getData(4 * 50)
        #self.printBlock(data)
        if self.makeOp2Debug:
            self.op2Debug.write('block3header\n')

        self.parseApproachCode(data)  # 3
        ## element type
        self.addDataParameter(data, 'elementType', 'i', 3, False)
        ## load set ID
        self.addDataParameter(data, 'loadSet', 'i', 8, False)
        ## format code
        self.addDataParameter(data, 'formatCode', 'i', 9, False)
        ## number of words per entry in record
        ## @note is this needed for this table ???
        self.addDataParameter(data, 'numWide', 'i', 10, False)
        ## stress/strain codes
        self.addDataParameter(data, 'sCode', 'i', 11, False)
        ## thermal flag; 1 for heat ransfer, 0 otherwise
        self.addDataParameter(data, 'thermal', 'i', 23, False)

        #print "loadset=%s formatCode=%s numWordsEntry=%s sCode=%s" % (self.loadSet,self.formatCode,self.numWide,self.sCode)
        #print "thermal(23)=%s elementType(3)=%s" % (self.thermal,self.elementType)

        ## assuming tCode=1
        if self.analysisCode == 1:   # statics / displacement / heat flux
            ## load set number
            self.addDataParameter(data, 'lsdvmn', 'i', 5,False)
            self.applyDataCodeValue('dataNames', ['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysisCode == 2:  # real eigenvalues
            ## mode number
            self.addDataParameter(data, 'mode', 'i', 5)
            ## real eigenvalue            self.addDataParameter(data, 'eign', 'f', 6, False)
            ## mode or cycle @todo confused on the type - F1???
            self.addDataParameter(data, 'modeCycle', 'f', 7, False)
            self.applyDataCodeValue('dataNames', ['mode', 'eigr', 'modeCycle'])
        #elif self.analysisCode==3: # differential stiffness
            #self.lsdvmn = self.getValues(data,'i',5) ## load set number
            #self.dataCode['lsdvmn'] = self.lsdvmn
        #elif self.analysisCode==4: # differential stiffness
        #    self.lsdvmn = self.getValues(data,'i',5) ## load set number
        elif self.analysisCode == 5:   # frequency
            ## frequency
            self.addDataParameter(data, 'freq', 'f', 5)
            self.applyDataCodeValue('dataNames', ['freq'])
        elif self.analysisCode == 6:  # transient
            ## time step
            self.addDataParameter(data, 'dt', 'f', 5)
            self.applyDataCodeValue('dataNames', ['dt'])
        elif self.analysisCode == 7:  # pre-buckling
            ## load set
            self.addDataParameter(data, 'lsdvmn', 'i', 5)
            self.applyDataCodeValue('dataNames', ['lsdvmn'])
        elif self.analysisCode == 8:  # post-buckling
            ## mode number
            self.addDataParameter(data, 'lsdvmn', 'i', 5)
            self.addDataParameter(data, 'eigr', 'f', 6, False)  # real eigenvalue
            self.applyDataCodeValue('dataNames', ['lsdvmn', 'eigr'])
        elif self.analysisCode == 9:  # complex eigenvalues
            ## mode number
            self.addDataParameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            self.addDataParameter(data, 'eigr', 'f', 6, False)
            ## imaginary eigenvalue            self.addDataParameter(data, 'eigi', 'f', 7, False)
            self.applyDataCodeValue('dataNames', ['mode', 'eigr', 'eigi'])
        elif self.analysisCode == 10:  # nonlinear statics
            ## load step
            self.addDataParameter(data, 'lftsfq', 'f', 5)
            self.applyDataCodeValue('dataNames', ['lftsfq'])
        elif self.analysisCode == 11:  # old geometric nonlinear statics
            ## load set number
            self.addDataParameter(data, 'lsdvmn', 'i', 5)
            self.applyDataCodeValue('dataNames', ['lsdvmn'])
        elif self.analysisCode == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## Time step ??? --> straight from DMAP
            self.addDataParameter(data, 'dt', 'f', 5)
            self.applyDataCodeValue('dataNames', ['dt'])
        else:
            raise RuntimeError('invalid analysisCode...analysisCode=%s' %
                               (self.analysisCode))
        # tCode=2
        #if self.analysisCode==2: # sort2
        #    self.lsdvmn = self.getValues(data,'i',5)

        self.readTitle()
        #print "n4 = ",self.n

    def parseStressCode(self):
        """
        sCode =  0 -> stressBits = [0,0,0,0,0]
        sCode =  1 -> stressBits = [0,0,0,0,1]
        sCode =  2 -> stressBits = [0,0,0,1,0]
        sCode =  3 -> stressBits = [0,0,0,1,1]
        etc.
        sCode = 32 -> stressBits = [1,1,1,1,1]

        stressBits[0] = 0 -> isMaxShear=True       isVonMises=False
        stressBits[0] = 1 -> isMaxShear=False      isVonMises=True

        stressBits[1] = 0 -> isStress=True         isStrain=False
        stressBits[2] = 0 -> isFiberCurvature=True isFiberDistance=False
        stressBits[3] = 0 -> duplicate of Bit[1] (stress/strain)
        stressBits[4] = 0 -> material coordinate system flag
        """
        bits = [0, 0, 0, 0, 0]

        sCode = self.sCode
        i = 4
        while sCode > 0:
            value = sCode % 2
            sCode = (sCode - value) // 2
            bits[i] = value
            #print "    *bit = ",value
            #print "    sCode = ",sCode
            i -= 1
        #bits.reverse()
        #print "stressBits = ",bits
        self.stressBits = bits
        self.dataCode['stressBits'] = self.stressBits

    def readTable_OES_4_Data(self, iTable):
        isTable4Done = False
        isBlockDone = False
        #print self.printSection(100)

        data = self.getData(16)
        #print self.printBlock(data) # on
        #print "16 block..."
        #self.printBlock(data)
        bufferWords, = unpack(b'i', data[4:8])
        #print "bufferWords = ",bufferWords
        if self.makeOp2Debug:
            self.op2Debug.write('bufferWords=|%s|\n' % (str(bufferWords)))

        #print "*********************"
        #bufferWords = self.getMarker() # 87 - buffer
        #print "OES4 bufferWords = ",bufferWords,bufferWords*4
        #self.verifyBufferSize(bufferWords)

        isBlockDone = not(bufferWords)
        #print "self.firstPass = ",self.firstPass

        ## table -4 is done, restarting table -3
        if self.isBufferDone:  # table is done when the buffer is done
            isTable4Done = True
            #print "exitA"
            return (isTable4Done, isBlockDone)
        if bufferWords == 0:
            #print "bufferWords 0 - done with Table4"
            isTable4Done = True
            #isBlockDone  = True
            self.printSection(40)
            #print "exitB"
            return (isTable4Done, isBlockDone)

        self.readElementTable()
        return isTable4Done, isBlockDone

    def readElementTable(self):
        #print "**self.readElementTable"
        #print "*elementType = ",self.elementType
        #print "op2.tell=%s n=%s" % (self.op2.tell(),self.n)

        self.rewind(4)
        self.data = self.readBlock()  # 348
        #print "len(self.data) = ",len(self.data)

        if self.makeOp2Debug:
            self.op2Debug.write('reading big data block\n')
        #print self.printBlock(self.data)

        #msg = 'elementType=%s -> %s' % (self.elementType,self.ElementType(self.elementType))
        self.parseStressCode()

        if not self.isValidSubcase():  # lets the user skip a certain subcase
            self.log.debug("***skipping table=%s iSubcase=%s" %
                           (self.tableName, self.iSubcase))
            print("***skipping table=%s iSubcase=%s" % (
                self.tableName, self.iSubcase))
            self.skipOES_Element()
        elif self.thermal == 0:
            # Stress / Strain
            self.dataCode['elementName'] = self.ElementType(self.elementType)
            if self.tableCode == 5 and self.isSort1():
                assert self.tableName in ['OES1', 'OES1X', 'OES1X1', 'OES1C', 'OESNLXR', 'OESNLXD', 'OESNL1X', 'OESCP', 'OESTRCP',
                                          'OSTR1X', 'OSTR1C'], '%s is not supported' % (self.tableName)
                self.readOES_Data()
            else:
                self.NotImplementedOrSkip(
                    'invalid atfsCode=%s' % (self.atfsCode))
        elif self.thermal == 1:
            self.OES_Thermal()
        else:
            raise RuntimeError('invalid thermal option...')

        #print self.obj

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
                       24: 3,          # CVISC   ## @todo:  ?????
                       33: 9,          # CQUAD4
                       34: None,       # CBAR   ## @todo: 19 stress; 10 for strain???
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

        Real = realMapper[self.elementType]
        Imag = imagMapper[self.elementType]
        Random = randomMapper[self.elementType]
        return (Real, Imag, Random)

    def readOES_Data(self):
        #msg = '%s-OES elementType=%-3s -> %-6s\n' % (self.tableName,self.elementType,self.ElementType(self.elementType))
        msg = ''
        
        readCase = True
        if self.iSubcase in self.expectedTimes and len(self.expectedTimes[self.iSubcase]) > 0:
            readCase = self.updateDtMap()

        if readCase == False:
            self.skipOES_Element()
            return

        #print 'self.elementType  = ',self.elementType
        (numWideReal, numWideImag, numWideRandom) = self.OES_StressStrainCode()
        if self.elementType in [1, 3, 10]:  # crod/ctube/conrod
            #if self.elementType==1:    self.dataCode['elementName'] = 'CROD'
            #if self.elementType==3:    self.dataCode['elementName'] = 'CTUBE'
            #if self.elementType==10:   self.dataCode['elementName'] = 'CONROD'
            if self.numWide == numWideReal:
                resultName = self.makeOES_Object(self.rodStress, RodStressObject, 'rodStress',
                                                 self.rodStrain, RodStrainObject, 'rodStrain')
                self.handleResultsBuffer3(self.OES_basicElement, resultName)
            elif self.numWide == numWideImag:
                resultName = self.makeOES_Object(self.rodStress, ComplexRodStressObject, 'rodStress',
                                                 self.rodStrain, ComplexRodStrainObject, 'rodStrain')
                self.handleResultsBuffer3(self.OES_Rod1_alt, resultName)
            else:
                self.NotImplementedOrSkip()

        elif self.elementType == 2:   # cbeam - these have the same numWide
            self.dataCode['elementName'] = 'CBEAM'
            #print self.codeInformation()
            #print "numWideReal=%s numWideImag=%s" % (numWideReal,numWideImag)
            if self.formatCode == 1:  # Real
            #if self.numWide==numWideReal:
                resultName = self.makeOES_Object(self.beamStress, BeamStressObject, 'beamStress',
                                                 self.beamStrain, BeamStrainObject, 'beamStrain')
                self.handleResultsBuffer3(self.OES_CBEAM_2, resultName)
            #elif self.formatCode in [2,3]: # Imag
            ##elif self.numWide==numWideImag:
            #    self.makeOES_Object(self.beamStress, ComplexBeamStressObject, 'beamStress',
            #                        self.beamStrain, ComplexBeamStrainObject, 'beamStrain')
            #    self.handleResultsBuffer3(self.OES_CBEAM_2_alt,resultName)
            else:
                self.NotImplementedOrSkip()

        elif self.elementType in [4]:  # cshear
            self.dataCode['elementName'] = 'CSHEAR'
            if self.numWide == numWideReal:
                resultName = self.makeOES_Object(self.shearStress, ShearStressObject, 'shearStres',
                                                 self.shearStrain, ShearStrainObject, 'shearStrain')
                self.handleResultsBuffer3(self.OES_basicElement, resultName)
            #elif self.numWide==numWideImag:
            #    resultName = self.makeOES_Object(self.shearStress, ComplexShearStressObject, 'shearStress',
            #                                     self.shearStrain, ComplexShearStrainObject, 'shearStrain')
            #    self.handleResultsBuffer3(self.OES_shear4_alt,resultName)
            else:
                self.NotImplementedOrSkip()
        elif self.elementType in [11, 12, 13]:   # celas1/celas2/celas3
            #print "    found celas2_12"
            #if   self.elementType==11: self.dataCode['elementName'] = 'CELAS1'
            #elif self.elementType==12: self.dataCode['elementName'] = 'CELAS2'
            #elif self.elementType==13: self.dataCode['elementName'] = 'CELAS3'
            #else:  raise Exception('not implemented error')

            if self.numWide == numWideReal:
                resultName = self.makeOES_Object(self.celasStress, CelasStressObject, 'celasStress',
                                                 self.celasStrain, CelasStrainObject, 'celasStrain')
                self.handleResultsBuffer3(self.OES_basicElement, resultName)
            elif self.numWide == numWideImag:
                resultName = self.makeOES_Object(self.celasStress, ComplexCelasStressObject, 'celasStress',
                                                 self.celasStrain, ComplexCelasStrainObject, 'celasStrain')
                self.handleResultsBuffer3(self.OES_Elas1_alt, resultName)
                #print self.obj
                #raise NotImplementedError('add printing to stress CELAS complex...')
            else:
                self.NotImplementedOrSkip()
        elif self.elementType == 34:   # cbar
            self.dataCode['elementName'] = 'CBAR'
            #print "numWideReal=%s numWideImag=%s" % (numWideReal,numWideImag)
            if self.numWide == numWideReal:
                resultName = self.makeOES_Object(self.barStress, BarStressObject, 'barStress',
                                                 self.barStrain, BarStrainObject, 'barStrain')
                self.handleResultsBuffer3(self.OES_CBAR_34, resultName)
            elif self.numWide == numWideImag:
                resultName = self.makeOES_Object(self.barStress, ComplexBarStressObject, 'barStress',
                                                 self.barStrain, ComplexBarStrainObject, 'barStrain')
                self.handleResultsBuffer3(self.OES_CBAR_34_alt, resultName)
            else:
                self.NotImplementedOrSkip()

        elif self.elementType == 33:  # cquad4_33
            if self.numWide == numWideReal:
                self.dataCode['elementName'] = 'CQUAD4'
                resultName = self.makeOES_Object(self.plateStress, PlateStressObject, 'plateStress',
                                                 self.plateStrain, PlateStrainObject, 'plateStrain')
                self.handleResultsBuffer3(self.OES_CQUAD4_33, resultName)
            elif self.numWide == numWideImag:
                resultName = self.makeOES_Object(self.plateStress, ComplexPlateStressObject, 'plateStress',
                                                 self.plateStrain, ComplexPlateStrainObject, 'plateStrain')
                self.handleResultsBuffer3(self.OES_CQUAD4_33_alt, resultName)
            else:
                self.NotImplementedOrSkip()

        elif self.elementType == 53:  # ctriax6
            self.dataCode['elementName'] = 'CTRIAX6'
            #print self.codeInformation()
            if self.numWide == numWideReal:
                resultName = self.makeOES_Object(self.ctriaxStress, TriaxStressObject, 'ctriaxStress',
                                                 self.ctriaxStrain, TriaxStrainObject, 'ctriaxStrain')
                self.handleResultsBuffer3(self.OES_CTRIAX6_53, resultName)
            #elif self.numWide==numWideImag:
            #    resultName = self.makeOES_Object(self.ctriaxStress,ComplexTriaxStressObject,'ctriaxStress',
            #                                     self.ctriaxStrain,ComplexTriaxStrainObject,'ctriaxStrain')
            #    self.handleResultsBuffer3(self.OES_CTRIAX6_53_alt,resultName)
            else:
                self.NotImplementedOrSkip()

        elif self.elementType == 74:  # ctria
            self.dataCode['elementName'] = 'CTRIA3'
            if self.numWide == numWideReal:
                resultName = self.makeOES_Object(self.plateStress, PlateStressObject, 'plateStress',
                                                 self.plateStrain, PlateStrainObject, 'plateStrain')
                self.handleResultsBuffer3(self.OES_CTRIA3_74, resultName)
            elif self.numWide == numWideImag:
                resultName = self.makeOES_Object(self.plateStress, ComplexPlateStressObject, 'plateStress',
                                                 self.plateStrain, ComplexPlateStrainObject, 'plateStrain')
                self.handleResultsBuffer3(self.OES_CTRIA3_74_alt, resultName)
            else:
                self.NotImplementedOrSkip()
        elif self.elementType in [64, 144, 70, 75, 82]:  # 64-cquad8/cquad4/70-ctriar/ctria6/cquadr
            #print "    found cquad_144"
            if     self.elementType == 64:
                self.dataCode['elementName'] = 'CQUAD8'
            elif   self.elementType == 144:
                self.dataCode['elementName'] = 'CQUAD4'
            elif   self.elementType == 70:
                self.dataCode['elementName'] = 'CTRIAR'
            elif   self.elementType == 75:
                self.dataCode['elementName'] = 'CTRIA6'
            elif   self.elementType == 82:
                self.dataCode['elementName'] = 'CQUADR'
            else:
                msg = 'card not implemented elementType=%s' % (
                    self.elementType)
                raise NotImplementedError(msg)
            if self.numWide == numWideReal:
                resultName = self.makeOES_Object(self.plateStress, PlateStressObject, 'plateStress',
                                                 self.plateStrain, PlateStrainObject, 'plateStrain')
                self.handleResultsBuffer3(self.OES_CQUAD4_144, resultName)
            elif self.numWide == numWideImag:
                resultName = self.makeOES_Object(self.plateStress, ComplexPlateStressObject, 'plateStress',
                                                 self.plateStrain, ComplexPlateStrainObject, 'plateStrain')
                self.handleResultsBuffer3(self.OES_CQUAD4_144_alt, resultName)
            else:
                self.NotImplementedOrSkip()

        elif self.elementType in [39, 67, 68]:   # ctetra/chexa/cpenta (linear)
            if self.numWide == numWideReal:
                resultName = self.makeOES_Object(self.solidStress, SolidStressObject, 'solidStress',
                                                 self.solidStrain, SolidStrainObject, 'solidStrain')
                self.handleResultsBuffer3(self.OES_CSOLID_67, resultName)
            #elif self.numWide==numWideImag:
            #    resultName = self.makeOES_Object(self.solidStress,ComplexSolidStressObject,'solidStress',
            #                                     self.solidStrain,ComplexSolidStrainObject,'solidStrain')
            #    self.handleResultsBuffer3(self.OES_CSOLID_67_alt,resultName)
            else:
                self.NotImplementedOrSkip()

        elif self.elementType in [85]:   # ctetra/chexa/cpenta (91,93)  (nonlinear)
            if self.numWide == numWideReal:
                resultName = self.makeOES_Object(self.solidStress, SolidStressObject, 'solidStress',
                                                 self.solidStrain, SolidStrainObject, 'solidStrain')
                self.handleResultsBuffer3(self.OES_CSOLID_85, resultName)
            #elif self.numWide==numWideImag:
            #    resultName = self.makeOES_Object(self.solidStress,ComplexSolidStressObject,'solidStress',
            #                                     self.solidStrain,ComplexSolidStrainObject,'solidStrain')
            #    self.handleResultsBuffer3(self.OES_CSOLID_85_alt,resultName)
            else:
                self.NotImplementedOrSkip()

        elif self.elementType in [87, 89, 92]:   # CTUBENL, RODNL, CONRODNL
            #print "    found RODNL_89"
            resultName = self.makeOES_Object(self.nonlinearRodStress, NonlinearRodObject, 'nonlinearRodStress',
                                             self.nonlinearRodStrain, NonlinearRodObject, 'nonlinearRodStrain')
            self.handleResultsBuffer3(self.OES_RODNL_89_92, resultName)

        elif self.elementType in [88, 90]:  # CTRIA3NL, CQUAD4NL
            #print "cquad4_90"
            if self.numWide == numWideReal:
                resultName = self.makeOES_Object(self.nonlinearPlateStress, NonlinearQuadObject, 'nonlinearPlateStress',
                                                 self.nonlinearPlateStrain, NonlinearQuadObject, 'nonlinearPlateStrain')
                self.handleResultsBuffer3(self.OES_CQUAD4NL_90, resultName)
            elif self.numWide == numWideImag:  # @todo switch to ComplexNonlinearPlateObject
                resultName = self.makeOES_Object(self.nonlinearPlateStress, NonlinearQuadObject, 'nonlinearPlateStress',
                                                 self.nonlinearPlateStrain, NonlinearQuadObject, 'nonlinearPlateStrain')
                self.handleResultsBuffer3(self.OES_CQUAD4NL_90_alt, resultName)
            else:
                self.NotImplementedOrSkip()

        #elif self.elementType in [91]: # CPENTANL
        #    #print "hexa_93"
        #    self.handleResultsBuffer3(self.OES_CPENTANL_91)
        #elif self.elementType in [93]: # CHEXANL
        #    #print "hexa_93"
        #    self.handleResultsBuffer3(self.OES_CHEXANL_93)

        elif self.elementType in [95, 96, 97, 98]:  # CQUAD4, CQUAD8, CTRIA3, CTRIA6 (composite)
            #print "    found a 95/96/97 or 98!"
            self.eid2 = None  # stores the previous elementID
            if self.numWide == numWideReal:
                resultName = self.makeOES_Object(self.compositePlateStress, CompositePlateStressObject, 'compositePlateStress',
                                                 self.compositePlateStrain, CompositePlateStrainObject, 'compositePlateStrain')
                self.handleResultsBuffer3(self.OES_CQUAD4_95, resultName)
            #elif self.numWide==numWideImag:
            #    resultName = self.makeOES_Object(self.compositePlateStress,ComplexCompositePlateStressObject,'compositePlateStress',
            #                                     self.compositePlateStrain,ComplexCompositePlateStrainObject,'compositePlateStrain')
            #    self.handleResultsBuffer3(self.OES_CQUAD4_95_alt,resultName)
            else:
                self.NotImplementedOrSkip()
            del self.eid2

        #elif self.elementType in [94]: # CBEAM (nonlinear)
            #print "    found a 94!"
            #self.eid2 = None # stores the previous elementID
            #self.makeOES_Object(self.beamStress,beamStressObject,
            #                    self.beamStrain,beamStrainObject)
            #self.handleResultsBuffer3(self.OES_CBEAM_94)
            #raise NotImplementedError('stoping at end of CBEAM_94')
            #del self.eid2

        elif self.elementType in [139]:   # QUAD4FD (hyperelastic)
            #print "    found QUAD4FD_139"
            if self.numWide == numWideReal:
                resultName = self.makeOES_Object(self.hyperelasticPlateStress, HyperelasticQuadObject, 'hyperelasticPlateStress',
                                                 self.hyperelasticPlateStrain, HyperelasticQuadObject, 'hyperelasticPlateStrain')
                self.handleResultsBuffer3(self.OES_QUAD4FD_139, resultName)
            #elif self.numWide==numWideImag:
            #    resultName = self.makeOES_Object(self.hyperelasticPlateStress,ComplexHyperelasticQuadObject,'hyperelasticPlateStress',
            #                                     self.hyperelasticPlateStrain,ComplexHyperelasticQuadObject,'hyperelasticPlateStrain')
            #    self.handleResultsBuffer3(self.OES_QUAD4FD_139,resultName)
            else:
                self.NotImplementedOrSkip()

        #elif self.elementType in [2,53,61,70,86,88,90,94,102,189,232,]:
            #elementType=53  -> TRIAX6  is not supported
            #elementType=61  -> DUM9    is not supported
            #elementType=70  -> TRIAR   is not supported
            #elementType=86  -> GAPNL   is not supported
            #elementType=88  -> TRIA3NL is not supported
            #elementType=90  -> QUAD4NL is not supported
            #elementType=94  -> BEAMNL  is not supported
            #elementType=102 -> BUSH    is not supported
            #elementType=189 -> VUQUAD  is not supported
            #elementType=232 -> QUADRLC is not supported
        #elif self.elementType in [100]:   # BARS
        #    self.makeOES_Object(self.barStress,barStressObject,
        #                        self.barStrain,barStrainObject)
        #    self.handleResultsBuffer3(self.OES_basicElement)
        #elif self.elementType in [75,89,90,92,93]:
        #    msg = '%s-OES format1_sort0 elementType=%-3s -> %s is not supported - fname=%s\n' % (self.tableName,self.elementType,self.ElementType(self.elementType),self.op2FileName)
        #    raise AddNewElementError(msg)
        else:
            msg = '%s-OES format%s elementType=%-3s -> %s is not supported - fname=%s\n' % (self.tableName, self.formatCode, self.elementType, self.ElementType(self.elementType), self.op2FileName.strip())
            self.log.debug(msg)
            self.skippedCardsFile.write(msg)
            self.NotImplementedOrSkip()
        ###
        #elif self.elementType == 1:    # crod     (done)
        #elif self.elementType == 2:    # cbeam    (done)
        #elif self.elementType == 3:    # ctube    (done)
        #elif self.elementType == 4:    # cshear   (done)
        #elif self.elementType == 10:   # conrod   (done)

        #elif self.elementType == 33:   # cquad4_33 (done)
        #elif self.elementType == 34:   # cbar      (done)
        #elif self.elementType == 35:   # cconeax
        #elif self.elementType == 38:   # cgap
        #elif self.elementType == 39:   # ctetra    (done)
        #elif self.elementType == 40:   # cbush1d

        #elif self.elementType == 47:   # caxif2
        #elif self.elementType == 48:   # caxif3
        #elif self.elementType == 49:   # caxif4
        #elif self.elementType == 50:   # cslot3
        #elif self.elementType == 51:   # cslot4
        #elif self.elementType == 53:   # ctriax6
        #elif self.elementType == 55:   # cdum3
        #elif self.elementType == 56:   # cdum4
        #elif self.elementType == 57:   # cdum5
        #elif self.elementType == 58:   # cdum6
        #elif self.elementType == 59:   # cdum7
        #elif self.elementType == 60:   # cdum8/crac2d
        #elif self.elementType == 64:   # cquad8   (done)
        #elif self.elementType == 67:   # chexa    (done)
        #elif self.elementType == 68:   # cpenta   (done)
        #elif self.elementType == 69:   # cbend
        #elif self.elementType == 70:   # ctriar
        #elif self.elementType == 74:   # ctria3  (done)
        #elif self.elementType == 75:   # ctria6  (done)
        #elif self.elementType == 82:   # cquadr  (done)
        #elif self.elementType == 95:   # composite CQUAD4    (done)
        #elif self.elementType == 96:   # composite CQUAD8    (done)
        #elif self.elementType == 97:   # composite CTRIA3    (done)
        #elif self.elementType == 98:   # composite CTRIA6    (done)
        #elif self.elementType == 100:  # cbar w/ cbarao or pload1

        #elif self.elementType == 102:  # cbush
        #elif self.elementType == 144:  # cquad_144 - corner stresses (done)

        # rods/bars/beams
        #elif self.elementType == 1:    # crod   (done)
        #elif self.elementType == 2:    # cbeam  (done)
        #elif self.elementType == 3:    # ctube  (done)
        #elif self.elementType == 10:   # conrod (done)
        #elif self.elementType == 34:   # cbar   (done)
        #springs
        #elif self.elementType == 11:   # celas1 (done)
        #elif self.elementType == 12:   # celas2 (done)
        #elif self.elementType == 13:   # celas3 (done)
        #plate
        #elif self.elementType == 33:   # cquad_33 (done)
        #elif self.elementType == 74:   # ctria3_74  - (done)
        #elif self.elementType == 75:   # ctria6_75  - (in progress)
        #elif self.elementType == 64:   # cquad8_64  - corner stresses (done)
        #elif self.elementType == 144:  # cquad4_144 - corner stresses (done)
        #solid (???)
        #elif self.elementType == 39:  # ctetra (done)
        #elif self.elementType == 67:  # chexa  (done)
        #elif self.elementType == 68:  # cpenta (done)
        # composite plate
        #elif self.elementType == 95: # CQUAD4 (done)
        #elif self.elementType == 96: # CQUAD8 (done)
        #elif self.elementType == 97: # CTRIA3 (done)
        #elif self.elementType == 98: # CTRIA6 (done)
        # nonlinear
        #elif self.elementType == 85:   # tetra  (nonlinear,not integrated)
        #elif self.elementType == 86:   # gap    (nonlinear)
        #elif self.elementType == 87:   # ctube  (nonlinear, done)
        #elif self.elementType == 88:   # tria3  (nonlinear) - same as quad4
        #elif self.elementType == 89:   # crod   (nonlinear, done)
        #elif self.elementType == 90:   # quad4  (nonlinear, done)
        #elif self.elementType == 91:   # cpenta (nonlinear)
        #elif self.elementType == 92:   # conrod (nonlinear, done)
        #elif self.elementType == 93:   # chexa  (nonlinear,not integrated)
        #elif self.elementType == 94:   # cbeam  (nonlinear)
        # acoustic
        #elif self.elementType == 76:   # chexa  (acoustic)
        #elif self.elementType == 77:   # cpenta (acoustic)
        #elif self.elementType == 78:   # ctetra (acoustic)
        #elif self.elementType == 101:  # caabsf (acoustic)
    def makeOES_Object(self, stress, stressObject, stressName,
                             strain, strainObject, strainName):
        """
        Creates a stress/strain object if necessary
        """
        if self.isStress():
            self.createTransientObject(stress, stressObject)
            return stressName
        else:
            self.createTransientObject(strain, strainObject)
            return strainName
        #print "loading",self.obj.__class__.__name__
        #return self.obj
