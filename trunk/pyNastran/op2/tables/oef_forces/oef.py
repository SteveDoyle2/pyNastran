from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys
from struct import unpack


from .realForces import RealForces
from .complexForces import ComplexForces


from .oef_forceObjects import (RealRodForce, RealCBeamForce, RealCShearForce,
                               RealSpringForce, RealDamperForce, RealViscForce,
                               RealPlateForce, RealConeAxForce, RealPLATE2Force,
                               RealCBAR100Force, RealCGAPForce, RealBendForce,
                               RealPentaPressureForce, RealCBUSHForce,
                               RealForce_VU_2D, RealCBARForce, RealForce_VU)
from .oef_complexForceObjects import (ComplexRodForce, ComplexCBeamForce,
                                      ComplexCShearForce, ComplexSpringForce,
                                      ComplexDamperForce, ComplexViscForce,
                                      ComplexPlateForce, ComplexPLATE2Force,
                                      ComplexBendForce,
                                      ComplexPentaPressureForce,
                                      ComplexCBUSHForce, ComplexForce_VU_2D,
                                      ComplexCBARForce, ComplexForce_VU)
from .thermal_elements import ThermalElements


class OEF(ThermalElements, RealForces, ComplexForces):
    """Table of element forces"""
    def readTable_OEF(self):
        table3 = self.readTable_OEF_3
        table4Data = self.readOEF_Data
        self.readResultsTable(table3, table4Data)
        self.deleteAttributes_OEF()

    def deleteAttributes_OEF(self):
        params = ['elementType', 'dLoadID', 'loadID', 'obj', 'markerStart', 'oCode',
                  'eigr', 'eigi', 'eign', 'mode', 'freq', 'time', 'thermal', ]
        self.deleteAttributes(params)
        #print self.obj

    def readTable_OEF_3(self, iTable):  # iTable=-3
        bufferWords = self.getMarker()
        #print "2-bufferWords = ",bufferWords,bufferWords*4,'\n'

        data = self.getData(4)
        bufferSize, = unpack('i', data)
        data = self.getData(4 * 50)
        #self.printBlock(data)

        aCode = self.getBlockIntEntry(data, 1)

        self.parseApproachCode(data)
        ## element type
        self.addDataParameter( data, 'elementType', 'i', 3, False)  
        ## dynamic load set ID/random code        self.addDataParameter(data, 'dLoadID', 'i', 8, False)
        ## format code
        self.addDataParameter( data, 'formatCode', 'i', 9, False)
        ## number of words per entry in record
        ## @note is this needed for this table ???
        self.addDataParameter(data, 'numWide', 'i', 10, False)
        ## undefined in DMAP...
        self.addDataParameter(data, 'oCode', 'i', 11, False)
        ## thermal flag; 1 for heat ransfer, 0 otherwise
        self.addDataParameter(data, 'thermal', 'i', 23, False)

        #print "dLoadID(8)=%s formatCode(9)=%s numwde(10)=%s oCode(11)=%s thermal(23)=%s" %(self.dLoadID,self.formatCode,self.numWide,self.oCode,self.thermal)
        #print "thermal(23)=%s elementType(3)=%s" %(self.thermal,self.elementType)

        #if not self.isSort1():
            #raise NotImplementedError('sort2...')

        ## assuming tCode=1
        if self.analysisCode == 1:   # statics
            self.addDataParameter(data, 'loadID', 'i', 5,
                                  False)  # load set ID number
            #self.applyDataCodeValue('dataNames',['lsdvmn'])
            self.applyDataCodeValue('dataNames', ['loadID'])
            self.setNullNonlinearFactor()
        elif self.analysisCode == 2:  # normal modes/buckling (real eigenvalues)
            ## mode number
            self.addDataParameter(data, 'mode', 'i', 5)
            ## eigenvalue
            self.addDataParameter(data, 'eign', 'f', 6, False)
            self.applyDataCodeValue('dataNames', ['mode', 'eigr', 'modeCycle'])
            #print "mode(5)=%s eigr(6)=%s modeCycle(7)=%s" %(self.mode,self.eigr,self.modeCycle)
        elif self.analysisCode == 3:  # differential stiffness 0
            ## load set ID number
            self.addDataParameter(data, 'loadID', 'i', 5)
            self.applyDataCodeValue('dataNames', ['loadID'])
        elif self.analysisCode == 4:  # differential stiffness 1
            ## load set ID number
            self.addDataParameter(data, 'loadID', 'i', 5)
            self.applyDataCodeValue('dataNames', ['loadID'])
        elif self.analysisCode == 5:   # frequency
            self.addDataParameter(data, 'freq', 'f', 5)  # frequency
            self.applyDataCodeValue('dataNames', ['freq'])
        elif self.analysisCode == 6:  # transient
            self.addDataParameter(data, 'time', 'f', 5)  # time step
            self.applyDataCodeValue('dataNames', ['time'])
        elif self.analysisCode == 7:  # pre-buckling
            ## load set ID number
            self.addDataParameter(data, 'loadID', 'i', 5)
            #self.applyDataCodeValue('dataNames',['lsdvmn'])
            self.applyDataCodeValue('dataNames', ['loadID'])
        elif self.analysisCode == 8:  # post-buckling
            ## load set ID number
            self.addDataParameter(data, 'loadID', 'i', 5)
            ## real eigenvalue
            self.addDataParameter( data, 'eigr', 'f', 6, False)
            self.applyDataCodeValue('dataNames', ['lsdvmn', 'eigr'])
            #print "loadID(5)=%s  eigr(6)=%s" %(self.loadID,self.eigr)
        elif self.analysisCode == 9:  # complex eigenvalues
            ## mode number
            self.addDataParameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            self.addDataParameter(data, 'eigr', 'f', 6, False)
            ## imaginary eigenvalue
            self.addDataParameter(data, 'eigi', 'f', 7, False)
            self.applyDataCodeValue('dataNames', ['mode', 'eigr', 'eigi'])
        elif self.analysisCode == 10:  # nonlinear statics
            ## load step
            self.addDataParameter(data, 'loadStep', 'f', 5)
            self.applyDataCodeValue('dataNames', ['loadStep'])
        elif self.analysisCode == 11:  # geometric nonlinear statics
            ## load set ID number
            self.addDataParameter(data, 'loadID', 'i', 5)
            self.applyDataCodeValue('dataNames', ['loadID'])
            #print "loadID(5)=%s" %(self.loadID)
        else:
            raise RuntimeError('invalid analysisCode...analysisCode=%s' % (str(self.analysisCode) + '\n' + self.codeInformation()))

        # tCode=2
        #if self.analysisCode==2: # sort2
        #    self.loadID = self.getValues(data,'i',5) ## load set ID number

        #print "*iSubcase=%s"%(self.iSubcase)
        #print "analysisCode=%s tableCode=%s thermal=%s" %(self.analysisCode,self.tableCode,self.thermal)
        #print self.codeInformation()

        #self.printBlock(data)
        #print '-'*80
        self.readTitle()

    def OEF_ForceCode(self):
        """
        Gets the numwide codes for the element to determine if
        the real or complex result should be found.
        The format and sort codes do not always give the right answer...
        """
        realMapper = {
            1: 3,    # CROD
            2: 1 + (10 - 1) * 11,  # CBEAM
            3: 3,    # CTUBE
            4: 17,   # CSHEAR
            10: 3,    # CONROD
            11: 2,    # CELAS1
            12: 2,    # CELAS2
            13: 2,    # CELAS3
            14: 2,    # CELAS4

            20: 2,    # CDAMP1
            21: 2,    # CDAMP2
            22: 2,    # CDAMP3
            23: 2,    # CDAMP4
            24: 3,    # CVISC
            33: 9,    # CQUAD4
            34: 9,    # CBAR
            35: 7,    # CCONEAX
            38: 9,    # CGAP
            40: 8,    # CBUSH1D ???
            64: 2 + (11 - 2) * 5,  # CQUAD8
            69: 1 + (8 - 1) * 2,  # CBEND
            70: 2 + (11 - 2) * 4,  # CTRIAR
            74: 9,    # CTRIA3
            75: 2 + (11 - 2) * 4,  # CTRIA6


            #76:  16,   # Acoustic Velocity/Pressure CHEXA ???
            76: None,  # dummy so it doesnt go into the real results
            77: 10,   # Acoustic Velocity/Pressure CPENTA
            78: 10,   # Acoustic Velocity/Pressure CTETRA

            82: 2 + (11 - 2) * 5,  # CQUADR
            95: 9,    # composite CQUAD4 ???
            96: 9,    # composite CQUAD8 ???
            97: 9,    # composite CTRIA3 ???
            98: 9,    # composite CTRIA6 ???
            100: 8,    # BARS
            102: 7,    # CBUSH
            144: 2 + (11 - 2) * 5,  # bilinear CQUAD4
            189: 6 + (19 - 6) * 4,  # VUQUAD
            190: 6 + (19 - 6) * 3,  # VUTRIA
            191: 4 + (12 - 4) * 2,  # VUBEAM
            200: 9,    # CWELD
            232: 9,    # composite CQUADR ???
            233: 9,    # composite TRIAR ???
            235: 9,    # punch CQUADR...numWide in DMAP is wrong...left out first entry...
            236: 8,    # punch CTRIAR
        }
        imagMapper = {
            1: 5,    # CROD
            2: 1 + (17 - 1) * 11,  # CBEAM
            3: 5,    # CTUBE
            4: 33,   # CSHEAR
            10: 5,    # CONROD

            11: 3,    # CELAS1
            12: 3,    # CELAS2
            13: 3,    # CELAS3
            14: 3,    # CELAS4

            20: 3,    # CDAMP1
            21: 3,    # CDAMP2
            22: 3,    # CDAMP3
            23: 3,    # CDAMP4
            24: 5,    # CVISC
            33: 17,   # CQUAD4
            34: 17,   # CBAR
            35: 7,    # CCONEAX # needed to not crash the code...
            38: 9,    # CGAP
            40: 8,    # CBUSH1D ???
            64: 2 + (19 - 2) * 5,  # CQUAD8
            69: 1 + (14 - 1) * 2,  # CBEND
            70: 2 + (19 - 2) * 4,  # CTRIAR
            74: 17,   # CTRIA3
            75: 2 + (19 - 2) * 4,  # CTRIA6

            76: 16,   # Acoustic Velocity/Pressure CHEXA_PR
            77: 16,   # Acoustic Velocity/Pressure CPENTA_PR
            78: 16,   # Acoustic Velocity/Pressure CTETRA_PR

            82: 2 + (19 - 2) * 5,  # CQUADR
            95: 9,    # composite CQUAD4 ???
            96: 9,    # composite CQUAD8 ???
            97: 9,    # composite CTRIA3 ???
            98: 9,    # composite CTRIA6 ???
            100: 14,   # BARS
            102: 13,   # CBUSH
            144: 2 + (19 - 2) * 5,  # bilinear CQUAD4
            189: 6 + (31 - 6) * 4,  # VUQUAD
            190: 6 + (31 - 6) * 3,  # VUTRIA
            191: 4 + (18 - 4) * 2,  # VUBEAM
            200: 17,   # CWELD
            232: 9,    # composite CQUADR ???
            233: 9,    # composite TRIAR ???
            235: 17,   # punch CQUADR...numWide in DMAP is wrong...left out first entry...
            236: 16,   # punch CTRIAR
        }

        Real = realMapper[self.elementType]
        Imag = imagMapper[self.elementType]
        return (Real, Imag)

    def readOEF_Data(self):
        """
        OEF1X -
        DOEF1 -
        """
        if self.tableCode == 4 and self.tableName in ['OEF1X', 'DOEF1']:  # Forces/Heat Flux
            assert self.tableName in ['OEF1X', 'DOEF1'], 'tableName=%s tableCode=%s' % (self.tableName, self.tableCode)
            self.readOEF_Data_table4()
        #elif self.tableName in ['OEFATO2','OEFCRM2','OEFPSD2','OEFRMS2','OEFNO2',]:
            #self.skipOES_Element() # skipping entire table
        else:
            self.NotImplementedOrSkip()
        ###

    def readOEF_Data_table4(self):  # Forces/Heat Flux
        if self.thermal in [0, 8]:
            self.readOEF_Forces()
        elif self.thermal == 1:
            self.readOEF_Thermal()
        else:
            self.NotImplementedOrSkip('thermal=%s' % (self.thermal))

        #if self.thermal==8:
        #    self.NotImplementedOrSkip('thermal=%s' %(self.thermal))

    def readOEF_Forces(self):
        #print self.codeInformation()
        try:
            (numWideReal, numWideImag) = self.OEF_ForceCode()
        except KeyError:
            self.NotImplementedOrSkip()

        if self.elementType in [1, 3, 10]:  # CROD,CTUBE,CONROD
            resultName = 'rodForces'
            if self.numWide == numWideReal:
                self.createTransientObject(self.rodForces, RealRodForce)
                self.handleResultsBuffer3(self.OEF_Rod, resultName)
            elif self.numWide == numWideImag:
                self.createTransientObject(self.rodForces, ComplexRodForce)
                self.handleResultsBuffer3(self.OEF_Rod_alt, resultName)
            else:
                self.NotImplementedOrSkip()

        elif self.elementType in [2]:  # CBEAM
            resultName = 'beamForces'
            if self.numWide == numWideReal:
                self.createTransientObject(self.beamForces, RealCBeamForce)
                self.handleResultsBuffer3(self.OEF_Beam, resultName)
            elif self.numWide == numWideImag:
                self.createTransientObject(self.beamForces, ComplexCBeamForce)
                self.handleResultsBuffer3(self.OEF_Beam_alt, resultName)
            else:
                self.NotImplementedOrSkip()

        elif self.elementType in [4]:  # CSHEAR
            resultName = 'shearForces'
            if self.numWide == numWideReal:
                self.createTransientObject(self.shearForces, RealCShearForce)
                self.handleResultsBuffer3(self.OEF_Shear, resultName)
            elif self.numWide == numWideImag:
                self.createTransientObject(
                    self.shearForces, ComplexCShearForce)
                self.handleResultsBuffer3(self.OEF_Shear_alt, resultName)
            else:
                self.NotImplementedOrSkip()

        elif self.elementType in [11, 12, 13, 14]:  # CELAS1,CELAS2,CELAS3,CELAS4
            resultName = 'springForces'
            if self.numWide == numWideReal:  # @todo is this correct or is DMAP wrong (CELAS1)???
                self.createTransientObject(self.springForces, RealSpringForce)
                self.handleResultsBuffer3(self.OEF_Spring, resultName)
            elif self.numWide == numWideImag:
                self.createTransientObject(
                    self.springForces, ComplexSpringForce)
                self.handleResultsBuffer3(self.OEF_Spring_alt, resultName)
            else:
                self.NotImplementedOrSkip()

        elif self.elementType in [20, 21, 22, 23]:  # CDAMP1,CDAMP2,CDAMP3,CDAMP4
            resultName = 'damperForces'
            if self.numWide == numWideReal:
                self.createTransientObject(self.damperForces, RealDamperForce)
                self.handleResultsBuffer3(self.OEF_Spring,
                                          resultName)  # same reader as for springs
            elif self.numWide == numWideImag:
                self.createTransientObject(
                    self.damperForces, ComplexDamperForce)
                self.handleResultsBuffer3(self.OEF_Spring_alt, resultName)
            else:
                self.NotImplementedOrSkip()

        elif self.elementType in [24]:  # CVISC
            resultName = 'viscForces'
            if self.numWide == numWideReal:
                self.createTransientObject(self.viscForces, RealViscForce)
                self.handleResultsBuffer3(self.OEF_CVisc, resultName)
            elif self.numWide == numWideImag:
                self.createTransientObject(self.viscForces, ComplexViscForce)
                self.handleResultsBuffer3(self.OEF_CVisc_alt, resultName)
            else:
                self.NotImplementedOrSkip()

        elif self.elementType in [33, 74, 235]:  # CQUAD4,CTRIA3,CQUADR
            resultName = 'plateForces'
            if self.numWide == numWideReal:
                #print self.codeInformation()
                self.createTransientObject(self.plateForces, RealPlateForce)
                self.handleResultsBuffer3(self.OEF_Plate, resultName)
            elif self.numWide == numWideImag:
                self.createTransientObject(self.plateForces, ComplexPlateForce)
                self.handleResultsBuffer3(self.OEF_Plate_alt, resultName)
            else:
                self.NotImplementedOrSkip()

        #elif self.elementType in [235]: # CQUADR
            #if self.numWide == numWideReal:
                #self.OEF_Plate()
            #elif self.numWide == numWideImag:
                #self.OEF_Plate_alt()
            #else:
                #raise NotImplementedError(self.codeInformation())
        elif self.elementType in [35]:  # CCONEAX
            resultName = 'coneAxForces'
            if self.numWide == numWideReal:
                self.createTransientObject(self.coneAxForces, RealConeAxForce)
                self.handleResultsBuffer3(self.OEF_ConeAx, resultName)
            else:
                self.NotImplementedOrSkip()
        elif self.elementType in [64, 70, 75, 82, 144]:  # CQUAD8,CTRIAR,CTRIA6,CQUADR,CQUAD4-bilinear
            resultName = 'plateForces2'
            if self.numWide == numWideReal:
                self.createTransientObject(self.plateForces2, RealPLATE2Force)
                self.handleResultsBuffer3(self.OEF_Plate2, resultName)
            elif self.numWide == numWideImag:
                self.createTransientObject(
                    self.plateForces2, ComplexPLATE2Force)
                self.handleResultsBuffer3(self.OEF_Plate2_alt, resultName)
            else:
                self.NotImplementedOrSkip()
        elif self.elementType in [34]:  # CBAR
            resultName = 'barForces'
            if self.numWide == numWideReal:
                self.createTransientObject(self.barForces, RealCBARForce)
                self.handleResultsBuffer3(self.OEF_CBar, resultName)
            elif self.numWide == numWideImag:
                self.createTransientObject(self.barForces, ComplexCBARForce)
                self.handleResultsBuffer3(self.OEF_CBar_alt, resultName)
            else:
                self.NotImplementedOrSkip()
        elif self.elementType in [100]:  # CBAR
            resultName = 'barForces'
            if self.numWide == numWideReal:
                self.createTransientObject(self.bar100Forces, RealCBAR100Force)
                self.handleResultsBuffer3(self.OEF_CBar100, resultName)
            elif self.numWide == numWideImag:
                self.handleResultsBuffer3(self.OEF_CBar100_alt, resultName)
            else:
                self.NotImplementedOrSkip()
        elif self.elementType in [38]:  # CGAP
            resultName = 'gapForces'
            if self.numWide == numWideReal:
                self.createTransientObject(self.gapForces, RealCGAPForce)
                self.handleResultsBuffer3(self.OEF_CGap, resultName)
            elif self.numWide == numWideImag:
                self.handleResultsBuffer3(self.OEF_CGap_alt, resultName)
            else:
                self.NotImplementedOrSkip()
        elif self.elementType in [69]:  # CBEND
            resultName = 'bendForces'
            if self.numWide == numWideReal:
                self.createTransientObject(self.bendForces, RealBendForce)
                self.handleResultsBuffer3(self.OEF_Bend, resultName)
            elif self.numWide == numWideImag:
                self.createTransientObject(self.bendForces, ComplexBendForce)
                self.handleResultsBuffer3(self.OEF_Bend_alt, resultName)
            else:
                self.NotImplementedOrSkip()
        elif self.elementType in [76, 77, 78]:  # CHEXA_PR,PENTA_PR,CTETRA_PR
            resultName = 'solidPressureForces'
            if self.numWide == numWideReal:
                self.createTransientObject(self.solidPressureForces,
                                           RealPentaPressureForce)
                self.handleResultsBuffer3(self.OEF_PentaPressure, resultName)
            elif self.numWide == numWideImag:
                self.createTransientObject(self.solidPressureForces,
                                           ComplexPentaPressureForce)
                self.handleResultsBuffer3(
                    self.OEF_PentaPressure_alt, resultName)
            else:
                self.NotImplementedOrSkip()
        #elif self.elementType in [95,96,97,98]: # composite CQUAD4,CQUAD8,CTRIA3,CTRIA6
        #    resultName = '???' ## @todo why is there no object...
        #    if self.numWide == numWideReal:
        #        self.handleResultsBuffer3(self.OEF_CompositePlate,resultName)
        #    elif self.numWide == numWideImag:
        #        self.handleResultsBuffer3(self.OEF_CompositePlate_alt)
        #    else:
        #        self.NotImplementedOrSkip()
        elif self.elementType in [102]:  # CBUSH
            resultName = 'bushForces'
            if self.numWide == numWideReal:
                self.createTransientObject(self.bushForces, RealCBUSHForce)
                self.handleResultsBuffer3(self.OEF_CBush, resultName)
            elif self.numWide == numWideImag:
                self.createTransientObject(self.bushForces, ComplexCBUSHForce)
                self.handleResultsBuffer3(self.OEF_CBush_alt, resultName)
            else:
                self.NotImplementedOrSkip()
        elif self.elementType in [189, 190]:  # VUQUAD,VUTRIA
            resultName = 'force_VU_2D'
            if self.numWide == numWideReal:
                self.createTransientObject(self.force_VU_2D, RealForce_VU_2D)
                self.handleResultsBuffer3(self.OEF_Force_VUTRIA, resultName)
            elif self.numWide == numWideImag:
                self.createTransientObject(
                    self.force_VU_2D, ComplexForce_VU_2D)
                self.handleResultsBuffer3(
                    self.OEF_Force_VUTRIA_alt, resultName)
            else:
                self.NotImplementedOrSkip()
        elif self.elementType in [191]:  # VUBEAM
            resultName = 'force_VU'
            if self.numWide == numWideReal:
                self.createTransientObject(self.force_VU, RealForce_VU)
                self.handleResultsBuffer3(self.OEF_Force_VU, resultName)
            elif self.numWide == numWideImag:
                self.createTransientObject(self.force_VU, ComplexForce_VU)
                self.handleResultsBuffer3(self.OEF_Force_VU_alt, resultName)
            else:
                self.NotImplementedOrSkip()
        else:
            self.NotImplementedOrSkip()
        ###
