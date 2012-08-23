import sys
from struct import unpack

from pyNastran.op2.tables.oes_stressStrain.real.elementsStressStrain import RealElementsStressStrain

class OESNLXR(RealElementsStressStrain):  ##todo real or complex?? see r767
    """Table of stresses/strains"""

    def readTable_OESNLXR(self):
        table3 = self.readTable_OESNLXR_3
        table4Data = self.readTable_OES_4_Data
        self.readResultsTable(
            table3, table4Data, flag=1)  # flag=1 defines old style
        self.deleteAttributes_OES()

    def deleteAttributes_OESNLXR(self):
        params = ['sCode', 'elementType', 'obj', 'markerStart', 'loadSet', 'formatCode', 'sCode', 'thermal',
                  'lsdvmn', 'mode', 'eign', 'modeCycle', 'freq', 'mode', 'eigr', 'eigi', 'dt']
        self.deleteAttributes(params)

    def readTable_OESNLXR_3(self, iTable):
        #print "*iTable3 = ",iTable
        #if 0:
            #markers = self.readMarkers([0,2])
            #print "markers=%s" %(markers)
            #block = self.readBlock()
            #print "block = ",block
            #markers = self.readMarkers([-1,7])
            #print "markers=%s" %(markers)
            #print self.printSection(200)

        bufferWords = self.getBufferWords()

        data = self.getData(4)
        bufferSize, = unpack('i', data)
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

        #print "loadset=%s formatCode=%s numWordsEntry=%s sCode=%s" %(self.loadSet,self.formatCode,self.numWide,self.sCode)
        #print "thermal(23)=%s elementType(3)=%s" %(self.thermal,self.elementType)

        ## assuming tCode=1
        if self.analysisCode == 1:   # statics / displacement / heat flux
            ## load set number
            self.addDataParameter(data, 'lsdvmn', 'i', 5, False)
        elif self.analysisCode == 2:  # real eigenvalues
            ## mode number
            self.addDataParameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            self.addDataParameter(data, 'eign', 'f', 6, False)
            ## mode or cycle
            ## @todo confused on the type - F1???
            self.addDataParameter(data, 'modeCycle', 'f', 7, False)
        #elif self.analysisCode==3: # differential stiffness
        #    ## load set number
        #    self.lsdvmn = self.getValues(data,'i',5)
        #elif self.analysisCode==4: # differential stiffness
        #    self.lsdvmn = self.getValues(data,'i',5) ## load set number

        elif self.analysisCode == 5:   # frequency
            ## frequency
            self.addDataParameter(data, 'freq', 'f', 5)
        elif self.analysisCode == 6:  # transient
            ## time step
            self.addDataParameter(data, 'dt', 'f', 5)
            #print "DT(5)=%s" %(self.dt)
        elif self.analysisCode == 7:  # pre-buckling
            ## load set
            self.addDataParameter(data, 'lsdvmn', 'i', 5)
            #print "LSDVMN(5)=%s" %(self.lsdvmn)
        elif self.analysisCode == 8:  # post-buckling
            ## mode number
            self.addDataParameter(data, 'lsdvmn', 'i', 5)
            ## real eigenvalue
            self.addDataParameter(data, 'eigr', 'f', 6, False)
            #print "LSDVMN(5)=%s  EIGR(6)=%s" %(self.lsdvmn,self.eigr)
        elif self.analysisCode == 9:  # complex eigenvalues
            ## mode number
            self.addDataParameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            self.addDataParameter(data, 'eigr', 'f', 6, False)
            ## imaginary eigenvalue
            self.addDataParameter(data, 'eigi', 'f', 7, False)
            #print "mode(5)=%s  eigr(6)=%s  eigi(7)=%s" %(self.mode,self.eigr,self.eigi)
        elif self.analysisCode == 10:  # nonlinear statics
            ## load step
            self.addDataParameter(data, 'lftsfq', 'f', 5)
            #print "LFTSFQ(5) = %s" %(self.lftsfq)
        elif self.analysisCode == 11:  # old geometric nonlinear statics
            ## load set number
            self.addDataParameter(data, 'lsdvmn', 'i', 5)
            #print "LSDVMN(5)=%s" %(self.lsdvmn)
        elif self.analysisCode == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## Time step ??? --> straight from DMAP
            self.addDataParameter(data, 'dt', 'f', 5)
        else:
            raise RuntimeError('Invalid Analysis Code: analysisCode=%s' % (self.analysisCode))
        ###
        # tCode=2
        #if self.analysisCode==2: # sort2
        #    self.lsdvmn = self.getValues(data,'i',5)

        self.readTitle()
        #print "n4 = ",self.n

    def readTable_OESNLXR_4_Data(self, iTable):
        isTable4Done = False
        isBlockDone = False
        #print self.printSection(100)

        data = self.getData(16)
        #print self.printBlock(data) # on
        #print "16 block..."
        #self.printBlock(data)
        #self.printBlock(data)
        bufferWords, = unpack('i', data[4:8])
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
            return isTable4Done, isBlockDone
        if bufferWords == 0:
            #print "bufferWords 0 - done with Table4"
            isTable4Done = True
            #isBlockDone  = True
            self.printSection(40)
            #print "exitB"
            return isTable4Done, isBlockDone

        self.readOESNLXR_ElementTable()
        return isTable4Done, isBlockDone

    def readOESNLXR_ElementTable(self):
        #print "**self.readElementTable"
        #print "*elementType = ",self.elementType
        #print "op2.tell=%s n=%s" %(self.op2.tell(),self.n)

        self.rewind(4)
        self.data = self.readBlock()  # 348
        #print "len(self.data) = ",len(self.data)

        if self.makeOp2Debug:
            self.op2Debug.write('reading big data block\n')
        #print self.printBlock(self.data)

        #msg = 'elementType=%s -> %s' %(self.elementType,self.ElementType(self.elementType))
        tfsCode = [self.tableCode, self.formatCode, self.sortCode]

        self.parseStressCode()

        if not self.isValidSubcase():  # lets the user skip a certain subcase
            self.log.debug("***skipping table=%s iSubcase=%s" %
                           (self.tableName, self.iSubcase))
            self.skipOES_Element()
        elif self.thermal == 0:
            # Stress / Strain
            self.dataCode['elementName'] = self.ElementType(self.elementType)
            if   tfsCode == [5, 1, 0]:
                self.readOESNLXR_Data_format1_sort0()
            else:
                raise RuntimeError('invalid atfsCode=%s' % (self.atfsCode))
                self.skipOES_Element()
                pass
        else:
            raise RuntimeError('invalid thermal option...')
        ###

        #print self.obj

    def readOESNLXR_Data_format1_sort0(self):
        raise NotImplementedError('code this...')
        #msg = 'OES elementType=%-3s -> %-6s\n' %(self.elementType,self.ElementType(self.elementType))
        msg = ''
        
        if self.elementType in [1, 3, 10]:  # crod/ctube/conrod
            self.makeOES_Object(self.rodStress, rodStressObject,
                                self.rodStrain, rodStrainObject)
            self.basicElement()
        elif self.elementType == 2:   # cbeam
            #print "    found cbeam_2"
            self.dataCode['elementName'] = 'CBEAM'
            self.makeOES_Object(self.beamStress, beamStressObject,
                                self.beamStrain, beamStrainObject)
            self.CBEAM_2()

        elif self.elementType in [4]:  # cshear
            #print "    found crod_1"
            self.dataCode['elementName'] = 'CSHEAR'
            self.makeOES_Object(self.shearStress, shearStressObject,
                                self.shearStrain, shearStrainObject)
            self.basicElement()
        elif self.elementType in [11, 12, 13]:   # celas1/celas2/celas3
            self.makeOES_Object(self.celasStress, celasStressObject,
                                self.celasStrain, celasStrainObject)
            self.basicElement()
        elif self.elementType == 34:   # cbar
            #print "    found cbar_34"
            self.dataCode['elementName'] = 'CBAR'
            self.makeOES_Object(self.barStress, barStressObject,
                                self.barStrain, barStrainObject)
            self.CBAR_34()

        elif self.elementType == 33:  # cquad4_33
            #print "    found cquad_33"
            self.dataCode['elementName'] = 'CQUAD4'
            self.makeOES_Object(self.plateStress, plateStressObject,
                                self.plateStrain, plateStrainObject)
            self.CQUAD4_33()
        elif self.elementType == 74:  # ctria
            #print "    found ctria_74"
            self.dataCode['elementName'] = 'CTRIA3'
            self.makeOES_Object(self.plateStress, plateStressObject,
                                self.plateStrain, plateStrainObject)
            self.CTRIA3_74()  # ctria3
        elif self.elementType in [64, 144, 70, 75]:  # cquad8/cquad4/ctriar/ctria6
            #print "    found cquad_144"
            if     self.elementType == 64:
                self.dataCode['elementName'] = 'CQUAD8'
            elif   self.elementType == 144:
                self.dataCode['elementName'] = 'CQUAD4'
            elif   self.elementType == 70:
                self.dataCode['elementName'] = 'CTRIAR'
            elif   self.elementType == 75:
                self.dataCode['elementName'] = 'CTRIA6'
            else:
                raise NotImplementedError(self.elementType)
            self.makeOES_Object(self.plateStress, plateStressObject,
                                self.plateStrain, plateStrainObject)
            self.CQUAD4_144()

        elif self.elementType in [39, 67, 68]:   # ctetra/chexa/cpenta
            #print "    found ctetra_39 / hexa_67 / cpenta_68"
            self.makeOES_Object(self.solidStress, solidStressObject,
                                self.solidStrain, solidStrainObject)
            self.CSOLID_67()

        elif self.elementType in [85]:   # ctetra/chexa/cpenta (91,93)
            #print "    found ctetra_85 / hexa_93 / cpenta_91"
            self.makeOES_Object(self.solidStress, solidStressObject,
                                self.solidStrain, solidStrainObject)
            self.CSOLID_85()

        elif self.elementType in [89, 92]:   # RODNL
            #print "    found RODNL_89"
            self.makeOES_Object(self.rodStress, nonlinearRodObject,
                                self.rodStrain, nonlinearRodObject)
            self.RODNL_89_92()

        #elif self.elementType in [91]: # CPENTANL
        #    #print "hexa_93"
        #    self.CPENTANL_91()
        #elif self.elementType in [93]: # CHEXANL
        #    #print "hexa_93"
        #    self.CHEXANL_93()

        elif self.elementType in [95, 96, 97, 98]:  # CQUAD4, CQUAD8, CTRIA3, CTRIA6 (composite)
            #print "    found a 95/96/97 or 98!"
            self.eid2 = None  # stores the previous elementID
            self.makeOES_Object(
                self.compositePlateStress, compositePlateStressObject,
                self.compositePlateStrain, compositePlateStrainObject)
            self.CQUAD4_95()
            del self.eid2

        elif self.elementType in [94]:  # CBEAM (nonlinear)
            print "    found a 94!"
            #self.eid2 = None # stores the previous elementID
            self.makeOES_Object(self.beamStress, beamStressObject,
                                self.beamStrain, beamStrainObject)
            self.CBEAM_94()
            #sys.exit('stoping at end of CBEAM_94')
            #del self.eid2

        elif self.elementType in [139]:   # QUAD4FD
            #print "    found QUAD4FD_139"
            self.makeOES_Object(self.plateStress, nonlinearQuadObject,
                                self.plateStrain, nonlinearQuadObject)
            self.QUAD4FD_139()

        else:
            #self.printBlock(self.data[0:100])
            self.skipOES_Element()
            msg = 'OES format1_sort0 elementType=%-3s -> %s is not supported - fname=%s\n' % (self.elementType, self.ElementType(self.elementType), self.op2FileName)
            self.log.debug(msg)
            #msg = 'OES format1_sort0 elementType=%-3s -> %s is not supported' %(self.elementType,self.ElementType(self.elementType))
            #raise RuntimeError(msg)
            self.skippedCardsFile.write(msg)
        ###
