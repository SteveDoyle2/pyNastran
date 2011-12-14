import sys
import copy
from struct import unpack

from pyNastran.op2.op2Errors import *
from elementsStressStrain import ElementsStressStrain
from oes_rods   import rodStressObject,rodStrainObject
from oes_bars   import barStressObject,barStrainObject
from oes_beams  import beamStressObject,beamStrainObject

from oes_solids import solidStressObject,solidStrainObject
from oes_plates import plateStressObject,plateStrainObject
from oes_compositePlates import compositePlateStressObject,compositePlateStrainObject
from oes_springs import celasStressObject,celasStrainObject

class OES(ElementsStressStrain):
    """Table of stresses/strains"""

    def readTable_OES1(self):
        self.tableName = 'OES'
        table3 = self.readTable_OES_3
        table4Data = self.readTable_OES_4_Data
        self.readResultsTable(table3,table4Data,flag=1) # flag=1 defines old style
        self.deleteAttributes_OES()

    def deleteAttributes_OES(self):
        params = ['sCode','elementType','obj','markerStart','loadSet','formatCode','sCode','thermal',
                  'lsdvmn','mode','eign','modeCycle','freq','mode','eigr','eigi','dt','nonlinearFactor']
        self.deleteAttributes(params)

    def readTable_OES_3(self,iTable):
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
        #print "data = ".data
        bufferSize, = unpack('i',data)
        if self.makeOp2Debug:
            self.op2Debug.write('bufferSize=|%s|\n' %(str(bufferSize)))

        data = self.getData(4*50)
        self.printBlock(data) # on
        if self.makeOp2Debug:
            self.op2Debug.write('block3header\n')

        self.elementType = self.parseApproachCode(data) # 3
        self.loadSet    = self.getValues(data,'i',8)
        self.formatCode = self.getValues(data,'i',9)
        self.numWide    = self.getValues(data,'i',10)
        self.sCode      = self.getValues(data,'i',11)
        self.thermal    = self.getValues(data,'i',23) # 1 is heat transfer, 0 otherwise

        self.dataCode = {'analysisCode': self.analysisCode,'deviceCode':self.deviceCode,
                         'loadSet':self.loadSet,'formatCode':self.formatCode,'sortCode':self.sortCode,
                         'numWide': self.numWide,'sCode':self.sCode,
                         'thermal': self.thermal,'elementType':self.elementType}
        print "self.dataCode = ",self.dataCode
        #print "loadset=%s formatCode=%s numWordsEntry=%s sCode=%s" %(self.loadSet,self.formatCode,self.numWide,self.sCode)
        #print "thermal(23)=%s elementType(3)=%s" %(self.thermal,self.elementType)

        ## assuming tCode=1
        self.nonlinearFactor = None
        if self.analysisCode==1:   # statics / displacement / heat flux
            self.lsdvmn = self.getValues(data,'i',5) ## load set number
            self.dataCode['lsdvmn'] = self.lsdvmn
            self.dataCode['name'] = 'lsdvmn'
        elif self.analysisCode==2: # real eigenvalues
            self.mode      = self.getValues(data,'i',5) ## mode number
            self.eign      = self.getValues(data,'f',6) ## real eigenvalue
            self.modeCycle = self.getValues(data,'f',7) ## mode or cycle @todo confused on the type ???
            self.dataCode['mode'] = self.mode
            self.dataCode['eign'] = self.eign
            self.dataCode['modeCycle'] = self.modeCycle
            self.dataCode['name'] = 'mode'
            self.nonlinearFactor = self.mode
        #elif self.analysisCode==3: # differential stiffness
        #    self.lsdvmn = self.getValues(data,'i',5) ## load set number
        #    self.nonlinearFactor = self.lsdvmn
        #elif self.analysisCode==4: # differential stiffness
        #    self.lsdvmn = self.getValues(data,'i',5) ## load set number
        #    self.nonlinearFactor = self.lsdvmn

        elif self.analysisCode==5:   # frequency
            self.freq = self.getValues(data,'f',5) ## frequency
            self.dataCode['freq'] = self.freq
            self.dataCode['name'] = 'freq'
            self.nonlinearFactor = self.freq
        elif self.analysisCode==6: # transient
            self.dt = self.getValues(data,'f',5) ## time step
            self.dataCode['dt'] = self.dt
            self.dataCode['name'] = 'dt'
            
            print "DT(5)=%s" %(self.dt)
            self.nonlinearFactor = self.dt
        elif self.analysisCode==7: # pre-buckling
            self.lsdvmn = self.getValues(data,'i',5) ## load set
            self.dataCode['lsdvmn'] = self.lsdvmn
            self.nonlinearFactor = self.lsdvmn
            print "LSDVMN(5)=%s" %(self.lsdvmn)
        elif self.analysisCode==8: # post-buckling
            self.lsdvmn = self.getValues(data,'i',5) ## mode number
            self.eigr   = self.getValues(data,'f',6) ## real eigenvalue
            self.dataCode['lsdvmn'] = self.lsdvmn
            self.dataCode['eigr']   = self.eigr
            self.nonlinearFactor = self.lsdvmn
            print "LSDVMN(5)=%s  EIGR(6)=%s" %(self.lsdvmn,self.eigr)
        elif self.analysisCode==9: # complex eigenvalues
            self.mode   = self.getValues(data,'i',5) ## mode
            self.eigr   = self.getValues(data,'f',6) ## real eigenvalue
            self.eigi   = self.getValues(data,'f',7) ## imaginary eigenvalue
            self.nonlinearFactor = self.mode
            self.dataCode['mode'] = self.mode
            self.dataCode['eigr'] = self.eigr
            self.dataCode['eigi'] = self.eigi
            self.dataCode['name'] = 'mode'
            print "mode(5)=%s  eigr(6)=%s  eigi(7)=%s" %(self.mode,self.eigr,self.eigi)
        elif self.analysisCode==10: # nonlinear statics
            self.lftsfq = self.getValues(data,'f',5) ## load step
            self.dataCode['lftsfq'] = self.lftsfq
            self.dataCode['name'] = 'lftsfq'
            self.nonlinearFactor = self.lftsfq
            print "LFTSFQ(5) = %s" %(self.lftsfq)
        elif self.analysisCode==11: # old geometric nonlinear statics
            self.lsdvmn = self.getValues(data,'i',5)
            self.dataCode['lsdvmn'] = self.lsdvmn
            self.dataCode['name'] = 'lsdvmn'
            self.nonlinearFactor = self.lsdvmn
            print "LSDVMN(5)=%s" %(self.lsdvmn)
        elif self.analysisCode==12: # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            #self.lsdvmn = self.getValues(data,'i',5)
            self.dt = self.getValues(data,'f',5)  ## Time step ??? --> straight from DMAP
            self.dataCode['dt'] = self.dt
            self.dataCode['name'] = 'dt'
            self.nonlinearFactor = self.dt
            print "LSDVMN(5)=%s" %(self.lsdvmn)
        else:
            raise InvalidAnalysisCodeError('invalid analysisCode...analysisCode=%s' %(self.analysisCode))
        self.dataCode['nonlinearFactor'] = self.nonlinearFactor
        # tCode=2
        #if self.analysisCode==2: # sort2
        #    self.lsdvmn = self.getValues(data,'i',5)

        print self.codeInformation()
        self.readTitle()
        #print "n4 = ",self.n

    def parseStressCode(self):
        bits = [0,0,0,0,0]
        
        sCode = self.sCode
        i=4
        while sCode>0:
            value = sCode%2
            sCode = (sCode - value)/2
            bits[i] = value
            #print "    *bit = ",value
            #print "    sCode = ",sCode
            i-=1
        #bits.reverse()
        print "bits = ",bits
        self.stressBits = bits
        self.dataCode['stressBits'] = self.stressBits

    def readTable_OES_4_Data(self,iTable):
        print "-------------"
        print "**self.readTable_OES_4_Data"
        isTable4Done = False
        isBlockDone  = False
        #print self.printSection(100)

        data = self.getData(16)
        #print self.printBlock(data) # on
        #print "16 block..."
        #self.printBlock(data)
        #self.printBlock(data)
        bufferWords, = unpack('i',data[4:8])
        print "bufferWords = ",bufferWords
        if self.makeOp2Debug:
            self.op2Debug.write('bufferWords=|%s|\n' %(str(bufferWords)))

        #print "*********************"
        #bufferWords = self.getMarker() # 87 - buffer
        #print "OES4 bufferWords = ",bufferWords,bufferWords*4
        #self.verifyBufferSize(bufferWords)
        
        isBlockDone = not(bufferWords)
        #print "self.firstPass = ",self.firstPass
        
        ## table -4 is done, restarting table -3
        if self.isBufferDone: # table is done when the buffer is done
            isTable4Done = True
            #print "exitA"
            return isTable4Done,isBlockDone
        if bufferWords==0:
            #print "bufferWords 0 - done with Table4"
            isTable4Done = True
            #isBlockDone  = True
            self.printSection(40)
            #print "exitB"
            return isTable4Done,isBlockDone

        self.readElementTable()
        return isTable4Done,isBlockDone

    def readElementTable(self):
        #print "**self.readElementTable"
        print "*elementType = ",self.elementType
        #print "op2.tell=%s n=%s" %(self.op2.tell(),self.n)
        
        self.rewind(4)
        self.data = self.readBlock()  # 348
        #print "len(self.data) = ",len(self.data)

        if self.makeOp2Debug:
            self.op2Debug.write('reading big data block\n')
        #print "self.n = ",self.n
        #print 'reading big data block'
        #print self.printBlock(self.data)

        #msg = 'elementType=%s -> %s' %(self.elementType,self.ElementType(self.elementType))
        tfsCode = [self.tableCode,self.formatCode,self.sortCode]

        self.parseStressCode()
        if self.thermal==0:
            # Stress / Strain
            if tfsCode==[5,1,0]:
                self.readOES1_Data_format1_sort0()
            elif tfsCode==[5,1,1]:
                self.readOES1_Data_format1_sort1()
            elif tfsCode==[5,2,0]:
                self.readOES1_Data_format2_sort0()
            elif tfsCode==[5,2,1]:
                self.readOES1_Data_format2_sort1()
            elif tfsCode==[5,3,0]:
            
                self.readOES1_Data_format3_sort0()
            elif tfsCode==[5,3,1]:
                self.readOES1_Data_format3_sort1()
            elif tfsCode==[5,3,2]:
                self.readOES1_Data_format3_sort2()
            else:
                raise Exception('invalid atfsCode=%s' %(self.atfsCode))
            ###
        elif self.thermal==1:
            self.OES_Thermal(None)
            #raise Exception('thermal stress')
        else:
            raise Exception('invalid thermal option...')
        ###
        #print self.obj

    def readOES1_Data_format1_sort1(self):
        #msg = 'OES format3_sort2 elementType=%-3s -> %-6s is not supported - fname=%s\n' %(self.elementType,self.ElementType(self.elementType),self.op2FileName)
        msg = ''
        if 0:
            pass
        else:
            #self.printBlock(self.data[0:100])
            self.skipOES_Element()
            msg = 'OES format1_sort1 elementType=%-3s -> %-6s is not supported - fname=%s\n' %(self.elementType,self.ElementType(self.elementType),self.op2FileName)
            print msg
            #raise RuntimeError(msg)
        ###
        self.skippedCardsFile.write(msg)

    def readOES1_Data_format2_sort0(self):
        #msg = 'OES format2_sort0 elementType=%-3s -> %-6s is not supported - fname=%s\n' %(self.elementType,self.ElementType(self.elementType),self.op2FileName)
        msg = ''
        if 0:
            pass
        else:
            #self.printBlock(self.data[0:100])
            self.skipOES_Element()
            msg = 'OES format2_sort0 elementType=%-3s -> %-6s is not supported - fname=%s\n' %(self.elementType,self.ElementType(self.elementType),self.op2FileName)
            print msg
            #raise RuntimeError(msg)
        ###
        self.skippedCardsFile.write(msg)

    def readOES1_Data_format2_sort1(self):
        #msg = 'OES format2_sort1 elementType=%-3s -> %-6s is not supported - fname=%s\n' %(self.elementType,self.ElementType(self.elementType),self.op2FileName)
        msg = ''
        if self.elementType==1: # crod
            #print "    found crod_1"
            self.makeOES_Object(self.rodStress,rodStressObject,
                                               self.rodStrain,rodStrainObject)
            self.basicElement()
        #elif self.elementType in [10,11,12,33,74]:
        #    raise AddNewElementError('add element=%s' %(self.elementType))
        else:
            #self.printBlock(self.data[0:100])
            self.skipOES_Element()
            msg = 'OES format2_sort1 elementType=%-3s -> %-6s is not supported - fname=%s\n' %(self.elementType,self.ElementType(self.elementType),self.op2FileName)
            print msg
            #raise RuntimeError(msg)
        ###
        self.skippedCardsFile.write(msg)

    def readOES1_Data_format3_sort0(self):
        msg = 'OES format3_sort0 elementType=%-3s -> %-6s is not supported - fname=%s\n' %(self.elementType,self.ElementType(self.elementType),self.op2FileName)
        #msg = ''
        if 0:
            pass
        else:
            #self.printBlock(self.data[0:100])
            self.skipOES_Element()
            msg = 'OES format3_sort0 elementType=%-3s -> %-6s is not supported - fname=%s\n' %(self.elementType,self.ElementType(self.elementType),self.op2FileName)
            print msg
            #raise RuntimeError(msg)
        ###
        self.skippedCardsFile.write(msg)

    def readOES1_Data_format3_sort1(self):
        #msg = 'OES format3_sort1 elementType=%-3s -> %-6s is not supported - fname=%s\n' %(self.elementType,self.ElementType(self.elementType),self.op2FileName)
        msg = ''
        if 0:
            pass
        else:
            #self.printBlock(self.data[0:100])
            self.skipOES_Element()
            msg = 'OES format3_sort1 elementType=%-3s -> %-6s is not supported - fname=%s\n' %(self.elementType,self.ElementType(self.elementType),self.op2FileName)
            print msg
            #raise RuntimeError(msg)
        ###
        self.skippedCardsFile.write(msg)

    def readOES1_Data_format3_sort2(self):
        #msg = 'OES format3_sort2 elementType=%-3s -> %-6s is not supported - fname=%s\n' %(self.elementType,self.ElementType(self.elementType),self.op2FileName)
        msg = ''
        if 0:
            pass
        else:
            #self.printBlock(self.data[0:100])
            self.skipOES_Element()
            msg = 'OES format3_sort2 elementType=%-3s -> %-6s is not supported - fname=%s\n' %(self.elementType,self.ElementType(self.elementType),self.op2FileName)
            print msg
            #raise RuntimeError(msg)
        ###
        self.skippedCardsFile.write(msg)

    def readOES1_Data_format1_sort0(self):
        #msg = 'OES elementType=%-3s -> %-6s\n' %(self.elementType,self.ElementType(self.elementType))
        msg = ''
        assert self.analysisCode in [1,6,10],'self.atfsCode=%s' %(self.atfsCode)
        if self.elementType==1: # crod
            print "    found crod_1"
            self.dataCode['ElementName'] = 'CROD'
            self.makeOES_Object(self.rodStress,rodStressObject,
                                self.rodStrain,rodStrainObject)
            self.basicElement()
        elif self.elementType == 2:   # cbeam
            #print "    found cbeam_2"
            self.dataCode['ElementName'] = 'CBEAM'
            self.makeOES_Object(self.beamStress,beamStressObject,
                                self.beamStrain,beamStrainObject)
            self.CBEAM_2()
        elif self.elementType == 10:   # conrod
            #print "    found conrod_10"
            self.dataCode['ElementName'] = 'CONROD'
            self.makeOES_Object(self.conrodStress,conrodStressObject,
                                self.conrodStrain,conrodStrainObject)
            self.CONROD_10()

        elif self.elementType == 11:   # celas1
            #print "    found celas2_12"
            self.dataCode['ElementName'] = 'CELAS1'
            self.makeOES_Object(self.celasStress,celasStressObject,
                                self.celasStrain,celasStrainObject)
            self.basicElement()
        elif self.elementType == 12:   # celas2
            #print "    found celas2_12"
            self.dataCode['ElementName'] = 'CELAS2'
            self.makeOES_Object(self.celasStress,celasStressObject,
                                self.celasStrain,celasStrainObject)
            self.basicElement()
        elif self.elementType == 34:   # cbar
            #print "    found cbar_34"
            self.dataCode['ElementName'] = 'CBAR'
            self.makeOES_Object(self.barStress,barStressObject,
                                self.barStrain,barStrainObject)
            self.CBAR_34()

        elif self.elementType==33: # cquad4_33
            self.stopCode = True
            #print "    found cquad_33"
            self.dataCode['ElementName'] = 'CQUAD4'
            self.makeOES_Object(self.plateStress,plateStressObject,
                                self.plateStrain,plateStrainObject)
            self.CQUAD4_33()
        elif self.elementType==74:  # ctria
            #print "    found ctria_74"
            self.dataCode['ElementName'] = 'CTRIA3'
            self.makeOES_Object(self.plateStress,plateStressObject,
                                self.plateStrain,plateStrainObject)
            self.CTRIA3_74() # ctria3
        elif self.elementType==144: # cquad4
            self.stopCode = True
            #print "    found cquad_144"
            self.dataCode['ElementName'] = 'CQUAD4'
            self.makeOES_Object(self.plateStress,plateStressObject,
                                self.plateStrain,plateStrainObject)
            self.CQUAD4_144()

        elif self.elementType in [39,67,68]:   # ctetra/chexa/cpenta
            #print "    found ctetra_39 / hexa_67 / cpenta_68"
            self.makeOES_Object(self.solidStress,solidStressObject,
                                self.solidStrain,solidStrainObject)
            self.CSOLID_67()

        elif self.elementType in [85]:   # ctetra/chexa/cpenta (91,93)
            #print "    found ctetra_85 / hexa_93 / cpenta_91"
            self.makeOES_Object(self.solidStress,solidStressObject,
                                self.solidStrain,solidStrainObject)
            self.CSOLID_85()
        elif self.elementType in [91]: # CPENTANL
            #print "hexa_93"
            self.CPENTANL_91()
        elif self.elementType in [93]: # CHEXANL
            #print "hexa_93"
            self.CHEXANL_93()

        elif self.elementType in [95,96,97,98]: # CQUAD4, CQUAD8, CTRIA3, CTRIA6 (composite)
            #print "    found a 95/96/97 or 98!"
            self.eid2 = None # stores the previous elementID
            self.makeOES_Object(self.compositePlateStress,compositePlateStressObject,
                                self.compositePlateStrain,compositePlateStrainObject)
            self.CQUAD4_95()
            del self.eid2
        #elif self.elementType in [2,53,61,70,86,88,90,94,102,189,232,]:
        #    self.skipOES_Element()
            #elementType=2   -> BEAM    is not supported
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
        else:
            #self.printBlock(self.data[0:100])
            self.skipOES_Element()
            msg = 'OES format1_sort0 elementType=%-3s -> %s is not supported - fname=%s\n' %(self.elementType,self.ElementType(self.elementType),self.op2FileName)
            print msg
            #raise RuntimeError(msg)
        ###
        self.skippedCardsFile.write(msg)
        #elif self.elementType == 1:    # crod     (done)
        #elif self.elementType == 2:    # cbeam    (untested)
        #elif self.elementType == 3:    # ctube
        #elif self.elementType == 4:    # cshear
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
        #elif self.elementType == 64:   # cquad8
        #elif self.elementType == 67:   # chexa    (done)
        #elif self.elementType == 68:   # cpenta   (done)
        #elif self.elementType == 69:   # cbend
        #elif self.elementType == 70:   # ctriar
        #elif self.elementType == 74:   # ctria3  (done)
        #elif self.elementType == 75:   # ctria6
        #elif self.elementType == 82:   # cquadr
        #elif self.elementType == 95: # CQUAD4    (done)
        #elif self.elementType == 96: # CQUAD8    (done)
        #elif self.elementType == 97: # CTRIA3    (done)
        #elif self.elementType == 98: # CTRIA6    (done)
        #elif self.elementType == 100:  # cbar w/ cbarao or pload1

        #elif self.elementType == 102:  # cbush
        #elif self.elementType == 144:  # cquad_144 - corner stresses (done)


        # rods/bars/beams
        #elif self.elementType == 1:    # crod   (done)
        #elif self.elementType == 2:    # cbeam  (untested)
        #elif self.elementType == 3:    # ctube
        #elif self.elementType == 10:   # conrod (done)
        #elif self.elementType == 34:   # cbar   (done)

        #springs
        #elif self.elementType == 11:   # celas1
        #elif self.elementType == 12:   # celas2 (not integrated)
        #elif self.elementType == 13:   # celas3
        #elif self.elementType == 14:   # celas4
        
        #plate
        #elif self.elementType == 33:   # cquad_33 (done)

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
        #elif self.elementType == 85:   # tetra  (nonlinear,done)
        #elif self.elementType == 86:   # gap    (nonlinear)
        #elif self.elementType == 87:   # ctube  (nonlinear)
        #elif self.elementType == 88:   # tria3  (nonlinear) - same as quad4
        #elif self.elementType == 89:   # crod   (nonlinear)
        #elif self.elementType == 90:   # quad4  (nonlinear)
        #elif self.elementType == 91:   # cpenta (nonlinear)
        #elif self.elementType == 92:   # conrod (nonlinear)
        #elif self.elementType == 93:   # chexa  (nonlinear,done)
        #elif self.elementType == 94:   # cbeam  (nonlinear)
        
        # acoustic
        #elif self.elementType == 76:   # chexa  (acoustic)
        #elif self.elementType == 77:   # cpenta (acoustic)
        #elif self.elementType == 78:   # ctetra (acoustic)
        #elif self.elementType == 101:  # caabsf (acoustic)

    def isStress(self):
        if self.stressBits[1]==0:
            return True
        return False

    def makeOES_Object(self,stress,stressObject,strain,strainObject):
        """
        Creates a stress/strain object if necessary
        """
        if self.isStress():
            self.createTransientObject(stress,stressObject)
        else:
            self.createTransientObject(strain,strainObject)
        ###
        return self.obj
