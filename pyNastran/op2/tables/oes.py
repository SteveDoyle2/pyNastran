import sys
import copy
from struct import unpack
from pyNastran.op2.resultObjects.oes_Objects import (
    rodStressObject,rodStrainObject,
    barStressObject,barStrainObject,
    plateStressObject,plateStrainObject,
    solidStressObject,solidStrainObject,
    compositePlateStressObject,compositePlateStrainObject)

class OES(object):
    """Table of stresses/strains"""

    def readTable_OES1(self):
        table3 = self.readTable_OES_3
        table4Data = self.readTable_OES_4_Data
        self.readResultsTable(table3,table4Data,flag=1) # flag=1 defines old style
        self.deleteAttributes_OES()

    def deleteAttributes_OES(self):
        params = ['sCode','elementType','obj','markerStart','loadSet','formatCode','sCode','thermal',
                  'lsdvmn','mode','eign','modeCycle','freq','mode','eigr','eigi','dt','nonlinearFactor']
        self.deleteAttributes(params)

    def readTable_OES_3(self,iTable):
        #print "iTable3 = ",iTable
        bufferWords = self.getBufferWords()
        
        data = self.getData(4)
        bufferSize, = unpack('i',data)
        if self.makeOp2Debug:
            self.op2Debug.write('bufferSize=|%s|\n' %(str(bufferSize)))

        data = self.getData(4*50)
        #self.printBlock(data) # on
        if self.makeOp2Debug:
            self.op2Debug.write('block3header\n')

        self.elementType = self.parseApproachCode(data) # 3
        self.loadSet    = self.getValues(data,'i',8)
        self.formatCode = self.getValues(data,'i',9)
        self.numWide    = self.getValues(data,'i',10)
        self.sCode      = self.getValues(data,'i',11)
        self.thermal    = self.getValues(data,'i',23) # 1 is heat transfer, 0 otherwise

        self.dataCode = {'analysisCode': self.approachCode,'deviceCode':self.deviceCode,
                         'loadSet':self.loadSet,'formatCode':self.formatCode,
                         'numWide': self.numWide,'sCode':self.sCode,
                         'thermal': self.thermal}

        print "loadset=%s formatCode=%s numWordsEntry=%s sCode=%s" %(self.loadSet,self.formatCode,self.numWide,self.sCode)
        print "thermal(23)=%s elementType(3)=%s" %(self.thermal,self.elementType)

        ## assuming tCode=1
        self.nonlinearFactor = None
        if self.approachCode==1:   # statics / displacement / heat flux
            self.lsdvmn = self.getValues(data,'i',5) ## load set number
        elif self.approachCode==2: # real eigenvalues
            self.mode      = self.getValues(data,'i',5) ## mode number
            self.eign      = self.getValues(data,'f',6) ## real eigenvalue
            self.modeCycle = self.getValues(data,'f',7) ## mode or cycle @todo confused on the type ???
            self.nonlinearFactor = self.mode
        #elif self.approachCode==3: # differential stiffness
        #    self.lsdvmn = self.getValues(data,'i',5) ## load set number
        #    self.nonlinearFactor = self.lsdvmn
        #elif self.approachCode==4: # differential stiffness
        #    self.lsdvmn = self.getValues(data,'i',5) ## load set number
        #    self.nonlinearFactor = self.lsdvmn

        elif self.approachCode==5:   # frequency
            self.freq = self.getValues(data,'f',5) ## frequency
            self.nonlinearFactor = self.freq
        elif self.approachCode==6: # transient
            self.dt = self.getValues(data,'f',5) ## time step
            print "DT(5)=%s" %(self.dt)
            self.nonlinearFactor = self.dt
        elif self.approachCode==7: # pre-buckling
            self.lsdvmn = self.getValues(data,'i',5) ## load set
            self.nonlinearFactor = self.lsdvmn
            print "LSDVMN(5)=%s" %(self.lsdvmn)
        elif self.approachCode==8: # post-buckling
            self.lsdvmn = self.getValues(data,'i',5) ## mode number
            self.eigr   = self.getValues(data,'f',6) ## real eigenvalue
            self.nonlinearFactor = self.lsdvmn
            print "LSDVMN(5)=%s  EIGR(6)=%s" %(self.lsdvmn,self.eigr)
        elif self.approachCode==9: # complex eigenvalues
            self.mode   = self.getValues(data,'i',5) ## mode
            self.eigr   = self.getValues(data,'f',6) ## real eigenvalue
            self.eigi   = self.getValues(data,'f',7) ## imaginary eigenvalue
            print "mode(5)=%s  eigr(6)=%s  eigi(7)=%s" %(self.mode,self.eigr,self.eigi)
            self.nonlinearFactor = self.mode
        elif self.approachCode==10: # nonlinear statics
            self.lftsfq = self.getValues(data,'f',5) ## load step
            self.nonlinearFactor = self.lftsfq
            print "LFTSFQ(5) = %s" %(self.lftsfq)
        elif self.approachCode==11: # old geometric nonlinear statics
            self.lsdvmn = self.getValues(data,'i',5)
            self.nonlinearFactor = self.lsdvmn
            print "LSDVMN(5)=%s" %(self.lsdvmn)
        elif self.approachCode==12: # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            #self.lsdvmn = self.getValues(data,'i',5)
            self.dt = self.getValues(data,'f',5)  ## Time step ??? --> straight from DMAP
            self.nonlinearFactor = self.dt
            print "LSDVMN(5)=%s" %(self.lsdvmn)
        else:
            raise RuntimeError('invalid approach code...approachCode=%s' %(self.approachCode))
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
            print "    *bit = ",value
            #print "    sCode = ",sCode
            i-=1
        #bits.reverse()
        print "bits = ",bits
        return bits

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

        tfsCode = [self.tableCode,self.formatCode,self.sortCode]
        if self.thermal==0:
            if tfsCode==[5,1,0]:  # Stress
                self.readOES1_Data_format1_sort0()
            else:
                raise Exception('invalid atfsCode=%s' %(self.atfsCode))
            ###
        elif self.thermal==1:
            self.OES_Thermal(None)
            #raise Exception('thermal stress')
        else:
            raise Exception('invalid thermal option...')
        ###
        #print stressStrainObj


    def readOES1_Data_format1_sort0(self):
        if self.elementType==1: # crod
            #print "    found crod_1"
            stressStrainObj = self.instantiateRodObject()
            self.CROD_1(stressStrainObj)
        #elif self.elementType == 2:   # cbeam
        #    print "    found cbeam_2"
        #    #stressStrainObj = self.instantiateBeamObject()
        #    stressStrainObj = None
        #    self.CBEAM_2(stressStrainObj)
        #elif self.elementType == 10:   # conrod
            print "    found conrod_10"
            #stressStrainObj = self.instantiateConrodObject()
            stressStrainObj = None
            self.CONROD_10(stressStrainObj)

        elif self.elementType == 12:   # celas2
            #print "    found celas2_12"
            #stressStrainObj = self.instantiateCelasObject()
            stressStrainObj = None
            self.CELAS2_12(stressStrainObj)
        elif self.elementType == 34:   # cbar
            #print "    found cbar_34"
            stressStrainObj = self.instantiateBarObject()
            self.CBAR_34(stressStrainObj)

        elif self.elementType==33: # cquad4_33
            self.stopCode = True
            #print "    found cquad_33"
            stressStrainObj = self.instantiatePlateObject()
            self.CQUAD4_33(stressStrainObj)
        elif self.elementType==74:  # ctria
            #print "    found ctria_74"
            stressStrainObj = self.instantiatePlateObject()
            self.CTRIA3_74(stressStrainObj) # ctria3
        elif self.elementType==144: # cquad4
            self.stopCode = True
            #print "    found cquad_144"
            stressStrainObj = self.instantiatePlateObject()
            self.CQUAD4_144(stressStrainObj)

        elif self.elementType in [39,67,68]:   # ctetra/chexa/cpenta
            #print "    found ctetra_39 / hexa_67 / cpenta_68"
            stressStrainObj = self.instantiateSolidObject()
            self.CSOLID_67(stressStrainObj)

        elif self.elementType in [85]:   # ctetra/chexa/cpenta (91,93)
            #print "    found ctetra_85 / hexa_93 / cpenta_91"
            stressStrainObj = self.instantiateSolidObject()
            self.CSOLID_85(stressStrainObj)
        elif self.elementType in [93]: # CHEXANL
            stressStrainObj = None
            #print "hexa_93"
            self.CHEXANL_93(stressStrainObj)

        elif self.elementType in [95,96,97,98]: # CQUAD4, CQUAD8, CTRIA3, CTRIA6 (composite)
            #print "    found a 95/96/97 or 98!"
            self.eid2 = None # stores the previous elementID
            stressStrainObj = self.instantiateCompositePlateObject()
            self.CQUAD4_95(stressStrainObj)
            del self.eid2
        else:
            self.printBlock(self.data[0:100])
            self.skipOES_Element(None)
            msg = 'elementType=%s -> %s is not supported' %(self.elementType,self.ElementType(self.elementType))
            self.developerFile.write(msg+'\n')
            print msg
            #raise RuntimeError(msg)
        ###
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


    def instantiateRodObject(self): # 1 (CROD)
        """
        Creates a stress/strain object if necessary
        """
        #print "starting a ROD stress/strain object"
        #print "    iSubcase = ",self.iSubcase
        #print "    sCode    = ",self.sCode
        #bits = self.parseStressCode()
        #if (self.iSubcase not in self.rodStress) and (self.iSubcase not in self.rodStrain):
        #    print "making new subcase..."

        if self.sCode in [0,1]:
            self.createTransientObject(self.rodStress,rodStressObject)
        elif self.sCode in [10,11]:
            self.createTransientObject(self.rodStrain,rodStrainObject)
        else:
            raise Exception('invalid sCode...sCode=|%s|' %(self.sCode))
        ###
        return self.obj

    def instantiateBarObject(self): # 34 (CBAR)
        """
        Creates a stress/strain object if necessary
        """
        #print "starting a BAR stress/strain object"
        #print "    iSubcase = %s" %(self.iSubcase)
        #print "    sCode    = %s" %(self.sCode)
        #bits = self.parseStressCode()
        #if (self.iSubcase not in self.barStress) and (self.iSubcase not in self.barStrain):
        #    print "making new subcase..."

        if self.sCode in [0,1]:
            self.createTransientObject(self.barStress,barStressObject)
        elif self.sCode in [10,11]:
            self.createTransientObject(self.barStrain,barStrainObject)
        else:
            raise Exception('invalid sCode...sCode=|%s|' %(self.sCode))
        ###
        return self.obj

    def instantiateCompositePlateObject(self):  # 95,96,97,98 (CQUAD4, CQUAD8, CTRIA3, CTRIA6)
        """
        Creates a stress/strain object if necessary
        """
        #print "starting a COMPOSITE PLATE stress/strain object"
        #print "    iSubcase = %s" %(self.iSubcase)
        #print "    sCode    = %s" %(self.sCode)
        #bits = self.parseStressCode()
        #if (self.iSubcase not in self.compositePlateStress) and (self.iSubcase not in self.compositePlateStrain):
        #    print "making new subcase..."

        if self.sCode in [0,1,16,17]: # stress
            self.createTransientObject(self.compositePlateStress,compositePlateStressObject)
        elif self.sCode in [10,11,13,14,15,26,27,30,31]: # strain
            self.createTransientObject(self.compositePlateStrain,compositePlateStrainObject)
        else:
            raise Exception('invalid sCode...sCode=|%s|' %(self.sCode))
        ###
        return self.obj

    def instantiatePlateObject(self): # 74, 39, 144 (CTRIA3, CQUAD4)
        """
        Creates a stress/strain object if necessary
        """
        #print "starting a PLATE stress/strain object"
        #print "    iSubcase = %s" %(self.iSubcase)
        #print "    sCode    = %s" %(self.sCode)
        #bits = self.parseStressCode()
        #if (self.iSubcase not in self.plateStress) and (self.iSubcase not in self.plateStrain):
        #    print "making new subcase..."

        if self.sCode in [0,1]:
            self.createTransientObject(self.plateStress,plateStressObject)
        elif self.sCode in [10,11,15]:
            self.createTransientObject(self.plateStrain,plateStrainObject)
        else:
            raise Exception('invalid sCode...sCode=|%s|' %(self.sCode))
        ###
        print type(self.obj)
        return self.obj

    def instantiateSolidObject(self): # 39, 67, 68 (CTETRA, CPENTA, CHEXA)
        """
        Creates a stress/strain object if necessary
        """
        #print "starting a SOLID stress/strain object"
        #print "    iSubcase = %s" %(self.iSubcase)
        #print "    sCode    = %s" %(self.sCode)
        #bits = self.parseStressCode()
        #if (self.iSubcase not in self.solidStress) and (self.iSubcase not in self.solidStrain):
        #    print "making new subcase..."

        if self.sCode in [0,1]:
            self.createTransientObject(self.solidStress,solidStressObject)
        elif self.sCode in [10,11]:
            self.createTransientObject(self.solidStrain,solidStrainObject)
        else:
            raise Exception('invalid sCode...sCode=|%s|' %(self.sCode))
        ###
        return self.obj

