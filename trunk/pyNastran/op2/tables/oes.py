import sys
import copy
from struct import unpack
from pyNastran.op2.resultObjects.oes_Objects import (
    rodStressObject,rodStrainObject,
    barStressObject,barStrainObject,
    plateStressObject,plateStrainObject,
    solidStressObject,
    compositePlateStressObject,compositePlateStrainObject)

class OES(object):
    """Table of stresses/strains"""
    def readTable_OES1(self):
        if self.makeOp2Debug:
            self.op2Debug.write("***OES table***\n")
        tableName = self.readTableName(rewind=False) # OES1X1
        print "tableName = |%r|" %(tableName)
        self.tableInit(tableName)

        self.readMarkers([-1,7])
        data = self.readBlock()
        #print self.printBlock(data)

        self.readMarkers([-2,1,0,7])
        word = self.readStringBlock()  # OES1
        print "OESword = |%r|" %(word)

        iTable = -3
        markerA=None; markerB=None
        self.readMarkers([iTable,1,0],'OES') # 146
        while [markerA,markerB]!=[0,2]:
            #self.markerStart = self.n
            print "reading Table 3...iTable=%s" %(iTable)
            self.readTable_OES_3(iTable)
            #if isDone:
            #    self.n = self.markerStart
            isBlockDone = self.readTable_OES_4(iTable-1)
            
            if self.sCode==11:
                pass
                #print self.plateStrain[self.iSubcase]
            else:
                pass
                #print self.barStress[self.iSubcase]
                #print self.plateStress[self.iSubcase]
                #print self.solidStress[self.iSubcase]
            iTable -= 2

            if isBlockDone:
                #print "***"
                print "iTable = ",iTable
                #print "$$$$"
                #self.n = self.markerStart
                #self.op2.seek(self.n)
                break
            ###

            n = self.n
            self.readMarkers([iTable,1,0],'OES')
            #print "i read the markers!!!"
        ###
        self.readMarkers([iTable,1,0],'OES')
        #self.printSection(100)
        
        self.deleteAttributes_OES()
        if self.makeOp2Debug:
            self.op2Debug.write("***end of OES table***\n")

    def deleteAttributes_OES(self):
        params = ['sCode','elementType','obj','markerStart','loadSet','formatCode','sCode','thermal',
                  'lsdvmn','mode','eign','modeCycle','freq','mode','eigr','eigi','dt']
        self.deleteAttributes(params)

    def readTable_OES_3(self,iTable):
        bufferWords = self.getBufferWords()
        
        data = self.getData(4)
        bufferSize, = unpack('i',data)
        if self.makeOp2Debug:
            self.op2Debug.write('bufferSize=|%s|\n' %(str(bufferSize)))

        data = self.getData(4*50)
        #self.printBlock(data) # on
        if self.makeOp2Debug:
            self.op2Debug.write('block3header\n')

        self.parseApproachCode(data)
        
        self.elementType = self.getValues(data,'i',3)
        self.loadSet    = self.getValues(data,'i',8)
        self.formatCode = self.getValues(data,'i',9)
        self.numWide    = self.getValues(data,'i',10)
        self.sCode      = self.getValues(data,'i',11)
        self.thermal    = self.getValues(data,'i',23) # 1 is heat transfer, 0 otherwise
        assert self.thermal==0

        print "loadset=%s formatCode=%s numWordsEntry=%s sCode=%s" %(self.loadSet,self.formatCode,self.numWide,self.sCode)
        print "self.thermal=%s" %(self.thermal)

        ## assuming tCode=1
        if self.approachCode==1:   # statics / displacement / heat flux
            self.lsdvmn = self.getValues(data,'i',5) ## load set number
        elif self.approachCode==2: # real eigenvalues
            self.mode      = self.getValues(data,'i',5) ## mode number
            self.eign      = self.getValues(data,'f',6) ## real eigenvalue
            self.modeCycle = self.getValues(data,'f',7) ## mode or cycle @todo confused on the type ???
        elif self.approachCode==3: # differential stiffness
            self.lsdvmn = self.getValues(data,'i',5) ## load set number
        elif self.approachCode==4: # differential stiffness
            self.lsdvmn = self.getValues(data,'i',5) ## load set number
        elif self.approachCode==5:   # frequency
            self.freq = self.getValues(data,'f',5) ## frequency

        elif self.approachCode==6: # transient
            self.dt = self.getValues(data,'f',5) ## time step
            print "DT(5)=%s" %(self.dt)
        elif self.approachCode==7: # pre-buckling
            self.lsdvmn = self.getValues(data,'i',5) ## load set
            print "LSDVMN(5)=%s" %(self.lsdvmn)
        elif self.approachCode==8: # post-buckling
            self.lsdvmn = self.getValues(data,'i',5) ## mode number
            self.eigr   = self.getValues(data,'f',6) ## real eigenvalue
            print "LSDVMN(5)=%s  EIGR(6)=%s" %(self.lsdvmn,self.eigr)
        elif self.approachCode==9: # complex eigenvalues
            self.mode   = self.getValues(data,'i',5) ## mode
            self.eigr   = self.getValues(data,'f',6) ## real eigenvalue
            self.eigi   = self.getValues(data,'f',7) ## imaginary eigenvalue
            print "LFTSFQ(5)=%s  EIGR(6)=%s  EIGI(7)=%s" %(self.lftsfq,self.eigr,self.eigi)
        elif self.approachCode==10: # nonlinear statics
            self.lftsfq = self.getValues(data,'f',5) ## load step
            print "LFTSFQ(5) = %s" %(self.lftsfq)
        elif self.approachCode==11: # old geometric nonlinear statics
            self.lsdvmn = self.getValues(data,'i',5)
            print "LSDVMN(5)=%s" %(self.lsdvmn)
        elif self.approachCode==12: # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            #self.lsdvmn = self.getValues(data,'i',5)
            self.dt = self.getValues(data,'f',5)  ## Time step ??? --> straight from DMAP
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
        
    def readTable_OES_4(self,iTable):
        #self.readMarkers([iTable,1,0])
        markerA = 4
        
        j = 0
        while markerA>None:
            self.markerStart = copy.deepcopy(self.n)
            #self.printSection(180)
            self.readMarkers([iTable,1,0])
            print "starting OES table 4..."
            isTable4Done,isBlockDone = self.readTable_OES_4_Data()
            if isTable4Done:
                print "done with OES4"
                self.n = self.markerStart
                self.op2.seek(self.n)
                break
            print "finished reading stress table..."
            markerA = self.getMarker('A')
            self.n-=12
            self.op2.seek(self.n)

            self.n = self.op2.tell()
            print "***markerA = ",markerA
            #sys.exit('check...end of table issueB...')
            #self.printSection(140)

            #print self.plateStress[self.iSubcase]
            
            iTable-=1
            #if j>10000:
            #    sys.exit('check...')
            j+=1
            #print "isBlockDone = ",isBlockDone
            #print "*******restarting table 4**********"
        print "isBlockDone = ",isBlockDone
        return isBlockDone

    def readTable_OES_4_Data(self):
        isTable4Done = False
        isBlockDone = False
        #self.printSection(100)

        data = self.getData(16)
        #print self.printBlock(data) # on
        #print "16 block..."
        #self.printBlock(data)
        #self.printBlock(data)
        bufferWords, = unpack('i',data[4:8])
        if self.makeOp2Debug:
            self.op2Debug.write('bufferWords=|%s|\n' %(str(bufferWords)))

        print "*********************"
        #bufferWords = self.getMarker() # 87 - buffer
        print "OES4 bufferWords = ",bufferWords,bufferWords*4
        #self.verifyBufferSize(bufferWords)
        
        isBlockDone = not(bufferWords)

        if bufferWords==146:  # table -4 is done, restarting table -3
            isTable4Done = True
            return isTable4Done,isBlockDone
        elif bufferWords==0:
            print "bufferWords 0 - done with Table4"
            isTable4Done = True
            isBlockDone = True
            return isTable4Done,isBlockDone

        print "*elementType = ",self.elementType
        #print "op2.tell=%s n=%s" %(self.op2.tell(),self.n)
        
        self.data = self.getData(bufferWords*4)
        #print "bufferWords = ",bufferWords
        if self.makeOp2Debug:
            self.op2Debug.write('reading big data block\n')
        #print "self.n = ",self.n
        #print 'reading big data block'
        #print self.printBlock(self.data)

        #stressStrainObj = self.instatiateStressStrainObject()
        self.readElementTable()
        return isTable4Done,isBlockDone

    def readElementTable(self):
        if self.elementType==1: # crod
            #print "    found crod_1"
            stressStrainObj = self.instantiateRodObject()
            self.CROD_1(stressStrainObj)
        #elif self.elementType == 2:   # cbeam
        #    print "    found cbeam_2"
        #    stressStrainObj = self.instantiateBeamObject()
        #    self.CBEAM_2(stressStrainObj)
        #elif self.elementType == 10:   # conrod
        #    print "    found cbeam_2"
        #    stressStrainObj = self.instantiateConrodObject()
        #    self.CONROD_10(stressStrainObj)
        elif self.elementType == 34:   # cbar
            #print "    found cbar_34"
            stressStrainObj = self.instantiateBarObject()
            self.CBAR_34(stressStrainObj)

        elif self.elementType==33: # cquad4_33
            #print "    found cquad_33"
            stressStrainObj = self.instantiatePlateObject()
            self.CQUAD4_33(stressStrainObj)
        elif self.elementType==74:  # ctria
            #print "    found ctria_74"
            stressStrainObj = self.instantiatePlateObject()
            self.CTRIA3_74(stressStrainObj) # ctria3
        elif self.elementType==144: # cquad4
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

        elif self.elementType in [95,96,97,98]: # CQUAD4, CQUAD8, CTRIA3, CTRIA6 (composite)
            #print "    found a 95/96/97 or 98!"
            stressStrainObj = self.instantiateCompositePlateObject()
            self.CQUAD4_95(stressStrainObj)
        else:
            self.printBlock(self.data[0:100])
            msg = 'elementType=%s -> %s is not supported' %(self.elementType,self.ElementType(self.elementType))
            raise RuntimeError(msg)

        #print stressStrainObj

        #elif self.elementType == 1:    # crod     (done)
        #elif self.elementType == 2:    # cbeam    (untested)
        #elif self.elementType == 3:    # ctube
        #elif self.elementType == 4:    # cshear
        #elif self.elementType == 10:   # conrod

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
        #elif self.elementType == 10:   # conrod
        #elif self.elementType == 34:   # cbar   (done)

        #springs
        #elif self.elementType == 11:   # celas1
        #elif self.elementType == 12:   # celas2
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
        #elif self.elementType == 85:   # tetra  (nonlinear)
        #elif self.elementType == 86:   # gap    (nonlinear)
        #elif self.elementType == 87:   # ctube  (nonlinear)
        #elif self.elementType == 88:   # tria3  (nonlinear) - same as quad4
        #elif self.elementType == 89:   # crod   (nonlinear)
        #elif self.elementType == 90:   # quad4  (nonlinear)
        #elif self.elementType == 91:   # cpenta (nonlinear)
        #elif self.elementType == 92:   # conrod (nonlinear)
        #elif self.elementType == 93:   # chexa  (nonlinear)
        #elif self.elementType == 94:   # cbeam  (nonlinear)
        
        # acoustic
        #elif self.elementType == 76:   # chexa  (acoustic)
        #elif self.elementType == 77:   # cpenta (acoustic)
        #elif self.elementType == 78:   # ctetra (acoustic)
        #elif self.elementType == 101:  # caabsf (acoustic)


    def instantiateRodObject(self): # 1 (CROD)
        """
        Creates a stress/strain object if necessary
        @todo I dont like the double return blocks, but it'll work...
        """
        print "starting a ROD stress/strain object"
        print "    iSubcase = ",self.iSubcase
        print "    sCode    = ",self.sCode
        #bits = self.parseStressCode()
        if (self.iSubcase not in self.rodStress) and (self.iSubcase not in self.rodStrain):
            print "making new subcase..."

        if self.sCode in [0,1]:
            #assert self.tableCode==0,'only REAL stress/strain is supported...tableCode=%s' %(self.tableCode)
            if self.iSubcase not in self.rodStress:
                self.rodStress[self.iSubcase] = rodStressObject(self.iSubcase)
            return self.rodStress[self.iSubcase]

        elif self.sCode in [10,11]:
            #assert self.tableCode==0,'only REAL stress/strain is supported...tableCode=%s' %(self.tableCode)
            if self.iSubcase not in self.barStrain:
                self.rodStrain[self.iSubcase] = rodStrainObject(self.iSubcase)
            return self.rodStrain[self.iSubcase]
        else:
            raise Exception('invalid sCode...sCode=|%s|' %(self.sCode))
        ###

    def instantiateBarObject(self): # 34 (CBAR)
        """
        Creates a stress/strain object if necessary
        @todo I dont like the double return blocks, but it'll work...
        """
        print "starting a BAR stress/strain object"
        print "    iSubcase = %s" %(self.iSubcase)
        print "    sCode    = %s" %(self.sCode)
        #bits = self.parseStressCode()
        if (self.iSubcase not in self.barStress) and (self.iSubcase not in self.barStrain):
            print "making new subcase..."

        if self.sCode in [0,1]:
            #assert self.tableCode==0,'only REAL stress/strain is supported...tableCode=%s' %(self.tableCode)
            if self.iSubcase not in self.barStress:
                self.barStress[self.iSubcase] = barStressObject(self.iSubcase)
            return self.barStress[self.iSubcase]

        elif self.sCode in [10,11]:
            #assert self.tableCode==0,'only REAL stress/strain is supported...tableCode=%s' %(self.tableCode)
            if self.iSubcase not in self.barStrain:
                self.barStrain[self.iSubcase] = barStrainObject(self.iSubcase)
            return self.barStrain[self.iSubcase]
        else:
            raise Exception('invalid sCode...sCode=|%s|' %(self.sCode))
        ###

    def instantiateCompositePlateObject(self):  # 95,96,97,98 (CQUAD4, CQUAD8, CTRIA3, CTRIA6)
        """
        Creates a stress/strain object if necessary
        @todo I dont like the double return blocks, but it'll work...
        """
        print "starting a COMPOSITE PLATE stress/strain object"
        print "    iSubcase = %s" %(self.iSubcase)
        print "    sCode    = %s" %(self.sCode)
        #bits = self.parseStressCode()
        if (self.iSubcase not in self.compositePlateStress) and (self.iSubcase not in self.compositePlateStrain):
            print "making new subcase..."

        if self.sCode in [0,1,16,17]: # stress
            #assert self.sortCode==0,'only REAL stress/strain is supported...tableCode=%s' %(self.tableCode)
            if self.iSubcase not in self.compositePlateStress:
                self.compositePlateStress[self.iSubcase] = compositePlateStressObject(self.iSubcase)
            return self.compositePlateStress[self.iSubcase]

        elif self.sCode in [10,11,13,14,15,26,27,30,31]: # strain
            #assert self.sortCode==0,'only REAL stress/strain is supported...tableCode=%s' %(self.tableCode)
            if self.iSubcase not in self.plateStrain:
                self.compositePlateStrain[self.iSubcase] = compositePlateStrainObject(self.iSubcase)
            return self.compositePlateStrain[self.iSubcase]
        else:
            raise Exception('invalid sCode...sCode=|%s|' %(self.sCode))
        ###

    def instantiatePlateObject(self): # 74, 144 (CTRIA3, CQUAD4)
        """
        Creates a stress/strain object if necessary
        @todo I dont like the double return blocks, but it'll work...
        """
        print "starting a PLATE stress/strain object"
        print "    iSubcase = %s" %(self.iSubcase)
        print "    sCode    = %s" %(self.sCode)
        #bits = self.parseStressCode()
        if (self.iSubcase not in self.plateStress) and (self.iSubcase not in self.plateStrain):
            print "making new subcase..."

        if self.sCode in [0,1]:
            #assert self.tableCode==0,'only REAL stress/strain is supported...tableCode=%s' %(self.tableCode)
            if self.iSubcase not in self.plateStress:
                self.plateStress[self.iSubcase] = plateStressObject(self.iSubcase)
            return self.plateStress[self.iSubcase]

        elif self.sCode in [10,11,15]:
            #assert self.tableCode==0,'only REAL stress/strain is supported...tableCode=%s' %(self.tableCode)
            if self.iSubcase not in self.plateStrain:
                self.plateStrain[self.iSubcase] = plateStrainObject(self.iSubcase)
            return self.plateStrain[self.iSubcase]
        else:
            raise Exception('invalid sCode...sCode=|%s|' %(self.sCode))
        ###

    def instantiateSolidObject(self): # 39, 67, 68 (CTETRA, CPENTA, CHEXA)
        """
        Creates a stress/strain object if necessary
        @todo I dont like the double return blocks, but it'll work...
        """
        print "starting a SOLID stress/strain object"
        print "    iSubcase = %s" %(self.iSubcase)
        print "    sCode    = %s" %(self.sCode)
        #bits = self.parseStressCode()
        if (self.iSubcase not in self.solidStress) and (self.iSubcase not in self.solidStrain):
            print "making new subcase..."

        if self.sCode in [0,1]:
            #assert self.tableCode==0,'only REAL stress/strain is supported...tableCode=%s' %(self.tableCode)
            if self.iSubcase not in self.solidStress:
                print "created new solidObject for iSubcase=%s" %(self.iSubcase)
                self.solidStress[self.iSubcase] = solidStressObject(self.iSubcase)
            return self.solidStress[self.iSubcase]

        elif self.sCode in [10,11]:
            #assert self.tableCode==0,'only REAL stress/strain is supported...tableCode=%s' %(self.tableCode)
            if self.iSubcase not in self.solidStrain:
                self.solidStrain[self.iSubcase] = solidStrainObject(self.iSubcase)
            return self.solidStrain[self.iSubcase]
        else:
            raise Exception('invalid sCode...sCode=|%s|' %(self.sCode))
        ###

