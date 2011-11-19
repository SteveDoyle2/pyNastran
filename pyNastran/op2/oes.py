import sys
import copy
from struct import unpack
from op2_Objects import (rodStressObject,rodStrainObject,
                        barStressObject,
                        plateStressObject,plateStrainObject,
                        solidStressObject)

class OES(object):
    def readTable_OES1X1(self):
        tableName = self.readTableName(rewind=False) # OES1X1
        print "OES table = |%r|" %(tableName)

        self.readMarkers([-1,7])
        #print "****",self.op2.tell()
        data = self.readBlock()
        #self.printBlock(data)
        #print "****",self.op2.tell()
        #assert self.op2.tell()==4880

        self.readMarkers([-2,1,0,7])
        word = self.readStringBlock()  # OES1
        print "OESword = |%r|" %(word)

        iTable = -3
        markerA=None; markerB=None
        self.readMarkers([iTable,1,0]) # 146
        while [markerA,markerB]!=[0,2]:
            #self.markerStart = self.n
            print "reading Table 3...iTable=%s" %(iTable)
            self.readTable_OES_3(iTable)
            #if isDone:
            #    self.n = self.markerStart
            isBlockDone = self.readTable_OES_4(iTable-1)
            
            if self.sCode==11:
                print self.plateStrain[self.iSubcase]
            else:
                pass
                #print self.barStress[self.iSubcase]
                #print self.plateStress[self.iSubcase]
                #print self.solidStress[self.iSubcase]
            iTable -= 2

            n = self.n

            if isBlockDone:
                print "***"
                print "iTable = ",iTable
                print "$$$$"
                #self.n = self.markerStart
                #self.op2.seek(self.n)
                break
            ###

            self.readMarkers([iTable,1,0],'OES')
            print "i read the markers!!!"
            #markerB = self.getMarker('B')
            #self.n-=24
            #self.op2.seek(self.n)
            #print "markerA=%s markerB=%s" %(markerA,markerB)
        ###
        self.readMarkers([iTable,1,0],'OES')
        #print "tell5 = ",self.op2.tell()
        #self.printSection(100)
        #self.readMarkers([0,])
        #self.readMarker() # 0 or 2
        
        
        print "exiting stress/strain table..."
        #self.printSection(180)
        #self.readIntBlock()

        #word = self.readTableName(rewind=True) # OES1X1
        #print "OES table = |%r|" %(word)

        #print "end tell = ",self.op2.tell()


    def getBufferWords(self):
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4
        assert bufferWords >0
        return bufferWords

    def readTable_OES_3(self,iTable):
        bufferWords = self.getBufferWords()
        
        data = self.getData(4)
        bufferSize, = unpack('i',data)
        data = self.getData(4*50)

        nWide = self.getBlockIntEntry(data,10)
        #print "nWide = ",nWide
        thermal = self.getBlockIntEntry(data,21)

        (tCode,self.elementType,self.iSubcase) = self.parseApproachCode(data)
        data = data[16:]
        
        (word5,word6,word7) = unpack('iii',data[:12]) # depends on analysisCode,tCode
        print "word5=%s word6=%s word7=%s" %(word5,word6,word7)
        data = data[12:]

        (loadset,fcode,numWordsEntry,self.sCode) = unpack('iiii',data[:16])
        print "loadset=%s fcode=%s numWordsEntry=%s sCode=%s" %(loadset,fcode,numWordsEntry,self.sCode)
        print "thermal=%s" %(thermal)
        data = data[16:]
        assert thermal==0 # 1 is heat transfer

       
        word = self.readString(4*(63+33)) # title, subtitle, and label
        Title    = word[0:128]
        Subtitle = word[128:256]
        Label    = word[256:]
        print "Title    %s |%s|" %(len(Title   ),Title.strip())
        print "Subtitle %s |%s|" %(len(Subtitle),Subtitle.strip())
        print "Label    %s |%s|" %(len(Label   ),Label.strip())
        self.readHollerith()
        
        print "n4 = ",self.n

        #print "titleSubtitleLabel = |%s|" %(word)
        
    def isStatics(self):
        if self.approachCode==1 and self.tableCode==1:
            return True
        return False

    def isTransient(self):
        if self.approachCode==6 and self.tableCode==1:
            return True
        return False

    def isNonlinearStatics(self):
        if self.approachCode==10 and self.tableCode==1:
            return True
        return False

    def isStress(self):
        if self.sCode==1:
            return True
        return False

    def isStrain(self):
        if self.sCode==0:
            return True
        return False
    
    def isStaticStress(self):
        if self.isStatics() and self.isStress():
            return True
        return False

    def isStaticStrain(self):
        if self.isStatics() and self.isStrain():
            return True
        return False

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
        
    def verifyBufferSize(self,bufferWords):
        assert bufferWords>0,self.printSection(220)

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
            #self.printSection(140)

            #print self.plateStress[self.iSubcase]
            
            iTable-=1
            if j>10000:
                sys.exit('check...')
            j+=1
            print "isBlockDone = ",isBlockDone
            
        print "isBlockDone = ",isBlockDone
        return isBlockDone

    def readTable_OES_4_Data(self):
        isTable4Done = False
        isBlockDone = False
        #self.printSection(100)

        data = self.getData(16)
        #print "16 block..."
        #self.printBlock(data)
        #self.printBlock(data)
        bufferWords, = unpack('i',data[4:8])

        print "*********************"
        #bufferWords = self.getMarker() # 87 - buffer
        print "OES4 bufferWords = ",bufferWords,bufferWords*4
        #self.verifyBufferSize(bufferWords)
        
        isBlockDone = not(bufferWords)

        if bufferWords==146:  # table -4 is done, restarting table -3
            isTable4Done = True
            return isTable4Done,isBlockDone
        elif bufferWords==0:
            #print "bufferWords 0"
            isTable4Done = True
            isBlockDone = True
            return isTable4Done,isBlockDone

        print "*elementType = ",self.elementType
        
        #print "op2.tell=%s n=%s" %(self.op2.tell(),self.n)
        
        self.data = self.getData(bufferWords*4)
        #print "self.n = ",self.n
        #self.printBlock(data)

        #stressStrainObj = self.instatiateStressStrainObject()
        self.readElementTable()

        return isTable4Done,isBlockDone

    def readElementTable(self):
        if self.elementType==1: # crod
            print "    found crod_1"
            stressStrainObj = self.instantiateRodObject()
            self.CROD_1(stressStrainObj)

        elif self.elementType == 34:   # cbar
            print "    found cbar_34"
            stressStrainObj = self.instantiateBarObject()
            self.CBAR_34(stressStrainObj)

        elif self.elementType==144: # cquad4
            print "    found cquad_144"
            stressStrainObj = self.instantiatePlateObject()
            self.CQUAD4_144(stressStrainObj)
        elif self.elementType==74:  # ctria
            print "    found ctria_74"
            stressStrainObj = self.instantiatePlateObject()
            self.CTRIA3_74(stressStrainObj) # ctria3


        #elif self.elementType==39:
        #    print "    found ctetra_39"
        #    stressStrainObj = self.instantiateSolidObject()
        #    self.CTETRA_39(stressStrainObj)
        elif self.elementType in [39,67,68]:   # ctetra/chexa/cpenta
            print "    found hexa_67 / cpenta_68"
            stressStrainObj = self.instantiateSolidObject()
            self.CHEXA_67(stressStrainObj)
        else:
            self.printSection(100)
            msg = 'elementType=%s -> %s is not supported' %(self.elementType,self.ElementType(self.elementType))
            raise RuntimeError(msg)


        #elif self.elementType == 4:    # cshear

        #elif self.elementType == 33:   # cquad4_33
        #elif self.elementType == 34:   # cbar
        #elif self.elementType == 35:   # cconeax
        #elif self.elementType == 38:   # cgap
        #elif self.elementType == 39:   # ctetra
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
        #elif self.elementType == 67:   # chexa
        #elif self.elementType == 68:   # cpenta
        #elif self.elementType == 69:   # cbend
        #elif self.elementType == 70:   # ctriar
        #elif self.elementType == 74:   # ctria3
        #elif self.elementType == 75:   # ctria6
        #elif self.elementType == 82:   # cquadr
        #elif self.elementType == 100:  # cbar w/ cbarao or pload1

        #elif self.elementType == 102:  # cbush

        #elif self.elementType == 144:  # cquad_144 - corner stresses

        # rods
        #if   self.elementType == 1:    # crod
        #elif self.elementType == 2:    # cbeam
        #elif self.elementType == 3:    # ctube
        #elif self.elementType == 10:   # conrod

        #springs
        #elif self.elementType == 11:   # celas1
        #elif self.elementType == 12:   # celas2
        #elif self.elementType == 13:   # celas3
        #elif self.elementType == 14:   # celas4
        
        #plate
        #elif self.elementType == 74:   # ctria3
        #elif self.elementType == 144:  # cquad_144
        #elif self.elementType == 33:   # cquad_33

        #solid (???)
        #elif self.elementType == 39:  # ctetra
        #elif self.elementType == 67:  # chexa
        #elif self.elementType == 68:  # cpenta

        # composite plate
        #elif self.elementType == 94:   # quad4 (composite)
        #elif self.elementType == 95:   # quad8 (composite)
        #elif self.elementType == 97:   # tria3 (composite) - same as quad4 composite
        #elif self.elementType == 98:   # tria6 (composite) - same as quad8 composite ??? said quad4

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
        
        # acoustic
        #elif self.elementType == 76:   # chexa  (acoustic)
        #elif self.elementType == 77:   # cpenta (acoustic)
        #elif self.elementType == 78:   # ctetra (acoustic)
        #elif self.elementType == 101:  # caabsf (acoustic)


    def instantiateStressStrainObject(self):
        if self.elementType==34:
            return self.instatiateBarObject()
        elif self.elementType==67 or self.elementType==68:
            return self.instatiateSolidObject()
        elif self.elementType==144:
            return self.instatiatePlateObject()
        else:
            msg = 'elementType=%s -> %s is not supported' %(self.elementType,self.ElementType(self.elementType))
            raise Exception(msg)
        ###

    def instantiateRodObject(self):
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

        if (self.sCode==0 or self.sCode==1):
            #assert self.tableCode==0,'only REAL stress/strain is supported...tableCode=%s' %(self.tableCode)
            if self.iSubcase not in self.rodStress:
                self.rodStress[self.iSubcase] = rodStressObject(self.iSubcase)
            return self.rodStress[self.iSubcase]

        elif (self.sCode==10 or self.sCode==11):
            #assert self.tableCode==0,'only REAL stress/strain is supported...tableCode=%s' %(self.tableCode)
            if self.iSubcase not in self.barStrain:
                self.rodStrain[self.iSubcase] = rodStrainObject(self.iSubcase)
            return self.rodStrain[self.iSubcase]
        else:
            raise Exception('invalid sCode...sCode=|%s|' %(self.sCode))
        ###

    def instantiateBarObject(self):
        """
        Creates a stress/strain object if necessary
        @todo I dont like the double return blocks, but it'll work...
        """
        print "starting a BAR stress/strain object"
        print "    iSubcase = ",self.iSubcase
        print "    sCode    = ",self.sCode
        #bits = self.parseStressCode()
        if (self.iSubcase not in self.barStress) and (self.iSubcase not in self.barStrain):
            print "making new subcase..."

        if (self.sCode==0 or self.sCode==1):
            #assert self.tableCode==0,'only REAL stress/strain is supported...tableCode=%s' %(self.tableCode)
            if self.iSubcase not in self.barStress:
                self.barStress[self.iSubcase] = barStressObject(self.iSubcase)
            return self.barStress[self.iSubcase]

        elif (self.sCode==10 or self.sCode==11):
            #assert self.tableCode==0,'only REAL stress/strain is supported...tableCode=%s' %(self.tableCode)
            if self.iSubcase not in self.barStrain:
                self.barStrain[self.iSubcase] = barStrainObject(self.iSubcase)
            return self.barStrain[self.iSubcase]
        else:
            raise Exception('invalid sCode...sCode=|%s|' %(self.sCode))
        ###

    def instantiatePlateObject(self):
        """
        Creates a stress/strain object if necessary
        @todo I dont like the double return blocks, but it'll work...
        """
        print "starting a PLATE stress/strain object"
        print "    iSubcase = ",self.iSubcase
        print "    sCode    = ",self.sCode
        #bits = self.parseStressCode()
        if (self.iSubcase not in self.plateStress) and (self.iSubcase not in self.plateStrain):
            print "making new subcase..."

        if (self.sCode==0 or self.sCode==1):
            #assert self.tableCode==0,'only REAL stress/strain is supported...tableCode=%s' %(self.tableCode)
            if self.iSubcase not in self.plateStress:
                self.plateStress[self.iSubcase] = plateStressObject(self.iSubcase)
            return self.plateStress[self.iSubcase]

        elif (self.sCode==10 or self.sCode==11):
            #assert self.tableCode==0,'only REAL stress/strain is supported...tableCode=%s' %(self.tableCode)
            if self.iSubcase not in self.plateStrain:
                self.plateStrain[self.iSubcase] = plateStrainObject(self.iSubcase)
            return self.plateStrain[self.iSubcase]
        else:
            raise Exception('invalid sCode...sCode=|%s|' %(self.sCode))
        ###

    def instantiateSolidObject(self):
        """
        Creates a stress/strain object if necessary
        @todo I dont like the double return blocks, but it'll work...
        """
        print "starting a SOLID stress/strain object"
        print "    iSubcase = ",self.iSubcase
        print "    sCode    = ",self.sCode
        #bits = self.parseStressCode()
        if (self.iSubcase not in self.solidStress) and (self.iSubcase not in self.solidStrain):
            print "making new subcase..."

        if (self.sCode==0 or self.sCode==1):
            #assert self.tableCode==0,'only REAL stress/strain is supported...tableCode=%s' %(self.tableCode)
            if self.iSubcase not in self.solidStress:
                print "created new solidObject for iSubcase=%s" %(self.iSubcase)
                self.solidStress[self.iSubcase] = solidStressObject(self.iSubcase)
            return self.solidStress[self.iSubcase]

        elif (self.sCode==10 or self.sCode==11):
            #assert self.tableCode==0,'only REAL stress/strain is supported...tableCode=%s' %(self.tableCode)
            if self.iSubcase not in self.solidStrain:
                self.solidStrain[self.iSubcase] = solidStrainObject(self.iSubcase)
            return self.solidStrain[self.iSubcase]
        else:
            raise Exception('invalid sCode...sCode=|%s|' %(self.sCode))
        ###
