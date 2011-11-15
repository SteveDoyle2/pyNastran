import sys
from struct import unpack
from op2_Objects import stressObject

class OES(object):
    def readTable_OES1X1(self):
        word = self.readTableName(rewind=False) # OES1X1
        print "word = |%r|" %(word)

        self.readMarkers([-1,7])
        print "****",self.op2.tell()
        data = self.readBlock()
        #self.printBlock(data)
        print "****",self.op2.tell()
        #assert self.op2.tell()==4880

        self.readMarkers([-2,1,0,7])
        word = self.readStringBlock()  # OES1
        print "word = |%r|" %(word)

        iTable = -3
        markerA=None; markerB=None
        while [markerA,markerB]!=[0,2]:
            self.readTable_OES_3(iTable)
            self.readTable_OES_4(iTable-1)
            print self.stress
            iTable -= 2

            n = self.n
            self.readMarkers([iTable,1,0],'OUGV1')
            
            #markerB = self.getMarker('B')
            #self.n-=24
            #self.op2.seek(self.n)
            #print "markerA=%s markerB=%s" %(markerA,markerB)
            self.printSection(120)


        self.readMarkers([-5,1,0,])
        #print "tell5 = ",self.op2.tell()
        self.printSection(100)
        self.readMarkers([0,0,])
        #print "end tell = ",self.op2.tell()

    def readTable_OES_3(self,iTable):
        self.readMarkers([iTable,1,0]) # 146
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4
        
        data = self.getData(4)
        bufferSize, = unpack('i',data)
        data = self.getData(4*51)

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

       
        word = self.readString(4*(63+32)) # subcase and label
        self.readHollerith()
        
        print "n4 = ",self.n

        print "word* = |%s|" %(word)
        
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

    def readTable_OES_4(self,iTable):
        self.readMarkers([iTable,1,0])
        markerA = 4
        
        j = 0
        while markerA>None:
            print "starting OES table 4..."
            self.readTable_OES_4_Data()
            print "finished reading stress table..."
            markerA = self.getMarker('A')
            self.n-=12
            self.op2.seek(self.n)
            
            self.n = self.op2.tell()
            print "***markerA = ",markerA
            
            if j>1:
                sys.exit('check...')
            j+=1


    def readTable_OES_4_Data(self):
        #self.printSection(100)

        assert self.op2.tell()==self.n,'tell=%s n=%s' %(self.op2.tell(),self.n)
        data = self.getData(16)
        assert self.op2.tell()==self.n,'tell=%s n=%s' %(self.op2.tell(),self.n)
        #self.printBlock(data)
        bufferWords, = unpack('i',data[4:8])

        print "*********************"
        #bufferWords = self.getMarker() # 87 - buffer
        print "bufferWords = ",bufferWords,bufferWords*4
        print "*elementType = ",self.elementType
        
        print "op2.tell=%s n=%s" %(self.op2.tell(),self.n)
        
        assert self.op2.tell()==self.n,'tell=%s n=%s' %(self.op2.tell(),self.n)

        self.data = self.getData(bufferWords*4)
        assert self.op2.tell()==self.n,'tell=%s n=%s' %(self.op2.tell(),self.n)
        print "self.n = ",self.n
        #self.printBlock(data)

        self.stress = stressObject(self.iSubcase)
        if self.elementType==144:
            self.CQUAD4_144(self.stress)  # 144
            print "found cquad_144"
        elif self.elementType==74:
            print "found ctria_74"
            self.CTRIA3_74(self.stress)  # 74
        elif self.elementType==39:
            self.CTETRA_39(self.stress)  # 39
        else:
            self.printSection(100)
            raise RuntimeError('elementType=%s -> %s is not supported' %(self.elementType,self.ElementType(self.elementType)))

        assert self.op2.tell()==self.n,'tell=%s n=%s' %(self.op2.tell(),self.n)
        # rods
        #if   self.elementType == 1:    # crod
        #elif self.elementType == 2:    # cbeam
        #elif self.elementType == 34:   # cbar

        #plate
        #elif self.elementType == 74:   # ctria3
        #elif self.elementType == 144:  # cquad_144
        #elif self.elementType == 33:   # cquad_33

        #solid
        #elif self.elementType == 39:  # ctetra
        #elif self.elementType == 67:  # chexa
        #elif self.elementType == 68:  # cpenta

        # composite plate
        #elif self.elementType == 95:  # quad4
        #elif self.elementType == 96:  # quad8
        #elif self.elementType == 97:  # tria3
        #elif self.elementType == 98:  # tria6
        

