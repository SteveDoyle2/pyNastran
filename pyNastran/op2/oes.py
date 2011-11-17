import sys
import copy
from struct import unpack
from op2_Objects import stressObject,strainObject

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
            isDone = self.readTable_OES_4(iTable-1)
            
            if self.sCode==11:
                print self.strain[self.iSubcase]
            else:
                print self.stress[self.iSubcase]
            iTable -= 2

            n = self.n

            if isDone:
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

    def readTable_OES_4(self,iTable):
        #self.readMarkers([iTable,1,0])
        markerA = 4
        
        j = 0
        while markerA>None:
            self.markerStart = copy.deepcopy(self.n)
            #self.printSection(180)
            self.readMarkers([iTable,1,0])
            print "starting OES table 4..."
            isDone,isOesDone = self.readTable_OES_4_Data()
            if isDone:
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

            #print self.stress[self.iSubcase]
            
            iTable-=1
            if j>10000:
            #    sys.exit('check...')
            #j+=1
            print "isOesDone = ",isOesDone
            
        print "isOesDone = ",isOesDone
        return isOesDone
            

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
        
    def instatiateStressStrainObject(self):
        """
        Creates a stress/strain object if necessary
        @todo I dont like the double return blocks, but it'll work...
        """
        print "starting a stress/strain object"
        print "    iSubcase = ",self.iSubcase
        print "    sCode    = ",self.sCode
        #bits = self.parseStressCode()
        if (self.iSubcase not in self.stress) and (self.iSubcase not in self.strain):
            print "making new subcase..."

        if (self.sCode==0 or self.sCode==1):
            if self.iSubcase not in self.stress:
                self.stress[self.iSubcase] = stressObject(self.iSubcase)
            return self.stress[self.iSubcase]

        elif (self.sCode==10 or self.sCode==11):
            if self.iSubcase not in self.strain:
                self.strain[self.iSubcase] = strainObject(self.iSubcase)
            return self.strain[self.iSubcase]
        else:
            raise Exception('invalid sCode...sCode=|%s|' %(self.sCode))
        ###

    def verifyBufferSize(self,bufferWords):
        assert bufferWords>0,self.printSection(220)

    def readTable_OES_4_Data(self):
        isDone = False
        isOesDone = False
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
        
        isOesDone = not(bufferWords)

        if bufferWords==146:  # table -4 is done, restarting table -3
            isDone = True
            return isDone,isOesDone
        elif bufferWords==0:
            #print "bufferWords 0"
            isDone = True
            isOesDone = True
            return isDone,isOesDone
        print "*elementType = ",self.elementType
        
        #print "op2.tell=%s n=%s" %(self.op2.tell(),self.n)
        
        
        self.data = self.getData(bufferWords*4)
        #print "self.n = ",self.n
        #self.printBlock(data)

        stressStrainObj = self.instatiateStressStrainObject()
        if self.elementType==144:
            print "    found cquad_144"
            self.CQUAD4_144(stressStrainObj)  # 144
        elif self.elementType==74:
            print "    found ctria_74"
            self.CTRIA3_74(stressStrainObj)  # 74
        elif self.elementType==39:
            print "    found ctetra_39"
            self.CTETRA_39(stressStrainObj)  # 39
        else:
            self.printSection(100)
            raise RuntimeError('elementType=%s -> %s is not supported' %(self.elementType,self.ElementType(self.elementType)))

        #assert self.op2.tell()==self.n,'tell=%s n=%s' %(self.op2.tell(),self.n)
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
        return isDone,isOesDone

