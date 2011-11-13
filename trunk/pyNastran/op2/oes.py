import sys

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

        self.readMarkers([-3,1,0]) # 146
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4
        
        data = self.getData(4)
        bufferSize, = unpack('i',data)
        data = self.getData(4*51)

        nWide = self.getBlockIntEntry(data,10)
        #print "nWide = ",nWide
        thermal = self.getBlockIntEntry(data,21)

        (analysisCode,deviceCode,tCode,elementType,iSubcase) = self.parseAnalysisCode(data)
        data = data[16:]
        
        (word5,word6,word7) = unpack('iii',data[:12]) # depends on analysisCode,tCode
        print "word5=%s word6=%s word7=%s" %(word5,word6,word7)
        data = data[12:]

        (loadset,fcode,numWordsEntry,sCode) = unpack('iiii',data[:16])
        print "loadset=%s fcode=%s numWordsEntry=%s sCode=%s" %(loadset,fcode,numWordsEntry,sCode)
        print "thermal=%s" %(thermal)
        data = data[16:]

       
        word = self.readString(4*(63+32)) # subcase and label
        self.readHollerith()
        
        print "n4 = ",self.n

        print "word* = |%s|" %(word)
        self.readMarkers([-4,1,0])
        #self.printSection(100)

        data = self.getData(16)
        #self.printBlock(data)
        bufferWords, = unpack('i',data[4:8])

        print "*********************"
        #bufferWords = self.getMarker() # 87 - buffer
        print "bufferWords = ",bufferWords,bufferWords*4
        print "*elementType = ",elementType
        
        print "op2.tell=%s n=%s" %(self.op2.tell(),self.n)
        #assert self.op2.tell()==5656

        data = self.getData(bufferWords*4)
        #self.printBlock(data)

        stress = stressObject(1) # @todo dummy for now...
        if elementType==144:
            self.CQUAD4_144(data,stress)  # 144
            print "found cquad)144"
        elif elementType==74:
            print "found ctria_74"
            self.CTRIA3_74(data,stress)  # 74
        elif elementType==39:
            self.CTETRA_39(data,stress)  # 39
        else:
            raise RuntimeError('elementType=%s -> %s is not supported' %(elementType,self.ElementType(elementType)))

        print stress

        self.readMarkers([-5,1,0,])
        #print "tell5 = ",self.op2.tell()
        self.readMarkers([0,0,])
        #print "end tell = ",self.op2.tell()

