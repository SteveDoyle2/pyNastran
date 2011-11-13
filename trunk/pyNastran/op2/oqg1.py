import sys

class OQG1(object):
    def readTable_OQG1(self):
        ## OQG1
        word = self.readTableName(rewind=False) # OQG1
        print "word = |%r|" %(word)

        self.readMarkers([-1,7])
        ints = self.readIntBlock()
        print "*ints = ",ints

        self.readMarkers([-2,1,0]) # 7
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4
        ints = self.readIntBlock()
        print "*ints = ",ints

        self.readMarkers([-3,1,0])
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4

        data = self.getData(4)
        bufferSize, = unpack('i',data)
        print "bufferSize = ",bufferSize
        data = self.getData(4*50)
        aCode = self.getBlockIntEntry(data,1)
        print "aCode = ",aCode
        (analysisCode,deviceCode,tableCode,three,subcase) = self.parseAnalysisCode(data)


        word = self.readString(384)
        print "word = |%s|" %(word)
        self.readHollerith()
        
        self.readMarkers([-4,1,0])
        wordCount = self.getMarker()
        data = self.readBlock()
        #self.printBlock(data)

        iSubcase = 1 ## @todo temporary
        spcForcesObj = spcForcesObject(iSubcase)
        self.readScalars(deviceCode,data,spcForcesObj)

        self.readMarkers([-5,1,0])
        #print str(spcForcesObj)
    
