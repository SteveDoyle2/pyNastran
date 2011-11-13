import sys

class OUGV1(object):
    def readTable_OUGV1(self):
        ## OUGV1
        word = self.readTableName(rewind=False) # OUGV1
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
        print "bufferWords = ",bufferWords,bufferWords*4,'\n'

        data = self.getData(4)
        bufferSize, = unpack('i',data)
        print "bufferSize = ",bufferSize
        data = self.getData(4*50)
        aCode = self.getBlockIntEntry(data,1)
        print "aCode = ",aCode
        (analysisCode,deviceCode,tableCode,three,subcase) = self.parseAnalysisCode(data)
        #self.printBlock(data)

        word = self.readString(384)
        print "word = |%s|" %(word)
        self.readHollerith()

        self.readMarkers([-4,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()
        
        iSubcase = 1 ## @todo temporary
        dispObj = displacementObject(iSubcase)
        self.readScalars(deviceCode,data,dispObj)
        #print str(dispObj)

        self.readMarkers([-5,1,0])
        return
        #self.printSection(100)
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-6,1,0])
        bufferWords = self.getMarker()
        #self.printSection(100)
        data = self.readBlock()
        #sys.exit('ougv1 - stop')
        self.readMarkers([-7,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-8,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-9,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-10,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-11,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-12,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-13,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-14,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-15,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-16,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-17,1,0])
        bufferWords = self.getMarker()
        bufferWords = self.getMarker()
        data = self.readBlock()



        self.readMarkers([-1,7])
        ints = self.readIntBlock()

        self.readMarkers([-2,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-3,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-4,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()


        self.readMarkers([-5,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()


        self.readMarkers([-6,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()


        self.readMarkers([-7,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()


        self.readMarkers([-8,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-9,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-10,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-11,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()


        self.readMarkers([-12,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-13,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()


        self.readMarkers([-14,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-15,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-16,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-17,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-18,1,0])
        bufferWords = self.getMarker()
        eid = unpack('i',data[0:4])

        data = self.readBlock()
        print "---------"
        self.printBlock(data)
        
        
        

        #self.readMarkers([-1,7])
        #bufferWords = self.getMarker()
        #data = self.readBlock()
        #self.printSection(160)


        #assert self.op2.tell()==4780,self.op2.tell()
        #sys.exit('end of displacements')
        #self.readTable_OEF1X()
        sys.exit('end of ougv1')

    def readTable_OEF1X(self):
        ## OEF1X
        word = self.readTableName(rewind=False) # OEF1X
        print "word = |%r|" %(word)

        self.readMarkers([-1,7])
        ints = self.readIntBlock()
        print "*ints = ",ints

