import os
import sys
import struct
from struct import unpack

class GeometryTables(object):

    def readTable_Geom1(self):
        word = self.readTableName(rewind=False) # GEOM1
        #print "*word = |%r|" %(word)

        self.readMarkers([-1,7])
        fields = self.readIntBlock()
        print "fields = ",fields

        self.readMarkers([-2,1,0]) # 2
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4
        word = self.readStringBlock()

        self.readMarkers([-3,1,0])
        marker = self.getMarker() # 35,315,1571
        print "marker = ",marker
        #self.printTableCode(marker)
        
        ints = self.readIntBlock()
        #print "*ints = ",ints, len(ints)

        #while ints:  ## @todo is this correct???
        #    coord1 = ints[:6]
        #    ints = ints[6:]
        #    print "coord1 = ",coord1

        self.readMarkers([-4,1,0])  #3
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4
        ints = self.readIntBlock()
        print "*ints = ",ints

        self.readMarkers([-5,1,0])
        #assert self.op2.tell()==584,self.op2.tell()

    def readTable_Geom2(self):
        word = self.readTableName(rewind=False) # GEOM2
        print "word = |%r|" %(word)

        self.readMarkers([-1,7])
        ints = self.readIntBlock()
        print "*ints = ",ints
        self.readMarkers([-2,1,0]) #2
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4
        word = self.readStringBlock()
        #print "word = |%r|" %(word)

        self.readMarkers([-3,1,0]) # 17
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4

        #marker = self.getMarker() # 17
        #print "marker = ",marker
        #self.printTableCode(marker)

        ints = self.readIntBlock()
        print "*ints = ",ints
        self.readMarkers([-4,1,0]) # 3
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4
        ints = self.readIntBlock()
        print "*ints = ",ints
        self.readMarkers([-5,1,0])
        #assert self.op2.tell()==976,self.op2.tell()

    def readTable_Geom3(self):
        ## GEOM3
        word = self.readTableName(rewind=False) # GEOM3
        print "word = |%r|" %(word)

        self.readMarkers([-1,7])
        ints = self.readIntBlock()
        print "*ints = ",ints

        self.readMarkers([-2,1,0])
        bufferWords = self.getMarker() # 2
        print "bufferWords = ",bufferWords,bufferWords*4
        word = self.readStringBlock()
        print "word = |%r|" %(word)

        self.readMarkers([-3,1,0]) # 24
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4
        self.skip(4*26)

        self.readMarkers([-4,1,0]) # 9
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4
        self.skip(4*11)

        self.readMarkers([-5,1,0]) # 3
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4
        self.skip(4*5)

        self.readMarkers([-6,1,0])
        assert self.op2.tell()==1488,self.op2.tell()
        

    def readTable_Geom4(self):
        # GEOM4
        word = self.readTableName(rewind=False) # GEOM4

        self.readMarkers([-1,7])
        ints = self.readIntBlock()
        print "*ints = ",ints

        self.readMarkers([-2,1,0]) # 2
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4
        word = self.readStringBlock()
        print "word = |%r|" %(word)

        self.readMarkers([-3,1,0]) # 9
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4
        ints = self.readIntBlock()
        print "*ints = ",ints

        self.readMarkers([-4,1,0]) # 6
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4
        ints = self.readIntBlock()
        print "*ints = ",ints

        self.readMarkers([-5,1,0]) # 3
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4
        ints = self.readIntBlock()
        print "*ints = ",ints

        self.readMarkers([-6,1,0])
        print "------------"

    def readTable_EPT(self):
        # EPT
        word = self.readTableName(rewind=False) # EPT
        print "word = |%r|" %(word)

        self.readMarkers([-1,7])
        ints = self.readIntBlock()
        print "*ints = ",ints

        self.readMarkers([-2,1,0]) # 2
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4
        print "------------"
        word = self.readStringBlock()
        print "word = |%r|" %(word)

        self.readMarkers([-3,1,0]) # 14
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4
        self.skip(4*16)

        self.readMarkers([-4,1,0]) # 3
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4
        ints = self.readIntBlock()
        print "*ints = ",ints

        print "------------"
        self.readMarkers([-5,1,0])

    def readTable_MPTS(self):
        ## MPTS
        word = self.readTableName(rewind=False) # MPTS
        print "word = |%r|" %(word)

        self.readMarkers([-1,7])
        ints = self.readIntBlock()
        print "*ints = ",ints

        self.readMarkers([-2,1,0]) # 2
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4

        print "------------"
        word = self.readStringBlock()
        print "word = |%r|" %(word)

        self.readMarkers([-3,1,0]) # 15
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4
        ints = self.readIntBlock()
        print "*ints = ",ints

        self.readMarkers([-4,1,0]) # 3
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4
        ints = self.readIntBlock()
        print "*ints = ",ints

        self.readMarkers([-5,1,0])
        print "------------"
        assert self.op2.tell()==2692,self.op2.tell()


