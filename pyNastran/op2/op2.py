from fortranFile import FortranFile
import os
import sys
import struct
from struct import unpack

class Op2(FortranFile):
    def __init__(self,infileName): 
        self.infilename = infileName
    
    def read(self):
        self.op2 = open(self.infilename,'rb')
 
        self.n = self.op2.tell()
        print "self.n = ",self.n
        self.readMarkers([3])
        ints = self.readIntBlock()
        print "*ints = ",ints
        self.readMarkers([7])

        word = self.readStringBlock() # Nastran Fort Tape ID Code - 
        print "word = |%r|" %(word)

        self.readMarkers([2])

        self.skip(4*4)
        
        self.readMarkers([-1,0,2])

        self.skip(4)
        hname = self.readString(8)
        print "hname = |%s|" %(hname)
        self.skip(4)

        self.readMarkers([-1,7])
        
        fields = self.readIntBlock()
        print "fields = ",fields

        self.readMarkers([-2,1,0,2])
        word = self.readStringBlock()
        self.readMarkers([-3,1,0,35])
        
        ## end geom1
        self.skip(4*69)
        
        ## start geom2
        word = self.readStringBlock()
        #print "word = |%r|" %(word)
        self.readMarkers([-1,7,])
        ints = self.readIntBlock()
        print "*ints = ",ints
        self.readMarkers([-2,1,0,2])
        word = self.readStringBlock()
        #print "word = |%r|" %(word)

        self.readMarkers([-3,1,0])
        self.skip(4*39)
        print "------------"

        ## GEOM3
        self.readMarkers([-5,1,0,0,2])
        word = self.readStringBlock()
        print "word = |%r|" %(word)
        self.readMarkers([-1,7])
        ints = self.readIntBlock()
        print "*ints = ",ints
        self.readMarkers([-2,1,0,2])
        word = self.readStringBlock()
        print "word = |%r|" %(word)

        self.readMarkers([-3,1,0,24])
        self.skip(4*26)

        self.readMarkers([-4,1,0,9])
        self.skip(4*11)

        self.readMarkers([-5,1,0,3])
        self.skip(4*5)

        print "------------"
        # GEOM4
        self.startTable([-6,1,0,0,2])
        self.readMarkers([-1,7])
        ints = self.readIntBlock()
        print "*ints = ",ints
        self.readMarkers([-2,1,0,2])
        word = self.readStringBlock()
        print "word = |%r|" %(word)

        self.readMarkers([-3,1,0,9])
        ints = self.readIntBlock()
        print "*ints = ",ints
        self.readMarkers([-4,1,0,6])
        ints = self.readIntBlock()
        print "*ints = ",ints
        self.readMarkers([-5,1,0,3])
        ints = self.readIntBlock()
        print "*ints = ",ints
        self.readMarkers([-6,1,0,0,2])
        print "------------"

        # EPT
        word = self.readStringBlock()
        print "word = |%r|" %(word)
        self.readMarkers([-1,7])
        ints = self.readIntBlock()
        print "*ints = ",ints
        self.readMarkers([-2,1,0,2])
        print "------------"
        word = self.readStringBlock()
        print "word = |%r|" %(word)
        self.readMarkers([-3,1,0,14])
        self.skip(4*16)
        self.readMarkers([-4,1,0,3])
        ints = self.readIntBlock()
        print "*ints = ",ints

        print "------------"
        self.readMarkers([-5,1,0,0,2,])
        word = self.readStringBlock()
        print "word = |%r|" %(word)
        self.readMarkers([-1,7])
        ints = self.readIntBlock()
        print "*ints = ",ints
        self.readMarkers([-2,1,0,2])

        print "------------"
        word = self.readStringBlock()
        print "word = |%r|" %(word)
        self.readMarkers([-3,1,0,15])
        ints = self.readIntBlock()
        print "*ints = ",ints
        self.readMarkers([-4,1,0,3])
        ints = self.readIntBlock()
        print "*ints = ",ints
        self.readMarkers([-5,1,0,0,2])
        print "------------"

        self.readTable_OQG1()
        self.readTable_OES1X1()
        
        
    def readTable_OQG1(self):
        ## OQG1
        word = self.readStringBlock()
        print "word = |%r|" %(word)
        self.readMarkers([-1,7])
        ints = self.readIntBlock()
        print "*ints = ",ints
        self.readMarkers([-2,1,0,7])
        ints = self.readIntBlock()
        print "*ints = ",ints
        self.readMarkers([-3,1,0,146])
        
        self.skip(4*51)
        word = self.readString(384)
        print "word = |%s|" %(word)
        #data = self.op2.read(384)
        #self.n += 384
        #print "word = |%s|" %(''.join(data))
        #self.skip(0)
        self.skip(4)  # weird hollerith
        
        self.readMarkers([-4,1,0,32])
        self.scan(4*34)

        self.readMarkers([-5,1,0,0,2])
        #word = self.readStringBlock()
        #print "word = |%r|" %(word)
        self.scan(16)
        self.readMarkers([-1,7])
        ints = self.readIntBlock()
        print "*ints = ",ints
        self.readMarkers([-2,1,0,7])
        ints = self.readIntBlock()
        print "*ints = ",ints
        self.readMarkers([-3,1,0,146])
        self.skip(4*51)
        word = self.readString(384)
        print "word = |%s|" %(word)
        self.skip(4)  # weird hollerith
        self.readMarkers([-4,1,0,32])
        self.skip(4*34)
        self.readMarkers([-5,1,0,0,2])

    def readTable_OES1X1(self):
        # OES1X1
        word = self.readStringBlock() # OES1X1
        print "word = |%r|" %(word)
        self.readMarkers([-1,7])
        ints = self.readIntBlock()
        print "*ints = ",ints
        self.readMarkers([-2,1,0,7])
        word = self.readStringBlock()  # OES1
        print "word = |%r|" %(word)
        self.readMarkers([-3,1,0,146])
        self.skip(4*51)
        word = self.readString(384)
        self.skip(4)  # weird hollerith
        print "word* = |%s|" %(word)
        self.readMarkers([-4,1,0,87])
        self.skip(16)

        floats = self.readFloats(4*85)
        print "*floats = ",floats,'\n'

        self.readMarkers([-5,1,0,])
        self.readMarkers([0,0,])

        
        

if __name__=='__main__':
    op2 = Op2()
    #op2.read('quad4.op2')
    op2.read('tria3.op2')

    print "done..."