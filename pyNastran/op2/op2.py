from fortranFile import FortranFile
from op2Codes import Op2Codes
import os
import sys
import struct
from struct import unpack

class Op2(FortranFile,Op2Codes):
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
        self.readTable_Geom1()
        self.readTable_Geom2()
        self.readTable_Geom3()
        self.readTable_Geom4()
        self.readTable_OQG1()
        self.readTable_OES1X1()
        
        
    def readTable_Geom1(self):
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

        self.readMarkers([-3,1,0])
        marker = self.getMarker() # 35,1571
        print "marker = ",marker
        self.printTableCode(marker)
        ints = self.readIntBlock()
        print "*ints = ",ints, len(ints)

        while ints:
            coord1 = ints[:6]
            ints = ints[6:]
            print "coord1 = ",coord1

        self.readMarkers([-4,1,0,3])
        ints = self.readIntBlock()
        print "*ints = ",ints

        self.readMarkers([-5,1,0,0,2])

    def readTable_Geom2(self):
        word = self.readStringBlock()
        print "word = |%r|" %(word)
        self.readMarkers([-1,7])
        ints = self.readIntBlock()
        print "*ints = ",ints
        self.readMarkers([-2,1,0,2])
        word = self.readStringBlock()
        #print "word = |%r|" %(word)

        self.readMarkers([-3,1,0,17])

        #marker = self.getMarker() # 17
        #print "marker = ",marker
        #self.printTableCode(marker)

        ints = self.readIntBlock()
        print "*ints = ",ints
        self.readMarkers([-4,1,0,3])
        ints = self.readIntBlock()
        print "*ints = ",ints

    def readTable_Geom3(self):
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

    def readTable_Geom4(self):
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
        self.readHollerith()
        
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
        self.readHollerith()

        self.readMarkers([-4,1,0,32])
        self.skip(4*34)
        self.readMarkers([-5,1,0,0,2])

    def readTable_OES1X1(self):
        word = self.readStringBlock() # OES1X1
        print "word = |%r|" %(word)
        self.readMarkers([-1,7])
        print "****",self.op2.tell()
        data = self.readBlock()
        #self.printBlock(data)

        #print "len(block) = ",len(data)
        #(aCode,tCode,elType,iSubcase) = unpack('iiii',data[:16])
        #print "aCode=%s tCode=%s elType=%s iSubcase=%s" %(aCode,tCode,elType,iSubcase)
        #data = data[16:]
        
        #word5,word6,word7 = unpack('iii',data[:12]) # depends on aCode,tCode
        #data = data[12:]
        
        #(loadset,fcode,numWordsEntry,sCode) = unpack('iiii',data[:16])
        #print "loadset=%s fcode=%s numWordsEntry=%s sCode=%s" %(loadset,fcode,numWordsEntry,sCode)
        #data = data[16:]
        
        #self.printBlock(data)
        
        print "****",self.op2.tell()
        assert self.op2.tell()==4880
        self.readMarkers([-2,1,0,7])
        self.printSection(100)
        word = self.readStringBlock()  # OES1
        print "word = |%r|" %(word)
        self.readMarkers([-3,1,0,146])
        #sys.exit('oes')
        self.skip(4*51)
        word = self.readString(384)
        self.readHollerith()
        print "word* = |%s|" %(word)
        self.readMarkers([-4,1,0,87])
        self.skip(16)

        floats = self.readFloats(4*85)
        print "*floats = ",floats,'\n'

        self.readMarkers([-5,1,0,])
        self.readMarkers([0,0,])

if __name__=='__main__':
    op2 = Op2('quad4.op2')
    #op2 = Op2('tria3.op2')
    op2.read()

    print "done..."