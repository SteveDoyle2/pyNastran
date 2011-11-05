#from fortranFile import FortranFile
import struct
from struct import unpack

class FortranFile(object):
    def __init__(self):
        self.endian = '>'
    
    def readHeader(self):
        data = self.op2.read(12)
        ints = self.getInts(data)
        print "header ints = ",ints
        self.n += 12
        return ints[1]

    def getStrings(self,data):
        n = len(data)
        sFormat = 's'*n
        strings = unpack(sFormat,data)
        return strings

    def readInts(self,nInts):
        nData = 4*nInts
        #print "nData = ",nData
        data = self.op2.read(nData)

        iFormat = 'i'*nInts
        ints = unpack(iFormat,data)
        self.n+=nData
        return ints

    def getInts(self,data):
        n = len(data)
        nInts = (n/4)
        iFormat = 'i'*nInts
        ints = unpack(iFormat,data[:nInts*4])
        return ints

    def getFloats(self,data):
        n = len(data)
        nFloats = (n/4)
        iFormat = 'f'*nFloats
        ints = unpack(iFormat,data[:nFloats*4])
        return ints

    def printSection(self,nBytes):
        """
        prints data, but doesnt move the cursor
        """
        data = self.op2.read(nBytes)
        ints    = self.getInts(data)
        floats  = self.getFloats(data)
        strings = self.getStrings(data)
        print "ints    = ",ints
        print "floats  = ",floats
        print "strings = |%r|" %(''.join(strings))
        self.op2.seek(self.n)
    
    def readString(self,nData):
        data = self.op2.read(nData)
        string = ''.join(self.getStrings(data))
        self.n += nData
        return string

    def skip(self,n):
        #print "\n--SKIP--"
        #print "tell = ",self.op2.tell()
        #print "n = ",n
        #print "self.n = ",self.n
        self.n += n
        self.op2.seek(self.n)
        #print "*tell = ",self.op2.tell()
        #print "\n"
        
class Op2(FortranFile):
    def __init__(self): 
        self.infilename = 'quad4.op2'
    
    def read(self):
        self.op2 = open(self.infilename,'rb')
 
        self.n = self.op2.tell()
        print "self.n = ",self.n
        self.readMarkers([3])
        self.skip(36) # geometry header???
        
        #self.printSection(28) # nastran fort tape id code
        self.skip(28) # name
        self.skip(4)

        #word = self.readStringBlock()
        #print "word = |%r|" %(word)

        self.readMarkers([2])

        self.skip(16)
        
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
        self.skip(140)
        self.skip(136)
        
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
        self.skip(156)
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
        self.skip(104)

        self.readMarkers([-4,1,0,9])
        self.skip(44)

        self.readMarkers([-5,1,0,3])
        self.skip(20)

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
        self.skip(64)
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
        
        self.skip(204)
        data = self.op2.read(384)
        words = self.getStrings(data)
        print "word = |%s|" %(''.join(words))
        #self.skip(0)
        data = self.op2.read(4)  # weird hollerith
        self.n+=4
        
        self.readMarkers([-4,1,0,32])
        self.scan(136)

        self.readMarkers([-5,1,0,0,2])
        #word = self.readStringBlock()
        #print "word = |%r|" %(word)
        self.scan(16)
        self.readMarkers([-1,7])

        self.printSection(200) # 
        
        


    def scan(self,n):
        data = self.op2.read(n)
        self.n+=n

    def startTable(self,markers):
        self.readMarkers(markers)
        word = self.readStringBlock()
        return word

    def getTableCode(self):
        tableCode = self.readHeader()
        return tableCode
        
    def readMarkers(self,markers):
        for marker in markers:
            tableCode = self.readHeader()
            assert marker==tableCode
        ###

    def readStringBlock(self):
        data = self.op2.read(4)
        nValues, = unpack('i',data)
        
        self.n+=4
        word = self.readString(nValues)
        print "word = |%s|" %(word)
        self.skip(4)
        return word

    def readIntBlock(self):
        data = self.op2.read(4)
        self.n+=4

        nValues, = unpack('i',data)
        #print "nValues = ",nValues/4
        ints = self.readInts(nValues/4)
        self.skip(4)
        #print "ints = ",ints
        return ints

if __name__=='__main__':
    op2 = Op2()
    op2.read()
