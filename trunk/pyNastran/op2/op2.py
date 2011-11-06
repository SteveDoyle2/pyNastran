from fortranFile import FortranFile
from op2Codes import Op2Codes
import os
import sys
import struct
from struct import unpack
from op2_Objects import *

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
        assert self.op2.tell()==584,self.op2.tell()

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
        assert self.op2.tell()==976,self.op2.tell()

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


class Op2(FortranFile,Op2Codes,GeometryTables):
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
        self.readMarkers([-1])

        tableName = self.readTableName(rewind=True)
        print "tableName = |%r|" %(tableName)
        (isAnotherTable) = self.skipNextTable()

        #self.readTable_Geom1()
        assert self.op2.tell()==584,self.op2.tell()
        tableName = self.readTableName(rewind=True)

        self.readTable_Geom2()
        self.readTable_Geom3()
        self.readTable_Geom4()
        self.readTable_OQG1()
        self.readTable_OES1X1()
        print "end of oes1x1"

        #self.printSection(4*51+12)
        

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

        self.readMarkers([-5,1,0,0,2])
        print str(spcForcesObj)

        word = self.readStringBlock()  # OUGV1
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

        data = self.getData(4*51)
        #self.printBlock(data)

        #self.skip(4*51)
        word = self.readString(384)
        print "word = |%s|" %(word)
        self.readHollerith()

        self.readMarkers([-4,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()
        
        iSubcase = 1 ## @todo temporary
        dispObj = displacementObject(iSubcase)
        self.readScalars(deviceCode,data,dispObj)
        print str(dispObj)

        self.readMarkers([-5,1,0])
        assert self.op2.tell()==4780,self.op2.tell()
        #sys.exit('end of displacements')

    def readScalars(self,deviceCode,data,scalarObject):
        while data:
            (gridDevice,gridType,dx,dy,dz,rx,ry,rz) = unpack('iiffffff',data[0:32])
            #print "gridDevice = ",gridDevice
            #print "deviceCode = ",deviceCode
            grid = (gridDevice-deviceCode)/10
            print "grid=%g dx=%g dy=%g dz=%g" %(grid,dx,dy,dz)
            scalarObject.add(grid,dx,dy,dz,rx,ry,rz)
            data = data[32:]
        ###

    def readTable_OES1X1(self):
        word = self.readTableName(rewind=False) # OES1X1
        print "word = |%r|" %(word)

        self.readMarkers([-1,7])
        print "****",self.op2.tell()
        data = self.readBlock()
        #self.printBlock(data)
        print "****",self.op2.tell()
        assert self.op2.tell()==4880

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
        assert self.op2.tell()==5656

        data = self.getData(bufferWords*4)
        #self.printBlock(data)
        if elementType==144:
            self.CQUAD4(data)  # 144
        else:
            raise RuntimeError('elementType=%s is not supported' %(elmentType))

        self.readMarkers([-5,1,0,])
        #print "tell5 = ",self.op2.tell()
        self.readMarkers([0,0,])
        #print "end tell = ",self.op2.tell()

    def CQUAD4(self,data):
        """
        GRID-ID  DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR MINOR,VONMISES
        """
        nNodes = 5 # centroid + 4 corner points
        #self.printSection(20)
        #term = data[0:4] CEN/
        #data = data[4:]
        print "*****"
        while data:
            for nodeID in range(nNodes):   #nodes pts
                if nodeID==0:
                    (eid,_,_,_,_) = struct.unpack("issss",data[0:8])
                    data = data[8:]
                    eid = (eid - 1) / 10

                eData = data[0:4*17]
                data  = data[4*17: ]
                out = unpack('iffffffffffffffff',eData[0:68])
                (grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1,
                      fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2,) = out
                #print "eid=%i grid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" %(eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
                #print "               fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"          %(fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
                #print "len(data) = ",len(data)
            ###
            #sys.exit('asdf')
        self.skip(4)
        ###

    def parseAnalysisCode(self,data):
        #self.printBlock(data)
        (aCode,tCode,elementType,iSubcase) = unpack('iiii',data[:16])
        deviceCode   = aCode%10
        approachCode = (aCode-deviceCode)/10
        print "aCode=%s analysisCode=%s deviceCode=%s tCode=%s elementType=%s iSubcase=%s" %(aCode,approachCode,deviceCode,tCode,elementType,iSubcase)
        return (approachCode,deviceCode,tCode,elementType,iSubcase)
