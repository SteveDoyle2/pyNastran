from fortranFile import FortranFile
from op2Codes import Op2Codes
import os
import sys
import struct
from struct import unpack

from op2_Objects import *
from geometryTables import GeometryTables
from elementsStressStrain import ElementsStressStrain

class Op2(FortranFile,Op2Codes,GeometryTables,ElementsStressStrain):
    def __init__(self,infileName): 
        self.infilename = infileName
        self.tablesToRead = ['GEOM1','GEOM2','GEOM3','GEOM4','OQG1','OUGV1','OES1X1']  # 'OUGV1','GEOM1','GEOM2'
        ## GEOM1 & GEOM2 are skippable on simple problems...hmmm
    
    def readTapeCode(self):
        self.printSection(500)
        sys.exit('a')
        self.readMarkers([3])
        ints = self.readIntBlock()
        #print "*ints = ",ints
        self.readMarkers([7])

        word = self.readStringBlock() # Nastran Fort Tape ID Code - 
        #print "word = |%r|" %(word)

        self.readMarkers([2])
        ints = self.readIntBlock()
        #print "*ints = ",ints

        self.readMarkers([-1])

        #data = self.getData(60)
        #self.printBlock(data)

    def read(self):
        self.op2 = open(self.infilename,'rb')
        
        self.n = self.op2.tell()
        print "self.n = ",self.n
        self.readTapeCode()

        isAnotherTable = True
        while isAnotherTable:
            tableName = self.readTableName(rewind=True)
            print "tableName = |%r|" %(tableName)
 
            if tableName in self.tablesToRead:
                if tableName=='GEOM1': # nodes,coords,etc.
                    self.readTable_Geom1()
                elif tableName=='GEOM2': # elements
                    self.readTable_Geom2()
                elif tableName=='GEOM3': # static/thermal loads
                    self.readTable_Geom3()
                elif tableName=='GEOM4': # constraints
                    self.readTable_Geom3()

                elif tableName=='EPT':   # element properties
                    self.readTable_EPT()
                elif tableName=='MPTS':  # material properties
                    self.readTable_MPTS()

                elif tableName=='OQG1':  # spc forces
                    self.readTable_OQG1()
                elif tableName=='OUGV1': # displacements/velocity/acceleration
                    self.readTable_OUGV1()
                elif tableName=='OES1X1': # stress
                    self.readTable_OES1X1()
                else:
                    raise Exception('unhandled tableName=|%s|' %(tableName))
                print "---isAnotherTable---"
                (isAnotherTable) = self.hasMoreTables()
            else:
                (isAnotherTable) = self.skipNextTable()
                continue
            print "*** finished tableName = |%r|" %(tableName)
            ###
        ###
        #

        #tableName = self.readTableName(rewind=True)
        #self.readTable_Geom2()
        #self.readTable_Geom3()
        #self.readTable_Geom4()
        #self.readTable_OQG1()
        #self.readTable_OES1X1()
        print "---end of all tables---"

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

        self.readMarkers([-5,1,0])
        #print str(spcForcesObj)

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
        #assert self.op2.tell()==4780,self.op2.tell()
        #sys.exit('end of displacements')

    def readScalars(self,deviceCode,data,scalarObject):
        while data:
            (gridDevice,gridType,dx,dy,dz,rx,ry,rz) = unpack('iiffffff',data[0:32])
            #print "gridDevice = ",gridDevice
            #print "deviceCode = ",deviceCode
            grid = (gridDevice-deviceCode)/10
            #print "grid=%g dx=%g dy=%g dz=%g" %(grid,dx,dy,dz)
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
        if elementType==144:
            self.CQUAD4_144(data)  # 144
        elif elementType==74:
            self.CTRIA3_74(data)  # 74
        elif elementType==39:
            self.CTETRA_39(data)  # 39
        else:
            raise RuntimeError('elementType=%s -> %s is not supported' %(elementType,self.ElementType(elementType)))

        self.readMarkers([-5,1,0,])
        #print "tell5 = ",self.op2.tell()
        self.readMarkers([0,0,])
        #print "end tell = ",self.op2.tell()

    def parseAnalysisCode(self,data):
        #self.printBlock(data)
        (aCode,tCode,elementType,iSubcase) = unpack('iiii',data[:16])
        deviceCode   = aCode%10
        approachCode = (aCode-deviceCode)/10
        print "aCode=%s analysisCode=%s deviceCode=%s tCode=%s elementType=%s iSubcase=%s" %(aCode,approachCode,deviceCode,tCode,elementType,iSubcase)
        return (approachCode,deviceCode,tCode,elementType,iSubcase)
