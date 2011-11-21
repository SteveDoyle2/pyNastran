import sys
from struct import unpack
from op2_Objects import spcForcesObject

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

        iTable = -3
        self.readTable_OQG1_3_4(iTable)
        iTable -= 2

        self.readMarkers([-5,1,0])
        #print str(self.spcForces[self.iSubcase])

        #self.printSection(80)
        #sys.exit('end of oqg1')
    
    
    def readTable_OQG1_3_4(self,iTable):
        self.readMarkers([iTable,1,0])  # iTable=-3
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4

        data = self.getData(4)
        bufferSize, = unpack('i',data)
        print "bufferSize = ",bufferSize
        data = self.getData(4*50)
        aCode = self.getBlockIntEntry(data,1)
        print "aCode = ",aCode
        (three) = self.parseApproachCode(data)


        word = self.readString(384)
        print "word = |%s|" %(word)
        self.readHollerith()
        
        self.readMarkers([iTable-1,1,0])  # iTable=-4
        wordCount = self.getMarker()
        data = self.readBlock()
        #self.printBlock(data)

        self.spcForces[self.iSubcase] = spcForcesObject(self.iSubcase)
        self.readScalars(data,self.spcForces[self.iSubcase])
        
    def readScalars(self,data,scalarObject):
        while data:
            #print "len(data) = ",len(data)
            #self.printBlock(data[32:])
            (gridDevice,gridType,dx,dy,dz,rx,ry,rz) = unpack('iiffffff',data[0:32])
            #print "gridDevice = ",gridDevice
            #print "deviceCode = ",deviceCode
            grid = (gridDevice-self.deviceCode)/10
            #print "grid=%g dx=%g dy=%g dz=%g rx=%g ry=%g rz=%g" %(grid,dx,dy,dz,rx,ry,rz)
            scalarObject.add(grid,dx,dy,dz,rx,ry,rz)
            data = data[32:]
        ###
        

