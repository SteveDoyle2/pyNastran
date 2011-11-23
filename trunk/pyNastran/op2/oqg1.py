import sys
from struct import unpack
from pyNastran.op2.resultObjects.op2_Objects import spcForcesObject

class OQG1(object):
    def readTable_OQG1(self):
        ## OQG1
        word = self.readTableName(rewind=False) # OQG1
        self.tableInit(word)
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
        self.data = self.readBlock()
        self.printBlock(self.data)

        self.spcForces[self.iSubcase] = spcForcesObject(self.iSubcase)
        self.readScalars(self.spcForces[self.iSubcase])
        
    def readScalars(self,scalarObject):
        data = self.data
        deviceCode = self.deviceCode
        #print type(scalarObject)
        while data:
            #print "self.numWide = ",self.numWide
            #print "len(data) = ",len(data)
            #self.printBlock(data[32:])
            out = unpack('iiffffff',data[0:32])
            (gridDevice,gridType,dx,dy,dz,rx,ry,rz) = out
            self.op2Debug.write('%s\n' %(str(out)))
            #print "gridDevice = ",gridDevice
            #print "deviceCode = ",deviceCode
            grid = (gridDevice-deviceCode)/10
            #print "grid=%g dx=%g dy=%g dz=%g rx=%g ry=%g rz=%g" %(grid,dx,dy,dz,rx,ry,rz)
            scalarObject.add(grid,gridType,dx,dy,dz,rx,ry,rz)
            data = data[32:]
        ###
        self.data = data
        self.handleScalarBuffer(self.readScalars,scalarObject)

    def handleScalarBuffer(self,func,stress):
        """
        works by knowing that:
        the end of an unbuffered table has a
            - [4]
        the end of an table with a buffer has a
            - [4,4,x,4] where x is the next buffer size, which may have another buffer
        the end of the final buffer block has
            - nothing!
        
        The code knows that the large buffer is the default size and the
        only way there will be a smaller buffer is if there are no more
        buffers.  So, the op2 is shifted by 1 word (4 bytes) to account for
        this end shift.  An extra marker value is read, but no big deal.
        Beyond that it's just appending some data to the binary string
        and calling the function that's passed in
        """
        #print stress
        print "len(data) = ",len(self.data)
        #if marker[0]==4:
        #    print "found a 4 - end of unbuffered table"

        if len(self.data)>0:
            print "*******len(self.data)=%s...assuming a buffer block" %(len(self.data))
            markers = self.readHeader()
            data = self.readBlock()
            #if len(data)<marker:
            #    self.goto(self.n-4) # handles last buffer not having an extra 4
            self.data += data
            func(stress)
        ###
