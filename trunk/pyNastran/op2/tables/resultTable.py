#import sys
#from struct import unpack

from pyNastran.op2.tables.ougv1 import OUGV1
from pyNastran.op2.tables.oqg1  import OQG1
from pyNastran.op2.tables.oes   import OES
from pyNastran.op2.tables.oef   import OEF
from pyNastran.op2.tables.ogp   import OGP

class ResultTable(OQG1,OUGV1,OEF,OGP,OES):

    def readResultsTable(self,table3,table4Data):
        tableName = self.readTableName(rewind=False) # OEF
        self.tableInit(tableName)
        #print "tableName = |%r|" %(tableName)

        self.readMarkers([-1,7],tableName)
        ints = self.readIntBlock()
        #print "*ints = ",ints

        self.readMarkers([-2,1,0],tableName) # 7
        bufferWords = self.getMarker()
        #print "1-bufferWords = ",bufferWords,bufferWords*4
        ints = self.readIntBlock()
        #print "*ints = ",ints
        
        markerA = -4
        markerB = 0

        iTable=-3
        self.readMarkers([iTable,1,0],tableName)
        while [markerA,markerB]!=[0,2]:
            table3(iTable)
            isBlockDone = self.readTable4(table4Data,iTable-1)
            iTable -= 2

            if isBlockDone:
                #print "iTable = ",iTable
                #self.n = self.markerStart
                #self.op2.seek(self.n)
                break
            ###
            n = self.n
            self.readMarkers([iTable,1,0],tableName)
            #print "i read the markers!!!"
   
        ###
        self.readMarkers([iTable,1,0],tableName)
        #print str(self.obj)
        if self.makeOp2Debug:
            self.op2Debug.write("***end of %s table***\n" %(tableName))

    def handleResultsBuffer(self,func,stress,debug=False):
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
        #print "len(data) = ",len(self.data)
        #if marker[0]==4:
        #    print "found a 4 - end of unbuffered table"

        if debug:
            print self.printSection(120)

        nOld = self.n
        markers = self.readHeader()
        #print "markers = ",markers
        if markers<0:
            self.goto(nOld)
            #print self.printSection(120)
            #sys.exit('found a marker')
            #print 'found a marker'

        else:
            #print "*******len(self.data)=%s...assuming a buffer block" %(len(self.data))
            #markers = self.readHeader()
            #print "markers = ",markers
            data = self.readBlock()
            #if len(data)<marker:
            #    self.goto(self.n-4) # handles last buffer not having an extra 4
            self.data += data
            func(stress)
        ###

    def readScalars8(self,scalarObject):
        data = self.data
        deviceCode = self.deviceCode
        #print type(scalarObject)
        while len(data)>=32:
            #print "self.numWide = ",self.numWide
            #print "len(data) = ",len(data)
            #self.printBlock(data[32:])
            msg = 'len(data)=%s\n'%(len(data))
            assert len(data)>=32,msg+self.printSection(120)
            out = unpack('iiffffff',data[0:32])
            (gridDevice,gridType,dx,dy,dz,rx,ry,rz) = out
            if self.makeOp2Debug:
                self.op2Debug.write('%s\n' %(str(out)))
            #print "gridDevice = ",gridDevice
            #print "deviceCode = ",deviceCode
            grid = (gridDevice-deviceCode)/10
            #if grid<100:
            #    print "grid=%-3i dx=%g dy=%g dz=%g rx=%g ry=%g rz=%g" %(grid,dx,dy,dz,rx,ry,rz)
            scalarObject.add(grid,gridType,dx,dy,dz,rx,ry,rz)
            data = data[32:]
        ###
        self.data = data
        #print self.printSection(200)
        self.handleResultsBuffer(self.readScalars8,scalarObject,debug=False)

    def readScalars14(self,scalarObject):
        data = self.data
        deviceCode = self.deviceCode
        #print type(scalarObject)
        while len(data)>=56:
            #print "self.numWide = ",self.numWide
            #print "len(data) = ",len(data)
            #self.printBlock(data[56:])
            msg = 'len(data)=%s\n'%(len(data))
            assert len(data)>=56,msg+self.printSection(120)
            out = unpack('iiffffffffffff',data[0:56])
            (gridDevice,gridType,dx, dy, dz, rx, ry, rz,
                                 dxi,dyi,dzi,rxi,ryi,rzi) = out
            if self.makeOp2Debug:
                self.op2Debug.write('%s\n' %(str(out)))
            #print "gridDevice = ",gridDevice
            #print "deviceCode = ",deviceCode
            grid = (gridDevice-deviceCode)/10
            #if grid<100:
            #   print "grid=%-7i dx=%.2g dy=%g dz=%g rx=%g ry=%g rz=%g" %(grid,dx,dy,dz,rx,ry,rz)
            #scalarObject.add(grid,gridType,dx, dy, dz, rx, ry, rz,
            #                               dxi,dyi,dzi,rxi,ryi,rzi)
            data = data[56:]
        ###
        self.data = data
        #print self.printSection(200)
        self.handleResultsBuffer(self.readScalars14,scalarObject,debug=False)
    