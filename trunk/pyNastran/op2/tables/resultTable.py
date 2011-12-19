import sys
import copy
from struct import unpack

from pyNastran.op2.op2Errors     import *
from pyNastran.op2.tables.oug.ougv1  import OUGV1
from pyNastran.op2.tables.oes.oes    import OES
from pyNastran.op2.tables.oqg.oqg1   import OQG1
from pyNastran.op2.tables.oef.oef    import OEF
from pyNastran.op2.tables.ogp.ogp    import OGP
from pyNastran.op2.tables.oee.oee    import OEE
#from pyNastran.op2.tables.hisadd import HISADD - combined with R1TAB for now
from pyNastran.op2.tables.r1tab  import R1TAB
from pyNastran.op2.tables.destab import DESTAB


class ResultTable(OQG1,OUGV1,OEF,OGP,OES,OEE,R1TAB,DESTAB):

    def readTableA_DUMMY(self):
        """reads a dummy geometry table"""
        self.iTableMap = {
                         }
        self.readRecordTable('DUMMY')

    def readTableB_DUMMY(self):
        """reads a dummy results table"""
        self.tableName = 'DUMMY'
        table3     = self.readTable_DUMMY_3
        table4Data = self.readDUMMY_Data
        self.readResultsTable(table3,table4Data)
        self.deleteAttributes_OGP()

    def readTable_DUMMY_3(self,iTable):
        """sets dummy parameters"""
        self.analysisCode = None
        self.tableCode    = None
        self.formatCode   = None
        self.sortCode     = None

    def readDUMMY_Data(self):
        """creates a dummy object and skips the results"""
        self.obj = None
        self.readOES_Element()

    def createTransientObject(self,storageObj,classObj):
        """
        Creates a transient object (or None if the subcase should be skippied).
        @param storageObj  the dictionary to store the object in (e.g. self.bars)
        @param classObj    the class object to instantiate
        @note dt can also be loadStep depending on the class
        """
        #print "create Transient Object"
        #print "***NF = ",self.nonlinearFactor
        #print "DC = ",self.dataCode
        
        if self.iSubcase in storageObj:
            #print "updating dt..."
            self.obj = storageObj[self.iSubcase]
            
            try:
                self.obj.updateDt(self.dataCode,self.nonlinearFactor)
            except:
                #try:
                #    print "objName = ",self.obj.name()
                #except:
                #    print "objName = ",self.obj
                raise
            ###
        else:
            self.obj = classObj(self.dataCode,self.iSubcase,self.nonlinearFactor)
        storageObj[self.iSubcase] = self.obj
        ###

    def readResultsTable(self,table3,table4Data,flag=0):
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
            self.isBufferDone = False
            #print self.printSection(140)
            #print "reading iTable3=%s" %(iTable)
            #self.obj = None

            ## the results object
            self.obj = None
            ## dt/loadFactor/frequency/loadStep value (or None for static)
            self.nonlinearFactor = None
            self.dataCode = {}
            table3(iTable)
            ## developer parameter - Analysis/Table/Format/Sort Codes
            self.atfsCode = [self.analysisCode,self.tableCode,self.formatCode,self.sortCode]
            #print "self.tellA = ",self.op2.tell()
            
            ## ???
            self.isMarker = False

            isBlockDone = self.readTable4(table4Data,flag,iTable-1)
            #self.firstPass = False

            #print "self.tellB = ",self.op2.tell()
            iTable -= 2
            #print "isBlockDone = ",isBlockDone
            #sys.exit('stopping')
            if isBlockDone:
                #print "iTable = ",iTable
                #self.n = self.markerStart
                #self.op2.seek(self.n)
                break
            ###
            n = self.n
            #print self.printSection(100)
            self.readMarkers([iTable,1,0],tableName)
            self.log.debug("")
            #print "i read the markers!!!"
   
        ###
        nOld = self.op2.tell()
        #try:
        self.readMarkers([iTable,1,0],tableName)
        #except InvalidMarkersError:
        #    self.goto(nOld)
            #print self.printBlock(self.data)
            #print self.printSection(100)
            #markerZero = self.getMarker()
            #assert markerZero==0
            #self.goto(nOld)
            #print "finished markerZero"
            #return
            
        
        #print str(self.obj)
        if self.makeOp2Debug:
            self.op2Debug.write("***end of %s table***\n" %(tableName))

    def readTable4(self,table4Data,flag,iTable):
        #self.readMarkers([iTable,1,0])
        markerA = 4
        
        while markerA is not None:
            self.markerStart = copy.deepcopy(self.n)
            #self.printSection(180)
            self.readMarkers([iTable,1,0])
            #print "starting OEF table 4..."
            if flag:
                isTable4Done,isBlockDone = table4Data(iTable)
            else:
                isTable4Done,isBlockDone = self.readTable4DataSetup(table4Data,iTable)
            if isTable4Done:
                #print "done with OEF4"
                self.n = self.markerStart
                self.op2.seek(self.n)
                break
            #print "finished reading oef table..."
            markerA = self.getMarker('A')
            self.n-=12
            self.op2.seek(self.n)
            
            self.n = self.op2.tell()
            #print "***markerA = ",markerA
            
            iTable-=1
            #print "isBlockDone = ",isBlockDone
        ###    
        #print "isBlockDone = ",isBlockDone
        return isBlockDone

    def readTable4DataSetup(self,table4Data,iTable): # iTable=-4
        isTable4Done = False
        isBlockDone  = False

        bufferWords = self.getMarker('OUGV1') ## @todo replace with table name
        #print "bufferWords = ",bufferWords
        #print len(bufferWords)
        self.data = self.readBlock()
        #self.printBlock(data)

        if bufferWords==146:  # table -4 is done, restarting table -3
            isTable4Done = True
            return isTable4Done,isBlockDone
        elif bufferWords==0:
            #print "bufferWords 0 - done with Table4"
            isTable4Done = True
            isBlockDone = True
            return isTable4Done,isBlockDone

        isBlockDone = not(bufferWords)

        if self.isValidSubcase(): # lets the user skip a certain subcase
            table4Data()
        else:
            self.log.debug("***skipping table=%s iSubcase=%s" %(self.tableName,self.iSubcase))
            self.skipOES_Element()
        ###
        return (isTable4Done,isBlockDone)

    def handleResultsBuffer(self,func,debug=False):
        """
        works by knowing that:
        the end of an unbuffered table has a
            - [4]
        the end of an table with a buffer has a
            - [4,4,x,4] where x is the next buffer size, which may have another buffer
        the end of the final buffer block has
            - nothing!

        @param self the object pointer
        @param func the function to recursively call (the function that called this)
        @param debug developer debug

        @note
            The code knows that the large buffer is the default size and the
            only way there will be a smaller buffer is if there are no more
            buffers.  So, the op2 is shifted by 1 word (4 bytes) to account for
            this end shift.  An extra marker value is read, but no big deal.
            Beyond that it's just appending some data to the binary string
            and calling the function that's passed in

        @warning
            Dont modify this without LOTS of testing.
            It's a VERY important function
        """
        #print self.obj
        #print "len(data) = ",len(self.data)
        #if marker[0]==4:
        #    self.log.debug("found a 4 - end of unbuffered table")

        #if debug:
        #    self.log.debug(self.printSection(120))
        
        nOld = self.n
        #try:
        markers = self.readHeader()
        #except AssertionError:  # end of table - poor catch
        #    self.goto(nOld)
        #   self.log.debug(self.printSection(120))
        #    return

        #print "markers = ",markers
        #print self.printSection(160)
        
        if markers<0:  # not a buffer, the table may be done
            self.goto(nOld)
            if markers is not None and markers%2==1:
                self.isBufferDone = True

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
            func()
        ###

    def readScalars4(self,debug=False):
        """
        reads 4 values "ifff" and puts them into the result object
        """
        data = self.data
        deviceCode = self.deviceCode
        #print type(scalarObject)

        n = 0
        nEntries = len(data)//16
        for i in range(nEntries):
            eData = data[n:n+16]
            #print "self.numWide = ",self.numWide
            #print "len(data) = ",len(data)
            #self.printBlock(data[16:])
            #msg = 'len(data)=%s\n'%(len(data))
            #assert len(data)>=16,msg+self.printSection(120)
            out  = unpack('ifff',eData)
            a,b,c,d,E,F,G = unpack('ssssfff',eData)
            #print "abcd=|%s|" %(a+b+c+d)
            (gridDevice,dx,dy,dz) = out
            if self.makeOp2Debug:
                self.op2Debug.write('%s\n' %(str(out)))
            #print "gridDevice = ",gridDevice
            #print "deviceCode = ",deviceCode
            grid = (gridDevice-deviceCode) // 10
            #if grid<100:
            if debug:
                self.log.debug("grid=%-3i dx=%g dy=%g dz=%g" %(grid,dx,dy,dz))
            #print "grid=%g dx=%g dy=%g dz=%g" %(grid,dx,dy,dz)
            self.obj.add(grid,dx,dy,dz)
            n+=16
        ###
        self.data = data[n:]
        #print self.printSection(200)
        self.handleResultsBuffer(self.readScalars4,debug=False)

    def readScalars4o(self,debug=False):
        """
        reads len(strFormat) values and puts it into the result object
        the "o" in readScalars4o means "out" b/c it creates an out tuple
        instead of 4 values like readScalars4
        @note
            nTotal is the number of bytes
            strFormat = 'iiii'
        """
        data = self.data
        #print type(self.obj)
        (nTotal,strFormat) = self.obj.getLength()
        n = 0
        nEntries = len(data)//nTotal
        for i in range(nEntries):
            eData = data[n:n+nTotal]
            out  = unpack(strFormat,eData)
            if debug:
                self.log.debug("*out = %s" %(out))
            self.obj.add(out)
            n+=nTotal
        ###
        self.data = data[n:]
        #print self.printSection(200)
        self.handleResultsBuffer(self.readScalars4o,debug=False)

    def readScalarsX(self,strFormat,nTotal,debug=False):
        """
        unused...
        """
        data = self.data
        deviceCode = self.deviceCode
        #print type(self.boj)
        
        n = 0
        nEntries = len(data)//nTotal
        for i in range(nEntries):
            eData = data[n:n+nTotal]
            #print self.printBlock(data[n:n+nTotal])
            out = unpack(strFormat,eData)
            #print "Xout = ",out
            self.obj.add(out)
            n+=nTotal
        ###
        self.data = data[n:]
        self.handleResultsBuffer(self.readScalarsX,strFormat,nTotal,debug=False)

    #def readScalars8(self,debug=False):
    #    self.readScalarsX(self,'iiffffff',32,debug)

    def readScalars8(self,debug=False):
        """
        @see readScalars4
        """
        data = self.data
        deviceCode = self.deviceCode
        #print type(scalarObject)
        
        n = 0
        nEntries = len(data)//32
        for i in range(nEntries):
            #if debug:
            #    self.log.debug(self.printBlock(self.data[n:n+64]))
            #print self.printBlock(self.data[n:n+64])
            eData = data[n:n+32]
            #print "self.numWide = ",self.numWide
            #print "len(data) = ",len(data)
            #print self.printBlock(data[n:n+60])
            out = unpack('iiffffff',eData)
            (gridDevice,gridType,dx,dy,dz,rx,ry,rz) = out
            #if self.makeOp2Debug:
                #self.op2Debug.write('%s\n' %(str(out)))
                
            #print "gridDevice = ",gridDevice
            #print "deviceCode = ",deviceCode
            grid = (gridDevice-deviceCode) // 10
            #if grid<100:
            #print "grid=%-3s type=%s dx=%g dy=%g dz=%g rx=%g ry=%g rz=%g" %(grid,gridType,dx,dy,dz,rx,ry,rz)
            if debug:
                self.log.debug("grid=%-3i type=%s dx=%g dy=%g dz=%g rx=%g ry=%g rz=%g" %(grid,gridType,dx,dy,dz,rx,ry,rz))
                self.log.debug(self.printBlock(self.data[n:n+64]))
                sys.stdout.flush()
            self.obj.add(grid,gridType,dx,dy,dz,rx,ry,rz)
            n+=32
        ###
        self.data = data[n:]
        #print self.printSection(200)
        self.handleResultsBuffer(self.readScalars8,debug=False)

    #def readScalarsF8(self,debug=False):
    #    self.readScalars(self,'fiffffff',32,debug)

    def readScalarsF8(self,debug=False):
        """
        @see readScalars8
        F is Float
        """
        data = self.data
        deviceCode = self.deviceCode
        #print type(scalarObject)
        
        n = 0
        nEntries = len(data)//32
        #print "len(data) = ",len(data)
        for i in range(nEntries):
            #print self.printBlock(self.data[n:n+64])
            eData = data[n:n+32]
            #print "self.numWide = ",self.numWide
            #print "len(data) = ",len(data)
            self.printBlock(data[n:n+60])
            out = unpack('fiffffff',eData)
            (freq,gridType,dx,dy,dz,rx,ry,rz) = out
            if self.makeOp2Debug:
                self.op2Debug.write('%s\n' %(str(out)))
            #print "gridDevice = ",gridDevice
            #print "deviceCode = ",deviceCode
            #if grid<100:
            if debug:
                self.log.debug("freq=%-3s dx=%g dy=%g dz=%g rx=%g ry=%g rz=%g" %(freq,dx,dy,dz,rx,ry,rz))
            self.obj.add(freq,gridType,dx,dy,dz,rx,ry,rz)
            n+=32
        ###
        self.data = data[n:]
        #print self.printSection(200)
        self.handleResultsBuffer(self.readScalarsF8,debug=False)

    #def readScalars14(self,debug=False):
    #    self.readScalarsX(self,'iiffffffffffff',56,debug)

    def readScalars14(self,debug=True):
        """
        @see readScalars8
        """
        data = self.data
        deviceCode = self.deviceCode
        #print type(scalarObject)
        #print "objName = ",self.obj.name()
        n = 0
        nEntries = len(data)//56
        #print "len(data) = ",len(data)
        #print "nEntries = %s" %(nEntries)
        for i in range(nEntries):
            eData = data[n:n+56]
            #print self.printBlock(self.data[n:n+64])
            #print "self.numWide = ",self.numWide
            #print "len(data) = ",len(data)
            #self.printBlock(data[56:])
            #msg = 'len(data)=%s\n'%(len(data))
            #assert len(data)>=56,msg+self.printSection(120)
            out = unpack('iiffffffffffff',eData)
            (gridDevice,gridType,dx, dy, dz, rx, ry, rz,
                                 dxi,dyi,dzi,rxi,ryi,rzi) = out
            if self.makeOp2Debug:
                self.op2Debug.write('%s\n' %(str(out)))
            #print "gridDevice = ",gridDevice
            #print "deviceCode = ",deviceCode
            grid = (gridDevice-deviceCode) // 10
            #if grid<100:
            if debug:
                self.log.debug("grid=%-7i dx=%.2g dy=%g dz=%g rx=%g ry=%g rz=%g" %(grid,dx,dy,dz,rx,ry,rz))
            self.obj.add(grid,gridType,dx, dy, dz, rx, ry, rz,
                                           dxi,dyi,dzi,rxi,ryi,rzi)
            n+=56
        ###
        self.data = data[n:]
        #print self.printSection(200)
        self.handleResultsBuffer(self.readScalars14,debug=False)

    #def readScalarsF14(self,debug=False):
    #    self.readScalarsX(self,'fiffffffffffff',56,debug)

    def readScalarsF14(self,debug=False):
        """
        @see readScalars8
        F is Float
        """
        data = self.data
        deviceCode = self.deviceCode
        #print type(scalarObject)

        n = 0
        nEntries = len(data)//56
        for i in range(nEntries):
            eData = data[n:n+56]
            #print self.printBlock(self.data[n:n+64])
            #print "self.numWide = ",self.numWide
            #print "len(data) = ",len(data)
            #self.printBlock(data[56:])
            #msg = 'len(data)=%s\n'%(len(data))
            #assert len(data)>=56,msg+self.printSection(120)
            out = unpack('fiffffffffffff',eData)
            (freq,gridType,dx, dy, dz, rx, ry, rz,
                           dxi,dyi,dzi,rxi,ryi,rzi) = out
            if self.makeOp2Debug:
                self.op2Debug.write('%s\n' %(str(out)))
            #print "gridDevice = ",gridDevice
            #print "deviceCode = ",deviceCode
            #if grid<100:
            if debug:
                self.log.debug("gridType=%s freq=%-7i dx=%.2g dy=%g dz=%g rx=%g ry=%g rz=%g" %(gridType,freq,dx,dy,dz,rx,ry,rz))
            self.obj.add(freq,gridType,dx, dy, dz, rx, ry, rz,
                                       dxi,dyi,dzi,rxi,ryi,rzi)
            n+=56
        ###
        self.data = data[n:]
        #print self.printSection(200)
        self.handleResultsBuffer(self.readScalars14,debug=False)
    