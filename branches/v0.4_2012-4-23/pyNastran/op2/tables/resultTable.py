## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
import sys
import copy
from struct import unpack

from pyNastran.op2.op2Errors     import *
from pyNastran.op2.tables.oug.oug  import OUG
from pyNastran.op2.tables.oes_stressStrain.oes import OES

from pyNastran.op2.tables.oqg_constraintForces.oqg   import OQG
from pyNastran.op2.tables.oef_forces.oef import OEF
from pyNastran.op2.tables.opg_appliedLoads.opg import OPG
from pyNastran.op2.tables.oee_energy.oee import OEE
from pyNastran.op2.tables.ogf_gridPointForces.ogf import OGF
#from pyNastran.op2.tables.hisadd import HISADD - combined with R1TAB for now
from pyNastran.op2.tables.r1tab  import R1TAB
from pyNastran.op2.tables.destab import DESTAB
from pyNastran.op2.tables.lama_eigenvalues.lama import LAMA



class ResultTable(OQG,OUG,OEF,OPG,OES,OEE,OGF,R1TAB,DESTAB,LAMA):

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
        self.deleteAttributes_OPG()

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

    def createTransientObject(self,storageObj,classObj,debug=False):
        """
        Creates a transient object (or None if the subcase should be skippied).
        @param storageObj  the dictionary to store the object in (e.g. self.bars)
        @param classObj    the class object to instantiate
        @note dt can also be loadStep depending on the class
        """
        if debug:
            print "create Transient Object"
            print "***NF = ",self.nonlinearFactor
            #print "DC = ",self.dataCode
        
        if self.iSubcase in storageObj:
            #print "updating dt..."
            self.obj = storageObj[self.iSubcase]
            #print "obj = ",self.obj.__class__.__name__
            #print self.obj.writeF06(['',''],'PAGE ',1)[0]
            
            try:
                self.obj.updateDt(self.dataCode,self.nonlinearFactor)
            except:
                #try:
                    #print "objName = ",self.obj.name()
                #except:
                    #print "objName = ",self.obj
                raise
            ###
        else:
            self.obj = classObj(self.dataCode,self.iSubcase,self.nonlinearFactor)
            #print "obj2 = ",self.obj.__class__.__name__
        storageObj[self.iSubcase] = self.obj
        ###

    def readResultsTable(self,table3,table4Data,flag=0):
        self.dtMap = {}
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

        exitFast = False
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

            n = self.op2.tell()
            marker = self.getMarker()
            self.goto(n)
            if marker!=146:
                print "marker = ",marker
                exitFast = True
                break

            table3(iTable)
            self.dataCode['tableName'] = self.tableName
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
            #self.log.debug("")
            #print "i read the markers!!!"
   
        ###
        nOld = self.op2.tell()
        #try:
        if not(exitFast):
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
        del self.dtMap

    def readTable4(self,table4Data,flag,iTable):
        """loops over repeated table -4s"""
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
        """checks to see if table 4 is done, loads the data, and handles skipping"""
        isTable4Done = False
        isBlockDone  = False

        bufferWords = self.getMarker(self.tableName)
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

    def updateDtMap(self):
        """
        Interfaces with setTransientTimes(times) to limit the amount
        of output that is in the transient results.  This allows
        for models with 1000s of time steps that would otherwise
        crash with memory errors to run.  Every result may be extracted
        if the OP2 is read multiple times.
        
        While not ideal, this function prevents having to drastically
        change the code to support large models, which would
        make the OP2 reader not as useful for reasonably sized models.
        
        The code works by taking the user-provided array of values
        for a given subcase and b/c it's an array can be subtracted
        from the current value of dt.  Then the minimum absolute value
        is found and that is mapped back to the original value.  If a
        closer value to the target user-specified value is found, the
        previous result will be deleted, which keeps memory usage down.
        The new result will then be read.  If a dt is not closer than
        an existing value to a given case, that case will not be read.
        """
        #numArray = array([1.,2.])
        iSubcase = self.iSubcase
        numArray = self.expectedTimes[iSubcase]
        #nums = [0.9,1.11,  1.89,2.1]
        num = self.obj.getTransients()[-1] # they're sorted so the last value is the current dt

        readCase = True
        delta = numArray-num
        absDelta = list(abs(delta))
        closest = min(absDelta)
        iclose = absDelta.index(closest)
        actualValue = numArray[iclose]

        if iSubcase not in self.dtMap:
            self.dtMap[iSubcase] = {}
        if iclose in self.dtMap[iSubcase]:
            v1 = self.dtMap[iSubcase][iclose]
            vact = getCloseNum(v1,num,actualValue)

            if vact!=self.dtMap[iSubcase][iclose]:
                del self.dtMap[iSubcase][iclose]
                self.obj.deleteTransient(v1)
                #print "num=%s closest=%s iclose=%s" %(num,actualValue,iclose)
                #print "***deleted v1=%s num=%s vact=%s actual=%s" %(v1,num,vact,actualValue)
                self.dtMap[iSubcase][iclose] = vact
                readCase = True
                #print self.dtMap
                #print "A"
            else: # cleanup previous creation of empty dt case (happened in updateDt outside this function)
                readCase = False
                #print self.dtMap
                #print "num=%s closest=%s iclose=%s" %(num,actualValue,iclose)
                #print "B"
                self.obj.deleteTransient(num)
            ###
        else: # read case
            self.dtMap[iSubcase][iclose] = num
            readCase = True
            #print "num=%s closest=%s iclose=%s" %(num,actualValue,iclose)
            #print self.dtMap
            #print "C"
        ###
        #print "delta = ",delta,'\n'
        #print "readCase = ",readCase
        #if num>=0.14:
            #print self.obj.getTransients()
            #sys.exit('OUG !!!')
            
        return readCase

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

    def readMappedScalarsOut(self,debug=False):
        readCase = True
        #print "isSort1() = ",self.isSort1()
        if self.iSubcase in self.expectedTimes and len(self.expectedTimes[self.iSubcase])>0:
            readCase = self.updateDtMap()
        
        if self.obj and readCase and self.isSort1():
            self.readScalarsOut(debug=False)
        else:
            self.skipOES_Element()
        ###

    def readScalarsOut(self,debug=False):
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
        #print  "strFormat = ",strFormat
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
        self.handleResultsBuffer(self.readScalarsOut,debug=False)

    def readScalarsX(self,strFormat,nTotal,debug=False):
        """
        unused...
        """
        assert debug==True or debug==False
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
        assert debug==True or debug==False
        data = self.data
        deviceCode = self.deviceCode
        #print type(scalarObject)
        
        n = 0
        nEntries = len(data)//32
        if debug:
            self.log.debug('calling readScalars8 debug')
            print "self.obj = ",self.obj
            
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
                #self.log.debug(self.printBlock(self.data[n:n+64]))
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
        assert debug==True or debug==False
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

    def readScalars14(self,debug=False):
        """
        @see readScalars8
        """
        assert debug==True or debug==False
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
                self.log.debug("grid=%-7i dx=%.2g dy=%g dz=%g rx=%g ry=%g rz=%g" %(grid,dx, dy, dz, rx, ry, rz))
                self.log.debug("     %-7s dx=%.2g dy=%g dz=%g rx=%g ry=%g rz=%g" %('',  dxi,dyi,dzi,rxi,ryi,rzi))
            self.obj.add([grid,gridType,dx, dy, dz, rx, ry, rz,
                                        dxi,dyi,dzi,rxi,ryi,rzi])
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
        assert debug==True or debug==False
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

def getCloseNum(v1,v2,closePoint):
    numList = [v1,v2]
    delta = array([v1,v2])-closePoint
    #print "**delta=%s" %(delta)
    absDelta = list(abs(delta))
    closest = min(absDelta)
    iclose = absDelta.index(closest)
    actualValue = numList[iclose]
    return actualValue

    