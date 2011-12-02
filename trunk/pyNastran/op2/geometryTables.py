import os
import sys
import struct
from struct import unpack

from pyNastran.op2.op2Errors import *
from pyNastran.op2.tables.geom1 import Geometry1
from pyNastran.op2.tables.geom2 import Geometry2
from pyNastran.op2.tables.geom3 import Geometry3
from pyNastran.op2.tables.geom4 import Geometry4
from pyNastran.op2.tables.ept import EPT
from pyNastran.op2.tables.mpt import MPT
#from pyNastran.op2.tables.dynamics import DYNAMICS

class GeomObj(object):
    def __init__(self):
        pass
    def geomFunc(self,data):
        pass

class GeometryTables(Geometry1,Geometry2,Geometry3,Geometry4,EPT,MPT):

    def readRecordTable(self,expectedTableName):
        """
        @note assumes self.iTableMap has already been set
        """
        tableName = self.readTableName(rewind=False) # GEOM1
        self.tableInit(tableName)
        #print "*tableName = |%r|" %(tableName)

        self.readMarkers([-1,7])
        fields = self.readIntBlock()
        #print "fields = ",fields

        self.readMarkers([-2,1,0]) # 2
        bufferWords = self.getMarker()
        #print "bufferWords = ",bufferWords,bufferWords*4
        word = self.readStringBlock()

        iTable = -3
        while 1:  ## @todo could this cause an infinite loop...i dont this so...
            (tableName,isNextTable,isNextSubTable,isFileDone) = self.readGeomSubTable(iTable)
        
            if self.checkForNextTable() or isFileDone:
                #sys.exit('end of geom1')
                return
            iTable -= 1
        ###
        sys.exit('end of %s-this should never happen...' %(expectedTableName))

    def checkForNextTable(self):
        foundTable = False
        #print "---checking---"
        word = self.readTableName(rewind=True,debug=False,stopOnFailure=False)
        if word != None:
            foundTable = True
        #print '---checked---'
        #print "geomWord = ",word
        return foundTable

    def checkFileDone(self,n):
        isFileDone = False
        #print "tell = ",self.op2.tell()
        
        nOld = self.op2.tell()
        try:
            #print self.printSection(60)
            self.readMarkers([n,1,0])
            markerA = self.getMarker()
            markerB = self.getMarker()
            #print "markerA=%s markerB=%s" %(markerA,markerB)
            #self.readMarkers([0,0])
            #print "subtable :) = ",foundSubTable
            if [markerA,markerB]==[0,0]:
                isFileDone = True
            ###
        except:
            pass
        ###
        self.n = nOld
        self.op2.seek(self.n)
        #print "isFileDone = ",isFileDone
        return isFileDone

    def checkForNextSubTable(self,n):
        foundSubTable = False
        #print "tell = ",self.op2.tell()
        
        try:
            nOld = self.op2.tell()
            self.readMarkers([n,1,0])
            foundSubTable=True
            #print "subtable :) = ",foundSubTable
        except InvalidMarkersError:
            foundSubTable = False
        ###
        self.n = nOld
        self.op2.seek(self.n)

        return foundSubTable

    def readGeomSubTable(self,iTable):
        i=0
        isNextTable=False
        isNextSubTable=False
        self.readMarkers([iTable,1,0])
        #print self.iTableMap

        tableName = self.readTableName(rewind=True,stopOnFailure=False)
        if tableName:
            #print "**tableName = |%r|" %(tableName)
            return tableName,isNextTable,isNextSubTable,False

        data = ''
        isTableActive=False
        while isNextSubTable==False and isNextTable==False:
            #print self.printSection(200)
            marker = self.getMarker()
            #print "marker = ",marker
            if marker<0:
                msg = 'marker is less than 0...'
                raise Exception(msg)
            data += self.readBlock()
            if not isTableActive:
                tableType = unpack('iii',data[:12])
                data = data[12:]

            #print "iTable=%s lenGeomData=%s" %(iTable,len(data))

            if tableType in self.iTableMap:
                #print "reading  iTable=%-3s with tableType=%s" %(iTable,tableType)
                self.iTableMap[tableType](data)
            else:
                print "skipping iTable=%-3s with tableType=%s" %(iTable,tableType)

            #self.op2Debug.write('ints = %s\n' %(str(ints)))

            isNextTable = self.checkForNextTable()
            isNextSubTable = self.checkForNextSubTable(iTable-1)
            isFileDone = self.checkFileDone(iTable-1)
            if isFileDone:
                isNextTable=True
            
            #print "i=%s tell=%s isNextTable=%s isNextSubTable=%s" %(i,self.op2.tell(),isNextTable,isNextSubTable)
            #if i==13:
            #    sys.exit('stopA')
            i+=1
            isTableActive=True
        ### while
        
        #print "exiting the geom sub table"
        return (tableName,isNextTable,isNextSubTable,isFileDone)

    def readTable_DYNAMICS(self):
        self.iTableMap = {}
        self.readRecordTable('DYNAMICS')
