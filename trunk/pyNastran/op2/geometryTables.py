import os
import sys
import struct
from struct import unpack

from pyNastran.op2.op2Errors import *
from pyNastran.op2.geom1 import Geometry1
from pyNastran.op2.geom2 import Geometry2
#from pyNastran.op2.geom3 import Geometry3
#from pyNastran.op2.geom4 import Geometry4

class GeomObj(object):
    def __init__(self):
        pass
    def geomFunc(self,data):
        pass

class GeometryTables(Geometry1,Geometry2):

    def readRecordTable(self,expectedTableName):
        """
        @note assumes self.iTableMap has already been set
        """
        tableName = self.readTableName(rewind=False) # GEOM1
        self.tableInit(tableName)
        print "*tableName = |%r|" %(tableName)

        self.readMarkers([-1,7])
        fields = self.readIntBlock()
        #print "fields = ",fields

        self.readMarkers([-2,1,0]) # 2
        bufferWords = self.getMarker()
        #print "bufferWords = ",bufferWords,bufferWords*4
        word = self.readStringBlock()

        iTable = -3
        while 1:  ## @todo could this cause an infinite loop...i dont this so...
            (tableName,isNextTable,isNextSubTable) = self.readGeomSubTable(iTable)
        
            if self.checkForNextTable():
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
            print "**tableName = |%r|" %(tableName)
            return tableName,isNextTable,isNextSubTable

        data = ''
        isTableActive=False
        while isNextSubTable==False and isNextTable==False:
            #print self.printSection(100)
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
            #print "i=%s tell=%s isNextTable=%s isNextSubTable=%s" %(i,self.op2.tell(),isNextTable,isNextSubTable)
            #if i==13:
            #    sys.exit('stopA')
            i+=1
            isTableActive=True
        ### while
        
        #print "exiting the geom sub table"
        return (tableName,isNextTable,isNextSubTable)

    def readTable_Geom3(self):
        self.iTableMap = {}
        self.readRecordTable('GEOM3')

    def readTable_Geom4(self):
        self.iTableMap = {}
        self.readRecordTable('GEOM4')

    def readTable_EPT(self):
        self.iTableMap = {}
        self.readRecordTable('EPT')

    def readTable_MPTS(self):
        self.iTableMap = {}
        self.readRecordTable('MPTS')

    def readTable_DYNAMICS(self):
        self.iTableMap = {}
        self.readRecordTable('DYNAMICS')
