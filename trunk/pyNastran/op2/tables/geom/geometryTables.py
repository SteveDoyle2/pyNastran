import os
import sys
import struct
from struct import unpack

from pyNastran.op2.op2Errors import *
from pyNastran.op2.tables.geom.geom1 import Geometry1
from pyNastran.op2.tables.geom.geom2 import Geometry2
from pyNastran.op2.tables.geom.geom3 import Geometry3
from pyNastran.op2.tables.geom.geom4 import Geometry4
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
        if self.makeGeom==False:
            self.iTableMap = {}

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
                if not(tableType[0]==tableType[1]==tableType[2]):
                    if self.makeGeom:
                        print "skipping %s iTable=%-3s with tableType=%s" %(self.tableName,iTable,tableType)
                    ###
                ###
            ###
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

    def readTable_PCOMPTS(self):
        self.iTableMap = {
                         }
        self.readRecordTable('PCOMPTS')

    def readTable_DUMMY_GEOM(self,tableName):
        self.readRecordTable(tableName)

    def readTable_DIT(self):
        self.iTableMap = {
                            (1105,11,133): self.readFake,
                            (105,  1, 93): self.readFake,
                            (15, 21, 162): self.readFake,
                            (56, 26, 303): self.readFake,
                            (3105,31, 97): self.readFake,
                         }
        self.readRecordTable('DIT')

    def readTable_DYNAMICS(self):
        self.iTableMap = {
                            (37, 18, 183): self.readFake,
                            (57,   5,123): self.readFake,
                            (107,  1, 86): self.readFake,
                            (207,  2, 87): self.readFake,
                            (307,  3, 85): self.readFake,
                            (308,  8,348): self.readFake,
                            (707,  7,124): self.readFake,
                            (1007,10,125): self.readFake,
                            (1307,13,126): self.readFake,
                            (3107,31,127): self.readFake,
                            (5107,51,131): self.readFake,
                            (5207,52,132): self.readFake,
                            (6207,62,136): self.readFake,
                            (6607,66,137): self.readFake,
                            (7107,71,138): self.readFake,
                            (7207,72,139): self.readFake,
                            (8307,83,142): self.readFake,
                            (2107,21,195): self.readFake,
                            (2207,22,196): self.readFake,
                            
                            
                         }
        self.readRecordTable('DYNAMICS')
