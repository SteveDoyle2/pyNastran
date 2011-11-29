import os
import sys
import struct
from struct import unpack

from pyNastran.op2.op2Errors import *
from pyNastran.op2.geom1 import Geometry1

class GeomObj(object):
    def __init__(self):
        pass
    def geomFunc(self,data):
        pass

class GeometryTables(Geometry1):

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

#-----
# GEOM2

    def readCQUAD4(self,data):
        """
        (2958,51,177)    - the marker for Record 69
        (13900,139,9989) - the marker for Record 70
        """
        print "reading CQUAD4"
        while len(data)>=56: # 14*4
            eData = data[:56]
            data  = data[56:]
            (eid,pid,n1,n2,n3,n4,theta,zoffs,blank,tflag,t1,t2,t3,t4) = unpack('iiiiiiffiiffff',eData)
            dataInit = [eid,pid,n1,n2,n3,n4,theta,zoffs,blank,tflag,t1,t2,t3,t4]
            CQUAD4(None,dataInit)
        ###

    def readCROD(self,data):
        """
        (3001,30,48)    - the marker for Record 93
        """
        print "reading CROD"
        while len(data)>=16: # 4*4
            eData = data[:16]
            data  = data[16:]
            (eid,pid,n1,n2) = unpack('iiii',eData)
            dataInit = [eid,pid,n1,n2]
            CROD(None,dataInit)
        ###

    def readCTRIA3(self,data):
        """
        (5959,59,282)    - the marker for Record 93
        """
        print "reading CTRIA3"
        while len(data)>=48: # 12*4
            eData = data[:48]
            data  = data[48:]
            (eid,pid,n1,n2,n3,theta,zoffs,blank,tflag,t1,t2,t3) = unpack('iiiiiffiifff',eData)
            dataInit = [eid,pid,n1,n2,n3,n4,theta,zoffs,blank,tflag,t1,t2,t3,t4]
            CQUAD4(None,dataInit)
        ###

    def readCTUBE(self,data):
        """
        (3701,37,49)    - the marker for Record 103
        """
        print "reading CTUBE"
        while len(data)>=16: # 4*4
            eData = data[:16]
            data  = data[16:]
            (eid,pid,n1,n2) = unpack('iiii',eData)
            dataInit = [eid,pid,n1,n2]
            CTUBE(None,dataInit)
        ###

    def readSPOINT(self,data):
        """
        (5551,49,105)    - the marker for Record 118
        """
        print "reading SPOINT"
        while len(data)>=4: # 4*4
            eData = data[:4]
            data  = data[4:]
            (nid) = unpack('i',eData)
            SPOINT(None,[nid])
        ###

    def readTable_Geom2(self):
        self.iTableMap = {}
        self.readRecordTable('GEOM2')

# GEOM2
#-----
# GEOM3

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
