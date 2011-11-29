import os
import sys
import struct
from struct import unpack

#from pyNastran.op2.op2Errors import *
from pyNastran.bdf.cards.elementsShell import CTRIA3,CQUAD4


class Geometry2(object):

    def readTable_Geom2(self):
        self.iTableMap = {
                           (2958,51,177):    self.readCQUAD4, # record 69
                           (13900,139,9989): self.readCQUAD4, # record 70
                           (3001,30,48):     self.readCROD,   # record 93
                           (5959,59,282):    self.readCTRIA3, # record 93
                           (3701,37,49):     self.readCTUBE,  # record 103
                           (5551,49,105):    self.readSPOINT, # record 118
                         }
        self.readRecordTable('GEOM2')

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
            elem = CQUAD4(None,dataInit)
            self.addElement(elem)
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
            #elem = CROD(None,dataInit)
            #self.addElement(elem)
        ###

    def readCTRIA3(self,data):
        """
        (5959,59,282)    - the marker for Record 93
        """
        print "reading CTRIA3"
        while len(data)>=52: # 13*4
            eData = data[:52]
            data  = data[52:]
            (eid,pid,n1,n2,n3,theta,zoffs,blank1,blank2,tflag,t1,t2,t3) = unpack('iiiiiffiiifff',eData)
            #print "eid=%s pid=%s n1=%s n2=%s n3=%s theta=%s zoffs=%s blank1=%s blank2=%s tflag=%s t1=%s t2=%s t3=%s" %(eid,pid,n1,n2,n3,theta,zoffs,blank1,blank2,tflag,t1,t2,t3)
            dataInit = [eid,pid,n1,n2,n3,theta,zoffs,blank1,tflag,t1,t2,t3]
            elem = CTRIA3(None,dataInit)
            self.addElement(elem)
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
            #elem = CTUBE(None,dataInit)
            self.addElement(elem)
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
            #node = SPOINT(None,[nid])
            #self.addNode(node)
        ###

