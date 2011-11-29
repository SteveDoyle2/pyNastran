import os
import sys
import struct
from struct import unpack

#from pyNastran.op2.op2Errors import *
#from pyNastran.bdf.cards.nodes import GRID


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
            #CQUAD4(None,dataInit)
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
            #CROD(None,dataInit)
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
            dataInit = [eid,pid,n1,n2,n3,theta,zoffs,blank,tflag,t1,t2,t3]
            #CQUAD4(None,dataInit)
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
            #CTUBE(None,dataInit)
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
            #SPOINT(None,[nid])
        ###

