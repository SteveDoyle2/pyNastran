import os
import sys
import struct
from struct import unpack

#from pyNastran.op2.op2Errors import *
from pyNastran.bdf.cards.loads import FORCE,FORCE1,FORCE2

class Geometry3(object):

    def readTable_Geom3(self):
        self.iTableMap = {
                         (4201,42,18): self.readForce,
                         (4001,40,20): self.readForce1,
                         (4101,41,22): self.readForce2,
                         }
        self.readRecordTable('GEOM3')

# ACCEL
# ACCEL1

    def readForce(self,data):
        """
        (4201,42,18) - the marker for Record 3
        """
        print "reading FORCE"
        while len(data)>=28: # 7*4
            eData = data[:28]
            data  = data[28:]
            (sid,g,cid,f,n1,n2,n3) = unpack('iiiffff',eData)

            load = FORCE(None,[sid,g,cid,f,n1,n2,n3])
            self.addLoad(load)
        ###

    def readForce1(self,data):
        """
        (4001,40,20) - the marker for Record 4
        """
        print "reading FORCE1"
        while len(data)>=16: # 4*4
            eData = data[:16]
            data  = data[16:]
            (sid,g,f,n1,n2) = unpack('iifii',eData)

            load = FORCE1(None,[sid,g,f,n1,n2])
            self.addLoad(load)
        ###

    def readForce2(self,data):
        """
        (4101,41,22) - the marker for Record 5
        """
        print "reading FORCE2"
        while len(data)>=28: # 7*4
            eData = data[:28]
            data  = data[28:]
            (sid,g,f,n1,n2,n3,n4) = unpack('iifiiii',eData)

            load = FORCE2(None,[sid,g,f,n1,n2,n3,n4])
            self.addLoad(load)
        ###

# GMLOAD
# GRAV
# LOAD
# LOADCYH
# LOADCYN
# LOADCYT
# LSEQ
# MOMENT
# MOMENT1
# MOMENT2
# PLOAD
# PLOAD1
# PLOAD2
# PLOAD3
# PLOAD4
# PLOADX
# PLOADX1
# PRESAX
# QBDY1
# QBDY2
# QBDY3
# QHBDY
# QVECT
# QVOL
# RFORCE
# SLOAD
# TEMP
# TEMPD
# TEMPEST
# TEMPF
# TEMP1C
# TEMPP1
# TEMPP2
# TEMPP3
# TEMPRB
# PFACE
# PEDGE

#geom3
