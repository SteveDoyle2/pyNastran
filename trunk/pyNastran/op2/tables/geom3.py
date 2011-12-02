import os
import sys
import struct
from struct import unpack

#from pyNastran.op2.op2Errors import *
from pyNastran.bdf.cards.aero  import GRAV
from pyNastran.bdf.cards.loads import FORCE,FORCE1,FORCE2,LOAD

class Geometry3(object):

    def readTable_Geom3(self):
        self.iTableMap = {
                         (4201,42,18): self.readForce,  # record 3
                         (4001,40,20): self.readForce1, # record 4
                         (4101,41,22): self.readForce2, # record 5
                         (4401,44,26): self.readGrav,   # record 7
                         (4551,61,84): self.readLoad,   # record 8

                         #(2109,21,414): self.read
                         #(5641,65,98):  self.read
                         #(5701,57,27):  self.read
                         #(5641,65,98):  self.read

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

    def readGrav(self,data):
        """
        GRAV(4401,44,26) - the marker for Record 7
        @todo add object
        """
        print "reading GRAV"
        while len(data)>=42: # 7*4
            eData = data[:42]
            data  = data[42:]
            out = unpack('iiffffi',eData)
            (sid,cid,a,n1,n2,n3,mb) = out
            grav = GRAV(None,out)
            self.addGrav(grav)
        ###

    def readLoad(self,data):
        """
        (4551, 61, 84) - the marker for Record 8
        @todo add object
        """
        print "reading LOAD"
        while len(data)>=16: # 4*4
            eData = data[:16]
            data  = data[16:]
            (sid,s,si,l1) = unpack('iffi',eData)
            Si=[si]; L1=[l1]
            #print Si,L1
            while 1:
                eData = data[:8]
                data  = data[8:]
                #print "len(eData) = ",len(eData)
                (si,l1) = unpack('fi',eData)
                siTest, = unpack('i',eData[0:4])
                #print si,siTest,l1
                #print type(si)
                
                if [siTest,l1] == [-1,-1]:
                    break
                Si.append(si); L1.append(l1)
                #print Si,L1
            
            dataIn = [sid,s,Si,L1]
            load = LOAD(None,dataIn)
            self.addLoad(load)
        ###

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
