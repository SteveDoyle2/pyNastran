import os
import sys
import struct
from struct import unpack

#from pyNastran.op2.op2Errors import *
from pyNastran.bdf.cards.aero    import GRAV
from pyNastran.bdf.cards.loads   import *
from pyNastran.bdf.cards.thermal import *

class Geometry3(object):

    def readTable_Geom3(self):
        self.iTableMap = {
                         (4201,42,18): self.readFORCE,  # record 3
                         (4001,40,20): self.readFORCE1, # record 4
                         (4101,41,22): self.readFORCE2, # record 5
                         (4401,44,26): self.readGRAV,   # record 7
                         (4551,61,84): self.readLOAD,   # record 8

                         #(7209,72,299): self.readPLOAD4,# record 20 - buggy
                         (4509,45,239): self.readQBDY1, # record 24
                         (4909,49,240): self.readQBDY2, # record 25
                         (2109,21,414): self.readQBDY3, # record 26
                         #(5641,65,98):  self.readTEMP,  # record 32 - buggy
                         (5701,57,27):  self.readTEMPD, # record 33
                         (4801,48,19):  self.readFake,
                         (4701,47,23):  self.readFake,
                         (6909,69,198): self.readFake,
                         (8409,84,204): self.readFake,
                         (8109,81,201): self.readFake,
                         (6802,68,199): self.readFake,
                         (3709,37,331): self.readFake,
                         (5401,54,25):  self.readFake,
                         (7309,73,351): self.readFake,
                         (3609,36,188): self.readFake,
                         (5101,51,24):  self.readFake,

                         (7209,72,299): self.readFake, # PLOAD4
                         (5641,65,98):  self.readFake, # TEMP
                         #(5408,54,261)

                         }
        self.readRecordTable('GEOM3')


# ACCEL
# ACCEL1

    def readFORCE(self,data):
        """
        FORCE(4201,42,18) - the marker for Record 3
        """
        print "reading FORCE"
        while len(data)>=28: # 7*4
            eData = data[:28]
            data  = data[28:]
            (sid,g,cid,f,n1,n2,n3) = unpack('iiiffff',eData)

            load = FORCE(None,[sid,g,cid,f,n1,n2,n3])
            self.addLoad(load)
        ###

    def readFORCE1(self,data):
        """
        FORCE1(4001,40,20) - the marker for Record 4
        """
        print "reading FORCE1"
        n=0
        nEntries = len(data)/20
        for i in range(0,nEntries):
            eData = data[n:n+20]  # 5*4
            (sid,g,f,n1,n2) = unpack('iifii',eData)

            load = FORCE1(None,[sid,g,f,n1,n2])
            self.addLoad(load)
            n+=20
        ###
        data = data[n:]

    def readFORCE2(self,data):
        """
        FORCE2(4101,41,22) - the marker for Record 5
        """
        print "reading FORCE2"
        n=0
        nEntries = len(data)/28
        for i in range(0,nEntries):
            eData = data[n:n+28]  # 7*4
            (sid,g,f,n1,n2,n3,n4) = unpack('iifiiii',eData)

            load = FORCE2(None,[sid,g,f,n1,n2,n3,n4])
            self.addLoad(load)
            n+=28
        ###
        data = data[n:]

# GMLOAD

    def readGRAV(self,data):
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

    def readLOAD(self,data):
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

    def readPLOAD4(self,data): ## inconsistent with DMAP
        """
        PLOAD4(7209,72,299) - the marker for Record 20
        """
        print "reading PLOAD4"
        while len(data)>=48: # 13*4
            eData = data[:48]
            data  = data[48:]
                         #iiffffiiifffi   ssssssssssssssss
            out = unpack('iiffffiiifff',eData)
            (sid,eid,p1,p2,p3,p4,g1,g34,cid,n1,n2,n3) = out
            #s1,s2,s3,s4,s5,s6,s7,s8,L1,L2,L3,L4,L5,L6,L7,L8
            #sdrlA = s1+s2+s3+s4
            #sdrlB = s5+s6+s7+s8
            #ldirA = L1+L2+L3+L4
            #ldirB = L5+L6+L7+L8
            sdrlA = None
            sdrlB = None
            ldirA = None
            ldirB = None
            load = PLOAD4(None,[sid,eid,[p1,p2,p3,p4],g1,g34,cid,[n1,n2,n3],sdrlA,sdrlB,ldirA,ldirB])
            self.addLoad(load)
        ###

# PLOADX
# PLOADX1
# PRESAX
# QBDY1
# QBDY2


    def readQBDY1(self,data):
        """
        QBDY1(4509,45,239) - the marker for Record 24
        """
        print "reading QBDY1"
        while len(data)>=12: # 3*4
            eData = data[:12]
            data  = data[12:]
            out = unpack('ifi',eData)
            (sid,q0,eid) = out
            load = QBDY1(None,out)
            self.addThermalLoad(load)
        ###

    def readQBDY2(self,data):
        """
        QBDY2(4909,49,240) - the marker for Record 25
        """
        print "reading QBDY2"
        while len(data)>=40: # 10*4
            eData = data[:40]
            data  = data[40:]
            out = unpack('iiffffffff',eData)
            (sid,eid,q1,q2,q3,q4,q5,q6,q7,q8) = out
            load = QBDY2(None,out)
            self.addThermalLoad(load)
        ###

    def readQBDY3(self,data):
        """
        QBDY3(2109,21,414) - the marker for Record 26
        """
        print "reading QBDY3"
        while len(data)>=16: # 4*4
            eData = data[:16]
            data  = data[16:]
            out = unpack('ifii',eData)
            (sid,q0,cntrlnd,eid) = out
            load = QBDY3(None,out)
            self.addThermalLoad(load)
        ###

    def readTEMP(self,data):
        """
        TEMP(5701,57,27) - the marker for Record 32
        """
        print "reading TEMP"
        while len(data)>=12: # 3*4
            eData = data[:12]
            data  = data[12:]
            out = unpack('iif',eData)
            (sid,g,T) = out
            load = TEMP(None,out)
            self.addThermalLoad(load)
        ###

    def readTEMPD(self,data):
        """
        TEMPD(5641,65,98) - the marker for Record 33
        @todo add object
        """
        print "reading TEMPD"
        while len(data)>=8: # 4*4
            eData = data[:8]
            data  = data[8:]
            out = unpack('if',eData)
            (sid,T) = out
            load = TEMPD(None,out)
            #self.addThermalLoad(load)
        ###

# QHBDY
# QVECT
# QVOL
# RFORCE
# SLOAD
# TEMP(5701,57,27) # 32
# TEMPD(5641,65,98) # 33
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
