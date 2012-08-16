import sys
from struct import unpack

from pyNastran.bdf.cards.loads.staticLoads import (FORCE, FORCE1, FORCE2, GRAV,
                                                   LOAD, PLOAD1, PLOAD2,
                                                   PLOAD4)  # PLOAD3,
from pyNastran.bdf.cards.thermal.loads import QBDY1, QBDY2, QBDY3, TEMP, TEMPD


class Geometry3(object):

    def readTable_Geom3(self):
        self.iTableMap = {
            (4201, 42, 18): self.readFORCE,   # record 3
            (4001, 40, 20): self.readFORCE1,  # record 4
            (4101, 41, 22): self.readFORCE2,  # record 5
            (4401, 44, 26): self.readGRAV,    # record 7 - buggy
            (4551, 61, 84): self.readLOAD,    # record 8
            (3709, 37, 331): self.readLOADCYH,  # record 9 - not done
            (3609, 36, 188): self.readLSEQ,    # record 12 - not done
            (4801, 48, 19): self.readMOMENT,  # record 13 - not tested
            (4601, 46, 21): self.readMOMENT1,  # record 14 - not tested
            (4701, 47, 23): self.readMOMENT2,  # record 15 - not tested
            (5101, 51, 24): self.readPLOAD,   # record 16 - not done
            #(6909,69,198): self.readPLOAD1,  # record 17 - buggy
            #(6802,68,199): self.readPLOAD2,  # record 18 - buggy
            (7109, 81, 255): self.readPLOAD3,  # record 19 - not done
            #(7209,72,299): self.readPLOAD4,  # record 20 - buggy - g1/g3/g4
            (7309, 73, 351): self.readPLOADX1,  # record 22
            (4509, 45, 239): self.readQBDY1,   # record 24
            (4909, 49, 240): self.readQBDY2,   # record 25
            (2109, 21, 414): self.readQBDY3,   # record 26
            (5509, 55, 190): self.readRFORCE,  # record 30 - not done
            (5401, 54, 25): self.readSLOAD,   # record 31 - not done
            (5701, 57, 27): self.readTEMP,    # record 32
            (5641, 65, 98): self.readTEMPD,   # record 33
            (8409, 84, 204): self.readTEMPRB,  # record 40 - not done
            (8109, 81, 201): self.readTEMPP1,  # record 37 - not done
            (8209, 82, 202): self.readTEMPP2,  # record 38 - not done
            (8309, 83, 203): self.readTEMPP3,  # record 39 - not done
            #(8409,84,204): self.readTEMPP4,  # record 40 - not done

        }
        self.readRecordTable('GEOM3')

    def readTable_Geom3S(self):
        self.readTable_Geom3()
        #self.iTableMap = {
        #                 }
        #self.readRecordTable('GEOM3S')

# ACCEL
# ACCEL1

    def readFORCE(self, data):
        """
        FORCE(4201,42,18) - the marker for Record 3
        """
        #print "reading FORCE"
        n = 0
        nEntries = len(data) // 28  # 7*4
        for i in xrange(nEntries):
            eData = data[n:n + 28]
            (sid, g, cid, f, n1, n2, n3) = unpack('iiiffff', eData)
            load = FORCE(None, [sid, g, cid, f, n1, n2, n3])
            self.addLoad(load)
            n += 28
        ###
        data = data[n:]

    def readFORCE1(self, data):
        """
        FORCE1(4001,40,20) - the marker for Record 4
        """
        #print "reading FORCE1"
        n = 0
        nEntries = len(data) // 20  # 5*4
        for i in xrange(nEntries):
            eData = data[n:n + 20]
            (sid, g, f, n1, n2) = unpack('iifii', eData)

            load = FORCE1(None, [sid, g, f, n1, n2])
            self.addLoad(load)
            n += 20
        ###
        data = data[n:]

    def readFORCE2(self, data):
        """
        FORCE2(4101,41,22) - the marker for Record 5
        """
        #print "reading FORCE2"
        n = 0
        nEntries = len(data) // 28  # 7*4
        for i in xrange(nEntries):
            eData = data[n:n + 28]
            (sid, g, f, n1, n2, n3, n4) = unpack('iifiiii', eData)

            load = FORCE2(None, [sid, g, f, n1, n2, n3, n4])
            self.addLoad(load)
            n += 28
        ###
        data = data[n:]

# GMLOAD

    def readGRAV(self, data):
        """
        GRAV(4401,44,26) - the marker for Record 7
        """
        #print "reading GRAV"
        n = 0
        nEntries = len(data) // 28  # 7*4
        for i in xrange(nEntries):
            eData = data[n:n + 28]
            out = unpack('iiffffi', eData)
            (sid, cid, a, n1, n2, n3, mb) = out
            grav = GRAV(None, out)
            self.addLoad(grav)
            n += 28
        ###
        data = data[n:]

    def readLOAD(self, data):
        """
        (4551, 61, 84) - the marker for Record 8
        @todo add object
        """
        #print "reading LOAD"
        while len(data) >= 16:  # 4*4
            eData = data[:16]
            data = data[16:]
            (sid, s, si, l1) = unpack('iffi', eData)
            Si = [si]
            L1 = [l1]
            #print Si,L1
            while 1:
                eData = data[:8]
                data = data[8:]
                #print "len(eData) = ",len(eData)
                (si, l1) = unpack('fi', eData)
                siTest, = unpack('i', eData[0:4])
                #print si,siTest,l1
                #print type(si)

                if [siTest, l1] == [-1, -1]:
                    break
                Si.append(si)
                L1.append(l1)
                #print Si,L1

            dataIn = [sid, s, Si, L1]
            load = LOAD(None, dataIn)
            self.addLoad(load)
        ###

    def readLOADCYH(self, data):
        self.skippedCardsFile.write('skipping LOADCYG in GEOM3\n')

# LOADCYN
# LOADCYT

    def readLSEQ(self, data):
        self.skippedCardsFile.write('skipping LSEQ in GEOM3\n')

    def readMOMENT(self, data):
        """
        MOMENT(4801,48,19) - the marker for Record 13
        """
        #print "reading MOMENT"
        n = 0
        nEntries = len(data) // 28  # 7*4
        for i in xrange(nEntries):
            eData = data[n:n + 28]
            out = unpack('iiiffff', eData)
            (sid, g, cid, m, n1, n2, n3) = out

            load = FORCE1(None, out)
            self.addLoad(load)
            n += 28
        ###
        data = data[n:]

    def readMOMENT1(self, data):
        """
        MOMENT1(4601,46,21) - the marker for Record 14
        """
        #print "reading MOMENT1"
        n = 0
        nEntries = len(data) // 20  # 5*4
        for i in xrange(nEntries):
            eData = data[n:n + 20]
            out = unpack('iifii', eData)
            (sid, g, m, n1, n2) = out
            load = FORCE1(None, out)
            self.addLoad(load)
            n += 20
        ###
        data = data[n:]

    def readMOMENT2(self, data):
        """
        MOMENT2(4701,47,23) - the marker for Record 15
        """
        #print "reading MOMENT2"
        n = 0
        nEntries = len(data) // 28  # 7*4
        for i in xrange(nEntries):
            eData = data[n:n + 28]
            out = unpack('iifiiii', eData)
            (sid, g, m, n1, n2, n3, n4) = out

            load = FORCE1(None, out)
            self.addLoad(load)
            n += 28
        ###
        data = data[n:]

    def readPLOAD(self, data):
        pass

    def readPLOAD1(self, data):
        """
        PLOAD2(6802,68,199) - the marker for Record 17
        """
        #print "reading PLOAD1"
        n = 0
        nEntries = len(data) // 32  # 8*4
        for i in xrange(nEntries):
            eData = data[n:n + 32]
            out = unpack('iiiiffff', eData)
            (sid, eid, Type, scale, x1, p1, x2, p2) = out
            #print "PLOAD1 = ",out
            load = PLOAD1(None, out)
            self.addLoad(load)
            n += 32
        ###
        data = data[n:]

    def readPLOAD2(self, data):
        """
        PLOAD2(6802,68,199) - the marker for Record 18
        """
        #print "reading PLOAD2"
        n = 0
        nEntries = len(data) // 12  # 3*4
        for i in xrange(nEntries):
            eData = data[n:n + 12]
            out = unpack('ifi', eData)
            (sid, p, eid) = out
            load = PLOAD2(None, out)
            self.addLoad(load)
            n += 12
        ###
        data = data[n:]

    def readPLOAD3(self, data):
        """
        PLOAD3(7109,71,255) - the marker for Record 19
        """
        #print "reading PLOAD3"
        n = 0
        nEntries = len(data) // 20  # 5*4
        for i in xrange(nEntries):
            eData = data[n:n + 20]
            out = unpack('if3i', eData)
            (sid, p, eid, n1, n2) = out
            load = PLOAD3(None, out)
            self.addLoad(load)
            n += 20
        ###
        data = data[n:]

    def readPLOAD4(self, data):  ## inconsistent with DMAP
        """
        PLOAD4(7209,72,299) - the marker for Record 20
        """
        #print "reading PLOAD4"
        n = 0
        nEntries = len(data) // 48  # 13*4
        for i in range(nEntries):
            eData = data[n:n + 48]
                         #iiffffiiifffi   ssssssssssssssss
            out = unpack('2i4f3i3f', eData)
            (sid, eid, p1, p2, p3, p4, g1, g34, cid, n1, n2, n3) = out
            #s1,s2,s3,s4,s5,s6,s7,s8,L1,L2,L3,L4,L5,L6,L7,L8
            #sdrlA = s1+s2+s3+s4
            #sdrlB = s5+s6+s7+s8
            #ldirA = L1+L2+L3+L4
            #ldirB = L5+L6+L7+L8
            sdrlA = None
            sdrlB = None
            ldirA = None
            ldirB = None
            load = PLOAD4(None, [sid, eid, [p1, p2, p3, p4], g1, g34,
                                 cid, [n1, n2, n3], sdrlA, sdrlB, ldirA, ldirB])
            self.addLoad(load)
            n += 48
        ###
        data = data[n:]

# PLOADX - obsolete
    def readPLOADX1(self, data):
        self.skippedCardsFile.write('skipping PLOADX1 in GEOM3\n')

# PRESAX

    def readQBDY1(self, data):
        """
        QBDY1(4509,45,239) - the marker for Record 24
        """
        #print "reading QBDY1"
        n = 0
        nEntries = len(data) // 12  # 3*4
        for i in range(nEntries):
            eData = data[n:n + 12]
            out = unpack('ifi', eData)
            (sid, q0, eid) = out
            load = QBDY1(None, out)
            self.addThermalLoad(load)
            n += 12
        ###
        data = data[n:]

    def readQBDY2(self, data):
        """
        QBDY2(4909,49,240) - the marker for Record 25
        """
        #print "reading QBDY2"
        n = 0
        nEntries = len(data) // 40  # 10*4
        for i in range(nEntries):
            eData = data[n:n + 40]
            out = unpack('iiffffffff', eData)
            (sid, eid, q1, q2, q3, q4, q5, q6, q7, q8) = out
            load = QBDY2(None, out)
            self.addThermalLoad(load)
            n += 40
        ###
        data = data[n:]

    def readQBDY3(self, data):
        """
        QBDY3(2109,21,414) - the marker for Record 26
        """
        #print "reading QBDY3"
        n = 0
        nEntries = len(data) // 16  # 4*4
        for i in range(nEntries):
            eData = data[n:n + 16]
            out = unpack('ifii', eData)
            (sid, q0, cntrlnd, eid) = out
            load = QBDY3(None, out)
            self.addThermalLoad(load)
            n += 16
        ###
        data = data[n:]

    def readTEMP(self, data):
        """
        TEMP(5701,57,27) - the marker for Record 32
        @warning buggy
        """
        #print "reading TEMP"
        n = 0
        nEntries = len(data) // 12  # 3*4
        for i in range(nEntries):
            eData = data[n:n + 12]
            out = unpack('iif', eData)
            (sid, g, T) = out
            if g < 10000000:
                load = TEMP(None, out)
                self.addThermalLoad(load)
            else:
                self.log.debug('TEMP = %s' % (out))
            n += 12
        ###
        data = data[n:]

    def readTEMPD(self, data):
        """
        TEMPD(5641,65,98) - the marker for Record 33
        @todo add object
        """
        #print "reading TEMPD"
        n = 0
        nEntries = len(data) // 8  # 2*4
        for i in range(nEntries):
            eData = data[n:n + 8]
            out = unpack('if', eData)
            (sid, T) = out
            load = TEMPD(None, out)
            #self.addThermalLoad(load)
            n += 8
        ###
        data = data[n:]

# QHBDY
# QVECT
# QVOL

    def readRFORCE(self, data):
        self.skippedCardsFile.write('skipping RFORCE in GEOM3\n')

    def readSLOAD(self, data):
        self.skippedCardsFile.write('skipping SLOAD in GEOM3\n')

# TEMP(5701,57,27) # 32
# TEMPD(5641,65,98) # 33
# TEMPEST
# TEMPF
# TEMP1C

    def readTEMPP1(self, data):
        self.skippedCardsFile.write('skipping TEMPP1 in GEOM3\n')

    def readTEMPP2(self, data):
        self.skippedCardsFile.write('skipping TEMPP2 in GEOM3\n')

    def readTEMPP3(self, data):
        self.skippedCardsFile.write('skipping TEMPP3 in GEOM3\n')

    def readTEMPRB(self, data):
        self.skippedCardsFile.write('skipping TEMPRB in GEOM3\n')

# PFACE
# PEDGE
