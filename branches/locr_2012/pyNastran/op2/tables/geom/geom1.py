import sys
from struct import unpack

from pyNastran.bdf.cards.nodes import GRID
from pyNastran.bdf.cards.coordinateSystems import (CORD1R, CORD1C, CORD1S,
                                                   CORD2R, CORD2C, CORD2S,
                                                   CORD3G)


class Geometry1(object):

    def readTable_Geom1(self):
        self.iTableMap = {
            (1701, 17, 6): self.readCord1C,  # record 1
            (1801, 18, 5): self.readCord1R,  # record 2
            (1901, 19, 7): self.readCord1S,  # record 3
            (2001, 20, 9): self.readCord2C,  # record 4
            (2101, 21, 8): self.readCord2R,  # record 5
            (2201, 22, 10): self.readCord2S,  # record 6
            #(14301,143,651): self.readCord3G, # record 7
            (4501, 45, 1): self.readGrid,   # record 17 - slow, but works
            (5301, 53, 4): self.readSEQGP,  # record 27 - not done

            (1101, 11, 66): self.readFake,  # record
            (3901, 39, 50): self.readFake,  # record
            (2201, 22, 10): self.readFake,  # record
            (6101, 61, 388): self.readFake,  # record
        }
        self.readRecordTable('GEOM1')

    def readTable_Geom1S(self):
        self.readTable_Geom1()
        #self.iTableMap = {
        #                 }
        #self.readRecordTable('GEOM1S')

    def readTable_Geom1N(self):
        self.iTableMap = {
        }
        self.readRecordTable('GEOM1N')

    def readCord1C(self, data):
        """
        (1701,17,6) - the marker for Record 1
        """
        #print "reading CORD1C"
        n = 0
        nEntries = len(data) // 24
        for i in xrange(nEntries):
            eData = data[n:n + 24]  # 6*4
            (cid, one, two, g1, g2, g3) = unpack(b'iiiiii', eData)
            dataIn = [cid, g1, g2, g3]
            coord = CORD1C(None, None, dataIn)
            self.addCoord(coord)
            n += 24
        data = data[n:]

    def readCord1R(self, data):
        """
        (1801,18,5) - the marker for Record 2
        """
        #print "reading CORD1R"
        n = 0
        nEntries = len(data) // 24
        for i in xrange(nEntries):
            eData = data[n:n + 24]  # 6*4
            (cid, one, one, g1, g2, g3) = unpack(b'iiiiii', eData)
            dataIn = [cid, g1, g2, g3]
            coord = CORD1R(None, None, dataIn)
            self.addCoord(coord)
            n += 24
        data = data[n:]

    def readCord1S(self, data):
        """
        (1901,19,7) - the marker for Record 3
        """
        #print "reading CORD1S"
        n = 0
        nEntries = len(data) // 24
        for i in xrange(nEntries):
            eData = data[n:n + 24]  # 6*4
            (cid, three, one, g1, g2, g3) = unpack(b'iiiiii', eData)
            dataIn = [cid, g1, g2, g3]
            coord = CORD1S(None, dataIn)
            self.addCoord(coord, allowOverwrites=True)
            n += 24
        data = data[n:]

    def readCord2C(self, data):
        """
        (2001,20,9) - the marker for Record 4
        """
        #print "reading CORD2C"
        n = 0
        nEntries = len(data) // 52
        for i in xrange(nEntries):
            eData = data[n:n + 52]  # 13*4
            (cid, two, two, rid, a1, a2, a3, b1, b2, b3, c1,
                c2, c3) = unpack(b'4i9f', eData)
            #print "cid=%s two=%s two=%s rid=%s a1=%s a2=%s a3=%s b1=%s b2=%s b3=%s c1=%s c2=%s c3=%s" %(cid,two,two,rid,a1,a2,a3,b1,b2,b3,c1,c2,c3)
            dataIn = [cid, rid, a1, a2, a3, b1, b2, b3, c1, c2, c3]
            coord = CORD2C(None, dataIn)
            self.addCoord(coord, allowOverwrites=True)
            n += 52
        data = data[n:]

    def readCord2R(self, data):
        """
        (2101,21,8) - the marker for Record 5
        """
        #print "reading CORD2R"
        n = 0
        nEntries = len(data) // 52
        for i in xrange(nEntries):
            eData = data[n:n + 52]  # 13*4
            (cid, one, two, rid, a1, a2, a3, b1, b2, b3, c1,
                c2, c3) = unpack(b'4i9f', eData)
            dataIn = [cid, rid, a1, a2, a3, b1, b2, b3, c1, c2, c3]
            #print "cid=%s one=%s two=%s rid=%s a1=%s a2=%s a3=%s b1=%s b2=%s b3=%s c1=%s c2=%s c3=%s" %(cid,one,two,rid,a1,a2,a3,b1,b2,b3,c1,c2,c3)
            coord = CORD2R(None, dataIn)
            self.addCoord(coord, allowOverwrites=True)
            n += 52
        data = data[n:]

    def readCord2S(self, data):
        """
        (2201,22,10) - the marker for Record 6
        """
        #print "reading CORD2S"
        n = 0
        nEntries = len(data) // 52
        for i in xrange(nEntries):
            eData = data[n:n + 52]  # 13*4
            (cid, sixty5, eight, rid, a1, a2, a3, b1, b2, b3,
                c1, c2, c3) = unpack(b'4i9f', eData)
            #print "cid=%s sixty5=%s eight=%s rid=%s a1=%s a2=%s a3=%s b1=%s b2=%s b3=%s c1=%s c2=%s c3=%s" %(cid,sixty5,eight,rid,a1,a2,a3,b1,b2,b3,c1,c2,c3)
            dataIn = [cid, rid, a1, a2, a3, b1, b2, b3, c1, c2, c3]
            coord = CORD2S(dataIn)
            self.addCoord(coord, allowOverwrites=True)
            n += 52
        data = data[n:]

    def readCord3G(self, data):
        """
        (14301,143,651) - the marker for Record 7
        @todo isnt this a CORD3G, not a CORD3R ???
        """
        #print "reading CORD3G"
        n = 0
        nEntries = len(data) // 16
        for i in xrange(nEntries):
            eData = data[n:n + 16]  # 4*4
            (cid, n1, n2, n3) = unpack(b'iiii', eData)
            dataIn = [cid, n1, n2, n3]
            coord = CORD3G(None, dataIn)
            self.addCoord(coord, allowOverwrites=True)
            n += 16
        data = data[n:]

    def readGrid(self, data):  # 21.8 sec, 18.9
        """(4501,45,1) - the marker for Record 17"""
        #print "reading GRID"
        n = 0
        nEntries = len(data) // 32
        for i in xrange(nEntries):
            eData = data[n:n + 32]
            out = unpack(b'iifffiii', eData)

            (nID, cp, x1, x2, x3, cd, ps, seid) = out
            if cd >= 0 and nID < 10000000:
                node = GRID(None, out)
                self.addNode(node)
            else:
                self.log.debug("*nID=%s cp=%s x1=%-5.2f x2=%-5.2f x3=%-5.2f cd=%-2s ps=%s seid=%s" % (nID, cp, x1, x2, x3, cd, ps, seid))
            #print str(grid)[:-1]
            n += 32
        data = data[n:]
        #assert len(data)==0,'len(data)!=0   len(data)=%s' %(len(data))
        #print "len(data) = ",len(data)

    def readSEQGP(self, data):
        """(5301,53,4) - the marker for Record 27"""
        self.skippedCardsFile.write('skipping SEQGP in GEOM1\n')
