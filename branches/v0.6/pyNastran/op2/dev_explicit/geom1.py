#pylint: disable=C0301,C0103,W0612,R0914,C0326
import StringIO
from struct import unpack, Struct

from pyNastran.bdf.cards.nodes import GRID
from pyNastran.bdf.cards.coordinateSystems import (CORD1R, CORD1C, CORD1S,
                                                   CORD2R, CORD2C, CORD2S,
                                                   CORD3G)

class GEOM1(object):
    def add_node(self, node, allowOverwrites=True):
        raise RuntimeError('this should be overwritten')
    def add_coord(self, coord, allowOverwrites=True):
        raise RuntimeError('this should be overwritten')
    def readFake(self, data, n):
        raise RuntimeError('this should be overwritten')

    def __init__(self):
        self.skippedCardsFile = StringIO.StringIO()
        self.card_count = {}
        self._geom1_map = {
            (1701,  17,  6): ['CORD1C', self.readCord1C],  # record 1
            (1801,  18,  5): ['CORD1R', self.readCord1R],  # record 2
            (1901,  19,  7): ['CORD1S', self.readCord1S],  # record 3
            (2001,  20,  9): ['CORD2C', self.readCord2C],  # record 4
            (2101,  21,  8): ['CORD2R', self.readCord2R],  # record 5
            (2201,  22, 10): ['CORD2S', self.readCord2S],  # record 6
            (14301,143,651): ['CORD3G', self.readCord3G],  # record 7
            (4501,  45,  1): ['GRID',   self.readGrid],    # record 17 - slow, but works
            (5301,  53,  4): ['SEQGP',  self.readSEQGP],   # record 27 - not done

            (1101,  11,  66): ['', self.readFake],  # record
            (2201,  22,  10): ['', self.readFake],  # record
            (3901,  39,  50): ['', self.readFake],  # record
            (5501, 55, 297): ['', self.readFake],  # record
            (5601, 56, 296): ['', self.readFake],  # record
            (6001,  60, 377): ['', self.readFake],  # record
            (6101,  61, 388): ['', self.readFake],  # record
            (6601, 66, 392): ['', self.readFake],  # record
            (7902,  79, 302): ['', self.readFake],  # record
            (13301, 133, 509): ['', self.readFake],  # record
            (1227, 12, 462): ['', self.readFake],  # record
            (6201, 62, 389): ['', self.readFake],  # record
            (10101, 101, 394): ['', self.readFake],  # record
            (6401, 64, 402): ['', self.readFake],    # record
            (1427, 14, 465): ['', self.readFake],  # record
            (5678, 71, 475): ['', self.readFake],  # record
        }
    def readCord1C(self, data, n):
        """
        (1701,17,6) - the marker for Record 1
        """
        s = Struct(b'6i')
        nEntries = (len(data) - n) // 24
        for i in xrange(nEntries):
            eData = data[n:n + 24]  # 6*4
            (cid, one, two, g1, g2, g3) = s.unpack(eData)
            assert one in [1, 2], one
            assert two in [1, 2], two
            dataIn = [cid, g1, g2, g3]
            coord = CORD1C(None, None, dataIn)
            self.add_coord(coord)
            n += 24
        self.card_count['CORD1C'] = nEntries
        return n

    def readCord1R(self, data, n):
        """
        (1801,18,5) - the marker for Record 2
        """
        s = Struct(b'6i')
        nEntries = (len(data) - n) // 24
        for i in xrange(nEntries):
            eData = data[n:n + 24]  # 6*4
            (cid, one1, one2, g1, g2, g3) = s.unpack(eData)
            assert one1 == 1, one1
            assert one2 == 1, one2
            dataIn = [cid, g1, g2, g3]
            coord = CORD1R(None, None, dataIn)
            self.add_coord(coord)
            n += 24
        self.card_count['CORD1R'] = nEntries
        return n

    def readCord1S(self, data, n):
        """
        (1901,19,7) - the marker for Record 3
        """
        s = Struct(b'6i')
        nEntries = (len(data) - n) // 24
        for i in xrange(nEntries):
            edata = data[n:n + 24]  # 6*4
            (cid, three, one, g1, g2, g3) = s.unpack(edata)
            assert three == 3, three
            assert one == 1, one
            dataIn = [cid, g1, g2, g3]
            coord = CORD1S(None, dataIn)
            self.add_coord(coord, allowOverwrites=True)
            n += 24
        self.card_count['CORD1S'] = nEntries
        return n

    def readCord2C(self, data, n):
        """
        (2001,20,9) - the marker for Record 4
        """
        s = Struct(b'4i9f')
        nEntries = (len(data) - n) // 52
        for i in xrange(nEntries):
            eData = data[n:n + 52]  # 13*4
            (cid, two1, two2, rid, a1, a2, a3, b1, b2, b3, c1, c2, c3) = s.unpack(eData)
            assert two1 == 2, two1
            assert two2 == 2, two2
            dataIn = [cid, rid, a1, a2, a3, b1, b2, b3, c1, c2, c3]
            coord = CORD2C(None, dataIn)
            self.add_coord(coord, allowOverwrites=True)
            n += 52
        self.card_count['CORD2C'] = nEntries
        return n

    def readCord2R(self, data, n):
        """
        (2101,21,8) - the marker for Record 5
        """
        nEntries = (len(data) - n) // 52
        for i in xrange(nEntries):
            eData = data[n:n + 52]  # 13*4
            (cid, one, two, rid, a1, a2, a3, b1, b2, b3, c1,
                c2, c3) = unpack(b'4i9f', eData)
            assert one == 1, one
            assert two == 2, two
            dataIn = [cid, rid, a1, a2, a3, b1, b2, b3, c1, c2, c3]
            #print "cid=%s one=%s two=%s rid=%s a1=%s a2=%s a3=%s b1=%s b2=%s b3=%s c1=%s c2=%s c3=%s" %(cid,one,two,rid,a1,a2,a3,b1,b2,b3,c1,c2,c3)
            coord = CORD2R(None, dataIn)
            self.add_coord(coord, allowOverwrites=True)
            n += 52
        self.card_count['CORD2R'] = nEntries
        return n

    def readCord2S(self, data, n):
        """
        (2201,22,10) - the marker for Record 6
        """
        s = Struct(b'4i9f')
        nEntries = (len(data) - n) // 52
        for i in xrange(nEntries):
            eData = data[n:n + 52]  # 13*4
            (cid, sixty5, eight, rid, a1, a2, a3, b1, b2, b3,
                c1, c2, c3) = s.unpack(eData)
            dataIn = [cid, rid, a1, a2, a3, b1, b2, b3, c1, c2, c3]
            coord = CORD2S(dataIn)
            self.add_coord(coord, allowOverwrites=True)
            n += 52
        self.card_count['CORD2S'] = nEntries
        return n

    def readCord3G(self, data, n):
        """
        (14301,143,651) - the marker for Record 7
        .. todo:: isnt this a CORD3G, not a CORD3R ???
        """
        s = Struct(b'4i')
        nEntries = (len(data) - n) // 16
        for i in xrange(nEntries):
            eData = data[n:n + 16]  # 4*4
            (cid, n1, n2, n3) = s.unpack(eData)
            dataIn = [cid, n1, n2, n3]
            coord = CORD3G(None, dataIn)
            self.add_coord(coord, allowOverwrites=True)
            n += 16
        self.card_count['CORD3G'] = nEntries
        return n

    def readGrid(self, data, n):  # 21.8 sec, 18.9
        """(4501,45,1) - the marker for Record 17"""
        s = Struct(b'ii3f3i')
        ntotal = 32
        nentries = (len(data) - n) // ntotal
        for i in xrange(nentries):
            edata = data[n:n + 32]
            out = s.unpack(edata)
            (nID, cp, x1, x2, x3, cd, ps, seid) = out
            if cd >= 0 and nID < 10000000:
                node = GRID(None, out)
                self.add_node(node)
            else:
                self.log.debug("*nID=%s cp=%s x1=%-5.2f x2=%-5.2f x3=%-5.2f cd=%-2s ps=%s seid=%s" % (nID, cp, x1, x2, x3, cd, ps, seid))
            n += ntotal
        self.card_count['GRID'] = nentries
        return n

    def readSEQGP(self, data, n):
        """(5301,53,4) - the marker for Record 27"""
        self.skippedCardsFile.write('skipping SEQGP in GEOM1\n')
        return n
