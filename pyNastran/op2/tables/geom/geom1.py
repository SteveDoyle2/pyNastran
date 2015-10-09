#pylint: disable=C0301,C0103,W0612,R0914,C0326
from six.moves import range
from struct import unpack, Struct

from pyNastran.bdf.cards.nodes import GRID
from pyNastran.bdf.cards.coordinateSystems import (CORD1R, CORD1C, CORD1S,
                                                   CORD2R, CORD2C, CORD2S,
                                                   CORD3G)

class GEOM1(object):
    def _is_same_fields(self, fields1, fields2):
        for (field1, field2) in zip(fields1, fields2):
            if not is_same(field1, field2):
                return False
        return True

    def add_node(self, node, allowOverwrites=False):
        """GRDSET creates duplicate nodes...what about duplicate nodes?"""
        key = node.nid
        assert key > 0, 'nid=%s node=%s' % (key, node)
        if key in self.nodes:
            fields1 = self.nodes[key].raw_fields()
            fields2 = node.raw_fields()
            #grid, nid, cp, x1, x2, x3, cd, ps, seid
            for i, (v1, v2) in enumerate(zip(fields1, fields2)):
                if v1 != v2:
                    print('i=%s v1=%r v2=%r fields1=%s\nfields2=%s' % (
                          i, v1, v2, fields1, fields2))
        else:
            self._type_to_id_map[node.type].append(key)
        self.nodes[key] = node

    def add_coord(self, coord, allowOverwrites=True):
        raise RuntimeError('this should be overwritten')
    def _readFake(self, data, n):
        return len(data)

    def _read_geom1_4(self, data):
        return self._read_geom_4(self._geom1_map, data)

    def __init__(self):
        self.card_count = {}
        self._geom1_map = {
            (1701,  17,  6): ['CORD1C', self._readCord1C],  # record 1
            (1801,  18,  5): ['CORD1R', self._readCord1R],  # record 2
            (1901,  19,  7): ['CORD1S', self._readCord1S],  # record 3
            (2001,  20,  9): ['CORD2C', self._readCord2C],  # record 4
            (2101,  21,  8): ['CORD2R', self._readCord2R],  # record 5
            (2201,  22, 10): ['CORD2S', self._readCord2S],  # record 6
            (14301,143,651): ['CORD3G', self._readCord3G],  # record 7

            (4501,  45,  1): ['GRID',   self._readGrid],    # record 17
            (5301,  53,  4): ['SEQGP',  self._readSEQGP],   # record 27 - not done

            (2301,  23, 304): ['CSUPER',  self._readFake],  # record 8
            (5501,  55, 297): ['CSUPEXT', self._readFake],  # record 9
            (1627,  16, 463): ['EXTRN',   self._readFake],  # record 10
            (6101,  61, 388): ['FEEDGE',  self._readFake],  # record 11
            (6601,  66, 392): ['GMCURVE', self._readFake],  # record 12
            (6201,  62, 389): ['FEFACE',  self._readFake],  # record 13
            (6001,  60, 377): ['POINT',   self._readFake],  # record 14
            (10101,101, 394): ['GMSURF',  self._readFake],  # record 15
            (6401,  64, 402): ['GMCORD',  self._readFake],  # record 16
            # 17 - GRID  (above)
            (1527, 15, 466): ['SEBNDRY', self._readFake],  # record 18
            (1427, 14, 465): ['SEBULK',  self._readFake],  # record 19
            (427,   4, 453): ['SECONCT', self._readFake],  # record 20

            (7902, 79, 302): ['SEELT',   self._readFake],  # record 21
            (527,  72, 454): ['SEEXCLD', self._readFake],  # record 22
            (1027, 10, 459): ['SELABEL', self._readFake],  # record 23
            (827,   8, 457): ['SELOC',   self._readFake],  # record 24
            (927,   9, 458): ['SEMPLN',  self._readFake],  # record 25
            (1327, 13, 464): ['SENQSET', self._readFake],  # record 26
            # 27 - SEQGP (above)
            (5401, 54, 305): ['SEQSEP',  self._readFake],  # record 28
            (5601, 56, 296): ['SESET',   self._readFake],  # record 29
            (1227, 12, 462): ['SETREE',  self._readFake],  # record 30
            (5678, 71, 475): ['SNORM',   self._readFake],  # record 31
            (5701, 57, 323): ['CSUPER1', self._readFake],  # record 32

            (5801,   58, 324): ['SUPUP', self._readFake],  # record 33 - CSUPUP in NX; SUPUP in MSC
            (14101, 141, 403): ['SWLDPRM', self._readFake],  # record 34

            (1101,   11,  66): ['', self._readFake],  # record
            (2201,   22,  10): ['', self._readFake],  # record
            (3901,   39,  50): ['', self._readFake],  # record
            (13301, 133, 509): ['', self._readFake],  # record
            (1127,   11, 461) : ['SELOAD', self._readFake],  # record NX
            #(4501, 45, 1120001) : ['', self._readFake],  # record
        }

    def _readCord1C(self, data, n):
        """
        (1701,17,6) - the marker for Record 1
        """
        s = Struct(self._endian + b'6i')
        nentries = (len(data) - n) // 24
        for i in range(nentries):
            eData = data[n:n + 24]  # 6*4
            out = s.unpack(eData)
            (cid, one, two, g1, g2, g3) = out
            assert one in [1, 2], one
            assert two in [1, 2], two
            self.binary_debug.write('  CORD1C=%s\n' % str(out))
            dataIn = [cid, g1, g2, g3]
            coord = CORD1C(None, None, dataIn)
            self.add_coord(coord)
            n += 24
        self._increase_card_count('CORD1C', nentries)
        return n

    def _readCord1R(self, data, n):
        """
        (1801,18,5) - the marker for Record 2
        """
        s = Struct(self._endian + b'6i')
        nentries = (len(data) - n) // 24
        for i in range(nentries):
            eData = data[n:n + 24]  # 6*4
            out = s.unpack(eData)
            (cid, one1, one2, g1, g2, g3) = out
            self.binary_debug.write('  CORD1R=%s\n' % str(out))
            assert one1 == 1, one1
            assert one2 == 1, one2
            dataIn = [cid, g1, g2, g3]
            coord = CORD1R(None, None, dataIn)
            self.add_coord(coord)
            n += 24
        self._increase_card_count('CORD1R', nentries)
        return n

    def _readCord1S(self, data, n):
        """
        (1901,19,7) - the marker for Record 3
        """
        s = Struct(self._endian + b'6i')
        nentries = (len(data) - n) // 24
        for i in range(nentries):
            edata = data[n:n + 24]  # 6*4
            out = s.unpack(edata)
            (cid, three, one, g1, g2, g3) = out
            self.binary_debug.write('  CORD1S=%s\n' % str(out))
            assert three == 3, three
            assert one == 1, one
            data_in = [cid, g1, g2, g3]
            coord = CORD1S(None, None, data_in)
            self.add_coord(coord, allowOverwrites=True)
            n += 24
        self._increase_card_count('CORD1S', nentries)
        return n

    def _readCord2C(self, data, n):
        """
        (2001,20,9) - the marker for Record 4
        """
        s = Struct(self._endian + b'4i9f')
        nentries = (len(data) - n) // 52
        for i in range(nentries):
            eData = data[n:n + 52]  # 13*4
            out = s.unpack(eData)
            (cid, two1, two2, rid, a1, a2, a3, b1, b2, b3, c1, c2, c3) = out
            assert two1 == 2, two1
            assert two2 == 2, two2
            dataIn = [cid, rid, a1, a2, a3, b1, b2, b3, c1, c2, c3]
            coord = CORD2C(None, dataIn)
            self.binary_debug.write('  CORD2C=%s\n' % str(out))
            self.add_coord(coord, allowOverwrites=True)
            n += 52
        self._increase_card_count('CORD2C', nentries)
        return n

    def _readCord2R(self, data, n):
        """
        (2101,21,8) - the marker for Record 5
        """
        nentries = (len(data) - n) // 52
        for i in range(nentries):
            eData = data[n:n + 52]  # 13*4
            (cid, one, two, rid, a1, a2, a3, b1, b2, b3, c1,
                c2, c3) = unpack(self._endian + b'4i9f', eData)
            assert one == 1, one
            assert two == 2, two
            dataIn = [cid, rid, a1, a2, a3, b1, b2, b3, c1, c2, c3]
            #print("cid=%s rid=%s a1=%s a2=%s a3=%s b1=%s b2=%s b3=%s c1=%s c2=%s c3=%s" %(cid,rid,a1,a2,a3,b1,b2,b3,c1,c2,c3))
            self.binary_debug.write('  CORD2R=%s\n' % dataIn)
            coord = CORD2R(None, dataIn)
            self.add_coord(coord, allowOverwrites=True)
            n += 52
        self._increase_card_count('CORD2R', nentries)
        return n

    def _readCord2S(self, data, n):
        """
        (2201,22,10) - the marker for Record 6
        """
        s = Struct(self._endian + b'4i9f')
        nentries = (len(data) - n) // 52
        for i in range(nentries):
            eData = data[n:n + 52]  # 13*4
            out = s.unpack(eData)
            (cid, sixty5, eight, rid, a1, a2, a3, b1, b2, b3, c1, c2, c3) = out
            dataIn = [cid, rid, a1, a2, a3, b1, b2, b3, c1, c2, c3]
            self.binary_debug.write('  CORD2S=%s\n' % str(out))
            coord = CORD2S(dataIn)
            self.add_coord(coord, allowOverwrites=True)
            n += 52
        self._increase_card_count('CORD2S', nentries)
        return n

    def _readCord3G(self, data, n):
        """
        (14301,143,651) - the marker for Record 7
        .. todo:: isnt this a CORD3G, not a CORD3R ???
        """
        s = Struct(self._endian + b'4i')
        nentries = (len(data) - n) // 16
        for i in range(nentries):
            eData = data[n:n + 16]  # 4*4
            out = s.unpack(eData)
            (cid, n1, n2, n3) = out
            coord = CORD3G(None, out)
            self.binary_debug.write('  CORD3G=%s\n' % str(out))
            self.add_coord(coord, allowOverwrites=True)
            n += 16
        self._increase_card_count('CORD3G', nentries)
        return n

    def _readGrid(self, data, n):  # 21.8 sec, 18.9
        """(4501,45,1) - the marker for Record 17"""
        s = Struct(self._endian + b'ii3f3i')
        ntotal = 32
        nentries = (len(data) - n) // ntotal
        self._increase_card_count('GRID', nentries)
        for i in range(nentries):
            edata = data[n:n + 32]
            out = s.unpack(edata)
            (nID, cp, x1, x2, x3, cd, ps, seid) = out
            self.binary_debug.write('  GRID=%s\n' % str(out))
            if cd >= 0 and nID < 10000000:
                node = GRID(None, out)
                self.add_node(node)
            else:
                self.log.debug("*nID=%s cp=%s x1=%-5.2f x2=%-5.2f x3=%-5.2f cd=%-2s ps=%s seid=%s" % (nID, cp, x1, x2, x3, cd, ps, seid))
            n += ntotal
        return n

    def _readSEQGP(self, data, n):
        """(5301,53,4) - the marker for Record 27"""
        self.log.debug('skipping SEQGP in GEOM1\n')
        return len(data)
