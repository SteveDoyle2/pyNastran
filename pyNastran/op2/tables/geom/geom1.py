"""
defines readers for BDF objects in the OP2 GEOM1/GEOM1S table
"""
#pylint: disable=C0301,C0103,W0612,R0914,C0326
from __future__ import annotations
from struct import Struct
from collections import defaultdict
from typing import Union, TYPE_CHECKING

import numpy as np

from pyNastran.bdf.cards.nodes import GRID, POINT, SEQGP
#from pyNastran.bdf.cards.parametric.geometry import FEFACE

from pyNastran.bdf.cards.coordinate_systems import (
    CORD1R, CORD1C, CORD1S,
    CORD2R, CORD2C, CORD2S,
    CORD3G)
from pyNastran.bdf.cards.elements.damper import CVISC
#from pyNastran.bdf.cards.elements.mass import CMASS2
from pyNastran.op2.op2_interface.op2_reader import mapfmt, reshape_bytes_block
from .utils import get_minus1_start_end
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2_geom import OP2Geom

class GEOM1:
    """defines methods for reading op2 nodes/coords"""

    def _read_fake(self, data: bytes, n: int) -> int:
        return self.op2._read_fake(data, n)

    def read_geom1_4(self, data: bytes, ndata: int):
        return self.op2._read_geom_4(self.geom1_map, data, ndata)

    def __init__(self, op2: OP2Geom):
        self.op2 = op2
        geom2 = self.op2.reader_geom2
        self.geom1_map = {
            (1701, 17, 6): ['CORD1C', self._read_cord1c],    # record 1
            (1801, 18, 5): ['CORD1R', self._read_cord1r],    # record 2
            (1901, 19, 7): ['CORD1S', self._read_cord1s],    # record 3

            (2001, 20, 9): ['CORD2C', self._read_cord2c],    # record 4
            (2101, 21, 8): ['CORD2R', self._read_cord2r],    # record 5
            (2201, 22, 10): ['CORD2S', self._read_cord2s],   # record 6

            (14301,143,651): ['CORD3G', self._read_cord3g],  # record 7

            (4501,  45,  1): ['GRID',   self._read_grid],    # record 17
            (5301,  53,  4): ['SEQGP',  self._read_seqgp],   # record 27

            (2301,  23, 304): ['CSUPER',  self._read_fake],  # record 8
            (5501,  55, 297): ['CSUPEXT', self._read_fake],  # record 9
            (1627,  16, 463): ['EXTRN',   self._read_extrn],  # record 10
            (6101,  61, 388): ['FEEDGE',  self._read_feedge],  # record 11
            (6601,  66, 392): ['GMCURVE', self._read_gmcurv],  # record 12
            (6201,  62, 389): ['FEFACE',  self._read_feface],  # record 13
            (6001,  60, 377): ['POINT',   self._read_point],  # record 14
            (10101,101, 394): ['GMSURF',  self._read_gmsurf],  # record 15
            (6401,  64, 402): ['GMCORD',  self._read_gmcord],  # record 16
            # 17 - GRID  (above)
            (1527, 15, 466): ['SEBNDRY', self._read_fake],  # record 18
            (1427, 14, 465): ['SEBULK',  self._read_sebulk],  # record 19 - superelements/see103q4.op2
            (427,   4, 453): ['SECONCT', self._read_seconct],  # record 20

            (7902, 79, 302): ['SEELT',   self._read_seelt],  # record 21
            (527,  72, 454): ['SEEXCLD', self._read_fake],  # record 22
            (1027, 10, 459): ['SELABEL', self._read_selabel],  # record 23 - superelements/see103q4.op2
            (827,   8, 457): ['SELOC',   self._read_seloc],  # record 24 - superelements/see103q4.op2
            (927,   9, 458): ['SEMPLN',  self._read_sempln],  # record 25 - superelements/see103q4.op2
            (1327, 13, 464): ['SENQSET', self._read_fake],  # record 26
            # 27 - SEQGP (above)
            (5401, 54, 305): ['SEQSEP',  self._read_fake],  # record 28
            (5601, 56, 296): ['SESET',   self._read_seset],  # record 29
            (1227, 12, 462): ['SETREE',  self._read_fake],  # record 30
            (5678, 71, 475): ['SNORM',   self._read_snorm],  # record 31
            (5701, 57, 323): ['CSUPER1', self._read_fake],  # record 32

            (5801,   58, 324): ['SUPUP', self._read_fake],  # record 33 - CSUPUP in NX; SUPUP in MSC
            (14101, 141, 403): ['SWLDPRM', self._read_fake],  # record 34

            (1101,   11,  66): ['CMASS2', geom2._read_cmass2],  # record
            (3901,   39,  50): ['CVISC', self._read_cvisc],  # record
            (13301, 133, 509): ['', self._read_fake],  # record
            (1127,   11, 461) : ['SELOAD', self._read_fake],  # record NX
            (4501, 45, 1120001): ['GRID', self._read_grid_maybe],  # record ???; test_ibulk
            (4501, 45, 810001): ['GRID', self._read_grid],


            #F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_consolid31.op2
            (2001, 20, 2220009): ['CORD2C?', self._read_cord2cx],

            # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_boltsold01d.op2
            (2101, 21, 2220008) : ['CORD2R?', self._read_cord2rx],
            (2001, 20, 1310009) : ['CORD2C-NX', self._read_cord2c_nx],

            (2101, 21, 1310008) : ['CORD2R?', self._read_fake],
            (501, 5, 43) : ['CORDx?', self._read_cord3g],
            (6591, 65, 677) : ['ATVBULK', self._read_fake],
            (2201, 22, 2220010) : ['CORDx?', self._read_fake],
            (1209, 96, 665) : ['IMPERF', self._read_fake],

            # nx
            #(707, 7, 124) :  ['EPOINT', self._read_epoint],  # record 12
        }

    #def _read_fake_c(self, data: bytes, n: int) -> int:
        #"""(2101, 21, 2220008)

          #ints    = (20, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          #1.875, 0, 0, 0, 0, 0,
          #1.75, 0, -1.75,
          #50, 1, 2, 0, 0,
          #2.565, 0, 2.565, 0,
          #3.390625, 0,
          #2.5625, 0, 2.5625, 0,
          #1079590912, 0, 1076232192, 0, 1076101120, 0, 1079574528)
          #floats  = (20, 1, 2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          #1.875, 0.0, 0.0, 0.0, 0.0, 0.0,
          #1.75, 0.0, -1.75,
          #50, 1, 2, 0.0,
          #0.0, 2.5625, 0.0,
          #2.5625, 0.0, 3.390625,
          #0.0, 2.5625, 0.0,
          #2.5625, 0.0, 3.39453125,
          #0.0, 2.59375, 0.0,
          #2.5625, 0.0, 3.390625)
        #"""
        #self.show_data(data[n:])
        #return aaa

    #def _read_new(self, data: bytes, n: int) -> int:
        #"""
        #C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\e402conm1_02.op2
        #ndata = 104:
            #ints    = (64, 0,   1, 0,   2, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 1072693248, -2409647, -1126432769, 0, 0, -2409647, -1126432769, 0, -1074790400, 46723, 1022656512)
            #floats  = (64, 0.0, 1, 0.0, 2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.875, nan, -0.02685546688735485, 0.0, 0.0, nan, -0.02685546688735485, 0.0, -1.875, 6.547286814864843e-41, 0.02984619140625)
            #doubles (float64) = (64, 1, 2, 0.0, 0.0, 0.0, 0.0, 1.0, -1.554312234e-15, 0.0, -1.554312234e-15, -1.0, 4.551914401e-15)
            #long long (int64) = (64, 1, 2, 0,   0,   0,   0,   1.0, -4837991899705164975, 0, -4837991899705164975, -4616189618054758400, 4392276274081478275)
        #"""
        #self.show_data(data[n:], 'qsdfi')

    def _read_seelt(self, data: bytes, n: int) -> int:
        """
        Record – SEELT(7902,79,302)
        Word Name Type Description
        1 SEID I Superelement identification number
        2 EID  I Element identification number
        Word 2 repeats until End of Record

        data = (
            5, 24, -1,
            66, 662001, 662003, 662007, 662019, -1,
            67, 672001, 672008, 672015, -1,
            68, 682001, 682019, -1)

        """
        op2 = self.op2
        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        istart, iend = get_minus1_start_end(ints)
        ncards = 0

        size = op2.size
        #di = 0
        for (i0, i1) in zip(istart, iend):
            assert ints[i1] == -1, ints[i1]
            seid, *eids = ints[i0:i1]

            seelt = op2.add_seelt(seid, eids)
            str(seelt)
            n += (i1 - i0 + 1) * size
            ncards += 1
            op2.card_count['SEELT'] = ncards
        return n

    def _read_seset(self, data: bytes, n: int) -> int:
        """
        Record 30 -- SESET(5601,56,296)

        Word Name Type Description
        1 SEID I Superelement identification number
        2 G    I Grid or scalar point identification number
        Word 2 repeats until End of Record

        data = (
            1, 33, 34, 37, 38, -1,
            1, 93, -98, -1,
            7, 1, -8, -1)

        SESET SEID G1 G2 G3 G4 G5 G6 G7

        This card is straight nonsense...
        When you have nodes [-8, 1], it means:
          1, THRU, 8

        The subsequent case can be [-50, -8, 1, 45].
        Well, we saw the 1,THRU,8 already, so we filter that out and:
          [-50, 45]

        is leftover, so we turn that into:
          1, THRU, 8
          45, THRU, 50

        """
        op2 = self.op2
        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        istart, iend = get_minus1_start_end(ints)

        ncards = 0
        size = op2.size
        set_ids_helper = defaultdict(set)
        set_ids = defaultdict(set)
        di = 0
        for (i0, i1) in zip(istart, iend):
            assert ints[i1] == -1, ints[i1]
            seid, *nids = ints[i0:i1]
            #print('*', seid, nids)
            if min(nids) > 0:
                set_ids_helper[seid].update(set(nids))
                set_ids[seid].update(set(nids))
            #if seid not in set_ids_helper and min(nids) > 0:
                #set_ids_helper[seid].update(set(nids))
                #set_ids[seid].update(set(nids))
                #continue
            else:
                #print(seid, nids, set_ids_helper[seid])
                diff = list(set(nids).difference(set_ids_helper[seid]))
                diff.sort()
                assert len(diff) == 2, diff

                # -10, 1  -> 1,THRU,10
                max_value, min_value = diff
                nids_new = list(range(min_value, -max_value+1))
                assert min(nids_new) > 0, nids_new
                set_ids[seid].update(set(nids_new))

                #print(seid, set_ids[seid])
                #print(seid, nids, diff, nids_new)
                set_ids_helper[seid].update(set(nids))
            n += (i1 - i0 + 1) * size
            di += (i1 - i0 + 1)
            assert i1 +1 == di, f'di={di} i1+1={i1+1} nints={len(ints)}'
        assert i1+1 == len(ints), f'i1+1={i1+1} nints={len(ints)}'

        ncards += len(set_ids)
        for seid, nids in sorted(set_ids.items()):
            nids = list(nids)
            nids.sort()
            assert min(nids) > 0, nids
            seset = op2.add_seset(seid, nids)
            seset.write_card_16()
            #print(seset)
        op2.card_count['SESET'] = ncards
        assert n == len(data), f'factor={op2.factor} size={size} n={n} ndata={len(data)}'
        return n

    def _read_snorm(self, data: bytes, n: int) -> int:
        """
        Record – SNORM(5678,71,475)
        Word Name Type Description
        1 GID I Grid point identification number
        2 CID I Coordinate system identification number
        3 N1 RS Normal component in direction 1 of CID
        4 N2 RS Normal component in direction 2 of CID
        5 N3 RS Normal component in direction 3 of CID

        ints    = (-1059, 101000001, 0, 0, 0)
        floats  = (nan, 2.5040420265274087e-35, 0.0, 0.0, 0.0)

        """
        op2 = self.op2
        structi = Struct(op2._endian + b'2i 3f')
        ntotal = 20 * op2.factor
        ndatai = len(data) - n
        nentries = ndatai // ntotal
        assert nentries > 0, nentries
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            (nid, cid, n1, n2, n3) = out
            normal = [n1, n2, n3]
            n += ntotal
            if nid < 0:
                op2.log.warning(f'geom skipping SNORM nid={nid} cid={cid} normal={normal}')
                continue
            snorm = op2.add_snorm(nid, normal, cid=cid)
            snorm.write_card_16()
        return n

    def _read_seconct(self, data: bytes, n: int) -> int:
        """
        Record – SECONCT(427,4,453)

        Word Name Type Description
        1 SEIDA I Superelement A identification number
        2 SEIDB I Superelement B identification number
        3 TOL  RS Location tolerance
        4 LOC   I Coincident location check option: yes=1 or no=2
        5 UNDEF(4) None
        9  GA I Grid point identification number in SEIDA
        10 GB I Grid point identification number in SEIDB
        Words 9 through 10 repeat until (-1,-1) occurs

        """
        op2 = self.op2
        #n0 = n
        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        floats = np.frombuffer(data[n:], op2.fdtype8).copy()
        iminus1 = np.where(ints == -1)[0]
        iminus1_start = iminus1[::2]
        iminus1_end = iminus1[1::2]

        ncards = 0
        istart = [0] + list(iminus1_end + 1)
        iend = iminus1_start
        size = op2.size
        for (i0, i1) in zip(istart, iend):
            assert ints[i1] == -1, ints[i1]
            seid_a, seid_b = ints[i0:i0+2]
            tol = floats[i0+2]
            loc_int = ints[i0+3]
            gab = ints[i0+4:i1]
            nrows = len(gab) // 2
            gab = gab.reshape(nrows, 2)
            nodes_a = gab[:, 0]
            nodes_b = gab[:, 1]

            # yes=1 or no=2
            if loc_int == 1:
                loc = 'YES'
            elif loc_int == 2:
                loc = 'NO'
            else:
                raise NotImplementedError(loc_int)

            seconct = op2.add_seconct(seid_a, seid_b, tol, loc,
                                      nodes_a, nodes_b)
            #print(seconct)
            str(seconct)
            n += (i1 - i0 + 2) * size
            ncards += 1
        op2.card_count['SECONCT'] = ncards
        return n

    def _read_cord2c_nx(self, data: bytes, n: int) -> int:
        """
        doubles (float64) = (6, 2, 2, 0.0, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 1.0])
        long long (int64) = (6, 2, 2, 0, 0, 0, 0, 0, 0, 4607182418800017408, 4607182418800017408, 0, 4607182418800017408)
        """
        #assert len(d)
        #self.op2.show_data(data[n:], types='ifsqd', endian=None, force=False)
        #raise RuntimeError('is this a cord2s?')
        op2 = self.op2
        structi = Struct(op2._endian + b'qqqq 9d')
        ntotal = 52 * op2.factor
        ndatai = len(data) - n
        nentries = ndatai // ntotal
        assert nentries > 0, nentries
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            (cid, one, two, rid, a1, a2, a3, b1, b2, b3, c1, c2, c3) = out
            assert (one, two) == (2, 2) # CORD-2-C
            origin = [a1, a2, a3]
            zaxis = [b1, b2, b3]
            xzplane = [c1, c2, c3]
            coord = op2.add_cord2c(cid, origin, zaxis, xzplane)
            str(coord)
            n += ntotal
        op2.to_nx('; because CORD2C-NX was found')
        return n

    def _read_cord2cx(self, data: bytes, n: int) -> int:
        """
        I think this is a ROTOR coordinate system...

        data = (2001, 20, 2220009,
               100002, 2, 2, 100001,
                   0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0,
                   0.0, 1.875, 0.0,
                   0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0,
                   0.0, 0.0, 1.875)
        """
        op2 = self.op2
        n = self._read_cordx(data, n, cord_type=2, cord_n=2)
        return n

        #self.show_data(data, types='ifs')
        ndatai = len(data) - n
        s = Struct(op2._endian + b'4i 9d')
        # CORD2C    100002  100001     0.0     0.0     0.0     1.0     0.0     0.0+
        # +            0.0     0.0     1.0
        # 100002 2 2 100001
        # [0.0, 0.0, 0.0] [1.0, 0.0, 0.0] [0.0, 0.0, 1.0]
        assert ndatai == 88, ndatai
        #print(len(data[n:]))
        out = s.unpack(data[n:])
        (cid, two_a, two_b, rid,
         a1, a2, a3, b1, b2, b3, c1, c2, c3) = out
        origin = [a1, a2, a3]
        zaxis = [b1, b2, b3]
        xzplane = [c1, c2, c3]
        assert (two_a, two_b) == (2, 2), (two_a, two_b)
        #print(cid, rid)
        #print(a, b, c)
        coord = op2.add_cord2c(cid, origin, zaxis, xzplane, rid=rid,
                               setup=True, comment='')
        print(coord)
        return len(data)

    def _read_cord2rx(self, data: bytes, n: int) -> int:
        n = self._read_cordx(data, n, cord_type=1, cord_n=2)
        return n

    def _read_cordx(self, data: bytes, n: int, cord_type: int, cord_n: int) -> int:
        """
        (2101, 21, 2220008)
        CORD2R  4               0.      0.      0.      0.      -1.     0.      +
        +       0.      0.      -1.
        $  Integrated Coordinate System (spring)
        CORD2R  5               0.      0.      0.      0.      -1.     0.      +
        +       0.      0.      -1.
        $  Integrated Coordinate System (dashpot)
        CORD2R  6               0.      0.      0.      0.      -1.     0.      +
        +       0.      0.      -1.

                   ?  ?  ?
                   i  i  i  [i/f]*10                         f      [i/f] * 7               f
        ints    = (4, 1, 2, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], -1.875, [0, 0, 0, 0, 0, 0, 0], -1.875,
                   5, 1, 2, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], -1.875, [0, 0, 0, 0, 0, 0, 0], -1.875,
                   6, 1, 2, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], -1.875, [0, 0, 0, 0, 0, 0, 0], -1.875)
        floats  = (4, 1, 2, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], -1.875, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.875,
                   5, 1, 2, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], -1.875, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.875,
                   6, 1, 2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.875, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.875)

        """
        op2 = self.op2
        ntotal = 88 * op2.factor  # 22*4
        #structi = Struct(mapfmt(op2._endian + b'iiq 9d', op2.size))
        structi = Struct(mapfmt(op2._endian + b'iiii 9d', op2.size))
        ndatai = len(data) - n
        nentries = ndatai // ntotal
        assert nentries > 0
        assert ndatai % ntotal == 0, ndatai
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            (cid, one, two, rid, a1, a2, a3, b1, b2, b3, c1, c2, c3) = out
            assert cid > 0, out
            assert one == cord_type, (one, out)
            assert two == cord_n, (two, out)
            origin = [a1, a2, a3]
            zaxis = [b1, b2, b3]
            xzplane = [c1, c2, c3]
            if (two, one) == (2, 1):
                coord = op2.add_cord2r(
                    cid, origin, zaxis, xzplane, rid=rid,
                    setup=True, comment='')
            elif (two, one) == (2, 2):
                coord = op2.add_cord2c(
                    cid, origin, zaxis, xzplane, rid=rid,
                    setup=True, comment='')
            else:
                raise RuntimeError((two, one))
            print(coord)
            n += ntotal
        return n

    def _read_cord1c(self, data: bytes, n: int) -> int:
        """
        (1701,17,6) - the marker for Record 1
        """
        op2 = self.op2
        ntotal = 24 * op2.factor  # 6*4
        struct_6i = Struct(mapfmt(op2._endian + b'6i', op2.size))
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = struct_6i.unpack(edata)
            (cid, one, two, g1, g2, g3) = out
            assert one in [1, 2], one
            assert two in [1, 2], two
            if op2.is_debug_file:
                op2.binary_debug.write('  CORD1C=%s\n' % str(out))
            data_in = [cid, g1, g2, g3]
            coord = CORD1C.add_op2_data(data_in)
            op2._add_methods._add_coord_object(coord)
            n += ntotal
        op2.increase_card_count('CORD1C', nentries)
        return n

    def _read_cord1r(self, data: bytes, n: int) -> int:
        """
        (1801,18,5) - the marker for Record 2
        """
        op2 = self.op2
        ntotal = 24 * op2.factor  # 6*4
        struct_6i = Struct(mapfmt(op2._endian + b'6i', op2.size))
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = struct_6i.unpack(edata)
            (cid, one1, one2, g1, g2, g3) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  CORD1R=%s\n' % str(out))
            assert one1 == 1, one1
            assert one2 == 1, one2
            data_in = [cid, g1, g2, g3]
            coord = CORD1R.add_op2_data(data_in)
            op2._add_methods._add_coord_object(coord)
            n += ntotal
        op2.increase_card_count('CORD1R', nentries)
        return n

    def _read_cord1s(self, data: bytes, n: int) -> int:
        """
        (1901,19,7) - the marker for Record 3
        """
        op2 = self.op2
        ntotal = 24 * op2.factor  # 6*4
        struct_6i = Struct(mapfmt(op2._endian + b'6i', op2.size))
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = struct_6i.unpack(edata)
            (cid, three, one, g1, g2, g3) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  CORD1S=%s\n' % str(out))
            assert three == 3, three
            assert one == 1, one
            data_in = [cid, g1, g2, g3]
            coord = CORD1S.add_op2_data(data_in)
            op2._add_methods._add_coord_object(coord, allow_overwrites=False)
            n += ntotal
        op2.increase_card_count('CORD1S', nentries)
        return n

    def _read_cord2c(self, data: bytes, n: int) -> int:
        """
        (2001,20,9) - the marker for Record 4
        """
        n = self._read_cord2x(data, n, 'CORD2C', CORD2C, (2, 2))
        return n

    def _read_cord2r(self, data: bytes, n: int) -> int:
        """
        (2101,21,8) - the marker for Record 5
        """
        n = self._read_cord2x(data, n, 'CORD2R', CORD2R, (1, 2))
        return n

    def _read_cord2s(self, data: bytes, n: int) -> int:
        """
        (2201,22,10) - the marker for Record 6
        """
        n = self._read_cord2x(data, n, 'CORD2S', CORD2S, (3, 2))
        return n

    def _read_cord2x(self, data: bytes, n: int, card_name: str, card_obj,
                     coord_flag: tuple[int, int]) -> int:
        op2 = self.op2
        if op2.table_name == b'GEOM1N' and op2.factor == 1:
            try:
                n2, coords = self._read_cord2x_22(data, n, card_name, card_obj, coord_flag)
            except Exception:
                n2, coords = self._read_cord2x_13(data, n, card_name, card_obj, coord_flag)
        else:
            n2, coords = self._read_cord2x_13(data, n, card_name, card_obj, coord_flag)

        ncoords = len(coords)
        assert n is not None
        for coord in coords:
            op2._add_methods._add_coord_object(coord, allow_overwrites=False)
        op2.card_count[card_name] = ncoords
        return n2

    def _read_cord2x_22(self, data: bytes, n: int,
                        coord_name: str,
                        coord_cls: Union[CORD2R, CORD2C, CORD2S],
                        flags: tuple[int, int]) -> int:
        """
        (2101,21,8) - CORD2R
        (2201,22,10) - CORD2S

        """
        op2 = self.op2
        ntotal = 88 * op2.factor # 22*4
        s = Struct(op2._endian + b'4i9d')
        ndatai = len(data) - n
        nentries = ndatai // ntotal
        assert ndatai % ntotal == 0
        coords = []
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            (cid, sixty5, eight, rid, a1, a2, a3, b1, b2, b3, c1, c2, c3) = out
            data_in = [cid, rid, a1, a2, a3, b1, b2, b3, c1, c2, c3]
            assert (sixty5, eight) == flags, f'(sixty5,eight)={(sixty5, eight)} flags={flags}'
            if op2.is_debug_file:
                op2.binary_debug.write(f'  {coord_name}={out}\n')
            coord = coord_cls.add_op2_data(data_in)
            coords.append(coord)
            n += ntotal
        return n, coords

    def _read_cord2x_13(self, data: bytes, n: int,
                        coord_name: str,
                        coord_cls: Union[CORD2R, CORD2C, CORD2S],
                        flags: tuple[int, int]) -> int:
        op2 = self.op2
        ntotal = 52 * op2.factor # 13*4
        s = Struct(mapfmt(op2._endian + b'4i9f', op2.size))
        ndatai = len(data) - n
        nentries = ndatai // ntotal
        assert ndatai % ntotal == 0
        coords = []
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            (cid, sixty5, eight, rid, a1, a2, a3, b1, b2, b3, c1, c2, c3) = out
            data_in = [cid, rid, a1, a2, a3, b1, b2, b3, c1, c2, c3]
            if op2.is_debug_file:
                op2.binary_debug.write(f'  {coord_name}={out}\n')
            coord = coord_cls.add_op2_data(data_in)
            coords.append(coord)
            n += ntotal
        return n, coords

    def _read_cord3g(self, data: bytes, n: int) -> int:
        """
        (14301,143,651) - the marker for Record 7
        .. todo:: isn't this a CORD3G, not a CORD3R ???
        """
        op2 = self.op2
        ntotal = 16 * op2.factor
        struct_4i = Struct(op2._endian + b'4i')
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]  # 4*4
            out = struct_4i.unpack(edata)
            #(cid, n1, n2, n3) = out
            coord = CORD3G.add_op2_data(out)
            if op2.is_debug_file:
                op2.binary_debug.write('  CORD3G=%s\n' % str(out))
            op2._add_methods._add_coord_object(coord, allow_overwrites=False)
            n += ntotal
        op2.increase_card_count('CORD3G', nentries)
        return n

    def _read_grid_maybe(self, data: bytes, n: int) -> int:  # pragma: no cover
        """
        (4501, 45, 1120001) - the marker for Record 17
        this is a GRID card with double vales for xyz
        """
        op2 = self.op2
        ntotal = 44 * op2.factor
        structi = Struct(mapfmt(op2._endian + b'2i 3d 3i', op2.size))
        nentries = (len(data) - n) // ntotal
        leftover = (len(data) - n) % ntotal
        assert leftover == 0, f'ndata={len(data)-n} leftover={leftover}'
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            (nid, cp, x1, x2, x3, cd, ps, seid) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  GRID=%s\n' % str(out))

            # cd can be < 0
            if ps == 0:
                ps = ''
            node = GRID(nid, np.array([x1, x2, x3]), cp, cd, ps, seid)
            op2._type_to_id_map['GRID'].append(nid)
            op2.nodes[nid] = node
            #op2.log.debug(f'  nid={nid} cp={cp} x=[{x1:g}, {x2:g}, {x3:g}] cd={cd} ps={ps} seid={seid}')

            n += ntotal
        op2.increase_card_count('GRID', nentries)
        return n

    def _read_cord3g(self, data: bytes, n: int) -> int:
        """
        Record – CORD3G(501,5,43)
        Word Name Type Description
        1 CID           I Coordinate system identification number
        2 METHOD(2) CHAR4 Methods
        4 FORM(2)   CHAR4 Forms
        6 THETAID(3)    I Identification number for DEQATN or TABLE
        9 CIDREF        I Coordinate system identification number
        """
        op2 = self.op2
        assert op2.size == 4, op2.size
        structi = Struct(op2._endian + b'i 8s 8s 4i')
        ntotal = 36 * op2.factor
        ndatai = len(data) - n
        nentries = ndatai // ntotal
        assert nentries > 0, nentries
        assert ndatai % ntotal == 0, f'ndatai={ndatai} ntotal={ntotal} leftover={ndatai % ntotal}'
        #grids = {}
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            (cid, method_bytes, form_bytes, theta1, theta2, theta3, cid_ref) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  CORD3G=%s\n' % str(out))
            method = method_bytes.decode('ascii').strip()
            method_es = method[0]
            method_int = int(method[1:])

            form = form_bytes.decode('ascii').strip()

            thetas = [theta1, theta2, theta3]
            coord = op2.add_cord3g(cid,
                                   method_es, method_int, form,
                                   thetas, cid_ref, comment='')
            str(coord)
            n += ntotal
        return n

    def _read_grid(self, data: bytes, n: int) -> int:
        """(4501,45,1) - the marker for Record 17"""
        op2 = self.op2
        if op2.table_name == b'GEOM1N' and op2.factor == 1:
            try:
                n, grids = self._read_grid_11(data, n)
            except Exception:
                n, grids = self._read_grid_8(data, n)
        else:
            n, grids = self._read_grid_8(data, n)

        ngrids = len(grids)
        assert n is not None
        for nid, grid in grids.items():
            op2.nodes[nid] = grid
            op2._type_to_id_map['GRID'].append(nid)

        op2.card_count['GRID'] = ngrids
        return n

    def _read_grid_8(self, data: bytes, n: int) -> tuple[int, dict[int, GRID]]:  # 21.8 sec, 18.9
        """(4501,45,1) - the marker for Record 17"""
        op2 = self.op2
        structi = Struct(mapfmt(op2._endian + b'ii 3f 3i', op2.size))
        ntotal = 32 * op2.factor
        ndatai = len(data) - n
        nentries = ndatai // ntotal
        assert nentries > 0, nentries
        assert ndatai % ntotal == 0, f'ndatai={ndatai} ntotal={ntotal} leftover={ndatai % ntotal}'
        grids = {}
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            (nid, cp, x1, x2, x3, cd, ps, seid) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  GRID=%s\n' % str(out))
            #if nid < 10000000:
            # cd can be < 0
            if ps == 0:
                ps = ''
            node = GRID(nid, np.array([x1, x2, x3]), cp, cd, ps, seid)
            #op2._type_to_id_map['GRID'].append(nid)
            #self.nodes[nid] = node
            grids[nid] = node
            #if nid in self.nodes:
                #self.reject_lines.append(str(node))
            #else:
            #self.nodes[nid] = node
            #self.add_node(node)
            #else:
                #op2.log.warning('*nid=%s cp=%s x1=%-5.2f x2=%-5.2f x3=%-5.2f cd=%-2s ps=%s '
                                 #'seid=%s' % (nid, cp, x1, x2, x3, cd, ps, seid))
                #node = GRID(nid, np.array([x1, x2, x3]), cp, cd, ps, seid)
                #self.rejects.append(str(node))
                #nfailed += 1
            n += ntotal
        #op2.increase_card_count('GRID', nentries - nfailed)
        return n, grids

    def _read_grid_11(self, data: bytes, n: int) -> tuple[int, dict[int, GRID]]:  # 21.8 sec, 18.9
        """(4501,45,1) - the marker for Record 17"""
        op2 = self.op2
        ntotal = 44
        structi = Struct(op2._endian + b'ii 3d 3i')

        ndatai = len(data) - n
        nentries = ndatai // ntotal
        #print(f'len(data)={len(data)} ndatai={ndatai} ntotal={ntotal} nentries={nentries}')
        assert ndatai % ntotal == 0, f'len(data)={len(data)} ndatai={ndatai} ntotal={ntotal} nentries={nentries}'
        assert nentries > 0, f'len(data)={len(data)} ndatai={ndatai} ntotal={ntotal} nentries={nentries}'
        nfailed = 0
        grids = {}
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            #cp, x1, x2, x3, cd, ps, seid
            (nid, cp, x1, x2, x3, cd, ps, seid) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  GRID=%s\n' % str(out))
            #print(f'nid={nid}, cp={cp} x=({x1}, {x2}, {x3}), cd={cd} ps={ps}, seid={seid}')
            assert cd >= 0, f'nid={nid}, cp={cp} x=({x1}, {x2}, {x3}), cd={cd} ps={ps}, seid={seid}'
            assert seid == 0, f'nid={nid}, cp={cp} x=({x1}, {x2}, {x3}), cd={cd} ps={ps}, seid={seid}'
            if nid < 10000000:
                # cd can be < 0
                if ps == 0:
                    ps = ''
                node = GRID(nid, np.array([x1, x2, x3]), cp, cd, ps, seid)
                #print(node)
                #op2._type_to_id_map['GRID'].append(nid)
                #self.nodes[nid] = node
                grids[nid] = node
                #if nid in self.nodes:
                    #self.reject_lines.append(str(node))
                #else:
                #self.nodes[nid] = node
                #self.add_node(node)
            else:
                #op2.log.warning('*nid=%s cp=%s x1=%-5.2f x2=%-5.2f x3=%-5.2f cd=%-2s ps=%s '
                                 #'seid=%s' % (nid, cp, x1, x2, x3, cd, ps, seid))
                #node = GRID(nid, np.array([x1, x2, x3]), cp, cd, ps, seid)
                #self.rejects.append(str(node))
                nfailed += 1
            n += ntotal
        return n, grids

    def _read_seqgp(self, data: bytes, n: int) -> int:
        """(5301,53,4) - the marker for Record 27"""
        op2 = self.op2
        ntotal = 8 * op2.factor
        struct_2i = Struct(op2._endian + b'2i')
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]  # 2*4
            out = struct_2i.unpack(edata)
            # (nid, seid) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  SEQGP=%s\n' % str(out))
            seqgp = SEQGP.add_op2_data(out)
            op2._add_methods._add_seqgp_object(seqgp)
            n += ntotal
        op2.increase_card_count('SEQGP', nentries)
        return n

    def _read_point(self, data: bytes, n: int) -> int:
        """
        POINT(6001,60,377)
        """
        op2 = self.op2
        s = Struct(op2._endian + b'2i3f')
        ntotal = 20 * op2.factor
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]  # 5*4
            out = s.unpack(edata)
            # (nid, cid, x1, x2, x3) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  POINT=%s\n' % str(out))
            point = POINT.add_op2_data(out)
            op2._add_methods._add_point_object(point)
            n += ntotal
        op2.increase_card_count('POINT', nentries)
        return n

    #def _read_cmass2(self, data: bytes, n: int) -> int:
        #struct_i4fi = Struct(op2._endian + b'if4i')
        #nentries = (len(data) - n) // 24
        #for unused_i in range(nentries):
            #edata = data[n:n + 24]  # 6*4
            #out = struct_i4fi.unpack(edata)
            ## (eid, mass, g1, g2, c1, c2) = out
            #if op2.is_debug_file:
                #op2.binary_debug.write('  CMASS2=%s\n' % str(out))
            #op2.log.debug('  CMASS2=%s\n' % str(out))
            #element = CMASS2.add_op2_data(out)
            #print(element)
            #self.add_op2_mass(element)
            #n += 24
        #print(self.elements)
        #op2.increase_card_count('CMASS2', nentries)
        #return n
        #return len(data)

    def _read_cvisc(self, data: bytes, n: int) -> int:
        """CVISC(3901,39,50) - the marker for Record 105"""
        op2 = self.op2
        struct_4i = Struct(op2._endian + b'4i')
        ntotal = 16 * op2.factor  # 4*4
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = struct_4i.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CVISC=%s\n' % str(out))
            # (eid, pid, n1, n2) = out
            element = CVISC.add_op2_data(out)
            op2.add_op2_element(element)
            n += ntotal
        op2.increase_card_count('CVISC', nentries)
        return n


    def _read_extrn(self, data: bytes, n: int) -> int:
        """
        Record - EXTRN(1627,16,463)
        Word Name Type Description
        1 GID I Grid point identification numbers to connect external SE
        2 C   I Component numbers
        Words 1 through 2 repeat until (-1,-1) occurs

        (1627, 16, 463,
         1, 123456,
         2, 123456,
         3, 123456,
         4, 123456,
         100001, 1,
         100002, 1,
         100003, 1,
         100004, 1,
        -1, -1)
        """
        op2 = self.op2
        #Partitioned External Superelement Connection
        #Defines a boundary connection for an external superelement.
        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        iminus1 = np.where(ints == -1)[0]
        iminus1_start = iminus1[::2]
        iminus1_end = iminus1[1::2]

        #ncards = 0
        istart = [0] + list(iminus1_end + 1)
        iend = iminus1_start
        #size = op2.size
        for (i0, i1) in zip(istart, iend):
            assert ints[i1] == -1, ints[i1]
            nids = ints[i0:i1:2].tolist()
            comps = ints[i0+1:i1:2].tolist()
            #print(nids)
            #print(comps)
            for c in comps:
                assert c in [1, 2, 3, 4, 5, 6, 123, 123456], f'c={c}; nids={nids} comps={comps}'
            op2.add_aset(nids, comps)
            #self.add_extrn(nids, comps)
            assert len(nids) == len(comps)
        return len(data)

    def _read_feedge(self, data: bytes, n: int) -> int:
        """
        (2901, 29, 9601)

        Word Name Type Description
        1 EDGEID     I Edge identification number
        2 GRID1      I Identification number of end GRID 1
        3 GRID2      I Identification number of end GRID 2
        4 CID        I Coordinate system identification number
        5 GEOMIN CHAR4 Type of referencing entry: "GMCURV" or "POINT"
        6 GEOMID1    I Identification number of a POINT or GMCURV entry
        7 GEOMID2    I Identification number of a POINT or GMCURV entry
        """
        op2 = self.op2
        # C:\NASA\m4\formats\git\examples\move_tpl\phsflux4.op2
        #(200000002, 3, 1002, 6, 12, 0, 0)
        # FEEDGE EDGEID GRID1 GRID2 CIDBC GEOMIN ID1 ID2
        #FEEDGE    1002    6     12
        #self.show_data(data[12:])
        ntotal = 28 * op2.factor  # 7*4
        if op2.size == 4:
            s = Struct(op2._endian + b'4i 4s 2i') #expected
        else:
            s = Struct(op2._endian + b'4q 8s 2q') #expected
        #s = Struct(op2._endian + b'7i')
        nelements = (len(data) - n)// ntotal
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            #print(out)
            edge_id, n1, n2, cid, geomin, geom1, geom2 = out # expected

            if op2.is_debug_file:
                op2.binary_debug.write('  FEEDGE=%s\n' % str(out))

            geomin = geomin.rstrip()
            if geomin == b'POIN':
                geomin_str = 'POINT'
            elif geomin == b'GMCU':
                geomin_str = 'GMCURV'
            else:  # pragma: no cover
                raise RuntimeError(geomin)

            if cid == -1:
                cid = None
            unused_elem = op2.add_feedge(edge_id, [n1, n2], cid, [geom1, geom2],
                                         geomin=geomin_str)
            n += ntotal
        op2.card_count['FEEDGE'] = nelements
        return n

    def _read_gmcurv(self, data: bytes, n: int) -> int:
        """
        Word Name Type Description
        1 CURVID       I Curve identification number
        2 GROUP(2) CHAR4 Group of curves/surfaces to which this curve belongs
        4 CIDIN        I Coordinate system identification number for the geometry
        5 CIDBC        I Coordinate system identification number for the constraints
        6 DATA     CHAR4 Geometry evaluator specific data
        """
        op2 = self.op2
        size = op2.size
        if size == 4:
            struct_i = op2.struct_i
            structi = Struct(b'i 8s ii')
        else:
            struct_i = op2.struct_q
            structi = Struct(b'q 16s qq')

        ntotal1 = 20 * op2.factor
        ntotal2 = 64 * op2.factor
        while n < len(data):
            datab = data[n:n+ntotal1]
            curve_id, group_bytes, cid_in, cid_bc = structi.unpack(datab)
            if size == 8:
                group_bytes = reshape_bytes_block(group_bytes)
            group = group_bytes.decode('latin1').rstrip()
            #print(curve_id, group, cid_in, cid_bc)
            assert group in ['MSCGRP0', 'MSCGRP1', 'MSCGRP2'], f'GMCURV: curve_id={curve_id} group={group!r} cid_in={cid_in} cid_bc={cid_bc}'
            n += ntotal1

            databi_bytes = data[n:n+size]
            n += size
            databi = data[n:n+size]
            datab_int, = struct_i.unpack(databi)
            n += size
            while datab_int != -1:
                databi_bytes += databi
                databi = data[n:n+size]
                datab_int, = struct_i.unpack(databi)
                n += size
            datai = databi_bytes.decode('latin1').rstrip()

            data_split = ['        %s\n' % datai[i:i+ntotal2].strip() for i in range(0, len(datai), ntotal2)]
            op2.add_gmcurv(curve_id, group, data_split, cid_in=cid_in, cid_bc=cid_bc)
            #print(datai)

        #ints = np.frombuffer(data[n:], dtype=self.idtype).copy()
        #iminus1 = np.where(ints == -1)[0].tolist()
        #i0 = 0
        #for iminus1i in iminus1:
            #curve_id = ints[i0]
            #cid_in, cid_bc = ints[i0+3:i0+5]
            #s0 = n + 4
            #s1 = s0 + 8
            #group = data[s0:s1].decode('latin1').rstrip()
            #print(curve_id, group, cid_in, cid_bc)
            #assert group in ['MSCGRP1', 'MSCGRP2'], f'GMCURV: curve_id={curve_id} group={group!r} cid_in={cid_in} cid_bc={cid_bc}'

            #s2 = s1 + 8
            #s3 = 12 + iminus1i * 4
            #datai = data[s2:s3].decode('latin1').rstrip()
            #print('datai = %r' % datai)
            #i0 = iminus1i + 1
            ## n = s3 + 4
            #n = 12+(iminus1i + 1)*4
            #print('-----------------')
        #return len(data)
        return n

    def _read_feface(self, data: bytes, n: int) -> int:
        """
        Word Name Type Description
        1 FACEID    I Face identification number
        2 GRID1     I Identification number of end GRID 1
        3 GRID2     I Identification number of end GRID 2
        4 GRID3     I Identification number of end GRID 3
        5 GRID4     I Identification number of end GRID 4
        6 CIDBC     I Coordinate system identification number for the constraints
        7 SURFID(2) I Alternate method used to specify the geometry
        """
        op2 = self.op2
        ntotal = 32 * op2.factor
        structi = Struct(mapfmt(op2._endian + b'8i', op2.size))
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            (face_id, n1, n2, n3, n4, cid, surf_id1, surf_id2) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  FEFACE=%s\n' % str(out))

            nodes = [n1, n2, n3, n4]
            surf_ids = [surf_id1, surf_id2]
            unused_feface = op2.add_feface(face_id, nodes, cid, surf_ids)
            n += ntotal
        op2.increase_card_count('FEFACE', nentries)
        return n

    def _read_gmsurf(self, data: bytes, n: int) -> int:
        op2 = self.op2
        op2.log.info('geom skipping GMSURF in GEOM1')
        return len(data)

    def _read_gmcord(self, data: bytes, n: int) -> int:
        op2 = self.op2
        op2.log.info('geom skipping GMCORD in GEOM1')
        return len(data)

    def _read_sebulk(self, data: bytes, n: int) -> int:
        """
        Record 18 -- SEBULK(1427,14,465)

        Word Name Type Description
        1 SEID    I Superelement identification number
        2 TYPE    I Superelement type
        3 RSEID   I Reference superelement identification number
        4 METHOD  I Boundary point search method: 1=automatic or 2=manual
        5 TOL    RS Location tolerance
        6 LOC     I Coincident location check option: yes=1 or no=2
        7 MEDIA   I Media format of boundary data of external SE
        8 UNIT    I FORTRAN unit number of OP2 and OP4 input of external SE
        """
        op2 = self.op2
        ntotal = 32 * op2.factor # 4*8
        nentries = (len(data) - n) // ntotal
        structi = Struct(mapfmt(op2._endian + b'4if3i', op2.size))

        superelement_type_int_to_superelement_type = {
            1 : 'PRIMARY',
            5 : 'EXTOP2',
            6 : 'MIRROR',
            7 : 'FRFOP2',
        }
        loc_int_to_loc = {
            1 : 'YES',
            2 : 'NO',
        }
        method_int_to_method = {
            1: 'AUTO',
            2: 'MANUAL',
        }
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]  # 4*8
            out = structi.unpack(edata)
            (seid, superelement_type_int, rseid, method_int, tol, loc_int, media, unit) = out
            try:
                superelement_type = superelement_type_int_to_superelement_type[superelement_type_int]
            except KeyError:  # pragma: no cover
                raise NotImplementedError(f'superelement_type={superelement_type_int} not in [PRIMARY, EXTOP2, MIRROR, FRFOP2]')

            try:
                loc = loc_int_to_loc[loc_int]
            except KeyError:  # pragma: no cover
                raise NotImplementedError(f'loc={loc_int} not in [YES, NO]')

            try:
                method = method_int_to_method[method_int]
            except KeyError:  # pragma: no cover
                raise NotImplementedError(f'method={method_int} not in [AUTO, MANUAL]')

            if op2.is_debug_file:
                op2.binary_debug.write('  SEBULK=%s\n' % str(out))
            #media,
            sebulk = op2.add_sebulk(seid, superelement_type, rseid,
                                    method=method, tol=tol, loc=loc, unitno=unit)
            sebulk.validate()
            n += 32 * op2.factor
        op2.increase_card_count('SEBULK', nentries)
        return n

    def _read_seloc(self, data: bytes, n: int) -> int:
        """
        Record 23 -- SELOC(827,8,457)

        Word Name Type Description
        1 SEID I Superelement identification number
        2 GA1  I Grid point 1 identification number in SEID
        3 GA2  I Grid point 2 identification number in SEID
        4 GA3  I Grid point 3 identification number in SEID
        5 GB1  I Grid point 1 identification number in the main Bulk Data
        6 GB2  I Grid point 2 identification number in the main Bulk Data
        7 GB3  I Grid point 3 identification number in the main Bulk Data
        """
        op2 = self.op2
        ntotal = 28 * op2.factor
        structi = Struct(op2._endian + b'7i')
        nentries = (len(data) - n) // ntotal # 4*7
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            (seid, ga1, ga2, ga3, gb1, gb2, gb3) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  SELOC=%s\n' % str(out))
            op2.add_seloc(seid, [ga1, ga2, ga3], [gb1, gb2, gb3])
            n += ntotal
        op2.increase_card_count('SELOC', nentries)
        return n

    def _read_sempln(self, data: bytes, n: int) -> int:
        """
        Record 24 -- SEMPLN(927,9,458)

        1 SEID     I Superelement identification number
        2 MIRRTYPE I Mirror type

        MIRRTYPE=1 Plane
        3 G1       I    Grid point 1 identification number in the main Bulk Data
        4 G2       I    Grid point 2 identification number in the main Bulk Data
        5 G3       I    Grid point 3 identification number in the main Bulk Data
        6 UNDEF(2) none Not Defined

        MIRRTYPE=2 Normal
        3 G    I Grid point identification number in the main Bulk Data
        4 CID  I Coordinate system identification number
        5 N1  RS Normal component in direction 1 of CID
        6 N2  RS Normal component in direction 2 of CID
        7 N3  RS Normal component in direction 3 of CID
        """
        op2 = self.op2
        struct2i = Struct(op2._endian + b'2i') # 8
        struct5i = Struct(op2._endian + b'5i') # 20
        #struct2i_3f = Struct(op2._endian + b'2i3f') # 20

        nentries = (len(data) - n) // 28 # 4*7
        for unused_i in range(nentries):
            edata1 = data[n:n + 8]  # 4*2
            edata2 = data[n+8:n + 28]  # 4*7
            out = struct2i.unpack(edata1)
            (seid, mirror_type) = out
            if mirror_type == 1:
                g1, g2, g3, unused_junk1, unused_junk2 = struct5i.unpack(edata2)
                op2.add_sempln(seid, g1, g2, g3)
            else:
                raise NotImplementedError(mirror_type)
            if op2.is_debug_file:
                op2.binary_debug.write('  SEMPLN=%s\n' % str(out))
            n += 28
        op2.increase_card_count('SEMPLN', nentries)
        return n

    def _read_selabel(self, data: bytes, n: int) -> int:
        """
        Record 22 -- SELABEL(1027,10,459)

        Word Name Type Description
        1 SEID I Superelement identification number
        2 LABEL(14) CHAR4 Label associated with superelement SEID

        """
        op2 = self.op2
        ntotal = 60 * op2.factor
        structi = Struct(op2._endian + b'i14s') # 18
        structi = Struct(op2._endian + b'i56s') # 60
        nentries = (len(data) - n) // ntotal # 4+18
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            (seid, label) = out
            label = label.decode(op2._encoding).rstrip()
            if op2.is_debug_file:
                op2.binary_debug.write('  SELABEL=%s\n' % str(out))
            selabel = op2.add_selabel(seid, label)
            selabel.validate()
            n += ntotal
        op2.increase_card_count('SELABEL', nentries)
        return n
