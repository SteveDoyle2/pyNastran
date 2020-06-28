"""
defines readers for BDF objects in the OP2 GEOM1/GEOM1S table
"""
#pylint: disable=C0301,C0103,W0612,R0914,C0326
from struct import Struct
from typing import Tuple, Union

import numpy as np

from pyNastran.bdf.cards.nodes import GRID, POINT, SEQGP
#from pyNastran.bdf.cards.parametric.geometry import FEFACE

from pyNastran.bdf.cards.coordinate_systems import (
    CORD1R, CORD1C, CORD1S,
    CORD2R, CORD2C, CORD2S,
    CORD3G)
from pyNastran.bdf.cards.elements.damper import CVISC
#from pyNastran.bdf.cards.elements.mass import CMASS2
from pyNastran.op2.tables.geom.geom_common import GeomCommon
from pyNastran.op2.op2_interface.op2_reader import mapfmt, reshape_bytes_block


class GEOM1(GeomCommon):
    """defines methods for reading op2 nodes/coords"""

    #def _add_node_object(self, node, allow_overwrites=False):
        #"""GRDSET creates duplicate nodes...what about duplicate nodes?"""
        #key = node.nid
        #assert key > 0, 'nid=%s node=%s' % (key, node)
        #if key in self.nodes:
            #fields1 = self.nodes[key].raw_fields()
            #fields2 = node.raw_fields()
            ##grid, nid, cp, x1, x2, x3, cd, ps, seid
            #for i, (v1, v2) in enumerate(zip(fields1, fields2)):
                #if v1 != v2:
                    #self.log.info('i=%s v1=%r v2=%r fields1=%s\nfields2=%s' % (
                        #i, v1, v2, fields1, fields2))
        #else:
            #self._type_to_id_map[node.type].append(key)
        #self.nodes[key] = node

    def _read_geom1_4(self, data, ndata):
        return self._read_geom_4(self._geom1_map, data, ndata)

    def __init__(self):
        GeomCommon.__init__(self)
        self._geom1_map = {
            (1701, 17, 6): ['CORD1C', self._read_cord1c],    # record 1
            (1801, 18, 5): ['CORD1R', self._read_cord1r],    # record 2
            (1901, 19, 7): ['CORD1S', self._read_cord1s],    # record 3
            (2001, 20, 9): ['CORD2C', self._read_cord2c],    # record 4

            #F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_consolid31.op2
            (2001, 20, 2220009): ['GRIDx?', self._read_fake],

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
            (427,   4, 453): ['SECONCT', self._read_fake],  # record 20

            (7902, 79, 302): ['SEELT',   self._read_fake],  # record 21
            (527,  72, 454): ['SEEXCLD', self._read_fake],  # record 22
            (1027, 10, 459): ['SELABEL', self._read_selabel],  # record 23 - superelements/see103q4.op2
            (827,   8, 457): ['SELOC',   self._read_seloc],  # record 24 - superelements/see103q4.op2
            (927,   9, 458): ['SEMPLN',  self._read_sempln],  # record 25 - superelements/see103q4.op2
            (1327, 13, 464): ['SENQSET', self._read_fake],  # record 26
            # 27 - SEQGP (above)
            (5401, 54, 305): ['SEQSEP',  self._read_fake],  # record 28
            (5601, 56, 296): ['SESET',   self._read_fake],  # record 29
            (1227, 12, 462): ['SETREE',  self._read_fake],  # record 30
            (5678, 71, 475): ['SNORM',   self._read_fake],  # record 31
            (5701, 57, 323): ['CSUPER1', self._read_fake],  # record 32

            (5801,   58, 324): ['SUPUP', self._read_fake],  # record 33 - CSUPUP in NX; SUPUP in MSC
            (14101, 141, 403): ['SWLDPRM', self._read_fake],  # record 34

            (1101,   11,  66): ['CMASS2', self._read_cmass2],  # record
            (3901,   39,  50): ['CVISC', self._read_cvisc],  # record
            (13301, 133, 509): ['', self._read_fake],  # record
            (1127,   11, 461) : ['SELOAD', self._read_fake],  # record NX
            (4501, 45, 1120001): ['GRID/BCT?/BOLT?', self._read_grid_maybe],  # record ???; test_ibulk

            # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_boltsold01d.op2
            (2101, 21, 2220008) : ['CORDx?', self._read_fake],

            # nx
            #(707, 7, 124) :  ['EPOINT', self._read_epoint],  # record 12
        }

    def _read_cord1c(self, data: bytes, n: int) -> int:
        """
        (1701,17,6) - the marker for Record 1
        """
        ntotal = 24 * self.factor  # 6*4
        struct_6i = Struct(mapfmt(self._endian + b'6i', self.size))
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = struct_6i.unpack(edata)
            (cid, one, two, g1, g2, g3) = out
            assert one in [1, 2], one
            assert two in [1, 2], two
            if self.is_debug_file:
                self.binary_debug.write('  CORD1C=%s\n' % str(out))
            data_in = [cid, g1, g2, g3]
            coord = CORD1C.add_op2_data(data_in)
            self._add_coord_object(coord)
            n += ntotal
        self.increase_card_count('CORD1C', nentries)
        return n

    def _read_cord1r(self, data: bytes, n: int) -> int:
        """
        (1801,18,5) - the marker for Record 2
        """
        ntotal = 24 * self.factor  # 6*4
        struct_6i = Struct(mapfmt(self._endian + b'6i', self.size))
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = struct_6i.unpack(edata)
            (cid, one1, one2, g1, g2, g3) = out
            if self.is_debug_file:
                self.binary_debug.write('  CORD1R=%s\n' % str(out))
            assert one1 == 1, one1
            assert one2 == 1, one2
            data_in = [cid, g1, g2, g3]
            coord = CORD1R.add_op2_data(data_in)
            self._add_coord_object(coord)
            n += ntotal
        self.increase_card_count('CORD1R', nentries)
        return n

    def _read_cord1s(self, data: bytes, n: int) -> int:
        """
        (1901,19,7) - the marker for Record 3
        """
        ntotal = 24 * self.factor  # 6*4
        struct_6i = Struct(mapfmt(self._endian + b'6i', self.size))
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = struct_6i.unpack(edata)
            (cid, three, one, g1, g2, g3) = out
            if self.is_debug_file:
                self.binary_debug.write('  CORD1S=%s\n' % str(out))
            assert three == 3, three
            assert one == 1, one
            data_in = [cid, g1, g2, g3]
            coord = CORD1S.add_op2_data(data_in)
            self._add_coord_object(coord, allow_overwrites=False)
            n += ntotal
        self.increase_card_count('CORD1S', nentries)
        return n

    def _read_cord2c(self, data: bytes, n: int) -> int:
        """
        (2001,20,9) - the marker for Record 4
        """
        if self.table_name == b'GEOM1N' and self.factor == 1:
            n2 = self._read_cord2x_22(data, n, 'CORD2C', CORD2C, (2, 2))
        else:
            n2 = self._read_cord2x_13(data, n, 'CORD2C', CORD2C, (2, 2))
        return n2

    def _read_cord2r(self, data: bytes, n: int) -> int:
        """
        (2101,21,8) - the marker for Record 5
        """
        if self.table_name == b'GEOM1N' and self.factor == 1:
            n2 = self._read_cord2x_22(data, n, 'CORD2R', CORD2R, (1, 2))
        else:
            n2 = self._read_cord2x_13(data, n, 'CORD2R', CORD2R, (1, 2))
        return n2

    def _read_cord2s(self, data: bytes, n: int) -> int:
        """
        (2201,22,10) - the marker for Record 6
        """
        if self.table_name == b'GEOM1N' and self.factor == 1:
            n2 = self._read_cord2x_22(data, n, 'CORD2S', CORD2S, (3, 2))
        else:
            n2 = self._read_cord2x_13(data, n, 'CORD2S', CORD2S, (3, 2))
        return n2

    def _read_cord2x_22(self, data: bytes, n: int,
                        coord_name: str,
                        coord_cls: Union[CORD2R, CORD2C, CORD2S],
                        flags: Tuple[int, int]) -> int:
        """
        (2201,22,10) - the marker for Record 6
        """
        ntotal = 88 # 22*4
        s = Struct(self._endian + b'4i9d')
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            (cid, sixty5, eight, rid, a1, a2, a3, b1, b2, b3, c1, c2, c3) = out
            data_in = [cid, rid, a1, a2, a3, b1, b2, b3, c1, c2, c3]
            assert (sixty5, eight) == flags, f'(sixty5,eight)={(sixty5, eight)} flags={flags}'
            if self.is_debug_file:
                self.binary_debug.write(f'  {coord_name}={out}\n')
            coord = coord_cls.add_op2_data(data_in)
            self._add_coord_object(coord, allow_overwrites=False)
            n += ntotal
        self.increase_card_count(coord_name, nentries)
        return n

    def _read_cord2x_13(self, data: bytes, n: int,
                        coord_name: str,
                        coord_cls: Union[CORD2R, CORD2C, CORD2S],
                        flags: Tuple[int, int]) -> int:
        ntotal = 52 * self.factor # 13*4
        s = Struct(mapfmt(self._endian + b'4i9f', self.size))
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            (cid, sixty5, eight, rid, a1, a2, a3, b1, b2, b3, c1, c2, c3) = out
            data_in = [cid, rid, a1, a2, a3, b1, b2, b3, c1, c2, c3]
            if self.is_debug_file:
                self.binary_debug.write(f'  {coord_name}={out}\n')
            coord = coord_cls.add_op2_data(data_in)
            self._add_coord_object(coord, allow_overwrites=False)
            n += ntotal
        self.increase_card_count(coord_name, nentries)
        return n

    def _read_cord3g(self, data: bytes, n: int) -> int:
        """
        (14301,143,651) - the marker for Record 7
        .. todo:: isnt this a CORD3G, not a CORD3R ???
        """
        struct_4i = Struct(self._endian + b'4i')
        nentries = (len(data) - n) // 16
        for unused_i in range(nentries):
            edata = data[n:n + 16]  # 4*4
            out = struct_4i.unpack(edata)
            #(cid, n1, n2, n3) = out
            coord = CORD3G.add_op2_data(out)
            if self.is_debug_file:
                self.binary_debug.write('  CORD3G=%s\n' % str(out))
            self._add_coord_object(coord, allow_overwrites=False)
            n += 16
        self.increase_card_count('CORD3G', nentries)
        return n

    def _read_grid_maybe(self, data: bytes, n: int) -> int:  # pragma: no cover
        """(4501, 45, 1120001) - the marker for Record 17"""
        return len(data)
        #nfields = (len(data) - n) // 4
        # nfields = 3 * 11 * 17 * 71

        # it's not 11, 17...
        #self.show_data(data[12:], types='if')

        # 2i: correct
        #   id, 0
        # i: ???
        # i/f: ???
        # f: correct
        # i: ???
        # f: correct
        # 5i: ???
        #                                    ? ?   ?
        structi = Struct(self._endian + b'2i i f i f 5i') # 11...decent

        #structi = Struct(self._endian + b'17i') # 17...not a chance

        # i: id
        # i/f: ???
        # 3f
        #structi = Struct(self._endian + b'i i 3f 28i') # 33...better?...still not right

        #structi = Struct(self._endian + b'51i') # 17*3..nope
        #structi = Struct(self._endian + b'71i') # 71...nope
        #structi = Struct(self._endian + b'187i') # 11*17...

        ntotal = 4 * 11

        nentries = (len(data) - n) // ntotal
        leftover = (len(data) - n) % ntotal
        assert leftover == 0, f'ndata={len(data)-n} leftover={leftover}'
        nfailed = 0
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            self.log.debug(out)
            n += ntotal
            continue
            (nid, cp, x1, x2, x3, cd, ps, seid) = out
            if self.is_debug_file:
                self.binary_debug.write('  GRID=%s\n' % str(out))
            if nid < 10000000:
                # cd can be < 0
                if ps == 0:
                    ps = ''
                node = GRID(nid, np.array([x1, x2, x3]), cp, cd, ps, seid)
                self._type_to_id_map['GRID'].append(nid)
                self.nodes[nid] = node
                #if nid in self.nodes:
                    #self.reject_lines.append(str(node))
                #else:
                #self.nodes[nid] = node
                #self.add_node(node)
            else:
                #self.log.warning('*nid=%s cp=%s x1=%-5.2f x2=%-5.2f x3=%-5.2f cd=%-2s ps=%s '
                                 #'seid=%s' % (nid, cp, x1, x2, x3, cd, ps, seid))
                #node = GRID(nid, np.array([x1, x2, x3]), cp, cd, ps, seid)
                #self.rejects.append(str(node))
                nfailed += 1
            n += ntotal
        self.increase_card_count('GRID', nentries - nfailed)
        return n

    def _read_grid(self, data: bytes, n: int) -> int:
        """(4501,45,1) - the marker for Record 17"""
        if self.table_name == b'GEOM1N' and self.factor == 1:
            n2 = self._read_grid_11(data, n)
        else:
            n2 = self._read_grid_8(data, n)
        return n2

    def _read_grid_8(self, data: bytes, n: int) -> int:  # 21.8 sec, 18.9
        """(4501,45,1) - the marker for Record 17"""
        structi = Struct(mapfmt(self._endian + b'ii3f3i', self.size))
        ntotal = 32 * self.factor
        ndatai = len(data) - n
        nentries = ndatai // ntotal
        assert nentries > 0, nentries
        nfailed = 0
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            (nid, cp, x1, x2, x3, cd, ps, seid) = out
            if self.is_debug_file:
                self.binary_debug.write('  GRID=%s\n' % str(out))
            #if nid < 10000000:
            # cd can be < 0
            if ps == 0:
                ps = ''
            node = GRID(nid, np.array([x1, x2, x3]), cp, cd, ps, seid)
            self._type_to_id_map['GRID'].append(nid)
            self.nodes[nid] = node
            #if nid in self.nodes:
                #self.reject_lines.append(str(node))
            #else:
            #self.nodes[nid] = node
            #self.add_node(node)
            #else:
                #self.log.warning('*nid=%s cp=%s x1=%-5.2f x2=%-5.2f x3=%-5.2f cd=%-2s ps=%s '
                                 #'seid=%s' % (nid, cp, x1, x2, x3, cd, ps, seid))
                #node = GRID(nid, np.array([x1, x2, x3]), cp, cd, ps, seid)
                #self.rejects.append(str(node))
                #nfailed += 1
            n += ntotal
        assert ndatai % ntotal == 0, f'ndatai={ndatai} ntotal={ntotal} leftover={ndatai % ntotal}'
        self.increase_card_count('GRID', nentries - nfailed)
        return n

    def _read_grid_11(self, data: bytes, n: int) -> int:  # 21.8 sec, 18.9
        """(4501,45,1) - the marker for Record 17"""
        ntotal = 44
        structi = Struct(self._endian + b'ii 3d 2ii')

        ndatai = len(data) - n
        nentries = ndatai // ntotal
        assert ndatai % ntotal == 0
        assert nentries > 0, nentries
        nfailed = 0
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            #cp, x1, x2, x3, cd, ps, seid
            (nid, cp, x1, x2, x3, cd, ps, seid) = out
            if self.is_debug_file:
                self.binary_debug.write('  GRID=%s\n' % str(out))
            #print(f'nid={nid}, cp={cp} x1={x1}, x2={x2}, x3={x3}, cd={cd} ps={ps}, seid={seid}')
            assert cd >= 0, f'nid={nid}, cp={cp} x1({x1}, {x2}, {x3}), cd={cd} ps={ps}, seid={seid}'
            assert seid == 0, f'nid={nid}, cp={cp} x1({x1}, {x2}, {x3}), cd={cd} ps={ps}, seid={seid}'
            if nid < 10000000:
                # cd can be < 0
                if ps == 0:
                    ps = ''
                node = GRID(nid, np.array([x1, x2, x3]), cp, cd, ps, seid)
                self._type_to_id_map['GRID'].append(nid)
                self.nodes[nid] = node
                #if nid in self.nodes:
                    #self.reject_lines.append(str(node))
                #else:
                #self.nodes[nid] = node
                #self.add_node(node)
            else:
                #self.log.warning('*nid=%s cp=%s x1=%-5.2f x2=%-5.2f x3=%-5.2f cd=%-2s ps=%s '
                                 #'seid=%s' % (nid, cp, x1, x2, x3, cd, ps, seid))
                #node = GRID(nid, np.array([x1, x2, x3]), cp, cd, ps, seid)
                #self.rejects.append(str(node))
                nfailed += 1
            n += ntotal
        assert ndatai % ntotal == 0, f'ndatai={ndatai} ntotal={ntotal} leftover={ndatai % ntotal}'
        self.increase_card_count('GRID', nentries - nfailed)
        return n

    def _read_seqgp(self, data: bytes, n: int) -> int:
        """(5301,53,4) - the marker for Record 27"""
        struct_2i = Struct(self._endian + b'2i')
        nentries = (len(data) - n) // 8
        for unused_i in range(nentries):
            edata = data[n:n + 8]  # 2*4
            out = struct_2i.unpack(edata)
            # (nid, seid) = out
            if self.is_debug_file:
                self.binary_debug.write('  SEQGP=%s\n' % str(out))
            seqgp = SEQGP.add_op2_data(out)
            self._add_seqgp_object(seqgp)
            n += 8
        self.increase_card_count('SEQGP', nentries)
        return n

    def _read_point(self, data: bytes, n: int) -> int:
        """
        POINT(6001,60,377)
        """
        s = Struct(self._endian + b'2i3f')
        ntotal = 20
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]  # 5*4
            out = s.unpack(edata)
            # (nid, cid, x1, x2, x3) = out
            if self.is_debug_file:
                self.binary_debug.write('  POINT=%s\n' % str(out))
            point = POINT.add_op2_data(out)
            self._add_point_object(point)
            n += ntotal
        self.increase_card_count('POINT', nentries)
        return n

    #def _read_cmass2(self, data: bytes, n: int) -> int:
        #struct_i4fi = Struct(self._endian + b'if4i')
        #nentries = (len(data) - n) // 24
        #for unused_i in range(nentries):
            #edata = data[n:n + 24]  # 6*4
            #out = struct_i4fi.unpack(edata)
            ## (eid, mass, g1, g2, c1, c2) = out
            #if self.is_debug_file:
                #self.binary_debug.write('  CMASS2=%s\n' % str(out))
            #self.log.debug('  CMASS2=%s\n' % str(out))
            #element = CMASS2.add_op2_data(out)
            #print(element)
            #self.add_op2_mass(element)
            #n += 24
        #print(self.elements)
        #self.increase_card_count('CMASS2', nentries)
        #return n
        #return len(data)

    def _read_cvisc(self, data: bytes, n: int) -> int:
        """CVISC(3901,39,50) - the marker for Record 105"""
        struct_4i = Struct(self._endian + b'4i')
        ntotal = 16  # 4*4
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = struct_4i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CVISC=%s\n' % str(out))
            # (eid, pid, n1, n2) = out
            element = CVISC.add_op2_data(out)
            self.add_op2_element(element)
            n += ntotal
        self.increase_card_count('CVISC', nentries)
        return n


    def _read_extrn(self, data: bytes, n: int) -> int:
        """
        Record - EXTRN(1627,16,463)
        Word Name Type Description
        1 GID I Grid point identification numbers to connect external SE
        2 C   I Component numbers
        Words 1 through 2 repeat until (-1,-1) occurs

        (1, 123456, 2, 123456, 3, 123456, 4, 123456, 100001, 1, 100002, 1, 100003, 1, 100004, 1,
        -1, -1)
        """
        self.log.info('skipping EXTRN in GEOM1')
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
        # C:\NASA\m4\formats\git\examples\move_tpl\phsflux4.op2
        #(200000002, 3, 1002, 6, 12, 0, 0)
        # FEEDGE EDGEID GRID1 GRID2 CIDBC GEOMIN ID1 ID2
        #FEEDGE    1002    6     12
        #self.show_data(data[12:])
        ntotal = 28 * self.factor  # 7*4
        if self.size == 4:
            s = Struct(self._endian + b'4i 4s 2i') #expected
        else:
            s = Struct(self._endian + b'4q 8s 2q') #expected
        #s = Struct(self._endian + b'7i')
        nelements = (len(data) - n)// ntotal
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            #print(out)
            edge_id, n1, n2, cid, geomin, geom1, geom2 = out # expected

            if self.is_debug_file:
                self.binary_debug.write('  FEEDGE=%s\n' % str(out))

            geomin = geomin.rstrip()
            if geomin == b'POIN':
                geomin_str = 'POINT'
            elif geomin == b'GMCU':
                geomin_str = 'GMCURV'
            else:  # pragma: no cover
                raise RuntimeError(geomin)

            if cid == -1:
                cid = None
            unused_elem = self.add_feedge(edge_id, [n1, n2], cid, [geom1, geom2], geomin=geomin_str)
            n += ntotal
        self.card_count['FEEDGE'] = nelements
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
        size = self.size
        if size == 4:
            struct_i = self.struct_i
            structi = Struct(b'i 8s ii')
        else:
            struct_i = self.struct_q
            structi = Struct(b'q 16s qq')

        ntotal1 = 20 * self.factor
        ntotal2 = 64 * self.factor
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
            self.add_gmcurv(curve_id, group, data_split, cid_in=cid_in, cid_bc=cid_bc)
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
        structi = Struct(self._endian + b'8i')
        nentries = (len(data) - n) // 32
        for unused_i in range(nentries):
            edata = data[n:n + 32]
            out = structi.unpack(edata)
            (face_id, n1, n2, n3, n4, cid, surf_id1, surf_id2) = out
            if self.is_debug_file:
                self.binary_debug.write('  FEFACE=%s\n' % str(out))

            nodes = [n1, n2, n3, n4]
            surf_ids = [surf_id1, surf_id2]
            feface = self.add_feface(face_id, nodes, cid, surf_ids)
            n += 32
        self.increase_card_count('FEFACE', nentries)
        return n

    def _read_gmsurf(self, data: bytes, n: int) -> int:
        self.log.info('skipping GMSURF in GEOM1')
        return len(data)

    def _read_gmcord(self, data: bytes, n: int) -> int:
        self.log.info('skipping GMCORD in GEOM1')
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
        ntotal = 32 * self.factor # 4*8
        nentries = (len(data) - n) // ntotal
        structi = Struct(mapfmt(self._endian + b'4if3i', self.size))

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

            if self.is_debug_file:
                self.binary_debug.write('  SEBULK=%s\n' % str(out))
            #media,
            sebulk = self.add_sebulk(seid, superelement_type, rseid,
                                     method=method, tol=tol, loc=loc, unitno=unit)
            sebulk.validate()
            n += 32 * self.factor
        self.increase_card_count('SEBULK', nentries)
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
        structi = Struct(self._endian + b'7i')
        nentries = (len(data) - n) // 28 # 4*7
        for unused_i in range(nentries):
            edata = data[n:n + 28]
            out = structi.unpack(edata)
            (seid, ga1, ga2, ga3, gb1, gb2, gb3) = out
            if self.is_debug_file:
                self.binary_debug.write('  SELOC=%s\n' % str(out))
            self.add_seloc(seid, [ga1, ga2, ga3], [gb1, gb2, gb3])
            n += 28
        self.increase_card_count('SELOC', nentries)
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
        struct2i = Struct(self._endian + b'2i') # 8
        struct5i = Struct(self._endian + b'5i') # 20
        #struct2i_3f = Struct(self._endian + b'2i3f') # 20

        nentries = (len(data) - n) // 28 # 4*7
        for unused_i in range(nentries):
            edata1 = data[n:n + 8]  # 4*2
            edata2 = data[n+8:n + 28]  # 4*7
            out = struct2i.unpack(edata1)
            (seid, mirror_type) = out
            if mirror_type == 1:
                g1, g2, g3, unused_junk1, unused_junk2 = struct5i.unpack(edata2)
                self.add_sempln(seid, g1, g2, g3)
            else:
                raise NotImplementedError(mirror_type)
            if self.is_debug_file:
                self.binary_debug.write('  SEMPLN=%s\n' % str(out))
            n += 28
        self.increase_card_count('SEMPLN', nentries)
        return n

    def _read_selabel(self, data: bytes, n: int) -> int:
        """
        Record 22 -- SELABEL(1027,10,459)

        Word Name Type Description
        1 SEID I Superelement identification number
        2 LABEL(14) CHAR4 Label associated with superelement SEID
        """
        structi = Struct(self._endian + b'i14s') # 18
        structi = Struct(self._endian + b'i56s') # 60
        nentries = (len(data) - n) // 60 # 4+18
        for unused_i in range(nentries):
            edata = data[n:n + 60]
            out = structi.unpack(edata)
            (seid, label) = out
            label = label.decode(self._encoding).rstrip()
            if self.is_debug_file:
                self.binary_debug.write('  SELABEL=%s\n' % str(out))
            selabel = self.add_selabel(seid, label)
            selabel.validate()
            n += 60
        self.increase_card_count('SELABEL', nentries)
        return n
