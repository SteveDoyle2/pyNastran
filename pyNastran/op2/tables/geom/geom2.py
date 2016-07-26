# pylint: disable=W0612,C0103,C0302,W0613,C0111,R0914,R0201
from struct import unpack, Struct
from six import b
from six.moves import range

from pyNastran.bdf.cards.elements.elements import CGAP, PLOTEL
from pyNastran.bdf.cards.elements.damper import (CDAMP1, CDAMP2, CDAMP3,
                                                 CDAMP4, CDAMP5, CVISC)
from pyNastran.bdf.cards.elements.springs import CELAS1, CELAS2, CELAS3, CELAS4
from pyNastran.bdf.cards.elements.shell import (CTRIA3, CQUAD4, CTRIA6,
                                                CQUADR, CQUAD8, CQUAD, CQUADX,
                                                CSHEAR)
from pyNastran.bdf.cards.elements.rods import CROD, CTUBE, CONROD
from pyNastran.bdf.cards.elements.bars import CBAR
from pyNastran.bdf.cards.elements.beam import CBEAM
from pyNastran.bdf.cards.elements.mass import (CONM1, CONM2, CMASS1, CMASS2,
                                               CMASS3, CMASS4)
from pyNastran.bdf.cards.elements.solid import (CTETRA4, CPYRAM5, CPENTA6, CHEXA8,
                                                CTETRA10, CPYRAM13, CPENTA15, CHEXA20,)
from pyNastran.bdf.cards.thermal.thermal import CHBDYG, CONV, CHBDYP, CHBDYE, CONVM
from pyNastran.bdf.cards.nodes import SPOINTs
from pyNastran.op2.tables.geom.geom_common import GeomCommon
from pyNastran.bdf.cards.elements.bush import CBUSH

class GEOM2(GeomCommon):
    """defines methods for reading op2 elements"""

    def _read_geom2_4(self, data, ndata):
        return self._read_geom_4(self._geom2_map, data, ndata)

    def __init__(self):
        GeomCommon.__init__(self)
        self._geom2_map = {
            # per dmap-nx-10.pdf
            # BEAMAERO(2601,26,0)
            # CAABSF(2708,27,59)
            # CAXIF2(2108,21,224)
            # CAXIF3(2208,22,225)
            # CAXIF4(2308,23,226)
            (2408, 24, 180): ['CBAR', self._read_cbar],         # record 8
            (4001, 40, 275): ['CBARAO', self._read_cbarao],     # record 9  - not done
            (5408, 54, 261): ['CBEAM', self._read_cbeam],       # record 10
            (11401, 114, 9016): ['CBEAMP', self._read_cbeamp],  # record 11 - not done
            (4601, 46, 298): ['CBEND', self._read_cbend],     # record 12 - not done
            (2608, 26, 60): ['CBUSH', self._read_cbush],      # record 13
            (5608, 56, 218): ['CBUSH1D', self._read_cbush1d], # record 14 - not done
            (2315, 23, 146): ['CCONE', self._read_ccone],     # record 15 - not done
            (201, 2, 69): ['CDAMP1', self._read_cdamp1],      # record 16
            (301, 3, 70): ['CDAMP2', self._read_cdamp2],      # record 17
            (401, 4, 71): ['CDAMP3', self._read_cdamp3],      # record 18
            (501, 5, 72): ['CDAMP4', self._read_cdamp4],      # record 19
            (10608, 106, 404): ['CDAMPS', self._read_cdamp5], # record 20

            (601, 6, 73): ['CELAS1', self._read_celas1],      # record 29
            (701, 7, 74): ['CELAS2', self._read_celas2],      # record 30
            (801, 8, 75): ['CELAS3', self._read_celas3],      # record 31
            (901, 9, 76): ['CELAS4', self._read_celas4],      # record 32
            # record 33
            (9801, 98, 506): ['CFAST', self._read_cfast],     # record 34 - not done
            # record 34
            # record 35
            (8515, 85, 209): ['CFLUID2', self._read_cfluid2],  # record 35 - not done
            (8615, 86, 210): ['CFLUID3', self._read_cfluid3],  # record 36 - not done
            (8715, 87, 211): ['CFLUID4', self._read_cfluid4],  # record 37 - not done
            (1908, 19, 104): ['CGAP', self._read_cgap],        # record 39 - buggy
            # record 40
            # record 41
            # record 42
            (10808, 108, 406): ['CHBDYG', self._read_chbdyg], # record 43
            (10908, 109, 407): ['CHBDYP', self._read_chbdyp], # record 44 - not done
            (7308, 73, 253): ['CHEXA', self._read_chexa],     # record 45
            (1001, 10, 65): ['CMASS1', self._read_cmass1],    # record 51
            (1101, 11, 66): ['CMASS2', self._read_cmass2],    # record 52
            (1201, 12, 67): ['CMASS3', self._read_cmass3],    # record 53
            (1301, 13, 68): ['CMASS4', self._read_cmass4],    # record 54
            (2508, 25, 0): ['CMFREE', self._read_cmfree],     # record 55 - not done
            (1401, 14, 63): ['CONM1', self._read_conm1],      # record 56 - not done
            (1501, 15, 64): ['CONM2', self._read_conm2],      # record 57
            (1601, 16, 47): ['CONROD', self._read_conrod],    # record 58
            (12701, 127, 408): ['CONV', self._read_conv],     # record 59 - not tested
            (8908, 89, 422): ['CONVM', self._read_convm],     # record 60 - not tested
            # record 61
            (4108, 41, 280): ['CPENTA', self._read_cpenta],   # record 62
            # record 63
            # record 64
            # record 65
            # record 66
            # record 67
            (9108, 91, 507): ['CQUAD', self._read_cquad],       # record 68 - not tested
            (2958, 51, 177): ['CQUAD4', self._read_cquad4],     # record 69 - maybe buggy on theta/Mcsid field
            (13900, 139, 9989): ['CQUAD4', self._read_cquad4],  # record 70 - maybe buggy on theta/Mcsid field
            (4701, 47, 326): ['CQUAD8', self._read_cquad8],     # record 71 - maybe buggy on theta/Mcsid field
            # record 72
            # record 73
            (8009, 80, 367): ['CQUADR', self._read_cquadr],   # record 74 - not tested
            (9008, 90, 508): ['CQUADX', self._read_cquadx],   # record 75 - not tested
            # record 76
            # record 77
            # record 78
            # record 79
            (3001, 30, 48): ['CROD', self._read_crod],        # record 80
            # record 81
            # record 82
            # record 83
            # record 84
            # record 85
            (12201, 122, 9013): ['CTETP', self._read_ctetrap],  # record 86 - not done
            (5508, 55, 217): ['CTETRA', self._read_ctetra],     # record 87
            # record 88
            # record 89
            # record 90
            # record 91
            # record 92
            (5959, 59, 282): ['CTRIA3', self._read_ctria3],   # record 93 - maybe buggy on theta/Mcsid field
            # record 94
            (4801, 48, 327): ['CTRIA6', self._read_ctria6],   # record 95 - buggy
            # record 96
            # record 97
            (9200, 92, 385): ['CTRIAR', self._read_ctriar],   # record 98  - not done
            # record 99
            (6108, 61, 107): ['CTRIAX6', self._read_ctriax6], # record 100 - not done
            # record 101
            # record 102
            (3701, 37, 49): ['CTUBE', self._read_ctube],      # record 103
            (3901, 39, 50): ['CVISC', self._read_cvisc],      # record 104 - not done
            # record 105
            # record 106
            # record 107
            # record 108
            # record 109
            # record 110
            # record 111
            # record 112
            # record 113
            (5201, 52, 11):   ['PLOTEL', self._read_plotel],    # record 114 - not done
            # record 115
            # record 116
            # record 117
            (5551, 49, 105): ['SPOINT', self._read_spoint],     # record 118
            (11601, 116, 9942): ['VUBEAM', self._read_vubeam],  # record 119 - not done
            (2108, 21, 224): ['CAXIF2', self._read_fake],
            (3101, 31, 61): ['CSHEAR', self._read_fake],
            (4301, 43, 28): ['GENEL', self._read_fake],
            (5601, 56, 296): ['SESET', self._read_fake],
            (6808, 68, 114): ['CDUM8', self._read_fake],
            (6908, 69, 115): ['CDUM9', self._read_fake],
            (7409, 74, 9991): ['CHEXPR', self._read_fake],
            (7509, 75, 9992): ['CPENPR', self._read_fake],
            (7609, 76, 9993): ['CTETPR', self._read_fake],
            (8100, 81, 381): ['CHACAB', self._read_fake],
            (8200, 82, 383): ['CHACBR', self._read_fake],
            (8308, 83, 405): ['CHBDYE', self._read_chbdye],
            (11201, 112, 9940): ['VUQUAD4', self._read_fake],
            (12801, 128, 417): ['RADBC', self._read_fake],
            (2708, 27, 59): ['CAABSF', self._read_fake],
            (3201, 32, 478): ['GMBNDC', self._read_fake],
            (13900, 139, 9984): ['', self._read_fake],
            (14000, 140, 9990): ['', self._read_fake],
            (16000, 160, 9988): ['', self._read_fake],
            (16100, 161, 9986): ['', self._read_fake],
            (16300, 163, 9989): ['', self._read_fake],
            (16700, 167, 9981): ['', self._read_fake],
            (16800, 168, 9978): ['', self._read_fake],
            (16500, 165, 9987): ['', self._read_fake],
            (5008, 50, 258): ['', self._read_fake],
            (16400, 164, 9983) : ['', self._read_fake],
            (11000, 110, 6667): ['', self._read_fake],
            (12301, 123, 9921): ['', self._read_fake],
            (12401, 124, 9922): ['', self._read_fake],
            (12600, 126, 6661): ['', self._read_fake],
            (14700, 147, 6662): ['', self._read_fake],
            (7309, 73, 0): ['', self._read_fake],
            (17200, 172, 6663): ['', self._read_fake],
            (17300, 173, 6664): ['', self._read_fake],
            (11501, 115, 9941): ['', self._read_fake],    # record
            (12501, 125, 9923): ['', self._read_fake],    # record
            (3401, 34, 9600): ['', self._read_fake],    # record
            (7701, 77, 8881): ['', self._read_fake],  # record
            (2901, 29, 9601): ['', self._read_fake],  # record
            (16600, 166, 9985) : ['', self._read_fake],  # record
            (16200, 162, 9982) : ['', self._read_fake],  # record
            (16900, 169, 9977) : ['', self._read_fake],  # record
            (1701, 17, 980) : ['', self._read_fake],  # record
            (1801, 18, 986) : ['', self._read_fake],  # record
            (17200, 172, 1000) : ['CPYRAM', self._read_fake],  # record
            (23500, 235, 6662) : ['', self._read_fake],  # record
            (23800, 238, 6665) : ['', self._read_fake],  # record
            (23900, 239, 6666) : ['', self._read_fake],  # record

            (2208, 22, 225): ['CAXIF3', self._read_fake],  # record 7
            (4408, 44, 227): ['CSLOT3', self._read_fake],  # record 85
            (4508, 45, 228): ['CSLOT4', self._read_fake],  # record 86-R
            (12901, 129, 482): ['GMBNDS', self._read_fake],  # record 112
            (7801, 78, 8883): ['SINT', self._read_fake],  # record 118
            (8801, 88, 984) : ['CPLSTS3', self._read_fake],  # record
            (8401, 84, 985) : ['CPLSTS4', self._read_fake],  # record
            (17000, 170, 9980): ['CQDX4FD', self._read_fake],  # record
            (17100, 171, 9979): ['CQDX9FD', self._read_fake],  # record

            (1976, 1, 1996) : ['', self._read_fake],  # record
            (6120, 1, 60434) : ['', self._read_fake],  # record
            (2024, 1001, 2024) : ['', self._read_fake],  # record
            (801, 1, 572) : ['', self._read_fake],  # record

            (5701, 57, 981) : ['CPLSTN4', self._read_fake],  # record
            (5801, 58, 982) : ['CPLSTN6', self._read_fake],  # record
            (1801, 18, 986) : ['CPLSTS6', self._read_fake],  # record
            (6111, 61, 996) : ['CTRAX3', self._read_fake],  # record
            (6112, 61, 997) : ['CQUADX4', self._read_fake],  # record
            (6113, 61, 998) : ['CTRAX6', self._read_fake],  # record
            (6114, 61, 999) : ['CQUADX8', self._read_fake],  # record
            (3501, 35, 1) : ['', self._read_fake],  # record
            (1001, 100, 10000) : ['', self._read_fake],  # record
            (1118, 1, 1874) : ['', self._read_fake],  # record

            # NX specific
            (17200, 172, 1000) : ['CPYRAM', self._read_cpyram],
            (25700,257,9948) : ['CPYRA5FD', self._read_cpyram],
            (25800,258,9947) : ['CPYRA13F', self._read_cpyram],
            (7909, 79, 9946) : ['CPYRAMPR', self._read_cpyram],
        }

    def add_element(self, elem, allow_overwrites=True):
        raise RuntimeError('this should be overwritten')

    def add_op2_element(self, elem):
        if elem.eid <= 0:
            self.log.debug(elem)
            raise ValueError(elem)
            #return

        if elem.eid > 100000000:
            raise RuntimeError('bad parsing...')

        if elem.type in ['CTRIA6', 'CQUAD8']:
            for nid in elem.nodes:
                if nid == -1:
                    nid = None
        else:
            for nid in elem.nodes:
                if nid == -1:
                    assert nid > 0, elem
        self.add_element(elem, allow_overwrites=True)
        #print(str(elem)[:-1])

# 1-AEROQ4 (???)
# AEROT3   (???)
# 1-BEAMAERO (1701,17,0)
# 2-CAABSF (2708,27,59)
# 3-CAXIF2 (2108,21,224)
# 4-CAXIF3 (2208,22,225)
# 5-CAXIF4 (2308,23,226)

    def _read_cbar(self, data, n):
        """
        CBAR(2408,24,180) - the marker for Record 8
        """
        nelements = (len(data) - n) // 64
        for i in range(nelements):
            edata = data[n:n + 64]  # 16*4
            f, = self.struct_i.unpack(edata[28:32])
            if f == 0:
                out = unpack(b(self._endian + '4i3f3i6f'), edata)
                (eid, pid, ga, gb, x1, x2, x3, f, pa, pb,
                 w1a, w2a, w3a, w1b, w2b, w3b) = out
                data_in = [[eid, pid, ga, gb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                           [f, x1, x2, x3]]
            elif f == 1:
                out = unpack(b(self._endian + '4i3f3i6f'), edata)
                (eid, pid, ga, gb, x1, x2, x3, f, pa, pb,
                 w1a, w2a, w3a, w1b, w2b, w3b) = out
                data_in = [[eid, pid, ga, gb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                           [f, x1, x2, x3]]
            elif f == 2:
                out = unpack(b(self._endian + '7if2i6f'), edata)
                (eid, pid, ga, gb, g0, junk, junk, f, pa,
                 pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                data_in = [[eid, pid, ga, gb, pa, pb, w1a,
                            w2a, w3a, w1b, w2b, w3b], [f, g0]]
            else:
                raise RuntimeError('invalid f value...f=%s' % (f))
            elem = CBAR.add_op2_data(data_in)
            self.add_op2_element(elem)
            n += 64
        self.card_count['CBAR'] = nelements
        return n

    def _read_cbarao(self, data, n):
        """
        CBARAO(4001,40,275) - the marker for Record 9
        """
        self.log.debug('skipping CBARAO in GEOM2\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping CBARAO in GEOM2\n')
        return len(data)

    def _read_cbeam(self, data, n):
        """
        CBEAM(5408,54,261) - the marker for Record 10
        """
        nelements = (len(data) - n) // 72
        for i in range(nelements):
            edata = data[n:n + 72]  # 18*4
            fe, = self.struct_i.unpack(edata[40:44])

            # per DMAP: F = FE bit-wise AND with 3
            f = fe & 3
            if f == 0:  # basic cid
                out = unpack(b(self._endian + '6i3f3i6f'), edata)
                (eid, pid, ga, gb, sa, sb, x1, x2, x3, fe, pa,
                 pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                #self.log.info('CBEAM: eid=%s fe=%s f=%s; basic cid' % (eid, fe, f))
                data_in = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                           [f, x1, x2, x3]]
            elif f == 1:  # global cid
                out = unpack(b(self._endian + '6i3f3i6f'), edata)
                (eid, pid, ga, gb, sa, sb, x1, x2, x3, fe, pa,
                 pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                #self.log.info('CBEAM: eid=%s fe=%s f=%s; global cid' % (eid, fe, f))
                data_in = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                           [f, x1, x2, x3]]
            elif f == 2:  # grid option
                out = unpack(b(self._endian + '12i6f'), edata)
                (eid, pid, ga, gb, sa, sb, g0, xx, xx, fe, pa,
                 pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                #self.log.info('CBEAM: eid=%s fe=%s f=%s; grid option' % (eid, fe, f))
                data_in = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                           [f, g0]]
            else:
                raise RuntimeError('invalid f value...f=%r' % f)
            elem = CBEAM.add_op2_data(data_in, f)
            self.add_op2_element(elem)
            n += 72
        self.card_count['CBEAM'] = nelements
        return n

    def _read_cbeamp(self, data, n):
        """
        CBEAMP(11401,114,9016) - the marker for Record 11
        """
        self.log.debug('skipping CBEAMP in GEOM2\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping CBEAMP in GEOM2\n')
        return len(data)

    def _read_cbend(self, data, n):
        """
        CBEND(4601,46,298) - the marker for Record 12
        """
        self.log.debug('skipping CBEND in GEOM2\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping CBEND in GEOM2\n')
        return len(data)

    def _read_cbush(self, data, n):
        """
        CBUSH(2608,26,60) - the marker for Record 13
        """
        nelements = (len(data) - n) // 56
        struct_obj1 = Struct(b(self._endian + '4i iii i ifi3f'))
        struct_obj2 = Struct(b(self._endian + '4i fff i ifi3f'))
        for i in range(nelements):
            edata = data[n:n + 56]  # 14*4
            out = struct_obj1.unpack(edata)
            eid, pid, ga, gb, five, six, seven, f, cid, s, ocid, s1, s2, s3 = out
            si = [s1, s2, s3]
            if f == -1: # Use Element CID below for orientation
                x = [None, None, None]
                g0 = None
            elif f in [0, 1]:
                # 5, 6, 7, f
                #0:4 4:8 8:12 12:16 16:20 20:24 24:28 28:32 32:36
                #0   1   2    3     4     5     6     7     8
                #x1, x2, x3, f2 = unpack('3f i', edata[20:36])
                out = struct_obj2.unpack(edata)
                eid, pid, ga, gb, x1, x2, x3, f2, cid, s, ocid, s1, s2, s3 = out

                assert f == f2, 'f=%s f2=%s' % (f, f2)
                x = [x1, x2, x3]
                g0 = None
            elif f == 2:
                x = [None, None, None]
                g0 = five
            else:
                raise RuntimeError('invalid f value...f=%r' % f)
            if cid == -1:
                cid = None
            data_in = [[eid, pid, ga, gb, cid, s, ocid, si], x, g0]

            elem = CBUSH.add_op2_data(data_in, f)
            self.add_op2_element(elem)
            n += 56
        self.card_count['CBUSH'] = nelements
        return n

    def _read_cbush1d(self, data, n):
        """
        CBUSH1D(5608,56,218) - the marker for Record 14
        """
        self.log.debug('skipping CBUSH1D in GEOM2\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping CBUSH1D in GEOM2\n')
        return len(data)

    def _read_ccone(self, data, n):
        """
        CCONE(2315,23,0) - the marker for Record 15
        """
        self.log.debug('skipping CCONE in GEOM2\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping CCONE in GEOM2\n')
        return len(data)

    def _read_cdamp1(self, data, n):
        """
        CDAMP1(201,2,69) - the marker for Record 16
        """
        nelements = (len(data) - n) // 24
        for i in range(nelements):
            edata = data[n:n + 24]  # 6*4
            out = unpack(b(self._endian + '6i'), edata)
            if self.is_debug_file:
                self.binary_debug.write('  CDAMP1=%s\n' % str(out))
            (eid, pid, g1, g2, c1, c2) = out
            elem = CDAMP1.add_op2_data(out)
            self.add_op2_element(elem)
            n += 24
        self.card_count['CDAMP1'] = nelements
        return n

    def _read_cdamp2(self, data, n):
        """
        CDAMP2(301,3,70) - the marker for Record 17
        """
        nelements = (len(data) - n) // 24
        for i in range(nelements):
            edata = data[n:n + 24]  # 6*4
            out = unpack(b(self._endian + 'if4i'), edata)
            if self.is_debug_file:
                self.binary_debug.write('  CDAMP2=%s\n' % str(out))
            (eid, bdamp, g1, g2, c1, c2) = out
            elem = CDAMP2.add_op2_data(out)
            self.add_op2_element(elem)
            n += 24
        self.card_count['CDAMP2'] = nelements
        return n

    def _read_cdamp3(self, data, n):
        """
        CDAMP3(401,4,71) - the marker for Record 18
        """
        s = Struct(b(self._endian + '4i'))
        nelements = (len(data) - n) // 16
        for i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CDAMP3=%s\n' % str(out))
            (eid, pid, s1, s2) = out
            elem = CDAMP3.add_op2_data(out)
            self.add_op2_element(elem)
            n += 16
        self.card_count['CDAMP3'] = nelements
        return n

    def _read_cdamp4(self, data, n):
        """
        CDAMP4(501,5,72) - the marker for Record 19
        """
        s = Struct(b(self._endian + 'ifii'))
        nelements = (len(data) - n) // 16
        for i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CDAMP4=%s\n' % str(out))
            (eid, bdamp, s1, s2) = out
            elem = CDAMP4.add_op2_data(out)
            self.add_op2_element(elem)
            n += 16
        self.card_count['CDAMP4'] = nelements
        return n

    def _read_cdamp5(self, data, n):
        """
        CDAMP5(10608,106,404) - the marker for Record 20
        """
        s = Struct(b(self._endian + '4i'))
        nelements = (len(data) - n) // 16
        for i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CDAMP5=%s\n' % str(out))
            (eid, pid, s1, s2) = out
            elem = CDAMP5.add_op2_data(out)
            self.add_op2_element(elem)
            n += 16
        self.card_count['CDAMP5'] = nelements
        return n

# CDUM2
# CDUM3
# CDUM4
# CDUM5
# CDUM6
# CDUM7
# CDUM8
# CDUM9

    def _read_celas1(self, data, n):
        """
        CELAS1(601,6,73) - the marker for Record 29
        """
        ntotal = 24  # 6*4
        s = Struct(b(self._endian + '6i'))
        nelements = (len(data) - n) // ntotal
        for i in range(nelements):
            edata = data[n:n+24]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CELAS1=%s\n' % str(out))
            (eid, pid, g1, g2, c1, c2) = out
            elem = CELAS1.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CELAS1'] = nelements
        return n

    def _read_celas2(self, data, n):
        """
        CELAS2(701,7,74) - the marker for Record 30
        """
        s1 = Struct(b(self._endian + 'if4iff'))
        ntotal = 32
        nelements = (len(data) - n) // ntotal
        for i in range(nelements):
            edata = data[n:n+32]
            out = s1.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CELAS2=%s\n' % str(out))
            (eid, k, g1, g2, c1, c2, ge, s) = out
            elem = CELAS2.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CELAS2'] = nelements
        return n

    def _read_celas3(self, data, n):
        """
        CELAS3(801,8,75) - the marker for Record 31
        """
        ntotal = 16  # 4*4
        s = Struct(b(self._endian + '4i'))
        nelements = (len(data) - n) // ntotal
        for i in range(nelements):
            edata = data[n:n+16]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CELAS3=%s\n' % str(out))
            (eid, pid, s1, s2) = out
            elem = CELAS3.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CELAS3'] = nelements
        return n

    def _read_celas4(self, data, n):
        """
        CELAS4(901,9,76) - the marker for Record 32
        """
        s = Struct(b(self._endian + 'ifii'))
        nelements = (len(data) - n) // 16
        for i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CELAS4=%s\n' % str(out))
            (eid, k, s1, s2) = out
            elem = CELAS4.add_op2_data(out)
            self.add_op2_element(elem)
            n += 16
        self.card_count['CELAS4'] = nelements
        return n

    def _read_cfast(self, data, n):
        """
        CFAST(9801,98,506) - the marker for Record ???
        """
        self.log.debug('skipping CFAST in GEOM2\n')
        return len(data)

# CFASTP

    def _read_cfluid2(self, data, n):
        """
        CFLUID2(8515,85,209) - the marker for Record 35
        """
        self.log.debug('skipping CFLUID2 in GEOM2\n')
        return len(data)

    def _read_cfluid3(self, data, n):
        """
        CFLUID3(8615,86,210) - the marker for Record 36
        """
        self.log.debug('skipping CFLUID3 in GEOM2\n')
        return len(data)

    def _read_cfluid4(self, data, n):
        """
        CFLUID4(8715,87,211) - the marker for Record 37
        """
        self.log.debug('skipping CFLUID4 in GEOM2\n')
        return len(data)

# CINT

    def _read_cgap(self, data, n):
        """
        CGAP(1908,19,104) - the marker for Record 39
        """
        s1 = Struct(b(self._endian + '4i3fii'))
        nelements = (len(data) - n) // 36
        for i in range(nelements):
            edata = data[n:n + 36]  # 9*4
            out = s1.unpack(edata)
            (eid, pid, ga, gb, x1, x2, x3, f, cid) = out  # f=0,1
            g0 = None
            f2, = self.struct_i.unpack(edata[28:32])
            assert f == f2, 'f=%s f2=%s' % (f, f2)
            if f == 2:
                g0 = self.struct_i.unpack(edata[16:20])
                x1 = None
                x2 = None
                x3 = None
            else:
                assert f == 1, 'CGAP - f=%r f2=%r' % (f, f2)
                assert f2 == 1, 'CGAP - f=%r f2=%r' % (f, f2)
                #raise NotImplementedError('CGAP - f=%r f2=%r' % (f, f2))
            data_in = [eid, pid, ga, gb, g0, x1, x2, x3, cid]
            elem = CGAP.add_op2_data(data_in)
            self.add_op2_element(elem)
            n += 36
        self.card_count['CGAP'] = nelements
        return n

# CHACAB
# CHACBR

    def _read_chbdye(self, data, n):
        """
        CHBDYE(8308,83,405) - the marker for Record ???
        """
        ntotal = 28  # 7*4
        s = Struct(b(self._endian + '7i'))
        nelements = (len(data) - n) // ntotal
        for i in range(nelements):
            edata = data[n:n+28]
            out = s.unpack(edata)
            (eid, eid2, side, iviewf, iviewb, radmidf, radmidb) = out
            if self.is_debug_file:
                self.binary_debug.write('  CHBDYE=%s\n' % str(out))
            #self.log.debug('  CHBDYE=%s' % str(out))
            data_in = [eid, eid2, side, iviewf, iviewb, radmidf, radmidb]
            elem = CHBDYE.add_op2_data(data_in)
            self.add_thermal_element(elem)
            n += ntotal
        self.card_count['CHBDYE'] = nelements
        return n

    def _read_chbdyg(self, data, n):
        """
        CHBDYG(10808,108,406) - the marker for Record 43
        """
        ntotal = 64  # 16*4
        s = Struct(b(self._endian + '16i'))
        nelements = (len(data) - n) // ntotal
        for i in range(nelements):
            edata = data[n:n+64]
            out = s.unpack(edata)
            (eid, blank, Type, iviewf, iviewb, radmidf, radmidb, blank2,
             g1, g2, g3, g4, g5, g6, g7, g8) = out
            if self.is_debug_file:
                self.binary_debug.write('  CHBDYG=%s\n' % str(out))
            #self.log.debug('  CHBDYG=%s' % str(out))
            data_in = [eid, Type, iviewf, iviewb, radmidf, radmidb,
                       g1, g2, g3, g4, g5, g6, g7, g8]
            elem = CHBDYG.add_op2_data(data_in)
            self.add_thermal_element(elem)
            n += ntotal
        self.card_count['CHBDYG'] = nelements
        return n

    def _read_chbdyp(self, data, n):
        """
        CHBDYP(10908,109,407)
        """
        ntotal = 60  # 16*4
        s = Struct(b(self._endian + '12i 3f'))
        nelements = (len(data) - n) // ntotal
        for i in range(nelements):
            edata = data[n:n+60]
            out = s.unpack(edata)
            (eid, pid, Type, iviewf, iviewb, g1, g2, g0, radmidf, radmidb,
             dislin, ce, e1, e2, e3) = out
            if self.is_debug_file:
                self.binary_debug.write('  CHBDYP=%s\n' % str(out))
            #self.log.debug('  CHBDYP=%s' % str(out))
            data_in = [eid, pid, Type, iviewf, iviewb, g1, g2, g0, radmidf, radmidb,
             dislin, ce, e1, e2, e3]
            elem = CHBDYP.add_op2_data(data_in)
            self.add_thermal_element(elem)
            n += ntotal
        self.card_count['CHBDYP'] = nelements
        return n

    def _read_chexa(self, data, n):
        """
        CHEXA(7308,73,253) - the marker for Record 45
        """
        s = Struct(b(self._endian + '22i'))
        ntotal = 88  # 22*4
        nelements = (len(data) - n) // ntotal
        for i in range(nelements):
            edata = data[n:n+88]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CHEXA=%s\n' % str(out))
            (eid, pid, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10,
             g11, g12, g13, g14, g15, g16, g17, g18, g19, g20) = out

            data_in = [eid, pid, g1, g2, g3, g4, g5, g6, g7, g8, ]
            big_nodes = [g9, g10, g11, g12, g13, g14, g15, g16,
                         g17, g18, g19, g20]
            if sum(big_nodes) > 0:
                elem = CHEXA20.add_op2_data(data_in + big_nodes)
            else:
                elem = CHEXA8.add_op2_data(data_in)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CHEXA'] = nelements
        return n

# CHEXA20F
# CHEXAFD
# CHEXAL
# CHEXP
# CHEXPR

    def _read_cmass1(self, data, n):
        """
        CMASS1(1001,10,65) - the marker for Record 51
        """
        s = Struct(b(self._endian + '6i'))
        nelements = (len(data) - n) // 24
        for i in range(nelements):
            edata = data[n:n + 24]  # 6*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CMASS1=%s\n' % str(out))
            #(eid, pid, g1, g2, c1, c2) = out
            elem = CMASS1.add_op2_data(out)
            self.add_mass(elem)
            n += 24
        self.card_count['CMASS1'] = nelements
        return n

    def _read_cmass2(self, data, n):
        """
        CMASS2(1101,11,66) - the marker for Record 52
        """
        s = Struct(b(self._endian + 'if4i'))
        nelements = (len(data) - n) // 24
        for i in range(nelements):
            edata = data[n:n + 24]  # 6*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CMASS2=%s\n' % str(out))
            #(eid, m, g1, g2, c1, c2) = out
            elem = CMASS2.add_op2_data(out)
            self.add_mass(elem)
            n += 24
        self.card_count['CMASS2'] = nelements
        return n

    def _read_cmass3(self, data, n):
        """
        CMASS3(1201,12,67) - the marker for Record 53
        """
        s = Struct(b(self._endian + '4i'))
        nelements = (len(data) - n) // 16
        for i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CMASS3=%s\n' % str(out))
            #(eid, pid, s1, s2) = out
            elem = CMASS3.add_op2_data(out)
            self.add_mass(elem)
            n += 16
        self.card_count['CMASS3'] = nelements
        return n

    def _read_cmass4(self, data, n):
        """
        CMASS4(1301,13,68) - the marker for Record 54
        """
        nelements = (len(data) - n) // 16
        s = Struct(b(self._endian + 'ifii'))
        for i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = s.unpack(edata)
            #(eid, m,s 1, s2) = out
            elem = CMASS4.add_op2_data(out)
            self.add_mass(elem)
            n += 16
        self.card_count['CMASS4'] = nelements
        return n

    def _read_cmfree(self, data, n):
        """
        CMFREE(2508,25,0) - the marker for Record 55
        """
        self.log.debug('skipping CMFREE in GEOM2\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping CMFREE in GEOM2\n')
        return len(data)

    def _read_conm1(self, data, n):
        """
        CONM1(1401,14,63) - the marker for Record 56
        """
        s = Struct(b(self._endian + '3i21f'))
        nelements = (len(data) - n) // 96
        for i in range(nelements):
            edata = data[n:n + 96]  # 24*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CONM1=%s\n' % str(out))
            (eid, g, cid, m1, m2a, m2b, m3a, m3b, m3c, m4a, m4b, m4c, m4d,
             m5a, m5b, m5c, m5d, m5e, m6a, m6b, m6c, m6d, m6e, m6f) = out
            elem = CONM1.add_op2_data(out)
            self.add_mass(elem)
            n += 96
        self.card_count['CONM1'] = nelements
        return n

    def _read_conm2(self, data, n):
        """
        CONM2(1501,15,64) - the marker for Record 57
        """
        ntotal = 52  # 13*4
        s = Struct(b(self._endian + '3i10f'))
        nelements = (len(data) - n) // ntotal
        for i in range(nelements):
            edata = data[n:n+52]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CONM2=%s\n' % str(out))
            (eid, g, cid, m, x1, x2, x3, i1, i2a, i2b, i3a, i3b, i3c) = out
            elem = CONM2.add_op2_data(out)
            self.add_mass(elem)
            n += ntotal
        self.card_count['CONM2'] = nelements
        return n

    def _read_conrod(self, data, n):
        """
        CONROD(1601,16,47) - the marker for Record 58
        """
        ntotal = 32  # 8*4
        s = Struct(b(self._endian + '4i4f'))
        nelements = (len(data) - n) // ntotal
        for i in range(nelements):
            edata = data[n:n+32]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CONROD=%s\n' % str(out))
            (eid, n1, n2, mid, a, j, c, nsm) = out
            elem = CONROD.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CONROD'] = nelements
        return n

    def _read_conv(self, data, n):
        """
        The CONV card is different between MSC and NX Nastran.
        The MSC version is 8 fields longer.
        """
        n0 = n
        if self.is_nx:
            try:
                n, elements = self._read_conv_nx(data, n)
            except AssertionError:
                n, elements = self._read_conv_msc(data, n0)
        else:
            try:
                n, elements = self._read_conv_msc(data, n)
            except AssertionError:
                n, elements = self._read_conv_nx(data, n0)

        nelements = len(elements)
        for elem in elements:
            self.add_thermal_BC(elem, elem.eid)

        self.card_count['CONV'] = nelements
        return n

    def _read_dual_card(self, data, n, nx_read, msc_read, card_name, add_method):
        """
        generalization of multi read methods (MSC, NX)
        """
        n0 = n
        if self.is_nx:
            try:
                n, elements = nx_read(data, n)
            except AssertionError:
                #raise
                n, elements = msc_read(data, n0)
        else:
            try:
                n, elements = msc_read(data, n)
            except AssertionError:
                #raise
                n, elements = nx_read(data, n0)

        nelements = len(elements)
        for elem in elements:
            add_method(elem)

        self.card_count[card_name] = nelements
        return n

    def _read_conv_nx(self, data, n):
        """
        CONV(12701,127,408) - the marker for Record 59
        """
        ntotal = 48  # 12*4
        s = Struct(b(self._endian + '4i 8i'))
        nelements = (len(data) - n) // ntotal
        elements = []
        for i in range(nelements):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CONV=%s; len=%s\n' % (str(out), len(out)))
            (eid, pcon_id, flmnd, cntrlnd,
             ta1, ta2, ta3, ta4, ta5, ta6, ta7, ta8) = out
            assert eid > 0, out
            #assert eid > 0, out

            #ta = [ta1, ta2, ta3, ta5, ta6, ta7, ta8]
            weights = [None] * 8
            data_in = [eid, pcon_id, flmnd, cntrlnd,
                       [ta1, ta2, ta3, ta4, ta5, ta6, ta7, ta8],
                       weights]

            elem = CONV.add_op2_data(data_in)
            elements.append(elem)
            n += ntotal
        return n, elements

    def _read_conv_msc(self, data, n):
        """
        CONV(12701,127,408) - the marker for Record 60
        """
        ntotal = 80  # 20*4
        s = Struct(b(self._endian + '12i 8f'))
        nelements = (len(data) - n) // ntotal
        elements = []
        for i in range(nelements):
            edata = data[n:n+80]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CONV=%s; len=%s\n' % (str(out), len(out)))
            (eid, pcon_id, flmnd, cntrlnd,
             ta1, ta2, ta3, ta4, ta5, ta6, ta7, ta8,
             wt1, wt2, wt3, wt4, wt5, wt6, wt7, wt8) = out
            assert eid > 0, out
            data_in = [eid, pcon_id, flmnd, cntrlnd,
                       [ta1, ta2, ta3, ta5, ta6, ta7, ta8],
                       [wt1, wt2, wt3, wt5, wt6, wt7, wt8]]
            elem = CONV.add_op2_data(data_in)
            elements.append(elem)
            n += ntotal
        return n, elements

    def _read_convm(self, data, n):
        """
        CONVM(8908,89,422) - the marker for Record 60

        TODO: MSC has 7 fields in QRG (6 defined in DMAP)
              NX has 6 fields in QRG (6 defined in DMAP)
              MSC has extra MDOT field.
              I think it's 6...
        """
        #return len(data)
        ntotal = 24  # 7*4
        s = Struct(b(self._endian + '6i'))
        nelements = (len(data) - n) // ntotal
        for i in range(nelements):
            edata = data[n:n+24]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CONVM=%s\n' % str(out))
            (eid, pconID, flmnd, cntrlnd, ta1, ta2) = out
            assert eid > 0, out  # TODO: I'm not sure that this really has 7 fields...
            mdot = 0.
            data_in = [eid, pconID, flmnd, cntrlnd, ta1, ta2, mdot]
            elem = CONVM.add_op2_data(data_in)
            self.add_thermal_BC(elem, elem.eid)
            n += ntotal
        self.card_count['CONVM'] = nelements
        return n

    def _read_cpyram(self, data, n):
        """
        CPYRAM(17200,172,1000) - the marker for Record ???

        Specific to NX Nastran
        """
        s = Struct(b(self._endian + '16i'))
        nelements = (len(data) - n) // 64
        for i in range(nelements):
            edata = data[n:n + 64]  # 15*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CPENTA=%s\n' % str(out))
            (eid, pid, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10,
             g11, g12, g13, _g14) = out

            data_in = [eid, pid, g1, g2, g3, g4, g5]
            big_nodes = [g6, g7, g8, g9, g10, g11, g12, g13]
            if sum(big_nodes) > 0:
                elem = CPYRAM13.add_op2_data(data_in + big_nodes)
            else:
                elem = CPYRAM5.add_op2_data(data_in)
            self.add_op2_element(elem)
            n += 64
        self.card_count['CPENTA'] = nelements
        return n

# CPENP

    def _read_cpenta(self, data, n):
        """
        CPENTA(4108,41,280)      - the marker for Record 63
        CPENPR(7509,75,9992)     - the marker for Record 64
        CPENT15F(16500,165,9999) - the marker for Record 65
        CPENT6FD(16000,160,9999) - the marker for Record 66
        """
        s = Struct(b(self._endian + '17i'))
        nelements = (len(data) - n) // 68
        for i in range(nelements):
            edata = data[n:n + 68]  # 17*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CPENTA=%s\n' % str(out))
            (eid, pid, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10,
             g11, g12, g13, g14, g15) = out

            data_in = [eid, pid, g1, g2, g3, g4, g5, g6]
            big_nodes = [g7, g8, g9, g10, g11, g12, g13, g14, g15]
            if sum(big_nodes) > 0:
                elem = CPENTA15.add_op2_data(data_in + big_nodes)
            else:
                elem = CPENTA6.add_op2_data(data_in)
            self.add_op2_element(elem)
            n += 68
        self.card_count['CPENTA'] = nelements
        return n

# CQDX4FD
# CQDX9FD - same as CQDX4FD

    def _read_cquad(self, data, n):
        """
        CQUAD(9108,91,507)  - the marker for Record 69
        """
        return self.run_cquad(data, n, CQUAD)

    def run_cquad(self, data, n, element):
        """common method for CQUAD, CQUADX"""
        s = Struct(b(self._endian + '11i'))
        nelements = (len(data) - n) // 44  # 11*4
        if self.is_debug_file:
            self.binary_debug.write('ndata=%s\n' % (nelements * 44))
        for i in range(nelements):
            edata = data[n:n + 44]
            out = s.unpack(edata)
            (eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, n9) = out
            if self.is_debug_file:
                self.binary_debug.write('  %s=%s\n' % (element.type, str(out)))
            #print("eid=%s pid=%s n1=%s n2=%s n3=%s n4=%s theta=%s zoffs=%s tflag=%s t1=%s t2=%s t3=%s t4=%s" % (
                #eid, pid, n1, n2, n3, n4, theta, zoffs, tflag, t1, t2, t3, t4))
            #data_init = [eid, pid, n1, n2, n3, n4, theta, zoffs, tflag, t1, t2, t3, t4]
            data = [eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, n9]
            elem = element.add_op2_data(data)
            self.add_op2_element(elem)
            n += 44
        self.card_count[element.type] = nelements
        return n

    def _read_cquad4(self, data, n):
        """
        CQUAD4(2958,51,177)    - the marker for Record 70
        CQUAD4(13900,139,9989) - the marker for Record 71
        """
        return self.run_cquad4(data, n, CQUAD4)

    def run_cquad4(self, data, n, element):
        """
        common method for CQUAD4, CQUADR
        """
        nelements = (len(data) - n) // 56
        s = Struct(b(self._endian + '6iffii4f'))
        if self.is_debug_file:
            self.binary_debug.write('ndata=%s\n' % (nelements * 44))
        for i in range(nelements):
            edata = data[n:n + 56]  # 14*4
            out = s.unpack(edata)
            (eid, pid, n1, n2, n3, n4, theta, zoffs, blank, tflag,
             t1, t2, t3, t4) = out
            if self.is_debug_file:
                self.binary_debug.write('  %s=%s\n' % (element.type, str(out)))
            #print("eid=%s pid=%s n1=%s n2=%s n3=%s n4=%s theta=%s zoffs=%s blank=%s tflag=%s t1=%s t2=%s t3=%s t4=%s" % (
                #eid, pid, n1, n2, n3, n4, theta, zoffs, blank, tflag, t1, t2, t3, t4))
            data_init = [
                eid, pid, n1, n2, n3, n4, theta, zoffs,
                tflag, t1, t2, t3, t4]
            elem = element.add_op2_data(data_init)
            self.add_op2_element(elem)
            n += 56
        self.card_count[element.type] = nelements
        return n

# CQUAD4FD

    def _read_cquad8(self, data, n):
        """
        CQUAD8(4701,47,326)  - the marker for Record 72
        .. warning:: inconsistent with dmap manual
        """
        #return n
        nelements = (len(data) - n) // 68  # 17*4
        s = Struct(b(self._endian + '10i 6f i'))
        for i in range(nelements):
            edata = data[n:n + 68]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CQUAD8=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, t1, t2,
             t3, t4, theta, zoffs, tflag) = out
            #print("eid=%s pid=%s n1=%s n2=%s n3=%s n4=%s theta=%s zoffs=%s tflag=%s t1=%s t2=%s t3=%s t4=%s" %
                  #(eid, pid, n1, n2, n3, n4, theta, zoffs, tflag, t1, t2, t3, t4))
            #data_init = [eid,pid,n1,n2,n3,n4,theta,zoffs,tflag,t1,t2,t3,t4]
            elem = CQUAD8.add_op2_data(out)
            self.add_op2_element(elem)
            n += 68
        self.card_count['CQUAD8'] = nelements
        return n

# CQUAD9FD
# CQUADP

    def _read_cquadr(self, data, n):
        """
        CQUADR(8009,80,367)  - the marker for Record 75
        """
        return self.run_cquad4(data, n, CQUADR)

    def _read_cquadx(self, data, n):
        """
        CQUADX(9008,90,508)  - the marker for Record 76
        """
        return self.run_cquad4(data, n, CQUADX)

# CRBAR
# CRBE1
# CRBE3
# CRJOINT

    def _read_crod(self, data, n):
        """
        CROD(3001,30,48)    - the marker for Record 81
        """
        s = Struct(b(self._endian + '4i'))
        nelements = (len(data) - n) // 16  # 4*4
        for i in range(nelements):
            edata = data[n:n + 16]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CROD=%s\n' % str(out))
            (eid, pid, n1, n2) = out
            elem = CROD.add_op2_data(out)
            self.add_op2_element(elem)
            n += 16
        self.card_count['CROD'] = nelements
        return n

# CRROD
# CSEAM

    def _read_cshear(self, data, n):
        """
        CSHEAR(3101,31,61)    - the marker for Record 84
        """
        s = Struct(b(self._endian + '6i'))
        nelements = (len(data) - n) // 24  # 6*4
        for i in range(nelements):
            edata = data[n:n + 24]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CSHEAR=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4) = out
            elem = CSHEAR.add_op2_data(out)
            self.add_op2_element(elem)
            n += 24
        self.card_count['CSHEAR'] = nelements
        return n

# CSLOT3
# CSLOT4

    def _read_ctetrap(self, data, n):
        """
        CTETP(12201,122,9013)    - the marker for Record 87
        .. todo:: needs work
        """
        nelements = (len(data) - n) // 108  # 27*4
        s = Struct(b(self._endian + '27i'))
        for i in range(nelements):
            edata = data[n:n+108]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CTETP=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12,
             f1, f2, f3, f4, b1, ee1, ee2, ee3, ee4) = out
            #print("out = ",out)
            e = [e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12]
            f = [f1, f2, f3, f4]
            ee = [ee1, ee2, ee3, ee4]

            #print("e  = ",e)
            #print("f  = ",f)
            #print("b1  = ",b1)
            #print("ee = ",ee)
            data_in = [eid, pid, n1, n2, n2, n3, n4]
            elem = CTETRA4.add_op2_data(data_in)
            self.add_op2_element(elem)

    def _read_ctetra(self, data, n):
        """
        CTETRA(5508,55,217)      - the marker for Record 88
        CTETPR(7609,76,9993)     - the marker for Record 89
        CTETR10F(16600,166,9999) - the marker for Record 90
        CTETR4FD(16100,161,9999) - the marker for Record 91
        """
        s = Struct(b(self._endian + '12i'))
        nelements = (len(data) - n)// 48  # 12*4
        for i in range(nelements):
            edata = data[n:n + 48]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CTETRA=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10) = out
            #print("out = ",out)

            data_in = [eid, pid, n1, n2, n3, n4]
            big_nodes = [n5, n6, n7, n8, n9, n10]
            if sum(big_nodes) > 0:
                elem = CTETRA10.add_op2_data(data_in + big_nodes)
            else:
                elem = CTETRA4.add_op2_data(data_in)
            self.add_op2_element(elem)
            n += 48
        self.card_count['CTETRA'] = nelements
        return n

# CTQUAD - 92
# CTTRIA - 93

    def _read_ctria3(self, data, n):
        """
        CTRIA3(5959,59,282)    - the marker for Record 94
        """
        ntotal = 52  # 13*4
        s = Struct(b(self._endian + '5iff3i3f'))
        nelements = (len(data) - n)// 52  # 13*4
        for i in range(nelements):
            edata = data[n:n+52]
            out = s.unpack(edata)
            #print("eid=%s pid=%s n1=%s n2=%s n3=%s theta=%s zoffs=%s blank1=%s blank2=%s tflag=%s t1=%s t2=%s t3=%s" %
                  #(eid, pid, n1, n2, n3, theta, zoffs, blank1, blank2, tflag, t1, t2, t3))
            (eid, pid, n1, n2, n3, theta, zoffs, blank1,
             blank2, tflag, t1, t2, t3) = out
            if self.is_debug_file:
                self.binary_debug.write('  CTRIA3=%s\n' % str(out))
            data_in = [eid, pid, n1, n2, n3, theta, zoffs, tflag, t1, t2, t3]
            elem = CTRIA3.add_op2_data(data_in)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CTRIA3'] = nelements
        return n


# CTRIAFD - 95

    def _read_ctria6(self, data, n):
        """
        CTRIA6(4801,48,327)    - the marker for Record 96
        .. warning:: inconsistent with dmap manual
        """
        s = Struct(b(self._endian + '8i 5f i'))
        nelements = (len(data) - n) // 56  # 14*4
        for i in range(nelements):
            edata = data[n:n + 56]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CTRIA6=%s\n' % str(out))
            #print("eid=%s pid=%s n1=%s n2=%s n3=%s theta=%s zoffs=%s blank1=%s blank2=%s tflag=%s t1=%s t2=%s t3=%s" %
                  #(eid, pid, n1, n2, n3, theta, zoffs, blank1, blank2, tflag, t1, t2, t3))
            (eid, pid, n1, n2, n3, n4, n5, n6, theta, zoffs, t1, t2, t3, tflag) = out
            elem = CTRIA6.add_op2_data(out)
            self.add_op2_element(elem)
            n += 56
        self.card_count['CTRIA6'] = nelements
        return n

# CTRIA6FD
# CTRIAP

    def _read_ctriar(self, data, n):  # 98
        """
        CTRIAR(9200,92,385)    - the marker for Record 99
        """
        ntotal = 52  # 13*4
        s = Struct(b(self._endian + '5iff3i3f'))
        nelements = (len(data) - n)// 52  # 13*4
        for i in range(nelements):
            edata = data[n:n+52]
            out = s.unpack(edata)
            #print("eid=%s pid=%s n1=%s n2=%s n3=%s theta=%s zoffs=%s blank1=%s blank2=%s tflag=%s t1=%s t2=%s t3=%s" %
                  #(eid, pid, n1, n2, n3, theta, zoffs, blank1, blank2, tflag, t1, t2, t3))
            (eid, pid, n1, n2, n3, theta, zoffs, blank1,
             blank2, tflag, t1, t2, t3) = out
            if self.is_debug_file:
                self.binary_debug.write('  CTRIAR=%s\n' % str(out))
            data_in = [eid, pid, n1, n2, n3, theta, zoffs, tflag, t1, t2, t3]
            elem = CTRIA3.add_op2_data(data_in)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CTRIA3'] = nelements
        return n

# CTRIAX - 100

    def _read_ctriax6(self, data, n):  # 101
        self.log.debug('skipping CTRIAX6 in GEOM2\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping CTRIAX6 in GEOM2\n')
        return len(data)

# CTRIX3FD - 102
# CTRIX6FD - 103

    def _read_ctube(self, data, n):
        """
        CTUBE(3701,37,49)    - the marker for Record 104
        """
        s = Struct(b(self._endian + '4i'))
        nelements = (len(data) - n) // 16
        for i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CTUBE=%s\n' % str(out))
            (eid, pid, n1, n2) = out
            elem = CTUBE.add_op2_data(out)
            self.add_op2_element(elem)
            n += 16
        self.card_count['CTUBE'] = nelements
        return n

    def _read_cvisc(self, data, n):
        """CVISC(3901,39,50) - the marker for Record 105"""
        s = Struct(b(self._endian + '4i'))
        nelements = (len(data) - n) // 16
        for i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CVISC=%s\n' % str(out))
            #(eid,pid,n1,n2) = out
            elem = CVISC.add_op2_data(out)
            self.add_op2_element(elem)
            n += 16
        self.card_count['CVISC'] = nelements
        return n

    def _read_cweld(self, data, n):
        """
        CWELD(11701,117,559) - Record 106
        same as CFAST
        """
        self.log.debug('skipping CWELD in GEOM2\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping CWELD in GEOM2\n')
        return len(data)

    def _read_cweldc(self, data, n):  # 107
        self.log.debug('skipping CWELDC in GEOM2\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping CWELDC in GEOM2\n')
        return len(data)

    def _read_cweldg(self, data, n):  # 108
        self.log.debug('skipping CWELDG in GEOM2\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping CWELDG in GEOM2\n')
        return len(data)

# TDOO: above are checked by DMAP...
#-------------------------------
# CWSEAM
# GENEL
# GMDNDC
# GMBNDS
# GMINTC
# GMINTS
    def _read_plotel(self, data, n):  # 114
        s = Struct(b(self._endian + '3i'))
        ntotal = 12
        nelements = (len(data) - n) // ntotal
        for i in range(nelements):
            edata = data[n:n + ntotal]  # 4*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  PLOTEL=%s\n' % str(out))
            #(eid,n1,n2) = out
            elem = PLOTEL.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['PLOTEL'] = nelements
        return n
# RADBC
# RADINT
# SINT

    def add_spoint(self, spooint):
        raise RuntimeError('this should be overwritten by the BDF')

    def _read_spoint(self, data, n):
        """
        (5551,49,105)    - the marker for Record 118
        """
        npoints = (len(data) - n) // 4
        fmt = b(self._endian + '%ii' % npoints)
        nids = unpack(fmt, data[n:])
        if self.is_debug_file:
            self.binary_debug.write('SPOINT=%s\n' % str(nids))
        spoint = SPOINTs.add_op2_data(list(nids))
        self.add_spoint(spoint)
        self.card_count['SPOINT'] = npoints
        return n

    def _read_vubeam(self, data, n):  # 119
        self.log.debug('skipping VUBEAM in GEOM2\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping VUBEAM in GEOM2\n')
        return len(data)

# VUHEXA
# VUQUAD4
# VUPENTA
# VUTETRA
# VUTRIA
# VUBEAM
# VUHEXA
# VUQUAD4
# WELDP
