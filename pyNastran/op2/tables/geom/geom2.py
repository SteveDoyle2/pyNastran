"""
defines readers for BDF objects in the OP2 GEOM2/GEOM2S table
"""
# pylint: disable=C0103
from struct import Struct
import numpy as np

from pyNastran.bdf.cards.elements.elements import CGAP, PLOTEL
from pyNastran.bdf.cards.elements.damper import (CDAMP1, CDAMP2, CDAMP3,
                                                 CDAMP4, CDAMP5, CVISC)
from pyNastran.bdf.cards.elements.springs import CELAS1, CELAS2, CELAS3, CELAS4
from pyNastran.bdf.cards.elements.axisymmetric_shells import CQUADX, CTRIAX6
from pyNastran.bdf.cards.elements.shell import (CTRIA3, CQUAD4, CTRIA6,
                                                CQUADR, CTRIAR,
                                                CQUAD8, CQUAD,
                                                CSHEAR)
from pyNastran.bdf.cards.elements.rods import CROD, CTUBE, CONROD
from pyNastran.bdf.cards.elements.bars import CBAR, CBEND
from pyNastran.bdf.cards.elements.beam import CBEAM
from pyNastran.bdf.cards.elements.mass import (CONM1, CONM2, CMASS1, CMASS2,
                                               CMASS3, CMASS4)
from pyNastran.bdf.cards.elements.solid import (CTETRA4, CPYRAM5, CPENTA6, CHEXA8,
                                                CTETRA10, CPYRAM13, CPENTA15, CHEXA20,)
from pyNastran.bdf.cards.thermal.thermal import CHBDYG, CONV, CHBDYP, CHBDYE, CONVM
from pyNastran.bdf.cards.nodes import SPOINTs
from pyNastran.bdf.cards.elements.bush import CBUSH

from pyNastran.op2.errors import MixedVersionCard
from pyNastran.op2.tables.geom.geom_common import GeomCommon


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

            # MSC
            (8515, 85, 0): ['CFLUID2', self._read_cfluid2],  # record 36 - not done
            (8615, 86, 0): ['CFLUID3', self._read_cfluid3],  # record 37 - not done
            (8715, 87, 0): ['CFLUID4', self._read_cfluid4],  # record 38 - not done
            (7701, 77, 8881): ['CINT', self._read_fake],     # record 39 - not done
            (1908, 19, 104): ['CGAP', self._read_cgap],       # record 40 - buggy
            (8100, 81, 381): ['CHACAB', self._read_fake],     # record 41 - not done
            (8200, 82, 383): ['CHACBR', self._read_fake],     # record 42 - not done
            (8308, 83, 405): ['CHBDYE', self._read_fake],     # record 43 - not done
            (10808, 108, 406): ['CHBDYG', self._read_chbdyg], # record 44
            (10908, 109, 407): ['CHBDYP', self._read_chbdyp], # record 45
            (7308, 73, 253): ['CHEXA', self._read_chexa],     # record 46
            # CHEXA20F record 47
            # CHEXAFD record 48
            # CHEXAL record 49
            (12001, 120, 9011): ['CHEXP', self._read_chexp],    # record 50
            # CHEXPR record 51
            (1001, 10, 65): ['CMASS1', self._read_cmass1],    # record 52
            (1101, 11, 66): ['CMASS2', self._read_cmass2],    # record 53
            (1201, 12, 67): ['CMASS3', self._read_cmass3],    # record 54
            (1301, 13, 68): ['CMASS4', self._read_cmass4],    # record 55
            (2508, 25, 0): ['CMFREE', self._read_cmfree],     # record 56 - not done
            (1401, 14, 63): ['CONM1', self._read_conm1],      # record 57
            (1501, 15, 64): ['CONM2', self._read_conm2],      # record 58
            (1601, 16, 47): ['CONROD', self._read_conrod],    # record 59
            (12701, 127, 408): ['CONV', self._read_conv],     # record 60 - not tested
            (8908, 89, 422): ['CONVM', self._read_convm],     # record 61 - not tested
            (12101, 121, 9012) : ['CPENPR', self._read_fake],  # record 62
            (4108, 41, 280): ['CPENTA', self._read_cpenta],   # record 63
            (7509, 75, 9992) : ['CPENPR', self._read_fake],  # record 64
            (16500, 165, 9999) : ['CPENT15F', self._read_fake],  # record 65
            (16000, 160, 9999) : ['CPENT6FD', self._read_fake],  # record 66
            (17000, 170, 9999) : ['CQDX4FD', self._read_fake],  # record 67
            (17100, 171, 9999) : ['CQDX9FD', self._read_fake],  # record 68
            (9108, 91, 507): ['CQUAD', self._read_cquad],       # record 69 - not tested
            (2958, 51, 177): ['CQUAD4', self._read_cquad4],     # record 70
            (13900, 139, 9989): ['CQUAD4', self._read_cquad4],  # record 71
            (4701, 47, 326): ['CQUAD8', self._read_cquad8],     # record 72
            (16400, 164, 9999) : ['CQUAD9FD', self._read_fake],  # record 73
            (11101, 111, 9014) : ['CQUADP', self._read_fake],  # record 74
            (8009, 80, 367): ['CQUADR', self._read_cquadr],   # record 75 - not tested
            (9008, 90, 508): ['CQUADX', self._read_cquadx],   # record 76 - not tested
            (14700, 147, 6662) : ['CRBAR', self._read_fake],  # record 77
            (17300, 173, 6664) : ['CRBE1', self._read_fake],  # record 78
            (17200, 172, 6663) : ['CRBE3', self._read_fake],  # record 79
            (11000, 110, 6667) : ['CRJOINT', self._read_fake],  # record 80
            (3001, 30, 48): ['CROD', self._read_crod],        # record 81
            (12600, 126, 6661) : ['CRROD', self._read_fake],  #  record 82
            (13801, 138, 570) : ['CSEAM', self._read_fake],  #  record 83
            (3101, 31, 61): ['CSHEAR', self._read_cshear],        # record 84
            (4408, 44, 227) : ['CSLOT3', self._read_fake],  # record 85
            (4508, 45, 228) : ['CSLOT4', self._read_fake],  # record 86
            (12201, 122, 9013): ['CTETP', self._read_ctetrap],  # record 87 - not done
            (5508, 55, 217): ['CTETRA', self._read_ctetra],     # record 88
            (7609, 76, 9993) : ['CTETPR', self._read_fake],  # record 89
            (16600, 166, 9999) : ['CTETR10F', self._read_fake],  # record 90
            (16100, 161, 9999) : ['CTETR4FD', self._read_fake],  # record 91
            (14801, 148, 643) : ['CTQUAD', self._read_fake],  # record 92
            (14901, 149, 644) : ['CTTRIA', self._read_fake],  # record 93
            (5959, 59, 282): ['CTRIA3', self._read_ctria3],   # record 94
            (16200, 162, 9999) : ['CTRIA3FD', self._read_fake],  # record 95
            (4801, 48, 327): ['CTRIA6', self._read_ctria6],   # record 96 - buggy
            # : ['RADBC', self._read_fake],  record 97
            # : ['RADBC', self._read_fake],  record 98
            (9200, 92, 385): ['CTRIAR', self._read_ctriar],   # record 99  - not done
            # : ['RADBC', self._read_fake],  record 100
            (6108, 61, 107): ['CTRIAX6', self._read_ctriax6], # record 101 - not done
            (16700, 167, 9999) : ['CTRIA6FD', self._read_fake],  # record 102
            (11301, 113, 9015) : ['CTRIAP', self._read_fake],  # record 103
            (3701, 37, 49): ['CTUBE', self._read_ctube],      # record 104
            (3901, 39, 50): ['CVISC', self._read_cvisc],      # record 105 - not done
            (11701, 117, 559) : ['CWELD', self._read_fake],  # record 106
            (13501, 135, 564) : ['CWELDC', self._read_fake],  # record 107
            (13601, 136, 562) : ['CWELDG', self._read_fake],  # record 108
            # : ['RADBC', self._read_fake],  record 109
            # : ['RADBC', self._read_fake],  record 110
            # : ['RADBC', self._read_fake],  record 111
            # : ['RADBC', self._read_fake],  record 112
            # : ['RADBC', self._read_fake],  record 113
            # : ['RADBC', self._read_fake],  record 114
            (5201, 52, 11) :   ['PLOTEL', self._read_plotel],    # record 115 - not done
            (12801, 128, 417) : ['RADBC', self._read_fake],  # record 116
            (15501, 155, 634) : ['RADINT', self._read_fake], # record 117
            (7801, 78, 8883) : ['SINT', self._read_fake],    # record 118
            (5551, 49, 105) : ['SPOINT', self._read_spoint],     # record 119
            (11601, 116, 9942) : ['VUBEAM', self._read_vubeam],  # record 120 - not done
            # record 121
            # record 121
            # record 121
            # record 121
            # record 121
            # record 121
            # record 121
            (2108, 21, 224): ['CAXIF2', self._read_fake],
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
            (12801, 128, 417): ['RADBC', self._read_radbc],
            (2708, 27, 59): ['CAABSF', self._read_fake],
            (3201, 32, 478): ['GMBNDC', self._read_fake],
            (13900, 139, 9984): ['CQUAD4FD', self._read_fake],
            (14000, 140, 9990): ['CHEXAFD', self._read_fake],
            (16000, 160, 9988): ['CPENTA6FD', self._read_fake],
            (16100, 161, 9986): ['CTETRAFD', self._read_fake],
            (16300, 163, 9989): ['CHEXA20F', self._read_fake],
            (16700, 167, 9981): ['CTRIA6FD', self._read_fake],
            (16800, 168, 9978): ['CTRIAX3FD', self._read_fake],
            (16500, 165, 9987): ['CPENT15F', self._read_fake],
            (5008, 50, 258): ['', self._read_fake],
            (16400, 164, 9983) : ['CQUAD9FD', self._read_fake],
            (11000, 110, 6667): ['', self._read_fake],
            (12301, 123, 9921): ['', self._read_fake],
            (12401, 124, 9922): ['FEFACE/PVAL?', self._read_fake],
            (12600, 126, 6661): ['', self._read_fake],
            (14700, 147, 6662): ['', self._read_fake],
            (7309, 73, 0): ['', self._read_fake],
            (11501, 115, 9941): ['', self._read_fake],    # record
            (12501, 125, 9923): ['', self._read_fake],    # record
            (3401, 34, 9600): ['', self._read_fake],    # record
            (7701, 77, 8881): ['', self._read_fake],  # record
            (2901, 29, 9601): ['', self._read_fake],  # record
            (16600, 166, 9985) : ['', self._read_fake],  # record
            (16200, 162, 9982) : ['', self._read_fake],  # record
            (16900, 169, 9977) : ['', self._read_fake],  # record
            (1701, 17, 980) : ['CPLSTN3', self._read_fake],  # record
            (1801, 18, 986) : ['', self._read_fake],  # record
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
            (6111, 61, 996) : ['CTRAX3', self._read_trax3],  # record
            (6112, 61, 997) : ['CQUADX4', self._read_cquadx4],  # record
            (6113, 61, 998) : ['CTRAX6', self._read_ctrax6],  # record
            (6114, 61, 999) : ['CQUADX8', self._read_cquadx8],  # record
            (3501, 35, 1) : ['', self._read_fake],  # record
            (1001, 100, 10000) : ['', self._read_fake],  # record
            (1118, 1, 1874) : ['', self._read_fake],  # record

            # NX specific
            (17200, 172, 1000) : ['CPYRAM', self._read_cpyram],
            (25700, 257, 9948) : ['CPYRA5FD', self._read_cpyram],
            (25800, 258, 9947) : ['CPYRA13F', self._read_cpyram],
            (7909, 79, 9946) : ['CPYRAMPR', self._read_cpyram],
            (7201, 72, 983) : ['CPLSTN8', self._read_fake],
            (11701, 117, 559) : ['CWELD', self._read_fake],
            (13501, 135, 564) : ['CWELDC', self._read_fake],
            (3601, 36, 987) : ['CPLSTS8', self._read_fake],
            (13701, 137, 565) : ['CWELDP', self._read_fake],
        }

    def add_op2_element(self, elem):
        """checks that eids are positive and that -1 node ids become None"""
        if elem.eid <= 0:
            self.log.debug(elem)
            raise ValueError(elem)
            #return

        #if elem.eid > 100000000:
            #raise RuntimeError('bad parsing...elem:\n%s' % elem)

        if elem.type in ['CTRIA6', 'CQUAD8']:
            for nid in elem.nodes:
                if nid == -1:
                    nid = None
        else:
            for nid in elem.nodes:
                if nid == -1:
                    assert nid > 0, elem
        self._add_element_object(elem, allow_overwrites=False)
        #print(str(elem)[:-1])

# 1-AEROQ4 (???)
# AEROT3   (???)
# 1-BEAMAERO (1701,17,0)
# 2-CAABSF (2708,27,59)
# 3-CAXIF2 (2108,21,224)
# 4-CAXIF3 (2208,22,225)
# 5-CAXIF4 (2308,23,226)

    def _read_cbar(self, data, n):
        r"""
        CBAR(2408,24,180) - the marker for Record 8

        MSC/NX
        Word Name Type Description
        1 EID    I  Element identification number
        2 PID    I  Property identification number
        3 GA     I  Grid point identification number at end A
        4 GB     I  Grid point identification number at end B

        F=0* XYZ option -- basic coordinate system
           5 X1 RS  T1 component of orientation vector from GA
           6 X2 RS  T2 component of orientation vector from GA
           7 X3 RS  T3 component of orientation vector from GA
           8 FE  I  Orientation vector flag (encoded)
        F=1* XYZ option -- global coordinate system
           5 X1 RS  T1 component of orientation vector from GA
           6 X2 RS  T2 component of orientation vector from GA
           7 X3 RS  T3 component of orientation vector from GA
           8 FE  I   Orientation vector flag (encoded)
        F=2* Grid option
           5 GO I Grid point identification number at end of orientation vector
           6 UNDEF(2) none Not used
           8 FE I Orientation vector flag (encoded)
        *F = FE bit-wise AND with 3
        End F

        9  PA   I Pin flags for end A
        10 PB   I Pin flags for end B
        11 W1A RS T1 component of offset vector from GA
        12 W2A RS T2 component of offset vector from GA
        13 W3A RS T3 component of offset vector from GA
        14 W1B RS T1 component of offset vector from GB
        15 W2B RS T2 component of offset vector from GB
        16 W3B RS T3 component of offset vector from GB
        F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_sebload1.op2
        """
        nelements = (len(data) - n) // 64
        s1 = Struct(self._endian + b'4i3f3i6f')
        #s2 = Struct(self._endian + b'4i3f3i6f')
        s2 = s1
        s3 = Struct(self._endian + b'7ii2i6f')
        for unused_i in range(nelements):
            edata = data[n:n + 64]  # 16*4
            fe, = self.struct_i.unpack(edata[28:32])
            # per DMAP: F = FE bit-wise AND with 3
            f = fe & 3
            if f == 0:
                out = s1.unpack(edata)
                (eid, pid, ga, gb, x1, x2, x3, _f, pa, pb,
                 w1a, w2a, w3a, w1b, w2b, w3b) = out
                data_in = [[eid, pid, ga, gb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                           [f, x1, x2, x3]]
            elif f == 1:
                out = s2.unpack(edata)
                (eid, pid, ga, gb, x1, x2, x3, _f, pa, pb,
                 w1a, w2a, w3a, w1b, w2b, w3b) = out
                data_in = [[eid, pid, ga, gb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                           [f, x1, x2, x3]]
            elif f == 2:
                out = s3.unpack(edata)
                (eid, pid, ga, gb, g0, unused_junk1, unused_junk2, _f, pa,
                 pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                data_in = [[eid, pid, ga, gb, pa, pb, w1a,
                            w2a, w3a, w1b, w2b, w3b], [f, g0]]
            else:
                raise RuntimeError('invalid f value...f=%s' % (f))
            elem = CBAR.add_op2_data(data_in)
            assert f == fe, 'f=%s type(f)=%s fe=%s\n%s' % (f, type(f), fe, elem)

            self.add_op2_element(elem)
            n += 64
        self.card_count['CBAR'] = nelements
        return n

    def _read_cbarao(self, data, n):
        """
        CBARAO(4001,40,275) - the marker for Record 9

        1 EID   I Element identification number
        2 SCALE I Scale of Xi values
        3 X1 RS 1st intermediate station for data recovery
        4 X2 RS 2nd intermediate station for data recovery
        5 X3 RS 3rd intermediate station for data recovery
        6 X4 RS 4th intermediate station for data recovery
        7 X5 RS 5th intermediate station for data recovery
        8 X6 RS 6th intermediate station for data recovery
        9 UNDEF none Not used
        """
        #self.log.info('skipping CBARAO in GEOM2')
        #if self.is_debug_file:
            #self.binary_debug.write('skipping CBARAO in GEOM2\n')
        #return len(data)
        nelements = (len(data) - n) // 36
        s = Struct(self._endian + b'2i7f')
        for unused_i in range(nelements):
            edata = data[n:n + 36]  # 9*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CBARAO=%s\n' % str(out))
            (eid, scale, x1, x2, x3, x4, x5, x6, unused_null) = out
            if scale == 2:
                scale = 'FR'
            else:
                NotImplementedError('CBARAO scale=%r; 2=FR' % scale)
            x = [x1, x2, x3, x4, x5, x6]
            self.add_cbarao(eid, scale, x, comment='')
            n += 36
        self.card_count['CBARAO'] = nelements
        return n

    def _read_cbeam(self, data, n):
        """
        CBEAM(5408,54,261) - the marker for Record 10
        """
        nelements = (len(data) - n) // 72
        s1 = Struct(self._endian + b'6i3f3i6f')
        s3 = Struct(self._endian + b'12i6f')
        for unused_i in range(nelements):
            edata = data[n:n + 72]  # 18*4
            fe, = self.struct_i.unpack(edata[40:44])

            # per DMAP: F = FE bit-wise AND with 3
            f = fe & 3
            if f == 0:  # basic cid
                out = s1.unpack(edata)
                (eid, pid, ga, gb, sa, sb, x1, x2, x3, fe,
                 pa, pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                #self.log.info('CBEAM: eid=%s fe=%s f=%s; basic cid' % (eid, fe, f))

                data_in = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                           [f, x1, x2, x3]]
            elif f == 1:  # global cid
                # CBEAM    89616   5       384720  384521  0.      0.     -1.
                out = s1.unpack(edata)
                (eid, pid, ga, gb, sa, sb, x1, x2, x3, fe,
                 pa, pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                #self.log.info('CBEAM: eid=%s fe=%s f=%s; global cid' % (eid, fe, f))
                data_in = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                           [f, x1, x2, x3]]
            elif f == 2:  # grid option
                out = s3.unpack(edata)
                (eid, pid, ga, gb, sa, sb, g0, xxa, xxb, fe,
                 pa, pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                #self.log.info('CBEAM: eid=%s fe=%s f=%s; grid option '
                              #'(g0=%s xxa=%s xxb=%s)' % (eid, fe, f, g0, xxa, xxb))
                if g0 <= 0 or g0 >= 100000000 or xxa != 0 or xxb != 0:
                    # Nastran set this wrong...MasterModelTaxi
                    #CBEAM    621614  2672    900380  900379 .197266 -.978394.0600586
                    #    6


                    f = 1
                    out = s1.unpack(edata)
                    (eid, pid, ga, gb, sa, sb, x1, x2, x3, fe, pa,
                     pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                    #self.log.info('CBEAM: eid=%s fe=%s f=%s; global cid' % (eid, fe, f))
                    data_in = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                               [f, x1, x2, x3]]
                    #self.log.info('   (x1=%s x2=%s x3=%s)' % (x1, x2, x3))
                else:
                    data_in = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                               [f, g0]]
            else:
                raise RuntimeError('invalid f value...f=%r' % f)
            if self.is_debug_file:
                self.binary_debug.write('  CBEAM eid=%s f=%s fe=%s %s\n' % (
                    eid, f, fe, str(data_in)))

            elem = CBEAM.add_op2_data(data_in, f)
            self.add_op2_element(elem)
            n += 72
        self.card_count['CBEAM'] = nelements
        return n

    def _read_cbeamp(self, data, n):
        """
        CBEAMP(11401,114,9016) - the marker for Record 11
        """
        self.log.info('skipping CBEAMP in GEOM2')
        return len(data)

    def _read_cbend(self, data, n):
        """
        CBEND(4601,46,298) - the marker for Record 12

        1 EID I Element identification number
        2 PID I Property identification number
        3 GA  I Grid point End A identification number
        4 GB  I Grid point End B identification number

        F = 0 Z
          5 X1 RS T1 component of orientation vector from GA
          6 X2 RS T2 component of orientation vector from GA
          7 X3 RS T3 component of orientation vector from GA
        8 F    I     Orientation vector flag = 0
        F = 1 XYZ option - global cooridnate system
          5 X1 RS T1 component of orientation vector from GA
          6 X2 RS T2 component of orientation vector from GA
          7 X3 RS T3 component of orientation vector from GA
          8 F   I    Orientation vector flag = 1
        F = 2 Grid option
          5 GO       I Grid point ID at end of orientation vector
          6 UNDEF(2)    None
          8 F        I Orientation vector flag = 2
        End F
        9 UNDEF(4) None
        13 GEOM I Element geometry option
        """
        ntotal = 52 # 4*13
        nentries = (len(data) - n) // ntotal
        fstruc = Struct(self._endian + b'4i 3f 6i')
        istruc = Struct(self._endian + b'4i 3i 6i')

        for unused_i in range(nentries):
            edata = data[n:n + 52]  # 13*4
            fe, = self.struct_i.unpack(edata[28:32])
            # per DMAP: F = FE bit-wise AND with 3
            f = fe & 3
            if f == 0:
                out = fstruc.unpack(edata)
                (eid, pid, ga, gb, x1, x2, x3, fe,
                 unused_dunnoa, unused_dunnob, unused_dunnoc, unused_dunnod, geom) = out
                data_in = [[eid, pid, ga, gb, geom],
                           [f, x1, x2, x3]]
            elif f == 1:
                out = fstruc.unpack(edata)
                (eid, pid, ga, gb, x1, x2, x3, fe,
                 unused_dunnoa, unused_dunnob, unused_dunnoc, unused_dunnod, geom) = out
                data_in = [[eid, pid, ga, gb, geom],
                           [f, x1, x2, x3]]
            elif f == 2:
                out = istruc.unpack(edata)
                (eid, pid, ga, gb, g0, unused_junk1, unused_junk2, fe,
                 unused_dunnoa, unused_dunnob, unused_dunnoc, unused_dunnod, geom) = out
                data_in = [[eid, pid, ga, gb, geom],
                           [f, g0]]
            else:
                raise RuntimeError('invalid f value...f=%s' % (f))
            elem = CBEND.add_op2_data(data_in)
            elem.validate()
            assert f == fe, 'f=%s type(f)=%s fe=%s\n%s' % (f, type(f), fe, elem)

            self.add_op2_element(elem)
            n += 52
        self.increase_card_count('CBEND', nentries)
        return n

    def _read_cbush(self, data, n):
        """
        CBUSH(2608,26,60) - the marker for Record 13
        """
        nelements = (len(data) - n) // 56
        struct_obj1 = Struct(self._endian + b'4i iii i ifi3f')
        struct_obj2 = Struct(self._endian + b'4i fff i ifi3f')
        for unused_i in range(nelements):
            edata = data[n:n + 56]  # 14*4
            out = struct_obj1.unpack(edata)
            eid, pid, ga, gb, five, unused_sixi, unused_seven, f, cid, s, ocid, s1, s2, s3 = out
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

        1 EID  I Element identification number
        2 PID  I Property identification number
        3 G(2) I Grid point identification numbers
        5 CID  I Coordinate system identification number
        6 UNDEF(3) none

        """
        ntotal = 32 # 4*8
        nelements = (len(data) - n) // ntotal
        struct_6i = Struct(self._endian + b'8i')
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_6i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CBUSH1D=%s\n' % str(out))
            (eid, pid, g1, g2, cid, unused_a, unused_b, unused_c) = out
            if cid == -1:
                cid = None
            self.add_cbush1d(eid, pid, [g1, g2], cid=cid)
            n += ntotal
        self.card_count['CBUSH1D'] = nelements
        return n
        #self.log.info('skipping CBUSH1D in GEOM2')
        #if self.is_debug_file:
            #self.binary_debug.write('skipping CBUSH1D in GEOM2\n')
        #return len(data)

    def _read_ccone(self, data, n):
        """
        CCONE(2315,23,0) - the marker for Record 15
        """
        self.log.info('skipping CCONE in GEOM2')
        if self.is_debug_file:
            self.binary_debug.write('skipping CCONE in GEOM2\n')
        return len(data)

    def _read_cdamp1(self, data, n):
        """
        CDAMP1(201,2,69) - the marker for Record 16
        """
        nelements = (len(data) - n) // 24
        struct_6i = Struct(self._endian + b'6i')
        for unused_i in range(nelements):
            edata = data[n:n + 24]  # 6*4
            out = struct_6i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CDAMP1=%s\n' % str(out))
            #(eid, pid, g1, g2, c1, c2) = out
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
        s = Struct(self._endian + b'if4i')
        for unused_i in range(nelements):
            edata = data[n:n + 24]  # 6*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CDAMP2=%s\n' % str(out))
            #(eid, bdamp, g1, g2, c1, c2) = out
            elem = CDAMP2.add_op2_data(out)
            self.add_op2_element(elem)
            n += 24
        self.card_count['CDAMP2'] = nelements
        return n

    def _read_cdamp3(self, data, n):
        """
        CDAMP3(401,4,71) - the marker for Record 18
        """
        struct_4i = Struct(self._endian + b'4i')
        nelements = (len(data) - n) // 16
        for unused_i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = struct_4i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CDAMP3=%s\n' % str(out))
            #(eid, pid, s1, s2) = out
            elem = CDAMP3.add_op2_data(out)
            self.add_op2_element(elem)
            n += 16
        self.card_count['CDAMP3'] = nelements
        return n

    def _read_cdamp4(self, data, n):
        """
        CDAMP4(501,5,72) - the marker for Record 19
        """
        s = Struct(self._endian + b'ifii')
        nelements = (len(data) - n) // 16
        for unused_i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CDAMP4=%s\n' % str(out))
            #(eid, bdamp, s1, s2) = out
            elem = CDAMP4.add_op2_data(out)
            self.add_op2_element(elem)
            n += 16
        self.card_count['CDAMP4'] = nelements
        return n

    def _read_cdamp5(self, data, n):
        """
        CDAMP5(10608,106,404) - the marker for Record 20
        """
        s = Struct(self._endian + b'4i')
        nelements = (len(data) - n) // 16
        for unused_i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CDAMP5=%s\n' % str(out))
            #(eid, pid, s1, s2) = out
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
        struct_4i = Struct(self._endian + b'6i')
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n+24]
            out = struct_4i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CELAS1=%s\n' % str(out))
            #(eid, pid, g1, g2, c1, c2) = out
            elem = CELAS1.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CELAS1'] = nelements
        return n

    def _read_celas2(self, data, n):
        """
        CELAS2(701,7,74) - the marker for Record 30
        """
        s1 = Struct(self._endian + b'if4iff')
        ntotal = 32
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n+32]
            out = s1.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CELAS2=%s\n' % str(out))
            #(eid, k, g1, g2, c1, c2, ge, s) = out
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
        struct_4i = Struct(self._endian + b'4i')
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n+16]
            out = struct_4i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CELAS3=%s\n' % str(out))
            #(eid, pid, s1, s2) = out
            elem = CELAS3.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CELAS3'] = nelements
        return n

    def _read_celas4(self, data, n):
        """
        CELAS4(901,9,76) - the marker for Record 32
        """
        s = Struct(self._endian + b'ifii')
        nelements = (len(data) - n) // 16
        for unused_i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CELAS4=%s\n' % str(out))
            #(eid, k, s1, s2) = out
            elem = CELAS4.add_op2_data(out)
            self.add_op2_element(elem)
            n += 16
        self.card_count['CELAS4'] = nelements
        return n

    def _read_cfast(self, data, n):
        """
        CFAST(9801,98,506) - the marker for Record ???
        """
        self.log.info('skipping CFAST in GEOM2')
        return len(data)

# CFASTP

    def _read_cfluid2(self, data, n):
        """
        CFLUID2(8515,85,209) - the marker for Record 35

        1 EID       I Element identification number
        2 IDF1      I RINGFL point 1 identification number
        3 IDF2      I RINGFL point 2 identification number
        4 RHO      RS Mass density
        5 B        RS Bulk modulus
        6 HARMINDX  I Harmonic index
        """
        s = Struct(self._endian + b'3i2fi')
        nelements = (len(data) - n) // 24
        for unused_i in range(nelements):
            edata = data[n:n + 24]  # 6*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CFLUID2=%s\n' % str(out))
            eid, idf1, idf2, rho, bi, harmonic = out
            self.add_cfluid2(eid, [idf1, idf2], rho, bi, harmonic)
            n += 24
        self.card_count['CFLUID2'] = nelements
        return n

    def _read_cfluid3(self, data, n):
        """
        CFLUID3(8615,86,210) - the marker for Record 36

        1 EID       I Element identification number
        2 IDF1      I RINGFL point 1 identification number
        3 IDF2      I RINGFL point 2 identification number
        4 IDF3      I RINGFL point 3 identification number
        5 RHO      RS Mass density
        6 B        RS Bulk modulus
        7 HARMINDX  I Harmonic index
        """
        s = Struct(self._endian + b'4i2fi')
        nelements = (len(data) - n) // 28
        for unused_i in range(nelements):
            edata = data[n:n + 28]  # 7*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CFLUID3=%s\n' % str(out))
            eid, idf1, idf2, idf3, rho, b, harmonic = out
            self.add_cfluid3(eid, [idf1, idf2, idf3], rho, b, harmonic)
            n += 28
        self.card_count['CFLUID3'] = nelements
        return n

    def _read_cfluid4(self, data, n):
        """
        CFLUID4(8715,87,211) - the marker for Record 37

        1 EID       I Element identification number
        2 IDF1      I RINGFL point 1 identification number
        3 IDF2      I RINGFL point 2 identification number
        4 IDF3      I RINGFL point 3 identification number
        5 IDF4      I RINGFL point 4 identification number
        6 RHO      RS Mass density
        7 B        RS Bulk modulus
        8 HARMINDX  I Harmonic index
        """
        s = Struct(self._endian + b'5i2fi')
        nelements = (len(data) - n) // 32
        for unused_i in range(nelements):
            edata = data[n:n + 32]  # 8*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CFLUID4=%s\n' % str(out))
            eid, idf1, idf2, idf3, idf4, rho, bi, harmonic = out
            self.add_cfluid4(eid, [idf1, idf2, idf3, idf4], rho, bi, harmonic)
            n += 32
        self.card_count['CFLUID4'] = nelements
        return n

# CINT

    def _read_cgap(self, data, n):
        """
        CGAP(1908,19,104) - the marker for Record 39
        """
        s1 = Struct(self._endian + b'4i3fii')
        nelements = (len(data) - n) // 36
        for unused_i in range(nelements):
            edata = data[n:n + 36]  # 9*4
            out = s1.unpack(edata)
            (eid, pid, ga, gb, x1, x2, x3, f, cid) = out  # f=0,1
            g0 = None
            f2, = self.struct_i.unpack(edata[28:32])
            assert f == f2, 'f=%s f2=%s' % (f, f2)
            if f == 2:
                g0, = self.struct_i.unpack(edata[16:20])
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
        s = Struct(self._endian + b'7i')
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n+28]
            out = s.unpack(edata)
            (eid, eid2, side, iviewf, iviewb, radmidf, radmidb) = out
            if self.is_debug_file:
                self.binary_debug.write('  CHBDYE=%s\n' % str(out))
            #self.log.debug('  CHBDYE=%s' % str(out))
            data_in = [eid, eid2, side, iviewf, iviewb, radmidf, radmidb]
            elem = CHBDYE.add_op2_data(data_in)
            self._add_thermal_element_object(elem)
            n += ntotal
        self.card_count['CHBDYE'] = nelements
        return n

    def _read_chbdyg(self, data, n):
        """
        CHBDYG(10808,108,406) - the marker for Record 43
        """
        ntotal = 64  # 16*4
        s = Struct(self._endian + b'16i')
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n+64]
            out = s.unpack(edata)
            (eid, unused_blank, Type, iviewf, iviewb, radmidf, radmidb, unused_blank2,
             g1, g2, g3, g4, g5, g6, g7, g8) = out
            if self.is_debug_file:
                self.binary_debug.write('  CHBDYG=%s\n' % str(out))
            #self.log.debug('  CHBDYG=%s' % str(out))
            data_in = [eid, Type, iviewf, iviewb, radmidf, radmidb,
                       g1, g2, g3, g4, g5, g6, g7, g8]
            elem = CHBDYG.add_op2_data(data_in)
            self._add_thermal_element_object(elem)
            n += ntotal
        self.card_count['CHBDYG'] = nelements
        return n

    def _read_chbdyp(self, data, n):
        """
        CHBDYP(10908,109,407)
        """
        ntotal = 60  # 16*4
        s = Struct(self._endian + b'12i 3f')
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
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
            self._add_thermal_element_object_safe(elem)
            n += ntotal
        self.card_count['CHBDYP'] = nelements
        return n

    def _add_thermal_element_object_safe(self, obj):
        if obj.eid in self.elements:
            self.reject_lines.append(obj.write_card(size=16))
        else:
            self._add_element_object(obj)
        #raise RuntimeError('this should be overwritten by the BDF class')

    def _read_chexa(self, data, n):
        """
        CHEXA(7308,73,253) - the marker for Record 45
        """
        s = Struct(self._endian + b'22i')
        ntotal = 88  # 22*4
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
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
    def _read_chexp(self, data, n):
        """
        CHEXP(12001,120,9011) - the marker for Record 50
        """
        self.log.info('skipping CHEXP in GEOM2')
        if self.is_debug_file:
            self.binary_debug.write('skipping CHEXP in GEOM2\n')
        return len(data)

    def _read_cmass1(self, data, n):
        """
        CMASS1(1001,10,65) - the marker for Record 51
        """
        struct_6i = Struct(self._endian + b'6i')
        nelements = (len(data) - n) // 24
        for unused_i in range(nelements):
            edata = data[n:n + 24]  # 6*4
            out = struct_6i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CMASS1=%s\n' % str(out))
            #(eid, pid, g1, g2, c1, c2) = out
            elem = CMASS1.add_op2_data(out)
            self._add_mass_object(elem)
            n += 24
        self.card_count['CMASS1'] = nelements
        return n

    def _read_cmass2(self, data, n):
        """
        CMASS2(1101,11,66) - the marker for Record 52
        """
        s = Struct(self._endian + b'if4i')
        nelements = (len(data) - n) // 24
        for unused_i in range(nelements):
            edata = data[n:n + 24]  # 6*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CMASS2=%s\n' % str(out))
            #(eid, m, g1, g2, c1, c2) = out
            elem = CMASS2.add_op2_data(out)
            self._add_mass_object(elem)
            n += 24
        self.card_count['CMASS2'] = nelements
        return n

    def _read_cmass3(self, data, n):
        """
        CMASS3(1201,12,67) - the marker for Record 53
        """
        struct_4i = Struct(self._endian + b'4i')
        nelements = (len(data) - n) // 16
        for unused_i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = struct_4i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CMASS3=%s\n' % str(out))
            #(eid, pid, s1, s2) = out
            elem = CMASS3.add_op2_data(out)
            self._add_mass_object(elem)
            n += 16
        self.card_count['CMASS3'] = nelements
        return n

    def _read_cmass4(self, data, n):
        """
        CMASS4(1301,13,68) - the marker for Record 54
        """
        nelements = (len(data) - n) // 16
        struct_if2i = Struct(self._endian + b'ifii')
        for unused_i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = struct_if2i.unpack(edata)
            #(eid, m,s 1, s2) = out
            elem = CMASS4.add_op2_data(out)
            self._add_mass_object(elem)
            n += 16
        self.card_count['CMASS4'] = nelements
        return n

    def _read_cmfree(self, data, n):
        """
        CMFREE(2508,25,0) - the marker for Record 55

        1 EID  I Element identification number
        2   S  I
        3  S2  I
        4   Y RS
        5   N  I
        """
        assert n == 12, n
        nelements = (len(data) - n) // 20
        assert (len(data) - n) % 20 == 0
        struct_3ifi = Struct(self._endian + b'3ifi')
        for unused_i in range(nelements):
            edata = data[n:n + 20]  # 5*4
            out = struct_3ifi.unpack(edata)
            eid, s, s2, y, ncm = out
            self.add_cmfree(eid, s, s2, y, ncm)
            n += 20
        self.card_count['CMFREE'] = nelements
        return n

    def _read_conm1(self, data, n):
        """
        CONM1(1401,14,63) - the marker for Record 56
        """
        s = Struct(self._endian + b'3i21f')
        nelements = (len(data) - n) // 96
        for unused_i in range(nelements):
            edata = data[n:n + 96]  # 24*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CONM1=%s\n' % str(out))
            #(eid, g, cid, m1, m2a, m2b, m3a, m3b, m3c, m4a, m4b, m4c, m4d,
             #m5a, m5b, m5c, m5d, m5e, m6a, m6b, m6c, m6d, m6e, m6f) = out
            elem = CONM1.add_op2_data(out)
            self._add_mass_object(elem)
            n += 96
        self.card_count['CONM1'] = nelements
        return n

    def _read_conm2(self, data, n):
        """
        CONM2(1501,15,64) - the marker for Record 57
        """
        ntotal = 52  # 13*4
        s = Struct(self._endian + b'3i10f')
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n+52]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CONM2=%s\n' % str(out))
            #(eid, g, cid, m, x1, x2, x3, i1, i2a, i2b, i3a, i3b, i3c) = out
            elem = CONM2.add_op2_data(out)
            self._add_mass_object(elem)
            n += ntotal
        self.card_count['CONM2'] = nelements
        return n

    def _read_conrod(self, data, n):
        """
        CONROD(1601,16,47) - the marker for Record 58
        """
        ntotal = 32  # 8*4
        s = Struct(self._endian + b'4i4f')
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n+32]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CONROD=%s\n' % str(out))
            #(eid, n1, n2, mid, a, j, c, nsm) = out
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
            except (AssertionError, MixedVersionCard):
                n, elements = self._read_conv_msc(data, n0)
        else:
            try:
                n, elements = self._read_conv_msc(data, n)
            except (AssertionError, MixedVersionCard):
                n, elements = self._read_conv_nx(data, n0)

        nelements = len(elements)
        for elem in elements:
            self._add_thermal_bc_object(elem, elem.eid)
        self.card_count['CONV'] = nelements
        return n

    def _read_split_card(self, data, n, read1, read2, card_name, add_method):
        """
        generalization of multi read methods for different
        versions of MSC Nastran
        """
        n0 = n
        try:
            n, elements = read1(data, n)
        except AssertionError:
            self.log.info('AssertionError...try again reading %r' % card_name)
            n, elements = read2(data, n0)

        nelements = len(elements)
        for elem in elements:
            add_method(elem)
        self.card_count[card_name] = nelements
        return n

    def _read_dual_card(self, data, n, nx_read, msc_read, card_name, add_method):
        """
        generalization of multi read methods (MSC, NX)
        """
        n0 = n
        if self.is_nx:
            try:
                n, elements = nx_read(data, n)
            except (AssertionError, MixedVersionCard):
                #raise
                n, elements = msc_read(data, n0)
        else:
            try:
                n, elements = msc_read(data, n)
            except (AssertionError, MixedVersionCard):
                #raise
                n, elements = nx_read(data, n0)

        nelements = len(elements)
        assert n is not None
        for elem in elements:
            add_method(elem)
        self.card_count[card_name] = nelements
        return n

    def _read_conv_nx(self, data, n):
        """
        CONV(12701,127,408) - the marker for Record 59
        """
        ntotal = 48  # 12*4
        s = Struct(self._endian + b'4i 8i')
        nelements = (len(data) - n) // ntotal
        elements = []
        for unused_i in range(nelements):
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
        s = Struct(self._endian + b'12i 8f')
        nelements = (len(data) - n) // ntotal
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n+80]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CONV=%s; len=%s\n' % (str(out), len(out)))
            (eid, pcon_id, flmnd, cntrlnd,
             # TODO: why is ta4 and wt4 unused?
             ta1, ta2, ta3, unused_ta4, ta5, ta6, ta7, ta8,
             wt1, wt2, wt3, unused_wt4, wt5, wt6, wt7, wt8) = out
            assert eid > 0, out
            data_in = [eid, pcon_id, flmnd, cntrlnd,
                       [ta1, ta2, ta3, ta5, ta6, ta7, ta8],
                       [wt1, wt2, wt3, wt5, wt6, wt7, wt8]]
            elem = CONV.add_op2_data(data_in)
            elements.append(elem)
            n += ntotal
        return n, elements

    #def _read_convm(self, data, n):
        #"""CONVM"""
        #n = self._read_dual_card(data, n, self._read_convm_nx, self._read_convm_msc,
                                 #'CONVM', self._add_thermal_bc_object)
        #return n

    def _read_convm(self, data, n):
        """
        CONVM(8908,89,422) - the marker for Record 60

        MSC
        1 EID I Element identification number
        2 PCONID I Convection property identification number
        3 FLMND I Point for film convection fluid property temperature
        4 CNTMDOT I Control point used for controlling mass flow.
        5 TA I Ambient points used for convection
        Word 5 repeats 2 times

        NX
        1 EID I Element identification number
        2 PCONID I Convection property identification number
        3 FLMND I Point for film convection fluid property temperature
        4 CNTMDOT I Control point used for controlling mass flow.
        5 TA I Ambient points used for convection
        Word 5 repeats 2 times

        [110, 200, 0, 50000, 99999, 99999, 1.0,
        111, 200, 0, 50000, 99999, 99999, 1.0,
        112, 200, 0, 50000, 99999, 99999, 1.0,
        113, 200, 0, 50000, 99999, 99999, 1.0,
        114, 200, 0, 50000, 99999, 99999, 1.0,
        115, 200, 0, 50000, 99999, 99999, 1.0,
        116, 200, 0, 50000, 99999, 99999, 1.0,
        117, 200, 0, 50000, 99999, 99999, 1.0,
        118, 200, 0, 50000, 99999, 99999, 1.0,
        119, 200, 0, 50000, 99999, 99999, 1.0,
        130, 200, 0, 50000, 99999, 99999, 1.0,
        131, 200, 0, 50000, 99999, 99999, 1.0,
        132, 200, 0, 50000, 99999, 99999, 1.0,
        133, 200, 0, 50000, 99999, 99999, 1.0,
        134, 200, 0, 50000, 99999, 99999, 1.0,
        135, 200, 0, 50000, 99999, 99999, 1.0,
        136, 200, 0, 50000, 99999, 99999, 1.0,
        137, 200, 0, 50000, 99999, 99999, 1.0,
        138, 200, 0, 50000, 99999, 99999, 1.0,
        139, 200, 0, 50000, 99999, 99999, 1.0,
        150, 200, 0, 50000, 99999, 99999, 1.0,
        151, 200, 0, 50000, 99999, 99999, 1.0,
        152, 200, 0, 50000, 99999, 99999, 1.0,
        153, 200, 0, 50000, 99999, 99999, 1.0,
        154, 200, 0, 50000, 99999, 99999, 1.0,
        155, 200, 0, 50000, 99999, 99999, 1.0,
        156, 200, 0, 50000, 99999, 99999, 1.0,
        157, 200, 0, 50000, 99999, 99999, 1.0,
        158, 200, 0, 50000, 99999, 99999, 1.0,
        159, 200, 0, 50000, 99999, 99999, 1.0,
        170, 200, 0, 50000, 99999, 99999, 1.0,
        171, 200, 0, 50000, 99999, 99999, 1.0,
        172, 200, 0, 50000, 99999, 99999, 1.0,
        173, 200, 0, 50000, 99999, 99999, 1.0,
        174, 200, 0, 50000, 99999, 99999, 1.0,
        175, 200, 0, 50000, 99999, 99999, 1.0,
        176, 200, 0, 50000, 99999, 99999, 1.0,
        177, 200, 0, 50000, 99999, 99999, 1.0,
        178, 200, 0, 50000, 99999, 99999, 1.0,
        179, 200, 0, 50000, 99999, 99999, 1.0,
        190, 200, 0, 50000, 99999, 99999, 1.0,
        191, 200, 0, 50000, 99999, 99999, 1.0,
        192, 200, 0, 50000, 99999, 99999, 1.0,
        193, 200, 0, 50000, 99999, 99999, 1.0,
        194, 200, 0, 50000, 99999, 99999, 1.0,
        195, 200, 0, 50000, 99999, 99999, 1.0,
        196, 200, 0, 50000, 99999, 99999, 1.0,
        197, 200, 0, 50000, 99999, 99999, 1.0,
        198, 200, 0, 50000, 99999, 99999, 1.0,
        199, 200, 0, 50000, 99999, 99999, 1.0]
        """
        #C:\Users\sdoyle\Dropbox\move_tpl\ht15330.op2
        ntotal = 24  # 7*4
        struct_6i = Struct(self._endian + b'6i')
        ndata = len(data)
        nelements = (ndata - n) // ntotal
        if (ndata - n) % ntotal != 0:
            msg = 'CONVM error; ndata-n=%s ntotal=%s ndata-n/ntotal=%s' % (
                ndata-n, ntotal, (ndata-n)/float(ntotal))
            self.log.error(msg)
            return n + ndata

        for unused_i in range(nelements):
            edata = data[n:n+24]
            out = struct_6i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CONVM=%s\n' % str(out))
            (eid, pcon_id, flmnd, cntrlnd, ta1, ta2) = out
            if eid <= 0:
                self.show_data(data, 'if')
                # TODO: I'm not sure that this really has 7 fields...
                raise RuntimeError('eid=%s < 0' % eid)
            mdot = 0.
            data_in = [eid, pcon_id, flmnd, cntrlnd, ta1, ta2, mdot]
            elem = CONVM.add_op2_data(data_in)
            self._add_thermal_bc_object(elem, elem.eid)
            n += ntotal
        self.card_count['CONVM'] = nelements
        return n

    def _read_cpyram(self, data, n):
        """
        CPYRAM(17200,172,1000) - the marker for Record ???

        Specific to NX Nastran
        """
        struct_16i = Struct(self._endian + b'16i')
        nelements = (len(data) - n) // 64
        for unused_i in range(nelements):
            edata = data[n:n + 64]  # 15*4
            out = struct_16i.unpack(edata)
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
        s = Struct(self._endian + b'17i')
        nelements = (len(data) - n) // 68
        for unused_i in range(nelements):
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
        s = Struct(self._endian + b'11i')
        nelements = (len(data) - n) // 44  # 11*4
        if self.is_debug_file:
            self.binary_debug.write('ndata=%s\n' % (nelements * 44))
        for unused_i in range(nelements):
            edata = data[n:n + 44]
            out = s.unpack(edata)
            (eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, n9) = out
            if self.is_debug_file:
                self.binary_debug.write('  %s=%s\n' % (element.type, str(out)))
            #print('CQUAD eid=%s pid=%s n1=%s n2=%s n3=%s n4=%s n5=%s n6=%s n7=%s n8=%s' % (
                #eid, pid, n1, n2, n3, n4, n5, n6, n7, n8))
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
        s = Struct(self._endian + b'6iffii4f')
        if self.is_debug_file:
            self.binary_debug.write('ndata=%s\n' % (nelements * 44))

        for unused_i in range(nelements):
            edata = data[n:n + 56]  # 14*4
            out = s.unpack(edata)
            (eid, pid, n1, n2, n3, n4, theta, zoffs, unused_blank, tflag,
             t1, t2, t3, t4) = out
            if self.is_debug_file:
                self.binary_debug.write('  %s=%s\n' % (element.type, str(out)))

            theta_mcid = convert_theta_to_mcid(theta)
            #print('eid=%s pid=%s n1=%s n2=%s n3=%s n4=%s theta=%s zoffs=%s '
                  #'blank=%s tflag=%s t1=%s t2=%s t3=%s t4=%s' % (
                      #eid, pid, n1, n2, n3, n4, theta, zoffs,
                      #blank, tflag, t1, t2, t3, t4))

            data_init = [
                eid, pid, n1, n2, n3, n4, theta_mcid, zoffs,
                tflag, t1, t2, t3, t4]
            elem = element.add_op2_data(data_init)
            self.add_op2_element(elem)
            n += 56
        #if stop:
            #raise RuntimeError('theta is too large...make the quad wrong')
        self.card_count[element.type] = nelements
        return n

# CQUAD4FD

    def _read_cquad8(self, data, n):
        """common method for reading CQUAD8s"""
        n = self._read_split_card(data, n,
                                  self._read_cquad8_current, self._read_cquad8_v2001,
                                  'CQUAD8', self.add_op2_element)
        return n

    def _read_cquad8_current(self, data, n):
        """
        CQUAD8(4701,47,326)  - the marker for Record 72
        .. warning:: inconsistent with dmap manual

        1 EID     I Element identification number
        2 PID     I Property identification number
        3 G(8)    I Grid point identification numbers of
                    connection points
        11 T(4)  RS Membrane thickness of element at grid
                    points
        15 THETA RS Material property orientation angle or
                    coordinate system identification number
        16 ZOFFS RS Offset from the surface of grid points
                    reference plane
        17 TFLAG  I Relative thickness flag
        """
        nelements = (len(data) - n) // 68  # 17*4
        s = Struct(self._endian + b'10i 6f i')
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n + 68]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CQUAD8=%s\n' % str(out))
            #(eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, t1, t2,
             #t3, t4, theta, zoffs, tflag) = out
            tflag = out[-1]
            #self.log.info('cquad8 tflag = %s' % tflag)
            assert isinstance(tflag, int), tflag
            assert tflag in [-1, 0, 1], tflag
            #print('eid=%s pid=%s n1=%s n2=%s n3=%s n4=%s theta=%s zoffs=%s '
                  #'tflag=%s t1=%s t2=%s t3=%s t4=%s' % (
                      #eid, pid, n1, n2, n3, n4, theta, zoffs,
                      #tflag, t1, t2, t3, t4))
            #data_init = [eid,pid,n1,n2,n3,n4,theta,zoffs,tflag,t1,t2,t3,t4]
            elem = CQUAD8.add_op2_data(out)
            elements.append(elem)
            n += 68
        return n, elements

    def _read_cquad8_v2001(self, data, n):
        """
        CQUAD8(4701,47,326)  - the marker for Record 72

        1 EID     I Element identification number
        2 PID     I Property identification number
        3 G(8)    I Grid point identification numbers of
                    connection points
        11 T(4)  RS Membrane thickness of element at grid
                    points
        15 THETA RS Material property orientation angle or
                    coordinate system identification number
        16 ZOFFS RS Offset from the surface of grid points
                    reference plane
        """
        #self.show_data(data, types='if')
        nelements = (len(data) - n) // 64  # 16*4
        s = Struct(self._endian + b'10i 6f')
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n + 64]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CQUAD8=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, t1, t2,
             t3, t4, theta, zoffs) = out
            tflag = None
            out = (eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, t1, t2,
                   t3, t4, theta, zoffs, tflag)
            #print('eid=%s pid=%s n1=%s n2=%s n3=%s n4=%s theta=%s zoffs=%s '
                  #'tflag=%s t1=%s t2=%s t3=%s t4=%s' % (
                      #eid, pid, n1, n2, n3, n4, theta, zoffs, tflag, t1, t2, t3, t4))
            #data_init = [eid,pid,n1,n2,n3,n4,theta,zoffs,tflag,t1,t2,t3,t4]
            elem = CQUAD8.add_op2_data(out)
            elements.append(elem)
            self.add_op2_element(elem)
            n += 64
        return n, elements

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
        struct_4i = Struct(self._endian + b'4i')
        nelements = (len(data) - n) // 16  # 4*4
        #is_long_ids = False
        for unused_i in range(nelements):
            edata = data[n:n + 16]
            out = struct_4i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CROD=%s\n' % str(out))
            #(eid, pid, n1, n2) = out
            #if n1 > 100000000 or n2 > 100000000:
                #is_long_ids = True
            elem = CROD.add_op2_data(out)
            self.add_op2_element(elem)
            n += 16
        #self._is_long_ids = is_long_ids
        self.card_count['CROD'] = nelements
        return n

# CRROD
# CSEAM

    def _read_cshear(self, data, n):
        """
        CSHEAR(3101,31,61)    - the marker for Record 84
        """
        struct_6i = Struct(self._endian + b'6i')
        nelements = (len(data) - n) // 24  # 6*4
        for unused_i in range(nelements):
            edata = data[n:n + 24]
            out = struct_6i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CSHEAR=%s\n' % str(out))
            #(eid, pid, n1, n2, n3, n4) = out
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
        struct_27i = Struct(self._endian + b'27i')
        for unused_i in range(nelements):
            edata = data[n:n+108]
            out = struct_27i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CTETP=%s\n' % str(out))

            eid, pid, n1, n2, n3, n4 = out[:6]
            #(eid, pid, n1, n2, n3, n4,
             #e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12,
             #f1, f2, f3, f4, b1, ee1, ee2, ee3, ee4) = out
            #print("out = ",out)
            #e = [e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12]
            #f = [f1, f2, f3, f4]
            #ee = [ee1, ee2, ee3, ee4]

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
        s = Struct(self._endian + b'12i')
        nelements = (len(data) - n) // 48  # 12*4
        for unused_i in range(nelements):
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
        s = Struct(self._endian + b'5iff3i3f')
        nelements = (len(data) - n)// 52  # 13*4
        for unused_i in range(nelements):
            edata = data[n:n+52]
            out = s.unpack(edata)
            #print('eid=%s pid=%s n1=%s n2=%s n3=%s theta=%s zoffs=%s '
                  #'blank1=%s blank2=%s tflag=%s t1=%s t2=%s t3=%s' % (
                      #eid, pid, n1, n2, n3, theta, zoffs,
                      #blank1, blank2, tflag, t1, t2, t3))
            (eid, pid, n1, n2, n3, theta, zoffs, unused_blank1,
             unused_blank2, tflag, t1, t2, t3) = out
            if self.is_debug_file:
                self.binary_debug.write('  CTRIA3=%s\n' % str(out))

            theta_mcid = convert_theta_to_mcid(theta)
            data_in = [eid, pid, n1, n2, n3, theta_mcid, zoffs, tflag, t1, t2, t3]
            elem = CTRIA3.add_op2_data(data_in)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CTRIA3'] = nelements
        return n


# CTRIAFD - 95

    def _read_ctria6(self, data, n):
        """
        common method for reading CTRIA6

        CTRIA6(4801,48,327) # MSC 2005 - GEOM201
        Word Name Type Description
        1  EID    I Element identification number
        2  PID    I Property identification number
        3  G(6)   I Grid point identification numbers of connection points
        9 THETA  RS Material property orientation angle or coordinate system identification number
        10 ZOFFS RS Offset from the surface of grid points reference plane
        11 T(3)  RS Membrane thickness of element at grid points

        Record 90 -- CTRIA6(4801,48,327) # MSC 2005 - GEOM2
        CTRIA6(4801,48,327)
        Word Name Type Description
        1 EID     I Element identification number
        2 PID     I Property identification number
        3 G(6)    I Grid point identification numbers of connection points
        9 THETA  RS Material property orientation angle or coordinate system identification number
        10 ZOFFS RS Offset from the surface of grid points reference plane
        11 T(3)  RS Membrane thickness of element at grid points
        14 TFLAG  I Relative thickness flag
        """
        n = self._read_split_card(data, n,
                                  self._read_ctria6_current, self._read_ctria6_v2001,
                                  'CTRIA6', self.add_op2_element)
        return n

    def _read_ctria6_current(self, data, n):
        """
        CTRIA6(4801,48,327) - the marker for Record 96

        Record 90 -- CTRIA6(4801,48,327) # MSC 2005 - GEOM2
        Word Name Type Description
        1 EID     I Element identification number
        2 PID     I Property identification number
        3 G(6)    I Grid point identification numbers of connection points
        9 THETA  RS Material property orientation angle or coordinate system identification number
        10 ZOFFS RS Offset from the surface of grid points reference plane
        11 T(3)  RS Membrane thickness of element at grid points
        14 TFLAG  I Relative thickness flag
        """
        s = Struct(self._endian + b'8i 5f i')
        nelements = (len(data) - n) // 56  # 14*4
        assert (len(data) - n) % 56 == 0
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n + 56]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CTRIA6=%s\n' % str(out))
            #print('eid=%s pid=%s n1=%s n2=%s n3=%s theta=%s zoffs=%s '
                  #'blank1=%s blank2=%s tflag=%s t1=%s t2=%s t3=%s' % (
                      #eid, pid, n1, n2, n3, theta, zoffs,
                      #blank1, blank2, tflag, t1, t2, t3))
            #(eid, pid, n1, n2, n3, n4, n5, n6, theta, zoffs, t1, t2, t3, tflag) = out
            tflag = out[-1]
            #self.log.info('ctria6 tflag = %s' % tflag)
            elem = CTRIA6.add_op2_data(out)
            self.add_op2_element(elem)
            assert tflag in [-1, 0, 1], tflag
            elements.append(elem)
            n += 56
        return n, elements

    def _read_ctria6_v2001(self, data, n):
        """
        CTRIA6(4801,48,327) - the marker for Record 96
        """
        s = Struct(self._endian + b'8i 5f')
        nelements = (len(data) - n) // 52  # 13*4
        assert (len(data) - n) % 52 == 0
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n + 52]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CTRIA6=%s\n' % str(out))
            #print('eid=%s pid=%s n1=%s n2=%s n3=%s theta=%s zoffs=%s '
                  #'blank1=%s blank2=%s tflag=%s t1=%s t2=%s t3=%s' % (
                      #eid, pid, n1, n2, n3, theta, zoffs,
                      #blank1, blank2, tflag, t1, t2, t3))
            (eid, pid, n1, n2, n3, n4, n5, n6, theta, zoffs, t1, t2, t3) = out
            out = (eid, pid, n1, n2, n3, n4, n5, n6, theta, zoffs, t1, t2, t3, 0)
            elem = CTRIA6.add_op2_data(out)
            elements.append(elem)
            n += 52
        return n, elements

# CTRIA6FD
# CTRIAP

    def _read_ctriar(self, data, n):
        """
        CTRIAR(9200,92,385)    - the marker for Record 99
        """
        ntotal = 52  # 13*4
        s = Struct(self._endian + b'5iff3i3f')
        nelements = (len(data) - n)// 52  # 13*4
        for unused_i in range(nelements):
            edata = data[n:n+52]
            out = s.unpack(edata)
            #print('eid=%s pid=%s n1=%s n2=%s n3=%s theta=%s zoffs=%s '
                  #'blank1=%s blank2=%s tflag=%s t1=%s t2=%s t3=%s' % (
                      #eid, pid, n1, n2, n3, theta, zoffs,
                      #blank1, blank2, tflag, t1, t2, t3))
            (eid, pid, n1, n2, n3, theta, zoffs, unused_blank1,
             unused_blank2, tflag, t1, t2, t3) = out
            if self.is_debug_file:
                self.binary_debug.write('  CTRIAR=%s\n' % str(out))
            data_in = [eid, pid, n1, n2, n3, theta, zoffs, tflag, t1, t2, t3]
            elem = CTRIAR.add_op2_data(data_in)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CTRIAR'] = nelements
        return n

    def _read_ctriax(self, data, n): # 100
        self.log.info('skipping CTRIAX in GEOM2')
        return len(data)

    def _read_ctriax6(self, data, n):  # 101
        """(6108, 61, 107)"""
        ntotal = 44  # 11*4
        nentries = (len(data) - n) // ntotal
        struc = Struct(self._endian + b'8i f ii')
        for unused_i in range(nentries):
            edata = data[n:n + 44]
            out = struc.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CTRIAX6=%s\n' % str(out))
            elem = CTRIAX6.add_op2_data(out)
            self.add_op2_element(elem)
            n += 44
        self.card_count['CTRIAX6'] = nentries
        return n

# CTRIX3FD - 102
# CTRIX6FD - 103

    def _read_ctube(self, data, n):
        """
        CTUBE(3701,37,49) - the marker for Record 104
        """
        struct_4i = Struct(self._endian + b'4i')
        nelements = (len(data) - n) // 16
        for unused_i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = struct_4i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CTUBE=%s\n' % str(out))
            #(eid, pid, n1, n2) = out
            elem = CTUBE.add_op2_data(out)
            self.add_op2_element(elem)
            n += 16
        self.card_count['CTUBE'] = nelements
        return n

    def _read_cvisc(self, data, n):
        """CVISC(3901,39,50) - the marker for Record 105"""
        struct_4i = Struct(self._endian + b'4i')
        nelements = (len(data) - n) // 16
        for unused_i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = struct_4i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CVISC=%s\n' % str(out))
            #(eid, pid, n1, n2) = out
            element = CVISC.add_op2_data(out)
            self.add_op2_element(element)
            n += 16
        self.card_count['CVISC'] = nelements
        return n

    def _read_cweld(self, data, n):
        """
        CWELD(11701,117,559) - Record 106
        same as CFAST
        """
        self.log.info('skipping CWELD in GEOM2')
        if self.is_debug_file:
            self.binary_debug.write('skipping CWELD in GEOM2\n')
        return len(data)

    def _read_cweldc(self, data, n):  # 107
        self.log.info('skipping CWELDC in GEOM2')
        if self.is_debug_file:
            self.binary_debug.write('skipping CWELDC in GEOM2\n')
        return len(data)

    def _read_cweldg(self, data, n):  # 108
        self.log.info('skipping CWELDG in GEOM2')
        if self.is_debug_file:
            self.binary_debug.write('skipping CWELDG in GEOM2\n')
        return len(data)

# TDOO: above are checked by DMAP...
#-------------------------------
# CWSEAM
    def _read_genel(self, data, n):
        self.log.info('skipping GENEL in GEOM2')
        return len(data)
# GMDNDC
# GMBNDS
# GMINTC
# GMINTS

    def _read_plotel(self, data, n):  # 114
        """(5201, 52, 11)"""
        struct_3i = Struct(self._endian + b'3i')
        ntotal = 12
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]  # 4*4
            out = struct_3i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  PLOTEL=%s\n' % str(out))
            #(eid,n1,n2) = out
            elem = PLOTEL.add_op2_data(out)
            self._add_plotel_object(elem)
            n += ntotal
        self.card_count['PLOTEL'] = nelements
        return n

    def _read_radbc(self, data, n):
        self.log.info('skipping RADBC in GEOM2')
        return len(data)

# RADINT
# SINT

    def _read_spoint(self, data, n):
        """
        (5551,49,105)    - the marker for Record 118
        """
        npoints = (len(data) - n) // 4
        nids = np.frombuffer(data[n:], self.idtype).tolist()
        if self.is_debug_file:
            self.binary_debug.write('SPOINT=%s\n' % nids)
        spoint = SPOINTs.add_op2_data(nids)
        self._add_spoint_object(spoint)
        self.card_count['SPOINT'] = npoints
        return n

    def _read_vubeam(self, data, n):  # 119
        self.log.info('skipping VUBEAM in GEOM2')
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

    def _read_trax3(self, data, n):
        self.log.info('skipping CTRAX3 in GEOM2')
        return len(data)

    def _read_cquadx4(self, data, n):
        self.log.info('skipping CQUADX4 in GEOM2')
        return len(data)

    def _read_ctrax6(self, data, n):
        self.log.info('skipping CTRAX6 in GEOM2')
        return len(data)

    def _read_cquadx8(self, data, n):
        self.log.info('skipping CQUADX8 in GEOM2')
        return len(data)

def convert_theta_to_mcid(theta):
    """odd function..."""
    # sort of guessed at this number...it seems reasonable-ish
    if theta > 511.:
        # per DMAP...you couldn't make a new record number?
        # theta = 512. * (cid + 1)
        # theta/512 = cid + 1
        # cid = theta/512. - 1
        #
        cid_float = theta / 512. - 1
        cid = int(cid_float)
        assert np.allclose(cid, cid_float), 'theta=%s cid=%s cid_float=%s' % (theta, cid, cid_float)
        theta = cid
    return theta
