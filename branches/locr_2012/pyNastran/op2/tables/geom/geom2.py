# pylint: disable=E1101,C0103,R0902,R0904,R0914
import sys
from struct import unpack

from pyNastran.bdf.cards.elements.elements import CGAP
from pyNastran.bdf.cards.elements.damper import (CDAMP1, CDAMP2, CDAMP3,
                                                 CDAMP4, CDAMP5, CVISC)
from pyNastran.bdf.cards.elements.springs import CELAS1, CELAS2, CELAS3, CELAS4
from pyNastran.bdf.cards.elements.shell import (CTRIA3, CQUAD4, CTRIA6,
                                                CQUADR, CQUAD8, CQUAD, CQUADX,
                                                CSHEAR)
from pyNastran.bdf.cards.elements.bars import (CROD, CBAR, CTUBE, CONROD,
                                               CBEAM)
from pyNastran.bdf.cards.elements.mass import (CONM1, CONM2, CMASS1, CMASS2,
                                               CMASS3, CMASS4)
from pyNastran.bdf.cards.elements.solid import (CTETRA4, CTETRA10, CPENTA6,
                                                CPENTA15, CHEXA8, CHEXA20)
from pyNastran.bdf.cards.thermal.thermal import CHBDYG, CONV  # , CONVM, CHBDYP
from pyNastran.bdf.cards.nodes import SPOINTs


class Geometry2(object):
    def readTable_Geom2(self):
        self.iTableMap = {
            (2408, 24, 180): self.readCBAR,    # record 8
            (4001, 40, 275): self.readCBARAO,  # record 9  - not done
            (5408, 54, 261): self.readCBEAM,   # record 10
            (11401, 114, 9016): self.readCBEAMP,  # record 11 - not done
            (4601, 46, 298): self.readCBEND,   # record 12 - not done
            (5608, 56, 218): self.readCBUSH1D,  # record 14 - not done
            (2315, 23, 146): self.readCCONE,   # record 15 - not done
            (201, 2, 69): self.readCDAMP1,  # record 16
            (301, 3, 70): self.readCDAMP2,  # record 17
            (401, 4, 71): self.readCDAMP3,  # record 18
            (501, 5, 72): self.readCDAMP4,  # record 19
            (10608, 106, 404): self.readCDAMP5,  # record 20
            (601, 6, 73): self.readCELAS1,  # record 29
            (701, 7, 74): self.readCELAS2,  # record 30
            (801, 8, 75): self.readCELAS3,  # record 31
            (901, 9, 76): self.readCELAS4,  # record 32
            #(8515,85,209):    self.readCFLUID2, # record 35 - not done
            #(8615,86,210):    self.readCFLUID3, # record 36 - not done
            #(8715,87,211):    self.readCFLUID4, # record 37 - not done
            (1908, 19, 104): self.readCGAP,    # record 39 - buggy

            (10808, 108, 406): self.readCHBDYG,   # record 43
            (10908, 109, 407): self.readCHBDYP,   # record 44 - not done
            (7308, 73, 253): self.readCHEXA,    # record 45
            (1001, 10, 65): self.readCMASS1,   # record 51
            (1101, 11, 66): self.readCMASS2,   # record 52
            (1201, 12, 67): self.readCMASS3,   # record 53
            (1301, 13, 68): self.readCMASS4,   # record 54
            (2508, 25, 0): self.readCMFREE,   # record 55 - not done
            (1401, 14, 63): self.readCONM1,    # record 56 - not done
            (1501, 15, 64): self.readCONM2,    # record 57
            (1601, 16, 47): self.readCONROD,   # record 58
            (12701, 127, 408): self.readCONV,     # record 59 - not tested
            (8908, 89, 422): self.readCONVM,    # record 60 - not tested
            (4108, 41, 280): self.readCPENTA,   # record 62

            (9108, 91, 507): self.readCQUAD,    # record 68 - not tested
            (2958, 51, 177): self.readCQUAD4,   # record 69 - maybe buggy on theta/Mcsid field
            (13900, 139, 9989): self.readCQUAD4,   # record 70 - maybe buggy on theta/Mcsid field
            (4701, 47, 326): self.readCQUAD8,   # record 71 - maybe buggy on theta/Mcsid field
            (8009, 80, 367): self.readCQUADR,   # record 74 - not tested
            (9008, 90, 508): self.readCQUADX,   # record 75 - not tested

            (3001, 30, 48): self.readCROD,     # record 80
            #(12201,122,9013): self.readCTETP,    # record 86 - not done
            (5508, 55, 217): self.readCTETRA,   # record 87
            (5959, 59, 282): self.readCTRIA3,   # record 93 - maybe buggy on theta/Mcsid field
            (4801, 48, 327): self.readCTRIA6,   # record 95 - buggy
            (9200, 92, 385): self.readCTRIAR,   # record 98  - not done
            (6108, 61, 107): self.readCTRIAX6,  # record 100 - not done
            (3701, 37, 49): self.readCTUBE,    # record 103
            (3901, 39, 50): self.readCVISC,   # record 104 - not done
            #(5201,52,11):      self.readPLOTEL,  # record 114 - not done
            (5551, 49, 105): self.readSPOINT,   # record 118
            #(11601,116,9942):  self.readVUBEAM,  # record 119 - not done
            #(2608, 26, 60)
        }
        self.readRecordTable('GEOM2')

    def readTable_Geom2S(self):
        self.readTable_Geom2()
        #self.iTableMap = {
        #                 }
        #self.readRecordTable('GEOM2S')

    def readFake(self, data):
        pass

    def addOp2Element(self, elem):
        self.addElement(elem, allowOverwrites=True)

# 1-AEROQ4 (???)
# AEROT3   (???)
# 1-BEAMAERO (1701,17,0)
# 2-CAABSF (2708,27,59)
# 3-CAXIF2 (2108,21,224)
# 4-CAXIF3 (2208,22,225)
# 5-CAXIF4 (2308,23,226)

    def readCBAR(self, data):
        """
        CBAR(2408,24,180) - the marker for Record 8
        """
        #print "reading CBAR"
        n = 0
        nEntries = len(data) // 64
        for i in xrange(nEntries):
            eData = data[n:n + 64]  # 16*4
            f, = unpack(b'i', eData[28:32])
            #print "len(eData) = %s" %(len(eData))
            if   f == 0:
                out = unpack(b'iiiifffiiiffffff', eData)
                (eid, pid, ga, gb, x1, x2, x3, f, pa, pb,
                    w1a, w2a, w3a, w1b, w2b, w3b) = out
                dataIn = [[eid, pid, ga, gb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                          [f, x1, x2, x3]]
            elif f == 1:
                out = unpack(b'iiiifffiiiffffff', eData)
                (eid, pid, ga, gb, x1, x2, x3, f, pa, pb,
                    w1a, w2a, w3a, w1b, w2b, w3b) = out
                dataIn = [[eid, pid, ga, gb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                          [f, x1, x2, x3]]
            elif f == 2:
                out = unpack(b'iiiiiiifiiffffff', eData)
                (eid, pid, ga, gb, g0, junk, junk, f, pa,
                    pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                dataIn = [[eid, pid, ga, gb, pa, pb, w1a,
                           w2a, w3a, w1b, w2b, w3b], [f, g0]]
            else:
                raise RuntimeError('invalid f value...f=%s' % (f))
            elem = CBAR(None, dataIn)
            self.addOp2Element(elem)
            n += 64
        data = data[n:]

    def readCBARAO(self, data):
        """
        CBARAO(4001,40,275) - the marker for Record 9
        """
        self.skippedCardsFile.write('skipping CBARAO in GEOM2\n')

    def readCBEAM(self, data):
        """
        CBEAM(5408,54,261) - the marker for Record 10
        """
        #print "reading CBEAM"
        n = 0
        nEntries = len(data) // 72
        for i in xrange(nEntries):
            eData = data[n:n + 72]  # 18*4
            f, = unpack(b'i', eData[40:44])
            #print "f = ",f
            #print "len(eData) = %s" %(len(eData))
            if   f == 0:  # basic cid
                out = unpack(b'iiiiiifffiiiffffff', eData)
                (eid, pid, ga, gb, sa, sb, x1, x2, x3, f, pa,
                    pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                dataIn = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                          [f, x1, x2, x3]]
            elif f == 1:  # global cid
                out = unpack(b'iiiiiifffiiiffffff', eData)
                (eid, pid, ga, gb, sa, sb, x1, x2, x3, f, pa,
                    pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                dataIn = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                          [f, x1, x2, x3]]
            elif f == 2:  # grid option
                out = unpack(b'iiiiiiiiiiiiffffff', eData)
                (eid, pid, ga, gb, sa, sb, g0, xx, xx, f, pa,
                    pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                dataIn = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                          [f, g0]]
            else:
                raise RuntimeError('invalid f value...f=%s' % (f))
            elem = CBEAM(None, dataIn)
            self.addOp2Element(elem)
            n += 72
        data = data[n:]

    def readCBEAMP(self, data):
        """
        CBEAMP(11401,114,9016) - the marker for Record 11
        """
        self.skippedCardsFile.write('skipping CBEAMP in GEOM2\n')

    def readCBEND(self, data):
        """
        CBEND(4601,46,298) - the marker for Record 12
        """
        self.skippedCardsFile.write('skipping CBEND in GEOM2\n')

    def readCBUSH(self, data):
        """
        CBUSH(2608,26,60) - the marker for Record 13
        """
        self.skippedCardsFile.write('skipping CBUSH in GEOM2\n')

    def readCBUSH1D(self, data):
        """
        CBUSH1D(5608,56,218) - the marker for Record 14
        """
        self.skippedCardsFile.write('skipping CBUSH1D in GEOM2\n')

    def readCCONE(self, data):
        """
        CCONE(2315,23,0) - the marker for Record 15
        """
        self.skippedCardsFile.write('skipping CCONE in GEOM2\n')

    def readCDAMP1(self, data):
        """
        CDAMP1(201,2,69) - the marker for Record 16
        """
        #print "reading CDAMP1"
        n = 0
        nEntries = len(data) // 24
        for i in xrange(nEntries):
            eData = data[n:n + 24]  # 6*4
            out = unpack(b'iiiiii', eData)
            (eid, pid, g1, g2, c1, c2) = out
            elem = CDAMP1(None, out)
            self.addOp2Element(elem)
            n += 24
        data = data[n:]

    def readCDAMP2(self, data):
        """
        CDAMP2(301,3,70) - the marker for Record 17
        """
        #print "reading CDAMP2"
        n = 0
        nEntries = len(data) // 24
        for i in xrange(nEntries):
            eData = data[n:n + 24]  # 6*4
            out = unpack(b'ifiiii', eData)
            (eid, b, g1, g2, c1, c2) = out
            elem = CDAMP2(None, out)
            self.addOp2Element(elem)
            n += 24
        data = data[n:]

    def readCDAMP3(self, data):
        """
        CDAMP3(401,4,71) - the marker for Record 18
        """
        #print "reading CDAMP3"
        n = 0
        nEntries = len(data) // 16
        for i in xrange(nEntries):
            eData = data[n:n + 16]  # 4*4
            out = unpack(b'iiii', eData)
            (eid, pid, s1, s2) = out
            elem = CDAMP3(None, out)
            self.addOp2Element(elem)
            n += 16
        data = data[n:]

    def readCDAMP4(self, data):
        """
        CDAMP4(501,5,72) - the marker for Record 19
        """
        #print "reading CDAMP4"
        n = 0
        nEntries = len(data) // 16
        for i in xrange(nEntries):
            eData = data[n:n + 16]  # 4*4
            out = unpack(b'ifii', eData)
            (eid, b, s1, s2) = out
            elem = CDAMP4(None, out)
            self.addOp2Element(elem)
            n += 16
        data = data[n:]

    def readCDAMP5(self, data):
        """
        CDAMP5(10608,106,404) - the marker for Record 20
        """
        #print "reading CDAMP5"
        n = 0
        nEntries = len(data) // 16
        for i in xrange(nEntries):
            eData = data[n:n + 16]  # 4*4
            out = unpack(b'iiii', eData)
            (eid, pid, s1, s2) = out
            elem = CDAMP5(None, out)
            self.addOp2Element(elem)
            n += 16
        data = data[n:]

# CDUM2
# CDUM3
# CDUM4
# CDUM5
# CDUM6
# CDUM7
# CDUM8
# CDUM9

    def readCELAS1(self, data):
        """
        CELAS1(601,6,73) - the marker for Record 29
        """
        #print "reading CELAS1"
        while len(data) >= 24:  # 6*4
            eData = data[:24]
            data = data[24:]
            out = unpack(b'iiiiii', eData)
            (eid, pid, g1, g2, c1, c2) = out
            elem = CELAS1(None, out)
            self.addOp2Element(elem)

    def readCELAS2(self, data):
        """
        CELAS2(701,7,74) - the marker for Record 30
        """
        #print "reading CELAS2"
        while len(data) >= 32:  # 8*4
            eData = data[:32]
            data = data[32:]
            out = unpack(b'ifiiiiff', eData)
            (eid, k, g1, g2, c1, c2, ge, s) = out
            elem = CELAS2(None, out)
            self.addOp2Element(elem)

    def readCELAS3(self, data):
        """
        CELAS3(801,8,75) - the marker for Record 31
        """
        #print "reading CELAS3"
        while len(data) >= 16:  # 4*4
            eData = data[:16]
            data = data[16:]
            out = unpack(b'iiii', eData)
            (eid, pid, s1, s2) = out
            elem = CELAS3(None, out)
            self.addOp2Element(elem)

    def readCELAS4(self, data):
        """
        CELAS4(901,9,76) - the marker for Record 32
        """
        #print "reading CELAS4"
        n = 0
        nEntries = len(data) // 16
        for i in xrange(nEntries):
            eData = data[n:n + 16]  # 4*4
            out = unpack(b'ifii', eData)
            (eid, k, s1, s2) = out
            elem = CELAS4(None, out)
            self.addOp2Element(elem)
            n += 16
        data = data[n:]

# CFAST
# CFASTP
# CFLUID2
# CFLUID3
# CFLUID4
# CINT

    def readCGAP(self, data):
        """
        CGAP(1908,19,104) - the marker for Record 39
        """
        #print "reading CGAP"
        n = 0
        nEntries = len(data) // 36
        for i in xrange(nEntries):
            eData = data[n:n + 36]  # 9*4
            out = unpack(b'iiiifffii', eData)
            (eid, pid, ga, gb, x1, x2, x3, f, cid) = out  # f=0,1
            g0 = None
            f2, = unpack(b'i', eData[28:32])
            assert f == f2, 'f=%s f2=%s' % (f, f2)
            if f == 2:
                g0 = unpack(b'i', eData[16:20])
                x1 = None
                x2 = None
                x3 = None

            dataIn = [eid, pid, ga, gb, g0, x1, x2, x3, cid]
            elem = CGAP(None, dataIn)
            self.addOp2Element(elem)
            n += 36
        data = data[n:]

# CHACAB
# CHACBR
# CHBDYE
# CHBDYG

    def readCHBDYG(self, data):
        """
        CHBDYG(10808,108,406) - the marker for Record 43
        """
        return
        #print "reading CHBDYG"
        while len(data) >= 64:  # 16*4
            eData = data[:64]
            data = data[64:]
            (eid, blank, Type, iviewf, iviewb, radmidf, radmidb, blank2,
             g1, g2, g3, g4, g5, g6, g7, g8) = unpack(b'iiiiiiiiiiiiiiii', eData)
            dataIn = [eid, Type, iviewf, iviewb, radmidf, radmidb,
                      g1, g2, g3, g4, g5, g6, g7, g8]
            elem = CHBDYG(None, dataIn)
            self.addOp2Element(elem)

    def readCHBDYP(self, data):
        self.skippedCardsFile.write('skipping CHBDYP in GEOM2\n')

    def readCHEXA(self, data):
        """
        CHEXA(7308,73,253) - the marker for Record 45
        """
        #print "reading CHEXA"
        while len(data) >= 88:  # 22*4
            eData = data[:88]
            data = data[88:]
            out = unpack(b'22i', eData)
            (eid, pid, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10,
             g11, g12, g13, g14, g15, g16, g17, g18, g19, g20) = out

            dataIn = [eid, pid, g1, g2, g3, g4, g5, g6, g7, g8, ]
            bigNodes = [g9, g10, g11, g12, g13, g14, g15, g16,
                        g17, g18, g19, g20]
            if sum(bigNodes) > 0:
                elem = CHEXA20(None, dataIn + bigNodes)
            else:
                elem = CHEXA8(None, dataIn)
            self.addOp2Element(elem)

# CHEXA20F
# CHEXAFD
# CHEXAL
# CHEXP
# CHEXPR

    def readCMASS1(self, data):
        """
        CMASS1(1001,10,65) - the marker for Record 51
        """
        #print "reading CMASS1"
        n = 0
        nEntries = len(data) // 24
        for i in xrange(nEntries):
            eData = data[n:n + 24]  # 6*4
            out = unpack(b'iiiiii', eData)
            #(eid,pid,g1,g2,c1,c2) = out
            elem = CMASS1(None, out)
            self.addOp2Element(elem)
            n += 24
        data = data[n:]

    def readCMASS2(self, data):
        """
        CMASS2(1101,11,66) - the marker for Record 52
        """
        #print "reading CMASS2"
        n = 0
        nEntries = len(data) // 24
        for i in xrange(nEntries):
            eData = data[n:n + 24]  # 6*4
            out = unpack(b'ifiiii', eData)
            #(eid,m,g1,g2,c1,c2) = out
            elem = CMASS2(None, out)
            self.addOp2Element(elem)
            n += 24
        data = data[n:]

    def readCMASS3(self, data):
        """
        CMASS3(1201,12,67) - the marker for Record 53
        """
        #print "reading CMASS3"
        n = 0
        nEntries = len(data) // 16
        for i in xrange(nEntries):
            eData = data[n:n + 16]  # 4*4
            out = unpack(b'iiii', eData)
            #(eid,pid,s1,s2) = out
            elem = CMASS3(None, out)
            self.addOp2Element(elem)
            n += 16
        data = data[n:]

    def readCMASS4(self, data):
        """
        CMASS4(1301,13,68) - the marker for Record 54
        """
        #print "reading CMASS4"
        n = 0
        nEntries = len(data) // 16
        for i in xrange(nEntries):
            eData = data[n:n + 16]  # 4*4
            out = unpack(b'ifii', eData)
            #(eid,m,s1,s2) = out
            elem = CMASS4(None, out)
            self.addOp2Element(elem)
            n += 16
        data = data[n:]

    def readCMFREE(self, data):
        """
        CMFREE(2508,25,0) - the marker for Record 55
        """
        self.skippedCardsFile.write('skipping CMFREE in GEOM2\n')

    def readCONM1(self, data):
        """
        CONM1(1401,14,63) - the marker for Record 56
        """
        #print "reading CONM1"
        n = 0
        nEntries = len(data) // 96
        for i in xrange(nEntries):
            eData = data[n:n + 96]  # 24*4
            out = unpack(b'iii21f', eData)
            (eid, g, cid, m1, m2a, m2b, m3a, m3b, m3c, m4a, m4b, m4c, m4d,
             m5a, m5b, m5c, m5d, m5e, m6a, m6b, m6c, m6d, m6e, m6f) = out
            elem = CONM1(None, out)
            self.addOp2Element(elem)
            n += 96
        data = data[n:]

    def readCONM2(self, data):
        """
        CONM2(1501,15,64) - the marker for Record 57
        """
        #print "reading CONM2"
        while len(data) >= 52:  # 13*4
            eData = data[:52]
            data = data[52:]
            out = unpack(b'iiiffffffffff', eData)
            (eid, g, cid, m, x1, x2, x3, i1, i2a, i2b, i3a, i3b, i3c) = out
            elem = CONM2(None, out)
            self.addOp2Element(elem)

    def readCONROD(self, data):
        """
        CONROD(1601,16,47) - the marker for Record 58
        """
        #print "reading CONROD"
        while len(data) >= 32:  # 8*4
            eData = data[:32]
            data = data[32:]
            out = unpack(b'iiiiffff', eData)
            (eid, n1, n2, mid, a, j, c, nsm) = out
            elem = CONROD(None, out)
            self.addOp2Element(elem)

    def readCONV(self, data):
        """
        CONV(12701,127,408) - the marker for Record 59
        """
        #print "reading CONV"
        return
        while len(data) >= 80:  # 20*4
            eData = data[:80]
            data = data[80:]
            out = unpack(b'iiiiiiiiiiiiffffffff', eData)
            (eid, pconID, flmnd, cntrlnd,
             ta1, ta2, ta3, ta5, ta6, ta7, ta8,
             wt1, wt2, wt3, wt5, wt6, wt7, wt8) = out
            dataIn = [eid, pconID, flmnd, cntrlnd,
                      [ta1, ta2, ta3, ta5, ta6, ta7, ta8],
                      [wt1, wt2, wt3, wt5, wt6, wt7, wt8]]
            elem = CONV(None, dataIn)
            self.addOp2Element(elem)

    def readCONVM(self, data):
        """
        CONVM(8908,89,422) - the marker for Record 60
        """
        #print "reading CONVM"
        return
        while len(data) >= 28:  # 7*4
            eData = data[:28]
            data = data[28:]
            out = unpack(b'iiiiiii', eData)
            (eid, pconID, flmnd, cntrlnd,
             [ta1, ta2, ta3]) = out
            dataIn = [eid, pconID, flmnd, cntrlnd,
                      [ta1, ta2, ta3]]
            elem = CONVM(None, dataIn)
            self.addOp2Element(elem)

# CPENP

    def readCPENTA(self, data):
        """
        CPENTA(4108,41,280) - the marker for Record 62
        """
        #print "reading CPENTA"
        n = 0
        nEntries = len(data) // 68
        for i in xrange(nEntries):
            eData = data[n:n + 68]  # 17*4
            out = unpack(b'iiiiiiiiiiiiiiiii', eData)
            (eid, pid, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10,
             g11, g12, g13, g14, g15) = out

            dataIn = [eid, pid, g1, g2, g3, g4, g5, g6]
            bigNodes = [g7, g8, g9, g10, g11, g12, g13, g14, g15]
            if sum(bigNodes) > 0:
                elem = CPENTA15(None, dataIn + bigNodes)
            else:
                elem = CPENTA6(None, dataIn)
            self.addOp2Element(elem)
            n += 68
        data = data[n:]

# CPENPR
# CPENT15F
# CPENT6FD
# CQDX4FD
# CQDX9FD

    def readCQUAD(self, data):
        """
        CQUAD(9108,91,507)  - the marker for Record 68
        """
        #print "reading CQUAD"
        self.runCQUAD(data, CQUAD)

    def runCQUAD(self, data, Element):
        """common method for CQUAD, CQUADX"""
        n = 0
        nEntries = len(data) // 44  # 11*4
        for i in xrange(nEntries):
            eData = data[n:n + 44]
            (eid, pid, n1, n2, n3, n4, n5, n6, n7, n8,
                n9) = unpack(b'iiiiiiiiiii', eData)
            #print "eid=%s pid=%s n1=%s n2=%s n3=%s n4=%s theta=%s zoffs=%s tflag=%s t1=%s t2=%s t3=%s t4=%s" %(eid,pid,n1,n2,n3,n4,theta,zoffs,tflag,t1,t2,t3,t4)
            #dataInit = [eid,pid,n1,n2,n3,n4,theta,zoffs,tflag,t1,t2,t3,t4]
            data = [eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, n9]
            elem = Element(None, data)
            self.addOp2Element(elem)
            n += 44
        data = data[n:]

    def readCQUAD4(self, data):
        """
        CQUAD4(2958,51,177)    - the marker for Record 69
        CQUAD4(13900,139,9989) - the marker for Record 70
        """
        #print "reading CQUAD4"
        self.runCQUAD4(data, CQUAD4)

    def runCQUAD4(self, data, Element):
        """
        common method for CQUAD4, CQUADR
        """
        n = 0
        nEntries = len(data) // 56
        for i in xrange(nEntries):
            eData = data[n:n + 56]  # 14*4
            (eid, pid, n1, n2, n3, n4, theta, zoffs, blank, tflag,
                t1, t2, t3, t4) = unpack(b'iiiiiiffiiffff', eData)
            #print "eid=%s pid=%s n1=%s n2=%s n3=%s n4=%s theta=%s zoffs=%s blank=%s tflag=%s t1=%s t2=%s t3=%s t4=%s" %(eid,pid,n1,n2,n3,n4,theta,zoffs,blank,tflag,t1,t2,t3,t4)
            dataInit = [eid, pid, n1, n2, n3, n4, theta, zoffs,
                        tflag, t1, t2, t3, t4]
            elem = Element(None, dataInit)
            self.addOp2Element(elem)
            n += 56
        data = data[n:]

# CQUAD4FD

    def readCQUAD8(self, data):
        """
        CQUAD8(4701,47,326)  - the marker for Record 71
        @warning inconsistent with dmap manual
        """
        #print "reading CQUAD8"
        return
        n = 0
        nEntries = len(data) // 64  # 17*4
        for i in xrange(nEntries):
            eData = data[n:n + 64]
            out = unpack(b'iiiiiiiiiifffffi', eData)
            (eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, t1, t2,
                t3, t4, theta, tflag) = out
            #print "eid=%s pid=%s n1=%s n2=%s n3=%s n4=%s theta=%s zoffs=%s tflag=%s t1=%s t2=%s t3=%s t4=%s" %(eid,pid,n1,n2,n3,n4,theta,zoffs,tflag,t1,t2,t3,t4)
            #dataInit = [eid,pid,n1,n2,n3,n4,theta,zoffs,tflag,t1,t2,t3,t4]
            elem = CQUAD8(None, out)
            self.addOp2Element(elem)
            n += 64
        data = data[n:]

# CQUAD9FD
# CQUADP
    def readCQUADR(self, data):
        """
        CQUADR(8009,80,367)  - the marker for Record 74
        """
        #print "reading CQUADR"
        self.runCQUAD4(data, CQUADR)

    def readCQUADX(self, data):
        """
        CQUADX(9008,90,508)  - the marker for Record 75
        """
        #print "reading CQUADX"
        self.runCQUAD4(data, CQUADX)

# CRBAR
# CRBE1
# CRBE3
# CRJOINT

    def readCROD(self, data):
        """
        CROD(3001,30,48)    - the marker for Record 80
        """
        #print "reading CROD"
        n = 0
        nEntries = len(data) // 16  # 4*4
        for i in xrange(nEntries):
            eData = data[n:n + 16]
            out = unpack(b'iiii', eData)
            (eid, pid, n1, n2) = out
            elem = CROD(None, out)
            self.addOp2Element(elem)
            n += 16
        data = data[n:]

# CRROD
# CSEAM

    def readCSHEAR(self, data):
        """
        CSHEAR(3101,31,61)    - the marker for Record 83
        """
        #print "reading CSHEAR"
        n = 0
        nEntries = len(data) // 24  # 6*4
        for i in xrange(nEntries):
            eData = data[n:n + 24]
            out = unpack(b'6i', eData)
            (eid, pid, n1, n2, n3, n4) = out
            elem = CSHEAR(None, out)
            self.addOp2Element(elem)
            n += 24
        data = data[n:]

# CSLOT3
# CSLOT4

    def readCTETP(self, data):
        """
        CTETP(12201,122,9013)    - the marker for Record 86
        @todo create object
        """
        #print "reading CTETP"
        raise NotImplementedError('needs work...')
        while len(data) >= 108:  # 27*4
            eData = data[:108]
            data = data[108:]
            out = unpack(b'27i', eData)
            (
                eid, pid, n1, n2, n3, n4, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12,
                f1, f2, f3, f4, b1, ee1, ee2, ee3, ee4) = out
            #print "out = ",out
            e = [e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12]
            f = [f1, f2, f3, f4]
            ee = [ee1, ee2, ee3, ee4]

            #print "e  = ",e
            #print "f  = ",f
            #print "b1  = ",b1
            #print "ee = ",ee
            dataIn = [eid, pid, n1, n2, n2, n3, n4]
            elem = CTETRA4(None, dataIn)
            self.addOp2Element(elem)

    def readCTETRA(self, data):
        """
        CTETRA(5508,55,217)    - the marker for Record 87
        """
        #print "reading CTETRA"
        n = 0
        nEntries = len(data) // 48  # 12*4
        for i in xrange(nEntries):
            eData = data[n:n + 48]
            out = unpack(b'iiiiiiiiiiii', eData)
            (eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10) = out
            #print "out = ",out

            dataIn = [eid, pid, n1, n2, n3, n4]
            bigNodes = [n5, n6, n7, n8, n9, n10]
            if sum(bigNodes) > 0:
                elem = CTETRA10(None, dataIn + bigNodes)
            else:
                elem = CTETRA4(None, dataIn)
            self.addOp2Element(elem)
            n += 48
        data = data[n:]

# CTETPR
# CTETR10F
# CTETR4FD
# CTQUAD
# CTTRIA

    def readCTRIA3(self, data):
        """
        CTRIA3(5959,59,282)    - the marker for Record 93
        """
        #print "reading CTRIA3"
        while len(data) >= 52:  # 13*4
            eData = data[:52]
            data = data[52:]
            out = unpack(b'iiiiiffiiifff', eData)
            #print "eid=%s pid=%s n1=%s n2=%s n3=%s theta=%s zoffs=%s blank1=%s blank2=%s tflag=%s t1=%s t2=%s t3=%s" %(eid,pid,n1,n2,n3,theta,zoffs,blank1,blank2,tflag,t1,t2,t3)
            (eid, pid, n1, n2, n3, theta, zoffs, blank1,
                blank2, tflag, t1, t2, t3) = out
            dataIn = [eid, pid, n1, n2, n3, theta, zoffs, tflag, t1, t2, t3]
            elem = CTRIA3(None, dataIn)
            self.addOp2Element(elem)

# CTRIAFD

    def readCTRIA6(self, data):
        """
        CTRIA6(4801,48,327)    - the marker for Record 95
        @warning inconsistent with dmap manual
        """
        #print "reading CTRIA6"
        return
        n = 0
        nEntries = len(data) // 52  # 13*4
        for i in xrange(nEntries):
            eData = data[n:n + 52]
            out = unpack(b'iiiiiiiiffffi', eData)
            #print "eid=%s pid=%s n1=%s n2=%s n3=%s theta=%s zoffs=%s blank1=%s blank2=%s tflag=%s t1=%s t2=%s t3=%s" %(eid,pid,n1,n2,n3,theta,zoffs,blank1,blank2,tflag,t1,t2,t3)
            (eid, pid, n1, n2, n3, n4, n5, n6, theta, t1, t2, t3, tflag) = out
            elem = CTRIA6(None, out)
            self.addOp2Element(elem)
            n += 52
        data = data[n:]

# CTRIA6FD
# CTRIAP

    def readCTRIAR(self, data):  # 98
        self.skippedCardsFile.write('skipping CTRIAR in GEOM2\n')

# CTRIAX

    def readCTRIAX6(self, data):  # 100
        self.skippedCardsFile.write('skipping CTRIAX6 in GEOM2\n')

# CTRIX3FD
# CTRIX6FD

    def readCTUBE(self, data):
        """
        CTUBE(3701,37,49)    - the marker for Record 103
        """
        #print "reading CTUBE"
        n = 0
        nEntries = len(data) // 16
        for i in xrange(nEntries):
            eData = data[n:n + 16]  # 4*4
            out = unpack(b'4i', eData)
            (eid, pid, n1, n2) = out
            elem = CTUBE(None, out)
            self.addOp2Element(elem)
            n += 16
        data = data[n:]

    def readCVISC(self, data):
        """CVISC(3901,39, 50) - the marker for Record 104"""
        #print "reading CVISC"
        n = 0
        nEntries = len(data) // 16
        for i in xrange(nEntries):
            eData = data[n:n + 16]  # 4*4
            out = unpack(b'4i', eData)
            #(eid,pid,n1,n2) = out
            elem = CVISC(None, out)
            self.addOp2Element(elem)
            n += 16
        data = data[n:]

    def readCWELD(self, data):  # 105
        self.skippedCardsFile.write('skipping CWELD in GEOM2\n')

    def readCWELDC(self, data):  # 106
        self.skippedCardsFile.write('skipping CWELDC in GEOM2\n')

    def readCWELDG(self, data):  # 107
        self.skippedCardsFile.write('skipping CWELDG in GEOM2\n')

# CWSEAM
# GENEL
# GMDNDC
# GMBNDS
# GMINTC
# GMINTS
# PLOTEL
# RADBC
# RADINT
# SINT

    def readSPOINT(self, data):
        """
        (5551,49,105)    - the marker for Record 118
        """
        while len(data) >= 4:  # 4*4
            eData = data[:4]
            data = data[4:]
            (nid) = unpack(b'i', eData)
            spoint = SPOINTs(None, [nid])
            self.addSPoint(spoint)

# VUBEAM
# VUHEXA
# VUQUAD4
# VUPENTA
# VUTETRA
# VUTRIA
# VUBEAM
# VUHEXA
# VUQUAD4
# WELDP
