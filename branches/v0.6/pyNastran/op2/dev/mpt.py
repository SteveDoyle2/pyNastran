#pylint: disable=C0111,C0103,C0301,W0612,W0613,R0914,R0201
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import StringIO
from struct import unpack, Struct

from pyNastran.bdf.cards.materials import (CREEP, MAT1, MAT2, MAT3, MAT4, MAT5,
                                           MAT8, MAT9, MAT10, MATHP)
from pyNastran.bdf.cards.material_deps import (MATT1, MATS1)
from pyNastran.bdf.cards.dynamic import NLPARM, TSTEP, TSTEPNL


class MPT(object):
    def add_TSTEPNL(self, card, allowOverwrites=True):
        raise RuntimeError('this should be overwritten')
    def add_NLPARM(self, card, allowOverwrites=True):
        raise RuntimeError('this should be overwritten')
    def add_material_dependence(self, material, allowOverwrites=True):
        raise RuntimeError('this should be overwritten')
    def add_creep_material(self, material, allowOverwrites=True):
        raise RuntimeError('this should be overwritten')
    def add_structural_material(self, material, allowOverwrites=True):
        raise RuntimeError('this should be overwritten')
    def add_thermal_material(self, material, allowOverwrites=True):
        raise RuntimeError('this should be overwritten')

    def __init__(self):
        self.skippedCardsFile = StringIO.StringIO()
        self.card_count = {}
        self.bigMaterials = {}
        self._mpt_map = {
            (1003, 10, 245): ['CREEP', self.readCREEP],  # record 1
            ( 103,  1,  77): ['MAT1', self.readMAT1],    # record 2
            ( 203,  2,  78): ['MAT2', self.readMAT2],    # record 3
            (1403, 14, 122): ['MAT3', self.readMAT3],    # record 4
            (2103, 21, 234): ['MAT4', self.readMAT4],    # record 5
            (2203, 22, 235): ['MAT5', self.readMAT5],    # record 6
            (2503, 25, 288): ['MAT8', self.readMAT8],    # record 7
            (2603, 26, 300): ['MAT9', self.readMAT9],    # record 8 - buggy
            (2801, 28, 365): ['MAT10', self.readMAT10],  # record 9
            (4506, 45, 374): ['MATHP', self.readMATHP],  # record 11
            (503,  5,   90): ['MATS1', self.readMATS1],  # record 12
            (703,   7,  91): ['MATT1',   self.readMATT1],   # record 13 - not done
            (803,   8, 102): ['MATT2',   self.readMATT2],   # record 14 - not done
            (1503, 14, 189): ['MATT3',   self.readMATT3],   # record 15 - not done
            (2303, 23, 237): ['MATT4',   self.readMATT4],   # record 16 - not done
            (2403, 24, 238): ['MATT5',   self.readMATT5],   # record 17 - not done
            (2703, 27, 301): ['MATT9',   self.readMATT9],   # record 19 - not done
            (8802, 88, 413): ['RADM',    self.readRADM],    # record 25 - not done
            # record 26
            (3003, 30, 286): ['NLPARM',  self.readNLPARM],  # record 27
            (3104, 32, 350): ['NLPCI',   self.readNLPCI],   # record 28
            (3103, 31, 337): ['TSTEPNL', self.readTSTEPNL], # record 29

            (903, 9, 336): ['', self.readFake],
            (8902, 89, 423): ['', self.readFake],
            (9002, 90, 410): ['', self.readFake],

        }

    def addOp2Material(self, mat):
        self.add_structural_material(mat, allowOverwrites=True)
        #print(str(mat)[:-1])

    def readCREEP(self, data, n):
        """
        CREEP(1003,10,245) - record 1
        """
        #print "reading CREEP"
        nmaterials = (len(data) - n) // 64
        s = Struct(b'i2f4ifi7f')
        for i in xrange(nmaterials):
            edata = data[n:n+64]
            out = s.unpack(edata)
            (mid, T0, exp, form, tidkp, tidcp, tidcs, thresh,
                Type, ag1, ag2, ag3, ag4, ag5, ag6, ag7) = out
            self.add_creep_material(CREEP(None, out), allowOverwrites=True)
            n += 64
        self.card_count['CREEP'] = nmaterials
        return n

    def readMAT1(self, data, n):
        """
        MAT1(103,1,77) - record 2
        """
        #print "reading MAT1"
        ntotal = 48  # 12*4
        s = Struct(b'i10fi')
        nmaterials = (len(data) - n) // ntotal
        for i in xrange(nmaterials):
            eData = data[n:n+48]
            out = s.unpack(eData)
            (mid, E, G, nu, rho, A, TRef, ge, St, Sc, Ss, mcsid) = out
            self.addOp2Material(MAT1(None, out))
            n += ntotal
        self.card_count['MAT1'] = nmaterials
        return n

    def readMAT2(self, data, n):
        """
        MAT2(203,2,78) - record 3
        """
        #print "reading MAT2"
        ntotal = 68  # 17*4
        s = Struct(b'i15fi')
        nmaterials = (len(data) - n) // ntotal
        for i in xrange(nmaterials):
            edata = data[n:n+68]
            out = s.unpack(edata)
            (mid, g1, g2, g3, g4, g5, g6, rho, aj1, aj2, aj3,
                TRef, ge, St, Sc, Ss, mcsid) = out
            #print "MAT2 = ",out
            mat = MAT2(None, out)

            if mid > 1e8 or mid < 0:  # just a checker for out of range materials
                self.bigMaterials[mid] = mat
            else:
                self.addOp2Material(mat)
            n += ntotal
        self.card_count['MAT2'] = nmaterials
        return n

    def readMAT3(self, data, n):
        """
        MAT3(1403,14,122) - record 4
        """
        #print "reading MAT3"
        s = Struct(b'i8fi5fi')
        nmaterials = (len(data) - n) // 64
        for i in xrange(nmaterials):
            out = s.unpack(data[n:n+64])
            (mid, ex, eth, ez, nuxth, nuthz, nuzx, rho, gzx,
                blank, ax, ath, az, TRef, ge, blank) = out
            mat = MAT3(None, [mid, ex, eth, ez, nuxth, nuthz,
                              nuzx, rho, gzx, ax, ath, az, TRef, ge])
            self.addOp2Material(mat)
            n += 64
        self.card_count['MAT3'] = nmaterials

    def readMAT4(self, data, n):
        """
        MAT4(2103,21,234) - record 5
        """
        #print "reading MAT4"
        s = Struct(b'i10f')
        nmaterials = (len(data) - n) // 40
        for i in xrange(nmaterials):
            out = s.unpack(data[n:n+44])
            (mid, k, cp, rho, h, mu, hgen, refenth, tch, tdelta, qlat) = out
            self.add_thermal_material(MAT4(None, out), allowOverwrites=True)
            n += 44
        self.card_count['MAT4'] = nmaterials

    def readMAT5(self, data, n):
        """
        MAT5(2203,22,235) - record 6
        """
        #print "reading MAT5"
        s = Struct(b'i9f')
        nmaterials = (len(data) - n) // 40
        for i in xrange(nmaterials):
            out = s.unpack(data[n:n+40])
            (mid, k1, k2, k3, k4, k5, k6, cp, rho, hgen) = out
            self.add_thermal_material(MAT5(None, out), allowOverwrites=True)
            n += 40
        self.card_count['MAT5'] = nmaterials

    def readMAT8(self, data, n):
        """
        MAT8(2503,25,288) - record 7
        """
        #print "reading MAT8"
        s = Struct(b'i18f')
        nmaterials = (len(data) - n) // 76
        for i in xrange(nmaterials):
            out = s.unpack(data[n:n+76])
            (mid, E1, E2, nu12, G12, G1z, G2z, rho, a1, a2,
             TRef, Xt, Xc, Yt, Yc, S, ge, f12, strn) = out
            self.addOp2Material(MAT8(None, out))
            n += 76
        self.card_count['MAT8'] = nmaterials

    def readMAT9(self, data, n):
        """
        MAT9(2603,26,300) - record 9
        .. todo:: buggy
        """
        #print "reading MAT9"
        s = Struct(b'22i9f4i')
        nmaterials = (len(data) - n) // 140
        for i in xrange(nmaterials):
            out = s.unpack(data[n:n+140])
            (mid, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12, g13, g14, g15, g16, g17, g18, g19, g20, g21,
             rho, a1, a2, a3, a4, a5, a6, TRef, ge, blank1, blank2, blank3, blank4) = out
            dataIn = [mid, [g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12, g13, g14, g15, g16, g17, g18, g19, g20, g21],
                      rho, [a1, a2, a3, a4, a5, a6],
                      TRef, ge]
            #print "dataIn = ",dataIn
            self.addOp2Material(MAT9(None, dataIn))
            n += 140
        self.card_count['MAT9'] = nmaterials

    def readMAT10(self, data, n):
        """
        MAT10(2801,28,365) - record 9
        """
        #print "reading MAT10"
        ntotal = 44  # 5*4
        s = Struct(b'i4f')
        nmaterials = (len(data) - n) // ntotal
        for i in xrange(nmaterials):
            edata = data[n:n+20]
            out = s.unpack(edata)
            (mid, bulk, rho, c, ge) = out
            self.addOp2Material(MAT10(None, out))
        self.card_count['MAT10'] = nmaterials

# MAT11 - unused

    def readMATHP(self, data, n):
        """MATHP(4506,45,374) - Record 11"""
        #print "reading MATHP"
        nmaterials = 0
        s1 = Struct(b'i7f3i23fi')
        s2 = Struct(b'8i')
        n2 = n
        while n2 < n:
            eData = data[n:n+140]
            n += 140
            out1 = s1.unpack(eData)
            (mid, a10, a01, d1, rho, alpha, tref, ge, sf, na, nd, kp,
             a20, a11, a02, d2,
             a30, a21, a12, a03, d3,
             a40, a31, a22, a13, a04, d4,
             a50, a41, a32, a23, a14, a05, d5,
             continueFlag) = out1
            dataIn = [out1]

            if continueFlag:
                eData = data[n:n+32]  # 7*4
                n += 32
                out2 = s2.unpack(eData)
                (tab1, tab2, tab3, tab4, x1, x2, x3, tab5) = out2
                data.append(out2)
            self.addOp2Material(MATHP(None, dataIn))
            nmaterials += 1
        self.card_count['MATHP'] = nmaterials

    def readMATS1(self, data, n):
        """
        MATS1(503,5,90) - record 12
        """
        #print "reading MATS1"
        ntotal = 44  # 11*4
        s = Struct(b'3ifiiff3i')
        nmaterials = (len(data) - n) // ntotal
        for i in xrange(nmaterials):
            edata = data[n:n+44]
            out = s.unpack(edata)
            (mid, tid, Type, h, yf, hr, limit1, limit2, a, b, c) = out
            dataIn = [mid, tid, Type, h, yf, hr, limit1, limit2]
            self.add_material_dependence(MATS1(None, dataIn), allowOverwrites=True)
        self.card_count['MATS1'] = nmaterials
        return n

    def readMATT1(self, data, n):
        self.skippedCardsFile.write('skipping MATT1 in MPT\n')
        return n

    def readMATT2(self, data, n):
        self.skippedCardsFile.write('skipping MATT2 in MPT\n')
        return n

    def readMATT3(self, data, n):
        self.skippedCardsFile.write('skipping MATT3 in MPT\n')
        return n

    def readMATT4(self, data, n):
        self.skippedCardsFile.write('skipping MATT4 in MPT\n')
        return n

    def readMATT5(self, data, n):
        self.skippedCardsFile.write('skipping MATT5 in MPT\n')
        return n

# MATT8 - unused
    def readMATT9(self, data, n):
        self.skippedCardsFile.write('skipping MATT9 in MPT\n')
        return n

# MBOLT
# MBOLTUS
# MSTACK
# NLAUTO
# RADBND

    def readRADM(self, data, n):
        """
        RADM(8802,88,413) - record 25
        .. todo:: add object
        """
        #print "reading RADM"
        return
        s = Struct(b'i', )
        while len(data) >= 4:  # 1*4
            eData = data[:4]
            data = data[4:]
            number, = s.unpack(eData)

            iFormat = 'if%if' % (number + 1)
            eDataLen = len(strings) * 4

            eData = data[:eDataLen]
            data = data[eDataLen:]
            iFormat = bytes(iFormat)
            pack = list(unpack(iFormat, eData))
            packs = []

            while data:
                eData = data[:eDataLen]
                data = data[eDataLen:]
                pack = list(unpack(iFormat, eData))
                packs.append(pack)

            #mat = RADM(None, packs)
            #self.addOp2Material(mat)

# RADMT

    def readNLPARM(self, data, n):
        """
        NLPARM(3003,30,286) - record 27
        """
        #print "reading NLPARM"
        ntotal = 76  # 19*4
        s = Struct(b'iif5i3f3iffiff')
        nentries = (len(data) - n) // ntotal
        for i in xrange(nentries):
            edata = data[n:n+76]
            out = s.unpack(edata)
            #(sid,ninc,dt,kMethod,kStep,maxIter,conv,intOut,epsU,epsP,epsW,
            # maxDiv,maxQn,maxLs,fStress,lsTol,maxBisect,maxR,rTolB) = out
            self.add_NLPARM(NLPARM(None, out))
            n += ntotal
        self.card_count['NLPARM'] = nentries
        return n

    def readNLPCI(self, data, n):
        self.skippedCardsFile.write('skipping NLPCI in MPT\n')
        return n

    def readTSTEPNL(self, data, n):
        """
        TSTEPNL(3103,31,337) - record 29
        """
        #print "reading TSTEPNL"
        ntotal = 88  # 19*4
        s = Struct(b'iif5i3f3if3i4f')
        nentries = (len(data) - n) // ntotal
        for i in xrange(nentries):
            edata = data[n:n+88]
            out = s.unpack(edata)
            #(sid,ndt,dt,no,kMethod,kStep,maxIter,conv,epsU,epsP,epsW,
            # maxDiv,maxQn,maxLs,fStress,lsTol,maxBisect,adjust,mStep,rb,maxR,uTol,rTolB) = out
            self.add_TSTEPNL(TSTEPNL(None, out))
            n += ntotal
        self.card_count['TSTEPNL'] = nentries
        return n
