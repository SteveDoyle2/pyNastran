from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import sys
from struct import unpack

from pyNastran.bdf.cards.materials import (CREEP, MAT1, MAT2, MAT3, MAT4, MAT5,
                                           MAT8, MAT9, MAT10, MATS1, MATHP)
from pyNastran.bdf.cards.dynamic import NLPARM, TSTEP, TSTEPNL


class MPT(object):

    def readTable_MPTS(self):
        self.tableName = 'MPT'
        self.bigMaterials = {}
        self.iTableMap = {
            (1003, 10, 245): self.readCREEP,  # record 1
            (103, 1, 77): self.readMAT1,  # record 2
            (203, 2, 78): self.readMAT2,  # record 3
            (1403, 14, 122): self.readMAT3,  # record 4
            (2103, 21, 234): self.readMAT4,  # record 5
            (2203, 22, 235): self.readMAT5,  # record 6
            (2503, 25, 288): self.readMAT8,  # record 7
            (2603, 26, 300): self.readMAT9,  # record 8 - buggy
            (2801, 28, 365): self.readMAT10,  # record 9
            (4506, 45, 374): self.readMATHP,  # record 11
            (503, 5, 90): self.readMATS1,  # record 12

            (3003, 30, 286): self.readNLPARM,  # record 27
            (3104, 32, 350): self.readNLPCI,   # record 28
            (3103, 31, 337): self.readTSTEPNL,  # record 29
            (703, 7, 91): self.readMATT1,   # record 13 - not done
            (803, 8, 102): self.readMATT2,   # record 14 - not done
            (1503, 14, 189): self.readMATT3,   # record 15 - not done
            (2303, 23, 237): self.readMATT4,   # record 16 - not done
            (2403, 24, 238): self.readMATT5,   # record 17 - not done
            (2703, 27, 301): self.readMATT9,   # record 19 - not done
            #(8802,88,413): self.readRADM,    # record 25 - not done

        }

        self.readRecordTable('MPTS')

    def addOp2Material(self, mat):
        self.addMaterial(mat, allowOverwrites=True)

    def readCREEP(self, data):
        """
        CREEP(1003,10,245) - record 1
        """
        #print "reading CREEP"
        while len(data) >= 64:  # 16*4
            eData = data[:64]
            data = data[64:]
            out = unpack(b'iffiiiififffffff', eData)
            (mid, T0, exp, form, tidkp, tidcp, tidcs, thresh,
                Type, ag1, ag2, ag3, ag4, ag5, ag6, ag7) = out
            mat = CREEP(None, out)
            self.addCreepMaterial(mat, allowOverwrites=True)

    def readMAT1(self, data):
        """
        MAT1(103,1,77) - record 2
        """
        #print "reading MAT1"
        while len(data) >= 48:  # 12*4
            eData = data[:48]
            data = data[48:]
            out = unpack(b'iffffffffffi', eData)
            (mid, E, G, nu, rho, A, TRef, ge, St, Sc, Ss, mcsid) = out
            mat = MAT1(None, out)
            self.addOp2Material(mat)

    def readMAT2(self, data):
        """
        MAT2(203,2,78) - record 3
        """
        #print "reading MAT2"
        while len(data) >= 68:  # 17*4
            eData = data[:68]
            data = data[68:]
            out = unpack(b'ifffffffffffffffi', eData)
            (mid, g1, g2, g3, g4, g5, g6, rho, aj1, aj2, aj3,
                TRef, ge, St, Sc, Ss, mcsid) = out
            #print "MAT2 = ",out
            mat = MAT2(None, out)

            if mid > 1e8 or mid < 0:  # just a checker for out of range materials
                self.bigMaterials[mid] = mat
            else:
                self.addOp2Material(mat)

    def readMAT3(self, data):
        """
        MAT3(1403,14,122) - record 4
        """
        #print "reading MAT3"
        while len(data) >= 64:  # 16*4
            eData = data[:64]
            data = data[64:]
            out = unpack(b'iffffffffifffffi', eData)
            (mid, ex, eth, ez, nuxth, nuthz, nuzx, rho, gzx,
                blank, ax, ath, az, TRef, ge, blank) = out
            mat = MAT3(None, [mid, ex, eth, ez, nuxth, nuthz,
                              nuzx, rho, gzx, ax, ath, az, TRef, ge])
            self.addOp2Material(mat)

    def readMAT4(self, data):
        """
        MAT4(2103,21,234) - record 5
        """
        #print "reading MAT4"
        while len(data) >= 44:  # 11*4
            eData = data[:44]
            data = data[44:]
            out = unpack(b'i10f', eData)
            (mid, k, cp, rho, h, mu, hgen, refenth, tch, tdelta, qlat) = out
            mat = MAT4(None, out)
            self.addThermalMaterial(mat, allowOverwrites=True)

    def readMAT5(self, data):
        """
        MAT5(2203,22,235) - record 6
        """
        #print "reading MAT5"
        while len(data) >= 40:  # 10*4
            eData = data[:40]
            data = data[40:]
            out = unpack(b'ifffffffff', eData)
            (mid, k1, k2, k3, k4, k5, k6, cp, rho, hgen) = out
            mat = MAT5(None, out)
            self.addThermalMaterial(mat, allowOverwrites=True)

    def readMAT8(self, data):
        """
        MAT8(2503,25,288) - record 7
        """
        #print "reading MAT8"
        while len(data) >= 76:  # 19*4
            eData = data[:76]
            data = data[76:]
            out = unpack(b'i18f', eData)
            (mid, E1, E2, nu12, G12, G1z, G2z, rho, a1, a2,
                TRef, Xt, Xc, Yt, Yc, S, ge, f12, strn) = out
            mat = MAT8(None, out)
            self.addOp2Material(mat)

    def readMAT9(self, data):
        """
        MAT9(2603,26,300) - record 9
        @todo buggy
        """
        #print "reading MAT9"
        while len(data) >= 140:  # 35*4
            eData = data[:140]
            data = data[140:]
            out = unpack(b'iiiiiiiiiiiiiiiiiiiiiifffffffffiiii', eData)

            (
                mid, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12, g13, g14, g15, g16, g17, g18, g19, g20, g21,
                rho, a1, a2, a3, a4, a5, a6, TRef, ge, blank1, blank2, blank3, blank4) = out
            dataIn = [mid, [g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12, g13, g14, g15, g16, g17, g18, g19, g20, g21],
                      rho, [a1, a2, a3, a4, a5, a6],
                      TRef, ge]
            #print "dataIn = ",dataIn
            mat = MAT9(None, dataIn)
            self.addOp2Material(mat)

    def readMAT10(self, data):
        """
        MAT10(2801,28,365) - record 9
        """
        #print "reading MAT10"
        while len(data) >= 20:  # 5*4
            eData = data[:20]
            data = data[20:]
            out = unpack(b'iffff', eData)
            (mid, bulk, rho, c, ge) = out
            mat = MAT10(None, out)
            self.addOp2Material(mat)

# MAT11 - unused

    def readMATHP(self, data):
        """MATHP(4506,45,374) - Record 11"""
        #print "reading MATHP"
        while len(data) >= 140:  # 35*4
            eData = data[:140]
            data = data[140:]
            out1 = unpack(b'ifffffffiiifffffffffffffffffffffffi', eData)
            (mid, a10, a01, d1, rho, alpha, tref, ge, sf, na, nd, kp,
             a20, a11, a02, d2,
             a30, a21, a12, a03, d3,
             a40, a31, a22, a13, a04, d4,
             a50, a41, a32, a23, a14, a05, d5,
             continueFlag) = out1
            dataIn = [out1]

            if continueFlag:
                eData = data[:32]  # 7*4
                data = data[32:]
                out2 = unpack(b'iiiiiiii', eData)
                (tab1, tab2, tab3, tab4, x1, x2, x3, tab5) = out2
                data.append(out2)
            mat = MATHP(None, dataIn)
            self.addOp2Material(mat)

    def readMATS1(self, data):
        """
        MATS1(503,5,90) - record 12
        """
        #print "reading MATS1"
        while len(data) >= 44:  # 11*4
            eData = data[:44]
            data = data[44:]
            out = unpack(b'iiifiiffiii', eData)
            (mid, tid, Type, h, yf, hr, limit1, limit2, a, b, c) = out
            dataIn = [mid, tid, Type, h, yf, hr, limit1, limit2]
            mat = MATS1(None, dataIn)
            self.addMaterialDependence(mat, allowOverwrites=True)

    def readMATT1(self, data):
        self.skippedCardsFile.write('skipping MATT1 in MPT\n')

    def readMATT2(self, data):
        self.skippedCardsFile.write('skipping MATT2 in MPT\n')

    def readMATT3(self, data):
        self.skippedCardsFile.write('skipping MATT3 in MPT\n')

    def readMATT4(self, data):
        self.skippedCardsFile.write('skipping MATT4 in MPT\n')

    def readMATT5(self, data):
        self.skippedCardsFile.write('skipping MATT5 in MPT\n')

# MATT8 - unused

    def readMATT9(self, data):
        self.skippedCardsFile.write('skipping MATT9 in MPT\n')

# MBOLT
# MBOLTUS
# MSTACK
# NLAUTO
# RADBND

    def readRADM(self, data):
        """
        RADM(8802,88,413) - record 25
        @todo add object
        """
        #print "reading RADM"
        return
        while len(data) >= 4:  # 1*4
            eData = data[:4]
            data = data[4:]
            number, = unpack(b'i', eData)

            iFormat = 'if' + str(number) + 'f'
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

            #mat = RADM(None,packs)
            #self.addOp2Material(mat)

# RADMT

    def readNLPARM(self, data):
        """
        NLPARM(3003,30,286) - record 27
        """
        #print "reading NLPARM"
        while len(data) >= 76:  # 19*4
            eData = data[:76]
            data = data[76:]
            out = unpack(b'iifiiiiifffiiiffiff', eData)
            #(sid,ninc,dt,kMethod,kStep,maxIter,conv,intOut,epsU,epsP,epsW,
            # maxDiv,maxQn,maxLs,fStress,lsTol,maxBisect,maxR,rTolB) = out
            nlparm = NLPARM(None, out)
            self.addNLParm(nlparm)

    def readNLPCI(self, data):
        self.skippedCardsFile.write('skipping NLPCI in MPT\n')

    def readTSTEPNL(self, data):
        """
        TSTEPNL(3103,31,337) - record 29
        """
        #print "reading TSTEPNL"
        while len(data) >= 88:  # 19*4
            eData = data[:88]
            data = data[88:]
            out = unpack(b'iifiiiiifffiiifiiiffff', eData)
            #(sid,ndt,dt,no,kMethod,kStep,maxIter,conv,epsU,epsP,epsW,
            # maxDiv,maxQn,maxLs,fStress,lsTol,maxBisect,adjust,mStep,rb,maxR,uTol,rTolB) = out
            tstepnl = TSTEPNL(None, out)
            self.addTSTEPNL(tstepnl)
