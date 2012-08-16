import sys
from struct import unpack

#from pyNastran.bdf.cards.constraints import SPC,SPCADD
from pyNastran.bdf.cards.elements.rigid import RBE2
from pyNastran.bdf.cards.constraints import SUPORT, SPCD, SPC


class Geometry4(object):

    def readTable_Geom4(self):
        self.iTableMap = {
            (5561, 76, 215): self.readASET,    # record 1  - not done
            (5571, 77, 216): self.readASET1,   # record 2  - not done
            (10200, 102, 473): self.readBNDGRID,  # record 3  - not done
            (1510, 15, 328): self.readCYAX,    # record 8  - not done
            (5210, 52, 257): self.readCYJOIN,  # record 9  - not done
            (1710, 17, 330): self.readCYSYM,   # record 11 - not done
            (4901, 49, 17): self.readMPC,     # record 16 - not done
            (4891, 60, 83): self.readMPCADD,  # record 17 - not done
            (4951, 63, 92): self.readOMIT1,   # record 19 - not done
            (610, 6, 316): self.readQSET1,   # record 21 - not done
            (6601, 66, 292): self.readRBAR,    # record 22 - not done
            (6801, 68, 294): self.readRBE1,    # record 23 - not done
            (6901, 69, 295): self.readRBE2,    # record 24 - buggy
            (7101, 71, 187): self.readRBE3,    # record 25 - not done
            (6501, 65, 291): self.readRROD,    # record 30 - not done
            (7001, 70, 186): self.readRSPLINE,  # record 31 - not done
            (7201, 72, 398): self.readRSSCON,  # record 32 - not done
            (1210, 12, 322): self.readSEQSET1,  # record 40 - not done
            (5501, 55, 16): self.readSPC,     # record 44 - buggy
            (5481, 58, 12): self.readSPC1,    # record 45 - not done
            (5491, 59, 13): self.readSPCADD,  # record 46 - not done
            (5110, 51, 256): self.readSPCD,    # record 47 - buggy
            (5601, 56, 14): self.readSUPORT,  # record 59 - not done
            (10100, 101, 472): self.readSUPORT1,  # record 60 - not done
            #(4901,49,420017):self.readFake,    # record
            #(5561,76,0):     self.readFake,    # record
            #(5110,51,256):   self.readFake,    # record
            #(610,6,0):       self.readFake,    # record
            #(5501,55,620016):self.readFake,    # record

        }
        self.readRecordTable('GEOM4')

    def readTable_Geom4S(self):
        self.readTable_Geom4()
        #self.iTableMap = {
        #                 }
        #self.readRecordTable('GEOM4S')

    def readASET(self, data):
        """ASET(5561,76,215) - Record 1"""
        self.skippedCardsFile.write('skipping ASET in GEOM4\n')

    def readASET1(self, data):
        """ASET1(5571,77,216) - Record 2"""
        self.skippedCardsFile.write('skipping ASET1 in GEOM4\n')

    def readBNDGRID(self, data):
        """BNDGRID(10200,102,473) - Record 3 """
        self.skippedCardsFile.write('skipping BNDGRID in GEOM4\n')

# BSET
# BSET1
# CSET
# CSET1

    def readCYAX(self, data):
        """CYAX(1510,15,328) - Record 8 """
        self.skippedCardsFile.write('skipping CYAX in GEOM4\n')

    def readCYJOIN(self, data):
        """CYJOIN(5210,52,257) - Record 9 """
        self.skippedCardsFile.write('skipping CYJOIN in GEOM4\n')

# CYSUP
    def readCYSYM(self, data):
        """CYSYM(1710,17,330) - Record 11"""
        self.skippedCardsFile.write('skipping CYSYM in GEOM4\n')

# EGENDT
# GMBC
# GMSPC
    def readMPC(self, data):
        """MPC(4901,49,17) - Record 16"""
        self.skippedCardsFile.write('skipping MPC in GEOM4\n')

    def readMPCADD(self, data):
        """MPCADD(4891,60,83) - Record 17"""
        self.skippedCardsFile.write('skipping MPCADD in GEOM4\n')

    def readOMIT1(self, data):
        """OMIT1(4951,63,92) - Record 19"""
        self.skippedCardsFile.write('skipping OMIT1 in GEOM4\n')

    def readQSET1(self, data):
        """QSET1(610, 6, 316) - Record 21"""
        self.skippedCardsFile.write('skipping QSET1 in GEOM4\n')

    def readRBAR(self, data):
        """RBAR(6601,66,292) - Record 22"""
        self.skippedCardsFile.write('skipping RBAR in GEOM4\n')

    def readRBE1(self, data):
        """RBE1(6801,68,294) - Record 23"""
        self.skippedCardsFile.write('skipping RBE1 in GEOM4\n')

    def readRBE2(self, data):
        """RBE2(6901,69,295) - Record 24"""
        self.skippedCardsFile.write('skipping RBE2 in GEOM4\n')
        return
        #n=0
        #nData = len(data)  # 5*4
        if 1:
            eData = data[:12]
            (eid, gn, cm, gm) = unpack(b'iiii', eData)

            eData = data[12:-4]
            nGm = len(eData) // 4
            iFormat = 'i' * nGm
            iFormat = bytes(iFormat)
            Gm = list(unpack(iFormat, eData))
            alpha, = unpack(b'f', data[-4:])
        ###
        elem = RBE2(None, [eid, gn, cm, Gm, alpha])
        self.addRigidElement(elem)
        data = data[-1:]

    def readRBE3(self, data):
        """RBE3(7101,71,187) - Record 25"""
        self.skippedCardsFile.write('skipping RBE3 in GEOM4\n')

# RBJOINT
# RBJSTIF
# RELEASE
# RPNOM
    def readRROD(self, data):
        """RROD(6501,65,291) - Record 30"""
        self.skippedCardsFile.write('skipping RROD in GEOM4\n')

    def readRSPLINE(self, data):
        """RSPLINE(7001,70,186) - Record 31"""
        self.skippedCardsFile.write('skipping RSPLINE in GEOM4\n')

    def readRSSCON(self, data):
        """RSSCON(7201,72,398) - Record 32"""
        self.skippedCardsFile.write('skipping RSSCON in GEOM4\n')

# RTRPLT
# RWELD
# SEBSET
# SEBSET1
# SECSET
# SECSET1
# SEQSET

    def readSEQSET1(self, data):
        """SEQSET1(1210,12,322) - Record 40"""
        self.skippedCardsFile.write('skipping SEQSET1 in GEOM4\n')

# SESUP
# SEUSET
# SEUSET1

    def readSPC(self, data):
        """SPC(5501,55,16) - Record 44"""
        #self.skippedCardsFile.write('skipping SPC in GEOM4\n')
        n = 0
        nEntries = len(data) // 20  # 5*4
        for i in xrange(nEntries):
            eData = data[n:n + 20]
            (sid, ID, c, xxx, dx) = unpack(b'iiiif', eData)

            constraint = SPC(None, [sid, ID, c, dx])
            self.addConstraint_SPC(constraint)
            n += 20
        ###
        data = data[n:]

    def readSPC1(self, data):
        """SPC1(5481,58,12) - Record 45"""
        self.skippedCardsFile.write('skipping SPC1 in GEOM4\n')
        return
        n = 0
        nEntries = len(data) // 20  # 5*4
        for i in xrange(nEntries):
            eData = data[n:n + 20]
            (sid, c, thruFlag) = unpack(b'iifii', eData)

            constraint = SPC1(None, [sid, g, f, n1, n2])
            self.addConstraint_SPC(constraint)
            n += 20
        ###
        data = data[n:]

    def readSPCADD(self, data):
        """SPCADD(5491,59,13) - Record 46"""
        self.skippedCardsFile.write('skipping SPCADD in GEOM4\n')

    def readSPCD(self, data):
        """SPCD(5110,51,256) - Record 47"""
        #self.skippedCardsFile.write('skipping SPCD in GEOM4\n')
        n = 0
        nEntries = len(data) // 20  # 5*4
        for i in xrange(nEntries):
            eData = data[n:n + 20]
            (sid, ID, c, xxx, dx) = unpack(b'iiiif', eData)

            constraint = SPCD(None, [sid, ID, c, dx])
            self.addConstraint_SPC(constraint)
            n += 20
        ###
        data = data[n:]

    def readSPCDE(self, data):
        self.skippedCardsFile.write('skipping SPCDE in GEOM4\n')

    def readSPCDF(self, data):
        self.skippedCardsFile.write('skipping SPCDF in GEOM4\n')

    def readSPCDG(self, data):
        self.skippedCardsFile.write('skipping SPCDG in GEOM4\n')

    def readSPCE(self, data):
        self.skippedCardsFile.write('skipping SPCE in GEOM4\n')

    def readSPCEB(self, data):
        self.skippedCardsFile.write('skipping SPCEB in GEOM4\n')

    def readSPCF(self, data):
        self.skippedCardsFile.write('skipping SPCF in GEOM4\n')

    def readSPCFB(self, data):
        self.skippedCardsFile.write('skipping SPCFB in GEOM4\n')

    def readSPCGB(self, data):
        self.skippedCardsFile.write('skipping SPCGB in GEOM4\n')

# SPCGRID
# SPCOFF
# SPCOFF1

    def readSUPORT(self, data):
        """SUPORT(5601,56, 14) - Record 59"""
        #self.skippedCardsFile.write('skipping SUPORT in GEOM4\n')
        n = 0
        nEntries = len(data) // 8  # 2*4
        for i in xrange(nEntries):
            eData = data[n:n + 8]
            (sid, c) = unpack(b'ii', eData)

            suport = SUPORT(None, [sid, c])
            self.addSuport(suport)
            n += 8
        ###
        data = data[n:]

    def readSUPORT1(self, data):
        """SUPORT1(10100,101,472) - Record 60"""
        self.skippedCardsFile.write('skipping SUPORT1 in GEOM4\n')

# TEMPBC

    def readUSET(self, data):
        self.skippedCardsFile.write('skipping USET in GEOM4\n')

    def readUSET1(self, data):
        self.skippedCardsFile.write('skipping USET1 in GEOM4\n')
