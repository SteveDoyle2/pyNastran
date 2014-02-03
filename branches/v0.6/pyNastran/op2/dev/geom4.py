#pylint: disable=C0301,C0111,C0103,W0613
import StringIO
from struct import unpack, Struct

#from pyNastran.bdf.cards.constraints import SPC,SPCADD
from pyNastran.bdf.cards.elements.rigid import RBE2
from pyNastran.bdf.cards.constraints import SUPORT, SPCD, SPC, SPC1


class GEOM4(object):
    def add_constraint_SPC(self, constraint):
        raise RuntimeError('this should be overwritten')
    def add_rigid_element(self, constraint):
        raise RuntimeError('this should be overwritten')
    def add_suport(self, constraint):
        raise RuntimeError('this should be overwritten')
    def readFake(self, data, n):
        raise RuntimeError('this should be overwritten')

    def __init__(self):
        self.skippedCardsFile = StringIO.StringIO()
        self.card_count = {}
        self._geom4_map = {
            (5561,   76, 215): ['ASET', self.readASET],      # record 1  - not done
            (5571,   77, 216): ['ASET1', self.readASET1],    # record 2  - not done
            (10200, 102, 473): ['BNDGRID', self.readBNDGRID],# record 3  - not done
            (1510,   15, 328): ['CYAX', self.readCYAX],      # record 8  - not done
            (5210,   52, 257): ['CYJOIN', self.readCYJOIN],  # record 9  - not done
            (1710,   17, 330): ['CYSYM', self.readCYSYM],    # record 11 - not done
            (4901,   49,  17): ['MPC', self.readMPC],        # record 16 - not done
            (4891,   60,  83): ['MPCADD', self.readMPCADD],  # record 17 - not done
            (4951,   63,  92): ['OMIT1', self.readOMIT1],    # record 19 - not done

            (610,     6, 316): ['QSET1', self.readQSET1],    # record 21 - not done
            (6601,   66, 292): ['RBAR', self.readRBAR],      # record 22 - not done
            (6801,   68, 294): ['RBE1', self.readRBE1],      # record 23 - not done
            (6901,   69, 295): ['RBE2', self.readRBE2],      # record 24 - buggy
            (7101,   71, 187): ['RBE3', self.readRBE3],      # record 25 - not done

            (6501,   65, 291): ['RROD', self.readRROD],        # record 30 - not done
            (7001,   70, 186): ['RSPLINE', self.readRSPLINE],  # record 31 - not done
            (7201,   72, 398): ['RSSCON', self.readRSSCON],    # record 32 - not done
            (1210,   12, 322): ['SEQSET1', self.readSEQSET1],  # record 40 - not done
            (5501,   55,  16): ['SPC', self.readSPC],          # record 44 - buggy
            (5481,   58,  12): ['SPC1', self.readSPC1],        # record 45 - not done
            (5491,   59,  13): ['SPCADD', self.readSPCADD],    # record 46 - not done
            (5110,   51, 256): ['SPCD', self.readSPCD],        # record 47 - buggy
            (5601,   56,  14): ['SUPORT', self.readSUPORT],    # record 59 - not done
            (10100, 101, 472): ['SUPORT1', self.readSUPORT1],  # record 60 - not done

            (1310, 13,    247): ['', self.readFake],    # record
            (4901, 49, 420017): ['', self.readFake],    # record
            (5561, 76,      0): ['', self.readFake],     # record
            (5110, 51,    256): ['', self.readFake],     # record
            (610,   6,      0): ['', self.readFake],     # record
            (5110, 51, 620256): ['', self.readFake],    # record
            (5501, 55, 620016): ['', self.readFake],    # record
            (5001, 50, 15): ['', self.readFake],    # record
            (2010, 20, 193): ['', self.readFake],    # record
            (410, 4, 0): ['', self.readFake],    # record
            (110, 1, 584): ['', self.readFake],    # record
            (210, 2, 585): ['', self.readFake],    # record
            (6210, 62, 344): ['', self.readFake],    # record
            (510, 5, 315): ['', self.readFake],    # record
            (6701, 67, 293): ['', self.readFake],    # record
            (8801, 88, 9022): ['', self.readFake],    # record
            (9001, 90, 9024): ['', self.readFake],    # record
            (9901, 99, 80): ['', self.readFake],  # record
            (1010, 10, 320): ['', self.readFake],  # record
            (9801, 98, 79): ['', self.readFake],  # record
        }

    def readASET(self, data, n):
        """ASET(5561,76,215) - Record 1"""
        self.skippedCardsFile.write('skipping ASET in GEOM4\n')
        return n

    def readASET1(self, data, n):
        """ASET1(5571,77,216) - Record 2"""
        self.skippedCardsFile.write('skipping ASET1 in GEOM4\n')
        return n

    def readBNDGRID(self, data, n):
        """BNDGRID(10200,102,473) - Record 3 """
        self.skippedCardsFile.write('skipping BNDGRID in GEOM4\n')
        return n

# BSET
# BSET1
# CSET
# CSET1

    def readCYAX(self, data, n):
        """CYAX(1510,15,328) - Record 8 """
        self.skippedCardsFile.write('skipping CYAX in GEOM4\n')
        return n

    def readCYJOIN(self, data, n):
        """CYJOIN(5210,52,257) - Record 9 """
        self.skippedCardsFile.write('skipping CYJOIN in GEOM4\n')
        return n

# CYSUP
    def readCYSYM(self, data, n):
        """CYSYM(1710,17,330) - Record 11"""
        self.skippedCardsFile.write('skipping CYSYM in GEOM4\n')
        return n

# EGENDT
# GMBC
# GMSPC
    def readMPC(self, data, n):
        """MPC(4901,49,17) - Record 16"""
        self.skippedCardsFile.write('skipping MPC in GEOM4\n')
        return n

    def readMPCADD(self, data, n):
        """MPCADD(4891,60,83) - Record 17"""
        self.skippedCardsFile.write('skipping MPCADD in GEOM4\n')
        return n

    def readOMIT1(self, data, n):
        """OMIT1(4951,63,92) - Record 19"""
        self.skippedCardsFile.write('skipping OMIT1 in GEOM4\n')
        return n

    def readQSET1(self, data, n):
        """QSET1(610, 6, 316) - Record 21"""
        self.skippedCardsFile.write('skipping QSET1 in GEOM4\n')
        return n

    def readRBAR(self, data, n):
        """RBAR(6601,66,292) - Record 22"""
        self.skippedCardsFile.write('skipping RBAR in GEOM4\n')
        return n

    def readRBE1(self, data, n):
        """RBE1(6801,68,294) - Record 23"""
        self.skippedCardsFile.write('skipping RBE1 in GEOM4\n')
        return n

    def readRBE2(self, data, n):
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
        elem = RBE2(None, [eid, gn, cm, Gm, alpha])
        self.add_rigid_element(elem)
        data = data[-1:]

    def readRBE3(self, data, n):
        """RBE3(7101,71,187) - Record 25"""
        self.skippedCardsFile.write('skipping RBE3 in GEOM4\n')
        return n

# RBJOINT
# RBJSTIF
# RELEASE
# RPNOM
    def readRROD(self, data, n):
        """RROD(6501,65,291) - Record 30"""
        self.skippedCardsFile.write('skipping RROD in GEOM4\n')
        return n

    def readRSPLINE(self, data, n):
        """RSPLINE(7001,70,186) - Record 31"""
        self.skippedCardsFile.write('skipping RSPLINE in GEOM4\n')
        return n

    def readRSSCON(self, data, n):
        """RSSCON(7201,72,398) - Record 32"""
        self.skippedCardsFile.write('skipping RSSCON in GEOM4\n')
        return n

# RTRPLT
# RWELD
# SEBSET
# SEBSET1
# SECSET
# SECSET1
# SEQSET

    def readSEQSET1(self, data, n):
        """SEQSET1(1210,12,322) - Record 40"""
        self.skippedCardsFile.write('skipping SEQSET1 in GEOM4\n')
        return n

# SESUP
# SEUSET
# SEUSET1

    def readSPC(self, data, n):
        """SPC(5501,55,16) - Record 44"""
        #self.skippedCardsFile.write('skipping SPC in GEOM4\n')
        n = 0
        nEntries = len(data) // 20  # 5*4
        for i in xrange(nEntries):
            eData = data[n:n + 20]
            (sid, ID, c, xxx, dx) = unpack(b'iiiif', eData)

            constraint = SPC(None, [sid, ID, c, dx])
            self.add_constraint_SPC(constraint)
            n += 20
        return n

    def readSPC1(self, data, n):
        """SPC1(5481,58,12) - Record 45"""
        self.skippedCardsFile.write('skipping SPC1 in GEOM4\n')
        #return n
        n2 = n
        #nentries = (len(data) - n - 12) // 4  # 5*4
        nentries = 0
        while n2 < n:
            eData = data[n:n+12]
            n += 12
            (sid, g, thru_flag, n1) = unpack('4i', eData)
            eData = data[n:n + 12]

            nids = [n1]
            if thru_flag == 0:  # repeat 4 to end
                nnodes = (len(data) - n) // 4
                nodes = unpack(b'%ii' % nnodes, data[n:])
                nids += list(nodes)
                n += 4 * nentries
            else:
                n2 = unpack(b'i', data[n:n+4])
                n += 4
                nids.append(n2)
            nentries += 1
            constraint = SPC1(None, [sid, g, nids])
            self.add_constraint_SPC(constraint)
        self.card_count['SPC1'] = nentries
        return n

    def readSPCADD(self, data, n):
        """SPCADD(5491,59,13) - Record 46"""
        self.skippedCardsFile.write('skipping SPCADD in GEOM4\n')
        return n

    def readSPCD(self, data, n):
        """SPCD(5110,51,256) - Record 47"""
        #self.skippedCardsFile.write('skipping SPCD in GEOM4\n')
        n = 0
        s = Struct(b'4if')
        nEntries = len(data) // 20  # 5*4
        for i in xrange(nEntries):
            eData = data[n:n + 20]
            (sid, ID, c, xxx, dx) = s.unpack(eData)

            constraint = SPCD(None, [sid, ID, c, dx])
            self.add_constraint_SPC(constraint)
            n += 20
        self.card_count['SPCD'] = nEntries
        return n

    def readSPCDE(self, data, n):
        self.skippedCardsFile.write('skipping SPCDE in GEOM4\n')
        return n

    def readSPCDF(self, data, n):
        self.skippedCardsFile.write('skipping SPCDF in GEOM4\n')
        return n

    def readSPCDG(self, data, n):
        self.skippedCardsFile.write('skipping SPCDG in GEOM4\n')
        return n

    def readSPCE(self, data, n):
        self.skippedCardsFile.write('skipping SPCE in GEOM4\n')
        return n

    def readSPCEB(self, data, n):
        self.skippedCardsFile.write('skipping SPCEB in GEOM4\n')
        return n

    def readSPCF(self, data, n):
        self.skippedCardsFile.write('skipping SPCF in GEOM4\n')
        return n

    def readSPCFB(self, data, n):
        self.skippedCardsFile.write('skipping SPCFB in GEOM4\n')
        return n

    def readSPCGB(self, data, n):
        self.skippedCardsFile.write('skipping SPCGB in GEOM4\n')
        return n

# SPCGRID
# SPCOFF
# SPCOFF1

    def readSUPORT(self, data, n):
        """SUPORT(5601,56, 14) - Record 59"""
        #self.skippedCardsFile.write('skipping SUPORT in GEOM4\n')
        n = 0
        nEntries = len(data) // 8  # 2*4
        s = Struct(b'2i')
        for i in xrange(nEntries):
            self.add_suport(SUPORT(None, list(s.unpack(data[n:n + 8])))) # extracts [sid, c]
            n += 8
        return n

    def readSUPORT1(self, data, n):
        """SUPORT1(10100,101,472) - Record 60"""
        self.skippedCardsFile.write('skipping SUPORT1 in GEOM4\n')
        return n

# TEMPBC

    def readUSET(self, data, n):
        self.skippedCardsFile.write('skipping USET in GEOM4\n')
        return n

    def readUSET1(self, data, n):
        self.skippedCardsFile.write('skipping USET1 in GEOM4\n')
        return n
