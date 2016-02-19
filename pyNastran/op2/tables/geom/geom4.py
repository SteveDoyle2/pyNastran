#pylint: disable=C0301,C0111,C0103,W0613
from six import b
from six.moves import range
from struct import unpack, Struct

#from pyNastran.bdf.cards.constraints import SPC,SPCADD
from pyNastran.bdf.cards.elements.rigid import RBE2
from pyNastran.bdf.cards.constraints import SUPORT, SPC, SPC1
from pyNastran.bdf.cards.loads.loads import SPCD


class GEOM4(object):
    def add_constraint_SPC(self, constraint):
        raise RuntimeError('this should be overwritten')
    def add_rigid_element(self, constraint):
        raise RuntimeError('this should be overwritten')
    def add_suport(self, constraint):
        raise RuntimeError('this should be overwritten')

    def _read_geom4_4(self, data, ndata):
        return self._read_geom_4(self._geom4_map, data, ndata)

    def __init__(self):
        self.card_count = {}
        self._geom4_map = {
            (5561,   76, 215): ['ASET', self._readASET],      # record 1  - not done
            (5571,   77, 216): ['ASET1', self._readASET1],    # record 2  - not done
            (10200, 102, 473): ['BNDGRID', self._readBNDGRID],# record 3  - not done
            (1510,   15, 328): ['CYAX', self._readCYAX],      # record 8  - not done
            (5210,   52, 257): ['CYJOIN', self._readCYJOIN],  # record 9  - not done
            (1710,   17, 330): ['CYSYM', self._readCYSYM],    # record 11 - not done
            (4901,   49,  17): ['MPC', self._readMPC],        # record 16 - not done
            (4891,   60,  83): ['MPCADD', self._readMPCADD],  # record 17 - not done
            (4951,   63,  92): ['OMIT1', self._readOMIT1],    # record 19 - not done

            (610,     6, 316): ['QSET1', self._readQSET1],    # record 21 - not done
            (6601,   66, 292): ['RBAR', self._readRBAR],      # record 22 - not done
            (6801,   68, 294): ['RBE1', self._readRBE1],      # record 23 - not done
            (6901,   69, 295): ['RBE2', self._readRBE2],      # record 24 - buggy
            (7101,   71, 187): ['RBE3', self._readRBE3],      # record 25 - not done

            (6501,   65, 291): ['RROD', self._readRROD],        # record 30 - not done
            (7001,   70, 186): ['RSPLINE', self._readRSPLINE],  # record 31 - not done
            (7201,   72, 398): ['RSSCON', self._readRSSCON],    # record 32 - not done
            (1210,   12, 322): ['SEQSET1', self._readSEQSET1],  # record 40 - not done
            (5501,   55,  16): ['SPC', self._readSPC],          # record 44 - buggy
            (5481,   58,  12): ['SPC1', self._readSPC1],        # record 45 - not done
            (5491,   59,  13): ['SPCADD', self._readSPCADD],    # record 46 - not done
            (5110,   51, 256): ['SPCD', self._readSPCD],        # record 47 - buggy
            (5601,   56,  14): ['SUPORT', self._readSUPORT],    # record 59 - not done
            (10100, 101, 472): ['SUPORT1', self._readSUPORT1],  # record 60 - not done
            (2010,   20, 193) : ['USET', self._readUSET], # Record 62 -- USET(2010,20,193)

            (1310, 13,    247): ['', self._read_fake],    # record
            (4901, 49, 420017): ['', self._read_fake],    # record
            (5561, 76,      0): ['', self._read_fake],     # record
            (5110, 51,    256): ['', self._read_fake],     # record
            (610,   6,      0): ['', self._read_fake],     # record
            (5110, 51, 620256): ['', self._read_fake],    # record
            (5501, 55, 620016): ['', self._read_fake],    # record
            (5001, 50, 15): ['', self._read_fake],    # record
            (410, 4, 0): ['', self._read_fake],    # record
            (110, 1, 584): ['', self._read_fake],    # record
            (210, 2, 585): ['', self._read_fake],    # record
            (6210, 62, 344): ['', self._read_fake],    # record
            (510, 5, 315): ['', self._read_fake],    # record
            (6701, 67, 293): ['', self._read_fake],    # record
            (8801, 88, 9022): ['', self._read_fake],    # record
            (9001, 90, 9024): ['', self._read_fake],    # record
            (9901, 99, 80): ['', self._read_fake],  # record
            (1010, 10, 320): ['', self._read_fake],  # record
            (9801, 98, 79): ['', self._read_fake],  # record
            (12001, 120, 601) : ['', self._read_fake],  # record
            (2110, 21, 194) : ['', self._read_fake],  # record
            (310, 3, 586) : ['', self._read_fake],  # record
        }

    def _readASET(self, data, n):
        """ASET(5561,76,215) - Record 1"""
        self.log.debug('skipping ASET in GEOM4\n')
        return n

    def _readASET1(self, data, n):
        """ASET1(5571,77,216) - Record 2"""
        self.log.debug('skipping ASET1 in GEOM4\n')
        return n

    def _readBNDGRID(self, data, n):
        """BNDGRID(10200,102,473) - Record 3 """
        self.log.debug('skipping BNDGRID in GEOM4\n')
        return n

# BSET
# BSET1
# CSET
# CSET1

    def _readCYAX(self, data, n):
        """CYAX(1510,15,328) - Record 8 """
        self.log.debug('skipping CYAX in GEOM4\n')
        return n

    def _readCYJOIN(self, data, n):
        """CYJOIN(5210,52,257) - Record 9 """
        self.log.debug('skipping CYJOIN in GEOM4\n')
        return n

# CYSUP
    def _readCYSYM(self, data, n):
        """CYSYM(1710,17,330) - Record 11"""
        self.log.debug('skipping CYSYM in GEOM4\n')
        return n

# EGENDT
# GMBC
# GMSPC
    def _readMPC(self, data, n):
        """MPC(4901,49,17) - Record 16"""
        self.log.debug('skipping MPC in GEOM4\n')
        return n

    def _readMPCADD(self, data, n):
        """MPCADD(4891,60,83) - Record 17"""
        self.log.debug('skipping MPCADD in GEOM4\n')
        return n

    def _readOMIT1(self, data, n):
        """OMIT1(4951,63,92) - Record 19"""
        self.log.debug('skipping OMIT1 in GEOM4\n')
        return n

    def _readQSET1(self, data, n):
        """QSET1(610, 6, 316) - Record 21"""
        self.log.debug('skipping QSET1 in GEOM4\n')
        return n

    def _readRBAR(self, data, n):
        """RBAR(6601,66,292) - Record 22"""
        self.log.debug('skipping RBAR in GEOM4\n')
        return n

    def _readRBE1(self, data, n):
        """RBE1(6801,68,294) - Record 23"""
        self.log.debug('skipping RBE1 in GEOM4\n')
        return n

    def _readRBE2(self, data, n):
        """RBE2(6901,69,295) - Record 24"""
        self.log.debug('skipping RBE2 in GEOM4\n')
        return n
        #n=0
        #nData = len(data)  # 5*4
        if 1:
            edata = data[:12]
            (eid, gn, cm, gm) = unpack(b(self._endian + '4i'), edata)

            edata = data[12:-4]
            nGm = len(edata) // 4
            iFormat = 'i' * nGm
            iFormat = bytes(iFormat)
            Gm = list(unpack(iFormat, edata))
            alpha, = unpack(b(self._endian + 'f'), data[-4:])
        elem = RBE2(None, [eid, gn, cm, Gm, alpha])
        self.add_rigid_element(elem)
        data = data[-1:]

    def _readRBE3(self, data, n):
        """RBE3(7101,71,187) - Record 25"""
        self.log.debug('skipping RBE3 in GEOM4\n')
        return n

# RBJOINT
# RBJSTIF
# RELEASE
# RPNOM
    def _readRROD(self, data, n):
        """RROD(6501,65,291) - Record 30"""
        self.log.debug('skipping RROD in GEOM4\n')
        return n

    def _readRSPLINE(self, data, n):
        """RSPLINE(7001,70,186) - Record 31"""
        self.log.debug('skipping RSPLINE in GEOM4\n')
        return n

    def _readRSSCON(self, data, n):
        """RSSCON(7201,72,398) - Record 32"""
        self.log.debug('skipping RSSCON in GEOM4\n')
        return n

# RTRPLT
# RWELD
# SEBSET
# SEBSET1
# SECSET
# SECSET1
# SEQSET

    def _readSEQSET1(self, data, n):
        """SEQSET1(1210,12,322) - Record 40"""
        self.log.debug('skipping SEQSET1 in GEOM4\n')
        return n

# SESUP
# SEUSET
# SEUSET1

    def _readSPC(self, data, n):
        """SPC(5501,55,16) - Record 44"""
        #self.log.debug('skipping SPC in GEOM4\n')
        n = 0
        nentries = len(data) // 20  # 5*4
        for i in range(nentries):
            edata = data[n:n + 20]
            (sid, ID, c, xxx, dx) = unpack(b(self._endian + 'iiiif'), edata)

            constraint = SPC(None, [sid, ID, c, dx])
            self.add_constraint_SPC(constraint)
            n += 20
        return n

    def _readSPC1(self, data, n):
        """SPC1(5481,58,12) - Record 45"""
        self.log.debug('skipping SPC1 in GEOM4\n')
        #return n
        n2 = n
        #nentries = (len(data) - n - 12) // 4  # 5*4
        nentries = 0
        while n2 < n:
            edata = data[n:n+12]
            n += 12
            out = unpack('4i', edata)
            (sid, g, thru_flag, n1) = out
            if self.is_debug_file:
                self.binary_debug.write('  SPC1=%s\n' % str(out))
            edata = data[n:n + 12]

            nids = [n1]
            if thru_flag == 0:  # repeat 4 to end
                nnodes = (len(data) - n) // 4
                nodes = unpack(b(self._endian + '%ii' % nnodes), data[n:])
                nids += list(nodes)
                n += 4 * nentries
            else:
                n2 = self.struct_i.unpack(data[n:n+4])
                n += 4
                nids.append(n2)
            if self.is_debug_file:
                self.binary_debug.write('   nids=%s\n' % str(nids[1:]))
            nentries += 1
            constraint = SPC1(None, [sid, g, nids])
            self.add_constraint_SPC(constraint)
        self.card_count['SPC1'] = nentries
        return n

    def _readSPCADD(self, data, n):
        """SPCADD(5491,59,13) - Record 46"""
        self.log.debug('skipping SPCADD in GEOM4\n')
        return n

    def _readSPCD(self, data, n):
        """SPCD(5110,51,256) - Record 47"""
        #self.log.debug('skipping SPCD in GEOM4\n')
        n = 0
        s = Struct(b(self._endian + '4if'))
        nentries = len(data) // 20  # 5*4
        for i in range(nentries):
            edata = data[n:n + 20]
            out = s.unpack(edata)
            (sid, ID, c, xxx, dx) = out
            if self.is_debug_file:
                self.binary_debug.write('  SPCD=%s\n' % str(out))

            constraint = SPCD(None, [sid, ID, c, dx])
            self.add_constraint_SPC(constraint)
            n += 20
        self.card_count['SPCD'] = nentries
        return n

    def _readSPCDE(self, data, n):
        self.log.debug('skipping SPCDE in GEOM4\n')
        return n

    def _readSPCDF(self, data, n):
        self.log.debug('skipping SPCDF in GEOM4\n')
        return n

    def _readSPCDG(self, data, n):
        self.log.debug('skipping SPCDG in GEOM4\n')
        return n

    def _readSPCE(self, data, n):
        self.log.debug('skipping SPCE in GEOM4\n')
        return n

    def _readSPCEB(self, data, n):
        self.log.debug('skipping SPCEB in GEOM4\n')
        return n

    def _readSPCF(self, data, n):
        self.log.debug('skipping SPCF in GEOM4\n')
        return n

    def _readSPCFB(self, data, n):
        self.log.debug('skipping SPCFB in GEOM4\n')
        return n

    def _readSPCGB(self, data, n):
        self.log.debug('skipping SPCGB in GEOM4\n')
        return n

# SPCGRID
# SPCOFF
# SPCOFF1

    def _readSUPORT(self, data, n):
        """SUPORT(5601,56, 14) - Record 59"""
        #self.log.debug('skipping SUPORT in GEOM4\n')
        n = 0
        nentries = len(data) // 8  # 2*4
        s = Struct(b(self._endian + '2i'))
        for i in range(nentries):
            data_in = list(s.unpack(data[n:n + 8]))
            suport = SUPORT.add_op2_data(data_in)
            self.add_suport(suport) # extracts [sid, c]
            n += 8
        return n

    def _readSUPORT1(self, data, n):
        """SUPORT1(10100,101,472) - Record 60"""
        self.log.debug('skipping SUPORT1 in GEOM4\n')
        return n

# TEMPBC

    def _readUSET(self, data, n):
        """USET(2010,20,193) - Record 62"""
        self.log.debug('skipping USET in GEOM4\n')
        return n

    def _readUSET1(self, data, n):
        self.log.debug('skipping USET1 in GEOM4\n')
        return n
