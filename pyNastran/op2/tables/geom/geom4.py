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
            (5561,   76, 215): ['ASET', self._read_aset],      # record 1  - not done
            (5571,   77, 216): ['ASET1', self._read_aset1],    # record 2  - not done
            (10200, 102, 473): ['BNDGRID', self._read_bndgrid],# record 3  - not done
            (1510,   15, 328): ['CYAX', self._read_cyax],      # record 8  - not done
            (5210,   52, 257): ['CYJOIN', self._read_cyjoin],  # record 9  - not done
            (1710,   17, 330): ['CYSYM', self._read_cysym],    # record 11 - not done
            (4901,   49,  17): ['MPC', self._read_mpc],        # record 16 - not done
            (4891,   60,  83): ['MPCADD', self._read_mpcadd],  # record 17 - not done
            (4951,   63,  92): ['OMIT1', self._read_omit1],    # record 19 - not done

            (610,     6, 316): ['QSET1', self._read_qset1],    # record 21 - not done
            (6601,   66, 292): ['RBAR', self._read_rbar],      # record 22 - not done
            (6801,   68, 294): ['RBE1', self._read_rbe1],      # record 23 - not done
            (6901,   69, 295): ['RBE2', self._read_rbe2],      # record 24 - buggy
            (7101,   71, 187): ['RBE3', self._read_rbe3],      # record 25 - not done

            (6501,   65, 291): ['RROD', self._read_rrod],        # record 30 - not done
            (7001,   70, 186): ['RSPLINE', self._read_rspline],  # record 31 - not done
            (7201,   72, 398): ['RSSCON', self._read_rsscon],    # record 32 - not done
            (1210,   12, 322): ['SEQSET1', self._read_seqset1],  # record 40 - not done
            (5501,   55,  16): ['SPC', self._read_spc],          # record 44 - buggy
            (5481,   58,  12): ['SPC1', self._read_spc1],        # record 45 - not done
            (5491,   59,  13): ['SPCADD', self._read_spcadd],    # record 46 - not done
            (5110,   51, 256): ['SPCD', self._read_spcd],        # record 47 - buggy
            (5601,   56,  14): ['SUPORT', self._read_suport],    # record 59 - not done
            (10100, 101, 472): ['SUPORT1', self._read_suport1],  # record 60 - not done
            (2010,   20, 193) : ['USET', self._read_uset], # Record 62 -- USET(2010,20,193)

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

    def _read_aset(self, data, n):
        """ASET(5561,76,215) - Record 1"""
        self.log.debug('skipping ASET in GEOM4\n')
        return n

    def _read_aset1(self, data, n):
        """ASET1(5571,77,216) - Record 2"""
        self.log.debug('skipping ASET1 in GEOM4\n')
        return n

    def _read_bndgrid(self, data, n):
        """BNDGRID(10200,102,473) - Record 3 """
        self.log.debug('skipping BNDGRID in GEOM4\n')
        return n

# BSET
# BSET1
# CSET
# CSET1

    def _read_cyax(self, data, n):
        """CYAX(1510,15,328) - Record 8 """
        self.log.debug('skipping CYAX in GEOM4\n')
        return n

    def _read_cyjoin(self, data, n):
        """CYJOIN(5210,52,257) - Record 9 """
        self.log.debug('skipping CYJOIN in GEOM4\n')
        return n

# CYSUP
    def _read_cysym(self, data, n):
        """CYSYM(1710,17,330) - Record 11"""
        self.log.debug('skipping CYSYM in GEOM4\n')
        return n

# EGENDT
# GMBC
# GMSPC
    def _read_mpc(self, data, n):
        """MPC(4901,49,17) - Record 16"""
        self.log.debug('skipping MPC in GEOM4\n')
        return n

    def _read_mpcadd(self, data, n):
        """MPCADD(4891,60,83) - Record 17"""
        self.log.debug('skipping MPCADD in GEOM4\n')
        return n

    def _read_omit1(self, data, n):
        """OMIT1(4951,63,92) - Record 19"""
        self.log.debug('skipping OMIT1 in GEOM4\n')
        return n

    def _read_qset1(self, data, n):
        """QSET1(610, 6, 316) - Record 21"""
        self.log.debug('skipping QSET1 in GEOM4\n')
        return n

    def _read_rbar(self, data, n):
        """RBAR(6601,66,292) - Record 22"""
        self.log.debug('skipping RBAR in GEOM4\n')
        return n

    def _read_rbe1(self, data, n):
        """RBE1(6801,68,294) - Record 23"""
        self.log.debug('skipping RBE1 in GEOM4\n')
        return n

    def _read_rbe2(self, data, n):
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

    def _read_rbe3(self, data, n):
        """RBE3(7101,71,187) - Record 25"""
        self.log.debug('skipping RBE3 in GEOM4\n')
        return n

# RBJOINT
# RBJSTIF
# RELEASE
# RPNOM
    def _read_rrod(self, data, n):
        """RROD(6501,65,291) - Record 30"""
        self.log.debug('skipping RROD in GEOM4\n')
        return n

    def _read_rspline(self, data, n):
        """RSPLINE(7001,70,186) - Record 31"""
        self.log.debug('skipping RSPLINE in GEOM4\n')
        return n

    def _read_rsscon(self, data, n):
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

    def _read_seqset1(self, data, n):
        """SEQSET1(1210,12,322) - Record 40"""
        self.log.debug('skipping SEQSET1 in GEOM4\n')
        return n

# SESUP
# SEUSET
# SEUSET1

    def _read_spc(self, data, n):
        """SPC(5501,55,16) - Record 44"""
        #self.log.debug('skipping SPC in GEOM4\n')
        n = 0
        nentries = len(data) // 20  # 5*4
        for i in range(nentries):
            edata = data[n:n + 20]
            (sid, ID, c, xxx, dx) = unpack(b(self._endian + 'iiiif'), edata)

            constraint = SPC.add_op2_data([sid, ID, c, dx])
            self.add_constraint_SPC(constraint)
            n += 20
        return n

    def _read_spc1(self, data, n):
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

    def _read_spcadd(self, data, n):
        """SPCADD(5491,59,13) - Record 46"""
        self.log.debug('skipping SPCADD in GEOM4\n')
        return n

    def _read_spcd(self, data, n):
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

    def _read_spcde(self, data, n):
        self.log.debug('skipping SPCDE in GEOM4\n')
        return n

    def _read_spcf(self, data, n):
        self.log.debug('skipping SPCDF in GEOM4\n')
        return n

    def _read_spcdg(self, data, n):
        self.log.debug('skipping SPCDG in GEOM4\n')
        return n

    def _read_spce(self, data, n):
        self.log.debug('skipping SPCE in GEOM4\n')
        return n

    def _readSPCEB(self, data, n):
        self.log.debug('skipping SPCEB in GEOM4\n')
        return n

    def _read_spcf(self, data, n):
        self.log.debug('skipping SPCF in GEOM4\n')
        return n

    def _read_spcfb(self, data, n):
        self.log.debug('skipping SPCFB in GEOM4\n')
        return n

    def _read_spcgb(self, data, n):
        self.log.debug('skipping SPCGB in GEOM4\n')
        return n

# SPCGRID
# SPCOFF
# SPCOFF1

    def _read_suport(self, data, n):
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

    def _read_suport1(self, data, n):
        """SUPORT1(10100,101,472) - Record 60"""
        self.log.debug('skipping SUPORT1 in GEOM4\n')
        return n

# TEMPBC

    def _read_uset(self, data, n):
        """USET(2010,20,193) - Record 62"""
        self.log.debug('skipping USET in GEOM4\n')
        return n

    def _read_uset1(self, data, n):
        self.log.debug('skipping USET1 in GEOM4\n')
        return n
