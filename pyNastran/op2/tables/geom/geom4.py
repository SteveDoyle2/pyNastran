#pylint: disable=C0301,C0111,C0103,W0613
from __future__ import print_function
from struct import unpack, Struct
from six import b
from six.moves import range

#from pyNastran.bdf.cards.constraints import SPC,SPCADD
from pyNastran.bdf.cards.elements.rigid import RBE2
from pyNastran.bdf.cards.bdf_sets import (
    ASET, ASET1, QSET, QSET1, USET, USET1, SEQSET, SEQSET1)
from pyNastran.bdf.cards.loads.loads import SPCD
from pyNastran.op2.tables.geom.geom_common import GeomCommon
from pyNastran.bdf.cards.constraints import SUPORT1, SUPORT, SPC, SPC1, MPC
    #SPCADD, SPCAX, MPCADD, SESUP, GMSPC

class GEOM4(GeomCommon):
    """defines methods for reading op2 constraints"""

    def _read_geom4_4(self, data, ndata):
        return self._read_geom_4(self._geom4_map, data, ndata)

    def __init__(self):
        GeomCommon.__init__(self)
        self._geom4_map = {
            (5561, 76, 215): ['ASET', self._read_aset],          # record 1
            (5571, 77, 216): ['ASET1', self._read_aset1],        # record 2
            (10200, 102, 473): ['BNDGRID', self._read_bndgrid],  # record 3  - not done

            (110, 1, 311): ['BSET', self._read_fake],           # record 5  - not done
            (410, 4, 314): ['BSET1', self._read_fake],          # record 6  - not done
            (310, 3, 313): ['CSET', self._read_fake],           # record 7  - not done
            (210, 2, 312): ['CSET1', self._read_fake],          # record 8  - not done

            (1510, 15, 328): ['CYAX', self._read_cyax],          # record 9  - not done
            (5210, 52, 257): ['CYJOIN', self._read_cyjoin],      # record 10 - not done
            (1610, 16, 329) : ['CYSUP', self._read_fake],        # record 11 - not done
            (1710, 17, 330): ['CYSYM', self._read_cysym],        # record 12 - not done
            (8801, 88, 9022) : ['EGENDT', self._read_fake],      # record 13 - not done
            (9001, 90, 9024): ['FCENDT', self._read_fake],       # record 14 - not done
            (8001, 80, 395): ['GMBC', self._read_fake],          # record 15 - not done
            (7801, 78, 393): ['GMSPC', self._read_fake],         # record 16 - not done
            #: ['', self._read_fake],


            (4901, 49, 17) : ['MPC', self._read_mpc],             # record 17 - not done
            (4891, 60, 83) : ['MPCADD', self._read_mpcadd],       # record 18 - not done
            (5001, 50, 15) : ['OMIT', self._read_fake],           # record 19 - not done
            (4951, 63, 92) : ['OMIT1', self._read_omit1],         # record 20 - not done
            (510, 5, 315) : ['QSET', self._read_qset],            # record 21
            (610, 6, 316) : ['QSET1', self._read_qset1],          # record 22

            (6601, 66, 292) : ['RBAR', self._read_rbar],          # record 23 - not done
            (6801, 68, 294) : ['RBE1', self._read_rbe1],          # record 24 - not done
            (6901, 69, 295) : ['RBE2', self._read_rbe2],          # record 25 - buggy
            (7101, 71, 187) : ['RBE3', self._read_rbe3],          # record 26 - not done
            (14201, 142, 652) : ['RBJOINT', self._read_fake],     # record 27 - not done
            (14301, 143, 653) : ['RBJSTIF', self._read_fake],     # record 28 - not done
            (1310, 13, 247) : ['RELEASE', self._read_fake],       # record 29 - not done
            (14101, 141, 640): ['RPNOM', self._read_fake],          # record 30 - not done
            (6501, 65, 291): ['RROD', self._read_rrod],          # record 31 - not done
            (7001, 70, 186): ['RSPLINE', self._read_rspline],    # record 32 - not done
            (7201, 72, 398): ['RSSCON', self._read_rsscon],      # record 33 - not done
            #: ['', self._read_fake],
            #: ['', self._read_fake],
            #: ['', self._read_fake],
            (1110, 11, 321): ['SEQSET', self._read_seqset],      # record 40
            (1210, 12, 322): ['SEQSET1', self._read_seqset1],    # record 41
            (5110, 51, 256): ['SPCD', self._read_spcd],          # record 48 - buggy

            # these ones are not fully marked...
            (5501, 55, 16): ['SPC', self._read_spc],             # record 44 - buggy
            (5481, 58, 12): ['SPC1', self._read_spc1],           # record 45 - not done
            (5491, 59, 13): ['SPCADD', self._read_spcadd],       # record 46 - not done
            (5601, 56, 14): ['SUPORT', self._read_suport],       # record 59 - not done
            (10100, 101, 472): ['SUPORT1', self._read_suport1],  # record 60 - not done
            (2010, 20, 193) : ['USET', self._read_uset],         # Record 62 -- USET(2010,20,193)

            (1310, 13, 247): ['RELEASE', self._read_fake],       # record
            (6210, 62, 344): ['SPCOFF1', self._read_fake],    # record
            (510, 5, 315): ['QSET', self._read_fake],    # record
            (2110, 21, 194) : ['USET1', self._read_fake],  # record
            (1010, 10, 320): ['SECSET1', self._read_fake],  # record
            (5001, 50, 15): ['OMIT', self._read_fake],    # record 22

            (4901, 49, 420017): ['', self._read_fake],    # record
            (5561, 76, 0): ['', self._read_fake],         # record
            (610, 6, 0): ['', self._read_fake],           # record
            (5110, 51, 620256): ['', self._read_fake],    # record
            (5501, 55, 620016): ['', self._read_fake],    # record
            (410, 4, 0): ['', self._read_fake],    # record
            (6701, 67, 293): ['RTRPLT', self._read_fake],    # record 34
            (8801, 88, 9022): ['', self._read_fake],    # record
            (9001, 90, 9024): ['', self._read_fake],    # record
            (9801, 98, 79): ['', self._read_fake],  # record
            (9901, 99, 80): ['', self._read_fake],  # record
            (12001, 120, 601) : ['', self._read_fake],  # record

            # GEOM4705 - pre MSC 2001
            (110, 1, 584): ['BNDFIX', self._read_fake],    # record 3
            (210, 2, 585): ['BNDFIX1', self._read_fake],    # record 4
            (310, 3, 586) : ['BNDFREE', self._read_fake],  # record 5
        }

    def _read_aset(self, data, n):
        """ASET(5561,76,215) - Record 1"""
        return self._read_xset(data, n, 'ASET', ASET, self.add_ASET)

    def _read_qset(self, data, n):
        """QSET(610, 6, 316) - Record 21"""
        return self._read_xset(data, n, 'QSET', QSET, self.add_QSET)

    def _read_xset(self, data, n, card_name, cls, add_method):
        """common method for ASET, QSET"""
        self.log.debug('skipping %s in GEOM4\n' % card_name)
        return len(data)
        s = Struct(b(self._endian + '2i'))
        ntotal = 8
        nelements = (len(data) - n) // ntotal
        for i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  %s=%s\n' % (card_name, str(out)))
            #(id, component) = out
            elem = cls.add_op2_data(out)
            self.add_method(elem)
            n += ntotal
            self._increase_card_count(card_name, 1)
        return n

    def _read_aset1(self, data, n):
        """ASET1(5571,77,216) - Record 22"""
        #self.log.debug('skipping ASET1 in GEOM4\n')
        #return len(data)
        return self._read_xset1(data, n, 'ASET1', ASET1, self.add_ASET)

    def _read_xset1(self, data, n, card_name, cls, add_method, debug=False):
        """common method for ASET1, QSET1"""
        self.log.debug('skipping %s in GEOM4\n' % card_name)
        return len(data)
        ndata = len(data)
        nfields = (ndata - n) // 4
        fmt = '%ii' % nfields
        out = unpack(b(self._endian + fmt), data[n:])
        if self.is_debug_file:
            self.binary_debug.write('  %s=%s\n' % (card_name, str(out)))
        if debug:
            self.log.info('  %s=%s\n' % (card_name, str(out)))
        assert -1 not in out, '  %s=%s\n' % (card_name, str(out))
        card = cls.add_op2_data(out)
        add_method(card)
        self._increase_card_count(card_name, 1)
        return ndata


    def _read_bndgrid(self, data, n):
        """BNDGRID(10200,102,473) - Record 3 """
        self.log.debug('skipping BNDGRID in GEOM4\n')
        return len(data)

# BSET
# BSET1
# CSET
# CSET1

    def _read_cyax(self, data, n):
        """CYAX(1510,15,328) - Record 8 """
        self.log.debug('skipping CYAX in GEOM4\n')
        return len(data)

    def _read_cyjoin(self, data, n):
        """CYJOIN(5210,52,257) - Record 9 """
        self.log.debug('skipping CYJOIN in GEOM4\n')
        return len(data)

# CYSUP
    def _read_cysym(self, data, n):
        """CYSYM(1710,17,330) - Record 11"""
        self.log.debug('skipping CYSYM in GEOM4\n')
        return len(data)

# EGENDT
# GMBC
# GMSPC
    def _read_mpc(self, data, n):
        """MPC(4901,49,17) - Record 16"""
        self.log.debug('skipping MPC in GEOM4\n')
        return len(data)
        #self.show_data(data)
        ndata = len(data)
        nfields = (ndata-n) // 4
        print('nfields = %s' % nfields)
        datan = data[n:]
        ints = unpack(b(self._endian + '%ii' % nfields), datan)
        floats = unpack(b(self._endian + '%if' % nfields), datan)
        i = 0
        nentries = 0
        mpc = []
        j = 0
        while i < nfields:
            if ints[i] == -1:
                assert ints[i+1] == -1, ints
                assert ints[i+2] == -1, ints
                mpci = MPC.add_op2_data(mpc)
                self.add_suport(mpci) # extracts [sid, nid, c]
                print(mpc)
                nentries += 1
                if self.is_debug_file:
                    self.binary_debug.write('  MPC=%s\n' % str(mpc))
                mpc = []
                j = 0
                i += 2
                continue
            if j == 0:
                mpc = [ints[i], ints[i+1], ints[i+2], floats[i+3]]
                i = 4
            i += 1
            #print(suport)
            assert -1 not in mpc, mpc

        MPC.add_op2_data(data)
        aaa

    def _read_mpcadd(self, data, n):
        """MPCADD(4891,60,83) - Record 17"""
        self.log.debug('skipping MPCADD in GEOM4\n')
        return len(data)

    def _read_omit1(self, data, n):
        """OMIT1(4951,63,92) - Record 19"""
        self.log.debug('skipping OMIT1 in GEOM4\n')
        return len(data)

    def _read_qset1(self, data, n):
        """QSET1(610,6,316) - Record 22"""
        #self.log.debug('skipping QSET1 in GEOM4\n')
        #return len(data)
        return self._read_xset1(data, n, 'QSET1', QSET1, self.add_QSET)

    def _read_rbar(self, data, n):
        """RBAR(6601,66,292) - Record 22"""
        self.log.debug('skipping RBAR in GEOM4\n')
        return len(data)

    def _read_rbe1(self, data, n):
        """RBE1(6801,68,294) - Record 23"""
        self.log.debug('skipping RBE1 in GEOM4\n')
        return len(data)

    def _read_rbe2(self, data, n):
        """RBE2(6901,69,295) - Record 24"""
        self.log.debug('skipping RBE2 in GEOM4\n')
        return len(data)
        #n=0
        #ndata = len(data)  # 5*4
        if 1:
            edata = data[:12]
            (eid, gn, cm, gm) = unpack(b(self._endian + '4i'), edata)

            edata = data[12:-4]
            nGm = len(edata) // 4
            iformat = bytes('i' * nGm)
            Gm = list(unpack(iformat, edata))
            alpha, = unpack(b(self._endian + 'f'), data[-4:])
        elem = RBE2.add_op2_data([eid, gn, cm, Gm, alpha])
        self.add_rigid_element(elem)
        data = data[-1:]

    def _read_rbe3(self, data, n):
        """RBE3(7101,71,187) - Record 25"""
        self.log.debug('skipping RBE3 in GEOM4\n')
        return len(data)

# RBJOINT
# RBJSTIF
# RELEASE
# RPNOM
    def _read_rrod(self, data, n):
        """RROD(6501,65,291) - Record 30"""
        self.log.debug('skipping RROD in GEOM4\n')
        return len(data)

    def _read_rspline(self, data, n):
        """RSPLINE(7001,70,186) - Record 31"""
        self.log.debug('skipping RSPLINE in GEOM4\n')
        return len(data)

    def _read_rsscon(self, data, n):
        """RSSCON(7201,72,398) - Record 32"""
        self.log.debug('skipping RSSCON in GEOM4\n')
        return len(data)

# RTRPLT
# RWELD
# SEBSET
# SEBSET1
# SECSET
# SECSET1
# SEQSET

    def _read_seqset(self, data, n):
        """SEQSET(1110,11,321) - Record 40"""
        self.log.debug('skipping SEQSET in GEOM4\n')
        return len(data)
        #return self._read_xset(data, n, 'SEQSET', SEQSET, self.add_SEQSET)

    def _read_seqset1(self, data, n):
        """SEQSET1(1210,12,322) - Record 41"""
        #self.log.debug('skipping SEQSET1 in GEOM4\n')
        #return len(data)
        return self._read_xset1(data, n, 'SEQSET1', SEQSET1, self.add_SEQSET, debug=True)

# SESUP
# SEUSET
# SEUSET1

    def _read_spc(self, data, n):
        """SPC(5501,55,16) - Record 44"""
        #self.log.debug('skipping SPC in GEOM4\n')
        #n = 0
        ntotal = 16
        nentries = (len(data) - n) // ntotal
        #nentries = len(data) // 16  # 4*4
        for i in range(nentries):
            edata = data[n:n + 16]
            (sid, ID, c, dx) = unpack(b(self._endian + 'iiif'), edata)
            if self.is_debug_file:
                self.log.debug('SPC sid=%s id=%s c=%s dx=%s' % (sid, ID, c, dx))
            constraint = SPC.add_op2_data([sid, ID, c, dx])
            self.add_constraint_SPC(constraint)
            n += 16
        return n

    #def _read_spcd(self, data, n):
        #"""SPCD(5110,51,256) - Record 44 - dmap2005"""
        #n = 0
        #nentries = len(data) // 20  # 5*4
        #for i in range(nentries):
            #edata = data[n:n + 20]
            #(sid, ID, c, xxx, dx) = unpack(b(self._endian + 'iiiif'), edata)
            #load = SPCD.add_op2_data(sid, ID, c, dx)
            #self.add_load(load)
            #n += 20
        #return n

    def _read_spc1(self, data, n):
        """SPC1(5481,58,12) - Record 45"""
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
            elif thru_flag == 1:
                n2 = self.struct_i.unpack(data[n:n+4])
                n += 4
                nids.append(n2)
            else:
                raise NotImplementedError('SPC1; thru_flag=%s' % thru_flag)
            if self.is_debug_file:
                self.binary_debug.write('   nids=%s\n' % str(nids[1:]))
            self.log.debug('   nids=%s\n' % str(nids[1:]))
            nentries += 1
            constraint = SPC1.add_op2_data([sid, g, nids])
            self.add_constraint_SPC(constraint)
        self.card_count['SPC1'] = nentries
        return n

    def _read_spcadd(self, data, n):
        """SPCADD(5491,59,13) - Record 46"""
        self.log.debug('skipping SPCADD in GEOM4\n')
        return len(data)

    def _read_spcd(self, data, n):
        """common method for reading SPCDs"""
        n = self._read_dual_card(data, n, self._read_spcd_nx, self._read_spcd_msc,
                                 'SPCD', self.add_load)
        return n

    def _read_spcd_nx(self, data, n):
        """SPCD(5110,51,256) - NX specific"""
        s = Struct(b(self._endian + '3if'))
        ntotal = 16 # 4*4
        nentries = (len(data) - n) // ntotal
        assert nentries > 0, nentries
        assert (len(data) - n) % ntotal == 0
        loads = []
        for i in range(nentries):
            edata = data[n:n + ntotal]
            #self.show_data(edata)
            out = s.unpack(edata)
            (sid, ID, c, dx) = out
            #print(out)
            if self.is_debug_file:
                self.binary_debug.write('  SPCD=%s\n' % str(out))
            constraint = SPCD.add_op2_data([sid, ID, c, dx])
            loads.append(constraint)
            n += ntotal
        return n, loads

    def _read_spcd_msc(self, data, n):
        """SPCD(5110,51,256) - MSC specific - Record 47"""
        s = Struct(b(self._endian + '4ifi'))
        ntotal = 20 # 5*4
        nentries = (len(data) - n) // ntotal
        assert nentries > 0, nentries
        assert (len(data) - n) % ntotal == 0
        loads = []
        for i in range(nentries):
            edata = data[n:n + ntotal]
            #self.show_data(edata)
            out = s.unpack(edata)
            (sid, ID, c, xxx, dx) = out
            #print(out)
            if self.is_debug_file:
                self.binary_debug.write('  SPCD=%s\n' % str(out))
            constraint = SPCD.add_op2_data([sid, ID, c, dx])
            loads.append(constraint)
            n += ntotal
        return n, loads

    def _read_spcde(self, data, n):
        self.log.debug('skipping SPCDE in GEOM4\n')
        return len(data)

    def _read_spcf(self, data, n):
        self.log.debug('skipping SPCDF in GEOM4\n')
        return len(data)

    def _read_spcdg(self, data, n):
        self.log.debug('skipping SPCDG in GEOM4\n')
        return len(data)

    def _read_spce(self, data, n):
        self.log.debug('skipping SPCE in GEOM4\n')
        return len(data)

    def _readSPCEB(self, data, n):
        self.log.debug('skipping SPCEB in GEOM4\n')
        return len(data)

    def _read_spcfb(self, data, n):
        self.log.debug('skipping SPCFB in GEOM4\n')
        return len(data)

    def _read_spcgb(self, data, n):
        self.log.debug('skipping SPCGB in GEOM4\n')
        return len(data)

# SPCGRID
# SPCOFF
# SPCOFF1

    def _read_suport(self, data, n):
        """SUPORT(5601,56, 14) - Record 59"""
        nentries = (len(data) - n) // 8 # 2*4
        s = Struct(b(self._endian + '2i'))
        for i in range(nentries):
            out = list(s.unpack(data[n:n + 8]))
            if self.is_debug_file:
                self.binary_debug.write('  SUPORT=%s\n' % str(out))
                #self.log.info(out)
            suport = SUPORT.add_op2_data(out)
            self.add_suport(suport) # extracts [sid, c]
            n += 8
        return n

    def _read_suport1(self, data, n):
        """SUPORT1(10100,101,472) - Record 60"""
        nfields = (len(data) - n) // 4 - 2
        out = unpack(b(self._endian + '%ii' % nfields), data[n:n+nfields*4])

        i = 0
        nsuports = 0
        suport = []
        while i < len(out):
            if out[i] == -1:
                assert out[i+1] == -1, out
                suporti = SUPORT1.add_op2_data(suport)
                self.add_suport(suporti) # extracts [sid, nid, c]
                #print(suporti)
                nsuports += 1
                if self.is_debug_file:
                    self.binary_debug.write('  SUPORT1=%s\n' % str(suport))
                suport = []
                i += 2
                continue
            suport.append(out[i])
            i += 1
            #print(suport)
            assert -1 not in suport, suport

        if self.is_debug_file:
            self.binary_debug.write('  SUPORT1=%s\n' % str(suport))
        suporti = SUPORT1.add_op2_data(suport)
        self.add_suport(suporti) # extracts [sid, nid, c]
        nsuports += 1
        self.card_count['SUPOT1'] = nsuports

        assert n+nfields*4+8 == len(data), 'a=%s b=%s' % (n+nfields*4+8, len(data))
        return len(data)

# TEMPBC

    def _read_uset(self, data, n):
        """USET(2010,20,193) - Record 63"""
        return self._read_xset(data, n, 'USET', USET, self.add_USET)

    def _read_uset1(self, data, n):
        """USET1(2110,21,194) - Record 65"""
        return self._read_xset1(data, n, 'USET1', USET1, self.add_USET)
