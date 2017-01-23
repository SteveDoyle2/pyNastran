#pylint: disable=C0111,C0103,W0613
from __future__ import print_function
from struct import unpack, Struct
from six import b, integer_types
from six.moves import range

from pyNastran.bdf.cards.elements.rigid import RBAR, RBE2
from pyNastran.bdf.cards.bdf_sets import (
    ASET, ASET1, QSET, QSET1, USET, USET1, SEQSET1 # SEQSET
)
from pyNastran.bdf.cards.loads.loads import SPCD
from pyNastran.op2.tables.geom.geom_common import GeomCommon
from pyNastran.bdf.cards.constraints import (
    SUPORT1, SUPORT, SPC, SPC1, SPCADD,
    MPC, #SPCAX, MPCADD, SESUP, GMSPC
)

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
            (2010, 20, 193) : ['USET', self._read_uset],         # Record 62

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
        return self._read_xset(data, n, 'ASET', ASET, self._add_aset_object)

    def _read_qset(self, data, n):
        """QSET(610, 6, 316) - Record 21"""
        return self._read_xset(data, n, 'QSET', QSET, self._add_qset_object)

    def _read_xset(self, data, n, card_name, cls, add_method):
        """common method for ASET, QSET"""
        self.log.info('skipping %s in GEOM4\n' % card_name)
        return len(data)
        #s = Struct(b(self._endian + '2i'))
        #ntotal = 8
        #nelements = (len(data) - n) // ntotal
        #for i in range(nelements):
            #edata = data[n:n + ntotal]
            #out = s.unpack(edata)
            #if self.is_debug_file:
                #self.binary_debug.write('  %s=%s\n' % (card_name, str(out)))
            ##(id, component) = out
            #elem = cls.add_op2_data(out)
            #self.add_method(elem)
            #n += ntotal
            #self._increase_card_count(card_name, 1)
        #return n

    def _read_aset1(self, data, n):
        """ASET1(5571,77,216) - Record 22"""
        #self.log.info('skipping ASET1 in GEOM4\n')
        #return len(data)
        return self._read_xset1(data, n, 'ASET1', ASET1, self._add_aset_object)

    def _read_xset1(self, data, n, card_name, cls, add_method, debug=False):
        """common method for ASET1, QSET1"""
        self.log.info('skipping %s in GEOM4\n' % card_name)
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
        self.log.info('skipping BNDGRID in GEOM4\n')
        return len(data)

# BSET
# BSET1
# CSET
# CSET1

    def _read_cyax(self, data, n):
        """CYAX(1510,15,328) - Record 8 """
        self.log.info('skipping CYAX in GEOM4\n')
        return len(data)

    def _read_cyjoin(self, data, n):
        """CYJOIN(5210,52,257) - Record 9 """
        self.log.info('skipping CYJOIN in GEOM4\n')
        return len(data)

# CYSUP
    def _read_cysym(self, data, n):
        """CYSYM(1710,17,330) - Record 11"""
        self.log.info('skipping CYSYM in GEOM4\n')
        return len(data)

# EGENDT
# GMBC
# GMSPC
    def _read_mpc(self, data, n):
        """MPC(4901,49,17) - Record 16"""
        self.log.info('skipping MPC in GEOM4\n')
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
        mpc_data = []
        j = 0
        self.show_data(data, 'if')
        while i < nfields:
            if ints[i] == -1:
                # starting a new configuration
                assert ints[i+1] == -1, ints
                assert ints[i+2] == -1, ints
                mpci = MPC.add_op2_data(mpc_data)
                self._add_constraint_mpc_object(mpci) # extracts [sid, nid, c]
                print(mpc_data)
                nentries += 1
                if self.is_debug_file:
                    self.binary_debug.write('  MPC=%s\n' % str(mpc))
                mpc = []
                j = 0
                i += 2
                continue
            if j == 0:
                mpc_data = [ints[i], ints[i+1], ints[i+2], floats[i+3]]
                i += 4
            i += 1
            print(i, nfields)
            for val in mpc_data:
                if isinstance(val, integer_types):
                    assert val != -1, mpc_data

        mpc = MPC.add_op2_data(mpc_data)


    def _read_mpcadd(self, data, n):
        """MPCADD(4891,60,83) - Record 17"""
        self.log.info('skipping MPCADD in GEOM4\n')
        #mpcadd
        return len(data)

    def _read_omit1(self, data, n):
        """OMIT1(4951,63,92) - Record 19"""
        self.log.info('skipping OMIT1 in GEOM4\n')
        return len(data)

    def _read_qset1(self, data, n):
        """QSET1(610,6,316) - Record 22"""
        #self.log.info('skipping QSET1 in GEOM4\n')
        #return len(data)
        return self._read_xset1(data, n, 'QSET1', QSET1, self._add_qset_object)

    def _read_rbar(self, data, n):
        """RBAR(6601,66,292) - Record 22"""
        n = self._read_dual_card(data, n, self._read_rbar_nx, self._read_rbar_nx,
                                 'RBAR', self._add_rigid_element_object)
        return n

    def _read_rbar_nx(self, data, n):
        """RBAR(6601,66,292) - Record 22 - NX version"""
        # nx
        s = Struct(b(self._endian + '7i'))
        ntotal = 28
        nelements = (len(data) - n) // ntotal
        elems = []
        for i in range(nelements):
            edata = data[n:n + ntotal]  # 8*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  RBAR=%s\n' % str(out))
            (eid, ga, gb, cna, cnb, cma, cmb) = out
            out = list(out)
            out.append(0.)
            elem = RBAR.add_op2_data(out)
            elems.append(elem)
            n += ntotal
        return n, elems

    def _read_rbar_msc(self, data, n):
        """RBAR(6601,66,292) - Record 22 - MSC version"""
        s = Struct(b(self._endian + '7if'))
        ntotal = 32
        nelements = (len(data) - n) // ntotal
        elems = []
        for i in range(nelements):
            edata = data[n:n + ntotal]  # 8*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  RBAR=%s\n' % str(out))
            #(eid, ga, gb, cna, cnb, cma, cmb, alpha) = out
            elem = RBAR.add_op2_data(out)
            elems.append(elem)
            n += ntotal
        return n, elems

    def _read_rbe1(self, data, n):
        """RBE1(6801,68,294) - Record 23"""
        self.log.info('skipping RBE1 in GEOM4\n')
        return len(data)

    def _read_rbe2(self, data, n):
        """RBE2(6901,69,295) - Record 24"""
        s_nx = Struct(b(self._endian + '3i f 3i'))
        s_msc = Struct(b(self._endian + '5i'))
        struct_i = Struct(b(self._endian + 'i'))
        nelements = 0

        while n < len(data):
            # (eid, gn, cm, gm, ..., alpha)
            out = s_msc.unpack(data[n:n+20])
            eid, gn, cm, gm1, gm2 = out
            n += 20

            Gmi = [gm1, gm2]
            while gm2 != -1:
                gm2, = struct_i.unpack(data[n:n+4])
                Gmi.append(gm2)
                n += 4
            Gmi = [gmi for gmi in Gmi if gmi != -1]

            ## TODO: according to the MSC/NX manual, alpha should be here,
            ##       but it's not...
            alpha = 0.

            if self.is_debug_file:
                self.binary_debug.write('  RBE2=%s\n' % str(out))

            out = (eid, gn, cm, Gmi, alpha)
            elem = RBE2.add_op2_data(out)
            self._add_rigid_element_object(elem)
            nelements += 1
        self.card_count['RBE2'] = nelements
        return n

    def _read_rbe3(self, data, n):
        """RBE3(7101,71,187) - Record 25"""
        self.log.info('skipping RBE3 in GEOM4\n')
        return len(data)

# RBJOINT
# RBJSTIF
# RELEASE
# RPNOM
    def _read_rrod(self, data, n):
        """RROD(6501,65,291) - Record 30"""
        self.log.info('skipping RROD in GEOM4\n')
        return len(data)

    def _read_rspline(self, data, n):
        """RSPLINE(7001,70,186) - Record 31"""
        self.log.info('skipping RSPLINE in GEOM4\n')
        return len(data)

    def _read_rsscon(self, data, n):
        """RSSCON(7201,72,398) - Record 32"""
        self.log.info('skipping RSSCON in GEOM4\n')
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
        self.log.info('skipping SEQSET in GEOM4\n')
        return len(data)
        #return self._read_xset(data, n, 'SEQSET', SEQSET, self.add_SEQSET)

    def _read_seqset1(self, data, n):
        """SEQSET1(1210,12,322) - Record 41"""
        #self.log.info('skipping SEQSET1 in GEOM4\n')
        #return len(data)
        return self._read_xset1(data, n, 'SEQSET1', SEQSET1, self._add_seqset_object, debug=True)

# SESUP
# SEUSET
# SEUSET1

    def _read_spc(self, data, n):
        """SPC(5501,55,16) - Record 44"""
        #self.log.info('skipping SPC in GEOM4\n')
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
            self._add_constraint_spc_object(constraint)
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
            #self._add_load_object(load)
            #n += 20
        #return n

    def _read_spc1(self, data, n):
        """
        SPC1(5481,58,12) - Record 45

        odd case = (
            # sid, comp, thru_flag
            12, 12345, 0, 110039, 110040, 110041, 110042, 110043, 110044, 110045,
                          110046, 110047, 110048, 110049, -1,
            # sid, comp, thru_flag, ???
            12, 12345, 0, -1)
        """
        nentries = 0
        nints = (len(data) - n) // 4
        idata = unpack('%s%ii' % (self._endian, nints), data[n:])
        i = 0
        nidata = len(idata)
        while i < nidata:
            sid, comp, thru_flag = idata[i:i+3]
            i += 3
            if thru_flag == 0:  # repeat 4 to end
                nid = idata[i]
                nids = [nid]
                i += 1
                if i == nidata:
                    break
                while idata[i] != -1:
                    nid = idata[i]
                    nids.append(nid)
                    i += 1
                i += 1
            elif thru_flag == 1:
                n1, n2 = idata[i:i+2]
                nids = list(range(n1, n2+1))
                i += 2
            else:
                raise NotImplementedError('SPC1; thru_flag=%s' % thru_flag)

            if self.is_debug_file:
                self.binary_debug.write('SPC1: sid=%s comp=%s thru_flag=%s' % (
                    sid, comp, thru_flag))
                self.binary_debug.write('   nids=%s\n' % str(nids))
            #print('SPC1: sid=%s comp=%s thru_flag=%s' % (
            #    sid, comp, thru_flag))
            #print('   nids=%s\n' % str(nids))
            in_data = [sid, comp, nids]

            constraint = SPC1.add_op2_data(in_data)
            self._add_constraint_spc_object(constraint)
        self.card_count['SPC1'] = nentries
        return len(data)

    def _read_spcadd(self, data, n):
        """SPCADD(5491,59,13) - Record 46"""
        nentries = (len(data) - n) // 4
        datai = unpack('%s%si' % (self._endian, nentries), data[n:])
        if self.is_debug_file:
            self.binary_debug.write('  SPCADD - %s' % str(datai))
        #spcadd_id = datai[0]
        #values = list(datai[1:-1])
        assert datai[-1] == -1, datai
        #print('spcadd_id=%s values=%s' % (spcadd_id, values))

        constraint = SPCADD.add_op2_data(datai)
        self._add_constraint_spc_object(constraint)
        self._increase_card_count('SPCADD', count_num=1)
        return n

    def _read_spcd(self, data, n):
        """common method for reading SPCDs"""
        n = self._read_dual_card(data, n, self._read_spcd_nx, self._read_spcd_msc,
                                 'SPCD', self._add_load_object)
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
            out = s.unpack(edata)
            (sid, ID, c, xxx, dx) = out
            if self.is_debug_file:
                self.binary_debug.write('  SPCD=%s\n' % str(out))
            constraint = SPCD.add_op2_data([sid, ID, c, dx])
            loads.append(constraint)
            n += ntotal
        return n, loads

    def _read_spcde(self, data, n):
        self.log.info('skipping SPCDE in GEOM4\n')
        return len(data)

    def _read_spcf(self, data, n):
        self.log.info('skipping SPCDF in GEOM4\n')
        return len(data)

    def _read_spcdg(self, data, n):
        self.log.info('skipping SPCDG in GEOM4\n')
        return len(data)

    def _read_spce(self, data, n):
        self.log.info('skipping SPCE in GEOM4\n')
        return len(data)

    def _read_spceb(self, data, n):
        self.log.info('skipping SPCEB in GEOM4\n')
        return len(data)

    def _read_spcfb(self, data, n):
        self.log.info('skipping SPCFB in GEOM4\n')
        return len(data)

    def _read_spcgb(self, data, n):
        self.log.info('skipping SPCGB in GEOM4\n')
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
            self._add_suport_object(suport) # extracts [sid, c]
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
                #self.log.info(suporti)
                self._add_suport_object(suporti) # extracts [sid, nid, c]
                nsuports += 1
                if self.is_debug_file:
                    self.binary_debug.write('  SUPORT1=%s\n' % str(suport))
                suport = []
                i += 2
                continue
            suport.append(out[i])
            i += 1
            assert -1 not in suport, suport

        if self.is_debug_file:
            self.binary_debug.write('  SUPORT1=%s\n' % str(suport))

        suporti = SUPORT1.add_op2_data(suport)
        self._add_suport_object(suporti) # extracts [sid, nid, c]
        nsuports += 1
        self.card_count['SUPOT1'] = nsuports

        assert n+nfields*4+8 == len(data), 'a=%s b=%s' % (n+nfields*4+8, len(data))
        return len(data)

# TEMPBC

    def _read_uset(self, data, n):
        """USET(2010,20,193) - Record 63"""
        return self._read_xset(data, n, 'USET', USET, self._add_uset_object)

    def _read_uset1(self, data, n):
        """USET1(2110,21,194) - Record 65"""
        return self._read_xset1(data, n, 'USET1', USET1, self._add_uset_object)
