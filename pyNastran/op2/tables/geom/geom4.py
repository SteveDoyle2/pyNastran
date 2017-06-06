"""
defines readers for BDF objects in the OP2 GEOM4/GEOM4S table
"""
#pylint: disable=C0111,C0103
from __future__ import print_function
from struct import unpack, Struct
from six import b, integer_types
from six.moves import range
import numpy as np

from pyNastran.bdf.cards.elements.rigid import RBAR, RBE2, RBE3, RROD
from pyNastran.bdf.cards.bdf_sets import (
    ASET, ASET1, BSET, BSET1, CSET, CSET1, QSET, QSET1, USET, USET1, SEQSET1 # SEQSET
)
from pyNastran.bdf.cards.loads.loads import SPCD
from pyNastran.op2.tables.geom.geom_common import GeomCommon
from pyNastran.bdf.cards.constraints import (
    SUPORT1, SUPORT,
    SPC, SPC1, SPCADD, SPCOFF, SPCOFF1,
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

            (110, 1, 311): ['BSET', self._read_bset],            # record 5  - not done
            (410, 4, 314): ['BSET1', self._read_bset1],          # record 6  - not done
            (310, 3, 313): ['CSET', self._read_cset],            # record 7  - not done
            (210, 2, 312): ['CSET1', self._read_cset1],          # record 8  - not done

            (1510, 15, 328): ['CYAX', self._read_cyax],          # record 9  - not done
            (5210, 52, 257): ['CYJOIN', self._read_cyjoin],      # record 10 - not done
            (1610, 16, 329) : ['CYSUP', self._read_cysup],       # record 11 - not done
            (1710, 17, 330): ['CYSYM', self._read_cysym],        # record 12 - not done
            (8801, 88, 9022) : ['EGENDT', self._read_egendt],    # record 13 - not done (NX)
            (9001, 90, 9024): ['FCENDT', self._read_fcendt],     # record 14 - not done (NX)
            (8001, 80, 395): ['GMBC', self._read_gmbc],          # record 15 - not done
            (7801, 78, 393): ['GMSPC', self._read_gmspc],        # record 16 - not done
            #: ['', self._read_fake],


            (4901, 49, 17) : ['MPC', self._read_mpc],             # record 17 - not done
            (4891, 60, 83) : ['MPCADD', self._read_mpcadd],       # record 18 - not done
            (5001, 50, 15) : ['OMIT', self._read_omit],           # record 19 - not done
            (4951, 63, 92) : ['OMIT1', self._read_omit1],         # record 20 - not done
            (510, 5, 315) : ['QSET', self._read_qset],            # record 21
            (610, 6, 316) : ['QSET1', self._read_qset1],          # record 22

            (6601, 66, 292) : ['RBAR', self._read_rbar],          # record 23 - not done
            (6801, 68, 294) : ['RBE1', self._read_rbe1],          # record 24 - not done
            (6901, 69, 295) : ['RBE2', self._read_rbe2],          # record 25 - buggy
            (7101, 71, 187) : ['RBE3', self._read_rbe3],          # record 26 - not done
            (14201, 142, 652) : ['RBJOINT', self._read_rbjoint],  # record 27 - not done
            (14301, 143, 653) : ['RBJSTIF', self._read_rbjstif],  # record 28 - not done
            (1310, 13, 247) : ['RELEASE', self._read_release],    # record 29 - not done
            (14101, 141, 640): ['RPNOM', self._read_rpnom],       # record 30 - not done
            (6501, 65, 291): ['RROD', self._read_rrod],           # record 31 - not done
            (7001, 70, 186): ['RSPLINE', self._read_rspline],     # record 32 - not done
            (7201, 72, 398): ['RSSCON', self._read_rsscon],       # record 33 - not done
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

            (6210, 62, 344): ['SPCOFF1', self._read_spcoff1],    # record
            (2110, 21, 194) : ['USET1', self._read_uset1],  # record
            (1010, 10, 320): ['SECSET1', self._read_secset1],  # record

            (4901, 49, 420017): ['', self._read_fake],    # record
            (5561, 76, 0): ['', self._read_fake],         # record
            (610, 6, 0): ['', self._read_fake],           # record
            (5110, 51, 620256): ['', self._read_fake],    # record
            (5501, 55, 620016): ['', self._read_fake],    # record
            (410, 4, 0): ['', self._read_fake],    # record
            (6701, 67, 293): ['RTRPLT', self._read_rtrplt],    # record 34
            (9801, 98, 79): ['', self._read_fake],  # record
            (9901, 99, 80): ['', self._read_fake],  # record
            (12001, 120, 601) : ['BLTMPC', self._read_bltmpc],  # record (NX)

            # GEOM4705 - pre MSC 2001
            (110, 1, 584): ['BNDFIX', self._read_bndfix],    # record 3 (NX)
            (210, 2, 585): ['BNDFIX1', self._read_bndfix1],    # record 4 (NX)
            (310, 3, 586) : ['BNDFREE', self._read_bndfree],  # record 5 (NX)
        }

    def _read_aset(self, data, n):
        """ASET(5561,76,215) - Record 1"""
        return self._read_xset(data, n, 'ASET', ASET, self._add_aset_object)

    def _read_qset(self, data, n):
        """QSET(610, 6, 316) - Record 21"""
        return self._read_xset(data, n, 'QSET', QSET, self._add_qset_object)

    def _read_aset1(self, data, n):
        """
        ASET1(5571,77,216) - Record 22

        ASET1=(5, 0, 4, 10, -1,
               12345, 0, 1, 2, 3, -1,
               12345, 0, 8, 9)
        """
        return self._read_xset1(data, n, 'ASET1', ASET1, self._add_aset_object)

    def _read_xset(self, data, n, card_name, cls, add_method):
        """common method for ASET, QSET; not USET"""
        s = Struct(b(self._endian + '2i'))
        #self.show_data(data, types='ifs')
        ntotal = 8
        nelements = (len(data) - n) // ntotal
        for i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  %s=%s\n' % (card_name, str(out)))
            #(id, component) = out
            set_obj = cls.add_op2_data(out)
            add_method(set_obj)
            n += ntotal
            self._increase_card_count(card_name, 1)
        return n

    def _read_xset1(self, data, n, card_name, cls, add_method, debug=False):
        r"""
        common method for ASET1, QSET1; not USET1

        F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_cntlmtlboltld02.op2
        [123   0   6  10  -1]
        """
        ndata = len(data)
        #nfields = (ndata - n) // 4
        #fmt = '%ii' % nfields
        out = np.fromstring(data[n:], self.idtype)
        #print(out)
        izero = np.where(out == -1)[0]
        if len(izero) == 0:
            card = cls.add_op2_data(out)
            add_method(card)
            self._increase_card_count(card_name, 1)
        else:
            i = np.hstack([[0], izero[:-1]+1])
            j = np.hstack([izero[:-1], -1])
            #print(i, j)
            for ii, jj in zip(i, j):
                outi = out[ii:jj]
                #print(outi)
                assert -1 not in outi, outi
                if self.is_debug_file:
                    self.binary_debug.write('  %s=%s\n' % (card_name, str(out)))
                card = cls.add_op2_data(outi)
                add_method(card)
            self._increase_card_count(card_name, len(i))
        return ndata

    def _add_superset_card(self, cls, card_name, add_method, out):
        """helper method for ``_read_superxset1``"""
        #print('out =', out)
        seid = out[0]
        components = out[1]
        thru_flag = out[2]
        if thru_flag == 0:
            nids = out[3:]
            thru_check = False
        else:
            nids = list(range(out[3], out[4]+1))
            thru_check = True

        in_data = [seid, components, nids]
        card = cls.add_op2_data(in_data)
        add_method(card)
        self._increase_card_count(card_name, 1)
        if thru_check and len(out) > 5:
            card = out[5:]
            #print('out[5:] =', out[5:])
            self._add_superset_card(cls, card_name, add_method, out[5:])

    def _read_superxset1(self, data, n, card_name, cls, add_method, debug=False):
        r"""
        common method for ASET1, QSET1; not USET1

        F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_cntlmtlboltld02.op2
        [123   0   6  10  -1]
        [  1   0   1 101 112
           2   0   1 113 124]
        """
        ndata = len(data)
        #nfields = (ndata - n) // 4
        #fmt = '%ii' % nfields
        out = np.fromstring(data[n:], self.idtype)
        #print(out)
        iminus1 = np.where(out == -1)[0]
        if len(iminus1) == 0:
            self._add_superset_card(cls, card_name, add_method, out)
        else:
            i = np.hstack([[0], iminus1[:-1]+1])
            j = np.hstack([iminus1[:-1], -1])
            #print(i, j)
            for ii, jj in zip(i, j):
                outi = out[ii:jj]
                #print(outi)
                assert -1 not in outi, outi
                if self.is_debug_file:
                    self.binary_debug.write('  %s=%s\n' % (card_name, str(out)))

                self._add_superset_card(cls, card_name, add_method, out)

                #seid = data[0]
                #components = data[1]
                #thru_flag = outi[2]
                #if thru_flag == 0:
                    #nids = outi[3:]
                #else:
                    #assert len(outi) == 5, outi
                    #nids = list(range(outi[3], outi[4]+1))

                #in_data = [seid, components, nids]

                #card = cls.add_op2_data(in_data)
                #add_method(card)
            #self._increase_card_count(card_name, len(i))
        return ndata

    def _read_bndgrid(self, data, n):
        """BNDGRID(10200,102,473) - Record 3 """
        self.log.info('skipping BNDGRID in GEOM4\n')
        return len(data)

    def _read_bset(self, data, n):
        return self._read_xset(data, n, 'BSET', BSET, self._add_bset_object)

    def _read_bset1(self, data, n):
        return self._read_xset1(data, n, 'BSET1', BSET1, self._add_bset_object)

    def _read_cset(self, data, n):
        return self._read_xset(data, n, 'CSET', CSET, self._add_cset_object)

    def _read_cset1(self, data, n):
        return self._read_xset1(data, n, 'CSET1', CSET1, self._add_cset_object)

    def _read_cyax(self, data, n):
        """CYAX(1510,15,328) - Record 8 """
        self.log.info('skipping CYAX in GEOM4\n')
        return len(data)

    def _read_cyjoin(self, data, n):
        """CYJOIN(5210,52,257) - Record 9 """
        self.log.info('skipping CYJOIN in GEOM4\n')
        return len(data)

    def _read_cysup(self, data, n):
        self.log.info('skipping CYSUP in GEOM4\n')
        return len(data)

    def _read_cysym(self, data, n):
        """CYSYM(1710,17,330) - Record 11"""
        self.log.info('skipping CYSYM in GEOM4\n')
        return len(data)

    def _read_egendt(self, data, n):
        self.log.info('skipping EGENDT in GEOM4\n')
        return len(data)

    def _read_fcendt(self, data, n):
        self.log.info('skipping FCENDT in GEOM4\n')
        return len(data)

    def _read_gmbc(self, data, n):
        self.log.info('skipping GMBC in GEOM4\n')
        return len(data)

    def _read_gmspc(self, data, n):
        self.log.info('skipping GMSPC in GEOM4\n')
        return len(data)

    def _read_mpc(self, data, n):
        """MPC(4901,49,17) - Record 16"""
        self.log.info('skipping MPC in GEOM4\n')
        return len(data)
        #self.show_data(data)
        ndata = len(data)
        nfields = (ndata-n) // 4
        #print('nfields = %s' % nfields)
        datan = data[n:]
        ints = unpack(b(self._endian + '%ii' % nfields), datan)
        floats = unpack(b(self._endian + '%if' % nfields), datan)
        i = 0
        nentries = 0
        mpc_data = []
        j = 0
        #self.show_data(data, 'if')
        while i < nfields:
            if ints[i] == -1:
                # starting a new configuration
                assert ints[i+1] == -1, ints
                assert ints[i+2] == -1, ints
                mpci = MPC.add_op2_data(mpc_data)
                self._add_constraint_mpc_object(mpci) # extracts [sid, nid, c]
                #print(mpc_data)
                nentries += 1
                if self.is_debug_file:
                    self.binary_debug.write('  MPC=%s\n' % str(mpc_data))
                mpc = []
                j = 0
                i += 2
                continue
            if j == 0:
                mpc_data = [ints[i], ints[i+1], ints[i+2], floats[i+3]]
                i += 4
            i += 1
            #print(i, nfields)
            for val in mpc_data:
                if isinstance(val, integer_types):
                    assert val != -1, mpc_data
        mpc = MPC.add_op2_data(mpc_data)
        return len(data)

    def _read_mpcadd(self, data, n):
        """MPCADD(4891,60,83) - Record 17"""
        self.log.info('skipping MPCADD in GEOM4\n')
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
        idata = np.fromstring(data[n:], self.idtype)
        fdata = np.fromstring(data[n:], self.fdtype)
        rbe3s = read_rbe3s_from_idata_fdata(self, idata, fdata)
        return n

    def _read_rbjoint(self, data, n):
        self.log.info('skipping RBJOINT in GEOM4\n')
        return len(data)

    def _read_rbjstif(self, data, n):
        self.log.info('skipping RBJSTIF in GEOM4\n')
        return len(data)

    def _read_release(self, data, n):
        self.log.info('skipping RELEASE in GEOM4\n')
        return len(data)

    def _read_rpnom(self, data, n):
        self.log.info('skipping RPNOM in GEOM4\n')
        return len(data)

    def _read_rrod(self, data, n):
        """common method for reading RROD"""
        n = self._read_dual_card(data, n, self._read_rrod_nx, self._read_rrod_msc,
                                 'RROD', self._add_rigid_element_object)
        return n

    def _read_rrod_nx(self, data, n):
        """RROD(6501,65,291) - Record 30"""
        s = Struct(b(self._endian + '5i'))
        #self.show_data(data)
        ntotal = 20
        nelements = (len(data) - n) // ntotal
        elements = []
        for i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  RROD=%s\n' % str(out))
            (eid, ga, gb, cma, cmb) = out
            out = (eid, ga, gb, cma, cmb, 0.0) # alpha
            elem = RROD.add_op2_data(out)
            elements.append(elem)
            n += ntotal
        return n, elements

    def _read_rrod_msc(self, data, n):
        """RROD(6501,65,291) - Record 30"""
        s = Struct(b(self._endian + '5if'))
        #self.show_data(data)
        ntotal = 24
        nelements = (len(data) - n) // ntotal
        elements = []
        for i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  RROD=%s\n' % str(out))
            #(eid, ga, gb, cma, cmb, alpha) = out
            elem = RROD.add_op2_data(out)
            elements.append(elem)
            n += ntotal
        return n, elements

    def _read_rspline(self, data, n):
        """RSPLINE(7001,70,186) - Record 31"""
        self.log.info('skipping RSPLINE in GEOM4\n')
        return len(data)

    def _read_rsscon(self, data, n):
        """RSSCON(7201,72,398) - Record 32"""
        self.log.info('skipping RSSCON in GEOM4\n')
        return len(data)

    def _read_rweld(self, data, n):
        self.log.info('skipping RWELD in GEOM4\n')
        return len(data)

    def _read_sebset(self, data, n):
        self.log.info('skipping SEBSET in GEOM4\n')
        return len(data)

    def _read_sebset1(self, data, n):
        self.log.info('skipping SEBSET1 in GEOM4\n')
        return len(data)

    def _read_secset(self, data, n):
        self.log.info('skipping SECSET in GEOM4\n')
        return len(data)

    def _read_secset1(self, data, n):
        self.log.info('skipping SECSET1 in GEOM4\n')
        return len(data)

    def _read_seqset(self, data, n):
        """SEQSET(1110,11,321) - Record 40"""
        self.log.info('skipping SEQSET in GEOM4\n')
        return len(data)
        #return self._read_xset(data, n, 'SEQSET', SEQSET, self.add_SEQSET)

    def _read_seqset1(self, data, n):
        """
        SEQSET1(1210,12,322) - Record 41

        SEQSET1=(1, 0, 0, 700, -1,
                 2, 0, 0, 200, -1,
                 3, 0, 0, 300, -1,
                 4, 0, 0, 400, -1,
                 5, 0, 0, 500, -1,
                 6, 0, 0, 600)
        SEQSET1=[1   0   1 101 112
                 2   0   1 113 124]
        """
        #self.log.info('skipping SEQSET1 in GEOM4\n')
        #return len(data)
        return self._read_superxset1(data, n, 'SEQSET1', SEQSET1, self._add_seqset_object,
                                     debug=True)

    def _read_sesup(self, data, n):
        self.log.info('skipping SESUP in GEOM4\n')
        return len(data)

    def _read_seuset(self, data, n):
        self.log.info('skipping SEUSET in GEOM4\n')
        return len(data)

    def _read_seuset1(self, data, n):
        self.log.info('skipping SEUSET1 in GEOM4\n')
        return len(data)

    def _read_spcoff(self, data, n):
        """SPCOFF(5501,55,16) - Record 44"""
        ntotal = 16
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n + 16]
            (sid, ID, c, dx) = unpack(b(self._endian + 'iiif'), edata)
            if self.is_debug_file:
                self.log.debug('SPCOFF sid=%s id=%s c=%s dx=%s' % (sid, ID, c, dx))
            constraint = SPCOFF.add_op2_data([sid, ID, c, dx])
            self._add_constraint_spcoff_object(constraint)
            n += 16
        return n

    def _read_spc(self, data, n):
        """SPC(5501,55,16) - Record 44"""
        ntotal = 16
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n + 16]
            (sid, ID, c, dx) = unpack(b(self._endian + 'iiif'), edata)
            if self.is_debug_file:
                self.log.debug('SPC sid=%s id=%s c=%s dx=%s' % (sid, ID, c, dx))
            constraint = SPC.add_op2_data([sid, ID, c, dx])
            self._add_constraint_spc_object(constraint)
            n += 16
        return n

    def _read_spcoff1(self, data, n):
        """
        SPCOFF1(6210, 62, 344) - Record
        see SPC1

        Record - SPC1(5481,58,12)
        Word Name Type Description
        1 SID I Set identification number  <------ removing...
        2 C I Component numbers
        3 THRUFLAG I Thru range flag
        THRUFLAG=0 No
        4 ID I Grid or scalar point identification number
        Word 4 repeats until End of Record
        THRUFLAG=1 Yes
        4 ID1 I First grid or scalar point identification number
        5 ID2 I Second grid or scalar point identification number
        End THRUFLAG

        Word Name Type Description
        1 C I Component numbers
        2 THRUFLAG I Thru range flag
        THRUFLAG=0 No
        3 ID I Grid or scalar point identification number
        Word 3 repeats until End of Record
        THRUFLAG=1 Yes
        3 ID1 I First grid or scalar point identification number
        4 ID2 I Second grid or scalar point identification number
        End THRUFLAG
        """
        nentries = 0
        nints = (len(data) - n) // 4
        idata = np.fromstring(data[n:], self.idtype)
        if not idata[-1] == -1:
            idata = np.hstack([idata, -1])
        iminus1 = np.where(idata == -1)[0]
        assert len(iminus1) > 0, idata

        i = np.hstack([[0], iminus1[:-1]+1])
        j = np.hstack([iminus1[:-1], -1])
        for ii, jj in zip(i, j):
            outi = idata[ii:jj]
            self._add_spcoff1_card(outi)
        return len(data)

    def _add_spcoff1_card(self, out):
        """helper method for ``_read_spcoff1``"""
        components, thru_flag = out[:2]
        if thru_flag == 0:  # repeat 4 to end
            nids = out[2:].tolist()
            thru_check = False
        elif thru_flag == 1:
            n1 = out[2]
            n2 = out[3]
            nids = list(range(n1, n2+1))
            thru_check = True
        else:
            raise NotImplementedError('SPCOFF1; thru_flag=%s' % thru_flag)

        assert -1 not in out, out.tolist()
        if self.is_debug_file:
            self.binary_debug.write('SPCOFF1: components=%s thru_flag=%s' % (
                components, thru_flag))
            self.binary_debug.write('   nids=%s\n' % str(nids))
        if len(nids) == 0:
            #self.log.warning('skipping SPC1 because its empty...%s' % out)
            return
        in_data = [components, nids]
        constraint = SPCOFF1.add_op2_data(in_data)
        self._add_constraint_spcoff_object(constraint)
        self._increase_card_count('SPCOFF1', 1)
        if thru_check and len(out) > 5:
            card = out[5:]
            self._add_spcoff1_card(out[5:])

    def _read_spc1(self, data, n):
        r"""
        SPC1(5481,58,12) - Record 45

        odd case = (
            # sid, comp, thru_flag
            12, 12345, 0, 110039, 110040, 110041, 110042, 110043, 110044, 110045,
                          110046, 110047, 110048, 110049, -1,
            # sid, comp, thru_flag, ???
            12, 12345, 0, -1)

        F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_acmsnnns.op2

        [1, 123456, 0, 31, 35, 39, 43, 47, 48, 53, 63, 64, 69, 70, 71, 72, -1,
         3, 456, 1, 1]
        TestOP2.test_op2_solid_bending_02_geom

        [123456, 456, 1, 5, 13,
         123456, 123456, 0, 22, 23, 24, 25, -1]
        TestOP2.test_op2_solid_shell_bar_01_geom
        """
        nentries = 0
        nints = (len(data) - n) // 4
        idata = np.fromstring(data[n:], self.idtype)
        if not idata[-1] == -1:
            idata = np.hstack([idata, -1])
        iminus1 = np.where(idata == -1)[0]
        assert len(iminus1) > 0, idata

        i = np.hstack([[0], iminus1[:-1]+1])
        j = np.hstack([iminus1[:-1], -1])
        for ii, jj in zip(i, j):
            outi = idata[ii:jj]
            self._add_spc1_card(outi)
        return len(data)

    def _add_spc1_card(self, out):
        """helper method for ``_read_spc1``"""
        sid, components = out[:2]
        thru_flag = out[2]
        if thru_flag == 0:  # repeat 4 to end
            nids = out[3:].tolist()
            thru_check = False
        elif thru_flag == 1:
            n1 = out[3]
            n2 = out[4]
            nids = list(range(n1, n2+1))
            thru_check = True
        else:
            raise NotImplementedError('SPC1; thru_flag=%s' % thru_flag)

        assert -1 not in out, out.tolist()
        if self.is_debug_file:
            self.binary_debug.write('SPC1: sid=%s components=%s thru_flag=%s' % (
                sid, components, thru_flag))
            self.binary_debug.write('   nids=%s\n' % str(nids))
        if len(nids) == 0:
            #self.log.warning('skipping SPC1 because its empty...%s' % out)
            return
        in_data = [sid, components, nids]
        constraint = SPC1.add_op2_data(in_data)
        self._add_constraint_spc_object(constraint)
        self._increase_card_count('SPC1', 1)
        if thru_check and len(out) > 5:
            card = out[5:]
            self._add_spc1_card(out[5:])

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

    def _read_spcgrid(self, data, n):
        self.log.info('skipping SPCGRID in GEOM4\n')
        return len(data)

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

    def _read_tempbc(self, data, n):
        self.log.info('skipping TEMPBC in GEOM4\n')
        return len(data)

    def _read_uset(self, data, n):
        """
        USET(2010,20,193) - Record 63
        (sid, nid, comp), ...
        """
        s = Struct(b(self._endian + '3i'))
        ntotal = 12
        #self.show_data(data, types='is')
        nelements = (len(data) - n) // ntotal
        for i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  USET=%s\n' % str(out))
            #(sid, id, component) = out
            set_obj = USET.add_op2_data(out)
            self._add_uset_object(set_obj)
            n += ntotal
        self._increase_card_count('USET', len(self.usets))
        return n

    def _read_uset1(self, data, n):
        """USET1(2110,21,194) - Record 65

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
        #print('idata = %s' % idata)
        nidata = len(idata)
        while i < nidata:
            sname = data[n+i*(4) : n+(i+1)*4]
            sname_str = unpack('4s', sname)
            #print('sname_str = %r' % sname_str)
            comp, thru_flag = idata[i+1:i+3]
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
                raise NotImplementedError('USET1; thru_flag=%s' % thru_flag)

            if self.is_debug_file:
                self.binary_debug.write('USET1: sname=%s comp=%s thru_flag=%s' % (
                    sname_str, comp, thru_flag))
                self.binary_debug.write('   nids=%s\n' % str(nids))
            #print('SPC1: sid=%s comp=%s thru_flag=%s' % (
            #    sid, comp, thru_flag))
            #print('   nids=%s\n' % str(nids))
            in_data = [sname_str, comp, nids]

            constraint = USET1.add_op2_data(in_data)
            self._add_uset_object(constraint)
        self.card_count['USET1'] = nentries
        return len(data)


    def _read_omit(self, data, n):
        self.log.info('skipping OMIT in GEOM4\n')
        return len(data)

    def _read_rtrplt(self, data, n):
        self.log.info('skipping RTRPLT in GEOM4\n')
        return len(data)

    def _read_bndfix(self, data, n):
        self.log.info('skipping BNDFIX in GEOM4\n')
        return len(data)

    def _read_bndfix1(self, data, n):
        self.log.info('skipping BNDFIX1 in GEOM4\n')
        return len(data)

    def _read_bndfree(self, data, n):
        self.log.info('skipping BNDFREE in GEOM4\n')
        return len(data)

    def _read_bltmpc(self, data, n):
        self.log.info('skipping BLTMPC in GEOM4\n')
        return len(data)

def read_rbe3s_from_idata_fdata(self, idata, fdata):
    """
    1 EID   I Element identification number
    2 REFG  I Reference grid point identification number
    3 REFC  I Component numbers at the reference grid point
    4 WT1  RS Weighting factor for components of motion at G
    5 C     I Component numbers
    6 G     I Grid point identification number

    Word 6 repeats until -1 occurs
    Words 4 through 6 repeat until -2 occurs

    7 GM    I Grid point identification number for dependent degrees-of-freedom
    8 CM    I Component numbers of dependent degrees-of-freedom

    Words 7 through 8 repeat until -3 occurs

    data = [99           99 123456 1.0    123    44    45  48  49  -1    -3]
    data = [61           71 123456 1.0    123    70    75  77      -1    -3
            62           71 123456 1.0    123    58    59  72      -1    -3]
    data = [1001100 1001100 123456 1.0 123456 10011 10002          -1 -2 -3
            1002500 1002500 123456 1.0 123456 10025 10020          -1 -2 -3]
            eid     refg    refc   wt  c      g     ...
    """
    rbe3s = []
    #iminus1 = np.where(idata == -1)[0]
    #iminus2 = np.where(idata == -2)[0]
    iminus3 = np.where(idata == -3)[0]
    #assert len(iminus1) == 1, idata
    #assert len(iminus2) == 1, idata
    #assert len(iminus3) == 1, idata
    i = np.hstack([[0], iminus3[:-1]+1])
    j = np.hstack([iminus3[:-1], len(idata)])

    #print('idata =', idata)
    for ii, jj in zip(i, j):

        eid, refg, refc, dummy, c, g = idata[ii:ii+6]
        wt = fdata[ii+3]
        weights = [wt]
        comps = [c]
        gijs = [g]
        #print('eid=%s refg=%s refc=%s wt=%s c=%s g=%s' % (
            #eid, refg, refc, wt, c, g))

        idatai = idata[ii:jj]
        iminus2 = np.where(idatai == -2)[0]
        if len(iminus2):
            self.log.info('skipping RBE3 in GEOM4\n')
            p = np.hstack([[6], iminus2[:-1]+1])
            q = np.hstack([iminus2[:-1], len(idatai)-1])
            #print('p=%s q=%s' % (p, q))

            # -2 loop (gij repeats until -2)
            for pi, qi in zip(p, q):
                break
                #print('p=%s q=%s' % (p, q))
                #idataii = idatai[pi:qi]
                #print(idataii)
                #iminus1 = np.where(idataii == -1)[0]
                #r = np.hstack([[6], iminus1[:-1]+1])
                #s = np.hstack([iminus1[:-1], len(idatai)-1])

                #ri = r[0]
                #si = s[0]
                #idataiii = idataii[ri:si]
                #print('idatai ', idataiii)

                # -1 loop ((wt, c, gij) repeats until -1)
                #aaa
                #print()
            #gijs += endi
            #Gijs = [gijs]
        else:
            assert len(iminus2) == 0, '\nii:jj=%s\nall=%s' % (idatai, idata)
            iminus1 = np.where(idatai == -1)[0]
            r = np.hstack([[6], iminus1[:-1]+1])
            s = np.hstack([iminus1[:-1], len(idatai)])
            ri = r[0]
            si = s[0]
            endi = idatai[ri:ri].tolist()
            #print('endi1 =', endi)
            for ri, si in zip(r[1:], s[1:]):
                wt = fdata[ii+ri]
                g = idata[ii+ri+1]
                c = idata[ii+ri+2]
                weights = [wt]
                comps = [c]
                gijs = [g]
                endi = idatai[ri+2:si].tolist()
                gijs += idatai[ri:si].tolist()
                #print('wt=%s g=%s c=%s gijs=%s' % (
                    #wt, g, c, str(gijs)))

            if ii + s[-1] != jj:
                msg = 'gm/cm is not supported; ii+q[0]+1=%s ii=%s jj=%s' % (ii+s[-1], ii, jj)
                raise RuntimeError(msg)

            #print('idatai =', idatai)
            #print('iminus1 =', iminus1)
            #print('----------')
            gmi = []
            cmi = []
            alpha = 0.0
            in_data = [eid, refg, refc, weights, comps, gijs,
                       gmi, cmi, alpha]
            rbe3 = RBE3.add_op2_data(in_data)
            self._add_rigid_element_object(rbe3)
            rbe3s.append(rbe3)
    return rbe3s
