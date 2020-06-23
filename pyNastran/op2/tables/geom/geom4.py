"""
defines readers for BDF objects in the OP2 GEOM4/GEOM4S table
"""
#pylint: disable=C0111,C0103,C1801
from struct import unpack, Struct
import numpy as np

from pyNastran.bdf.cards.elements.rigid import RBAR, RBE2, RBE3, RROD
from pyNastran.bdf.cards.bdf_sets import (
    ASET, ASET1, BSET, BSET1, CSET, CSET1, QSET, QSET1, USET, USET1, SEQSET1,
    OMIT1, # SEQSET
)
from pyNastran.op2.errors import MixedVersionCard
from pyNastran.op2.op2_interface.op2_reader import mapfmt, reshape_bytes_block
from pyNastran.op2.tables.geom.geom_common import GeomCommon
from pyNastran.bdf.cards.loads.loads import SPCD
from pyNastran.bdf.cards.constraints import (
    SUPORT1, SUPORT,
    SPC, SPC1, SPCADD, SPCOFF, SPCOFF1,
    MPC, MPCADD, #SPCAX, SESUP, GMSPC
)

class GEOM4(GeomCommon):
    """defines methods for reading op2 constraints"""

    def _read_geom4_4(self, data, ndata):
        """reads the GEOM4/GEOM4OLD table"""
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


            #(4901, 49, 420017): ['', self._read_fake],    # record
            (4901, 49, 420017) : ['MPC', self._read_mpc2],  # this theoretically shouldn't exist
            (4901, 49, 17) : ['MPC', self._read_mpc],             # record 17
            (4891, 60, 83) : ['MPCADD', self._read_mpcadd],       # record 18
            (5001, 50, 15) : ['OMIT', self._read_omit],           # record 19 - not done
            (4951, 63, 92) : ['OMIT1', self._read_omit1],         # record 20
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

            (5561, 76, 0): ['PLOTEL/SESET/SEQSET1?', self._read_fake],         # record
            #(5561, 76, 0): ['PLOTEL/SESET/SEQSET1?', self._read_seqset1b],         # record
            (610, 6, 0): ['SESET/SEQSET1?', self._read_fake],           # record
            (5110, 51, 620256): ['SPCD?', self._read_fake],    # record
            (5501, 55, 620016): ['SPC/SPC1?', self._read_fake],    # record
            (410, 4, 0): ['', self._read_fake],    # record
            (6701, 67, 293): ['RTRPLT', self._read_rtrplt],    # record 34
            (9801, 98, 79): ['', self._read_fake],  # record
            (9901, 99, 80): ['', self._read_fake],  # record
            (12001, 120, 601) : ['BLTMPC', self._read_bltmpc],  # record (NX)

            # GEOM4705 - pre MSC 2001
            (110, 1, 584): ['BNDFIX', self._read_bndfix],    # record 3 (NX)
            (210, 2, 585): ['BNDFIX1', self._read_bndfix1],    # record 4 (NX)
            (310, 3, 586) : ['BNDFREE', self._read_bndfree],  # record 5 (NX)

            (9801, 98, 609) : ['RVDOF', self._read_fake],
            (9901, 99, 610) : ['RVDOF1', self._read_fake],
            (11901, 119, 561) : ['RWELD', self._read_fake],
            (5571, 77, 0) : ['', self._read_fake],

            # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_sdr_s111se.op2
            (210, 2, 0) : ['', self._read_fake],
            (810, 8, 318) : ['SESET?', self._read_fake],
        }

    def _read_aset(self, data: bytes, n: int) -> int:
        """ASET(5561,76,215) - Record 1"""
        return self._read_xset(data, n, 'ASET', ASET, self._add_aset_object)

    def _read_qset(self, data: bytes, n: int) -> int:
        """QSET(610, 6, 316) - Record 21"""
        return self._read_xset(data, n, 'QSET', QSET, self._add_qset_object)

    def _read_aset1(self, data: bytes, n: int) -> int:
        """
        ASET1(5571,77,216) - Record 22

        ASET1=(5, 0, 4, 10, -1,
               12345, 0, 1, 2, 3, -1,
               12345, 0, 8, 9)
        """
        return self._read_xset1(data, n, 'ASET1', ASET1, self._add_aset_object)

    def _read_xset(self, data, n, card_name, cls, add_method):
        """common method for ASET, QSET; not USET

        Word Name Type Description
        1 ID I Grid or scalar point identification number
        2  C I Component numbers
        """
        struct_2i = Struct(self._endian + b'2i')
        #self.show_data(data, types='ifs')
        ntotal = 8
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_2i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  %s=%s\n' % (card_name, str(out)))
            #(id, component) = out
            set_obj = cls.add_op2_data(out)
            add_method(set_obj)
            n += ntotal
            self.increase_card_count(card_name, 1)
        return n

    def _read_xset1(self, data, n, card_name, cls, add_method, debug=False):
        r"""
        common method for ASET1, QSET1; not USET1

        F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_cntlmtlboltld02.op2
        [123   0   6  10  -1]

        C:\Users\sdoyle\Dropbox\move_tpl\cc439e.op2
        # this should be 4 cards
        [3, 1, 1, 8,
         3, 1, 10, 16,
         3, 1, 18, 24,
         3, 1, 26, 40]

         C:\Users\sdoyle\Dropbox\move_tpl\dogdr.op2
        [0, 1, 10001, 10050,
         3, 0, 101, 102, -1,
         3, 0, 201, 202, -1,
         3, 0, 301, 302, -1,
         3, 0, 1001, 1002, -1,
         3, 0, 2001, 2002, -1,
         3, 0, 3001, 3002, -1]
        """
        def break_by_thru_type(data):
            """helper for ``read_xset1``"""
            i = 0
            packs = []
            while i < len(data):
                #print('data[i:] = ', data[i:])
                if data[i+1] == 1:
                    pack = data[i:i+4]
                    #print('pack1 = %s' % pack)
                    packs.append(pack)
                    i += 4
                    continue

                i1 = i
                i += 3
                while data[i] != -1:
                    i += 1
                #print('pack2', data[i1:i])
                pack = data[i1:i]
                packs.append(pack)

                # get rid of the trailing -1
                i += 1
            return packs

        ndata = len(data)
        out = np.frombuffer(data[n:], self.idtype8).copy()

        #print(out)
        #izero = np.where(out == -1)[0]
        nentries = 0

        packs = break_by_thru_type(out)
        for pack in packs:
            card = cls.add_op2_data(pack)
            add_method(card)
            nentries += 1

        self.increase_card_count(card_name, nentries)
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
        assert -1 not in nids, (seid, components, nids.tolist())
        card = cls.add_op2_data(in_data)
        add_method(card)
        self.increase_card_count(card_name, 1)
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
        out = np.frombuffer(data[n:], self.idtype8).copy()
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
                assert -1 not in outi, outi
                if self.is_debug_file:
                    self.binary_debug.write('  %s=%s\n' % (card_name, str(out)))

                self._add_superset_card(cls, card_name, add_method, outi)

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
            #self.increase_card_count(card_name, len(i))
        return ndata

    def _read_bndgrid(self, data: bytes, n: int) -> int:
        """BNDGRID(10200,102,473) - Record 3 """
        self.log.info('skipping BNDGRID in GEOM4')
        return len(data)

    def _read_bset(self, data: bytes, n: int) -> int:
        return self._read_xset(data, n, 'BSET', BSET, self._add_bset_object)

    def _read_bset1(self, data: bytes, n: int) -> int:
        return self._read_xset1(data, n, 'BSET1', BSET1, self._add_bset_object)

    def _read_cset(self, data: bytes, n: int) -> int:
        return self._read_xset(data, n, 'CSET', CSET, self._add_cset_object)

    def _read_cset1(self, data: bytes, n: int) -> int:
        return self._read_xset1(data, n, 'CSET1', CSET1, self._add_cset_object)

    def _read_cyax(self, data: bytes, n: int) -> int:
        """CYAX(1510,15,328) - Record 8 """
        self.log.info('skipping CYAX in GEOM4')
        return len(data)

    def _read_cyjoin(self, data: bytes, n: int) -> int:
        """CYJOIN(5210,52,257) - Record 9 """
        self.log.info('skipping CYJOIN in GEOM4')
        return len(data)

    def _read_cysup(self, data: bytes, n: int) -> int:
        self.log.info('skipping CYSUP in GEOM4')
        return len(data)

    def _read_cysym(self, data: bytes, n: int) -> int:
        """CYSYM(1710,17,330) - Record 11"""
        self.log.info('skipping CYSYM in GEOM4')
        return len(data)

    def _read_egendt(self, data: bytes, n: int) -> int:
        self.log.info('skipping EGENDT in GEOM4')
        return len(data)

    def _read_fcendt(self, data: bytes, n: int) -> int:
        self.log.info('skipping FCENDT in GEOM4')
        return len(data)

    def _read_gmbc(self, data: bytes, n: int) -> int:
        self.log.info('skipping GMBC in GEOM4')
        return len(data)

    def _read_gmspc(self, data: bytes, n: int) -> int:
        self.log.info('skipping GMSPC in GEOM4')
        return len(data)

    def _read_mpc2(self, data: bytes, n: int) -> int:
        """MPC(4901,49,420017) - Record 16"""
        self.log.info('skipping MPC? in GEOM4')
        return len(data)

    def _read_mpc(self, data: bytes, n: int) -> int:
        """MPC(4901,49,17) - Record 16"""
        ndata = len(data)
        nfields = (ndata - n) // self.size
        datan = data[n:]
        ints = unpack(mapfmt(self._endian + b'%ii' % nfields, self.size), datan)
        floats = unpack(mapfmt(self._endian + b'%if' % nfields, self.size), datan)

        i = 0
        nentries = 0
        while i < nfields:
            sid, grid, comp = ints[i:i+3]
            coeff = floats[i+3]
            mpc_data = [sid, grid, comp, coeff]
            nodes = [grid]
            components = [comp]
            coefficients = [coeff]

            intsi = ints[i+4:i+7]
            assert len(intsi) == 3, intsi
            while intsi != (-1, -1, -1):
                gridi, compi, coeffi = ints[i+4], ints[i+5], floats[i+6]
                mpc_data.extend([gridi, compi, coeffi])
                nodes.append(gridi)
                components.append(compi)
                coefficients.append(coeffi)
                i += 3
                intsi = ints[i+4:i+7]
            mpc_data.extend([-1, -1, -1])
            i += 7 # 3 + 4 from (-1,-1,-1) and (sid,grid,comp,coeff)
            if self.is_debug_file:
                self.binary_debug.write('  MPC=%s\n' % str(mpc_data))
            mpci = MPC.add_op2_data((sid, nodes, components, coefficients))
            self._add_constraint_mpc_object(mpci)

            nentries += 1
        self.increase_card_count('MPC', nentries)
        return len(data)

    def _read_mpcadd(self, data: bytes, n: int) -> int:
        """
        MPCADD(4891,60,83) - Record 17
        """
        datai = np.frombuffer(data[n:], self.idtype8).copy()
        _read_spcadd_mpcadd(self, 'MPCADD', datai)
        return len(data)

    def _read_omit1(self, data: bytes, n: int) -> int:
        """OMIT1(4951,63,92) - Record 19"""
        return self._read_xset1(data, n, 'OMIT1', OMIT1, self._add_omit_object)

    def _read_qset1(self, data: bytes, n: int) -> int:
        """QSET1(610,6,316) - Record 22"""
        return self._read_xset1(data, n, 'QSET1', QSET1, self._add_qset_object)

    def _read_rbar_nx(self, data, n):
        """
        RBAR(6601,66,292) - Record 22 - NX version

        1 EID I Element identification number
        2 GA  I Grid point A identification number
        3 GB  I Grid point B identification number
        4 CNA I Component numbers of independent degrees-of-freedom at end A
        5 CNB I Component numbers of independent degrees-of-freedom at end B
        6 CMA I Component numbers of dependent degrees-of-freedom at end A
        7 CMB I Component numbers of dependent degrees-of-freedom at end B

        50501, 50001, 52125, 123456, 0, 0, 654321,
        50502, 50002, 52126, 123456, 0, 0, 654321,
        50503, 50003, 52127, 123456, 0, 0, 654321,
        """
        s = Struct(mapfmt(self._endian + b'7i', self.size))
        ntotal = 28 * self.factor
        nelements = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        elems = []
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]  # 8*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  RBAR NX=%s\n' % str(out))
            (eid, unused_ga, unused_gb, unused_cna, unused_cnb, unused_cma, unused_cmb) = out
            assert eid > 0, out
            out = list(out)
            out.append(0.)
            elem = RBAR.add_op2_data(out)
            elems.append(elem)
            #if self.is_debug_file:
                #self.binary_debug.write('	eid	ga	gb	cna	cnb	cma	cmb	alpha\n')
                #self.binary_debug.write(str(elem))
            n += ntotal
        self.to_nx()
        return n, elems

    def _read_rbar_msc(self, data, n):
        """RBAR(6601,66,292) - Record 22 - MSC version"""
        s = Struct(mapfmt(self._endian + b'7if', self.size))
        ntotal = 32 * self.factor  # 8*4
        nelements = (len(data) - n) // ntotal
        if not (len(data) - n) % ntotal == 0:
            raise MixedVersionCard('failed reading as MSC')
        elems = []
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  RBAR MSC=%s\n' % str(out))
            #(eid, ga, gb, cna, cnb, cma, cmb, alpha) = out
            assert out[0] > 0, out
            elem = RBAR.add_op2_data(out)
            elems.append(elem)
            n += ntotal
        return n, elems

    def _read_rbe1(self, data: bytes, n: int) -> int:
        """
        RBE1(6801,68,294) - Record 23

        MSC/NX
        Word Name Type Description
        1 EID I Element identification number
        2 GN  I Grid point identification number for independent degrees-of-freedom
        3 CN  I Component numbers of independent degrees-of-freedom
        Words 2 through 3 repeat until (-2,-2) occurs

        4 GM  I Grid point identification number for dependent degrees-of-freedom
        5 CM  I Component numbers of dependent degreesof-freedom
        Words 4 through 5 repeat until (-1,-1) occurs

        6 ALPHA RS Thermal expansion coefficient
        7 UNDEF none Not used
        """
        # TODO: neither reader or writer considers alpha; no current examples
        idata = np.frombuffer(data[n:], self.idtype8).copy()
        i = 0
        nelements = 0
        nfields = len(idata)
        while i < nfields:
            eid, gn, cn = idata[i:i+3]
            gni, cni = idata[i+3:i+5]
            Gni = [gn]
            Cni = [cn]
            while (gni, cni) != (-2, -2):
                #print(eid, gn, cn, (gni, cni))
                Gni.append(gni)
                Cni.append(cni)
                i += 2
                gni, cni = idata[i+3:i+5]
                #print((gni, cni))
                #print('----')
            i += 2

            gmi, cmi = idata[i+3:i+5]
            #print("gmi,cmi=", gmi, cmi)
            Gmi = []
            Cmi = []
            while (gmi, cmi) != (-1, -1):
                Gmi.append(gmi)
                Cmi.append(cmi)
                i += 2
                gmi, cmi = idata[i+3:i+5]
            i += 5
            #print(idata[i+3:])
            #idata
            #print(idata[i:])
            self.add_rbe1(eid, Gni, Cni, Gmi, Cmi, alpha=0.)

            nelements += 1
        self.card_count['RBE1'] = nelements
        return len(data)

    def _read_rbe2(self, data: bytes, n: int) -> int:
        """
        RBE2(6901,69,295) - Record 24

        Word Name Type Description
        1 EID I Element identification number
        2  GN I Grid point identification number for independent degrees-of-freedom
        3  CM I Component numbers of dependent degrees of-freedom
        4  GM I Grid point identification number for dependent degrees-of-freedom
        Word 4 repeats until End of Record
        5 ALPHA RS Thermal expansion coefficient

        ::

          data = (1, 1, 123456, 10000, -1, 0.0,
                  2, 2, 123456, 20000, -1, 0.0,
                  3, 3, 12345,  30000, 30001, 30002, 30003, 30004, 30005, -1, 0.0,
                  4, 4, 123,    40000, 40001, 40010, 40011, 40020, 40021, -1, 0.0,
                  5, 5, 123,    50000, 50001, 50010, 50011, 50020, 50021, -1, 0.0)
        (1,  35,  2, 33, 34, 36, 37, 133, 134, 135, 136, 137, -1, 0.0,
         3,  3,   2,  1,  2,  4,  5, 101, 102, 103, 104, 105, -1, 0.0,
         5,  8,   2,  6,  7,  9, 10, 106, 107, 108, 109, 110, -1, 0.0,
         6,  13,  2, 11, 12, 14, 15, 111, 112, 113, 114, 115, -1, 0.0
         0,  9,  30,  2, 28, 29, 31, 32,  128, 129, 130, 131, 132, -1, 0.0,
         10, 25,  2, 23, 24, 26, 27, 123, 124, 125, 126, 127, -1, 0.0)

        idata = [
            10101, 10101, 123456, 1, 2, -1,
            10102, 10102, 123456, 3, 4, -1,
            10103, 10103, 123456, 5, -1,
        ]
        """
        idata = np.frombuffer(data[n:], self.idtype8).copy()
        iminus1 = np.where(idata == -1)[0]
        if idata[-1] == -1:
            is_alpha = False
            i = np.hstack([[0], iminus1[:-1]+1])
            j = np.hstack([iminus1[:-1], len(idata)])
        else:
            is_alpha = True
            i = np.hstack([[0], iminus1[:-1]+2])
            fdata = np.frombuffer(data[n:], self.fdtype8).copy()
            j = np.hstack([iminus1[:-1]+1, len(idata)-1])
        #print('is_alpha=%s' % is_alpha)
        #print('i=%s' % i)
        #print('j=%s' % j)
        #print('idata=%s' % idata.tolist())
        #print(fdata, len(fdata))
        nelements = len(j)
        if is_alpha:
            for ii, jj in zip(i, j):
                eid, gn, cm = idata[ii:ii + 3]
                gm = idata[ii+3:jj-1].tolist()
                #print('eid=%s gn=%s cm=%s gm=%s' % (eid, gn, cm, gm))
                #alpha = fdata[jj]
                alpha = fdata[jj]
                #print('eid=%s gn=%s cm=%s gm=%s alpha=%s' % (eid, gn, cm, gm, alpha))

                out = (eid, gn, cm, gm, alpha)
                if self.is_debug_file:
                    self.binary_debug.write('  RBE2=%s\n' % str(out))
                #print('  RBE2=%s\n' % str(out))
                elem = RBE2.add_op2_data(out)
                self._add_op2_rigid_element(elem)
        else:
            alpha = 0.0
            for ii, jj in zip(i, j):
                #eid, gn, cm, gm1, gm2 = idata[ii:ii + 5]
                eid, gn, cm = idata[ii:ii + 3]
                gm = idata[ii+3:jj].tolist()
                if -1 in gm:
                    gm = gm[:-1]
                assert -1 not in gm, 'eid=%s gn=%s cm=%s gm=%s' % (eid, gn, cm, gm)
                #print('eid=%s gn=%s cm=%s gm=%s' % (eid, gn, cm, gm))

                out = (eid, gn, cm, gm, alpha)
                if self.is_debug_file:
                    self.binary_debug.write('  RBE2=%s\n' % str(out))
                #print('  RBE2=%s\n' % str(out))
                elem = RBE2.add_op2_data(out)
                self._add_op2_rigid_element(elem)
        self.card_count['RBE2'] = nelements
        return len(data)

        #while n < len(data):
            ## (eid, gn, cm, gm, ..., alpha)
            #out = s_msc.unpack(data[n:n+20])
            #eid, gn, cm, gm1, gm2 = out
            #n += 20

            #Gmi = [gm1, gm2]
            #while gm2 != -1:
                #gm2, = struct_i.unpack(data[n:n+4])
                #Gmi.append(gm2)
                #n += 4
            #Gmi = [gmi for gmi in Gmi if gmi != -1]

            ### TODO: according to the MSC/NX manual, alpha should be here,
            ###       but it's not...
            #alpha = 0.

            #if self.is_debug_file:
                #self.binary_debug.write('  RBE2=%s\n' % str(out))
            #print('  RBE2=%s\n' % str(out))
            #out = (eid, gn, cm, Gmi, alpha)
            #elem = RBE2.add_op2_data(out)
            #self._add_op2_rigid_element(elem)
            #nelements += 1
        #self.card_count['RBE2'] = nelements
        #return n

    def _read_rbe3(self, data: bytes, n: int) -> int:
        """RBE3(7101,71,187) - Record 25"""
        #self.show_data(data[n+80:], 'ifs')
        idata = np.frombuffer(data[n:], self.idtype8).copy()
        fdata = np.frombuffer(data[n:], self.fdtype8).copy()
        read_rbe3s_from_idata_fdata(self, idata, fdata)
        return len(data)

    def _read_rbjoint(self, data: bytes, n: int) -> int:
        self.log.info('skipping RBJOINT in GEOM4')
        return len(data)

    def _read_rbjstif(self, data: bytes, n: int) -> int:
        self.log.info('skipping RBJSTIF in GEOM4')
        return len(data)

    def _read_release(self, data: bytes, n: int) -> int:
        """
        Record - RELEASE(1310,13,247)

        Word Name Type Description
        1 SEID     I Superelement identification number
        2 C        I Component numbers
        3 THRUFLAG I Thru range flag
        THRUFLAG=0 No
            4 ID I Grid or scalar point identification number
            Word 4 repeats until End of Record
        THRUFLAG=1 Yes
            4 ID1 I First grid or scalar point identification number
            5 ID2
        """
        #[1310, 13, 247,
         #10, 456, 0, 10, -1,
         #20, 456, 0, 11, -1]
        from pyNastran.bdf.field_writer_16 import print_card_16
        ints = np.frombuffer(data[n:], self.idtype).copy()
        nfields = len(ints)
        i = 0
        while i < nfields:
            seid = ints[i]
            comp = ints[i + 1]
            thru_flag = ints[i + 2]
            i += 3
            if thru_flag == 0:
                value = ints[i]
                values = []
                i += 1
                while value != -1:
                    values.append(value)
                    value = ints[i]
                    i += 1
                assert len(values) > 0, 'seid=%s comp=%s thru_flag=%s values=%s' % (seid, comp, thru_flag, values)
                fields = ['RELEASE', seid, comp] + values
                self.reject_lines.append(print_card_16(fields))
            else:
                raise NotImplementedError(thru_flag)
        return len(data)

    def _read_rpnom(self, data: bytes, n: int) -> int:
        self.log.info('skipping RPNOM in GEOM4')
        return len(data)

    def _read_rrod(self, data: bytes, n: int) -> int:
        """common method for reading RROD"""
        n = self._read_dual_card(data, n, self._read_rrod_nx, self._read_rrod_msc,
                                 'RROD', self._add_op2_rigid_element)
        return n

    def _read_rrod_nx(self, data, n):
        """RROD(6501,65,291) - Record 30"""
        struct_5i = Struct(mapfmt(self._endian + b'5i', self.size))
        ntotal = 20 * self.factor
        nelements = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_5i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  RROD=%s\n' % str(out))
            (eid, ga, gb, cma, cmb) = out
            assert eid > 0, out
            out = (eid, ga, gb, cma, cmb, 0.0) # alpha
            elem = RROD.add_op2_data(out)
            elements.append(elem)
            n += ntotal
        self.to_nx()
        return n, elements

    def _read_rrod_msc(self, data, n):
        """RROD(6501,65,291) - Record 30"""
        s = Struct(mapfmt(self._endian + b'5if', self.size))
        ntotal = 24 * self.factor
        nelements = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  RROD=%s\n' % str(out))
            #(eid, ga, gb, cma, cmb, alpha) = out
            assert out[0] > 0, out
            elem = RROD.add_op2_data(out)
            elements.append(elem)
            n += ntotal
        return n, elements

    def _read_rspline(self, data: bytes, n: int) -> int:
        """RSPLINE(7001,70,186) - Record 31"""
        self.log.info('skipping RSPLINE in GEOM4')
        return len(data)

    def _read_rsscon(self, data: bytes, n: int) -> int:
        """RSSCON(7201,72,398) - Record 32"""
        self.log.info('skipping RSSCON in GEOM4')
        return len(data)

    def _read_rweld(self, data: bytes, n: int) -> int:
        self.log.info('skipping RWELD in GEOM4')
        return len(data)

    def _read_sebset(self, data: bytes, n: int) -> int:
        self.log.info('skipping SEBSET in GEOM4')
        return len(data)

    def _read_sebset1(self, data: bytes, n: int) -> int:
        self.log.info('skipping SEBSET1 in GEOM4')
        return len(data)

    def _read_secset(self, data: bytes, n: int) -> int:
        self.log.info('skipping SECSET in GEOM4')
        return len(data)

    def _read_secset1(self, data: bytes, n: int) -> int:
        self.log.info('skipping SECSET1 in GEOM4')
        return len(data)

    def _read_seqset(self, data: bytes, n: int) -> int:
        """SEQSET(1110,11,321) - Record 40"""
        self.log.info('skipping SEQSET in GEOM4')
        return len(data)
        #return self._read_xset(data, n, 'SEQSET', SEQSET, self.add_SEQSET)

    def _read_seqset1(self, data: bytes, n: int) -> int:
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

        data = (
             1, 0, 0, 700, -1,
             2, 0, 0, 200, -1,
             3, 0, 0, 300, -1,
             4, 0, 0, 400, -1,
             5, 0, 0, 500, -1,
             6, 0, 0, 600, -1)

        """
        return self._read_superxset1(data, n, 'SEQSET1', SEQSET1, self._add_seqset_object,
                                     debug=True)

    def _read_sesup(self, data: bytes, n: int) -> int:
        self.log.info('skipping SESUP in GEOM4')
        return len(data)

    def _read_seuset(self, data: bytes, n: int) -> int:
        self.log.info('skipping SEUSET in GEOM4')
        return len(data)

    def _read_seuset1(self, data: bytes, n: int) -> int:
        self.log.info('skipping SEUSET1 in GEOM4')
        return len(data)

    def _read_spcoff(self, data: bytes, n: int) -> int:
        """SPCOFF(5501,55,16) - Record 44"""
        ntotal = 16
        nentries = (len(data) - n) // ntotal
        struct_3if = Struct(self._endian + b'iiif')
        for unused_i in range(nentries):
            edata = data[n:n + 16]
            (sid, ID, c, dx) = struct_3if.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('SPCOFF sid=%s id=%s c=%s dx=%s\n' % (sid, ID, c, dx))
            constraint = SPCOFF.add_op2_data([sid, ID, c, dx])
            self._add_constraint_spcoff_object(constraint)
            n += 16
        return n

    def _read_spc(self, data: bytes, n: int) -> int:
        """common method for reading SPCs"""
        n = self._read_dual_card(data, n, self._read_spc_nx, self._read_spc_msc,
                                 'SPC', self._add_constraint_spc_object)
        return n

    def _read_spc_msc(self, data, n):
        """SPC(5501,55,16) - Record 44

        1 SID   I    Set identification number
        2 ID    I    Grid or scalar point identification number
        3 C     I    Component numbers
        4 UNDEF none Not used
        5 D     RX   Enforced displacement

        """
        ntotal = 20
        nentries = (len(data) - n) // ntotal
        assert nentries > 0, nentries
        assert (len(data) - n) % ntotal == 0
        #self.show_data(data, types='if')

        constraints = []
        struc = Struct(self._endian + b'iiiif')
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            (sid, nid, comp, xxx, dx) = struc.unpack(edata)
            assert xxx == 0, xxx
            if self.is_debug_file:
                self.binary_debug.write('SPC-MSC sid=%s id=%s comp=%s dx=%s\n' % (
                    sid, nid, comp, dx))
            assert comp != 7, 'SPC-MSC sid=%s id=%s comp=%s dx=%s\n' % (sid, nid, comp, dx)
            constraint = SPC.add_op2_data([sid, nid, comp, dx])
            constraints.append(constraint)
            n += ntotal
        return n, constraints

    def _read_spc_nx(self, data, n):
        """SPC(5501,55,16) - Record 44

        1 SID I  Set identification number
        2 ID  I  Grid or scalar point identification number
        3 C   I  Component numbers
        4 D   RS Enforced displacement

        """
        msg = ''
        ntotal = 16 * self.factor
        ndatai = len(data) - n
        nentries = ndatai // ntotal
        assert nentries > 0, nentries
        assert ndatai % ntotal == 0
        #self.show_data(data, types='if')

        struc = Struct(mapfmt(self._endian + b'iiif', self.size))
        constraints = []
        def check_component(component, msg):
            scomponent = str(component)
            sorted_components = list(set(scomponent))
            ssorted_components = ''.join(sorted_components)
            if component in [11, 22, 33, 44, 66]:
                # 11 : C:\Users\sdoyle\Dropbox\move_tpl\ifct23.op2
                # 22 : C:\Users\sdoyle\Dropbox\move_tpl\ifcpi44.op2
                # 33 : C:\Users\sdoyle\Dropbox\move_tpl\ifcq11.op2
                # 44 : C:\Users\sdoyle\Dropbox\move_tpl\pshp54.op2
                # 66 : C:\Users\sdoyle\Dropbox\move_tpl\ifsh12p.op2
                pass
            else:
                msg2 = msg + 'scomponent=%r sorted_components=%r' % (scomponent, ssorted_components)
                assert len(scomponent) == len(ssorted_components), msg2
            for componenti in scomponent:
                # 8 : C:\Users\sdoyle\Dropbox\move_tpl\beamp13.op2
                # 9 : C:\Users\sdoyle\Dropbox\move_tpl\ifcq11r.op2
                assert componenti in '0123456789', msg

        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            (sid, nid, comp, dx) = struc.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('SPC-NX sid=%s nid=%s comp=%s dx=%s\n' % (
                    sid, nid, comp, dx))

            msg = 'SPC-NX sid=%s nid=%s comp=%s dx=%s\n' % (sid, nid, comp, dx)
            check_component(comp, msg)
            if nid > 100000000:
                self._is_long_ids = True
            constraint = SPC.add_op2_data([sid, nid, comp, dx])
            constraints.append(constraint)
            #msg += '  SPC-NX sid=%s nid=%s comp=%s dx=%s\n' % (sid, nid, comp, dx)
            #else:
                #msg += '  SPC-NX sid=%s nid=%s comp=%s dx=%s\n' % (sid, nid, comp, dx)
            n += ntotal
        #if msg:
            #self.log.warning('Invalid Node IDs; skipping\n' + msg)
        #self.to_nx()
        return n, constraints

    def _read_spcoff1(self, data: bytes, n: int) -> int:
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
        #nentries = 0
        #nints = (len(data) - n) // 4
        idata = np.frombuffer(data[n:], self.idtype).copy()
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
        if max(nids) > 100000000:
            self._is_long_ids = True
        in_data = [components, nids]
        constraint = SPCOFF1.add_op2_data(in_data)
        self._add_constraint_spcoff_object(constraint)
        self.increase_card_count('SPCOFF1', 1)
        if thru_check and len(out) > 5:
            #card = out[5:]
            self._add_spcoff1_card(out[5:])

    def _read_spc1(self, data: bytes, n: int) -> int:
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
        #nentries = 0
        #nints = (len(data) - n) // 4
        idata = np.frombuffer(data[n:], self.idtype8).copy()
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
        self.increase_card_count('SPC1', 1)
        if thru_check and len(out) > 5:
            #card = out[5:]
            self._add_spc1_card(out[5:])

    def _read_spcadd(self, data: bytes, n: int) -> int:
        """SPCADD(5491,59,13) - Record 46"""
        #nentries = (len(data) - n) // 4
        datai = np.frombuffer(data[n:], self.idtype8).copy()
        _read_spcadd_mpcadd(self, 'SPCADD', datai)
        return len(data)

    def _read_spcd(self, data: bytes, n: int) -> int:
        """common method for reading SPCDs"""
        n = self._read_dual_card(data, n, self._read_spcd_nx, self._read_spcd_msc,
                                 'SPCD', self._add_load_object)
        return n

    def _read_spcd_nx(self, data, n):
        """SPCD(5110,51,256) - NX specific"""
        struct_3if = Struct(self._endian + b'3if')
        ntotal = 16 # 4*4
        nentries = (len(data) - n) // ntotal
        assert nentries > 0, nentries
        assert (len(data) - n) % ntotal == 0
        constraints = []
        is_long_ids = False
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = struct_3if.unpack(edata)
            (sid, node_id, c, dx) = out
            if self.is_debug_file:
                self.binary_debug.write('  SPCD-NX=%s\n' % str(out))

            if node_id > 100000000:
                is_long_ids = True
            constraint = SPCD.add_op2_data([sid, node_id, c, dx])
            constraints.append(constraint)
            n += ntotal
        self._is_long_ids = is_long_ids
        self.to_nx()
        return n, constraints

    def _read_spcd_msc(self, data, n):
        """
        SPCD(5110,51,256) - MSC specific - Record 47

        Word Name Type Description
        1 SID   I    Superelement identification number
        2 ID    I    Grid or scalar point identification number
        3 C     I    Component numbers
        4 UNDEF none Not used
        5 D     RX   Enforced displacement

        """
        struct_4if = Struct(self._endian + b'4if')
        ntotal = 20 # 5*4
        nentries = (len(data) - n) // ntotal
        assert nentries > 0, nentries
        assert (len(data) - n) % ntotal == 0
        constraints = []
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = struct_4if.unpack(edata)
            (sid, ID, c, xxx, dx) = out
            assert xxx == 0, xxx

            if self.is_debug_file:
                self.binary_debug.write('  SPCD-MSC=%s\n' % str(out))
            constraint = SPCD.add_op2_data([sid, ID, c, dx])
            constraints.append(constraint)
            n += ntotal
        return n, constraints

    def _read_spcde(self, data: bytes, n: int) -> int:
        self.log.info('skipping SPCDE in GEOM4')
        return len(data)

    def _read_spcf(self, data: bytes, n: int) -> int:
        self.log.info('skipping SPCDF in GEOM4')
        return len(data)

    def _read_spcdg(self, data: bytes, n: int) -> int:
        self.log.info('skipping SPCDG in GEOM4')
        return len(data)

    def _read_spce(self, data: bytes, n: int) -> int:
        self.log.info('skipping SPCE in GEOM4')
        return len(data)

    def _read_spceb(self, data: bytes, n: int) -> int:
        self.log.info('skipping SPCEB in GEOM4')
        return len(data)

    def _read_spcfb(self, data: bytes, n: int) -> int:
        self.log.info('skipping SPCFB in GEOM4')
        return len(data)

    def _read_spcgb(self, data: bytes, n: int) -> int:
        self.log.info('skipping SPCGB in GEOM4')
        return len(data)

    def _read_spcgrid(self, data: bytes, n: int) -> int:
        self.log.info('skipping SPCGRID in GEOM4')
        return len(data)

    def _read_suport(self, data: bytes, n: int) -> int:
        """SUPORT(5601,56, 14) - Record 59"""
        nentries = (len(data) - n) // 8 # 2*4
        struct_2i = Struct(self._endian + b'2i')
        for unused_i in range(nentries):
            out = list(struct_2i.unpack(data[n:n + 8]))
            if self.is_debug_file:
                self.binary_debug.write('  SUPORT=%s\n' % str(out))
                #self.log.info(out)
            suport = SUPORT.add_op2_data(out)
            self._add_suport_object(suport) # extracts [sid, c]
            n += 8
        return n

    def _read_suport1(self, data: bytes, n: int) -> int:
        """SUPORT1(10100,101,472) - Record 60"""
        nfields = (len(data) - n) // 4 - 2
        out = unpack(self._endian + b'%ii' % nfields, data[n:n+nfields*4])

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

    def _read_tempbc(self, data: bytes, n: int) -> int:
        self.log.info('skipping TEMPBC in GEOM4')
        return len(data)

    def _read_uset(self, data: bytes, n: int) -> int:
        """
        USET(2010,20,193) - Record 63
        (sid, nid, comp), ...

        """
        struct_3i = Struct(mapfmt(self._endian + b'3i', self.size))
        ntotal = 12 * self.factor
        #self.show_data(data, types='is')
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_3i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  USET=%s\n' % str(out))
            #(sid, id, component) = out
            set_obj = USET.add_op2_data(out)
            self._add_uset_object(set_obj)
            n += ntotal
        self.increase_card_count('USET', len(self.usets))
        return n

    def _read_seqset1b(self, data: bytes, n: int) -> int:  # pragma: no cover
        """

        SEQSET1(1210,12,322)
        SESET, 0, 120
        SEQSET1, 11, 0, 2001, THRU, 2008

        (5561, 76, 0, 120, -1,
         2001, -1,
         2002, -1,
         2003, -1,
         2004, -1,
         2005, -1,
         2006, -1,
         2007, -1,
         2008, -1)

        """
        #C:\NASA\m4\formats\git\examples\move_tpl\fsp11j.op2
        self.show_data(data)
        #sss

    def _read_uset1(self, data: bytes, n: int) -> int:
        """USET1(2110,21,194) - Record 65

        odd case = (
            # sid, comp, thru_flag
            12, 12345, 0, 110039, 110040, 110041, 110042, 110043, 110044, 110045,
                          110046, 110047, 110048, 110049, -1,
            # sid, comp, thru_flag, ???
            12, 12345, 0, -1)

        data =(2110, 21, 194,
               2.0, 123456, 0, 44, 45, 48, 49, -1)

        """
        assert self.factor == 1, self.factor
        nentries = 0
        size = self.size
        ints = np.frombuffer(data[n:], self.idtype8).copy()
        floats = np.frombuffer(data[n:], self.fdtype8).copy()
        i = 0
        #print('ints = %s' % ints)
        nidata = len(ints)
        while i < nidata:
            sname = data[n+i*size : n+(i+1)*size]
            sname_str = 'U%i' % int(floats[i])
            #print('sname_str = %r' % sname_str)
            #print('sname_str2 = %r' % sname_str2)

            comp, thru_flag = ints[i+1:i+3]
            comp_str = str(comp)
            for comp_stri in comp_str:
                assert comp_stri in '123456', comp_str
            i += 3
            assert thru_flag in [0, 1, 2], ints[:i+5]
            if thru_flag == 0:  # repeat 4 to end
                nid = ints[i]
                nids = [nid]
                i += 1
                if i == nidata:
                    break
                while ints[i] != -1:
                    nid = ints[i]
                    nids.append(nid)
                    i += 1
                i += 1
            elif thru_flag == 1:
                n1, n2 = ints[i:i+2]
                nids = list(range(n1, n2+1))
                i += 2
            else:
                raise NotImplementedError('USET1; thru_flag=%s' % thru_flag)

            if self.is_debug_file:
                self.binary_debug.write('USET1: sname=%s comp=%s thru_flag=%s' % (
                    sname_str, comp, thru_flag))
                self.binary_debug.write('   nids=%s\n' % str(nids))
            #print('USET1: sid=%s comp=%s thru_flag=%s' % (
            #    sid, comp, thru_flag))
            #print('   nids=%s\n' % str(nids))
            in_data = [sname_str, comp_str, nids]
            #print(in_data)

            constraint = USET1.add_op2_data(in_data)
            self._add_uset_object(constraint)
        self.card_count['USET1'] = nentries
        return len(data)


    def _read_omit(self, data: bytes, n: int) -> int:
        self.log.info('skipping OMIT in GEOM4')
        return len(data)

    def _read_rtrplt(self, data: bytes, n: int) -> int:
        self.log.info('skipping RTRPLT in GEOM4')
        return len(data)

    def _read_bndfix(self, data: bytes, n: int) -> int:
        self.log.info('skipping BNDFIX in GEOM4')
        return len(data)

    def _read_bndfix1(self, data: bytes, n: int) -> int:
        self.log.info('skipping BNDFIX1 in GEOM4')
        return len(data)

    def _read_bndfree(self, data: bytes, n: int) -> int:
        self.log.info('skipping BNDFREE in GEOM4')
        return len(data)

    def _read_bltmpc(self, data: bytes, n: int) -> int:
        self.log.info('skipping BLTMPC in GEOM4')
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
    data = [107, 1, 123456, 1.0, 1234, 10600, 10601, -1, -2, -3, 0/0.0,
            207, 2, 123456, 1.0, 1234, 20600, 20601, -1, -2, -3, 0/0.0,
            307, 3, 123456, 1.0, 1234, 30600, 30601, -1, -2, -3, 0/0.0]]

    data = [
    407, 4, 123, 1.0, 123, 41201, 41210, 41212, 41221, -1,
    -0.25, 123, 41200, 41202, 41220, 41222, -1, -2, -3, 0,
    408, 4, 456, 1.0, 123, 41201, 41210, 41212, 41221, -1,
    1.0, 123, 41200, 41202, 41220, 41222, -1, -2, -3, 0)
    ]
    RBE3    407             4       123     1.0     123     41201   41210
            41212   41221   -.25    123     41200   41202   41220   41222
    RBE3    408             4       456     1.0     123     41201   41210
            41212   41221   +.25    123     41200   41202   41220   41222

    """
    # C:\Users\sdoyle\Dropbox\move_tpl\ecbep1.op2
    rbe3s = []

    #iminus1 = np.where(idata == -1)[0]
    #iminus2 = np.where(idata == -2)[0]
    #assert len(iminus1) == 1, idata
    #assert len(iminus2) == 1, idata
    #assert len(iminus3) == 1, idata

    # -3 is the flag for not quite the end of the card, but pretty close...well depending on alpha
    # if alpha exists, we correct for it, so i/j define the start/end of the cards
    iminus3 = np.where(idata == -3)[0]
    if idata[-1] == -3:
        is_alpha = False
        i = np.hstack([[0], iminus3[:-1]+1])
    else:
        is_alpha = True
        i = np.hstack([[0], iminus3[:-1]+2])
    j = np.hstack([iminus3[:-1], len(idata)])

    #print('idata = %s' % idata.tolist())
    #print('fdata = %s' % fdata.tolist())
    #print('i = %s' % i.tolist())
    #print('j = %s' % j.tolist())
    #print('is_alpha = %s' % is_alpha)
    assert len(i) == len(j)
    for ii, jj in zip(i, j):

        idatai = idata[ii:jj]
        eid, refg, refc, dummy, unused_comp, unused_grid = idatai[:6]
        unused_weight = fdata[ii+3]

        #print('eid=%s refgrid=%s refc=%s weight=%s comp=%s grid=%s' % (
            #eid, refg, refc, weight, comp, grid))

        ii, weights, comps, gijs = fill_rbe3_wt_comp_gijs(ii, jj, idata, fdata)
        ii, gmi, cmi = _get_rbe3_um(ii, jj, idata, fdata)
        #print(idata[ii:jj].tolist())
        assert len(gijs) > 0, gijs

        if is_alpha:
            alpha = fdata[ii]
            ii += 1
        else:
            alpha = 0.0
        #print('alpha = %s' % alpha)
        in_data = [eid, refg, refc, weights, comps, gijs,
                   gmi, cmi, alpha]
        if self.is_debug_file:
            self.binary_debug.write('  RBE3=%s\n' % str(in_data))
        #print('rbe3 =', in_data)
        rbe3 = RBE3.add_op2_data(in_data)
        #print(rbe3.rstrip())

        self._add_op2_rigid_element(rbe3)
        rbe3s.append(rbe3)
        #print('--------------------------------------')
    return rbe3s

def _get_rbe3_um(i, unused_j, idata, unused_fdata):
    """helper for ``read_rbe3s_from_idata_fdata``"""
    gmi = []
    cmi = []
    if idata[i] == -2:
        i += 1
    if idata[i] != -3:
        # gm/cm
        #print('idata[i+6]*', idata[i+6])
        #print('idata[i+7]*', idata[i+7])
        #print('end = ', idata[i+6:j+1].tolist())

        while idata[i] != -3:
            gm = idata[i]
            cm = idata[i+1]
            #print('  gm=%s cm=%s' % (gm, cm))
            gmi.append(gm)
            cmi.append(cm)
            i += 2

    if idata[i] == -3:
        i += 1
    return i, gmi, cmi

def get_minus_2_index(idata):
    """helper for ``get_minus_2_index``"""
    #print('idata =', idata)
    i = np.where((idata == -2) | (idata == -3))[0]
    if len(i) == 0:
        return len(idata)
    return i[0]

def fill_rbe3_wt_comp_gijs(i, j, idata, fdata):
    """helper for ``read_rbe3s_from_idata_fdata``"""
    weights = []
    comps = []
    grids = []

    i2 = i + get_minus_2_index(idata[i:j])
    #print('i=%s i2=%s' % (i, i2))
    iold = -1
    i += 3
    while i < i2:
        #print('  idata[i:i2]=%s' % idata[i:i2].tolist())
        if iold == i:
            raise RuntimeError('infinite loop in the rbe3...')
        iold = i
        weight = fdata[i]
        comp = idata[i + 1]
        grid = idata[i + 2]
        assert grid > 0, grid
        weights.append(weight)
        comps.append(comp)
        gijs = [grid]
        grids.append(gijs)
        #print('weight=%s comp=%s grid=%s' % (
            #weight, comp, grid))
        #if idata[i + 6] in [-2, -3]:
            #break

        i += 3
        while idata[i] > 0:
            grid = idata[i]
            assert grid > 0, grid
            #print('  gridi =', grid)
            #print('  i+3=%s gridi=%s extra=%s' % (i+3, grid, idata[i+3:j+1].tolist()))
            gijs.append(grid)
            i += 1
        #print('i = ', i)
        #print('idata[i:] =', idata[i:].tolist())

        # this is the end of the -1 section (wt, comp, Gijs), but we have UM fields at the end
        if idata[i] == -3:
            #print('breaking -3...')
            break
        #assert idata[i] == -1, idata[i:j].tolist()

        # -1 defines the end of the Gij block
        #if idata[i] != -1:
            #continue
        #if idata[i] != -1 and 0:
            #msg = (
                #'Expected -1 in %i slot of idata\n'
                #'RBE3: eid=%s, refg=%s, refc=%s\n'
                #'weights=%s comps=%s\n'
                #'weight=%s comp=%s gijs=%s\n'
                #'idatai=%s\n'
                #'idata[ii+6:jj] = %s' % (i, eid, refg, refc,
                                         #weights, comps,
                                         #weight, comp, gijs,
                                         #idatai.tolist(), idata[i:j].tolist())
            #)
            #raise RuntimeError(msg)

        i += 1 # -1
        #if idata[i+3] == -2:
            #i += 1 # -2
            #break
        #print('  weight=%s comp=%s gijs=%s' % (weight, comp, gijs))
        #print('-------------------------------')
    #print('weights=%s comps=%s gijs=%s' % (weights, comps, grids))
    assert len(weights) > 0, weights
    return i, weights, comps, grids

def _read_spcadd_mpcadd(model, card_name, datai):
    """
    reads a SPCADD/MPCADD card

    Word Name Type Description
    1 SID I Set identification number
    2 S   I Set identification number
    Word 2 repeats until End of Record

    Parameters
    ----------
    model : OP2Geom()
        the model to store the data in
    card_name : str
        SPCADD or MPCADD
    datai : (n, ) int ndarray
        the data array; cannot be a List[int]
        [2  1 10 -1]
        [3  1 -1]

    """
    if model.is_debug_file:
        model.binary_debug.write('  %s - %s' % (card_name, str(datai)))
    iend = np.where(datai == -1)[0]
    if len(datai) == 3:
        dataii = datai[:-1]
        if card_name == 'MPCADD':
            constraint = MPCADD.add_op2_data(dataii)
            model._add_constraint_mpcadd_object(constraint)
        else:
            constraint = SPCADD.add_op2_data(dataii)
            model._add_constraint_spcadd_object(constraint)
        model.increase_card_count(card_name, count_num=1)
        return

    i0 = 0
    count_num = len(iend)
    for iendi in iend:
        dataii = datai[i0:iendi]
        #print('dataii = ', dataii)
        i0 = iendi + 1

        assert -1 not in dataii, dataii

        #print('%r %s' % (card_name, dataii))
        if card_name == 'MPCADD':
            constraint = MPCADD.add_op2_data(dataii)
            model._add_constraint_mpcadd_object(constraint)
        else:
            constraint = SPCADD.add_op2_data(dataii)
            model._add_constraint_spcadd_object(constraint)
    model.increase_card_count(card_name, count_num=count_num)
