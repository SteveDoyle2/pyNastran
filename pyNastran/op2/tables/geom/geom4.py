"""
defines readers for BDF objects in the OP2 GEOM4/GEOM4S table
"""
#pylint: disable=C0111,C0103,C1801
from __future__ import annotations
from struct import unpack, Struct
from typing import Union, Type, TYPE_CHECKING
import numpy as np

from pyNastran.bdf.cards.elements.rigid import RBAR, RBE2, RBE3, RROD
from pyNastran.bdf.cards.bdf_sets import (
    ASET, ASET1, BSET, BSET1, CSET, CSET1, QSET, QSET1, USET, USET1, SEQSET1,
    OMIT, OMIT1, # SEQSET
)
from pyNastran.op2.errors import MixedVersionCard
from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.tables.geom.geom_common import GeomCommon
from pyNastran.bdf.cards.loads.loads import SPCD
from pyNastran.bdf.cards.constraints import (
    SUPORT1, SUPORT,
    SPC, SPC1, SPCADD, SPCOFF, SPCOFF1,
    MPC, MPCADD, #SPCAX, SESUP, GMSPC
)
from pyNastran.bdf.cards.optimization import DCONADD
from pyNastran.bdf.field_writer_16 import print_card_16
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2_geom import OP2Geom

class GEOM4(GeomCommon):
    """defines methods for reading op2 constraints"""

    @property
    def size(self) -> int:
        return self.op2.size
    @property
    def factor(self) -> int:
        return self.op2.factor

    def read_fake(self, data: bytes, n: int) -> int:
        return self.op2._read_fake(data, n)

    def read_stop(self, data: bytes, n: int) -> int:
        return self.op2.reader_geom1.read_stop(data, n)

    def read_geom4_4(self, data: bytes, ndata: int) -> int:
        """reads the GEOM4/GEOM4OLD table"""
        return self.op2._read_geom_4(self.geom4_map, data, ndata)

    def __init__(self, op2: OP2Geom):
        GeomCommon.__init__(self)
        self.op2 = op2
        self.geom4_map = {
            (5561, 76, 215): ['ASET', self.read_aset],          # record 1
            (5571, 77, 216): ['ASET1', self.read_aset1],        # record 2
            (10200, 102, 473): ['BNDGRID', self.read_bndgrid],  # record 3  - not done

            (110, 1, 311): ['BSET', self.read_bset],            # record 5  - not done
            (210, 2, 312): ['BSET1', self.read_bset1],          # record 6  - not done
            (310, 3, 313): ['CSET', self.read_cset],            # record 7  - not done
            (410, 4, 314): ['CSET1', self.read_cset1],          # record 8  - not done

            (1510, 15, 328): ['CYAX', self.read_cyax],          # record 9  - not done
            (5210, 52, 257): ['CYJOIN', self.read_cyjoin],      # record 10 - not done
            (1610, 16, 329) : ['CYSUP', self.read_cysup],       # record 11 - not done
            (1710, 17, 330): ['CYSYM', self.read_cysym],        # record 12 - not done
            (8801, 88, 9022) : ['EGENDT', self.read_egendt],    # record 13 - not done (NX)
            (9001, 90, 9024): ['FCENDT', self.read_fcendt],     # record 14 - not done (NX)
            (8001, 80, 395): ['GMBC', self.read_gmbc],          # record 15 - not done
            (7801, 78, 393): ['GMSPC', self.read_gmspc],        # record 16 - not done
            #: ['', self.read_fake],


            #(4901, 49, 420017): ['', self.read_fake],    # record
            (4901, 49, 420017) : ['MPC', self.read_mpc2],  # this theoretically shouldn't exist
            (4901, 49, 17) : ['MPC', self.read_mpc],             # record 17
            (4891, 60, 83) : ['MPCADD', self.read_mpcadd],       # record 18
            (5001, 50, 15) : ['OMIT', self.read_omit],           # record 19 - not done
            (4951, 63, 92) : ['OMIT1', self.read_omit1],         # record 20
            (510, 5, 315) : ['QSET', self.read_qset],            # record 21
            (610, 6, 316) : ['QSET1', self.read_qset1],          # record 22

            (6601, 66, 292) : ['RBAR', self.read_rbar],          # record 23 - not done
            (6801, 68, 294) : ['RBE1', self.read_rbe1],          # record 24 - not done
            (6901, 69, 295) : ['RBE2', self.read_rbe2],          # record 25 - buggy
            (7101, 71, 187) : ['RBE3', self.read_rbe3],          # record 26 - not done
            (14201, 142, 652) : ['RBJOINT', self.read_rbjoint],  # record 27 - not done
            (14301, 143, 653) : ['RBJSTIF', self.read_rbjstif],  # record 28 - not done
            (1310, 13, 247) : ['RELEASE', self.read_release],    # record 29 - not done
            (14101, 141, 640): ['RPNOM', self.read_rpnom],       # record 30 - not done
            (6501, 65, 291): ['RROD', self.read_rrod],           # record 31 - not done
            (7001, 70, 186): ['RSPLINE', self.read_rspline],     # record 32 - not done
            (7201, 72, 398): ['RSSCON', self.read_rsscon],       # record 33 - not done
            #: ['', self.read_fake],
            #: ['', self.read_fake],
            #: ['', self.read_fake],
            (1110, 11, 321): ['SEQSET', self.read_seqset],      # record 40
            (1210, 12, 322): ['SEQSET1', self.read_seqset1],    # record 41
            (5110, 51, 256): ['SPCD', self.read_spcd],          # record 48 - buggy

            # these ones are not fully marked...
            (5501, 55, 16): ['SPC', self.read_spc],             # record 44 - buggy
            (5481, 58, 12): ['SPC1', self.read_spc1],           # record 45 - not done
            (5491, 59, 13): ['SPCADD', self.read_spcadd],       # record 46 - not done
            (5601, 56, 14): ['SUPORT', self.read_suport],       # record 59 - not done
            (10100, 101, 472): ['SUPORT1', self.read_suport1],  # record 60 - not done
            (2010, 20, 193) : ['USET', self.read_uset],         # Record 62

            (6210, 62, 344): ['SPCOFF1', self.read_spcoff1],    # record
            (2110, 21, 194) : ['USET1', self.read_uset1],  # record
            (1010, 10, 320): ['SECSET1', self.read_secset1],  # record
            (910,   9, 319): ['SECSET', self.read_secset],  # record
            (710,   7, 317): ['SEBSET', self.read_sebset],  # record
            (810,   8, 318): ['SEBSET1', self.read_sebset1],  # record

            (1810, 18, 334): ['SEUSET', self.read_seuset],  # record
            (1910, 19, 335): ['SEUSET1', self.read_seuset1],  # record

            (5561, 76, 0): ['PLOTEL/SESET/SEQSET1?', self.read_fake],         # record
            #(5561, 76, 0): ['PLOTEL/SESET/SEQSET1?', self.read_seqset1b],         # record
            (610, 6, 0): ['SESET/SEQSET1?', self.read_seseta],           # record
            (5110, 51, 620256): ['SPCD?', self.read_fake],    # record
            (5501, 55, 620016): ['SPC', self.read_spcb],    # record
            (410, 4, 0): ['', self.read_fake],    # record
            (6701, 67, 293): ['RTRPLT', self.read_rtrplt],    # record 34
            (9801, 98, 79): ['', self.read_fake],  # record
            (9901, 99, 80): ['', self.read_fake],  # record
            (12001, 120, 601) : ['BLTMPC', self.read_bltmpc],  # record (NX)

            # GEOM4705 - pre MSC 2001
            (110, 1, 584): ['BNDFIX', self.read_bndfix],    # record 3 (NX)
            (210, 2, 585): ['BNDFIX1', self.read_bndfix1],    # record 4 (NX)
            (310, 3, 586) : ['BNDFREE', self.read_bndfree],  # record 5 (NX)

            (9801, 98, 609) : ['RVDOF', self.read_fake],
            (9901, 99, 610) : ['RVDOF1', self.read_fake],
            (11901, 119, 561) : ['RWELD', self.read_fake],
            (5571, 77, 0) : ['', self.read_fake],

            # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_sdr_s111se.op2
            (210, 2, 0) : ['', self.read_fake],
            (810, 8, 318) : ['SESET?', self.read_fake],

            (8420, 84, 641) : ['CYCADD', self.read_fake],
            (8510, 85, 643) : ['CYCAXIS', self.read_fake],
            (8220, 82, 640) : ['CYCSET', self.read_fake],

            # MSC
            (16010, 160, 934) : ['???', self.read_fake],
            (15010, 150, 841) : ['???', self.read_fake],
            (12000, 120, 255) : ['???', self.read_fake],

            # ???
            (4901, 49, 310017) : ['MPC-3', self.read_mpc],
            (4901, 49, 320017) : ['???', self.read_fake],
            (12101, 121, 696) : ['JCON', self.read_fake],
            #(4901, 49, 320017) : ['???', self.read_fake],
            #(4901, 49, 320017) : ['???', self.read_fake],

            (12100, 121, 259) : ['RBAR1', self.read_fake],
            (12200, 122, 254) : ['RTRPLT1', self.read_fake],
            (14410, 144, 320700) : ['MPCY', self.read_fake],
            (14910, 149, 820) : ['ACCSSPT', self.read_fake],
        }

    def read_seseta(self, data: bytes, n: int) -> int:  # pragma: no cover
        r"""
        (610, 6, 0,
         12, 0, 1, -1)
        C:\MSC.Software\msc_nastran_runs\see10195.op2
        """
        self.op2.log.info(f'geom skipping ??? in {self.op2.table_name}; ndata={len(data)-n}')
        self.op2.show_data(data, types='ifs')
        return len(data)

    def read_mpc3(self, data: bytes, n: int) -> int:  # pragma: no cover
        """
        MPC(4901, 49, 310017)

        SOL 402
        MPC*                   1            3746               3     1.000000000+
        *                   3786               3    -1.000000000
        data = (
            1, 3746, 3, 1.0, 3786, 3, -1.0, -1, -1, -1,
            1, 3747, 3, 1.0, 3786, 3, -1.0, -1, -1, -1,
            1, 3748, 3, 1.0, 3786, 3, -1.0, -1, -1, -1,
            1, 3749, 3, 1.0, 3786, 3, -1.0, -1, -1, -1,
            1, 3750, 3, 1.0, 3786, 3, -1.0, -1, -1, -1,
            1, 3776, 3, 1.0, 3786, 3, -1.0, -1, -1, -1,
            1, 3777, 3, 1.0, 3786, 3, -1.0, -1, -1, -1,
            1, 3778, 3, 1.0, 3786, 3, -1.0, -1, -1, -1,
            1, 3779, 3, 1.0, 3786, 3, -1.0, -1, -1, -1,
            1, 3780, 3, 1.0, 3786, 3, -1.0, -1, -1, -1,
            1, 3864, 3, 1.0, 3904, 3, -1.0, -1, -1, -1,
            1, 3865, 3, 1.0, 3904, 3, -1.0, -1, -1, -1,
            1, 3866, 3, 1.0, 3904, 3, -1.0, -1, -1, -1,
            1, 3867, 3, 1.0, 3904, 3, -1.0, -1, -1, -1,
            1, 3868, 3, 1.0, 3904, 3, -1.0, -1, -1, -1,
            1, 3894, 3, 1.0, 3904, 3, -1.0, -1, -1, -1,
            1, 3895, 3, 1.0, 3904, 3, -1.0, -1, -1, -1,
            1, 3896, 3, 1.0, 3904, 3, -1.0, -1, -1, -1,
            1, 3897, 3, 1.0, 3904, 3, -1.0, -1, -1, -1,
            1, 3898, 3, 1.0, 3904, 3, -1.0, -1, -1, -1)
        """
        self.op2.show_data(data[n:], types='qds')
        dd
        return len(data)

    def read_aset(self, data: bytes, n: int) -> int:
        """ASET(5561,76,215) - Record 1"""
        return self._read_xset(data, n, 'ASET', ASET, self.op2._add_methods._add_aset_object)

    def read_qset(self, data: bytes, n: int) -> int:
        """QSET(610, 6, 316) - Record 21"""
        return self._read_xset(data, n, 'QSET', QSET, self.op2._add_methods._add_qset_object)

    def read_aset1(self, data: bytes, n: int) -> int:
        """
        ASET1(5571,77,216) - Record 22

        ASET1=(5, 0, 4, 10, -1,
               12345, 0, 1, 2, 3, -1,
               12345, 0, 8, 9)
        """
        return self._read_xset1(data, n, 'ASET1', ASET1, self.op2._add_methods._add_aset_object)

    def _read_xset(self, data: bytes, n: int, card_name: str,
                   cls: Type[Union[ASET, BSET, CSET, QSET, OMIT]], add_method) -> int:
        """common method for ASET, QSET; not USET

        Word Name Type Description
        1 ID I Grid or scalar point identification number
        2  C I Component numbers
        """
        op2: OP2Geom = self.op2
        ints = np.frombuffer(data, dtype=op2.idtype8)[3:]
        assert len(ints) % 2 == 0, len(ints)
        grids = ints[::2]
        components = ints[1::2]
        ncards = len(grids)
        for out in zip(grids, components):
            if op2.is_debug_file:
                op2.binary_debug.write('  %s=%s\n' % (card_name, str(out)))
            #(id, component) = out
            set_obj = cls.add_op2_data(out)
            add_method(set_obj)
        op2.increase_card_count(card_name, ncards)
        return len(data)

    def _read_xset1(self, data: bytes, n: int, card_name: str,
                    cls: Type[Union[ASET1, BSET1, CSET1, QSET1, OMIT1]], add_method, debug: bool=False):
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
        op2: OP2Geom = self.op2
        ndata = len(data)
        out = np.frombuffer(data[n:], op2.idtype8).copy()

        #print(out)
        #izero = np.where(out == -1)[0]
        nentries = 0
        packs = xset1_break_by_thru_type(out)
        for pack in packs:
            card = cls.add_op2_data(pack)
            add_method(card)
            nentries += 1

        op2.increase_card_count(card_name, nentries)
        return ndata

    def _add_superset_card(self, cls, card_name: str, add_method, out):
        """helper method for ``_read_superxset1``"""
        op2: OP2Geom = self.op2
        #print('out =', out)
        seid = out[0]
        components = out[1]
        thru_flag = out[2]
        if thru_flag == 0:
            nids = out[3:].tolist()
            thru_check = False
        else:
            nids = list(range(out[3], out[4]+1))
            thru_check = True

        in_data = [seid, components, nids]
        assert -1 not in nids, (seid, components, nids)
        card = cls.add_op2_data(in_data)
        add_method(card)
        op2.increase_card_count(card_name, 1)
        if thru_check and len(out) > 5:
            card = out[5:]
            #print('out[5:] =', out[5:])
            self._add_superset_card(cls, card_name, add_method, out[5:])

    def _read_superxset1(self, data: bytes, n: int, card_name: str, cls, add_method, debug: bool=False):
        r"""
        common method for ASET1, QSET1; not USET1

        F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_cntlmtlboltld02.op2
        [123   0   6  10  -1]
        [  1   0   1 101 112
           2   0   1 113 124]
        """
        op2: OP2Geom = self.op2
        ndata = len(data)
        #nfields = (ndata - n) // 4
        #fmt = '%ii' % nfields
        out = np.frombuffer(data[n:], op2.idtype8).copy()
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
                if op2.is_debug_file:
                    op2.binary_debug.write('  %s=%s\n' % (card_name, str(out)))

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
            #op2.increase_card_count(card_name, len(i))
        return ndata

    def read_bndgrid(self, data: bytes, n: int) -> int:
        """
        BNDGRID(10200,102,473) - Record 3

        Specifies a list of grid point identification numbers on design
        boundaries or surfaces for shape optimization (SOL 200).

        Word Name Type Description
        1 GPI I Shape boundary grid point identification number
        1 ID  I Grid or scalar point identification number
        2 C   I Component numbers
        ^The cake is a lie

        BNDGRID C       GP1     GP2     ...
        BNDGRID 123     1       2       3       4       5       6       7
                71      72      73      74      75      76      77
        BNDGRID 123     8       14      15      21      22      28      29
                35      36      42      43      49      50      56      57
                63      64      70
        ints = (123, 0, 1, 2, 3, 4, 5, 6, 7, 71, 72, 73, 74, 75, 76, 77, -1,
                123, 0, 8, 14, 15, 21, 22, 28, 29, 35, 36, 42, 43, 49, 50, 56, 57, 63, 64, 70, -1)
        """

        op2: OP2Geom = self.op2
        op2.show_data(data[n:])
        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        nvalues = len(ints)

        i = 0
        value = 0
        while value != 0:
            components = ints[i]
            assert ints[i+1] == 0, ints[i:]
            value = 0
            i += 2
            values = []
            while value != -1:
                value = ints[i]
                i += 1
                if value == -1:
                    break
                values.append(value)
            components_str = str(components)
            op2.add_bndgrid(components_str, values)
        return len(data)

    def read_bset(self, data: bytes, n: int) -> int:
        return self._read_xset(data, n, 'BSET', BSET, self.op2._add_methods._add_bset_object)

    def read_bset1(self, data: bytes, n: int) -> int:
        return self._read_xset1(data, n, 'BSET1', BSET1, self.op2._add_methods._add_bset_object)

    def read_cset(self, data: bytes, n: int) -> int:
        return self._read_xset(data, n, 'CSET', CSET, self.op2._add_methods._add_cset_object)

    def read_cset1(self, data: bytes, n: int) -> int:
        return self._read_xset1(data, n, 'CSET1', CSET1, self.op2._add_methods._add_cset_object)

    def read_cyax(self, data: bytes, n: int) -> int:
        """CYAX(1510,15,328) - Record 8 """
        self.op2.log.info('geom skipping CYAX in GEOM4')
        return len(data)

    def read_cyjoin(self, data: bytes, n: int) -> int:
        """CYJOIN(5210,52,257) - Record 9 """
        self.op2.log.info('geom skipping CYJOIN in GEOM4')
        return len(data)

    def read_cysup(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping CYSUP in GEOM4')
        return len(data)

    def read_cysym(self, data: bytes, n: int) -> int:
        """CYSYM(1710,17,330) - Record 11"""
        self.op2.log.info('geom skipping CYSYM in GEOM4')
        return len(data)

    def read_egendt(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping EGENDT in GEOM4')
        return len(data)

    def read_fcendt(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping FCENDT in GEOM4')
        return len(data)

    def read_gmbc(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping GMBC in GEOM4')
        return len(data)

    def read_gmspc(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping GMSPC in GEOM4')
        return len(data)

    def read_mpc2(self, data: bytes, n: int) -> int:
        """MPC(4901,49,420017) - Record 16"""
        self.op2.log.info('geom skipping MPC? in GEOM4')
        return len(data)

    def read_mpc(self, data: bytes, n: int) -> int:
        """MPC(4901,49,17) - Record 16"""
        op2: OP2Geom = self.op2
        ndata = len(data)
        nfields = (ndata - n) // self.size
        datan = data[n:]
        ints = unpack(mapfmt(op2._endian + b'%ii' % nfields, self.size), datan)
        floats = unpack(mapfmt(op2._endian + b'%if' % nfields, self.size), datan)

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
            if op2.is_debug_file:
                op2.binary_debug.write('  MPC=%s\n' % str(mpc_data))
            mpci = MPC.add_op2_data((sid, nodes, components, coefficients))
            op2._add_methods._add_constraint_mpc_object(mpci)

            nentries += 1
        op2.increase_card_count('MPC', nentries)
        return len(data)

    def read_mpcadd(self, data: bytes, n: int) -> int:
        """
        MPCADD(4891,60,83) - Record 17
        """
        op2: OP2Geom = self.op2
        datai = np.frombuffer(data[n:], op2.idtype8).copy()
        _read_spcadd_mpcadd(op2, 'MPCADD', datai)
        return len(data)

    def read_omit1(self, data: bytes, n: int) -> int:
        """OMIT1(4951,63,92) - Record 19"""
        return self._read_xset1(data, n, 'OMIT1', OMIT1, self.op2._add_methods._add_omit_object)

    def read_qset1(self, data: bytes, n: int) -> int:
        """QSET1(610,6,316) - Record 22"""
        return self._read_xset1(data, n, 'QSET1', QSET1, self.op2._add_methods._add_qset_object)

    def read_rbar(self, data: bytes, n: int) -> int:
        """RBAR(6601,66,292) - Record 22"""
        card_name = 'RBAR'
        card_obj = RBAR
        methods = {
            28 : self._read_rbar_nx_28,
            32 : self._read_rbar_msc_32,
            36 : self._read_rbar_msc_36,
        }
        #try:
        #n = self._read_rbar_msc_36(data, n)
        n = self.op2.reader_geom2._read_double_card(
            card_name, card_obj,
            self.op2.reader_geom3._add_op2_rigid_element,
            methods, data, n)
        #except DoubleCardError:
            #raise
        return n
        #n = self.op2.reader_geom2._read_dual_card(
            #data, n, self._read_rbar_nx, self._read_rbar_msc,
            #'RBAR', self.op2.reader_geom3._add_op2_rigid_element)
        return n

    def _read_rbar_nx_28(self, card_obj, data: bytes, n: int) -> tuple[int, list[RBAR]]:
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
        op2: OP2Geom = self.op2
        ntotal = 28 * self.factor
        nelements = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        elems: list[RBAR] = []
        s = Struct(mapfmt(op2._endian + b'7i', self.size))
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]  # 8*4
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  RBAR NX=%s\n' % str(out))
            (eid, unused_ga, unused_gb, unused_cna, unused_cnb, unused_cma, unused_cmb) = out
            assert eid > 0, out
            out_list = list(out)
            out_list.append(0.)
            elem = RBAR.add_op2_data(out_list)
            elems.append(elem)
            #if op2.is_debug_file:
                #op2.binary_debug.write('	eid	ga	gb	cna	cnb	cma	cmb	alpha\n')
                #op2.binary_debug.write(str(elem))
            n += ntotal
        #op2.to_nx(' because RBAR-NX was found')
        return n, elems

    def _read_rbar_msc_32(self, card_obj, data: bytes, n: int) -> tuple[int, list[RBAR]]:
        """RBAR(6601,66,292) - Record 22 - MSC version"""
        op2: OP2Geom = self.op2
        ntotal = 32 * self.factor  # 8*4
        nelements = (len(data) - n) // ntotal
        if not (len(data) - n) % ntotal == 0:
            raise MixedVersionCard('failed reading as MSC')
        elems = []
        s = Struct(mapfmt(op2._endian + b'7if', self.size))
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  RBAR MSC=%s\n' % str(out))
            #(eid, ga, gb, cna, cnb, cma, cmb, alpha) = out
            assert out[0] > 0, out
            elem = RBAR.add_op2_data(out)
            elems.append(elem)
            n += ntotal
        return n, elems


    def _read_rbar_msc_36(self, card_obj, data: bytes, n: int) -> tuple[int, list[RBAR]]:
        """RBAR(6601,66,292) - Record 22 - MSC version

        datai  = (392, 757, 758, 123456, 0,   0,   123456, 0,   0)
        dataf  = (392, 757, 758, 123456, 0.0, 0.0, 123456, 0.0, 0.0)
        """
        op2: OP2Geom = self.op2
        log = op2.log
        ntotal = 36 * self.factor  # 9*4
        nelements = (len(data) - n) // ntotal
        if not (len(data) - n) % ntotal == 0:
            raise MixedVersionCard('failed reading as MSC')
        elems = []
        s = Struct(mapfmt(op2._endian + b'7i 2f', self.size))
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  RBAR MSC=%s\n' % str(out))
            #(eid, ga, gb, cna, cnb, cma, cmb, alpha, tref) = out
            assert out[0] > 0, out
            if out[-1] != 0.0:
                log.warning(f'new RBE3 tref? field={out[-1]} (assuming float)')
            #print(out)
            elem = RBAR.add_op2_data(out)
            elems.append(elem)
            n += ntotal
        assert n == len(data)
        return n, elems

    def read_rbe1(self, data: bytes, n: int) -> int:
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
        op2: OP2Geom = self.op2
        # TODO: neither reader or writer considers alpha; no current examples
        idata = np.frombuffer(data[n:], op2.idtype8).copy()
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
            op2.add_rbe1(eid, Gni, Cni, Gmi, Cmi, alpha=0.)

            nelements += 1
        op2.card_count['RBE1'] = nelements
        return len(data)

    def read_rbe2(self, data: bytes, n: int) -> int:
        """
        RBE2(6901,69,295) - Record 24
        """
        op2: OP2Geom = self.op2
        idata = np.frombuffer(data, op2.idtype8, offset=n).copy()
        fdata = np.frombuffer(data, op2.fdtype8, offset=n).copy()
        read_rbe2s_from_idata_fdata(op2, idata, fdata)
        return len(data)

    def read_rbe3(self, data: bytes, n: int) -> int:
        """RBE3(7101,71,187) - Record 25"""
        #self.show_data(data[n+80:], 'ifs')
        op2: OP2Geom = self.op2
        idata = np.frombuffer(data, op2.idtype8, offset=n).copy()
        fdata = np.frombuffer(data, op2.fdtype8, offset=n).copy()
        read_rbe3s_from_idata_fdata(op2, idata, fdata)
        return len(data)

    def read_rbjoint(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping RBJOINT in GEOM4')
        return len(data)

    def read_rbjstif(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping RBJSTIF in GEOM4')
        return len(data)

    def read_release(self, data: bytes, n: int) -> int:
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
        op2: OP2Geom = self.op2
        #[1310, 13, 247,
         #10, 456, 0, 10, -1,
         #20, 456, 0, 11, -1]
        ints = np.frombuffer(data[n:], op2.idtype).copy()

        #i, cards = ints_to_secset1s('RELEASE', ints)
        #for (seid, comp, values) in cards:
            ##print('SECSET1', seid, comp, values)
            #assert len(values) > 0, 'seid=%s comp=%s thru_flag=%s values=%s' % (seid, comp, thru_flag, values)
            #fields = ['RELEASE', seid, comp] + values
            #op2.reject_lines.append(print_card_16(fields))

        nfields = len(ints)
        i = 0
        cards = []
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
                cards.append((seid, comp, values))
            elif thru_flag == -1:
                cards.append((seid, comp, ['ALL']))
            else:
                raise NotImplementedError(thru_flag)

        for (seid, comp, values) in cards:
            fields = ['RELEASE', seid, comp] + values
            op2.reject_lines.append(print_card_16(fields))
        return len(data)

    def read_rpnom(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping RPNOM in GEOM4')
        return len(data)

    def read_rrod(self, data: bytes, n: int) -> int:
        """common method for reading RROD"""
        op2: OP2Geom = self.op2
        n = op2.reader_geom2._read_dual_card(
            data, n,
            self._read_rrod_nx, self._read_rrod_msc,
            'RROD', self.op2.reader_geom3._add_op2_rigid_element)
        return n

    def _read_rrod_nx(self, data, n):
        """RROD(6501,65,291) - Record 30"""
        op2: OP2Geom = self.op2
        struct_5i = Struct(mapfmt(op2._endian + b'5i', self.size))
        ntotal = 20 * self.factor
        nelements = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_5i.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  RROD=%s\n' % str(out))
            (eid, ga, gb, cma, cmb) = out
            assert eid > 0, out
            out = (eid, ga, gb, cma, cmb, 0.0) # alpha
            elem = RROD.add_op2_data(out)
            elements.append(elem)
            n += ntotal
        op2.to_nx(' because RROD-NX was found')
        return n, elements

    def _read_rrod_msc(self, data, n):
        """RROD(6501,65,291) - Record 30"""
        op2: OP2Geom = self.op2
        s = Struct(mapfmt(op2._endian + b'5if', self.size))
        ntotal = 24 * self.factor
        nelements = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  RROD=%s\n' % str(out))
            #(eid, ga, gb, cma, cmb, alpha) = out
            assert out[0] > 0, out
            elem = RROD.add_op2_data(out)
            elements.append(elem)
            n += ntotal
        return n, elements

    def read_rspline(self, data: bytes, n: int) -> int:
        """RSPLINE(7001,70,186) - Record 31"""
        self.op2.log.info('geom skipping RSPLINE in GEOM4')
        return len(data)

    def read_rsscon(self, data: bytes, n: int) -> int:
        """RSSCON(7201,72,398) - Record 32"""
        self.op2.log.info('geom skipping RSSCON in GEOM4')
        return len(data)

    def read_rweld(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping RWELD in GEOM4')
        return len(data)

    def read_secset(self, data: bytes, n: int) -> int:
        """
        Word Name Type Description
        1 SEID I Superelement identification number
        2 ID I Grid or scalar point identification number
        3 C I Component numbers
        """
        op2: OP2Geom = self.op2
        ints = np.frombuffer(data[n:], dtype=op2.idtype8)
        nints = len(ints)
        assert nints % 3 == 0, nints
        ints = ints.reshape(nints // 3, 3)
        for seid, nid, component in ints:
            op2.add_secset(seid, [nid], [str(component)])
        return len(data)

    def read_sebset(self, data: bytes, n: int) -> int:
        op2: OP2Geom = self.op2
        ints = np.frombuffer(data[n:], dtype=op2.idtype8)
        nints = len(ints)
        assert nints % 3 == 0, nints
        ints = ints.reshape(nints // 3, 3)
        for seid, nid, component in ints:
            op2.add_sebset(seid, [nid], [str(component)])
        return len(data)

    def read_sebset1(self, data: bytes, n: int) -> int:
        #asdf
        #self.op2.log.info('geom skipping SEBSET1 in GEOM4')
        op2: OP2Geom = self.op2
        ints = np.frombuffer(data[n:], dtype=op2.idtype8)
        i, cards = ints_to_secset1s('SEBSET1', ints)
        for (seid, comp, values) in cards:
            #print('SEBSET1', seid, comp, values)
            op2.add_sebset1(seid, values, comp)
        return len(data)

    def read_secset1(self, data: bytes, n: int) -> int:
        """
        (1010, 10, 320,
         103, 123456, 0, 1, 101, -1)

        SECSET
        Word Name Type Description
        1 SEID I Superelement identification number
        2 ID I Grid or scalar point identification number
        3 C I Component numbers

        Word Name Type Description
        1 SEID I Superelement identification number
        2 C    I Component numbers
        3 THRUFLAG I Thru range flag
        THRUFLAG=0 No
           4 ID I Grid or scalar point identification number
           Word 4 repeats until End of Record
        THRUFLAG=1 Yes
           4 ID1 I First grid or scalar point identification number
           5 ID2 I Second grid or scalar point identification number
        End THRUFLAG
        """
        op2: OP2Geom = self.op2
        ints = np.frombuffer(data[n:], dtype=op2.idtype8)
        i, cards = ints_to_secset1s('SECSET1', ints)
        for (seid, comp, values) in cards:
            #print('SECSET1', seid, comp, values)
            op2.add_secset1(seid, values, str(comp))
        return len(data)

    def read_seqset(self, data: bytes, n: int) -> int:
        """SEQSET(1110,11,321) - Record 40"""
        op2: OP2Geom = self.op2
        ints = np.frombuffer(data[n:], dtype=op2.idtype8)
        nints = len(ints)
        assert nints % 3 == 0, nints
        ints = ints.reshape(nints // 3, 3)
        for seid, nid, component in ints:
            op2.add_seqset(seid, [nid], [str(component)])
        return len(data)
        #return self._read_xset(data, n, 'SEQSET', SEQSET, self.add_SEQSET)

    def read_seqset1(self, data: bytes, n: int) -> int:
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
        nbytes = self._read_superxset1(
            data, n,
            'SEQSET1', SEQSET1, self.op2._add_methods._add_seqset_object,
            debug=True)
        #for seqset in self.se_qsets:
            #print(seqset)
        return nbytes

    def read_sesup(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping SESUP in GEOM4')
        return len(data)

    def read_seuset(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping SEUSET in GEOM4')
        asdf
        return len(data)

    def read_seuset1(self, data: bytes, n: int) -> int:
        op2: OP2Geom = self.op2
        ints = np.frombuffer(data[n:], dtype=op2.idtype8)
        i, cards = ints_to_secset1s('SEUSET1', ints)
        for (seid, comp, values) in cards:
            op2.add_seuset1(seid, values, comp)
        #self.op2.log.info('geom skipping SEUSET1 in GEOM4')
        return len(data)

    def read_spcoff(self, data: bytes, n: int) -> int:
        """SPCOFF(5501,55,16) - Record 44"""
        op2: OP2Geom = self.op2
        ntotal = 16 * self.factor
        nentries = (len(data) - n) // ntotal
        struct_3if = Struct(op2._endian + b'iiif')
        for unused_i in range(nentries):
            edata = data[n:n + 16]
            (sid, ID, c, dx) = struct_3if.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('SPCOFF sid=%s id=%s c=%s dx=%s\n' % (sid, ID, c, dx))
            constraint = SPCOFF.add_op2_data([sid, ID, c, dx])
            op2._add_methods._add_constraint_spcoff_object(constraint)
            n += 16
        return n

    def read_spc(self, data: bytes, n: int) -> int:
        """common method for reading SPCs"""
        op2: OP2Geom = self.op2
        n = op2.reader_geom2._read_dual_card(
            data, n,
            self._read_spc_nx, self._read_spc_msc,
            'SPC', op2._add_methods._add_constraint_spc_object)
        return n

    def _read_spc_msc(self, data, n: int):
        """SPC(5501,55,16) - Record 44

        1 SID   I    Set identification number
        2 ID    I    Grid or scalar point identification number
        3 C     I    Component numbers
        4 UNDEF none Not used
        5 D     RX   Enforced displacement

        """
        op2: OP2Geom = self.op2
        #log = op2.log
        #log.debug('read_spc_mpc')
        ntotal = 20 * self.factor
        nentries = (len(data) - n) // ntotal
        assert nentries > 0, nentries
        assert (len(data) - n) % ntotal == 0
        #self.show_data(data, types='if')

        constraints = []
        struc = Struct(op2._endian + b'iiiif')
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            (sid, nid, comp, xxx, dx) = struc.unpack(edata)
            assert xxx == 0, xxx
            if op2.is_debug_file:
                op2.binary_debug.write('SPC-MSC sid=%s id=%s comp=%s dx=%s\n' % (
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
        op2: OP2Geom = self.op2
        #op2.log.debug('read_spc_nx')
        msg = ''
        ntotal = 16 * self.factor
        ndatai = len(data) - n
        nentries = ndatai // ntotal
        assert nentries > 0, nentries
        assert ndatai % ntotal == 0
        #self.show_data(data, types='if')

        struc = Struct(mapfmt(op2._endian + b'iiif', self.size))
        constraints = []

        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            (sid, nid, comp, dx) = struc.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('SPC-NX sid=%s nid=%s comp=%s dx=%s\n' % (
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
            #op2.log.warning('Invalid Node IDs; skipping\n' + msg)
        #op2.to_nx()
        return n, constraints

    def read_spcoff1(self, data: bytes, n: int) -> int:
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
        idata = np.frombuffer(data[n:], self.op2.idtype).copy()
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
        op2: OP2Geom = self.op2
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
        if op2.is_debug_file:
            op2.binary_debug.write('SPCOFF1: components=%s thru_flag=%s' % (
                components, thru_flag))
            op2.binary_debug.write('   nids=%s\n' % str(nids))
        if len(nids) == 0:
            #op2.log.warning('geom skipping SPC1 because its empty...%s' % out)
            return
        if max(nids) > 100000000:
            self._is_long_ids = True
        in_data = [components, nids]
        constraint = SPCOFF1.add_op2_data(in_data)
        op2._add_methods._add_constraint_spcoff_object(constraint)
        op2.increase_card_count('SPCOFF1', 1)
        if thru_check and len(out) > 5:
            #card = out[5:]
            self._add_spcoff1_card(out[5:])

    def read_spc1(self, data: bytes, n: int) -> int:
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
        idata = np.frombuffer(data[n:], self.op2.idtype8).copy()
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

    def read_spcb(self, data: bytes, n: int) -> int:
        r"""
        SPC(5501, 55, 620016) - Record 45 - MSC

        1     SID   I    Set identification number
        2     ID    I    Grid or scalar point identification number
        3     C     I    Component numbers
        4/5/6 UNDEF none Not used
        4/5/6 D     RX   Enforced displacement

        odd case = (
            # sid, id, comp, [undef, undef, undef]...one of the undef is dx...
            1, 12, 123456, 0, 0, 0,
            1, 13, 123456, 0, 0, 0,
            1, 14, 123456, 0, 0, 0,
            1, 15, 123456, 0, 0, 0,
            1, 16, 123456, 0, 0, 0,
            ....
        )
        SPC     203     10000   3               10002   3       .01

        """
        assert self.size == 4, self.size
        op2: OP2Geom = self.op2
        idata = np.frombuffer(data[n:], 'int32').copy()
        nints = len(idata)
        nrows = nints // 6
        idata = idata.reshape(nrows, 6)[:, :3]
        enforceds = np.frombuffer(data[n:], 'float64').copy()[2::3]

        for idatai, enforced in zip(idata, enforceds):
            sid, idi, comp = idatai
            op2.add_spc(sid, [idi], [str(comp)], [enforced])
        return len(data)

    def _add_spc1_card(self, out):
        """helper method for ``_read_spc1``"""
        op2: OP2Geom = self.op2
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
        if op2.is_debug_file:
            op2.binary_debug.write('SPC1: sid=%s components=%s thru_flag=%s' % (
                sid, components, thru_flag))
            op2.binary_debug.write('   nids=%s\n' % str(nids))
        if len(nids) == 0:
            #op2.log.warning('geom skipping SPC1 because its empty...%s' % out)
            return
        in_data = [sid, components, nids]
        constraint = SPC1.add_op2_data(in_data)
        op2._add_methods._add_constraint_spc_object(constraint)
        op2.increase_card_count('SPC1', 1)
        if thru_check and len(out) > 5:
            #card = out[5:]
            self._add_spc1_card(out[5:])

    def read_spcadd(self, data: bytes, n: int) -> int:
        """SPCADD(5491,59,13) - Record 46"""
        #nentries = (len(data) - n) // 4
        op2: OP2Geom = self.op2
        datai = np.frombuffer(data[n:], op2.idtype8).copy()
        _read_spcadd_mpcadd(op2, 'SPCADD', datai)
        return len(data)

    def read_spcd(self, data: bytes, n: int) -> int:
        """common method for reading SPCDs"""
        op2: OP2Geom = self.op2
        n = op2.reader_geom2._read_dual_card(
            data, n,
            self._read_spcd_nx, self._read_spcd_msc,
            'SPCD', op2._add_methods._add_load_object)
        return n

    def _read_spcd_nx(self, data, n):
        """SPCD(5110,51,256) - NX specific"""
        op2: OP2Geom = self.op2
        struct_3if = Struct(op2._endian + b'3if')
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
            if op2.is_debug_file:
                op2.binary_debug.write('  SPCD-NX=%s\n' % str(out))

            if node_id > 100000000:
                is_long_ids = True
            constraint = SPCD.add_op2_data([sid, node_id, c, dx])
            constraints.append(constraint)
            n += ntotal
        self._is_long_ids = is_long_ids
        #op2.to_nx(' because SPCD-NX was found')
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
        op2: OP2Geom = self.op2
        ntotal = 20 * self.factor  # 5*4
        nentries = (len(data) - n) // ntotal
        assert nentries > 0, nentries
        assert (len(data) - n) % ntotal == 0
        constraints = []
        struct_4if = Struct(op2._endian + b'4if')
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = struct_4if.unpack(edata)
            (sid, ID, c, xxx, dx) = out
            assert xxx == 0, xxx

            if op2.is_debug_file:
                op2.binary_debug.write('  SPCD-MSC=%s\n' % str(out))
            constraint = SPCD.add_op2_data([sid, ID, c, dx])
            constraints.append(constraint)
            n += ntotal
        return n, constraints

    def read_spcde(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping SPCDE in GEOM4')
        return len(data)

    def read_spcf(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping SPCDF in GEOM4')
        return len(data)

    def read_spcdg(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping SPCDG in GEOM4')
        return len(data)

    def read_spce(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping SPCE in GEOM4')
        return len(data)

    def read_spceb(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping SPCEB in GEOM4')
        return len(data)

    def read_spcfb(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping SPCFB in GEOM4')
        return len(data)

    def read_spcgb(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping SPCGB in GEOM4')
        return len(data)

    def read_spcgrid(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping SPCGRID in GEOM4')
        return len(data)

    def read_suport(self, data: bytes, n: int) -> int:
        """SUPORT(5601,56, 14) - Record 59"""
        op2: OP2Geom = self.op2
        ntotal = 8 * self.factor  # 2*4
        nentries = (len(data) - n) // ntotal
        struct_2i = Struct(op2._endian + mapfmt(b'2i', self.size))
        for unused_i in range(nentries):
            out = list(struct_2i.unpack(data[n:n + ntotal]))
            if op2.is_debug_file:
                op2.binary_debug.write('  SUPORT=%s\n' % str(out))
                #op2.log.info(out)
            suport = SUPORT.add_op2_data(out)
            op2._add_methods._add_suport_object(suport) # extracts [sid, c]
            n += ntotal
        return n

    def read_suport1(self, data: bytes, n: int) -> int:
        """SUPORT1(10100,101,472) - Record 60"""
        op2: OP2Geom = self.op2
        nfields = (len(data) - n) // 4 - 2
        out = unpack(op2._endian + b'%ii' % nfields, data[n:n+nfields*4])

        i = 0
        nsuports = 0
        suport: list[int] = []
        while i < len(out):
            if out[i] == -1:
                assert out[i+1] == -1, out
                suporti = SUPORT1.add_op2_data(suport)
                #op2.log.info(suporti)
                op2._add_methods._add_suport_object(suporti) # extracts [sid, nid, c]
                nsuports += 1
                if op2.is_debug_file:
                    op2.binary_debug.write('  SUPORT1=%s\n' % str(suport))
                suport = []
                i += 2
                continue
            suport.append(out[i])
            i += 1
            assert -1 not in suport, suport

        if op2.is_debug_file:
            op2.binary_debug.write('  SUPORT1=%s\n' % str(suport))

        suporti = SUPORT1.add_op2_data(suport)
        op2._add_methods._add_suport_object(suporti) # extracts [sid, nid, c]
        nsuports += 1
        op2.card_count['SUPORT1'] = nsuports
        assert n+nfields*4+8 == len(data), 'a=%s b=%s' % (n+nfields*4+8, len(data))
        return len(data)

    def read_tempbc(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping TEMPBC in GEOM4')
        return len(data)

    def read_uset(self, data: bytes, n: int) -> int:
        """
        USET(2010,20,193) - Record 63
        (sid, nid, comp), ...

        """
        op2: OP2Geom = self.op2
        struct_3i = Struct(mapfmt(op2._endian + b'3i', self.size))
        ntotal = 12 * self.factor
        #self.show_data(data, types='is')
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_3i.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  USET=%s\n' % str(out))
            #(sid, id, component) = out
            set_obj = USET.add_op2_data(out)
            op2._add_methods._add_uset_object(set_obj)
            n += ntotal
        op2.increase_card_count('USET', len(op2.usets))
        return n

    def read_seqset1b(self, data: bytes, n: int) -> int:  # pragma: no cover
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
        self.op2.show_data(data)
        return len(data)

    def read_uset1(self, data: bytes, n: int) -> int:
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
        op2: OP2Geom = self.op2
        assert self.factor == 1, self.factor
        nentries = 0
        size = self.size
        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        floats = np.frombuffer(data[n:], op2.fdtype8).copy()
        i = 0
        #print('ints = %s' % ints)
        nidata = len(ints)
        while i < nidata:
            unused_sname = data[n+i*size : n+(i+1)*size]
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

            if op2.is_debug_file:
                op2.binary_debug.write('USET1: sname=%s comp=%s thru_flag=%s' % (
                    sname_str, comp, thru_flag))
                op2.binary_debug.write('   nids=%s\n' % str(nids))
            #print('USET1: sid=%s comp=%s thru_flag=%s' % (
            #    sid, comp, thru_flag))
            #print('   nids=%s\n' % str(nids))
            in_data = [sname_str, comp_str, nids]
            #print(in_data)

            constraint = USET1.add_op2_data(in_data)
            op2._add_methods._add_uset_object(constraint)
        op2.card_count['USET1'] = nentries
        return len(data)


    def read_omit(self, data: bytes, n: int) -> int:
        return self._read_xset(data, n, 'OMIT', OMIT, self.op2._add_methods._add_omit_object)

    def read_rtrplt(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping RTRPLT in GEOM4')
        return len(data)

    def read_bndfix(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping BNDFIX in GEOM4')
        return len(data)

    def read_bndfix1(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping BNDFIX1 in GEOM4')
        return len(data)

    def read_bndfree(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping BNDFREE in GEOM4')
        return len(data)

    def read_bltmpc(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping BLTMPC in GEOM4')
        return len(data)

def read_rbe2s_from_idata_fdata(op2: OP2Geom, idata, fdata) -> list[RBE2]:
    """
    Word Name Type Description
    1 EID I Element identification number
    2  GN I Grid point identification number for independent degrees-of-freedom
    3  CM I Component numbers of dependent degrees of-freedom
    4  GM I Grid point identification number for dependent degrees-of-freedom
    Word 4 repeats until End of Record
    5 ALPHA RS Thermal expansion coefficient  (optional)
    6 TREF  RS Reference temperature          (optional if alpha exists)

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

    idata = [
        132, 253, 123456,
            4, 5, 6, 7, 8, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 57, 58, 59, 60,
            -1, 0,
        133, 254, 123456,
            126, 127, 128, 129, 138, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214,
            -1, 0]
    i=[ 0, 30]
    j=[29, 59]

    idata = [
        132, 253, 123456,
            4, 5, 6, 7, 8, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 57, 58, 59, 60,
            -1, 0, 0,
        133, 254, 123456,
            126, 127, 128, 129, 138, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214,
            -1, 0, 0,
    ]
    i=[ 0, 31]
    j=[29, 61]

    idata = [
        201, 5001, 123456, 201, -1,
        221, 5002, 123456, 221, -1,
        5001, 5001, 123, 1001, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1118, 1119, 1120, 1121, 1122, 1123, 1124, 1125, 1126, 1127, 1128, 1129, 1226, 1227, 1228, 1229, 1230, 1231, 1232, 1233, 1234, 1235, 1236, 1237, -1, 5002, 5002, 123, 2097, 2106, 2107, 2108, 2109, 2110, 2111, 2112, 2113, 2114, 2115, 2116, 2117, 2214, 2215, 2216, 2217, 2218, 2219, 2220, 2221, 2222, 2223, 2224, 2225, 2322, 2323, 2324, 2325, 2326, 2327, 2328, 2329, 2330, 2331, 2332, 2333, -1,
        22222, 4455, 123456, 4467, 4443, -1
    ]
    """
    rbe2s = []
    geom3 = op2.reader_geom3
    iminus1 = np.where(idata == -1)[0]
    is_alpha = False
    is_alpha_tref = False
    if idata[-1] == -1:
        is_alpha = False
        i = np.hstack([[0], iminus1[:-1]+1])
        j = np.hstack([iminus1[:-1], len(idata)])  # the index for the last integer
    else:
        if idata[-2] == -1:
            is_alpha = True
            i = np.hstack([[0], iminus1[:-1]+2])
            j = np.hstack([iminus1[:-1]+1, len(idata)-1])
        else:
            is_alpha_tref = True
            i = np.hstack([[0], iminus1[:-1]+3])
            j = np.hstack([iminus1[:-1]+1, len(idata)-2])
            assert is_alpha_tref is True, is_alpha_tref
        assert -1 not in idata[j], idata[j]
    assert -1 not in idata[i], idata[i]
    #print(is_alpha, is_alpha_tref)
    #print('is_alpha=%s' % is_alpha)
    #print('i=%s' % i)
    #print('j=%s' % j)
    #print('idata=%s' % idata.tolist())
    #print(fdata, len(fdata))
    nelements = len(j)
    alpha = 0.0
    for ii, jj in zip(i, j):
        eid, gn, cm = idata[ii:ii + 3]
        if is_alpha or is_alpha_tref:
            gm = idata[ii+3:jj-1].tolist()
            #print('eid=%s gn=%s cm=%s gm=%s' % (eid, gn, cm, gm))
            alpha = fdata[jj]
            if is_alpha_tref:
                tref = fdata[jj+1]
        else:
            gm = idata[ii+3:jj].tolist()
            if -1 in gm:
                gm = gm[:-1]
        assert -1 not in gm, f'eid={eid} gn={gn} cm={cm} gm={gm}'
        #assert len(gm) > 0, gm

        out = [eid, gn, cm, gm, alpha]
        if is_alpha_tref:
            out.append(tref)
        if op2.is_debug_file:
            op2.binary_debug.write('  RBE2=%s\n' % str(out))
        #print('rbe2 =', out)
        elem = RBE2.add_op2_data(out)
        geom3._add_op2_rigid_element(elem)
        rbe2s.append(elem)
    op2.card_count['RBE2'] = nelements
    return rbe2s

def read_rbe3s_from_idata_fdata(op2: OP2Geom, idata, fdata) -> list[RBE3]:
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

    idata = [130, 251, 123456, 1.0, 123, 12, 13, 14, 15, 16, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89,
             -1, -2, -3, 0.0, 0.0,
             131, 252, 123456, 1.0, 123, 129, 130, 131, 132, 133, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 182, 183, 184, 185,
             -1, -2, -3, 0.0, 0.0]

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

    is_alpha = False
    is_alpha_tref = False
    if idata[-1] == -3:
        i = np.hstack([[0], iminus3[:-1]+1])
    else:
        if idata[-2] == -3:
            is_alpha = True
            i = np.hstack([[0], iminus3[:-1]+2])
        else:
            is_alpha_tref = True
            assert idata[-3] == -3, idata
            i = np.hstack([[0], iminus3[:-1]+3])
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

        if is_alpha_tref:
            alpha = fdata[ii]
            tref = fdata[ii+1]
            if tref != 0.0:
                op2.log.warning(f'new RBE3 field={tref} (assuming float)')
            ii += 2
        elif is_alpha:
            alpha = fdata[ii]
            ii += 1
        else:
            alpha = 0.0
        #print('alpha = %s' % alpha)
        in_data = [eid, refg, refc, weights, comps, gijs,
                   gmi, cmi, alpha]
        if is_alpha_tref:
            in_data.append(tref)
        if op2.is_debug_file:
            op2.binary_debug.write('  RBE3=%s\n' % str(in_data))
        #print('rbe3 =', in_data)
        rbe3 = RBE3.add_op2_data(in_data)
        #print(rbe3.rstrip())
        #if eid in self.rigid_elements:
            #old_rigid = self.rigid_elements[eid]
            #op2.log.warning(f'geom skipping RBE3 eid={eid} because its duplicated\n{rbe3} by:\n{old_rigid}')
            #continue
        op2.reader_geom3._add_op2_rigid_element(rbe3)
        rbe3s.append(rbe3)
        #print('--------------------------------------')
    return rbe3s

def _get_rbe3_um(i: int, unused_j, idata, unused_fdata):
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

def get_minus_2_index(idata) -> int:
    """helper for ``get_minus_2_index``"""
    #print('idata =', idata)
    i = np.where((idata == -2) | (idata == -3))[0]
    if len(i) == 0:
        return len(idata)
    return i[0]

def fill_rbe3_wt_comp_gijs(i: int, j: int,
                           idata, fdata) -> tuple[int, list[float], list[int], list[list[int]]]:
    """helper for ``read_rbe3s_from_idata_fdata``"""
    weights: list[float] = []
    comps: list[int] = []
    grids: list[list[int]] = []

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

def _read_spcadd_mpcadd(model: OP2Geom, card_name: str, datai: np.ndarray) -> None:
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
        the data array; cannot be a list[int]
        [2  1 10 -1]
        [3  1 -1]

    """
    add_methods = model._add_methods
    if model.is_debug_file:
        model.binary_debug.write('  %s - %s' % (card_name, str(datai)))
    iend = np.where(datai == -1)[0]
    if len(datai) == 3:
        dataii = datai[:-1]
        if card_name == 'MPCADD':
            constraint = MPCADD.add_op2_data(dataii)
            add_methods._add_constraint_mpcadd_object(constraint)
        elif card_name == 'SPCADD':
            constraint = SPCADD.add_op2_data(dataii)
            add_methods._add_constraint_spcadd_object(constraint)
        elif card_name == 'DCONADD':
            constraint = DCONADD.add_op2_data(dataii)
            add_methods._add_dconstr_object(constraint)
        else:  # pragma: no cover
            raise NotImplementedError(card_name)

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
            add_methods._add_constraint_mpcadd_object(constraint)
        elif card_name == 'SPCADD':
            constraint = SPCADD.add_op2_data(dataii)
            add_methods._add_constraint_spcadd_object(constraint)
        elif card_name == 'DCONADD':
            constraint = DCONADD.add_op2_data(dataii)
            add_methods._add_dconstr_object(constraint)
        else:  # pragma: no cover
            raise NotImplementedError(card_name)
    model.increase_card_count(card_name, count_num=count_num)

def check_component(component: int, msg: str) -> None:
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



def ints_to_secset1s(card_name: str, ints: np.ndarray) -> tuple[int, list[tuple[int, int, list[int]]]]:
    """
    [61,  123456,    1, 610101, 610124]
    [100,    123,    0,     41,     42,  -1,
     200,    123,    0,     35,     36,  -1]
     [ 1, 123456,    1,   1001,   1006,
       1, 123456,    0,   1066,     -1]

    """
    iword = 1
    i = 0
    cards = []
    while i < len(ints):
        #print(i, iword, ints[i])

        if iword == 1:
            seid = ints[i]
        elif iword == 2:
            comp = ints[i]
        elif iword == 3:
            thru_flag = ints[i]
        elif iword == 4:
            if thru_flag == 0:
                value = ints[i]
                values = []
                while value != -1:
                    values.append(value)
                    i += 1
                    value = ints[i]
                #print('SECSET1', seid, comp, thru_flag, values)
                #print('A', (seid, comp, values))
                cards.append((seid, comp, values))
                iword = 0
            elif thru_flag == 1:
                value0 = ints[i]
                value1 = ints[i+1]
                i += 1
                values = list(range(value0, value1+1, 1))
                #print('B', (seid, comp, values))
                cards.append((seid, comp, values))
                iword = 0
            else:
                raise NotImplementedError(f'{card_name} thru_flag={thru_flag}')
        else:
            raise NotImplementedError(f'{card_name} iword={iword}')
        i += 1
        iword += 1
    i -= 1
    assert len(ints) >= i, f'nints={len(ints)} i={i}'
    return i, cards

def xset1_break_by_thru_type(data: np.ndarray) -> list[np.ndarray]:
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
