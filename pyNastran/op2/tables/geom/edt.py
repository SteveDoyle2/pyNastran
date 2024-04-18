"""
defines readers for BDF objects in the OP2 EDT/EDTS table
"""
from __future__ import annotations
from struct import Struct
from typing import Any, TYPE_CHECKING

import numpy as np
from pyNastran.bdf.cards.aero.aero import (
    #AECOMP, AECOMPL, AEFACT, AELINK, AELIST, AEPARM, AESURF, AESURFS,
    CAERO1, # CAERO2, CAERO3, CAERO4, CAERO5,
    PAERO1, # PAERO2, PAERO3, PAERO4, PAERO5,
    MONPNT1, MONPNT2, MONPNT3, MONDSP1,
    SPLINE1, SPLINE2, # SPLINE3,
    SPLINE4, SPLINE5)

from pyNastran.op2.errors import DoubleCardError
from pyNastran.op2.op2_interface.op2_reader import mapfmt, reshape_bytes_block_size # reshape_bytes_block,
from pyNastran.bdf.cards.elements.acoustic import ACMODL
from pyNastran.bdf.cards.aero.dynamic_loads import AERO, FLUTTER
from .utils import get_minus1_start_end
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2_geom import OP2Geom

class EDT:
    """defines methods for reading aero and element deformations"""

    @property
    def size(self) -> int:
        return self.op2.size
    @property
    def factor(self) -> int:
        return self.op2.factor

    def _read_fake(self, data: bytes, n: int) -> int:
        return self.op2._read_fake(data, n)

    def read_stop(self, data: bytes, n: int) -> int:
        return self.op2.reader_geom1.read_stop(data, n)

    def read_edt_4(self, data: bytes, ndata: int):
        """
        3.21 EDT
        Aero and element deformations.

        """
        return self.op2._read_geom_4(self.edt_map, data, ndata)

    def __init__(self, op2: OP2Geom):
        self.op2 = op2

        # F:\Program Files\Siemens\NXNastran\nxn10p1\nxn10p1\nast\tpl\fsw_eng.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_boltld04i.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_eliter17.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_weld01i.op2
        # F:\work\pyNastran\examples\Dropbox\move_tpl\ac10901a_new.op2
        self.edt_map = {
            (5201, 52, 373) : ['ACMODL', self.read_acmodl],
            (6301, 63, 397) : ['ADAPT', self._read_fake],
            (7801, 78, 582) : ['AECOMP', self.read_aecomp],
            (7901, 79, 583) : ['AECOMPL', self.read_aecompl],
            (7301, 73, 574) : ['AEDW', self._read_fake],
            (4002, 40, 273) : ['AEFACT', self.read_aefact],
            (7501, 75, 576) : ['AEFORCE', self.read_aeforce],
            (2602, 26, 386) : ['AELINK', self.read_aelink],
            (2302, 23, 341) : ['AELIST', self.read_aelist],
            (7001, 70, 571) : ['AEPARM', self.read_aeparm],
            (7401, 74, 575) : ['AEPRESS', self.read_aepress],
            (3202, 32, 265) : ['AERO', self.read_aero],
            (2202, 22, 340) : ['AEROS', self.read_aeros],
            (2102, 21, 339) : ['AESTAT', self.read_aestat],
            (2002, 20, 338) : ['AESURF', self.read_aesurf],
            (7701, 77, 581) : ['AESURFS', self.read_aesurfs],
            (3002, 30, 263) : ['CAERO1', self.read_caero1],
            (4301, 43, 167) : ['CAERO2', self.read_caero2],
            (4401, 44, 168) : ['CAERO3', self.read_caero3],
            (4501, 45, 169) : ['CAERO4', self.read_caero4],
            (5001, 50, 175) : ['CAERO5', self.read_caero5],
            (6201, 62, 143) : ['CLOAD', self._read_fake],
            (6401, 64, 307) : ['CSSCHD', self.read_csschd],
            (104, 1, 81) : ['DEFORM', self.read_deform],
            (2702, 27, 387) : ['DIVERG', self.read_diverg],
            (4102, 41, 274) : ['FLFACT', self.read_flfact],
            (3902, 39, 272) : ['FLUTTER', self.read_flutter],
            (17400, 174, 616) : ['GROUP', self.read_group],
            (3802, 38, 271) : ['MKAERO1', self.read_mkaero1],
            (3702, 37, 270) : ['MKAERO2', self.read_mkaero2],
            (7601, 76, 577) : ['MONPNT1', self.read_monpnt1],
            (3102, 31, 264) : ['PAERO1', self.read_paero1],
            (4601, 46, 170) : ['PAERO2', self.read_paero2],
            (4701, 47, 171) : ['PAERO3', self.read_paero3],
            (4801, 48, 172) : ['PAERO4', self.read_paero4],
            (5101, 51, 176) : ['PAERO5', self.read_paero5],
            (5301, 53, 378) : ['PANEL', self.read_panel],
            (3502, 35, 268) : ['SET1', self.read_set1],
            (3602, 36, 269) : ['SET2', self.read_set2],
            (4302, 43, 607) : ['SET3', self.read_set3],
            (3302, 33, 266) : ['SPLINE1', self.read_spline1],
            (3402, 34, 267) : ['SPLINE2', self.read_spline2],
            (4901, 49, 173) : ['SPLINE3', self.read_spline3],
            (6501, 65, 308) : ['SPLINE4', self.read_spline4],
            (6601, 66, 309) : ['SPLINE5', self.read_spline5],
            (2402, 24, 342) : ['TRIM', self.read_trim],
            (7201, 72, 573) : ['UXVEC', self.read_uxvec],
            (7108, 822, 51) : ['BOLT', self._read_fake],
            (7108, 71, 251) : ['???', self._read_fake],
            (5808, 58, 220) : ['ITER', self._read_fake],
            (14000, 140, 568) : ['SWLDPRM', self._read_fake],
            (11001, 110, 581) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],

            (7108, 82, 251): ['BOLT', self._read_fake],

            # MSC
            #(1247, 12, 667): ['MONPNT2', self._read_monpnt2],
            (11204, 112, 821): ['ERPPNL', self._read_fake],
            (8001, 80, 511): ['SET3', self.read_set3],
            (9400, 94, 641): ['MDLPRM', self.read_mdlprm],
            (11004, 110, 1820_720): ['HADACRI', self._read_fake],
            (8804, 88, 628): ['MONDSP1', self.read_mondsp1],

            (10904, 109, 719): ['HADAPTL', self._read_fake],
            (8204, 82, 621): ['MONPNT2', self.read_monpnt2],
            (8304, 83, 622): ['MONPNT3', self.read_monpnt3],
            (1247, 12, 667): ['MONPNT2', self.read_monpnt2], # nx
            #(8001, 80, 511): ['???', self._read_fake],
            #(8001, 80, 511): ['???', self._read_fake],
            #(8001, 80, 511): ['???', self._read_fake],
            #(8001, 80, 511): ['???', self._read_fake],
            #(8001, 80, 511): ['???', self._read_fake],
            (9004, 90, 619): ['MASSSET', self.read_massset],
            (10104, 101, 683): ['SPBLND1', self.read_sblnd1],
            (10004, 100, 682): ['SPLINRB', self._read_fake],
            (10404, 104, 686): ['???', self.read_stop],
            (12804, 128, 945): ['???', self.read_stop],

        }
    def read_sblnd1(self, data: bytes, n: int) -> int:
        """
        24 - SPBLND1(10104,101,683)
        Word Name Type Description
        1 EID           I Identification of Blended Spline
        2 EID1          I Ident. of 1st Spline to be Blended
        3 EID2          I Ident. of 2nd Spline to be Blended
        4 OPTION(2) CHAR4 Blending Option: WAVG, LIN, CUB
        6 W1           RS Weight for 1st Spline
        7 GID           I Aerodynamic Reference Grid
        8 D1           RS Blending Depth of 1st Spline
        9 D2           RS Blending Depth of 2nd Spline
        10 X1          RS X1 Component of Direction Vector
        11 X2          RS X2 Component of Direction Vector
        12 X3          RS X3 Component of Direction Vector
        13 CID          I Coordinate System Identification Number
        """
        op2: OP2Geom = self.op2
        ntotal = 13 * self.size # 4*13
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        assert self.factor == 1, self.factor
        structi = Struct(op2._endian + b'3i 8s fi 5fi')
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)

            eid, eid1, eid2, option_bytes, w1, aero_grid, d1, d2, x1, x2, x3, cid = out
            option = reshape_bytes_block_size(option_bytes, size=self.size)
            x = [x1, x2, x3]
            spline_blend = op2.add_sblnd1(eid, eid1, eid2, option_bytes, w1, aero_grid, d1, d2, x, cid)
            str(spline_blend)
            n += ntotal
        return n

    def read_massset(self, data: bytes, n: int) -> int:
        """
        51 - MASSSET(9004,90,602)
        Word Name Type Description
        1 MID     I massset id
        2 S      RS Overall scale factor
        3 SI     RS massid scale factor
        4 MASSID  I massid
        Words 3 through 4 repeat until (-1,-1) occurs
        """
        op2: OP2Geom = self.op2
        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        floats = np.frombuffer(data[n:], op2.fdtype8).copy()
        iminus1 = np.where(ints == -1)[0]
        iminus1_start = iminus1[::2]
        iminus1_end = iminus1[1::2]

        #ncards = 0
        istart = [0] + list(iminus1_end + 1)
        iend = iminus1_start
        for (i0, i1) in zip(istart, iend):
            assert ints[i1] == -1, ints[i1]
            mid = ints[i0]
            scale = floats[i0+1]
            scales = floats[i0+2:i1:2].tolist()
            massids = ints[i0+3:i1:2].tolist()
            assert len(scales) == len(massids)
            assert len(scales) > 0, scales
            op2.add_massset(mid, scale, scales, massids)
        return len(data)

    def read_aeforce(self, data: bytes, n: int) -> int:
        """Word Name Type Description
        1 MACH     RS
        2 SYMXZ(2) CHAR4
        4 SYMXY(2) CHAR4
        6 UXID     I
        7 MESH(2)  CHAR4
        9 FORCE    I
        10 DMIK(2) CHAR4
        12 PERQ(2) CHAR4
        """
        op2: OP2Geom = self.op2
        ntotal = 52 * self.factor # 4*13
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        assert self.factor == 1, self.factor
        structi = Struct(op2._endian + b'f 8s 8s i 8s i 8s 8s')
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)

            mach, sym_xz_bytes, sym_xy_bytes, ux_id, mesh_bytes, force, dmik_bytes, perq_bytes = out
            sym_xz = reshape_bytes_block_size(sym_xz_bytes, size=self.size)
            sym_xy = reshape_bytes_block_size(sym_xy_bytes, size=self.size)
            mesh = reshape_bytes_block_size(mesh_bytes, size=self.size)
            dmik = reshape_bytes_block_size(dmik_bytes, size=self.size)
            perq = reshape_bytes_block_size(perq_bytes, size=self.size)
            assert isinstance(mesh, str), mesh
            assert isinstance(dmik, str), dmik
            assert isinstance(perq, str), perq

            aeforce = op2.add_aeforce(mach, sym_xz, sym_xy, ux_id, mesh, force, dmik, perq)
            str(aeforce)
            n += ntotal
        return n

    def read_aepress(self, data: bytes, n: int) -> int:
        """
        Parametric pressure loading for aerodynamics.
        Word Name Type Description
        1 MACH     RS    Mach number
        2 SYMXZ(2) CHAR4 Character string for identifying symmetry of the
                         force vector. Allowable values are SYMM, ASYMM, and ANTI
        4 SYMXY(2) CHAR4 Character string for identifying symmetry of the
                         force vector. Allowable values are SYMM, ASYMM, and ANTI
        6 UXID     I     The identification number of a UXVEC entry
        7 DMIJ(2)  CHAR4 The name of a DMI or DMIJ entry that defines the pressure
                         per unit dynamic pressure
        9 DMIJI(2) CHAR4 The name of a DMIJI entry that defines the CAERO2
                         interference element downwashes
        """
        op2: OP2Geom = self.op2
        ntotal = 40 * self.factor # 4*10
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        assert self.factor == 1, self.factor
        structi = Struct(op2._endian + b'f 8s 8s i 8s 8s')
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)

            mach, sym_xz_bytes, sym_xy_bytes, ux_id, dmij_bytes, dmiji_bytes= out
            sym_xz = reshape_bytes_block_size(sym_xz_bytes, size=self.size)
            sym_xy = reshape_bytes_block_size(sym_xy_bytes, size=self.size)
            dmij = reshape_bytes_block_size(dmij_bytes, size=self.size)
            dmiji = reshape_bytes_block_size(dmiji_bytes, size=self.size)
            assert isinstance(dmij, str), dmij
            assert isinstance(dmiji, str), dmiji
            aepress = op2.add_aepress(mach, sym_xz, sym_xy, ux_id, dmij, dmiji)
            str(aepress)
            #print(mach, sym_xz, sym_xy, ux_id, dmij, dmiji)
            n += ntotal
        return n

    def read_mkaero2(self, data: bytes, n: int) -> int:
        mkaero2x

    def read_csschd(self, data: bytes, n: int) -> int:
        """
        (6401, 64, 307,
         10, 510, 10, 20, 30,
         20, 510, 100, 20, 300)

        csschd	10	510	10	20	30
        csschd	20	510	100	20	300
        """
        ints = np.frombuffer(data[n:], dtype=self.op2.idtype8)
        nints = len(ints)
        assert nints % 5 == 0
        ints = ints.reshape(nints//5, 5)
        for intsi in ints:
            sid, aesid, lschd, lalpha, lmach = intsi
            self.op2.add_csschd(sid, aesid, lschd, lalpha=lalpha, lmach=lmach)
        return len(data)

    def read_diverg(self, data: bytes, n: int) -> int:
        """
        Record – DIVERG(2702,27,387)
        Divergence analysis data.
        Word Name Type Description
        1 SID     I    Unique set identification number
        2 NROOT   I    Number of divergence roots to output
        3 M       RS   Mach number
        Word 3 repeats until -1 occurs
        """
        op2: OP2Geom = self.op2
        ints = np.frombuffer(data[n:], op2.idtype).copy()
        floats = np.frombuffer(data[n:], op2.fdtype).copy()
        istart, iend = get_minus1_start_end(ints)

        for (i0, i1) in zip(istart, iend):
            sid, nroots = ints[i0:i0+2]
            machs = floats[i0+2:i1]
            #print(sid, nroots, machs)
            assert ints[i1] == -1, ints[i1]
            diverg = op2.add_diverg(sid, nroots, machs)
            str(diverg)
        return len(data)

    def read_flfact(self, data: bytes, n: int) -> int:
        """
        data = (1, 0.206, -1,
                2, 1.3, -1,
                3, 14400.0, 15600.0, 16800.0, 18000.0, 19200.0, 20400.0, -1)
        """
        op2: OP2Geom = self.op2
        ints = np.frombuffer(data[n:], op2.idtype).copy()
        floats = np.frombuffer(data[n:], op2.fdtype).copy()
        istart, iend = get_minus1_start_end(ints)

        for (i0, i1) in zip(istart, iend):
            sid = ints[i0]
            factors = floats[i0+1:i1]
            assert ints[i1] == -1, ints[i1]
            flfact = op2.add_flfact(sid, factors)
            str(flfact)
        return len(data)

    def read_mkaero1(self, data: bytes, n: int) -> int:
        """
        (3802, 38, 271)
        Kinda brilliant way to write the card.  Weird to parse though.

        data = (1.3, -1, -1, -1, -1, -1, -1, -1,
                0.03, 0.04, 0.05, -1, -1, -1, -1, -1)
        """
        op2: OP2Geom = self.op2
        #assert len(data) == 76, len(data)
        nvalues = (len(data) - n) // 4
        nrows = nvalues // 16
        assert nrows > 0, nrows
        ints = np.frombuffer(data[12:], dtype=op2.idtype).reshape(nrows, 16)
        floats = np.frombuffer(data[12:], dtype=op2.fdtype).reshape(nrows, 16)
        irows, icols = np.where(ints != -1)
        uirows = np.unique(irows)
        for irow in uirows:
            iaero = np.where(irows == irow)[0]
            ifloats = icols[iaero]
            imachsi = np.where(ifloats < 8)[0]
            ikfreqsi = np.where(ifloats >= 8)[0]
            imachs = ifloats[imachsi]
            ikfreqs = ifloats[ikfreqsi]

            machs = floats[irow, imachs]
            kfreqs = floats[irow, ikfreqs]
            mkaero1 = op2.add_mkaero1(machs, kfreqs)
            str(mkaero1)
        return len(data)

    def read_group(self, data: bytes, n: int) -> int:
        """
        GROUP(17400,174,616) - NX specific

        1 GID          I Group identification number
        2 NDESC(C)     I Length of group description
        3 GDESC(2) CHAR4 Group description
        Word 3 repeats NDESC times

        NDESC+3 GTYPE I Group type
        -2 = Meta data
        -3 = Property identification numbers
        -4 = Grid identification numbers
        -5 = Element identification numbers
        GTYPE = -2 Meta data
        NDESC+4 NMETA        I Length of meta data (includes -1 terminator)
        NDESC+5 MDESC(2) CHAR4 Meta data
        Word NDESC+5 repeats NMETA times

        GTYPE = -3 Property identification numbers
        NDESC+5
        +NMETA
        ID I Property identification numbers
        > 0 for ID
        = 0 for THRU
        = -6 for BY
        = -7 for ALL
        Word NDESC+5+NMETA repeats until -1 occurs
        GTYPE = -4 Grid identification numbers
        NDESC+5+NMETA:
        ID I Grid identification numbers
        > 0 for ID
        = 0 for THRU
        = -6 for BY
        = -7 for ALL
        Word NDESC+5+NMETA repeats until -1 occurs
        GTYPE = -5 Element identification numbers
        NDESC+5
        +NMETA
        ID I Element identification numbers
        > 0 for ID
        = 0 for THRU
        = -6 for BY
        = -7 for ALL
        Word NDESC+5+NMETA repeats until -1 occurs

        (
            17400, 174, 616,
            6, 0,
                -2, 1, -1,
                -4, 1, 0, 440, -1,
            -1
        )
        (
            17400, 174, 616,
            55, 0,
                -5, 90011, -1,
                -1,
             65, 0,
                -5, 90012, -1,
                -1,
             75, 0,
                -5 90013, -1,
                -1)

        GROUP 10 Assembly AA4
        META 100 RPM
        META Optionally continue the meta data
        GRID 1 2 3 4 5 6 7 8
        GRID 10 THRU 20
        GRID 100 THRU 200
        GRID 341 THRU 360 BY 2
        ELEM 30 THRU 40
        PROP ALL
          strings = (b'o\x00\x00\x00\x05\x00\x00\x00THIS IS GROUP 111   \xfe\xff\xff\xff\x05\x00\x00\x00THIS IS METADATA\xff\xff\xff\xff\xfb\xff\xff\xff\x01\x00\x00\x00\x00\x00\x00\x00\n\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff',)
          ints    = (111, 5, 'THIS IS GROUP 111   ', -2, 5, 'THIS IS METADATA', -1, -5, 1, 0, 10, -1, -1)
          floats  = (111, 5, 'THIS IS GROUP 111   ', -2, 5, 'THIS IS METADATA', -1, -5, 1, 0.0, 10, -1, -1)

        # double
        FLEXIBLE SLIDER(1)
          doubles (float64) = (1, 5, 6.03e-154, 6.08e-154, 6.01-154, 6.01e-154, 6.04e-154,
          -4, 14, -1, -1, 2, 6, 23e-154, 6.08e-154, 6.82e-154, 6.9e-154, 6.3e-154, 6.0e-154,
          -5, 1, 0.0, 10, -1, -1)
          long long (int64) = (1, 5, 234137606, 231102857, 23148032, 955957572, 231888857,
          -4, 14, -1, -1, 2, 6, 2413766, 742102857, 231216032, 23997572, 23192817, 23453545,
          -5, 1, 0,   10, -1, -1)
        """
        op2: OP2Geom = self.op2
        #print('reading group')
        #assert self.factor == 1, self.factor
        nentries = 0
        ints = np.frombuffer(data[n:], dtype=op2.idtype8)
        if self.factor == 1:
            strs = np.frombuffer(data[n:], dtype='|S4')
        else:
            op2.show_data(data[n:], types='qds')
            strs = np.frombuffer(data[n:], dtype='|S8')
        size = self.size
        #print(ints)
        #print(strs)

        i = 0
        #minus1_count = 0
        #minus1 = np.where(ints == -1)[0]
        ndata = len(data)
        while n < ndata:
            #1 GID          I Group identification number
            #2 NDESC(C)     I Length of group description
            #3 GDESC(2) CHAR4 Group description
            #Word 3 repeats NDESC times
            group_id, ndesc = ints[i:i+2]
            i += 2
            n += 8
            group_desc = reshape_bytes_block_size(b''.join(strs[i:i+ndesc]), size=size)
            #if self.factor == 1:
                #group_desc = ''.join(stri.decode('latin1') for stri in strs[i:i+ndesc]).strip()
            #else:
                #group_desc_bytes = reshape_bytes_block(b''.join(strs[i:i+ndesc]))
                #group_desc = group_desc_bytes.decode('latin1').rstrip()
            i += ndesc
            n += 4 * ndesc

            #------------------------------
            #gtype, nmeta, mdesc
            gtype = ints[i]
            #i += 1
            #n += 4
            op2.log.debug(f'group_id={group_id} ndesc={ndesc} group_desc={group_desc!r}; gtype={gtype!r}')

            data_dict = {
                'meta': '',
                'property': [],
                'grid': [],
                'element': [],
            }
            #i += 1
            #n += 4

            while n < ndata:
                #Group type
                #-2 = Meta data
                #-3 = Property identification numbers
                #-4 = Grid identification numbers
                #-5 = Element identification numbers
                #print(f'-----gtype={gtype}----')
                #print(ints[i:])
                if gtype == -1:
                    # end of card
                    i += 1
                    n += size
                    break
                elif gtype == -2:
                    assert ints[i] == -2, ints[i]
                    i += 1
                    n += size

                    # meta-data
                    nmeta = ints[i]
                    assert nmeta >= 0, nmeta
                    #i += 1
                    #n += 4
                    #print(i, nmeta)
                    #print(strs[i:i+nmeta-1])
                    #self.show_data(data[i*4:])
                    istop = i+nmeta-1
                    assert istop > i, f'i={i} nmeta={nmeta}'
                    #print('strs[i:istop] =', strs[i:istop])
                    #print('istop =', istop, ints[istop])
                    #meta_desc = ''.join(stri.decode('latin1') for stri in strs[i:istop])
                    datai = data[n+(i+1)*size:n+istop*size]
                    meta_desc = datai.decode('latin1')
                    data_dict['meta'] = meta_desc
                    i += nmeta + 1
                    n += size * (nmeta + 1)
                    #print(f'  gtype={gtype} nmeta={nmeta} meta_desc={meta_desc!r}')
                    #iminus1 = minus1[minus1_count+2]
                    #print('ints: ', ints[i:iminus1].tolist())
                    #minus1_count += 1
                elif gtype == -3:
                    assert ints[i] == -3, ints[i]
                    i, n, props = _read_group_elem_prop_nids(ints, i, n, size)
                    data_dict['property'].append(props)
                elif gtype == -4:
                    assert ints[i] == -4, ints[i]
                    i, n, grids = _read_group_elem_prop_nids(ints, i, n, size)
                    data_dict['grid'].append(grids)
                elif gtype == -5:
                    assert ints[i] == -5, ints[i]
                    i += 1
                    n += size

                    #print(f'gtype=5 (eids); ints[{i}]={ints[i]}')
                    #self.show_data(data[12:], types='ifs')
                    #print('data', ints[i:].tolist())
                    #GTYPE = -5 Element identification numbers
                    #NDESC+5+NMETA
                    #ID I Element identification numbers
                    #> 0 for ID
                    #= 0 for THRU
                    #= -6 for BY
                    #= -7 for ALL
                    #Word NDESC+5+NMETA repeats until -1 occurs
                    #print('ints[i:] =', ints[i:])
                    assert ints[i:][0] > 0, ints[i:]
                    for j, nj in enumerate(ints[i:]):
                        if nj == -1:
                            break
                    eids_array = ints[i:i+j].tolist()
                    eids2 = _expand_vals(eids_array)
                    # print(f'  eids1 = {eids_array}')
                    #print(f'  eids2 = {eids2}')
                    assert 'THRU' != eids2[0], eids2
                    assert 'BY' != eids2[0], eids2
                    assert 'ALL' != eids2[0], eids2
                    data_dict['element'].append(eids2)
                    nstop = len(eids_array) + 1
                    i += nstop
                    n += nstop * self.size
                else:
                    raise NotImplementedError(gtype)
                gtype = ints[i]
                assert gtype <= -1, ints[i]
                #print(f'***gtype={gtype} (ndata-n)={(ndata-n)}')
                #print('---------------')
                #if gtype == -1 and (ndata - n) == 4:
                    #print('break')
                    #minus1_count += 1
                    #i += 1
                    #n += 4
                    #break
            #grid=1 ndesc=4 group_desc='GROUP(1)_ELAR'
            # $ROUP         ID         DESC
            # GROUP          1Group(1)_elar                                           +
            # $           TYPE     ID1  "THRU"     ID2
            # +           ELEM      21    THRU      36
            #print(data_dict)
            #op2.add_group(group_id, group_desc, data_dict)
                #i += 1
                #n += 4
            #assert ints[i] == -1, ints[i:]

            meta = data_dict['meta']
            nodes = data_dict['grid']
            elements = data_dict['element']
            properties = data_dict['property']
            op2.add_group(group_id, nodes, elements, properties)
            # self.log.warning(f'geom skipping GROUP in {self.table_name}')
            nentries += 1

        assert n == len(data), f'n={n} ndata={len(data)}'
        op2.increase_card_count('GROUP', nentries)
        assert nentries > 0, nentries
        return n

    def read_aero(self, data: bytes, n: int) -> int:
        """
        (3202, 32, 265)

        Word Name Type Description
        1 ACSID     I
        2 VELOCITY RS
        3 REFC     RS
        4 RHOREF   RS
        5 SYMXZ     I
        6 SYMXY     I

        """
        op2: OP2Geom = self.op2
        assert len(data) == 36, len(data)
        struct = Struct(op2._endian + b'i 3f 2i')
        out = struct.unpack(data[n:])
        acsid, velocity, cref, rho_ref, sym_xz, sym_xy = out
        aero = AERO(velocity, cref, rho_ref, acsid=acsid, sym_xz=sym_xz, sym_xy=sym_xy)
        op2._add_methods._add_aero_object(aero, allow_overwrites=True)
        #op2.add_aero(velocity, cref, rho_ref,
                     #acsid=acsid, sym_xz=sym_xz, sym_xy=sym_xy)
        n = 36
        return n

    def read_aeros(self, data: bytes, n: int) -> int:
        """
        AEROS(2202, 22, 340)

        AEROS   0       100     36.     360.    12960.
        data = (0, 100, 36.0, 360.0, 12960.0, 0, 0)

        """
        op2: OP2Geom = self.op2
        assert len(data) == 40 * self.factor, len(data)
        struct = Struct(mapfmt(op2._endian + b'2i 3f 2i', self.size))
        out = struct.unpack(data[n:])
        acsid, rcsid, cref, bref, sref, sym_xz, sym_xy = out
        op2.add_aeros(cref, bref, sref,
                      acsid=acsid, rcsid=rcsid,
                      sym_xz=sym_xz, sym_xy=sym_xy)
        n = 40 * self.factor
        return n

    def read_deform(self, data: bytes, n: int) -> int:
        """
        (104, 1, 81)
        NX 2019.2

        Word Name Type Description
        1 SID I Deformation set identification number
        2 EID I Element number
        3 D RS Deformation

        """
        op2: OP2Geom = self.op2
        ntotal = 12 * self.factor # 4*3
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        structi = Struct(mapfmt(op2._endian + b'2i f', self.size))
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            sid, eid, deformation = out
            deform = op2.add_deform(sid, eid, deformation)
            str(deform)
            n += ntotal
        return n

    def read_aeparm(self, data: bytes, n: int) -> int:
        """
        MSC 2020

        Word Name Type Description

        data = (601, 'PLOAD   ', 'LBS.    ')

        """
        op2: OP2Geom = self.op2
        ntotal = 20 * self.factor # 4*5
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        structi = Struct(mapfmt(op2._endian + b'i 8s 8s', self.size))
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            aeparm_id, label_bytes, units_bytes = structi.unpack(edata)
            label = label_bytes.decode(op2._encoding).rstrip()
            units = units_bytes.decode(op2._encoding).rstrip()
            op2.add_aeparm(aeparm_id, label, units, comment='')
            n += ntotal
        return n

    def read_uxvec(self, data: bytes, n: int) -> int:
        """
        MSC 2020

        Word Name Type Description
        1 ID           I Control vector identification number
        2 LABEL(2) CHAR4 Controller name
        4 UX          RS Magnitude of aerodynamic extra point dof
        Words 2 thru 4 repeat until (-1,-1) occurs

        data  = (5001, 'PLOAD   ', 1.0, 'INTERCPT', 0.0, -1, -1)
        """
        op2: OP2Geom = self.op2
        ntotal1 = 4 * self.factor # 4*1
        ntotal_end = 8 * self.factor # 4*2
        ntotal2 = 12 * self.factor # 4*3
        #ndatai = len(data) - n
        #ncards = ndatai // ntotal
        #assert ndatai % ntotal == 0
        struct1 = Struct(mapfmt(op2._endian + b'i', self.size))
        struct2 = Struct(mapfmt(op2._endian + b'8s f', self.size))
        struct_end = Struct(mapfmt(op2._endian + b'2i', self.size))
        #for unused_i in range(ncards):
            #edata = data[n:n + ntotal]
            #aeparm_id, label_bytes, units_bytes = structi.unpack(edata)
            #label = label_bytes.decode(op2._encoding).rstrip()
            #units = units_bytes.decode(op2._encoding).rstrip()
            #op2.add_uxvec()
            #n += ntotal

        while n < len(data):
            edata1 = data[n:n+ntotal1]
            idi, = struct1.unpack(edata1)
            n += ntotal1
            labels = []
            uxs = []

            edata_end = data[n:n+ntotal_end]
            while struct_end.unpack(edata_end) != (-1, -1):
                edata2 = data[n:n+ntotal2]
                label_bytes, ux = struct2.unpack(edata2)
                label = reshape_bytes_block_size(label_bytes, self.size)
                labels.append(label)
                uxs.append(ux)
                n += ntotal2
                edata_end = data[n:n+ntotal_end]
            n += ntotal_end
            uxvec = op2.add_uxvec(idi, labels, uxs)
            str(uxvec)
        return n # len(data)

    def read_caero1(self, data: bytes, n: int) -> int:
        """
        (3002, 30, 263)
        MSC 2018.2

        Word Name Type Description
        1 EID    I
        2 PID    I
        3 CP     I
        4 NSPAN  I
        5 NCHORD I
        6 LSPAN  I
        7 LCHORD I
        8 IGID   I
        9 X1    RS
        10 Y1   RS
        11 Z1   RS
        12 X12  RS
        13 X4   RS
        14 Y4   RS
        15 Z4   RS
        16 X43  RS

        CAERO1  100001  100001  0       10                      24      1
        99.2956521.45381-11.654442.85999101.8387122.6196-2.6930832.70996

        data = (100001, 100001, 0, 10, 0, 0, 24, 1,
                99.3, 21.45, -11.65, 42.86, 101.8387, 122.62, -2.69, 32.71)

        """
        op2: OP2Geom = self.op2
        add_methods = op2._add_methods
        ntotal = 64 * self.factor # 4*16
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        structi = Struct(mapfmt(op2._endian + b'8i 8f', self.size))
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            eid, pid, cp, nspan, nchord, lspan, lchord, igid, x1, y1, z1, x12, x4, y4, z4, x43 = out
            caero = CAERO1(
                eid, pid, igid,
                [x1, y1, z1], x12,
                [x4, y4, z4], x43, cp=cp,
                nspan=nspan, lspan=lspan, nchord=nchord, lchord=lchord,
                comment='')
            add_methods._add_caero_object(caero, allow_overwrites=True)
            #op2.add_caero1(eid, pid, igid,
                           #[x1, y1, z1], x12,
                           #[x4, y4, z4], x43,
                           #cp=cp,
                           #nspan=nspan, lspan=lspan,
                           #nchord=nchord, lchord=lchord)
            n += ntotal
        return n

    def read_caero2(self, data: bytes, n: int) -> int:
        """
        MSC 2018.2

        Word Name Type Description
        1 EID          I
        2 PID          I
        3 CP           I
        4 NSB          I
        5 NINT         I
        6 LSB          I
        7 LINT         I
        8 IGID         I
        9 X1          RS
        10 Y1         RS
        11 Z1         RS
        12 X12        RS
        13 UNDEF(4) none

        data = (54000, 4020, 0, 8, 8, 0, 0, 1, -5.0, 0, 0, 40.0, 0, 0, 0, 0),
        """
        op2: OP2Geom = self.op2
        ntotal = 64 * self.factor # 4*16
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        structi = Struct(mapfmt(op2._endian + b'8i 4f 4i', self.size))
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            eid, pid, cp, nsb, nint, lsb, lint, igroup, x1, y1, z1, x12, zero1, zero2, zero3, zero4 = out
            assert min(zero1, zero2, zero3, zero4) == max(zero1, zero2, zero3, zero4)
            p1 = [x1, y1, z1]
            caero2 = op2.add_caero2(eid, pid, igroup, p1, x12,
                                    cp=cp,
                                    nsb=nsb, nint=nint,
                                    lsb=lsb, lint=lint)
            str(caero2)
            n += ntotal
        return n

    def read_caero3(self, data: bytes, n: int) -> int:
        """
        Aerodynamic panel element configuration.

        Word Name Type Description
        1 EID    I Element identification number
        2 PID    I Property identification number of a PAERO3 entry
        3 CP     I Coordinate system for locating points 1 and 4
        4 LISTW  I Identification number of an AEFACT entry that lists
                   coordinate pairs for structural interpolation of the wing
        5 LISTC1 I Identification number of an AEFACT entry that lists
                   coordinate pairs for control surfaces
        6 LISTC2 I Identification number of an AEFACT entry that lists
                   coordinate pairs for control surfaces
        7 UNDEF(2) None
        9  X1   RS X-coordinate of point 1 in coordinate system CP
        10 Y1   RS Y-coordinate of point 1 in coordinate system CP
        11 Z1   RS Z-coordinate of point 1 in coordinate system CP
        12 X12  RS Edge chord length in aerodynamic coordinate system
        13 X4   RS X-coordinate of point 4 in coordinate system CP
        14 Y4   RS Y-coordinate of point 4 in coordinate system CP
        15 Z4   RS Z-coordinate of point 4 in coordinate system CP
        16 X43  RS Edge chord length in aerodynamic coordinate system
        """
        op2: OP2Geom = self.op2
        ntotal = 64 * self.factor # 4*16
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        structi = Struct(mapfmt(op2._endian + b'6i 2i 8f', self.size))
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            eid, pid, cp, list_w, list_c1, list_c2, zero1, zero2, x1, y1, z1, x12, x4, y4, z4, x43 = out
            #eid, pid, cp, nspan, lspan, zero1, zero2, zero3, x1, y1, z1, x12, x4, y4, z4, x43 = out
            assert min(zero1, zero2) == max(zero1, zero2)
            p1 = [x1, y1, z1]
            p4 = [x4, y4, z4]
            #4 LISTW  I Identification number of an AEFACT entry that lists
                       #coordinate pairs for structural interpolation of the wing
            #5 LISTC1 I Identification number of an AEFACT entry that lists
                       #coordinate pairs for control surfaces
            #6 LISTC2 I Identification number of an AEFACT entry that lists
                       #coordinate pairs for control surfaces
            if list_w == 0:
                list_w = None
            if list_c1 == 0:
                list_c1 = None
            if list_c2 == 0:
                list_c2 = None
            caero3 = op2.add_caero3(
                eid, pid, list_w, p1, x12, p4, x43,
                cp=cp, list_c1=list_c1, list_c2=list_c2, comment='')
            str(caero3)
            #print(caero3)
            n += ntotal
        return n

    def read_caero4(self, data: bytes, n: int) -> int:
        """
        Word Name Type Description
        1 EID   I Element identification number
        2 PID   I Property identification number of a PAERO4 entry
        3 CP    I Coordinate system for locating points 1 and 4
        4 NSPAN I Number of strips
        5 LSPAN I Identification number of an AEFACT entry
                  containing a list of division points for strips
        6 UNDEF(3) None
        9  X1  RS X-coordinate of point 1 in coordinate system CP
        10 Y1  RS Y-coordinate of point 1 in coordinate system CP
        11 Z1  RS Z-coordinate of point 1 in coordinate system CP
        12 X12 RS Edge chord length in aerodynamic coordinate system
        13 X4  RS X-coordinate of point 4 in coordinate system CP
        14 Y4  RS Y-coordinate of point 4 in coordinate system CP
        15 Z4  RS Z-coordinate of point 4 in coordinate system CP
        16 X43 RS Edge chord length in aerodynamic coordinate system
        """
        op2: OP2Geom = self.op2
        ntotal = 64 * self.factor # 4*16
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        structi = Struct(mapfmt(op2._endian + b'5i 3i 8f', self.size))
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            eid, pid, cp, nspan, lspan, zero1, zero2, zero3, x1, y1, z1, x12, x4, y4, z4, x43 = out
            assert min(zero1, zero2, zero3) == max(zero1, zero2, zero3)
            p1 = [x1, y1, z1]
            p4 = [x4, y4, z4]
            caero4 = op2.add_caero4(eid, pid, p1, x12, p4, x43,
                                    cp=cp, nspan=nspan, lspan=lspan, comment='')
            str(caero4)
            #print(caero4)
            n += ntotal
        return n

    def read_caero5(self, data: bytes, n: int) -> int:
        """
        MSC 2018.2

        Word Name Type Description
        1 EID        I
        2 PID        I
        3 CP         I
        4 NSPAN      I
        5 LSPAN      I
        6 NTHRY      I
        7 NTHICK     I
        8 UNDEF   none
        9 X1        RS
        10 Y1       RS
        11 Z1       RS
        12 X12      RS
        13 X4       RS
        14 Y4       RS
        15 Z4       RS
        16 X43      RS

        """
        op2: OP2Geom = self.op2
        ntotal = 64 # 4*16
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        structi = Struct(op2._endian + b'8i 8f')
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            eid, pid, cp, nspan, lspan, ntheory, nthick, undef, x1, y1, z1, x12, x4, y4, z4, x43 = out
            p1 = [x1, y1, z1]
            p4 = [x4, y4, z4]
            caero5 = op2.add_caero5(eid, pid,
                   p1, x12,
                   p4, x43,
                   cp=cp,
                   nspan=nspan, lspan=lspan,
                   ntheory=ntheory,
                   nthick=nthick)
            str(caero5)
            n += ntotal
        return n

    def read_paero1(self, data: bytes, n: int) -> int:
        r"""
        (3102, 31, 264)
        MSC 2018.2

        Word Name Type Description
        1 PID      I
        2 B1       I
        3 B2       I
        4 B3       I
        5 B4       I
        6 B5       I
        7 B6       I
        8 UNDEF none

        PAERO1  100001

        data = (100001, 100001, 0, 10, 0, 0, 24, 1,
                99.3, 21.45, -11.65, 42.86, 101.8387, 122.62, -2.69, 32.71)

        C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\adb144_2.op2
        PAERO1      1000   74000   74510   84610
        """
        op2: OP2Geom = self.op2
        ntotal = 32 * self.factor # 4 * 8
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        structi = Struct(mapfmt(op2._endian + b'8i', self.size))
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            pid, b1, b2, b3, b4, b5, b6, empty = out
            caero_body_ids = []
            for body in [b1, b2, b3, b4, b5, b6, empty]:
                if body != 0:
                    caero_body_ids.append(body)
            paero1 = PAERO1(pid, caero_body_ids=caero_body_ids, comment='')
            op2._add_methods._add_paero_object(paero1, allow_overwrites=True)
            #paero1 = op2.add_paero1(pid, caero_body_ids=caero_body_ids)
            str(paero1)
            #if caero_body_ids:
                #self.log.warning(str(paero1))
            n += ntotal
        return n

    def read_paero2(self, data: bytes, n: int) -> int:
        """
        MSC 2018.2

        Word Name Type Description
        1 PID        I
        2 ORIENT CHAR4
        3 UNDEF   none (orient carryover)
        4 WIDTH     RS
        5 AR        RS
        6 LRSB       I
        7 LRIB       I
        8 LTH1       I
        9 LTH2       I
        10 THI1      I
        11 THN1      I
        12 THI2      I
        13 THN2      I
        14 THI3      I
        15 THN3      I

        PAERO2  100001

        data = (100001, 100001, 0, 10, 0, 0, 24, 1,
                99.3, 21.45, -11.65, 42.86, 101.8387, 122.62, -2.69, 32.71)

        """
        op2: OP2Geom = self.op2
        ntotal = 60 * self.factor # 4 * 15
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        if self.size == 4:
            structi = Struct(op2._endian + b'i4si 2f 10i')
        else:
            structi = Struct(op2._endian + b'q8sq 2d 10q')
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            (pid, orient_bytes, undef, width, ar, lrsb, lrib, lth1, lth2,
             thi1, thn1,
             thi2, thn2,
             thi3, thn3) = out
            lth = [lth1, lth2]
            thi = [thi1, thi2, thi3]
            thn = [thn1, thn2, thn3]
            orient = reshape_bytes_block_size(orient_bytes, self.size)
            assert isinstance(orient, str), orient
            paero2 = op2.add_paero2(pid, orient, width, ar,
                                    thi, thn,
                                    lrsb=lrsb,
                                    lrib=lrib,
                                    lth=lth)
            n += ntotal
            str(paero2)
        return n

    def read_paero3(self, data: bytes, n: int) -> int:
        """
        NX 2019.2
        Record – PAERO3(4701,47,171)
        Aerodynamic panel property.
        1 PID  I Property identification number
        2 NBOX I Number of Mach boxes in the flow direction
        3 FLAG I
        FLAG = 0
        4 UNDEF None
        5 X5 RS X-coordinate of point 5 in the aerodynamic coordinate system defining the cranks and control surface geometry
        6 Y5 RS Y-coordinate of point 5 in the aerodynamic coordinate system defining the cranks and control surface geometry
        7 X6 RS X-coordinate of point 6 in the aerodynamic coordinate system defining the cranks and control surface geometry
        8 Y6 RS Y-coordinate of point 6 in the aerodynamic coordinate system defining the cranks and control surface geometry
        FLAG = 1
        4 UNDEF None
        5 X5 RS X-coordinate of point 5 in the aerodynamic coordinate system defining the cranks and control surface geometry
        6 Y5 RS Y-coordinate of point 5 in the aerodynamic coordinate system defining the cranks and control surface geometry
        7 X6 RS X-coordinate of point 6 in the aerodynamic coordinate system defining the cranks and control surface geometry
        8 Y6 RS Y-coordinate of point 6 in the aerodynamic coordinate system defining the cranks and control surface geometry
        9 X7 RS X-coordinate of point 7 in the aerodynamic coordinate system defining the cranks and control surface geometry
        10 Y7 RS Y-coordinate of point 7 in the aerodynamic coordinate system defining the cranks and control surface geometry
        11 X8 RS X-coordinate of point 8 in the aerodynamic coordinate system defining the cranks and control surface geometry
        12 Y8 RS Y-coordinate of point 8 in the aerodynamic coordinate system defining the cranks and control surface geometry
        13 X9 RS X-coordinate of point 9 in the aerodynamic coordinate system defining the cranks and control surface geometry
        14 Y9 RS Y-coordinate of point 9 in the aerodynamic coordinate system defining the cranks and control surface geometry
        15 X10 RS X-coordinate of point 10 in the aerodynamic coordinate system defining the cranks and control surface geometry
        16 Y10 RS Y-coordinate of point 10 in the aerodynamic coordinate system defining the cranks and control surface geometry
        FLAG = 2
        4 UNDEF None
        5 X5 RS X-coordinate of point 5 in the aerodynamic coordinate system defining the cranks and control surface geometry
        6 Y5 RS Y-coordinate of point 5 in the aerodynamic coordinate system defining the cranks and control surface geometry
        7 X6 RS X-coordinate of point 6 in the aerodynamic coordinate system defining the cranks and control surface geometry
        8 Y6 RS Y-coordinate of point 6 in the aerodynamic coordinate system defining the cranks and control surface geometry
        9 X7 RS X-coordinate of point 7 in the aerodynamic coordinate system defining the cranks and control surface geometry
        10 Y7 RS Y-coordinate of point 7 in the aerodynamic coordinate system defining the cranks and control surface geometry
        11 X8 RS X-coordinate of point 8 in the aerodynamic coordinate system defining the cranks and control surface geometry
        12 Y8 RS Y-coordinate of point 8 in the aerodynamic coordinate system defining the cranks and control surface geometry
        13 X9 RS X-coordinate of point 9 in the aerodynamic coordinate system defining the cranks and control surface geometry
        14 Y9 RS Y-coordinate of point 9 in the aerodynamic coordinate system defining the cranks and control surface geometry
        15 X10 RS X-coordinate of point 10 in the aerodynamic coordinate system defining the cranks and control surface geometry
        16 Y10 RS Y-coordinate of point 10 in the aerodynamic coordinate system defining the cranks and control surface geometry
        17 X11 RS X-coordinate of point 11 in the aerodynamic coordinate system defining the cranks and control surface geometry
        18 Y11 RS Y-coordinate of point 11 in the aerodynamic coordinate system defining the cranks and control surface geometry
        19 X12 RS X-coordinate of point 12 in the aerodynamic coordinate system defining the cranks and control surface geometry
        20 Y12 RS Y-coordinate of point 12 in the aerodynamic coordinate system defining the cranks and control surface geometry

        data = ()
        ints    = (4701, 47, 171,
                   1, 11, 0, 4, 0, 0, 0, 0)
        floats  = (4701, 47, 171,
                   1, 11, 0.0, 4, 0.0, 0.0, 0.0, 0.0)

        """
        op2: OP2Geom = self.op2
        #op2.show_data(data, types='ifs', endian=None, force=False)
        #asdf
        #ntotal = 60 * self.factor # 4 * 15
        #ndatai = len(data) - n
        #ncards = ndatai // ntotal
        #assert ndatai % ntotal == 0
        if self.size == 4:
            struct_init = Struct(op2._endian + b'3i')
            struct0 = Struct(op2._endian + b'i 4f')
            struct1 = Struct(op2._endian + b'i 12f')
        else:
            struct_init = Struct(op2._endian + b'3q')
            struct0 = Struct(op2._endian + b'q 4d')
            struct1 = Struct(op2._endian + b'q 12d')
        n_init = 3 * self.size
        n0 = 5 * self.size
        n1 = 13 * self.size

        #ncards = 0
        while n < len(data):
            #op2.show_data(data[n:])
            data_init = data[n:n+n_init]
            pid, nbox, flag = struct_init.unpack(data_init)
            n += n_init
            #print(pid, nbox, flag)
            #if pid == 0 and nbox == 0 and flag == 0:
                #continue
            ncontrol_surfaces = 0
            x = []
            y = []
            if flag == 0:
                data0 = data[n:n+n0]
                assert len(data0) == n0, (len(data0), n0)
                #4 UNDEF None
                #5 X5 RS X-coordinate of point 5 in the aerodynamic coordinate system defining the cranks and control surface geometry
                #6 Y5 RS Y-coordinate of point 5 in the aerodynamic coordinate system defining the cranks and control surface geometry
                #7 X6 RS X-coordinate of point 6 in the aerodynamic coordinate system defining the cranks and control surface geometry
                #8 Y6 RS Y-coordinate of point 6 in the aerodynamic coordinate system defining the cranks and control surface geometry
                dummy, x5, y5, x6, y6 = struct0.unpack(data0)
                x.extend([x5, x6])
                y.extend([y5, y6])
                #print('  ', dummy, x5, y5, x6, y6)
                paero3 = op2.add_paero3(pid, nbox, ncontrol_surfaces, x, y)
                str(paero3)
                n += n0
            elif flag == 1:
                data1 = data[n:n+n1]
                assert len(data1) == n1, (len(data1), n1)
                #4 UNDEF None
                #5 X5 RS X-coordinate of point 5 in the aerodynamic coordinate system defining the cranks and control surface geometry
                #6 Y5 RS Y-coordinate of point 5 in the aerodynamic coordinate system defining the cranks and control surface geometry
                #7 X6 RS X-coordinate of point 6 in the aerodynamic coordinate system defining the cranks and control surface geometry
                #8 Y6 RS Y-coordinate of point 6 in the aerodynamic coordinate system defining the cranks and control surface geometry
                #9 X7 RS X-coordinate of point 7 in the aerodynamic coordinate system defining the cranks and control surface geometry
                #10 Y7 RS Y-coordinate of point 7 in the aerodynamic coordinate system defining the cranks and control surface geometry
                #11 X8 RS X-coordinate of point 8 in the aerodynamic coordinate system defining the cranks and control surface geometry
                #12 Y8 RS Y-coordinate of point 8 in the aerodynamic coordinate system defining the cranks and control surface geometry
                #13 X9 RS X-coordinate of point 9 in the aerodynamic coordinate system defining the cranks and control surface geometry
                #14 Y9 RS Y-coordinate of point 9 in the aerodynamic coordinate system defining the cranks and control surface geometry
                #15 X10 RS X-coordinate of point 10 in the aerodynamic coordinate system defining the cranks and control surface geometry
                #16 Y10 RS Y-coordinate of point 10 in the aerodynamic coordinate system defining the cranks and control surface geometry
                dummy, x5, y5, x6, y6, x7, y7, x8, y8, x9, y9, x10, y10 = struct1.unpack(data1)
                x.extend([x5, x6, x7, x8, x9, x10])
                y.extend([y5, y6, y7, y8, y9, y10])
                paero3 = op2.add_paero3(pid, nbox, ncontrol_surfaces, x, y)
                str(paero3)
                n += n1
            else:
                raise NotImplementedError(flag)
            #ncards += 1
        return n

    def read_paero4(self, data: bytes, n: int) -> int:
        """
        67 - PAERO4(4801,48,172)
        1 PID   I
        2 CLA   I
        3 LCLA  I
        4 CIRC  I
        5 LCIRC I
        6 DOCI   RS
        7 CAOCI  RS
        8 GAPOCI RS
        Words 6 through 8 repeat until End of Record
        """
        op2: OP2Geom = self.op2
        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        floats = np.frombuffer(data[n:], op2.fdtype8).copy()
        istart, iend = get_minus1_start_end(ints)

        for (i0, i1) in zip(istart, iend):
            pid, cla, lcla, circ, lcirc = ints[i0:i0+5]
            doci_caoci_gapoci = floats[i0+5:i1]
            nrow = len(doci_caoci_gapoci) // 3
            doci_caoci_gapoci = doci_caoci_gapoci.reshape(nrow, 3)
            assert ints[i1] == -1, ints[i1]
            docs = doci_caoci_gapoci[:, 0].tolist()
            caocs = doci_caoci_gapoci[:, 1].tolist()
            gapocs = doci_caoci_gapoci[:, 2].tolist()
            paero4 = op2.add_paero4(
                pid, docs, caocs, gapocs,
                cla=cla, lcla=lcla, circ=circ, lcirc=lcirc)
            str(paero4)
        return len(data)

    def read_paero5(self, data: bytes, n: int) -> int:
        """
        MSC 2018.2

        Word Name Type Description
        1 PID    I
        2 NALPHA I
        3 LALPHA I
        4 NXIS   I
        5 LXIS   I
        6 NTAUS  I
        7 LTAUS  I
        8 CAOCI RS
        Word 8 repeats until End of Record

        """
        op2: OP2Geom = self.op2
        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        floats = np.frombuffer(data[n:], op2.fdtype8).copy()
        istart, iend = get_minus1_start_end(ints)

        for (i0, i1) in zip(istart, iend):
            pid, nalpha, lalpha, nxis, lxis, ntaus, ltaus = ints[i0:i0+7]
            caoci = floats[i0+7:i1]
            assert ints[i1] == -1, ints[i1]
            paero5 = op2.add_paero5(
                pid, caoci,
                nalpha=nalpha, lalpha=lalpha,
                nxis=nxis, lxis=lxis,
                ntaus=ntaus, ltaus=ltaus)
            str(paero5)
        return len(data)

    def read_panel(self, data: bytes, n: int) -> int:
        """
        MSC 2018.2

        Word Name Type Description
        1 NAME(2) CHAR4
        3 SETID I
        Words 1 through 3 repeat until End of Record

        ('PANEL1', 1, -1)
        """
        op2: OP2Geom = self.op2
        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        istart, iend = get_minus1_start_end(ints)

        for (i0, i1) in zip(istart, iend):
            assert ints[i1] == -1, ints[i1]
            names = []
            set_ids = []
            while i0 < i1:
                name_bytes = data[n:n+8]
                name = reshape_bytes_block_size(name_bytes, self.size)
                set_id = ints[i0+2]
                names.append(name)
                set_ids.append(set_id)
                n += 12
                i0 += 3
            panel = op2.add_panel(names, set_ids)
            str(panel)
        return len(data)

    def read_acmodl(self, data: bytes, n: int) -> int:
        """Reads the ACMODL card"""
        op2: OP2Geom = self.op2
        card_name = 'ACMODL'
        card_obj = ACMODL
        methods = {
            72 : self._read_acmodl_nx_72,
            64 : self._read_acmodl_msc_64,
        }
        try:
            n = op2.reader_geom2._read_double_card(
                card_name, card_obj,
                op2._add_methods._add_acmodl_object,
                methods, data, n)
        except DoubleCardError:
            raise
        return n

        #n = self._read_dual_card(data, n, self._read_acmodl_nx, self._read_acmodl_msc,
                                 #'ACMODL', op2._add_methods._add_acmodl_object)
        ##return self._read_acmodl_msc(data, n)
        #return n


    def _read_acmodl_nx_72(self, card, data: bytes, n: int) -> tuple[int, list[ACMODL]]:
        """
        NX 2019.2 - 72 bytes
        Word Name Type Description
        1 INTER(2)      CHAR4 IDENT or DIFF method specification
        3 INFO(2)       CHAR4 Allowable values are ALL, ELEMENTS, PID, SET3, and NONE
        5 FSET              I Fluid set ID
        6 SSET              I Structure set ID
        7 NORML            RS Outward normal search distance to detect fluid-structure interface
        8 METHOD(2)     CHAR4 Interface calculation method
        10 OVPLANG         RS Angular tolerance in degrees used to decide
                              whether a fluid free face and a structural face
                              are overlapping
        11 SRCHUNIT(2)  CHAR4 Search unit
        13 INTOL           RS Inward normal search distance to detect
                              fluid-structure interface
        14 AREAOPT          I Area option
        15 SKNEPS          RS SKNEPS option. Only used when AREAOPT = 0
        16 INTORD           I Integration order
        17 CTYPE(2)     CHAR4 Coupling type (STRONG, WEAK, WEAKINT, or WEAKEXT)

        strings = (b'IDENT   NONE   ,                     \x00\x00pBREL     \xcd\xccL>\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00STRONG  ',)
        ints    = ('IDENT', 'NONE', 0,     0, 1.0e-4, AS, 60.0, 541869394, 538976288,                           0.2, 0, 0.0, 0, STRONG)
        floats  = ('IDENT', 'NONE', 0.0, 0.0, 1.0e-4, AS, 60.0, 1.7302408291410174e-19, 1.3563156426940112e-19, 0.2, 0, 0.0, 0, STRONG)

        MSC 2018.2
        | ACMODL | INTER | INFOR   |   FSET   | SSET | NORMAL | METHOD |  SKNEPS | DSKNEPS  |
        |        | INTOL | ALLSSET | SRCHUNIT |      |        |        |         |          |
        # NX 2019.2
        | ACMODL |       |  INFOR  |   FSET   | SSET | NORMAL |        | OVLPANG | SRCHUNIT |
        |        | INTOL | AREAOP  |          |      |  CTYPE |        |         |          |
        ACMODL,IDENT,,,,1.0-4

        """
        op2: OP2Geom = self.op2
        ntotal = 72 *  self.factor # 4 * 8
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0, ndatai % ntotal
        #structi = Struct(op2._endian + b'4i f 8s 8s 3i f') # msc
        if self.size == 4:
            structi = Struct(op2._endian + b'8s 8s 2i f 8s f 8s f ifi 8s')
        else:
            structi = Struct(op2._endian + b'16s 16s 2q d 16s d 16s d qdq 16s')

        acmodls = []
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            #Type of structure-fluid interface. (Character = IDENT or DIFF;
            #inter: good
            #olvpang: good
            #search_unit_bytes: good
            #ctype: good
            #area_op: good
            #sk_neps/olvpang;  Default=60.0
            inter_bytes, infor_bytes, fset, sset, normal, method_bytes, olvpang, search_unit_bytes, intol, area_op, sk_neps, intord, ctype_bytes = out
            inter = reshape_bytes_block_size(inter_bytes, self.size)
            infor = reshape_bytes_block_size(infor_bytes, self.size)
            method = reshape_bytes_block_size(method_bytes, self.size)
            search_unit = reshape_bytes_block_size(search_unit_bytes, self.size)
            ctype = reshape_bytes_block_size(ctype_bytes, self.size)

            assert inter in ['IDEN', 'IDENT', 'DIFF'], inter
            assert ctype in ['STRONG', 'WEAK', 'WEAKINT', 'WEAKEXT'], ctype
            assert method in ['AS'], method
            assert area_op in [0, 1], area_op
            #If CTYPE = STRONG
            #If CTYPE = WEAK
            #print(f'inter={inter!r} infor={infor!r} fset={fset} sset={sset} normal={normal:g} method={method} olvpang={olvpang} search_unit={search_unit!r}\n'
                  #f'intol={intol:g} area_op={area_op} sk_neps={sk_neps} intord={intord} ctype={ctype!r}')

            #If SRCHUNIT = ABS, then the model units are absolute.
            #If SRCHUNIT = REL, then the relative model units are based on element size.
            assert search_unit in ['ABS', 'REL'], search_unit
            #INTOL Inward normal sea
            # set2 = op2.add_set2(sid, macro, sp1, sp2, ch1, ch2, zmax, zmin)
            acmodl = ACMODL(infor, fset, sset,
                            normal=normal, olvpang=olvpang,
                            search_unit=search_unit, intol=intol,
                            area_op=area_op,
                            ctype=ctype,
                            method='BW',
                            sk_neps=sk_neps,
                            #dsk_neps=0.75,
                            #all_set='NO',
                            inter=inter,
                            nastran_version='nx')
            #print(acmodl)
            str(acmodl)
            acmodls.append(acmodl)
            n += ntotal
        return n, acmodls

    def _read_acmodl_msc_64(self, card, data: bytes, n: int) ->  tuple[int, list[ACMODL]]:
        """
        MSC 2018.2 - 64 bytes
        Word Name Type Description


        1 INTER(2)     CHAR4 Type of structure-fluid interface:
                             "IDENT","DIFF"-def"DIFF"
        3 INFOR(2)     CHAR4 If INTER="DIFF", defines the type of list:
                             "GRIDS","ELEMENTS","NONE"-def"NONE"
        5 FSET             I SET1 ID for fluid elements or grids ID list, FSET>0 or blank
        6 SSET             I SET1 ID for struc elements or grids ID list, SSET>0 or blank
        7 NORMAL          RS Fluid normal tolerance - def 1.0
        8 METHOD(2)    CHAR4 Method -def" "
        10 SKNEPS         RS Fluid skin growth tolereance - def 0.75
        11 DSKNEPS        RS Fluid secondary skin growth tolerance - def 0.75
        12 INTOL          RS Tolerance of inward normal - def 0.5
        13 ALLSSET(2)  CHAR4 coupled all structure given by SSET if "YES" - def"NO"
        15 SRCHUNIT(2) CHAR4 Search Units:"ABS","REL"-def"REL"

        ndata = 64: # MSC
          strings = (b'IDENT   NONE    \x00\x00\x00\x00\x00\x00\x00\x00\x17\xb7\xd18        \x00\x00\x00?\x00\x00@?\x00\x00\x00?NO      REL     ',)
          ints    = (IDENT, NONE, 0,   0,   953267991, 538976288, 538976288,                                       0.5, 0.75, 0.5, NO, REL)
          floats  = (IDENT, NONE, 0.0, 0.0, 9.999999747378752e-05, 1.3563156426940112e-19, 1.3563156426940112e-19, 0.5, 0.75, 0.5, NO, REL)

        MSC 2018.2
        | ACMODL | INTER | INFOR   |   FSET   | SSET | NORMAL | METHOD |  SKNEPS | DSKNEPS  |
        |        | INTOL | ALLSSET | SRCHUNIT |      |        |        |         |          |

        """
        op2: OP2Geom = self.op2
        ntotal = 64 *  self.factor # 4 * 8
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0, ndatai % ntotal
        #structi = Struct(op2._endian + b'4i f 8s 8s 3i f') # msc
        if self.size == 4:
            structi = Struct(op2._endian + b'8s 8s 2if 8s 3f 8s 8s')
        else:
            raise NotImplementedError(('ACMODL-MSC', self.size))
            structi = Struct(op2._endian + b'16s 16s 2q d 16s d 16s d qdq 16s')

        acmodls = []
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            #Type of structure-fluid interface. (Character = IDENT or DIFF;
            #inter: good
            #olvpang: good
            #search_unit_bytes: good
            #ctype: good
            #area_op: good
            #sk_neps/olvpang;  Default=60.0
            inter_bytes, infor_bytes, fset, sset, normal, method_bytes, sk_neps, dsk_neps, intol, all_sset, search_unit_bytes = out
            area_op= None
            olvpang = None
            #inter_bytes, infor_bytes, fset, sset, normal, method_bytes, olvpang, search_unit_bytes, intol, area_op, sk_neps, intord, ctype_bytes = out
            inter = reshape_bytes_block_size(inter_bytes, self.size)
            infor = reshape_bytes_block_size(infor_bytes, self.size)
            method = reshape_bytes_block_size(method_bytes, self.size)
            search_unit = reshape_bytes_block_size(search_unit_bytes, self.size)
            #ctype = reshape_bytes_block_size(ctype_bytes, self.size)

            assert inter in {'IDENT', 'DIFF'}, inter
            #assert ctype in ['STRONG', 'WEAK', 'WEAKINT', 'WEAKEXT'], ctype
            assert method in {'CP', 'BW', ''}, method
            #assert area_op in [0, 1], area_op
            #If CTYPE = STRONG
            #If CTYPE = WEAK
            #print(f'inter={inter!r} infor={infor!r} fset={fset} sset={sset} normal={normal:g} method={method} olvpang={olvpang} search_unit={search_unit!r}\n'
                  #f'intol={intol:g} area_op={area_op} sk_neps={sk_neps} intord={intord} ctype={ctype!r}')

            #If SRCHUNIT = ABS, then the model units are absolute.
            #If SRCHUNIT = REL, then the relative model units are based on element size.
            assert search_unit in ['ABS', 'REL'], search_unit
            #INTOL Inward normal sea
            # set2 = op2.add_set2(sid, macro, sp1, sp2, ch1, ch2, zmax, zmin)
            acmodl = ACMODL(infor, fset, sset,
                            normal=normal, olvpang=olvpang,
                            search_unit=search_unit, intol=intol,
                            area_op=area_op,
                            #ctype=ctype,
                            method='BW',
                            sk_neps=sk_neps,
                            #dsk_neps=0.75,
                            #all_set='NO',
                            inter=inter,
                            nastran_version='msc')
            #print(acmodl)
            str(acmodl)
            acmodls.append(acmodl)
            n += ntotal
        return n, acmodls

    def read_aelist(self, data: bytes, n: int) -> int:
        """
        MSC 2018.2

        Word Name Type Description
        1 SID I
        2 E I
        Word 2 repeats until End of Record

        """
        op2: OP2Geom = self.op2
        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        floats = np.frombuffer(data[n:], op2.fdtype8).copy()
        istart, iend = get_minus1_start_end(ints)

        for (i0, i1) in zip(istart, iend):
            sid = ints[i0]
            elements = floats[i0+1:i1].tolist()
            assert ints[i1] == -1, ints[i1]
            op2.add_aelist(sid, elements)
            #n += ntotal
        return len(data)

    def read_set1(self, data: bytes, n: int) -> int:
        """
        SET1: (3502, 35, 268)
        MSC 2018.2

        Word Name Type Description
        1 SID I
        2 G1  I Grid ID or -2 when SKIN is specified
        Word 2 repeats until End of Record

        """
        op2: OP2Geom = self.op2
        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        istart, iend = get_minus1_start_end(ints)

        for (i0, i1) in zip(istart, iend):
            sid = ints[i0]
            elements = ints[i0+1:i1].tolist()
            assert -2 not in elements, elements
            assert ints[i1] == -1, ints[i1]
            op2.add_set1(sid, elements, is_skin=False)
        return len(data)

    def read_set2(self, data: bytes, n: int) -> int:
        """
        SET2
        MSC 2018.2
        NX 2019.2

        Word Name Type Description
        1 SID I
        2 G1  I Grid ID or -2 when SKIN is specified
        Word 2 repeats until End of Record

        Record 71 - SET2(3602,36,269)
        Word Name Type Description
        1 SID   I  Unique identification number
        2 MACRO I  Identification number of an aerodynamic macro element
        3 SP1   RS Lower span division point defining the prism containing the set
        4 SP2   RS Higher span division point defining the prism containing the set
        5 CH1   RS Lower chord division point defining the prism containing the set
        6 CH2   RS Higher chord division point defining the prism containing the set
        7 ZMAX  RS Z-coordinate of the top of the prism containing the set
        8 ZMIN  RS Z-coordinate of the bottom of the prism containing the set.

        data = (3602, 36, 269,
                200, 101, -0.10, 1.10, -0.10, 1.0, 1.0, -0.10)
        """
        op2: OP2Geom = self.op2
        #op2.show_data(data, types='if', endian=None, force=False)
        ntotal = 32 * self.factor # 4 * 8
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0, (ndatai, ntotal)
        #structi = Struct(op2._endian + b'4i f 8s 8s 3i f') # msc
        structi = Struct(op2._endian + b'2i 6f')
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            sid, macro, sp1, sp2, ch1, ch2, zmax, zmin = out
            set2 = op2.add_set2(sid, macro, sp1, sp2, ch1, ch2, zmax, zmin)
            str(set2)
            n += ntotal
        op2.to_nx(' because SET2 was found')
        return n

    def read_set3(self, data: bytes, n: int) -> int:
        """
        MSC 2018.2
        Word Name Type Description
        1 SID I
        2 DES I Set description:
            1=ELEM
            2=GRID
            3=POINT
            4=PROP
            5=RBEin
            6=RBEex
        3 ID1 I IDs of Grids, Elements, Points or Properties.
        4 "ID1 THRU ID2" format will be EXPANDED into explicit IDs.
        Words 3 through 4 repeat until End of Record

        data = (1, 1, 190, ..., 238, -1,
                2, 1, 71, ..., 189, -1,
                4, 1, 309, ..., ..., 378, -1)
        """
        op2: OP2Geom = self.op2
        # this is setup for NX
        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        istart, iend = get_minus1_start_end(ints)

        for (i0, i1) in zip(istart, iend):
            sid = ints[i0]
            desc_int = ints[i0+1]
            elements = ints[i0+2:i1].tolist()
            if desc_int == 1:
                desc = 'ELEM'
            elif desc_int == 2:
                desc = 'GRID'
            elif desc_int == 3:
                desc = 'POINT'
            elif desc_int == 4:
                desc = 'PROP'
            elif desc_int == 5:
                desc = 'RBEin'
            elif desc_int == 6:
                desc = 'RBEex'
            else:
                raise NotImplementedError(desc_int)

            assert min(elements) > 0, elements
            assert ints[i1] == -1, ints[i1]
            set3 = op2.add_set3(sid, desc, elements)
            str(set3)
        return len(data)

    def read_aelink(self, data: bytes, n: int) -> int:
        """
        MSC 2018.2

        Word Name Type Description
        1 ID I
        2 LABLD(2) CHAR4
        4 LABLI(2) CHAR4
        6 C1 RS
        Words 4 through 6 repeat until (-1,-1,-1) occurs

        """
        op2: OP2Geom = self.op2
        struct1 = Struct(op2._endian + b'i8s')
        struct2 =Struct(op2._endian + b'8sf')
        struct_end = Struct(op2._endian + b'3i')
        ntotal = 12
        while n < len(data):
            edata = data[n:n+ntotal]
            aelink_id, label_bytes = struct1.unpack(edata)
            #print(f'  {aelink_id}, {label_bytes}')
            if aelink_id == 0:
                aelink_id = 'ALWAYS'
            #assert aelink_id > 0, aelink_id
            label = reshape_bytes_block_size(label_bytes, self.size)
            n += ntotal
            linking_coefficents = []
            independent_labels = []

            edata = data[n:n+ntotal]
            while struct_end.unpack(edata) != (-1, -1, -1):
                ind_label_bytes, coeff = struct2.unpack(edata)
                ind_label = reshape_bytes_block_size(ind_label_bytes, self.size)
                independent_labels.append(ind_label)
                assert isinstance(ind_label, str), ind_label
                linking_coefficents.append(coeff)
                #print(f'  {ind_label_bytes}, {coeff}')
                n += ntotal
                edata = data[n:n+ntotal]
            n += ntotal
            #print('  (-1, -1, -1)\n')
            assert isinstance(label, str), label
            aelink = op2.add_aelink(aelink_id, label,
                                     independent_labels, linking_coefficents)
            str(aelink)
        return len(data)

    def read_aecomp(self, data: bytes, n: int) -> int:
        """
        MSC 2018.2

        Word Name Type Description
        1 NAME(2)     CHAR4
        3 LISTTYPE(2) CHAR4
        5 LISTID          I
        Word 5 repeats until End of Record

        """
        op2: OP2Geom = self.op2
        ints = np.frombuffer(data[n:], op2.idtype).copy()
        #floats = np.frombuffer(data[n:], op2.fdtype).copy()
        istart, iend = get_minus1_start_end(ints)

        for (i0, i1) in zip(istart, iend):
            name_bytes = data[n+i0*4:n+i0*4+8]
            list_type_bytes = data[n+i0*4+8:n+i0*4+16]
            lists = ints[i0+4:i1].tolist()
            assert ints[i1] == -1, ints[i1]
            name = name_bytes.rstrip().decode('ascii')
            list_type = list_type_bytes.rstrip().decode('ascii')
            #print(name, list_type, lists)
            aecomp = op2.add_aecomp(name, list_type, lists)
            str(aecomp)
        """
        Word Name Type Description
        1 NAME(2) CHAR4
        3 LABEL(2) CHAR4
        Words 3 through 4 repeat until (-1,-1) occurs
        """
        return len(data)

    def read_aecompl(self, data: bytes, n: int) -> int:
        """
        MSC 2018.2

        Word Name Type Description
        1 NAME(2)  CHAR4
        3 LABEL(2) CHAR4
        Words 3 through 4 repeat until (-1,-1) occurs
        """
        op2: OP2Geom = self.op2
        struct1 = Struct(op2._endian + b'8s')
        struct2 = Struct(op2._endian + b'8s')
        struct_end = Struct(op2._endian + b'2i')
        ntotal = 8
        while n < len(data):
            edata = data[n:n+ntotal]
            name_bytes, = struct1.unpack(edata)
            name = name_bytes.decode('latin1').rstrip()
            n += ntotal
            labels = []

            edata = data[n:n+ntotal]
            while struct_end.unpack(edata) != (-1, -1):
                label_bytes, = struct2.unpack(edata)
                label = reshape_bytes_block_size(label_bytes, self.size)
                labels.append(label)
                n += ntotal
                edata = data[n:n+ntotal]
            n += ntotal
            aecompl = op2.add_aecompl(name, labels)
            str(aecompl)
        return len(data)

    def read_spline1(self, data: bytes, n: int) -> int:
        """reads the SPLINE1 card"""
        n = self._read_spline1_nx(data, n)
        return n

    def _read_spline1_nx(self, data: bytes, n: int) -> int:
        """
        Word Name Type Description
        1 EID   I
        2 CAERO I
        3 BOX1  I
        4 BOX2  I
        5 SETG  I
        6 DZ    RS
        7 METHOD(2) CHAR4 Method: IPS|TPS|FPS
        9 USAGE(2) CHAR4 Usage flag: FORCE|DISP|BOTH
        11 NELEM I Number of elements for FPS on x-axis
        12 MELEM I Number of elements for FPS on y-axis

        """
        op2: OP2Geom = self.op2
        ntotal = 48 # 4 * 12
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        #structi = Struct(op2._endian + b'4i f 8s 8s 3i f') # msc
        structi = Struct(op2._endian + b'5if 8s 8s 2i')
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            eid, caero, box1, box2, setg, dz, method_bytes, usage_bytes, nelements, melements = out
            method = method_bytes.rstrip().decode('ascii')
            usage = usage_bytes.rstrip().decode('ascii')
            #if nelements == 0:
                #nelements = None
            #if melements == 0:
                #melements = None

            spline1 = SPLINE1(eid, caero, box1, box2, setg,
                              dz=dz, method=method,
                            usage=usage, nelements=nelements,
                            melements=melements)
            op2._add_methods._add_spline_object(spline1, allow_overwrites=True)
            #spline1 = op2.add_spline1(eid, caero, box1, box2, setg,
                                       #dz=dz, method=method,
                                       #usage=usage, nelements=nelements,
                                       #melements=melements)
            str(spline1)
            n += ntotal
        #op2.to_nx()
        return n

    def read_spline2(self, data: bytes, n: int) -> int:
        """
        Reads the SPLINE2 card

        Word Name Type Description
        1 EID   I
        2 CAERO I
        3 ID1   I
        4 ID2   I
        5 SETG  I
        6 DZ    RS
        7 DTOR  RS
        8 CID   I
        9 DTHX  RS
        10 DTHY RS
        11 USAGE(2) CHAR4 Usage flag: FORCE|DISP|BOTH

        """
        op2: OP2Geom = self.op2
        ntotal = 48 * self.factor # 4 * 12
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        #assert ndatai % ntotal == 0
        if self.size == 4:
            structi = Struct(op2._endian + b'5i 2f i 2f 8s')
        else:
            structi = Struct(op2._endian + b'5q 2d q 2d 16s')

        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            eid, caero, id1, id2, setg, dz, dtor, cid, dthx, dthy, usage_bytes = out
            usage = usage_bytes.rstrip().decode('latin1')
            spline2 = SPLINE2(
                eid, caero,
                id1, id2, setg,
                dz=dz, dtor=dtor, cid=cid,
                dthx=dthx, dthy=dthy,
                usage=usage)
            op2._add_methods._add_spline_object(spline2, allow_overwrites=True)

            #spline2 = op2.add_spline2(
                #eid, caero,
                #id1, id2, setg,
                #dz=dz, dtor=dtor, cid=cid,
                #dthx=dthx, dthy=dthy,
                #usage=usage)
            str(spline2)
            n += ntotal
        #n = self._read_spline2_nx(data, n)
        return n

    def read_spline3(self, data: bytes, n: int) -> int:
        """reads the SPLINE3 card"""
        spline3
        #n = self._read_spline2_nx(data, n)
        return n

    def read_spline4(self, data: bytes, n: int) -> int:
        """reads the SPLINE4 card"""
        op2: OP2Geom = self.op2
        card_name = 'SPLINE4'
        card_obj = SPLINE4
        methods = {
            44 : self._read_spline4_nx_44,
            52 : self._read_spline4_msc_52,
        }
        try:
            n = op2.reader_geom2._read_double_card(
                card_name, card_obj,
                op2._add_methods._add_spline_object,
                methods, data, n)
        except DoubleCardError:
            raise
        return n
        #n = self._read_spline4_nx(data, n)
        #return n

    def _read_spline4_nx_44(self, spline: SPLINE4, data: bytes, n: int) -> tuple[int, SPLINE4]:
        """
        MSC 2018.2

        Word Name Type Description
        1 EID           I Spline element Identification
        2 CAERO         I Component Identifification
        3 AELIST        I AELIST Id for boxes
        4 SETG          I SETi Id for grids
        5 DZ           RS Smoothing Parameter
        6 METHOD(2) CHAR4 Method: IPS|TPS|FPS
        8 USAGE(2)  CHAR4 Usage flag: FORCE|DISP|BOTH
        10 NELEM        I Number of elements for FPS on x-axis
        11 MELEM        I Number of elements for FPS on y-axis
        12 FTYPE        I Radial interpolation function for METHOD=RIS  (not in NX)
        13 RCORE       RS Radius of radial interpolation function      (not in NX)

        """
        op2: OP2Geom = self.op2
        # 792/4 = 198
        # 198 = 2 * 99 = 2 * 11 * 9
        ntotal = 44 # 4 * 11
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        #structi = Struct(op2._endian + b'4i f 8s 8s 3i f') # msc
        structi = Struct(op2._endian + b'4i f 8s 8s 2i')
        splines = []
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            #eid, caero, aelist, setg, dz, method_bytes, usage_bytes, nelements, melements, ftype, rcore = out  # msc
            eid, caero, aelist, setg, dz, method_bytes, usage_bytes, nelements, melements = out
            method = method_bytes.rstrip().decode('ascii')
            usage = usage_bytes.rstrip().decode('ascii')
            spline = SPLINE4(eid, caero, aelist, setg,
                             dz, method, usage,
                             nelements, melements)
            str(spline)
            splines.append(spline)
            n += ntotal
        op2.to_nx(' because SPLINE4-NX was found')
        return n, splines

    def _read_spline4_msc_52(self, spline: SPLINE4, data: bytes, n: int) -> tuple[int, SPLINE4]:
        """
        MSC 2018.2

        Word Name Type Description
        1 EID           I Spline element Identification
        2 CAERO         I Component Identifification
        3 AELIST        I AELIST Id for boxes
        4 SETG          I SETi Id for grids
        5 DZ           RS Smoothing Parameter
        6 METHOD(2) CHAR4 Method: IPS|TPS|FPS
        8 USAGE(2)  CHAR4 Usage flag: FORCE|DISP|BOTH
        10 NELEM        I Number of elements for FPS on x-axis
        11 MELEM        I Number of elements for FPS on y-axis
        12 FTYPE        I Radial interpolation function for METHOD=RIS  (not in NX)
        13 RCORE       RS Radius of radial interpolation function      (not in NX)

        """
        op2: OP2Geom = self.op2
        # 792/4 = 198
        # 198 = 2 * 99 = 2 * 11 * 9
        ntotal = 52 # 4 * 13
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        #structi = Struct(op2._endian + b'4i f 8s 8s 3i f') # msc
        structi = Struct(op2._endian + b'4i f 8s 8s 2i if')
        splines = []
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            eid, caero, aelist, setg, dz, method_bytes, usage_bytes, nelements, melements, ftype, rcore = out  # msc
            #print(eid, caero, aelist, setg, dz, method_bytes, usage_bytes, nelements, melements, ftype, rcore)
            method = method_bytes.rstrip().decode('ascii')
            usage = usage_bytes.rstrip().decode('ascii')
            if ftype == 2:
                ftype = 'WF2'
            assert isinstance(ftype, str), ftype
            spline = SPLINE4(eid, caero, aelist, setg,
                             dz, method, usage,
                             nelements, melements,
                             ftype=ftype, rcore=rcore)
            str(spline)
            splines.append(spline)
            n += ntotal
        op2.to_msc(' because SPLINE4-MSC was found')
        return n, splines

    def read_spline5(self, data: bytes, n: int) -> int:
        """reads the SPLINE5 card"""
        op2: OP2Geom = self.op2
        card_name = 'SPLINE5'
        card_obj = SPLINE5
        methods = {
            60 : self._read_spline5_nx_60,
            68 : self._read_spline5_msc_68,
        }
        try:
            n = op2.reader_geom2._read_double_card(
                card_name, card_obj, op2._add_methods._add_spline_object,
                methods, data, n)
        except DoubleCardError:
            raise
        return n

    def _read_spline5_msc_68(self, spline: SPLINE5, data: bytes, n: int) -> tuple[int, list[SPLINE5]]:
        """
        reads the SPLINE5 card

        Word Name Type Description
        1 EID            I Spline element Identification
        2 CAERO          I Component Identifification
        3 AELIST         I AELIST Id for boxes
        4 SETG           I SETi Id for grids
        5 DZ            RS Smoothing Parameter
        6 DTORXY        RS Flexibility ratio in XY Plane
        7 CID            I Coordinate Sys. Id. for Beam CS
        8 DTHX          RS Smoothing/Attachment Flags for X rotations
        9 DTHY          RS Smoothing/Attachment Flags for Y rotations
        10 DTHZ         RS Smoothing/Attachment Flags for Z rotations
        11 USAGE(2)  CHAR4 Usage flag: FORCE|DISP|BOTH
        13 METHOD(2) CHAR4 Method: IPS|TPS|FPS|RIS
        15 DTORZY       RS Flexibility ratio in ZY Plane
        16 FTYPE         I Radial interpolation function for METHOD=RIS (not in NX)
        17 RCORE        RS Radius of radial interpolation function     (not in NX)

        """
        op2: OP2Geom = self.op2
        ntotal = 68 * self.factor # 4 * 17
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0, ndatai % ntotal
        if self.size == 4:
            structi = Struct(op2._endian + b'4i 2f i 3f 8s8s fif')
        else:
            asdf
            #structi = Struct(op2._endian + b'5q 2d q 2d 16s')

        splines = []
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            #ftype, rcore
            eid, caero, aelist, setg, dz, dtorxy, cid, dthx, dthy, dthz, usage_bytes, method_bytes, dtorzy, ftype_int, rcore = out
            method = reshape_bytes_block_size(method_bytes, self.size)
            usage = reshape_bytes_block_size(usage_bytes, self.size)
            #print(f'eid={eid} caero={caero} aelist={aelist} setg={setg} dz={dz} '
                  #f'dtorxy={dtorxy} cid={cid} dthx={dthx} dthy={dthy} dthz={dthz} '
                  #f'usage={usage!r} method={method!r} dtorzy={dtorzy} ftype={ftype} rcore={rcore}')
            assert method in ['IPS','TPS','FPS','RIS', 'BEAM'], method
            assert usage in ['FORCE','DISP','BOTH'], usage
            thx = dthx
            thy = dthy
            dtor = dtorzy
            if ftype_int == 2:
                ftype = 'WF2'
            spline = SPLINE5(
                eid, caero, aelist, setg, thx, thy,
                dz=dz, dtor=dtor, cid=cid,
                usage=usage, method=method,
                ftype=ftype, rcore=rcore,  # not in NX
            )
            spline.validate()
            str(spline)
            splines.append(spline)
            n += ntotal
        op2.to_msc(' because SPLINE5-MSC was found')
        return n, splines

    def _read_spline5_nx_60(self, spline: SPLINE5, data: bytes, n: int) -> tuple[int, list[SPLINE5]]:
        """
        reads the SPLINE5 card

        Word Name Type Description
        1 EID            I Spline element Identification
        2 CAERO          I Component Identifification
        3 AELIST         I AELIST Id for boxes
        4 SETG           I SETi Id for grids
        5 DZ            RS Smoothing Parameter
        6 DTORXY        RS Flexibility ratio in XY Plane
        7 CID            I Coordinate Sys. Id. for Beam CS
        8 DTHX          RS Smoothing/Attachment Flags for X rotations
        9 DTHY          RS Smoothing/Attachment Flags for Y rotations
        10 DTHZ         RS Smoothing/Attachment Flags for Z rotations
        11 USAGE(2)  CHAR4 Usage flag: FORCE|DISP|BOTH
        13 METHOD(2) CHAR4 Method: IPS|TPS|FPS|RIS
        15 DTORZY       RS Flexibility ratio in ZY Plane
        16 FTYPE         I Radial interpolation function for METHOD=RIS (not in NX?)
        17 RCORE        RS Radius of radial interpolation function     (not in NX?)

        """
        op2: OP2Geom = self.op2
        ntotal = 60 * self.factor # 4 * 12
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0, ndatai % ntotal
        if self.size == 4:
            structi = Struct(op2._endian + b'4i 2f i 3f 8s8s f')
        else:
            asdf
            #structi = Struct(op2._endian + b'5q 2d q 2d 16s')

        splines = []
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            #ftype, rcore
            eid, caero, aelist, setg, dz, dtorxy, cid, dthx, dthy, dthz, usage_bytes, method_bytes, dtorzy = out
            method = reshape_bytes_block_size(method_bytes, self.size)
            usage = reshape_bytes_block_size(usage_bytes, self.size)
            #print(f'eid={eid} caero={caero} aelist={aelist} setg={setg} dz={dz} dtorxy={dtorxy} cid={cid} dthx={dthx} dthy={dthy} dthz={dthz} usage={usage!r} method={method!r} dtorzy={dtorzy}')
            assert method in ['IPS','TPS','FPS','RIS'], method
            assert usage in ['FORCE','DISP','BOTH'], usage
            thx = dthx
            thy = dthy
            dtor = dtorzy
            spline = SPLINE5(
                eid, caero, aelist, setg, thx, thy,
                dz=dz, dtor=dtor, cid=cid,
                usage=usage, method=method,
                #ftype=ftype, rcore=rcore,  # not in NX
            )
            str(spline)
            splines.append(spline)
            n += ntotal
        op2.to_nx(' because SPLINE5-NX was found')
        return n, splines

    def read_monpnt1(self, data: bytes, n: int) -> int:
        """Reads the MONPNT1 card"""
        op2: OP2Geom = self.op2
        card_name = 'MONPNT1'
        card_obj = MONPNT1
        methods = {
            92 : self._read_monpnt1_nx_92,
            96 : self._read_monpnt1_96,
        }
        try:
            n = op2.reader_geom2._read_double_card(
                card_name, card_obj,
                op2._add_methods._add_monpnt_object,
                methods, data, n)
        except DoubleCardError:
            raise
        return n

    def _read_monpnt1_nx_92(self, monpnt1: MONPNT1, data: bytes, n: int) -> tuple[int, list[MONPNT1]]:
        """
        MSC 2018.2

        Word Name Type Description
        1 NAME(2)   CHAR4
        3 LABEL(14) CHAR4
        17 AXES         I
        18 COMP(2)  CHAR4
        20 CP           I
        21 X           RS
        22 Y           RS
        23 Z           RS

        """
        op2: OP2Geom = self.op2
        #ntotal = 4 * 24 # 4 * 24
        ntotal = 92 # 4 * 23
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        #structi = Struct(op2._endian + b'8s 56s i 8s i 3f i')  # msc
        structi = Struct(op2._endian + b'8s 56s i 8s i 3f')  # nx
        monpnt1s = []
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            #name_bytes, label_bytes, axes, comp_bytes, cp, x, y, z, cd = out
            name_bytes, label_bytes, axes, aecomp_name_bytes, cp, x, y, z = out
            name = reshape_bytes_block_size(name_bytes, self.size)
            label = reshape_bytes_block_size(label_bytes, self.size)
            aecomp_name = reshape_bytes_block_size(aecomp_name_bytes, self.size)
            xyz = [x, y, z]
            monpnt1 = MONPNT1(name, label, axes, aecomp_name,
                              xyz, cp=cp)
            str(monpnt1)
            n += ntotal
            monpnt1s.append(monpnt1)
        op2.to_nx(' because MONPNT1-NX was found')
        return n, monpnt1s

    def _read_monpnt1_96(self, monpnt1: MONPNT1, data: bytes, n: int) -> tuple[int, list[MONPNT1]]:
        """
        MSC 2018.2

        Word Name Type Description
        1 NAME(2)   CHAR4
        3 LABEL(14) CHAR4
        17 AXES         I
        18 COMP(2)  CHAR4
        20 CP           I
        21 X           RS
        22 Y           RS
        23 Z           RS
        24 CD           I

        """
        op2: OP2Geom = self.op2
        #ntotal = 4 * 24 # 4 * 24
        ntotal = 96 # 4 * 23
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        structi = Struct(op2._endian + b'8s 56s i 8s i 3f i')  # msc
        #structi = Struct(op2._endian + b'8s 56s i 8s i 3f')  # nx
        monpnt1s = []
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            name_bytes, label_bytes, axes, aecomp_name_bytes, cp, x, y, z, cd = out
            #name_bytes, label_bytes, axes, aecomp_name_bytes, cp, x, y, z = out
            name = reshape_bytes_block_size(name_bytes, self.size)
            label = reshape_bytes_block_size(label_bytes, self.size)
            aecomp_name = reshape_bytes_block_size(aecomp_name_bytes, self.size)
            xyz = [x, y, z]
            monpnt1 = MONPNT1(name, label, axes, aecomp_name,
                              xyz, cp=cp, cd=cd)
            str(monpnt1)
            n += ntotal
            monpnt1s.append(monpnt1)
        #op2.to_nx(' because MONPNT1-NX was found')
        return n, monpnt1s

    def read_monpnt2(self, data: bytes, n: int) -> int:
        """
        Record 59 - MONPNT2(8204,82,621)
        Word Name Type Description
        1 NAME(2)   CHAR4
        3 LABEL(14) CHAR4
        17 TABLE(2) CHAR4
        19 ELTYP(2) CHAR4
        21 ITEM(2)  CHAR4
        23 EID      I

        NX?
        Words 17 thru 23 repeat until -1 occurs
        Words 1 thru 23 repeat until (-2,-2) occurs
        """
        op2: OP2Geom = self.op2
        ntotal = 92 * self.factor # 4 * 23
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        structi = Struct(op2._endian + b'8s 56s 8s8s8s i')  # msc
        monpnts = []
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            name_bytes, label_bytes, table_bytes, eltype_bytes, item_bytes, eid = out
            name = reshape_bytes_block_size(name_bytes, self.size)
            label = reshape_bytes_block_size(label_bytes, self.size)
            table = reshape_bytes_block_size(table_bytes, self.size)
            Type = reshape_bytes_block_size(eltype_bytes, self.size)
            nddl_item = reshape_bytes_block_size(item_bytes, self.size)
            monpnt = MONPNT2(name, label, table, Type, nddl_item, eid, comment='')
            op2._add_methods._add_monpnt_object(monpnt)
            str(monpnt)
            #print(monpnt)
            n += ntotal
            monpnts.append(monpnt)
        #op2.to_nx(' because MONPNT3-NX was found')
        return n # , monpnt1s

    def read_monpnt3(self, data: bytes, n: int) -> int:
        """
        Record 60 - MONPNT3(8304,83,622)
        Word Name    Type   Description
        1  NAME(2)   CHAR4
        3  LABEL(14) CHAR4
        17 AXES      I      Axes to compute
        18 GRIDSET   I      GPF Grid Set
        19 ELEMSET   I      GPF Elem Set
        20 CID       I      Coord system x,y,z input in
        21 X         RS
        22 Y         RS
        23 Z         RS
        24 XFLAG     I      Exclude forces from class
        25 CD        I

        """
        op2: OP2Geom = self.op2
        ntotal = 100 * self.factor # 4 * 25
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        structi = Struct(op2._endian + b'8s 56s 4i 3f 2i')  # msc
        monpnts = []

        #XFLAG Exclusion flag. Exclude the indicated Grid Point Force types from summation at the
        #monitor point. Default = blank (no type excluded). See Remark 4.
        #S SPCforces
        #M MPC forces
        #A, L, or P applied loads
        #D dmig’s (and any other type not described above) at the monitored point.

        # A = L = P
        # 2 ^ 4 = 16

        # official: 0, 2, 4, 8, 16, 28, 30
        # guess:    6, 24, 26
        xflag_map = {
            0: None,
            #1: 'A',?
            2: 'S',
            #3: 'SA',?
            4: 'M',
            #5: 'MA',
            6: 'MS',
            #7: 'MSA',?
            8: 'A', # A = L = P

            16: 'D',

            24: 'DP',
            26: 'SDP',
            28: 'MAD',
            30: 'SMAD',

            # 32: C?
            62: 'SMADC',
        }
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            name_bytes, label_bytes, axes, grid_set, elem_set, cp, x, y, z, xflag, cd = out
            name = reshape_bytes_block_size(name_bytes, self.size)
            label = reshape_bytes_block_size(label_bytes, self.size)
            xyz = [x, y, z]
            try:
                xflag_str = xflag_map[xflag]
            except Exception:
                raise RuntimeError((name, label, xflag))
            monpnt = MONPNT3(name, label, axes, grid_set, elem_set, xyz,
                              cp=cp, cd=cd, xflag=xflag_str, comment='')
            op2._add_methods._add_monpnt_object(monpnt)
            str(monpnt)
            #print(monpnt)
            n += ntotal
            monpnts.append(monpnt)
        #op2.to_nx(' because MONPNT3-NX was found')
        return n # , monpnt1s

    def read_mondsp1(self, data: bytes, n: int) -> int:
        """
        Record 56 - MONDSP1(8804,88,628)
        Word Name Type Description
        1  NAME(2)   CHAR4
        3  LABEL(14) CHAR4
        17 AXES    I
        18 COMP(2) CHAR4
        20 CP      I
        21 X       RS
        22 Y       RS
        23 Z       RS
        24 CD      I
        25 INDDOF  I
        """
        op2: OP2Geom = self.op2
        ntotal = 100 * self.factor # 4 * 25
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        structi = Struct(op2._endian + b'8s 56s i8s i 3f 2i')  # msc
        monpnts = []
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            name_bytes, label_bytes, axes, aecomp_name_bytes, cp, x, y, z, cd, ind_dof = out
            name = reshape_bytes_block_size(name_bytes, self.size)
            label = reshape_bytes_block_size(label_bytes, self.size)
            aecomp_name = reshape_bytes_block_size(aecomp_name_bytes, self.size)
            xyz = [x, y, z]
            monpnt = MONDSP1(name, label, axes, aecomp_name, xyz, cp=cp, cd=cd,
                             ind_dof='123', comment='')
            op2._add_methods._add_monpnt_object(monpnt)
            str(monpnt)
            n += ntotal
            monpnts.append(monpnt)
        #op2.to_nx(' because MONPNT3-NX was found')
        return n # , monpnt1s

    def read_mdlprm(self, data: bytes, n: int) -> int:
        """
        Word Name Type Description
        1 NAME(2) CHAR4 User defined parameter
        3 VALUE I parameter value

          strings = (b'NSGRDS4 \x14\x00\x00\x00PEXTS4  \x00\x00\x00\x00SPBLNDX \xcd\xcc\xcc\xcc',)
          ints    = (b'NSGRDS4 ', 20, b'PEXTS4  ', 0,   b'SPBLNDX ', -858993459)
          floats  = (b'NSGRDS4 ', 20, b'PEXTS4  ', 0.0, b'SPBLNDX ', -107374184.0)
          MDLPRM, nsgrds4, 20, pexts4, 50., spblndx, 3.1
        """
        op2: OP2Geom = self.op2
        ntotal = 12 * self.factor # 4 * 3
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        structi = Struct(op2._endian + b'8s i')  # msc
        structf = Struct(op2._endian + b'8s f')  # msc
        data_dict = {}

        MDLPRM_FLOAT_KEYS_1 = {
            'DBCTOLE', 'DELELAS', 'DELFAST', 'DELMASS', 'DELSEAM', 'DELWELD',
            'PEXTS4', 'PIVTHRSH', 'SPBLNDX'}
        float_names = {('-%8s' % name).encode('ascii') for name in MDLPRM_FLOAT_KEYS_1}
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            name_bytes, value = out
            if name_bytes in float_names:
                name_bytes, value = structf.unpack(edata)
            name = reshape_bytes_block_size(name_bytes, self.size)

            if name == 'SHEARP':
                if value == 2:
                    value = 'HARDER'
                else:
                    raise NotImplementedError((name, value))
            elif name == 'OFFDEF':
                if value == 8:
                    value = 'ELMOFF'
                elif value == 65:
                    value = 'NODIFF'
                elif value == 96:
                    value = 'NOMASS'
                elif value in [128, 192]:
                    value = 'LROFF'
                else:
                    raise NotImplementedError((name, value))

            data_dict[name] = value
            n += ntotal

        if 'SPBLNDX' in data_dict:
            raise RuntimeError(f'SPBLNDX exists and has the wrong value...{data_dict}')

        if op2.mdlprm is not None:
            return n
        op2.add_mdlprm(data_dict)
        op2.to_msc(' because MDLPRM-MSC was found')
        #monpnt = MONDSP1(name, label, axes, aecomp_name, xyz, cp=cp, cd=cd,
                         #ind_dof='123', comment='')
        #op2._add_methods._add_monpnt_object(monpnt)
        return n # , monpnt1s

    def read_aestat(self, data: bytes, n: int) -> int:
        """
        MSC 2018.2

        Word Name Type Description
        1 ID I
        2 LABEL(2) CHAR4

        """
        op2: OP2Geom = self.op2
        ntotal = 12 * self.factor # 4 * 8
        if self.size == 4:
            structi = Struct(op2._endian + b'i 8s')
        else:
            structi = Struct(op2._endian + b'q 16s')
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0

        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            aestat_id, label_bytes = out
            label = reshape_bytes_block_size(label_bytes, self.size)
            aestat = op2.add_aestat(aestat_id, label)
            str(aestat)
            n += ntotal
        return n

    def read_flutter(self, data: bytes, n: int) -> int:
        r"""
        (3902, 39, 272)
        MSC 2018.2

        Word Name Type Description
        1 SID           I
        2 METHOD(2) CHAR4
        4 DENS          I
        5 MACH          I
        6 RFREQ         I
        7 IMETH(2)  CHAR4
        SFLG=0 (std)
          9 NEIGN  I  nvalue
        SFLG=1 (sweep)
          9 FMAX  RS maximum frequency
        End SFLG
        10 EPR  RS
        11 SFLG  I SWEEP FLAG
        Words 1 through max repeat until End of Record

        NX:
                  sid method,     d, m, k, imethod,  neign, epr, end
          data = (30, PK,         1, 2, 3, L,         3, 0.001, -1)       # NX
          data = (30, KE, '    ', 1, 2, 3, L, '    ', 3, 0.001, 0.0, -1)  # MSC

        C:\MSC.Software\msc_nastran_runs\pkswep.op2
                 sid method,    d, m, k, imethod,   neign/fmax, epr,   sweep_flag
        data  = (3, 'PKS     ', 1, 2, 3, 'L       ', 5.0,       0.01,  1, -1,
                 4, 'K       ', 1, 2, 4, 'L       ', 3,         0.001, 0, -1)
        data  = (3, 'PKS     ', 1, 2, 3, 'L       ', 5.0,       0.01,  1, -1,
                 4, 'K       ', 1, 2, 4, 'L       ', 3,         0.001, 0, -1)
        """
        op2: OP2Geom = self.op2
        ints = np.frombuffer(data[n:], op2.idtype).copy()
        floats = np.frombuffer(data[n:], op2.fdtype).copy()
        istart, iend = get_minus1_start_end(ints)

        for (i0, i1) in zip(istart, iend):
            sid = ints[i0]
            assert ints[i1] == -1, ints[i1]
            method_bytes = data[n+i0*4+4:n+i0*4+12]

            density = ints[i0+3]
            mach = ints[i0+4]
            reduced_freq_velocity = ints[i0+5]
            imethod_bytes = data[n+i0*4+24:n+i0*4+32]
            nvalue = ints[i0+8]  # nvalue
            fmax = None
            if ints[i0+10] == -1:
                epsilon = floats[i0+9]
                assert ints[i0+10] == -1, ints[i0:i1]
                op2.to_nx(' because FLUTTER was found')
            else:
                # msc
                op2.to_msc(' because FLUTTER was found')
                sweep_flag = ints[i0+10]
                assert ints[i0+11] == -1, ints[i0:i1+1]
                if sweep_flag == 0:
                    nvalue = ints[i0+8]
                elif sweep_flag == 1:
                    nvalue = None
                    fmax = floats[i0+8]
                else:
                    raise RuntimeError(sweep_flag)
                epsilon = floats[i0+9]

                #if fmax is None:
                    #op2.log.debug(f'FLUTTER: 9: nvalue={nvalue} 10: epsilon={epsilon:g}')
                #else:
                    #op2.log.debug(f'FLUTTER: 9: fmax={fmax}; 10: epsilon={epsilon:g}')
            method = method_bytes.rstrip().decode('ascii')
            imethod = imethod_bytes.rstrip().decode('ascii')
            flutter = FLUTTER(
                sid, method,
                density, mach, reduced_freq_velocity,
                imethod=imethod, # 'L'
                nvalue=nvalue,
                epsilon=epsilon,
                omax=fmax,
                validate=True)
            op2._add_methods._add_flutter_object(flutter, allow_overwrites=True)
            #op2.add_flutter(sid, method,
                            #density, mach, reduced_freq_velocity,
                            #imethod=imethod, # 'L'
                            #nvalue=nvalue,
                            #epsilon=epsilon,
                            #omax=fmax,
                            #validate=True)
        return len(data)
        #ntotal = 12 # 4 * 8
        #ndatai = len(data) - n
        #ncards = ndatai // ntotal
        #assert ndatai % ntotal == 0
        #structi = Struct(op2._endian + b'i 8s')
        #for unused_i in range(ncards):
            #edata = data[n:n + ntotal]
            #out = structi.unpack(edata)
            #aestat_id, label = out
            #label = label.rstrip().decode('latin1')
            #op2.add_aestat(aestat_id, label)
            #n += ntotal
        #return n

    def read_trim(self, data: bytes, n: int) -> int:
        """
        (2402, 24, 342)
        MSC 2018.2

        Word Name Type Description
        1 ID           I
        2 MACH        RS
        3 Q           RS
        4 AEQRATIO    RS

        5 LABEL(2) CHAR4
        7 UX          RS
        Words 5 through 7 repeat until (-1,-1,-1) occurs

        """
        op2: OP2Geom = self.op2
        ntotal1 = 16 * self.factor # 4 * 4
        ntotal2 = 12 * self.factor # 4 * 3
        #ndatai = len(data) - n
        #ncards = ndatai // ntotal
        if self.size == 4:
            struct1 = Struct(op2._endian + b'i 3f')
            struct2 = Struct(op2._endian + b'8sf')
            struct_end = Struct(op2._endian + b'3i')
        else:
            struct1 = Struct(op2._endian + b'q 3d')
            struct2 = Struct(op2._endian + b'16sd')
            struct_end = Struct(op2._endian + b'3q')
        while n < len(data):
            edata = data[n:n+ntotal1]
            trim_id, mach, q, aeqr = struct1.unpack(edata)
            n += ntotal1
            labels = []
            uxs = []

            edata = data[n:n+ntotal2]
            while struct_end.unpack(edata) != (-1, -1, -1):
                label_bytes, ux = struct2.unpack(edata)
                label = reshape_bytes_block_size(label_bytes, self.size)
                labels.append(label)
                uxs.append(ux)
                n += ntotal2
                edata = data[n:n+ntotal2]
            n += ntotal2
            trim = op2.add_trim(trim_id, mach, q, labels, uxs, aeqr=aeqr, trim_type=1)
            str(trim)
        return n

    def read_aesurf(self, data: bytes, n: int) -> int:
        """
        MSC 2018.2

        Word Name Type Description
        1 ID           I Identification of an aerodynamic trim variable degree of freedom >0, no default
        2 LABEL(2) CHAR4 Control Surface (CS) name, no default
        4 CID1         I IDentification of a rectangular Coordinate system with y-axis that defines the hinge line of the CS component
        5 ALID1        I IDentification of an AELIST bulk data entry that identifies all aerodynamic elements that make up the CS comp
        6 CID2         I IDentification of a rectangular Coordinate system with y-axis that defines the hinge line of the CS component
        7 ALID2        I IDentification of an AELIST bulk data entry that identifies all aerodynamic elements that make up the CS comp
        8 EFF         RS Control surface EFFectiveness, default=1.0
        9 LDW          I =0 create a linear down wash, >0 no linear downwash
        10 CREFC      RS Reference Chord Length for the CS > 0.0, default = 1.0
        11 CREFS      RS Reference area for the CS > 0.0, default = 1.0
        12 PLLIM      RS Lower deflection   Limit for the control surface in radians, default=no limit
        13 PULIM      RS Upper deflection   Limit for the control surface in radians, default=no limit
        14 HMLLIM     RS Lower Hinge Moment Limit for the control surface in force-length units, default=no limit
        15 HMULIM     RS Upper Hinge Moment Limit for the control surface in force-length units, default=no limit
        16 TQLLIM     RS Lower deflection   Limit for the control surface as fct(q), >0, default=no limit
        17 TQULIM     RS Upper deflection   Limit for the control surface as fct(q), >0, default=no limit

        """
        op2: OP2Geom = self.op2
        ntotal = 68 *  self.factor # 4 * 17
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0

        ntotali = -24 * self.factor
        if self.size == 4:
            struct2 = Struct(op2._endian + b'6f')
            struct1 = Struct(op2._endian + b'i 8s 4i fi 2f 4s4s 4s4s 4s4s')
            nan = b'    '
        else:
            struct2 = Struct(op2._endian + b'6d')
            struct1 = Struct(op2._endian + b'q 16s 4q dq 2d 8s8s 8s8s 8s8s')
            nan = b'        '

        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out1 = struct1.unpack(edata)
            (aesurf_id, label_bytes,
             cid1, aelist_id1,
             cid2, aelist_id2,
             eff, ldw_int, crefc, crefs,
             pllim, pulim, hmllim, hmulim, tqllim, tqulim) = out1
            #print(out1)
            label = reshape_bytes_block_size(label_bytes, self.size)

            if ldw_int == 0:
                ldw = 'LDW'
            elif ldw_int == 1:
                ldw = 'NOLDW'
            else:
                raise NotImplementedError(ldw_int)
            assert isinstance(ldw, str), ldw

            if (pllim, pulim, hmllim, hmulim, tqllim, tqulim) == (nan, nan, nan, nan, nan, nan):
                pllim = None
                pulim = None
                hmllim = None
                hmulim = None
                tqllim = None
                tqulim = None
            else:
                pllim2, pulim2, hmllim2, hmulim2, tqllim2, tqulim2 = struct2.unpack(edata[ntotali:])
                pllim = pllim2 if pllim != nan else None
                pulim = pulim2 if pulim != nan else None

                hmllim = hmllim2 if hmllim != nan else None
                hmulim = hmulim if hmulim != nan else None

                tqllim = tqllim2 if tqllim != nan else None
                tqulim = tqulim2 if tqulim != nan else None
                #print('pllim, pulim, hmllim, hmulim, tqllim, tqulim', pllim, pulim, hmllim, hmulim, tqllim, tqulim)

            op2.add_aesurf(aesurf_id, label, cid1, aelist_id1, cid2=None, aelist_id2=None,
                           eff=eff, ldw=ldw, crefc=crefc, crefs=crefs,
                           #pllim=-np.pi/2., pulim=np.pi/2.,
                           pllim=pllim, pulim=pulim,
                           hmllim=hmllim, hmulim=hmulim, # hinge moment lower/upper limits
                           tqllim=tqllim, tqulim=tqulim)
            n += ntotal
        return n

    def read_aesurfs(self, data: bytes, n: int) -> int:
        """
        Word Name Type Description
        1 ID       I     Identification of an aerodynamic trim variable degree
                         of freedom >0, no default
        2 LABEL(2) CHAR4 Control Surface (CS) name, no default
        4 LIST1    I     Identification of a SET1 that contains the grids ids
                         associated with this control surface
        5 LIST2    I     Identification of a SET1 that contains the grids ids
                         associated with this control surface
        """
        op2: OP2Geom = self.op2
        ntotal = 20 *  self.factor # 4 * 5
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        if self.size == 4:
            struct1 = Struct(op2._endian + b'i 8s 2i')
        else:
            struct1 = Struct(op2._endian + b'i 16s 2i')

        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out1 = struct1.unpack(edata)

            aesid, label_bytes, list1, list2 = out1
            label = reshape_bytes_block_size(label_bytes, self.size)

            aesurfs = op2.add_aesurfs(aesid, label, list1, list2)
            str(aesurfs)
            n += ntotal
        return n

    def read_aefact(self, data: bytes, n: int) -> int:
        """
        MSC 2018.2

        Word Name Type Description
        1 SID I
        2 D RS
        Word 2 repeats until End of Record

        (1, 0.0, 0.1, 0.2, 1.0, -1,
         2, 0.0, 0.1, 0.2, 0.5, 1.0, -1,
        )

        """
        op2: OP2Geom = self.op2
        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        floats = np.frombuffer(data[n:], op2.fdtype8).copy()
        istart, iend = get_minus1_start_end(ints)

        for (i0, i1) in zip(istart, iend):
            #if i0 == i1:
                #continue
            sid = ints[i0]
            fractions = floats[i0+1:i1]
            assert ints[i1] == -1, ints[i1]
            aefact = op2.add_aefact(sid, fractions)
            #print(aefact)
            str(aefact)
            #n += ntotal
        return len(data)


def _expand_vals(grids):
    grids2 = []
    for val in grids:
        #> 0 for ID
        #= 0 for THRU
        #= -6 for BY
        #= -7 for ALL
        if val > 0:
            pass
        elif val == 0:
            val = 'THRU'
        elif val == -6:
            val = 'BY'
        elif val == -7:
            val = 'ALL'
        else:
            raise NotImplementedError(f'val={val} data={grids}')
        grids2.append(val)
    return grids2


def _read_group_elem_prop_nids(ints, i, n, size) -> tuple[int, int, Any]:
    """helper for _read_group"""
    i += 1
    n += size

    # grids
    #iminus1 = minus1[minus1_count] # + 1
    #print(ints[iminus1:])
    #grids = ints[i+1:iminus1].tolist()
    #print('ints[i:]', ints[i:])
    assert ints[i:][0] > 0, ints[i:]
    for j, nj in enumerate(ints[i:]):
        if nj == -1:
            break
    grids = ints[i+1:i+j].tolist()
    print('grids', grids)
    grids2 = _expand_vals(grids)
    print(f'  grids = {grids2}')
    assert 'THRU' != grids2[0]
    assert 'BY' != grids2[0]
    assert 'ALL' != grids2[0]

    #minus1_count += 1

    nstop = len(grids) + 2
    i += nstop
    n += nstop * size
    #i = iminus1
    #n = iminus1 * 4
    return grids2
