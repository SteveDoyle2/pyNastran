"""
defines readers for BDF objects in the OP2 EDT/EDTS table
"""
from struct import Struct
from typing import Tuple, Any

import numpy as np
from pyNastran.op2.tables.geom.geom_common import GeomCommon
from pyNastran.op2.op2_interface.op2_reader import mapfmt, reshape_bytes_block


class EDT(GeomCommon):
    """defines methods for reading aero and element deformations"""

    def _read_edt_4(self, data, ndata):
        """
        3.21 EDT
        Aero and element deformations.

        """
        return self._read_geom_4(self._edt_map, data, ndata)

    def __init__(self):
        GeomCommon.__init__(self)

        # F:\Program Files\Siemens\NXNastran\nxn10p1\nxn10p1\nast\tpl\fsw_eng.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_boltld04i.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_eliter17.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_weld01i.op2
        # F:\work\pyNastran\examples\Dropbox\move_tpl\ac10901a_new.op2
        self._edt_map = {
            (5201, 52, 373) : ['ACMODL', self._read_acmodl],
            (6301, 63, 397) : ['ADAPT', self._read_fake],
            (7801, 78, 582) : ['AECOMP', self._read_aecomp],
            (7901, 79, 583) : ['AECOMPL', self._read_aecompl],
            (7301, 73, 574) : ['AEDW', self._read_fake],
            (4002, 40, 273) : ['AEFACT', self._read_aefact],
            (7501, 75, 576) : ['AEFORCE', self._read_aeforce],
            (2602, 26, 386) : ['AELINK', self._read_aelink],
            (2302, 23, 341) : ['AELIST', self._read_aelist],
            (7001, 70, 571) : ['AEPARM', self._read_fake],
            (7401, 74, 575) : ['AEPRESS', self._read_aepress],
            (3202, 32, 265) : ['AERO', self._read_aero],
            (2202, 22, 340) : ['AEROS', self._read_aeros],
            (2102, 21, 339) : ['AESTAT', self._read_aestat],
            (2002, 20, 338) : ['AESURF', self._read_aesurf],
            (7701, 77, 581) : ['AESURFS', self._read_fake],
            (3002, 30, 263) : ['CAERO1', self._read_caero1],
            (4301, 43, 167) : ['CAERO2', self._read_caero2],
            (4401, 44, 168) : ['CAERO3', self._read_caero3],
            (4501, 45, 169) : ['CAERO4', self._read_caero4],
            (5001, 50, 175) : ['CAERO5', self._read_caero5],
            (6201, 62, 143) : ['CLOAD', self._read_fake],
            (6401, 64, 307) : ['CSSCHD', self._read_fake],
            (104, 1, 81) : ['DEFORM', self._read_deform],
            (2702, 27, 387) : ['DIVERG', self._read_fake],
            (4102, 41, 274) : ['FLFACT', self._read_flfact],
            (3902, 39, 272) : ['FLUTTER', self._read_flutter],
            (17400, 174, 616) : ['GROUP', self._read_group],
            (3802, 38, 271) : ['MKAERO1', self._read_mkaero1],
            (3702, 37, 270) : ['MKAERO2', self._read_mkaero2],
            (7601, 76, 577) : ['MONPNT1', self._read_monpnt1],
            (3102, 31, 264) : ['PAERO1', self._read_paero1],
            (4601, 46, 170) : ['PAERO2', self._read_paero2],
            (4701, 47, 171) : ['PAERO3', self._read_fake],
            (4801, 48, 172) : ['PAERO4', self._read_fake],
            (5101, 51, 176) : ['PAERO5', self._read_paero5],
            (5301, 53, 378) : ['PANEL', self._read_panel],
            (3502, 35, 268) : ['SET1', self._read_set1],
            (3602, 36, 269) : ['SET2', self._read_set2],
            (4302, 43, 607) : ['SET3', self._read_set3],
            (3302, 33, 266) : ['SPLINE1', self._read_spline1],
            (3402, 34, 267) : ['SPLINE2', self._read_spline2],
            (4901, 49, 173) : ['SPLINE3', self._read_spline3],
            (6501, 65, 308) : ['SPLINE4', self._read_spline4],
            (6601, 66, 309) : ['SPLINE5', self._read_spline5],
            (2402, 24, 342) : ['TRIM', self._read_trim],
            (7201, 72, 573) : ['UXVEC', self._read_fake],
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
        }
    def _read_aeforce(self, data: bytes, n: int) -> int:
        aeforcex
    def _read_aepress(self, data: bytes, n: int) -> int:
        aepressx
    def _read_mkaero2(self, data: bytes, n: int) -> int:
        mkaero2x


    def _read_flfact(self, data: bytes, n: int) -> int:
        """
        data = (1, 0.206, -1,
                2, 1.3, -1,
                3, 14400.0, 15600.0, 16800.0, 18000.0, 19200.0, 20400.0, -1)
        """
        ints = np.frombuffer(data[n:], self.idtype).copy()
        floats = np.frombuffer(data[n:], self.fdtype).copy()
        iminus1 = np.where(ints == -1)[0]

        istart = [0] + list(iminus1[:-1] + 1)
        iend = iminus1

        for (i0, i1) in zip(istart, iend):
            sid = ints[i0]
            factors = floats[i0+1:i1]
            assert ints[i1] == -1, ints[i1]
            flfact = self.add_flfact(sid, factors)
            str(flfact)
        return len(data)

    def _read_mkaero1(self, data: bytes, n: int) -> int:
        """
        (3802, 38, 271)
        Kinda brilliant way to write the card.  Weird to parse though.

        data = (1.3, -1, -1, -1, -1, -1, -1, -1,
                0.03, 0.04, 0.05, -1, -1, -1, -1, -1)
        """
        assert len(data) == 76, len(data)
        nvalues = (len(data) - n) // 4
        nrows = nvalues // 16
        assert nrows > 0, nrows
        ints = np.frombuffer(data[12:], dtype=self.idtype).reshape(nrows, 16)
        floats = np.frombuffer(data[12:], dtype=self.fdtype).reshape(nrows, 16)
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
            mkaero1 = self.add_mkaero1(machs, kfreqs)
            str(mkaero1)
        return len(data)

    def _read_group(self, data: bytes, n: int) -> int:
        """
        GROUP(17400,174,616)

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
        #print('reading group')
        #assert self.factor == 1, self.factor
        nentries = 0
        ints = np.frombuffer(data[n:], dtype=self.idtype8)
        if self.factor == 1:
            strs = np.frombuffer(data[n:], dtype='|S4')
        else:
            self.show_data(data[n:], types='qds')
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

            if self.factor == 1:
                group_desc = ''.join(stri.decode('latin1') for stri in strs[i:i+ndesc]).strip()
            else:
                group_desc_bytes = reshape_bytes_block(b''.join(strs[i:i+ndesc]))
                group_desc = group_desc_bytes.decode('latin1').rstrip()
            i += ndesc
            n += 4 * ndesc

            #------------------------------
            #gtype, nmeta, mdesc
            gtype = ints[i]
            #i += 1
            #n += 4
            print(f'group_id={group_id} ndesc={ndesc} group_desc={group_desc!r}; gtype={gtype!r}')

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
                print(f'-----gtype={gtype}----')
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
                    print(f'  gtype={gtype} nmeta={nmeta} meta_desc={meta_desc!r}')
                    #iminus1 = minus1[minus1_count+2]
                    #print('ints: ', ints[i:iminus1].tolist())
                    #minus1_count += 1
                elif gtype == -3:
                    assert ints[i] == -3, ints[i]
                    i, n, props = _read_group_elem_prop_nids(ints, i, n, size)
                    data_dict['prop'].append(props)
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
                    print(f'  eids2 = {eids2}')
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
            print(data_dict)
            #self.add_group(group_id, group_desc, data_dict)
                #i += 1
                #n += 4
            #assert ints[i] == -1, ints[i:]
            self.log.warning(f'skipping GROUP in {self.table_name}')
            nentries += 1

        assert n == len(data), f'n={n} ndata={len(data)}'
        self.increase_card_count('GROUP', nentries)
        assert nentries > 0, nentries
        return n

    def _read_aero(self, data: bytes, n: int) -> int:
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
        assert len(data) == 36, len(data)
        struct = Struct(self._endian + b'i 3f 2i')
        out = struct.unpack(data[n:])
        acsid, velocity, cref, rho_ref, sym_xz, sym_xy = out
        self.add_aero(velocity, cref, rho_ref,
                      acsid=acsid, sym_xz=sym_xz, sym_xy=sym_xy)
        n = 36
        return n

    def _read_aeros(self, data: bytes, n: int) -> int:
        """
        AEROS(2202, 22, 340)

        AEROS   0       100     36.     360.    12960.
        data = (0, 100, 36.0, 360.0, 12960.0, 0, 0)

        """
        assert len(data) == 40 * self.factor, len(data)
        struct = Struct(mapfmt(self._endian + b'2i 3f 2i', self.size))
        out = struct.unpack(data[n:])
        acsid, rcsid, cref, bref, sref, sym_xz, sym_xy = out
        self.add_aeros(cref, bref, sref,
                       acsid=acsid, rcsid=rcsid,
                       sym_xz=sym_xz, sym_xy=sym_xy)
        n = 40 * self.factor
        return n

    def _read_deform(self, data: bytes, n: int) -> int:
        """
        (104, 1, 81)
        NX 2019.2

        Word Name Type Description
        1 SID I Deformation set identification number
        2 EID I Element number
        3 D RS Deformation

        """
        ntotal = 12 * self.factor # 4*3
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        structi = Struct(mapfmt(self._endian + b'2i f', self.size))
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            sid, eid, deformation = out
            deform = self.add_deform(sid, eid, deformation)
            str(deform)
            n += ntotal
        return n

    def _read_caero1(self, data: bytes, n: int) -> int:
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
        ntotal = 64 * self.factor # 4*16
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        structi = Struct(mapfmt(self._endian + b'8i 8f', self.size))
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            eid, pid, cp, nspan, nchord, lspan, lchord, igid, x1, y1, z1, x12, x4, y4, z4, x43 = out
            self.add_caero1(eid, pid, igid,
                            [x1, y1, z1], x12,
                            [x4, y4, z4], x43,
                            cp=cp,
                            nspan=nspan, lspan=lspan,
                            nchord=nchord, lchord=lchord)
            n += ntotal
        return n

    def _read_caero2(self, data: bytes, n: int) -> int:
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
        ntotal = 64 * self.factor # 4*16
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        structi = Struct(mapfmt(self._endian + b'8i 4f 4i', self.size))
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            eid, pid, cp, nsb, nint, lsb, lint, igroup, x1, y1, z1, x12, zero1, zero2, zero3, zero4 = out
            assert min(zero1, zero2, zero3, zero4) == max(zero1, zero2, zero3, zero4)
            p1 = [x1, y1, z1]
            caero2 = self.add_caero2(eid, pid, igroup, p1, x12,
                                     cp=cp,
                                     nsb=nsb, nint=nint,
                                     lsb=lsb, lint=lint)
            str(caero2)
            n += ntotal
        return n

    def _read_caero3(self, data: bytes, n: int) -> int:
        caero3x
    def _read_caero4(self, data: bytes, n: int) -> int:
        caero4x

    def _read_caero5(self, data: bytes, n: int) -> int:
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
        ntotal = 64 # 4*16
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        structi = Struct(self._endian + b'8i 8f')
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            eid, pid, cp, nspan, lspan, ntheory, nthick, undef, x1, y1, z1, x12, x4, y4, z4, x43 = out
            p1 = [x1, y1, z1]
            p4 = [x4, y4, z4]
            caero5 = self.add_caero5(eid, pid,
                   p1, x12,
                   p4, x43,
                   cp=cp,
                   nspan=nspan, lspan=lspan,
                   ntheory=ntheory,
                   nthick=nthick)
            str(caero5)
            n += ntotal
        return n

    def _read_paero1(self, data: bytes, n: int) -> int:
        """
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

        """
        ntotal = 32 * self.factor # 4 * 8
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        structi = Struct(mapfmt(self._endian + b'8i', self.size))
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            pid, b1, b2, b3, b4, b5, b6, empty = out
            caero_body_ids = []
            for body in [b1, b2, b3, b4, b5, b6, empty]:
                if body != 0:
                    caero_body_ids.append(body)
            paero1 = self.add_paero1(pid, caero_body_ids=caero_body_ids)
            if caero_body_ids:
                self.log.warning(str(paero1))
            n += ntotal
        return n

    def _read_paero2(self, data: bytes, n: int) -> int:
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

        PAERO1  100001

        data = (100001, 100001, 0, 10, 0, 0, 24, 1,
                99.3, 21.45, -11.65, 42.86, 101.8387, 122.62, -2.69, 32.71)

        """
        ntotal = 60 * self.factor # 4 * 15
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        if self.size == 4:
            structi = Struct(self._endian + b'i4si 2f 10i')
        else:
            structi = Struct(self._endian + b'q8sq 2d 10q')
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
            orient = cast_string(orient_bytes, self.size, encoding='latin1')
            paero2 = self.add_paero2(pid, orient, width, ar,
                                     thi, thn,
                                     lrsb=lrsb,
                                     lrib=lrib,
                                     lth=lth)
            n += ntotal
            str(paero2)
        return n

    def _read_paero5(self, data: bytes, n: int) -> int:
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
        ints = np.frombuffer(data[n:], self.idtype8).copy()
        floats = np.frombuffer(data[n:], self.fdtype8).copy()
        iminus1 = np.where(ints == -1)[0]

        istart = [0] + list(iminus1[:-1] + 1)
        iend = iminus1

        for (i0, i1) in zip(istart, iend):
            pid, nalpha, lalpha, nxis, lxis, ntaus, ltaus = ints[i0:i0+7]
            caoci = floats[i0+7:i1]
            assert ints[i1] == -1, ints[i1]
            paero5 = self.add_paero5(
                pid, caoci,
                nalpha=nalpha, lalpha=lalpha,
                nxis=nxis, lxis=lxis,
                ntaus=ntaus, ltaus=ltaus)
            str(paero5)
        return len(data)

    def _read_panel(self, data: bytes, n: int) -> int:
        """
        MSC 2018.2

        Word Name Type Description
        1 NAME(2) CHAR4
        3 SETID I
        Words 1 through 3 repeat until End of Record

        ('PANEL1', 1, -1)
        """
        ints = np.frombuffer(data[n:], self.idtype8).copy()
        iminus1 = np.where(ints == -1)[0]

        istart = [0] + list(iminus1[:-1] + 1)
        iend = iminus1

        for (i0, i1) in zip(istart, iend):
            assert ints[i1] == -1, ints[i1]
            names = []
            set_ids = []
            while i0 < i1:
                name_bytes = data[n:n+8]
                name = cast_string(name_bytes, self.size, encoding='latin1')
                set_id = ints[i0+2]
                names.append(name)
                set_ids.append(set_id)
                n += 12
                i0 += 3
            panel = self.add_panel(names, set_ids)
            str(panel)
        return len(data)

    def _read_acmodl(self, data: bytes, n: int) -> int:
        """
        MSC 2018.2
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

        NX 2019.2
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
        ntotal = 72 *  self.factor # 4 * 8
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0, ndatai % ntotal
        #structi = Struct(self._endian + b'4i f 8s 8s 3i f') # msc
        if self.size == 4:
            structi = Struct(self._endian + b'8s 8s 2i f 8s f 8s f ifi 8s')
        else:
            structi = Struct(self._endian + b'16s 16s 2q d 16s d 16s d qdq 16s')
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
            inter = cast_string(inter_bytes, self.size, encoding='latin1')
            infor = cast_string(infor_bytes, self.size, encoding='latin1')
            method = cast_string(method_bytes, self.size, encoding='latin1')
            search_unit = cast_string(search_unit_bytes, self.size, encoding='latin1')
            ctype = cast_string(ctype_bytes, self.size, encoding='latin1')

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
            # set2 = self.add_set2(sid, macro, sp1, sp2, ch1, ch2, zmax, zmin)
            acmodl = self.add_acmodl(infor, fset, sset,
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
            n += ntotal
        return n

    def _read_aelist(self, data: bytes, n: int) -> int:
        """
        MSC 2018.2

        Word Name Type Description
        1 SID I
        2 E I
        Word 2 repeats until End of Record

        """
        #self.show_data(data[12:], types='if')
        ints = np.frombuffer(data[n:], self.idtype8).copy()
        floats = np.frombuffer(data[n:], self.fdtype8).copy()
        iminus1 = np.where(ints == -1)[0]

        istart = [0] + list(iminus1[:-1] + 1)
        iend = iminus1

        for (i0, i1) in zip(istart, iend):
            sid = ints[i0]
            elements = floats[i0+1:i1]
            assert ints[i1] == -1, ints[i1]
            self.add_aelist(sid, elements)
            #n += ntotal
        return len(data)

    def _read_set1(self, data: bytes, n: int) -> int:
        """
        SET1: (3502, 35, 268)
        MSC 2018.2

        Word Name Type Description
        1 SID I
        2 G1  I Grid ID or -2 when SKIN is specified
        Word 2 repeats until End of Record

        """
        ints = np.frombuffer(data[n:], self.idtype8).copy()
        iminus1 = np.where(ints == -1)[0]

        istart = [0] + list(iminus1[:-1] + 1)
        iend = iminus1

        for (i0, i1) in zip(istart, iend):
            sid = ints[i0]
            elements = ints[i0+1:i1].tolist()
            assert -2 not in elements, elements
            assert ints[i1] == -1, ints[i1]
            self.add_set1(sid, elements, is_skin=False)
        return len(data)

    def _read_set2(self, data: bytes, n: int) -> int:
        """
        SET2
        MSC 2018.2

        Word Name Type Description
        1 SID I
        2 G1  I Grid ID or -2 when SKIN is specified
        Word 2 repeats until End of Record

        Record 71 - SET2(3602,36,269)
        Word Name Type Description
        1 SID I
        2 MACRO I
        3 SP1 RS
        4 SP2 RS
        5 CH1 RS
        6 CH2 RS
        7 ZMAX RS
        8 ZMIN RS

        """
        #self.show_data(data)
        ntotal = 32 # 4 * 8
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        #structi = Struct(self._endian + b'4i f 8s 8s 3i f') # msc
        structi = Struct(self._endian + b'2i 6f')
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            sid, macro, sp1, sp2, ch1, ch2, zmax, zmin = out
            set2 = self.add_set2(sid, macro, sp1, sp2, ch1, ch2, zmax, zmin)
            str(set2)
            n += ntotal
        self.to_nx()
        return n

    def _read_set3(self, data: bytes, n: int) -> int:

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
        # this is setup for NX
        ints = np.frombuffer(data[n:], self.idtype8).copy()
        iminus1 = np.where(ints == -1)[0]

        istart = [0] + list(iminus1[:-1] + 1)
        iend = iminus1

        for (i0, i1) in zip(istart, iend):
            sid = ints[i0]
            desc_int = ints[i0+1]
            elements = ints[i0+2:i1].tolist()
            if desc_int == 1:
                desc = 'ELEM'
            elif desc_int == 2:
                desc = 'GRID'
            elif desc_int == 3:
                desc = 'PROP'
            else:
                raise NotImplementedError(desc_int)

            assert min(elements) > 0, elements
            assert ints[i1] == -1, ints[i1]
            set3 = self.add_set3(sid, desc, elements)
            str(set3)
        return len(data)

    def _read_aelink(self, data: bytes, n: int) -> int:
        """
        MSC 2018.2

        Word Name Type Description
        1 ID I
        2 LABLD(2) CHAR4
        4 LABLI(2) CHAR4
        6 C1 RS
        Words 4 through 6 repeat until (-1,-1,-1) occurs

        """
        struct1 = Struct(self._endian + b'i8s')
        struct2 =Struct(self._endian + b'8sf')
        struct_end = Struct(self._endian + b'3i')
        ntotal = 12
        while n < len(data):
            edata = data[n:n+ntotal]
            aelink_id, label_bytes = struct1.unpack(edata)
            label = cast_string(label_bytes, self.size, encoding='latin1')
            n += ntotal
            linking_coefficents = []
            independent_labels = []

            edata = data[n:n+ntotal]
            while struct_end.unpack(edata) != (-1, -1, -1):
                ind_label_bytes, coeff = struct2.unpack(edata)
                ind_label = cast_string(ind_label_bytes, self.size, encoding='latin1')
                independent_labels.append(ind_label)
                linking_coefficents.append(coeff)
                n += ntotal
                edata = data[n:n+ntotal]
            n += ntotal
            self.add_aelink(aelink_id, label,
                            independent_labels, linking_coefficents)
        return len(data)

    def _read_aecomp(self, data: bytes, n: int) -> int:
        """
        MSC 2018.2

        Word Name Type Description
        1 NAME(2)     CHAR4
        3 LISTTYPE(2) CHAR4
        5 LISTID          I
        Word 5 repeats until End of Record

        """
        ints = np.frombuffer(data[n:], self.idtype).copy()
        #floats = np.frombuffer(data[n:], self.fdtype).copy()
        iminus1 = np.where(ints == -1)[0]

        istart = [0] + list(iminus1[:-1] + 1)
        iend = iminus1

        for (i0, i1) in zip(istart, iend):
            name_bytes = data[n+i0*4:n+i0*4+8]
            list_type_bytes = data[n+i0*4+8:n+i0*4+16]
            lists = ints[i0+4:i1].tolist()
            assert ints[i1] == -1, ints[i1]
            name = name_bytes.rstrip().decode('ascii')
            list_type = list_type_bytes.rstrip().decode('ascii')
            #print(name, list_type, lists)
            aecomp = self.add_aecomp(name, list_type, lists)
            str(aecomp)
        """
        Word Name Type Description
        1 NAME(2) CHAR4
        3 LABEL(2) CHAR4
        Words 3 through 4 repeat until (-1,-1) occurs
        """
        return len(data)

    def _read_aecompl(self, data: bytes, n: int) -> int:
        """
        MSC 2018.2

        Word Name Type Description
        1 NAME(2)  CHAR4
        3 LABEL(2) CHAR4
        Words 3 through 4 repeat until (-1,-1) occurs
        """
        struct1 = Struct(self._endian + b'8s')
        struct2 = Struct(self._endian + b'8s')
        struct_end = Struct(self._endian + b'2i')
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
                label = cast_string(label_bytes, self.size, encoding='latin1')
                labels.append(label)
                n += ntotal
                edata = data[n:n+ntotal]
            n += ntotal
            aecompl = self.add_aecompl(name, labels)
            str(aecompl)
        return len(data)

    def _read_spline1(self, data: bytes, n: int) -> int:
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
        ntotal = 48 # 4 * 12
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        #structi = Struct(self._endian + b'4i f 8s 8s 3i f') # msc
        structi = Struct(self._endian + b'5if 8s 8s 2i')
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
            self.add_spline1(eid, caero, box1, box2, setg,
                    dz=dz, method=method,
                    usage=usage, nelements=nelements,
                    melements=melements)
            n += ntotal
        #self.to_nx()
        return n

    def _read_spline2(self, data: bytes, n: int) -> int:
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
        ntotal = 48 * self.factor # 4 * 12
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        #assert ndatai % ntotal == 0
        if self.size == 4:
            structi = Struct(self._endian + b'5i 2f i 2f 8s')
        else:
            structi = Struct(self._endian + b'5q 2d q 2d 16s')

        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            eid, caero, id1, id2, setg, dz, dtor, cid, dthx, dthy, usage_bytes = out
            usage = usage_bytes.rstrip().decode('latin1')
            spline2 = self.add_spline2(
                eid, caero,
                id1, id2, setg,
                dz=dz, dtor=dtor, cid=cid,
                dthx=dthx, dthy=dthy,
                usage=usage)
            str(spline2)
            n += ntotal
        #n = self._read_spline2_nx(data, n)
        return n

    def _read_spline3(self, data: bytes, n: int) -> int:
        """reads the SPLINE3 card"""
        spline3
        #n = self._read_spline2_nx(data, n)
        return n
    def _read_spline4(self, data: bytes, n: int) -> int:
        """reads the SPLINE4 card"""
        n = self._read_spline4_nx(data, n)
        return n

    def _read_spline4_nx(self, data: bytes, n: int) -> int:
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
        12 FTYPE        I Radial interpolation funtion fir METHOD=RIS  (not in NX)
        13 RCORE       RS Radius of radial interpolation function      (not in NX)

        """
        # 792/4 = 198
        # 198 = 2 * 99 = 2 * 11 * 9
        ntotal = 4 * 11 # 4 * 13
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        #structi = Struct(self._endian + b'4i f 8s 8s 3i f') # msc
        structi = Struct(self._endian + b'4i f 8s 8s 2i')
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            #eid, caero, aelist, setg, dz, method_bytes, usage_bytes, nelements, melements, ftype, rcore = out  # msc
            eid, caero, aelist, setg, dz, method_bytes, usage_bytes, nelements, melements = out
            method = method_bytes.rstrip().decode('ascii')
            usage = usage_bytes.rstrip().decode('ascii')
            self.add_spline4(eid, caero, aelist, setg,
                             dz, method, usage,
                             nelements, melements)
            n += ntotal
        self.to_nx()
        return n

    def _read_spline5(self, data: bytes, n: int) -> int:
        """reads the SPLINE6 card"""
        n = self._read_spline5_nx(data, n)
        return n

    def _read_spline5_nx(self, data: bytes, n: int) -> int:
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
        16 FTYPE         I Radial interpolation funtion fir METHOD=RIS (not in NX?)
        17 RCORE        RS Radius of radial interpolation function     (not in NX?)

        """

        ntotal = 60 * self.factor # 4 * 12
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0, ndatai % ntotal
        if self.size == 4:
            structi = Struct(self._endian + b'4i 2f i 3f 8s8s f')
        else:
            asdf
            #structi = Struct(self._endian + b'5q 2d q 2d 16s')

        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            #ftype, rcore
            eid, caero, aelist, setg, dz, dtorxy, cid, dthx, dthy, dthz, usage_bytes, method_bytes, dtorzy = out
            method = cast_string(method_bytes, self.size, encoding='latin1')
            usage = cast_string(usage_bytes, self.size, encoding='latin1')
            #print(f'eid={eid} caero={caero} aelist={aelist} setg={setg} dz={dz} dtorxy={dtorxy} cid={cid} dthx={dthx} dthy={dthy} dthz={dthz} usage={usage!r} method={method!r} dtorzy={dtorzy}')
            assert method in ['IPS','TPS','FPS','RIS'], method
            assert usage in ['FORCE','DISP','BOTH'], usage
            thx = dthx
            thy = dthy
            dtor = dtorzy
            spline5 = self.add_spline5(
                eid, caero, aelist, setg, thx, thy,
                dz=dz, dtor=dtor, cid=cid,
                usage=usage, method=method,
                #ftype=ftype, rcore=rcore,  # not in NX
            )
            str(spline5)
            n += ntotal
        self.to_nx()
        return n

    def _read_monpnt1(self, data: bytes, n: int) -> int:
        """Reads the MONPNT1 card"""
        n = self._read_monpnt1_nx(data, n)
        return n

    def _read_monpnt1_nx(self, data: bytes, n: int) -> int:
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
        24 CD           I  (not in NX)

        """
        #ntotal = 4 * 24 # 4 * 24
        ntotal = 4 * 23 # 4 * 23
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        #structi = Struct(self._endian + b'8s 56s i 8s i 3f i')  # msc
        structi = Struct(self._endian + b'8s 56s i 8s i 3f')  # nx
        #monpnt1s = []
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            #name_bytes, label_bytes, axes, comp_bytes, cp, x, y, z, cd = out
            name_bytes, label_bytes, axes, aecomp_name_bytes, cp, x, y, z = out
            name = cast_string(name_bytes, self.size, encoding='latin1')
            label = cast_string(label_bytes, self.size, encoding='latin1')
            aecomp_name = cast_string(aecomp_name_bytes, self.size, encoding='latin1')
            xyz = [x, y, z]
            monpnt1 = self.add_monpnt1(name, label, axes, aecomp_name,
                                       xyz, cp=cp)
            str(monpnt1)
            n += ntotal
            #monpnt1s.append(monpnt1)
        self.to_nx()
        return n

    def _read_aestat(self, data: bytes, n: int) -> int:
        """
        MSC 2018.2

        Word Name Type Description
        1 ID I
        2 LABEL(2) CHAR4

        """
        ntotal = 12 * self.factor # 4 * 8
        if self.size == 4:
            structi = Struct(self._endian + b'i 8s')
        else:
            structi = Struct(self._endian + b'q 16s')
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0

        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            aestat_id, label_bytes = out
            label = cast_string(label_bytes, self.size, encoding='latin1')
            self.add_aestat(aestat_id, label)
            n += ntotal
        return n

    def _read_flutter(self, data: bytes, n: int) -> int:
        """
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
          10 EPR  RS
          11 SFLG  I SWEEP FLAG
        SFLG=1 (sweep)
          9 FMAX  RS maximum frequency
          10 EPR  RS
          11 SFLG  I SWEEP FLAG
        End SFLG
        Words 1 through max repeat until End of Record

        NX:
          data = (30, PK, 1, 2, 3, L, 3, 0.001, -1)
        """
        ints = np.frombuffer(data[n:], self.idtype).copy()
        floats = np.frombuffer(data[n:], self.fdtype).copy()
        iminus1 = np.where(ints == -1)[0]

        istart = [0] + list(iminus1[:-1] + 1)
        iend = iminus1

        for (i0, i1) in zip(istart, iend):
            sid = ints[i0]
            assert ints[i1] == -1, ints[i1]
            method_bytes = data[n+i0*4+4:n+i0*4+12]
            density = ints[i0+3]
            mach = ints[i0+4]
            reduced_freq_velocity = ints[i0+5]
            imethod_bytes = data[n+i0*4+24:n+i0*4+32]
            nvalue = ints[i0+8]
            epsilon = floats[i0+9]
            assert ints[i0+10] == -1, ints[i0:i1]

            method = method_bytes.rstrip().decode('ascii')
            imethod = imethod_bytes.rstrip().decode('ascii')
            self.add_flutter(sid, method,
                             density, mach, reduced_freq_velocity,
                             imethod=imethod, # 'L'
                             nvalue=nvalue,
                             epsilon=epsilon)
        self.to_nx()
        return len(data)
        #ntotal = 12 # 4 * 8
        #ndatai = len(data) - n
        #ncards = ndatai // ntotal
        #assert ndatai % ntotal == 0
        #structi = Struct(self._endian + b'i 8s')
        #for unused_i in range(ncards):
            #edata = data[n:n + ntotal]
            #out = structi.unpack(edata)
            #aestat_id, label = out
            #label = label.rstrip().decode('latin1')
            #self.add_aestat(aestat_id, label)
            #n += ntotal
        #return n

    def _read_trim(self, data: bytes, n: int) -> int:
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
        ntotal1 = 16 * self.factor # 4 * 4
        ntotal2 = 12 * self.factor # 4 * 3
        #ndatai = len(data) - n
        #ncards = ndatai // ntotal
        if self.size == 4:
            struct1 = Struct(self._endian + b'i 3f')
            struct2 = Struct(self._endian + b'8sf')
            struct_end = Struct(self._endian + b'3i')
        else:
            struct1 = Struct(self._endian + b'q 3d')
            struct2 = Struct(self._endian + b'16sd')
            struct_end = Struct(self._endian + b'3q')
        while n < len(data):
            edata = data[n:n+ntotal1]
            trim_id, mach, q, aeqr = struct1.unpack(edata)
            n += ntotal1
            labels = []
            uxs = []

            edata = data[n:n+ntotal2]
            while struct_end.unpack(edata) != (-1, -1, -1):
                label_bytes, ux = struct2.unpack(edata)
                label = cast_string(label_bytes, self.size, encoding='latin1')
                labels.append(label)
                uxs.append(ux)
                n += ntotal2
                edata = data[n:n+ntotal2]
            n += ntotal2
            trim = self.add_trim(trim_id, mach, q, labels, uxs, aeqr=aeqr, trim_type=1)
            str(trim)
        return n

    def _read_aesurf(self, data: bytes, n: int) -> int:
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
        10 CREFC      RS REFerence Chord Length for the CS > 0.0, default = 1.0
        11 CREFS      RS REFerence area for the CS > 0.0, default = 1.0
        12 PLLIM      RS Lower deflection   Limit for the control surface in radians, default=no limit
        13 PULIM      RS Upper deflection   Limit for the control surface in radians, default=no limit
        14 HMLLIM     RS Lower Hinge Moment Limit for the control surface in force-length units, default=no limit
        15 HMULIM     RS Upper Hinge Moment Limit for the control surface in force-length units, default=no limit
        16 TQLLIM     RS Lower deflection   Limit for the control surface as fct(q), >0, default=no limit
        17 TQULIM     RS Upper deflection   Limit for the control surface as fct(q), >0, default=no limit

        """
        ntotal = 68 *  self.factor # 4 * 17
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        if self.size == 4:
            struct1 = Struct(self._endian + b'i 8s 4i fi 8f')
            #struct2 = Struct(self._endian + b'i 8s 4i fi 4f4i')
        else:
            struct1 = Struct(self._endian + b'q 16s 4q dq 8d')
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out1 = struct1.unpack(edata)
            #out2 = struct2.unpack(edata)
            (aesurf_id, label_bytes, cid1, alid1, cid2, alid2, eff, ldw_int, crefc, crefs,
             pllim, pulim, hmllim, hmulim, tqllim, tqulim) = out1

            label = cast_string(label_bytes, self.size, encoding='latin1')

            if ldw_int == 0:
                ldw = 'LDW'
            elif ldw_int == 1:
                ldw = 'NOLDW'
            else:
                raise NotImplementedError(ldw_int)

            assert isinstance(ldw, str), ldw
            # TODO: incorrect...too strict
            # hmllim, hmulim, tqllim, tqulim

            # this is super janky
            if hmllim == 1.3563156426940112e-19:
                hmllim = None
            if hmulim == 1.3563156426940112e-19:
                hmulim = None
            if tqllim == 1.3563156426940112e-19:
                tqllim = None
            if tqulim == 1.3563156426940112e-19:
                tqulim = None
            #print(pllim, pulim, hmllim, hmulim, tqllim, tqulim)
            self.add_aesurf(aesurf_id, label, cid1, alid1, cid2=None, alid2=None,
                            eff=eff, ldw=ldw, crefc=crefc, crefs=crefs,
                            #pllim=-np.pi/2., pulim=np.pi/2.,
                            pllim=pllim, pulim=pulim,
                            hmllim=hmllim, hmulim=hmulim, # hinge moment lower/upper limits
                            tqllim=tqllim, tqulim=tqulim)
            n += ntotal
        return n

    def _read_aefact(self, data: bytes, n: int) -> int:
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
        #self.show_data(data[12:], types='if')
        ints = np.frombuffer(data[n:], self.idtype8).copy()
        floats = np.frombuffer(data[n:], self.fdtype8).copy()
        iminus1 = np.where(ints == -1)[0]

        istart = [0] + list(iminus1[:-1] + 1)
        iend = iminus1

        for (i0, i1) in zip(istart, iend):
            #if i0 == i1:
                #continue
            sid = ints[i0]
            fractions = floats[i0+1:i1]
            assert ints[i1] == -1, ints[i1]
            aefact = self.add_aefact(sid, fractions)
            #print(aefact)
            str(aefact)
            #n += ntotal
        return len(data)

def cast_string(name_bytes: bytes, size: int, encoding='latin1') -> str:
    if size == 4:
        name = name_bytes.decode(encoding).rstrip()
    else:
        name = reshape_bytes_block(name_bytes).decode(encoding).rstrip()
    return name


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


def _read_group_elem_prop_nids(ints, i, n, size) -> Tuple[int, int, Any]:
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
