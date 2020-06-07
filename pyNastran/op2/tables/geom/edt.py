"""
defines readers for BDF objects in the OP2 EDT/EDTS table
"""
from struct import Struct
import numpy as np
from pyNastran.op2.tables.geom.geom_common import GeomCommon


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
            (5201, 52, 373) : ['ACMODL', self._read_fake],
            (6301, 63, 397) : ['ADAPT', self._read_fake],
            (7801, 78, 582) : ['AECOMP', self._read_aecomp],
            (7901, 79, 583) : ['AECOMPL', self._read_aecompl],
            (7301, 73, 574) : ['AEDW', self._read_fake],
            (4002, 40, 273) : ['AEFACT', self._read_aefact],
            (7501, 75, 576) : ['AEFORCE', self._read_fake],
            (2602, 26, 386) : ['AELINK', self._read_aelink],
            (2302, 23, 341) : ['AELIST', self._read_aelist],
            (7001, 70, 571) : ['AEPARM', self._read_fake],
            (7401, 74, 575) : ['AEPRESS', self._read_fake],
            (3202, 32, 265) : ['AERO', self._read_fake],
            (2202, 22, 340) : ['AEROS', self._read_aeros],
            (2102, 21, 339) : ['AESTAT', self._read_aestat],
            (2002, 20, 338) : ['AESURF', self._read_aesurf],
            (7701, 77, 581) : ['AESURFS', self._read_fake],
            (3002, 30, 263) : ['CAERO1', self._read_caero1],
            (4301, 43, 167) : ['CAERO2', self._read_fake],
            (4401, 44, 168) : ['CAERO3', self._read_fake],
            (4501, 45, 169) : ['CAERO4', self._read_fake],
            (5001, 50, 175) : ['CAERO5', self._read_fake],
            (6201, 62, 143) : ['CLOAD', self._read_fake],
            (6401, 64, 307) : ['CSSCHD', self._read_fake],
            (104, 1, 81) : ['DEFORM', self._read_fake],
            (2702, 27, 387) : ['DIVERG', self._read_fake],
            (4102, 41, 274) : ['FLFACT', self._read_fake],
            (3902, 39, 272) : ['FLUTTER', self._read_fake],
            (17400, 174, 616) : ['GROUP', self._read_group],
            (3802, 38, 271) : ['MKAERO1', self._read_fake],
            (3702, 37, 270) : ['MKAERO2', self._read_fake],
            (7601, 76, 577) : ['MONPNT1', self._read_monpnt1],
            (3102, 31, 264) : ['PAERO1', self._read_paero1],
            (4601, 46, 170) : ['PAERO2', self._read_fake],
            (4701, 47, 171) : ['PAERO3', self._read_fake],
            (4801, 48, 172) : ['PAERO4', self._read_fake],
            (5101, 51, 176) : ['PAERO5', self._read_fake],
            (5301, 53, 378) : ['PANEL', self._read_fake],
            (3502, 35, 268) : ['SET1', self._read_set1],
            (3602, 36, 269) : ['SET2', self._read_fake],
            (4302, 43, 607) : ['SET3', self._read_fake],
            (3302, 33, 266) : ['SPLINE1', self._read_fake],
            (3402, 34, 267) : ['SPLINE2', self._read_fake],
            (4901, 49, 173) : ['SPLINE3', self._read_fake],
            (6501, 65, 308) : ['SPLINE4', self._read_spline4],
            (6601, 66, 309) : ['SPLINE5', self._read_fake],
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
        """
        nentries = 0
        ints = np.frombuffer(data[12:], dtype=self.idtype)
        strs = np.frombuffer(data[12:], dtype='|S4')
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
            grid, ndesc = ints[i:i+2]
            i += 2
            n += 8

            gdesc = ''.join(stri.decode('latin1') for stri in strs[i:i+ndesc]).strip()
            i += ndesc
            n += 4 * ndesc

            #------------------------------
            #gtype, nmeta, mdesc
            gtype = ints[i]
            i += 1
            n += 4
            print(f'grid={grid} ndesc={ndesc} gdesc={gdesc!r}; gtype={gtype!r}')

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

            while n < ndata:
                if gtype == -1:
                    # end of card
                    i += 1
                    n += 4
                    break
                elif gtype == -2:
                    # meta-data
                    nmeta = ints[i]
                    i += 1
                    n += 4
                    #print(i, nmeta)
                    #print(strs[i:i+nmeta-1])
                    mdesc = ''.join(stri.decode('latin1') for stri in strs[i:i+nmeta-1])
                    i += nmeta
                    n += 4 * nmeta
                    print(f'gtype={gtype} nmeta={nmeta} mdesc={mdesc!r}')
                    #iminus1 = minus1[minus1_count+2]
                    #print('ints: ', ints[i:iminus1].tolist())
                    #minus1_count += 1
                elif gtype == -4:
                    # grids
                    #iminus1 = minus1[minus1_count] # + 1
                    #print(ints[iminus1:])
                    #grids = ints[i+1:iminus1].tolist()

                    for j, nj in enumerate(ints[i:]):
                        if nj == -1:
                            break
                    grids = ints[i+1:i+j].tolist()
                    grids2 = _expand_vals(grids)
                    print(f'grids = {grids2}')
                    #minus1_count += 1

                    nstop = len(grids) + 2
                    i += nstop
                    n += nstop * 4
                    #i = iminus1
                    #n = iminus1 * 4
                elif gtype == -5:
                    #print('data', ints[i:].tolist())
                    #GTYPE = -5 Element identification numbers
                    #NDESC+5+NMETA
                    #ID I Element identification numbers
                    #> 0 for ID
                    #= 0 for THRU
                    #= -6 for BY
                    #= -7 for ALL
                    #Word NDESC+5+NMETA repeats until -1 occurs
                    for j, nj in enumerate(ints[i:]):
                        if nj == -1:
                            break
                    eids = ints[i+1:i+j].tolist()
                    eids2 = _expand_vals(eids)
                    print(f'eids = {eids2}')
                    nstop = len(eids) + 2
                    i += nstop
                    n += nstop * 4
                else:
                    raise NotImplementedError(gtype)
                gtype = ints[i]
                print(f'***gtype={gtype} (ndata-n)={(ndata-n)}')
                if gtype == -1 and (ndata - n) == 4:
                    print('break')
                    #minus1_count += 1
                    i += 1
                    n += 4
                    break
                #i += 1
                #n += 4
            #assert ints[i] == -1, ints[i:]
            self.log.warning(f'skipping GROUP in {self.table_name}')
            #self.add_rgyro(sid, asynci, refrot, unit, speed_low, speed_high, speed)
        self.increase_card_count('GROUP', nentries)
        return n

    def _read_aeros(self, data: bytes, n: int) -> int:
        """
        AEROS   0       100     36.     360.    12960.
        data = (0, 100, 36.0, 360.0, 12960.0, 0, 0)
        """
        assert len(data) == 40, len(data)
        struct = Struct(self._endian + b'2i 3f 2i')
        out = struct.unpack(data[12:])
        acsid, rcsid, cref, bref, sref, sym_xz, sym_xy = out
        self.add_aeros(cref, bref, sref, acsid=0, rcsid=0, sym_xz=0,
                       sym_xy=0, comment='aeros')
        n = 40
        return n

    def _read_caero1(self, data: bytes, n: int) -> int:
        """
        MSC 2018.2

        Word Name Type Description
        1 EID I
        2 PID I
        3 CP I
        4 NSPAN I
        5 NCHORD I
        6 LSPAN I
        7 LCHORD I
        8 IGID I
        9 X1 RS
        10 Y1 RS
        11 Z1 RS
        12 X12 RS
        13 X4 RS
        14 Y4 RS
        15 Z4 RS
        16 X43 RS

        CAERO1  100001  100001  0       10                      24      1
        99.2956521.45381-11.654442.85999101.8387122.6196-2.6930832.70996

        data = (100001, 100001, 0, 10, 0, 0, 24, 1,
                99.3, 21.45, -11.65, 42.86, 101.8387, 122.62, -2.69, 32.71)

        """
        ntotal = 4 * 16
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        structi = Struct(self._endian + b'8i 8f')
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

    def _read_paero1(self, data: bytes, n: int) -> int:
        """
        MSC 2018.2

        Word Name Type Description
        1 PID I
        2 B1 I
        3 B2 I
        4 B3 I
        5 B4 I
        6 B5 I
        7 B6 I
        8 UNDEF none

        PAERO1  100001

        data = (100001, 100001, 0, 10, 0, 0, 24, 1,
                99.3, 21.45, -11.65, 42.86, 101.8387, 122.62, -2.69, 32.71)

        """
        ntotal = 32 # 4 * 8
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        structi = Struct(self._endian + b'8i')
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            pid, b1, b2, b3, b4, b5, b6, empty = out
            caero_body_ids = []
            for body in [b1, b2, b3, b4, b5, b6, empty]:
                if body != 0:
                    caero_body_ids.append(body)
                    self.log.warning(out)
            self.add_paero1(pid, caero_body_ids=caero_body_ids, comment='')
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
        ints = np.frombuffer(data[n:], self.idtype).copy()
        floats = np.frombuffer(data[n:], self.fdtype).copy()
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
        MSC 2018.2

        Word Name Type Description
        1 SID I
        2 G1  I Grid ID or -2 when SKIN is specified
        Word 2 repeats until End of Record

        """
        ints = np.frombuffer(data[n:], self.idtype).copy()
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
            label = label_bytes.decode('latin1').rstrip()
            n += ntotal
            linking_coefficents = []
            independent_labels = []

            edata = data[n:n+ntotal]
            while struct_end.unpack(edata) != (-1, -1, -1):
                ind_label_bytes, coeff = struct2.unpack(edata)
                ind_label = ind_label_bytes.decode('latin1').rstrip()
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
            sid = ints[i0]
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
                label = label_bytes.decode('latin1').rstrip()
                labels.append(label)
                n += ntotal
                edata = data[n:n+ntotal]
            n += ntotal
            aecompl = self.add_aecompl(name, labels)
            str(aecompl)
        return len(data)

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
            name_bytes, label_bytes, axes, comp_bytes, cp, x, y, z = out
            name = name_bytes.rstrip().decode('ascii')
            label = label_bytes.rstrip().decode('ascii')
            aecomp_name = comp_bytes.rstrip().decode('ascii')
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
        ntotal = 12 # 4 * 8
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        structi = Struct(self._endian + b'i 8s')
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = structi.unpack(edata)
            aestat_id, label = out
            label = label.rstrip().decode('ascii')
            self.add_aestat(aestat_id, label)
            n += ntotal
        return n

    def _read_trim(self, data: bytes, n: int) -> int:
        """
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
        ntotal1 = 16 # 4 * 4
        ntotal2 = 12 # 4 * 3
        #ndatai = len(data) - n
        #ncards = ndatai // ntotal
        struct1 = Struct(self._endian + b'i 3f')
        struct2 =Struct(self._endian + b'8sf')
        struct_end = Struct(self._endian + b'3i')
        while n < len(data):
            edata = data[n:n+ntotal1]
            trim_id, mach, q, aeqr = struct1.unpack(edata)
            n += ntotal1
            labels = []
            uxs = []

            edata = data[n:n+ntotal2]
            while struct_end.unpack(edata) != (-1, -1, -1):
                label_bytes, ux = struct2.unpack(edata)
                label = label_bytes.decode('latin1').rstrip()
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
        ntotal = 68 # 4 * 17
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        struct1 = Struct(self._endian + b'i 8s 4i fi 8f')
        struct2 = Struct(self._endian + b'i 8s 4i fi 4f4i')
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out1 = struct1.unpack(edata)
            out2 = struct2.unpack(edata)
            (aesurf_id, label_bytes, cid1, alid1, cid2, alid2, eff, ldw_int, crefc, crefs,
             pllim, pulim, hmllim, hmulim, tqllim, tqulim) = out1
            label = label_bytes.decode('ascii').rstrip()
            if ldw_int == 0:
                ldw = 'LDW'
            elif ldw_int == 1:
                ldw = 'NOLDW'
            else:
                raise NotImplementedError(ldw_int)

            assert isinstance(ldw, str), ldw
            # TODO: incorrect...too strict
            # hmllim, hmulim, tqllim, tqulim

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
        ints = np.frombuffer(data[n:], self.idtype).copy()
        floats = np.frombuffer(data[n:], self.fdtype).copy()
        iminus1 = np.where(ints == -1)[0]

        istart = [0] + list(iminus1[:-1] + 1)
        iend = iminus1

        for (i0, i1) in zip(istart, iend):
            sid = ints[i0]
            fractions = floats[i0+1:i1]
            assert ints[i1] == -1, ints[i1]
            self.add_aefact(sid, fractions)
            #n += ntotal
        return len(data)
