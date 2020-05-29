"""
defines readers for BDF objects in the OP2 EDT/EDTS table
"""
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
            (7801, 78, 582) : ['AECOMP', self._read_fake],
            (7901, 79, 583) : ['AECOMPL', self._read_fake],
            (7301, 73, 574) : ['AEDW', self._read_fake],
            (4002, 40, 273) : ['AEFACT', self._read_fake],
            (7501, 75, 576) : ['AEFORCE', self._read_fake],
            (2602, 26, 386) : ['AELINK', self._read_fake],
            (2302, 23, 341) : ['AELIST', self._read_fake],
            (7001, 70, 571) : ['AEPARM', self._read_fake],
            (7401, 74, 575) : ['AEPRESS', self._read_fake],
            (3202, 32, 265) : ['AERO', self._read_fake],
            (2202, 22, 340) : ['AEROS', self._read_fake],
            (2102, 21, 339) : ['AESTAT', self._read_fake],
            (2002, 20, 338) : ['AESURF', self._read_fake],
            (7701, 77, 581) : ['AESURFS', self._read_fake],
            (3002, 30, 263) : ['CAERO1', self._read_fake],
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
            (7601, 76, 577) : ['MONPNT1', self._read_fake],
            (3102, 31, 264) : ['PAERO1', self._read_fake],
            (4601, 46, 170) : ['PAERO2', self._read_fake],
            (4701, 47, 171) : ['PAERO3', self._read_fake],
            (4801, 48, 172) : ['PAERO4', self._read_fake],
            (5101, 51, 176) : ['PAERO5', self._read_fake],
            (5301, 53, 378) : ['PANEL', self._read_fake],
            (3502, 35, 268) : ['SET1', self._read_fake],
            (3602, 36, 269) : ['SET2', self._read_fake],
            (4302, 43, 607) : ['SET3', self._read_fake],
            (3302, 33, 266) : ['SPLINE1', self._read_fake],
            (3402, 34, 267) : ['SPLINE2', self._read_fake],
            (4901, 49, 173) : ['SPLINE3', self._read_fake],
            (6501, 65, 308) : ['SPLINE4', self._read_fake],
            (6601, 66, 309) : ['SPLINE5', self._read_fake],
            (2402, 24, 342) : ['TRIM', self._read_fake],
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

