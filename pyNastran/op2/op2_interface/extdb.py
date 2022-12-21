from __future__ import annotations
from struct import Struct
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger
    from pyNastran.op2.op2 import OP2


def _read_extdb_extdb(self, name_str: str, data: bytes, endian: bytes, idtype: str) -> None:
    #TODO: needs work...
    # ints    = (0, 6, 2, 2, 0, 0, 1018, 1, 1, 1, 1, 725010254, 1099302303, -1, -1)
    # strings = (b'r\x00\x00\x00\x01\x00\x00\x00x\x00\x00\x00EXTDB   \x00\x00\x00\x00\x06\x00\x00\x00\x02\x00\x00\x00\x02\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xfa\x03\x00\x00\x01\x00\x00\x00\x01\x00\x00\x00\x01\x00\x00\x00\x01\x00\x00\x00N\xc76+\x9f\x05\x86A\xff\xff\xff\xff\xff\xff\xff\xff\x01\x00\x00\x00\x02\x00\x00\x00\x01\x00\x00\x00\x01\x00\x00\x00\xee\x92\x91Z',)
    # ints    = (114, 1, 120, EXTDB, 0, 6, 2, 2, 0, 0, 1018, 1, 1, 1, 1, 725010254, 1099302303, -1, -1, 1, 2, 1, 1, 1519489774)
    # floats  = (114, 1, 120, EXTDB, 0.0, 6, 2, 2, 0.0, 0.0, 1.4265218366826638e-42, 1, 1, 1, 1, 6.493597977039189e-13, 16.752744674682617, nan, nan, 1, 2, 1, 1, 2.048771126145843e+16)
    #self.log.debug('start _read_extdb_extdb')
    size = self.size
    factor = self.factor
    if size == 4:
        struct2i = Struct(b'2i')
        struct2id = Struct(b'2i d')
        struct7i = Struct(b'7i')
        ntotal2 = 16
    else:
        struct2i = Struct(b'2q')
        struct2id = Struct(b'2q d')
        struct7i = Struct(b'7q')
        ntotal2 = 24

    name_bytes = data[12*factor:20*factor]
    unused_name = name_bytes.decode('latin1').rstrip()
    #print(name)
    #self.show_data(data[20:4000], types='if', force=True)
    ints = np.frombuffer(data[20*factor:], dtype=idtype)
    iminus1 = np.where(ints == -1)[0]
    i2 = np.where(iminus1[:-1] + 1 == iminus1[1:])[0]
    #print(iminus1.tolist())
    #print(i2)
    istart = [0] + (iminus1[i2] + 1).tolist()# [:-1]
    iend = iminus1[i2]
    #print(istart, len(istart))
    #print(iend, len(iend))

    n = 20 * factor
    for i0, i1 in zip(istart, iend):
        edata = data[n+i0*size:n+i1*size]
        if len(edata) == 0:
            break

        n2 = 0
        flag, niii = Struct(mapfmt(self._endian + b'ii', size)).unpack(edata[:8*factor])

        #print(n, i0, i1)
        # inputs = (0, 6, 2, 2, 0, 0,
        #           1018, 1, 1, 1, 1, 725010254, 1099302303, -1, -1
        # ints = (  0, 9, 2, 2, 0, 0, 1, 1, 1, # ???
        #          1, 1, 0, 1.875,
        #          1, 2, 0, 1.875,
        #          1, 3, 0, 1.875,
        #          1, 4, 0, 1.875,
        #          1, 5, 0, 1.875,
        #          1, 6, 0, 1.875,

        #          2, 1, 0, 1.875,
        #          2, 2, 0, 1.875,
        #          2, 3, 0, 1.875,
        #          2, 4, 0, 1.875,
        #          2, 5, 0, 1.875,
        #          2, 6, 0, 1.875,

        #          3, 1, 0, 1.875,
        #          3, 2, 0, 1.875,
        #          3, 3, 0, 1.875,
        #          3, 4, 0, 1.875,
        #          3, 5, 0, 1.875,
        #          3, 6, 0, 1.875,

        #          4, 1, 0, 1.875,
        #          4, 2, 0, 1.875,
        #          4, 3, 0, 1.875,
        #          4, 4, 0, 1.875,
        #          4, 5, 0, 1.875,
        #          4, 6, 0, 1.875,
        #          100001, 0, 0, 1.875,
        #          100002, 0, 0, 1.875,
        #          100003, 0, 0, 1.875,
        #          100004, 0, 0, 1.875)
        #self.show_data(data, types='ifqsd', force=True)
        #print(len(edata))

        # 0, 6, 2, 2, 0, 0,
        # 1018...

        # DMIG    KAAX           0       6       2       0                     243

        #(1, 2,
         #1, 1, 37772923, 1094840434,
         #1, 2, 2351332.017589388, 9779310.02896148)

        if flag == 0:  #  was counter
            out = struct7i.unpack(edata[:7*size])
            zero_a, ifo, tin, tout, polar, zero_b, ncols = out
            #self.log.info(f'flag0: ifo={ifo} tin={tin} tout={tout} polar={polar} ncols={ncols}')
            #, zero_a, ifo, tin, tout, polar, zero_b, ncols
            assert zero_a == 0, out
            assert zero_b == 0, out
            assert polar == 0, polar
            assert ifo in [6, 9], ifo
            #ifo : int
            #    matrix shape
            #    4=Lower Triangular
            #    5=Upper Triangular
            #    6=Symmetric
            #    8=Identity (m=nRows, n=m)

            #tin : int
            #    matrix input precision
            #    1=Real, Single Precision
            #    2=Real, Double Precision
            #    3=Complex, Single Precision
            #    4=Complex, Double Precision
            #tout : int
            #    matrix output precision
            #    0=same as tin
            #    1=Real, Single Precision
            #    2=Real, Double Precision
            #    3=Complex, Single Precision
            #    4=Complex, Double Precision
            tin_tout = tin, tout
            # TODO: this should only be tin or tout (probably tout), but until we have validation...
            if tin_tout == (2, 2):
                n2 = 7 * size
            elif tin_tout == (1, 1):
                n2 = 7 * size
            else:
                raise NotImplementedError(tin_tout)
        else:
            #self.log.debug(f'flag = {flag}')
            assert flag == -1, flag
            n2 = size
            #edata = edata[7*size:]

        len_edata = len(edata) - n2

        nvalues = len_edata // size  # assume doubles are 2

        if tin_tout == (2, 2):
            # 2 = (gi, ci)
            # 4 = (gj, cj, real_a, real_b) = (gj, cj, real)
            ngrid_components = (nvalues - 2) // 4  # 4 values
            assert (nvalues - 2) % 4 == 0
            #int_type = 'int32'
            #float_type = 'float64'
            #dint = 1
            #dfloat = 2
            #float_factor = 2
        elif tin_tout == (1, 1):
            # 2 = (gi, ci)
            # 3 = (gj, cj, real)
            ngrid_components = (nvalues - 2) // 3  # 3 values
            #int_type = 'int32'
            #float_type = 'float32'
            #dint = 1
            #dfloat = 1
            #float_factor = 1
        else:
            raise RuntimeError(tin_tout)
        #float_factor = dfloat // dint
        #print(float_type)

        GCj = []
        GCi = []
        reals = []
        edatai = edata[n2:n2+2*size]
        n2 += 2*size

        #1, 1,
        gj, cj = struct2i.unpack(edatai)
        gcj = gj, cj
        for ii in range(ngrid_components):
            # (1, 1, 11933743.780688906)
            edata2 = edata[n2:n2+ntotal2]
            # (1, 1, 0, 1.875)
            out = struct2id.unpack(edata2)
            #self.log.debug(' 16: ' + str(out))
            gi, ci, real = out
            gci = gi, ci
            GCj.append(gcj)
            GCi.append(gci)
            reals.append(real)
            n2 += ntotal2
            if ii > 0:
                #self.log.info(f'  breaking n2={n2}')
                break

        dmig = self.op2.add_dmig(name_str, ifo, tin, tout, polar, ncols, GCj, GCi,
                                 reals, Complex=None, comment='')
        if flag == 0:
            print(dmig)
        #print(dmig)
        str(dmig)
        #n += (i1 - i0 + 1) * 4
    #self.log.debug('end _read_extdb_extdb')
    return

def _read_extdb_geomx(self, data: bytes, endian: bytes, function_map):
    """this is literally a GEOMx table, but embedded in a EXTRN table"""
    size = self.size

    struct3i = Struct(mapfmt(endian + b'3i', size))

    factor = self.factor
    n = 12 * factor
    code = struct3i.unpack(data[:n])
    if code in function_map:
        name, func = function_map[code]
        self.log.debug(f'code = {code} -> {name}')
        func(data, n)
    else:
        self.log.error(f'code = {code}')
        raise RuntimeError(code)
    return


def _read_extdb_phip(self, marker: int,
                     data: bytes, endian: bytes,
                     int_type: str, idtype: str):
    #   ints    = (1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 9, 0,
    #              10, 0, 11, 0, 12, 0, 13, 0, 14, 0, 15, 0, 16, 0, 17, 0, 18, 0, 19, 0,
    #              20, 0, 21, 0, 22, 0, 23, 0, 24, 0, 25, 0, 26, 0, 27, 0, 28, 0)
    factor = self.factor
    if marker == -3:
        self.log.warning('showing for marker=-3')
        self.show_data(data, types=int_type)
    elif marker == -4:
        #data = (
        #    'TYPE  IDCOMP ROW    TYPE  IDCOMP ROW',
        #    1, 1, 6, 1, 0,
        #    1, 2, 6, 7, 0,
        # ...
        #    1, 168, 6, 1003, 0,
        #    1, 169, 6, 1009, 0,
        #    2, 100001, 6, 1015, 0,
        #    2, 100002, 6, 1021, 0,
        #    2, 100003, 6, 1027, 0,
        #    2, 100004, 6, 1033, 0)
        word = data[:40*factor].decode('latin1')
        print(f'word = {word!r}')
        self.log.info(f'itable = {marker}')

        ints = np.frombuffer(data[40*factor:], dtype=idtype)
        nints = len(ints)
        ints2 = ints.reshape(nints//5, 5).copy()
        print(ints2)
    elif marker == -5:
        # 33      1       -2.500002.500000
        # 33      2       -2.500002.500000
        # 33      3       -2.500002.500000
        # 33      4       -2.500002.500000
        # 33      5       -2.500002.500000
        # 33      6       -2.500002.500000
        # 33      7       -2.500002.500000
        # 33      8       -2.500002.500000
        # 33      9       -2.500002.500000
        # 33      10      -2.500002.500000
        # 33      11      -2.500002.500000
        # 33      12      -2.500002.500000
        # 33      13      -2.500002.500000
        # 33      14      -2.500002.500000
        # 33      15      -2.500002.500000
        # 33      16      -2.500002.500000
        # 33      17      -2.500002.500000
        # 33      18      -2.500002.500000
        # 33      19      -2.500002.500000
        # 33      20      -2.500002.500000
        # 33      21      -2.500002.500000
        # 33      22      -2.500002.500000
        # 33      23      -2.500002.500000
        # 33      24      -2.500002.500000
        # 33      25      -2.500002.500000
        # 33      26      -2.500002.500000
        # 33      27      -2.500002.500000
        # 33      28      -2.500002.500000
        # 33      29      -2.500002.500000
        # 33      30      -2.500002.500000
        # 33      31      -2.500002.500000
        # 33      32      -2.500002.500000
        # 33      33      -2.500002.500000
        # 33      34      -2.500002.500000
        # 33      35      -2.500002.500000
        # 33      36      -2.500002.500000
        # 33      37      -2.500002.500000
        # 33      38      -2.500002.500000
        # 33      39      -2.500002.500000
        # 33      40      -2.500002.500000
        # 33      41      -2.500002.500000
        # 33      42      -2.500002.500000
        # 33      43      -2.500002.500000
        # 33      44      -2.500002.500000
        # 33      45      -2.500002.500000
        # 33      46      -2.500002.500000
        # 33      47      -2.500002.500000
        # 33      48      -2.500002.500000
        # 33      49      -2.500002.500000
        # 33      50      -2.500002.500000
        # 33      51      -2.500002.500000
        # 33      52      -2.500002.500000
        # 33      53      -2.500002.500000
        # 33      54      -2.500002.500000
        # 33      55      -2.500002.500000
        # 33      56      -2.500002.500000
        # 33      57      -2.500002.500000
        # 33      58      -2.500002.500000
        # 33      59      -2.500002.500000
        # 33      60      -2.500002.500000
        # 33      61      -2.500002.500000
        # 33      62      -2.500002.500000
        # 33      63      -2.500002.500000
        # 33      64      -2.500002.500000
        # 33      65      -2.500002.500000
        # 33      66      -2.500002.500000
        # 33      67      -2.500002.500000
        # 33      68      -2.500002.500000
        # 33      69      -2.500002.500000
        # 33      70      -2.500002.500000
        # 33      71      -2.500002.500000
        # 33      72      -2.500002.500000
        # 33      73      -2.500002.500000
        # 33      74      -2.500002.500000
        # 33      75      -2.500002.500000
        # 33      76      -2.500002.500000
        # 33      77      -2.500002.500000
        # 33      78      -2.500002.500000
        # 33      79      -2.500002.500000
        # 33      80      -2.500002.500000
        # 33      81      -2.500002.500000
        # 33      82      -2.500002.500000
        # 33      83      -2.500002.500000
        # 33      84      -2.500002.500000
        # 33      85      -2.500002.500000
        # 33      86      -2.500002.500000
        # 33      87      -2.500002.500000
        # 33      88      -2.500002.500000
        # 33      89      -2.500002.500000
        # 33      90      -2.500002.500000
        # 33      91      -2.500002.500000
        # 33      92      -2.500002.500000
        # 33      93      -2.500002.500000
        # 33      94      -2.500002.500000
        # 33      95      -2.500002.500000
        # 33      96      -2.500002.500000
        # 33      97      -2.500002.500000
        # 33      98      -2.500002.500000
        # 33      99      -2.500002.500000
        # 33      100     -2.500002.500000
        # 33      101     -2.500002.500000
        # 33      102     -2.500002.500000
        # 33      103     -2.500002.500000
        # 33      104     -2.500002.500000
        # 33      105     -2.500002.500000
        # 33      106     -2.500002.500000
        # 33      107     -2.500002.500000
        # 33      108     -2.500002.500000
        # 33      109     -2.500002.500000
        # 33      110     -2.500002.500000
        # 33      111     -2.500002.500000
        # 33      112     -2.500002.500000
        # 33      113     -2.500002.500000
        # 33      114     -2.500002.500000
        # 33      115     -2.500002.500000
        # 33      116     -2.500002.500000
        # 33      117     -2.500002.500000
        # 33      118     -2.500002.500000
        # 33      119     -2.500002.500000
        # 33      120     -2.500002.500000
        # 33      121     -2.500002.500000
        # 33      122     -2.500002.500000
        # 33      123     -2.500002.500000
        # 33      124     -2.500002.500000
        # 33      125     -2.500002.500000
        # 33      126     -2.500002.500000
        # 33      127     -2.500002.500000
        # 33      128     -2.500002.500000
        # 33      129     -2.500002.500000
        # 33      130     -2.500002.500000
        # 33      131     -2.500002.500000
        # 33      132     -2.500002.500000
        # 33      133     -2.500002.500000
        # 33      134     -2.500002.500000
        # 33      135     -2.500002.500000
        # 33      136     -2.500002.500000
        # 33      137     -2.500002.500000
        # 33      138     -2.500002.500000
        # 33      139     -2.500002.500000
        # 33      140     -2.500002.500000
        # 33      141     -2.500002.500000
        # 33      142     -2.500002.500000
        # 33      143     -2.500002.500000
        # 33      144     -2.500002.500000
        ndata = len(data)
        ntotali = 32
        nrows = ndata // ntotali
        #self.show_data(data, types='s', force=True)
        for irow in range(nrows):
            out = data[irow*ntotali:(irow+1)*ntotali]
            if irow > 5:
                self.log.info(f'  breaking irow={irow}')
                break
            print(out)
    elif marker == -5:
        self.log.warning('showing for marker=-5')
        self.show_data(data, types='ifs', force=True)
        aaa
    return
