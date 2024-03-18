from __future__ import annotations
from struct import Struct
from typing import TYPE_CHECKING

import numpy as np

from pyNastran.op2.op2_interface.utils import reshape_bytes_block

#from pyNastran.op2.op2_interface.read_extdb import (
    #read_extdb, read_tug1, read_mef1)
#from pyNastran.op2.op2_interface.read_external_superelement import read_cmodext

from pyNastran.op2.op2_interface.utils_matpool import (
    #read_matpool_dmig,
    read_matpool_dmig_4, read_matpool_dmig_8)

#from pyNastran.op2.op2_interface.read_matrix import (
    #read_matrix_mat)
if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger
    from pyNastran.op2.op2_interface.op2_reader import OP2Reader
    from pyNastran.op2.op2 import OP2


def read_matrix_matpool(op2_reader: OP2Reader) -> None:
    """
    Reads a MATPOOL matrix

    MATPOOL matrices are always sparse

    +------+-----------------+
    | Form | Meaning         |
    +======+=================+
    |  1   | Square          |
    |  2   | Rectangular     |
    |  6   | Symmetric       |
    |  9   | Pseudo identity |
    +------+-----------------+


    Record 2  - BNDFL(9614,96,0)
    Record 3  - DMIAX(214,2,221)
    Record 4  - DMIG(114,1,120)
    Record 5  - DMIJ(514,5,578)
    Record 6  - DMIJI(614,6,579)
    Record 7  - DMIK(714,7,580)
    Record 8  - ELIST(314,3,279)
    Record 9  - MFLUID(414,4,284)
    Record 10 - RADCAV(2509,25,418)
    Record 11 - RADSET(8602,86,421)
    Record 12 - RADLST(2014,20,243)
    Record 13 - RADMTX(3014,30,244)
    """
    #print('-------------------------------------')
    op2: OP2 = op2_reader.op2
    read_matpool = op2.read_matpool

    table_name = op2_reader._read_table_name(rewind=False, stop_on_failure=True)
    utable_name = table_name.decode('utf-8')
    #print(utable_name)
    op2_reader.read_markers([-1])

    # (104, 32768, 0, 0, 0, 0, 0)
    data = op2_reader._read_record()
    #op2_reader.show_data(data, types='ifs', endian=None, force=False)

    op2_reader.read_3_markers([-2, 1, 0])

    # MATPOOL
    data = op2_reader._read_record()
    #op2_reader.show_data(data, types='ifs', endian=None, force=False)
    if op2_reader.size == 8:
        data = reshape_bytes_block(data)
    #op2_reader.show_data(data)
    ndata = len(data)
    if ndata == 8:
        table_name2, = op2.struct_8s.unpack(data)
        utable_name = table_name2.decode('utf-8').strip()
        #assert utable_name == utable_name2, f'utable_name={utable_name} utable_name2={utable_name2}'

    op2_reader.read_markers([-3, 1, 0])
    data = op2_reader._read_record()

    if op2_reader.size == 4:
        struct_3i = op2.struct_3i
    else:
        struct_3i = op2.struct_3q
    itable = -4
    while 1:
        #op2_reader.show_data(data[:12*op2_reader.factor], types='ifsqd')
        code = struct_3i.unpack(data[:12*op2_reader.factor])
        #op2_reader.show(36)

        #if utable_name == 'DELTAK':
            #op2_reader.show_data(data)

        #nvalues = len(data) // 4
        assert len(data) % 4 == 0, len(data) / 4

        if op2.read_mode == 1 and read_matpool:
            read_matpool_result(op2_reader, code, op2, data, utable_name)

        op2_reader.read_3_markers([itable, 1, 0])
        #op2_reader.log.debug(f'  read [{itable},1,0]')
        expected_marker = itable - 1
        data, ndatas = op2_reader.read_long_block(expected_marker)
        if ndatas == 0:
            itable -= 1
            #op2_reader.log.debug(f'  read [{itable},1,0]')
            op2_reader.read_3_markers([itable, 1, 0])
            break
        elif op2.read_mode == 1 and read_matpool:
            #op2_reader.show_data(data, types='ifs', endian=None, force=False)
            code = struct_3i.unpack(data[:12*op2_reader.factor])
            read_matpool_result(op2_reader, code, op2, data, utable_name)
            #op2_reader.log.info('showing data...')
            #data, ndatas = op2_reader.read_long_block(expected_marker)
        itable -= 1

    #op2_reader.show(100, types='ifsqd')
    op2_reader.read_markers([0])
    #raise RuntimeError('failed on _read_matpool_matrix')

def read_matpool_result(op2_reader: OP2Reader,
                        code: tuple[int, int, int],
                        op2: OP2,
                        data: bytes, utable_name: str) -> None:
    if code == (114, 1, 120):
        op2_reader.log.debug(f'  code = {code}')
        #raise NotImplementedError('read_matpool_dmig')
        try:
            if op2_reader.size == 4:
                read_matpool_dmig_4(op2, data, utable_name, debug=False)
            else:
                read_matpool_dmig_8(op2, data, utable_name, debug=False)
        except Exception as excep:
            op2_reader.log.error(str(excep))
            op2_reader.log.warning('  skipping MATPOOL-DMIG')
            #raise
    elif code == (314, 3, 279):
        # geom
        read_matpool_elist(op2_reader, op2, data, utable_name, debug=False)
    elif code == (414, 4, 284):
        # geom
        read_matpool_mfluid(op2_reader, op2, data, utable_name, debug=False)
    elif code == (2509, 25, 418):
        # C:\NASA\m4\formats\git\examples\pyNastran_examples\demo_sort2_post_m2\hd15305.op2
        op2_reader.log.warning('  skipping MATPOOL-RADCAV')
    elif code == (3014, 30, 244):
        # C:\NASA\m4\formats\git\examples\pyNastran_examples\demo_sort2_post_m2\hd15305.op2
        op2_reader.log.warning('  skipping MATPOOL-RADMTX')
    elif code == (8602, 86, 421):
        # C:\NASA\m4\formats\git\examples\pyNastran_examples\demo_sort2_post_m2\hd15305.op2
        op2_reader.log.warning('  skipping MATPOOL-RADSET')
    elif code == (2014, 20, 243):
        # C:\NASA\m4\formats\git\examples\pyNastran_examples\demo_sort2_post_m2\hd15306.op2
        op2_reader.log.warning('  skipping MATPOOL-RADLST')
    elif code == (9614, 96, 0):
        # some axisymmetric matrix
        read_matpool_bndfl(op2_reader, op2, data, utable_name, debug=False)
    else:
        print(f'  code = {code}')
        op2_reader.show_data(data, types='ifs', endian=None, force=False)
        raise NotImplementedError(code)

def read_matpool_mfluid(op2_reader: OP2Reader,
                        op2: OP2,
                        data: bytes, utable_name: str,
                        debug: bool=False) -> None:
    r"""
    Word Name Type Description
    1 SID       I
    2 CID       I
    3 ZFR      RS
    4 RHO      RS
    5 ELIST1    I
    6 ELIST2    I
    7 PLANE1    I
    8 PLANE2    I
    9 RMAX     RS
    10 FMEXACT RS

    C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\cms01.op2
    """
    endian = op2_reader._endian
    assert len(data) == 12 + 40, len(data)
    structi = Struct(endian + b'2i 2f 4i 2f')
    (unused_sid, unused_cid, unused_zfr, unused_rho, unused_elist1, unused_elist2,
     unused_plane1, unused_plane2, unused_rmax, unused_fmexact) = structi.unpack(data[12:])
    op2_reader.log.warning('skipping MATPOOL-MFLUID table')

def read_matpool_elist(op2_reader: OP2Reader,
                       op2: OP2,
                       data: bytes, utable_name: str,
                       debug: bool=False) -> None:
    """
    Word Name Type Description
    1 LID I
    2 E1 I
    Word 2 repeats until End of Record
    """
    n = 12 * op2_reader.factor
    datai = data[n:]
    ints = np.frombuffer(datai, dtype=op2.idtype8).copy()
    assert ints[-1] == 0, ints
    elist_id = ints[0]
    element_ids = ints[1:-1].tolist()
    if not hasattr(op2_reader, 'elist'):
        op2_reader.elist = {}
    op2_reader.elist[elist_id] = element_ids

def read_matpool_bndfl(op2_reader: OP2Reader,
                       op2: OP2,
                       data: bytes, utable_name: str,
                       debug: bool=False) -> None:
    r"""
    Word Name Type Description
    1 CSF      I
    2 G        RS
    3 RHO      RS
    4 B        RS
    5 NOSYM    I
    6 M        I
    7 S1       I
    8 S2       I
    9 NHARM(C) I
    10 NI      I
    Word 10 repeats NHARM times
    11 IDFL I
    12 R    RS
    13 Z    RS
    14 L    RS
    15 C    RS
    16 S    RS
    17 RHOI RS
    18 G    I
    19 PHI  RS
    Words 18 through 19 repeat until (-1,-1) occurs
    Words 11 through 19 repeat until End of Record

    ints = [7, 26, 27, 43, 44, 60, 61]
    C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1072x.op2
    """
    endian = op2_reader._endian
    struct1 = Struct(endian + b'i 3f 5i')
    #          c   g    rho   b   nosym m s1  s2 nharm n1, n2
    #ints    = (2, 32.2, 0.03, 0.0, 0,   4, 1, -1, 2,    4, 8,
    #           # idfl r    z     l    c    s    rho
    #           2,     8.0, 10.0, 2.5, 1.0, 0,   0.03,
    #           # g  phi
    #           3, 0,
    #           4, 30.0,
    #           5, 60.0,
    #           6, 90.0,
    #           -1, -1,
    #
    #           # idfl  r    z    l      c     s    rho
    #           8,      8.0, 5.0, 5.0, 0.93, -0.35, 0.03,
    #           # g phi
    #           9, 0,
    #           10, 30.0,
    #           11, 60.0,
    #           12, 90.0,
    #           -1, -1,
    #           ...)

    size = op2_reader.size
    factor = op2_reader.factor
    log: SimpleLogger = op2_reader.log
    ndata = len(data)
    log.warning('skipping MATPOOL-BNDFL table')
    n = 12 * factor
    datai = data[n:]
    ints = np.frombuffer(datai, dtype=op2.idtype8).copy()
    floats = np.frombuffer(datai, dtype=op2.fdtype8).copy()
    iminus1 = np.where(ints == -1)[0]
    #print(iminus1)

    iharm = 9
    nharm = ints[iharm]
    #idfl = ints[iharm + nharm]
    #print(f'nharm = {nharm}')

    #b_idfl = (iharm + nharm) * size + 12
    #op2_reader.show_data(datai, types='if')
    #op2_reader.show_data(datai[:b_idfl], types='ifs')

    ntotal1 = 9 * size
    #print('i =', i)

    edata = data[n:n+ntotal1]
    #print('len(edata)', len(edata))
    c, g, rho, b, nosym, m, s1, s2, nharm = struct1.unpack(edata)
    log.debug(f'  c={c} g={g:g} rho={rho:g} nosym={nosym} m={m} s1={s1} s2={s2} nharm={nharm}')
    n += ntotal1
    i = 9

    ntotal2 = nharm * size
    #ni = Struct(b'%di' % nharm).unpack(data[n:n+ntotal2])
    ni = ints[i:i+nharm]
    #print(f'n={n} i={i} ni={ni}')
    n += ntotal2
    i += nharm

    while n < ndata:
        #print('------------------------------------')
        #print(ints[i-2:i+10])
        #print(floats[i-2:i+10])
        # 11 IDFL I
        # 12 R    RS
        # 13 Z    RS
        # 14 L    RS
        # 15 C    RS
        # 16 S    RS
        # 17 RHOI RS
        idfl = ints[i]
        r, z, l, c, s, rhoi = floats[i+1:i+7]
        log.debug(f'    idfl={idfl} r={r} z={z} L={l:g} c={c:g} s={s:g} rhoi={rhoi:g}')
        i += 7
        n += 7 * size
        assert idfl < 1000

        ints2 = ints[i+1::2]
        ints1 = ints[i::2][:len(ints2)]
        iminus1i = np.where((ints1 == -1) & (ints2 == -1))[0][0] - 1

        #intsi = ints[i:i+iminus1i+4:2]  # good
        intsi = ints1[:iminus1i+1]
        #floatsi = floats[i+1::2][:iminus1i+1]  #good
        floatsi = floats[i+1:i+2*(iminus1i+1):2]
        log.debug(f'    {intsi} {floatsi}')
        assert intsi.min() >= 1, intsi
        assert intsi.max() <= 1000, intsi
        #n += (iminus1i + 4) * size# ???
        #i += iminus1i + 1
        i += len(intsi) * 2 + 2 # ???
        n2 = (i + 3) * size
        n = n2
        #print(f'n={n} i={i} -> n2={n2}')
    return

