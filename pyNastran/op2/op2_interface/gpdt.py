from __future__ import annotations
from struct import Struct
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger
    from pyNastran.op2.op2 import OP2


def _get_gpdt_nnodes_numwide(size: int, ndata: int,
                             header_ints,
                             log: SimpleLogger) -> tuple[int, int]:  # pragma: no cover
    """
    size=4; ndata=1120 [102  40   0   0   0   0   0] -> nnodes=40
      1120=(7*4)*4*10 -> 7 words; nnodes=40
      1120=28*40
    """
    unused_table_id, nnodes_nbytes, nwords, *other = header_ints
    #print("header_ints =", header_ints)
    nvalues = ndata // 4
    #nnodes_nbytes,
    if nwords == 0:
        # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\c402pen12f.op2
        # size=8 ndata=1568 nnodes=28
        #    1568 = 2* 4^2 * 7^2
        #
        # [102  18   0   0   0   0   0]
        # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\rbarthm1.op2
        # size=4; ndata=168
        #   168 = 4 * 2 * 3*7
        #if self.op2.table_name in [b'GPDT', b'GPDTS']:
            #is_nodes = True
            ##self.size = 4 168
            #self.log.warning(f'size={self.size} ndata={ndata}')
        #else:  # pragma: no cover
        try:
            is_nodes, numwide, nnodes = _get_gpdt_numwide_from_nodes_null_nwords(size, nnodes_nbytes, ndata)
        except Exception:
            is_nodes = True
            numwide = 0
            raise
        log.debug(f'ndata={ndata} numwide={numwide} nnodes={nnodes}; is_nodes={is_nodes}')

        if is_nodes:
            nnodes = nnodes_nbytes
            #nnodes = ndata // nwords // 4
            #assert nnodes == 3, nnodes
            numwide = nvalues // nnodes
        if numwide not in [7, 10, 14]:
            if size == 8:
                numwide = 14
            else:
                numwide = 7
            nnodes = ndata // (4 * numwide)
            assert ndata % (4 * numwide) == 0
            if numwide not in [7, 10, 14]:
                raise RuntimeError(f'numwide={numwide} must be 7, 10, or 14')
    else:
        nnodes = nnodes_nbytes
        numwide = nvalues // nnodes
        assert ndata % 4 == 0

    #print(f'nnodes={nnodes} numwide={numwide} error={nvalues % nnodes}')
    return nnodes, numwide


def _get_gpdt_numwide_from_nodes_null_nwords(size: bool, nbytes: int,  # pragma: no cover
                                             ndata: int) -> tuple[bool, int]:
    """
    size=4 ndata=392 nnodes_nbytes=28
      392 = 4 * 2. * 7^2
    """
    is_nodes = False
    if nbytes % 2: # or (size == 4 and nnodes_nbytes not in [28, 40]):
        is_nodes = True
        #self.log.warning('is_nodes = True')
        numwide = 0
        nnodes = 0
        return is_nodes, numwide, nnodes

    nnodes = ndata // nbytes
    #print('size =', size, ndata, nbytes)
    assert ndata % nbytes == 0

    nvalues = ndata // 4
    assert ndata % 4 == 0
    numwide = nvalues // nnodes
    assert nvalues % nnodes == 0, f'size={size} ndata={ndata} nvalues={nvalues} nnodes={nnodes}'
    #if (self.size == 4 and nnodes_nbytes in [28, 40]):
    if size == 4:
        assert ndata == numwide * nnodes * 4, f'size=4 ndata={ndata} numwide={numwide} nnodes={nnodes}'
        if numwide not in [7, 10]:
            is_nodes = True
    elif size == 8:
        assert ndata == numwide * nnodes * 4, f'size=8 ndata={ndata} numwide={numwide} nnodes={nnodes}'
        assert numwide == 14, numwide
        #if numwide not in [7, 10]:
            #is_nodes = True
    return is_nodes, numwide, nnodes

def _read_gpdt_8_14(op2: OP2, data: bytes, nnodes: int):
    """
    (nid, cp, x, y, z, cd, 0)
    """
    i = 0
    ntotal = 7 * 8  #  6 ints, 3 floats
    # (nid, cp, x, y, z, cd, 0)
    structi = Struct(op2._endian + b'qq 3d qq')

    for j in range(nnodes):
        edata = data[i:i+ntotal]
        out = structi.unpack(edata)
        #print(out)
        i += ntotal

    #   0   1  2  3  4  5   6
    # nid, cp, x, y, z, cd, ps
    iints = [0, 1, 5, 6]
    ifloats = [2, 3, 4]

    ints = np.frombuffer(data, op2.idtype8).reshape(nnodes, 7).copy()[:, iints]
    floats = np.frombuffer(data, op2.idtype8).reshape(nnodes, 7).copy()[:, ifloats]
    nid_cp_cd_ps = ints
    xyz = floats
    return nid_cp_cd_ps, xyz


def _read_gpdt_4_7(op2: OP2, data: bytes, nnodes: int):
    """
    (nid, cp, x, y, z, cd, 0)
    """
    i = 0
    ntotal = 7 * 4  #  6 ints, 3 floats
    # (nid, cp, x, y, z, cd, 0)
    structi = Struct(op2._endian + b'ii 3f ii')

    #ntotal = 16
    #structi = Struct(self._endian + b'')
    #self.show_data(data, types='if')
    for j in range(nnodes):
        edata = data[i:i+ntotal]
        #self.show_data(edata, types='if')
        out = structi.unpack(edata)
        #print(out)
        nid, zero_a, x, y, z, cd, zero_b = out
        outs = f'nid={nid} zero_a={zero_a} xyz=({x},{y},{z}) cd={cd} zero_b={zero_b}'
        assert nid > 0, outs
        #assert zero_a == 0, (nid, zero_a, x, y, z, cd, zero_b)
        #assert zero_b == 0, (nid, zero_a, x, y, z, cd, zero_b)
        i += ntotal

    #   0   1  2  3  4  5   6
    # nid, cp, x, y, z, cd, ps
    iints = [0, 1, 5, 6]
    ifloats = [2, 3, 4]

    ints = np.frombuffer(data, op2.idtype).reshape(nnodes, 7).copy()[:, iints]
    floats = np.frombuffer(data, op2.idtype).reshape(nnodes, 7).copy()[:, ifloats]
    nid_cp_cd_ps = ints
    xyz = floats
    return nid_cp_cd_ps, xyz

def _read_gpdt_4_10(op2: OP2, data: bytes, nnodes: int):
    """
    (nid, 0, x, y, z, cd, 0)
    """
    i = 0
    ntotal = 40  #  6 ints, 3 floats
    # (nid, 0, x, y, z, cd, 0)
    structi = Struct(op2._endian + b'2i 3d 2i')
    for j in range(nnodes):
        edata = data[i:i+ntotal]
        #self.show_data(edata, types='ifqd')
        out = structi.unpack(edata)
        nid, zero_a, x, y, z, cd, zero_b = out
        outs = f'nid={nid} zero_a={zero_a} xyz=({x},{y},{z}) cd={cd} zero_b={zero_b}'
        #assert zero_a == 0, (nid, zero_a, x, y, z, cd, zero_b)
        #assert zero_b == 0, (nid, zero_a, x, y, z, cd, zero_b)
        i += ntotal
        #print(outs)
        assert nid > 0, outs

    #   0   1    2  3  4  5  6  7  8   9
    # nid, zero, x, _, y, _, z, _, cd, zero
    iints = [0, 1, 8, 9]

    # 0  1  2  3  4
    # a, x, y, z, b
    ifloats = [1, 2, 3]

    ints = np.frombuffer(data, op2.idtype).reshape(nnodes, 10).copy()[:, iints]
    floats = np.frombuffer(data, 'float64').reshape(nnodes, 5).copy()[:, ifloats]
    nid_cp_cd_ps = ints
    xyz = floats
    return nid_cp_cd_ps, xyz


def _get_gpdt_nnodes2(ndata: int, header_ints: list[int], size: int):
    unused_table_id, nnodes_nbytes, nwords, *other = header_ints
    nvalues = ndata // 4
    assert nvalues > 0
    assert ndata % 4 == 0
    try:
        # assume nodes
        nnodes = nnodes_nbytes
        numwide = nvalues // nnodes
        assert nvalues % nnodes == 0
        if size == 4:
            assert numwide in [7, 10], numwide
        else:
            assert numwide == 14, numwide
    except AssertionError:
        # calculate the bytes
        if size == 4:
            numwide = 7
        elif numwide == 8:
            numwide = 14
        nnodes = nvalues // numwide
    assert ndata == nnodes * numwide * 4
    return nnodes, numwide
