# pylint: disable=C0301,W0201
from __future__ import annotations
from struct import Struct
from typing import Union, Any, TYPE_CHECKING

import numpy as np
from pyNastran.op2.op2_helper import polar_to_real_imag
from pyNastran.op2.op2_interface.op2_reader import mapfmt

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2

def read_real_table_static(op2: OP2, obj: Any, flag:str,
                           data: bytes, nnodes: int, ntotal: int) -> int:
    n = 0
    dt = np.nan
    fmt = mapfmt(op2._endian + b'2i6f', op2.size)
    s = Struct(fmt)
    for unused_inode in range(nnodes):
        out = s.unpack(data[n:n+ntotal])
        (eid_device, grid_type, tx, ty, tz, rx, ry, rz) = out
        #print(out)
        eid = eid_device // 10
        if op2.is_debug_file:
            op2.binary_debug.write('  %s=%i; %s\n' % (flag, eid, str(out)))
        obj.add_sort1(dt, eid, grid_type, tx, ty, tz, rx, ry, rz)
        n += ntotal
    return n

def read_real_table_sort1(op2: OP2, obj: Any, dt: Union[int, float], flag: str,
                          data: bytes, nnodes: int, ntotal: int) -> int:
    n = 0
    assert nnodes > 0, nnodes
    fmt = mapfmt(op2._endian + b'2i6f', op2.size)
    s = Struct(fmt)
    for unused_inode in range(nnodes):
        out = s.unpack(data[n:n+ntotal])
        (eid_device, grid_type, tx, ty, tz, rx, ry, rz) = out
        #print(out)
        assert grid_type != 1065353216, out  # caused by an op2 writer bug with int64 numbers being downcast directly to float32
        eid = eid_device // 10
        if op2.is_debug_file:
            op2.binary_debug.write('  %s=%i; %s\n' % (flag, eid, str(out)))
        obj.add_sort1(dt, eid, grid_type, tx, ty, tz, rx, ry, rz)
        n += ntotal
    return n

def read_complex_table_sort1_mag(op2: OP2, obj: Any, dt: Union[int, float], flag: str,
                                 data: bytes, nnodes: int, ntotal: int) -> int:
    n = 0
    fmt = mapfmt(op2._endian + b'2i12f', op2.size)
    s = Struct(fmt)
    for unused_inode in range(nnodes):
        out = s.unpack(data[n:n+ntotal])
        (eid_device, grid_type, txr, tyr, tzr, rxr, ryr, rzr,
         txi, tyi, tzi, rxi, ryi, rzi) = out
        eid = eid_device // 10
        if op2.is_debug_file:
            op2.binary_debug.write('  %s=%i %s\n' % (flag, eid, str(out)))
        tx = polar_to_real_imag(txr, txi)
        ty = polar_to_real_imag(tyr, tyi)
        tz = polar_to_real_imag(tzr, tzi)
        rx = polar_to_real_imag(rxr, rxi)
        ry = polar_to_real_imag(ryr, ryi)
        rz = polar_to_real_imag(rzr, rzi)
        obj.add_sort1(dt, eid, grid_type, tx, ty, tz, rx, ry, rz)
        n += ntotal
    return n

def read_complex_table_sort1_imag(op2: OP2, obj: Any, dt: Union[int, float], flag: str,
                                  data: bytes, nnodes: int, ntotal: int) -> int:
    n = 0
    fmt = mapfmt(op2._endian + b'2i12f', op2.size)
    s = Struct(fmt)

    assert op2.obj is not None
    assert nnodes > 0
    for unused_inode in range(nnodes):
        out = s.unpack(data[n:n+ntotal])
        (eid_device, grid_type, txr, tyr, tzr, rxr, ryr, rzr,
         txi, tyi, tzi, rxi, ryi, rzi) = out
        eid = eid_device // 10
        if op2.is_debug_file:
            op2.binary_debug.write('  %s=%i %s\n' % (flag, eid, str(out)))
        tx = complex(txr, txi)
        ty = complex(tyr, tyi)
        tz = complex(tzr, tzi)
        rx = complex(rxr, rxi)
        ry = complex(ryr, ryi)
        rz = complex(rzr, rzi)
        obj.add_sort1(dt, eid, grid_type, tx, ty, tz, rx, ry, rz)
        n += ntotal
    return n

def read_complex_table_sort2_imag(op2: OP2, obj: Any, node_id: int,
                                  flag: str, flag_type: str,
                                  data: bytes, nnodes: int, ntotal: int) -> int:
    n = 0
    #ntotal = 56  # 14 * 4
    fmt = mapfmt(op2._endian + op2._analysis_code_fmt + b'i12f', op2.size)
    s = Struct(fmt)
    assert op2.obj is not None
    assert nnodes > 0
    #assert ndata % ntotal == 0

    binary_debug_fmt = '  %s=%s %%s\n' % (flag, flag_type)

    for unused_inode in range(nnodes):
        edata = data[n:n+ntotal]
        out = s.unpack(edata)

        (freq, grid_type, txr, tyr, tzr, rxr, ryr, rzr,
         txi, tyi, tzi, rxi, ryi, rzi) = out

        if op2.is_debug_file:
            op2.binary_debug.write(binary_debug_fmt % (freq, str(out)))
        tx = complex(txr, txi)
        ty = complex(tyr, tyi)
        tz = complex(tzr, tzi)
        rx = complex(rxr, rxi)
        ry = complex(ryr, ryi)
        rz = complex(rzr, rzi)
        obj.add_sort2(freq, node_id, grid_type, tx, ty, tz, rx, ry, rz)
        n += ntotal
    return n

def read_complex_table_sort2_mag(op2: OP2, obj: Any, node_id: int,
                                 flag: str, flag_type: str,
                                 data: bytes, nnodes: int, ntotal: int) -> int:
    n = 0
    s = Struct(mapfmt(op2._endian + op2._analysis_code_fmt + b'i12f', op2.size))
    binary_debug_fmt = '  %s=%s %%s\n' % (flag, flag_type)
    for unused_inode in range(nnodes):
        edata = data[n:n+ntotal]
        out = s.unpack(edata)
        (freq, grid_type, txr, tyr, tzr, rxr, ryr, rzr,
         txi, tyi, tzi, rxi, ryi, rzi) = out

        if op2.is_debug_file:
            op2.binary_debug.write(binary_debug_fmt % (freq, str(out)))
        tx = polar_to_real_imag(txr, txi)
        ty = polar_to_real_imag(tyr, tyi)
        tz = polar_to_real_imag(tzr, tzi)
        rx = polar_to_real_imag(rxr, rxi)
        ry = polar_to_real_imag(ryr, ryi)
        rz = polar_to_real_imag(rzr, rzi)
        obj.add_sort2(freq, node_id, grid_type, tx, ty, tz, rx, ry, rz)
        n += ntotal
    return n

def read_real_table_sort2(self: OP2, obj: Any, flag: str, nid: int,
                          data: bytes, nnodes: int, ntotal: int) -> int:
    n = 0
    assert nnodes > 0
    fmt = mapfmt(self._endian + self._analysis_code_fmt + b'i6f', self.size)
    structi = Struct(fmt)

    #psds = ('CRM2', 'NO2', 'PSD2', 'RMS2')
    #print('sort_method=%s' % self.sort_method)
    #if self.table_name_str.endswith(psds):
    for unused_inode in range(nnodes):
        edata = data[n:n+ntotal]
        out = structi.unpack(edata)
        (dt, grid_type, tx, ty, tz, rx, ry, rz) = out
        #print(out)
        if self.is_debug_file:
            self.binary_debug.write(
                    f'  nid={nid} {flag}={dt} ({type(dt)}); {str(out)}\n')
        obj.add_sort2(dt, nid, grid_type, tx, ty, tz, rx, ry, rz)
        n += ntotal
    return n
