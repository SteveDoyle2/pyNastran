"""验证 EXTRN _read_extrn 对非标准分量值的处理。

NASTRAN ASET 的 C 字段合法值是 DOF 1-6 的任意排序组合（共 63 个 + 0），
但原始代码只允许 [1, 2, 3, 4, 5, 6, 123, 123456] 共 8 个值，
遇到 456、1234 等合法值就 assert 崩溃，导致整个 OP2 解析中断。
"""
import sys
import struct
import numpy as np
from pathlib import Path

REPO = Path(__file__).parent
sys.path.insert(0, str(REPO))

from pyNastran.op2.tables.geom.geom1 import GEOM1


class MockLog:
    def warning(self, msg):
        print(f"  [WARNING] {msg}")

    def info(self, msg):
        pass

    def debug(self, msg):
        pass


class _MockGeom2:
    def read_cmass2(self, *a, **k):
        pass


class MockOP2:
    idtype8 = np.dtype('<i4')
    fdtype8 = np.dtype('<f4')
    size = 4
    factor = 1
    _endian = b'<'
    log = MockLog()
    asets = {}
    reader_geom2 = _MockGeom2()

    def add_aset(self, nids, comps):
        print(f"  add_aset called: nids={nids}, comps={comps}")
        self.asets[id(tuple(nids))] = (nids, comps)


def test_extrn_component_values():
    """测试 _read_extrn 对各种分量值的处理。"""
    geom1 = GEOM1(MockOP2())

    # 构造 EXTRN 二进制数据: (GID, C) pairs, terminated by (-1, -1)
    # 包含原始白名单外的合法值: 456, 1234, 0
    data = struct.pack('<10i',
                       1, 456,      # 非标准但合法 (DOF 4,5,6)
                       2, 1234,     # 非标准但合法 (DOF 1,2,3,4)
                       3, 123456,   # 标准
                       4, 0,        # 0 = 无约束
                       -1, -1)      # 终止符
    n = 0

    try:
        result = geom1._read_extrn(data, n)
        print(f"PASS: _read_extrn completed without crash, returned {result}")
        print(f"  asets registered: {len(geom1.op2.asets)}")
        return True
    except AssertionError as e:
        print(f"FAIL: AssertionError on non-standard component: {e}")
        return False


if __name__ == '__main__':
    print("=== Testing EXTRN with non-standard component values ===")
    ok = test_extrn_component_values()
    sys.exit(0 if ok else 1)
