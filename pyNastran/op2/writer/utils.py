from typing import Any
import numpy as np

def op2_stringify(values: list) -> str:  # pragma: no cover
    out = []
    for value in values:
        if isinstance(value, (float, np.float32, np.float64)):
            out.append('%.3f' % float(value))
        elif isinstance(value, (int, np.int32, np.int64)):
            out.append(int(value))
        elif isinstance(value, (bytes, str)):
            out.append(value)
        else:
            raise NotImplementedError(value)
    return str(out)

def fix_table3_types(table3, size: int=4) -> list[Any]:
    assert size == 4, size
    table3_new = []
    n = 0
    for v in table3:
        if isinstance(v, (int, float, np.int32, np.float32)):
            n += 4
        elif isinstance(v, str):
            n += len(v)
        elif isinstance(v, np.int64):
            v = v.astype('int32')
            n += 4
        elif isinstance(v, np.float64):
            v = v.astype('float32')
            n += 4
        else:
            assert isinstance(v, bytes), v
            #print('write_table_3', v)
            n += len(v)
        table3_new.append(v)
    assert n == 584, n
    return table3_new
