from typing import List, Any
import numpy as np

def fix_table3_types(table3, size: int=4) -> List[Any]:
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
            #print('write_table_3', v)
            n += len(v)
        table3_new.append(v)
    assert n == 584, n
    return table3_new
