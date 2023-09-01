import numpy as np

def hstack_msg(mylist, msg: str, min_size=0) -> np.ndarray:
    if isinstance(mylist, list) and len(mylist) == 0:
        raise ValueError(f'empty list; {msg}')

    try:
        stacked = np.hstack(mylist)
    except ValueError:
        if len(mylist) == 0:
            raise ValueError(f'empty list; {msg}')
        raise
    if len(stacked) == 0:
        raise ValueError(f'empty list; {msg}')
    return stacked

def cast_int_array(list_ints: list[int]) -> np.ndarray:
    try:
        return np.array(list_ints, dtype='int32')
    except OverflowError:
        return np.array(list_ints, dtype='int64')
