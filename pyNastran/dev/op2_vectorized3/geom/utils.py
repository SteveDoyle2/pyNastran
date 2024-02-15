from typing import Tuple
import numpy as np

def get_minus1_start_end(ints: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    iminus1 = np.where(ints == -1)[0]
    istart = [0] + list(iminus1[:-1] + 1)
    iend = iminus1
    return istart, iend

def get_ints(data: bytes, n: int,
             nrows: int, ncols: int,
             size: int=4, endian: bytes=b'<') -> tuple[int, np.ndarray]:
    if size == 4:
        ints = np.frombuffer(data, dtype='int32', offset=n).reshape(nrows, ncols)
    else:
        #intsi = np.frombuffer(data, dtype='int64', offset=n).tolist()
        ints = np.frombuffer(data, dtype='int64', offset=n).reshape(nrows, ncols)
    if endian == b'>':
        ints = ints.newbyteorder('>')
    n = len(data)
    return n, ints

def get_ints_floats(data: bytes, n: int,
                    nrows: int, ncols: int,
                    size: int=4, endian: bytes=b'<') -> tuple[int, np.ndarray, np.ndarray]:
    if size == 4:
        ints = np.frombuffer(data, dtype='int32', offset=n).reshape(nrows, ncols)
        floats = np.frombuffer(data, dtype='float32', offset=n).reshape(nrows, ncols)
    else:
        #intsi = np.frombuffer(data, dtype='int64', offset=n)
        #print(intsi.tolist(), len(intsi), nrows, ncols)
        ints = np.frombuffer(data, dtype='int64', offset=n).reshape(nrows, ncols)
        floats = np.frombuffer(data, dtype='float64', offset=n).reshape(nrows, ncols)
    if endian == b'>':
        ints = ints.newbyteorder('>')
        floats = floats.newbyteorder('>')
    n = len(data)
    return n, ints, floats

def get_ints_strings(data: bytes, n: int,
                     nrows: int, ncols: int,
                     size: int=4, endian: str='<') -> tuple[int, np.ndarray, np.ndarray]:
    if size == 4:
        ints = np.frombuffer(data, dtype='int32', offset=n).reshape(nrows, ncols)
        strings = np.frombuffer(data, dtype='<S4', offset=n).reshape(nrows, ncols)
    else:
        ints = np.frombuffer(data, dtype='int64', offset=n).reshape(nrows, ncols)
        strings = np.frombuffer(data, dtype='<S4', offset=n).reshape(nrows, ncols*2)[:, ::2]
    if endian == b'>':
        ints = ints.newbyteorder('>')
        strings = strings.newbyteorder('>')
    n = len(data)
    return n, ints, strings

def get_ints_floats_strings(data: bytes, n: int,
                            nrows: int, ncols: int,
                            size: int=4, endian: str='<') -> tuple[int, np.ndarray,
                                                                   np.ndarray, np.ndarray]:
    if size == 4:
        ints = np.frombuffer(data, dtype='int32', offset=n).reshape(nrows, ncols)
        floats = np.frombuffer(data, dtype='float32', offset=n).reshape(nrows, ncols)
        strings = np.frombuffer(data, dtype='<S4', offset=n).reshape(nrows, ncols)
    else:
        ints = np.frombuffer(data, dtype='int64', offset=n).reshape(nrows, ncols)
        floats = np.frombuffer(data, dtype='float64', offset=n).reshape(nrows, ncols)
        strings = np.frombuffer(data, dtype='<S4', offset=n).reshape(nrows, ncols*2)[:, ::2]
    if endian == b'>':
        ints = ints.newbyteorder('>')
        floats = floats.newbyteorder('>')
        strings = strings.newbyteorder('>')
    n = len(data)
    return n, ints, floats, strings
