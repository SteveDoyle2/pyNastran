import os
from typing import BinaryIO
from pyNastran.utils import PathLike, print_bad_path
import struct

def read_zona_aic(aic_filename: PathLike):
    assert os.path.exists(aic_filename), print_bad_path(aic_filename)
    with open(aic_filename, "rb") as aic_file:
        data = aic_file.read()
    #show_n(data, 1000)
    #asdf

def show_n(data: bytes, n: int):
    n = min(n, len(data))
    print(len(data))
    sfmt = f'{n}s'
    datan = data[:n]
    sdata = struct.unpack(sfmt, datan)
    nints = n // 4
    ifmt = f'{nints}i'
    ffmt = f'{nints}f'
    print(ifmt)
    idata = struct.unpack(ifmt, datan)
    fdata = struct.unpack(ffmt, datan)
    print(f'strings = {sdata}')
    print(f'ints    = {idata}')
    print(f'floats  = {fdata}')
