import numpy as np
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from typing import Callable

MAX_8_CHAR_INT = 99_999_999

def print_card_8_comment(fields: list[str]) -> str:
    msg = '$%-7s' % fields[0]
    for i, word in enumerate(fields[1:]):
        msg += '%8s' % word
        if i > 0 and i % 8 == 0:
            msg += '\n$       '
    return msg.rstrip('$ \n') + '\n'

def print_card_16_comment(fields: list[str]) -> str:
    msg = '$%-7s' % fields[0]
    for i, word in enumerate(fields[1:]):
        msg += '%16s' % word
        if i > 0 and i % 4 == 0:
            msg += '\n$       '
    return msg.rstrip('$ \n') + '\n'

def get_print_card(size: int, max_int: int) -> Callable:
    if size == 16 or max_int > MAX_8_CHAR_INT:
        print_card = print_card_16
    else:
        print_card = print_card_8
    return print_card

def get_print_card_size(size: int, max_int: int) -> tuple[Callable, int]:
    if size == 16 or max_int > MAX_8_CHAR_INT:
        print_card = print_card_16
        size = 16
    else:
        print_card = print_card_8
    return print_card, size

#def array_8(ndarray: np.ndarray) -> np.ndarray:
    #str_array = ndarray.astype('|U8')
    #return str_array

def array_str(ndarray: np.ndarray, size: int=8) -> np.ndarray:
    if size == 8:
        str_array = ndarray.astype('|U8')
    else:
        str_array = ndarray.astype('|U16')
    return str_array

def array_default_int(ndarray: np.ndarray, default: int=0, size: int=8) -> np.ndarray:
    idefault = np.where(ndarray == default)
    if size == 8:
        str_array = ndarray.astype('|U8')
    else:
        str_array = ndarray.astype('|U16')
    str_array[idefault] = ''
    return str_array

def array_default_str(ndarray: np.ndarray, default: str='', size: int=8) -> np.ndarray:
    idefault = np.where(ndarray == default)
    val0 = ndarray.ravel()[0]
    assert isinstance(val0, str), f'array is not made of unicode strings (e.g., "|U8"); dtype="{ndarray.dtype}"'
    if size == 8:
        str_array = ndarray.astype('|U8')
    else:
        str_array = ndarray.astype('|U16')
    str_array[idefault] = ''
    return str_array

#def array_default_8(ndarray: np.ndarray, default: str='0') -> np.ndarray:
    #str_array = ndarray.astype('|U8')
    #idefault = np.where(str_array == default)[0]
    #str_array[idefault] = ''
    #return str_array

#def array_default_16(ndarray: np.ndarray, default: str='0') -> np.ndarray:
    #str_array = ndarray.astype('|16')
    #idefault = np.where(str_array == default)[0]
    #str_array[idefault] = ''
    #return str_array
