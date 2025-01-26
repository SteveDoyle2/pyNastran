import numpy as np
from pyNastran.bdf.field_writer_8 import print_card_8, print_float_8
from pyNastran.bdf.field_writer_16 import print_card_16, print_float_16
from pyNastran.bdf.field_writer_double import print_scientific_double
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

#def get_print_card(size: int, max_int: int) -> Callable:
    #if size == 16 or max_int > MAX_8_CHAR_INT:
        #print_card = print_card_16
    #else:
        #print_card = print_card_8
    #return print_card

def update_field_size(max_int: int, size: int) -> int:
    if max_int > MAX_8_CHAR_INT:
        size = 16
    return size

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
    assert ndarray.dtype.name not in {'float32', 'float64'}, ndarray.dtype.name
    if size == 8:
        str_array = ndarray.astype('|U8')
    else:
        str_array = ndarray.astype('|U16')
    return str_array

def array_default_int(ndarray: np.ndarray, default: int=0, size: int=8) -> np.ndarray:
    assert ndarray.dtype.name not in {'float32', 'float64'}, f'name={ndarray.dtype.name!r}'
    idefault = np.where(ndarray == default)
    if size == 8:
        str_array = ndarray.astype('|U8')
    else:
        str_array = ndarray.astype('|U16')
    str_array[idefault] = ''
    return str_array

def array_float(ndarray: np.ndarray, size: int=8, is_double: bool=False,
                nan_check: bool=True) -> np.ndarray:
    """setup the nan values and fill in the holes"""
    assert ndarray.dtype.name in {'float32', 'float64'}, ndarray.dtype.name
    if nan_check and np.any(np.isnan(ndarray)):
        raise RuntimeError('nans found')

    if size == 8:
        str_array = np.zeros(ndarray.shape, dtype='|U8')
        print_float = print_float_8
    else:
        str_array = np.zeros(ndarray.shape, dtype='|U16')
        print_float = print_float_16
        if is_double:
            print_float = print_scientific_double

    if ndarray.ndim == 1:
        for i, value in enumerate(ndarray):
            str_array[i] = print_float(value)
    elif ndarray.ndim == 2:
        for i, values in enumerate(ndarray):
            for j, value in enumerate(values):
                str_array[i, j] = print_float(value)
    else:
        raise NotImplementedError(ndarray.shape)
    return str_array

def array_float_nan(ndarray: np.ndarray, size: int=8, is_double: bool=False) -> np.ndarray:
    """setup the nan values and fill in the holes"""
    assert ndarray.dtype.name in {'float32', 'float64'}, ndarray.dtype.name
    inan = np.isnan(ndarray)
    ivalue = ~inan
    if size == 8:
        str_array = np.zeros(ndarray.shape, dtype='|U8')
        print_float = print_float_8
    else:
        str_array = np.zeros(ndarray.shape, dtype='|U16')
        print_float = print_float_16
        if is_double:
            print_float = print_scientific_double

    if ndarray.ndim == 1:
        i = np.where(ivalue)[0]
        for ii in i:
            str_array[ii] = print_float(ndarray[ii])
    elif ndarray.ndim == 2:
        i, j = np.where(ivalue)
        for ij in zip(i, j):
            str_array[ij] = print_float(ndarray[ij])
    else:
        raise NotImplementedError(ndarray.shape)

    str_array[inan] = ''
    return str_array

def array_default_float(ndarray: np.ndarray, default: float=0.,
                        size: int=8, is_double: bool=False) -> np.ndarray:
    """setup the nan values and fill in the holes"""
    assert ndarray.dtype.name in {'float32', 'float64'}, ndarray.dtype.name
    if np.any(np.isnan(ndarray)):
        raise RuntimeError('nans found')

    idefault = (ndarray == default)
    ivalue = ~idefault
    if size == 8:
        str_array = np.zeros(ndarray.shape, dtype='|U8')
        print_float = print_float_8
    else:
        str_array = np.zeros(ndarray.shape, dtype='|U16')
        print_float = print_float_16
        if is_double:
            print_float = print_scientific_double

    if ndarray.ndim == 1:
        i = np.where(ivalue)[0]
        for ii in i:
            str_array[ii] = print_float(ndarray[ii])
    elif ndarray.ndim == 2:
        i, j = np.where(ivalue)
        for ij in zip(i, j):
            str_array[ij] = print_float(ndarray[ij])
    else:
        raise NotImplementedError(ndarray.shape)

    str_array[idefault] = ''
    return str_array

def array_default_float_nan(ndarray: np.ndarray, default=0.,
                            size: int=8, is_double: bool=False) -> np.ndarray:
    """setup the nan values and fill in the holes"""
    assert ndarray.dtype.name in {'float32', 'float64'}, ndarray.dtype.name
    idefault = (np.isnan(ndarray) | (ndarray == default))
    ivalue = ~idefault
    if size == 8:
        str_array = np.zeros(ndarray.shape, dtype='|U8')
        print_float = print_float_8
    else:
        str_array = np.zeros(ndarray.shape, dtype='|U16')
        print_float = print_float_16
        if is_double:
            print_float = print_scientific_double

    if ndarray.ndim == 1:
        i = np.where(ivalue)[0]
        for ii in i:
            str_array[ii] = print_float(ndarray[ii])
    elif ndarray.ndim == 2:
        i, j = np.where(ivalue)
        for ij in zip(i, j):
            str_array[ij] = print_float(ndarray[ij])
    else:
        raise NotImplementedError(ndarray.shape)

    str_array[idefault] = ''
    return str_array

def array_default_floats(ndarray: np.ndarray, defaults: np.ndarray,
                         size: int=8, is_double: float=False,
                         nan_check: bool=True) -> np.ndarray:
    assert ndarray.dtype.name in {'float32', 'float64'}, ndarray.dtype.name
    values_str = array_float(ndarray, size=size, is_double=is_double, nan_check=nan_check)
    defaults_str = array_float(defaults, size=size, is_double=is_double, nan_check=True)
    idefault = (values_str == defaults_str)
    values_str[idefault] = ''
    return values_str

def array_default_str(ndarray: np.ndarray, default: str='', size: int=8) -> np.ndarray:
    assert np.dtype.name not in {'float32', 'float64'}, np.dtype.name
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
