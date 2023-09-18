import sys
import warnings
from typing import Union

import numpy as np
from numpy import float32, float64
from pyNastran.dev.bdf_vectorized3.bdf_interface.fast_float_print import get_float_format # print_float_8,
from pyNastran.bdf.field_writer_8 import print_float_8
from pyNastran.bdf.bdf_interface.assign_type import double_from_str

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

def print_card_8(fields: list[Union[int, float, str, None]]) -> str:
    """
    Prints a nastran-style card with 8-character width fields.

    Parameters
    ----------
    fields : list[int/float/str/None]
        all the fields in the BDF card (no trailing Nones)

    Returns
    -------
    card : str
        string representation of the card in small field format

    .. note:: An internal field value of None or '' will be treated as
              a blank field
    .. note:: A small field format follows the  8-8-8-8-8-8-8-8 = 80
              format where the first 8 is the card name or
              blank (continuation).  The last 8-character field indicates
              an optional continuation, but because it's a left-justified
              unnecessary field, print_card doesn't use it.

    .. code-block:: python

       >>> fields = ['DUMMY', 1, 2, 3, None, 4, 5, 6, 7, 8.]
       >>> print_card_8(fields)
       DUMMY          1       2       3               4       5       6       7
                     8.

    """
    try:
        out = '%-8s' % fields[0]
    except Exception:
        warnings.warn("ERROR!  fields=%s" % fields)
        sys.stdout.flush()
        raise

    for i in range(1, len(fields)):
        field = fields[i]
        try:
            out += print_field_8(field)
        except Exception:
            warnings.warn("bad fields = %s" % fields)
            raise
        if i % 8 == 0:  # allow 1+8 fields per line
            out = out.rstrip(' ')
            if out[-1] == '\n':  # empty line
                out += '+'
            out += '\n        '
    out = out.rstrip(' \n+') + '\n'  # removes blank lines at the end of cards
    return out

def print_field_8(value: Union[int, float, str, None]) -> str:
    """
    Prints a 8-character width field

    Parameters
    ----------
    value : int/float/str
        the value to print

    Returns
    -------
    field : str
        an 8-character string

    """
    assert not isinstance(value, list), value

    if isinstance(value, int):
        field = '%8i' % value
    elif isinstance(value, (float, float32, float64)):
        field_old = '%8s' % print_float_8(value)
        field_new = get_float_format(value)
        assert len(field_new) == 8, f'{field_new!r}; n={len(field_new)}'

        try:
            value_old = double_from_str(field_old.strip())
        except:
            raise SyntaxError(f'value={value}; field_old={field_old!r}')

        try:
            value_new = double_from_str(field_new.strip())
        except:
            raise SyntaxError(f'value={value}; field_new={field_new!r}')
        if not value_old == value_new:
            print(f'value={value}; field_old={field_old!r} field_new={field_new!r}')

        #assert value_old == value_new, f'value={value}; field_old={field_old!r} field_new={field_new!r}'
        field = field_new

    elif value is None:
        field = '        '
    else:
        field = '%8s' % value

    if len(field) != 8:
        nchars = len(field)
        msg = f'field={field!r} is not 8 characters long; its {nchars}...raw_value={value!r}'
        raise RuntimeError(msg)
    return field
