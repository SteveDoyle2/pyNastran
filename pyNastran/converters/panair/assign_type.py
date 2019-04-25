"""
defines methods for reading panair values:
 - double(value, name)
 - integer(value, name)
 - integer_or_blank(value, name, default=None)
 - double_or_blank(value, name, default=None)

defines methods for writing panair values:

 - fortran_value(value)
"""
from __future__ import  print_function
from typing import List, Union, Optional

def double(value, name):
    # type: (str, str) -> float
    """casts to an float value"""
    if isinstance(value, float):
        return value
    fvalue = float(value)
    return fvalue

def integer(value, name):
    # type: (str, str) -> int
    """casts to an integer value"""
    if isinstance(value, int):
        return value
    value = value
    fvalue = float(value)
    if not fvalue.is_integer():
        raise RuntimeError('%s=%r is not an integer' % (name, fvalue))
    return int(fvalue)

def fortran_value(value):
    # type: (float) -> str
    return "%8.4E" % value

def integer_or_blank(value, name, default=None):
    # type: (str, str, Optional[Union[float, int]]) -> Optional[Union[float, int]]
    value = value.strip()
    if not value:
        return default

    fvalue = float(value)
    if not fvalue.is_integer():
        raise RuntimeError('%s=%r is not an integer' % (name, fvalue))
    return int(fvalue)

def double_or_blank(value, name, default=None):
    # type: (str, str, Optional[float]) -> Optional[float]
    value = value.strip()
    if not value:
        return default
    try:
        fvalue = float(value)
    except ValueError:
        raise SyntaxError('%s=%r is not a float' % (name, value))
    return fvalue
