"""
defines:
 - expand_thru
 - expand_thru_by
 - expand_thru_exclude
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from typing import  List, Union
from six import string_types

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.bdf_interface.assign_type import interpret_value


def expand_thru(fields, set_fields=True, sort_fields=False):
    # type: (List[str], bool, bool) -> List[int]
    """
    Expands a list of values of the form [1,5,THRU,9,13]
    to be [1,5,6,7,8,9,13]

    Parameters
    ----------
    fields : List[int/str]
        the fields to expand
    set_fields : bool; default=True
        Should the fields be converted to a set and then back to a list?
        This is useful for [2, 'THRU' 5, 1]
    sort_fields : bool; default=False
        Should the fields be sorted at the end?

    """
    # ..todo:  should this be removed...is the field capitalized when read in?
    if isinstance(fields, integer_types):
        return [fields]
    #elif isinstance(fields[0], integer_types):  # don't use this [1, 'THRU', 10]
        #return fields
    elif len(fields) == 1:
        return [int(fields[0])]

    fields = [field.upper()
              if isinstance(field, string_types) else field for field in fields]

    out = []
    nfields = len(fields)
    i = 0
    while i < nfields:
        if isinstance(fields[i], string_types) and fields[i] == 'THRU':
            istart = int(fields[i - 1])
            iend = int(fields[i + 1])

            # adding 1 to iend for the range offset
            for j in range(istart+1, iend + 1):
                out.append(j)
            i += 2
        else:
            out.append(int(fields[i]))
            i += 1

    if set_fields:
        out = list(set(out))
    if sort_fields:
        out.sort()
    return out


def expand_thru_by(fields, set_fields=True, sort_fields=True,
                   require_int=True, allow_blanks=False):
    # type: (List[str], bool, bool, bool, bool) -> List[int]
    """
    Expands a list of values of the form [1,5,THRU,9,BY,2,13]
    to be [1,5,7,9,13]

    Parameters
    ----------
    fields : List[int/str]
        the fields to expand
    set_fields : bool; default=True
        Should the fields be converted to a set and then back to a list
        to remove duplicates?
        This is useful for [2, 'THRU' 5, 1]
    sort_fields : bool; default=False
        Should the fields be sorted at the end?
    require_int : bool; default=True
        True : all data must be integers
        False : floats are allowed (e.g., DDVAL)
    allow_blanks : bool; default=Fals
        True : blank/Nones are ignored (e.g., NSM1/NSML1)
        False : crash

    .. todo:: not tested

    Notes
    -----
    used for QBDY3 and what else ???

    """
    if require_int:
        func = int
    else:
        func = interpret_value

    # ..todo:  should this be removed...is the field capitalized when read in?
    fields = [field.upper()
              if isinstance(field, string_types) else field for field in fields]

    if len(fields) == 1:
        return [func(fields[0])]
    out = []
    nfields = len(fields)
    i = 0
    by = 1
    while i < nfields:
        #print('fields[i]=%r' % fields[i])
        is_blank = (
            allow_blanks and (
                (isinstance(fields[i], string_types) and fields[i].strip() == '') or
                fields[i] is None)
        )
        if is_blank:
            #print('blank=%s' % fields[i])
            i += 1
            continue
        if fields[i] == 'THRU':
            by = 1
            by_case = False
            if i + 2 < nfields and fields[i + 2] == 'BY':
                by = func(fields[i + 3])
            else:
                by = 1
                by_case = True
            min_value = func(fields[i - 1])
            max_value = func(fields[i + 1])
            max_range = int((max_value - min_value) // by + 1)  # max range value

            for j in range(0, max_range):  # +1 is to include final point
                value = min_value + by * j
                out.append(value)
            out.append(max_value)

            if by_case:  # null/standard case
                # A thru B
                i += 2
            else:     # BY case
                # A thru B by C
                i += 4
        else:
            out.append(func(fields[i]))
            i += 1

    if set_fields:
        out = list(set(out))
    if sort_fields:
        out.sort()
    return out


def expand_thru_exclude(fields):
    # type: (List[str]) -> List[int]
    """
    Expands a list of values of the form [1,5,THRU,11,EXCEPT,7,8,13]
    to be [1,5,6,9,10,11,13]

    .. warning:: hasn't been tested

    """
    # ..todo:  should this be removed...is the field capitalized when read in?
    isfields = [interpret_value(field.upper())
                if isinstance(field, string_types) else field
                for field in fields]  # type: List[Union[str,int]]

    fields_out = []  # type: List[int]
    nfields = len(isfields)
    for i in range(nfields):
        #print('fields[%i] = %r' % (i, isfields[i]))
        if isfields[i] == 'THRU':
            sorted_list = []
            for j in range(isfields[i - 1], isfields[i + 1]):
                sorted_list.append(isfields[j])

        elif isfields[i] == 'EXCLUDE':
            stored_set = set(sorted_list)
            while isfields[i] < max(sorted_list):
                stored_set.remove(isfields[i])
            sorted_list = list(stored_set)
        else:
            if sorted_list:
                fields_out += sorted_list
            fields_out.append(isfields[i])
    return fields_out
