"""
defines:
 - expand_thru
 - expand_thru_by
 - expand_thru_exclude

"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from typing import  List, Union

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
              if isinstance(field, str) else field for field in fields]

    out = []
    nfields = len(fields)
    i = 0
    while i < nfields:
        if isinstance(fields[i], str) and fields[i] == 'THRU':
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
              if isinstance(field, str) else field for field in fields]

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
                (isinstance(fields[i], str) and fields[i].strip() == '') or
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


def split_comma_space(datai):
    if ',' in datai:
        sline1 = [slinei.strip() for slinei in datai.split(',')]
        sline = []
        for slinei in sline1:
            if ' ' in slinei:
                sline.extend(slinei.split())
            else:
                sline.append(slinei)
    else:
        sline = datai.split()
    return sline

def setup_data(data_in):
    data = []
    for datai in data_in:
        if isinstance(datai, int):
            data.append(datai)
        else:
            datai = datai.upper()
            sline = split_comma_space(datai)
            data.extend(sline)

    stype = 'int'
    for datai in  data:
        #print(datai, any(char.isalpha() for char in datai))
        if isinstance(datai, int) or datai in ['THRU', 'EXCEPT', 'BY'] or isinteger(datai):
            continue
        #elif datai.isn
        elif '.' in datai:
            stype = 'float'
            break
        elif '/' in datai or any(char.isalpha() for char in datai):
            stype = 'str'
            break
        else:
            raise RuntimeError('datai=%r data=%s' % (datai, data))

    #if stype in ['float', 'str']:
        #print(data)
    return data, stype

def isinteger(astring):
    """ Is the given string an integer? """
    try:
        int(astring)
    except ValueError:
        return False
    else:
        return True

def get_except(data, i, ndata, end_value):
    removed = set()
    while i < len(data):
        value = data[i]
        #print('  exclude?', i, value)
        if isinstance(value, int):
            ivalue = value
        else:
            ivalue = int(value)
        #print('  ivalue=%s > end_value=%s' %  (ivalue, end_value))
        if ivalue > end_value:
            #print('    break')
            break
        removed.add(ivalue)
        i += 1
    #print('removed...', removed)
    #print('*datai =', data[i])
    return i, removed

def expand(data_in):
    """new expand method

    Odd behavior:
     - 1 THRU 10 EXCEPT 4,5,6,11 results in [1,2,3,7,8,9,10,11]
       the EXCEPT is active while the value is less than the THRU value of 10

    What happens when:
     - a value is excluded and later added or
     - a value is included and later excluded
    """
    data, stype = setup_data(data_in)
    if len(data) == 1 and stype == 'str':
        return data

    #print('***************************')
    #print('data =', data)
    if stype in ['float', 'str']:
        raise NotImplementedError(data_in)

    assert stype == 'int', data
    ndata = len(data) - 1
    out = []
    removed_set = set()

    i = 0
    is_thru = False
    while i < len(data):
        value = data[i]
        #print('value =', value)
        if isinstance(value, int):
            out.append(value)
            i += 1
            continue

        if is_thru:
            is_thru = False
            continue_after_range = False

            out2 = []
            #print('*', i, value)
            value0 = out.pop()
            end_value = int(value)
            i += 1

            by_value = 1
            if i < ndata:
                next_svalue = data[i]
                #print('found by?', next_svalue)
                if next_svalue.isdigit():  # catches integers
                    continue_after_range = True
                elif next_svalue == 'BY':
                    i += 1
                    next_svalue2 = data[i]
                    if isinstance(next_svalue2, int):
                        by_value = next_svalue2
                    else:
                        next_svalue2 = next_svalue2.upper()
                        #print('next_svalue2 =', next_svalue2)
                        by_value = int(next_svalue2)
                    #print('by_value =', by_value)
                elif next_svalue == 'EXCEPT':
                    i, removed_seti = get_except(data, i+1, ndata, end_value)
                    removed_set.update(removed_seti)
                    continue_after_range = True
                else:
                    raise RuntimeError('next_svalue=%s data=%s' % (next_svalue, data))
            rangei = range(value0, end_value + 1, by_value)
            #print('range(%r, %r, %r) = %s' % (value0, end_value + 1, by_value, rangei))
            out2.extend(rangei)
            if end_value != out2[-1]:
                out2.append(end_value)

            if continue_after_range:
                #print('continue_after_range')
                out.extend(out2)
                continue

            i += 1
            if i < ndata:
                exclude_svalue = data[i]
                #print('found exclude?', exclude_svalue)
                if exclude_svalue == 'EXCEPT':
                    i, removed = get_except(data, i, ndata, end_value)
                    removed_set.update(removed_seti)
                    asdf
                    continue
                else:
                    not_except
            else:
                i -= 1

            out.extend(out2)
            #print(i, value0, 'THRU', end_value)
            #print('------------')
            del by_value
            i += 1
        else:
            #print('add', i, value)
            assert is_thru is False, data
            if value == 'THRU':
                is_thru = True
            else:
                ivalue = int(value)
                out.append(ivalue)
            i += 1
    if len(removed_set):
        #print('data =', data)
        #print('out =', out)
        #print('removed_set =', removed_set)
        out_set = set(out) - removed_set
        out = list(out_set)
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
                if isinstance(field, str) else field
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
