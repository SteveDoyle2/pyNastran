"""
defines:
 - collapse_thru_by(fields, get_packs=False)
 - collapse_thru_by_float(fields)
 - collapse_thru(fields, nthru=None)
 - collapse_thru_packs(fields)
 - collapse_colon_packs(fields, thru_split=3)
 - condense(value_list)
 - build_thru_packs(packs, max_dv=1, thru_split=3)
 - build_thru(packs, max_dv=None, nthru=None)
 - build_thru_float(packs, max_dv=None)

"""
from collections import Counter
from typing import List, Union, Optional
import numpy as np


def collapse_thru_by(fields: List[int], get_packs: bool=False) -> List[List[int]]:
    """
    Parameters
    ----------
    fields : List[int]
        the list of fields to collapse
    get_packs : bool; default=False
        get the list of packs so "special" formatting can be done

    fields              packs
    [1, 2, 3...150]  -> [1, 150, 1]
    [1, 3, 5...150]  -> [1, 150, 2]

    """
    _check_sort_fields(fields)
    packs = condense(fields)
    if get_packs:
        return packs
    fields2 = build_thru(packs)
    #assert fields == expand_thru_by(fields2)  # why doesn't this work?
    return fields2


def collapse_thru_by_float(fields: List) -> List[Union[int, str]]:
    _check_sort_fields(fields)
    packs = condense(fields)
    fields2 = build_thru_float(packs)
    #assert fields == expand_thru_by(fields2)  # why doesn't this work?
    return fields2


def collapse_thru(fields, nthru=None):
    """
    Collapses fields into a set of packs

    Parameters
    ----------
    fields : list[int, int, ...]
        the list of integers to compress

    Returns
    -------
    packs = list[pack]
       pack = list[int first_val, int last_val, int_by]

    """
    _check_sort_fields(fields)
    packs = condense(fields)
    fields2 = build_thru(packs, max_dv=1) # , nthru=nthru
    if nthru is not None and Counter(fields2)['THRU'] > 2:
        return fields
    #assert fields == expand_thru_by(fields2), fields2  # why doesn't this work?
    return fields2

def _check_sort_fields(fields):
    if isinstance(fields, np.ndarray):
        # if an int array, we don't need to test it beyond that
        if not isinstance(fields[0], (np.int32, np.int64)):  # pragma: no cover
            raise NotImplementedError(f'fields={fields}; dtype={fields.dtype}')
        #'THRU' in fields
        #assert not np.any(''), fields
    elif isinstance(fields, (list, tuple)):
        assert 'THRU' not in fields, fields
    else:  # pragma: no cover
        raise NotImplementedError(f'fields={fields}; type={type(fields)}')
    fields.sort()

def collapse_thru_packs(fields):
    assert isinstance(fields, list), fields
    assert 'THRU' not in fields, fields
    packs = condense(fields)
    singles, doubles = build_thru_packs(packs, max_dv=1)

    #assert fields == expand_thru_by(fields2), fields2  # why doesn't this work?
    return singles, doubles


def collapse_colon_packs(fields, thru_split=3):
    """
    Applies colons (:) to packs to represent THRU and BY as is used by
    Patran.

    Parameters
    ----------
    fields : List[int]
        the values to collapse
    thru_split : int; default=3
        the length to not write THRU
        3 : [10, 11, 12] will write as '10 THRU 12'
        4 : [10, 11, 12] will write as '10 11 12'

    Returns
    -------
    singles : List[int]
        the list of singles
    doubles : List[pack]
        pack : List[varies]
            [3, :, 13]
            [3, :, 13, :, 5]
        the double packs

    # invalid
    SET1,4000, 1, 3, :, 10, 20, :, 30

    # valid
    SET1,4000, 1
    SET1,4000, 3,  :, 10
    SET1,4000, 20, :, 30

    # output
    singles = [1]
    doubles = [[3, ':', 10], [20, ':', 30]]

    """
    packs = condense(fields)
    singles, doubles = build_thru_packs(packs, max_dv=None, thru_split=thru_split)
    doubles2 = []
    for double in doubles:
        if len(double) == 3:
            double[1] = ':'
        elif len(double) == 5:
            double[1] = ':'
            double[3] = ':'
        else:
            raise RuntimeError(double)
        doubles2.append(double)
    return singles, doubles2


def condense(value_list):
    """
    Builds a list of packs (list of 3 values representing the first, last,
    and delta values for condensing a SET card.

    Parameters
    ----------
    value_list : List[int]
        list of values to pack

    Returns
    -------
    packs : List[pack]
       pack : List[id_low, id_high, delta_id]
           a list representation of the min/max/delta id values

    .. seealso:: build_thru

    """
    if len(value_list) == 0:
        return []
    if len(value_list) == 1:
        return [[value_list[0], value_list[0], 1]]
    value_list.sort()
    packs = []

    dv_old = None
    first_val = value_list[0]
    last_val = first_val

    for val in value_list[1:]:
        try:
            dv = val - last_val
        except TypeError:
            print("last_val=%r val=%r" % (last_val, val))
            print("value_list=%r" % value_list)
            raise

        # sets up the first item of the pack
        if dv_old is None:
            dv_old = dv

        # fill up the pack
        if dv_old == dv:
            last_val = val
        else:
            packs.append([first_val, last_val, dv_old])
            last_val = val
            dv_old = None
            first_val = val

    # fills the last pack
    if dv_old == dv:
        packs.append([first_val, val, dv])
    else:
        packs.append([first_val, val, dv_old])
    return packs


def build_thru_packs(packs, max_dv=1, thru_split=3):
    """
    Applies THRU and BY to packs to shorten output as Nastran does on
    cards like the SPOINT

    Parameters
    ----------
    packs : List[pack]
       pack : List[id_low, id_high, delta_id]
           a list representation of the min/max/delta id values
    max_dv : int; default=1
        maximum allowed delta between two values
    thru_split : int; default=3
        the length to not write THRU
        3 : [10, 11, 12] will write as '10 THRU 12'
        4 : [10, 11, 12] will write as '10 11 12'

    Returns
    -------
    singles : List[int]
        the list of singles
    doubles : List[pack]
        pack : List[varies]
            [3, THRU, 13]
            [3, THRU, 13, BY, 5]
        the double packs

    # invalid
    SET1,4000, 1, 3, THRU, 10, 20, THRU, 30

    # valid
    SET1,4000, 1
    SET1,4000, 3,  THRU, 10
    SET1,4000, 20, THRU, 30

    returns
      singles = [1]
      doubles = [[3, 'THRU', 10], [20, 'THRU', 30]]

    """
    singles = []
    doubles = []
    for (first_val, last_val, by) in packs:
        if first_val == last_val:
            singles.append(first_val)
        else:
            if by == 1:
                if last_val - first_val < thru_split: # dont make extra THRU cards
                    singlei = list(range(first_val, last_val + 1, 1))
                    singles += singlei
                else:
                    double = [first_val, 'THRU', last_val]
                    doubles.append(double)
            else:
                diff = last_val - first_val
                if max_dv == 1 or diff == by:
                    singlei = list(range(first_val, last_val + by, by))
                    singles += singlei
                else:
                    double = [first_val, 'THRU', last_val, 'BY', by]
                    doubles.append(double)
    return singles, doubles


def build_thru(packs, max_dv: Optional[int]=None, nthru: Optional[int]=None):
    """
    Takes a pack [1,7,2] and converts it into fields used by a SET card.
    The values correspond to the first value, last value, and delta in the
    list.  This means that [1,1001,2] represents 500 values.
    [1,1001,1] represents 1001 values and will be written as [1,THRU,1001]..

    Parameters
    ----------
    packs : List[pack]
        pack : List[int first, int last, int delta]
            the first, last, delta id values
    max_dv : int; default=None -> no limit
        defines the max allowable delta between two values
    nthru : int; default=None
        don't use this; it will crash

    Returns
    -------
    value : varies
        the value of the field

    """
    #singles = []
    fields = []
    if nthru is not None:
        raise NotImplementedError('nthru=%s' % nthru)
        #assert nthru == 1, nthru # others
        #assert nthru == 1, nthru
        #packs2 = []
        #nvalues = []
        #for (first_val, last_val, dv) in packs:
            #nvalue = (last_val - first_val + 1) // dv
            ##print('first=%s last=%s delta=%s dv=%s n=%s' % (
                ##first_val, last_val, last_val-first_val, dv, nvalue))
            #nvalues.append(nvalue)
        #i = nvalues.index(max(nvalues))
        ##print('nvalues =', nvalues, i)
        #packs = []

    for (first_val, last_val, dv) in packs:
        if first_val == last_val:
            fields.append(first_val)
        elif dv == 1:
            if last_val - first_val > 2:
                fields.extend([first_val, 'THRU', last_val])
            elif last_val - first_val == 2:
                # no point in writing 'A THRU A+2'
                fields.extend([first_val, first_val + 1, last_val])
            else:
                # no point in writing 'A THRU A+1'
                fields.append(first_val)
                fields.append(last_val)
        else:
            if max_dv is None:
                # no point in writing 'A THRU B BY C'
                if last_val - first_val > 4 * dv:
                    fields.extend([first_val, 'THRU', last_val, 'BY', dv])
                else:
                    fields += list(range(first_val, last_val + dv, dv))
            else:
                for v in range(first_val, last_val + dv, dv):
                    fields.append(v)
    return fields


def build_thru_float(packs: List[List[int]], max_dv: Optional[int]=None) -> List[Union[int, str]]:
    """
    Takes a pack [1,7,2] and converts it into fields used by a SET card.
    The values correspond to the first value, last value, and delta in the
    list.  This means that [1,1001,2] represents 500 values.
    [1,1001,1] represents 1001 values and will be written as [1,THRU,1001].

    Parameters
    ----------
    packs : List[pack]
        pack : List[first, last, delta]
        first, last, delta are integers
    max_dv : int; default=None -> no limit
        integer defining the max allowable delta between two values
        (default=None; no limit)

    """
    fields = []
    for (first_val, last_val, dv) in packs:
        if last_val - first_val > 4 * dv:
            fields.extend([first_val, 'THRU', last_val, 'BY', dv])
        else:
            nv = int(round((last_val - first_val) / dv)) + 1
            for i in range(nv):
                v = first_val + i * dv
                fields.append(v)
    return fields
