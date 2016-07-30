from __future__ import print_function


def collapse_thru_by(fields, get_packs=False):
    """
    :param fields:    the list of fields to collapse
    :param get_packs: get the list of packs so "special" formatting can be done

    fields              packs
    [1, 2, 3...150]  -> [1, 150, 1]
    [1, 3, 5...150]  -> [1, 150, 2]
    """
    assert 'THRU' not in fields, fields
    fields.sort()
    packs = condense(fields)
    if get_packs:
        return packs
    fields2 = build_thru(packs)
    #assert fields == expand_thru_by(fields2)  # why doesn't this work?
    return fields2


def collapse_thru_by_float(fields):
    assert 'THRU' not in fields, fields
    fields.sort()
    packs = condense(fields)
    fields2 = build_thru_float(packs)
    #assert fields == expand_thru_by(fields2)  # why doesn't this work?
    return fields2


def collapse_thru(fields):
    assert 'THRU' not in fields, fields
    fields.sort()
    packs = condense(fields)
    fields2 = build_thru(packs, max_dv=1)
    #assert fields == expand_thru_by(fields2), fields2  # why doesn't this work?
    return fields2


def collapse_thru_packs(fields):
    assert 'THRU' not in fields, fields
    fields.sort()
    packs = condense(fields)
    singles, doubles = build_thru_packs(packs, max_dv=1)

    #assert fields == expand_thru_by(fields2), fields2  # why doesn't this work?
    return singles, doubles


def collapse_colon_packs(fields):
    """
    # invalid
    SET1,4000, 1, 3, :, 10, 20, :, 30

    # valid
    SET1,4000, 1
    SET1,4000, 3,  :, 10
    SET1,4000, 20, :, 30

    Returns
    -------
      singles = [1]
      doubles = [[3, ':', 10], [20, ':', 30]]
    """
    fields.sort()
    packs = condense(fields)
    singles, doubles = build_thru_packs(packs, max_dv=None)
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


def build_thru_packs(packs, max_dv=1):
    """
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
                if last_val - first_val < 3: # dont make extra THRU cards
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


def build_thru(packs, max_dv=None):
    """
    Takes a pack [1,7,2] and converts it into fields used by a SET card.
    The values correspond to the first value, last value, and delta in the
    list.  This means that [1,1001,2] represents 500 values.
    [1,1001,1] represents 1001 values and will be written as [1,THRU,1001]..

    :param packs: list of packs (list of 3 values: [first, last, delta] )
    :param maxDV: integer defining the max allowable delta between two values
            (default=None; no limit)
    """
    fields = []
    for (first_val, last_val, dv) in packs:
        if first_val == last_val:
            fields.append(first_val)
        elif dv == 1:
            if last_val - first_val > 2:
                fields.append(first_val)
                fields.append('THRU')
                fields.append(last_val)
            elif last_val - first_val == 2:
                fields.append(first_val)
                fields.append(first_val + 1)
                fields.append(last_val)
            else:
                fields.append(first_val)
                fields.append(last_val)
        else:
            if max_dv is None:
                if last_val - first_val > 4 * dv:
                    fields.append(first_val)
                    fields.append('THRU')
                    fields.append(last_val)
                    fields.append('BY')
                    fields.append(dv)
                else:
                    fields += list(range(first_val, last_val + dv, dv))
            else:
                for v in range(first_val, last_val + dv, dv):
                    fields.append(v)
    return fields


def build_thru_float(packs, max_dv=None):
    """
    Takes a pack [1,7,2] and converts it into fields used by a SET card.
    The values correspond to the first value, last value, and delta in the
    list.  This means that [1,1001,2] represents 500 values.
    [1,1001,1] represents 1001 values and will be written as [1,THRU,1001].

    Parameters
    ----------
    packs : List[ List[int, int, int], ... ]
        list of packs
        pack : List[int, int, int]
            list of [first, last, delta]
    max_dv : int/None; default=None
        integer defining the max allowable delta between two values
        (default=None; no limit)
    """
    fields = []
    for (first_val, last_val, dv) in packs:
        if last_val - first_val > 4 * dv:
            fields.append(first_val)
            fields.append('THRU')
            fields.append(last_val)
            fields.append('BY')
            fields.append(dv)
        else:
            nv = int(round((last_val - first_val) / dv)) + 1
            for i in range(nv):
                v = first_val + i * dv
                fields.append(v)
    return fields
