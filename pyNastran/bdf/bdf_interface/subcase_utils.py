"""
defines:
 - expand_thru_case_control(set_value)
 - write_set(set_id, values, spaces='')
 - write_stress_type(key, options, value, spaces='')

"""
from typing import List, Optional, Union, Set, Any
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.collpase_card import collapse_thru_packs


def expand_thru_case_control(data_in):
    # type: (List[Union[int, float, str]]) -> List[int]
    """
    Expands a case control SET card

    Parameters
    ----------
    set_value : str???
        ???

    Returns
    -------
    values : List[int] / List[float]
        the values in the SET

    Examples
    --------
    set_value = ['1 THRU 5', '20 THRU 30 BY 2']
    >>> values = expand_thru_case_control(set_value)
    [1, 2, 3, 4, 5, 20, 22, 24, 26, 28, 30]

    Odd behavior:
     - 1 THRU 10 EXCEPT 4,5,6,11 results in [1,2,3,7,8,9,10,11]
       the EXCEPT is active while the value is less than the THRU value of 10
       and the values are in ascending order

    What happens when:
     - a value is excluded and later added or
     - a value is included and later excluded

    """
    data, stype = setup_data(data_in)
    if len(data) == 1 and stype == 'str':
        return data

    #print('***************************')
    #print('data =', data)
    if stype == 'int':
        out = [datai if isinstance(datai, int) else int(datai) for datai in data]
        out.sort()
        return out
    elif stype == 'float':
        return expand_float(data)
    elif stype == 'str':
        return data

    #-------------------------------------------------------------------
    #  has a THRU, BY, or EXCEPT

    assert stype == 'int_thru', data
    ndata = len(data) - 1
    out = []
    removed_set = set()

    i = 0
    is_thru = False
    while i < len(data):
        value = data[i]
        #print('---')
        #print('i=%s value=%r is_thru=%s' % (i, value, is_thru))

        if is_thru:
            #print('$$$ starting thru')
            is_thru = False
            continue_after_range = False

            out2 = []
            #print('*', i, value)
            value0 = out.pop()
            end_value = int(value)
            i += 1
            #print('*value0 =', value0)
            #print('*end_value =', i, end_value)

            by_value = 1
            if i <= ndata:
                next_svalue = data[i]
                #print('found by?', next_svalue)
                if isinstance(next_svalue, int) or next_svalue.isdigit():  # catches integers
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
            #else:
                # out of range

            rangei = range(value0, end_value + 1, by_value)
            #print('range(%r, %r, %r) = %s' % (value0, end_value + 1, by_value, rangei))
            out2.extend(rangei)
            if end_value != out2[-1]:
                out2.append(end_value)

            if continue_after_range:
                #print('continue_after_range; adding %s' % str(out2))
                out.extend(out2)
                continue

            #print('looking for except...')
            i += 1
            if i < ndata:
                exclude_svalue = data[i]
                #print('found exclude?', exclude_svalue)
                if exclude_svalue == 'EXCEPT':
                    i, removed_seti = get_except(data, i, ndata, end_value)
                    removed_set.update(removed_seti)
                    raise RuntimeError('need a test problem for EXCEPT')
                    #continue
                else:
                    raise RuntimeError('need a test problem for skipping EXCEPT; '
                                       'datai=%s' %  exclude_svalue)
            else:
                i -= 1

            out.extend(out2)
            #print(i, value0, 'THRU', end_value)
            #print('------------')
            del by_value
            i += 1

        elif isinstance(value, int):
            #print('value =', i, value)
            out.append(value)
            i += 1
            continue
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

def split_comma_space(datai):
    """normalizes the form of a SET into a comma separated list instead of being mixed"""
    if ',' in datai:
        # do I need a strip?
        sline = datai.replace(',', ' ').split()
    else:
        sline = datai.split()
    return sline

def setup_data(data_in):
    """helper method for ``expand_thru_case_control``"""
    data = []
    for datai in data_in:
        if isinstance(datai, (int, float)):
            data.append(datai)
        else:
            datai = datai.upper().replace('INCLUDE', ' ')
            sline = split_comma_space(datai)
            data.extend(sline)

    stype = 'int'
    for datai in  data:
        #print(datai, any(char.isalpha() for char in datai))
        if isinstance(datai, int) or isinteger(datai):
            continue
        elif datai in ['THRU', 'EXCEPT', 'BY']:
            stype = 'int_thru'
            break
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
    """Is the given string an integer?"""
    try:
        int(astring)
    except ValueError:
        return False
    else:
        return True

def get_except(data, i, ndata, end_value):
    """helper method for expand that gets the values until the end of an
    EXCEPT chain"""
    removed = []
    ivalue_old = 0
    while i < len(data):
        value = data[i]
        #print('  exclude?', i, value)
        if isinstance(value, int):
            ivalue = value
        elif value.isdigit():
            ivalue = int(value)
        elif value == 'THRU':
            ivalue_old = int(data[i-1])
            ivalue_new = int(data[i+1])
            rangei = range(ivalue_old, ivalue_new+1)
            removed.extend(rangei)
            i += 2
            continue
        else:
            raise NotImplementedError('data[%i] = %s; data=%s' % (i, value, data))
        #print('  ivalue=%s > end_value=%s' %  (ivalue, end_value))

        if ivalue > end_value and ivalue > ivalue_old:
            #print('    break')
            break
        removed.append(ivalue)
        ivalue_old = ivalue
        i += 1

    removed_set = set(removed)
    #print('removed...', removed)
    #print('*datai =', data[i])
    return i, removed_set


def expand_float(data):
    """helper method for ``expand_thru_case_control``"""
    out = []
    for datai in  data:
        out.append(float(datai))
    out.sort()
    return out


def write_stress_type(key, options, value, spaces=''):
    # type: (str, List[str], Optional[str], str) -> str
    """
    writes:
     - STRESS(SORT1) = ALL
     - GROUNDCHECK(PRINT,SET=(G,N,N+AUTOSPC,F,A),THRESH=1e-2,DATAREC=NO) = YES

    """
    msg = ''
    str_options = ','.join(options)
    #print("str_options = %r" % (str_options))
    #print("STRESS-type key=%s value=%s options=%s"
          #% (key, value, options))
    if value is None:
        val = ''
    else:
        val = ' = %s' % value
    if len(str_options) > 0:
        msgi = '%s(%s)%s\n' % (key, str_options, val)
        if len(msgi) > 64:
            msgi = '%s(' % key
            msg_done = ''
            i = 0
            while i < len(options):
                new_msgi = '%s,' % options[i]
                if (len(msgi) + len(new_msgi)) < 64:
                    msgi += new_msgi
                else:
                    msg_done += msgi + '\n'
                    msgi = spaces + new_msgi
                i += 1
            msg_done += msgi
            msgi = ''
            msg_done = msg_done.rstrip(' ,\n') + ')%s\n' % val
            assert len(msgi) < 68, 'len(msg)=%s; msg=\n%s' % (len(msgi), msgi)
            msg += msg_done
        else:
            assert len(msgi) < 68, 'len(msg)=%s; msg=\n%s' % (len(msgi), msgi)
            msg += msgi
    else:
        msgi = '%s%s\n' % (key, val)
        assert len(msgi) < 68, 'len(msg)=%s; msg=\n%s' % (len(msgi), msgi)
        msg += msgi
    msg = spaces + msg
    return msg


def write_set(set_id: int, values: List[int], spaces: str='') -> str:
    """
    writes
    SET 80 = 3926, 3927, 3928, 4141, 4142, 4143, 4356, 4357, 4358, 4571,
         4572, 4573, 3323 THRU 3462, 3464 THRU 3603, 3605 THRU 3683,
         3910 THRU 3921, 4125 THRU 4136, 4340 THRU 4351

    Parameters
    ----------
    set_id : int  / str?
        the Set ID
    values : List[int]
        the Set values
    spaces : str; default=''
        indentation

    Returns
    -------
    msg : str
       the string of the set

    Examples
    --------
    **Example 1**
    >>> set_id = 80
    >>> values = [1, 2, 3, 4, 5, 7]
    >>> set = write_set(set_id, values, spaces='')
    >>> set
    SET 80 = 1 THRU 5, 7

    **Example 2**
    >>> set_id = ''
    >>> values = ['ALL']
    >>> set = write_set(set_id, values, spaces='')
    >>> set
    SET = ALL

    """
    values.sort()
    starter = f'SET {set_id:d} = '
    msg2 = spaces + starter

    msg = ''
    nchars = len(msg2)
    is_valid = True
    for value in values:
        if not isinstance(value, integer_types):
            is_valid = False
            break

    if is_valid:
        singles, doubles = collapse_thru_packs(values)

        out_values = singles
        for double in doubles:
            assert len(double) == 3, double
            sdouble = '%i THRU %i' % (double[0], double[2])
            out_values.append(sdouble)
    else:
        out_values = values

    for out_value in out_values:
        new_string = '%s, ' % out_value
        if len(msg2 + new_string) > 70:
            msg += msg2 + '\n'
            msg2 = ' ' * nchars + new_string
        else:
            msg2 += new_string
    return msg + msg2.rstrip(' \n,') + '\n'
