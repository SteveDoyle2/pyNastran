"""
defines:
 - expand_thru_case_control(set_value)
 - write_set(set_id, values, spaces='')
 - write_stress_type(key, options, value, spaces='')
"""
from __future__ import print_function
from typing import List, Optional, Union, Set, Any
from six import string_types
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.collpase_card import collapse_thru_packs
from pyNastran.bdf.bdf_interface.assign_type import interpret_value

def expand_thru_int(set_value):   # pragma: no cover
    """
    27,35,25,41234,123,thru,134,9701,9901
    1,thru,9,by,2
    9,THRU,19,EXCEPT,12
    0.1 0.3 0.5 1.0 3.0 5.0 10.0 14.0
    """
    #packs = []
    assert '/' not in set_value, set_value
    values = ','.join(set_value.split()).split(',')
    nvalues = len(values)
    i = 0
    while i < nvalues:
        start = int(values[i])
        # A, THRU, B, BY, C
        if i < (nvalues - 2):
            values.append(start)
            i += 1
            continue

        thru_value = values[i + 1].upper()
        if thru_value != 'THRU':
            values.append(start)
            i += 1
            continue

        # there is a thru
        stop = int(values[i + 2])
        if i < (nvalues - 4):
            values += range(start, stop, 1)
            i += 3
            continue


        by_value = values[i + 3].upper()
        if by_value != 'BY':
            values += range(start, stop, 1)
            i += 3
            continue

        by = int(values[i+4])
        values += range(start, stop, by)
        values.append(by)
        i += 4


def expand_thru_case_control(set_value):
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
    **This hasn't been verified**
    #set_value = 'SET 88 = 1 THRU 5, 20 THRU 30 BY 2'
    set_value = ['1 THRU 5', '20 THRU 30 BY 2']
    >>> values = expand_thru_case_control(set_value)
    [1, 2, 3, 4, 5, 20, 22, 24, 26, 28, 30]
    """
    set_value2 = set()
    add_mode = True
    imin_int = 0
    imax_int = 0
    #print('set_value = %r' % set_value)
    for ivalue in set_value:
        #print('----------')
        if isinstance(ivalue, integer_types):
            assert add_mode is True, add_mode
            set_value2.add(ivalue)
            continue

        ivalue_str = ivalue.strip()
        #print('  ivalue=%r; type=%s add=%s' % (ivalue, type(ivalue), add_mode))
        if '/' in ivalue_str:
            set_value2.add(ivalue_str)
        else:
            #if 'ALL' in ivalue:
                #msg = ('ALL is not supported on CaseControlDeck '
                       #'SET card\nvalue=%r\nset=%r' % (ivalue, set_value))
                #raise RuntimeError(msg)
            #elif 'EXCEPT' in ivalue:
                #msg = ('EXCEPT is not supported on CaseControlDeck '
                       #'SET card\nvalue=%r\nset=%r' % (ivalue, set_value))
                #raise RuntimeError(msg)

            ivalue2 = interpret_value(ivalue_str, card=str(set_value)) #  type: Optional[Union[int, float, str]]
            if isinstance(ivalue2, integer_types):
                #print('  isdigit')
                #print('  imin=%s ivalue2=%s imax=%s' % (imin_int, ivalue2, imax_int))
                #ivalue = int(ivalue)
                if not imin_int < ivalue2 < imax_int:
                    add_mode = True
                    #print('  break...\n')
                if add_mode:
                    #print('  adding %s' % ivalue)
                    set_value2.add(ivalue2)
                else:
                    #print('  removing %s' % ivalue2)
                    set_value2.remove(ivalue2)
                    imin_int = int(ivalue2)
            elif isinstance(ivalue2, float):
                assert add_mode is True, add_mode
                set_value2.add(ivalue2)
            elif isinstance(ivalue2, string_types):
                #print('  not digit=%r' % set_value)
                if set_value == 'EXCLUDE':
                    msg = ('EXCLUDE is not supported on CaseControlDeck '
                           'SET card\n')
                    raise RuntimeError(msg)
                elif 'THRU' in ivalue2:
                    set_valuesi, imin_int, imax_int, add_mode = _expand_thru_case_control_string_thru(set_value, ivalue2, add_mode)
                    set_value2.update(set_valuesi)
                else:
                    assert add_mode is True, add_mode
                    #print('else...', ivalue)
                    set_value2.add(ivalue.strip())
            else:
                raise NotImplementedError(ivalue2)

    list_values = list(set_value2)
    try:
        list_values.sort()
    except TypeError:
        msg = 'sort error: list_values=%s set_value2=%s'  % (list_values, set_value2)
        raise TypeError(msg)

    #print('end of expansion = %s' % list_values)
    return list_values

def _expand_thru_case_control_string_thru(set_value, svalue, add_mode):
    # type: (Any, str, bool) -> Any
    """helper for ``expand_thru_case_control``"""
    ## type: (List[str], Any, bool) -> (Set[int], bool)
    set_values_out = set()  # type: Set[int]
    svalues = svalue.split()  # type: List[str]
    #print('  set_value=%s' %  set_value)
    #print('  svalues=%s' % svalues)
    if len(svalues) == 3:
        assert add_mode is True, add_mode
        imin, thru, imax = svalues
        assert thru == 'THRU', thru
        imin_int = int(imin)
        imax_int = int(imax)
        assert imax_int > imin_int, 'imin=%s imax=%s' % (imin_int, imax_int)
        for jthru in range(imin_int, imax_int + 1):
            set_values_out.add(jthru)

    elif len(svalues) == 4:
        imin, thru, imax, by_except = svalues
        imin_int = int(imin)
        imax_int = int(imax)
        assert imax_int > imin_int, 'imin=%s imax=%s' % (imin_int, imax_int)
        assert by_except == 'EXCEPT', by_except

        for jthru in range(imin_int, imax_int + 1):
            set_values_out.add(jthru)
        add_mode = False

    elif len(svalues) == 5:
        imin, thru, imax, by_except, increment_except = svalues
        imin_int = int(imin)
        imax_int = int(imax)
        assert imax_int > imin_int, 'imin=%s imax=%s' % (imin_int, imax_int)
        increment_except_int = int(increment_except)
        if by_except == 'BY':
            for jthru in range(imin_int, imax_int + 1, increment_except_int):
                set_values_out.add(jthru)
            add_mode = True
        elif by_except == 'EXCEPT':
            for jthru in range(imin_int, imax_int + 1):
                if jthru == increment_except_int:
                    continue
                set_values_out.add(jthru)
            add_mode = False
        else:
            raise RuntimeError(svalue)
    else:
        msg = ('expected data of the form: '
               '"10 THRU 20" or "10 THRU 20 BY 5"\n'
               'actual=%r; input=%s' % (svalue.strip(), set_value))
        raise RuntimeError(msg)
    #print('  set_values_out =',  set_values_out)
    #print('---')
    return set_values_out, imin_int, imax_int, add_mode


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


def write_set(set_id, values, spaces=''):
    # type: (List[int], List[int], str) -> str
    """
    writes
    SET 80 = 3926, 3927, 3928, 4141, 4142, 4143, 4356, 4357, 4358, 4571,
         4572, 4573, 3323 THRU 3462, 3464 THRU 3603, 3605 THRU 3683,
         3910 THRU 3921, 4125 THRU 4136, 4340 THRU 4351

    Parameters
    ----------
    values : List[int]
        the Set values
    options : int / str; default=''
        the Set ID
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
    starter = 'SET %s = ' % (set_id)
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
