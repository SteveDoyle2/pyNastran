from __future__ import print_function
from six import string_types, iteritems
from pyNastran.utils import integer_types
from pyNastran.bdf.cards.collpase_card import collapse_thru_packs
from pyNastran.bdf.bdf_interface.assign_type import interpret_value

def expand_thru_case_control(set_value):
    set_value2 = set([])
    add_mode = True
    imin = 0
    imax = 0
    #print('set_value = %r' % set_value)
    for ivalue in set_value:
        if isinstance(ivalue, integer_types):
            assert add_mode is True, add_mode
            set_value2.add(ivalue)
            continue
        ivalue = ivalue.strip()
        #print('  ivalue=%r; type=%s' % (ivalue, type(ivalue)))
        if '/' in ivalue:
            set_value2.add(ivalue)
        else:

            #if 'ALL' in ivalue:
                #msg = ('ALL is not supported on CaseControlDeck '
                       #'SET card\nvalue=%r\nset=%r' % (ivalue, set_value))
                #raise RuntimeError(msg)
            #elif 'EXCEPT' in ivalue:
                #msg = ('EXCEPT is not supported on CaseControlDeck '
                       #'SET card\nvalue=%r\nset=%r' % (ivalue, set_value))
                #raise RuntimeError(msg)

            ivalue = interpret_value(ivalue, card=str(set_value))
            if isinstance(ivalue, integer_types):
                #print('  isdigit')
                #ivalue = int(ivalue)
                if not imin < ivalue < imax:
                    add_mode = True
                    #print('  break...\n')
                if add_mode:
                    #print('  adding %s' % ivalue)
                    set_value2.add(ivalue)
                else:
                    #print('  removing %s' % ivalue)
                    set_value2.remove(ivalue)
                    imin = ivalue
            elif isinstance(ivalue, float):
                assert add_mode is True, add_mode
                set_value2.add(ivalue)
            elif isinstance(ivalue, string_types):
                #print('  not digit=%r' % set_value)
                if set_value is 'EXCLUDE':
                    msg = ('EXCLUDE is not supported on CaseControlDeck '
                           'SET card\n')
                    raise RuntimeError(msg)
                elif 'THRU' in ivalue:
                    svalues = ivalue.split()
                    #print('  svalues=%s' % svalues)
                    if len(svalues) == 3:
                        assert add_mode is True, add_mode
                        imin, thru, imax = svalues
                        assert thru == 'THRU', thru
                        imin = int(imin)
                        imax = int(imax)
                        assert imax > imin, 'imin=%s imax=%s' % (imin, imax)
                        for jthru in range(imin, imax + 1):
                            set_value2.add(jthru)

                    elif len(svalues) == 4:
                        imin, thru, imax, by_except = svalues
                        imin = int(imin)
                        imax = int(imax)
                        assert imax > imin, 'imin=%s imax=%s' % (imin, imax)
                        assert by_except == 'EXCEPT', by_except

                        for jthru in range(imin, imax + 1):
                            set_value2.add(jthru)
                        add_mode = False

                    elif len(svalues) == 5:
                        imin, thru, imax, by_except, increment_except = svalues
                        imin = int(imin)
                        imax = int(imax)
                        assert imax > imin, 'imin=%s imax=%s' % (imin, imax)
                        increment_except = int(increment_except)
                        if by_except == 'BY':
                            for jthru in range(imin, imax + 1, by_except):
                                set_value2.add(jthru)
                            add_mode = True
                        elif by_except == 'EXCEPT':
                            for jthru in range(imin, imax + 1):
                                if jthru == increment_except:
                                    continue
                                set_value2.add(jthru)
                            add_mode = False
                        else:
                            raise RuntimeError(ivalue)
                    else:
                        msg = ('expected data of the form: '
                               '"10 THRU 20" or "10 THRU 20 BY 5"\n'
                               'actual=%r; input=%s' % (ivalue.strip(), set_value))
                        raise RuntimeError(msg)
                else:
                    assert add_mode is True, add_mode
                    set_value2.add(ivalue)

    list_values = list(set_value2)
    try:
        list_values.sort()
    except TypeError:
        msg = 'sort error: list_values=%s'  % (list_values)
        raise TypeError(msg)

    #print('end of expansion = %s' % list_values)
    return list_values

def write_stress_type(key, options, value, spaces):
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
                option = options[i]
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


def write_set(value, options, spaces=''):
    """
    writes
    SET 80 = 3926, 3927, 3928, 4141, 4142, 4143, 4356, 4357, 4358, 4571,
         4572, 4573, 3323 THRU 3462, 3464 THRU 3603, 3605 THRU 3683,
         3910 THRU 3921, 4125 THRU 4136, 4340 THRU 4351

    Parameters
    ----------
    value : List[int]
        the Set values
    options : int
        the Set ID
    spaces : str; default=''
        indentation

    Returns
    -------
    msg : str
       the string of the set

    Example
    -------
    value = 80
    options = [1, 2, 3, 4, 5, 7]
    set = write_set(value, options, spaces='')
    print(set)
    >>> SET 80 = 1 THRU 5, 7
    """
    value.sort()
    starter = 'SET %s = ' % (options)
    msg2 = spaces + starter

    msg = ''
    nchars = len(msg2)
    is_valid = True
    for valuei in value:
        if not isinstance(valuei, integer_types):
            is_valid = False
            break

    if is_valid:
        singles, doubles = collapse_thru_packs(value)

        out_value = singles
        for double in doubles:
            assert len(double) == 3, double
            sdouble = '%i THRU %i' % (double[0], double[2])
            out_value.append(sdouble)
    else:
        out_value = value

    for i, out_valuei in enumerate(out_value):
        new_string = '%s, ' % out_valuei
        if len(msg2 + new_string) > 70:
            msg += msg2 + '\n'
            msg2 = ' ' * nchars + new_string
        else:
            msg2 += new_string
    return msg + msg2.rstrip(' \n,') + '\n'
