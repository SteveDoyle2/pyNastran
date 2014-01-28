import sys
from pyNastran.bdf.bdfInterface.BDF_Card import wipe_empty_fields

def print_scientific_16(value):
    """
    Prints a value in 16-character scientific notation.
    This is a sub-method and shouldnt typically be called

    .. seealso:: print_float_16 for a better method
    """
    pythonValue = '%16.14e' % value  # -1.e-2
    (svalue, sExponent) = pythonValue.strip().split('e')
    exponent = int(sExponent)  # removes 0s

    if abs(value) < 0.01:
        sign = '-'
    else:
        sign = '+'
    sExp2 = str(exponent).strip('-+')  # the exponent will be added later...
    value2 = float(svalue)

    lenSExp = len(sExp2) + 1  # the plus 1 is for the sign
    leftover = 16 - lenSExp

    if value < 0:
        format = "%%1.%sf" % (leftover - 3)
    else:
        format = "%%1.%sf" % (leftover - 2)

    svalue3 = format % value2
    svalue4 = svalue3.strip('0')
    field = "%16s" % (svalue4 + sign + sExp2)
    return field


def print_float_16(value, tol=0.):
    """
    Prints a float in nastran 16-character width syntax
    using the highest precision possbile.
    .. seealso:: print_float_8
    .. warning:: completely unimplemented & untested
    """
    if abs(value) <= tol:  # tol=1e-8
        field = "%16s" % ('0.')
    else:
        if value > 0.:  # positive, not perfect...

            if value < 5e-16:  #: .. todo:: does this work properly with tol
                field = print_scientific_16(value)
                return field
            elif value < 0.001:
                field = print_scientific_16(value)
                field2 = "%16.15f" % value  # small value
                field2 = field2.strip('0 ')

                field1 = field.replace('-', 'e-')

                if field2 == '.':
                    return print_scientific_16(value)
                if len(field2) <= 16 and float(field1) == float(field2):
                    field = field2
                    field = field.strip(' 0')
            elif value < 0.1:
                field = "%16.15f" % value
            elif value < 1.:
                field = "%16.15f" % value
            elif value < 10.:
                field = "%16.14f" % value
            elif value < 100.:
                field = "%16.13f" % value   # 10 < x < 100
            elif value < 1000.:
                field = "%16.12f" % value
            elif value < 10000.:
                field = "%16.11f" % value
            elif value < 100000.:
                field = "%16.10f" % value
            elif value < 1000000.:
                field = "%16.9f" % value
            elif value < 10000000.:
                field = "%16.8f" % value
            elif value < 100000000.:
                field = "%16.7f" % value
            elif value < 1000000000.:
                field = "%16.6f" % value
            elif value < 10000000000.:
                field = "%16.5f" % value
            elif value < 100000000000.:
                field = "%16.4f" % value
            elif value < 1000000000000.:
                field = "%16.3f" % value
            elif value < 10000000000000.:
                field = "%16.2f" % value
            elif value < 100000000000000.:
                field = "%16.1f" % value
            else:  # big value  # 123456789012345.
                field = "%16.1f" % value
                if field.index('.') < 16:
                    field = '%16.1f' % (round(value))
                    field = field[0:16]  # drop off the .1f
                    assert '.' != field[0], field
                else:
                    field = print_scientific_16(value)
                return field
        else:
            if value > -5e-15:  #: .. todo:: does this work properly with tol
                field = print_scientific_16(value)
                return field
            elif value > -0.01:  # -0.001
                field = print_scientific_16(value)
                field2 = "%16.14f" % value  # small value
                field2 = field2.strip('0 ')

                field1 = '-' + field.strip(' 0-').replace('-', 'e-')  # get rid of the first minus sign, add it on afterwards

                if len(field2) <= 16 and float(field1) == float(field2):
                    field = field2.rstrip(' 0')
                    field = field.replace('-0.', '-.')
            elif value > -0.1:
                field = "%16.14f" % value   # -0.01 >x>-0.1...should be 5 (maybe scientific...)
                field = field.replace('-0.', '-.')
            elif value > -1.:
                field = "%16.14f" % value   # -0.1  >x>-1.....should be 6, but the baseline 0 is kept...
                field = field.replace('-0.', '-.')
            elif value > -10.:
                field = "%16.14f" % value   # -1    >x>-10
            elif value > -100.:
                field = "%16.13f" % value  #       -1 > x >      -10
            elif value > -1000.:
                field = "%16.12f" % value  #      -10 > x >     -100
            elif value > -10000.:
                field = "%16.11f" % value  #     -100 > x >    -1000
            elif value > -100000.:
                field = "%16.10f" % value  #    -1000 > x >   -10000
            elif value > -1000000.:
                field = "%16.9f" % value   #  -10,000 > x > -100,000
            elif value > -10000000.:
                field = "%16.8f" % value   #           -100,000 > x >          -1,000,000
            elif value > -100000000.:
                field = "%16.7f" % value   #         -1,000,000 > x >         -10,000,000
            elif value > -1000000000.:
                field = "%16.6f" % value   #        -10,000,000 > x >        -100,000,000
            elif value > -10000000000.:
                field = "%16.5f" % value   #       -100,000,000 > x >      -1,000,000,000
            elif value > -100000000000.:
                field = "%16.4f" % value   #     -1,000,000,000 > x >     -10,000,000,000
            elif value > -1000000000000.:
                field = "%16.3f" % value   #    -10,000,000,000 > x >    -100,000,000,000
            elif value > -10000000000000.:
                field = "%16.2f" % value   #   -100,000,000,000 > x >  -1,000,000,000,000
            elif value > -100000000000000.:
                field = "%16.1f" % value   # -1,000,000,000,000 > x > -10,000,000,000,000
            else:
                field = "%16.1f" % value
                if field.index('.') < 16:
                    field = '%15s.' % (int(round(value, 0)))
                    assert '.' != field[0], field
                else:
                    field = print_scientific_16(value)
                return field
        field = field.strip(' 0')
        field = '%16s' % (field)

    assert len(field) == 16, ('value=%r field=%r is not 16 characters '
                              'long, its %s' % (value, field, len(field)))
    return field


def print_field_16(value):
    """
    Prints a single 16-character width field

    :param value:   the value to print
    :returns field: an 16-character string
    """
    if isinstance(value, int):
        field = "%16s" % value
    elif isinstance(value, float):
        field = print_float_16(value)
    elif value is None:
        field = "                "
    else:
        field = "%16s" % value
    if len(field) != 16:
        msg = 'field=|%s| is not 16 characters long...rawValue=|%s|' % (field,
                                                                        value)
        raise RuntimeError(msg)
    return field


def print_card_16(fields):
    """
    Prints a nastran-style card with 16-character width fields.

    :param fields: all the fields in the BDF card (no blanks)

    .. note:: A large field format follows the  8-16-16-16-16-8 = 80
     format where the first 8 is the card name or blank (continuation).
     The last 8-character field indicates an optional continuation,
     but because it's a left-justified unneccessary field,
     print_card doesnt use it.
    """
    fields = wipe_empty_fields(fields)
    nFieldsMain = len(fields) - 1  # chop off the card name
    nBDFLines = nFieldsMain // 8
    if nFieldsMain % 8 != 0:
        nBDFLines += 1
        nExtraFields = 8 * nBDFLines -  nFieldsMain
        fields += [None] * nExtraFields

    try:
        out = '%-8s' % (fields[0]+'*')
    except:
        print("ERROR!  fields=%s" % fields)
        sys.stdout.flush()
        raise

    for i in xrange(1, len(fields)):
        field = fields[i]
        try:
            out += print_field_16(field)
        except:
            print("bad fields = %s" % fields)
            raise
        if i % 4 == 0:  # allow 1+4 fields per line
            out = out.rstrip(' ')
            if out[-1] == '\n':  # empty line
                out += '*'
            out += '\n*       '
    out = out.rstrip(' *')  # removes one continuation star
    if not out.endswith('\n'):
        out += '\n'
    return out

if __name__ == '__main__':
    field = print_float_16(-55.1040257079)
    field = print_float_16(-55.1040257078872)
    field = print_float_16(-3.76948125497534)
