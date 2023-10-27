import math
from math import log10
from pyNastran.bdf.field_writer_8 import print_scientific_8, print_float_8
from pyNastran.bdf.bdf_interface.assign_type import double_from_str

plus_log_dict = {
    #-1. : "%%8.7f",
    #-1: '',
    0 : "%8.6f",  # 1.23
    1 : "%8.5f",
    2 : "%8.4f",
    3 : "%8.3f",
    4 : "%8.2f",
    5 : "%8.1f",
    6 : "%7.0f.",
}

def get_float_format(value: float):
    #value = 1000000.
    #value = 100000.0
    #value = 123456.78
    #value = 1.2e-10
    if value == 0.0:
        #print('block 0')
        return ' 0.0    '
    if value < 0.:
        abs_value = abs(value)
        # sign, decimal, 1
        if abs_value > 999_995:
            #print('A')
            # -123_456.
            field = print_scientific_8(value)
            assert len(field) == 8, f'{field!r}; n={len(field)}'
            return field
        elif abs_value < 5e-7:
            #print('B')
            field = print_scientific_8(value)
            assert len(field) == 8, f'{field!r}; n={len(field)}'
            return field
        elif value > -0.01:  # -0.001
            field = print_scientific_8(value)
            field2 = "%8.6f" % value  # small value
            field2 = field2.strip('0 ')

            # get rid of the first minus sign, add it on afterwards
            field1 = '-' + field.strip(' 0-').replace('-', 'e-')

            if len(field2) <= 8 and float(field1) == float(field2):
                field = field2.rstrip(' 0')
                field = field.replace('-0.', '-.')
                field = '%8s' % field
            assert len(field) == 8, f'{field!r}; n={len(field)}'
            return field
        elif abs_value < 1.:
            #print('C')
            field = f'{abs_value:7.6f}'
            field2 = '-' + field[1:]
            assert len(field2) == 8, f'{field2!r}; n={len(field2)}'
            return field2
        elif abs_value < 10000.:
            #print('C')
            magnitude = int(math.floor(math.log10(abs_value)))
            format_string = "-%7." + str(5 - magnitude) + "f"
            field = format_string % abs_value
            assert len(field) == 8, f'{field!r}; n={len(field)}'
            #print(f'value = {value!r}')
            return field

        magnitude = int(math.floor(math.log10(abs_value)))
        format_string = "-%6." + str(5 - magnitude) + "f."
        field = format_string % abs_value
        #print('D', value, magnitude, field, len(field))
        #if value < -99_995:
        #if value < 0.1:
            # -0.1 to -1.0
        assert len(field) == 8, f'{field!r}; n={len(field)}'
        return field
    if value < 5e-8:
        field = print_scientific_8(value)
        #print('scientific 1;', field)
        assert len(field) == 8, f'{field!r}; n={len(field)}'
        return field
    elif value < 0.001:
        field = print_scientific_8(value)
        field2 = "%8.7f" % value  # small value
        field2 = field2.strip('0 ')
        field1 = field.replace('-', 'e-')
        if field2 == '.':
            return print_scientific_8(value)
        if len(field2) <= 8 and float(field1) == float(field2):
            field = field2
            field = '%8s' % field.strip(' 0')
            assert len(field) == 8, f'{field!r}; n={len(field)}'
        #print('near 0.001;', field2, field1, field, len(field))
        assert len(field) == 8, f'{field!r}; n={len(field)}'
        return field
    elif value < 1.:
        field = "%8.7f" % value  # same as before...
        #print('near 1;', field, len(field))
    #elif value < 1.0:
        #field = print_scientific_8(value)
        #field2 = "%8.7f" % value  # small value
        #field2 = field2.strip('0 ')

        #field1 = field.replace('-', 'e-')
        #if field2 == '.':
            #return print_scientific_8(value)
        #if len(field2) <= 8 and float(field1) == float(field2):
            #field = field2
            #field = field.strip(' 0')
        field2 = '%8s' % field[1:9].rstrip('0')
        #print('field2 =', field2)
        assert len(field2) == 8, f'{field2!r}; n={len(field2)}'
        return field2
    elif value > 1000000.:
        field = "%8.1f" % value
        if field.index('.') < 8:
            field = '%8.1f' % round(value)
            field = field[0:8]
            #assert '.' != field[0], field
            #print('big;', field)
        else:
            field = print_scientific_8(value)
            #print('scientifc;', field)
        assert len(field) == 8, f'{field!r}; n={len(field)}'
        return field

    i = int(log10(value))
    #if i >= 6.:

    #if value < 5e-8:
        #field = print_scientific_8(value)
        #return field
    # elif value < 1.:
    #     field = "%8.7f" % value  # same as before...
    # elif value < 10.:
    #     field = "%8.6f" % value
    # elif value < 100.:
    #     field = "%8.5f" % value
    # elif value < 1000.:
    #     field = "%8.4f" % value
    # elif value < 10000.:
    #     field = "%8.3f" % value
    # elif value < 100000.:
    #     field = "%8.2f" % value
    # elif value < 1000000.:
    #     field = "%8.1f" % value

    fmt = plus_log_dict.get(i)
    field1 = fmt % value
    field = '%8s' % field1.rstrip('0')
    #print(f'log print; i={i} field={field}')
    assert len(field) == 8, field
    return field

def get_field_zach(x: float):  # pragma: no cover
    if x == 0:
        field_zach = ' 0.0    '
    if x > 0:
        if x < 1e-4: # small value
            [coefficient, power] = ('%8.11e' % x).strip().split('e')
            power = "-" + str(power[1:].lstrip('0'))
            coefficient = coefficient.rstrip('0')
            s = coefficient + power
            if len(s) > 8: # correction to scale
                s = coefficient[:-(len(s)-8)] + power
            # pad for aesthetics
            l = len(s)
            if l <= 7:
                s = ' ' + s + (7 - l)*' '
        elif x > 1e6: # large value
            s = '%8.1f' % x
            if len(s) < 8:
                s = '%8.1f' % round(x)
                s = s[0:8]
            else:
                [coefficient, power] = ('%8.11e' % x).strip().split('e')
                power = "+" + str(power[1:].lstrip('0'))
                coefficient = coefficient.rstrip('0')
                s = coefficient + power
                if len(s) > 8: # correction to scale
                    s = coefficient[:-(len(s)-8)] + power
                # pad for aesthetics
                l = len(s)
                if l <= 7:
                    s = ' ' + s + (7 - l)*' '
        else:
            if x < 1:
                s = ('%8.7f' % x).strip('0')
                if len(s) > 8: # correction to scale
                    s = coefficient[:-(len(s)-8)] + power
                # pad for aesthetics
                l = len(s)
                if l <= 7:
                    s = ' ' + s + (7 - l)*' '
            else:
                magnitude = int(math.floor(math.log10(abs(x))))
                format_string = "%8." + str(6 - magnitude) + "f"
                s = format_string % x
        field_zach = s
    else:
        field_zach = 0.
        raise RuntimeError(x)
    return field_zach

def compare(x: float, stop_on_error: bool=False):
    """ Get optimal way to print the number in 8 characters"""
    #print(x)
    field_new = get_float_format(x).strip()
    if x > 0.1:
        field_new = field_new.rstrip('0')
    field_old = print_float_8(x).strip()
    #assert field_new == field_old, f'field={field_new!r} old_field={field_old!r}'


    field_old = '%8s' % field_old
    #field_zach = get_field_zach(x)
    field_new = get_float_format(x)
    assert len(field_old) == 8
    assert len(field_new) == 8
    #assert len(field_zach) == 8
    field_old = field_old.strip()
    #field_zach = field_zach.strip()
    field_new = field_new.strip()
    value_old = double_from_str(field_old)
    value_new = double_from_str(field_new)
    if '.' not in field_new:
        # it must be a float; major error
        msg = f'x={x!r} old={field_old!r} new={field_new!r}'
        raise RuntimeError(msg)

    if value_old != value_new:
        # loss of precision
        #msg = f'x={x!r} old={field_old!r} field_zach={field_zach!r} new={field_new!r}'
        msg = f'x={x!r} old={field_old!r} new={field_new!r}'
        if stop_on_error:
            raise RuntimeError(msg)

def main(get_optimal_short_form_float):  # pragma: no cover
    #get_optimal_short_form_float = get_field_zach
    #get_optimal_short_form_float = print_float_8
    #get_optimal_short_form_float = get_float_format

    value=-0.8732806594999999 #; field_old='-.873281' field_new='-8.733-1'
    s = get_optimal_short_form_float(value)
    value=-0.6862676525499999 #; field_old='-.686268' field_new='-6.863-1'
    s = get_optimal_short_form_float(value)
    value=-0.647050147 #; field_old=' -.64705' field_new='-6.471-1'
    s = get_optimal_short_form_float(value)
    value=-0.9526646215 #; field_old='-.952665' field_new='-9.527-1'
    s = get_optimal_short_form_float(value)
    value=-0.6291313619 #; field_old='-.629131' field_new='-6.291-1'
    s = get_optimal_short_form_float(value)

    x = 2e-11
    s = get_optimal_short_form_float(x)
    x = 12e-11
    s = get_optimal_short_form_float(x)
    x = 12.1e-11
    s = get_optimal_short_form_float(x)
    x = 12.12e-11
    s = get_optimal_short_form_float(x)
    x = 12.123e-11
    s = get_optimal_short_form_float(x)
    x = 12.1234e-11
    s = get_optimal_short_form_float(x)
    x = 12e-12
    s = get_optimal_short_form_float(x)
    x = 12.1e-12
    s = get_optimal_short_form_float(x)
    x = 12.12e-12
    s = get_optimal_short_form_float(x)
    x = 12.123e-12
    s = get_optimal_short_form_float(x)
    x = 12.1234e-12
    s = get_optimal_short_form_float(x)
    x = 1.0e-5
    s = get_optimal_short_form_float(x)
    x = 1.0e-4
    s = get_optimal_short_form_float(x)
    x = 1.0e-3
    s = get_optimal_short_form_float(x)
    x = 1.0e-2
    s = get_optimal_short_form_float(x)
    x = 1.1e-5
    s = get_optimal_short_form_float(x)
    x = 1.1e-4
    s = get_optimal_short_form_float(x)
    x = 1.1e-3
    s = get_optimal_short_form_float(x)
    x = 1.1e-2
    s = get_optimal_short_form_float(x)
    x = 1.12e-5
    s = get_optimal_short_form_float(x)
    x = 1.12e-4
    s = get_optimal_short_form_float(x)
    x = 1.12e-3
    s = get_optimal_short_form_float(x)
    x = 1.12e-2
    s = get_optimal_short_form_float(x)
    x = 1.123e-5
    s = get_optimal_short_form_float(x)
    x = 1.123e-4
    s = get_optimal_short_form_float(x)
    x = 1.123e-3
    s = get_optimal_short_form_float(x)
    x = 1.123e-2
    s = get_optimal_short_form_float(x)
    x = 1.e5
    s = get_optimal_short_form_float(x)
    x = 1.e6
    s = get_optimal_short_form_float(x)
    x = 1.e7
    s = get_optimal_short_form_float(x)
    x = 1.e8
    s = get_optimal_short_form_float(x)
    x = 1.2e8
    s = get_optimal_short_form_float(x)
    x = 1.23e8
    s = get_optimal_short_form_float(x)
    x = 1.234e8
    s = get_optimal_short_form_float(x)
    x = 1.2345e8
    s = get_optimal_short_form_float(x)
    x = 1.23456e8
    s = get_optimal_short_form_float(x)
    #------------------------------------------------
    x = 12345678.
    s = get_optimal_short_form_float(x)
    x = 1234567.8
    s = get_optimal_short_form_float(x)
    x = 123456.78
    s = get_optimal_short_form_float(x)
    x = 12345.678
    s = get_optimal_short_form_float(x)
    x = 1234.5678
    s = get_optimal_short_form_float(x)
    x = 123.45678
    s = get_optimal_short_form_float(x)
    x = 12.345678
    s = get_optimal_short_form_float(x)
    x = 1.2345678
    s = get_optimal_short_form_float(x)
    x = .12345678
    s = get_optimal_short_form_float(x)
    x = .012345678
    s = get_optimal_short_form_float(x)
    x = .0012345678
    s = get_optimal_short_form_float(x)
    x = .00012345678
    s = get_optimal_short_form_float(x)
    x = .000012345678
    s = get_optimal_short_form_float(x)
    x = .0000012345678
    s = get_optimal_short_form_float(x)
    # }}}

def main_og():  # pragma: no cover
    main(print_float_8)
def main_zach():  # pragma: no cover
    main(get_field_zach)
def main_new():  # pragma: no cover
    main(get_float_format)

def run():  # pragma: no cover
    value = -0.8732806594999999
    out = get_float_format(value)
    print(value, out, len(out))
    return
    import timeit

    #main(get_field_zach)
    #2000000 loops, best of 5: 149 nsec per loop
    #2000000 loops, best of 5: 151 nsec per loop
    #2000000 loops, best of 5: 150 nsec per loop
    #2000000 loops, best of 5: 149 nsec per loop
    #2000000 loops, best of 5: 149 nsec per loop
    #2000000 loops, best of 5: 159 nsec per loop
    #2000000 loops, best of 5: 152 nsec per loop
    #2000000 loops, best of 5: 158 nsec per loop
    #2000000 loops, best of 5: 152 nsec per loop
    #2000000 loops, best of 5: 152 nsec per loop
    num = 10

    print('OG:')
    print(timeit.timeit('main_og()', number=num, setup="from __main__ import main_og"))

    print('zach:')
    print(timeit.timeit('main_zach()', number=num, setup="from __main__ import main_zach"))

    #get_optimal_short_form_float = get_field_zach
    #get_optimal_short_form_float = print_float_8
    #get_optimal_short_form_float = get_float_format

    print('new:')
    print(timeit.timeit('main_new()', number=num, setup="from __main__ import main_new"))
    x = 1
    #main()


if __name__ == '__main__':  # pragma: no cover
    run()
