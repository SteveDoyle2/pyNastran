"""
defines methods for reading panair values:
 - double(value, name)
 - integer(value, name)
 - integer_or_blank(value, name, default=None)
 - double_or_blank(value, name, default=None)

defines methods for writing panair values:

 - fortran_value(value)
"""
def double(value, name):
    """casts to an float value"""
    if isinstance(value, float):
        return value
    value = float(value)
    return value

def integer(value, name):
    """casts to an integer value"""
    if isinstance(value, int):
        return value
    value = value
    value = float(value)
    if not value.is_integer():
        raise RuntimeError('%s=%r is not an integer' % name, value)
    return int(value)

def fortran_value(value):
    return "%8.4E" % value

def integer_or_blank(value, name, default=None):
    value = value.strip()
    if not value:
        return default

    value = float(value)
    if not value.is_integer():
        raise RuntimeError('%s=%r is not an integer' % name, value)
    return int(value)

def double_or_blank(value, name, default=None):
    value = value.strip()
    if not value:
        return default
    try:
        value = float(value)
    except ValueError:
        raise SyntaxError('%s=%r is not a float' % (name, value))
    return value
