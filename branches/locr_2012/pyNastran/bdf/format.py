def field(card, n, fieldname):
    return integer_double_string_or_blank(card, n, fieldname, default=None)

def fields(card, fieldname, i, j):
    pass

def integer(card, n, fieldname):
    try:
        svalue = card[n]
    except IndexError:
        raise RuntimeError('%s (field # %s) on card must be an integer.\ncard=\n%s' % (fieldname, n, card) )

    try:
        return int(svalue)
    except:
        raise RuntimeError('%s = %s (field # %s) on card must be an integer.\ncard=\n%s' % (fieldname, svalue, n, card) )

def integer_or_blank(card, n, fieldname, default=None):
    try:
        svalue = card[n]
    except IndexError:
        return default

    if svalue.strip():
        try:
            return int(svalue)
        except:
            raise RuntimeError('%s = %s (field # %s) on card must be an integer.or blank.\ncard=\n%s' % (fieldname, svalue, n, card) )
    return default
    
def double(card, n, fieldname):
    try:
        svalue = card[n]
    except IndexError:
        raise RuntimeError('%s (field # %s) on card must be an float.\ncard=\n%s' % (fieldname, n, card) )

    try:  # 1.0, 1.0E+3, 1.0E-3
        value = float(svalue)
    except ValueError:
        try:
            svalue = svalue.upper()
            if 'D' in svalue:  # 1.0D+3, 1.0D-3
                svalue = svalue.replace('D','E')
                return float(value)

            # 1.0+3, 1.0-3
            sign = ''
            if '+' in svalue[0] or '-' in svalue[0]:
                svalue = svalue[1:]
                sign = svalue[0]
            if '+' in svalue:
                svalue = sign + svalue.replace('+','E+')
            elif '-' in svalue:
                svalue = sign + svalue.replace('-','E-')

            value = float(svalue)
        except ValueError:
            raise RuntimeError('%s = %s (field # %s) on card must be an float.\ncard=\n%s' % (fieldname, svalue, n, card) )
    return value

def double_or_blank(card, n, fieldname, default=None):
    try:
        svalue = card[n]
    except IndexError:
        return default

    if svalue:
        try:
            return double(card, n, fieldname)
        except:
            raise RuntimeError('%s = %s (field # %s) on card must be an float or blank.\ncard=\n%s' % (fieldname, svalue, n, card) )
    return default
    
def double_or_string(card, n, fieldname):
    try:
        svalue = card[n]
    except IndexError:
        raise RuntimeError('%s (field # %s) on card must be an integer.\ncard=\n%s' % (fieldname, n, card) )

    if '.' in svalue: # float
        try:
            return double(card, n, fieldname)
        except:
            raise RuntimeError('%s = %s (field # %s) on card must be an float or string (not blank).\ncard=\n%s' % (fieldname, svalue, n, card) )
    else: # string
        return svalue.strip()
    raise RuntimeError('%s = %s (field # %s) on card must be an float or string (not blank).\ncard=\n%s' % (fieldname, svalue, n, card) )
    
def double_string_or_blank(card, n, fieldname, default):
    try:
        svalue = card[n].strip()
    except IndexError:
        return default

    if '.' in svalue:
        try:
            return double(card, n, fieldname)
        except:
            raise RuntimeError('%s = %s (field # %s) on card must be an float, string or blank.\ncard=\n%s' % (fieldname, svalue, n, card) )
    else:
        return svalue
    return default
    
def integer_or_double(card, n, fieldname):
    try:
        svalue = card[n]
    except IndexError:
        raise RuntimeError('%s (field # %s) on card must be an integer or float.\ncard=\n%s' % (fieldname, n, card) )

    if '.' in svalue:  # float/exponent
        try:
            double(card, n, fieldname)
        except ValueError:
            raise RuntimeError('%s = %s (field # %s) on card must be a integer or a float.\ncard=\n%s' % (fieldname, svalue, n, card) )
    else:  # int
        try:
            value = int(svalue)
        except:           
            raise RuntimeError('%s = %s (field # %s) on card must be an integer or a float.\ncard=\n%s' % (fieldname, svalue, n, card) )
    return value

def integer_double_or_blank(card, n, fieldname, default=None):
    try:
        svalue = card[n].strip()
    except IndexError:
        return default

    if svalue:  # integer/float
        try:
            return integer_or_double(card, n, fieldname)
        except:
            raise RuntimeError('%s = %s (field # %s) on card must be an integer, float, or blank.\ncard=\n%s' % (fieldname, svalue, n, card) )
    return default
    
def integer_or_string(card, n, fieldname):
    try:
        svalue = card[n].strip()
    except IndexError:
        raise RuntimeError('%s (field # %s) on card must be an integer or string.\ncard=\n%s' % (fieldname, n, card) )

    if svalue.isdigit():  # int
        try:
            value = int(svalue)
        except ValueError:            
            raise RuntimeError('%s = %s (field # %s) on card must be an integer or string.\ncard=\n%s' % (fieldname, svalue, n, card) )
    else:  # string
        return svalue
    return value

def integer_string_or_default(card, n, fieldname, default=None):
    try:
        svalue = card[n]
    except IndexError:
        return default

    if svalue.strip():  # integer/string
        try:
            return integer_or_string(card, n, fieldname)
        except:
            raise RuntimeError('%s = %s (field # %s) on card must be an integer, string, or blank.\ncard=\n%s' % (fieldname, svalue, n, card) )
    return default

def integer_float_string(card, n, fieldname):
    try:
        svalue = card[n].strip()
    except IndexError:
        raise RuntimeError('%s (field # %s) on card must be an integer, float, or string.\ncard=\n%s' % (fieldname, n, card) )

    if svalue:  # integer/float/string
        if '.' in svalue:  # float
            value = double(card, n, fieldname)
        elif svalue.isdigit():  # int
            try:
                value = int(svalue)
            except ValueError:
                raise RuntimeError('%s = %s (field # %s) on card must be an integer, float, or string (not blank).\ncard=\n%s' % (fieldname, svalue, n, card) )
        else:
            value = svalue
        return value
    raise RuntimeError('%s = %s (field # %s) on card must be an integer, float, or string (not blank).\ncard=\n%s' % (fieldname, svalue, n, card) )

def integer_float_string_or_blank(card, n, fieldname, default=None):
    try:
        svalue = card[n].strip()
    except IndexError:
        raise RuntimeError('%s (field # %s) on card must be an integer, float, string, or blank.\ncard=\n%s' % (fieldname, n, card) )

    if svalue:  # integer/float/string
        if '.' in svalue:  # float
            value = double(card, n, fieldname)
        elif svalue.isdigit():  # int
            try:
                value = int(svalue)
            except ValueError:
                raise RuntimeError('%s = %s (field # %s) on card must be an integer, float, or string.\ncard=\n%s' % (fieldname, svalue, n, card) )
        else:
            value = svalue
        return value
    return default

def string(card, n, fieldname):
    svalue = card[n].strip()
    if svalue:  # string
        return svalue
    raise RuntimeError('%s = %s (field # %s) on card must be an string (not blank).\ncard=\n%s' % (fieldname, svalue, n, card) )

def string_or_blank(card, n, fieldname, default=None):
    svalue = card[n].strip()
    if svalue:  # string
        return svalue
    return default

# int                    - done
# int/blank              - done
# int/float              - done
# int/float/blank        - done
# int/float/string       - done
# int/float/string/blank - done
# int/string             - done
# int/string/blank       - done


# float              - done
# float/blank        - done
# float/string       - done
# float/string/blank - done

# string       - done
# string/blank - done

