def components(card, n, fieldname):
    svalue = card.field(n)
    try:
        value = int(svalue)
    except:
        Type = getType(svalue)
        raise SyntaxError('%s = %r (field #%s) on card must be an integer (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )

    svalue2 = str(value)
    svalue3 = ''.join(sorted(svalue2))
    for i,v in enumerate(svalue3):
        assert v not in svalue3[i+1:]
    return svalue3

def components_or_blank(card, n, fieldname, default=None):
    svalue = card.field(n)
    if svalue is None:
        return default
    elif isinstance(svalue, int):
        svalue = str(svalue)
    else:
        svalue = svalue.strip()

    if svalue:
        return components(card, n, fieldname)
    else:
        return default

def blank(card, n, fieldname, default=None):
    svalue = card.field(n)
    if svalue is None:
        return default
        
    if isinstance(svalue, str):
        svalue = svalue.strip()
        if len(svalue) == 0:
           return default
    Type = getType(svalue)
    raise SyntaxError('%s = %r (field #%s) on card must be blank (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )
    
def field(card, n, fieldname):
    return integer_double_string_or_blank(card, n, fieldname, default=None)

#def integer_double_string_or_blank(card, n, fieldname, default=None):
    #svalue = card.field(n)
    #if svalue is None:
        #return default
    #elif isinstance(svalue, int) or isinstance(svalue, float):
        #return svalue
    
    #svalue = svalue.strip()
    #if '.' in svalue:
        #return double(card, field, n)
    #elif svalue.isdigit():
        #return int(svalue)
    #elif len(svalue) == 0:
        #return default
    #else:
        #return svalue

def integer_double_string_or_blank(card, n, fieldname, default=None):
    #try:
    svalue = card.field(n)
    #except IndexError:
    #raise SyntaxError('%s (field #%s) on card must be an integer, float, string, or blank.\ncard=%s' % (fieldname, n, card) )

    if isinstance(svalue, int):
        return svalue
    elif svalue is None:
        return default
    elif isinstance(svalue, float):
        return svalue

    svalue = svalue.strip()
    if svalue:  # integer/float/string
        if '.' in svalue:  # float
            value = double(card, n, fieldname)
        elif svalue.isdigit():  # int
            try:
                value = int(svalue)
            except ValueError:
                Type = getType(svalue)
                raise SyntaxError('%s = %r (field #%s) on card must be an integer, float, or string (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )
        else:
            value = svalue
        return value
    return default

#def assert_int_bounded_range(card, n, fieldname, lower=None, upper=None):

def fields(f, card, fieldname, i, j=None):  ## @todo: improve fieldname
    fs = []
    if j is None:
        j = len(card)
    for ii in range(i,j):
        fs.append( f(card, ii, fieldname) )
    return fs

def fields_or_blank(f, card, fieldname, i, j=None, defaults=None):
    fs = []
    if j is None:
        j = len(card)
    
    assert j - i == len(defaults), 'j=%s i=%s j-i=%s len(defaults)=%s\ncard=%s' % (j,i,j-i, len(defaults), card)
    for ii, default in enumerate(defaults):
        fs.append( f(card, ii + i, fieldname + str(ii), default) )
    return fs

def integer(card, n, fieldname):
    try:
        svalue = card.field(n)
    except IndexError:
        raise SyntaxError('%s (field #%s) on card must be an integer.' % (fieldname, n) )

    if isinstance(svalue, float):
        Type = getType(svalue)
        raise SyntaxError('%s = %r (field #%s) on card must be an integer (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )
        
    try:
        return int(svalue)
    except:
        Type = getType(svalue)
        raise SyntaxError('%s = %r (field #%s) on card must be an integer (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )

def integer_or_blank(card, n, fieldname, default=None):
    #try:
    svalue = card.field(n)
    #except IndexError:
    #    return default

    if isinstance(svalue, int):
        return svalue
    elif svalue is None:
        return default
    elif isinstance(svalue, str) or isinstance(svalue, unicode):
        if len(svalue) == 0:
            return default
        try:
            return int(svalue)
        except:
            Type = getType(svalue)
            raise SyntaxError('%s = %r (field #%s) on card must be an integer or blank (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )

    Type = getType(svalue)
    raise SyntaxError('%s = %r (field #%s) on card must be an integer (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )
    #return default
    
def double(card, n, fieldname):
    try:
        svalue = card.field(n)
    except IndexError:
        raise SyntaxError('%s (field #%s) on card must be a float.\ncard=%s' % (fieldname, n, card) )
    
    if isinstance(svalue, float):
        return svalue
    elif isinstance(svalue, int):
        Type = getType(svalue)
        raise SyntaxError('%s = %r (field #%s) on card must be a float (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )
    elif svalue is None or len(svalue) == 0:  ## None
        Type = getType(svalue)
        raise SyntaxError('%s = %r (field #%s) on card must be a float (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )
    
    if svalue.isdigit(): # if only int
        Type = getType(int(svalue))
        raise SyntaxError('%s = %r (field #%s) on card must be a float (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )
    
    #svalue = svalue.strip()
    try:  # 1.0, 1.0E+3, 1.0E-3
        value = float(svalue)
    except TypeError:
        Type = getType(svalue)
        raise SyntaxError('%s = %r (field #%s) on card must be a float (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )
    except ValueError:  # 1D+3, 1D-3, 1-3
        try:
            svalue = svalue.upper()
            if 'D' in svalue:  # 1.0D+3, 1.0D-3
                svalue2 = svalue.replace('D','E')
                return float(svalue2)
            
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
            Type = getType(svalue)
            raise SyntaxError('%s = %r (field #%s) on card must be a float (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )
    return value

def double_or_blank(card, n, fieldname, default=None):
    #try:
    svalue = card.field(n)
    #except IndexError:
        #return default

    if isinstance(svalue, float):
        return svalue
    elif isinstance(svalue, str) or isinstance(svalue, unicode):
        try:
            return double(card, n, fieldname)
        except:
            Type = getType(svalue)
            raise SyntaxError('%s = %r (field #%s) on card must be a float or blank (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )
    return default
    
def double_or_string(card, n, fieldname):
    #try:
    svalue = card.field(n)
    #except IndexError:
        #raise SyntaxError('%s (field #%s) on card must be a float or string (not blank).\ncard=%s' % (fieldname, n, card) )

    if isinstance(svalue, float):
        return svalue
    elif svalue is None or isinstance(svalue, int):
        Type = getType(svalue)
        raise SyntaxError('%s = %r (field #%s) on card must be an float or string (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )
    elif isinstance(svalue, str) or isinstance(svalue, unicode):
        svalue = svalue.strip()

    if '.' in svalue: # float
        try:
            return double(card, n, fieldname)
        except:
            Type = getType(svalue)
            raise SyntaxError('%s = %r (field #%s) on card must be an float or string (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )
    elif svalue.isdigit(): # fail
        pass
    elif svalue: # string
        return svalue
    Type = getType(svalue)
    raise SyntaxError('%s = %r (field #%s) on card must be an float or string (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )


def double_string_or_blank(card, n, fieldname, default=None):
    """
    @param card BDF card as a list
    @param n field number
    @param fieldname name of field
    @param default the default value for the field (default=None)
    @return a double, string, or default value
    @raise SyntaxError if there is an invalid type
    """
    #try:
    svalue = card.field(n)
    #except IndexError:
        #return default
    if isinstance(svalue, float):
        return svalue
    elif svalue is None:
        return default
    elif isinstance(svalue, str) or isinstance(svalue, unicode):
        svalue = svalue.strip()
    elif isinstance(svalue, int):
        Type = getType(svalue)
        raise SyntaxError('%s = %r (field #%s) on card must be an float, string, or blank (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )

    if '.' in svalue:
        try:
            return double(card, n, fieldname)
        except:
            Type = getType(svalue)
            msg = '%s = %r (field #%s) on card must be a float, string or blank (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card)
            raise SyntaxError(msg)
    elif svalue.isdigit():
        Type = getType(svalue)
        raise SyntaxError('%s = %r (field #%s) on card must be a float, string or blank (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )
    elif svalue == '':
        return default
    return svalue

def integer_or_double(card, n, fieldname):
    """
    @param card a BDF_card object
    @param n the field index
    @param fieldname the name of the field
    @return the value with the proper type
    @raise SyntaxError if there's an invalid type
    """
    #try:
    svalue = card.field(n)
    #except IndexError:
        #raise SyntaxError('%s (field #%s) on card must be an integer or float.\ncard=%s' % (fieldname, n, card) )
    if isinstance(svalue, int) or isinstance(svalue, float):
        return svalue
    elif svalue is None:
        raise SyntaxError('%s (field #%s) on card must be an integer or float (not blank).\ncard=%s' % (fieldname, n, card) )

    if '.' in svalue:  # float/exponent
        try:
            value = double(card, n, fieldname)
        except ValueError:
            Type = getType(svalue)
            raise SyntaxError('%s = %r (field #%s) on card must be a integer or a float (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )
    else:  # int
        try:
            value = int(svalue)
        except:
            Type = getType(svalue)
            raise SyntaxError('%s = %r (field #%s) on card must be an integer or a float (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )
    return value

def integer_double_or_blank(card, n, fieldname, default=None):
    #try:
    svalue = card.field(n)
    #except IndexError:
    #    return default

    if isinstance(svalue, int) or isinstance(svalue, float):
        return svalue
    elif svalue is None:
        return default

    if svalue:  # integer/float
        try:
            return integer_or_double(card, n, fieldname)
        except:
            Type = getType(svalue)
            raise SyntaxError('%s = %r (field #%s) on card must be an integer, float, or blank (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )
    return default
    
def integer_or_string(card, n, fieldname):
    #try:
    svalue = card.field(n)
    #except IndexError:
        #raise SyntaxError('%s (field #%s) on card must be an integer or string.\ncard=%s' % (fieldname, n, card) )
    if isinstance(svalue, int):
        return svalue
    elif svalue is None:
        Type = getType(svalue)
        raise SyntaxError('%s = %r (field #%s) on card must be an integer or string (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )

    if svalue.isdigit():  # int
        try:
            value = int(svalue)
        except ValueError:
            Type = getType(svalue)
            raise SyntaxError('%s = %r (field #%s) on card must be an integer or string (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )
    elif isinstance(svalue, float):
        Type = getType(svalue)
        raise SyntaxError('%s = %r (field #%s) on card must be an integer or string (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )
    else:  # string
        return str(svalue)
    return value

def integer_string_or_blank(card, n, fieldname, default=None):
    #try:
    svalue = card.field(n)
    #except IndexError:
        #return default
    if isinstance(svalue, int):
        return svalue
    elif svalue is None:
        return default
    elif isinstance(svalue, float):
        Type = getType(svalue)
        raise SyntaxError('%s = %r (field #%s) on card must be an integer or string (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )

    if svalue.strip():  # integer/string
        try:
            return integer_or_string(card, n, fieldname)
        except:
            Type = getType(svalue)
            raise SyntaxError('%s = %r (field #%s) on card must be an integer, string, or blank (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )
    return default

def getType(value):
    if value is None:
        Type = 'blank'
    elif isinstance(value, int):
        Type = 'an integer'
    elif isinstance(value, float):
        Type = 'a double'
    elif isinstance(value, str) or isinstance(value, unicode):
        Type = 'a string'
    else:
        Type = str(type(value))
    return Type


def integer_double_or_string(card, n, fieldname):
    #try:
    svalue = card.field(n)
    #except IndexError:
        #raise SyntaxError('%s (field #%s) on card must be an integer, float, or string.\ncard=%s' % (fieldname, n, card) )
    if isinstance(svalue, int) or isinstance(svalue, float):
        return svalue
    
    svalue = str(svalue.strip())
    if svalue:  # integer/float/string
        if '.' in svalue:  # float
            value = double(card, n, fieldname)
        elif svalue.isdigit():  # int
            try:
                value = int(svalue)
            except ValueError:
                raise SyntaxError('%s = %r (field #%s) on card must be an integer, float, or string (not blank).\ncard=%s' % (fieldname, svalue, n, card) )
        else:
            value = svalue
        return value
    Type = getType(svalue)
    raise SyntaxError('%s = %r (field #%s) on card must be an integer, float, or string (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )


def string(card, n, fieldname):
    svalue = card.field(n)
    if isinstance(svalue, str) or isinstance(svalue, unicode):
        svalue = svalue.strip()
    else:
        Type = getType(svalue)
        raise SyntaxError('%s = %r (field #%s) on card must be an string (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )
    
    if svalue.isdigit() or '.' in svalue:
        value = integer_or_double(card, n, fieldname)
        Type = getType(value)
        raise SyntaxError('%s = %r (field #%s) on card must be an string with a character (not %s).\ncard=%s' % (fieldname, value, n, Type, card) )

    if svalue:  # string
        return str(svalue)
    Type = getType(svalue)
    msg = '%s = %r (field #%s) on card must be an string (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card)
    raise SyntaxError(msg)


def string_or_blank(card, n, fieldname, default=None):
    svalue = card.field(n)
    if svalue is None:
        return default
    elif isinstance(svalue, str) or isinstance(svalue, unicode):
        svalue = svalue.strip()
    else:
        Type = getType(svalue)
        raise SyntaxError('%s = %r (field #%s) on card must be an string (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )

    svalue = svalue.strip()
    if svalue.isdigit() or '.' in svalue:
        Type = getType(svalue)
        raise SyntaxError('%s = %r (field #%s) on card must be an string or blank (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )
        
    if svalue:  # string
        return str(svalue)
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

