## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
from pyNastran.utils import is_string
from .BDF_Card import BDFCard

def components(card, n, fieldname):
    """
    :param card:      BDF card as a list
    :param n:         field number
    :param fieldname: name of field
    """
    assert isinstance(card, BDFCard), type(card)
    assert isinstance(n, int), type(n)
    assert is_string(fieldname), type(fieldname)
    svalue = card.field(n)
    if isinstance(svalue, int):
        pass
    elif svalue is None or '.' in svalue:
        Type = getType(svalue)
        msg = '%s = %r (field #%s) on card must be an integer (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) 
        raise SyntaxError(msg)

    try:
        value = int(svalue)
    except:
        Type = getType(svalue)
        msg = '%s = %r (field #%s) on card must be an integer (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) 
        raise SyntaxError(msg)
    if value > 0  and is_string(svalue):
        if '0' in svalue:
            value2 = str(svalue).replace('0', '')
            msg = '%s = %r (field #%s) on card must contain 0 or %s (not both).\ncard=%s' % (fieldname, svalue, n, value2, card) 
            raise SyntaxError(msg)
    svalue2 = str(value)
    svalue3 = ''.join(sorted(svalue2))
    for i,v in enumerate(svalue3):
        if v not in '0123456':
            msg = '%s = %r (field #%s) on card contains an invalid component %r.\ncard=%s' % (fieldname, svalue, n, v, card) 
            raise SyntaxError(msg)
        if v in svalue3[i + 1:]:
            msg = '%s = %r (field #%s) on card must not contain duplicate entries.\ncard=%s' % (fieldname, svalue, n, card) 
            raise SyntaxError(msg)
    return svalue3

def components_or_blank(card, n, fieldname, default=None):
    """
    :param card:      BDF card as a list
    :param n:         field number
    :param fieldname: name of field
    :param default:   the default value for the field (default=None)
    """
    assert isinstance(card, BDFCard), type(card)
    assert isinstance(n, int), type(n)
    assert is_string(fieldname), type(fieldname)
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
    """
    :param card:      BDF card as a list
    :param n:         field number
    :param fieldname: name of field
    :param default:   the default value for the field (default=None)
    """
    assert isinstance(card, BDFCard), type(card)
    assert isinstance(n, int), type(n)
    assert is_string(fieldname), type(fieldname)
    svalue = card.field(n)
    if svalue is None:
        return default
        
    if is_string(svalue):
        svalue = svalue.strip()
        if len(svalue) == 0:
           return default
    Type = getType(svalue)
    raise SyntaxError('%s = %r (field #%s) on card must be blank (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )
    
def field(card, n, fieldname):
    """
    :param card:      BDF card as a list
    :param n:         field number
    :param fieldname: name of field
    """
    assert isinstance(card, BDFCard), type(card)
    assert isinstance(n, int), type(n)
    assert is_string(fieldname), type(fieldname)
    return integer_double_string_or_blank(card, n, fieldname, default=None)

def integer_double_string_or_blank(card, n, fieldname, default=None):
    """
    :param card:      BDF card as a list
    :param n:         field number
    :param fieldname: name of field
    :param default:   the default value for the field (default=None)
    """
    assert isinstance(card, BDFCard), type(card)
    assert isinstance(n, int), type(n)
    assert is_string(fieldname), type(fieldname)
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

def fields(f, card, fieldname, i, j=None):
    """
    .. todo:: improve fieldname
    """
    assert isinstance(card, BDFCard), type(card)
    assert is_string(fieldname), type(fieldname)
    fs = []
    if j is None:
        j = len(card)
    for ii in range(i,j):
        fs.append( f(card, ii, fieldname) )
    return fs

def fields_or_blank(f, card, fieldname, i, j=None, defaults=None):
    assert isinstance(card, BDFCard), type(card)
    assert is_string(fieldname), type(fieldname)
    fs = []
    if j is None:
        j = len(card)
    
    assert j - i == len(defaults), 'j=%s i=%s j-i=%s len(defaults)=%s\ncard=%s' % (j,i,j-i, len(defaults), card)
    for ii, default in enumerate(defaults):
        fs.append( f(card, ii + i, fieldname + str(ii), default) )
    return fs

def integer(card, n, fieldname):
    """
    :param card:      BDF card as a list
    :param n:         field number
    :param fieldname: name of field
    """
    assert isinstance(card, BDFCard), type(card)
    assert isinstance(n, int), type(n)
    assert is_string(fieldname), type(fieldname)
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
    """
    :param card:      BDF card as a list
    :param n:         field number
    :param fieldname: name of field
    :param default:   the default value for the field (default=None)
    """
    assert isinstance(card, BDFCard), type(card)
    assert isinstance(n, int), type(n)
    assert is_string(fieldname), type(fieldname)
    #try:
    svalue = card.field(n)
    #except IndexError:
    #    return default

    if isinstance(svalue, int):
        return svalue
    elif svalue is None:
        return default
    elif is_string(svalue):
        if len(svalue) == 0:
            return default
        elif '.' in svalue:
            Type = getType(svalue)
            raise SyntaxError('%s = %r (field #%s) on card must be an integer or blank (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )

        try:
            return int(svalue)
        except:
            Type = getType(svalue)
            raise SyntaxError('%s = %r (field #%s) on card must be an integer or blank (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )

    Type = getType(svalue)
    raise SyntaxError('%s = %r (field #%s) on card must be an integer (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )
    #return default
    
def double(card, n, fieldname):
    """
    :param card:      BDF card as a list
    :param n:         field number
    :param fieldname: name of field
    """
    assert isinstance(card, BDFCard), type(card)
    assert isinstance(n, int), type(n)
    assert is_string(fieldname), type(fieldname)
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
    """
    :param card:      BDF card as a list
    :param n:         field number
    :param fieldname: name of field
    :param default:   the default value for the field (default=None)
    """
    assert isinstance(card, BDFCard), type(card)
    assert isinstance(n, int), type(n)
    assert is_string(fieldname), type(fieldname)
    #try:
    svalue = card.field(n)
    #except IndexError:
        #return default
    
    if isinstance(svalue, float):
        return svalue
    elif isinstance(svalue, int):
        Type = getType(svalue)
        raise SyntaxError('%s = %r (field #%s) on card must be a float or blank (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )
    elif is_string(svalue):
        svalue = svalue.strip()
        if not svalue:
            return default
        try:
            return double(card, n, fieldname)
        except:
            Type = getType(svalue)
            raise SyntaxError('%s = %r (field #%s) on card must be a float or blank (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )
    return default
    
def double_or_string(card, n, fieldname):
    """
    :param card:      BDF card as a list
    :param n:         field number
    :param fieldname: name of field
    """
    assert isinstance(card, BDFCard), type(card)
    assert isinstance(n, int), type(n)
    assert is_string(fieldname), type(fieldname)
    #try:
    svalue = card.field(n)
    #except IndexError:
        #raise SyntaxError('%s (field #%s) on card must be a float or string (not blank).\ncard=%s' % (fieldname, n, card) )

    if isinstance(svalue, float):
        return svalue
    elif svalue is None or isinstance(svalue, int):
        Type = getType(svalue)
        raise SyntaxError('%s = %r (field #%s) on card must be an float or string (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )
    elif is_string(svalue):
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
    :param card:      BDF card as a list
    :param n:         field number
    :param fieldname: name of field
    :param default:   the default value for the field (default=None)
    :returns value:   a double, string, or default value
    :raises SyntaxError: if there is an invalid type
    """
    assert isinstance(card, BDFCard), type(card)
    assert isinstance(n, int), type(n)
    assert is_string(fieldname), type(fieldname)
    #try:
    svalue = card.field(n)
    #except IndexError:
        #return default
    if isinstance(svalue, float):
        return svalue
    elif svalue is None:
        return default
    elif is_string(svalue):
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
    :param card:      BDF card as a list
    :param n:         field number
    :param fieldname: name of field
    :returns value:   the value with the proper type
    :raises SyntaxError: if there's an invalid type
    """
    assert isinstance(card, BDFCard), type(card)
    assert isinstance(n, int), type(n)
    assert is_string(fieldname), type(fieldname)
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
    """
    :param card:      BDF card as a list
    :param n:         field number
    :param fieldname: name of field
    :param default:   the default value for the field (default=None)
    """
    assert isinstance(card, BDFCard), type(card)
    assert isinstance(n, int), type(n)
    assert is_string(fieldname), type(fieldname)
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
    """
    :param card:      BDF card as a list
    :param n:         field number
    :param fieldname: name of field
    :param default:   the default value for the field (default=None)
    """
    assert isinstance(card, BDFCard), type(card)
    assert isinstance(n, int), type(n)
    assert is_string(fieldname), type(fieldname)
    #try:
    svalue = card.field(n)
    #except IndexError:
        #raise SyntaxError('%s (field #%s) on card must be an integer or string.\ncard=%s' % (fieldname, n, card) )
    if isinstance(svalue, int):
        return svalue
    elif isinstance(svalue, float):
        Type = getType(svalue)
        raise SyntaxError('%s = %r (field #%s) on card must be an integer or string (not %s).\ncard=%s' % (fieldname, svalue, n, Type, card) )
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
    """
    :param card:      BDF card as a list
    :param n:         field number
    :param fieldname: name of field
    :param default:   the default value for the field (default=None)
    """
    assert isinstance(card, BDFCard), type(card)
    assert isinstance(n, int), type(n)
    assert is_string(fieldname), type(fieldname)
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
    """
    Get the type of the input value in a form that is clear.
    :param value: the value to get the type of
    """
    #print('Type value=%s' % value)
    try:
        value = interpret_value(value)
    except:
        pass
    if value is None:
        Type = 'blank'
    elif isinstance(value, int):
        Type = 'an integer'
    elif isinstance(value, float):
        Type = 'a double'
    elif is_string(value):
        Type = 'a string'
    else:
        Type = str(type(value))
    return Type


def integer_double_or_string(card, n, fieldname):
    """
    :param card:      BDF card as a list
    :param n:         field number
    :param fieldname: name of field
    """
    assert isinstance(card, BDFCard), type(card)
    assert isinstance(n, int), type(n)
    assert is_string(fieldname), type(fieldname)
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
    """
    :param card:      BDF card as a list
    :param n:         field number
    :param fieldname: name of field
    """
    assert isinstance(card, BDFCard), type(card)
    assert isinstance(n, int), type(n)
    assert is_string(fieldname), type(fieldname)
    svalue = card.field(n)
    if is_string(svalue):
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
    """
    :param card:      BDF card as a list
    :param n:         field number
    :param fieldname: name of field
    :param default:   the default value for the field (default=None)
    """
    assert isinstance(card, BDFCard), type(card)
    assert isinstance(n, int), type(n)
    assert is_string(fieldname), type(fieldname)
    svalue = card.field(n)
    if svalue is None:
        return default
    elif is_string(svalue):
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


def interpretValue(valueRaw, card='', debug=False):
    """
    .. seealso:: interpret_value
    .. deprecated: will be replaced in version 0.7 with interpret_value
    """
    warnings.warn('interpretValue has been deprecated; use '
                  'interpret_value', DeprecationWarning, stacklevel=2)
    return interpret_value(valueRaw, card, debug)


def interpret_value(valueRaw, card='', debug=False):
    """Converts a value from nastran format into python format."""
    #debug = True
    if debug:
        print("v1 = |%s|" % valueRaw)
    #lvalue = valueRaw.lstrip()
    try:
        valueIn = valueRaw.lstrip().rstrip(' *').upper()
    except AttributeError:  # it's already an int/float
        msg = 'valueRaw=%s type=%s' % (valueRaw, type(valueRaw))
        assert isinstance(valueRaw, int) or isinstance(valueRaw, float), msg
        return valueRaw

    if debug:
        pass
        #print("v2 = |%s|" % valueIn)
    if len(valueIn) == 0:
        if debug:
            print("BLANK!")
        return None

    if valueIn[0].isalpha():
        if debug:
            print("STRING!")
            print("valueStr = %s" % valueIn)
        return valueIn

    if '=' in valueIn or '(' in valueIn or '*' in valueRaw:
        if debug:
            print("=(! - special formatting")
        return valueRaw.strip()
    #valueIn = valueIn.upper()
    # int, float, string, exponent
    valuePositive = valueIn.strip('+-')
    if debug:
        print("isDigit = %s" % valuePositive.isdigit())
    if valuePositive.isdigit():
        if debug:
            print("INT!")
        return int(valueIn)
    try:
        value = float(valueIn)
        if debug:
            print("FLOAT!")
        return value
    except ValueError:
        pass

    #if('=' in valueIn or '(' in valueIn or ')' in valueIn):
    #    print("=()!")
    #    return valueIn

    # if there are non-floats/scientific notation -> string
    noED = list(set(valueIn) - set('ED 1234567890+-'))
    word = ''.join(noED)
    #print "word=|%s|" %word
    if word.isalpha():
        if debug:
            print("WORD!")
        return valueIn

    v0 = valueIn[0]
    if '-' == v0 or '+' == v0:
        valueLeft = valueIn[1:]  # truncate the sign for now
    else:
        v0 = '+'  # inplied positive value
        valueLeft = valueIn

    #print "valueIn = |%s|" % valueIn
    #print "v0 = |%s|" %v0
    if v0 == '-':
        vFactor = -1.
    elif v0 == '+' or v0.isdigit():
        vFactor = 1.
    else:
        msg = ('the only 2 cases for a float/scientific are +/- for v0...'
               'valueRaw=|%s| v0=|%s| card=%s' % (valueRaw, v0, card))
        raise SyntaxError(msg)

    # dont include the 1st character, find the exponent
    vm = valueIn.find('-', 1)
    vp = valueIn.find('+', 1)
    if vm > 0:
        sline = valueLeft.split('-')
        expFactor = -1.
    elif vp > 0:
        sline = valueLeft.split('+')
        expFactor = 1.
    else:
        msg = ("I thought this was in scientific notation, but i can't find "
               "the exponent sign...valueRaw=|%s| valueLeft=|%s| "
               "card=%s\nYou also might have mixed tabs/spaces/commas."
               % (valueRaw, valueLeft, card))
        raise SyntaxError(msg)

    try:
        s0 = vFactor * float(sline[0])
        s1 = expFactor * int(sline[1])
    except ValueError:
        msg = "vm=%s vp=%s valueRaw=|%s| sline=%s" % (vm, vp, valueRaw, sline)
        raise SyntaxError('cannot parse sline[0] into a float and sline[1] '
                          'into an integer\n%s\nYou HAVE mixed '
                          'tabs/spaces/commas!  Fix it!' % msg)

    value = s0 * 10 ** (s1)
    #print "valueOut = |%s|" %value

    if debug:
        print("SCIENTIFIC!")
    return value


def string_parser(stringIn):
    """not used"""
    typeCheck = ''
    n = 0
    for (i, s) in enumerate(stringIn):
        if s in "+-":
            state = '+'
        elif s == " ":
            state = ' '
        elif s == ".":
            state = '.'
        elif s in "eEdD":
            state = 'e'
        elif s.isdigit():
            state = '1'
        elif s.isalpha() or s in "()*/=]['\"":
            return 'string'  # string character
        else:
            msg = "s=|%r|" % s
            raise SyntaxError(msg)

        #print "s=%s stringIn[i-1]=%s" % (state, typeCheck[i-1])
        #print "i=%s s=%s typeCheck=%s" % (i, s, typeCheck)
        if i == 0:
            typeCheck += state
            n += 1
        elif typeCheck[n - 1] != state:
            typeCheck += state
            n += 1
        elif state in 'e .+':  # double e, space, dot, plus
            return 'string'

    if typeCheck == ' ':
        return None

    typeCheck = typeCheck.strip()
    if typeCheck in ['1', '+1']:  # integer
        return int(stringIn)

    elif typeCheck in ['1.', '1.1', '.1',  # float
                       '+1.', '+1.1', '+.1']:
        return float(stringIn)

    elif typeCheck in ['1.1e1', '1.1e+1', '1.e1', '1.e+1',  # python scientific
                       '+1.1e1', '+1.1e+1', '+1.e1', '+1.e+1',
                       '.1e1', '.1e+1', '+.1e1', '+.1e+1', ]:
        return float(stringIn)

    elif typeCheck in ['1+1', '+1+1', '.1+1', '+.1+1']:  # nastran scientific
        stringReversed = stringIn[::-1]
        i = stringReversed.index('+')
        lString = list(stringIn)
        lString.insert(-i - 1, 'e')
        #print "lString = ", lString
        out = ''.join(lString)
        print("out = %s" % out)
        return float(out)
    else:
        #print "string = ", stringIn
        #print "typeCheck = ", typeCheck
        #return 'string'
        return stringIn

    print("typeCheck = |%s|" % typeCheck)
    raise RuntimeError('error parsing a card...this should never happen...')

if __name__ == '__main__':
    print(string_parser('123'))
    print(string_parser('+123'))
    print(string_parser('.234'))
    print(string_parser('+.234'))
    print(string_parser('-.234'))
    print(string_parser('1+5'))
    print("abc = |%s|" % string_parser('abc'))
    print("eeg = |%s|" % string_parser('eeg'))
    #print("e1 = |%s|" % string_parser('\T'))
    print(string_parser('.e1'))