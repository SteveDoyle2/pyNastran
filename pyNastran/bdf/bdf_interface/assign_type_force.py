import warnings
from typing import Optional
from .assign_type import double, _get_dtype
from .bdf_card import BDFCard
from pyNastran.utils.numpy_utils import (
    integer_types, float_types)
from .assign_type import integer_double_or_blank

def force_integer(card: BDFCard, ifield: int, fieldname: str) -> int:
    """see ``integer``"""
    svalue = card.field(ifield)
    if isinstance(svalue, float_types):
        warnings.warn('%s = %r (field #%s) on card must be an integer (not a double).\n'
                      'card=%s' % (fieldname, svalue, ifield, card))
        return int(svalue)

    try:
        return int(svalue)
    except ValueError:
        # int('s')
        # int('5.0')
        if '.' in svalue:
            try:
                return _force_integer(svalue)
            except ValueError:
                pass
        dtype = _get_dtype(svalue)
        raise SyntaxError('%s = %r (field #%s) on card must be an integer (not %s).\n'
                          'card=%s' % (fieldname, svalue, ifield, dtype, card))

    except TypeError:
        dtype = _get_dtype(svalue)
        raise SyntaxError('%s = %r (field #%s) on card must be an integer (not %s).\n'
                          'card=%s' % (fieldname, svalue, ifield, dtype, card))

def _force_integer(svalue: str) -> int:
    sline = svalue.split('.')
    if len(sline) == 2:
        avalue = int(sline[0])
        if sline[1] == '':
            return avalue

        bvalue = int(sline[1])
        if bvalue != 0:
            raise ValueError()
        return avalue


def force_double(card: BDFCard, ifield: int, fieldname: str,
                 end: str='') -> float:
    """see ``double``"""
    svalue = card.field(ifield)

    if isinstance(svalue, float_types):
        return svalue
    elif isinstance(svalue, integer_types):
        dtype = _get_dtype(svalue)
        warnings.warn('%s = %r (field #%s) on card must be a float (not %s).\n'
                      'card=%s%s' % (fieldname, svalue, ifield, dtype, card, end))
        return float(svalue)
    elif svalue is None or len(svalue) == 0:  ## None
        dtype = _get_dtype(svalue)
        raise SyntaxError('%s = %r (field #%s) on card must be a float (not %s).\n'
                          'card=%s%s' % (fieldname, svalue, ifield, dtype, card, end))

    #if svalue.isdigit():  # 1, not +1, or -1
        ## if only int
        #raise SyntaxError('%s = %r (field #%s) on card must be a float (not an integer).\n'
                          #'card=%s' % (fieldname, svalue, ifield, card))

    try:
        # 1.0, 1.0E+3, 1.0E-3
        value = float(svalue)
    except TypeError:
        dtype = _get_dtype(svalue)
        raise SyntaxError('%s = %r (field #%s) on card must be a float (not %s).\n'
                          'card=%s%s' % (fieldname, svalue, ifield, dtype, card, end))
    except ValueError:
        # 1D+3, 1D-3, 1-3
        try:
            svalue = svalue.upper()
            if 'D' in svalue:
                # 1.0D+3, 1.0D-3
                svalue2 = svalue.replace('D', 'E')
                return float(svalue2)

            # 1.0+3, 1.0-3
            sign = ''
            if svalue[0] in ('+', '-'):
                sign = svalue[0]
                svalue = svalue[1:]
            if '+' in svalue:
                svalue = sign + svalue.replace('+', 'E+')
            elif '-' in svalue:
                svalue = sign + svalue.replace('-', 'E-')

            value = float(svalue)
        except ValueError:
            dtype = _get_dtype(svalue)
            raise SyntaxError('%s = %r (field #%s) on card must be a float (not %s).\n'
                              'card=%s%s' % (fieldname, svalue, ifield, dtype, card, end))
    return value

def force_integer_or_blank(card: BDFCard, ifield: int, fieldname: str,
                           default: Optional[int]=None) -> Optional[int]:
    """see ``integer_or_blank``"""
    svalue = card.field(ifield)

    if isinstance(svalue, integer_types):
        return svalue
    elif svalue is None:
        return default
    elif '.' in svalue:
        # float
        fvalue = force_double(card, ifield, fieldname)
        # TODO: warn if not a whole number
        return int(fvalue)
    elif isinstance(svalue, str):
        if len(svalue) == 0:
            return default
        elif '.' in svalue or '-' in svalue[1:] or '+' in svalue[1:]:
            dtype = _get_dtype(svalue)
            raise SyntaxError('%s = %r (field #%s) on card must be an integer or blank (not %s).\n'
                              'card=%s' % (fieldname, svalue, ifield, dtype, card))

        try:
            return int(svalue)
        except(ValueError, TypeError):
            dtype = _get_dtype(svalue)
            raise SyntaxError('%s = %r (field #%s) on card must be an integer or blank (not %s).\n'
                              'card=%s' % (fieldname, svalue, ifield, dtype, card))

    # float
    dtype = _get_dtype(svalue)
    raise SyntaxError('%s = %r (field #%s) on card must be an integer (not %s).\n'
                      'card=%s' % (fieldname, svalue, ifield, dtype, card))

def force_double_or_blank(card: BDFCard, ifield: int, fieldname: str,
                          default: Optional[float]=None):
    """see ``double_or_blank``"""
    svalue = card.field(ifield)

    if isinstance(svalue, float_types):
        return svalue
    elif isinstance(svalue, integer_types):
        fvalue = float(svalue)
        warnings.warn(f'{fieldname} = {svalue!r} (field #{ifield}) on card must be a float or blank (not an integer).\n'
                      f'card={card}')
        return fvalue
    elif isinstance(svalue, str):
        try:
            # if it casts as an integer, it's not typed right
            ivalue = int(svalue)
            fvalue = float(ivalue)
            warnings.warn('%s = %r (field #%s) on card must be a float or blank (not an integer) -> %s.\n'
                          'card=%s' % (fieldname, svalue, ifield, fvalue, card))
            return fvalue
        except Exception:
            svalue = svalue.strip().upper()
            if not svalue:
                out = default
            try:
                out = double(card, ifield, fieldname)
            except Exception:
                if svalue == '.':
                    return 0.
                dtype = _get_dtype(svalue)
                raise SyntaxError('%s = %r (field #%s) on card must be a float or blank (not %s).\n'
                                  'card=%s' % (fieldname, svalue, ifield, dtype, card))
            return out
    return default

def force_double_or_string(card: BDFCard, ifield: int, fieldname: str,
                           end: str=''):
    """see ``double_or_string``"""
    svalue = card.field(ifield)

    if isinstance(svalue, float_types):
        return svalue
    elif isinstance(svalue, integer_types):
        fvalue = float(svalue)
        warnings.warn('%s = %r (field #%s) on card must be a float or string (not an integer).\n'
                      'card=%s' % (fieldname, svalue, ifield, card))
        return fvalue

    elif isinstance(svalue, str):
        svalue_nospace = svalue.replace(' ', '')
        if len(svalue_nospace) == 0:
            warnings.warn('%s = %r (field #%s) on card must be a float or string (not an blank).\n'
                          'card=%s' % (fieldname, svalue_nospace, ifield, card))
            raise RuntimeError('no blanks allowed')

        assert len(card.card) > 0, card.card
        card.card[ifield] = svalue_nospace
        try:
            # float
            fvalue = force_double(card, ifield, fieldname)
            return fvalue
        except SyntaxError:
            pass

        svalue_upper = svalue.upper()
        if svalue_upper[0].isdigit() or '.' in svalue_upper or '+' in svalue_upper or '-' in svalue_upper or ' ' in svalue_upper:
            card.card[ifield] = svalue
            dtype = _get_dtype(svalue_upper)
            raise SyntaxError('%s = %r (field #%s) on card must be an float or '
                              'string (without a blank space; not %s).\n'
                              'card=%s%s' % (fieldname, svalue_upper, ifield, dtype, card, end))
        return str(svalue_upper)
        # raise NotImplementedError(svalue)
        #try:
            #ivalue = int(svalue)
            #fvalue = float(ivalue)
            #warnings.warn('%s = %r (field #%s) on card must be a float or blank (not an integer) -> %s.\n'
                          #'card=%s' % (fieldname, svalue, ifield, fvalue, card))
            #return fvalue
        #except Exception:
            #svalue = svalue.strip().upper()
            #if not svalue:
                #return default
            #try:
                #return double(card, ifield, fieldname)
            #except Exception:
                #if svalue == '.':
                    #return 0.
                #dtype = _get_dtype(svalue)
                #raise SyntaxError('%s = %r (field #%s) on card must be a float or blank (not %s).\n'
                                  #'card=%s' % (fieldname, svalue, ifield, dtype, card))
    else:
        dtype = _get_dtype(svalue)
        raise SyntaxError('%s = %r (field #%s) on card must be a float or string (not %s).\n'
                          'card=%s' % (fieldname, svalue, ifield, dtype, card))

    raise SyntaxError('%s = %r (field #%s) on card must be a float or string -> %s.\n'
                      'card=%s' % (fieldname, svalue, ifield, dtype, card))
    #return default

def lax_double_or_blank(card: BDFCard, ifield: int, fieldname: str,
                        default: Optional[float]=None,
                        end: str='') -> float:
    value = integer_double_or_blank(card, ifield, fieldname, default=default)
    if isinstance(value, int):
        value = float(value)
    return value

def parse_components(card: BDFCard, ifield: int, fieldname: str) -> str:
    #assert isinstance(ifield, int), type(ifield)
    #assert isinstance(fieldname, str), type(fieldname)
    svalue = card.field(ifield)

    if isinstance(svalue, integer_types):
        svalue = str(svalue)

    if svalue is None or '.' in svalue:
        dtype = _get_dtype(svalue)
        msg = ('%s = %r (field #%s) on card must be an integer (not %s).\n'
               'card=%s' % (fieldname, svalue, ifield, dtype, card))
        raise SyntaxError(msg)

    svalue = svalue.replace(' ', '')
    try:
        value = int(svalue)
    except ValueError:
        dtype = _get_dtype(svalue)
        msg = ('%s = %r (field #%s) on card must be an integer (not %s).\n'
               'card=%s' % (fieldname, svalue, ifield, dtype, card))
        raise SyntaxError(msg)
    if value > 0 and isinstance(svalue, str):
        if '0' in svalue:
            value2 = str(svalue).replace('0', '')
            msg = ('%s = %r (field #%s) on card must contain 0 or %s (not both).\n'
                   'card=%s' % (fieldname, svalue, ifield, value2, card))
            raise SyntaxError(msg)
    svalue2 = str(value)
    svalue3 = ''.join(sorted(svalue2))
    for i, component in enumerate(svalue3):
        if component not in '0123456':
            msg = ('%s = %r (field #%s) on card contains an invalid component %r.\n'
                   'card=%s' % (fieldname, svalue, ifield, component, card))
            raise SyntaxError(msg)
        if component in svalue3[i + 1:]:
            msg = ('%s = %r (field #%s) on card must not contain duplicate entries.\n'
                   'card=%s' % (fieldname, svalue, ifield, card))
            raise SyntaxError(msg)
    return svalue3

def parse_components_or_blank(card: BDFCard, ifield: int, fieldname: str, default: str='0') -> str:
    #assert isinstance(card, BDFCard), type(card)
    #assert isinstance(ifield, int), type(ifield)
    #assert isinstance(fieldname, str), type(fieldname)
    svalue = card.field(ifield)
    #if isinstance(svalue, integer_types):
        #pass
    if svalue is None:
        return default
    elif '.' in svalue:
        dtype = _get_dtype(svalue)
        msg = ('%s = %r (field #%s) on card must be an integer or blank (not %s).\n'
               'card=%s' % (fieldname, svalue, ifield, dtype, card))
        raise SyntaxError(msg)

    svalue = svalue.replace(' ', '')
    try:
        value = int(svalue)
    except ValueError:
        dtype = _get_dtype(svalue)
        msg = ('%s = %r (field #%s) on card must be an integer or blank (not %s).\n'
               'card=%s' % (fieldname, svalue, ifield, dtype, card))
        raise SyntaxError(msg)
    if value > 0 and isinstance(svalue, str):
        if '0' in svalue:
            value2 = str(svalue).replace('0', '')
            msg = ('%s = %r (field #%s) on card must contain 0 or %s (not both).\n'
                   'card=%s' % (fieldname, svalue, ifield, value2, card))
            raise SyntaxError(msg)
    svalue2 = str(value)
    svalue3 = ''.join(sorted(svalue2))
    for i, component in enumerate(svalue3):
        if component not in '0123456':
            msg = ('%s = %r (field #%s) on card contains an invalid component %r.\n'
                   'card=%s' % (fieldname, svalue, ifield, component, card))
            raise SyntaxError(msg)
        if component in svalue3[i + 1:]:
            msg = ('%s = %r (field #%s) on card must not contain duplicate entries.\n'
                   'card=%s' % (fieldname, svalue, ifield, card))
            raise SyntaxError(msg)
    return svalue3
