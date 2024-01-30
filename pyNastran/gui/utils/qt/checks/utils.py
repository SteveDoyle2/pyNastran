from functools import wraps
import numpy as np

def check_locale_float(text: str) -> tuple[float, bool]:
    """
    This is intended as a locale safe way to parse a float.
    It supports various types of floats.

    Example
    -------
    >>> _check_float('3.14')
    3.14, True
    >>> _check_float('2,557')
    2,557, True
    """
    value = np.nan
    is_valid = False
    try:
        value = float(text)
        is_valid = True
    except ValueError:
        #'2,557'
        if ',' not in text:  # not valid
            pass
        elif '.' in text and ',' in text:  # might be valid; 1,234.56  -> 1.234,56
            # not supported
            pass
        else:
            text2 = text.replace(',', '.')
            value = float(text2)
            is_valid = True
    except Exception:
        pass
    return value, is_valid


def is_ranged_value(value: float,
                    min_value=None,
                    max_value=None,
                    min_inclusive: bool=True,
                    max_inclusive: bool=True) -> bool:
    """
    Parameters
    ----------
    value : float
        float : the value as a float
    min_value / max_value : float / None
        float : the constraint is active
        None : the constraint is inactive
    min_inclusive / max_inclusive; bool; default=True
        flips [min_value, max_value] to:
          - (min_value, max_value)
          - [min_value, max_value)
          - (min_value, max_value]

    Returns
    -------
    is_ranged : bool
        is the value in range

    """
    is_ranged = True
    if min_value is not None:
        #print("min:")
        if min_inclusive:
            # ( is the goal
            if value < min_value:
                is_ranged = False
                #print('  min_exclusive, %s %s' % (value, min_value))
            #else:
                #print('  passed minA=%s' % value)
        else:
            # [ is the goal
            if value <= min_value:
                is_ranged = False
                #print('  min_inclusive, %s %s' % (value, min_value))
            #else:
                #print('  passed minB=%s' % value)
    #else:
        #print('no limit on min')

    if max_value is not None:
        #print("max:")
        if max_inclusive:
            # ] is the goal
            if value > max_value:
                #print('  max_exclusive, %s %s' % (value, max_value))
                is_ranged = False
            #else:
                #print('  passed maxA=%s' % value)
        else:
            # ) is the goal
            if value >= max_value:
                is_ranged = False
                #print('  max_inclusive, %s %s' % (value, max_value))
            #else:
                #print('  passed maxB=%s' % value)
    #else:
        #print('no limit on max')
    return is_ranged

def nocrash_str_bool(func):
    @wraps(func)
    def wrapper(_str):
        try:
            out = func(_str)
        except Exception:
            #print('dont crash...')
            out = (_str, False)
        return out
    return wrapper

@nocrash_str_bool
def check_format_str(text: str) -> tuple[str, bool]:
    """
    Checks a QLineEdit string formatter

    Parameters
    ----------
    text : str
        a QLineEdit containing a string formatter like:
        {'%s', '%i', '%f', '%g', '%.3f', '%e'}

    Returns
    -------
    text : str / None
        str : the validated text of the QLineEdit
        None : the format is invalid
    is_valid : bool
        The str/None flag to indicate if the string formatter is valid

    >>> check_format_str("")
    >>> check_format_str(".3E")
    >>> check_format_str(".3g")
    >>> check_format_str(".4f")
    >>> check_format_str("08,.1f")
    """

    text = text.strip()
    is_valid = True
    is_int_fmt = False

    # basic length checks
    if len(text) < 2:
        is_valid = False
    elif 's' in text.lower():
        is_valid = False
    elif '%' not in text[0]:
        is_valid = False
    elif text[-1].lower() in ['i', 'd']:
        is_int_fmt = True
    elif text[-1].lower() not in ['g', 'f', 'i', 'e', 'd']:
        is_valid = False

    # the float formatter handles ints/floats?
    if is_int_fmt:
        try:
            text % 1
            text % .2
            text % 1e3
            text % -5.
            text % -5
        except (ValueError, TypeError):
            is_valid = False
    else:
        text_no_percent = text[1:]
        try:
            format(1, text_no_percent)
            format(.2, text_no_percent)
            format(1e3, text_no_percent)
            format(-5., text_no_percent)
            format(-5, text_no_percent)
        except ValueError:
            is_valid = False
        except TypeError:
            is_valid = False

    if not is_valid:
        return text, is_valid

    # The float formatter isn't supposed to handle strings?
    # Doesn't this break %g?
    #     '%g' % 's'
    #     Traceback (most recent call last):
    #       Python Shell, prompt 3, line 1
    #     builtins.TypeError: must be real number, not str
    # apparently not...
    #
    # The point of this is to make sure that string formatting is invalid.
    # This function is for numbers!
    try:
        text % 's'
        is_valid = False
    except (ValueError, TypeError):
        pass
    return text, is_valid
