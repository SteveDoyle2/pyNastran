# pylint: disable=C0103,R0911,R0912,R0914,R0915
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import warnings


def interpretValue(valueRaw, card='', debug=False):
    """
    @see interpret_value
    @warning will be removed after v0.7 in favor of interpret_value
    """
    warnings.warn('interpretValue has been deprecated; use '
                  'interpret_value', DeprecationWarning, stacklevel=2)
    return interpret_value(valueRaw, card, debug)


def isLargeField(card):
    """
    @returns True if the card is in 16-character width fields
    """
    star_fields = ['*' in field for field in card]
    if any(star_fields):
        return True
    return False


def parse_csv(sline):
    """
    Parses a Nastran single-precision formatted CSV line
    """
    slineA = sline.split(',')[:9]
    sline2 = [''] * 9

    for (i, s) in enumerate(slineA):
        sline2[i] = s
    #sline = sline2
    #print "sline2 = ",sline2
    #sline = sline2
    #sline = sline.split(',')[0:9]  # doesnt fill all fields on line
    return sline2


def make_single_streamed_card(log, card, debug=False):
    """
    Takes a card that has been split b/c it's a multiline card
    and gets rid of the required blanks in it.
    """
    cardOut = []
    n = 0
    if debug:
        log.debug("card = %s" % card)
    for field in card:
        if n - 9 == 0:
            pass
        elif n - 10 == 0:
            n -= 10
        else:
            cardOut.append(field)
        n += 1
    #print "cardOut = ",cardOut
    return cardOut


def nastran_split(log, line, is_large_field, debug=False):
    """
    Splits a single BDF line into large field or small field format
    @param line
      the BDF line
    @param is_large_field
      a flag indicating small/large field format (True/False)
    @param
      debug extra developer debug

    @retval fields
      the 9 (small) or 5 (large) fields for the line

    @note CSV Format is handled by parse_csv
    @note tabs are handled prior to running this
    """
    if debug:
        print("is_large_field = %s" % (is_large_field))
    if is_large_field:
        fields = [line[0:8], line[8:24], line[24:40], line[40:56], line[56:72]]
    else:  # small field
        fields = [line[0:8], line[8:16], line[16:24], line[24:32],
                  line[32:40], line[40:48], line[48:56], line[56:64],
                  line[64:72]]
    #if debug:
    #    print("  fields = ", collapse(fields))

    fields2 = []
    for (i, rawField) in enumerate(fields):
        field = rawField.strip()
        #if debug:
        #if (i==9) and field=='' or field=='*' or field=='+':
            #print "skipping * or + or empty field"
        #    pass
        #else:
        if debug:
            log.debug("i=%s rawField=|%s| field=|%s|" % (i, rawField, field))
        fields2.append(field)
    return fields2


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