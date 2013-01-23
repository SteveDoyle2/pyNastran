# pylint: disable=E1101,C0103,R0902,R0904,R0914
"""
All damper properties are defined in this file.  This includes:
 *   PDAMP
 *   PDAMP5 (not implemented)
 *   PDAMPT
 *   PVISC

All bush properties are BushingProperty and Property objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys


class CardMethods(object):
    def __init__(self, nCardLinesMax=1000):
        self.nCardLinesMax = nCardLinesMax
        pass

    def parseDynamicSyntax(self, key):
        """
        Applies the dynamic syntax for %varName
        @param self the object pointer
        @param key the uppercased key
        @retval value the dynamic value defined by dictOfVars
        @note
          %varName is actually %VARNAME b/c of auto-uppercasing the string,
          so the setDynamicSyntax method uppercases the key prior to this step.
        @see setDynamicSyntax
        """
        #print "*** valueRaw.lstrip() = |%r|" %(valueRaw.lstrip())
        #key = key.lstrip('%%')
        key = key[1:]
        self.log.info("dynamic key = |%r|" % (key))
        #self.dictOfVars = {'P5':0.5,'ONEK':1000.}
        return self.dictOfVars[key]

    def get_value(self, valueRaw, card, debug=False):
        """Converts a value from nastran format into python format."""
        if debug:
            print("v1 = |%s|" % (valueRaw))
        lvalue = valueRaw.lstrip()
        if self._is_dynamic_syntax and '%' in lvalue[0:1]:
            return self.parseDynamicSyntax(valueRaw)
        valueIn = valueRaw.lstrip().rstrip(' *').upper()

        if debug:
            pass
            #print "v2 = |%s|" %(valueIn)
        if len(valueIn) == 0:
            if debug:
                print("BLANK!")
            return None

        if valueIn[0].isalpha():
            if debug:
                print("STRING!")
            return valueIn

        #print "valueIn = |%s|" %(valueIn)
        if ' ' in valueIn:
            msg = ('there are embedded blanks in the field (mixed '
                   'tabs/commas/spaces).\nvalueRaw=|%s| valueIn=|%s| card=%s'
                   % (valueRaw, valueIn, card))
            raise SyntaxError(msg)

        if '=' in valueIn or '(' in valueIn or '*' in valueRaw:
            if debug:
                print("=(! - special formatting")
            return valueRaw.strip()
        #valueIn = valueIn.upper()
        # int, float, string, exponent
        valuePositive = valueIn.strip('+-')
        if debug:
            print("valuePositive=|%r| isDigit=%s" % (valuePositive,
                                                     valuePositive.isdigit()))
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

        #print "valueIn = |%s|" %(valueIn)
        #print "v0 = |%s|" %v0
        if v0 == '-':
            vFactor = -1.
        elif v0 == '+' or v0.isdigit():
            vFactor = 1.
        else:
            msg = ('the only 2 cases for a float/scientific are +/- for v0...'
                   'valueRaw=|%s| v0=|%s| card=%s' % (valueRaw, v0, card))
            raise SyntaxError(msg)

        vm = valueIn.find('-', 1)
            # dont include the 1st character, find the exponent
        vp = valueIn.find('+', 1)
        if vm > 0:
            sline = valueLeft.split('-')
            expFactor = -1.
        elif vp > 0:
            sline = valueLeft.split('+')
            expFactor = 1.
        else:
            msg = ('thought this was in scientific notation, but there is no '
                   'exponent sign...valueRaw=|%s| valueLeft=|%s| card=%s\n'
                   'You also might have mixed tabs/spaces/commas or embedded '
                   'blanks in the field.' % (valueRaw, valueLeft, card))
            raise SyntaxError(msg)

        sline0 = sline[0].rstrip('Dd')
        sline1 = sline[1]
        try:
            s0 = vFactor * float(sline0)
            s1 = expFactor * int(sline1)
        except ValueError:
            msg = ("vm=%s vp=%s valueRaw=|%s| sline0=%s sline1=%s\ncard=%s"
                   % (vm, vp, valueRaw, sline0, sline1, card))
            msg2 = ('cannot parse sline0 into a float and sline1 into an '
                    'integer\n%s\nYou might have mixed tabs/spaces/commas!  '
                    'Fix it!\nfem=%s' % (msg, self.bdf_filename))
            raise SyntaxError(msg2)

        value = s0 * 10 ** (s1)
        #print "valueOut = |%s|" %value

        if debug:
            print("SCIENTIFIC!")
        return value

    def isLargeField(self, card):
        """returns True if the card is in 16-character width fields"""
        star_fields = ['*' in field for field in card]

        if any(star_fields):
            return True
        return False


def parse_csv(sline):
    slineA = sline.split(',')[:9]
    sline2 = [''] * 9
    for (i, s) in enumerate(slineA):
        #print i,s
        sline2[i] = s
    #sline = sline2
    #print "sline2 = ",sline2
    #sline = sline2
    #sline = sline.split(',')[0:9]  # doesnt fill all fields on line
    return sline2


def make_single_streamed_card(log, card, debug=False):
    """
    takes a card that has been split b/c it's a multiline card
    and gets rid of the required blanks in it.
    """
    cardOut = []
    n = 0
    if debug:
        log.debug("card = %s" % (card))
    for (i, field) in enumerate(card):
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
    #    print("  fields = ",collapse(fields))

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


def interpretValue(valueRaw, card='', debug=False):
    """converts a value from nastran format into python format."""
    #debug = True
    if debug:
        print("v1 = |%s|" % (valueRaw))
    #lvalue = valueRaw.lstrip()
    valueIn = valueRaw.lstrip().rstrip(' *').upper()

    if debug:
        pass
        #print "v2 = |%s|" %(valueIn)
    if len(valueIn) == 0:
        if debug:
            print("BLANK!")
        return None

    if valueIn[0].isalpha():
        if debug:
            print("STRING!")
            print("valueStr = %s" % (valueIn))
        return valueIn

    if '=' in valueIn or '(' in valueIn or '*' in valueRaw:
        if debug:
            print("=(! - special formatting")
        return valueRaw.strip()
    #valueIn = valueIn.upper()
    # int, float, string, exponent
    valuePositive = valueIn.strip('+-')
    if debug:
        print("isDigit = %s" % (valuePositive.isdigit()))
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

    #print "valueIn = |%s|" %(valueIn)
    #print "v0 = |%s|" %v0
    if v0 == '-':
        vFactor = -1.
    elif v0 == '+' or v0.isdigit():
        vFactor = 1.
    else:
        msg = ('the only 2 cases for a float/scientific are +/- for v0...'
               'valueRaw=|%s| v0=|%s| card=%s' % (valueRaw, v0, card))
        raise SyntaxError(msg)

    vm = valueIn.find(
        '-', 1)  # dont include the 1st character, find the exponent
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


def stringParser(stringIn):
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
            msg = "s=|%r|" % (s)
            raise SyntaxError(msg)

        #print "s=%s stringIn[i-1]=%s" % (state,typeCheck[i-1])
        #print "i=%s s=%s typeCheck=%s" % (i,s,typeCheck)
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
        #print "lString = ",lString
        out = ''.join(lString)
        print("out = %s" % (out))
        return float(out)
    else:
        #print "string = ",stringIn
        #print "typeCheck = ",typeCheck
        #return 'string'
        return stringIn

    print("typeCheck = |%s|" % (typeCheck))
    raise RuntimeError('error parsing a card...this should never happen...')

if __name__ == '__main__':
    print(stringParser('123'))
    print(stringParser('+123'))
    print(stringParser('.234'))
    print(stringParser('+.234'))
    print(stringParser('-.234'))
    print(stringParser('1+5'))
    print("abc = |%s|" % (stringParser('abc')))
    print("eeg = |%s|" % (stringParser('eeg')))
    #print("e1 = |%s|" %(stringParser('\T')))
    print(stringParser('.e1'))