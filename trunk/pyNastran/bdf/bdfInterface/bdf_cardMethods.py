# pylint: disable=E1101,C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys


class CardMethods(object):
    def __init__(self, nCardLinesMax=1000):
        self.nCardLinesMax = nCardLinesMax

    def _make_lines_pack(self, debug=False):
        emptyLines = 0
        if not self.linesPack:
            return ['']

        while len(self.linesPack[-1]) < self.nCardLinesMax:
            line = self.infilesPack[-1].readline()
            line = line.split('$')[0].rstrip('\n\r\t ')
            if('$' not in line and len(line) > 0):
                if debug:
                    print("line = |%r|" % (line))
                self.linesPack[-1].append(line)
            else:
                emptyLines += 1
            ###
            if emptyLines == 500:
                break
        return self.linesPack[-1]

    def update_card_lines(self, lines):
        """expands a card with tabs in it"""
        lines2 = []
        for line in lines:
            if '\t' in line:
                #raise SyntaxError('lines are ambiguous when there are tabs...'
                #                  'fix them...line=|%r|' %(line))
                if ',' in line:
                    #expandTabCommas(line2)
                    raise SyntaxError('tabs and commas in the same line are '
                                      'not supported...line=|%r|' % (line))
                line = line.expandtabs()
            ###
            lines2.append(line)
        return lines2

    def _get_card(self, debug=False):
        """gets a single unparsed card"""
        #debug = True

        linesPack = self._make_lines_pack(debug=debug)
        #if debug:
            #print '-------------------------------------------------'
            #print "pack = \n",'\n'.join(linesPack)

        # fix for unhandled card at end of deck
        #if linesPack[0][0] in ['+','*',' ']:
            #self.doneReading = True
            #return(None, None, None)

        if linesPack == []:
            self.close_file()
            linesPack = self._make_lines_pack(debug=debug)

        tempcard = [linesPack[0]]
        #if 'CQUAD4' in tempcard[0]:
        #    debug = True
        i = 1
        #if emptyLines<50:
        try:
            (i, tempcard) = self._get_multi_line_card(i, tempcard, debug=debug)
        except IndexError:
            #try:
            #    tempcard = self._get_multi_line_card(i, tempcard, debug=True)
            #except IndexError:
            #    print("workscard = %s" %(tempcard))
            #    print("")
            #    #raise
            #    pass
            #raise
            #print "error..."
            #print "line = ",self.lines[i]
            self.doneReading = True
            #raise
        ###
        #print "i = ",i

        #try:
        linesPack[:] = linesPack[i:]
        #except IndexError:
        #    self.lines = []

        #print "tempcard = ",''.join(tempcard)

        tempcard = self.update_card_lines(tempcard)
        upperCard = [line.upper() for line in tempcard]
        cardName = self._get_card_name(upperCard)
        #print "|%s|" %(cardName)
        #if cardName=='CTRIA3':
        if debug:
        #if 1:
            self.log.debug("cardName  = |%s|" % (cardName))
            self.log.debug("upperCard = |%s|" % (upperCard))
            self.log.debug("tempcard  = |%s|" % (tempcard))
            self.log.debug("-------\n")
        self._increaseCardCount(cardName)
        return (tempcard, upperCard, cardName)

    def _increaseCardCount(self, cardName):
        """
        Used for testing to check that the number of cards going in is the
        same as each time the model is read verifies proper writing of cards
        @warning
            this wont guarantee proper reading of cards, but will help
        """
        if cardName == '':  # stupid null case
            return

        if cardName in self.cardCount:
            self.cardCount[cardName] += 1
        else:
            self.cardCount[cardName] = 1

    def _get_multi_line_card(self, i, tempcard, isCSV=False, debug=False):
        iline = self.linesPack[-1][i].rstrip()
        #while iline=='':
        #    i+=1
        #    iline = self.lines[i].rstrip()

        sCardName = iline[0:8].rstrip()  # trying to find if it's blank...
        isNotDone = len(iline) > 0 and (iline.rstrip()[0] in ['*', '+', ',','\t']
                                        or sCardName == ''
                                        or not sCardName[0].isalpha())

    #debug = True
        if debug:
            print("_get_multi_line_card...i=%s" % (i))
            print("tempcard1 = %s" % (tempcard))

            self.log.debug("CRITERIA A")
            self.log.debug("  iline      = |%r|" % (iline))
            self.log.debug("  len(iline) = %-10s -> len(iline)>0         = %s" % ('|' + str(len(iline)) + '|', str(len(iline) > 0)))
            self.log.debug("  iline[0]   = %-10s -> line[0] in [*,+,','] = %s" % ('|' + iline[0] + '|', iline.strip()[0] in ['*', '+', ',']))
            self.log.debug("  sCardName  = %-10s -> name=''              = %s" % ('|' + sCardName + '|', sCardName == ''))
            self.log.debug("  iline = |%s|" % (iline))
            self.log.debug("isNotDone A = %s" % (isNotDone))

        while(isNotDone):
            if debug:
                print("not done...i=%s" % (i))
            tempcard.append(iline)
            i += 1
            #if debug:
            #    print(card)
            if 'ENDDATA' in iline:
                #print "********"
                self.log.debug('found ENDDATA')
                self.doneReading = True
                break
            if i + 1 == len(self.linesPack[-1]):
                if debug:
                    self.log.debug("breaking b/c empty pack???")
                break
            if i == len(self.linesPack[-1]):  # pre-catching the raise...
                self.doneReading = True
                break
            #ilineA = self.linesPack[-1]
            #print "iLineA = |%r|" %(ilineA)
            iline = self.linesPack[-1][i]
            #try:
            #    iline = self.linesPack[-1][i]
            #except IndexError:
            #    iline = ''

            sCardName = iline[0:8].strip()  # trying to find if it's blank...
            #slot0 = iline[0:8]
            #if '\t' in slot0:
            #    slot0 = slot0.expandtabs()
            #sCardName = slot0.strip()  # trying to find if it's blank...
            isNotDone = len(iline) > 0 and (iline.rstrip()[0] in ['*', '+', ',','\t']
                                            or sCardName == ''
                                            or not sCardName[0].isalpha())
            if debug:
                print(tempcard)
                self.log.debug("CRITERIA B")
                self.log.debug("  iline       = |%r|" % (iline))
                self.log.debug("  len(iline) = %-10s -> len(iline)>0         = %s" % ('|' + str(len(iline)) + '|', str(len(iline) > 0)))
                self.log.debug("  iline[0]   = %-10s -> line[0] in [*,+,','] = %s" % ('|' + iline[0] + '|', iline.strip()[0] in ['*', '+', ',']))
                self.log.debug("  sCardName  = %-10s -> name=''              = %s" % ('|' + sCardName + '|', sCardName == ''))
                self.log.debug("  isNotDone B = %s" % isNotDone)
        ###
        #if debug:
        #self.log.debug("tempcard2 = |%s|" %(tempcard))
            #print ""
        if debug:
            print("done...i=%s" % (i))
            print("")
        return (i, tempcard)

    def isLargeField(self, card):
        """returns True if the card is in 16-character width fields"""
        starField = ['*' in field for field in card]
        #print "starField = ",starField

        if any(starField):
            return True
        return False

    def processCard(self, tempcard, debug=False):
        """
        takes a list of strings and returns a list with the
        proper value in the fields of the list
        """
        #debug = True
        card = []
        #isLargeField = self.isLargeField(tempcard)
        #print "*** isLargeField = ",isLargeField

        #print "tempcard = ",tempcard
        for (i, line) in enumerate(tempcard):
            #print "line = ",line
            isLargeField = self.isLargeField(line)
            #print "i = ",i
            if debug:
                self.log.debug("  line  = |%r|" % (line))
            sline = line[0:72]
            if not(sline):
                break
            if debug:
                self.log.debug("  line2 = |%r|" % (sline))

            if ',' in sline:  # CSV - doesnt support large field CSV cards which I'd never used...
                try:
                    sline = parse_csv(sline)
                except:
                    print("cant parse sline=%s" % (sline))
                    raise
                #print "sline  = ",sline
                #self.log.debug("sline = %s" %(sline))
            else:  # standard
                sline = nastran_split(self.log, sline, isLargeField,
                                      debug=debug)
            #name = sline[0]
            #nFields = len(sline)
            #print "sline = ",sline

            for (fieldCounter, valueIn) in enumerate(sline):
                #if fieldCounter==8:
                    #print "**type(value) = ",type(value)
                    #break
                #if debug:
                    #print "type(value) = ",type(value)
                    #print ""
                if i > 0 and fieldCounter == 0:  # blank leading field
                    pass
                else:
                    #debug = True
                    try:
                        value = self.getValue(valueIn, sline, debug=debug)
                    except:
                        self.log.error("card = |%r|" % (card))
                        raise
                    card.append(value)
                    #print "fieldCounter=%s valueIn=%s value=%s type=%s" %(fieldCounter,valueIn,value,type(value))
                ###
            ###
            #print "cardEnd temp = ",card
        ###
        #print "cardOut&& = ",card
        #if debug:
            #self.log.debug("  sline2 = %s" %(card))
            #self.log.debug("  sline2 = %s" %(collapse(card)))
        #return make_single_streamed_card(self.log, card)
        return card

    def expandTabCommas(self, line):
        """
        The only valid tab/commas format in nastran is having the
        first field be a tab and the rest of the fields be separated by commas.
        @param self the object pointer
        @param line a BDF line
        """
        fields = []
        isWord = True
        field = ''
        i = 0
        for (i, letter) in enumerate(line):
            if letter not in ['\t', ',', ' ']:  # tab or comma
                if isWord:
                    field += letter
                else:
                    isWord = True
                    fields.append(field)
                    field = ''
                    field += letter
            elif letter == ' ' or letter == ',':
                isWord = False
                break
        #fields.append(field)
        sline = [field] + line[i:72].split(',')
        print("expandTabCommas = |%r|" % (sline))
        return fields

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

    def getValue(self, valueRaw, card, debug=False):
        """converts a value from nastran format into python format."""
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
            msg = ('there are embedded blanks in the field (mixed tabs/commas/spaces).\n'
                   'valueRaw=|%s| valueIn=|%s| card=%s' % (valueRaw, valueIn,
                                                           card))
            raise SyntaxError(msg)

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
                    'Fix it!' % (msg))
            raise SyntaxError(msg2)

        value = s0 * 10 ** (s1)
        #print "valueOut = |%s|" %value

        if debug:
            print("SCIENTIFIC!")
        return value


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
    #return collapse(cardOut)


def nastran_split(log, line, isLargeField, debug=False):
    """
    Splits a single BDF line into large field or small field format
    @param self the object pointer
    @param line the BDF line
    @param a flag indicating small/large field format (True/False)
    @param debug extra developer debug
    @retval fields the 9 (small) or 5 (large) fields for the line
    @note CSV Format is handled by parse_csv
    @note tabs are handled prior to running this
    """
    if debug:
        print("isLargeField = %s" % (isLargeField))
    if isLargeField:
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
        msg = 'the only 2 cases for a float/scientific are +/- for v0...valueRaw=|%s| v0=|%s| card=%s' % (valueRaw, v0, card)
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
        msg = 'thought this was in scientific notation, but i cant find the exponent sign...valueRaw=|%s| valueLeft=|%s| card=%s\nYou also might have mixed tabs/spaces/commas.' % (valueRaw, valueLeft, card)
        raise SyntaxError(msg)

    try:
        s0 = vFactor * float(sline[0])
        s1 = expFactor * int(sline[1])
    except ValueError:
        msg = "vm=%s vp=%s valueRaw=|%s| sline=%s" % (vm, vp, valueRaw, sline)
        raise SyntaxError('cannot parse sline[0] into a float and sline[1] into an integer\n%s\nYou HAVE mixed tabs/spaces/commas!  Fix it!' % (msg))

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
        ###

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
        ###
    ###
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
