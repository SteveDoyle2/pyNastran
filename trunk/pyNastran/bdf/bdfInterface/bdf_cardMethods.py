import sys
#from BDF_Card import collapse
from pyNastran.bdf.errors import *

class cardMethods(object):
    def __init__(self):
        pass

    def makeLinesPack(self,debug=False):
        emptyLines=0
        if not self.linesPack:
            return ['']

        while len(self.linesPack[-1])<40:
            line = self.infilesPack[-1].readline()
            line = line.split('$')[0].rstrip('\n\r\t ')
            if '\t' in line:
                raise TabCharacterError('lines are ambiguous when there are tabs...fix them...line=|%r|' %(line))
            if('$' not in line and len(line)>0):
                if debug:
                    print "line = |%r|" %(line)
                self.linesPack[-1].append(line)
            else:
                emptyLines += 1
            ###
            if emptyLines==50:
                break
        return self.linesPack[-1]

    def getCard(self,debug=False):
        """gets a single unparsed card"""
        #debug = True
        
        linesPack = self.makeLinesPack(debug=debug)
        tempcard = [linesPack[0]]
        #if 'CQUAD4' in tempcard[0]:
        #    debug = True
        i=1
        #if emptyLines<50:
        try:
            (i,tempcard) = self.getMultiLineCard(i,tempcard,debug=debug)
        except IndexError:
            #try:
            #    tempcard = self.getMultiLineCard(i,tempcard,debug=True)
            #except IndexError:
            #    print "workscard = ",tempcard
            #    print ""
            #    #raise
            #    pass
            pass
            #raise
            #print "error..."
            #print "line = ",self.lines[i]
            self.doneReading=True
            #raise
        ###
        
        #try:
        linesPack[:] = linesPack[i:]
        #except IndexError:
        #    self.lines = []
        
        #print "tempcard = ",''.join(tempcard)
        
        upperCard = [line.upper() for line in tempcard]
        cardName = self.getCardName(upperCard)
        #print "|%s|" %(cardName)
        #if cardName=='CTRIA3':
        if debug:
        #if 1:
            self.log.debug("cardName  = |%s|" %(cardName))
            self.log.debug("upperCard = |%s|" %(upperCard))
            self.log.debug("tempcard  = |%s|" %(tempcard))
            self.log.debug("-------\n")
        self._increaseCardCount(cardName)
        return (tempcard,upperCard,cardName)

    def _increaseCardCount(self,cardName):
        """
        used for testing to check that the number of cards going in is the same as each time the model is read
        verifies proper writing of cards
        @warning
            this wont guarantee proper reading of cards, but will help
        """
        if cardName=='': # stupid null case
            return

        if cardName in self.cardCount: # dict
            self.cardCount[cardName] += 1
        else:
            self.cardCount[cardName]  = 1
        ###

    def getMultiLineCard(self,i,tempcard,isCSV=False,debug=False):
        if debug:
            print "tempcard1 = ",tempcard
        iline = self.linesPack[-1][i].rstrip()
        #while iline=='':
        #    i+=1
        #    iline = self.lines[i].rstrip()

        sCardName = iline[0:8].strip()  # trying to find if it's blank...
        isNotDone = len(iline)>0 and (iline[0] in ['*','+',','] or sCardName=='')
        if debug:
            self.log.debug("  len(iline) = |%s|" %(len(iline)))
            self.log.debug("  iline[0]   = |%s|" %(iline[0]  ))
            self.log.debug("  sCardName  = |%s|" %(sCardName ))
            self.log.debug("  iline = |%s|" %(iline))
            print ""
            print "isNotDone A = %s" %(isNotDone)
        
        while(isNotDone):
            tempcard.append(iline)
            i+=1
            #if debug:
            #    print card
            if 'ENDDATA' in iline:
                print "********"
                self.doneReading=True
                break
            if i+1==len(self.linesPack[-1]):
                break
            iline = self.linesPack[-1][i]
            #try:
            #    iline = self.linesPack[-1][i]
            #except IndexError:
            #    iline = ''

            sCardName = iline[0:8].strip()  # trying to find if it's blank...
            isNotDone = len(iline)>0 and (iline[0] in ['*','+',','] or sCardName=='')
            if debug:
                self.log.debug("CRITERIA")
                self.log.debug("iline       = |%r|" %(iline))
                self.log.debug("sCardName   = %s" %sCardName)
                self.log.debug("len(iline)  = %s" %len(iline))
                self.log.debug("iline[0]    = %s" %iline[0])
                self.log.debug("isNotDone B = %s" %isNotDone)
        ###
        #if debug:
        #self.log.debug("tempcard2 = |%s|" %(tempcard))
            #print ""
        #sys.exit('asdf')
        return (i,tempcard)
    
    def nastranSplit(self,line,isLargeField,debug=False):
        if debug:
            print "isLargeField = %s" %(isLargeField)
        if isLargeField:
            #print "large"
            #print "line = |%s|" %(line)
            fields = [line[0:8],line[8:24],line[24:40],line[40:56],line[56:72]]
        else:
            #print "small"
            #fields = [line[0 :8 ],line[8 :16],line[16:24],line[24:32],line[32:40],
            #          line[40:48],line[48:56],line[56:64],line[64:72],line[72:80]]
            fields = [line[0 :8 ],line[8 :16],line[16:24],line[24:32],line[32:40],
                      line[40:48],line[48:56],line[56:64],line[64:72]]
        #if debug:
        #    print "  fields = ",collapse(fields)
        
        fields2 = []
        for i,rawField in enumerate(fields):
            field = rawField.strip()
            #if debug:
            #if (i==9) and field=='' or field=='*' or field=='+':
                #print "skipping * or + or empty field"
            #    pass
            #else:
            if debug:
                self.log.debug("i=%s rawField=|%s| field=|%s|" %(i,rawField,field))
            fields2.append(field)
        return fields2

    def isLargeField(self,card):
        """returns True if the card is in 16-character width fields"""
        starField = ['*' in field for field in card]
        #print "starField = ",starField
        
        if any(starField):
            return True
        return False

    def processCard(self,tempcard,debug=False):
        """
        takes a list of strings and returns a list with the 
        proper value in the fields of the list
        """
        #debug = True
        card = []
        #isLargeField = self.isLargeField(tempcard)
        #print "*** isLargeField = ",isLargeField

        #print "tempcard = ",tempcard
        for i,line in enumerate(tempcard):
            isLargeField = self.isLargeField(line)
            #print "i = ",i
            if debug:
                self.log.debug("  line  = |%r|" %(line))
            sline = line[0:73]
            if not(sline):
                break
            if debug:
                self.log.debug("  line2 = |%r|" %(sline))

            if ',' in sline:  #CSV
                sline = sline.split(',')[0:9]
                #self.log.debug("sline = %s" %(sline))
            else: # standard
                sline = self.nastranSplit(sline,isLargeField,debug=debug)
            #name = sline[0]
            #nFields = len(sline)
            #print "sline = ",sline
            
            for (fieldCounter,valueIn) in enumerate(sline):
                #if fieldCounter==8:
                #    print "**type(value) = ",type(value)
                #    break
                    #sys.exit(12131)
                #if debug:
                #    print "type(value) = ",type(value)
                #    print ""
                if i>0 and fieldCounter==0: # blank leading field
                    pass
                else:
                    value = self.getValue(valueIn,debug=debug)
                    card.append(value)
                    #print "fieldCounter=%s valueIn=%s value=%s type=%s" %(fieldCounter,valueIn,value,type(value))
                ###
            ###
            #print "cardEnd temp = ",card
        ###
        #print "cardOut&& = ",card
        if debug:
            self.log.debug("  sline2 = %s" %(card))
            #self.log.debug("  sline2 = %s" %(collapse(card)))
        #return self.makeSingleStreamedCard(card)
        return card
        
    def makeSingleStreamedCard(self,card,debug=False):
        """
        takes a card that has been split b/c it's a multiline card
        and gets rid of the required blanks in it.
        """
        cardOut = []
        n=0
        if debug:
            self.log.debug("card = %s" %(card))
        for i,field in enumerate(card):
            if n-9==0:
                pass
            elif n-10==0:
                n-=10
            else:
                cardOut.append(field)
            n+=1
        ###
        #print "cardOut = ",cardOut
        return cardOut
        #return collapse(cardOut)
        
    def getValue(self,valueRaw,debug=False):
        """converts a value from nastran format into python format."""
        if debug:
            print "v1 = |%s|" %(valueRaw)
        valueIn = valueRaw.lstrip().rstrip(' *').upper()
        
        if debug:
            pass
            #print "v2 = |%s|" %(valueIn)
        if len(valueIn)==0:
            if debug:
                print "BLANK!"
            return None

        if '=' in valueIn or '(' in valueIn or '*' in valueRaw:
            if debug:
                print "=(! - special formatting"
            return valueRaw.strip()
        #valueIn = valueIn.upper()
        # int, float, string, exponent
        valuePositive = valueIn.strip('+-')
        if debug:
            print "isDigit = ",valuePositive.isdigit()
        if valuePositive.isdigit():
            if debug:
                print "INT!"
            return int(valueIn)
        try:
            value = float(valueIn)
            if debug:
                print "FLOAT!"
            return value
        except:
             pass

        #if('=' in valueIn or '(' in valueIn or ')' in valueIn):
        #    print "=()!"
        #    return valueIn

        noED = list(set(valueIn)-set('ED 1234567890+-')) # if there are non-floats/scientific notation -> string
        word = ''.join(noED)
        #print "word=|%s|" %word
        if word.isalpha():
            if debug:
                print "WORD!"
            return valueIn

        v0 = valueIn[0]
        if '-'==v0 or '+'==v0:
            valueLeft = valueIn[1:] # truncate the sign for now
        else:
            v0 = '+' # inplied positive value
            valueLeft = valueIn

        #print "valueIn = |%s|" %(valueIn)
        #print "v0 = |%s|" %v0
        if v0=='-':
            vFactor=-1.
        elif v0=='+' or v0.isdigit():
            vFactor=1.
        else:
            msg = 'the only 2 cases for a float/scientific are +/- for v0...valueRaw=|%s| v0=|%s|' %(valueRaw,v0)
            raise FloatScientificParseError(msg)

        vm = valueIn.find('-',1) # dont include the 1st character, find the exponent
        vp = valueIn.find('+',1)
        if vm>0:
            sline = valueLeft.split('-')
            expFactor = -1.
        elif vp>0:
            sline = valueLeft.split('+')
            expFactor = 1.
        else:
            msg = 'thought this was in scientific notation, but i cant find the exponent sign...valueRaw=|%s| valueLeft=|%s|' %(valueRaw,valueLeft)
            raise ScientificParseError(msg)

        s0 = vFactor*float(sline[0])
        s1 = expFactor*int(sline[1])
        #except:
        #    print "vm=%s vp=%s valueRaw=|%s| sline=|%s|" %(vm,vp,valueRaw,sline)

        value = s0*10**(s1)
        #print "valueOut = |%s|" %value
        
        if debug:
            print "SCIENTIFIC!"
        return value
    

def stringParser(stringIn):
    typeCheck = ''
    n=0
    for i,s in enumerate(stringIn):
        if s in "+-":   state='+'
        elif s==" ":      state=' '
        elif s==".":      state='.'
        elif s in "eEdD": state='e'
        elif s.isdigit(): state='1'
        elif s.isalpha() or s in "()*/=]['\"": return 'string' # string character
        else:
            msg = "s=|%r|" %(s)
            raise Exception(msg)
        ###
        
        #print "s=%s stringIn[i-1]=%s" %(state,typeCheck[i-1])
        #print "i=%s s=%s typeCheck=%s" %(i,s,typeCheck)
        if i==0:
            typeCheck += state
            n+=1
        elif typeCheck[n-1]!=state:
            typeCheck += state
            n+=1
        elif state in 'e .+': # double e, space, dot, plus
            return 'string'
        ###
    ###
    if typeCheck==' ':
        return None

    typeCheck = typeCheck.strip()
    if typeCheck in ['1','+1']: # integer
        return int(stringIn)

    elif typeCheck in [ '1.', '1.1', '.1',  # float
                       '+1.','+1.1','+.1']:
        return float(stringIn)

    elif typeCheck in [ '1.1e1', '1.1e+1', '1.e1', '1.e+1', # python scientific
                       '+1.1e1','+1.1e+1','+1.e1','+1.e+1',
                       '.1e1','.1e+1','+.1e1','+.1e+1',]:
        return float(stringIn)

    elif typeCheck in ['1+1','+1+1','.1+1','+.1+1']: # nastran scientific
        stringReversed = stringIn[::-1]
        i = stringReversed.index('+')
        lString = list(stringIn)
        lString.insert(-i-1,'e')
        #print "lString = ",lString
        out = ''.join(lString)
        print "out = ",out
        return float(out)
    else:
        #print "string = ",stringIn
        #print "typeCheck = ",typeCheck
        #return 'string'
        return stringIn
    
    print "typeCheck = |%s|" %(typeCheck)
    raise Exception('this should never happen...')

if __name__=='__main__':
    print stringParser('123')
    print stringParser('+123')
    print stringParser('.234')
    print stringParser('+.234')
    print stringParser('-.234')
    print stringParser('1+5')
    print "abc = |%s|" %(stringParser('abc'))
    print "eeg = |%s|" %(stringParser('eeg'))
    #print "e1 = |%s|" %(stringParser('\T'))
    print stringParser('.e1')
 
