from pyNastran.bdf.bdfInterface.bdf_cardMethods import interpretValue
#from pyNastran.bdf.caseControlDeck import CaseControlDeck


def parseSetSline(listA):
    print "listA = ", listA
    listB = []
    for spot in listA:
        spot = spot.strip()
        print "spot = ", spot
        if spot == '':
            pass
        elif ' ' in spot:
            sline = spot.split(' ')
            print "sline = ", sline
            if sline[1] == 'THRU':
                if 'BY' in sline:
                    by = sline[4]
                else:
                    by = 1
                #print "BY = %s" %(by)
                vals = set([])
                startValue = interpretValue(sline[0])
                endValue = interpretValue(sline[2]) + 1
                for i in xrange(startValue, endValue, by):
                    vals.add(i)
                #print "vals = ",vals
                if 'EXCEPT' in sline:
                    iExcept = sline.index('EXCEPT')
                    #print "e = ",sline[iExcept]
                    excepted = int(sline[iExcept + 1])
                    vals.remove(excepted)
                vals = list(vals)
                listB += vals
                print "vals = ", vals
            else:
                print "sline = ", sline
        else:
            #print "spot = %s" %(spot)
            if '/' in spot:
                listB.append(spot)
            else:
                listB.append(interpretValue(spot))
    return listB


def parseSetType(i, line, lines, key, value):
    (key, ID) = key.split()
    key = key + ' ' + ID

    #print('SET-type key=%s ID=%s' %(key, ID))
    sline = value.strip(' ,').split(',')  # float/int values
    fivalues = parseSetSline(sline)

    ## @todo should be more efficient multiline reader...
    # read more lines....
    if line[-1].strip() == ',':
        i += 1
        #print "rawSETLine = |%r|" %(lines[i])
        while 1:
            if ',' == lines[i].strip()[-1]:
                sline = lines[i][:-1].strip().split(',')
                fivalues += parseSetSline(sline)
            else:  # last case
                sline = lines[i].strip().split(',')
                fivalues += parseSetSline(sline)
                #print "fivalues last = i=%s |%r|" %(i,lines[i])
                i += 1
                break
            i += 1
        ###
    ###
    #print "len(fivalues) = ",len(fivalues)
    value = fivalues

    options = int(ID)  # needed a place to put it...
    paramType = 'SET-type'
    return (i, key, value, options, paramType)


def _parseEntry(lines):
    i = 0
    options = []
    value = None
    key = None
    paramType = None

    line = lines[i]
    #print line
    #print "*****lines = ",lines

    equalsCount = 0
    for letter in line:
        if letter == '=':
            equalsCount += 1
    lineUpper = line.upper()

    if lineUpper.startswith('SUBCASE'):
        #print "line = |%r|" %(line)
        line2 = line.replace('=', '')
        sline = line2.split()
        if len(sline) != 2:
            msg = "trying to parse |%s|..." % (line)
            raise SyntaxError("Invalid Subcase: %s", (msg)) 
        (key, iSubcase) = sline
        #print "key=|%s| iSubcase=|%s|" %(key,iSubcase)
        value = int(iSubcase)
        #self.iSubcase = int(iSubcase)
        paramType = 'SUBCASE-type'
    elif (lineUpper.startswith('LABEL') or lineUpper.startswith('SUBTITLE')
          or lineUpper.startswith('TITLE')):
        eIndex = line.index('=')
        key = line[0:eIndex].strip()
        value = line[eIndex + 1:].strip()
        options = []
        paramType = 'STRING-type'
    elif equalsCount == 1:  # STRESS
        if '=' in line:
            (key, value) = line.strip().split('=')
        else:
            msg = 'expected item of form "name = value"   line=|%r|' % (
                line.strip())
            raise RuntimeError(msg)

        key = key.strip()
        #print value
        value = value.strip()
        #print("key=|%s| value=|%s|" %(key, value))
        paramType = 'STRESS-type'

        #if 'SET' in key:
            #(i,key,value,options,paramType) = parseSetType(i,line,lines,key,value)
        if '(' in key:  # comma may be in line - STRESS-type
            #paramType = 'STRESS-type'
            sline = key.strip(')').split('(')
            key = sline[0]
            options = sline[1].split(',')

            # handle TEMPERATURE(INITIAL) and TEMPERATURE(LOAD) cards
            if key == 'TEMPERATURE' or key == 'TEMP':
                key = 'TEMPERATURE(%s)' % (options[0])
                options = []
            #print "key=|%s| options=%s" %(key,options)

        elif ' ' in key and ',' in value:  # SET-type (with spaces, so SET 1 = 1, not SET = ALL)
            (i, key, value, options, paramType) = parseSetType(
                i, line, lines, key, value)
        elif ',' in value:  # STRESS-type; special TITLE = stuffA,stuffB
            #print 'A ??? line = ',line
            #raise RuntimeError(line)
            pass
        else:  # STRESS-type; TITLE = stuff
            #print 'B ??? line = ',line
            if ' ' in value:
                sline = value.split(' ')
                print "sline = ", sline
                value = parseSetSline(sline)
            else:
                value = interpretValue(value)
        ###
    ### = in line
    elif lineUpper.startswith('BEGIN'):  # begin bulk
        try:
            (key, value) = lineUpper.split(' ')
        except:
            msg = 'excepted "BEGIN BULK" found=|%r|' % (line)
            raise RuntimeError(msg)
        paramType = 'BEGIN_BULK-type'
    elif 'PARAM' in lineUpper:  # param
        sline = line.split(',')
        if len(sline) != 3:
            raise SyntaxError("Param Parse: trying to parse |%s|..." % (line))
        (key, value, options) = sline
        ###
        paramType = 'CSV-type'
    elif ' ' not in line:
        key = line.strip()
        value = line.strip()
        options = None
        paramType = 'KEY-type'
    else:
        msg = 'generic catch all...line=|%r|' % (line)
        #print 'C ??? line = ',line
        #raise RuntimeError(msg)
        key = ''
        value = line
        options = None
        paramType = 'KEY-type'
    ###
    i += 1
    #print "done with ",key
    return (i, key, value, options, paramType)

if __name__ == '__main__':
    #print _parseEntry(['SET =ALL'])
    #print _parseEntry(['SET 77=5'])
    print _parseEntry(['SET 88=5, 6, 7, 8, 9, 10 THRU 55 EXCEPT 15, 16, 77, 78, 79, 100 THRU 300'])
    #print _parseEntry(['SET 99=1 THRU 110'])
    #print _parseEntry(['SET 2001=M1,M2'])
    #print _parseEntry(['SET 101=1.0, 2.0, 3.0'])
    #deck = CaseControlDeck(lines)

    #print _parseEntry(['SET 105=1.009, 10.2, 13.4, 14.0, 15.0'])
    #print _parseEntry(['SET 1001=101/T1, 501/T3, 991/R3'])
