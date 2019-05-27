from pyNastran.bdf.bdf_interface.assign_type import interpret_value
#from pyNastran.bdf.caseControlDeck import CaseControlDeck


def parse_set_sline(list_a):
    print("list_a = ", list_a)
    list_b = []
    for spot in list_a:
        spot = spot.strip()
        print("spot = ", spot)
        if spot == '':
            pass
        elif ' ' in spot:
            sline = spot.split(' ')
            print("sline = ", sline)
            if sline[1] == 'THRU':
                if 'BY' in sline:
                    by = sline[4]
                else:
                    by = 1
                #print("BY = %s" % by)
                vals = set()
                start_value = interpret_value(sline[0])
                end_value = interpret_value(sline[2]) + 1
                for i in range(start_value, end_value, by):
                    vals.add(i)
                #print("vals = ", vals)
                if 'EXCEPT' in sline:
                    iexcept = sline.index('EXCEPT')
                    #print("e = ", sline[iexcept])
                    excepted = int(sline[iexcept + 1])
                    vals.remove(excepted)
                vals = list(vals)
                list_b += vals
                print("vals = ", vals)
            else:
                print("sline = ", sline)
        else:
            #print("spot = %s" % (spot))
            if '/' in spot:
                list_b.append(spot)
            else:
                list_b.append(interpret_value(spot))
    return list_b


def parse_set_type(i, line, lines, key, value):
    (key, ID) = key.split()
    key = key + ' ' + ID

    #print('SET-type key=%s ID=%s' %(key, ID))
    sline = value.strip(' ,').split(',')  # float/int values
    fivalues = parse_set_sline(sline)

    #: .. todo:: should be more efficient multiline reader...
    # read more lines....
    if line[-1].strip() == ',':
        i += 1
        #print("rawSETLine = %r" % lines[i])
        while 1:
            if ',' == lines[i].strip()[-1]:
                sline = lines[i][:-1].strip().split(',')
                fivalues += parse_set_sline(sline)
            else:  # last case
                sline = lines[i].strip().split(',')
                fivalues += parse_set_sline(sline)
                #print("fivalues last = i=%s |%r|" % (i, lines[i]))
                i += 1
                break
            i += 1

    #print("len(fivalues) = %s" % len(fivalues))
    value = fivalues

    options = int(ID)  # needed a place to put it...
    param_type = 'SET-type'
    return (i, key, value, options, param_type)


def _parse_entry(lines):
    i = 0
    options = []
    value = None
    key = None
    param_type = None

    line = lines[i]
    #print(line)
    #print("*****lines = ", lines)

    equals_count = 0
    for letter in line:
        if letter == '=':
            equals_count += 1
    line_upper = line.upper()

    if line_upper.startswith('SUBCASE'):
        #print("line = %r" % line)
        line2 = line.replace('=', '')
        sline = line2.split()
        if len(sline) != 2:
            msg = "trying to parse %r..." % line
            raise SyntaxError("Invalid Subcase: %s" % msg)
        (key, isubcase) = sline
        #print("key=%r isubcase=%r" % (key, isubcase))
        value = int(isubcase)
        #self. subcase = int(isubcase)
        param_type = 'SUBCASE-type'
    elif (line_upper.startswith(('LABEL', 'SUBTITLE')) or
          line_upper.startswith('TITLE')):
        eindex = line.index('=')
        key = line[0:eindex].strip()
        value = line[eindex + 1:].strip()
        options = []
        param_type = 'STRING-type'
    elif equals_count == 1:  # STRESS
        if '=' in line:
            (key, value) = line.strip().split('=')
        else:
            msg = 'expected item of form "name = value"   line=|%r|' % (
                line.strip())
            raise RuntimeError(msg)

        key = key.strip()
        #print(value)
        value = value.strip()
        #print("key=%r value=%r" % (key, value))
        param_type = 'STRESS-type'

        #if 'SET' in key:
            #(i, key, value, options, param_type) = parse_set_type(i,line,lines,key,value)
        if '(' in key:  # comma may be in line - STRESS-type
            #param_type = 'STRESS-type'
            sline = key.strip(')').split('(')
            key = sline[0]
            options = sline[1].split(',')

            # handle TEMPERATURE(INITIAL) and TEMPERATURE(LOAD) cards
            if key in ['TEMPERATURE', 'TEMP']:
                key = 'TEMPERATURE(%s)' % (options[0])
                options = []
            #print("key=%r options=%s" % (key, options))

        elif ' ' in key and ',' in value:
            # SET-type (with spaces, so SET 1 = 1, not SET = ALL)
            (i, key, value, options, param_type) = parse_set_type(
                i, line, lines, key, value)
        elif ',' in value:  # STRESS-type; special TITLE = stuffA,stuffB
            #print('A ??? line = ',line)
            #raise RuntimeError(line)
            pass
        else:  # STRESS-type; TITLE = stuff
            #print('B ??? line = ',line)
            if ' ' in value:
                sline = value.split(' ')
                print("sline = ", sline)
                value = parse_set_sline(sline)
            else:
                value = interpret_value(value)

    elif line_upper.startswith('BEGIN'):  # begin bulk
        try:
            (key, value) = line_upper.split(' ')
        except ValueError:
            msg = 'excepted "BEGIN BULK" found=%r' % (line)
            raise RuntimeError(msg)
        param_type = 'BEGIN_BULK-type'
    elif 'PARAM' in line_upper:  # param
        sline = line.split(',')
        if len(sline) != 3:
            raise SyntaxError("Param Parse: trying to parse %r..." % line)
        (key, value, options) = sline
        param_type = 'CSV-type'
    elif ' ' not in line:
        key = line.strip()
        value = line.strip()
        options = None
        param_type = 'KEY-type'
    else:
        msg = 'generic catch all...line=%r' % line
        #print('C ??? line = ', line)
        #raise RuntimeError(msg)
        key = ''
        value = line
        options = None
        param_type = 'KEY-type'
    i += 1
    #print("done with %s" % key)
    return (i, key, value, options, param_type)

if __name__ == '__main__':  # pragma: no cover
    pass
    #print(_parse_entry(['SET =ALL']))
    #print(_parse_entry(['SET 77=5']))
    #print(_parse_entry(['SET 88=5, 6, 7, 8, 9, 10 THRU 55 EXCEPT 15, 16, 77, 78, 79, 100 THRU 300']))
    #print(_parse_entry(['SET 99=1 THRU 110']))
    #print(_parse_entry(['SET 2001=M1,M2']))
    #print(_parse_entry(['SET 101=1.0, 2.0, 3.0']))
    #deck = CaseControlDeck(lines)

    #print(_parse_entry(['SET 105=1.009, 10.2, 13.4, 14.0, 15.0']))
    #print(_parse_entry(['SET 1001=101/T1, 501/T3, 991/R3']))
