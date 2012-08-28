def printScientific16(value):
    """
    Prints a value in 16-character scientific notation.
    This is a sub-method and shouldnt typically be called
    @warning not tested...
    """
    #print "scientific...%s" %(value)
    pythonValue = '%16.14e' % (value)  # -1.e-2
    #print "pythonValue = ",pythonValue
    svalue, sExponent = pythonValue.strip().split('e')
    exponent = int(sExponent)  # removes 0s
    #print "svalue=%s exponent=%s" %(svalue,exponent)

    if abs(value) < 0.01:
        sign = '-'
    else:
        sign = '+'
    sExp2 = str(exponent).strip('-+')  # the exponent will be added later...
    #print "sExp2 = ",sExp2

    value2 = float(svalue)
    #lenSValue = len(svalue)

    lenSExp = len(sExp2) + 1  # the plus 1 is for the sign
    leftover = 16 - lenSExp

    svalue2 = svalue.strip('0')

    if value < 0:
        #print "sExp2 = ",sExp2
        format = "%%1.%sf" % (leftover - 3)
    else:
        #print "greater..."
        format = "%%1.%sf" % (leftover - 2)

    #print "format = ",format
    svalue3 = format % (value2)
    #print "svalue3 = ",svalue3
    svalue4 = svalue3.strip('0')
    field = "%16s" % (svalue4 + sign + sExp2)
    return field


def printFloat16(value, tol=1e-8):
    """
    Prints a float in nastran 16-character width syntax
    using the highest precision possbile.
    @see printFloat8
    @warning completely unimplemented & untested
    """
    strVal = printFloat8(value, tol)
    return '%16s' % (strVal)

    #value = round(value,4)
    #print "float...%s" %value
    if abs(value) <= tol:  # tol=1e-8
        #print "below tol %s" %(value)
        field = "%16s" % ('0.')
    else:
        if value > 0.:  # positive, not perfect...
            #print "positive"

            if value < 5e-16:  ## @todo does this work properly with tol
                #print "scientific"
                field = printScientific16(value)
                return field
            elif value < 0.001:
                #print "A"
                if 1:
                    #print value
                    field = printScientific16(value)
                    field2 = "%16.15f" % (value)  # small value
                    field2 = field2.strip('0 ')

                    #if 'e' not in field:
                    field1 = field.replace('-', 'e-')

                    #print "value=|%s| field1=|%s| field2=|%s|" %(value,field,field2)
                    #print "same - ",float(field1)==float(field2)
                    if field2 == '.':
                        return "%16s" % (field)
                    if len(field2) <= 16 and float(field1) == float(field2):
                        field = field2
                        #print "*field = ",field
                        field = field.strip(' 0')

                        #print "AA"
                        #print "field  = ",field
                        #print "field1 = ",field1
                        #print "field2 = ",field2
                        #print ""
                    ###
                ###
                if 0:
                    field = "%16.7f" % (value)
                    #print "field = ",field
                    field = field.strip('0')
                    if len(field) < 16:
                        assert '.' == field[0], field
                    else:
                        field = printScientific16(value)
                        return field
                    ###
                    #print "field = ",field
                ###
            elif value < 0.1:
                #print "B*"
                field = "%16.15f" % (value)
                #field = field.strip('0 ')
                #print field
                #field = field[1:]
            elif value < 1.:
                field = "%16.15f" % (value)
            elif value < 10.:
                field = "%16.14f" % (value)
            elif value < 100.:
                field = "%16.13f" % (value)
            elif value < 1000.:
                field = "%16.12f" % (value)
            elif value < 10000.:
                field = "%16.11f" % (value)
            elif value < 100000.:
                field = "%16.10f" % (value)
            elif value < 1000000.:
                field = "%16.9f" % (value)
            elif value < 10000000.:
                field = "%16.8f" % (value)
            elif value < 100000000.:
                field = "%16.7f" % (value)
            elif value < 1000000000.:
                field = "%16.6f" % (value)
            elif value < 10000000000.:
                field = "%16.5f" % (value)
            elif value < 100000000000.:
                field = "%16.4f" % (value)
            elif value < 1000000000000.:
                field = "%16.3f" % (value)
            elif value < 10000000000000.:
                field = "%16.2f" % (value)
            elif value < 100000000000000.:
                field = "%16.1f" % (value)
            else:  # big value  # 123456789012345.
                #print "big"   # 100000000000000.
                field = "%16.1f" % (round(value))
                if field.index('.') < 16:
                    field = field[0:16]  # drop off the .1f
                    assert '.' != field[0], field
                else:
                    field = printScientific16(value)
                ###
                return field
            ###
        ###
        else:
            #print "negative"
            if value > -5e-15:  ## @todo does this work properly with tol
                #print "really small"
                field = printScientific16(value)
                return field
            elif value > -0.01:  # -0.001
                #print "tiny"
                field = printScientific16(value)
                field2 = "%16.14f" % (value)  # small value
                field2 = field2.strip('0 ')

                #if 'e' not in field:
                field1 = '-' + field.strip(' 0-').replace('-', 'e-')  # get rid of the first minus sign, add it on afterwards

                #print "value=%s field=%s field1=%s field2=%s" %(value,field[1:],field1,field2)
                #print "same - ",float(field1)==float(field2)
                if len(field2) <= 16 and float(field1) == float(field2):
                    field = field2
                    #print "*field = ",field
                    field = field.strip(' 0')

                    #print "AA"
                    #print "field  = ",field
                    #print "field1 = ",field1
                    #print "field2 = ",field2
                #print ""
                ###
            #elif value>-0.01:
            #    #print "A"
            #    field = "%16.16f" %(value)   # -0.001>x>-0.01..should be 4
            #    field = '-'+field[2:]
            elif value > -0.1:
                #print "B"
                field = "%16.14f" % (value)   # -0.01 >x>-0.1...should be 5 (maybe scientific...)
                field = field.replace('-0.', '-.')
            elif value > -1.:
                #print "C"
                field = "%16.14f" % (value)   # -0.1  >x>-1.....should be 6, but the baseline 0 is kept...
                field = field.replace('-0.', '-.')
            elif value > -10.:
                field = "%16.14f" % (value)   # -1    >x>-10
            elif value > -100.:
                field = "%16.13f" % (value)   # -10   >x>-100
            elif value > -1000.:
                field = "%16.12f" % (value)   # -100  >x>-1000
            elif value > -10000.:
                field = "%16.11f" % (value)   # -1000 >x>-10000
            elif value > -100000.:
                field = "%16.10f" % (value)   # -10000>x>-100000
            elif value > -1000000.:
                field = "%16.9f" % (value)
            elif value > -10000000.:
                field = "%16.8f" % (value)
            elif value > -100000000.:
                field = "%16.7f" % (value)
            elif value > -1000000000.:
                field = "%16.6f" % (value)
            elif value > -10000000000.:
                field = "%16.5f" % (value)
            elif value > -100000000000.:
                field = "%16.4f" % (value)
            elif value > -1000000000000.:
                field = "%16.3f" % (value)
            elif value > -10000000000000.:
                field = "%16.2f" % (value)
            elif value > -100000000000000.:
                field = "%16.1f" % (value)
            else:
                field = "%16.1f" % (value)
                if field.index('.') < 16:
                    field = field[0:16]
                    assert '.' != field[0], field
                else:
                    field = printScientific16(value)
                ###
                return field
            ###
        ###
        field = field.strip(' 0')
        field = '%16s' % (field)
    ###
    #print len(field)
    #print "value=|%s| field=|%s|\n" %(value,field)
    assert len(field) == 16, 'value=|%s| field=|%s| is not 16 characters long, its %s' % (value, field, len(field))
    return field


def printField16(value, tol=0.):
    """
    prints a single 16-character width field
    @param value the value to print
    @param tol the abs(tol) to consider value=0 (default=0.)
    @retval field an 16-character (tested) string
    """
    if isinstance(value, int):
        field = "%16s" % (value)
    elif isinstance(value, float):
        #print "float..."
        field = printFloat16(value)
    elif isinstance(value, NoneType):
        field = "                "
    else:
        field = "%16s" % (value)
    assert len(field) == 16, 'field=|%s| is not 16 characters long...rawValue=|%s|' % (field, value)
    return field


def printCard_16(fields, tol=0.):
    """
    Prints a nastran-style card with 16-character width fields.
    @param fields all the fields in the BDF card (no blanks)
    @param tol the abs(tol) to consider value=0 (default=0.)
    @note A large field format follows the  8-16-16-16-16-8 = 80
     format where the first 8 is the card name or blank (continuation).
     The last 8-character field indicates an optional continuation,
     but because it's a left-justified unneccessary field,
     printCard doesnt use it.
    """
    #print fields
    try:
        out = '%-8s' % (fields[0])
    except:
        print("ERROR!  fields=%s" % (fields))
        sys.stdout.flush()
        raise

    for i in xrange(1, len(fields)):
        field = fields[i]
        try:
            out += printField16(field)
            #print "|%r|" %(printField(field))
        except AssertionError:
            print("bad fields = ", fields)
            raise
        if i % 4 == 0:  # allow 1+4 fields per line
            #print "-------------------------"
            #print "out = ***\n%s***" %(out)
            #print "fields = ",fields[:i+1]
            out = out.rstrip()
            if out[-1] == '\n':  # empty line
                out += '*'
            out += '\n%8s' % ('')
    out = out.rstrip(' \n*') + '\n'  # removes blank lines at the end of cards
    return out
