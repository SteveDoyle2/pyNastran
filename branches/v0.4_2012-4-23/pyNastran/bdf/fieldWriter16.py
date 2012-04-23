def printScientific16(value):
    """
    Prints a value in 16-character scientific notation.
    This is a sub-method and shouldnt typically be called
    @warning not tested...
    """
    pythonValue = '%16.14e' %(value)
    #print "pythonValue = ",pythonValue
    svalue,sExponent = pythonValue.strip().split('e')
    exponent = int(sExponent) # removes 0s
    #print "svalue=%s exponent=%s" %(svalue,exponent)
    
    if abs(value)<0.01:
        sign = '-'
    else:
        sign = '+'
    sExp2 = str(exponent).strip('-+')  # the exponent will be added later...
    #print "sExp2 = ",sExp2

    value2 = float(svalue)
    #lenSValue = len(svalue)
    
    lenSExp   = len(sExp2)+1 # the plus 1 is for the sign
    leftover = 16-lenSExp
    
    svalue2 = svalue.strip('0')
    
    if value<0:
        #print "sExp2 = ",sExp2
        format = "%%1.%sf" %(leftover-3)
    else:
        #print "greater..."
        format = "%%1.%sf" %(leftover-2)

    #print "format = ",format
    svalue3 = format %(value2)
    #print "svalue3 = ",svalue3
    svalue4 = svalue3.strip('0')
    field = "%16s" %(svalue4 + sign +sExp2)
    return field

def printFloat16(value,tol=1e-8):
    """
    @todo preliminary print method uses 16-character width...
    @see printFloat8
    """
    strVal = printFloat8(value,tol)
    return '%16s' %(strVal)

def printField16(value,tol=1e-8):
    """
    prints a single 16-character width field
    @param value the value to print
    @param tol the abs(tol) to consider value=0 (default=1e-8)
    @retval field an 16-character (tested) string
    """
    if isinstance(value,int):
        field = "%16s" %(value)
    elif isinstance(value,float):
        #print "float..."
        field = printFloat16(value)
    elif isinstance(value,NoneType):
        field = "                "
    else:
        field = "%16s" %(value)
    assert len(field)==16,'field=|%s| is not 16 characters long...rawValue=|%s|' %(field,value)
    return field

def printCard_16(fields):
    """
    Prints a nastran-style card with 16-character width fields.
    
    @note A large field format follows the  8-16-16-16-16-8 = 80.
    format where the first 8 is the card name or blank (continuation).
    The last 8-character field indicates an optional continuation,
    but because it's a left-justified unneccessary field,
    printCard doesnt use it.
    """
    try:
        out = '%-8s' %(fields[0])
    except:
        print "ERROR!  fields=%s" %(fields)
        raise
    
    for i in range(1,len(fields)):
        field = fields[i]
        try:
            out += printField16(field)
            #print "|%r|" %(printField(field))
        except AssertionError:
            print "bad fields = ",fields
            raise
        if i%4==0: # allow 1+4 fields per line
            #print "-------------------------"
            #print "out = ***\n%s***" %(out)
            #print "fields = ",fields[:i+1]
            out = out.rstrip()
            out += '\n%8s' %('')
        ###
    ###
    out = out.rstrip()+'\n'
    return out

