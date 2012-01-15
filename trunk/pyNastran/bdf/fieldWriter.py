from types import NoneType
from numpy import allclose,isinf

def isSame(value1,value2):
    """
    checks to see if 2 values are the same
    @note this method is used by almost every card when printing
    """
    #print "value=%s default=%s" %(value1,value2)
    if isinstance(value1,str) or isinstance(value1,NoneType):
        if value1==value2:
            return True
        return False
    elif type(value1)==type(value2) and not isinf(value1) and allclose(value1,value2):
        #print "value=%s default=%s same=%s" %(value,default,True)
        return True
    #print "value=%s default=%s same=%s" %(value,default,False)
    return False

def setBlankIfDefault(value,default):
    """
    used when setting the output data of a card to clear default values
    @note this method is used by almost every card when printing
    """
    if isSame(value,default):
        return None
    return value

def setDefaultIfBlank(value,default):
    """
    used when initializing a card and the default value isnt set
    used on PBARL
    """
    if value is None or value=='':
        return default
    return value

def printScientific8(value):
    """
    Prints a value in 8-character scientific notation.
    This is a sub-method and shouldnt typically be called
    @see printFloat for a better method
    """
    #print "scientific...%s" %(value)
    pythonValue = '%8.6e' %(value)
    #print "pythonValue = ",pythonValue
    svalue,sExponent = pythonValue.strip().split('e')
    exponent = int(sExponent) # removes 0s
    #print "svalue=%s exponent=%s" %(svalue,exponent)
    
    if abs(value)<0.01:
        sign = '-'
    else:
        sign = '+'
    sExp2 = str(exponent).strip('-+')  # the exponent will be added later...

    value2 = float(svalue)
    #lenSValue = len(svalue)
    
    lenSExp   = len(sExp2)+1 # the plus 1 is for the sign
    leftover = 8-lenSExp
    
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
    field = "%8s" %(svalue4 + sign +sExp2)
    #print "fieldA = ",field

    #print "format=%s svalue4=%s sExp2=%s" %(format,svalue4,sExp2)
    #field = "%8s" %(svalue4 + sign +sExp2)
    #print "fieldB = ",field
    
    #if '+' in field and '-' in field:
    #print "scientific...value=%s field=%s" %(value, field)
    return field

def printScientific16(value):
    """
    Prints a value in 16-character scientific notation.
    This is a sub-method and shouldnt typically be called
    @warning not tested...
    """
    pythonValue = '%16.14e' %(value)
    svalue,sExponent = pythonValue.strip().split('e')
    exponent = int(sExponent) # removes 0s
    
    if abs(value)<0.01:
        sign = '-'
    else:
        sign = '+'
    sExp2 = str(exponent).strip('-+')  # the exponent will be added later...

    value2 = float(svalue)
    
    lenSExp   = len(sExp2)+1 # the plus 1 is for the sign
    leftover = 16-lenSExp
    
    svalue2 = svalue.strip('0')
    
    if value<0:
        format = "%%1.%sf" %(leftover-3)
    else:
        format = "%%1.%sf" %(leftover-2)

    svalue3 = format %(value2)
    svalue4 = svalue3.strip('0')
    field = "%16s" %(svalue4 + sign +sExp2)
    return field

def printFloat8(value,tol=1e-8):
    """
    Prints a float in nastran 8-character width syntax.
    using the highest precision possbile.
    @todo bad for small values...positive or negative...
    @warning hasnt really be tested for tolerancing
    """
    #value = round(value,4)
    #print "float...%s" %value
    if abs(value)<tol:  # tol=1e-8
        field = "%8s" %('0.')
    else:
        if value>0.:  # positive, not perfect...
            #print "positive"

            if value<5e-8:  ## @todo does this work properly with tol
                field = printScientific8(value)
            elif value<0.0001:
                #print "A"
                if 1:
                    #print value
                    field = printScientific(value)
                    field2 = "%8.7f" %(value) # small value
                    field2 = field2.strip('0 ')

                    #if 'e' not in field:
                    field1 = field.replace('-','e-')

                    #print "value=|%s| field1=|%s| field2=|%s|" %(value,field1,field2)
                    #print "same - ",float(field1)==float(field2)
                    if field2=='.':
                        return "%8s" %('0.')
                    if len(field2)<=8 and float(field1)==float(field2):
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
                    field = "%8.7f" %(value)
                    #print "field = ",field
                    field = field.strip('0')
                    if len(field)<8:
                        assert '.' == field[0],field
                    else:
                        field = printScientific8(value)
                    ###
                    print "field = ",field
                ###
            elif value<0.1:
                #print "B*"
                field = "%8.7f" %(value)
                #field = field.strip('0 ')
                #print field
                #field = field[1:]
            elif value<1.:       field = "%8.7f" %(value)
            elif value<10:       field = "%8.6f" %(value)
            elif value<100.:     field = "%8.5f" %(value)
            elif value<1000.:    field = "%8.4f" %(value)
            elif value<10000.:   field = "%8.3f" %(value)
            elif value<100000.:  field = "%8.2f" %(value)
            elif value<1000000.: field = "%8.1f" %(value)
            else: # big value
                #print "big"
                field = "%8.1f" %(value)
                if field.index('.')<8:
                    field = field[0:8]
                    assert '.' != field[0],field
                else:
                    field = printScientific8(value)
                ###
            ###
        ###
        else:
            #print "negative"
            if value>-5e-7:  ## @todo does this work properly with tol
                #print "really small"
                field = printScientific8(value)
            elif value>-0.01:  # -0.001
                #print "tiny"
                field = printScientific8(value)
                field2 = "%8.6f" %(value) # small value
                field2 = field2.strip('0 ')

                #if 'e' not in field:
                field1 = '-'+field.strip(' 0-').replace('-','e-') # get rid of the first minus sign, add it on afterwards

                #print "value=%s field=%s field1=%s field2=%s" %(value,field[1:],field1,field2)
                #print "same - ",float(field1)==float(field2)
                if len(field2)<=8 and float(field1)==float(field2):
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
            #    field = "%8.8f" %(value)   # -0.001>x>-0.01..should be 4
            #    field = '-'+field[2:]
            elif value>-0.1:
                #print "B"
                field = "%8.6f" %(value)   # -0.01 >x>-0.1...should be 5 (maybe scientific...)
                field = field.replace('-0.','-.')
            elif value>-1.:
                #print "C"
                field = "%8.6f" %(value)   # -0.1  >x>-1.....should be 6, but the baseline 0 is kept...
                field = field.replace('-0.','-.')
            elif value>-10.:   field = "%8.5f" %(value)   # -1    >x>-10
            elif value>-100:   field = "%8.4f" %(value)   # -10   >x>-100
            elif value>-1000:  field = "%8.3f" %(value)   # -100  >x>-1000
            elif value>-10000: field = "%8.2f" %(value)   # -1000 >x>-10000
            elif value>-100000:field = "%8.1f" %(value)   # -10000>x>-100000
            else:
                field = printScientific8(value)
                #field2 = 
            ###
        ###
        field = field.strip(' 0')
        field = '%8s' %(field)
    ###
    #print len(field)
    #print "value=|%s| field=|%s|\n" %(value,field)
    assert len(field)==8,'value=|%s| field=|%s| is not 8 characters long, its %s' %(value,field,len(field))
    return field

def printFloat16(value,tol=1e-8):
    """
    @todo preliminary print method uses 16-character width...
    @see printFloat8
    """
    strVal = printFloat8(value,tol)
    return '%16s' %(strVal)

def printField(value,tol=1e-8):
    """
    prints a single 8-character width field
    @param value the value to print
    @param tol the abs(tol) to consider value=0 (default=1e-8)
    @retval field an 8-character (tested) string
    """
    if isinstance(value,int):
        field = "%8s" %(value)
    elif isinstance(value,float):
        #print "float..."
        field = printFloat8(value)
    elif isinstance(value,NoneType):
        field = "        "
    else:
        field = "%8s" %(value)
    assert len(field)==8,'field=|%s| is not 8 characters long...rawValue=|%s|' %(field,value)
    return field

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

#def printCard(fields,size=8):
    #"""
    #prints a nastran-style card with 8 or 16-character width fields
    #@warning 8 or 16 is required, but 16 is not checked for
    #"""
    #if size==8:
    #    return self.printCard_8(fields)
    #else:
    #    return self.printCard_16(fields)
    ###

def printCard(fields):
    """
    Prints a nastran-style card with 8-character width fields.
    
    @note A small field format follows the  8-8-8-8-8-8-8-8 = 80
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
            out += printField(field)
        except AssertionError:
            print "bad fields = ",fields
            raise
        if i%8==0: # allow 1+8 fields per line
            out = out.rstrip()
            out += '\n%8s' %('')
        ###
    ###
    out = out.rstrip()+'\n'
    return out

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
        except AssertionError:
            print "bad fields = ",fields
            raise
        if i%4==0: # allow 1+4 fields per line
            out = out.rstrip()
            out += '\n%8s' %('')
        ###
    ###
    out = out.rstrip()+'\n'
    return out

def displayCard(fields):
    """
    Prints a cards fields in an easy to read/debug format.
    @param fields the fields to "print" to standard out.
    @note this method has been made obsolete by the maturation
    of card printing methods, but it's still useful.
    """
    for i,field in enumerate(fields):
        print "field[%s] = %s" %(i,field)
    ###

if __name__=='__main__':
    #printField(-0.021004)
    #print printField(-4.21704e-6)
    #print printField(4.21704e-6)
    #print printField(8.17272e-6)
    #print printField(10300000.0)
    #print printField(-10300000.0)
    #printField(-0.021004)
    
    field = printField(-.723476);   assert '-.723476' == field,'|%s|' %(field)
    field = printField(125000. );   assert ' 125000.' == field,'|%s|' %(field)
    field = printField(12500000.);  assert '  1.25+7' == field,'|%s|' %(field)
    field = printField(47.77267);   assert '47.77267' == field,'|%s|' %(field)
    field = printField(.001);       assert '    .001' == field,'|%s|' %(field)
    field = printField(.0000001);   assert '.0000001' == field,'|%s|' %(field)
    field = printField(-5.007e-3);  assert '-5.007-3' == field,'|%s|' %(field)
    field = printField(-0.0748662); assert '-.074866' == field,'|%s|' %(field)
    #print printField(12500000. )
    #print printField(47.77267)
    #print printField(.0000001)
    #print printField(-5.007e-3)
    

