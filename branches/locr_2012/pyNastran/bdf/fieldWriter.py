# pylint: disable=C0103,R0902,R0904,R0914,C0301
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import sys
from numpy import allclose, isinf


def is_same(value1, value2):
    """
    checks to see if 2 values are the same
    @note this method is used by almost every card when printing
    """
    if isinstance(value1, unicode) or value1 is None:
        return True if value1 == value2 else False
    return True if (value1 == value2 or type(value1) == type(value2) and
                    not isinf(value1) and allclose(value1, value2)) else False

def set_blank_if_default(value, default):
    """
    Used when setting the output data of a card to clear default values
    @param value
      the field value the may be set to None (blank) if value=default
    @param default
      the default value for the field
    @note
      this method is used by almost every card when printing
    """
    return None if is_same(value, default) else value


def set_default_if_blank(value, default):
    """
    used when initializing a card and the default value isnt set
    used on PBARL
    """
    return default if value is None or value == '' else value


def print_scientific_8(value):
    """
    Prints a value in 8-character scientific notation.
    This is a sub-method and shouldnt typically be called
    @see print_float for a better method
    """
    #print "scientific...%s" %(value)
    pythonValue = '%8.11e' % (value)
    #print "pythonValue = ",pythonValue
    (svalue, sExponent) = pythonValue.strip().split('e')
    exponent = int(sExponent)  # removes 0s
    #print "svalue=%s exponent=%s" %(svalue,exponent)

    sign = '-' if abs(value) < 0.01 else '+'
    
    sExp2 = str(exponent).strip('-+')  # the exponent will be added later...
    #print "sExp2 = ",sExp2

    value2 = float(svalue)
    #lenSValue = len(svalue)

    lenSExp = len(sExp2) + 1  # the plus 1 is for the sign
    leftover = 8 - lenSExp

    #svalue2 = svalue.strip('0')

    if value < 0:
        #print "sExp2 = ",sExp2
        Format = "%%1.%sf" % (leftover - 3)
    else:
        #print "greater..."
        Format = "%%1.%sf" % (leftover - 2)

    #print("Format = ",Format)
    svalue3 = Format % (value2)
    #print("svalue3 = ",svalue3)
    svalue4 = svalue3.strip('0')
    field = "%8s" % (svalue4 + sign + sExp2)
    #print("fieldA = ", field)

    #print("Format=%s svalue4=%s sExp2=%s" %(Format,svalue4,sExp2))
    #field = "%8s" %(svalue4 + sign +sExp2)
    #print("fieldB = ",field)

    #if '+' in field and '-' in field:
    #print("scientific...value=%s field=%s" %(value, field))
    return field


def print_float_8(value, tol=0.):
    """
    Prints a float in nastran 8-character width syntax using the
    highest precision possbile.
    @todo bad for small values...positive or negative...
    @warning hasnt really be tested for tolerancing
    """
    #value = round(value, 4)
    #print "float...%s" % value
    if abs(value) <= tol:  # tol=1e-8
        #print "below tol %s" %(value)
        field = "%8s" % ('0.')
    else:
        if value > 0.:  # positive, not perfect...
            #print("positive")

            if value < 5e-8:  ## @todo does this work properly with tol
                #print("scientific")
                field = print_scientific_8(value)
                return field
            elif value < 0.001:
                #print "A"
                #print value
                field = print_scientific_8(value)
                field2 = "%8.7f" % (value)  # small value
                field2 = field2.strip('0 ')

                #if 'e' not in field:
                field1 = field.replace('-', 'e-')

                #print("value=|%s| field1=|%s| field2=|%s|" %(value, field,
                #                                             field2))
                #print("same - ", float(field1)==float(field2))
                if field2 == '.':
                    return print_scientific_8(value)
                if len(field2) <= 8 and float(field1) == float(field2):
                    field = field2
                    #print("*field = ",field)
                    field = field.strip(' 0')

                    #print "AA"
                    #print "field  = ",field
                    #print "field1 = ",field1
                    #print "field2 = ",field2
                    #print ""
                ###
            elif value < 0.1:
                #print "B*"
                field = "%8.7f" % (value)
                #field = field.strip('0 ')
                #print field
                #field = field[1:]
            elif value < 1.:
                field = "%8.7f" % (value)
            elif value < 10.:
                field = "%8.6f" % (value)
            elif value < 100.:
                field = "%8.5f" % (value)
            elif value < 1000.:
                field = "%8.4f" % (value)
            elif value < 10000.:
                field = "%8.3f" % (value)
            elif value < 100000.:
                field = "%8.2f" % (value)
            elif value < 1000000.:
                field = "%8.1f" % (value)
            else:  # big value
                #print "big"
                field = "%8.1f" % (value)
                if field.index('.') < 8:
                    field = '%8.1f' % (round(value))
                    field = field[0:8]
                    #field = '%7s.' %(int(field))
                    assert '.' != field[0], field
                else:
                    field = print_scientific_8(value)
                ###
                return field
            ###
        ###
        else:
            #print "negative"
            if value > -5e-7:  ## @todo does this work properly with tol
                #print "really small"
                field = print_scientific_8(value)
                return field
            elif value > -0.01:  # -0.001
                #print "tiny"
                field = print_scientific_8(value)
                field2 = "%8.6f" % (value)  # small value
                field2 = field2.strip('0 ')

                # get rid of the first minus sign, add it on afterwards
                field1 = '-' + field.strip(' 0-').replace('-', 'e-')

                #print "value=%s field=%s field1=%s field2=%s" %(value,
                #                                   field[1:], field1,field2)
                #print "same - ",float(field1)==float(field2)
                if len(field2) <= 8 and float(field1) == float(field2):
                    field = field2.rstrip(' 0')
                    field = field.replace('-0.', '-.')

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
            elif value > -0.1:
                #print "B"
                # -0.01 >x>-0.1...should be 5 (maybe scientific...)
                field = "%8.6f" % (value)
                field = field.replace('-0.', '-.')
            elif value > -1.:
                #print "C"
                # -0.1  >x>-1.....should be 6, but the baseline 0 is kept...
                field = "%8.6f" % (value)
                field = field.replace('-0.', '-.')
            elif value > -10.:
                field = "%8.5f" % (value)   # -1    >x>-10
            elif value > -100.:
                field = "%8.4f" % (value)   # -10   >x>-100
            elif value > -1000.:
                field = "%8.3f" % (value)   # -100  >x>-1000
            elif value > -10000.:
                field = "%8.2f" % (value)   # -1000 >x>-10000
            elif value > -100000.:
                field = "%8.1f" % (value)   # -10000>x>-100000
            else:
                field = "%8.1f" % (value)
                if field.index('.') < 8:
                    field = '%7s.' % (int(round(value, 0)))
                    assert '.' != field[0], field
                else:
                    field = print_scientific_8(value)
                ###
                return field
            ###
        ###
        field = field.strip(' 0')
        field = '%8s' % (field)
    ###
    #print len(field)
    #print "value=|%s| field=|%s|\n" %(value, field)
    assert len(field) == 8, ('value=|%s| field=|%s| is not 8 characters '
                             'long, its %s' % (value, field, len(field)))
    return field


def print_field(value, tol=0.):
    """
    prints a single 8-character width field
    @param value the value to print
    @param tol the abs(tol) to consider value=0 (default=0.)
    @retval field an 8-character (tested) string
    """
    if isinstance(value, int):
        field = "%8s" % (value)
    elif isinstance(value, float):
        field = print_float_8(value)
    elif value is None:
        field = "        "
    else:
        field = "%8s" % (value)
    if len(field) != 8:
        msg = 'field=|%s| is not 8 characters long...rawValue=|%s|' % (field,
                                                                       value)
        raise RuntimeError(msg)
    return field

#def printCard(fields,size=8, tol=0.):
    #"""
    #prints a nastran-style card with 8 or 16-character width fields
    #@param fields all the fields in the BDF card (no blanks)
    #@param tol the abs(tol) to consider value=0 (default=0.)
    #@param size the width of a field (size=8 or 16)
    #@warning 8 or 16 is required, but 16 is not checked for
    #"""
    #if size==8:
    #    return self.printCard_8(fields)
    #else:
    #    return self.printCard_16(fields)
    ###


def printCard(fields, tol=0.):
    """
    Prints a nastran-style card with 8-character width fields.
    @param fields all the fields in the BDF card (no blanks)
    @param tol the abs(tol) to consider value=0 (default=0.)
    @note A small field format follows the  8-8-8-8-8-8-8-8 = 80
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
            out += print_field(field, tol=tol)
            #print "|%r|" %(printField(field))
        except:
            print("bad fields = %s" % (fields))
            raise
        if i % 8 == 0:  # allow 1+8 fields per line
            #print "-------------------------"
            #print "out = ***\n%s***" %(out)
            #print "fields = ",fields[:i+1]
            out = out.rstrip(' ')
            #print "out[-1] = |%r|" %(out[-1])
            if out[-1] == '\n':  # empty line
                out += '+'
            out += '\n        '
    #print "out = ",out
    out = out.rstrip(' \n+') + '\n'  # removes blank lines at the end of cards
    return out


def print_int_card(fields, tol=0.):
    """
    All fields (other than the first field must be integers.
    This is used to speed up SET cards.
    Prints a nastran-style card with 8-character width fields.
    @warning Blanks are not allowed!
    """
    try:
        out = '%-8s' % (fields[0])
    except:
        print("ERROR!  fields=%s" % (fields))
        sys.stdout.flush()
        raise

    for i in xrange(1, len(fields)):
        field = fields[i]
        try:
            out += "%8i" % (field)  # balks if you have None or string fields
        except:
            print("bad fields = %s" % (fields))
            raise
        if i % 8 == 0:  # allow 1+8 fields per line
            out = out.rstrip(' ')
            out += '\n        '
    out = out.rstrip(' \n+') + '\n'  # removes blank lines at the end of cards
    return out


if __name__ == '__main__':
    pass
