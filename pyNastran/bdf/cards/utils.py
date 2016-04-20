from __future__ import unicode_literals, print_function
from six import string_types

def build_table_lines(fields, nstart=1, nend=0):
    """
    Builds a table of the form:

    +--------+-------+--------+-------+--------+-------+-------+-------+
    | DESVAR | DVID1 | DVID2  | DVID3 |DVID4   | DVID5 | DVID6 | DVID7 |
    +--------+-------+--------+-------+--------+-------+-------+-------+
    |        | DVID8 | -etc.- |       |        |       |       |       |
    +--------+-------+--------+-------+--------+-------+-------+-------+
    |        |  UM   | VAL1   | VAL2  | -etc.- |       |       |       |
    +--------+-------+--------+-------+--------+-------+-------+-------+

    and pads the rest of the fields with None's (e.g. at the end of the
    DVID8 line).

    Parameters
    ----------
    fields: List[int/float/str]
        the fields to enter, including DESVAR
    nStart: int; default=1
        the number of blank fields at the start of the line
    nEnd : int; default=0)
        the number of blank fields at the end of the line

    Returns
    -------
    fields_out : List[str/None]
        the values ready for card printing

    .. note:: will be used for DVPREL2, RBE1, RBE3
    .. warning:: only works for small field format???
    """
    fields_out = []
    n = 8 - nstart - nend

    # pack all the fields into a list.  Only consider an entry as isolated
    for (i, field) in enumerate(fields):
        fields_out.append(field)
        if i > 0 and i % n == 0:  # beginning of line
            #pad = [None] * (i + j)
            #fields_out += pad
            fields_out += [None] * (nstart + nend)

    # make sure they're aren't any trailing None's (from a new line)
    fields_out = wipe_empty_fields(fields_out)

    # push the next key (aka next fields[0]) onto the next line
    nspaces = 8 - (len(fields_out)) % 8  # puts UM onto next line
    if nspaces < 8:
        fields_out += [None] * nspaces
    return fields_out


def wipe_empty_fields_off(fields):
    """
    Removes any trailing Nones from the card.
    Also converts empty strings to None.

    .. warning:: doesn't allows floats & ints.
    """
    nfields = len(fields)
    for i, field in enumerate(reversed(fields)):
        if field:
            break
    # print('i =', i)
    # print('fields =', fields[:nfields-i])

    # has the potential issue of not returning embedded Nones
    return [field.strip() if isinstance(field, string_types) else field for field in fields[:nfields-i]]
    # return [field.strip() if field.strip() else None for field in fields[:nfields-i] ] # fails on ints/Nones

def wipe_empty_fields(fields):
    i = len(fields) - 1
    # print('fields1', fields)
    while i > 0:
        field = fields[i]
        #print('i=%s %r; %s' % (i, field, type(field)))
        if isinstance(field, string_types):
            field_strip = field.strip()
            if field_strip:
                break
            else:
                #print('popping %r' % field_strip)
                fields.pop()
                i -= 1
        else:
            zafd

    while i > 0:
        field = fields[i]
        if isinstance(field, string_types):
            field_strip = field.strip()
            if not field_strip:
                field = field_strip
                i -= 1
            else:
                # return fields
                break
        fields[i] = field
    # print('fields2', fields)
    return fields

#def _wipe_empty_fields(card):
    #"""
    #For testing:
        #this method is:   wipe_empty_fields
        #the one above is: wipe_empty_fields_new

    #In general:
        #this method is     _wipe_empty_fields
        #the  one above is: wipe_empty_fields
    #"""
    #typed = wipe_empty_fields_typed(card)
    #untyped = wipe_empty_fields_new(card)  # the above method is wipe_empty_fields_new
    #untyped = wipe_empty_fields_new(card)  # the above method is wipe_empty_fields_new
    #if typed != untyped:
        #msg = 'typed   = %s\n' % typed
        #msg += 'untyped = %s\n' % untyped
        #print(msg)
    #return untyped

def wipe_empty_fields_typed(card):
    """
    Removes any trailing Nones from the card.
    Also converts empty strings to None.
    Allows floats & ints.

    Parameters
    ----------
    card : List[str]
        the fields on the card as a list

    Returns
    -------
    short_card : List[str]
        the card with no trailing blank fields
    """
    short_card = []
    for field in card:
        if isinstance(field, string_types):
            field = field.strip()
            if field == '':
                field = None
        short_card.append(field)

    imax = 0
    for i, field in enumerate(card):
        if short_card[i] is not None:
            imax = i
    out = short_card[:imax + 1]
    return out

wipe_empty_fields = wipe_empty_fields_typed
