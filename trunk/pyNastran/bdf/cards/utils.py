from six import string_types

def build_table_lines(fields, nStart=1, nEnd=0):
    """
    Builds a table of the form:

    +--------+-------+--------+-------+--------+-------+-------+-------+
    | DESVAR | DVID1 | DVID2  | DVID3 |DVID4   | DVID5 | DVID6 | DVID7 |
    +--------+-------+--------+-------+--------+-------+-------+-------+
    |        | DVID8 | -etc.- |       |        |
    +--------+-------+--------+-------+--------+
    |        |  UM   | VAL1   | VAL2  | -etc.- |
    +--------+-------+--------+-------+--------+

    and then pads the rest of the fields with None's

    :param fields: the fields to enter, including DESVAR
    :type fields:  list of values
    :param nStart: the number of blank fields at the start of the
                   line (default=1)
    :param nStart: int
    :param nEnd:   the number of blank fields at the end of the
                   line (default=0)
    :param nEnd:   int

    .. note:: will be used for DVPREL2, RBE1, RBE3
    .. warning:: only works for small field format???
    """
    fieldsOut = []
    n = 8 - nStart - nEnd

    # pack all the fields into a list.  Only consider an entry as isolated
    for (i, field) in enumerate(fields):
        fieldsOut.append(field)
        if i > 0 and i % n == 0:  # beginning of line
            #print("i=%s" %(i))
            #pad = [None]*(i+j)
            #fieldsOut += pad
            fieldsOut += [None] * (nStart + nEnd)

    # make sure they're aren't any trailing None's (from a new line)
    fieldsOut = wipe_empty_fields(fieldsOut)

    # push the next key (aka next fields[0]) onto the next line
    nSpaces = 8 - (len(fieldsOut)) % 8  # puts UM onto next line
    if nSpaces < 8:
        fieldsOut += [None] * nSpaces
    return fieldsOut


def wipe_empty_fields(card):
    """
    Removes an trailing Nones from the card.
    Also converts empty strings to None.

    :param card:        the fields on the card as a list
    :returns shortCard: the card with no trailing blank fields
    """
    cardB = []
    for field in card:
        if isinstance(field, string_types):
            field = field.strip()
            if field == '':
                field = None
        cardB.append(field)

    i = 0
    iMax = 0
    while i < len(card):
        if cardB[i] is not None:
            iMax = i
        i += 1
    return cardB[:iMax + 1]
