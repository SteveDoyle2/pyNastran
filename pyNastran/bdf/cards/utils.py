"""
defines:
 - fields_out = build_table_lines(fields, nstart=1, nend=0)
 - fields_out = wipe_empty_fields(card)

"""
from typing import List, Union


def build_table_lines(fields, nstart=1, nend=0):
    """
    Builds a table of the form:

    +--------+-------+--------+-------+--------+-------+-------+-------+
    | DESVAR | DVID1 | DVID2  | DVID3 |DVID4   | DVID5 | DVID6 | DVID7 |
    +--------+-------+--------+-------+--------+-------+-------+-------+
    |        | DVID8 |  etc.  |       |        |       |       |       |
    +--------+-------+--------+-------+--------+-------+-------+-------+
    |        |  UM   | VAL1   | VAL2  |  etc.  |       |       |       |
    +--------+-------+--------+-------+--------+-------+-------+-------+

    and pads the rest of the fields with None's (e.g. at the end of the
    DVID8 line).

    Parameters
    ----------
    fields: List[int/float/str]
        the fields to enter, including DESVAR
    nstart: int; default=1
        the number of blank fields at the start of the line
    nend : int; default=0)
        the number of blank fields at the end of the line

    Returns
    -------
    fields_out : List[str/None]
        the values ready for card printing

    .. note:: will be used for DVCREL2, DVMREL2, DVPREL2, RBE1, RBE3, DRESP2, DRESP3
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


def wipe_empty_fields(card: List[Union[str, int, float, None]]) -> List[Union[str, int, float, None]]:
    """
    Removes any trailing Nones from the card.
    Also converts empty strings to None.
    Allows floats & ints, but that's not the intended value,
    though it is ok (it's just extra work).

    Parameters
    ----------
    card : List[str]
        the fields on the card as a list

    Returns
    -------
    short_card : List[str]
        the card with no trailing blank fields

    """
    short_card = []  # type: List[Union[str, int, float, None]]
    for field in card:
        if isinstance(field, str):
            field = field.strip()
            if field == '':
                short_card.append(None)
            else:
                short_card.append(field)
        else:
            short_card.append(field)

    imax = 0
    for i, field in enumerate(card):
        if short_card[i] is not None:
            imax = i
    out = short_card[:imax + 1]
    return out
