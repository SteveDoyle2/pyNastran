"""
Defines various utilities for replicated BDF parsing including:
 - to_fields

"""
from typing import List, Optional

from pyNastran.bdf.bdf_interface.utils import expand_tabs
from pyNastran.bdf.cards.utils import wipe_empty_fields
from pyNastran.bdf.errors import ReplicationError


def to_fields_replication(card_lines: List[str]) -> List[Optional[str]]:
    """
    Converts a series of lines in a card into string versions of the field.
    Handles large, small, and CSV formatted cards.  Same as to_fields, but
    uses a different method to determine if it's large or small field format.

    Parameters
    ----------
    lines : List[str]
        the lines of the BDF card object

    Returns
    -------
    fields : List[str]
        the string formatted fields of the card

    .. warning:: this function is used by the reader and isn't intended
                 to be called by a separate process

    .. code-block:: python

       >>> card_lines = ['GRID,1,,1.0,2.0,3.0']
       >>> card_name = 'GRID'
       >>> fields = to_fields_replication(lines)
       >>> fields
       ['GRID', '1', '', '1.0', '2.0', '3.0']

    """
    #print('to_fields_replicationA =', card_lines)
    #line0 = card_lines[0]
    fields = []
    for iline, line in enumerate(card_lines):
        line = line.rstrip()

        if '\t' in line:
            line = expand_tabs(line)
        #print('  line = %r' % line)
        if ',' in line:
            assert '\t' not in line, '%r' % line
            sline = line.split(',')
            #print('sline = ', iline, sline)
            if iline == 0:
                card_name = sline[0]
                #assert card_name != '=', card_lines
                fields.append(card_name)

            if '*' in sline[0]:
                fields.extend(sline[1:5])
                for unused_i in range(5 - len(sline)):
                    fields.append('')
            else:
                fields.extend(sline[1:9])
                for unused_i in range(9 - len(sline)):
                    fields.append('')
        else:
            assert '\t' not in line, line
            if iline == 0:
                card_name = line[0:8]
                assert card_name != '=', card_lines
                fields.append(card_name)

            if '*' in card_name:
                fields.extend([line[8:24], line[24:40], line[40:56],
                               line[56:72]])
            else:
                fields.extend([line[8:16], line[16:24], line[24:32],
                               line[32:40], line[40:48], line[48:56], line[56:64],
                               line[64:72]])
    if '*' in card_name:
        raise ReplicationError(f'* found in unexpected position; {card_name!r}\nlines = {card_lines}')

    wiped_fields = wipe_empty_fields(fields)
    for field in fields:
        sfield = field.strip()
        #while '= ' in sfield:
            #sfield = field.replace('= ','=')
            #print('sfield=%r' % sfield)
        if ' ' in sfield:
            raise RuntimeError(f'field={sfield!r} has embedded blanks '
                               f'(mixed comma/space separated fields)\nfields={fields}')
    return wiped_fields

def get_nrepeats(field: str, old_card: List[str], new_card: List[str]) -> int:
    """=4, =(11)"""
    msg = f'field={field!r}; expected =(1), =2, ...\nold_card={old_card}\nnew_card={new_card}'
    assert field[0] == '=', msg
    assert '.' not in field, msg
    assert '*' not in field, msg
    if '(' in field or ')' in field:
        assert field[1] == '(', msg
        assert field[-1] == ')', msg
        fieldi = field[2:-1]
        assert '(' not in fieldi, msg
        assert ')' not in fieldi, msg
    else:
        assert '(' not in field, msg
        assert ')' not in field, msg
        fieldi = field[1:]
    nrepeats = int(fieldi)
    return nrepeats

def float_replication(field: str, old_field: str) -> float:
    """*4., *(11.5)"""
    msg = f'field={field!r}; expected *(1.), *2., ..., *11.'
    assert field[0] == '*', msg
    assert '.' in field, msg
    assert len(field) >= 3, msg
    if '(' in field:
        assert field[1] == '(', msg
        assert field[-1] == ')', msg
        fieldi = field[2:-1]
        assert '(' not in fieldi, msg
        assert ')' not in fieldi, msg
    else:
        assert '(' not in field, msg
        assert ')' not in field, msg
        fieldi = field[1:]
    nfloat = float(fieldi)
    field2 = nfloat + float(old_field)
    return field2

def int_replication(field: str, old_field: str) -> int:
    """*4, *(11)"""
    msg = f'field={field!r}; expected *(1), *2, ..., *11'
    assert field[0] == '*', msg
    assert '.' not in field, msg
    if '(' in field:
        assert field[1] == '(', msg
        assert field[-1] == ')', msg
        fieldi = field[2:-1]
        assert '(' not in fieldi, msg
        assert ')' not in fieldi, msg
    else:
        assert '(' not in field, msg
        assert ')' not in field, msg
        fieldi = field[1:]

    assert '*' not in fieldi, msg
    assert len(field) >= 2, msg
    #assert len(field) == 2, 'field=%r; expected *1, *2, *3, ..., *9' % field
    # integer
    nint = int(fieldi)
    field2 = nint + int(old_field)
    return field2

def _field(old_card: List[str], ifield: int) -> str:
    """helper for replication"""
    #if isinstance(old_card, list):
    #print(old_card, ifield)
    field2 = old_card[ifield]
    #else:
        #field2 = old_card.field(ifield)
    return field2

def repeat_cards(old_card: List[str], new_card: List[str]) -> List[List[str]]:
    """helper for replication"""
    card = []
    cards = []
    #print('*old_card = %s' % old_card)
    #print('*new_card = %s' % new_card)
    assert old_card != new_card
    for ifield, field in enumerate(new_card):
        if field is None:
            field2 = _field(old_card, ifield)
            #field2 = old_card.field(ifield)
            #print(' %i: %r -> %r' % (ifield, field, field2))
            #assert field2 is None, 'field=%s field2=%s' % (field, field2)
            card.append(field2)
            continue

        #1. Duplication of fields from the preceding entry is accomplished
        #   by coding the symbol =.
        #2. Duplication of all trailing fields from the preceding entry is
        #   accomplished by coding the symbol
        if field == '':
            field2 = field
        elif field == '=':
            field2 = _field(old_card, ifield)
        elif field == '==':
            # just append the remaining fields
            #print(' %s : extending %s' % (ifield, old_card[ifield:]))
            card.extend(old_card[ifield:])
            break

        elif '*' in field:
            # this is an increment, not multiplication...
            old_field = _field(old_card, ifield)
            #old_field = old_card.field(ifield)
            if '.' in field:
                assert old_field is not None, f'old_card:{old_card}\nnew_card:\n{new_card}'
                field2 = float_replication(field, old_field)
            else:
                assert old_field is not None, f'old_card:{old_card}\nnew_card:\n{new_card}'
                field2 = int_replication(field, old_field)
        else:
            msg = f'field={field!r}\nold_card={old_card}\nnew_card={new_card}'
            assert '(' not in field, msg
            assert '*' not in field, msg
            assert '=' not in field, msg
            field = field2
        #print(' %i: %r -> %r' % (ifield, field, field2))
        card.append(field2)
    #print(' appending %s' % card)
    cards.append(card)
    return cards
