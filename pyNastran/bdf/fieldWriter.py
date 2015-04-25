"""
Defines legacy import functions
"""
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double

def print_card(fields, size=8, is_double=False):
    """
    Prints a card in 8-character of 16-character Nastran format.

    :param fields:    all the fields in the BDF card (no trailing Nones)
    :param size:      8/16
    :param is_double: True/False
    :returns card: string representation of the card

    .. note:: be careful of using is_double on cards that aren't
              GRID or COORDx
    """
    if size == 8:
        return print_card_8(fields)
    elif is_double:
        return print_card_double(fields)
    return print_card_16(fields)
