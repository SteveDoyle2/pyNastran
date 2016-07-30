"""
Defines legacy import functions
"""
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double

def print_card(fields, size=8, is_double=False):
    """
    Prints a card in 8-character of 16-character Nastran format.

    Parameters
    ----------
    fields : List[int/float/str/None]
        all the fields in the BDF card (no trailing Nones)
    size : int; default=8
        the size of the field (8/16)
    is_double : bool; default=False
        is the card double precision?
        Double precision applies to specific cards and turns
        1.234E+5 into 1.234D+5.  Applies to GRID, CORDx only?

    Returns
    -------
    card : str
        string representation of the card

    .. note:: be careful of using is_double on cards that aren't
              GRID or CORDx
    """
    if size == 8:
        return print_card_8(fields)
    elif is_double:
        return print_card_double(fields)
    return print_card_16(fields)
