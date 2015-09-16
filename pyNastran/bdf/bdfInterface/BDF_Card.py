"""
Defines the BDFCard class that is passed into the various Nastran cards.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from pyNastran.bdf.cards.utils import wipe_empty_fields
from pyNastran.bdf.field_writer_16 import print_field_16
from six.moves import range


class BDFCard(object):
    """
    A BDFCard is a list that has a default value of None for fields out of
    range.
    """
    def __init__(self, card, has_none=True):
        if has_none:
            card = wipe_empty_fields([print_field_16(field) for field in card])
        self.card = card
        self.nfields = len(self.card)

    def pop(self):
        """card.pop()"""
        self.nfields -= 1
        return self.card.pop()

    def __setitem__(self, key, value):
        """card[4] = value"""
        self.card.__setitem__(key, value)

    def __getitem__(self, key):
        """print card[5]"""
        return self.card.__getitem__(key)

    def __getslice__(self, i, j):
        """card[1:10]"""
        return self.card.__getslice__(i, j)

    def __setslice__(self, i, j, sequence):
        """card[1:10] = 2"""
        self.card.__setslice__(i, j, sequence)

    def index(self, value):
        """card.index(value)"""
        return self.card.index(value)

    def __repr__(self):
        """
        Prints the card as a list

        Parameters
        ----------
        self : BDFCard()
            the object pointer

        Returns
        -------
        msg : str
            the string representation of the card
        """
        return '%r' % self.card

    def nFields(self):
        """
        Gets how many fields are on the card

        Parameters
        ----------
        self : BDFCard()
            the object pointer

        Returns
        -------
        nfields : int
            the number of fields on the card
        """
        self.deprecated('self.nFields()', 'self.nfields', '0.8')
        return self.nfields

    def __len__(self):
        """len(card)"""
        return self.nfields

    def fields(self, i=0, j=None, defaults=None):
        """
        Gets multiple fields on the card

        Parameters
        ----------
        self : BDFCard()
            the object pointer
        i : int > 0
            the ith field on the card (following list notation)
        j : int / None
            int : the jth field on the card
            None : last field on the card
        defaults : List[int/float/str]
            the default value for the field (as a list)
            len(defaults)=i-j-1

        Returns
        -------
        value : int/float/str
            the values on the ith-jth fields
        """
        if defaults is None:
            defaults = []
        if j is None:
            if self.nfields is None:
                return [None]
            j = self.nfields

        if defaults == []:
            defaults = [None] * (j - i + 1)
        out = []

        d = 0
        for n in range(i, j):
            value = self.field(n, defaults[d])
            out.append(value)
            d += 1
        return out

    def field(self, i, default=None):
        """
        Gets the ith field on the card

        Parameters
        ----------
        self : BDFCard()
            the object pointer
        i : int
            the ith field on the card (following list notation)
        default : int/float/str/None
            the default value for the field

        Returns
        -------
        value : int/float/str/None
            the value on the ith field
        """
        if i < self.nfields and self.card[i] is not None and self.card[i] is not '':
            return self.card[i]
        else:
            if default is None:
                return u''
            return default
