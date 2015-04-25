"""
Defines the BDFCard class that is passed into the various Nastran cards.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from pyNastran.bdf.cards.utils import wipe_empty_fields
from six.moves import range


class BDFCard(object):
    """
    A BDFCard is a list that has a default value of None for fields out of
    range.
    """
    def __init__(self, card=None, debug=False):
        self.debug = debug
        if card:
            self.card = wipe_empty_fields(card)
            self.nfields = len(self.card)
        else:
            self.card = None
            self.nfields = None

    def pop(self):
        """
        Pops the last value off
        """
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

        :param self:  the object pointer
        :returns msg: the string representation of the card
        """
        return str(self.card)

    def nFields(self):
        """
        Gets how many fields are on the card

        :param self:      the object pointer
        :returns nfields: the number of fields on the card
        """
        return self.nfields

    def __len__(self):
        """len(card)"""
        return self.nfields

    def fields(self, i=0, j=None, defaults=None):
        """
        Gets multiple fields on the card

        :param self:     the object pointer
        :param i:        the ith field on the card (following list notation)
        :type i:         integer >= 0
        :param j:        the jth field on the card (None means till the end
                         of the card)
        :type j:         integer or None (default=end of card)
        :param defaults: the default value for the field (as a list)
                         len(defaults)=i-j-1

        :returns value: the values on the ith-jth fields
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

        :param self:    the object pointer
        :param i:       the ith field on the card (following list notation)
        :type i:        integer
        :param default: the default value for the field

        :returns value: the value on the ith field
        """
        if i < self.nfields and self.card[i] is not None and self.card[i] is not '':
            return self.card[i]
        else:
            return default
