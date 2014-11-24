from pyNastran.bdf.fieldWriter import print_card_8
from pyNastran.bdf.fieldWriter16 import print_card_16
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, blank, string, string_or_blank)


class NLPCI(object):
    type = 'NLPCI'

    def __init__(self):
        pass

    def add(self, card=None, comment=''):
        if comment:
            self._comment = comment
        self.nlpci_id = integer(card, 1, 'nlparm_id')
        self.Type = string_or_blank(card, 2, 'Type', 'CRIS')
        self.minalr = double_or_blank(card, 3, 'minalr', 0.25)
        self.maxalr = double_or_blank(card, 4, 'maxalr', 4.0)
        self.scale = double_or_blank(card, 5, 'scale', 0.0)
        blank(card, 6, 'blank')
        self.desiter = double_or_blank(card, 7, 'minalr', 12)
        self.mxinc = integer_or_blank(card, 8, 'minalr', 20)

    #def raw_fields(self):
        #list_fields =
        #return list_fields

    def repr_fields(self):
        #minalr = set_blank_if_default(self.minalr, 0.25)
        return self.rawFields()

    def write_bdf(self, f, size=8):
        card = ['NLPCI', self.nlpci_id, self.Type, self.minalr,
                self.maxalr, self.scale, None, self.desiter, self.mxinc]
        if size == 8:
            f.write(print_card_8(card))
        else:
            f.write(print_card_16(card))
