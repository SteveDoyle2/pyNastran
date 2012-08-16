# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from pyNastran.bdf.cards.baseCard import BaseCard


class PARAM(BaseCard):
    type = 'PARAM'

    def __init__(self, card):
        self.key = card.field(1)
        self.value = card.field(2)

    #def isSameCard(self, param, debug=False):
        #fields1 = [self.key, self.value ]
        #fields2 = [param.key, param.value]
        #for (field1, field2) in izip(fields1, fields2):
        #    if not self.isSame(field1, field2):
        #        return False
        #    ###
        ####
        #return True

    def rawFields(self):
        fields = ['PARAM', self.key, self.value]
        return fields
