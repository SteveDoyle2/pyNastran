# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from pyNastran.bdf.cards.baseCard import BaseCard
from pyNastran.bdf.format import string, integer_double_or_string

class PARAM(BaseCard):
    type = 'PARAM'

    def __init__(self, card, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            self.key = string(card, 1, 'key')
            self.value = integer_double_or_string(card, 2, 'value')
        else:
            self.key = data[0]
            self.value = data[1]

    #def isSameCard(self, param, debug=False):
        #fields1 = [self.key, self.value ]
        #fields2 = [param.key, param.value]
        #for (field1, field2) in izip(fields1, fields2):
        #    if not self.isSame(field1, field2):
        #        return False
        #return True

    def rawFields(self):
        list_fields = ['PARAM', self.key, self.value]
        return list_fields
