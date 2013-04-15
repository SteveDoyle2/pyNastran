# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from pyNastran.bdf.cards.baseCard import BaseCard
from pyNastran.bdf.assign_type import (integer_or_blank,
    double_or_blank, string, string, string_or_blank,
    integer_double_or_string)

class PARAM(BaseCard):
    type = 'PARAM'

    def __init__(self, card, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            self.key = string(card, 1, 'key')
            n = 1
            if self.key == 'ACOUT':
                self.value = string_or_blank(card, 1, 'value', 'PEAK')
            elif self.key == 'ACOWEAK':
                self.value = string_or_blank(card, 1, 'value', 'NO')
            elif self.key == 'ACSYM':
                self.value = string_or_blank(card, 1, 'value', 'YES')
            elif self.key == 'ADJMETH':
                self.value = integer_or_blank(card, 1, 'value', 0)
            elif self.key == 'ADMPOST':
                self.value = string_or_blank(card, 1, 'value', 0)
            elif self.key == 'ADSTAT':
                self.value = string_or_blank(card, 1, 'value', 'YES')
            elif self.key in ['ALPHA1', 'ALPHA2']:
                self.value1 = double_or_blank(card, 1, 'value1', 0.0)
                self.value2 = double_or_blank(card, 1, 'value2', 0.0)
                n = 2
            elif self.key in ['ALPHA1FL', 'ALPHA2FL']:
                self.value1 = double_or_blank(card, 1, 'value1', 0.0)
                self.value2 = double_or_blank(card, 1, 'value2', 0.0)
                n = 2
            elif self.key in ['CB1', 'CB2', 'CK1', 'CK2', 'CK3', 'CM1', 'CM2', 'CP1', 'CP2']:
                self.value1 = double_or_blank(card, 1, 'value1', 1.0)
                self.value2 = double_or_blank(card, 1, 'value2', 0.0)
                n = 2
            else:
                self.value = integer_double_or_string(card, 2, 'value')

            if hasattr(self, 'values1'):
                self.values = [self.value1, self.value2]
            else:
                self.values = [self.value]

            if n == 1:
                assert len(card) == 3, 'len(PARAM card)=%i card=%r' % (len(card), card)
            else:
                assert len(card) == 4, 'len(PARAM card)=%i card=%r' % (len(card), card)
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
        list_fields = ['PARAM', self.key] + self.values
        return list_fields

    def reprFields(self):
        return self.rawFields()
