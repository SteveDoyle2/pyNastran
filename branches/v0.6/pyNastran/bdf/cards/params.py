## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from pyNastran.bdf.cards.baseCard import BaseCard
from pyNastran.bdf.bdfInterface.BDF_Card import BDFCard
from pyNastran.bdf.bdfInterface.assign_type import (integer_or_blank,
    double_or_blank, string, string_or_blank,
    integer_double_or_string)

class PARAM(BaseCard):
    type = 'PARAM'

    def __init__(self, card, data=None, comment=''):
        if comment:
            self._comment = comment
        if data:
            card = BDFCard(['PARAM'] + data)

        self.key = string(card, 1, 'key')
        n = 1
        value = None
        if self.key == 'ACOUT':
            value = string_or_blank(card, 2, 'value', 'PEAK')
        elif self.key == 'ACOWEAK':
            value = string_or_blank(card, 2, 'value', 'NO')
        elif self.key == 'ACSYM':
            value = string_or_blank(card, 2, 'value', 'YES')
        elif self.key == 'ADJMETH':
            value = integer_or_blank(card, 2, 'value', 0)
        elif self.key == 'ADMPOST':
            value = string_or_blank(card, 2, 'value', 0)
        elif self.key == 'ADSTAT':
            value = string_or_blank(card, 2, 'value', 'YES')
        elif self.key in ['ALPHA1', 'ALPHA2']:
            value1 = double_or_blank(card, 2, 'value1', 0.0)
            value2 = double_or_blank(card, 2, 'value2', 0.0)
            n = 2
        elif self.key in ['ALPHA1FL', 'ALPHA2FL']:
            value1 = double_or_blank(card, 2, 'value1', 0.0)
            value2 = double_or_blank(card, 3, 'value2', 0.0)
            n = 2
        elif self.key in ['CB1', 'CB2', 'CK1', 'CK2', 'CK3', 'CM1', 'CM2', 'CP1', 'CP2']:
            value1 = double_or_blank(card, 2, 'value1', 1.0)
            value2 = double_or_blank(card, 3, 'value2', 0.0)
            n = 2
        else:
            value = integer_double_or_string(card, 2, 'value')

        if value is None:
            self.values = [value1, value2]
        else:
            self.values = [value]

        if n == 1:
            if len(card) != 3:
                raise RuntimeError('len(PARAM card)=%i card=%r' % (len(card), card))
        else:
            if len(card) != 4:
                raise RuntimeError('len(PARAM card)=%i card=%r' % (len(card), card))

    def update_value(value1, value2=None):
        self.values = [value1]
        if value2 is not None:
            self.values.append(value2)
            
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
