# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from pyNastran.bdf.cards.baseCard import BaseCard
from pyNastran.bdf.bdfInterface.BDF_Card import BDFCard
from pyNastran.bdf.bdfInterface.assign_type import (integer_or_blank,
    double_or_blank, string, string_or_blank,
    integer_double_or_string)
from pyNastran.bdf.fieldWriter import print_card_8
from pyNastran.bdf.fieldWriter16 import print_card_16


class PARAM(BaseCard):
    type = 'PARAM'
    _field_map = {1: 'key'}

    def _update_field_helper(self, n, value):
        if n - 2 >= 0:
            try:
                self.values[n - 2] = value
            except IndexError:
                raise IndexError('Field %r=%r is an invalid %s entry for key=%r.' % (n, value, self.type, self.key.upper()))
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, card, data=None, comment=''):
        """
        Creates a PARAM card.

        :param self: the object
        :param card: BDFCard object
        :param data: list of PARAM entries not including 'PARAM';
                     intended to be used by OP2 Reader (default=None)
        :param comment: optional string (default='')
        """
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
        elif self.key == 'ADPCON':
            value = float_or_blank(card, 2, 'value', 1.0)
        #elif self.key == 'ADMPOST':
            #value = string_or_blank(card, 2, 'value', 0) ## TODO: 0 is not a string
        elif self.key == 'ADSDISC':
            value = float_or_blank(card, 2, 'value', 1e-8)
        elif self.key == 'AESMAXIT':
            value = integer_or_blank(card, 2, 'value', 15)
        elif self.key == 'AESMETH':
            value = string_or_blank(card, 2, 'value', 'SELECT')
            assert value in ['SELECT', 'AUTO', 'DIRECT', 'RITZ', 'ITER'], 'value=%s' % value
        elif self.key == 'AESTOL':
            value = ifloat_or_blank(card, 2, 'value', 1e-10)
        elif self.key == 'ADSTAT':
            value = string_or_blank(card, 2, 'value', 'YES')
        elif self.key in ['ALPHA1', 'ALPHA2', 'ALPHA1FL', 'ALPHA2FL']:  # check alpha1/alpha1FL
            value1 = double_or_blank(card, 2, 'value1', 0.0)
            value2 = double_or_blank(card, 3, 'value2', 0.0)
            n = 2
        elif self.key in ['CB1', 'CB2', 'CK1', 'CK2', 'CK3', 'CM1', 'CM2', 'CP1', 'CP2']:
            value1 = double_or_blank(card, 2, 'value1', 1.0)
            value2 = double_or_blank(card, 3, 'value2', 0.0)
            n = 2
        elif self.key == 'POST':
            value = integer_or_blank(card, 2, 'value', 1)
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

    def update_values(value1=None, value2=None):
        """
        Updates value1 and value2.  Performs type checking based on the PARAM
        type after setting any default value(s).

        :param self:   the PARAM object
        :param value1: the main value (default=None)
        :param value2: optional value (default=None)

        If you want to access the data directly, use:
        >>>  param = bdf.params['POST']
        >>> param.values[0] = -1  # value1
        >>> param.values[1] = 3   # value2
        >>>

        .. note::  Most PARAM cards only have one value.  Some have two.
        """
        if self.key == 'ACOUT':
            if value1 is None:
                value1 = 'PEAK'
            if not isinstance(value1, str):
                msg = 'key=%s value1=%r must be an string.' % (key, value1)
                raise TypeError(msg)

        elif self.key == 'ACOWEAK':
            if value1 is None:
                value1 = 'NO'
            if not isinstance(value1, str):
                msg = 'key=%s value1=%r must be an string.' % (key, value1)
                raise TypeError(msg)

        elif self.key == 'ACSYM':
            if value1 is None:
                value1 = 'YES'
            if not isinstance(value1, str):
                msg = 'key=%s value1=%r must be an string.' % (key, value1)
                raise TypeError(msg)

        elif self.key == 'ADJMETH':
            if value1 is None:
                value1 = 0
            if not isinstance(value1, int):
                msg = 'key=%s value1=%r must be an integer.' % (key, value1)
                raise TypeError(msg)

        #elif self.key == 'ADMPOST': ## TODO: 0 is not a string
            #value = string_or_blank(card, 2, 'value', 0)

        elif self.key == 'ADSTAT':
            if value1 is None:
                value1 = 'YES'
            if not isinstance(value1, str):
                msg = 'key=%s value1=%r must be an string.' % (key, value1)
                raise TypeError(msg)

        elif self.key in ['ALPHA1', 'ALPHA2', 'ALPHA1FL', 'ALPHA2FL']:
            if value1 is None:
                value1 = 0.0
            if value2 is None:
                value2 = 0.0
            if not isinstance(value1, float):
                msg = 'key=%s value1=%r must be an float.' % (key, value1)
                raise TypeError(msg)
            if isinstance(value2, float):
                msg = 'key=%s value2=%r must be an float.' % (key, value2)
                raise TypeError(msg)

        elif self.key in ['CB1', 'CB2', 'CK1', 'CK2', 'CK3', 'CM1', 'CM2', 'CP1', 'CP2']:
            if value1 is None:
                value1 = 1.0
            if value2 is None:
                value2 = 0.0
            if not isinstance(value1, float):
                msg = 'key=%s value1=%r must be an float.' % (key, value1)
                raise TypeError(msg)
            if isinstance(value2, float):
                msg = 'key=%s value2=%r must be an float.' % (key, value2)
                raise TypeError(msg)

        else:
            if not (isinstance(value1, int) or isinstance(value1, float) or
                    isinstance(value1, str)):
                msg = 'key=%s value1=%r must be an integer, float, or string.' % (key, value1)
                raise TypeError(msg)

        self.values = [value1]
        if value2 is not None:
            self.values.append(value2)

    #def isSameCard(self, param, debug=False):
        #fields1 = [self.key] + self.values
        #fields2 = [param.key] + param.values
        #for (field1, field2) in izip(fields1, fields2):
        #    if not self.isSame(field1, field2):
        #        return False
        #return True

    def rawFields(self):
        list_fields = ['PARAM', self.key] + self.values
        return list_fields

    def reprFields(self):
        return self.rawFields()

    def write_bdf(self, size, card_writer):
        card = self.rawFields()
        if size == 8:
            return self.comment() + print_card_8(card)  # works
        return self.comment() + print_card_16(card)
        #return self.comment() + card_writer(card)