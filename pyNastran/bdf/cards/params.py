# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import string_types, integer_types
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer_or_blank, double_or_blank, string, string_or_blank,
    integer_double_string_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16

INT_WORDS_1 = [
    'POST', 'OPPHIPA', 'OPPHIPB', 'GRDPNT', 'RPOSTS1', 'BAILOUT',
    'COUPMASS', 'CURV', 'INREL', 'MAXRATI', 'OG',
    'S1AM', 'S1M', 'DDRMM', 'MAXIT', 'PLTMSG', 'LGDISP', 'NLDISP',
    'OUNIT2M', 'RESCOMP', 'PDRMSG', 'LMODES', 'USETPRT', 'NOCOMPS',]
#float_words_1 = [
    #b'K6ROT', b'WTMASS', b'SNORM', b'PATVER', b'MAXRATIO', b'EPSHT',
    #b'SIGMA', b'TABS']
STR_WORDS_1 = [
    'POSTEXT', 'PRTMAXIM', 'AUTOSPC', 'OGEOM', 'PRGPST',
    'RESVEC', 'RESVINER', 'ALTRED', 'OGPS', 'OIBULK', 'OMACHPR',
    'UNITSYS', 'F56', 'OUGCORD', 'OGEM', 'EXTSEOUT',]
INT_STR_WORDS_1 = INT_WORDS_1 + STR_WORDS_1

SMALL_FIELD_PARAMS = [
    'ACOUT', 'ACOWEAK', 'ACSYM', 'ADJMETH', 'AESMAXIT', 'AESMETH', 'ADSTAT',
    'MAXLINES'] #+ INT_WORDS_1 + STR_WORDS_1


class PARAM(BaseCard):
    type = 'PARAM'
    _field_map = {1: 'key'}

    def _update_field_helper(self, n, value):
        if n - 2 >= 0:
            try:
                self.values[n - 2] = value
            except IndexError:
                msg = 'Field %r=%r is an invalid %s entry for key=%r.' % (
                    n, value, self.type, self.key.upper())
                raise IndexError(msg)
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    @classmethod
    def _init_from_empty(cls):
        key = 'POST'
        values = -1
        return PARAM(key, values, comment='')

    def __init__(self, key, values, comment=''):
        """
        Creates a PARAM card

        Parameters
        ----------
        key : str
            the name of the PARAM
        values : int/float/str/List
            varies depending on the type of PARAM
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        self.key = key
        if isinstance(values, (int, float, str)):
            values = [values]
        self.values = values

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PARAM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """

        key = string(card, 1, 'key')

        n = 1
        value = None
        if key == 'ACOUT':
            value = string_or_blank(card, 2, 'value', 'PEAK')
        elif key == 'ACOWEAK':
            value = string_or_blank(card, 2, 'value', 'NO')
        elif key == 'ACSYM':
            value = string_or_blank(card, 2, 'value', 'YES')
        elif key == 'ADJMETH':
            value = integer_or_blank(card, 2, 'value', 0)
        elif key == 'ADPCON':
            value = double_or_blank(card, 2, 'value', 1.0)
        #elif key == 'ADMPOST':
            #value = string_or_blank(card, 2, 'value', 0) ## TODO: 0 is not a string
        elif key == 'ADSDISC':
            value = double_or_blank(card, 2, 'value', 1e-8)
        elif key == 'AESMAXIT':
            value = integer_or_blank(card, 2, 'value', 15)
        elif key == 'AESMETH':
            value = string_or_blank(card, 2, 'value', 'SELECT')
            assert value in ['SELECT', 'AUTO', 'DIRECT', 'RITZ', 'ITER'], 'value=%s' % value
        elif key == 'AESTOL':
            value = double_or_blank(card, 2, 'value', 1e-10)
        elif key == 'ADSTAT':
            value = string_or_blank(card, 2, 'value', 'YES')
        elif key in ['ALPHA1', 'ALPHA2', 'ALPHA1FL', 'ALPHA2FL']:  # check alpha1/alpha1FL
            value1 = double_or_blank(card, 2, 'value1', 0.0)
            value2 = double_or_blank(card, 3, 'value2', 0.0)
            n = 2
        elif key in ['CB1', 'CB2', 'CK1', 'CK2', 'CK3', 'CK41', 'CK42',
                     'CM1', 'CM2', 'CP1', 'CP2']:
            value1 = double_or_blank(card, 2, 'value1', 1.0)
            value2 = double_or_blank(card, 3, 'value2', 0.0)
            n = 2
        elif key == 'POST':
            value = integer_or_blank(card, 2, 'value', 1)
        elif key == 'UNITSYS':
            value = string(card, 2, 'value')
        else:
            n = 2
            value1 = integer_double_string_or_blank(card, 2, 'value1')
            value2 = integer_double_string_or_blank(card, 3, 'value2')
            if value2 is None:
                value = value1
                n = 1

        if value is None:
            # n=2 or blank
            if isinstance(value1, string_types):
                assert ' ' not in value1, 'PARAM value1=%r' % value1
            if isinstance(value2, string_types):
                assert ' ' not in value2, 'PARAM value2=%r' % value2
            values = [value1, value2]
        else:
            # n=1
            if isinstance(value, string_types):
                assert ' ' not in value, 'PARAM value=%r' % value
            values = [value]

        if n == 1:
            assert len(card) <= 3, 'len(PARAM card)=%i card=%r' % (len(card), card)
        else:
            assert len(card) <= 4, 'len(PARAM card)=%i card=%r' % (len(card), card)
        return PARAM(key, values, comment=comment)

    def update_values(self, value1=None, value2=None):
        """
        Updates value1 and value2.  Performs type checking based on the PARAM
        type after setting any default value(s).

        Parameters
        ----------
        value1 : varies; default=None
            the main value
        value2 : varies; default=None
            optional value

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
            if not isinstance(value1, string_types):
                msg = 'key=%s value1=%r must be an string.' % (self.key, value1)
                raise TypeError(msg)

        elif self.key == 'ACOWEAK':
            if value1 is None:
                value1 = 'NO'
            if not isinstance(value1, string_types):
                msg = 'key=%s value1=%r must be an string.' % (self.key, value1)
                raise TypeError(msg)

        elif self.key == 'ACSYM':
            if value1 is None:
                value1 = 'YES'
            if not isinstance(value1, string_types):
                msg = 'key=%s value1=%r must be an string.' % (self.key, value1)
                raise TypeError(msg)

        elif self.key == 'ADJMETH':
            if value1 is None:
                value1 = 0
            if not isinstance(value1, integer_types):
                msg = 'key=%s value1=%r must be an integer.' % (self.key, value1)
                raise TypeError(msg)

        #elif self.key == 'ADMPOST': ## TODO: 0 is not a string
            #value = string_or_blank(card, 2, 'value', 0)

        elif self.key == 'ADSTAT':
            if value1 is None:
                value1 = 'YES'
            if not isinstance(value1, string_types):
                msg = 'key=%s value1=%r must be an string.' % (self.key, value1)
                raise TypeError(msg)

        elif self.key in ['ALPHA1', 'ALPHA2', 'ALPHA1FL', 'ALPHA2FL']:
            if value1 is None:
                value1 = 0.0
            if value2 is None:
                value2 = 0.0
            if not isinstance(value1, float):
                msg = 'key=%s value1=%r must be an float.' % (self.key, value1)
                raise TypeError(msg)
            if isinstance(value2, float):
                msg = 'key=%s value2=%r must be an float.' % (self.key, value2)
                raise TypeError(msg)

        elif self.key in ['CB1', 'CB2', 'CK1', 'CK2', 'CK3', 'CM1', 'CM2', 'CP1', 'CP2']:
            if value1 is None:
                value1 = 1.0
            if value2 is None:
                value2 = 0.0
            if not isinstance(value1, float):
                msg = 'key=%s value1=%r must be an float.' % (self.key, value1)
                raise TypeError(msg)
            if isinstance(value2, float):
                msg = 'key=%s value2=%r must be an float.' % (self.key, value2)
                raise TypeError(msg)

        else:
            if not isinstance(value1, (int, float, string_types)):
                msg = 'key=%s value1=%r must be an integer, float, or string.' % (self.key, value1)
                raise TypeError(msg)

        self.values = [value1]
        if value2 is not None:
            self.values.append(value2)

    def raw_fields(self):
        list_fields = ['PARAM', self.key] + self.values
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        if self.key in INT_STR_WORDS_1:
            return '%sPARAM   %8s%8s\n' % (
                self.comment, self.key, self.values[0])

        if size == 8 or self.key in SMALL_FIELD_PARAMS:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)
