# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from pyNastran.bdf.cards.utils import wipe_empty_fields
from six import string_types
from six.moves import range
import sys
import copy


class BDFCard(object):
    def __init__(self, card=None, debug=False):
        self.debug = debug
        if card:
            self.card = wipe_empty_fields(card)
            self.nfields = len(self.card)
        else:
            self.old_card = None
            self.card = None
            self.nfields = None

    def pop(self):
        """
        Pops the last value off
        """
        self.nfields -= 1
        return self.card.pop()

    def __setitem__(self, key, value):  # card[4] = value
        self.card.__setitem__(key, value)

    def __getitem__(self, key):  # print card[5]
        return self.card.__getitem__(key)

    def __getslice__(self, i, j):  # card[1:10]
        return self.card.__getslice__(i, j)

    def __setslice__(self, i, j, sequence):  # card[1:10] = 2
        self.card.__setslice__(i, j, sequence)

    def index(self, i):
        return self.card.index(i)

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
        if(i < self.nfields and self.card[i] is not None and self.card[i] is not ''):
            return self.card[i]
        else:
            return default

    def _replace_expression(self, field_new, field_old, replace_char='=',
                          replace_char2=''):
        """Used for nastran = format"""
        field_new = field_new.replace(replace_char, str(field_old) + replace_char2)
        typeOld = type(field_old)
        if isinstance(field_old, int) or isinstance(field_old, float):
            field_new = typeOld(eval(field_new))
        else:
            field_new = str(field_new)
        return field_new

    def _is_same_name(self):
        """Used for nastran = format"""
        if '=' in self.card[0]:
            return True
        return False

    def _apply_old_fields(self, card_count=0):
        """Used for nastran = format"""
        if not self._is_same_name():
            return

        self.cardCount = card_count
        stop = False
        self.card_text_old = self.card
        if self.nfields == 1 and self.old_card.nfields > 1:
            #self.nfields = self.oldCard.nfields
            #self._apply_old_fields(self, card_count=1)
            card_count = self.old_card.cardCount + card_count
            if self.debug:
                print("newCount = %s" % (card_count))
            self.card = copy.deepcopy(self.old_card.card_text_old)
            self.nfields = len(self.card)
            self.old_card.nfields = len(self.old_card.card)

            if self.debug:
                print("oldCard = %s" % (self.old_card))
                print("selfCard = %s" % (self.card))
            #stop = True

        fields_new = self.fields()
        fields_old = self.old_card.fields()

        max_length = max(self.nFields(), self.old_card.nFields())
        min_length = min(self.nFields(), self.old_card.nFields())

        card_built = [fields_old[0]]

        for i in range(1, min_length):
            field_old = fields_old[i]
            field_new = fields_new[i]

            a = "|%s|" % (field_new)
            if '*' in field_new:
                new_char = '+%s*' % (card_count + 1)
                field_new = self._replace_expression(field_new, field_old, '*', new_char)
            elif '/' in field_new:
                new_char = '-%s*' % (card_count + 1)
                field_new = self._replace_expression(field_new, field_old, '/', new_char)
            elif '==' == field_new:
                #break
                field_new = field_old
            elif '=' in field_new:
                if field_new == '=':
                    field_new = field_old
                elif field_new == '==':  # handle this in the max length section
                    pass
                else:  # replace = with str(expression)
                    field_new = self._replace_expression(field_new, field_old)

            elif '' == field_new:
                field_new = field_old
            else:
                b = "|%s|" % field_new
                c = "|%s|" % field_old
                print("i=%s field_start %-10s field_new %-10s field_old %-10s" %
                     (i, a, b, c))
                raise RuntimeError('unhandled case...')
            if self.debug:
                b = "|%s|" % field_new
                c = "|%s|" % field_old
                print("i=%s field_start %-10s field_new %-10s field_old %-10s" %
                     (i, a, b, c))
            card_built.append(field_new)
            i += 1

        if max_length < len(card_built):
            # the new card is longer than card_built
            for i in range(self.nfields, max_length):
                card_built.append(self.card[i])
        elif len(card_built) < self.old_card.nfields:
            # card_built is shorter than the old card
            for i in range(self.nfields, max_length):
                card_built.append(self.old_card.field(i))
        #else: # same length
            #pass

        if self.debug:
            print("cardBuilt = %s" % card_built)
        self.card = card_built
        self.nfields = len(self.card)

        if stop:
            sys.exit("stopping in _apply_old_fields")

    def get_old_field(self, i):
        """Used for nastran = format"""
        return self.old_card.field(i)
