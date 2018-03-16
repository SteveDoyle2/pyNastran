from __future__ import print_function, absolute_import
from six import iteritems, itervalues, iterkeys
from six.moves import range

import tables
import numpy as np

from .card_table import CardTable, TableDef
from ..data_helper import DataHelper


class Partition(object):
    def __init__(self, h5n, input):
        self._h5n = h5n
        self._input = input

        self.set1 = SET1(self._h5n, self)
        self.uset = USET(self._h5n, self)

    def path(self):
        return self._input.path() + ['PARTITION']

    def read(self):
        for key, item in iteritems(self.__dict__):
            if key.startswith('_'):
                continue
            try:
                item.read()
            except AttributeError:
                pass


########################################################################################################################


class SET1(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PARTITION/SET1/IDENTITY')

    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())

        result = {
            'IDENTITY': {
                'SID': [], 'GIDS_POS': [], 'GIDS_LEN': [], 'DOMAIN_ID': []
            },
            'GIDS': {
                'IDENTITY': {
                    'G1': []
                }
            },
            '_subtables': ['GIDS']
        }

        g1 = result['GIDS']['IDENTITY']['G1']
        identity = result['IDENTITY']
        sid = identity['SID']
        gids_pos = identity['GIDS_POS']
        gids_len = identity['GIDS_LEN']

        pos = 0
        for card_id in card_ids:
            card = cards[card_id]

            sid.append(card.sid)
            gids_pos.append(pos)
            _len = len(card.ids)
            gids_len.append(_len)
            pos += _len
            g1.extend(card.ids)

        return result


########################################################################################################################


class USET(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PARTITION/USET')

    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())

        result = {
            'IDENTITY': {
                'SNAME': [], 'ID': [], 'C': [], 'DOMAIN_ID': []
            },
        }

        identity = result['IDENTITY']
        sname = identity['SNAME']
        id_ = identity['ID']
        c = identity['C']

        for card_id in card_ids:
            card_list = cards[card_id]

            for card in card_list:
                for i in range(len(card.ids)):
                    sname.append(card.name)
                    id_.append(card.ids[i])
                    c.append(card.components[i])

        return result