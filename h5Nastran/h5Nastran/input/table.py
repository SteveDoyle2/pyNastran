from __future__ import print_function, absolute_import
from six import iteritems, itervalues, iterkeys
from six.moves import range

import tables
import numpy as np

from .card_table import CardTable, TableDef
from ..data_helper import DataHelper


class Table(object):
    def __init__(self, h5n, input):
        self._h5n = h5n
        self._input = input

        self.tabdmp1 = TABDMP1(self._h5n, self)
        self.tabled1 = TABLED1(self._h5n, self)
        self.tabled2 = TABLED2(self._h5n, self)
        self.tablem1 = TABLEM1(self._h5n, self)

    def path(self):
        return self._input.path() + ['TABLE']

    def read(self):
        for key, item in iteritems(self.__dict__):
            if key.startswith('_'):
                continue
            try:
                item.read()
            except AttributeError:
                pass


########################################################################################################################

# TODO: TABDMP1 is missing Type in msc spec... why?

class TABDMP1_DAMP(object):
    name = 'DAMP'
    path = '/NASTRAN/INPUT/TABLE/TABDMP1'
    dtype = [('F', '<f8', ()), ('G', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


class TABDMP1_SPEC(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/TABLE/TABDMP1'
    dtype = [('ID', '<i8', ()), ('TYPE', 'S4', ()), ('POS', '<i8', ()), ('LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = [TABDMP1_DAMP]


class TABDMP1(CardTable):
    table_def = TableDef.create(TABDMP1_SPEC, rename={'DAMP_POS': 'POS', 'DAMP_LEN': 'LEN'})

    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())

        damp = {'IDENTITY': {'F': [], 'G': []}}

        result = {'IDENTITY': {'ID': [], 'TYPE': [], 'POS': [], 'LEN': [], 'DOMAIN_ID': []},
                  'DAMP': damp,
                  '_subtables': ['DAMP']}

        f = damp['IDENTITY']['F']
        g = damp['IDENTITY']['G']
        identity = result['IDENTITY']
        id_ = identity['ID']
        type_ = identity['TYPE']
        pos = identity['POS']
        len_ = identity['LEN']

        _pos = 0
        for card_id in card_ids:
            card = cards[card_id]

            id_.append(card.tid)
            type_.append(card.Type)
            pos.append(_pos)
            _len = len(card.x)
            len_.append(_len)
            f += list(card.x)
            g += list(card.y)

        return result


########################################################################################################################


class TABLED1(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/TABLE/TABLED1/IDENTITY', rename={'XY_POS': 'POS', 'XY_LEN': 'LEN'})

    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())

        xy = {'IDENTITY': {'X': [], 'Y': []}}

        result = {'IDENTITY': {'ID': [], 'CODEX': [], 'CODEY': [], 'POS': [], 'LEN': [], 'DOMAIN_ID': []},
                  'XY': xy,
                  '_subtables': ['XY']}

        x = xy['IDENTITY']['X']
        y = xy['IDENTITY']['Y']
        identity = result['IDENTITY']
        id_ = identity['ID']
        codex = identity['CODEX']
        codey = identity['CODEY']
        pos = identity['POS']
        len_ = identity['LEN']

        # TODO: TABLED1 confirm codes are correct
        _xy_code = {'LINEAR': 1, 'LOG': 2, '': 1, None: 1}

        _pos = 0

        for card_id in card_ids:
            card = cards[card_id]

            id_.append(card.tid)
            codex.append(_xy_code[card.xaxis])
            codey.append(_xy_code[card.yaxis])
            pos.append(_pos)
            _len = len(card.x)
            len_.append(_len)
            _pos += _len
            x.extend(list(card.x))
            y.extend(list(card.y))

        return result

########################################################################################################################


class TABLED2(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/TABLE/TABLED2/IDENTITY', rename={'XY_POS': 'POS', 'XY_LEN': 'LEN'})

    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())

        xy = {'IDENTITY': {'X': [], 'Y': []}}

        result = {'IDENTITY': {'ID': [], 'X1': [], 'POS': [], 'LEN': [], 'DOMAIN_ID': []},
                  'XY': xy,
                  '_subtables': ['XY']}

        x = xy['IDENTITY']['X']
        y = xy['IDENTITY']['Y']
        identity = result['IDENTITY']
        id_ = identity['ID']
        x1 = identity['X1']
        pos = identity['POS']
        len_ = identity['LEN']

        _pos = 0
        for card_id in card_ids:
            card = cards[card_id]

            id_.append(card.tid)
            x1.append(card.x1)
            pos.append(_pos)
            _len = len(card.x)
            len_.append(_len)
            _pos += _len
            x.extend(list(card.x))
            y.extend(list(card.y))

        return result
    

########################################################################################################################


class TABLEM1(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/TABLE/TABLEM1/IDENTITY', rename={'XY_POS': 'POS', 'XY_LEN': 'LEN'})
    from_bdf = TABLED1.from_bdf
