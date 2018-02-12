from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range

import tables
import numpy as np

from .card_table import CardTable, TableDef, TableData


class Load(object):
    def __init__(self, h5n, input):
        self._h5n = h5n
        self._input = input

        self.force = FORCE(self._h5n, self)
        self.load = LOAD(self._h5n, self)
        self.moment = MOMENT(self._h5n, self)

    def path(self):
        return self._input.path() + ['LOAD']

    def read(self):
        for key, item in iteritems(self.__dict__):
            if key.startswith('_'):
                continue
            try:
                item.read()
            except AttributeError:
                pass

########################################################################################################################


class FORCE(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/LOAD/FORCE')

    @staticmethod
    def from_bdf(card):
        data = []
        for _ in card:
            data.append([_.sid, _.node, _.cid, _.mag, _.xyz])
        return TableData(list(sorted(data, key=lambda x: x[1])))

########################################################################################################################


class LOAD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/LOAD/LOAD/IDENTITY')

    @staticmethod
    def from_bdf(card):
        # TODO: rewrite all bdf cards this way
        # TODO: domain ids are superelement ids, need to address that
        
        # this is much more logical than the TableData and subdata stuff

        sfactors = {
            'IDENTITY': {'SI': [], 'LI': []}
        }
      
        data = {
            'IDENTITY': {'SID': [], 'S': [], 'SFACTORS_POS': [], 'SFACTORS_LEN': []},
            'SFACTORS': sfactors,
            '_subtables': ['SFACTORS']
        }
        
        identity = data['IDENTITY']
        
        sid = identity['SID']
        s = identity['S']
        pos = identity['SFACTORS_POS']
        _len = identity['SFACTORS_LEN']
        sfactors = sfactors['IDENTITY']
        si = sfactors['SI']
        li = sfactors['LI']

        cards = sorted(card, key=lambda x: x.sid)
        
        _pos = 0        
        for card in cards:
            sid.append(card.sid)
            s.append(card.scale)
            pos.append(_pos)
            data_len = len(card.scale_factors)
            _pos += data_len
            _len.append(data_len)
            si += list(card.scale_factors)
            li += list(card.load_ids)

        return data


########################################################################################################################


class MOMENT(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/LOAD/MOMENT')

    @staticmethod
    def from_bdf(card):
        data = []
        for _ in card:
            data.append([_.sid, _.node, _.cid, _.mag, _.xyz])
        return TableData(list(sorted(data, key=lambda x: x[1])))
