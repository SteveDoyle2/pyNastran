from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range

import tables
import numpy as np

from .card_table import CardTable, TableDef
from ..data_helper import DataHelper


class Load(object):
    def __init__(self, h5n, input):
        self._h5n = h5n
        self._input = input

        self.force = FORCE(self._h5n, self)
        self.load = LOAD(self._h5n, self)
        self.moment = MOMENT(self._h5n, self)
        self.pload4 = PLOAD4(self._h5n, self)

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

    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())

        result = {
            'IDENTITY': {
                'SID': [], 'G': [], 'CID': [], 'F': [], 'N': [], 'DOMAIN_ID': []
            }
        }

        identity = result['IDENTITY']
        sid = identity['SID']
        g = identity['G']
        cid = identity['CID']
        f = identity['F']
        n = identity['N']

        for card_id in card_ids:
            cards = sorted(cards[card_id], key=lambda x: x.node)

            for card in cards:
                sid.append(card.sid)
                g.append(card.node)
                cid.append(card.cid)
                f.append(card.mag)
                n.append(card.xyz)

        identity['N'] = np.array(n)

        return result


########################################################################################################################


class LOAD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/LOAD/LOAD/IDENTITY')

    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())

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

        _pos = 0

        for card_id in card_ids:
            card_list = sorted(cards[card_id], key=lambda x: x.sid)

            for card in card_list:
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
    
    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())

        result = {
            'IDENTITY': {
                'SID': [], 'G': [], 'CID': [], 'M': [], 'N': [], 'DOMAIN_ID': []
            }
        }

        identity = result['IDENTITY']
        sid = identity['SID']
        g = identity['G']
        cid = identity['CID']
        m = identity['M']
        n = identity['N']

        for card_id in card_ids:
            cards = sorted(cards[card_id], key=lambda x: x.node)

            for card in cards:
                sid.append(card.sid)
                g.append(card.node)
                cid.append(card.cid)
                m.append(card.mag)
                n.append(card.xyz)

        identity['N'] = np.array(n)

        return result


########################################################################################################################


class PLOAD4(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/LOAD/PLOAD4')

    @staticmethod
    def from_bdf(cards):
        data = {
            'IDENTITY': {
                'SID': [],
                'EID': [],
                'P': [],
                'G1': [],
                'G34': [],
                'CID': [],
                'N': [],
                'SORL': [],
                'LDIR': [],
                'DOMAIN_ID': []
            }
        }

        identity = data['IDENTITY']
        sid = identity['SID']
        eid = identity['EID']
        p = identity['P']
        g1 = identity['G1']
        g34 = identity['G34']
        cid = identity['CID']
        n = identity['N']
        sorl = identity['SORL']
        ldir = identity['LDIR']

        card_ids = sorted(cards.keys())

        for card_id in card_ids:
            card_list = cards[card_id]

            for card in card_list:
                eids = card.eids

                eid_len = len(eids)

                sid.extend([card.sid] * eid_len)
                eid.extend(eids)

                _p = list(card.pressures)
                if len(_p) == 3:
                    _p.append(0.)
                p.extend([_p] * eid_len)

                _g1 = card.g1
                if _g1 is None:
                    _g1 = DataHelper.default_int

                _g34 = card.g34
                if _g34 is None:
                    _g34 = DataHelper.default_int

                g1.extend([_g1] * eid_len)
                g34.extend([_g34] * eid_len)

                cid.extend([card.cid] * eid_len)
                n.extend([card.nvector] * eid_len)
                sorl.extend([card.surf_or_line] * eid_len)
                ldir.extend([card.line_load_dir] * eid_len)

        identity['P'] = np.array(p)

        return data
