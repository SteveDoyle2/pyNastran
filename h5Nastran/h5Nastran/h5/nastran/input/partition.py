from __future__ import print_function, absolute_import

from six.moves import range

from h5Nastran.defaults import Defaults
from h5Nastran.h5nastrannode import H5NastranNode
from .input_table import InputTable, TableDef


class Partition(H5NastranNode):
    def __init__(self, h5n, input):
        self._h5n = h5n
        self._input = input

        self.aesurf = AESURF(self._h5n, self)
        self.set1 = SET1(self._h5n, self)
        self.uset = USET(self._h5n, self)

    def path(self):
        return self._input.path() + ['PARTITION']


########################################################################################################################


class AESURF(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PARTITION/AESURF')

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        result = {
            'IDENTITY': {
                'ID': [], 'LABEL': [], 'CID1': [], 'ALID1': [], 'CID2': [], 'ALID2': [], 'EFF': [], 'LDW': [],
                'CREFC': [], 'CREFS': [], 'S_PLLIM': [], 'PLLIM': [], 'S_PULIM': [], 'PULIM': [],
                'S_HMLLIM': [], 'HMLLIM': [], 'S_HMULIM': [], 'HMULIM': [],'S_TQLLIM': [], 'TQLLIM': [],
                'S_TQULIM': [], 'TQULIM': [], 'DOMAIN_ID': []
            }
        }

        i = result['IDENTITY']
        id_ = i['ID']
        label = i['LABEL']
        cid1 = i['CID1']
        alid1 = i['ALID1']
        cid2 = i['CID2']
        alid2 = i['ALID2']
        eff = i['EFF']
        ldw = i['LDW']
        crefc = i['CREFC']
        crefs = i['CREFS']
        s_pllim = i['S_PLLIM']
        pllim = i['PLLIM']
        s_pulim = i['S_PULIM']
        pulim = i['PULIM']
        s_hmllim = i['S_HMLLIM']
        hmllim = i['HMLLIM']
        s_hmulim = i['S_HMULIM']
        hmulim = i['HMULIM']
        s_tqllim = i['S_TQLLIM']
        tqllim = i['TQLLIM']
        s_tqulim = i['S_TQULIM']
        tqulim = i['TQULIM']

        def _get_vals(val):
            if val is None:
                return 0, 0.
            else:
                return 1, val

        # TODO: AESURF - verify ldw is correct
        def _get_ldw(val):
            if val in ('', None):
                val = 'LDW'

            if val == 'LDW':
                return 1
            elif val == 'NOLDW':
                return 0
            else:
                raise ValueError(val)

        for card_id in card_ids:
            card = cards[card_id]

            id_.append(card.aesid)
            label.append(card.label)
            cid1.append(card.cid1)
            alid1.append(card.alid1)
            cid2.append(card.cid2 if card.cid2 is not None else Defaults.default_int)
            alid2.append(card.alid2 if card.alid2 is not None else Defaults.default_int)
            eff.append(card.eff)
            ldw.append(_get_ldw(card.ldw))
            crefc.append(card.crefc)
            crefs.append(card.crefs)

            _s, _v = _get_vals(card.pllim)
            s_pllim.append(_s)
            pllim.append(_v)

            _s, _v = _get_vals(card.pulim)
            s_pulim.append(_s)
            pulim.append(_v)

            _s, _v = _get_vals(card.hmllim)
            s_hmllim.append(_s)
            hmllim.append(_v)

            _s, _v = _get_vals(card.hmulim)
            s_hmulim.append(_s)
            hmulim.append(_v)

            _s, _v = _get_vals(card.tqllim)
            s_tqllim.append(_s)
            tqllim.append(_v)

            _s, _v = _get_vals(card.tqulim)
            s_tqulim.append(_s)
            tqulim.append(_v)

        return result


########################################################################################################################


class SET1(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PARTITION/SET1/IDENTITY')

    def from_bdf(self, cards):
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


class USET(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PARTITION/USET')

    def from_bdf(self, cards):
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