from __future__ import print_function, absolute_import

import numpy as np
from six.moves import range

from h5Nastran.defaults import Defaults
from h5Nastran.h5nastrannode import H5NastranNode
from .input_table import InputTable, TableDef


class Load(H5NastranNode):
    def __init__(self, h5n, input):
        self._h5n = h5n
        self._input = input

        self.darea = DAREA(self._h5n, self)
        self.dload = DLOAD(self._h5n, self)
        self.force = FORCE(self._h5n, self)
        self.grav = GRAV(self._h5n, self)
        self.load = LOAD(self._h5n, self)
        self.moment = MOMENT(self._h5n, self)
        self.pload4 = PLOAD4(self._h5n, self)
        self.rload1 = RLOAD1(self._h5n, self)
        self.rload2 = RLOAD2(self._h5n, self)

    def path(self):
        return self._input.path() + ['LOAD']


########################################################################################################################


class DAREA(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/LOAD/DAREA')

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        result = {
            'IDENTITY': {
                'SID': [],
                'P': [],
                'C': [],
                'A': []
            }
        }

        identity = result['IDENTITY']

        sid = identity['SID']
        p = identity['P']
        c = identity['C']
        a = identity['A']

        for card_id in card_ids:
            card = cards[card_id]
            _sid = card.sid

            node_ids = card.node_ids
            components = card.components
            scales = card.scales

            for j in range(len(node_ids)):
                sid.append(_sid)
                p.append(node_ids[j])
                c.append(components[j])
                a.append(scales[j])

        return result


########################################################################################################################


class DLOAD(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/LOAD/DLOAD/IDENTITY')

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        sl = {'IDENTITY': {'SI': [], 'LI': []}}

        result = {
            'IDENTITY': {
                'SID': [], 'S': [], 'SL_LEN': [], 'SL_POS': [], 'DOMAIN_ID': []
            },
            'SL': sl,
            '_subtables': ['SL']
        }

        si = sl['IDENTITY']['SI']
        li= sl['IDENTITY']['LI']

        identity = result['IDENTITY']
        sid = identity['SID']
        s = identity['S']
        sl_len = identity['SL_LEN']
        sl_pos = identity['SL_POS']

        _pos = 0
        for card_id in card_ids:
            card_list = cards[card_id]

            for card in card_list:
                sid.append(card.sid)
                s.append(card.scale)
                sl_pos.append(_pos)
                _len = len(card.scale_factors)
                sl_len.append(_len)
                _pos += _len
                si += card.scale_factors
                li += card.load_ids

        return result


########################################################################################################################


class FORCE(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/LOAD/FORCE')

    def from_bdf(self, cards):
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
            card_list = sorted(cards[card_id], key=lambda x: x.node)

            for card in card_list:
                sid.append(card.sid)
                g.append(card.node)
                cid.append(card.cid)
                f.append(card.mag)
                n.append(card.xyz)

        identity['N'] = np.array(n)

        return result


########################################################################################################################


class GRAV(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/LOAD/GRAV')

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        result = {
            'IDENTITY': {
                'SID': [], 'CID': [], 'A': [], 'N': [], 'MB': [], 'DOMAIN_ID': []
            }
        }

        identity = result['IDENTITY']
        sid = identity['SID']
        cid = identity['CID']
        a = identity['A']
        n = identity['N']
        mb = identity['MB']

        for card_id in card_ids:
            card_list = cards[card_id]

            for card in card_list:
                sid.append(card.sid)
                cid.append(card.cid)
                a.append(card.scale)
                n.append(card.N)
                mb.append(card.mb)

        identity['N'] = np.array(n)

        return result


########################################################################################################################

class LOAD(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/LOAD/LOAD/IDENTITY')

    def from_bdf(self, cards):
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


class MOMENT(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/LOAD/MOMENT')
    
    def from_bdf(self, cards):
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
            card_list = sorted(cards[card_id], key=lambda x: x.node)

            for card in card_list:
                sid.append(card.sid)
                g.append(card.node)
                cid.append(card.cid)
                m.append(card.mag)
                n.append(card.xyz)

        identity['N'] = np.array(n)

        return result


########################################################################################################################


class PLOAD4(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/LOAD/PLOAD4')

    def from_bdf(self, cards):
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
                    _g1 = Defaults.default_int

                _g34 = card.g34
                if _g34 is None:
                    _g34 = Defaults.default_int

                g1.extend([_g1] * eid_len)
                g34.extend([_g34] * eid_len)

                cid.extend([card.cid] * eid_len)
                n.extend([card.nvector] * eid_len)
                sorl.extend([card.surf_or_line] * eid_len)
                ldir.extend([card.line_load_dir] * eid_len)

        identity['P'] = np.array(p)

        return data


########################################################################################################################


class RLOAD1(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/LOAD/RLOAD1')

    def from_bdf(self, cards):
        data = {
            'IDENTITY': {
                'SID': [],
                'DAREA': [],
                'DELAY': [],
                'DPHASE': [],
                'TC': [],
                'TD': [],
                'TYPE': [],
                'T': [],
                'PH': [],
                'RC': [],
                'RD': [],
                'DOMAIN_ID': []
            }
        }

        identity = data['IDENTITY']
        sid = identity['SID']
        darea = identity['DAREA']
        delay = identity['DELAY']
        dphase = identity['DPHASE']
        tc = identity['TC']
        td = identity['TD']
        type_ = identity['TYPE']
        t = identity['T']
        ph = identity['PH']
        rc = identity['RC']
        rd = identity['RD']

        card_ids = sorted(cards.keys())

        def _get_type(val):
            if isinstance(val, int):
                return val
            if val in {None, '', 'L', 'LO', 'LOA', 'LOAD'}:
                return 0
            elif val in {'D', 'DI', 'DIS', 'DISP'}:
                return 1
            elif val in {'V', 'VE', 'VEL', 'VELO'}:
                return 2
            elif val in {'A', 'AC', 'ACC', 'ACCE'}:
                return 3
            else:
                raise ValueError

        def _get_vals(val):
            if val is None:
                return 0, 0.
            elif isinstance(val, int):
                return val, 0.
            elif isinstance(val, float):
                return 0, val
            else:
                raise ValueError

        for card_id in card_ids:
            card_list = cards[card_id]

            for card in card_list:
                sid.append(card.sid)
                darea.append(card.excite_id)
                _delay, _t = _get_vals(card.delay)
                _dphase, _ph = _get_vals(card.dphase)
                delay.append(_delay)
                dphase.append(_dphase)
                _tc, _rc = _get_vals(card.tc)
                _td, _rd = _get_vals(card.td)
                tc.append(_tc)
                td.append(_td)
                rc.append(_rc)
                rd.append(_rd)
                type_.append(_get_type(card.Type))
                t.append(_t)
                ph.append(_ph)

        return data


########################################################################################################################


class RLOAD2(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/LOAD/RLOAD2')

    def from_bdf(self, cards):
        data = {
            'IDENTITY': {
                'SID': [],
                'DAREA': [],
                'DELAY': [],
                'DPHASE': [],
                'TB': [],
                'TP': [],
                'TYPE': [],
                'T': [],
                'PH': [],
                'RB': [],
                'RP': [],
                'DOMAIN_ID': []
            }
        }

        identity = data['IDENTITY']
        sid = identity['SID']
        darea = identity['DAREA']
        delay = identity['DELAY']
        dphase = identity['DPHASE']
        tb = identity['TB']
        tp = identity['TP']
        type_ = identity['TYPE']
        t = identity['T']
        ph = identity['PH']
        rb = identity['RB']
        rp = identity['RP']

        card_ids = sorted(cards.keys())

        def _get_type(val):
            if isinstance(val, int):
                return val
            if val in {None, '', 'L', 'LO', 'LOA', 'LOAD'}:
                return 0
            elif val in {'D', 'DI', 'DIS', 'DISP'}:
                return 1
            elif val in {'V', 'VE', 'VEL', 'VELO'}:
                return 2
            elif val in {'A', 'AC', 'ACC', 'ACCE'}:
                return 3
            else:
                raise ValueError

        def _get_vals(val):
            if val is None:
                return 0, 0.
            elif isinstance(val, int):
                return val, 0.
            elif isinstance(val, float):
                return 0, val
            else:
                raise ValueError

        for card_id in card_ids:
            card_list = cards[card_id]

            for card in card_list:
                sid.append(card.sid)
                darea.append(card.excite_id)
                _delay, _t = _get_vals(card.delay)
                _dphase, _ph = _get_vals(card.dphase)
                delay.append(_delay)
                dphase.append(_dphase)
                _tb, _rb = _get_vals(card.tb)
                _tp, _rp = _get_vals(card.tp)
                tb.append(_tb)
                tp.append(_tp)
                rb.append(_rb)
                rp.append(_rp)
                type_.append(_get_type(card.Type))
                t.append(_t)
                ph.append(_ph)

        return data
