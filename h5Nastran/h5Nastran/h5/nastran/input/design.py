from __future__ import print_function, absolute_import

import numpy as np
from collections import OrderedDict
from six import iteritems

from h5Nastran.defaults import Defaults
from h5Nastran.h5nastrannode import H5NastranNode
from .input_table import InputTable, TableDef


class Design(H5NastranNode):
    def __init__(self, h5n, input):
        self._h5n = h5n
        self._input = input

        self.deqatn = DEQATN(self._h5n, self)
        self.desvar = DESVAR(self._h5n, self)
        self.dvprel1 = DVPREL1(self._h5n, self)

    def path(self):
        return self._input.path() + ['DESIGN']

    def to_bdf(self, bdf):
        for key, item in iteritems(self.__dict__):
            if key.startswith('_'):
                continue
            try:
                item.to_bdf(bdf)
            except NotImplementedError:
                pass


########################################################################################################################


class DEQATN_SPEC(object):
    name = 'DEQATN'
    path = '/NASTRAN/INPUT/DESIGN'
    dtype = [('ID', '<i8', ()), ('EQUATION', 'S256', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


class DEQATN(InputTable):
    table_def = TableDef.create(DEQATN_SPEC)

    def to_bdf(self, bdf):
        add_card = bdf.add_deqatn

        data = self.identity

        id_ = data['ID']
        eqs = data['EQUATION']

        _eqs = OrderedDict()
        for i in range(id_.size):
            _id = id_[i]
            _eq = eqs[i].decode()
            try:
                _eqs[_id].append(_eq)
            except KeyError:
                _eqs[_id] = [_eq]

        for eq_id, eqs in iteritems(_eqs):
            add_card(eq_id, eqs)

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        result = {
            'IDENTITY': {'ID': [], 'EQUATION': []}
        }

        i = result['IDENTITY']

        id_ = i['ID']
        eqs = i['EQUATION']

        for card_id in card_ids:
            card = cards[card_id]

            for eq in card.eqs:
                id_.append(card.equation_id)
                eqs.append(eq)

        return result


########################################################################################################################


class DESVAR(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/DESIGN/DESVAR')

    def to_bdf(self, bdf):
        add_card = bdf.add_desvar

        data = self.identity

        id_ = data['ID']
        label = data['LABEL']
        xinit = data['XINIT']
        xlb = data['XLB']
        xub = data['XUB']
        delx = data['DELX']
        dvid = data['DVID']

        to_dvid = self._h5n.defaults.to_value_int

        for i in range(id_.size):
            add_card(id_[i], label[i].decode(), xinit[i], xlb[i], xub[i], delx[i], dvid[i])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        id_ = data['ID']
        label = data['LABEL']
        xinit = data['XINIT']
        xlb = data['XLB']
        xub = data['XUB']
        delx = data['DELX']
        dvid = data['DVID']

        get_dvid = self._h5n.defaults.get_value_int

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            id_[i] = card.desvar_id
            label[i] = card.label
            xinit[i] = card.xinit
            xlb[i] = card.xlb
            xub[i] = card.xub
            delx[i] = card.delx
            dvid[i] = get_dvid(card.ddval)  # TODO: DESVAR use ddval for dvid?

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class DVPREL1(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/DESIGN/DVPREL1/IDENTITY',
                                rename={'RELATION_POS': 'START', 'RELATION_LEN': 'LEN'})

    def to_bdf(self, bdf):
        add_card = bdf.add_dvprel1

        identity = self.identity
        relation = self.relation

        dvid = relation['DVID']
        coef = relation['COEF']

        id_ = identity['ID']
        type_ = identity['TYPE']
        pid = identity['PID']
        fid = identity['FID']
        pmin = identity['PMIN']
        pmax = identity['PMAX']
        c0 = identity['C0']
        pname = identity['PNAME']
        start = identity['START']
        len_ = identity['LEN']

        defaults = self._h5n.defaults
        to_value_str = defaults.to_value_str
        to_value_int = defaults.to_value_int
        to_value_double = defaults.to_value_double

        for i in range(id_.size):
            oid = int(id_[i])
            prop_type = to_value_str(type_[i].decode())
            pid_ = int(pid[i])
            _pname = to_value_str(pname[i].decode())
            _fid = to_value_int(fid[i])
            if _pname is None:
                pname_fid = _fid
            else:
                pname_fid = _pname
            j1 = start[i]
            j2 = j1 + len_[i]
            dvids = dvid[j1:j2].tolist()
            coeffs = coef[j1:j2].tolist()
            p_min = to_value_double(pmin[i])
            p_max = to_value_double(pmax[i])
            c0_ = to_value_double(c0[i])

            add_card(oid, prop_type, pid_, pname_fid, dvids, coeffs, p_min, p_max, c0_)

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        
        relation = {'IDENTITY': {'DVID': [], 'COEF': []}}
        result = {'IDENTITY': {'ID': [], 'TYPE': [], 'PID': [], 'FID': [], 'PMIN': [], 'PMAX': [], 'C0': [], 'PNAME': [],
                               'START': [], 'LEN': [], 'DOMAIN_ID': []},
                  'RELATION': relation,
                  '_subtables': ['RELATION']}
        
        rel = relation['IDENTITY']
        dvid = rel['DVID']
        coef = rel['COEF']
        
        identity = result['IDENTITY']
        id_ = identity['ID']
        type_ = identity['TYPE']
        pid = identity['PID']
        fid = identity['FID']
        pmin = identity['PMIN']
        pmax = identity['PMAX']
        c0 = identity['C0']
        pname = identity['PNAME']
        start = identity['START']
        len_ = identity['LEN']

        def _get_val(val, default):
            if val in ('', None):
                val = default
            return val

        defaults = self._h5n.defaults

        _pos = 0
        for card_id in card_ids:
            card = cards[card_id]

            if isinstance(card.pname_fid, int):
                _fid = card.pname_fid
                _pname = ''
            elif isinstance(card.pname_fid, str):
                _fid = 0
                _pname = card.pname_fid
            else:
                _fid = defaults.default_int
                _pname = defaults.default_str

            id_.append(card.oid)
            type_.append(card.prop_type)
            pid.append(card.pid)
            fid.append(_fid)
            pmin.append(card.p_min)
            pmax.append(card.p_max)
            c0.append(card.c0)
            pname.append(_pname)
            start.append(_pos)
            _len = len(card.dvids)
            len_.append(_len)
            _pos += _len
            dvid += card.dvids
            coef += card.coeffs

        return result
