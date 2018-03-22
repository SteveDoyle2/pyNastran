from __future__ import print_function, absolute_import
from six import iteritems, itervalues, iterkeys
from six.moves import range

import tables
import numpy as np

from .input_table import InputTable, TableDef
from ..data_helper import DataHelper
from ..h5nastrannode import H5NastranNode


class Design(H5NastranNode):
    def __init__(self, h5n, input):
        self._h5n = h5n
        self._input = input

        self.desvar = DESVAR(self._h5n, self)
        self.dvprel1 = DVPREL1(self._h5n, self)

    def path(self):
        return self._input.path() + ['DESIGN']


########################################################################################################################


class DESVAR(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/DESIGN/DESVAR')

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

        def _get_val(val, default):
            if val in ('', None):
                val = default
            return val

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
            dvid[i] = _get_val(card.ddval, DataHelper.default_int)  # TODO: DESVAR use ddval for dvid?

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class DVPREL1(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/DESIGN/DVPREL1/IDENTITY',
                                rename={'RELATION_POS': 'START', 'RELATION_LEN': 'LEN'})

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
                _fid = DataHelper.default_int
                _pname = DataHelper.default_str

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
