from __future__ import print_function, absolute_import
from six import iteritems, itervalues, iterkeys
from six.moves import range

import tables
import numpy as np

from .card_table import CardTable, TableDef
from ..data_helper import DataHelper


class Parameter(object):
    def __init__(self, h5n, input):
        self._h5n = h5n
        self._input = input

        self.nlparm = NLPARM(self._h5n, self)

    def path(self):
        return self._input.path() + ['PARAMETER']

    def read(self):
        for key, item in iteritems(self.__dict__):
            if key.startswith('_'):
                continue
            try:
                item.read()
            except AttributeError:
                pass


########################################################################################################################


class NLPARM(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PARAMETER/NLPARM')

    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())

        data = np.empty(len(card_ids), dtype=cls.table_def.dtype)

        sid = data['SID']
        ninc = data['NINC']
        dt = data['DT']
        kmethod = data['KMETHOD']
        kstep = data['KSTEP']
        maxiter = data['MAXITER']
        conv = data['CONV']
        intout = data['INTOUT']
        epsu = data['EPSU']
        epsp = data['EPSP']
        epsw = data['EPSW']
        maxdiv = data['MAXDIV']
        maxqn = data['MAXQN']
        maxls = data['MAXLS']
        fstress = data['FSTRESS']
        lstol = data['LSTOL']
        maxbis = data['MAXBIS']
        maxr = data['MAXR']
        rtolb = data['RTOLB']
        miniter = data['MINITER']

        # TODO: check that kmethod is correct
        _kmethod = {'AUTO': 0, 'ITER': 1, 'SEMI': 2, 'FNT': 3, 'PFNT': 4, None: 0, '': 0}

        # TODO: check that conv is correct; what are the allowed combinations?
        _conv = {}

        # TODO: check that intout is correct
        _intout = {'YES': 1, 'NO': 0, 'ALL': 3, None: DataHelper.default_int, '': DataHelper.default_int}

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            sid[i] = card.nlparm_id
            ninc[i] = card.ninc if card.ninc is not None else DataHelper.default_int
            dt[i] = card.dt
            kmethod[i] = _kmethod[card.kmethod]
            kstep[i] = card.kstep
            maxiter[i] = card.max_iter
            conv[i] = DataHelper.unknown_int
            if isinstance(card.int_out, int):
                intout[i] = card.int_out
            else:
                intout[i] = _intout[card.int_out]
            epsu[i] = card.eps_u
            epsp[i] = card.eps_p
            epsw[i] = card.eps_w
            maxdiv[i] = card.max_div
            maxqn[i] = card.max_qn
            maxls[i] = card.max_ls
            fstress[i] = card.fstress
            lstol[i] = card.ls_tol
            maxbis[i] = card.max_bisect
            maxr[i] = card.max_r
            rtolb[i] = card.rtol_b
            try:
                miniter[i] = card.min_iter
            except AttributeError:
                miniter[i] = DataHelper.default_int

        result = {'IDENTITY': data}

        return result
