from __future__ import print_function, absolute_import

import numpy as np

from h5Nastran.defaults import Defaults
from h5Nastran.h5nastrannode import H5NastranNode
from .input_table import InputTable, TableDef


class Parameter(H5NastranNode):
    def __init__(self, h5n, input):
        self._h5n = h5n
        self._input = input

        self.flfact = FLFACT(self._h5n, self)
        self.flutter = FLUTTER(self._h5n, self)
        self.nlparm = NLPARM(self._h5n, self)

    def path(self):
        return self._input.path() + ['PARAMETER']


########################################################################################################################


class FLFACT(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PARAMETER/FLFACT/IDENTITY')

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        fis = {'IDENTITY': {'FI': []}}

        result = {
            'IDENTITY': {'SID': [], 'FIS_POS': [], 'FIS_LEN': [], 'DOMAIN_ID': []},
            'FIS': fis,
            '_subtables': ['FIS']
        }

        fi = fis['IDENTITY']['FI']
        identity = result['IDENTITY']
        sid = identity['SID']
        fis_pos = identity['FIS_POS']
        fis_len = identity['FIS_LEN']

        pos = 0
        for card_id in card_ids:
            card = cards[card_id]

            sid.append(card.sid)
            fis_pos.append(pos)
            _len = len(card.factors)
            fis_len.append(_len)
            pos += _len
            fi.extend(list(card.factors))

        return result


########################################################################################################################


class FLUTTER(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PARAMETER/FLUTTER')

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        # flutter has no domain_id in msc spec
        result = {
            'IDENTITY': {'SID': [], 'METHOD': [], 'DENS': [], 'MACH': [], 'RFREQ': [], 'IMETH': [], 'SFLG': [],
            'NVALUE': [], 'OMAX': [], 'EPS': []}
        }

        identity = result['IDENTITY']
        sid = identity['SID']
        method = identity['METHOD']
        dens = identity['DENS']
        mach = identity['MACH']
        rfreq = identity['RFREQ']
        imeth = identity['IMETH']
        sflg = identity['SFLG']
        nvalue = identity['NVALUE']
        omax = identity['OMAX']
        eps = identity['EPS']

        for card_id in card_ids:
            card = cards[card_id]

            sid.append(card.sid)
            method.append(card.method)
            dens.append(card.density)
            mach.append(card.mach)
            rfreq.append(card.reduced_freq_velocity)
            imeth.append(card.imethod)

            nvalue_, omax_ = card.nvalue, card.omax

            if nvalue_ is None:
                sflg_ = 1
                nvalue_ = 0
            elif omax_ is None:
                sflg_ = 0
                omax_ = 0.
            else:
                raise ValueError

            sflg.append(sflg_)
            nvalue.append(nvalue_)
            omax.append(omax_)
            eps.append(card.epsilon)

        return result


########################################################################################################################


class NLPARM(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PARAMETER/NLPARM')

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

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
        _intout = {'YES': 1, 'NO': 0, 'ALL': 3, None: Defaults.default_int, '': Defaults.default_int}

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            sid[i] = card.nlparm_id
            ninc[i] = card.ninc if card.ninc is not None else Defaults.default_int
            dt[i] = card.dt
            kmethod[i] = _kmethod[card.kmethod]
            kstep[i] = card.kstep
            maxiter[i] = card.max_iter
            conv[i] = Defaults.unknown_int
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
                miniter[i] = Defaults.default_int

        result = {'IDENTITY': data}

        return result
