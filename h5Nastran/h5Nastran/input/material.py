from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range

import tables
import numpy as np

from .card_table import CardTable, TableDef


class Material(object):
    def __init__(self, h5n, input):
        self._h5n = h5n
        self._input = input

        self.mat1 = MAT1(self._h5n, self)
        self.mat4 = MAT4(self._h5n, self)
        self.mat8 = MAT8(self._h5n, self)

    def path(self):
        return self._input.path() + ['MATERIAL']

    def read(self):
        for key, item in iteritems(self.__dict__):
            if key.startswith('_'):
                continue
            try:
                item.read()
            except AttributeError:
                pass

########################################################################################################################


class MAT1(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/MATERIAL/MAT1')

    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=cls.table_def.dtype)

        mid = data['MID']
        e = data['E']
        g = data['G']
        nu = data['NU']
        rho = data['RHO']
        a = data['A']
        tref = data['TREF']
        ge = data['GE']
        st = data['ST']
        sc = data['SC']
        ss = data['SS']
        mcsid = data['MCSID']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            mid[i] = card.mid
            e[i] = card.e
            g[i] = card.g
            nu[i] = card.nu
            rho[i] = card.rho
            a[i] = card.a
            tref[i] = card.tref
            ge[i] = card.ge
            st[i] = card.St
            sc[i] = card.Sc
            ss[i] = card.Ss
            mcsid[i] = card.mcsid

        result = {'IDENTITY': data}

        return result


########################################################################################################################


class MAT4(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/MATERIAL/MAT4')

    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=cls.table_def.dtype)

        mid = data['MID']
        k = data['K']
        cp = data['CP']
        rho = data['RHO']
        h = data['H']
        mu = data['MU']
        hgen = data['HGEN']
        refenth = data['REFENTH']
        tch = data['TCH']
        tdelta = data['TDELTA']
        qlat = data['QLAT']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            mid[i] = card.mid
            k[i] = card.k
            cp[i] = card.cp
            rho[i] = card.rho
            h[i] = card.H
            mu[i] = card.mu
            hgen[i] = card.hgen
            refenth[i] = card.ref_enthalpy
            tch[i] = card.tch
            tdelta[i] = card.tdelta
            qlat[i] = card.qlat

        result = {'IDENTITY': data}

        return result

########################################################################################################################

class MAT8(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/MATERIAL/MAT8')

    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=cls.table_def.dtype)

        mid = data['MID']
        e1 = data['E1']
        e2 = data['E2']
        nu12 = data['NU12']
        g1z = data['G1Z']
        g2z = data['G2Z']
        rho = data['RHO']
        a1 = data['A1']
        a2 = data['A2']
        tref = data['TREF']
        xt = data['XT']
        xc = data['XC']
        yt = data['YT']
        yc = data['YC']
        s = data['S']
        ge = data['GE']
        f12 = data['F12']
        strn = data['STRN']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            mid[i] = card.mid
            e1[i] = card.e11
            e2[i] = card.e22
            nu12[i] = card.nu12
            g1z[i] = card.g1z
            g2z[i] = card.g2z
            rho[i] = card.rho
            a1[i] = card.a1
            a2[i] = card.a2
            tref[i] = card.tref
            xt[i] = card.Xt
            xc[i] = card.Xc
            yt[i] = card.Yt
            yc[i] = card.Yc
            s[i] = card.S
            ge[i] = card.ge
            f12[i] = card.F12
            strn[i] = card.strn

        result = {'IDENTITY': data}

        return result

########################################################################################################################
