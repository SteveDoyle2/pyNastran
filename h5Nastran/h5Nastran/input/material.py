from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range

import tables
import numpy as np

from .card_table import CardTable, TableDef, TableData


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

    @staticmethod
    def from_bdf(card):
        data = [card.mid, card.e, card.g, card.nu, card.rho, card.a, card.tref, card.ge, card.St, card.Sc, card.Ss,
                card.mcsid]
        return TableData([data])

########################################################################################################################


class MAT4(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/MATERIAL/MAT4')

    @staticmethod
    def from_bdf(card):
        data = [card.mid, card.k, card.cp, card.rho, card.H, card.mu, card.hgen, card.ref_enthalpy, card.tch,
                card.tdelta, card.qlat]
        return TableData([data])

########################################################################################################################

class MAT8(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/MATERIAL/MAT8')

    @staticmethod
    def from_bdf(card):
        data = [card.mid, card.e11, card.e22, card.nu12, card.g12, card.g1z, card.g2z, card.rho, card.a1, card.a2,
                card.tref, card.Xt, card.Xc, card.Yt, card.Yc, card.S, card.ge, card.F12, card.strn]
        return TableData([data])

########################################################################################################################
