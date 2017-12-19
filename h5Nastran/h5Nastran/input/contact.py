from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range

import tables
import numpy as np

from .card_table import CardTable, TableDef


class Contact(object):
    def __init__(self, h5n, input):
        self._h5n = h5n
        self._input = input

        self.bcbdprp = BCBDPRP(self._h5n, self)
        self.bcbmrad = BCBMRAD(self._h5n, self)
        # self.bcbody = BCBODY(self._h5n, self)
        self.bcbody1 = BCBODY1(self._h5n, self)
        self.bcbzier = BCBZIER(self._h5n, self)
        self.bcmove = BCMOVE(self._h5n, self)
        self.bcnurb2 = BCNURB2(self._h5n, self)
        self.bcnurbs = BCNURBS(self._h5n, self)
        self.bconect = BCONECT(self._h5n, self)
        self.bconprg = BCONPRG(self._h5n, self)
        self.bconprp = BCONPRP(self._h5n, self)
        self.bconuds = BCONUDS(self._h5n, self)
        self.bcpara = BCPARA(self._h5n, self)
        self.bcpatch = BCPATCH(self._h5n, self)
        self.bcpflg = BCPFLG(self._h5n, self)
        self.bcprop = BCPROP(self._h5n, self)
        self.bcrgsrf = BCRGSRF(self._h5n, self)
        self.bcrigid = BCRIGID(self._h5n, self)
        self.bcsap = BCSCAP(self._h5n, self)
        self.bctabl1 = BCTABL1(self._h5n, self)
        self.bctable = BCTABLE(self._h5n, self)
        self.bctrim = BCTRIM(self._h5n, self)
        self.blseg = BLSEG(self._h5n, self)
        self.boutput = BOUTPUT(self._h5n, self)
        self.bsqueal = BSQUEAL(self._h5n, self)
        self.bsurf = BSURF(self._h5n, self)
        self.bsurf_old = BSURF_OLD(self._h5n, self)
        self.prjcon = PRJCON(self._h5n, self)
        self.unglue = UNGLUE(self._h5n, self)

    def path(self):
        return self._input.path() + ['CONTACT']


########################################################################################################################


class BCBDPRP(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONTACT/BCBDPRP')

########################################################################################################################


class BCBMRAD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONTACT/BCBMRAD/IDENTITY')

########################################################################################################################

# TODO: BCBODY - xml doesn't conform
# class BCBODY(CardTable):
#     table_def = TableDef.create('/NASTRAN/INPUT/CONTACT/BCBODY/IDENTITY')

########################################################################################################################


class BCBODY1(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONTACT/BCBODY1')

########################################################################################################################


class BCBZIER(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONTACT/BCBZIER/IDENTITY')

########################################################################################################################


class BCHANGE(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONTACT/BCHANGE/IDENTITY')

########################################################################################################################


class BCMOVE(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONTACT/BCMOVE')

########################################################################################################################


class BCNURB2(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONTACT/BCNURB2/IDENTITY')

########################################################################################################################


class BCNURBS(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONTACT/BCNURBS/IDENTITY')

########################################################################################################################


########################################################################################################################


########################################################################################################################


########################################################################################################################


########################################################################################################################


########################################################################################################################


########################################################################################################################


########################################################################################################################


########################################################################################################################


########################################################################################################################


########################################################################################################################
