from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range

import tables
import numpy as np

from .card_table import CardTable, TableDef


class Constraint(object):
    def __init__(self, h5n, input):
        self._h5n = h5n
        self._input = input

        self.aelink = AELINK(self._h5n, self)
        self.csschd = CSSCHD(self._h5n, self)
        self.cysup = CYSUP(self._h5n, self)
        self.deform = DEFORM(self._h5n, self)
        self.grdset = GRDSET(self._h5n, self)
        self.mpc = MPC(self._h5n, self)
        self.mpcadd = MPCADD(self._h5n, self)
        self.mpcax = MPCAX(self._h5n, self)
        self.mpcd = MPCD(self._h5n, self)
        self.mpcy = MPCY(self._h5n, self)
        self.rspline = RSPLINE(self._h5n, self)
        self.sesup = SESUP(self._h5n, self)
        self.spblnd1 = SPBLND1(self._h5n, self)
        self.spblnd2 = SPBLND2(self._h5n, self)
        self.spc = SPC(self._h5n, self)
        self.spc1_g = SPC1_G(self._h5n, self)
        self.spc1_thru = SPC1_THRU(self._h5n, self)
        self.spcadd = SPCADD(self._h5n, self)
        self.spcax = SPCAX(self._h5n, self)
        self.spcd = SPCD(self._h5n, self)
        self.spcoff = SPCOFF(self._h5n, self)
        self.spcoff1 = SPCOFF1(self._h5n, self)
        self.spcr = SPCR(self._h5n, self)
        self.spline1 = SPLINE1(self._h5n, self)
        self.spline2 = SPLINE2(self._h5n, self)
        self.spline3 = SPLINE3(self._h5n, self)
        self.spline4 = SPLINE4(self._h5n, self)
        self.spline5 = SPLINE5(self._h5n, self)
        self.spline6 = SPLINE6(self._h5n, self)
        self.spline7 = SPLINE7(self._h5n, self)
        self.splinex = SPLINEX(self._h5n, self)
        self.splinrb = SPLINRB(self._h5n, self)
        self.supax = SUPAX(self._h5n, self)
        self.suport = SUPORT(self._h5n, self)
        self.suport1 = SUPORT1(self._h5n, self)
        self.temp = TEMP(self._h5n, self)
        self.tempax = TEMPAX(self._h5n, self)
        self.tempb3 = TEMPB3(self._h5n, self)
        self.tempbc = TEMPBC(self._h5n, self)
        self.tempd = TEMPD(self._h5n, self)
        self.tempn1 = TEMPN1(self._h5n, self)
        self.tempp1 = TEMPP1(self._h5n, self)
        self.tempp2 = TEMPP2(self._h5n, self)
        self.tempp3 = TEMPP3(self._h5n, self)
        self.temprb = TEMPRB(self._h5n, self)
        self.trim = TRIM(self._h5n, self)
        self.trim2 = TRIM2(self._h5n, self)
        self.uxvec = UXVEC(self._h5n, self)

    def path(self):
        return self._input.path() + ['CONSTRAINT']

    def read(self):
        for key, item in iteritems(self.__dict__):
            if key.startswith('_'):
                continue
            try:
                item.read()
            except AttributeError:
                pass

########################################################################################################################


class AELINK(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/AELINK/IDENTITY')

########################################################################################################################


class CSSCHD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/CSSCHD')

########################################################################################################################


class CYSUP(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/CYSUP')

########################################################################################################################


class DEFORM(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/DEFORM')

########################################################################################################################


class GRDSET(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/GRDSET')

########################################################################################################################


class MPC(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/MPC/IDENTITY')

########################################################################################################################


class MPCADD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/MPCADD/IDENTITY')

########################################################################################################################


class MPCAX(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/MPCAX/IDENTITY')

########################################################################################################################


class MPCD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/MPCD')

########################################################################################################################


class MPCY(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/MPCY/IDENTITY')

########################################################################################################################


class RSPLINE(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/RSPLINE/IDENTITY')

########################################################################################################################


class SESUP(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SESUP')

########################################################################################################################


class SPBLND1(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPBLND1')

########################################################################################################################


class SPBLND2(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPBLND2')

########################################################################################################################


class SPC(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPC')

########################################################################################################################


class SPC1_G(CardTable):
    card_id = 'SPC1'
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPC1/SPC1_G/IDENTITY')
    
    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())

        g = {
            'IDENTITY': {'ID': []}
        }

        result = {
            'IDENTITY': {'SID': [], 'C': [], 'G_POS': [], 'G_LEN': [], 'DOMAIN_ID': []},
            'G': g,
            '_subtables': ['G']
        }

        identity = result['IDENTITY']

        sid = identity['SID']
        c = identity['C']
        g_pos = identity['G_POS']
        g_len = identity['G_LEN']
        id = g['IDENTITY']['ID']

        _pos = 0
        for card_id in card_ids:
            card_list = sorted(cards[card_id], key=lambda x: x.conid)
            for card in card_list:
                nids = []
                for nid in card.node_ids:
                    if nid is None:
                        continue
                    nids.append(nid)

                sid.append(card.conid)
                c.append(card.components)
                g_pos.append(_pos)
                _g_len = len(nids)
                _pos += _g_len
                g_len.append(_g_len)
                id += nids

        return result


########################################################################################################################


class SPC1_THRU(CardTable):
    card_id = None  # don't register this
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPC1/SPC1_THRU')

########################################################################################################################


class SPCADD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPCADD/IDENTITY')

########################################################################################################################


class SPCAX(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPCAX')

########################################################################################################################


class SPCD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPCD')

########################################################################################################################


class SPCOFF(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPCOFF')


########################################################################################################################


class SPCOFF1(CardTable):
    table_def = TableDef.create(
        '/NASTRAN/INPUT/CONSTRAINT/SPCOFF1/IDENTITY',
        rename={'GIS_POS': 'G_POS', 'GIS_LEN': 'G_LEN'}
    )

########################################################################################################################


class SPCR(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPCR')

########################################################################################################################


class SPLINE1(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPLINE1')

########################################################################################################################


class SPLINE2(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPLINE2')

########################################################################################################################


class SPLINE3(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPLINE3/IDENTITY')

########################################################################################################################


class SPLINE4(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPLINE4')

########################################################################################################################


class SPLINE5(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPLINE5')

########################################################################################################################


class SPLINE6(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPLINE6')

########################################################################################################################


class SPLINE7(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPLINE7')

########################################################################################################################


class SPLINEX(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPLINEX')

########################################################################################################################


class SPLINRB(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPLINRB/IDENTITY')

########################################################################################################################


class SUPAX(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SUPAX')

########################################################################################################################


class SUPORT(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SUPORT')

########################################################################################################################


class SUPORT1(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SUPORT1/IDENTITY')

########################################################################################################################


class TEMP(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/TEMP')

########################################################################################################################


class TEMPAX(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/TEMPAX')

########################################################################################################################


class TEMPB3(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/TEMPB3')

########################################################################################################################


class TEMPBC(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/TEMPBC')

########################################################################################################################


class TEMPD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/TEMPD')

########################################################################################################################


class TEMPN1(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/TEMPN1')

########################################################################################################################


class TEMPP1(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/TEMPP1')

########################################################################################################################


class TEMPP2(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/TEMPP2')

########################################################################################################################


class TEMPP3(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/TEMPP3')


########################################################################################################################


class TEMPRB(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/TEMPRB')

########################################################################################################################


class TRIM(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/TRIM/IDENTITY')

########################################################################################################################


class TRIM2(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/TRIM2/IDENTITY')

########################################################################################################################


class UXVEC(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/UXVEC/IDENTITY')

########################################################################################################################
