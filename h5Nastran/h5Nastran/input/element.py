from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range

import tables
import numpy as np

from .card_table import CardTable, TableDef, TableData


class Element(object):
    def __init__(self, h5n, input):
        self._h5n = h5n
        self._input = input

        self.aequad4 = AEQUAD4(self._h5n, self)
        self.aerod = AEROD(self._h5n, self)
        self.aeroq4 = AEROQ4(self._h5n, self)
        self.aerot3 = AEROT3(self._h5n, self)
        self.aetria3 = AETRIA3(self._h5n, self)
        self.beamaero = BEAMAERO(self._h5n, self)
        self.bolt = BOLT(self._h5n, self)
        self.caabsf = CAABSF(self._h5n, self)
        self.cacinf3 = CACINF3(self._h5n, self)
        self.cacinf4 = CACINF4(self._h5n, self)
        self.caero1 = CAERO1(self._h5n, self)
        self.caero2 = CAERO2(self._h5n, self)
        self.caero3 = CAERO3(self._h5n, self)
        self.caero4 = CAERO4(self._h5n, self)
        self.caero5 = CAERO5(self._h5n, self)
        self.caxif2 = CAXIF2(self._h5n, self)
        self.caxif3 = CAXIF3(self._h5n, self)
        self.caxif4 = CAXIF4(self._h5n, self)
        self.caxisym = CAXISYM(self._h5n, self)
        self.cbar = CBAR(self._h5n, self)
        self.cbeam = CBEAM(self._h5n, self)
        self.cbeam3 = CBEAM3(self._h5n, self)
        self.cbend = CBEND(self._h5n, self)
        self.cbush = CBUSH(self._h5n, self)
        self.cbush1d = CBUSH1D(self._h5n, self)
        self.cbush2d = CBUSH2D(self._h5n, self)
        self.cconeax = CCONEAX(self._h5n, self)
        self.cdamp1 = CDAMP1(self._h5n, self)
        self.cdamp2 = CDAMP2(self._h5n, self)
        self.cdamp3 = CDAMP3(self._h5n, self)
        self.cdamp4 = CDAMP4(self._h5n, self)
        self.cdamp5 = CDAMP5(self._h5n, self)
        self.celas1 = CELAS1(self._h5n, self)
        self.celas2 = CELAS2(self._h5n, self)
        self.celas3 = CELAS3(self._h5n, self)
        self.celas4 = CELAS4(self._h5n, self)
        self.cfast = CFAST(self._h5n, self)
        self.cfastp = CFASTP(self._h5n, self)
        self.cfluid2 = CFLUID2(self._h5n, self)
        self.cfluid3 = CFLUID3(self._h5n, self)
        self.cfluid4 = CFLUID4(self._h5n, self)
        self.cgap = CGAP(self._h5n, self)
        self.chacab = CHACAB(self._h5n, self)
        self.chacbr = CHACBR(self._h5n, self)
        self.chbdye = CHBDYE(self._h5n, self)
        self.chbdyg = CHBDYG(self._h5n, self)
        self.chbdyp = CHBDYP(self._h5n, self)
        self.chexa = CHEXA(self._h5n, self)
        self.chexal = CHEXAL(self._h5n, self)
        self.chexp = CHEXP(self._h5n, self)
        self.cifhex = CIFHEX(self._h5n, self)
        self.cifpent = CIFPENT(self._h5n, self)
        self.cifqdx = CIFQDX(self._h5n, self)
        self.cifquad = CIFQUAD(self._h5n, self)
        self.cmass1 = CMASS1(self._h5n, self)
        self.cmass2 = CMASS2(self._h5n, self)
        self.cmass3 = CMASS3(self._h5n, self)
        self.cmass4 = CMASS4(self._h5n, self)
        self.conm2 = CONM2(self._h5n, self)
        self.conrod = CONROD(self._h5n, self)
        self.contrlt = CONTRLT(self._h5n, self)
        self.cpenta = CPENTA(self._h5n, self)
        self.cqdx4fd = CQDX4FD(self._h5n, self)
        self.cqdx9fd = CQDX9FD(self._h5n, self)
        self.cquad = CQUAD(self._h5n, self)
        self.cquad4 = CQUAD4(self._h5n, self)
        self.cquad4fd = CQUAD4FD(self._h5n, self)
        self.cquad8 = CQUAD8(self._h5n, self)
        self.cquad9fd = CQUAD9FD(self._h5n, self)
        self.cquadr = CQUADR(self._h5n, self)
        self.cquadx = CQUADX(self._h5n, self)
        self.crbe1 = CRBE1(self._h5n, self)
        self.crod = CROD(self._h5n, self)
        self.cseam = CSEAM(self._h5n, self)
        self.cseamp = CSEAMP(self._h5n, self)
        self.cshear = CSHEAR(self._h5n, self)
        self.cslot3 = CSLOT3(self._h5n, self)
        self.cslot4 = CSLOT4(self._h5n, self)
        self.ctetra = CTETRA(self._h5n, self)
        self.ctria3 = CTRIA3(self._h5n, self)
        self.ctria3fd = CTRIA3FD(self._h5n, self)
        self.ctria6 = CTRIA6(self._h5n, self)
        self.ctria6fd = CTRIA6FD(self._h5n, self)
        self.ctriah = CTRIAH(self._h5n, self)
        self.ctriar = CTRIAR(self._h5n, self)
        self.ctriax = CTRIAX(self._h5n, self)
        self.ctriax6 = CTRIAX6(self._h5n, self)
        self.ctrix3fd = CTRIX3FD(self._h5n, self)
        self.ctrix6fd = CTRIX6FD(self._h5n, self)
        self.ctube = CTUBE(self._h5n, self)
        self.cvisc = CVISC(self._h5n, self)
        self.cweld = CWELD(self._h5n, self)
        self.cweldc = CWELDC(self._h5n, self)
        self.cweldp = CWELDP(self._h5n, self)
        # self.genel = GENEL(self._h5n, self)
        self.plotel = PLOTEL(self._h5n, self)
        self.prim1 = PRIM1(self._h5n, self)
        self.prim2 = PRIM2(self._h5n, self)
        self.prim3 = PRIM3(self._h5n, self)
        self.prim4 = PRIM4(self._h5n, self)
        self.prim5 = PRIM5(self._h5n, self)
        self.prim6 = PRIM6(self._h5n, self)
        self.prim7 = PRIM7(self._h5n, self)
        self.prim8 = PRIM8(self._h5n, self)
        self.radcol = RADCOL(self._h5n, self)
        self.rbar = RBAR(self._h5n, self)
        self.rbar1 = RBAR1(self._h5n, self)
        self.rbe1 = RBE1(self._h5n, self)
        self.rbe2 = RBE2(self._h5n, self)
        self.rbe2gs = RBE2GS(self._h5n, self)
        self.rbe3 = RBE3(self._h5n, self)
        self.ringax = RINGAX(self._h5n, self)
        self.rjoint = RJOINT(self._h5n, self)
        self.rrod = RROD(self._h5n, self)
        self.rsscon = RSSCON(self._h5n, self)
        self.rtrplt = RTRPLT(self._h5n, self)
        self.rtrplt1 = RTRPLT1(self._h5n, self)
        self.sectax = SECTAX(self._h5n, self)
        self.wetelme = WETELME(self._h5n, self)
        self.wetelmg = WETELMG(self._h5n, self)

    def path(self):
        return self._input.path() + ['ElEMENT']

    def read(self):
        for key, item in iteritems(self.__dict__):
            if key.startswith('_'):
                continue
            try:
                item.read()
            except AttributeError:
                pass

########################################################################################################################


class AEQUAD4(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/AEQUAD4')

########################################################################################################################


class AEROD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/AEROD')

########################################################################################################################


class AEROQ4(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/AEROQ4')

########################################################################################################################


class AEROT3(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/AEROT3')

########################################################################################################################


class AETRIA3(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/AETRIA3')

########################################################################################################################


class BEAMAERO(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/BEAMAERO')

########################################################################################################################


class BOLT(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/BOLT/IDENTITY')

########################################################################################################################


class CAABSF(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CAABSF')

########################################################################################################################


class CACINF3(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CACINF3')

########################################################################################################################


class CACINF4(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CACINF4')

########################################################################################################################


class CAERO1(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CAERO1')

########################################################################################################################


class CAERO2(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CAERO2')

########################################################################################################################


class CAERO3(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CAERO3')

########################################################################################################################


class CAERO4(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CAERO4')

########################################################################################################################


class CAERO5(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CAERO5')

########################################################################################################################


class CAXIF2(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CAXIF2')

########################################################################################################################


class CAXIF3(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CAXIF3')

########################################################################################################################


class CAXIF4(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CAXIF4')

########################################################################################################################


class CAXISYM(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CAXISYM')

########################################################################################################################


class CBAR(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CBAR')

    @staticmethod
    def from_bdf(card):
        x = card.x
        if x is None:
            x = [None, None, None]
            flag = 2
        else:
            x = x.tolist()
            offt = card.offt
            if offt == '':
                offt = 'GGG'
            if offt.startswith('G'):
                flag = 1
            else:
                flag = 0
        # TODO: check that flag is correct
        data = [card.eid, card.pid]
        data += card.node_ids
        data.append(flag)
        data += x
        data += [card.g0, card.pa, card.pb]  # card.offt not used...
        data += card.wa.tolist()
        data += card.wb.tolist()
        return TableData([data])

########################################################################################################################

# CBARAO not in spec

class CBARAO_SPEC(object):
    name = 'CBARAO'
    path = '/NASTRAN/INPUT/ELEMENT/CBARAO'
    dtype = [('EID', '<i8', ()), ('SCALE', 'S2', ()), ('X', '<f8', (6,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


class CBARAO(CardTable):
    table_def = TableDef.create(CBARAO_SPEC)

    @staticmethod
    def from_bdf(card):
        return TableData([[card.eid, card.scale, card.x]])


########################################################################################################################


class CBEAM(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CBEAM')

    @staticmethod
    def from_bdf(card):
        x = card.x
        if x is None:
            x = [None, None, None]
            flag = 2
        else:
            offt = card.offt
            if offt == '':
                offt = 'GGG'
            if offt.startswith('G'):
                flag = 1
            else:
                flag = 0
        # TODO: CBEAM flag check
        node_ids = card.node_ids
        data = [card.eid, card.pid, node_ids[0], node_ids[1], card.sa, card.sb, x, card.g0, flag, card.pa, card.pb,
                card.wa, card.wb]
        return TableData([data])

########################################################################################################################


class CBEAM3(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CBEAM3')

########################################################################################################################


class CBEND(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CBEND')

    @staticmethod
    def from_bdf(card):
        node_ids = card.node_ids
        # TODO: CBEND flag check
        x = card.x
        g0 = card.g0
        flag = 0
        return TableData([[card.eid, card.pid, node_ids[0], node_ids[1], flag, x[0], x[1], x[2], g0, card.geom]])

########################################################################################################################


class CBUSH(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CBUSH')

    @staticmethod
    def from_bdf(card):
        x = card.x
        if x is None:
            x = [None, None, None]
        cid = card.cid
        if cid not in ('', None):
            if cid == 0:
                if None in x:
                    flag = 0
                else:
                    flag = 1
            else:
                flag = -1
        else:
            if None in x:
                flag = 2
            else:
                flag = 1

        # TODO: CBUSH flag check

        node_ids = card.node_ids
        si = card.si
        data = [card.eid, card.pid, node_ids[0], node_ids[1], flag, x[0], x[1], x[2], card.g0, cid, card.s,
                card.ocid, si[0], si[1], si[2]]
        return TableData([data])

########################################################################################################################


class CBUSH1D(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CBUSH1D')

    @staticmethod
    def from_bdf(card):
        return TableData([[card.eid, card.pid, card.node_ids, card.cid]])

########################################################################################################################


class CBUSH2D(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CBUSH2D')

    @staticmethod
    def from_bdf(card):
        return TableData([[card.eid, card.pid, card.node_ids, card.cid, card.plane, card.sptid]])

########################################################################################################################


class CCONEAX(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CCONEAX')

    @staticmethod
    def from_bdf(card):
        r = card.ring_ids
        return TableData([[card.eid, card.pid, r[0], r[1]]])

########################################################################################################################


class CDAMP1(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CDAMP1')

########################################################################################################################


class CDAMP2(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CDAMP2')

########################################################################################################################


class CDAMP3(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CDAMP3')

########################################################################################################################


class CDAMP4(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CDAMP4')

########################################################################################################################


class CDAMP5(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CDAMP5')

########################################################################################################################


class CELAS1(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CELAS1')

########################################################################################################################


class CELAS2(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CELAS2')

########################################################################################################################


class CELAS3(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CELAS3')

########################################################################################################################


class CELAS4(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CELAS4')

########################################################################################################################


class CFAST(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CFAST')

########################################################################################################################


class CFASTP(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CFASTP')

########################################################################################################################


class CFLUID2(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CFLUID2')

########################################################################################################################


class CFLUID3(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CFLUID3')

########################################################################################################################


class CFLUID4(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CFLUID4')

########################################################################################################################


class CGAP(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CGAP')

    @staticmethod
    def from_bdf(card):
        nids = card.node_ids
        x = card.x
        flag = 0
        # TODO: CGAP flag check
        return [card.eid, card.pid, nids[0], nids[1], flag, x[0], x[1], x[2], card.g0, card.cid]

########################################################################################################################


class CHACAB(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CHACAB')

########################################################################################################################


class CHACBR(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CHACBR')

########################################################################################################################


class CHBDYE(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CHBDYE')

########################################################################################################################


class CHBDYG(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CHBDYG')

########################################################################################################################


class CHBDYP(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CHBDYP')

########################################################################################################################


class CHEXA(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CHEXA')

    @staticmethod
    def from_bdf(card):
        return TableData([[card.eid, card.pid, card.node_ids]])

########################################################################################################################


class CHEXAL(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CHEXAL')

########################################################################################################################


class CHEXP(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CHEXP')

########################################################################################################################


class CIFHEX(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CIFHEX')

########################################################################################################################


class CIFPENT(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CIFPENT')

########################################################################################################################


class CIFQDX(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CIFQDX')

########################################################################################################################


class CIFQUAD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CIFQUAD')

########################################################################################################################


class CMASS1(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CMASS1')

########################################################################################################################


class CMASS2(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CMASS2')

########################################################################################################################


class CMASS3(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CMASS3')

########################################################################################################################


class CMASS4(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CMASS4')

########################################################################################################################


class CONM2(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CONM2')

########################################################################################################################


class CONROD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CONROD')

########################################################################################################################


class CONTRLT(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CONTRLT')

########################################################################################################################


class CPENTA(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CPENTA')

########################################################################################################################


class CQDX4FD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CQDX4FD')

########################################################################################################################


class CQDX9FD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CQDX9FD')

########################################################################################################################


class CQUAD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CQUAD')

########################################################################################################################


class CQUAD4(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CQUAD4')

    @staticmethod
    def from_bdf(card):
        theta = card.theta_mcid
        if not isinstance(theta, float):
            mcid = theta
            theta = None
        else:
            mcid = None
        data = [card.eid, card.pid, card.node_ids, theta, card.zoffset, card.tflag,
                [card.T1, card.T2, card.T3, card.T4], mcid]
        return TableData([data])

########################################################################################################################


class CQUAD4FD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CQUAD4FD')

########################################################################################################################


class CQUAD8(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CQUAD8')

########################################################################################################################


class CQUAD9FD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CQUAD9FD')

########################################################################################################################


class CQUADR(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CQUADR')

########################################################################################################################


class CQUADX(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CQUADX')

########################################################################################################################


class CRBE1(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CRBE1/IDENTITY',
                                rename={'LAGMULTIPL_POS': 'LAGM_POS', 'LAGMULTIPL_LEN': 'LAGM_LEN'})

########################################################################################################################


class CROD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CROD')

    @staticmethod
    def from_bdf(card):
        return TableData([[card.eid, card.pid, card.node_ids]])

########################################################################################################################


class CSEAM(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CSEAM')

########################################################################################################################


class CSEAMP(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CSEAMP')

########################################################################################################################


class CSHEAR(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CSHEAR')

    @staticmethod
    def from_bdf(card):
        return TableData([[card.eid, card.pid, card.node_ids]])

########################################################################################################################


class CSLOT3(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CSLOT3')

########################################################################################################################


class CSLOT4(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CSLOT4')

########################################################################################################################


class CTETRA(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTETRA')

########################################################################################################################


class CTRIA3(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTRIA3')

    @staticmethod
    def from_bdf(card):
        theta = card.theta_mcid
        if not isinstance(theta, float):
            mcid = theta
            theta = None
        else:
            mcid = None
        data = [card.eid, card.pid, card.node_ids, theta, card.zoffset, card.tflag,
                [card.T1, card.T2, card.T3], mcid]
        return TableData([data])

########################################################################################################################


class CTRIA3FD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTRIA3FD')

########################################################################################################################


class CTRIA6(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTRIA6')

########################################################################################################################


class CTRIA6FD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTRIA6FD')

########################################################################################################################


class CTRIAH(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTRIAH')

########################################################################################################################


class CTRIAR(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTRIAR')

########################################################################################################################


class CTRIAX(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTRIAX')

########################################################################################################################


class CTRIAX6(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTRIAX6')

########################################################################################################################


class CTRIX3FD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTRIX3FD')

########################################################################################################################


class CTRIX6FD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTRIX6FD')

########################################################################################################################


class CTUBE(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTUBE')

########################################################################################################################


class CVISC(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CVISC')

########################################################################################################################


class CWELD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CWELD')

########################################################################################################################


class CWELDC(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CWELDC')

########################################################################################################################


class CWELDP(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CWELDP')

########################################################################################################################


# TODO: GENEL doesn't conform
# class GENEL(CardTable):
#     table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/GENEL/IDENTITY',
#                                 rename={'UDLIST_POS': 'UD_POS', 'UDLIST_LEN': 'UD_LEN', 'UILIST_POS': 'UI_POS',
#                                         'UILIST_LEN': 'UI_LEN'})

########################################################################################################################


class PLOTEL(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/PLOTEL')

########################################################################################################################


class PRIM1(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/PRIM1')

########################################################################################################################


class PRIM2(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/PRIM2')

########################################################################################################################


class PRIM3(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/PRIM3')

########################################################################################################################


class PRIM4(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/PRIM4')

########################################################################################################################


class PRIM5(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/PRIM5')

########################################################################################################################


class PRIM6(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/PRIM6')

########################################################################################################################


class PRIM7(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/PRIM7')

########################################################################################################################


class PRIM8(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/PRIM8')

########################################################################################################################


class RADCOL(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/RADCOL')

########################################################################################################################


class RBAR(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/RBAR')

########################################################################################################################


class RBAR1(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/RBAR1')

########################################################################################################################


class RBE1(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/RBE1/IDENTITY')

########################################################################################################################


class RBE2(CardTable):
    table_def = TableDef.create(
        '/NASTRAN/INPUT/ELEMENT/RBE2/RB',
        subtables=[TableDef.create('/NASTRAN/INPUT/ELEMENT/RBE2/GM')],
    )

    @staticmethod
    def from_bdf(card):
        data = TableData()
        data.data = [[card.eid, card.gn, card.cm, card.alpha]]
        data.subdata_len = np.array([[len(card.Gmi)]])

        gmi = np.array(card.Gmi)
        gmi = gmi.reshape((gmi.shape[0], 1))

        subdata = TableData()
        subdata.data = gmi

        data.subdata.append(subdata)

        try:
            data.validate()
        except AssertionError:
            raise AssertionError('RBE2')

        return data

########################################################################################################################


class RBE2GS(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/RBE2GS/IDENTITY')

########################################################################################################################


class RBE3(CardTable):
    table_def = TableDef.create(
        '/NASTRAN/INPUT/ELEMENT/RBE3/IDENTITY',
    subtables=[
        TableDef.create('/NASTRAN/INPUT/ELEMENT/RBE3/GM'),
        TableDef.create(
            '/NASTRAN/INPUT/ELEMENT/RBE3/WTCG',
            subtables=[TableDef.create('/NASTRAN/INPUT/ELEMENT/RBE3/G')]
        ),
    ]
    )

    g_len = 0

    @staticmethod
    def from_bdf(card):
        data = TableData()
        data.data = [[card.eid, card.refgrid, card.refc, card.alpha]]

        #################################

        Gmi = card.Gmi
        Cmi = card.Cmi

        gm = [(Gmi[i], Cmi[i]) for i in range(len(Gmi))]

        subdata = TableData(gm)

        data.subdata.append(subdata)
        gm_len = len(subdata.data)

        ################################

        weights = card.weights
        comps = card.comps

        wtcg = [(weights[i], comps[i]) for i in range(len(weights))]

        subdata = TableData(wtcg)

        Gijs = card.Gijs
        subdata.subdata_len = np.array([[len(Gijs[i])] for i in range(len(Gijs))])

        gijs = [Gijs[i] for i in range(len(Gijs))]

        all_gijs = []

        for gij in gijs:
            all_gijs.extend(gij)

        gijs = TableData(np.array(all_gijs).reshape((len(all_gijs), 1)))

        subdata.subdata = [gijs]

        data.subdata.append(subdata)
        data.subdata_len = np.array([[gm_len, len(subdata.data)]])

        try:
            data.validate()
        except AssertionError:
            raise AssertionError('RBE3')

        return data

########################################################################################################################


class RINGAX(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/RINGAX')

########################################################################################################################


class RJOINT(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/RJOINT')

########################################################################################################################


class RROD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/RROD')

########################################################################################################################


class RSSCON(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/RSSCON')

########################################################################################################################


class RTRPLT(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/RTRPLT')

########################################################################################################################


class RTRPLT1(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/RTRPLT1')

########################################################################################################################


class SECTAX(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/SECTAX')

########################################################################################################################


class WETELME(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/WETELME')

########################################################################################################################


class WETELMG(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/WETELMG/IDENTITY')

########################################################################################################################
