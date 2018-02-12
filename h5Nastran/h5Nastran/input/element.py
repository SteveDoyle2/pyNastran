from __future__ import print_function, absolute_import
from six import iteritems, itervalues, iterkeys
from six.moves import range

import numpy as np

from .card_table import CardTable, TableDef
from ..data_helper import DataHelper
from ._shell_element_info import QuadShell, TriaShell, shell_element_info_format


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

        self._shell_element_info = None
        self._shell_element_info_dict = None

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

    def write_shell_element_info(self, bdf, cards):
        if self._shell_element_info is not None:
            return self._shell_element_info

        tables = self.__dict__
        table_ids = self.__dict__.keys()

        element_info = None

        for table_id in table_ids:
            if table_id.startswith('_'):
                continue

            try:
                _element_info = tables[table_id].get_shell_info(bdf, cards)
            except AttributeError:
                continue

            if element_info is None:
                element_info = _element_info
            else:
                element_info = np.append(element_info, _element_info)

        element_info = element_info[element_info['EID'].argsort()]

        h5f = self._h5n.h5f
        table = h5f.create_table('/PRIVATE/INPUT', 'SHELL_ELEMENT_INFO', shell_element_info_format,
                                 'SHELL ELEMENT INFO', expectedrows=len(element_info), createparents=True)
        table.append(element_info)
        h5f.flush()

        self._shell_element_info = element_info

        return element_info

    def get_shell_element_info(self):
        if self._shell_element_info is None:
            self._shell_element_info = self._h5n.h5f.get_node('/PRIVATE/INPUT/SHELL_ELEMENT_INFO').read()
        return self._shell_element_info

    def get_shell_element_info_dict(self):
        if self._shell_element_info_dict is not None:
            return self._shell_element_info_dict

        element_shell_info = self.get_shell_element_info()

        result = {'CENTER': {}, 'THETA_RAD': {}}

        eid = element_shell_info['EID']
        center = element_shell_info['CENTER']
        theta = element_shell_info['THETA_RAD']

        _center = result['CENTER']
        _theta = result['THETA_RAD']

        for i in range(len(eid)):
            _center[eid[i]] = center[i]
            _theta[eid[i]] = theta[i]

        self._shell_element_info_dict = result

        return self._shell_element_info_dict


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

    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=cls.table_def.dtype)

        eid = data['EID']
        pid = data['PID']
        ga = data['GA']
        gb = data['GB']
        flag = data['FLAG']
        x1 = data['X1']
        x2 = data['X2']
        x3 = data['X3']
        go = data['GO']
        pa = data['PA']
        pb = data['PB']
        w1a = data['W1A']
        w2a = data['W2A']
        w3a = data['W3A']
        w1b = data['W1B']
        w2b = data['W2B']
        w3b = data['W3B']
        # domain_id = data['DOMAIN_ID']

        i = -1

        for card_id in card_ids:
            i += 1
            card = cards[card_id]
            x = card.x
            if x is None:
                x = [None, None, None]
            else:
                x = x.tolist()

            # MSC 2018 DMAP p. 209
            # flag = 0, basic coordinate system
            # flag = 1, global coordinate system
            # flag = 2, grid option
            # TODO: what does *F = FE bit-wise AND with 3 mean in DMAP guide?

            g0 = card.g0

            if x[0] is None:
                _flag = 2
            elif g0 in ('', None):
                g0 = DataHelper.default_int
                _flag = 0
            else:
                _flag = 1

            eid[i] = card.eid
            pid[i] = card.pid
            ga[i], gb[i] = card.node_ids
            flag[i] = _flag

            if x[0] is None:
                x1[i], x2[i], x3[i] = [DataHelper.default_double] * 3
            else:
                x1[i], x2[i], x3[i] = x

            go[i] = g0
            pa[i] = card.pa
            pb[i] = card.pb
            w1a[i], w2a[i], w3a[i] = card.wa.tolist()
            w1b[i], w2b[i], w3b[i] = card.wb.tolist()

        result = {
            'IDENTITY': data
        }

        return result


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

    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=cls.table_def.dtype)

        eid = data['EID']
        scale = data['SCALE']
        x = data['X']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]
            eid[i] = card.eid
            scale[i] = card.scale
            x[i] = card.x

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class CBEAM(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CBEAM')

    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=cls.table_def.dtype)

        eid = data['EID']
        pid = data['PID']
        ga = data['GA']
        gb = data['GB']
        sa = data['SA']
        sb = data['SB']
        x = data['X']
        g0 = data['G0']
        f = data['F']
        pa = data['PA']
        pb = data['PB']
        wa = data['WA']
        wb = data['WB']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            _x = card.x
            if _x is None:
                _x = [DataHelper.default_double, DataHelper.default_double, DataHelper.default_double]

            eid[i] = card.eid
            pid[i] = card.pid
            ga[i], gb[i] = card.node_ids
            sa[i] = card.sa
            sb[i] = card.sb
            x[i] = _x
            g0[i] = card.g0
            f[i] = DataHelper.unknown_int  # TODO: CBEAM flag
            pa[i] = card.pa
            pb[i] = card.pb
            wa[i] = card.wa
            wb[i] = card.wb

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class CBEAM3(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CBEAM3')


########################################################################################################################


class CBEND(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CBEND')

    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=cls.table_def.dtype)

        eid = data['EID']
        pid = data['PID']
        ga = data['GA']
        gb = data['GB']
        flag = data['FLAG']
        x1 = data['X1']
        x2 = data['X2']
        x3 = data['X3']
        go = data['GO']
        geom = data['GEOM']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            eid[i] = card.eid
            pid[i] = card.pid
            ga[i], gb[i] = card.node_ids
            flag[i] = DataHelper.unknown_int  # TODO: CBEND flag
            x1[i], x2[i], x3[i] = card.x
            go[i] = card.g0
            geom[i] = card.geom

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class CBUSH(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CBUSH')

    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=cls.table_def.dtype)

        eid = data['EID']
        pid = data['PID']
        ga = data['GA']
        gb = data['GB']
        flag = data['FLAG']
        x1 = data['X1']
        x2 = data['X2']
        x3 = data['X3']
        go = data['GO']
        cid = data['CID']
        s = data['S']
        ocid = data['OCID']
        s1 = data['S1']
        s2 = data['S2']
        s3 = data['S3']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            _x = card.x
            if _x is None:
                _x = [DataHelper.default_double, DataHelper.default_double, DataHelper.default_double]

            _g0 = card.g0
            if _g0 is None:
                _g0 = DataHelper.default_int

            eid[i] = card.eid
            pid[i] = card.pid
            ga[i], gb[i] = card.node_ids
            flag[i] = DataHelper.unknown_int  # TODO: CBUSH flag
            x1[i], x2[i], x3[i] = _x
            go[i] = _g0
            cid[i] = card.cid
            s[i] = card.s
            ocid[i] = card.ocid
            s1[i], s2[i], s3[i] = card.si

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class CBUSH1D(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CBUSH1D')

    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=cls.table_def.dtype)

        eid = data['EID']
        pid = data['PID']
        g = data['G']
        cid = data['CID']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            eid[i] = card.eid
            pid[i] = card.pid
            g[i] = card.node_ids
            cid[i] = card.cid

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class CBUSH2D(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CBUSH2D')

    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=cls.table_def.dtype)

        eid = data['EID']
        pid = data['PID']
        g = data['G']
        cid = data['CID']
        plane = data['PLANE']
        sptid = data['SPTID']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            eid[i] = card.eid
            pid[i] = card.pid
            g[i] = card.node_ids
            cid[i] = card.cid
            plane[i] = card.plane
            sptid[i] = card.sptid

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class CCONEAX(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CCONEAX')

    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=cls.table_def.dtype)

        eid = data['EID']
        pid = data['PID']
        ra = data['RA']
        rb = data['RB']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            eid[i] = card.eid
            pid[i] = card.pid
            ra[i], rb[i] = card.ring_ids

        result = {
            'IDENTITY': data
        }

        return result


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

    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=cls.table_def.dtype)

        eid = data['EID']
        pid = data['PID']
        ga = data['GA']
        gb = data['GB']
        f = data['F']
        x1 = data['X1']
        x2 = data['X2']
        x3 = data['X3']
        go = data['GO']
        cid = data['CID']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            _x = card.x
            if _x[0] is None:
                _x = [DataHelper.default_double] * 3

            _g0 = card.g0
            if _g0 is None:
                _g0 = DataHelper.default_int

            eid[i] = card.eid
            pid[i] = card.pid
            ga[i], gb[i] = card.node_ids
            f[i] = DataHelper.unknown_int  # TODO: CGAP flag
            x1[i], x2[i], x3[i] = _x
            go[i] = _g0
            cid[i] = card.cid

        result = {
            'IDENTITY': data
        }

        return result


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

    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=cls.table_def.dtype)

        eid = data['EID']
        pid = data['PID']
        g = data['G']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            node_ids = list(card.node_ids)
            diff = 20 - len(node_ids)
            if diff > 0:
                node_ids += [0] * diff

            eid[i] = card.eid
            pid[i] = card.pid
            g[i] = node_ids

        result = {
            'IDENTITY': data
        }

        return result


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


class CQUAD4(CardTable, QuadShell):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CQUAD4')

    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=cls.table_def.dtype)

        eid = data['EID']
        pid = data['PID']
        g = data['G']
        theta = data['THETA']
        zoffs = data['ZOFFS']
        tflag = data['TFLAG']
        t = data['T']
        mcid = data['MCID']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            theta_mcid = card.theta_mcid
            if not isinstance(theta_mcid, float):
                _mcid = theta_mcid
                _theta = DataHelper.default_double
            else:
                _mcid = DataHelper.default_int
                _theta = theta_mcid

            eid[i] = card.eid
            pid[i] = card.pid
            g[i] = card.node_ids
            theta[i] = _theta
            zoffs[i] = card.zoffset
            tflag[i] = card.tflag
            t[i] = [card.T1, card.T2, card.T3, card.T4]
            mcid[i] = _mcid

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class CQUAD4FD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CQUAD4FD')


########################################################################################################################


class CQUAD8(CardTable, QuadShell):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CQUAD8')


########################################################################################################################


class CQUAD9FD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CQUAD9FD')


########################################################################################################################


class CQUADR(CardTable, QuadShell):
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

    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=cls.table_def.dtype)

        eid = data['EID']
        pid = data['PID']
        g = data['G']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            eid[i] = card.eid
            pid[i] = card.pid
            g[i] = card.node_ids

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class CSEAM(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CSEAM')


########################################################################################################################


class CSEAMP(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CSEAMP')


########################################################################################################################


class CSHEAR(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CSHEAR')

    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=cls.table_def.dtype)

        eid = data['EID']
        pid = data['PID']
        g = data['G']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            eid[i] = card.eid
            pid[i] = card.pid
            g[i] = card.node_ids

        result = {
            'IDENTITY': data
        }

        return result


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


class CTRIA3(CardTable, TriaShell):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTRIA3')

    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=cls.table_def.dtype)

        eid = data['EID']
        pid = data['PID']
        g = data['G']
        theta = data['THETA']
        zoffs = data['ZOFFS']
        tflag = data['TFLAG']
        t = data['T']
        mcid = data['MCID']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            theta_mcid = card.theta_mcid
            if not isinstance(theta_mcid, float):
                _mcid = theta_mcid
                _theta = DataHelper.default_double
            else:
                _mcid = DataHelper.default_int
                _theta = theta_mcid

            eid[i] = card.eid
            pid[i] = card.pid
            g[i] = card.node_ids
            theta[i] = _theta
            zoffs[i] = card.zoffset
            tflag[i] = card.tflag
            t[i] = [card.T1, card.T2, card.T3]
            mcid[i] = _mcid

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class CTRIA3FD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTRIA3FD')


########################################################################################################################


class CTRIA6(CardTable, TriaShell):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTRIA6')


########################################################################################################################


class CTRIA6FD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTRIA6FD')


########################################################################################################################


class CTRIAH(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTRIAH')


########################################################################################################################


class CTRIAR(CardTable, TriaShell):
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

    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())

        gm = {
            'IDENTITY': {'ID': []}
        }

        result = {
            'IDENTITY': {'EID': [], 'GN': [], 'CM': [], 'GM_POS': [], 'GM_LEN': [], 'ALPHA': [],
                         'DOMAIN_ID': []},
            'GM': gm,
            '_subtables': ['GM']
        }

        identity = result['IDENTITY']

        eid = identity['EID']
        gn = identity['GN']
        cm = identity['CM']
        gm_pos = identity['GM_POS']
        gm_len = identity['GM_LEN']
        alpha = identity['ALPHA']
        id = gm['IDENTITY']['ID']

        _pos = 0
        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            eid.append(card.eid)
            gn.append(card.gn)
            cm.append(card.cm)
            alpha.append(card.alpha)
            gm_pos.append(_pos)
            _gm_len = len(card.Gmi)
            gm_len.append(_gm_len)
            _pos += _gm_len
            id += list(card.Gmi)

        return result


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

    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())

        g = {'IDENTITY': {'ID': []}}
        wtcg = {
            'IDENTITY': {'WT1': [], 'C': [], 'G_POS': [], 'G_LEN': []},
            'G': g,
            '_subtables': ['G']
        }
        gm = {'IDENTITY': {'GM': [], 'CM': []}}

        result = {'IDENTITY': {'EID': [], 'REFG': [], 'REFC': [], 'WTCG_POS': [],
                               'WTCG_LEN': [], 'GM_POS': [], 'GM_LEN': [], 'ALPHA': [],
                               'DOMAIN_ID': []},
                  'WTCG': wtcg, 'GM': gm,
                  '_subtables': ['GM', 'WTCG']
                  }

        id = g['IDENTITY']['ID']
        _gm = gm['IDENTITY']
        gm = _gm['GM']
        cm = _gm['CM']
        _wtcg = wtcg['IDENTITY']
        wt1 = _wtcg['WT1']
        c = _wtcg['C']
        g_pos = _wtcg['G_POS']
        g_len = _wtcg['G_LEN']
        identity = result['IDENTITY']
        eid = identity['EID']
        refg = identity['REFG']
        refc = identity['REFC']
        wtcg_pos = identity['WTCG_POS']
        wtcg_len = identity['WTCG_LEN']
        gm_pos = identity['GM_POS']
        gm_len = identity['GM_LEN']
        alpha = identity['ALPHA']

        _gm_pos = 0
        _wtcg_pos = 0
        _g_pos = 0

        for card_id in card_ids:
            card = cards[card_id]

            eid.append(card.eid)
            refg.append(card.refgrid)
            refc.append(card.refc)
            alpha.append(card.alpha)

            wtcg_pos.append(_wtcg_pos)
            _wtcg_len = len(card.weights)
            wtcg_len.append(_wtcg_len)
            _wtcg_pos += _wtcg_len

            gm_pos.append(_gm_pos)
            _gm_len = len(card.Gmi)
            gm_len.append(_gm_len)
            _gm_pos += _gm_len

            gm += list(card.Gmi)
            cm += list(card.Cmi)

            wt1.extend(list(card.weights))
            c.extend(list(card.comps))

            gijs = []
            for gij in card.Gijs:
                gijs += gij
                g_pos.append(_g_pos)
                _g_len = len(gij)
                g_len.append(_g_len)
                _g_pos += _g_len

            id += gijs

        return result


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
