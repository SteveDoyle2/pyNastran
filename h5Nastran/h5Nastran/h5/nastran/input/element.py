from __future__ import print_function, absolute_import

import numpy as np
from six import iteritems
from six.moves import range

from h5Nastran.defaults import Defaults
from h5Nastran.h5nastrannode import H5NastranNode
from ._shell_element_info import QuadShell, TriaShell, shell_element_info_format
from .input_table import InputTable, TableDef


class Element(H5NastranNode):
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
        table_paths = self._h5n.table_paths
        table = h5f.create_table(table_paths.shell_element_info_path, table_paths.shell_element_info_table,
                                 shell_element_info_format, 'SHELL ELEMENT INFO',
                                 expectedrows=len(element_info), createparents=True)
        table.append(element_info)
        h5f.flush()

        self._shell_element_info = element_info

        return element_info

    def get_shell_element_info(self):
        if self._shell_element_info is None:
            table_paths = self._h5n.table_paths
            self._shell_element_info = self._h5n.h5f.get_node(table_paths.shell_element_info).read()
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

    def to_bdf(self, bdf):
        for key, item in iteritems(self.__dict__):
            if key.startswith('_'):
                continue
            try:
                item.to_bdf(bdf)
            except NotImplementedError:
                pass


########################################################################################################################


class AEQUAD4(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/AEQUAD4')


########################################################################################################################


class AEROD(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/AEROD')


########################################################################################################################


class AEROQ4(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/AEROQ4')


########################################################################################################################


class AEROT3(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/AEROT3')


########################################################################################################################


class AETRIA3(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/AETRIA3')


########################################################################################################################


class BEAMAERO(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/BEAMAERO')


########################################################################################################################


class BOLT(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/BOLT/IDENTITY')


########################################################################################################################


class CAABSF(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CAABSF')


########################################################################################################################


class CACINF3(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CACINF3')


########################################################################################################################


class CACINF4(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CACINF4')


########################################################################################################################


class CAERO1(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CAERO1')
    
    def to_bdf(self, bdf):
        add_card = bdf.add_caero1
        
        data = self.identity
        
        eid = data['EID']
        pid = data['PID']
        cp = data['CP']
        nspan = data['NSPAN']
        nchord = data['NCHORD']
        lspan = data['LSPAN']
        lchord = data['LCHORD']
        igid = data['IGID']
        x1 = data['X1']
        y1 = data['Y1']
        z1 = data['Z1']
        x12 = data['X12']
        x4 = data['X4']
        y4 = data['Y4']
        z4 = data['Z4']
        x43 = data['X43']
        
        for i in range(eid.shape[0]):
            add_card(eid[i], pid[i], igid[i], np.array([x1[i], y1[i], z1[i]]), x12[i],
                     np.array([x4[i], y4[i], z4[i]]), x43[i], cp[i], nspan[i], lspan[i], nchord[i], lchord[i])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        eid = data['EID']
        pid = data['PID']
        cp = data['CP']
        nspan = data['NSPAN']
        nchord = data['NCHORD']
        lspan = data['LSPAN']
        lchord = data['LCHORD']
        igid = data['IGID']
        x1 = data['X1']
        y1 = data['Y1']
        z1 = data['Z1']
        x12 = data['X12']
        x4 = data['X4']
        y4 = data['Y4']
        z4 = data['Z4']
        x43 = data['X43']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            eid[i] = card.eid
            pid[i] = card.pid
            cp[i] = card.cp
            nspan[i] = card.nspan
            nchord[i] = card.nchord
            lspan[i] = card.lspan
            lchord[i] = card.lchord
            igid[i] = card.igid
            x1[i], y1[i], z1[i] = card.p1
            x12[i] = card.x12
            x4[i], y4[i], z4[i] = card.p4
            x43[i] = card.x43

        return {'IDENTITY': data}


########################################################################################################################


class CAERO2(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CAERO2')


########################################################################################################################


class CAERO3(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CAERO3')


########################################################################################################################


class CAERO4(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CAERO4')


########################################################################################################################


class CAERO5(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CAERO5')


########################################################################################################################


class CAXIF2(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CAXIF2')


########################################################################################################################


class CAXIF3(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CAXIF3')


########################################################################################################################


class CAXIF4(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CAXIF4')


########################################################################################################################


class CAXISYM(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CAXISYM')


########################################################################################################################


class CBAR(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CBAR')

    def to_bdf(self, bdf):
        add_card = bdf.add_cbar

        data = self.identity

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

        defaults = self._h5n.defaults
        to_int = defaults.to_value_int

        for i in range(eid.shape[0]):
            _flag = flag[i]
            if _flag == 0:
                offt = 'BGG'
            else:
                offt = 'GGG'
            add_card(eid[i], to_int(pid[i]), [ga[i], gb[i]], [x1[i], x2[i], x3[i]], to_int(go[i]), offt=offt,
                     pa=to_int(pa[i]), pb=to_int(pb[i]),
                     wa=[w1a[i], w2a[i], w3a[i]], wb=[w1b[i], w2b[i], w3b[i]])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

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
            offt = card.offt
            first_offt = offt[0]

            if x[0] is None:
                _flag = 2
            elif first_offt == 'B':
                g0 = Defaults.default_int
                _flag = 0
            elif first_offt == 'G':
                g0 = Defaults.default_int
                _flag = 1

            # TODO: offt is not stored, need to convert wa and wb to basic coordinate system

            eid[i] = card.eid
            pid[i] = card.pid
            ga[i], gb[i] = card.node_ids
            flag[i] = _flag

            if x[0] is None:
                x1[i], x2[i], x3[i] = [Defaults.default_double] * 3
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


class CBARAO(InputTable):
    table_def = TableDef.create(CBARAO_SPEC)

    def to_bdf(self, bdf):
        add_card = bdf.add_cbarao

        data = self.identity
        eid = data['EID']
        scale = data['SCALE']
        x = data['X']

        for i in range(eid.shape[0]):
            add_card(eid[i], scale[i], x[i])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

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


class CBEAM(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CBEAM')

    def to_bdf(self, bdf):
        add_card = bdf.add_cbeam

        data = self.identity

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

        # TODO: what about p-elements for bit?

        for i in range(eid.size):
            _flag = f[i]
            if _flag == 0:
                offt = 'BGG'
            else:
                offt = 'GGG'

            _x = x[i]
            if _x[0] == Defaults.default_double:
                _x = None

            _g0 = g0[i]

            if _g0 == Defaults.default_int:
                _g0 = None

            add_card(eid[i], pid[i], [ga[i], gb[i]], _x, _g0, offt, bit=None,
                     pa=pa[i], pb=pb[i], wa=wa[i], wb=wb[i], sa=sa[i], sb=sb[i])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

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

            offt = card.offt
            first_offt = offt[0]

            _x = card.x
            if _x is None:
                _x = [Defaults.default_double, Defaults.default_double, Defaults.default_double]
                _flag = 2
                _g0 = card.g0
            elif first_offt == 'B':
                _g0 = Defaults.default_int
                _flag = 0
            elif first_offt == 'G':
                _g0 = Defaults.default_int
                _flag = 1

            if _g0 is None:
                _g0 = Defaults.default_int

            # TODO: offt is not stored, need to convert wa and wb to basic coordinate system

            eid[i] = card.eid
            pid[i] = card.pid
            ga[i], gb[i] = card.node_ids
            sa[i] = card.sa
            sb[i] = card.sb
            x[i] = _x
            g0[i] = _g0
            f[i] = _flag
            pa[i] = card.pa
            pb[i] = card.pb
            wa[i] = card.wa
            wb[i] = card.wb

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class CBEAM3(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CBEAM3')


########################################################################################################################


class CBEND(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CBEND')

    def to_bdf(self, bdf):
        add_card = bdf.add_cbend

        data = self.identity
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

        for i in range(eid.size):
            _go = go[i]
            if _go == DataHelper.default_int:
                _go = None
            if x1[i] == DataHelper.default_double:
                x = None
            else:
                x = [x1[i], x2[i], x3[i]]

            add_card(eid[i], pid[i], [ga[i], gb[i]], _go, x, geom[i])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

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
            flag[i] = Defaults.unknown_int  # TODO: CBEND flag
            x1[i], x2[i], x3[i] = card.x
            go[i] = card.g0
            geom[i] = card.geom

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class CBUSH(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CBUSH')

    def to_bdf(self, bdf):
        add_card = bdf.add_cbush

        data = self.identity

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

        defaults = self._h5n.defaults

        for i in range(eid.size):
            _go = go[i]
            if _go == defaults.default_int:
                _go = None
            if x1[i] == defaults.default_double:
                x = None
            else:
                x = [x1[i], x2[i], x3[i]]
            if cid[i] == defaults.default_int:
                _cid = None
            else:
                _cid = cid[i]

            if s1[i] == defaults.default_double:
                _s = None
            else:
                _s = [s1[i], s2[i], s3[i]]

            add_card(eid[i], pid[i], [ga[i], gb[i]], x, _go, _cid, s[i], ocid[i], _s)

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

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
                _x = [Defaults.default_double, Defaults.default_double, Defaults.default_double]

            _g0 = card.g0
            if _g0 is None:
                _g0 = Defaults.default_int

            eid[i] = card.eid
            pid[i] = card.pid
            ga[i], gb[i] = card.node_ids
            flag[i] = Defaults.unknown_int  # TODO: CBUSH flag
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


class CBUSH1D(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CBUSH1D')
    
    def to_bdf(self, bdf):
        add_card = bdf.add_cbush1d
        
        data = self.identity
        
        eid = data['EID']
        pid = data['PID']
        g = data['G']
        cid = data['CID']
        
        for i in range(eid.size):
            add_card(eid[i], pid[i], g[i], cid[i])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

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


class CBUSH2D(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CBUSH2D')

    def to_bdf(self, bdf):
        add_card = bdf.add_cbush2d

        data = self.identity

        eid = data['EID']
        pid = data['PID']
        g = data['G']
        cid = data['CID']
        plane = data['PLANE']
        sptid = data['SPTID']

        for i in range(eid.size):
            add_card(eid[i], pid[i], g[i], cid[i], plane[i], sptid[i])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

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


class CCONEAX(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CCONEAX')

    def to_bdf(self, bdf):
        add_card = bdf.add_cconeax

        data = self.identity

        eid = data['EID']
        pid = data['PID']
        ra = data['RA']
        rb = data['RB']

        for i in range(eid.size):
            add_card(eid[i], pid[i], [ra[i], rb[i]])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

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


class CDAMP1(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CDAMP1')

    def to_bdf(self, bdf):
        add_card = bdf.add_cdamp1

        data = self.identity

        eid = data['EID']
        pid = data['PID']
        g1 = data['G1']
        g2 = data['G2']
        c1 = data['C1']
        c2 = data['C2']

        for i in range(eid.size):
            add_card(eid[i], pid[i], [g1[i], g2[i]], c1[i], c2[i])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        eid = data['EID']
        pid = data['PID']
        g1 = data['G1']
        g2 = data['G2']
        c1 = data['C1']
        c2 = data['C2']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            eid[i] = card.eid
            pid[i] = card.pid
            g1[i], g2[i] = card.node_ids
            c1[i] = card.c1
            c2[i] = card.c2

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class CDAMP2(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CDAMP2')

    def to_bdf(self, bdf):
        add_card = bdf.add_cdamp2

        data = self.identity

        eid = data['EID']
        b = data['B']
        g1 = data['G1']
        g2 = data['G2']
        c1 = data['C1']
        c2 = data['C2']

        for i in range(eid.size):
            add_card(eid[i], b[i], [g1[i], g2[i]], c1[i], c2[i])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)
    
        eid = data['EID']
        b = data['B']
        g1 = data['G1']
        g2 = data['G2']
        c1 = data['C1']
        c2 = data['C2']
    
        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]
    
            eid[i] = card.eid
            b[i] = card.b
            g1[i], g2[i] = card.node_ids
            c1[i] = card.c1
            c2[i] = card.c2
    
        result = {
            'IDENTITY': data
        }
    
        return result


########################################################################################################################


class CDAMP3(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CDAMP3')

    def to_bdf(self, bdf):
        add_card = bdf.add_cdamp3

        data = self.identity

        eid = data['EID']
        pid = data['PID']
        s1 = data['S1']
        s2 = data['S2']

        for i in range(eid.size):
            _s = [s1[i], s2[i]]
            add_card(eid[i], pid[i], _s)

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        eid = data['EID']
        pid = data['PID']
        s1 = data['S1']
        s2 = data['S2']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            nids = [nid if nid else 0 for nid in card.node_ids]

            eid[i] = card.eid
            pid[i] = card.pid
            s1[i], s2[i] = nids

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class CDAMP4(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CDAMP4')

    def to_bdf(self, bdf):
        add_card = bdf.add_cdamp4

        data = self.identity

        eid = data['EID']
        b = data['B']
        s1 = data['S1']
        s2 = data['S2']

        for i in range(eid.size):
            _s = [s1[i], s2[i]]
            add_card(eid[i], b[i], _s)

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)
    
        eid = data['EID']
        b = data['B']
        s1 = data['S1']
        s2 = data['S2']
    
        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]
    
            nids = [nid if nid else 0 for nid in card.node_ids]
    
            eid[i] = card.eid
            b[i] = card.b
            s1[i], s2[i] = nids
    
        result = {
            'IDENTITY': data
        }
    
        return result


########################################################################################################################


class CDAMP5(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CDAMP5')

    def to_bdf(self, bdf):
        add_card = bdf.add_cdamp5

        data = self.identity

        eid = data['EID']
        pid = data['PID']
        s1 = data['S1']
        s2 = data['S2']

        for i in range(eid.size):
            _s = [s1[i], s2[i]]
            add_card(eid[i], pid[i], _s)

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        eid = data['EID']
        pid = data['PID']
        s1 = data['S1']
        s2 = data['S2']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            nids = [nid if nid else 0 for nid in card.node_ids]

            eid[i] = card.eid
            pid[i] = card.pid
            s1[i], s2[i] = nids

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class CELAS1(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CELAS1')

    def to_bdf(self, bdf):
        add_card = bdf.add_celas1

        data = self.identity

        eid = data['EID']
        pid = data['PID']
        g1 = data['G1']
        g2 = data['G2']
        c1 = data['C1']
        c2 = data['C2']

        for i in range(eid.size):
            add_card(eid[i], pid[i], [g1[i], g2[i]], c1[i], c2[i])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        eid = data['EID']
        pid = data['PID']
        g1 = data['G1']
        g2 = data['G2']
        c1 = data['C1']
        c2 = data['C2']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            nids = [nid if nid is not None else 0 for nid in card.node_ids]

            eid[i] = card.eid
            pid[i] = card.pid
            g1[i], g2[i] = nids
            c1[i] = card.c1
            c2[i] = card.c2

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class CELAS2(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CELAS2')

    def to_bdf(self, bdf):
        add_card = bdf.add_celas2

        data = self.identity

        eid = data['EID']
        k = data['K']
        g1 = data['G1']
        g2 = data['G2']
        c1 = data['C1']
        c2 = data['C2']
        ge = data['GE']
        s = data['S']

        for i in range(eid.size):
            add_card(eid[i], k[i], [g1[i], g2[i]], c1[i], c2[i], ge[i], s[i])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        eid = data['EID']
        k = data['K']
        g1 = data['G1']
        g2 = data['G2']
        c1 = data['C1']
        c2 = data['C2']
        ge = data['GE']
        s = data['S']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            nids = [nid if nid is not None else 0 for nid in card.node_ids]

            eid[i] = card.eid
            k[i] = card.k
            g1[i], g2[i] = nids
            c1[i] = card.c1
            c2[i] = card.c2
            ge[i] = card.ge
            s[i] = card.s

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class CELAS3(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CELAS3')

    def to_bdf(self, bdf):
        add_card = bdf.add_celas3

        data = self.identity

        eid = data['EID']
        pid = data['PID']
        s1 = data['S1']
        s2 = data['S2']

        for i in range(eid.size):
            add_card(eid[i], pid[i], [s1[i], s2[i]])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        eid = data['EID']
        pid = data['PID']
        s1 = data['S1']
        s2 = data['S2']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            nids = [nid if nid is not None else 0 for nid in card.node_ids]

            eid[i] = card.eid
            pid[i] = card.pid
            s1[i], s2[i] = nids

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class CELAS4(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CELAS4')

    def to_bdf(self, bdf):
        add_card = bdf.add_celas4

        data = self.identity

        eid = data['EID']
        k = data['K']
        s1 = data['S1']
        s2 = data['S2']

        for i in range(eid.size):
            add_card(eid[i], k[i], [s1[i], s2[i]])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        eid = data['EID']
        k = data['K']
        s1 = data['S1']
        s2 = data['S2']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            nids = [nid if nid is not None else 0 for nid in card.node_ids]

            eid[i] = card.eid
            k[i] = card.k
            s1[i], s2[i] = nids

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class CFAST(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CFAST')


########################################################################################################################


class CFASTP(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CFASTP')


########################################################################################################################


class CFLUID2(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CFLUID2')


########################################################################################################################


class CFLUID3(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CFLUID3')


########################################################################################################################


class CFLUID4(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CFLUID4')


########################################################################################################################


class CGAP(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CGAP')
    
    def to_bdf(self, bdf):
        add_card = bdf.add_cgap
        
        data = self.identity

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

        defaults = self._h5n.defaults
        
        for i in range(eid.size):
            _x = defaults.to_list_double([x1[i], x2[i], x3[i]])
            _go = defaults.to_value_int(go[i])
            _cid = defaults.to_value_int(cid[i])

            add_card(eid[i], pid[i], [ga[i], gb[i]], _x, _go, _cid)

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

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
        
        defaults = self._h5n.defaults

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            eid[i] = card.eid
            pid[i] = card.pid
            ga[i], gb[i] = card.node_ids
            f[i] = defaults.unknown_int  # TODO: CGAP flag
            x1[i], x2[i], x3[i] = defaults.get_list_double(card.x)
            go[i] = defaults.get_value_int(card.g0)
            cid[i] = defaults.get_value_int(card.cid)

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class CHACAB(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CHACAB')


########################################################################################################################


class CHACBR(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CHACBR')


########################################################################################################################


class CHBDYE(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CHBDYE')


########################################################################################################################


class CHBDYG(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CHBDYG')


########################################################################################################################


class CHBDYP(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CHBDYP')


########################################################################################################################


class CHEXA(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CHEXA')

    def to_bdf(self, bdf):
        add_card = bdf.add_chexa

        data = self.identity

        eid = data['EID']
        pid = data['PID']
        g = data['G']

        for i in range(eid.size):
            nids = g[i].tolist()
            try:
                nids.remove(0)
            except ValueError:
                pass
            add_card(eid[i], pid[i], nids)

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

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


class CHEXAL(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CHEXAL')


########################################################################################################################


class CHEXP(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CHEXP')


########################################################################################################################


class CIFHEX(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CIFHEX')


########################################################################################################################


class CIFPENT(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CIFPENT')


########################################################################################################################


class CIFQDX(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CIFQDX')


########################################################################################################################


class CIFQUAD(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CIFQUAD')


########################################################################################################################


class CMASS1(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CMASS1')


########################################################################################################################


class CMASS2(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CMASS2')


########################################################################################################################


class CMASS3(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CMASS3')


########################################################################################################################


class CMASS4(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CMASS4')

    def to_bdf(self, bdf):
        add_card = bdf.add_cmass4

        data = self.identity

        eid = data['EID']
        m = data['M']
        s1 = data['S1']
        s2 = data['S2']

        for i in range(eid.size):
            add_card(eid[i], m[i], [s1[i], s2[i]])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        eid = data['EID']
        m = data['M']
        s1 = data['S1']
        s2 = data['S2']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            eid[i] = card.eid
            m[i] = card.mass
            s1[i], s2[i] = card.nodes

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class CONM2(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CONM2')

    def to_bdf(self, bdf):
        add_card = bdf.add_conm2

        data = self.identity

        eid = data['EID']
        g = data['G']
        cid = data['CID']
        m = data['M']
        x1 = data['X1']
        x2 = data['X2']
        x3 = data['X3']
        i1 = data['I1']
        i2 = data['I2']
        i3 = data['I3']

        for i in range(eid.size):
            x = [x1[i], x2[i], x3[i]]
            i_ = [i1[i]]
            i_.extend(i2[i].tolist())
            i_.extend(i3[i].tolist())
            add_card(eid[i], g[i], m[i], cid[i], x, i_)

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        eid = data['EID']
        g = data['G']
        cid = data['CID']
        m = data['M']
        x1 = data['X1']
        x2 = data['X2']
        x3 = data['X3']
        i1 = data['I1']
        i2 = data['I2']
        i3 = data['I3']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            eid[i] = card.eid
            g[i] = card.nid
            cid[i] = card.cid
            m[i] = card.mass
            x1[i], x2[i], x3[i] = card.X
            i1[i] = card.I[0]
            i2[i] = card.I[1:3]
            i3[i] = card.I[3:]

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class CONROD(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CONROD')

    def to_bdf(self, bdf):
        add_card = bdf.add_conrod

        data = self.identity

        eid = data['EID']
        g1 = data['G1']
        g2 = data['G2']
        mid = data['MID']
        a = data['A']
        j = data['J']
        c = data['C']
        nsm = data['NSM']

        for i in range(eid.size):
            add_card(eid[i], mid[i], [g1[i], g2[i]], a[i], j[i], c[i], nsm[i])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        eid = data['EID']
        g1 = data['G1']
        g2 = data['G2']
        mid = data['MID']
        a = data['A']
        j = data['J']
        c = data['C']
        nsm = data['NSM']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            eid[i] = card.eid
            g1[i], g2[i] = card.node_ids
            mid[i] = card.mid
            a[i] = card.A
            j[i] = card.j
            c[i] = card.c
            nsm[i] = card.nsm

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class CONTRLT(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CONTRLT')


########################################################################################################################


class CPENTA(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CPENTA')

    def to_bdf(self, bdf):
        add_card = bdf.add_cpenta

        data = self.identity

        eid = data['EID']
        pid = data['PID']
        g = data['G']

        for i in range(eid.size):
            nids = [_ for _ in g[i] if _ != 0]
            add_card(eid[i], pid[i], nids)

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        eid = data['EID']
        pid = data['PID']
        g = data['G']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            eid[i] = card.eid
            pid[i] = card.pid

            nids = [nid if nid is not None else 0 for nid in card.node_ids]

            diff_len = 15 - len(nids)
            if diff_len > 0:
                nids += [0] * diff_len

            g[i] = nids

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class CQDX4FD(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CQDX4FD')


########################################################################################################################


class CQDX9FD(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CQDX9FD')


########################################################################################################################


class CQUAD(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CQUAD')


########################################################################################################################


class CQUAD4(InputTable, QuadShell):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CQUAD4')

    def to_bdf(self, bdf):
        add_card = bdf.add_cquad4

        data = self.identity

        eid = data['EID']
        pid = data['PID']
        g = data['G']
        theta = data['THETA']
        zoffs = data['ZOFFS']
        tflag = data['TFLAG']
        t = data['T']
        mcid = data['MCID']

        def _get_val(val, default):
            if val == default:
                return None
            return val

        default_int = self._h5n.defaults.default_int
        default_double = self._h5n.defaults.default_double

        for i in range(len(eid)):
            _eid = eid[i]
            _pid = pid[i]
            _g = g[i].tolist()
            _theta = _get_val(theta[i], default_double)
            _zoffs = zoffs[i]
            _tflag = tflag[i]
            _t = t[i]
            _mcid = _get_val(mcid[i], default_int)

            if _theta is None:
                theta_mcid = _mcid
            else:
                theta_mcid = _theta

            add_card(_eid, _pid, _g, theta_mcid, _zoffs, _tflag, _t[0], _t[1], _t[2], _t[3])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        eid = data['EID']
        pid = data['PID']
        g = data['G']
        theta = data['THETA']
        zoffs = data['ZOFFS']
        tflag = data['TFLAG']
        t = data['T']
        mcid = data['MCID']

        default_int = self._h5n.defaults.default_int
        default_double = self._h5n.defaults.default_double

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            theta_mcid = card.theta_mcid
            if not isinstance(theta_mcid, float):
                _mcid = theta_mcid
                _theta = default_double
            else:
                _mcid = default_int
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


class CQUAD4FD(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CQUAD4FD')


########################################################################################################################


class CQUAD8(InputTable, QuadShell):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CQUAD8')

    def to_bdf(self, bdf):
        add_card = bdf.add_cquad8

        data = self.identity

        eid = data['EID']
        pid = data['PID']
        g = data['G']
        theta = data['THETA']
        zoffs = data['ZOFFS']
        tflag = data['TFLAG']
        t = data['T']
        mcid = data['MCID']

        defaults = self._h5n.defaults

        get_mcid = defaults.to_value_int
        get_theta = defaults.to_value_double

        for i in range(eid.size):
            _mcid = get_mcid(mcid[i])
            _theta = get_theta(theta[i])
            if _mcid is None:
                theta_mcid = _theta
            else:
                theta_mcid = _mcid
            _t = t[i]
            add_card(eid[i], pid[i], g[i], theta_mcid, zoffs[i], tflag[i], _t[0], _t[1], _t[2], _t[3])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

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
                _theta = Defaults.default_double
            else:
                _mcid = Defaults.default_int
                _theta = theta_mcid

            nids = [nid if nid is not None else 0 for nid in card.node_ids]
            if len(nids) < 8:
                nids += [0] * (8 - len(nids))

            eid[i] = card.eid
            pid[i] = card.pid
            g[i] = nids
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


class CQUAD9FD(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CQUAD9FD')


########################################################################################################################


class CQUADR(InputTable, QuadShell):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CQUADR')
    
    def to_bdf(self, bdf):
        add_card = bdf.add_cquadr
        
        data = self.identity
        
        eid = data['EID']
        pid = data['PID']
        g = data['G']
        theta = data['THETA']
        zoffs = data['ZOFFS']
        tflag = data['TFLAG']
        t = data['T']
        mcid = data['MCID']
        
        defaults = self._h5n.defaults

        get_mcid = defaults.to_value_int
        get_theta = defaults.to_value_double

        for i in range(eid.size):
            _mcid = get_mcid(mcid[i])
            _theta = get_theta(theta[i])
            if _mcid is None:
                theta_mcid = _theta
            else:
                theta_mcid = _mcid
            _t = t[i]
            add_card(eid[i], pid[i], g[i], theta_mcid, zoffs[i], tflag[i], _t[0], _t[1], _t[2], _t[3])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

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
                _theta = Defaults.default_double
            else:
                _mcid = Defaults.default_int
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


class CQUADX(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CQUADX')


########################################################################################################################


class CRBE1(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CRBE1/IDENTITY',
                                rename={'LAGMULTIPL_POS': 'LAGM_POS', 'LAGMULTIPL_LEN': 'LAGM_LEN'})


########################################################################################################################


class CROD(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CROD')

    def to_bdf(self, bdf):
        add_card = bdf.add_crod

        data = self.identity

        eid = data['EID']
        pid = data['PID']
        g = data['G']

        for i in range(eid.size):
            add_card(eid[i], pid[i], g[i])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

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


class CSEAM(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CSEAM')


########################################################################################################################


class CSEAMP(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CSEAMP')


########################################################################################################################


class CSHEAR(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CSHEAR')

    def to_bdf(self, bdf):
        add_card = bdf.add_cshear

        data = self.identity

        eid = data['EID']
        pid = data['PID']
        g = data['G']

        for i in range(eid.size):
            add_card(eid[i], pid[i], g[i])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

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


class CSLOT3(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CSLOT3')


########################################################################################################################


class CSLOT4(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CSLOT4')


########################################################################################################################


class CTETRA(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTETRA')

    def to_bdf(self, bdf):
        add_card = bdf.add_ctetra

        data = self.identity

        eid = data['EID']
        pid = data['PID']
        g = data['G']

        for i in range(eid.size):
            nids = [_ for _ in g[i] if _ != 0]
            add_card(eid[i], pid[i], nids)

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        eid = data['EID']
        pid = data['PID']
        g = data['G']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            eid[i] = card.eid
            pid[i] = card.pid

            nids = [nid if nid is not None else 0 for nid in card.node_ids]

            diff_len = 10 - len(nids)
            if diff_len > 0:
                nids += [0] * diff_len

            g[i] = nids

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class CTRIA3(InputTable, TriaShell):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTRIA3')

    def to_bdf(self, bdf):
        add_card = bdf.add_ctria3

        data = self.identity

        eid = data['EID']
        pid = data['PID']
        g = data['G']
        theta = data['THETA']
        zoffs = data['ZOFFS']
        tflag = data['TFLAG']
        t = data['T']
        mcid = data['MCID']

        def _get_val(val, default):
            if val == default:
                return None
            return val

        for i in range(len(eid)):
            _eid = eid[i]
            _pid = pid[i]
            _g = g[i].tolist()
            _theta = _get_val(theta[i], Defaults.default_double)
            _zoffs = zoffs[i]
            _tflag = tflag[i]
            _t = t[i]
            _mcid = _get_val(mcid[i], Defaults.default_int)

            if _theta is None:
                theta_mcid = _mcid
            else:
                theta_mcid = _theta

            add_card(_eid, _pid, _g, theta_mcid, _zoffs, _tflag, _t[0], _t[1], _t[2])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

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
                _theta = Defaults.default_double
            else:
                _mcid = Defaults.default_int
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


class CTRIA3FD(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTRIA3FD')


########################################################################################################################


class CTRIA6(InputTable, TriaShell):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTRIA6')

    def to_bdf(self, bdf):
        add_card = bdf.add_ctria6

        data = self.identity

        eid = data['EID']
        pid = data['PID']
        g = data['G']
        theta = data['THETA']
        zoffs = data['ZOFFS']
        tflag = data['TFLAG']
        t = data['T']
        mcid = data['MCID']

        defaults = self._h5n.defaults

        get_mcid = defaults.to_value_int
        get_theta = defaults.to_value_double

        for i in range(eid.size):
            _mcid = get_mcid(mcid[i])
            _theta = get_theta(theta[i])
            if _mcid is None:
                theta_mcid = _theta
            else:
                theta_mcid = _mcid
            _t = t[i]
            add_card(eid[i], pid[i], g[i], theta_mcid, zoffs[i], tflag[i], _t[0], _t[1], _t[2])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        eid = data['EID']
        pid = data['PID']
        g = data['G']
        theta = data['THETA']
        zoffs = data['ZOFFS']
        t = data['T']
        tflag = data['TFLAG']
        mcid = data['MCID']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            theta_mcid = card.theta_mcid
            if not isinstance(theta_mcid, float):
                _mcid = theta_mcid
                _theta = Defaults.default_double
            else:
                _mcid = Defaults.default_int
                _theta = theta_mcid

            nids = [nid if nid else 0 for nid in card.node_ids]

            eid[i] = card.eid
            pid[i] = card.pid
            g[i] = nids
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


class CTRIA6FD(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTRIA6FD')


########################################################################################################################


class CTRIAH(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTRIAH')


########################################################################################################################


class CTRIAR(InputTable, TriaShell):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTRIAR')

    def to_bdf(self, bdf):
        add_card = bdf.add_ctriar

        data = self.identity

        eid = data['EID']
        pid = data['PID']
        g = data['G']
        theta = data['THETA']
        zoffs = data['ZOFFS']
        tflag = data['TFLAG']
        t = data['T']
        mcid = data['MCID']

        defaults = self._h5n.defaults

        get_mcid = defaults.to_value_int
        get_theta = defaults.to_value_double

        for i in range(eid.size):
            _mcid = get_mcid(mcid[i])
            _theta = get_theta(theta[i])
            if _mcid is None:
                theta_mcid = _theta
            else:
                theta_mcid = _mcid
            _t = t[i]
            add_card(eid[i], pid[i], g[i], theta_mcid, zoffs[i], tflag[i], _t[0], _t[1], _t[2])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

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
                _theta = Defaults.default_double
            else:
                _mcid = Defaults.default_int
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


class CTRIAX(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTRIAX')


########################################################################################################################


class CTRIAX6(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTRIAX6')


########################################################################################################################


class CTRIX3FD(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTRIX3FD')


########################################################################################################################


class CTRIX6FD(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTRIX6FD')


########################################################################################################################


class CTUBE(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CTUBE')

    def to_bdf(self, bdf):
        add_card = bdf.add_ctube

        data = self.identity

        eid = data['EID']
        pid = data['PID']
        g = data['G']

        for i in range(eid.size):
            add_card(eid[i], pid[i], g[i])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        eid = data['EID']
        pid = data['PID']
        g = data['G']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            nids = [nid if nid else 0 for nid in card.node_ids]

            eid[i] = card.eid
            pid[i] = card.pid
            g[i] = nids

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class CVISC(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CVISC')

    def to_bdf(self, bdf):
        add_card = bdf.add_cvisc

        data = self.identity

        eid = data['EID']
        pid = data['PID']
        g = data['G']

        for i in range(eid.size):
            add_card(eid[i], pid[i], g[i])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        eid = data['EID']
        pid = data['PID']
        g = data['G']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            nids = [nid if nid else 0 for nid in card.node_ids]

            eid[i] = card.eid
            pid[i] = card.pid
            g[i] = nids

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class CWELD(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CWELD')


########################################################################################################################


class CWELDC(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CWELDC')


########################################################################################################################


class CWELDP(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/CWELDP')


########################################################################################################################


# TODO: GENEL doesn't conform
# class GENEL(CardTable):
#     table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/GENEL/IDENTITY',
#                                 rename={'UDLIST_POS': 'UD_POS', 'UDLIST_LEN': 'UD_LEN', 'UILIST_POS': 'UI_POS',
#                                         'UILIST_LEN': 'UI_LEN'})

########################################################################################################################


class PLOTEL(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/PLOTEL')

    def to_bdf(self, bdf):
        add_card = bdf.add_plotel

        data = self.identity

        eid = data['EID']
        g = data['G']

        for i in range(eid.size):
            add_card(eid[i], g[i])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        eid = data['EID']
        g = data['G']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            nids = [nid if nid else 0 for nid in card.node_ids]

            eid[i] = card.eid
            g[i] = nids

        result = {
            'IDENTITY': data
        }

        return result


########################################################################################################################


class PRIM1(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/PRIM1')


########################################################################################################################


class PRIM2(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/PRIM2')


########################################################################################################################


class PRIM3(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/PRIM3')


########################################################################################################################


class PRIM4(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/PRIM4')


########################################################################################################################


class PRIM5(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/PRIM5')


########################################################################################################################


class PRIM6(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/PRIM6')


########################################################################################################################


class PRIM7(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/PRIM7')


########################################################################################################################


class PRIM8(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/PRIM8')


########################################################################################################################


class RADCOL(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/RADCOL')


########################################################################################################################


class RBAR(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/RBAR')

    def to_bdf(self, bdf):
        add_card = bdf.add_rbar

        data = self.identity

        eid = data['EID']
        ga = data['GA']
        gb = data['GB']
        cna = data['CNA']
        cnb = data['CNB']
        cma = data['CMA']
        cmb = data['CMB']
        alpha = data['ALPHA']

        def _get_val(val):
            return val if val != 0 else None

        for i in range(eid.size):
            _cna = _get_val(cna[i])
            _cnb = _get_val(cnb[i])
            _cma = _get_val(cma[i])
            _cmb = _get_val(cmb[i])

            add_card(eid[i], [ga[i], gb[i]], _cna, _cnb, _cma, _cmb, alpha[i])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        eid = data['EID']
        ga = data['GA']
        gb = data['GB']
        cna = data['CNA']
        cnb = data['CNB']
        cma = data['CMA']
        cmb = data['CMB']
        alpha = data['ALPHA']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            eid[i] = card.eid
            ga[i] = card.ga
            gb[i] = card.gb
            cna[i] = card.cna if card.cna != '' else 0
            cnb[i] = card.cnb if card.cnb != '' else 0
            cma[i] = card.cma if card.cma != '' else 0
            cmb[i] = card.cmb if card.cmb != '' else 0
            alpha[i] = card.alpha

        return {'IDENTITY': data}


########################################################################################################################


class RBAR1(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/RBAR1')


########################################################################################################################


class RBE1(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/RBE1/IDENTITY')


########################################################################################################################


class RBE2(InputTable):
    table_def = TableDef.create(
        '/NASTRAN/INPUT/ELEMENT/RBE2/RB',
        subtables=[TableDef.create('/NASTRAN/INPUT/ELEMENT/RBE2/GM')],
    )

    def to_bdf(self, bdf):
        add_card = bdf.add_rbe2

        identity = self.identity
        gm = self.gm

        eid = identity['EID']
        gn = identity['GN']
        cm = identity['CM']
        gm_pos = identity['GM_POS']
        gm_len = identity['GM_LEN']
        alpha = identity['ALPHA']
        gm_id = gm['ID']

        j = 0
        for i in range(eid.size):
            gms = gm_id[j:j+gm_len[i]]
            j += gm_len[i]
            add_card(eid[i], gn[i], cm[i], gms, alpha[i])

    def from_bdf(self, cards):
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


class RBE2GS(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/RBE2GS/IDENTITY')


########################################################################################################################


class RBE3(InputTable):
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

    def to_bdf(self, bdf):
        add_card = bdf.add_rbe3

        identity = self.identity
        _wtcg = self.wtcg
        _gm = self.gm
        _g = self.g

        id_ = _g['ID']
        gm = _gm['GM']
        cm = _gm['CM']
        wt1 = _wtcg['WT1']
        c = _wtcg['C']
        g_pos = _wtcg['G_POS']
        g_len = _wtcg['G_LEN']
        eid = identity['EID']
        refg = identity['REFG']
        refc = identity['REFC']
        wtcg_pos = identity['WTCG_POS']
        wtcg_len = identity['WTCG_LEN']
        gm_pos = identity['GM_POS']
        gm_len = identity['GM_LEN']
        alpha = identity['ALPHA']

        for i in range(eid.size):
            i1 = gm_pos[i]
            i2 = i1 + gm_len[i]
            gmi = gm[i1:i2]
            cmi = cm[i1:i2]

            i1 = wtcg_pos[i]
            i2 = i1 + wtcg_len[i]
            weights = wt1[i1:i2]
            comps = c[i1:i2]

            gijs = []

            for j in range(i1, i2):
                j1 = g_pos[j]
                j2 = j1 + g_len[j]
                gijs.append(id_[j1:j2])

            add_card(eid[i], refg[i], refc[i], weights, comps, gijs, gmi, cmi, alpha[i])

    def from_bdf(self, cards):
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


class RINGAX(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/RINGAX')


########################################################################################################################


class RJOINT(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/RJOINT')


########################################################################################################################


class RROD(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/RROD')


########################################################################################################################


class RSSCON(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/RSSCON')


########################################################################################################################


class RTRPLT(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/RTRPLT')


########################################################################################################################


class RTRPLT1(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/RTRPLT1')


########################################################################################################################


class SECTAX(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/SECTAX')


########################################################################################################################


class WETELME(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/WETELME')


########################################################################################################################


class WETELMG(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/ELEMENT/WETELMG/IDENTITY')

########################################################################################################################
