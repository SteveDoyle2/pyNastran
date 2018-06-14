from __future__ import print_function, absolute_import
from six import iteritems
from collections import OrderedDict

from h5Nastran.h5nastrannode import H5NastranNode
from .input_table import InputTable, TableDef


class Constraint(H5NastranNode):
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

    def to_bdf(self, bdf):
        for key, item in iteritems(self.__dict__):
            if key.startswith('_'):
                continue
            try:
                item.to_bdf(bdf)
            except NotImplementedError:
                pass
    

########################################################################################################################


class AELINK(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/AELINK/IDENTITY')

########################################################################################################################


class CSSCHD(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/CSSCHD')

########################################################################################################################


class CYSUP(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/CYSUP')

########################################################################################################################


class DEFORM(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/DEFORM')

########################################################################################################################


class GRDSET(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/GRDSET')

########################################################################################################################


class MPC(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/MPC/IDENTITY')

    def to_bdf(self, bdf):
        add_card = bdf.add_mpc

        identity = self.identity
        gca = self.gca

        _g = gca['G']
        _c = gca['C']
        _a = gca['A']
        sid = identity['SID']
        g = identity['G']
        c = identity['C']
        a = identity['A']
        gca_pos = identity['GCA_POS']
        gca_len = identity['GCA_LEN']

        for i in range(sid.size):
            j1 = gca_pos[i]
            j2 = j1 + gca_len[i]

            nids = [g[i]] + _g[j1:j2].tolist()
            comps = [c[i]] + _c[j1:j2].tolist()
            enforced = [a[i]] + _a[j1:j2].tolist()

            add_card(sid[i], nids, comps, enforced)

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        gca = {'IDENTITY': {'G': [], 'C': [], 'A': []}}

        result = {
            'IDENTITY': {'SID': [], 'G': [], 'C': [], 'A': [], 'GCA_POS': [], 'GCA_LEN': [], 'DOMAIN_ID': []},
            'GCA': gca,
            '_subtables': ['GCA']
        }

        gca = gca['IDENTITY']

        _g = gca['G']
        _c = gca['C']
        _a = gca['A']
        identity = result['IDENTITY']
        sid = identity['SID']
        g = identity['G']
        c = identity['C']
        a = identity['A']
        gca_pos = identity['GCA_POS']
        gca_len = identity['GCA_LEN']

        pos = 0

        for card_id in card_ids:
            card_list = sorted(cards[card_id], key=lambda i: i.nodes[0])

            for card in card_list:
                sid.append(card.conid)
                nodes = list(card.nodes)
                g1 = nodes[0]
                nodes = nodes[1:]
                comps = list(card.components)
                c1 = comps[0]
                comps = comps[1:]
                coefs = list(card.coefficients)
                a1 = coefs[0]
                coefs = coefs[1:]
                g.append(g1)
                c.append(c1)
                a.append(a1)
                gca_pos.append(pos)
                _len = len(nodes)
                gca_len.append(_len)
                pos += _len
                _g.extend(nodes)
                _c.extend(comps)
                _a.extend(coefs)

        return result


########################################################################################################################


class MPCADD(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/MPCADD/IDENTITY')

########################################################################################################################


class MPCAX(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/MPCAX/IDENTITY')

########################################################################################################################


class MPCD(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/MPCD')

########################################################################################################################


class MPCY(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/MPCY/IDENTITY')

########################################################################################################################


class RSPLINE(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/RSPLINE/IDENTITY')

########################################################################################################################


class SESUP(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SESUP')

########################################################################################################################


class SPBLND1(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPBLND1')

########################################################################################################################


class SPBLND2(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPBLND2')

########################################################################################################################


class SPC(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPC')

    def to_bdf(self, bdf):
        add_card = bdf.add_spc

        identity = self.identity

        sid = identity['SID']
        g = identity['G']
        c = identity['C']
        d = identity['D']

        _spcs = OrderedDict()
        for i in range(sid.size):
            try:
                _spc = _spcs[sid[i]]
            except KeyError:
                _spc = _spcs[sid[i]] = [[], [], []]

            _gids, _comps, _enf = _spc

            _gids.append(g[i])
            _comps.append(c[i])
            _enf.append(d[i])

        for key, _spc in iteritems(_spcs):
            _gids, _comps, _enf = _spc
            add_card(key, _gids, _comps, _enf)

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        result = {
            'IDENTITY': {'SID': [], 'G': [], 'C': [], 'D': [], 'DOMAIN_ID': []},
        }

        identity = result['IDENTITY']

        sid = identity['SID']
        g = identity['G']
        c = identity['C']
        d = identity['D']

        for card_id in card_ids:
            card_list = sorted(cards[card_id], key=lambda x: x.conid)
            for card in card_list:
                sid += [card.conid] * len(card.gids)
                g += card.gids
                c += card.components
                d += card.enforced

        return result

########################################################################################################################


class SPC1_G(InputTable):
    card_id = 'SPC1'
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPC1/SPC1_G/IDENTITY')
    
    def from_bdf(self, cards):
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


class SPC1_THRU(InputTable):
    card_id = None  # don't register this
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPC1/SPC1_THRU')

########################################################################################################################


class SPCADD(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPCADD/IDENTITY')

    def to_bdf(self, bdf):
        add_card = bdf.add_spcadd

        identity = self.identity
        sid = identity['SID']
        s_pos = identity['S_POS']
        s_len = identity['S_LEN']
        s = self.s['S']

        for i in range(sid.size):
            j1 = s_pos[i]
            j2 = j1 + s_len[i]
            add_card(sid[i], s[j1:j2].tolist())

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
    
        s = {
            'IDENTITY': {'S': []}
        }
    
        result = {
            'IDENTITY': {'SID': [], 'S_POS': [], 'S_LEN': [], 'DOMAIN_ID': []},
            'S': s,
            '_subtables': ['S']
        }
    
        identity = result['IDENTITY']
    
        sid = identity['SID']
        s_pos = identity['S_POS']
        s_len = identity['S_LEN']
        s = s['IDENTITY']['S']
    
        _pos = 0
        for card_id in card_ids:
            card_list = sorted(cards[card_id], key=lambda x: x.conid)
            for card in card_list:
                sid.append(card.conid)
                s_pos.append(_pos)
                _s_len = len(card.sets)
                s_len.append(_s_len)
                s += card.sets
    
        return result

########################################################################################################################


class SPCAX(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPCAX')

########################################################################################################################


class SPCD(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPCD')

    def to_bdf(self, bdf):
        add_card = bdf.add_spcd

        data = self.identity

        sid = data['SID']
        g = data['G']
        c = data['C']
        d = data['D']

        _spcds = OrderedDict()
        for i in range(sid.size):
            try:
                _spcd = _spcds[sid[i]]
            except KeyError:
                _spcd = _spcds[sid[i]] = [[], [], []]

            _gids, _constr, _enf = _spcd

            _gids.append(g[i])
            _constr.append(c[i])
            _enf.append(d[i])

        for key, _spcd in iteritems(_spcds):
            _gids, _constr, _enf = _spcd
            add_card(key, _gids, _constr, _enf)

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        result = {
            'IDENTITY': {'SID': [], 'G': [], 'C': [], 'D': [], 'DOMAIN_ID': []}
        }

        i = result['IDENTITY']
        sid = i['SID']
        g = i['G']
        c = i['C']
        d = i['D']

        for card_id in card_ids:
            card_list = cards[card_id]

            for card in card_list:
                sid.extend([card.sid] * len(card.nodes))
                g.extend(card.nodes)
                c.extend([int(_) if _ is not None else DataHelper.default_int for _ in card.constraints])
                d.extend(card.enforced)

        return result

########################################################################################################################


class SPCOFF(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPCOFF')


########################################################################################################################


class SPCOFF1(InputTable):
    table_def = TableDef.create(
        '/NASTRAN/INPUT/CONSTRAINT/SPCOFF1/IDENTITY',
        rename={'GIS_POS': 'G_POS', 'GIS_LEN': 'G_LEN'}
    )

########################################################################################################################


class SPCR(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPCR')

########################################################################################################################


class SPLINE1(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPLINE1')

########################################################################################################################


class SPLINE2(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPLINE2')

    def to_bdf(self, bdf):
        add_card = bdf.add_spline2

        identity = self.identity

        eid = identity['EID']
        caero = identity['CAERO']
        id1 = identity['ID1']
        id2 = identity['ID2']
        setg = identity['SETG']
        dz = identity['DZ']
        dtor = identity['DTOR']
        cid = identity['CID']
        dthx = identity['DTHX']
        dthy = identity['DTHY']
        usage = identity['USAGE']

        for i in range(eid.size):
            add_card(eid[i], caero[i], id1[i], id2[i], setg[i], dz[i], dtor[i], cid[i], dthx[i], dthy[i], usage[i])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        result = {'IDENTITY': {'EID': [], 'CAERO': [], 'ID1': [], 'ID2': [], 'SETG': [], 'DZ': [], 'DTOR': [],
                               'CID': [], 'DTHX': [], 'DTHY': [], 'USAGE': [], 'DOMAIN_ID': []}}

        identity = result['IDENTITY']
        eid = identity['EID']
        caero = identity['CAERO']
        id1 = identity['ID1']
        id2 = identity['ID2']
        setg = identity['SETG']
        dz = identity['DZ']
        dtor = identity['DTOR']
        cid = identity['CID']
        dthx = identity['DTHX']
        dthy = identity['DTHY']
        usage = identity['USAGE']

        for card_id in card_ids:
            card = cards[card_id]

            eid.append(card.eid)
            caero.append(card.caero)
            id1.append(card.id1)
            id2.append(card.id2)
            setg.append(card.setg)
            dz.append(card.dz)
            dtor.append(card.dtor)
            cid.append(card.cid)
            dthx.append(card.dthx)
            dthy.append(card.dthy)
            usage.append(card.usage)

        return result


########################################################################################################################


class SPLINE3(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPLINE3/IDENTITY')

########################################################################################################################


class SPLINE4(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPLINE4')

########################################################################################################################


class SPLINE5(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPLINE5')

########################################################################################################################


class SPLINE6(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPLINE6')

########################################################################################################################


class SPLINE7(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPLINE7')

########################################################################################################################


class SPLINEX(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPLINEX')

########################################################################################################################


class SPLINRB(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SPLINRB/IDENTITY')

########################################################################################################################


class SUPAX(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SUPAX')

########################################################################################################################


class SUPORT(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SUPORT')

########################################################################################################################


class SUPORT1(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/SUPORT1/IDENTITY')

########################################################################################################################


class TEMP(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/TEMP')

########################################################################################################################


class TEMPAX(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/TEMPAX')

########################################################################################################################


class TEMPB3(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/TEMPB3')

########################################################################################################################


class TEMPBC(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/TEMPBC')

########################################################################################################################


class TEMPD(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/TEMPD')

    def to_bdf(self, bdf):
        add_card = bdf.add_tempd

        data = self.identity

        sid = data['SID']
        t = data['T']

        for i in range(sid.size):
            add_card(sid[i], t[i])

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        result = {'IDENTITY': {'SID': [], 'T': [], 'DOMAIN_ID': []}}
        i = result['IDENTITY']
        sid = i['SID']
        t = i['T']

        for card_id in card_ids:
            card = cards[card_id]

            sid.append(card.sid)
            t.append(card.temperature)

        return result


########################################################################################################################


class TEMPN1(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/TEMPN1')

########################################################################################################################


class TEMPP1(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/TEMPP1')

########################################################################################################################


class TEMPP2(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/TEMPP2')

########################################################################################################################


class TEMPP3(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/TEMPP3')


########################################################################################################################


class TEMPRB(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/TEMPRB')

########################################################################################################################


class TRIM(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/TRIM/IDENTITY')

########################################################################################################################


class TRIM2(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/TRIM2/IDENTITY')

########################################################################################################################


class UXVEC(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/CONSTRAINT/UXVEC/IDENTITY')

########################################################################################################################
