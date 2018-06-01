from __future__ import print_function, absolute_import

import numpy as np
from six import iterkeys

from h5Nastran.defaults import Defaults
from h5Nastran.h5nastrannode import H5NastranNode
from .input_table import InputTable, TableDef


class Property(H5NastranNode):
    def __init__(self, h5n, input):
        self._h5n = h5n
        self._input = input

        self.mfluid = MFLUID(self._h5n, self)
        self.nsm = NSM(self._h5n, self)
        self.nsm1 = NSM1(self._h5n, self)
        self.nsmadd = NSMADD(self._h5n, self)
        self.nsml = NSML(self._h5n, self)
        self.nsml1 = NSML1(self._h5n, self)
        self.paabsf = PAABSF(self._h5n, self)
        self.pacabs = PACABS(self._h5n, self)
        self.pacbar = PACBAR(self._h5n, self)
        self.pacinf = PACINF(self._h5n, self)
        self.paero1 = PAERO1(self._h5n, self)
        self.paero2 = PAERO2(self._h5n, self)
        self.paero3 = PAERO3(self._h5n, self)
        self.paero4 = PAERO4(self._h5n, self)
        self.paero5 = PAERO5(self._h5n, self)
        self.paxisym = PAXISYM(self._h5n, self)
        self.paxsymh = PAXSYMH(self._h5n, self)
        self.pbar = PBAR(self._h5n, self)
        self.pbarl = PBARL(self._h5n, self)
        self.pbarn1 = PBARN1(self._h5n, self)
        self.pbcomp = PBCOMP(self._h5n, self)
        self.pbeam = PBEAM(self._h5n, self)
        self.pbeam3 = PBEAM3(self._h5n, self)
        self.pbeaml = PBEAML(self._h5n, self)
        self.pbemn1 = PBEMN1(self._h5n, self)
        self.pbend = PBEND(self._h5n, self)
        # self.pbmsect = PBMSECT(self._h5n, self)
        # self.pbrsect = PBRSECT(self._h5n, self)
        self.pbush = PBUSH(self._h5n, self)
        self.pbush1d = PBUSH1D(self._h5n, self)
        # self.pbush2d = PBUSH2D(self._h5n, self)
        self.pbusht = PBUSHT(self._h5n, self)
        self.pcohe = PCOHE(self._h5n, self)
        self.pcomp = PCOMP(self._h5n, self)
        self.pcompf = PCOMPF(self._h5n, self)
        self.pcompg = PCOMPG(self._h5n, self)
        self.pcompls = PCOMPLS(self._h5n, self)
        self.pconeax = PCONEAX(self._h5n, self)
        self.pconv = PCONV(self._h5n, self)
        self.pconv1 = PCONV1(self._h5n, self)
        self.pconvm = PCONVM(self._h5n, self)
        self.pdamp = PDAMP(self._h5n, self)
        self.pdamp5 = PDAMP5(self._h5n, self)
        self.pdampt = PDAMPT(self._h5n, self)
        self.pelas = PELAS(self._h5n, self)
        self.pelast = PELAST(self._h5n, self)
        self.pfast = PFAST(self._h5n, self)
        self.pgap = PGAP(self._h5n, self)
        self.phbdy = PHBDY(self._h5n, self)
        self.plcomp = PLCOMP(self._h5n, self)
        self.plplane = PLPLANE(self._h5n, self)
        self.plsolid = PLSOLID(self._h5n, self)
        self.pmass = PMASS(self._h5n, self)
        self.prod = PROD(self._h5n, self)
        self.prodn1 = PRODN1(self._h5n, self)
        self.pseam = PSEAM(self._h5n, self)
        self.pshear = PSHEAR(self._h5n, self)
        self.pshearn = PSHEARN(self._h5n, self)
        self.pshell = PSHELL(self._h5n, self)
        self.pshln1 = PSHLN1(self._h5n, self)
        self.pshln2 = PSHLN2(self._h5n, self)
        self.psldn1 = PSLDN1(self._h5n, self)
        self.psolid = PSOLID(self._h5n, self)
        self.ptube = PTUBE(self._h5n, self)
        self.pvisc = PVISC(self._h5n, self)
        self.pweld = PWELD(self._h5n, self)
        self.snorm = SNORM(self._h5n, self)
        # self.vcct = VCCT(self._h5n, self)
        # self.viewex = VIEWEX(self._h5n, self)

    def path(self):
        return self._input.path() + ['PROPERTY']

########################################################################################################################


class MFLUID(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/MFLUID')

########################################################################################################################


class NSM(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/NSM/IDENTITY')

########################################################################################################################


class NSM1(InputTable):
    """
    <group name="NSM1">
        <dataset name="IDENTITY">
            <field name="SID" type="integer"/>
            <field name="TYPE" type="character" size="8"/>
            <field name="PROP" type="character" size="8" description="Name of nonstructural mass entry: NSM1 or NSML1"/>
            <field name="VALUE" type="double"/>
            <field name="ALL" type="integer"/>
            <field name="LIST_POS" type="integer"/>
            <field name="LIST_LEN" type="integer"/>
            <field name="THRU_POS" type="integer"/>
            <field name="THRU_LEN" type="integer"/>
            <field name="THRUBY_POS" type="integer"/>
            <field name="THRUBY_LEN" type="integer"/>
            <field name="DOMAIN_ID" type="integer"/>
        </dataset>
        <dataset name="IDLIST">
            <field name="ID" type="integer"/>
        </dataset>
        <dataset name="THRU">
            <field name="ID1" type="integer"/>
            <field name="ID2" type="integer"/>
        </dataset>
        <dataset name="THRU_BY">
            <field name="ID1" type="integer"/>
            <field name="ID2" type="integer"/>
            <field name="N" type="integer"/>
        </dataset>
    </group>
    """
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/NSM1/IDENTITY',
                                rename={'IDLIST_POS': 'LIST_POS', 'IDLIST_LEN': 'LIST_LEN', 'THRU_BY_POS': 'THRUBY_POS',
                                        'THRU_BY_LEN': 'THRUBY_LEN'})

########################################################################################################################


class NSMADD(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/NSMADD/IDENTITY')

########################################################################################################################


class NSML(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/NSML/IDENTITY')

########################################################################################################################


class NSML1(InputTable):
    """
    <group name="NSML1">
        <dataset name="IDENTITY">
            <field name="SID" type="integer"/>
            <field name="TYPE" type="character" size="8"/>
            <field name="VALUE" type="double"/>
            <field name="ALL" type="integer"/>
            <field name="LIST_POS" type="integer"/>
            <field name="LIST_LEN" type="integer"/>
            <field name="THRU_POS" type="integer"/>
            <field name="THRU_LEN" type="integer"/>
            <field name="THRUBY_POS" type="integer"/>
            <field name="THRUBY_LEN" type="integer"/>
            <field name="DOMAIN_ID" type="integer"/>
        </dataset>
        <dataset name="IDLIST">
            <field name="ID" type="integer"/>
        </dataset>
        <dataset name="THRU">
            <field name="ID1" type="integer"/>
            <field name="ID2" type="integer"/>
        </dataset>
        <dataset name="THRU_BY">
            <field name="ID1" type="integer"/>
            <field name="ID2" type="integer"/>
            <field name="N" type="integer"/>
        </dataset>
    </group>
    """
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/NSML1/IDENTITY',
                                rename={'IDLIST_POS': 'LIST_POS', 'IDLIST_LEN': 'LIST_LEN', 'THRU_BY_POS': 'THRUBY_POS',
                                        'THRU_BY_LEN': 'THRUBY_LEN'}
                                )

########################################################################################################################


class PAABSF(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PAABSF')

########################################################################################################################


class PACABS(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PACABS')

########################################################################################################################


class PACBAR(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PACBAR')

########################################################################################################################


class PACINF(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PACINF')

########################################################################################################################


class PAERO1(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PAERO1')

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        pid = data['PID']
        b1 = data['B1']
        b2 = data['B2']
        b3 = data['B3']
        b4 = data['B4']
        b5 = data['B5']
        b6 = data['B6']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            pid[i] = card.pid

            bi = list(card.Bi)

            diff_len = 6 - len(bi)

            if diff_len > 0:
                bi += [None] * diff_len

            bi = [_ if _ is not None else Defaults.default_int for _ in bi]

            b1[i], b2[i], b3[i], b4[i], b5[i], b6[i] = bi

        return {'IDENTITY': data}


########################################################################################################################


class PAERO2(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PAERO2')

########################################################################################################################


class PAERO3(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PAERO3')

########################################################################################################################


class PAERO4(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PAERO4/IDENTITY')

########################################################################################################################


class PAERO5(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PAERO5/IDENTITY')

########################################################################################################################


class PAXISYM(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PAXISYM')

########################################################################################################################


class PAXSYMH(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PAXSYMH')

########################################################################################################################


class PBAR(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PBAR')

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        pid = data['PID']
        mid = data['MID']
        a = data['A']
        i1 = data['I1']
        i2 = data['I2']
        j = data['J']
        nsm = data['NSM']
        data['FE'] = Defaults.default_double  # blank field
        c1 = data['C1']
        c2 = data['C2']
        d1 = data['D1']
        d2 = data['D2']
        e1 = data['E1']
        e2 = data['E2']
        f1 = data['F1']
        f2 = data['F2']
        k1 = data['K1']
        k2 = data['K2']
        i12 = data['I12']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            pid[i] = card.pid
            mid[i] = card.mid
            a[i] = card.A
            i1[i] = card.i1
            i2[i] = card.i2
            j[i] = card.j
            nsm[i] = card.nsm
            c1[i] = card.c1
            c2[i] = card.c2
            d1[i] = card.d1
            d2[i] = card.d2
            e1[i] = card.e1
            e2[i] = card.e2
            f1[i] = card.f1
            f2[i] = card.f2
            k1[i] = card.k1
            k2[i] = card.k2
            i12[i] = card.i12

        result = {'IDENTITY': data}

        return result


########################################################################################################################

# PBARL msc spec is missing NSM for some reason

class PBARL_INFO_SPEC(object):
    name = 'INFO'
    path = '/NASTRAN/INPUT/PROPERTY/PBARL'
    dtype = [('VALUE', '<f8', (),)]
    is_subtable = True
    same_as = None
    subtables = []


class PBARL_SPEC(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PROPERTY/PBARL'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('GROUP', 'S8', ()), ('TYPE', 'S8', ()), ('NSM', '<f8', ()),
            ('INFO_POS', '<i8', ()), ('INFO_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = [PBARL_INFO_SPEC]


class PBARL(InputTable):
    table_def = TableDef.create(PBARL_SPEC)
    
    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        
        info = {'IDENTITY': {'VALUE': []}}
        
        result = {'IDENTITY': {'PID': [], 'MID': [], 'GROUP': [], 'TYPE': [], 'NSM': [], 'INFO_POS': [],
                               'INFO_LEN': [], 'DOMAIN_ID': []},
                  'INFO': info,
                  '_subtables': ['INFO']}
        
        identity = result['IDENTITY']
        value = info['IDENTITY']['VALUE']
        pid = identity['PID']
        mid = identity['MID']
        group = identity['GROUP']
        type_ = identity['TYPE']
        nsm = identity['NSM']
        info_pos = identity['INFO_POS']
        info_len = identity['INFO_LEN']

        _pos = 0
        for card_id in card_ids:
            card = cards[card_id]
            
            pid.append(card.pid)
            mid.append(card.mid)
            group.append(card.group)
            type_.append(card.beam_type)
            nsm.append(card.nsm)
            info_pos.append(_pos)
            _info_len = len(card.dim)
            info_len.append(_info_len)
            _pos += _info_len
            value += list(card.dim)

        return result
        

########################################################################################################################


class PBARN1(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PBARN1')

########################################################################################################################


class PBCOMP(InputTable):
    """
    <group name="PBCOMP">
        <dataset name="IDENTITY">
            <field name="PID" type="integer"/>
            <field name="MID" type="integer"/>
            <field name="A" type="double"/>
            <field name="I1" type="double"/>
            <field name="I2" type="double"/>
            <field name="I12" type="double"/>
            <field name="J" type="double"/>
            <field name="NSM" type="double"/>
            <field name="K1" type="double"/>
            <field name="K2" type="double"/>
            <field name="M1" type="double"/>
            <field name="M2" type="double"/>
            <field name="N1" type="double"/>
            <field name="N2" type="double"/>
            <field name="NSECT" type="integer"/>
            <field name="POS" type="integer"/>
            <field name="LEN" type="integer"/>
            <field name="DOMAIN_ID" type="integer"/>
        </dataset>
        <dataset name="SECTION">
            <field name="Y" type="double"/>
            <field name="Z" type="double"/>
            <field name="C" type="double"/>
            <field name="MID" type="integer"/>
        </dataset>
    </group>
    """
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PBCOMP/IDENTITY',
                                rename={'SECTION_POS': 'POS', 'SECTION_LEN': 'LEN'})

########################################################################################################################


def _resize(arr, size):

    arr = list(arr)

    first = arr[0]
    last = arr[-1]

    del arr[0]
    try:
        del arr[-1]
    except IndexError:
        pass

    size -= 2

    arr_len = len(arr)
    diff_len = size - arr_len

    if diff_len == 0:
        return arr
    elif diff_len < 0:
        raise Exception

    mid_arr = arr + [None] * diff_len

    return [first] + mid_arr + [last]


class PBEAM(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PBEAM')

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        pid = data['PID']
        mid = data['MID']
        nsegs = data['NSEGS']
        data['CCF'][:] = Defaults.unknown_int
        data['CWELD'][:] = Defaults.unknown_int

        ######################
        so = data['SO']
        xxb = data['XXB']
        a = data['A']
        i1 = data['I1']
        i2 = data['I2']
        i12 = data['I12']
        j = data['J']
        nsm = data['NSM']
        c1 = data['C1']
        c2 = data['C2']
        d1 = data['D1']
        d2 = data['D2']
        e1 = data['E1']
        e2 = data['E2']
        f1 = data['F1']
        f2 = data['F2']
        ######################

        k1 = data['K1']
        k2 = data['K2']
        s1 = data['S1']
        s2 = data['S2']
        nsia = data['NSIA']
        nsib = data['NSIB']
        cwa = data['CWA']
        cwb = data['CWB']
        m1a = data['M1A']
        m2a = data['M2A']
        m1b = data['M1B']
        m2b = data['M2B']
        n1a = data['N1A']
        n2a = data['N2A']
        n1b = data['N1B']
        n2b = data['N2B']

        # TODO: PBEAM - verify so is correct

        _so = {
            '': Defaults.default_double,
            None: Defaults.default_double,
            'NO': 0.,
            'YES': 1.,
            'YESA': 2.
        }

        # TODO: The first and last stations (xxb = 0.0, 1.0 are in slots 0 and 10).
        #       Intermediate slots are 0.0 if they are not defined and nonzero
        #       if the data is meaningful.  The data coming from the PBEAM/PBEAML
        #       is sorted, but only uses as many fields as is necessary to fully
        #       define the card.
        #
        # TODO: PBEAM: verify that above comment has been implemented correctly regarding resizing of data

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            pid[i] = card.pid
            mid[i] = card.mid
            nsegs[i] = len(card.so)
            so[i] = _resize([_so[_] for _ in card.so], 11)
            xxb[i] = _resize(card.xxb, 11)
            a[i] = _resize(card.A, 11)
            i1[i] = _resize(card.i1, 11)
            i2[i] = _resize(card.i2, 11)
            i12[i] = _resize(card.i12, 11)
            j[i] = _resize(card.j, 11)
            nsm[i] = _resize(card.nsm, 11)
            c1[i] = _resize(card.c1, 11)
            c2[i] = _resize(card.c2, 11)
            d1[i] = _resize(card.d1, 11)
            d2[i] = _resize(card.d2, 11)
            e1[i] = _resize(card.e1, 11)
            e2[i] = _resize(card.e2, 11)
            f1[i] = _resize(card.f1, 11)
            f2[i] = _resize(card.f2, 11)
            k1[i] = card.k1
            k2[i] = card.k2
            s1[i] = card.s1
            s2[i] = card.s2
            nsia[i] = card.nsia
            nsib[i] = card.nsib
            cwa[i] = card.cwa
            cwb[i] = card.cwb
            m1a[i] = card.m1a
            m2a[i] = card.m2a
            m1b[i] = card.m1b
            m2b[i] = card.m2b
            n1a[i] = card.n1a
            n2a[i] = card.n2a
            n1b[i] = card.n1b
            n2b[i] = card.n2b

        result = {'IDENTITY': data}

        return result


########################################################################################################################


class PBEAM3(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PBEAM3')

########################################################################################################################


class PBEAML(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PBEAML/IDENTITY',
                                subtables=[
                                    TableDef.create('/NASTRAN/INPUT/PROPERTY/PBEAML/SECTION',
                                                    subtables=[
                                                        TableDef.create('/NASTRAN/INPUT/PROPERTY/PBEAML/DIMS')
                                                    ]
                                                    )
                                ]
                                )

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        dims = {'IDENTITY': {'DIM': []}}
        section = {'IDENTITY': {'SO': [], 'RDIST': [], 'DIMS_POS': [], 'DIMS_LEN': [], 'NSM': []},
                   'DIMS': dims,
                   '_subtables': ['DIMS']}
        result = {'IDENTITY': {'PID': [], 'MID': [], 'GROUP': [], 'TYPE': [],
                               'SECTION_POS': [], 'SECTION_LEN': [], 'DOMAIN_ID': []},
                  'SECTION': section,
                  '_subtables': ['SECTION']}

        section = section['IDENTITY']
        identity = result['IDENTITY']
        dim = dims['IDENTITY']['DIM']
        so = section['SO']
        rdist = section['RDIST']
        dims_pos = section['DIMS_POS']
        dims_len = section['DIMS_LEN']
        nsm = section['NSM']
        pid = identity['PID']
        mid = identity['MID']
        group = identity['GROUP']
        type_ = identity['TYPE']
        section_pos = identity['SECTION_POS']
        section_len = identity['SECTION_LEN']

        # TODO: PBEAML - verify so is correct

        _so = {
            '': Defaults.default_double,
            None: Defaults.default_double,
            'NO': 0.,
            'YES': 1.,
            'YESA': 2.
        }

        _section_pos = 0
        _dims_pos = 0

        for card_id in card_ids:
            card = cards[card_id]

            pid.append(card.pid)
            mid.append(card.mid)
            group.append(card.group)
            type_.append(card.beam_type)
            _section_len = len(card.so)
            section_pos.append(_section_pos)
            section_len.append(_section_len)
            so += [_so[_s] for _s in card.so]
            rdist += list(card.xxb)
            nsm += list(card.nsm)

            for _dim in card.dim:
                _dim_len = len(_dim)
                dims_pos.append(_dims_pos)
                _dims_pos += _dim_len
                dims_len.append(_dim_len)
                dim += list(_dim)

        return result


########################################################################################################################


class PBEMN1(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PBEMN1')

########################################################################################################################


class PBEND(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PBEND')

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        pid = data['PID']
        mid = data['MID']
        a = data['A']
        i1 = data['I1']
        i2 = data['I2']
        j = data['J']
        fsi = data['FSI']
        rm = data['RM']
        t = data['T']
        p = data['P']
        rb = data['RB']
        thetab = data['THETAB']
        c1 = data['C1']
        c2 = data['C2']
        d1 = data['D1']
        d2 = data['D2']
        e1 = data['E1']
        e2 = data['E2']
        f1 = data['F1']
        f2 = data['F2']
        k1 = data['K1']
        k2 = data['K2']
        nsm = data['NSM']
        rc = data['RC']
        zc = data['ZC']
        deltan = data['DELTAN']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            pid[i] = card.pid
            mid[i] = card.mid
            a[i] = card.A
            i1[i] = card.i1
            i2[i] = card.i2
            j[i] = card.j
            fsi[i] = card.fsi
            rm[i] = card.rm
            t[i] = card.t
            p[i] = card.p
            rb[i] = card.rb
            thetab[i] = card.theta_b
            c1[i] = card.c1
            c2[i] = card.c2
            d1[i] = card.d1
            d2[i] = card.d2
            e1[i] = card.e1
            e2[i] = card.e2
            f1[i] = card.f1
            f2[i] = card.f2
            k1[i] = card.k1
            k2[i] = card.k2
            nsm[i] = card.nsm
            rc[i] = card.rc
            zc[i] = card.zc
            deltan[i] = card.delta_n

        result = {'IDENTITY': data}

        return result

########################################################################################################################


# TODO: PBMSECT is complex
# class PBMSECT(CardTable):
#     table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PBMSECT/IDENTITY',
#                                 subtables=[
#                                     TableDef.create('/NASTRAN/INPUT/PROPERTY/PBMSECT/SECTION',
#                                                     subtables=[
#                                                         TableDef.create('/NASTRAN/INPUT/PROPERTY/PBMSECT/BRP')
#                                                     ])
#                                 ])

########################################################################################################################

# TODO: PBRSECT is complex
# class PBRSECT(CardTable):
#     table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PBRSECT/IDENTITY')

########################################################################################################################


class PBUSH(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PBUSH')

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        pid = data['PID']
        k = data['K']
        b = data['B']
        ge = data['GE']
        sa = data['SA']
        st = data['ST']
        ea = data['EA']
        et = data['ET']
        m = data['M']

        def _get_value(obj, attr, default):
            try:
                return getattr(obj, attr)
            except AttributeError:
                return default

        def _get_list(obj, attr, default):
            lst = list(getattr(obj, attr))
            if len(lst) == 0:
                return [default] * 6
            return lst

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            pid[i] = card.pid
            k[i] = _get_list(card, 'Ki', Defaults.default_double)
            b[i] = _get_list(card, 'Bi', Defaults.default_double)
            ge[i] = _get_list(card, 'GEi', Defaults.default_double)
            sa[i] = _get_value(card, 'sa', Defaults.default_double)
            st[i] = _get_value(card, 'st', Defaults.default_double)
            ea[i] = _get_value(card, 'ea', Defaults.default_double)
            et[i] = _get_value(card, 'et', Defaults.default_double)
            m[i] = _get_value(card, 'm', Defaults.default_double)

        result = {'IDENTITY': data}

        return result


########################################################################################################################


# TODO: PBUSH1D verify correctness
class PBUSH1D(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PBUSH1D')

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        pid = data['PID']
        k = data['K']
        c = data['C']
        m = data['M']
        alpha = data['ALPHA']
        sa = data['SA']
        ea = data['EA']
        typea = data['TYPEA']
        cvt = data['CVT']
        cvc = data['CVC']
        expvt = data['EXPVT']
        expvc = data['EXPVC']
        idtsu = data['IDTSU']
        idtcu = data['IDTCU']
        idtsud = data['IDTSUD']
        idcsud = data['IDCSUD']
        types = data['TYPES']
        idts = data['IDTS']
        idcs = data['IDCS']
        idtdu1 = data['IDTDU1']
        idcdu1 = data['IDCDU1']
        typed = data['TYPED']
        idtd1 = data['IDTD1']
        idtd2 = data['IDTD2']
        idtdv1 = data['IDTDV1']
        idcdv1 = data['IDCDV1']
        typeg = data['TYPEG']
        idtg = data['IDTG']
        idcg = data['IDCG']
        idtdu2 = data['IDTDU2']
        idcdu2 = data['IDCDU2']
        idtdv2 = data['IDTDV2']
        idcdv2 = data['IDCDV2']
        typef = data['TYPEF']
        idtf = data['IDTF']
        idcf = data['IDCF']
        ut = data['UT']
        uc = data['UC']

        default_double = Defaults.default_double
        default_int = Defaults.default_int

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            pid[i] = card.pid
            k[i] = card.k
            c[i] = card.c
            m[i] = card.m
            alpha[i] = default_double
            sa[i] = card.sa
            ea[i] = card.se

            shock_type = card.__dict__.get('shock_type', None)

            if shock_type is None:
                typea[i] = 0
                cvt[i] = default_double
                cvc[i] = default_double
                expvt[i] = default_double
                expvc[i] = default_double
                idtsu[i] = default_int
                idtcu[i] = default_int
                idtsud[i] = default_int
                idcsud[i] = default_int
            else:
                if shock_type == 'TABLE':
                    shock_type = 1
                elif shock_type == 'EQUAT':
                    shock_type = 2

                assert shock_type in (1, 2)

                typea[i] = shock_type
                cvt[i] = card.shock_cvt
                cvc[i] = card.shock_cvc
                expvt[i] = card.shock_exp_vt
                expvc[i] = card.shock_exp_vc

                if shock_type == 1:
                    idtsu[i] = card.shock_idts
                    idtcu[i] = default_int
                    idtsud[i] = default_int
                    idcsud[i] = default_int
                else:
                    itdsu[i] = card.idets
                    idtcu[i] = card.idecs
                    idtsud[i] = card.idetsd
                    idcsud[i] = card.idecsd

            spring_type = card.__dict__.get('spring_type', None)

            if spring_type is None:
                types[i] = 0
                idts[i] = default_int
                idcs[i] = default_int
                idtdu1[i] = default_int
                idcdu1[i] = default_int
            else:
                if spring_type == 'TABLE':
                    spring_type = 1
                elif spring_type == 'EQUAT':
                    spring_type = 2

                assert spring_type in (1, 2)

                types[i] = spring_type
                idts[i] = card.spring_idt
                idcs[i] = card.spring_idc
                idtdu1[i] = card.spring_idtdu
                idcdu1[i] = card.spring_idcdu

            damper_type = card.__dict__.get('damper_type', None)

            if damper_type is None:
                typed[i] = default_int
                idtd1[i] = default_int
                idtd2[i] = default_int
                idtdv1[i] = default_int
                idcdv1[i] = default_int
            else:
                if damper_type == 'TABLE':
                    damper_type = 1
                elif damper_type == 'EQUAT':
                    damper_type = 2

                assert damper_type in (1, 2)

                typed[i] = damper_type
                idtd1[i] = card.damper_idt
                idtd2[i] = card.damper_idc
                idtdv1[i] = card.damper_idtdv
                idcdv1[i] = card.damper_idcdv

            gener_idt = card.__dict__.get('gener_idt', None)

            if gener_idt is None:
                typeg[i] = 0
                idtg[i] = default_int
                idcg[i] = default_int
                idtdu2[i] = default_int
                idcdu2[i] = default_int
                idtdv2[i] = default_int
                idcdv2[i] = default_int
            else:
                typeg[i] = 2
                idtg[i] = card.gener_idt
                idcg[i] = card.gener_idc
                idtdu2[i] = card.gener_idtdu
                idcdu2[i] = card.gener_idcdu
                idtdv2[i] = card.gener_idtdv
                idcdv2[i] = card.gener_idcdv

            typef[i] = Defaults.unknown_int
            idtf[i] = Defaults.unknown_int
            idcf[i] = Defaults.unknown_int
            ut[i] = Defaults.unknown_double
            uc[i] = Defaults.unknown_double

        return {'IDENTITY': data}


########################################################################################################################

# TODO: PBUSH2D is complex
# class PBUSH2D(CardTable):
#     table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PBUSH2D/IDENTITY')

########################################################################################################################


class PBUSHT(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PBUSHT')

########################################################################################################################


class PCOHE(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PCOHE')

########################################################################################################################


class PCOMP(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PCOMP/IDENTITY')

    def from_bdf(self, cards):
        _ft = {
            None: Defaults.default_int,
            '': Defaults.default_int,
            'HILL': 1,
            'HOFF': 2,
            'TSAI': 3,
            'STRN': 4
        }
        
        # TODO: check that sout is correct
        _convert_sout = {'YES': 1, 'NO': 0}

        ply = {
            'IDENTITY': {'MID': [], 'T': [], 'THETA': [], 'SOUT': []}
        }

        data = {
            'IDENTITY': {'PID': [],
                         'NPLIES': [],
                         'Z0': [],
                         'NSM': [],
                         'SB': [],
                         'FT': [],
                         'TREF': [],
                         'GE': [],
                         'PLY_POS': [],
                         'PLY_LEN': []},
            'PLY': ply,
            '_subtables': ['PLY']
        }

        identity = data['IDENTITY']
        pid = identity['PID']
        nplies = identity['NPLIES']
        z0 = identity['Z0']
        nsm = identity['NSM']
        sb = identity['SB']
        ft = identity['FT']
        tref = identity['TREF']
        ge = identity['GE']
        ply_pos = identity['PLY_POS']
        ply_len = identity['PLY_LEN']

        ply = ply['IDENTITY']
        mid = ply['MID']
        t = ply['T']
        theta = ply['THETA']
        sout = ply['SOUT']

        card_ids = sorted(iterkeys(cards))

        _ply_pos = 0
        for card_id in card_ids:
            card = cards[card_id]

            _plies = len(card.material_ids)

            pid.append(card.pid)
            nplies.append(_plies)
            z0.append(round(card.z0, 15))
            nsm.append(card.nsm)
            sb.append(card.sb)
            ft.append(_ft[card.ft])
            tref.append(card.tref)
            ge.append(card.ge)
            ply_pos.append(_ply_pos)
            ply_len.append(_plies)
            _ply_pos += _plies

            mid.extend(list(card.material_ids))
            t.extend(list(card.thicknesses))
            theta.extend(card.thetas)
            sout.extend([_convert_sout.get(_, 0) for _ in card.souts])
            
        return data


########################################################################################################################


class PCOMPF(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PCOMPF/IDENTITY',
                                rename={'IDLIST_POS': 'LIST_POS', 'IDLIST_LEN': 'LIST_LEN'})

########################################################################################################################


class PCOMPG(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PCOMPG/IDENTITY')
    
    def from_bdf(self, cards):
        _ft = {
            None: Defaults.default_int,
            '': Defaults.default_int,
            'HILL': 1,
            'HOFF': 2,
            'TSAI': 3,
            'STRN': 4
        }

        # TODO: check that sout is correct
        _convert_sout = {'YES': 1, 'NO': 0}

        ply = {
            'IDENTITY': {'GPLYID': [], 'MID': [], 'THICK': [], 'THETA': [], 'SOUT': [], 'MIDMTX': [],
                         'VF': [], 'VV': [], 'CTEMP': [], 'MOIST': [], 'CRIT': [], 'NFTI': [], 'FTI': []}
        }

        result = {
            'IDENTITY': {'PID': [],
                         'NPLIES': [],
                         'Z0': [],
                         'NSM': [],
                         'SB': [],
                         'FT': [],
                         'TREF': [],
                         'GE': [],
                         'MICRO': [],
                         'PLY_POS': [],
                         'PLY_LEN': [],
                         'DOMAIN_ID': []
                         },
            'PLY': ply,
            '_subtables': ['PLY']
        }

        identity = result['IDENTITY']
        pid = identity['PID']
        nplies = identity['NPLIES']
        z0 = identity['Z0']
        nsm = identity['NSM']
        sb = identity['SB']
        ft = identity['FT']
        tref = identity['TREF']
        ge = identity['GE']
        micro = identity['MICRO']
        ply_pos = identity['PLY_POS']
        ply_len = identity['PLY_LEN']

        ply = ply['IDENTITY']
        gplyid = ply['GPLYID']
        mid = ply['MID']
        thick = ply['THICK']
        theta = ply['THETA']
        sout = ply['SOUT']
        midmtx = ply['MIDMTX']
        vf = ply['VF']
        vv = ply['VV']
        ctemp = ply['CTEMP']
        moist = ply['MOIST']
        crit = ply['CRIT']
        nfti = ply['NFTI']
        fti = ply['FTI']

        card_ids = sorted(iterkeys(cards))
        
        _pos = 0
        for card_id in card_ids:
            card = cards[card_id]
            
            pid.append(card.pid)
            n = len(card.thicknesses)
            nplies.append(n)
            z0.append(card.z0)
            nsm.append(card.nsm)
            sb.append(card.sb)
            ft.append(_ft[card.ft])
            tref.append(card.tref)
            ge.append(card.ge)
            micro.append(Defaults.unknown_str)
            ply_pos.append(_pos)
            ply_len.append(n)
            _pos += n

            gplyid += list(card.global_ply_ids)
            mid += list(card.mids)
            thick += list(card.thicknesses)
            theta += list(card.thetas)
            sout += [_convert_sout[_] for _ in card.souts]
            midmtx += [Defaults.unknown_int] * n
            vf += [Defaults.unknown_double] * n
            vv += [Defaults.unknown_double] * n
            ctemp += [Defaults.unknown_double] * n
            moist += [Defaults.unknown_double] * n
            crit += [Defaults.unknown_str] * n
            nfti += [Defaults.unknown_int] * n
            fti += [Defaults.unknown_str] * n
            
        return result
            

########################################################################################################################


class PCOMPLS(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PCOMPLS/IDENTITY')

########################################################################################################################


class PCONEAX(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PCONEAX')

########################################################################################################################


class PCONV(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PCONV')

########################################################################################################################


class PCONV1(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PCONV1')

########################################################################################################################


class PCONVM(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PCONVM')

########################################################################################################################


class PDAMP(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PDAMP')

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        pid = data['PID']
        b = data['B']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]
            pid[i] = card.pid
            b[i] = card.b

        return {'IDENTITY': data}

########################################################################################################################


class PDAMP5(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PDAMP5')

########################################################################################################################


class PDAMPT(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PDAMPT')

########################################################################################################################


class PELAS(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PELAS')

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)
    
        pid = data['PID']
        k = data['K']
        ge = data['GE']
        s = data['S']
    
        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]
    
            pid[i] = card.pid
            k[i] = card.k
            ge[i] = card.ge
            s[i] = card.s
    
        result = {
            'IDENTITY': data
        }
    
        return result

########################################################################################################################


class PELAST(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PELAST')

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())
        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        pid = data['PID']
        tkid = data['TKID']
        tgeid = data['TGEID']
        tknid = data['TKNID']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            pid[i] = card.pid
            tkid[i] = card.tkid
            tgeid[i] = card.tgeid
            tknid[i] = card.tknid

        result = {
            'IDENTITY': data
        }

        return result

########################################################################################################################


class PFAST(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PFAST')

########################################################################################################################


class PGAP(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PGAP')

########################################################################################################################


class PHBDY(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PHBDY')

########################################################################################################################


class PLCOMP(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PLCOMP/IDENTITY')

########################################################################################################################


class PLPLANE(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PLPLANE')

########################################################################################################################


class PLSOLID(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PLSOLID')

########################################################################################################################


class PMASS(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PMASS')

########################################################################################################################


class PROD(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PROD')

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        pid = data['PID']
        mid = data['MID']
        a = data['A']
        j = data['J']
        c = data['C']
        nsm = data['NSM']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            pid[i] = card.pid
            mid[i] = card.mid
            a[i] = card.A
            j[i] = card.j
            c[i] = card.c
            nsm[i] = card.nsm

        result = {'IDENTITY': data}

        return result


########################################################################################################################


class PRODN1(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PRODN1')

########################################################################################################################


class PSEAM(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PSEAM')

########################################################################################################################


class PSHEAR(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PSHEAR')

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        pid = data['PID']
        mid = data['MID']
        t = data['T']
        nsm = data['NSM']
        f1 = data['F1']
        f2 = data['F2']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            pid[i] = card.pid
            mid[i] = card.mid
            t[i] = card.t
            nsm[i] = card.nsm
            f1[i] = card.f1
            f2[i] = card.f2

        result = {'IDENTITY': data}

        return result

########################################################################################################################


class PSHEARN(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PSHEARN')

########################################################################################################################


class PSHELL(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PSHELL')

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        pid = data['PID']
        mid1 = data['MID1']
        t = data['T']
        mid2 = data['MID2']
        bk = data['BK']
        mid3 = data['MID3']
        ts = data['TS']
        nsm = data['NSM']
        z1 = data['Z1']
        z2 = data['Z2']
        mid4 = data['MID4']

        def _get_mid(val, default):
            if val is None:
                val = default
            return val

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            pid[i] = card.pid
            mid1[i] = card.mid1
            t[i] = card.t
            mid2[i] = _get_mid(card.mid2, Defaults.default_int)
            bk[i] = card.twelveIt3
            mid3[i] = _get_mid(card.mid3, Defaults.default_int)
            ts[i] = card.tst
            nsm[i] = card.nsm
            z1[i] = card.z1
            z2[i] = card.z2
            mid4[i] = _get_mid(card.mid4, Defaults.default_int)

        result = {'IDENTITY': data}

        return result

########################################################################################################################


class PSHLN1(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PSHLN1')

########################################################################################################################


class PSHLN2(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PSHLN2')

########################################################################################################################


class PSLDN1(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PSLDN1')

########################################################################################################################


class PSOLID(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PSOLID')

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        pid = data['PID']
        mid = data['MID']
        cordm = data['CORDM']
        in_ = data['IN']
        stress = data['STRESS']
        isop = data['ISOP']
        fctn = data['FCTN']
            
        _integ = {
            0: 0, 1: 1, 2: 2, 3: 3, 'BUBBLE': 0, 'GAUSS': 1, 'TWO': 2, 'THREE': 3, 
            '': Defaults.default_int, None: Defaults.default_int
        }
        
        _stress = {
            'GRID': Defaults.default_int, 'GAUSS': 1, '': Defaults.default_int, None: Defaults.default_int,
            1: 1
        }
        
        _isop = {0: 0, 1: 1, 'REDUCED': 0, 'FULL': 1, '': Defaults.default_int, None: Defaults.default_int}

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            pid[i] = card.pid
            mid[i] = card.mid
            cordm[i] = card.cordm
            in_[i] = _integ[card.integ]
            stress[i] = _stress[card.stress]
            isop[i] = _isop[card.isop]
            fctn[i] = card.fctn

        result = {'IDENTITY': data}

        return result

########################################################################################################################


# msc format missing OD2
class PTUBE_SPEC(object):
    name = 'PTUBE'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('OD', '<f8', ()), ('T', '<f8', ()), ('NSM', '<f8', ()),
             ('OD2', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


class PTUBE(InputTable):
    table_def = TableDef.create(PTUBE_SPEC)

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        pid = data['PID']
        mid = data['MID']
        od = data['OD']
        t = data['T']
        nsm = data['NSM']
        od2 = data['OD2']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            pid[i] = card.pid
            mid[i] = card.mid
            od[i] = card.OD1
            t[i] = card.t
            nsm[i] = card.nsm
            od2[i] = card.OD2

        result = {'IDENTITY': data}

        return result

########################################################################################################################


class PVISC(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PVISC')

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        pid = data['PID']
        ce = data['CE']
        cr = data['CR']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            pid[i] = card.pid
            ce[i] = card.ce
            cr[i] = card.cr

        return {'IDENTITY': data}


########################################################################################################################


class PWELD(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PWELD')

########################################################################################################################


class SNORM(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/SNORM')

########################################################################################################################

# TODO: VCCT - where and how is dataset SETID used?
# class VCCT(CardTable):
#     """
#     <group name="VCCT">
#         <dataset name="GRID">
#             <field name="GI" type="integer"/>
#         </dataset>
#         <dataset name="IDENTITY">
#             <field name="ID" type="integer"/>
#             <field name="IDCR" type="integer"/>
#             <field name="ITYPE" type="integer"/>
#             <field name="IGROW" type="integer"/>
#             <field name="INCM" type="integer"/>
#             <field name="METHOD" type="integer"/>
#             <field name="TIME" type="double"/>
#             <field name="IACT" type="integer"/>
#             <field name="CGI" type="double"/>
#             <field name="GC" type="double"/>
#             <field name="GTH" type="double"/>
#             <field name="C" type="double"/>
#             <field name="M" type="double"/>
#             <field name="GMIN" type="double"/>
#             <field name="GC2" type="double"/>
#             <field name="GC3" type="double"/>
#             <field name="TABCGI" type="integer"/>
#             <field name="TABGC" type="integer"/>
#             <field name="TABGTH" type="integer"/>
#             <field name="TABC" type="integer"/>
#             <field name="TABM" type="integer"/>
#             <field name="TABGMIN" type="integer"/>
#             <field name="TABGC2" type="integer"/>
#             <field name="TABGC3" type="integer"/>
#             <field name="GRID_POS" type="integer"/>
#             <field name="GRID_LEN" type="integer"/>
#             <field name="THBY_POS" type="integer"/>
#             <field name="THBY_LEN" type="integer"/>
#             <field name="DOMAIN_ID" type="integer"/>
#         </dataset>
#         <dataset name="SETID">
#             <field name="SET3ID" type="integer"/>
#         </dataset>
#         <dataset name="THBY">
#             <field name="G1" type="integer"/>
#             <field name="G2" type="integer"/>
#             <field name="GINC" type="integer"/>
#         </dataset>
#     </group>
#     """
#     table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/VCCT/IDENTITY')

########################################################################################################################

# TODO: VIEWEX doesn't conform to the subtable pos/len scheme
# class VIEWEX(CardTable):
#     table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/VIEWEX/IDENTITY')

########################################################################################################################
