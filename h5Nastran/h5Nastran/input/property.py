from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range

import tables
import numpy as np

from .card_table import CardTable, TableDef, TableData
from ..data_helper import DataHelper


class Property(object):
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

    def read(self):
        for key, item in iteritems(self.__dict__):
            if key.startswith('_'):
                continue
            try:
                item.read()
            except AttributeError:
                pass

########################################################################################################################


class MFLUID(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/MFLUID')

########################################################################################################################


class NSM(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/NSM/IDENTITY')

########################################################################################################################


class NSM1(CardTable):
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


class NSMADD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/NSMADD/IDENTITY')

########################################################################################################################


class NSML(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/NSML/IDENTITY')

########################################################################################################################


class NSML1(CardTable):
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


class PAABSF(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PAABSF')

########################################################################################################################


class PACABS(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PACABS')

########################################################################################################################


class PACBAR(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PACBAR')

########################################################################################################################


class PACINF(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PACINF')

########################################################################################################################


class PAERO1(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PAERO1')

########################################################################################################################


class PAERO2(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PAERO2')

########################################################################################################################


class PAERO3(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PAERO3')

########################################################################################################################


class PAERO4(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PAERO4/IDENTITY')

########################################################################################################################


class PAERO5(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PAERO5/IDENTITY')

########################################################################################################################


class PAXISYM(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PAXISYM')

########################################################################################################################


class PAXSYMH(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PAXSYMH')

########################################################################################################################


class PBAR(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PBAR')

    @staticmethod
    def from_bdf(card):
        fe = DataHelper.unknown_double  # TODO: what is FE???
        data = [card.pid, card.mid, card.A, card.i1, card.i2, card.j, card.nsm, fe, card.c1, card.c2, card.d1, card.d2,
                card.e1, card.e2, card.f1, card.f2, card.k1, card.k2, card.i12]
        return TableData([data])

########################################################################################################################


class PBARL(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PBARL/IDENTITY')

    @staticmethod
    def from_bdf(card):
        data = [card.pid, card.mid, card.group, card.beam_type, card.nsm]
        data = TableData([data])
        data.subdata_len = np.array([[len(card.dim)]])

        dim = np.array(card.dim)
        dim = dim.reshape((dim.shape[0], 1))

        data.subdata.append(TableData(dim))

        return data

########################################################################################################################


class PBARN1(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PBARN1')

########################################################################################################################


class PBCOMP(CardTable):
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
    del arr[-1]

    size -= 2

    arr_len = len(arr)
    diff_len = size - arr_len

    if diff_len == 0:
        return arr
    elif diff_len < 0:
        raise Exception

    mid_arr = arr + [None] * diff_len

    return [first] + mid_arr + [last]


class PBEAM(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PBEAM')

    @staticmethod
    def from_bdf(card):
        # TODO: PBEAM ccf, cweld
        # Constant cross-section flag: 1=yes and 0=no; 2 is also possible, but no idea
        ccf = DataHelper.unknown_int
        cweld = DataHelper.unknown_int

        nsegs = len(card.so)

        # TODO: PBEAM - verify so is correct

        _so = {
            '': DataHelper.default_double,
            None: DataHelper.default_double,
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

        so = [_so[_] for _ in card.so]
        so = _resize(so, 11)
        xxb = _resize(card.xxb, 11)
        A = _resize(card.A, 11)
        i1 = _resize(card.i1, 11)
        i2 = _resize(card.i2, 11)
        i12 = _resize(card.i12, 11)
        j = _resize(card.j, 11)
        nsm = _resize(card.nsm, 11)
        c1 = _resize(card.c1, 11)
        c2 = _resize(card.c2, 11)
        d1 = _resize(card.d1, 11)
        d2 = _resize(card.d2, 11)
        e1 = _resize(card.e1, 11)
        e2 = _resize(card.e2, 11)
        f1 = _resize(card.f1, 11)
        f2 = _resize(card.f2, 11)

        data = [
            card.pid, card.mid, nsegs, ccf, cweld,
            so, xxb, A, i1, i2, i12, j, nsm, c1, c2, d1, d2, e1, e2, f1, f2,
            card.k1, card.k2, card.s1, card.s2, card.nsia, card.nsib,
            card.cwa, card.cwb, card.m1a, card.m2a, card.m1b, card.m2b, card.n1a, card.n2a, card.n1b, card.n2b
        ]
        return TableData([data])

########################################################################################################################


class PBEAM3(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PBEAM3')

########################################################################################################################


class PBEAML(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PBEAML/IDENTITY',
                                subtables=[
                                    TableDef.create('/NASTRAN/INPUT/PROPERTY/PBEAML/SECTION',
                                                    subtables=[
                                                        TableDef.create('/NASTRAN/INPUT/PROPERTY/PBEAML/DIMS')
                                                    ]
                                                    )
                                ]
                                )

    @staticmethod
    def from_bdf(card):
        data = [card.pid, card.mid, card.group, card.beam_type]
        data = TableData([data])

        data.subdata_len = np.array([[len(card.so)]])

        # TODO: PBEAML - verify so is correct

        _so = {
            '': DataHelper.default_double,
            None: DataHelper.default_double,
            'NO': 0.,
            'YES': 1.,
            'YESA': 2.
        }

        so = [_so(_so) for _so in card.so]

        subdata = [[so[i], card.xxb[i], card.nsm[i]] for i in range(len(card.so))]
        subdata = TableData(subdata)
        subdata.subdata_len = np.array([[len(card.dim[i])] for i in range(len(card.dim))])

        data.subdata.append(subdata)

        dims = []
        for i in range(len(card.dim)):
            for j in range(len(card.dim[i])):
                dims.append([card.dim[i][j]])

        _subdata = TableData(dims)

        subdata.subdata.append(_subdata)

        return data

########################################################################################################################


class PBEMN1(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PBEMN1')

########################################################################################################################


class PBEND(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PBEND')

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


class PBUSH(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PBUSH')

    """
        <dataset name="PBUSH">
            <field name="PID" type="integer"/>
            <field name="K" type="double" size="6"/>
            <field name="B" type="double" size="6"/>
            <field name="GE" type="double" size="6"/>
            <field name="SA" type="double"/>
            <field name="ST" type="double"/>
            <field name="EA" type="double"/>
            <field name="ET" type="double"/>
            <field name="M" type="double"/>
            <field name="DOMAIN_ID" type="integer"/>
        </dataset>
    """

    @staticmethod
    def from_bdf(card):
        def _get_value(obj, attr):
            try:
                return getattr(obj, attr)
            except AttributeError:
                return None

        def _get_list(obj, attr):
            lst = list(getattr(obj, attr))
            if len(lst) == 0:
                return [None] * 6
            return lst

        pid = card.pid
        ki = _get_list(card, 'Ki')
        bi = _get_list(card, 'Bi')
        gei = _get_list(card, 'GEi')
        sa = _get_value(card, 'sa')
        st = _get_value(card, 'st')
        ea = _get_value(card, 'ea')
        et = _get_value(card, 'et')
        m = _get_value(card, 'm')

        return TableData([[pid, ki, bi, gei, sa, st, ea, et, m]])

########################################################################################################################


class PBUSH1D(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PBUSH1D')

########################################################################################################################

# TODO: PBUSH2D is complex
# class PBUSH2D(CardTable):
#     table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PBUSH2D/IDENTITY')

########################################################################################################################


class PBUSHT(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PBUSHT')

########################################################################################################################


class PCOHE(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PCOHE')

########################################################################################################################


class PCOMP(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PCOMP/IDENTITY')

    """
    <group name="PCOMP">
        <dataset name="IDENTITY">
            <field name="PID" type="integer" description="Property identification number"/>
            <field name="NPLIES" type="integer" description="Number of plies, negative if symmetric"/>
            <field name="Z0" type="double" description="Distance from the reference plane to the bottom surface"/>
            <field name="NSM" type="double" description="Nonstructural mass per unit area"/>
            <field name="SB" type="double" description="Allowable shear stress of the bonding material"/>
            <field name="FT" type="integer" description="Failure Theory: 0=blank, 1=HILL, 2=HOFF, 3=TSAI, 4=STRN"/>
            <field name="TREF" type="double" description="Reference temperature"/>
            <field name="GE" type="double" description="Damping coefficient"/>
            <field name="PLY_POS" type="integer" description="Start position of PLY data in the PLY dataset"/>
            <field name="PLY_LEN" type="integer" description="Length of PLY data in the PLY dataset"/>
            <field name="DOMAIN_ID" type="integer"/>
        </dataset>
        <dataset name="PLY">
            <field name="MID" type="integer"/>
            <field name="T" type="double"/>
            <field name="THETA" type="double"/>
            <field name="SOUT" type="integer"/>
        </dataset>
    </group>
    """

    @staticmethod
    def from_bdf(card):

        # TODO: PCOMP - why is SOUT an integer????

        _ft = {
            None: DataHelper.default_int,
            '': DataHelper.default_int,
            'HILL': 1,
            'HOFF': 2,
            'TSAI': 3,
            'STRN': 4
        }

        data = [card.pid, len(card.material_ids), round(card.z0, 15), card.nsm, card.sb, _ft[card.ft], card.tref, card.ge]
        data = TableData([data])

        data.subdata_len = np.array([[len(card.material_ids)]])

        def _convert_sout(_sout):
            return DataHelper.unknown_int

        subdata = TableData()
        subdata.data = [
            [card.material_ids[i], card.thicknesses[i], card.thetas[i], _convert_sout(card.souts[i])]
            for i in range(len(card.material_ids))
        ]

        data.subdata.append(subdata)

        return data

########################################################################################################################


class PCOMPF(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PCOMPF/IDENTITY',
                                rename={'IDLIST_POS': 'LIST_POS', 'IDLIST_LEN': 'LIST_LEN'})

########################################################################################################################


class PCOMPG(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PCOMPG/IDENTITY')

########################################################################################################################


class PCOMPLS(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PCOMPLS/IDENTITY')

########################################################################################################################


class PCONEAX(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PCONEAX')

########################################################################################################################


class PCONV(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PCONV')

########################################################################################################################


class PCONV1(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PCONV1')

########################################################################################################################


class PCONVM(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PCONVM')

########################################################################################################################


class PDAMP(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PDAMP')

########################################################################################################################


class PDAMP5(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PDAMP5')

########################################################################################################################


class PDAMPT(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PDAMPT')

########################################################################################################################


class PELAS(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PELAS')

########################################################################################################################


class PELAST(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PELAST')

########################################################################################################################


class PFAST(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PFAST')

########################################################################################################################


class PGAP(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PGAP')

########################################################################################################################


class PHBDY(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PHBDY')

########################################################################################################################


class PLCOMP(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PLCOMP/IDENTITY')

########################################################################################################################


class PLPLANE(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PLPLANE')

########################################################################################################################


class PLSOLID(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PLSOLID')

########################################################################################################################


class PMASS(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PMASS')

########################################################################################################################


class PROD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PROD')

    @staticmethod
    def from_bdf(card):
        data = [card.pid, card.mid, card.A, card.j, card.c, card.nsm]
        return TableData([data])

########################################################################################################################


class PRODN1(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PRODN1')

########################################################################################################################


class PSEAM(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PSEAM')

########################################################################################################################


class PSHEAR(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PSHEAR')

    @staticmethod
    def from_bdf(card):
        data = [card.pid, card.mid, card.t, card.nsm, card.f1, card.f2]
        return TableData([data])

########################################################################################################################


class PSHEARN(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PSHEARN')

########################################################################################################################


class PSHELL(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PSHELL')

    @staticmethod
    def from_bdf(card):
        data = [card.pid, card.mid1, card.t, card.mid2, card.twelveIt3, card.mid3, card.tst, card.nsm, card.z1,
                card.z2, card.mid4]
        return TableData([data])

########################################################################################################################


class PSHLN1(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PSHLN1')

########################################################################################################################


class PSHLN2(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PSHLN2')

########################################################################################################################


class PSLDN1(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PSLDN1')

########################################################################################################################


class PSOLID(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PSOLID')

########################################################################################################################


class PTUBE(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PTUBE')

########################################################################################################################


class PVISC(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PVISC')

########################################################################################################################


class PWELD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PWELD')

########################################################################################################################


class SNORM(CardTable):
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
