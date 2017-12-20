from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range

import tables
import numpy as np

from .card_table import CardTable, TableDef, TableData


class Property(object):
    def __init__(self, h5n, input):
        self._h5n = h5n
        self._input = input

        self.pbar = PBAR(self._h5n, self)
        self.pbarl = PBARL(self._h5n, self)
        self.pbeam = PBEAM(self._h5n, self)
        self.pbeaml = PBEAML(self._h5n, self)
        self.pbush = PBUSH(self._h5n, self)
        self.pcomp = PCOMP(self._h5n, self)
        self.prod = PROD(self._h5n, self)
        self.pshear = PSHEAR(self._h5n, self)
        self.pshell = PSHELL(self._h5n, self)

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


class PBAR(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PBAR')

    @staticmethod
    def from_bdf(card):
        fe = np.nan  # TODO: what is FE???
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


def _resize(arr, size):
    arr_len = len(arr)
    diff_len = size - arr_len

    if diff_len == 0:
        return arr
    elif diff_len < 0:
        raise Exception

    return list(arr) + [None] * diff_len


class PBEAM(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PBEAM')

    @staticmethod
    def from_bdf(card):
        # TODO: PBEAM ccf, cweld
        # Constant cross-section flag: 1=yes and 0=no; 2 is also possible, but no idea
        ccf = -10
        cweld = -10

        nsegs = len(card.so)

        # TODO: PBEAM - why is SO a double in the spec and not a string??? supposed to be YES, YESA, or NO.........

        # _so = {
        #     '': 0.,
        #     'NO': -1.,
        #     'YES': 1.,
        #     'YESA': 2.
        # }

        _so = {
            '': np.nan,
            'NO': np.nan,
            'YES': np.nan,
            'YESA': np.nan
        }

        # TODO: The first and last stations (xxb = 0.0, 1.0 are in slots 0 and 10).
        #       Intermediate slots are 0.0 if they are not defined and nonzero
        #       if the data is meaningful.  The data coming from the PBEAM/PBEAML
        #       is sorted, but only uses as many fields as is necessary to fully
        #       define the card.
        #
        # _so = {
        #       'YES' : 0.,
        #       'NO' : 1,
        # }
        #
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
        # TODO: PBEAML - why is SO a double????

        data = [card.pid, card.mid, card.group, card.beam_type]
        data = TableData([data])

        data.subdata_len = np.array([[len(card.so)]])

        def _convert_so(_so):
            return np.nan

        so = [_convert_so(_so) for _so in card.so]

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
            None: 0,
            '': 0,
            'HILL': 1,
            'HOFF': 2,
            'TSAI': 3,
            'STRN': 4
        }

        data = [card.pid, len(card.material_ids), round(card.z0, 15), card.nsm, card.sb, _ft[card.ft], card.tref, card.ge]
        data = TableData([data])

        data.subdata_len = np.array([[len(card.material_ids)]])

        def _convert_sout(_sout):
            return -10

        subdata = TableData()
        subdata.data = [
            [card.material_ids[i], card.thicknesses[i], card.thetas[i], _convert_sout(card.souts[i])]
            for i in range(len(card.material_ids))
        ]

        data.subdata.append(subdata)

        return data

########################################################################################################################


class PROD(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PROD')

    @staticmethod
    def from_bdf(card):
        data = [card.pid, card.mid, card.A, card.j, card.c, card.nsm]
        return TableData([data])

########################################################################################################################


class PSHEAR(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PSHEAR')

    @staticmethod
    def from_bdf(card):
        data = [card.pid, card.mid, card.t, card.nsm, card.f1, card.f2]
        return TableData([data])

########################################################################################################################


class PSHELL(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/PROPERTY/PSHELL')

    @staticmethod
    def from_bdf(card):
        data = [card.pid, card.mid1, card.t, card.mid2, card.twelveIt3, card.mid3, card.tst, card.nsm, card.z1,
                card.z2, card.mid4]
        return TableData([data])

########################################################################################################################
