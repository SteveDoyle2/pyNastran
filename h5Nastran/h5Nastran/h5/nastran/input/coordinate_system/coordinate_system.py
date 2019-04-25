from __future__ import print_function, absolute_import

import numpy as np
from six import iteritems
from six.moves import range

from h5Nastran.h5nastrannode import H5NastranNode
from .transformation import Transformation
from ...input.input_table import InputTable, TableDef


class CoordinateSystem(H5NastranNode):
    def __init__(self, h5n, input):
        self._h5n = h5n
        self._input = input

        self.cord1c = CORD1C(self._h5n, self)
        self.cord1r = CORD1R(self._h5n, self)
        self.cord1s = CORD1S(self._h5n, self)
        self.cord2c = CORD2C(self._h5n, self)
        self.cord2r = CORD2R(self._h5n, self)
        self.cord2s = CORD2S(self._h5n, self)
        self.cord3g = CORD3G(self._h5n, self)
        self.cord3r = CORD3R(self._h5n, self)
        # self.transformation = TRANSFORMATION(self._h5n, self)
        self.h5n_transformation = Transformation(self._h5n, self)

    def path(self):
        return self._input.path() + ['COORDINATE_SYSTEM']

    def update(self):
        self.h5n_transformation.update()
        
    def to_bdf(self, bdf):
        for key, item in iteritems(self.__dict__):
            if key.startswith('_') or key == 'h5n_transformation':
                continue
            try:
                item.to_bdf(bdf)
            except NotImplementedError:
                pass


########################################################################################################################


class CORD1C(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/COORDINATE_SYSTEM/CORD1C')

########################################################################################################################


class CORD1R(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/COORDINATE_SYSTEM/CORD1R')

########################################################################################################################

class CORD1S(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/COORDINATE_SYSTEM/CORD1S')

########################################################################################################################


class _CORD2(InputTable):

    add_card = None

    def to_bdf(self, bdf):
        add_card = getattr(bdf, self.add_card)

        data = self.identity

        cid = data['CID']
        rid = data['RID']
        a1 = data['A1']
        a2 = data['A2']
        a3 = data['A3']
        b1 = data['B1']
        b2 = data['B2']
        b3 = data['B3']
        c1 = data['C1']
        c2 = data['C2']
        c3 = data['C3']

        for i in range(len(cid)):
            _cid = cid[i]
            _rid = rid[i]
            origin = [a1[i], a2[i], a3[i]]
            zaxis = [b1[i], b2[i], b3[i]]
            xzplane = [c1[i], c2[i], c3[i]]
            add_card(_cid, _rid, origin, zaxis, xzplane)

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        data = np.empty(len(card_ids), dtype=self.table_def.dtype)

        cid = data['CID']
        rid = data['RID']
        a1 = data['A1']
        a2 = data['A2']
        a3 = data['A3']
        b1 = data['B1']
        b2 = data['B2']
        b3 = data['B3']
        c1 = data['C1']
        c2 = data['C2']
        c3 = data['C3']

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            cid[i] = card.cid
            rid[i] = card.rid
            a1[i], a2[i], a3[i] = card.e1
            b1[i], b2[i], b3[i] = card.e2
            c1[i], c2[i], c3[i] = card.e3

        result = {'IDENTITY': data}

        return result


########################################################################################################################


class CORD2C(_CORD2):
    table_def = TableDef.create('/NASTRAN/INPUT/COORDINATE_SYSTEM/CORD2C')
    add_card = 'add_cord2c'

########################################################################################################################


class CORD2R(_CORD2):
    table_def = TableDef.create('/NASTRAN/INPUT/COORDINATE_SYSTEM/CORD2R')
    add_card = 'add_cord2r'

########################################################################################################################


class CORD2S(_CORD2):
    table_def = TableDef.create('/NASTRAN/INPUT/COORDINATE_SYSTEM/CORD2S')
    add_card = 'add_cord2s'

########################################################################################################################


class CORD3G(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/COORDINATE_SYSTEM/CORD3G')

########################################################################################################################


class CORD3R(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/COORDINATE_SYSTEM/CORD3R')

########################################################################################################################


# class TRANSFORMATION(CardTable):
#     table_def = TableDef.create('/NASTRAN/INPUT/COORDINATE_SYSTEM/TRANSFORMATION/IDENTITY')

########################################################################################################################
