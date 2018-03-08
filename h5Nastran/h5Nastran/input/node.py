from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range

import tables
import numpy as np
from typing import Dict

from .card_table import CardTable, TableDef, DataHelper
from h5Nastran.utilities import ImmutableDict


class Node(object):
    def __init__(self, h5n, input):
        self._h5n = h5n
        self._input = input

        self.grid = GRID(self._h5n, self)

    def path(self):
        return self._input.path() + ['NODE']

    def read(self):
        for key, item in iteritems(self.__dict__):
            if key.startswith('_'):
                continue
            try:
                item.read()
            except AttributeError:
                pass

########################################################################################################################


class GRID(CardTable):
    table_def = TableDef.create('/NASTRAN/INPUT/NODE/GRID')

    def __init__(self, *args, **kwargs):
        super(GRID, self).__init__(*args, **kwargs)
        self._grid_dict = None  # type: ImmutableDict[int, np.ndarray]
        self._grid_in_basic_dict = None  # type: ImmutableDict[int, np.ndarray]
        self._grid_index = {}  # type: Dict[int, int]

    def read(self):
        super(GRID, self).read()
        self._grid_dict = None
        self._grid_in_basic_dict = None
        self._grid_index.clear()

    @classmethod
    def from_bdf(cls, cards):
        card_ids = sorted(cards.keys())

        data = np.empty(len(card_ids), dtype=cls.table_def.dtype)

        id_ = data['ID']
        cp = data['CP']
        x = data['X']
        cd = data['CD']
        ps = data['PS']
        seid = data['SEID']

        def _get_val(val, default):
            if val in ('', None):
                val = default
            return val

        i = -1
        for card_id in card_ids:
            i += 1
            card = cards[card_id]

            id_[i] = card.nid
            cp[i] = card.cp
            x[i] = card.xyz
            cd[i] = card.cd
            ps[i] = _get_val(card.ps, DataHelper.default_int)
            seid[i] = card.seid

        result = {'IDENTITY': data}

        return result

    def get_grid(self, nid):
        try:
            index = self._grid_index[nid]
        except KeyError:
            index = self._grid_index[nid] = np.searchsorted(self.grid['ID'], nid)

        return self.grid[index]

    def grid_dict(self):
        # type: (bool) -> ImmutableDict[int, np.ndarray]

        if self._grid_dict is not None:
            return self._grid_dict

        result = {}

        data = self.grid

        ids = data['ID']
        x = data['X']

        for i in range(len(ids)):
            nid = ids[i]
            arr = np.array(x[i])
            arr.setflags(write=False)
            result[nid] = arr

        self._grid_dict = ImmutableDict(result)

        return self._grid_dict

    def grid_in_basic_dict(self):
        # type: () -> ImmutableDict[int, np.ndarray]

        if self._grid_in_basic_dict is not None:
            return self._grid_in_basic_dict

        result = {}

        data = self.grid

        ids = data['ID']
        x = data['X']
        cp = data['CP']

        position_to_basic = self._h5n.input.coordinate_system.h5n_transformation.position_to_basic

        for i in range(len(ids)):
            nid = ids[i]
            arr = position_to_basic(x[i], cp[i])
            arr.setflags(write=False)
            result[nid] = arr
            
        self._grid_in_basic_dict = ImmutableDict(result)

        return self._grid_in_basic_dict

########################################################################################################################

