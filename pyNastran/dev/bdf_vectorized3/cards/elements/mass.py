from __future__ import annotations
from itertools import zip_longest
from typing import TYPE_CHECKING
import numpy as np
from pyNastran.bdf.field_writer_8 import print_card_8
#from pyNastran.bdf.field_writer_16 import print_card_16, print_scientific_16, print_field_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, integer_or_blank, double_or_blank,
)
from pyNastran.dev.bdf_vectorized3.cards.base_card import Element
from pyNastran.dev.bdf_vectorized3.cards.write_utils import array_str, array_default_int
#from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.utils import cast_int_array
if TYPE_CHECKING:
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.dev.bdf_vectorized3.bdf import BDF


class CONM1(Element):
    """
    Concentrated Mass Element Connection, General Form
    Defines a 6 x 6 symmetric mass matrix at a geometric grid point

    +--------+-----+-----+-----+-----+-----+-----+-----+-----+
    |    1   |  2  |  3  |  4  |  5  |  6  |  7  |  8  |  9  |
    +========+=====+=====+=====+=====+=====+=====+=====+=====+
    |  CONM1 | EID |  G  | CID | M11 | M21 | M22 | M31 | M32 |
    +--------+-----+-----+-----+-----+-----+-----+-----+-----+
    |        | M33 | M41 | M42 | M43 | M44 | M51 | M52 | M53 |
    +--------+-----+-----+-----+-----+-----+-----+-----+-----+
    |        | M54 | M55 | M61 | M62 | M63 | M64 | M65 | M66 |
    +--------+-----+-----+-----+-----+-----+-----+-----+-----+

    """
    def __init__(self, model: BDF):
        super().__init__(model)
        self.coord_id = np.array([], dtype='int32')
        self.node_id = np.array([], dtype='int32')
        self._mass = np.zeros((0, 6, 6), dtype='float64')

    def add_card(self, card: BDFCard, comment: str=''):
        m = np.zeros((6, 6), dtype='float64')
        eid = integer(card, 1, 'eid')
        nid = integer(card, 2, 'nid')
        cid = integer_or_blank(card, 3, 'cid', default=0)

        m[0, 0] = double_or_blank(card, 4, 'M11', default=0.)
        m[1, 0] = double_or_blank(card, 5, 'M21', default=0.)
        m[1, 1] = double_or_blank(card, 6, 'M22', default=0.)
        m[2, 0] = double_or_blank(card, 7, 'M31', default=0.)
        m[2, 1] = double_or_blank(card, 8, 'M32', default=0.)
        m[2, 2] = double_or_blank(card, 9, 'M33', default=0.)
        m[3, 0] = double_or_blank(card, 10, 'M41', default=0.)
        m[3, 1] = double_or_blank(card, 11, 'M42', default=0.)
        m[3, 2] = double_or_blank(card, 12, 'M43', default=0.)
        m[3, 3] = double_or_blank(card, 13, 'M44', default=0.)
        m[4, 0] = double_or_blank(card, 14, 'M51', default=0.)
        m[4, 1] = double_or_blank(card, 15, 'M52', default=0.)
        m[4, 2] = double_or_blank(card, 16, 'M53', default=0.)
        m[4, 3] = double_or_blank(card, 17, 'M54', default=0.)
        m[4, 4] = double_or_blank(card, 18, 'M55', default=0.)
        m[5, 0] = double_or_blank(card, 19, 'M61', default=0.)
        m[5, 1] = double_or_blank(card, 20, 'M62', default=0.)
        m[5, 2] = double_or_blank(card, 21, 'M63', default=0.)
        m[5, 3] = double_or_blank(card, 22, 'M64', default=0.)
        m[5, 4] = double_or_blank(card, 23, 'M65', default=0.)
        m[5, 5] = double_or_blank(card, 24, 'M66', default=0.)
        self.cards.append((eid, nid, cid, m, comment))
        self.n += 1

    def parse_cards(self):
        assert self.n >= 0, self.n
        if len(self.cards) == 0:
            return
        ncards = len(self.cards)
        assert ncards > 0, ncards
        self.element_id = np.zeros(ncards, dtype='int32')
        self._mass = np.zeros((ncards, 6, 6), dtype='float64')
        self.coord_id = np.zeros(ncards, dtype='int32')
        self.node_id = np.zeros(ncards, dtype='int32')
        for icard, card in enumerate(self.cards):
            eid, nid, cid, m, comment = card
            self.element_id[icard] = eid
            self.node_id[icard] = nid
            self.coord_id[icard] = cid
            self._mass[icard, :, :] = m
        self.cards = []

    def write(self, size: int=8) -> str:
        if len(self.element_id) == 0:
            return ''
        lines = []
        if size == 8:
            print_card = print_card_8

        element_ids = array_str(self.element_id, size=size)
        node_ids = array_str(self.node_id, size=size)
        coord_ids = array_default_int(self.coord_id, default=0, size=size)
        for eid, nid, cid, m in zip_longest(element_ids, node_ids, coord_ids, self._mass):
            list_fields = [
                'CONM1', eid, nid, cid, m[0, 0], m[1, 0], m[1, 1],
                m[2, 0], m[2, 1], m[2, 2], m[3, 0], m[3, 1], m[3, 2],
                m[3, 3], m[4, 0], m[4, 1], m[4, 2], m[4, 3], m[4, 4],
                m[5, 0], m[5, 1], m[5, 2], m[5, 3], m[5, 4], m[5, 5],
            ]
            lines.append(print_card(list_fields))
        return ''.join(lines)

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        cid = self.model.coord.coord_id
        geom_check(self,
                   missing,
                   node=(nid, self.node_id),
                   coord=(cid, self.coord_id))

    def mass(self) -> np.ndarray:
        mass = (self._mass[:, 0, 0] + self._mass[:, 1, 1] + self._mass[:, 2, 2]) / 3.
        return mass

    def centroid(self) -> np.ndarray:
        nid = self.model.grid.node_id
        xyz = self.model.grid.xyz_cid0()
        inode = np.searchsorted(nid, self.node_id)
        centroid = xyz[inode, :]
        return centroid
        #assert np.array_equal(nid[inode], self.node_id)
        #xyz = xyz[inode, :] + self.xyz_offset
        #return xyz


class CONM2(Element):
    """
    +-------+--------+-------+-------+---------+------+------+------+
    |   1   |    2   |    3  |   4   |    5    |  6   |  7   |   8  |
    +=======+========+=======+=======+=========+======+======+======+
    | CONM2 |   EID  |  NID  |  CID  |  MASS   |  X1  |  X2  |  X3  |
    +-------+--------+-------+-------+---------+------+------+------+
    |       |   I11  |  I21  |  I22  |   I31   |  I32 |  I33 |      |
    +-------+--------+-------+-------+---------+------+------+------+
    | CONM2 | 501274 | 11064 |       | 132.274 |      |      |      |
    +-------+--------+-------+-------+---------+------+------+------+

    """

    def add_card(self, card: BDFCard, comment: str=''):
        eid = integer(card, 1, 'eid')
        nid = integer(card, 2, 'nid')
        cid = integer_or_blank(card, 3, 'cid', default=0)
        mass = double_or_blank(card, 4, 'mass', default=0.)

        X = [
            double_or_blank(card, 5, 'x1', default=0.0),
            double_or_blank(card, 6, 'x2', default=0.0),
            double_or_blank(card, 7, 'x3', default=0.0),
        ]

        I = [
            double_or_blank(card, 9, 'I11', default=0.0),
            double_or_blank(card, 10, 'I21', default=0.0),
            double_or_blank(card, 11, 'I22', default=0.0),
            double_or_blank(card, 12, 'I31', default=0.0),
            double_or_blank(card, 13, 'I32', default=0.0),
            double_or_blank(card, 14, 'I33', default=0.0),
        ]
        assert len(card) <= 15, f'len(CONM2 card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, nid, cid, mass, X, I, comment))
        self.n += 1

    def __apply_slice__(self, elem: CONM2, i: np.ndarray):
        elem.element_id = self.element_id[i]
        elem._mass = self._mass[i]
        elem.coord_id = self.coord_id[i]
        elem.node_id = self.node_id[i]
        elem.xyz_offset = self.xyz_offset[i, :]
        elem.inertia = self.inertia[i, :]

    def parse_cards(self):
        assert self.n >= 0, self.n
        if len(self.cards) == 0:
            return
        ncards = len(self.cards)
        assert ncards > 0, ncards
        element_id = []
        mass = np.zeros(ncards, dtype='float64')
        coord_id = np.zeros(ncards, dtype='int32')
        node_id = []
        xyz_offset = np.zeros((ncards, 3), dtype='float64')
        #I11, I21, I22, I31, I32, I33 = I
        inertia = np.zeros((ncards, 6), dtype='float64')
        for icard, card in enumerate(self.cards):
            (eid, nid, cid, massi, X, I, comment) = card
            element_id.append(eid)
            node_id.append(nid)
            coord_id[icard] = cid
            mass[icard] = massi
            xyz_offset[icard, :] = X
            inertia[icard, :] = I
        self._save(element_id, mass, coord_id, node_id, xyz_offset, inertia)
        self.sort()
        self.cards = []

    def _save(self, element_id, mass, coord_id, node_id, xyz_offset, inertia):
        self.element_id = cast_int_array(element_id)
        self._mass = mass
        self.coord_id = coord_id
        self.node_id = cast_int_array(node_id)
        self.xyz_offset = xyz_offset
        #I11, I21, I22, I31, I32, I33 = I
        self.inertia = inertia

    def write(self, size: int=8) -> str:
        if len(self.element_id) == 0:
            return ''
        lines = []
        if size == 8:
            print_card = print_card_8

        is_xyz = self.xyz_offset.max() != 0. or self.xyz_offset.min() != 0.
        is_inertia = self.inertia.max() != 0. or self.inertia.min() != 0.
        element_str = array_str(self.element_id, size=size)
        node_str = array_str(self.node_id, size=size)
        coord_str = array_str(self.coord_id, size=size)
        if is_inertia:
            for eid, nid, cid, mass, xyz_offset, I in zip_longest(element_str, node_str, coord_str,
                                                                  self._mass, self.xyz_offset, self.inertia):
                list_fields = (['CONM2', eid, nid, cid, mass] +
                               list(xyz_offset) + [''] + list(I))
                lines.append(print_card(list_fields))
        elif is_xyz:
            for eid, nid, cid, mass, xyz_offset in zip_longest(element_str, node_str, coord_str,
                                                               self._mass, self.xyz_offset):
                list_fields = ['CONM2', eid, nid, cid, mass] + list(xyz_offset)
                lines.append(print_card(list_fields))
        else:
            for eid, nid, cid, mass in zip_longest(element_str, node_str, coord_str, self._mass):
                list_fields = ['CONM2', eid, nid, cid, mass]
                lines.append(print_card(list_fields))
        return ''.join(lines)

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        cid = self.model.coord.coord_id
        ucid = np.unique(self.coord_id)
        if ucid[0] == -1:
            ucid = ucid[1:]
            geom_check(self,
                       missing,
                       node=(nid, self.node_id),
                       coord=(cid, ucid))

    def mass(self) -> np.ndarray:
        return self._mass

    def centroid(self) -> np.ndarray:
        nid = self.model.grid.node_id
        xyz = self.model.grid.xyz_cid0()
        inode = np.searchsorted(nid, self.node_id)
        assert np.array_equal(nid[inode], self.node_id)
        centroid = xyz[inode, :] + self.xyz_offset

        # handle cid=-1
        #ucid = np.unique(self.coord_id)
        icoord = np.where(self.coord_id == -1)
        centroid[icoord, :] = self.xyz_offset[icoord, :]
        return centroid

    def center_of_mass(self) -> np.ndarray:
        return self.centroid()
