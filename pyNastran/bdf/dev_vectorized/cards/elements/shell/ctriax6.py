from six.moves import zip

import numpy as np
from numpy import array, zeros, arange, searchsorted, cross

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.assign_type import (integer, integer_or_blank,
    double_or_blank)


class CTRIAX6(object):
    type = 'CTRIAX6'
    def __init__(self, model):
        self.model = model
        self.n = 0
        self._cards = []

    def add(self, card, comment=''):
        self._cards.append(card)

    def build(self):
        cards = self._cards
        ncards = len(cards)
        self.n = ncards
        if ncards:
            float_fmt = self.model.float
            #: Element ID
            self.element_id = zeros(ncards, 'int32')
            #: Material ID
            self.material_id = zeros(ncards, 'int32')
            #: Node IDs
            self.node_ids = zeros((ncards, 6), 'int32')

            self.theta = zeros(ncards, float_fmt)

            for i, card in enumerate(cards):
                self.element_id[i] = integer(card, 1, 'element_id')

                self.material_id[i] = integer(card, 2, 'material_id')

                nids = [integer(card, 3, 'n1'),
                        integer_or_blank(card, 4, 'n2'),
                        integer(card, 5, 'n3'),
                        integer_or_blank(card, 6, 'n4'),
                        integer(card, 7, 'n5'),
                        integer_or_blank(card, 8, 'n6')]
                self.node_ids[i, :] = nids

                #: theta
                self.theta[i] = double_or_blank(card, 9, 'theta', 0.0)

                assert len(card) <= 10, 'len(CTRIAX6 card) = %i' % len(card)
            i = self.element_id.argsort()
            self.element_id = self.element_id[i]
            self.material_id = self.material_id[i]
            self.node_ids = self.node_ids[i, :]
            self.theta = self.theta[i]
            self._cards = []
            self._comments = []
        else:
            self.element_id = array([], dtype='int32')

    def get_index_by_element_id(self, element_id=None):
        if element_id is None:
            return arange(self.n)
        return searchsorted(element_id, self.element_id)

    def write_card(self, f, size=8, element_id=None):
        if self.n:
            i = self.get_index_by_element_id(element_id)
            Theta = [theta if theta != 0.0 else '' for theta in self.theta[i]]
            N0 = [n if n != 0 else '' for n in self.node_ids[i, 0]]
            N1 = [n if n != 0 else '' for n in self.node_ids[i, 1]]
            N2 = [n if n != 0 else '' for n in self.node_ids[i, 2]]
            N3 = [n if n != 0 else '' for n in self.node_ids[i, 3]]
            N4 = [n if n != 0 else '' for n in self.node_ids[i, 4]]
            N5 = [n if n != 0 else '' for n in self.node_ids[i, 5]]
            for (eid, mid, n0, n1, n2, n3, n4, n5, theta) in zip(self.element_id[i], self.material_id[i],
                    N0, N1, N2, N3, N4, N5, Theta):
                card = ['CTRIAX6', eid, mid, n0, n1, n2, n3, n4, n5, theta]
                if size == 8:
                    f.write(print_card_8(card))
                else:
                    f.write(print_card_16(card))

    def _verify(self):
        self.mass()
        self.area()
        self.normal()

    def rebuild(self):
        pass

    def mass(self, element_id=None, total=False, node_ids=None, grids_cid0=None):
        """
        Gets the mass of the CTRIAX6s on a total or per element basis.

        :param element_id: the elements to consider (default=None -> all)
        :param total: should the mass be summed (default=False)

        :param node_ids:   the GRIDs as an (N, )  NDARRAY (or None)
        :param grids_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        .. note:: If node_ids is None, the positions of all the GRID cards
                  must be calculated
        """
        mass, _area, _normal = self._mass_area_normal(element_id=element_id,
                                                      node_ids=node_ids, grids_cid0=grids_cid0,
                                                      calculate_mass=True, calculate_area=False,
                                                      calculate_normal=False)

        if total:
            return mass.sum()
        else:
            return mass

    def area(self, element_id=None, total=False, node_ids=None, grids_cid0=None):
        """
        Gets the area of the CTRIAX6s on a total or per element basis.

        :param element_id: the elements to consider (default=None -> all)
        :param total: should the area be summed (default=False)

        :param node_ids:   the GRIDs as an (N, )  NDARRAY (or None)
        :param grids_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        .. note:: If node_ids is None, the positions of all the GRID cards
                  must be calculated
        """
        _mass, area, _normal = self._mass_area_normal(element_id=element_id,
            node_ids=node_ids, grids_cid0=grids_cid0,
            calculate_mass=False, calculate_area=True,
            calculate_normal=False)

        if total:
            return area.sum()
        else:
            return area

    def normal(self, element_id=None, node_ids=None, grids_cid0=None):
        """
        Gets the normals of the CTRIAX6s on per element basis.

        :param element_id: the elements to consider (default=None -> all)

        :param node_ids:   the GRIDs as an (N, )  NDARRAY (or None)
        :param grids_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        .. note:: If node_ids is None, the positions of all the GRID cards
                  must be calculated
        """
        _mass, area, normal = self._mass_area_normal(element_id=element_id,
            node_ids=node_ids, grids_cid0=grids_cid0,
            calculate_mass=False, calculate_area=False,
            calculate_normal=True)

        if total:
            return area.sum()
        else:
            return area

    def _mass_area_normal(self, element_id=None, node_ids=None, grids_cid0=None,
                          calculate_mass=True, calculate_area=True,
                          calculate_normal=True):
        """
        Gets the mass, area, and normals of the CTRIAX6s on a per
        element basis.

        :param element_id: the elements to consider (default=None -> all)

        :param node_ids:   the GRIDs as an (N, )  NDARRAY (or None)
        :param grids_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        :param calculate_mass: should the mass be calculated (default=True)
        :param calculate_area: should the area be calculated (default=True)
        :param calculate_normal: should the normals be calculated (default=True)

        .. note:: If node_ids is None, the positions of all the GRID cards
                  must be calculated
        """
        if nodes_cid0 is None:
            node_ids = self.model.grid.node_ids
            grids_cid0 = self.model.grid.get_position_from_node_index()

        p1 = self._positions(grids_cid0, self.node_ids[:, 0])
        p2 = self._positions(grids_cid0, self.node_ids[:, 1])
        p3 = self._positions(grids_cid0, self.node_ids[:, 2])

        v12 = p2 - p1
        v13 = p3 - p1
        v123 = cross(v12, v13)
        normi = np.linalg.norm(v123)
        if calculate_normal or calculate_area:
            normal = v123 / normi
        if calculate_area:
            A = 0.5 * normi
        if calculate_mass:
            t = self.model.pid.get_thickness_by_property_id(self.pid)
            massi = A * t
        return massi, A, normal

    def _positions(self, nids_to_get, node_ids, grids_cid0):
        """
        Gets the positions of a list of nodes

        :param nids_to_get:  the node IDs to get as an NDARRAY
        :param node_ids:     the node IDs that contains all the nids_to_get
                             as an NDARRAY
        :param grids_cid_0:  the GRIDs as an (N, )  NDARRAY

        :returns grids2_cid_0 : the corresponding positins of the requested
                                GRIDs
        """
        grids2_cid_0 = grids_cid0[searchsorted(nids_to_get, node_ids), :]
        return grids2_cid_0

    def __repr__(self):
        return '<%s object; n=%s>' % (self.type, self.n)
