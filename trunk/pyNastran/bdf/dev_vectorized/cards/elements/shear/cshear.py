import cStringIO
from numpy import array, zeros, arange, concatenate, searchsorted, where, unique

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank)

from pyNastran.bdf.dev_vectorized.cards.elements.element import Element

class CSHEAR(Element):
    type = 'CSHEAR'
    op2_id = 4

    def __init__(self, model):
        Element.__init__(self, model)

    def build(self):
        cards = self._cards
        ncards = len(cards)
        self.n = ncards
        if ncards:
            #: Element ID
            self.element_id = zeros(ncards, 'int32')
            #: Property ID
            self.property_id = zeros(ncards, 'int32')
            #: Node IDs
            self.node_ids = zeros((ncards, 4), 'int32')

            self.zoffset = zeros(ncards, 'int32')
            self.t_flag = zeros(ncards, 'int32')
            self.thickness = zeros((ncards, 4), self.model.float)

            for i, card in enumerate(cards):
                self.element_id[i] = integer(card, 1, 'eid')

                self.property_id[i] = integer(card, 2, 'pid')

                self.node_ids[i] = [integer(card, 3, 'n1'),
                                    integer(card, 4, 'n2'),
                                    integer(card, 5, 'n3'),
                                    integer(card, 6, 'n4')]
                assert len(card) <= 7, 'len(CSHEAR card) = %i' % len(card)
            i = self.element_id.argsort()
            self.element_id = self.element_id[i]
            self.property_id = self.property_id[i]
            self.node_ids = self.node_ids[i, :]
            self._cards = []
            self._comments = []
        else:
            self.element_id = array([], dtype='int32')
            self.property_id = array([], dtype='int32')

    def write_bdf(self, f, size=8, eids=None):
        if self.n:
            if eids is None:
                i = arange(self.n)
            else:
                assert len(unique(eids))==len(eids), unique(eids)
                i = searchsorted(self.element_id, eids)

            for (eid, pid, n) in zip(self.element_id[i], self.property_id[i], self.node_ids[i]):
                card = ['CSHEAR', eid, pid, n[0], n[1], n[2], n[3]]
                f.write(print_card(card, size=size))


    def _verify(self, xref=True):
        self.mass()
        self.area()
        self.normal()

    def rebuild(self):
        raise NotImplementedError()

    def mass(self, eids=None, total=False, node_ids=None, grids_cid0=None):
        """
        Gets the mass of the CSHEARs on a total or per element basis.

        :param self: the CSHEAR object
        :param eids: the elements to consider (default=None -> all)
        :param total: should the mass be summed (default=False)

        :param node_ids:   the GRIDs as an (N, )  NDARRAY (or None)
        :param grids_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        ..note:: If node_ids is None, the positions of all the GRID cards
                 must be calculated
        """
        mass, _area, _normal = self._mass_area_normal(eids=eids,
            node_ids=node_ids, grids_cid0=grids_cid0,
            calculate_mass=True, calculate_area=False,
            calculate_normal=False)

        if total:
            return mass.sum()
        else:
            return mass

    def area(self, eids=None, total=False, node_ids=None, grids_cid0=None):
        """
        Gets the area of the CSHEARs on a total or per element basis.

        :param self: the CSHEAR object
        :param eids: the elements to consider (default=None -> all)
        :param total: should the area be summed (default=False)

        :param node_ids:   the GRIDs as an (N, )  NDARRAY (or None)
        :param grids_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        ..note:: If node_ids is None, the positions of all the GRID cards
                 must be calculated
        """
        _mass, area, _normal = self._mass_area_normal(eids=eids,
            node_ids=node_ids, grids_cid0=grids_cid0,
            calculate_mass=False, calculate_area=True,
            calculate_normal=False)

        if total:
            return area.sum()
        else:
            return area

    def normal(self, eids=None, node_ids=None, grids_cid0=None):
        """
        Gets the normals of the CSHEARs on per element basis.

        :param self: the CSHEAR object
        :param eids: the elements to consider (default=None -> all)

        :param node_ids:   the GRIDs as an (N, )  NDARRAY (or None)
        :param grids_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        ..note:: If node_ids is None, the positions of all the GRID cards
                 must be calculated
        """
        _mass, area, normal = self._mass_area_normal(eids=eids,
            node_ids=node_ids, grids_cid0=grids_cid0,
            calculate_mass=False, calculate_area=False,
            calculate_normal=True)

        if total:
            return area.sum()
        else:
            return area

    def _mass_area_normal(self, eids=None, node_ids=None, grids_cid0=None,
                          calculate_mass=True, calculate_area=True,
                          calculate_normal=True):
        """
        Gets the mass, area, and normals of the CSHEARs on a per
        element basis.

        :param self: the CSHEAR object
        :param eids: the elements to consider (default=None -> all)

        :param node_ids:   the GRIDs as an (N, )  NDARRAY (or None)
        :param grids_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        :param calculate_mass: should the mass be calculated (default=True)
        :param calculate_area: should the area be calculated (default=True)
        :param calculate_normal: should the normals be calculated (default=True)

        ..note:: If node_ids is None, the positions of all the GRID cards
                 must be calculated
        """
        if nodes_cid0 is None:
            node_ids = self.model.grid.node_ids
            grids_cid0 = self.model.grid.get_positions()

        p1 = self._positions(grids_cid0, self.node_ids[:, 0])
        p2 = self._positions(grids_cid0, self.node_ids[:, 1])
        p3 = self._positions(grids_cid0, self.node_ids[:, 2])
        p4 = self._positions(grids_cid0, self.node_ids[:, 3])

        v12 = p2 - p1
        v13 = p3 - p1
        v123 = cross(v12, v13)
        if calculate_normal or calculate_area:
            normal = v123 / n
        if calculate_area:
            A = 0.5 * n
        if calculate_mass:
            t = self.model.pid.get_thickness(self.pid)
            massi = A * t
        return massi

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