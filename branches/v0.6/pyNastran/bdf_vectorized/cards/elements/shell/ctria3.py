from itertools import izip

from numpy import array, zeros, arange, concatenate, searchsorted, where, unique

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank)


class CTRIA3(object):
    type = 'CTRIA3'
    def __init__(self, model):
        self.model = model
        self.n = 0
        self._cards = []
        self._comments = []

    def add(self, card, comment):
        self._cards.append(card)
        self._comments.append(comment)

    def build(self):
        cards = self._cards
        ncards = len(cards)
        self.n = ncards
        if ncards:
            float_fmt = self.model.float
            #: Element ID
            self.element_id = zeros(ncards, 'int32')
            #: Property ID
            self.property_id = zeros(ncards, 'int32')
            #: Node IDs
            self.node_ids = zeros((ncards, 3), 'int32')

            self.zoffset = zeros(ncards, 'int32')
            self.t_flag = zeros(ncards, 'int32')
            self.thickness = zeros((ncards, 3), float_fmt)

            for i, card in enumerate(cards):
                self.element_id[i] = integer(card, 1, 'elemeent_id')

                self.property_id[i] = integer(card, 2, 'property_id')

                self.node_ids[i] = [integer(card, 3, 'n1'),
                        integer(card, 4, 'n2'),
                        integer(card, 5, 'n3')]

                #self.thetaMcid = integer_double_or_blank(card, 6, 'thetaMcid', 0.0)
                #self.zOffset = double_or_blank(card, 7, 'zOffset', 0.0)
                blank(card, 8, 'blank')
                blank(card, 9, 'blank')

                #self.TFlag = integer_or_blank(card, 10, 'TFlag', 0)
                #self.T1 = double_or_blank(card, 11, 'T1', 1.0)
                #self.T2 = double_or_blank(card, 12, 'T2', 1.0)
                #self.T3 = double_or_blank(card, 13, 'T3', 1.0)
            i = self.element_id.argsort()
            self.element_id = self.element_id[i]
            self.property_id = self.property_id[i]
            self.node_ids = self.node_ids[i, :]
            self._cards = []
            self._comments = []

    def write_bdf(self, f, size=8, element_ids=None):
        if self.n:
            if element_ids is None:
                i = arange(self.n)
            else:
                assert len(unique(element_ids))==len(element_ids), unique(element_ids)
                i = searchsorted(self.element_id, element_ids)
            for (eid, pid, n) in izip(self.element_id[i], self.property_id[i], self.node_ids[i]):
                card = ['CTRIA3', eid, pid, n[0], n[1], n[2]]
                f.write(print_card(card, size=size))


    def _verify(self):
        self.mass()
        self.area()
        self.normal()

    def rebuild(self):
        pass

    def mass(self, element_ids=None, total=False, node_ids=None, grids_cid0=None):
        """
        Gets the mass of the CTRIA3s on a total or per element basis.

        :param self: the CTRIA3 object
        :param element_ids: the elements to consider (default=None -> all)
        :param total: should the mass be summed (default=False)

        :param node_ids:   the GRIDs as an (N, )  NDARRAY (or None)
        :param grids_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        ..note:: If node_ids is None, the positions of all the GRID cards
                 must be calculated
        """
        mass, _area, _normal = self._mass_area_normal(element_ids=element_ids,
            node_ids=node_ids, grids_cid0=grids_cid0,
            calculate_mass=True, calculate_area=False,
            calculate_normal=False)

        if total:
            return mass.sum()
        else:
            return mass

    def area(self, element_ids=None, total=False, node_ids=None, grids_cid0=None):
        """
        Gets the area of the CTRIA3s on a total or per element basis.

        :param self: the CTRIA3 object
        :param element_ids: the elements to consider (default=None -> all)
        :param total: should the area be summed (default=False)

        :param node_ids:   the GRIDs as an (N, )  NDARRAY (or None)
        :param grids_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        ..note:: If node_ids is None, the positions of all the GRID cards
                 must be calculated
        """
        _mass, area, _normal = self._mass_area_normal(element_ids=element_ids,
            node_ids=node_ids, grids_cid0=grids_cid0,
            calculate_mass=False, calculate_area=True,
            calculate_normal=False)

        if total:
            return area.sum()
        else:
            return area

    def normal(self, element_ids=None, node_ids=None, grids_cid0=None):
        """
        Gets the normals of the CTRIA3s on per element basis.

        :param self: the CTRIA3 object
        :param element_ids: the elements to consider (default=None -> all)

        :param node_ids:   the GRIDs as an (N, )  NDARRAY (or None)
        :param grids_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        ..note:: If node_ids is None, the positions of all the GRID cards
                 must be calculated
        """
        _mass, area, normal = self._mass_area_normal(element_ids=element_ids,
            node_ids=node_ids, grids_cid0=grids_cid0,
            calculate_mass=False, calculate_area=False,
            calculate_normal=True)

        if total:
            return area.sum()
        else:
            return area

    def _mass_area_normal(self, element_ids=None, node_ids=None, grids_cid0=None,
                          calculate_mass=True, calculate_area=True,
                          calculate_normal=True):
        """
        Gets the mass, area, and normals of the CTRIA3s on a per
        element basis.

        :param self: the CTRIA3 object
        :param element_ids: the elements to consider (default=None -> all)

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
            grids_cid0 = self.model.grid.position()

        p1 = self._positions(grids_cid0, self.node_ids[:, 0])
        p2 = self._positions(grids_cid0, self.node_ids[:, 1])
        p3 = self._positions(grids_cid0, self.node_ids[:, 2])

        v12 = p2 - p1
        v13 = p3 - p1
        v123 = cross(v12, v13)

        normal = v123 / n
        A = None
        massi = None
        if calculate_area or calculate_mass:
            A = 0.5 * n
        if calculate_mass:
            t = self.model.properties_shell.get_thickness(self.property_id)
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

    def displacement_stress(self):
        pass

    def stiffness(self, positions):
        # Mindlin-Reissner (thick plate)

        # Kirchoff-Love (thin plate)

        eid = self.element_id[i]
        pid = self.property_id[i]
        prop = self.get_property_by_index(i)
        n1, n2, n3 = self.node_ids[i, :]
        p1 = positions[n1]
        p2 = positions[n2]
        p3 = positions[n3]

        mat = prop.get_material(pid)
        E = mat.E()
        nu = mat.Nu()
        t = prop.get_thickness(pid)


        #====
        # bending
        Cb = [
            [1., nu, 0.],
            [nu, 1., 0.],
            [0., 0., (1.-nu)/2.],
        ]
        Cb *= E * h**3 / (12*(1-nu**2))

        #====
        # shear
        Cs = [
            [1., 0.],
            [0., 1.],
        ]
        Cs *= E * h * k / (2.*(1+nu))
