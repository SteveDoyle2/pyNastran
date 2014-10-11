import cStringIO
from itertools import izip

from numpy import (array, zeros, arange, concatenate, searchsorted, where,
                   unique, cross)
from numpy.linalg import norm

from pyNastran.bdf.dev_vectorized.cards.elements.shell.cquad4 import ShellElement

from pyNastran.bdf.fieldWriter import print_card_8
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank)

class CTRIA3(ShellElement):
    type = 'CTRIA3'
    def __init__(self, model):
        ShellElement.__init__(self, model)
        #self.model = model
        #self.n = 0
        #self._cards = []
        #self._comments = []

    #def add(self, card, comment):
        #self._cards.append(card)
        #self._comments.append(comment)

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
                f.write(print_card_8(card))

    def _verify(self):
        self.get_mass()
        self.get_area()
        self.get_normal()

    def rebuild(self):
        pass

    #def get_mass(self, element_ids=None, total=False, xyz_cid0=None):
        #"""
        #Gets the mass of the CTRIA3s on a total or per element basis.

        #:param self: the CTRIA3 object
        #:param element_ids: the elements to consider (default=None -> all)
        #:param total: should the mass be summed (default=False)

        #:param node_ids: the GRIDs as an (N, )  NDARRAY (or None)
        #:param xyz_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        #..note:: If node_ids is None, the positions of all the GRID cards
                 #must be calculated
        #"""
        #mass, _area, _normal = self._mass_area_normal(element_ids=element_ids,
            #xyz_cid0=xyz_cid0,
            #calculate_mass=True, calculate_area=False,
            #calculate_normal=False)

        #if total:
            #return mass.sum()
        #else:
            #return mass

    #def get_area(self, element_ids=None, total=False, xyz_cid0=None):
        #"""
        #Gets the area of the CTRIA3s on a total or per element basis.

        #:param self: the CTRIA3 object
        #:param element_ids: the elements to consider (default=None -> all)
        #:param total: should the area be summed (default=False)

        #:param node_ids:  the GRIDs as an (N, )  NDARRAY (or None)
        #:param xyz_cid0:  the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)
        #"""
        #_mass, area, _normal = self._mass_area_normal(element_ids=element_ids,
            #xyz_cid0=xyz_cid0,
            #calculate_mass=False, calculate_area=True,
            #calculate_normal=False)

        #if total:
            #return area.sum()
        #else:
            #return area

    #def get_normal(self, element_ids=None, xyz_cid0=None):
        #"""
        #Gets the normals of the CTRIA3s on per element basis.

        #:param self: the CTRIA3 object
        #:param element_ids: the elements to consider (default=None -> all)

        #:param node_ids:   the GRIDs as an (N, )  NDARRAY (or None)
        #:param grids_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        #..note:: If node_ids is None, the positions of all the GRID cards
                 #must be calculated
        #"""
        #_mass, area, normal = self._mass_area_normal(element_ids=element_ids,
            #grids_cid0=grids_cid0,
            #calculate_mass=False, calculate_area=False,
            #calculate_normal=True)

        #if total:
            #return area.sum()
        #else:
            #return area

    def _node_locations(self, xyz_cid0, i=None):
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.get_positions()
        if i is None:
            n1 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 0]), :]
            n2 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 1]), :]
            n3 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 2]), :]
        else:
            n1 = xyz_cid0[self.model.grid.index_map(self.node_ids[i, 0]), :]
            n2 = xyz_cid0[self.model.grid.index_map(self.node_ids[i, 1]), :]
            n3 = xyz_cid0[self.model.grid.index_map(self.node_ids[i, 2]), :]
        return n1, n2, n3

    def _mass_area_normal(self, element_ids=None, xyz_cid0=None,
                          calculate_mass=True, calculate_area=True,
                          calculate_normal=True):
        """
        Gets the mass, area, and normals of the CTRIA3s on a per
        element basis.

        :param self: the CTRIA3 object
        :param element_ids: the elements to consider (default=None -> all)

        :param node_ids: the GRIDs as an (N, )  NDARRAY (or None)
        :param xyz_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        :param calculate_mass: should the mass be calculated (default=True)
        :param calculate_area: should the area be calculated (default=True)
        :param calculate_normal: should the normals be calculated (default=True)
        """
        if element_ids is None:
            element_ids = self.element_id
            property_id = self.property_id
            i = None
        else:
            i = searchsorted(self.element_id, element_ids)
            property_id = self.property_id[i]

        n1, n2, n3 = self._node_locations(xyz_cid0)
        v12 = n2 - n1
        v13 = n3 - n1
        v123 = cross(v12, v13)
        n = norm(v123, axis=1)

        normal = v123 / n
        A = None
        massi = None
        if calculate_area or calculate_mass:
            A = 0.5 * n
        if calculate_mass:
            #t = self.model.properties_shell.get_thickness(property_id)  # PSHELL
            #nsm = self.model.properties_shell.get_nonstructural_mass(property_id)  # PSHELL
            #rho = self.model.properties_shell.get_density(property_id)  # MAT1
            mpa = self.model.properties_shell.get_mass_per_area(property_id)
            assert mpa is not None
            #massi = rho * A * t + nsm
            massi = mpa * A
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

    #def get_property_by_index(self, i):
        #pid = self.property_id[i]
        #return self.model.property_shell.get_property_by_index[pid]

    def get_stiffness(self, i, model, positions, index0s):
        # Mindlin-Reissner (thick plate)

        # Kirchoff-Love (thin plate)

        eid = self.element_id[i]
        pid = self.property_id[i]
        #prop = self.get_property_by_index(i)
        n1, n2, n3 = self.node_ids[i, :]
        p1 = positions[n1]
        p2 = positions[n2]
        p3 = positions[n3]

        #mat = prop.get_material(pid)
        mat = model.properties_shell.get_material(pid)
        E = mat.E()
        nu = mat.Nu()
        #t = prop.get_thickness(pid)
        t = model.properties_shell.get_thickness(pid)


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

    def get_centroid(self, element_ids=None, node_ids=None, xyz_cid0=None):
        if element_ids is None:
            element_ids = self.element_id
            i = None
        else:
            i = searchsorted(self.element_id, element_ids)
        n1, n2, n3 = self._node_locations(xyz_cid0, i)
        return (n1 + n2 + n3) / 3.

    def __getitem__(self, index):
        obj = CTRIA3(self.model)
        obj.n = len(index)
        #obj._cards = self._cards[index]
        #obj._comments = obj._comments[index]
        #obj.comments = obj.comments[index]
        obj.element_id = self.element_id[index]
        obj.property_id = self.property_id[index]
        obj.node_ids = self.node_ids[index, :]
        obj.zoffset = self.zoffset[index]
        obj.t_flag = self.t_flag[index]
        obj.thickness = self.thickness[index, :]
        return obj

    def __repr__(self):
        f = cStringIO.StringIO()
        f.write('<CTRIA3 object> n=%s\n' % self.n)
        self.write_bdf(f)
        return f.getvalue()