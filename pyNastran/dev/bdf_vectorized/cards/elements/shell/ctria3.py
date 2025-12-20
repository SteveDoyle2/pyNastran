from numpy import (array, zeros, arange, searchsorted,
                   unique, cross)
from numpy.linalg import norm  # type: ignore

from pyNastran.dev.bdf_vectorized.cards.elements.shell.shell_element import ShellElement

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.bdf_interface.assign_type import (integer, integer_or_blank,
    double_or_blank, blank)
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class CTRIA3(ShellElement):
    type = 'CTRIA3'
    def __init__(self, model):
        ShellElement.__init__(self, model)

    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            self.n = ncards
            float_fmt = self.model.float_fmt
            #: Element ID
            self.element_id = zeros(ncards, dtype='int32')
            #: Property ID
            self.property_id = zeros(ncards, dtype='int32')
            #: Node IDs
            self.node_ids = zeros((ncards, 3), dtype='int32')
            self.zoffset = zeros(ncards, dtype='int32')
            self.t_flag = zeros(ncards, dtype='int32')
            self.thickness = zeros((ncards, 3), dtype=float_fmt)
            #else:
                #self.element_id = array([], dtype='int32')
                #self.property_id = array([], dtype='int32')

    def add_card(self, card: BDFCard, comment: str=''):
        i = self.i
        self.element_id[i] = integer(card, 1, 'element_id')
        self.property_id[i] = integer(card, 2, 'property_id')
        self.node_ids[i] = [integer(card, 3, 'n1'),
                integer(card, 4, 'n2'),
                integer(card, 5, 'n3')]

        #self.theta_mcid = integer_double_or_blank(card, 6, 'theta_mcid', 0.0)
        #self.zoffset = double_or_blank(card, 7, 'zoffset', 0.0)
        blank(card, 8, 'blank')
        blank(card, 9, 'blank')

        self.t_flag[i] = integer_or_blank(card, 10, 'tflag', 0)
        self.thickness[i] = [
            double_or_blank(card, 11, 'T1', 1.0),
            double_or_blank(card, 12, 'T2', 1.0),
            double_or_blank(card, 13, 'T3', 1.0), ]
        self.i += 1

    def build(self):
        if self.n:
            i = self.element_id.argsort()
            self.element_id = self.element_id[i]
            self.property_id = self.property_id[i]
            self.node_ids = self.node_ids[i, :]
            self._cards = []
            self._comments = []
        else:
            self.element_id = array([], 'int32')
            self.property_id = array([], dtype='int32')

    def update(self, maps):
        """
        maps = {
            'node_id' : nid_map,
            'property' : pid_map,
        }
        """
        if self.n:
            eid_map = maps['element']
            nid_map = maps['node']
            pid_map = maps['property']
            for i, (eid, pid, nids) in enumerate(zip(self.element_id, self.property_id, self.node_ids)):
                #print(self.print_card(i))
                self.element_id[i] = eid_map[eid]
                self.property_id[i] = pid_map[pid]
                self.node_ids[i, 0] = nid_map[nids[0]]
                self.node_ids[i, 1] = nid_map[nids[1]]
                self.node_ids[i, 2] = nid_map[nids[2]]

    def write_card(self, bdf_file, element_id=None, size=8):
        if self.n:
            if element_id is None:
                i = arange(self.n)
            else:
                if isinstance(element_id, int):
                    element_id = [element_id]
                #assert len(unique(element_id)) == len(element_id), unique(element_id)
                i = searchsorted(self.element_id, element_id)
            for (eid, pid, n) in zip(self.element_id[i], self.property_id[i], self.node_ids[i]):
                card = ['CTRIA3', eid, pid, n[0], n[1], n[2]]
                bdf_file.write(print_card_8(card))

    def _verify(self):
        self.get_mass_by_element_id()
        self.get_area_by_element_id()
        self.get_normal_by_element_id()

    def rebuild(self):
        pass

    #def get_mass_by_element_id(self, element_ids=None, total=False, xyz_cid0=None):
        #"""
        #Gets the mass of the CTRIA3s on a total or per element basis.

        #Parameters
        #----------
        #element_id : (nelements, ) int ndarray; default=None -> all
            #the elements to consider
        #total : bool; default=False
            #should the mass be summed

        #:param node_ids: the GRIDs as an (N, )  NDARRAY (or None)
        #:param xyz_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        #.. note:: If node_ids is None, the positions of all the GRID cards
                  #must be calculated
        #"""
        #mass, _area, _normal = self._mass_area_normal(element_id=element_id,
            #xyz_cid0=xyz_cid0,
            #calculate_mass=True, calculate_area=False,
            #calculate_normal=False)

        #if total:
            #return mass.sum()
        #else:
            #return mass

    #def get_area_by_element_id(self, element_id=None, total=False, xyz_cid0=None):
        #"""
        #Gets the area of the CTRIA3s on a total or per element basis.

        #Parameters
        #----------
        #element_id : (nelements, ) int ndarray; default=None -> all
            #the elements to consider
        #:param total: should the area be summed (default=False)

        #:param node_ids:  the GRIDs as an (N, )  NDARRAY (or None)
        #:param xyz_cid0:  the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)
        #"""
        #_mass, area, _normal = self._mass_area_normal(element_id=element_id,
            #xyz_cid0=xyz_cid0,
            #calculate_mass=False, calculate_area=True,
            #calculate_normal=False)

        #if total:
            #return area.sum()
        #else:
            #return area

    #def get_normal_by_element_id(self, element_id=None, xyz_cid0=None):
        #"""
        #Gets the normals of the CTRIA3s on per element basis.

        #Parameters
        #----------
        #element_id : (nelements, ) int ndarray; default=None -> all
            #the elements to consider

        #:param node_ids:   the GRIDs as an (N, )  NDARRAY (or None)
        #:param grids_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        #.. note:: If node_ids is None, the positions of all the GRID cards
                  #must be calculated
        #"""
        #_mass, area, normal = self._mass_area_normal(element_id=element_id,
            #grids_cid0=grids_cid0,
            #calculate_mass=False, calculate_area=False,
            #calculate_normal=True)

        #if total:
            #return area.sum()
        #else:
            #return area

    def get_node_indicies(self, i=None):
        if i is None:
            i1 = self.model.grid.get_node_index_by_node_id(self.node_ids[:, 0])
            i2 = self.model.grid.get_node_index_by_node_id(self.node_ids[:, 1])
            i3 = self.model.grid.get_node_index_by_node_id(self.node_ids[:, 2])
        else:
            i1 = self.model.grid.get_node_index_by_node_id(self.node_ids[i, 0])
            i2 = self.model.grid.get_node_index_by_node_id(self.node_ids[i, 1])
            i3 = self.model.grid.get_node_index_by_node_id(self.node_ids[i, 2])
        return i1, i2, i3

    def _node_locations(self, xyz_cid0, i=None):
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.get_position_by_node_index()
        i1, i2, i3 = self.get_node_indicies(i)
        n1 = xyz_cid0[i1, :]
        n2 = xyz_cid0[i2, :]
        n3 = xyz_cid0[i3, :]
        return n1, n2, n3

    def _mass_area_normal(self, element_id=None, xyz_cid0=None,
                          calculate_mass=True, calculate_area=True,
                          calculate_normal=True):
        """
        Gets the mass, area, and normals of the CTRIA3s on a per
        element basis.

        Parameters
        ----------
        element_id : (nelements, ) int ndarray; default=None -> all
            the elements to consider

        :param node_ids: the GRIDs as an (N, )  NDARRAY (or None)
        :param xyz_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        calculate_mass : bool; default=True
            should the mass be calculated
        calculate_area : bool; default=True
            should the area be calculated
        calculate_normal : bool; default=True
            should the normals be calculated
        """
        if element_id is None:
            element_ids = self.element_id
            property_id = self.property_id
            i = None
        else:
            i = searchsorted(self.element_id, element_id)
            property_id = self.property_id[i]

        n1, n2, n3 = self._node_locations(xyz_cid0)
        if calculate_mass:
            calculate_area = True
        normal, A = _ctria3_normal_A(n1, n2, n3, calculate_area=calculate_area, normalize=True)

        massi = None
        if calculate_mass:
            #t = self.model.properties_shell.get_thickness_by_property_id(property_id)  # PSHELL
            #nsm = self.model.properties_shell.get_nonstructural_mass_by_property_id(property_id)  # PSHELL
            #rho = self.model.properties_shell.get_density_by_property_id(property_id)  # MAT1
            mpa = self.model.properties_shell.get_mass_per_area_by_property_id(property_id)
            assert mpa is not None
            #massi = rho * A * t + nsm
            massi = mpa * A
        return massi, A, normal

    def _positions(self, nids_to_get, node_ids, grids_cid0):
        """
        Gets the positions of a list of nodes

        Parameters
        ----------
        :param nids_to_get:  the node IDs to get as an NDARRAY
        :param node_ids:     the node IDs that contains all the nids_to_get
                             as an NDARRAY
        :param grids_cid_0:  the GRIDs as an (N, )  NDARRAY

        Returns
        -------
        grids2_cid_0 : (nnodes, 3) float ndarray
            the corresponding positions of the requested GRIDs
        """
        grids2_cid_0 = grids_cid0[searchsorted(nids_to_get, node_ids), :]
        return grids2_cid_0

    def displacement_stress(self):
        pass

    #def get_property_by_index(self, i):
        #pid = self.property_id[i]
        #return self.model.property_shell.get_property_by_index[pid]

    #def get_stiffness_matrices(self, model, positions, index0s):
        #out = []

        ## area coordinates
        #d2 = area
        #L1 = a1 + b1 * x + c1 * y / d2
        #L2 = a2 + b2 * x + c2 * y / d2
        #L3 = a3 + b3 * x + c3 * y / d2


        #for i in range(self.n):
            #K, dofs, nijv = self.get_stiffness_matrix(
                #i, model, positions, index0s)
            #out.append(K, dofs, nijv)
        ##self.add_stiffness(K, dofs, nijv)
        #return out

    def get_stiffness_matrix(self, i, model, positions, index0s):
        # Mindlin-Reissner (thick plate)

        # Kirchoff-Love (thin plate)

        eid = self.element_id[i]
        pid = self.property_id[i]
        #prop = self.get_property_by_index(i)
        n1, n2, n3 = self.node_ids[i, :]
        xyz1 = positions[n1]
        xyz2 = positions[n2]
        xyz3 = positions[n3]

        #mat = prop.get_material(pid)
        mat = model.properties_shell.get_material(pid)
        E = mat.E()
        nu = mat.Nu()
        #t = prop.get_thickness_by_property_id(pid)
        t = model.properties_shell.get_thickness_by_property_id(pid)

        # [k_bend] = [b.T] [D] [b] -> 5.16 - Przemieniecki
        # [b] - 5.139 Prz.
        # [D] - 2.28 - Prz.
        p, q, r = xyz1, xyz2, xyz3
        pq = q - p
        dpq = norm(pq)
        L = pq / dpq
        t = array([
            p[0] + L[0] * dpq,
            p[1] + L[1] * dpq,
            p[2] + L[2] * dpq,
        ])

        b = array([
            [y32,    0, -y32,   0,  y21,    0],
            [0,   -x32, 0,    x31,    0, -x21],
            [-x32, y32, x31, -y32, -x21,  y21],
        ], dtype='float64')

        D = E / (1 - nu ** 2) * array([
            [1, nu, 0],
            [nu, 1, 0],
            [0, 0, (1-nu)/2.],
        ], dtype='float64')
        k = np.cross(b.T, cross(D, b))

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

    def get_centroid_by_element_id(self, element_id=None, xyz_cid0=None):
        i = self.get_element_index_by_element_id(element_id)
        return self.get_centroid_by_element_index(i, xyz_cid0)

    def get_area_by_element_id(self, element_id=None, xyz_cid0=None):
        i = self.get_element_index_by_element_id(element_id)
        return self.get_area_by_element_index(i, xyz_cid0)

    def get_normal_by_element_id(self, element_id=None, xyz_cid0=None):
        i = self.get_element_index_by_element_id(element_id)
        return self.get_normal_by_element_index(i, xyz_cid0)

    def get_centroid_by_element_index(self, i=None, xyz_cid0=None):
        n1, n2, n3 = self._node_locations(xyz_cid0, i)
        return (n1 + n2 + n3) / 3.

    def get_area_by_element_index(self, i=None, xyz_cid0=None):
        n1, n2, n3 = self._node_locations(xyz_cid0, i)
        normal, A = _ctria3_normal_A(n1, n2, n3, calculate_area=True, normalize=False)
        return A

    def get_normal_by_element_index(self, i=None, xyz_cid0=None):
        n1, n2, n3 = self._node_locations(xyz_cid0, i)
        normal, A = _ctria3_normal_A(n1, n2, n3, calculate_area=False, normalize=True)
        return normal

    #def slice_by_index(self, i):
        #i = self._validate_slice(i)
        #obj = CTRIA3(self.model)
        #obj.n = len(i)
        ##obj._cards = self._cards[i]
        ##obj._comments = obj._comments[i]
        ##obj.comments = obj.comments[i]
        #obj.element_id = self.element_id[i]
        #obj.property_id = self.property_id[i]
        #obj.node_ids = self.node_ids[i, :]
        #obj.zoffset = self.zoffset[i]
        #obj.t_flag = self.t_flag[i]
        #obj.thickness = self.thickness[i, :]
        #return obj


def _ctria3_normal_A(n1, n2, n3, calculate_area=True, normalize=True):
    v12 = n2 - n1
    v13 = n3 - n1
    normal = np.cross(v12, v13)

    A = None
    if calculate_area:
        n = norm(normal, axis=1)
        A = 0.5 * n
    elif normalize:
        n = norm(normal, axis=1)

    if normalize:
        #n = len(_norm)
        #for i in range(n):
            #normal[i] /= n
        #print(normal.shape)
        #print(n.shape)
        nx = len(n)
        normal /= n.reshape(nx, 1)
    return normal, A
