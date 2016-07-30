from six.moves import zip, range
from numpy import array, zeros, arange, searchsorted, unique, cross
from numpy.linalg import norm

from pyNastran.bdf.dev_vectorized.cards.elements.shell.shell_element import ShellElement

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.bdf_interface.assign_type import (integer, integer_or_blank,
    double_or_blank)


class CQUAD4(ShellElement):
    type = 'CQUAD4'
    def __init__(self, model):
        ShellElement.__init__(self, model)

    def allocate(self, ncards):
        self.n = ncards
        float_fmt = self.model.float
        #: Element ID
        self.element_id = zeros(ncards, 'int32')
        #: Property ID
        self.property_id = zeros(ncards, 'int32')
        #: Node IDs
        self.node_ids = zeros((ncards, 4), 'int32')
        self.zoffset = zeros(ncards, 'int32')
        self.t_flag = zeros(ncards, 'int32')
        self.thickness = zeros((ncards, 4), float_fmt)

    def add(self, card, comment=''):
        i = self.i
        self.element_id[i] = integer(card, 1, 'eid')

        self.property_id[i] = integer(card, 2, 'pid')

        self.node_ids[i, :] = [
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3'),
            integer(card, 6, 'n4')
        ]

        #self.thetaMcid = integer_double_or_blank(card, 6, 'thetaMcid', 0.0)
        #self.zOffset = double_or_blank(card, 7, 'zOffset', 0.0)
        #blank(card, 8, 'blank')
        #blank(card, 9, 'blank')

        self.t_flag[i] = integer_or_blank(card, 10, 'TFlag', 0)
        self.thickness[i, :] = [
            double_or_blank(card, 11, 'T1', 1.0),
            double_or_blank(card, 12, 'T2', 1.0),
            double_or_blank(card, 13, 'T3', 1.0),
            double_or_blank(card, 14, 'T4', 1.0),
        ]
        self.i += 1

    def build(self):
        if self.n:
            i = self.element_id.argsort()
            self.element_id = self.element_id[i]
            self.property_id = self.property_id[i]
            self.node_ids = self.node_ids[i, :]
            self.thickness = self.thickness[i, :]
            self.t_flag = self.t_flag[i]
            assert self.node_ids.min() > 0
        else:
            self.element_id = array([], 'int32')
            self.property_id = array([], dtype='int32')

    def get_node_indicies(self, i=None):
        if i is None:
            i1 = self.model.grid.get_node_index_by_node_id(self.node_ids[:, 0])
            i2 = self.model.grid.get_node_index_by_node_id(self.node_ids[:, 1])
            i3 = self.model.grid.get_node_index_by_node_id(self.node_ids[:, 2])
            i4 = self.model.grid.get_node_index_by_node_id(self.node_ids[:, 3])
        else:
            i1 = self.model.grid.get_node_index_by_node_id(self.node_ids[i, 0])
            i2 = self.model.grid.get_node_index_by_node_id(self.node_ids[i, 1])
            i3 = self.model.grid.get_node_index_by_node_id(self.node_ids[i, 2])
            i4 = self.model.grid.get_node_index_by_node_id(self.node_ids[i, 3])
        return i1, i2, i3, i4

    #=========================================================================
    def _node_locations(self, xyz_cid0, i=None):
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.get_position_by_node_index()
        i1, i2, i3, i4 = self.get_node_indicies(i)
        n1 = xyz_cid0[i1, :]
        n2 = xyz_cid0[i2, :]
        n3 = xyz_cid0[i3, :]
        n4 = xyz_cid0[i4, :]
        return n1, n2, n3, n4

    def _mass_area_normal(self, element_id=None, node_ids=None, xyz_cid0=None,
                          calculate_mass=True, calculate_area=True,
                          calculate_normal=True):
        """
        Gets the mass, area, and normals of the CQUAD4s on a per
        element basis.

        :param element_id: the elements to consider (default=None -> all)

        :param xyz_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        :param calculate_mass: should the mass be calculated (default=True)
        :param calculate_area: should the area be calculated (default=True)
        :param calculate_normal: should the normals be calculated (default=True)

        .. note:: If node_ids is None, the positions of all the GRID cards
                  must be calculated
        """
        if element_id is None:
            element_id = self.element_id
            property_id = self.property_id
            i = None
        else:
            i = searchsorted(self.element_id, element_id)
            property_id = self.property_id[i]

        n1, n2, n3, n4 = self._node_locations(xyz_cid0, i)
        if calculate_mass:
            calculate_area = True
        normal, A = _cquad4_normal_A(n1, n2, n3, n4, calculate_area=calculate_area, normalize=True)

        massi = None
        if calculate_mass:
            #t = self.model.properties_shell.get_thickness_by_property_id(property_id)  # PSHELL
            #nsm = self.model.properties_shell.get_nonstructural_mass_by_property_id(property_id)  # PSHELL
            #rho = self.model.properties_shell.get_density_by_property_id(property_id)  # MAT1
            mpa = self.model.properties_shell.get_mass_per_area_by_property_id(property_id)
            assert mpa is not None
            #assert t is not None
            #assert nsm is not None
            #assert rho is not None
            #print(type(nsm)
            #print("A  =%s" % A)
            #print("rho=%s" % rho)
            #print("t  =%s" % t)
            #print("nsm=%s" % nsm)
            #massi = rho * A * t + nsm
            massi = mpa * A
        #print "massi =", massi
        return massi, A, normal

    def get_centroid_by_element_id(self, element_id=None, node_ids=None, xyz_cid0=None):
        if element_id is None:
            element_id = self.element_id
            i = None
        else:
            i = searchsorted(self.element_id, element_id)
        return self.get_centroid_by_element_index(i, xyz_cid0)

    def get_centroid_by_element_index(self, i=None, xyz_cid0=None):
        n1, n2, n3, n4 = self._node_locations(xyz_cid0, i)
        return (n1 + n2 + n3 + n4) / 4.

    def get_area_by_element_index(self, i=None, xyz_cid0=None):
        n1, n2, n3, n4 = self._node_locations(xyz_cid0, i)
        normal, A = _cquad4_normal_A(n1, n2, n3, n4, calculate_area=True, normalize=False)
        return A

    def get_normal_by_element_index(self, i=None, xyz_cid0=None):
        n1, n2, n3, n4 = self._node_locations(xyz_cid0, i)
        normal, A = _cquad4_normal_A(n1, n2, n3, n4, calculate_area=False, normalize=True)
        return normal

    #=========================================================================
    def write_card(self, f, size=8, element_id=None):
        if self.n:
            #print('    self.n = %s' % self.n)
            if element_id is None:
                i = arange(self.n)
            else:
                assert len(unique(element_id)) == len(element_id), unique(element_id)
                i = searchsorted(self.element_id, element_id)
            #if len(i) == 1:
                #i = i[0]
            #print('    iwrite = %s' % i)
            #print('    all_nodes = %s' % str(self.node_ids))

            for (eid, pid, n) in zip(self.element_id[i], self.property_id[i], self.node_ids[i]):
                #print('    n = %s' % n)
                card = ['CQUAD4', eid, pid, n[0], n[1], n[2], n[3]]
                f.write(print_card_8(card))

    def _verify(self, xref=True):
        self.get_mass_by_element_id()
        self.get_area_by_element_id()
        self.get_normal_by_element_id()

    def rebuild(self):
        raise NotImplementedError()

    def _positions(self, nids_to_get):
        """
        Gets the positions of a list of nodes

        :param nids_to_get:  the node IDs to get as an NDARRAY
        :param node_ids:     the node IDs that contains all the nids_to_get
                             as an NDARRAY
        :param grids_cid_0:  the GRIDs as an (N, )  NDARRAY

        :returns grids2_cid_0 : the corresponding positins of the requested
                                GRIDs
        """
        positions = self.model.grid.get_position_by_node_id(nids_to_get)
        #grids2_cid_0 = grids_cid0[searchsorted(node_ids, nids_to_get), :]
        #return grids2_cid_0
        return positions

    #def slice_by_index(self, i):
        #i = self._validate_slice(i)
        #obj = CQUAD4(self.model)
        #obj.n = len(i)
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        #obj.element_id = self.element_id[i]
        #obj.property_id = self.property_id[i]
        #obj.node_ids = self.node_ids[i, :]
        #obj.zoffset = self.zoffset[i]
        #obj.t_flag = self.t_flag[i]
        #obj.thickness = self.thickness[i, :]
        #return obj


def _cquad4_normal_A(n1, n2, n3, n4, calculate_area=True, normalize=True):
    v13 = n1 - n3
    v24 = n2 - n4
    normal = cross(v13, v24)

    A = None
    if calculate_area:
        n = norm(normal, axis=1)
        A = 0.5 * n
    elif normalize:
        n = norm(normal, axis=1)

    if normalize:
        #n = len(_norm)
        #for i in range(n):
            #normal[i] /= _norm[i]
        #normal /= n
        nx = len(n)
        normal /= n.reshape(nx, 1)

    return normal, A
