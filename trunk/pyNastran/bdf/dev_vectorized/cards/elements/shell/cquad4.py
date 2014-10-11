import cStringIO
from numpy import (array, zeros, arange, concatenate, searchsorted,
    where, unique, cross, dot)
from numpy.linalg import norm

from pyNastran.bdf.fieldWriter import print_card_8
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank)

class ShellElement(object):
    def __init__(self, model):
        self.model = model
        self.n = 0
        self._cards = []
        self._comments = []

    def add(self, card, comment):
        self._cards.append(card)
        self._comments.append(comment)

    def get_mass(self, element_ids=None, total=False, node_ids=None, xyz_cid0=None):
        """
        Gets the mass of the CQUAD4s on a total or per element basis.

        :param self: the CQUAD4 object
        :param element_ids: the elements to consider (default=None -> all)
        :param total: should the mass be summed (default=False)

        :param xyz_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        ..note:: If node_ids is None, the positions of all the GRID cards
                 must be calculated
        """
        mass, _area, _normal = self._mass_area_normal(element_ids=element_ids,
            xyz_cid0=xyz_cid0,
            calculate_mass=True, calculate_area=False,
            calculate_normal=False)

        if total:
            return mass.sum()
        else:
            #print('mass.shape = %s' % mass.shape)
            return mass

    def get_normal(self, element_ids=None, xyz_cid0=None):
        """
        Gets the normals of the CQUAD4s on per element basis.

        :param self: the CQUAD4 object
        :param element_ids: the elements to consider (default=None -> all)

        :param xyz_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        ..note:: If node_ids is None, the positions of all the GRID cards
                 must be calculated
        """
        _mass, area, normal = self._mass_area_normal(element_ids=element_ids,
            xyz_cid0=xyz_cid0,
            calculate_mass=False, calculate_area=False,
            calculate_normal=True)
        return normal

    def get_area(self, element_ids=None, total=False, xyz_cid0=None):
        """
        Gets the area of the CQUAD4s on a total or per element basis.

        :param self: the CQUAD4 object
        :param element_ids: the elements to consider (default=None -> all)
        :param total: should the area be summed (default=False)

        :param node_ids:   the GRIDs as an (N, )  NDARRAY (or None)
        :param grids_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        ..note:: If node_ids is None, the positions of all the GRID cards
                 must be calculated
        """
        _mass, area, _normal = self._mass_area_normal(element_ids=element_ids,
            xyz_cid0=xyz_cid0,
            calculate_mass=False, calculate_area=True,
            calculate_normal=False)
        if total:
            return area.sum()
        else:
            return area

    def get_thickness(self, element_ids=None):
        if element_ids is None:
            element_ids = self.element_id
            property_id = self.property_id
            i = None
        else:
            i = searchsorted(self.element_id, element_ids)
            property_id = self.property_id[i]
        #print 'element_ids =', element_ids
        #print 'property_ids =', property_id
        t = self.model.properties_shell.get_thickness(property_id)
        return t

    def get_nonstructural_mass(self, element_ids=None):
        if element_ids is None:
            element_ids = self.element_id
            property_id = self.property_id
            i = None
        else:
            i = searchsorted(self.element_id, element_ids)
            property_id = self.property_id[i]
        nsm = self.model.properties_shell.get_nonstructural_mass(property_id)
        return nsm

    def get_density(self, element_ids=None):
        if element_ids is None:
            element_ids = self.element_id
            property_id = self.property_id
            i = None
        else:
            i = searchsorted(self.element_id, element_ids)
            property_id = self.property_id[i]
        #print('density - element_ids = %s' % element_ids)
        density = self.model.properties_shell.get_density(property_id)
        #print('density_out = %s' % density)
        return density


class CQUAD4(ShellElement):
    type = 'CQUAD4'
    def __init__(self, model):
        ShellElement.__init__(self, model)

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
            self.node_ids = zeros((ncards, 4), 'int32')

            self.zoffset = zeros(ncards, 'int32')
            self.t_flag = zeros(ncards, 'int32')
            self.thickness = zeros((ncards, 4), float_fmt)

            for i, card in enumerate(cards):
                self.element_id[i] = integer(card, 1, 'eid')

                self.property_id[i] = integer(card, 2, 'pid')

                self.node_ids[i, :] = [integer(card, 3, 'n1'),
                                    integer(card, 4, 'n2'),
                                    integer(card, 5, 'n3'),
                                    integer(card, 6, 'n4')]

                #self.thetaMcid = integer_double_or_blank(card, 6, 'thetaMcid', 0.0)
                #self.zOffset = double_or_blank(card, 7, 'zOffset', 0.0)
                #blank(card, 8, 'blank')
                #blank(card, 9, 'blank')

                #self.TFlag = integer_or_blank(card, 10, 'TFlag', 0)
                #self.T1 = double_or_blank(card, 11, 'T1', 1.0)
                #self.T2 = double_or_blank(card, 12, 'T2', 1.0)
                #self.T3 = double_or_blank(card, 13, 'T3', 1.0)
            i = self.element_id.argsort()
            self.element_id = self.element_id[i]
            self.property_id = self.property_id[i]
            self.node_ids = self.node_ids[i, :]
            assert self.node_ids.min() > 0
            self._cards = []
            self._comments = []

    #=========================================================================
    def _node_locations(self, xyz_cid0, i=None):
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.get_positions()
        if i is None:
            n1 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 0]), :]
            n2 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 1]), :]
            n3 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 2]), :]
            n4 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 3]), :]
        else:
            n1 = xyz_cid0[self.model.grid.index_map(self.node_ids[i, 0]), :]
            n2 = xyz_cid0[self.model.grid.index_map(self.node_ids[i, 1]), :]
            n3 = xyz_cid0[self.model.grid.index_map(self.node_ids[i, 2]), :]
            n4 = xyz_cid0[self.model.grid.index_map(self.node_ids[i, 3]), :]
        return n1, n2, n3, n4

    def _mass_area_normal(self, element_ids=None, node_ids=None, xyz_cid0=None,
                          calculate_mass=True, calculate_area=True,
                          calculate_normal=True):
        """
        Gets the mass, area, and normals of the CQUAD4s on a per
        element basis.

        :param self: the CQUAD4 object
        :param element_ids: the elements to consider (default=None -> all)

        :param xyz_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        :param calculate_mass: should the mass be calculated (default=True)
        :param calculate_area: should the area be calculated (default=True)
        :param calculate_normal: should the normals be calculated (default=True)

        ..note:: If node_ids is None, the positions of all the GRID cards
                 must be calculated
        """
        if element_ids is None:
            element_ids = self.element_id
            property_id = self.property_id
            i = None
        else:
            i = searchsorted(self.element_id, element_ids)
            property_id = self.property_id[i]

        n1, n2, n3, n4 = self._node_locations(xyz_cid0, i)
        #print('n1 = %s' % n1)
        #print('n2 = %s' % n2)
        #print('n3 = %s' % n3)
        #print('n4 = %s' % n4)
        v13 = n1 - n3
        v24 = n2 - n4
        normal = cross(v13, v24)

        #if calculate_normal or calculate_area or calculate_mass:
        _norm = norm(normal, axis=1)

        massi = None
        A = None

        n = len(_norm)
        for i in xrange(n):
            normal[i] /= _norm[i]

        if calculate_area or calculate_mass:
            A = 0.5 * _norm
        #print "A =", A
        if calculate_mass:
            #t = self.model.properties_shell.get_thickness(property_id)  # PSHELL
            #nsm = self.model.properties_shell.get_nonstructural_mass(property_id)  # PSHELL
            #rho = self.model.properties_shell.get_density(property_id)  # MAT1
            mpa = self.model.properties_shell.get_mass_per_area(property_id)
            assert mpa is not None
            #assert t is not None
            #assert nsm is not None
            #assert rho is not None
            #print type(nsm)
            #print "nsm =", nsm
            #print("A  =%s" % A)
            #print("rho=%s" % rho)
            #print("t  =%s" % t)
            #print("nsm=%s" % nsm)
            #massi = rho * A * t + nsm
            massi = mpa * A
        #print "massi =", massi
        return massi, A, normal

    def get_centroid(self, element_ids=None, node_ids=None, xyz_cid0=None):
        if element_ids is None:
            element_ids = self.element_id
            i = None
        else:
            i = searchsorted(self.element_id, element_ids)
        n1, n2, n3, n4 = self._node_locations(xyz_cid0, i)
        return (n1 + n2 + n3 + n4) / 4.

    #=========================================================================
    def write_bdf(self, f, size=8, element_ids=None):
        if self.n:
            if element_ids is None:
                i = arange(self.n)
            else:
                assert len(unique(element_ids))==len(element_ids), unique(element_ids)
                i = searchsorted(self.element_id, element_ids)

            for (eid, pid, n) in zip(self.element_id[i], self.property_id[i], self.node_ids[i]):
                card = ['CQUAD4', eid, pid, n[0], n[1], n[2], n[3]]
                f.write(print_card_8(card))

    def _verify(self, xref=True):
        self.get_mass()
        self.get_area()
        self.get_normal()

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
        positions = self.model.grid.get_positions(nids_to_get)
        #grids2_cid_0 = grids_cid0[searchsorted(node_ids, nids_to_get), :]
        #return grids2_cid_0
        return positions

    def __getitem__(self, i):
        obj = CQUAD4(self.model)
        obj.n = len(i)
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        obj.element_id = self.element_id[i]
        obj.property_id = self.property_id[i]
        obj.node_ids = self.node_ids[i, :]
        obj.zoffset = self.zoffset[i]
        obj.t_flag = self.t_flag[i]
        obj.thickness = self.thickness[i, :]
        return obj

    def __repr__(self):
        f = cStringIO.StringIO()
        f.write('<QUAD4 object> n=%s\n' % self.n)
        self.write_bdf(f)
        return f.getvalue()