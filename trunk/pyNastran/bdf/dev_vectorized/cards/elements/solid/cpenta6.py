import cStringIO
from numpy import zeros, arange, unique, dot, cross, abs, searchsorted, array, where, asarray
from numpy.linalg import norm

from pyNastran.bdf.fieldWriter import print_card_8
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank, fields)


def area_centroid(n1, n2, n3):
    """
    Gets the area, :math:`A`, and centroid of a triangle.::

      1-----2
      |   /
      | /
      3
    """
    a = n1 - n2
    b = n2 - n3
    area = 0.5 * norm(cross(a, b), axis=1)
    centroid = (n1 + n2 + n3) / 3.

    return(area, centroid)


class CPENTA6(object):
    type = 'CPENTA6'
    op2_id = 62
    def __init__(self, model):
        """
        Defines the CPENTA6 object.

        :param self: the CPENTA6 object
        :param model: the BDF object
        """
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
            self.element_id = zeros(ncards, 'int32')
            self.property_id = zeros(ncards, 'int32')
            self.node_ids = zeros((ncards, 6), 'int32')

            comments = {}
            for i, card in enumerate(cards):
                comment = self._comments[i]
                eid = integer(card, 1, 'element_id')
                if comment:
                    self.comments[eid] = comment

                #: Element ID
                self.element_id[i] = eid
                #: Property ID
                self.property_id[i] = integer(card, 2, 'property_id')
                #: Node IDs
                nids = [
                    integer(card, 3, 'nid1'),
                    integer(card, 4, 'nid2'),
                    integer(card, 5, 'nid3'),
                    integer(card, 6, 'nid4'),
                    integer(card, 7, 'nid5'),
                    integer(card, 8, 'nid6'),
                ]
                assert 0 not in nids, '%s\n%s' % (nids, card)
                self.node_ids[i, :] = nids
                assert len(card) == 9, 'len(CPENTA6 card) = %i' % len(card)

            i = self.element_id.argsort()
            self.element_id = self.element_id[i]
            self.property_id = self.property_id[i]
            self.node_ids = self.node_ids[i, :]
            self._cards = []
            self._comments = []
        else:
            self.element_id = array([], dtype='int32')
            self.property_id = array([], dtype='int32')

    def _verify(self, xref=True):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.nodeIDs()
        assert isinstance(eid, int)
        assert isinstance(pid, int)
        for i,nid in enumerate(nids):
            assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)
        if xref:
            c = self.centroid()
            v = self.volume()
            assert isinstance(v, float)
            for i in range(3):
                assert isinstance(c[i], float)

    def _node_locations(self, xyz_cid0):
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.get_positions()
        n1 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 0]), :]
        n2 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 1]), :]
        n3 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 2]), :]
        n4 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 3]), :]
        n5 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 4]), :]
        n6 = xyz_cid0[self.model.grid.index_map(self.node_ids[:, 5]), :]
        return n1, n2, n3, n4, n5, n6

    def _get_area_centroid(self, element_ids, xyz_cid0):
        n1, n2, n3, n4, n5, n6 = self._node_locations(xyz_cid0)
        (A1, c1) = area_centroid(n1, n2, n3)
        (A2, c2) = area_centroid(n4, n5, n6)
        return (A1, A2, c1, c2)

    def get_volume(self, element_ids=None, xyz_cid0=None, total=False):
        """
        Gets the volume for one or more CPENTA6 elements.

        :param element_ids: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the volume be summed (default=False)

        ..note:: Volume for a CPENTA is the average area of two opposing faces
        times the length between the centroids of those points
        """
        (A1, A2, c1, c2) = self._get_area_centroid(element_ids, xyz_cid0)
        volume = (A1 + A2) / 2. * norm(c1 - c2, axis=1)
        if total:
            volume = abs(volume).sum()
        else:
            volume = abs(volume)
        return volume

    def get_centroid_volume(self, element_ids=None, xyz_cid0=None, total=False):
        """
        Gets the centroid and volume for one or more CPENTA6 elements.

        :param element_ids: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the volume be summed; centroid be averaged (default=False)

        ..see:: CPENTA6.get_volume() and CPENTA6.get_centroid() for more information.
        """
        if element_ids is None:
            element_ids = self.element_id
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.get_positions()

        (A1, A2, c1, c2) = self._area_centroid(element_ids, xyz_cid0)
        centroid = (c1 * A1 + c2 * A2) / (A1 + A2)
        volume = (A1 + A2) / 2. * norm(c1 - c2, axis=1)
        if total:
            centroid = centroid.mean()
            volume = abs(volume).sum()
        else:
            volume = abs(volume)
        assert volume.min() > 0.0, 'volume.min() = %f' % volume.min()
        return centroid, volume

    def get_centroid(self, element_ids=None, xyz_cid0=None, total=False):
        """
        Gets the centroid for one or more CPENTA6 elements.

        :param element_ids: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the centroid be averaged (default=False)
        """
        if element_ids is None:
            element_ids = self.element_id
        (A1, A2, c1, c2) = self._get_area_centroid(element_ids, xyz_cid0)
        centroid = (c1 * A1 + c2 * A2) / (A1 + A2)
        if total:
            centroid = centroid.mean(axis=0)
        return centroid

    def get_mass(self, element_ids=None, xyz_cid0=None, total=False):
        """
        Gets the mass for one or more CPENTA6 elements.

        :param element_ids: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the centroid be summed (default=False)
        """
        if element_ids is None:
            element_ids = self.element_id
        V = self.get_volume(element_ids, xyz_cid0)

        mid = self.model.properties_solid.get_mid(self.property_id)
        rho = self.model.materials.get_rho(mid)

        mass = V * rho
        if total:
            mass = mass.sum()
        return mass

    def get_face_nodes(self, nid, nid_opposite):
        raise NotImplementedError()
        nids = self.nodeIDs()[:4]
        indx = nids.index(nid_opposite)
        nids.pop(indx)
        return nids

    def get_mass_centroid_inertia(self, p=None, element_ids=None, xyz_cid0=None, total=False):
        """
        Calculates the mass, centroid, and (3, 3) moment of interia
        matrix.  Considers position, but not the (hopefully) small
        elemental term.

        :param p: the point to take the moment of inertia about (default=None -> origin)

        a  = integral(mu * (y^2 + z^2), dV)
        b  = integral(mu * (x^2 + z^2), dV)
        c  = integral(mu * (y^2 + y^2), dV)
        a' = integral(mu * (yz), dV)
        b' = integral(mu * (xz), dV)
        c' = integral(mu * (xy), dV)

        I = [ a  -b', -c']
            [-b'  b   -a']
            [-c' -a'   c ]

        Exact MOI for tetrahedron
        http://www.thescipub.com/abstract/?doi=jmssp.2005.8.11
        """
        if p is None:
            p = zeros(3, self.model.float)

        r = centroid - p  # 2D array - 1D array
        I = mass * r**2 # column vector * 2D array
        return mass, centroid, I

    def write_bdf(self, f, size=8, element_ids=None):
        if self.n:
            if element_ids is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.element_id, element_ids)
            for (eid, pid, n) in zip(self.element_id[i], self.property_id[i], self.node_ids[i, :]):
                card = ['CPENTA', eid, pid, n[0], n[1], n[2], n[3], n[4], n[5]]
                f.write(print_card_8(card))

    def get_density(self, element_ids=None):
        if element_ids is None:
            element_ids = self.element_id

        rho = []
        i = where(element_ids == self.element_id)[0]
        for pid in self.property_id[i]:
            rhoi = self.model.properties_solid.psolid.get_density(pid)
            rho += rhoi
        return rho

    def __getitem__(self, element_ids):
        i = searchsorted(self.element_id, element_ids)
        return self.slice_by_index(i)

    def slice_by_index(self, i):
        i = asarray(i)
        obj = CPENTA6(self.model)
        obj.n = len(i)
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        obj.element_id = self.element_id[i]
        obj.property_id = self.property_id[i]
        obj.node_ids = self.node_ids[i, :]
        return obj

    def __repr__(self):
        f = cStringIO.StringIO()
        f.write('<CPENTA6 object> n=%s\n' % self.n)
        self.write_bdf(f)
        return f.getvalue()