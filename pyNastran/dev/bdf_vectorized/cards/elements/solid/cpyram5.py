import numpy as np
from numpy import arange, cross, searchsorted, array
from numpy.linalg import norm  # type: ignore

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.bdf_interface.assign_type import integer

from pyNastran.dev.bdf_vectorized.cards.elements.solid.solid_element import SolidElement

def tri_area_centroid(n1, n2, n3):
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
    return area, centroid


class CPYRAM5(SolidElement):
    type = 'CPYRAM5'
    nnodes = 5
    def __init__(self, model):
        """
        Defines the CPYRAM5 object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        SolidElement.__init__(self, model)

    def add_card(self, card, comment=''):
        i = self.i
        eid = integer(card, 1, 'element_id')
        if comment:
            self.set_comment(eid, comment)

        #: Element ID
        self.element_id[i] = eid
        #: Property ID
        self.property_id[i] = integer(card, 2, 'property_id')
        #: Node IDs
        nids = array([
            integer(card, 3, 'node_id_1'),
            integer(card, 4, 'node_id_2'),
            integer(card, 5, 'node_id_3'),
            integer(card, 6, 'node_id_4'),
            integer(card, 7, 'node_id_5'),
        ], dtype='int32')
        assert 0 not in nids, '%s\n%s' % (nids, card)
        self.node_ids[i, :] = nids
        assert len(card) == 9, 'len(CPYRAM5 card) = %i\ncard=%s' % (len(card), card)
        self.i += 1

    def get_mass_matrix(self, i, model, positions, index0s):
        return M, dofs, nijv

    def get_stiffness_matrix(self, i, model, positions, index0s):
        return K, dofs, nijv

    def _verify(self, xref=True):
        eid = self.eid
        pid = self.Pid()
        nids = self.node_ids
        assert isinstance(eid, int)
        assert isinstance(pid, int)
        for i, nid in enumerate(nids):
            assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)
        if xref:
            c = self.centroid()
            v = self.volume()
            assert isinstance(v, float)
            for i in range(3):
                assert isinstance(c[i], float)

    def get_node_indicies(self, i=None):
        if i is None:
            i1 = self.model.grid.get_node_index_by_node_id(self.node_ids[:, 0])
            i2 = self.model.grid.get_node_index_by_node_id(self.node_ids[:, 1])
            i3 = self.model.grid.get_node_index_by_node_id(self.node_ids[:, 2])
            i4 = self.model.grid.get_node_index_by_node_id(self.node_ids[:, 3])
            i5 = self.model.grid.get_node_index_by_node_id(self.node_ids[:, 4])
        else:
            i1 = self.model.grid.get_node_index_by_node_id(self.node_ids[i, 0])
            i2 = self.model.grid.get_node_index_by_node_id(self.node_ids[i, 1])
            i3 = self.model.grid.get_node_index_by_node_id(self.node_ids[i, 2])
            i4 = self.model.grid.get_node_index_by_node_id(self.node_ids[i, 3])
            i5 = self.model.grid.get_node_index_by_node_id(self.node_ids[i, 4])
        return i1, i2, i3, i4, i5

    def _get_node_locations_by_index(self, i, xyz_cid0):
        """
        :param i:        None or an array of node IDs
        :param xyz_cid0: the node positions as a dictionary
        """
        grid = self.model.grid
        #get_node_index_by_node_id = self.model.grid.get_node_index_by_node_id
        node_ids = self.node_ids

        msg = ', which is required by %s' % self.type
        i1, i2, i3, i4, i5, i6 = self.get_node_indicies(i)
        n1 = xyz_cid0[i1, :]
        n2 = xyz_cid0[i2, :]
        n3 = xyz_cid0[i3, :]
        n4 = xyz_cid0[i4, :]
        n5 = xyz_cid0[i5, :]
        return n1, n2, n3, n4, n5

    def get_volume_by_element_id(self, element_id=None, xyz_cid0=None, total=False):
        """
        Gets the volume for one or more elements.

        Parameters
        ----------
        element_id : (nelements, ) int ndarray; default=None -> all
            the elements to consider
        xyz_cid0 : dict[int node_id] : (3, ) float ndarray xyz (default=None -> auto)
            the positions of the GRIDs in CID=0
        total : bool; default=False
            should the volume be summed

        .. note:: Volume for a CPENTA is the average area of two opposing faces
                  times the length between the centroids of those points
        """
        n1, n2, n3, n4, n5 = self._get_node_locations_by_element_id(element_id, xyz_cid0)
        (A1, c1) = tri_area_centroid(n1, n2, n3)
        (A2, c2) = tri_area_centroid(n4, n5, n6)
        volume = (A1 + A2) / 2. * norm(c1 - c2, axis=1)
        if total:
            volume = np.abs(volume).sum()
        else:
            volume = np.abs(volume)
        return volume

    def get_mass_by_element_id(self, element_id=None, xyz_cid0=None, total=False):
        """
        Gets the mass for one or more CTETRA elements.

        Parameters
        ----------
        element_id : (nelements, ) int ndarray; default=None -> all
            the elements to consider
        xyz_cid0 : dict[int node_id] : (3, ) float ndarray xyz (default=None -> auto)
            the positions of the GRIDs in CID=0
        total : bool; default=False
            should the centroid be summed
        """
        if element_id is None:
            element_id = self.element_id
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.get_position_by_node_index()

        V = self.get_volume_by_element_id(element_id, xyz_cid0)
        mid = self.model.properties_solid.get_material_id_by_property_id(self.property_id)
        rho = self.model.materials.get_density_by_material_id(mid)

        mass = V * rho
        if total:
            mass = mass.sum()
        return mass

    def get_centroid_volume_by_element_id(self, element_id=None, xyz_cid0=None, total=False):
        """
        Gets the centroid and volume for one or more elements.

        Parameters
        ----------
        element_id : (nelements, ) int ndarray; default=None -> all
            the elements to consider
        xyz_cid0 : dict[int node_id] : (3, ) float ndarray xyz (default=None -> auto)
            the positions of the GRIDs in CID=0
        :param total: should the volume be summed; centroid be averaged (default=False)

        .. seealso:: CPYRAM5.get_volume_by_element_id() and CPYRAM5.get_centroid_by_element_id
        """
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.get_position_from_node_index()

        n1, n2, n3, n4, n5 = self._get_node_locations_by_element_id(element_id, xyz_cid0)
        (A1, c1) = tri_area_centroid(n1, n2, n3)
        (A2, c2) = tri_area_centroid(n4, n5, n6)
        centroid = (c1 * A1 + c2 * A2) / (A1 + A2)
        volume = (A1 + A2) / 2. * norm(c1 - c2, axis=1)
        if total:
            centroid = centroid.mean()
            volume = np.abs(volume).sum()
        else:
            volume = np.abs(volume)
        assert volume.min() > 0.0, 'volume.min() = %f' % volume.min()
        return centroid, volume

    def get_centroid_by_element_id(self, element_id=None, xyz_cid0=None, total=False):
        """
        Gets the centroid for one or more elements.

        Parameters
        ----------
        element_id : (nelements, ) int ndarray; default=None -> all
            the elements to consider
        xyz_cid0 : dict[int node_id] : (3, ) float ndarray xyz (default=None -> auto)
            the positions of the GRIDs in CID=0
            the positions of the GRIDs in CID=0
        total : bool; default=False
            should the centroid be averaged
        """
        if element_id is None:
            element_id = self.element_id
        n1, n2, n3, n4, n5 = self._get_node_locations_by_element_id(element_id, xyz_cid0)
        (A1, c1) = tri_area_centroid(n1, n2, n3)
        (A2, c2) = tri_area_centroid(n4, n5, n6)
        centroid = (c1 * A1 + c2 * A2) / (A1 + A2)
        if total:
            centroid = centroid.mean(axis=0)
        return centroid

    def get_face_nodes(self, nid, nid_opposite):
        raise NotImplementedError()
        #nids = self.node_ids[:4]
        #indx = nids.index(nid_opposite)
        #nids.pop(indx)
        #return nids

    def write_card(self, bdf_file, size=8, element_id=None):
        if self.n:
            if element_id is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.element_id, element_id)

            for (eid, pid, n) in zip(self.element_id[i], self.property_id[i], self.node_ids[i, :]):
                if eid in self._comments:
                    bdf_file.write(self._comments[eid])
                card = ['CPYRAM5', eid, pid, n[0], n[1], n[2], n[3], n[4]]
                bdf_file.write(print_card_8(card))
