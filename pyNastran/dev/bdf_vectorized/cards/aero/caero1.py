import numpy as np
from numpy import zeros, arange, searchsorted, cross

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.assign_type import (integer, integer_or_blank,
    double_or_blank)
from pyNastran.dev.bdf_vectorized.cards.vectorized_card import VectorizedCard


class CAERO1(VectorizedCard):
    """
    Defines an aerodynamic macro element (panel) in terms of two leading edge
    locations and side chords. This is used for Doublet-Lattice theory for
    subsonic aerodynamics and the ZONA51 theory for supersonic aerodynamics.::

    +--------+-----+-----+----+-------+--------+--------+--------+------+
    |   1    |  2  |  3  | 4  |   5   |   6    |    7   |   8    |   9  |
    +--------+-----+-----+----+-------+--------+--------+--------+------+
    | CAERO1 | EID | PID | CP | NSPAN | NCHORD |  LSPAN | LCHORD | IGID |
    +--------+-----+-----+----+-------+--------+--------+--------+------+
    |        |  X1 | Y1  | Z1 | X12   | X4     | Y4     | Z4     | X43  |
    +--------+-----+-----+----+-------+--------+--------+--------+------+
    """
    type = 'CAERO1'
    def __init__(self, model):
        """
        ::

          1
          | \
          |   \
          |     \
          |      4
          |      |
          |      |
          2------3
        """
        VectorizedCard.__init__(self, model)

    def add_card(self, card, comment=''):
        i = self.i
        self.element_id[i] = integer(card, 1, 'element_id')
        self.property_id[i] = integer(card, 2, 'property_id')
        self.coord_id[i] = integer_or_blank(card, 3, 'cp', 0)

        self.nspan[i] = integer_or_blank(card, 4, 'nspan', 0)
        self.nchord[i] = integer_or_blank(card, 5, 'nchord', 0)

        #if self.nspan==0:
        self.lspan[i] = integer_or_blank(card, 6, 'lspan', 0)

        #if self.nchord==0:
        self.lchord[i] = integer_or_blank(card, 7, 'lchord', 0)

        self.igid[i] = integer(card, 8, 'igid')

        self.p1[i, :] = [double_or_blank(card, 9, 'x1', 0.0),
                         double_or_blank(card, 10, 'y1', 0.0),
                         double_or_blank(card, 11, 'z1', 0.0)]
        self.x12[i] = double_or_blank(card, 12, 'x12', 0.)

        self.p4[i, :] = [double_or_blank(card, 13, 'x4', 0.0),
                         double_or_blank(card, 14, 'y4', 0.0),
                         double_or_blank(card, 15, 'z4', 0.0)]
        self.x43[i] = double_or_blank(card, 16, 'x43', 0.)
        assert len(card) <= 17, 'len(CAERO1 card) = %i\ncard=%s' % (len(card), card)
        self.i += 1

    def allocate(self, ncards):
        self.n = ncards
        float_fmt = self.model.float_fmt

        #: Element ID
        self.element_id = zeros(ncards, 'int32')

        #: Property ID of a PAERO2
        self.property_id = zeros(ncards, 'int32')

        #: Coordinate system for locating point 1.
        self.coord_id = zeros(ncards, 'int32')

        #: Node IDs
        self.node_ids = zeros((ncards, 4), 'int32')

        #: Number of spanwise boxes; if a positive value is given NSPAN, equal divisions are
        #: assumed; if zero or blank, a list of division points is given at LSPAN, field 7.
        #: (Integer > 0)
        self.nspan = zeros(ncards, 'int32')

        #: Number of chordwise boxes; if a positive value is given NCHORD, equal divisions
        #: are assumed; if zero or blank, a list of division points is given at LCHORD, field 8.
        #: (Integer >= 0)
        self.nchord = zeros(ncards, 'int32')

        #: ID of an AEFACT entry containing a list of division points for spanwise boxes.
        #; Used only if NSPAN, field 5 is zero or blank. (Integer > 0)
        self.lspan = zeros(ncards, 'int32')

        #: ID of an AEFACT data entry containing a list of division points for chordwise boxes.
        #: Used only if NCHORD, field 6 is zero or blank. (Integer > 0)
        self.lchord = zeros(ncards, 'int32')

        #: Interference group identification; aerodynamic elements with different IGIDs are uncoupled.
        self.igid = zeros(ncards, 'int32')

        #: Location of points 1, in coordinate system CP. (Real)
        self.p1 = zeros((ncards, 3), float_fmt)

        #: Location of points 4, in coordinate system CP. (Real)
        self.p4 = zeros((ncards, 3), float_fmt)
        self.x12 = zeros(ncards, float_fmt)
        self.x43 = zeros(ncards, float_fmt)

    def build(self):
        if self.n:
            i = self.element_id.argsort()
            self.model.log.debug(i)
            self.element_id = self.element_id[i]
            self.property_id = self.property_id[i]
            self.coord_id = self.coord_id[i]
            self.p1 = self.p1[i, :]
            self.p4 = self.p4[i, :]
            self.x12 = self.x12[i]
            self.x43 = self.x43[i]
            self.nspan = self.nspan[i]
            self.nchord = self.nchord[i]
            self.lspan = self.lspan[i]
            self.lchord = self.lchord[i]
            self.igid = self.igid[i]

    def write_card(self, bdf_file, size=8, is_double=True, element_id=None):
        assert size in [8, 16], size
        assert is_double in [True, False], is_double
        if self.n:
            if element_id is None:
                i = arange(self.n)
            else:
                #assert len(unique(element_id))==len(element_id), unique(element_id)
                i = searchsorted(self.element_id, element_id)

            Cid = [cid if cid != 0 else '' for cid in self.coord_id[i]]
            Nspan = [nspan if nspan != 0 else '' for nspan in self.nspan[i]]
            Nchord = [nchord if nchord != 0 else '' for nchord in self.nchord[i]]

            Igid = self.igid[i]
            Lspan = [lspan if lspan != 0. else '' for lspan in self.lspan[i]]
            Lchord = [lchord if lchord != 0. else '' for lchord in self.lchord[i]]

            for (eid, pid, cid, nspan, nchord, lspan, lchord, igid,
                 p1, x12, p4, x43) in zip(self.element_id[i], self.property_id[i],
                 Cid, Nspan, Nchord, Lspan, Lchord, Igid,
                 self.p1[i, :], self.x12[i], self.p4[i, :], self.x43[i]):
                card = ['CAERO1', eid, pid, cid, nspan, nchord, lspan, lchord, igid,
                        p1[0], p1[1], p1[2], x12,
                        p4[0], p4[1], p4[2], x43,]
                if size == 8:
                    bdf_file.write(print_card_8(card))
                else:
                    bdf_file.write(print_card_16(card))

    def _verify(self, xref=True):
        self.get_area()
        self.get_normal()

    def rebuild(self):
        raise NotImplementedError()

    def get_area(self, element_id=None, total=False, node_ids=None, grids_cid0=None):
        """
        Gets the area of the CAERO1s on a total or per element basis.

        Parameters
        ----------
        element_id : (N, ) int ndarray; (default=None -> all)
            the elements to consider
        total : bool; default=False
           should the area be summed
        node_ids : (nnodes, ) int ndarray; (default=None -> computes)
            the GRID ids for grids_cid0
        grids_cid0 : (nnodes, 3) float ndarray; (default=None -> computes)
            the GRIDs in CORD2R=0

            Notes
        -----
        If node_ids is None, the positions of all the GRID cards
        must be calculated
        """
        area, _normal = self._area_normal(element_id=element_id,
            node_ids=node_ids, grids_cid0=grids_cid0,
            calculate_area=True,
            calculate_normal=False)

        if total:
            return area.sum()
        else:
            return area

    def get_normal(self, element_id=None, node_ids=None, grids_cid0=None):
        """
        Gets the normals of the CAERO1s on per element basis.

        Parameters
        ----------
        element_id : (N, ) int ndarray; (default=None -> all)
            the elements to consider
        node_ids : (nnodes, ) int ndarray; (None -> computes)
            the GRID ids for grids_cid0
        grids_cid0 : (nnodes, 3) float ndarray; (None -> computes)
            the GRIDs in CORD2R=0

        Notes
        -----
        If node_ids is None, the positions of all the GRID cards
        must be calculated.
        """
        area, normal = self._area_normal(element_id=element_id,
            node_ids=node_ids, grids_cid0=grids_cid0,
            calculate_area=False,
            calculate_normal=True)

        return normal
        #if total:
            #return area.sum()
        #else:
            #return area

    def _area_normal(self, element_id=None, node_ids=None, grids_cid0=None,
                          calculate_area=True,
                          calculate_normal=True):
        """
        Gets the area, and normals of the CAERO1s on a per
        element basis.

        Parameters
        ----------
        element_id : (N, ) int ndarray; (default=None -> all)
            the elements to consider
        node_ids : (nnodes, ) int ndarray; (None -> computes)
            the GRID ids for grids_cid0
        grids_cid0 : (nnodes, 3) float ndarray; (None -> computes)
            the GRIDs in CORD2R=0

        calculate_area : bool; default=True
            should the area be calculated
        calculate_normal : bool; default=True
            should the normals be calculated

        Notes
        -----
        If node_ids is None, the positions of all the GRID cards
        must be calculated.
        """
        if grids_cid0 is None:
            node_ids = self.model.grid.node_ids
            grids_cid0 = self.model.grid.position()

        p1 = self._positions(grids_cid0, self.node_ids[:, 0])
        p2 = self._positions(grids_cid0, self.node_ids[:, 1])
        p3 = self._positions(grids_cid0, self.node_ids[:, 2])
        p4 = self._positions(grids_cid0, self.node_ids[:, 3])

        v12 = p2 - p1
        v13 = p3 - p1
        v123 = cross(v12, v13)
        normi = np.linalg.norm(v123, axis=0)

        A = None
        normal = v123 / normi
        if calculate_area:
            A = 0.5 * normi
        return A, normal

    def _positions(self, nids_to_get, node_ids, grids_cid0):
        """
        Gets the positions of a list of nodes

        Parameters
        ----------
        nids_to_get : (n,) ndarray
            the node IDs to get
        node_ids : (n,) ndarray
            the node IDs that contains all the nids_to_get
        grids_cid0 : (n,3) ndarray
            the GRIDs in CID=0

        Returns
        -------
        grids2_cid0 : ???
            the corresponding positions of the requested GRIDs
        """
        grids2_cid_0 = grids_cid0[searchsorted(nids_to_get, node_ids), :]
        return grids2_cid_0
