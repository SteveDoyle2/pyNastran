from itertools import izip

from numpy import array, zeros, arange, concatenate, searchsorted, where, unique

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank)


class CAERO1(object):
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
            
            for i, card in enumerate(cards):
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

                self.p1[i, :] = [double_or_blank(card, 9,  'x1', 0.0),
                                 double_or_blank(card, 10, 'y1', 0.0),
                                 double_or_blank(card, 11, 'z1', 0.0)]
                self.x12[i] = double_or_blank(card, 12, 'x12', 0.)

                self.p4[i, :] = [double_or_blank(card, 13, 'x4', 0.0),
                                 double_or_blank(card, 14, 'y4', 0.0),
                                 double_or_blank(card, 15, 'z4', 0.0)]
                self.x43[i] = double_or_blank(card, 16, 'x43', 0.)
                assert len(card) <= 17, 'len(CAERO1 card) = %i' % len(card)
            i = self.element_id.argsort()
            self.element_id = self.element_id[i]
            self.property_id = self.property_id[i]
            self.coord_id = self.coord_id[i]
            self.p1 = self.p1[i, :]
            self.p4 = self.p4[i, :]
            self.x12 = self.x12[i, :]
            self.x43 = self.x43[i, :]
            self.nspan = self.nspan[i]
            self.nchord = self.nchord[i]
            self.lspan = self.lspan[i]
            self.lchord = self.lchord[i]
            self.igid = self.igid[i]
            self._cards = []
            self._comments = []

    def write_bdf(self, f, size=8, element_ids=None):
        if self.n:
            if element_ids is None:
                i = arange(self.n)
            else:
                assert len(unique(element_ids))==len(element_ids), unique(element_ids)
                i = searchsorted(self.element_id, element_ids)

            Cid    = [cid    if cid    !=0 else '' for cid    in self.coord_id[i]]
            Nspan  = [nspan  if nspan  !=0 else '' for nspan  in self.nspan[i]]
            Nchord = [nchord if nchord !=0 else '' for nchord in self.nchord[i]]

            Igid = self.igid[i]
            Lspan  = [lspan  if lspan  !=0. else '' for lspan  in self.lspan[i]]
            Lchord = [lchord if lchord !=0. else '' for lchord in self.lchord[i]]

            for (eid, pid, cid, nspan, nchord, lspan, lchord, igid,
                    p1, x12, p4, x43) in izip(self.element_id[i], self.property_id[i],
                    Cid, Nspan, Nchord, Lspan, Lchord, Igid,
                    self.p1[i, :], self.x12[i], self.p4[i, :], self.x43[i]):
                card = ['CAERO1', eid, pid, cid, nspan, nchord, lspan, lchord, igid,
                            p1[0], p1[1], p1[2], x12,
                            p4[0], p4[1], p4[2], x43,]
                f.write(print_card(card, size=size))

    def _verify(self, xref=True):
        self.area()
        self.normal()

    def rebuild(self):
        raise NotImplementedError()

    def area(self, element_ids=None, total=False, node_ids=None, grids_cid0=None):
        """
        Gets the area of the CAERO1s on a total or per element basis.
        
        :param self: the CAERO1 object
        :param element_ids: the elements to consider (default=None -> all)
        :param total: should the area be summed (default=False)

        :param node_ids:   the GRIDs as an (N, )  NDARRAY (or None)
        :param grids_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)
        
        ..note:: If node_ids is None, the positions of all the GRID cards
                 must be calculated
        """
        area, _normal = self._area_normal(element_ids=element_ids,
            node_ids=node_ids, grids_cid0=grids_cid0,
            calculate_area=True,
            calculate_normal=False)

        if total:
            return area.sum()
        else:
            return area

    def normal(self, element_ids=None, node_ids=None, grids_cid0=None):
        """
        Gets the normals of the CAERO1s on per element basis.
        
        :param self: the CAERO1 object
        :param element_ids: the elements to consider (default=None -> all)

        :param node_ids:   the GRIDs as an (N, )  NDARRAY (or None)
        :param grids_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)
        
        ..note:: If node_ids is None, the positions of all the GRID cards
                 must be calculated
        """
        area, normal = self._area_normal(element_ids=element_ids,
            node_ids=node_ids, grids_cid0=grids_cid0,
            calculate_area=False,
            calculate_normal=True)

        if total:
            return area.sum()
        else:
            return area

    def _area_normal(self, element_ids=None, node_ids=None, grids_cid0=None,
                          calculate_area=True,
                          calculate_normal=True):
        """
        Gets the area, and normals of the CAERO1s on a per
        element basis.
        
        :param self: the CAERO1 object
        :param element_ids: the elements to consider (default=None -> all)

        :param node_ids:   the GRIDs as an (N, )  NDARRAY (or None)
        :param grids_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

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
        p4 = self._positions(grids_cid0, self.node_ids[:, 3])
        
        v12 = p2 - p1
        v13 = p3 - p1
        v123 = cross(v12, v13)
        
        A = None
        normal = v123 / n
        if calculate_area:
            A = 0.5 * n
        return A, normal
    
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