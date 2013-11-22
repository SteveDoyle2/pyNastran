from numpy import array, zeros, concatenate, searchsorted, where, unique

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank)


class ElementsShell(object):
    def __init__(self, model):
        """
        Defines the ShellProperties object.

        :param self: the ShellProperties object
        :param model: the BDF object
        """
        self.model = model

        self.ctria3 = CTRIA3(self.model)
        #self.cquad4 = CQUAD4(self.model)
        #self.cquad8 = CQUAD8(self.model)

        self._ctria3 = []
        self._cquad4 = []
        self._cquad8 = []

        self._ctria3_comment = []
        self._cquad4_comment = []
        self._cquad8_comment = []

    def build(self):
        self.ctria3.build(self._ctria3)
        #self.cquad4.build(self._cquad4)
        #self.cquad8.build(self._cquad8)
        
        self._ctria3 = []
        self._cquad4 = []
        self._cquad8 = []

        self._ctria3_comment = []
        self._cquad4_comment = []
        self._cquad8_comment = []

        #eid = concatenate(pshell.pid, pcomp.pid)
        #unique_eids = unique(eid)
        #if unique_eids != len(eid):
        #    raise RuntimeError('There are duplicate CTRIA3/CQUAD4 IDs...')

    def rebuild(self):
        raise NotImplementedError()

    def add_ctria3(self, card, comment):
        self._ctria3.append(card)
        self._ctria3_comment.append(comment)

    #def add_pcomp(self, card, comment):
        #self._pcomp.append(card)
        #self._pcomp_comment.append(comment)

    #def add_pcomp(self, card, comment):
        #self._pshear.append(card)
        #self._pshear_comment.append(comment)

    def write_bdf(self, f, size=8, eids=None):
        f.write('$ELEMENTS\n')
        self.ctria3.write_bdf(f, size=size, eids=eids)
        #self.cquad4.write_bdf(f, size=size, eids=eids)
        #self.cquad8.write_bdf(f, size=size, eids=eids)

    def get_stats(self):
        msg = []
        types = [self.ctria3]
        for element in types:
            nele = len(element.eid)
            if nele:
                msg.append('  %-8s: %i' % (element.type, nele))
        return msg

    def _verify(self):
        self.ctria3._verify()
        #self.ctria6._verify()
        #self.cquad4._verify()
        #self.cquad8._verify()

class CTRIA3(object):
    type = 'CTRIA3'
    def __init__(self, model):
        self.model = model

    def build(self, cards):
        ncards = len(cards)
        #: Element ID
        self.eid = zeros(ncards, 'int32')
        #: Property ID
        self.pid = zeros(ncards, 'int32')
        #: Node IDs
        self.node_ids = zeros((ncards, 3), 'int32')
        
        self.zoffset = zeros(ncards, 'int32')
        self.t_flag = zeros(ncards, 'int32')
        self.thickness = zeros((ncards, 3), 'float32')
        
        for i, card in enumerate(cards):
            self.eid[i] = integer(card, 1, 'eid')
            
            self.pid[i] = integer(card, 2, 'pid')

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
        i = self.eid.argsort()
        self.eid = self.eid[i]
        self.pid = self.pid[i]
        self.node_ids = self.node_ids[i, :]

    def write_bdf(self, f, size=8, eids=None):
        if eids is None:
            for (eid, pid, n123) in zip(self.eid, self.pid, self.node_ids):
                card = ['CTRIA3', eid, pid, n123[0], n123[1], n123[2]]
                f.write(print_card(card, size=size))
        else:
            #eids.sort()
            #ie = eids.argsort()
            assert len(unique(eids))==len(eids), unique(eids)
            i = searchsorted(self.eid, eids)
            
            iis = []
            for ii, eid in enumerate(self.eid):
                if eid in eids:
                    iis.append(ii)
                #else:
                    #print 'eid', eid
                    #asdf
            #iis = array(iis)
            #print "i =", i
            #print "iis =", iis
            assert sum(i - iis) == 0
            #i = iis
            #print "ictria3 =", i, len(i), len(unique(i))
            #assert len(unique(i))==len(i), unique(i)
            #print "-------"
            for (eid, pid, n123) in zip(self.eid[i], self.pid[i], self.node_ids[i]):
                card = ['CTRIA3', eid, pid, n123[0], n123[1], n123[2]]
                f.write(print_card(card, size=size))


    def _verify(self):
        self.mass()
        self.area()
        self.normal()

    def rebuild(self):
        pass

    def mass(self, eids=None, total=False, node_ids=None, grids_cid0=None):
        """
        Gets the mass of the CTRIA3s on a total or per element basis.
        
        :param self: the CTRIA3 object
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
        Gets the area of the CTRIA3s on a total or per element basis.
        
        :param self: the CTRIA3 object
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
        Gets the normals of the CTRIA3s on per element basis.
        
        :param self: the CTRIA3 object
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
        Gets the mass, area, and normals of the CTRIA3s on a per
        element basis.
        
        :param self: the CTRIA3 object
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
            grids_cid0 = self.model.grid.position()

        p1 = self._positions(grids_cid0, self.node_ids[:, 0])
        p2 = self._positions(grids_cid0, self.node_ids[:, 1])
        p3 = self._positions(grids_cid0, self.node_ids[:, 2])
        
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