# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import sys
from numpy import array, ndarray, unique

from pyNastran.bdf.fieldWriter import set_blank_if_default, print_card
from pyNastran.bdf.cards.baseCard import BaseCard, expand_thru, collapse_thru
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, blank, integer_or_string, components_or_blank)


class GRIDs(object):
    def __init__(self):
        #: stores the GRID comments
        self._comments = {}

        #: Node ID
        self.nid = []

        #: Grid point coordinate system
        self.cp = []

        #: node location in local frame
        self.xyz = []

        #: Analysis coordinate system
        self.cd = []

        #: SPC constraint
        self.ps = []

        #: Superelement ID
        self.seid = []
        
        # mapping dictionary from grid ID to index
        self._nidmap = {}
        
        # internal counter
        self._nnodes = 0

    def comment(self, nid):
        id = self._nidmap[nid]
        try:
            return self._comments[id]
        except KeyError:
            return ''

    def finalize(self):
        # convert the data to arrays for slicing
        # and memory storage
        self.nid = array(self.nid, dtype='int32')
        self.cp = array(self.cp, dtype='int32')
        self.xyz = array(self.xyz, dtype='float32')
        self.cd = array(self.cd, dtype='int32')
        self.ps = array(self.ps, dtype='string')
        self.seid = array(self.seid, dtype='int32')
        
    def get_xyz(self, model, nid):
        id = self._nidmap[nid]
        return self.xyz[id]

    def get_position(self, model, nid, debug=False):
        id = self._nidmap[nid]
        msg = 'which is required by GRID %i' % nid
        cp = self.cp[id]
        Cp = model.Coord(cp, msg)
        
        xyz = self.xyz[id]
        p, matrix = Cp.transformToGlobal(xyz, debug=debug)
        return p

    def resolve_grids(self, model, cid=0):
        """
        Puts all nodes in a common coordinate system (mainly for cid testing)

        Parameters
        ----------
        self: the object pointer
        cid:  the cid to resolve the nodes to (default=0)
        
        Notes
        -----
        loses association with previous coordinate systems so to go
        back requires another fem
        """
        assert cid in model.coords, ('cannot resolve nodes to '
                                     'cid=|%s| b/c it doesnt exist' % cid)
        for (i, nid) in enumerate(self.nid):
            p = self.PositionWRT(model, nid, cid)
            self.UpdatePosition(p, nid, cid)

    def get_positions(self, model, nids):
        if not isinstance(nids, ndarray): # prevents long msg
            nids = array(nids)

        # dict -> list
        #print("nids =", nids)
        ids = [self._nidmap[nid] for nid in nids]
        #print("ids  =", ids)
        cps = self.cp[ids]
        unique_cps = unique(cps)
        
        #d = dict((key, value) for (key, value) in sequence) # <= Python 2.6
        #d = {key: value for (key, value) in sequence} Python 2.7+
        msg = 'which is required by GRIDs %s' % str(nids)
        Cps = {cid: model.Coord(Cp, msg) for cid in unique_cps}
        
        positions = zeros((nnodes, 3), dtype='float32')
        for (id, cid) in zip(ids, cps):
            Cp = Cps[cid]
            xyz = self.xyz[id]
            assert len(xyz) == 3, len(xyz)
            p, matrix = Cp.transformToGlobal(xyz, debug=debug)
            positions[id, :] = p[:]
        return positions
        
    def PositionWRT(self, model, nid, cid, debug=False):
        """
        Gets the point which started in some arbitrary local coordinate
        system and returns it in the desired coordinate system

        Parameters
        ----------
        self :  the NODEs object
        model : the BDF model object
        cid :   int
                the desired coordinate ID
        nid :   int
                the desired node ID
        debug : bool
                developer debug (default=False)
        
        Returns
        -------
        position : (3,) array
                   The location of node ID nid in coordinate system C
        """
        #print("nid", nid)
        id = self._nidmap[nid]
        if cid == self.cp[id]:
            return self.xyz[id]
        
        msg = ' which is required by GRID nid=%s' % self.nid
        coordA = model.Coord(cid, msg=msg)

        xyz = self.xyz[id]
        assert len(xyz) == 3, len(xyz)
        # converting the xyz point arbitrary->global
        p, matrixDum = coordA.transformToGlobal(xyz, debug=debug)
        #print "wrt = ",p

        #msg = ' which is required by %s nid=%s' % (self.type, self.nid)
        coordB = model.Coord(cid, msg=msg)

        # a matrix global->local matrix is found
        pdum, matrix = coordB.transformToGlobal(
            array([1., 0., 0]), debug=debug)
        p2 = coordB.transformToLocal(p, matrix, debug=debug)
        return p2

    def resolve_grids(self, model, cid=0):
        """
        Puts all nodes in a common coordinate system (mainly for cid testing)

        Parameters
        ----------
        self: the object pointer
        cid:  the cid to resolve the nodes to (default=0)
        
        Notes
        -----
        Loses association with previous coordinate systems so to go back
        requires another fem
        """
        assert cid in model.coords, ('cannot resolve nodes to '
                                     'cid=|%s| b/c it doesnt exist' % cid)
        for i, nid in enumerate(self.nid):
            p = self.PositionWRT(model, nid, cid)
            self.UpdatePosition(p, nid, cid)

    def unresolve_grids(self, model, model_old):
        """
        Puts all nodes back to original coordinate system.

        Parameters
        ----------
        self:      the object pointer
        model_old: the old model that hasnt lost it's connection to
                          the node cids
        
        Warning
        -------
        hasnt been tested well...
        """
        debug = False
        for (nid, node_old) in model_old.nodes.iteritems():
            coord = node_old.cp
            xyz = self.xyz[id]
            (p, matrix) = coord.transformToGlobal(xyz, debug=debug)
            p2 = coord.transformToLocal(p, matrix, debug=debug)
            self.UpdatePosition(p2, nid, cid)

    def UpdatePosition(self, xyz, nid, cid):
        #print("xyz =", xyz)
        id = self._nidmap[nid]
        self.xyz[id, :] = array(xyz)
        self.cp[id] = cid
        #assert cid == 0

    #def cross_reference(self, model):
        #pass

    def add_grid(self, card=None, data=None, comment=''):
        """
        if coming from a BDF object, card is used
        if coming from the OP2, data is used
        """
        if comment:
            self._comments[self._nnodes] = comment
        if card:
            nid = integer(card, 1, 'nid')

            cp = integer_or_blank(card, 2, 'cp', 0)
            
            x = double_or_blank(card, 3, 'x1', 0.)
            y = double_or_blank(card, 4, 'x2', 0.)
            z = double_or_blank(card, 5, 'x3', 0.)
            xyz = [x, y, z]

            cd = integer_or_blank(card, 6, 'cd', 0)
            ps = components_or_blank(card, 7, 'ps', '')
            seid = integer_or_blank(card, 8, 'seid', 0)
            assert len(card) <= 9, 'len(GRID card) = %i' % len(card)
        else:
            nid = data[0]
            cp = data[1]
            xyz = data[2:5]
            cd = data[5]
            ps = data[6]
            seid = data[7]
            if ps == 0:
                ps = ''
            assert len(xyz) == 3, '%s len=%s' % (xyz, len(xyz))

        assert nid > 0, 'nid=%s' % nid
        assert cp >= 0, 'cp=%s' % cp
        assert cd >= -1, 'cd=%s' % cd
        assert seid >= 0, 'seid=%s' % seid
        
        #print("nid", self.nid)
        #print("nid", self._nidmap)
        if nid in self._nidmap:
            asdf
        self._nidmap[nid] = self._nnodes
        self.nid.append(nid)
        self.cp.append(cp)
        self.xyz.append(xyz)
        self.cd.append(cd)
        self.ps.append(ps)
        self.seid.append(seid)

        self._nnodes += 1
        #print("adding...")
        #print("%r" % self.write(nids=nid))
    
    def write(self, f=None, nids=None):
        """
        Writes the BDF fields for the given set of node IDs.
        
        Parameters
        ----------
        f    : file object
        nids : list/array/tuple
               write a subset of nids (default=None -> all nodes)
        
        Returns
        -------
        out : string
              a message of the unwritten cards
        
        Notes
        -----
          nids are useful when running interactively

        Examples
        --------
        >>> f = open('f.bdf', 'w')
        >>> all_nids = [1, 2, 3, 4]
        >>> bdf.nodes.write(f)  # no output, but all_nids is written to f
        >>>
        >>> nids = [1, 4]
        >>> bdf.nodes.write(f, nids) # write a subset of nodes to f
        >>>
        >>> bdf.nodes.write(nids=nids)
        GRID   1  0  0.1  0.2  0.0
        GRID   4  0  0.7  0.8  0.0

        >>> bdf.nodes.write()
        GRID   1  0  0.1  0.2  0.0
        GRID   2  0  0.3  0.4  0.0
        GRID   3  0  0.5  0.6  0.0
        GRID   4  0  0.7  0.8  0.0
        """
        #d = dict((key, value) for (key, value) in sequence) # <= Python 2.6
        #d = {key: value for (key, value) in sequence} Python 2.7+
        if nids is None:
            ids = range(self._nnodes)
        elif isinstance(nids, int):
            ids = [self._nidmap[nids]]  # nids = 1 integer
        else:
            # dict -> list
            ids = [self._nidmap[nid] for nid in nids]
        sys.stdout.flush()

        # storing the flag of whether to print or not
        file_exists = f is not None
        
        msg = []
        size = 8
        for i in ids:
            nid = self.nid[i]
            cp = self.cp[i]
            #print("*i=%i nid=%i" %(i, nid))
            (x, y, z) = self.xyz[i]  #  [i, :]
            cd = self.cd[i]
            ps = self.ps[i]
            seid = self.seid[i]
            comment = self.comment(nid)

            fields = ['GRID', nid, cp, x, y, z, cd, ps, seid]
            msg.append(print_card(fields, size=size))

            if file_exists and i % 1000 == 0:  # reduce memory usage
                f.write(''.join(msg))
                msg = ['']
        if file_exists:
            f.write(''.join(msg))
        return ''.join(msg)

    def __repr__(self):
        msg =  '<GRIDs> object\n'
        msg += '  nnodes = %i' % self._nnodes
        return msg