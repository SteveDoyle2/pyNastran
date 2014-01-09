# pylint: disable=C0103,R0902,R0904,R0914
"""
All static loads are defined in this file.  This includes:

 * LOAD
 * GRAV
 * ACCEL1
 * FORCE1
 * FORCE2
 * MOMENT
 * PLOAD
 * PLOAD2
 * PLOAD3
 * PLOAD4
 * PLOADX1

"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from itertools import izip

from numpy import array, cross, allclose, ndarray
from numpy.linalg import norm

from pyNastran.bdf.cards.loads.loads import Load, LoadCombination
from pyNastran.bdf.fieldWriter import set_blank_if_default
from ..baseCard import BaseCard, expand_thru, expand_thru_by
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank,
    string, string_or_blank,
    integer_or_string, fields,
    integer_string_or_blank)
from pyNastran.bdf.fieldWriter import print_card_8
from pyNastran.bdf.fieldWriter16 import print_card_16


class LOAD(LoadCombination):
    type = 'LOAD'

    def __init__(self, card=None, data=None, comment=''):
        LoadCombination.__init__(self, card, data)
        if comment:
            self._comment = comment

    def getLoadIDs(self):
        """
        .. note:: requires a cross referenced load
        """
        load_IDs = []
        for loads in self.loadIDs:
            for load in loads:
                #if isinstance(load,int):
                    #load_IDs += [load]

                if isinstance(load, LOAD):
                    lid = load.lid
                    if isinstance(lid, list):
                        load_IDs += load.lid
                    else:  # int
                        load_IDs += load.getLoadIDs()
                elif (isinstance(load, Force) or isinstance(load, Moment) or
                      isinstance(load, PLOAD4) or isinstance(load, GRAV)):
                    load_IDs += [load.lid]
                else:
                    msg = ('The getLoadIDs method doesnt support %s cards.\n'
                           '%s' % (load.__class__.__name__, str(load)))
                    raise NotImplementedError(msg)

        load_IDs = list(set(load_IDs))
        #print "load_IDs = ",load_IDs
        return load_IDs

    def getLoadTypes(self):
        """
        .. note:: requires a cross referenced load
        """
        loadTypes = []
        for loads in self.loadIDs:
            for load in loads:
                if isinstance(load, LOAD):
                    lid = load.lid
                    if isinstance(lid, list):
                        loadTypes += load.type
                    else:  # int
                        loadTypes += [load.type] + load.getLoadTypes()
                elif (isinstance(load, Force) or isinstance(load, Moment) or
                      isinstance(load, PLOAD4) or isinstance(load, GRAV)):
                    loadTypes += [load.type]
                else:
                    raise NotImplementedError(load)

        loadTypes = list(set(loadTypes))
        #print "loadTypes = ",loadTypes
        return loadTypes

    def writeCalculixGRAV(self, gx, gy, gz):
        msg = '*DLOAD\n'
        msg += 'AllElements,GRAV,%s,%s,%s\n' % (gx, gy, gz)
        return msg

    def writeCodeAsterLoad(self, model, gridWord='node'):
        loadIDs = self.getLoadIDs()
        loadTypes = self.getLoadTypes()

        #msg = '# Loads\n'
        msg = ''
        (typesFound, forceLoads, momentLoads,
         forceConstraints, momentConstraints,
         gravityLoads) = self.organizeLoads(model)

        nids = []
        for nid in forceLoads:
            nids.append(nid)
        for nid in momentLoads:
            nids.append(nid)

        if nids:
            msg += '# typesFound = %s\n' % (list(typesFound))
            msg += '# loadIDs    = %s\n' % (loadIDs)
            msg += "load_bc=AFFE_CHAR_MECA(MODELE=modmod,\n"
            #msg += "                      DDL_IMPO=(_F(GROUP_MA='Lleft',\n"
            msg += "                       FORCE_NODALE=(\n"

        #CHAR=AFFE_CHAR_MECA(MODELE=MODE,
        #             FORCE_NODALE=(
        #                     _F(NOEUD='N1',
        #                        FZ=-500.0),)

        #print("nids = ",nids)
        spaces = "                                 "
        for nid in sorted(nids):  # ,load in sorted(forceLoads.iteritems())
            #print("nid = ",nid)
            msg += spaces + "_F(NOEUD='%s%s',\n" % (gridWord, nid)
            #print "load = ",load

            if nid in forceLoads:
                force = forceLoads[nid]
                if abs(force[0]) > 0.:
                    msg += spaces + "  FX=%s,\n" % (force[0])
                if abs(force[1]) > 0.:
                    msg += spaces + "  FY=%s,\n" % (force[1])
                if abs(force[2]) > 0.:
                    msg += spaces + "  FZ=%s,\n" % (force[2])

            if nid in momentLoads:
                moment = momentLoads[nid]
                if abs(moment[0]) > 0.:
                    msg += spaces + "  MX=%s,\n" % (moment[0])
                if abs(moment[1]) > 0.:
                    msg += spaces + "  MY=%s,\n" % (moment[1])
                if abs(moment[2]) > 0.:
                    msg += spaces + "  MZ=%s,\n" % (moment[2])
            msg = msg[:-2]
            msg += '),\n'
            # finish the load

            #if moment in
            #msg += "                                   DX=0.0,\n"
            #msg += "                                   DY=0.0,\n"
            #msg += "                                   DZ=0.0,),\n"
            #msg += "                                _F(GROUP_MA='Lright',\n"
            #msg += "                                   DZ=0.0,),),\n"
        msg = msg[:-2]
        msg += ');\n'

        for gravityLoad in gravityLoads:
            msg += 'CA_GRAVITY(%s);\n' % (str(gravityLoad))
        return (msg, loadIDs, loadTypes)

    def getReducedLoads(self):
        """
        Get all load objects in a simplified form, which means all
        scale factors are already applied and only base objects
        (no LOAD cards) will be returned.

        .. todo:: lots more object types to support
        """
        scale_factors = []
        loads = []
        load_scale = self.scale # global
        for (loadsPack, i_scale) in izip(self.loadIDs, self.scaleFactors):
            scale = i_scale * load_scale # actual scale = global * local
            for load in loadsPack:
                if (isinstance(load, Force) or isinstance(load, Moment) or
                    isinstance(load, PLOAD4) or isinstance(load, GRAV)):
                    loads.append(load)
                    scale_factors.append(scale) # local
                elif isinstance(load, LOAD):
                    load_data = load.getReducedLoads()
                    (reduced_scale_factors, reduced_loads) = load_data

                    loads += reduced_loads
                    scale_factors += [scale * j_scale for j_scale
                                      in reduced_scale_factors]
                else:
                    msg = ('%s isnt supported in getReducedLoads method'
                           % load.__class__.__name__)
                    raise NotImplementedError(msg)

        return (scale_factors, loads)

    def organizeLoads(self, model):
        """
        Figures out magnitudes of the loads to be applied to the various nodes.
        This includes figuring out scale factors.
        """
        forceLoads = {}  # spc enforced displacement (e.g. FORCE=0)
        momentLoads = {}
        forceConstraints = {}
        momentConstraints = {}
        gravityLoads = []
        #print("self.loadIDs = ",self.loadIDs)

        typesFound = set()
        (scaleFactors, loads) = self.getReducedLoads()

        for (scaleFactor, load) in izip(scaleFactors, loads):
            #print("*load = ",load)
            out = load.transformLoad()
            typesFound.add(load.__class__.__name__)
            if isinstance(load, Force):
                (isLoad, node, vector) = out
                if isLoad:  # load
                    if node not in forceLoads:
                        forceLoads[node] = vector * scaleFactor
                    else:
                        forceLoads[node] += vector * scaleFactor
                else:  # constraint
                    if node not in forceLoads:
                        forceConstraints[node] = vector * scaleFactor
                    else:
                        forceConstraints[node] += vector * scaleFactor

            elif isinstance(load, Moment):
                (isLoad, node, vector) = out
                if isLoad:  # load
                    if node not in momentLoads:
                        momentLoads[node] = vector * scaleFactor
                    else:
                        momentLoads[node] += vector * scaleFactor
                else:  # constraint
                    if node not in momentLoads:
                        momentConstraints[node] = vector * scaleFactor
                    else:
                        momentConstraints[node] += vector * scaleFactor

            elif isinstance(load, PLOAD4):
                (isLoad, nodes, vectors) = out
                for (nid, vector) in izip(nodes, vectors):
                    # not the same vector for all nodes
                    forceLoads[nid] = vector * scaleFactor

            elif isinstance(load, GRAV):
                #(grav) = out
                gravityLoads.append(out * scaleFactor)  # grav
            else:
                msg = '%s not supported' % (load.__class__.__name__)
                raise NotImplementedError(msg)

        return (typesFound, forceLoads, momentLoads, forceConstraints,
                momentConstraints, gravityLoads)

    def rawFields(self):
        list_fields = ['LOAD', self.sid, self.scale]
        for (scaleFactor, loadID) in izip(self.scaleFactors, self.loadIDs):
            list_fields += [scaleFactor, self.LoadID(loadID)]
        return list_fields

    def reprFields(self):
        return self.rawFields()

    def write_bdf(self, size, card_writer):
        card = self.rawFields()
        if size == 8:
            return self.comment() + print_card_8(card)
        else:
            return self.comment() + print_card_16(card)
        #return self.comment() + card_writer(card)


class GRAV(BaseCard):
    """
    Defines acceleration vectors for gravity or other acceleration loading.::

      GRAV SID CID A     N1  N2 N3    MB
      GRAV 1   3   32.2 0.0 0.0 -1.0
    """
    type = 'GRAV'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            #: Set identification number
            self.sid = integer(card, 1, 'sid')
            #: Coordinate system identification number.
            self.cid = integer_or_blank(card, 2, 'cid', 0)
            #: scale factor
            self.scale = double(card, 3, 'scale')
            #: Acceleration vector components measured in coordinate system CID
            self.N = array([double_or_blank(card, 4, 'N1', 0.0),
                            double_or_blank(card, 5, 'N2', 0.0),
                            double_or_blank(card, 6, 'N3', 0.0)])
            #: Indicates whether the CID coordinate system is defined in the
            #: main Bulk Data Section (MB = -1) or the partitioned superelement
            #: Bulk Data Section (MB = 0). Coordinate systems referenced in the
            #: main Bulk Data Section are considered stationary with respect to
            #: the assembly basic coordinate system. See Remark 10.
            #: (Integer; Default = 0)
            self.mb = integer_or_blank(card, 7, 'mb', 0)
            assert len(card) <= 8, 'len(GRAV card) = %i' % len(card)
        else:
            self.sid = data[0]
            self.cid = data[1]
            self.a = data[2]
            self.N = array(data[3:6])
            self.mb = data[6]
            self.scale = 1.
            assert len(data) == 7

        assert not allclose(max(abs(self.N)), 0.), ('GRAV N is a zero vector, '
                                                   'N=%s' % (str(self.N)))

    def getLoads(self):
        return [self]

    def organizeLoads(self, model):
        typesFound = [self.type]
        forceLoads = {}
        momentLoads = {}
        forceConstraints = {}
        momentConstraints = {}
        gravityLoad = self.transformLoad()
        return (typesFound, forceLoads, momentLoads,
                           forceConstraints, momentConstraints,
                           gravityLoad)

    def transformLoad(self):
        g = self.GravityVector()
        (g2, matrix) = self.cid.transformToGlobal(g)
        return (g2)

    #def writeCodeAster(self,mag):
        #p = self.GravityVector()
        #msg = 'GRAV([%s,%s,%s])' %(p)
        #return msg

    def cross_reference(self, model):
        #print("xref GRAV")
        self.cid = model.Coord(self.cid)

    def Cid(self):
        if isinstance(self.cid, int):
            return self.cid
        return self.cid.cid

    def GravityVector(self):
        """returns the gravity vector in absolute coordinates"""
        (p, matrix) = self.cid.transformToGlobal(self.N)
        return self.scale * p

    def rawFields(self):
        N = list(self.N)
        list_fields = ['GRAV', self.sid, self.Cid(), self.scale] + N + [
                       self.mb]
        return list_fields

    def reprFields(self):
        N = []
        for n in self.N:
            N.append(set_blank_if_default(n, 0.0))

        mb = set_blank_if_default(self.mb, 0)
        list_fields = ['GRAV', self.sid, self.Cid(), self.scale] + N + [mb]
        return list_fields

    def write_bdf(self, size, card_writer):
        card = self.rawFields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class ACCEL1(BaseCard):
    """
    Acceleration Load

    Defines static acceleration loads at individual GRID points.
    """
    type = 'ACCEL1'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            #: Load set identification number (Integer>0)
            self.sid = integer(card, 1, 'sid')

            #: Coordinate system identification number. (Integer>0: Default=0)
            self.cid = integer_or_blank(card, 2, 'cid', 0)

            #: Acceleration vector scale factor. (Real)
            self.scale = double(card, 3, 'scale')

            #: Components of the acceleration vector measured in coordinate system
            #: CID. (Real; at least one Ni != 0)
            self.N = array([double_or_blank(card, 4, 'N1', 0.0),
                            double_or_blank(card, 5, 'N2', 0.0),
                            double_or_blank(card, 6, 'N3', 0.0)])
            assert max(abs(self.N)) > 0.

            nodes = fields(integer_or_string, card, 'node', i=9, j=len(card))
        else:
            raise NotImplementedError(data)

        #: nodes to apply the acceleration to
        self.nodes = expand_thru_by(nodes)

    def cross_reference(self, model):
        self.cid = model.Coord(self.cid)
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True)

    def Cid(self):
        if isinstance(self.cid, int):
            return self.cid
        return self.cid.cid

    def nodeIDs(self, nodes=None):  # this function comes from BaseCard.py
        """returns nodeIDs for repr functions"""
        if not nodes:
            nodes = self.nodes
        if isinstance(nodes[0], int):
            nodeIDs = [node for node in nodes]
        else:
            nodeIDs = [node.nid for node in nodes]
        assert 0 not in nodeIDs, 'nodeIDs = %s' % (nodeIDs)
        return nodeIDs

    def rawFields(self):
        list_fields = ['ACCEL1', self.sid, self.Cid(), self.scale,
                  self.N[0], self.N[1], self.N[2], None, None
                  ] + self.nodeIDs()
        return list_fields

    def write_bdf(self, size, card_writer):
        card = self.rawFields()
        return self.comment() + card_writer(card)


class Force(Load):
    """Generic class for all Forces"""
    type = 'Force'

    def __init__(self, card, data):
        Load.__init__(self, card, data)

    def normalize(self):
        """
        adjust the vector to a unit length
        scale up the magnitude of the vector
        """
        if self.mag != 0.0:  # enforced displacement
            normXYZ = norm(self.xyz)
            #mag = self.mag*normXYZ
            self.mag *= normXYZ
            try:
                self.xyz = self.xyz / normXYZ
            except FloatingPointError:
                msg = 'xyz = %s\n' % self.xyz
                msg += 'normXYZ = %s\n' % normXYZ
                msg += 'card =\n%s' % str(self)
                raise FloatingPointError(msg)

    def transformLoad(self):
        #print("self.xyz = ",self.xyz)
        (xyz, matrix) = self.cid.transformToGlobal(self.xyz)
        if self.mag > 0.:
            #print("mag=%s xyz=%s" % (self.mag, xyz))
            return (True, self.node, self.mag * xyz)  # load
        return (False, self.node, xyz)  # enforced displacement

    def getLoads(self):
        return [self]

    def F(self):
        return self.xyz * self.mag

    def getReducedLoads(self):
        scaleFactors = [1.]
        loads = self.F()
        return(scaleFactors, loads)

    def organizeLoads(self, model):
        (scaleFactors, forceLoads) = self.getReducedLoads()

        typesFound = [self.type]
        momentLoads = {}
        forceConstraints = {}
        momentConstraints = {}
        gravityLoads = []
        return (typesFound, forceLoads, momentLoads,
                forceConstraints, momentConstraints,
                gravityLoads)

    def write_bdf(self, size, card_writer):
        card = self.rawFields()
        return self.comment() + card_writer(card)


class Moment(Load):
    """Generic class for all Moments"""
    type = 'Moment'

    def __init__(self, card, data):
        Load.__init__(self, card, data)

    def normalize(self):
        """
        adjust the vector to a unit length
        scale up the magnitude of the vector
        """
        if self.mag != 0.0:  # enforced displacement
            normXYZ = norm(self.xyz)
            #mag = self.mag*normXYZ
            self.mag *= normXYZ
            self.xyz = self.xyz / normXYZ

    def transformLoad(self):
        #print("self.xyz = ",self.xyz)
        (xyz, matrix) = self.cid.transformToGlobal(self.xyz)
        if self.mag > 0.:
            #print("mag=%s xyz=%s" % (self.mag, xyz))
            return (True, self.node, self.mag * xyz)  # load
        return (False, self.node, xyz)  # enforced displacement

    def getLoads(self):
        return [self]

    def getReducedLoads(self):
        scaleFactors = [1.]
        loads = { self.node: self.M() }
        print(loads)
        return(scaleFactors, loads)

    def organizeLoads(self, model):
        (scaleFactors, momentLoads) = self.getReducedLoads()

        typesFound = [self.type]
        forceLoads = {}
        forceConstraints = {}
        momentConstraints = {}
        gravityLoads = []
        return (typesFound, forceLoads, momentLoads,
                forceConstraints, momentConstraints,
                gravityLoads)

    def M(self):
        return self.xyz * self.mag

    def write_bdf(self, size, card_writer):
        card = self.rawFields()
        return self.comment() + card_writer(card)


class FORCE(Force):
    type = 'FORCE'

    def __init__(self, card=None, data=None, comment=''):
        """
        ::

          FORCE          3       1            100.      0.      0.      1.
        """
        Force.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.sid = integer(card, 1, 'sid')
            self.node = integer(card, 2, 'node')
            self.cid = integer_or_blank(card, 3, 'cid', 0)
            self.mag = double(card, 4, 'mag')
            xyz = array([double_or_blank(card, 5, 'X1', 0.0),
                         double_or_blank(card, 6, 'X2', 0.0),
                         double_or_blank(card, 7, 'X3', 0.0)])
            assert len(card) <= 8, 'len(FORCE card) = %i' % len(card)
        else:
            self.sid = data[0]
            self.node = data[1]
            self.cid = data[2]
            self.mag = data[3]
            xyz = data[4:7]

        assert len(xyz) == 3, 'xyz=%s' % (xyz)
        self.xyz = array(xyz)

    def Cid(self):
        if isinstance(self.cid, int):
            return self.cid
        return self.cid.cid

    def F(self):
        return {self.node: self.mag * self.xyz}

    #def nodeID(self):
        #return self.node

    def cross_reference(self, model):
        """
        .. todo:: cross reference and fix repr function
        """
        self.cid = model.Coord(self.cid)

    def rawFields(self):
        list_fields = ['FORCE', self.sid, self.node, self.Cid(),
                       self.mag] + list(self.xyz)
        return list_fields

    def reprFields(self):
        cid = set_blank_if_default(self.Cid(), 0)
        list_fields = ['FORCE', self.sid, self.node, cid,
                       self.mag] + list(self.xyz)
        return list_fields

    def write_bdf(self, size, card_writer):
        card = self.rawFields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)
        #return self.comment() + card_writer(card)


class FORCE1(Force):
    """
    Defines a static concentrated force at a grid point by specification of a
    magnitude and two grid points that determine the direction.
    """
    type = 'FORCE1'

    def __init__(self, card=None, data=None, comment=''):
        Force.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.sid = integer(card, 1, 'sid')
            self.node = integer(card, 2, 'node')
            self.mag = double(card, 3, 'mag')
            self.g1 = integer(card, 4, 'g1')
            self.g2 = integer(card, 5, 'g2')
            assert len(card) == 6, 'len(FORCE1 card) = %i' % len(card)
        else:
            self.sid = data[0]
            self.node = data[1]
            self.mag = data[2]
            self.g1 = data[3]
            self.g2 = data[4]

    def cross_reference(self, model):
        """
        .. todo:: cross reference and fix repr function
        """
        self.node = model.Node(self.node)
        self.g1 = model.Node(self.g1)
        self.g2 = model.Node(self.g2)
        self.xyz = self.g2.Position() - self.g1.Position()
        self.normalize()

    def G1(self):
        if isinstance(self.g1, int) or isinstance(self.g1, float):
            return self.g1
        return self.g1.nid

    def G2(self):
        if isinstance(self.g2, int) or isinstance(self.g1, float):
            return self.g2
        return self.g2.nid

    def NodeID(self):
        if isinstance(self.node, int):
            return self.node
        return self.node.nid

    def rawFields(self):
        (node, g1, g2) = self.nodeIDs([self.node, self.g1, self.g2])
        list_fields = ['FORCE1', self.sid, node, self.mag, g1, g2]
        return list_fields

    def reprFields(self):
        return self.rawFields()

    def write_bdf(self, size, card_writer):
        card = self.rawFields()
        return self.comment() + card_writer(card)


class FORCE2(Force):
    """
    Defines a static concentrated force at a grid point by specification of a
    magnitude and four grid points that determine the direction.
    """
    type = 'FORCE2'

    def __init__(self, card=None, data=None, comment=''):
        """
        ::

          FORCE2 SID G F G1 G2 G3 G4
        """
        Force.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.sid = integer(card, 1, 'sid')
            self.node = integer(card, 2, 'node')
            self.mag = double(card, 3, 'mag')
            self.g1 = integer(card, 4, 'g1')
            self.g2 = integer(card, 5, 'g2')
            self.g3 = integer(card, 6, 'g3')
            self.g4 = integer(card, 7, 'g4')
            assert len(card) == 8, 'len(FORCE2 card) = %i' % len(card)
        else:
            self.sid = data[0]
            self.node = data[1]
            self.mag = data[2]
            self.g1 = data[3]
            self.g2 = data[4]
            self.g3 = data[5]
            self.g4 = data[6]

    def cross_reference(self, model):
        """
        .. todo:: cross reference and fix repr function
        """
        self.node = model.Node(self.node)

        v12 = model.Node(self.g2).Position() - model.Node(self.g1).Position()
        v34 = model.Node(self.g4).Position() - model.Node(self.g3).Position()
        v12 /= norm(v12)
        v34 /= norm(v34)
        self.xyz = cross(v12, v34)
        self.normalize()

    def NodeID(self):
        if isinstance(self.node, int):
            return self.node
        return self.node.nid

    def G1(self):
        if isinstance(self.g1, int):
            return self.g1
        return self.g1.nid

    def G2(self):
        if isinstance(self.g2, int):
            return self.g2
        return self.g2.nid

    def G3(self):
        if isinstance(self.g3, int):
            return self.g3
        return self.g3.nid

    def G4(self):
        if isinstance(self.g4, int):
            return self.g4
        return self.g4.nid

    def rawFields(self):
        (node, g1, g2, g3, g4) = self.nodeIDs([self.node, self.g1, self.g2, self.g3, self.g4])
        list_fields = ['FORCE2', self.sid, node, self.mag, g1, g2, g3, g4]
        return list_fields

    def reprFields(self):
        return self.rawFields()

    def write_bdf(self, size, card_writer):
        card = self.rawFields()
        return self.comment() + card_writer(card)


class MOMENT(Moment):
    type = 'MOMENT'

    def __init__(self, card=None, data=None, comment=''):
        """
        Defines a static concentrated moment at a grid point by specifying a
        scale factor and a vector that determines the direction.::

          MOMENT SID G CID M    N1  N2  N3
          MOMENT 2   5   6 2.9 0.0 1.0 0.0
        """
        Moment.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.sid = integer(card, 1, 'sid')
            self.node = integer(card, 2, 'node')
            self.cid = integer_or_blank(card, 3, 'cid', 0)
            self.mag = double(card, 4, 'mag')

            xyz = array([double_or_blank(card, 5, 'X1', 0.0),
                         double_or_blank(card, 6, 'X2', 0.0),
                         double_or_blank(card, 7, 'X3', 0.0)])
            assert len(card) <= 8, 'len(MOMENT card) = %i' % len(card)
        else:
            raise NotImplementedError(data)
        assert len(xyz) == 3, 'xyz=%s' % xyz
        self.xyz = xyz

    def Cid(self):
        if isinstance(self.cid, int):
            return self.cid
        return self.cid.cid

    def cross_reference(self, model):
        """
        .. todo:: cross reference and fix repr function
        """
        pass

    def rawFields(self):
        list_fields = ['MOMENT', self.sid, self.node, self.Cid(),
                  self.mag] + list(self.xyz)
        return list_fields

    def reprFields(self):
        cid = set_blank_if_default(self.Cid(), 0)
        list_fields = ['MOMENT', self.sid, self.node, cid,
                  self.mag] + list(self.xyz)
        return list_fields

    def write_bdf(self, size, card_writer):
        card = self.rawFields()
        return self.comment() + card_writer(card)


class MOMENT1(Moment):
    type = 'MOMENT1'

    def __init__(self, card=None, data=None, comment=''):
        """
        Defines a static concentrated moment at a grid point by specifying a
        magnitude and two grid points that determine the direction.::

          MOMENT1 SID G M G1 G2
        """
        Moment.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.sid = integer(card, 1, 'sid')
            self.node = integer(card, 2, 'node')
            self.mag = double(card, 3, 'mag')
            self.g1 = integer(card, 4, 'g1')
            self.g2 = integer(card, 5, 'g2')
            assert len(card) == 6, 'len(MOMENT1 card) = %i' % len(card)
        else:
            self.sid = data[0]
            self.node = data[1]
            self.mag = data[2]
            self.g1 = data[3]
            self.g2 = data[4]
            self.g3 = data[5]
            self.g4 = data[6]
            xyz = data[7:10]
            raise NotImplementedError('MOMENT1 is probably wrong')

        #assert len(xyz) == 3, 'xyz=%s' % (xyz)
        #self.xyz = array(xyz)
        self.xyz = None

    def cross_reference(self, model):
        """
        .. todo:: cross reference and fix repr function
        """
        self.node = model.Node(self.node)
        self.xyz = model.Node(self.g2).Position() - model.Node(self.g1).Position()
        self.normalize()

    def get_node_id(self):
        if isinstance(self.node, int):
            return self.node
        return self.node.nid

    def G1(self):
        if isinstance(self.g1, int):
            return self.g1
        return self.g1.nid

    def G2(self):
        if isinstance(self.g2, int):
            return self.g2
        return self.g2.nid

    def rawFields(self):
        node = self.get_node_id()
        g1 = self.G1()
        g2 = self.G2()
        list_fields = ['MOMENT1', self.sid, node, self.mag, g1, g2]
        return list_fields

    def reprFields(self):
        return self.rawFields()

    def write_bdf(self, size, card_writer):
        card = self.rawFields()
        return self.comment() + card_writer(card)


class MOMENT2(Moment):
    type = 'MOMENT2'

    def __init__(self, card=None, data=None, comment=''):
        """
        Defines a static concentrated moment at a grid point by specification
        of a magnitude and four grid points that determine the direction.::

          MOMENT2 SID G M G1 G2 G3 G4
        """
        Moment.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.sid = integer(card, 1, 'sid')
            self.node = integer(card, 2, 'node')
            self.mag = double(card, 3, 'mag')
            self.g1 = integer(card, 4, 'g1')
            self.g2 = integer(card, 5, 'g2')
            self.g3 = integer(card, 6, 'g3')
            self.g4 = integer(card, 7, 'g4')
            assert len(card) <= 8, 'len(MOMENT2 card) = %i' % len(card)
        else:
            self.sid = data[0]
            self.node = data[1]
            self.mag = data[2]
            self.g1 = data[3]
            self.g2 = data[4]
            self.g3 = data[5]
            self.g4 = data[6]
            xyz = data[7:10]
            self.xyz = array(xyz)
            assert len(xyz) == 3, 'xyz=%s' % xyz

    def cross_reference(self, model):
        """
        .. todo:: cross reference and fix repr function
        """
        self.node = model.Node(self.node)
        self.g1 = model.Node(self.g1)
        self.g2 = model.Node(self.g2)
        self.g3 = model.Node(self.g3)
        self.g4 = model.Node(self.g4)

        v12 = self.g2.Position() - self.g1.Position()
        v34 = self.g4.Position() - self.g3.Position()
        v12 = v12 / norm(v12)
        v34 = v34 / norm(v34)
        self.xyz = cross(v12, v34)

    def NodeID(self):
        if isinstance(self.node, int):
            return self.node
        return self.node.nid

    def G1(self):
        if isinstance(self.g1, int):
            return self.g1
        return self.g1.nid

    def G2(self):
        if isinstance(self.g2, int):
            return self.g2
        return self.g2.nid

    def G3(self):
        if isinstance(self.g3, int):
            return self.g3
        return self.g3.nid

    def G4(self):
        if isinstance(self.g4, int):
            return self.g4
        return self.g4.nid

    def rawFields(self):
        (node, g1, g2, g3, g4) = self.nodeIDs(nodes=[self.node, self.g1, self.g2,
                                               self.g3, self.g4])
        assert isinstance(g1, int), g1
        list_fields = ['MOMENT2', self.sid, node, self.mag, g1, g2, g3, g4]
        return list_fields

    def reprFields(self):
        return self.rawFields()

    def write_bdf(self, size, card_writer):
        card = self.rawFields()
        return self.comment() + card_writer(card)


class PLOAD(Load):
    type = 'PLOAD'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            self.sid = integer(card, 1, 'sid')
            self.p = double(card, 2, 'p')
            nodes = [integer(card, 3, 'n1'),
                     integer(card, 4, 'n2'),
                     integer(card, 5, 'n3')]
            n4 = integer_or_blank(card, 6, 'n5', 0)
            if n4:
                nodes.append(n4)
            #self.nodes = self._wipeEmptyFields(nodes)
            self.nodes = nodes
            assert len(card) <= 7, 'len(PLOAD card) = %i' % len(card)
        else:
            self.sid = data[0]
            self.p = data[1]
            self.nodes = data[2:]
            print("PLOAD = %s" % data)
            raise NotImplementedError('PLOAD')
        assert len(self.nodes) in [3, 4], 'nodes=%s' % self.nodes

    def cross_reference(self, model):
        """
        .. todo:: cross reference and fix repr function
        """
        pass

    def getLoads(self):
        return [self]

    def rawFields(self):
        list_fields = ['PLOAD', self.sid, self.p] + self.nodeIDs()
        return list_fields

    def reprFields(self):
        return self.rawFields()

    def write_bdf(self, size, card_writer):
        card = self.rawFields()
        return self.comment() + card_writer(card)


class PLOAD1(Load):
    type = 'PLOAD1'
    validTypes = ['FX', 'FY', 'FZ', 'FXE', 'FYE', 'FZE',
                  'MX', 'MY', 'MZ', 'MXE', 'MYE', 'MZE']
    validScales = ['LE', 'FR', 'LEPR', 'FRPR'] # LE: length-based; FR: fractional; PR:projected

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            self.sid = integer(card, 1, 'sid')
            self.eid = integer(card, 2, 'eid')
            self.Type = string(card, 3,  'Type ("%s")' % '",  "'.join(self.validTypes) )
            self.scale = string(card, 4, 'scale ("%s")' % '", "'.join(self.validScales) )
            self.x1 = double(card, 5, 'x1')
            self.p1 = double(card, 6, 'p1')
            self.x2 = double_or_blank(card, 7, 'x2', self.x1)
            self.p2 = double_or_blank(card, 8, 'p2', self.p1)
            assert len(card) <= 9, 'len(PLOAD1 card) = %i' % len(card)
        else:
            self.sid = data[0]
            self.eid = data[1]
            self.Type = data[2]
            self.scale = data[3]
            self.x1 = data[4]
            self.p1 = data[5]
            self.x2 = data[6]
            self.p2 = data[7]
        if self.Type not in self.validTypes:
            msg = '%s is an invalid type on the PLOAD1 card' % self.Type
            raise RuntimeError(msg)

        assert 0.0 <= self.x1 <= self.x2
        if self.scale in ['FR', 'FRPR']:
            assert self.x1 <= 1.0, 'x1=%r' % self.x1
            assert self.x2 <= 1.0, 'x2=%r' % self.x2
        assert self.scale in self.validScales, '%s is an invalid scale on the PLOAD1 card' % (self.scale)


    def cross_reference(self, model):
        """
        .. todo:: cross reference and fix repr function
        """
        self.eid = model.elements[self.eid]

    def transformLoad(self):
        p1 = self.eid.ga.Position()
        p2 = self.eid.gb.Position()

        g0 = self.eid.g0_vector
        #if not isinstance(g0, ndarray):
        #    g0 = g0.Position()

        x = p2 - p1
        y = p1 - g0
        z = cross(x, y)
        A = [x, y, z]
        #print("x =", x)
        #print("y =", y)
        #print("z =", z)
        #g = self.GravityVector()
        return A
        #(g2, matrix) = self.cid.transformToGlobal(A)
        #return (g2)

    def getReducedLoads(self):
        """
        Get all load objects in a simplified form, which means all
        scale factors are already applied and only base objects
        (no LOAD cards) will be returned.

        .. todo:: lots more object types to support
        """
        scale_factors = [1.0]
        loads = [self]
        return scale_factors, loads

    def organizeLoads(self, model):
        """
        Figures out magnitudes of the loads to be applied to the various nodes.
        This includes figuring out scale factors.
        """
        forceLoads = {}  # spc enforced displacement (e.g. FORCE=0)
        momentLoads = {}
        forceConstraints = {}
        momentConstraints = {}
        gravityLoads = []
        #print("self.loadIDs = ",self.loadIDs)

        typesFound = set()
        (scaleFactors, loads) = self.getReducedLoads()

        for (scaleFactor, load) in izip(scaleFactors, loads):
            #print("*load = ",load)
            out = load.transformLoad()
            typesFound.add(load.__class__.__name__)

            if isinstance(load, PLOAD1): # CBAR/CBEAM
                element = load.eid
                (ga, gb) = element.nodeIDs()
                load_type = load.Type

                scale = load.scale
                eType = element.type

                if load_type in ['FX', 'FY', 'FZ']:
                    p1 = element.ga.Position()
                    p2 = element.gb.Position()
                    r = p2 - p1

                if load_type == 'FX':
                    Fv = array([1., 0., 0.])
                elif load_type == 'FY':
                    Fv = array([0., 0., 1.])
                elif load_type == 'FZ':
                    Fv = array([0., 0., 1.])

                elif load_type == 'MX':
                    Fv = array([0., 0., 0.])
                    Mv = array([1., 0., 0.])
                elif load_type == 'MY':
                    Fv = array([0., 0., 0.])
                    Mv = array([0., 1., 0.])
                elif load_type == 'MZ':
                    Fv = array([0., 0., 0.])
                    Mv = array([0., 0., 1.])
                # FXE, FYE, FZE, MXE, MYE, MZE
                else:
                    raise NotImplementedError(load_type)

                p1 = load.p1
                p2 = load.p2
                if scale == 'FR':
                    x1 = load.x1
                    x2 = load.x2
                elif scale == 'LE':
                    L = element.Length()
                    x1 = load.x1 / L
                    x2 = load.x2 / L
                else:
                    raise NotImplementedError('scale=%s is not supported.  Use "FR", "LE".')

                assert x1 <= x2, '---load---\n%sx1=%r must be less than x2=%r' % (repr(self), self.x1, self.x2)
                if  x1 == x2:
                    msg = 'Point loads are not supported on...\n%sTry setting x1=%r very close to x2=%r and\n' % (repr(self), self.x1, self.x2)
                    msg += 'scaling p1=%r and p2=%r by x2-x1 (for "FR") and (x2-x1)/L (for "LE").' % (self.p1, self.p2)
                    raise NotImplementedError(msg)
                    if p1 != p2:
                        msg = 'p1=%r must be equal to p2=%r for x1=x2=%r'  %(self.p1, self.p2, self.x1)
                        raise RuntimeError(msg)

                dx = x2 - x1
                m = (p2 - p1) / dx
                #dx * (x2 + x1) = (x2^2-x1^2)
                dx2 = x2**2 - x1**2
                dx3 = x2**3 - x1**3
                dx4 = x2**4 - x1**4
                #F = (p1 - m * x1) * dx + m * dx2 / 2.
                F = p1 * dx + m * (dx2 / 2. - x1 * dx)

                F /= 2.

                if eType in ['CBAR', 'CBEAM']:
                    if load_type in ['FX', 'FY', 'FZ']:
                        M = p1 * dx3 / 6. + m * (dx4 / 12. - x1 * dx3 / 6.)
                        Mv = M * cross(r, Fv) / 2. # divide by 2 for 2 nodes

                        Fv *= F
                        #Mv = M
                        forceLoads[ga] = Fv
                        forceLoads[gb] = Fv
                        momentLoads[ga] = Mv
                        momentLoads[gb] = Mv
                    elif load_type in ['MX', 'MY', 'MZ']:
                        # these are really moments
                        Mv *= F
                        momentLoads[ga] = Mv
                        momentLoads[gb] = Mv
                    else:
                        raise NotImplementedError(load_type)
                else:
                    raise NotImplementedError(eType)

            else:
                msg = '%s not supported' % (load.__class__.__name__)
                raise NotImplementedError(msg)
        return (typesFound, forceLoads, momentLoads, forceConstraints,
                momentConstraints, gravityLoads)

    def getLoads(self):
        return [self]

    def Eid(self):
        if isinstance(self.eid, int):
            return self.eid
        return self.eid.eid

    def rawFields(self):
        list_fields = ['PLOAD1', self.sid, self.Eid(), self.Type, self.scale,
                  self.x1, self.p1, self.x2, self.p2]
        return list_fields

    def reprFields(self):
        return self.rawFields()

    def write_bdf(self, size, card_writer):
        card = self.rawFields()
        return self.comment() + card_writer(card)


class PLOAD2(Load):
    type = 'PLOAD2'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            self.sid = integer(card, 1, 'sid')
            self.p = double(card, 2, 'p')

            if integer_string_or_blank(card, 4, 'THRU') == 'THRU':
                e1 = integer(card, 3, 'Element1')
                e2 = integer(card, 5, 'Element1')
                eids = [i for i in xrange(e1, e2 + 1)]
                assert len(card) == 6, 'len(PLOAD2 card) = %i' % len(card)
            else:
                eids = fields(integer, card, 'eid', i=3, j=len(card))
            self.eids = eids
        else:
            self.sid = data[0]
            self.p = data[1]
            self.eids = list(data[2:])
            #print "PLOAD2 = ",data

    def cross_reference(self, model):
        """
        .. todo:: cross reference and fix repr function
        """
        pass

    def getLoads(self):
        return [self]

    def rawFields(self):
        list_fields = ['PLOAD2', self.sid, self.p]
        if len(self.eids) > 6:
            list_fields += [self.eids[0], 'THRU', self.eids[-1]]
        else:
            list_fields += self.eids
        return list_fields

    def reprFields(self):
        return self.rawFields()

    def write_bdf(self, size, card_writer):
        card = self.rawFields()
        return self.comment() + card_writer(card)


class PLOAD4(Load):
    type = 'PLOAD4'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            self.sid = integer(card, 1, 'sid')
            self.eid = integer(card, 2, 'eid')
            p1 = double_or_blank(card, 3, 'p1', 0.0)
            p = [p1,
                 double_or_blank(card, 4, 'p2', p1),
                 double_or_blank(card, 5, 'p3', p1),
                 double_or_blank(card, 6, 'p4', p1)]
            self.pressures = p

            self.eids = [self.eid]
            if (integer_string_or_blank(card, 7, 'g1/THRU') == 'THRU' and
                integer_or_blank(card, 8, 'eid2')):  # plates
                eid2 = integer(card, 8, 'eid2')
                if eid2:
                    self.eids = expand_thru([self.eid, 'THRU', eid2],
                                            set_fields=False, sort_fields=False)
                self.g1 = None
                self.g34 = None
            else:
                #: used for CPENTA, CHEXA
                self.eids = [self.eid]
                #: used for solid element only
                self.g1 = integer_or_blank(card, 7, 'g1')
                #: g3/g4 - different depending on CHEXA/CPENTA or CTETRA
                self.g34 = integer_or_blank(card, 8, 'g34')

            #: Coordinate system identification number. See Remark 2.
            #: (Integer >= 0;Default=0)
            self.cid = integer_or_blank(card, 9, 'cid', 0)
            #print "PLOAD4 cid = ",self.cid
            self.NVector = array([double_or_blank(card, 10, 'N1', 0.0),
                                  double_or_blank(card, 11, 'N2', 0.0),
                                  double_or_blank(card, 12, 'N3', 0.0)])
            self.sorl = string_or_blank(card, 13, 'sorl', 'SURF')
            self.ldir = string_or_blank(card, 14, 'ldir', 'NORM')
            assert len(card) <= 15, 'len(PLOAD4 card) = %i' % len(card)
        else:
            #print "PLOAD4 = ",data
            self.sid = data[0]
            self.eid = data[1]
            self.pressures = data[2]

            self.g1 = data[3]
            self.g34 = data[4]
            self.cid = data[5]
            self.NVector = data[6]

            self.sorl = data[7]
            #self.ldir = data[8]
            #assert len(data)==8

            self.g1 = self.g1
            self.g34 = self.g34
            self.eids = [self.eid]

    def getLoads(self):
        return [self]

    def transformLoad(self):
        """
        .. warning:: sorl=SURF is supported (not LINE)
        .. warning:: ldir=NORM is supported (not X,Y,Z)
        """
        if  self.sorl != 'SURF':
            msg = 'Only surface loads are supported.  required_sorl=SURF.  actual=%s' % self.sorl
            raise RuntimeError(msg)
        if self.ldir != 'NORM':
            msg = 'Only normal loads are supported.   required_ldir=NORM.  actual=%s' % self.ldir
            raise RuntimeError(msg)
        if len(self.eids) != 1:
            msg = 'Only one load may be defined on each PLOAD4.  nLoads=%s\n%s' % (len(self.eids), str(self))
            raise RuntimeError(msg)

        if self.g1 and self.g34:  # solid elements
            nid = self.g1.nid
            nidOpposite = self.g34.nid
            (faceNodeIDs, Area) = self.eid.getFaceNodesAndArea(self, nid,
                                                               nidOpposite)
        else:
            faceNodeIDs = self.eid.nodeIDs()
            Area = self.eid.Area()
        n = len(faceNodeIDs)

        vector = array(self.eid.Normal())
        vectors = []
        for (nid, p) in izip(faceNodeIDs, self.pressures):
            #: .. warning:: only supports normal pressures
            vectors.append(vector * p * Area / n)  # Force_i

        isLoad = None
        return (isLoad, faceNodeIDs, vectors)

    def Cid(self):
        if isinstance(self.cid, int):
            return self.cid
        return self.cid.cid

    def cross_reference(self, model):
        self.eid = model.Element(self.eid)
        self.cid = model.Coord(self.cid)
        if self.g1:
            self.g1 = model.Node(self.g1)
        if self.g34:
            self.g34 = model.Node(self.g34)
        if self.eids:
            self.eids = model.Elements(self.eids)

    def Eid(self, eid):
        if isinstance(eid, int):
            return eid
        return eid.eid

    def getElementIDs(self, eid=None):
        if eid:
            return self.Eid(eid)
        eids = []
        for element in self.eids:
            eids.append(self.Eid(element))
        return eids

    def rawFields(self):
        eid = self.Eid(self.eid)
        cid = set_blank_if_default(self.Cid(), 0)
        sorl = set_blank_if_default(self.sorl, 'SURF')
        ldir = set_blank_if_default(self.ldir, 'NORM')
        p1 = self.pressures[0]
        p2 = set_blank_if_default(self.pressures[1], p1)
        p3 = set_blank_if_default(self.pressures[2], p1)
        p4 = set_blank_if_default(self.pressures[3], p1)
        list_fields = ['PLOAD4', self.sid, eid, self.pressures[0], p2, p3, p4]

        #print "g3=|%s| g4=%s eids=|%s|" %(self.g3,self.g4,self.eids)
        if self.g1 is not None:  # is it a SOLID element
            (g1, g34) = self.nodeIDs([self.g1, self.g34])
            list_fields.append(g1)
            list_fields.append(g34)
        else:
            #print "eids = %s" %(self.eids)
            if len(self.eids) > 1:
                #print("self.eids = %s" %(self.eids))
                try:
                    list_fields.append('THRU')
                    eid = self.eids[-1]
                except:
                    print("g1  = %s" % self.g1)
                    print("g34 = %s" % self.g34)
                    print("self.eids = %s" % self.eids)
                    raise
                list_fields.append(self.getElementIDs(eid))
            else:
                list_fields += [None, None]
        list_fields.append(cid)

        n1 = set_blank_if_default(self.NVector[0], 0.0)
        n2 = set_blank_if_default(self.NVector[1], 0.0)
        n3 = set_blank_if_default(self.NVector[2], 0.0)
        list_fields += [n1, n2, n3]
        list_fields.append(sorl)
        list_fields.append(ldir)
        return list_fields

    def reprFields(self):
        return self.rawFields()

    def write_bdf(self, size, card_writer):
        card = self.rawFields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)
        #return self.comment() + card_writer(card)


class PLOADX1(Load):
    type = 'PLOADX1'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            self.sid = integer(card, 1, 'sid')
            self.eid = integer(card, 2, 'eid')
            self.pa = double(card, 3, 'pa')
            self.pb = double_or_blank(card, 4, 'pb', self.pa)
            self.ga = integer(card, 5, 'ga')
            self.gb = integer(card, 6, 'gb')
            self.theta = double_or_blank(card, 7, 'theta', 0.)
            assert len(card) <= 8, 'len(PLOADX1 card) = %i' % len(card)
        else:
            self.sid = data[0]
            print("PLOADX1 = %s" % data)
            raise NotImplementedError(data)

    def cross_reference(self, model):
        #self.eid = model.Element(self.eid)
        #self.ga = model.Node(self.ga)
        #self.gb = model.Node(self.gb)
        pass

    def getLoads(self):
        return [self]

    def rawFields(self):
        list_fields = ['PLOADX1', self.sid, self.eid, self.pa, self.pb,
                  self.ga, self.gb, self.theta]
        return list_fields

    def reprFields(self):
        return self.rawFields()

    def write_bdf(self, size, card_writer):
        card = self.rawFields()
        return self.comment() + card_writer(card)