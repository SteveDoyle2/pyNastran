# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from itertools import izip

from numpy import array, cross, allclose
from numpy.linalg import norm

from pyNastran.bdf.cards.loads.loads import Load, LoadCombination
from pyNastran.bdf.fieldWriter import set_blank_if_default
from ..baseCard import BaseCard, expand_thru, expand_thru_by


class LOAD(LoadCombination):
    type = 'LOAD'

    def __init__(self, card=None, data=None):
        LoadCombination.__init__(self, card, data)

    def getLoadIDs(self):
        """
        @note requires a cross referenced load
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

                ###
        ###
        load_IDs = list(set(load_IDs))
        #print "load_IDs = ",load_IDs
        return load_IDs

    def getLoadTypes(self):
        """
        @note requires a cross referenced load
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
                    raise RuntimeError(load)
                ###
        ###
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
        Get all load objects in a simplified form,
        which means all scale factors are already applied and
        only base objects (no LOAD cards) will be returned.
        @todo lots more object types to support
        """
        scaleFactors = []
        loads = []
        scale = self.scale
        for (loadsPack, scaleFactorI) in izip(self.loadIDs, self.scaleFactors):
            scale2 = scaleFactorI * scale
            for load in loadsPack:
                if (isinstance(load, Force) or isinstance(load, Moment) or
                    isinstance(load, PLOAD4) or isinstance(load, GRAV)):
                    loads.append(load)
                    scaleFactors.append(scale2)
                elif isinstance(load, LOAD):
                    (scaleFactorsi, loadsi) = load.getReducedLoads()
                    loads += loadsi
                    scaleFactors += [scale2 * scalei for scalei in scaleFactorsi]
                else:
                    msg = ('%s isnt supported in getReducedLoads method'
                        % (load.__class__.__name__))
                    raise NotImplementedError(msg)
                ###
            ###
        ###
        return (scaleFactors, loads)

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
                    ###
                else:  # constraint
                    if node not in forceLoads:
                        forceConstraints[node] = vector * scaleFactor
                    else:
                        forceConstraints[node] += vector * scaleFactor
                    ###
                ###
            elif isinstance(load, Moment):
                (isLoad, node, vector) = out
                if isLoad:  # load
                    if node not in momentLoads:
                        momentLoads[node] = vector * scaleFactor
                    else:
                        momentLoads[node] += vector * scaleFactor
                    ###
                else:  # constraint
                    if node not in momentLoads:
                        momentConstraints[node] = vector * scaleFactor
                    else:
                        momentConstraints[node] += vector * scaleFactor
                    ###
                ###
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
            ###
        ###
        return (typesFound, forceLoads, momentLoads, forceConstraints,
                momentConstraints, gravityLoads)

    def rawFields(self):
        fields = ['LOAD', self.sid, self.scale]
        for (scaleFactor, loadID) in izip(self.scaleFactors, self.loadIDs):
            fields += [scaleFactor, self.LoadID(loadID)]
        return fields

    def reprFields(self):
        return self.rawFields()


class GRAV(BaseCard):
    """
    Defines acceleration vectors for gravity or other acceleration loading
    GRAV SID CID A     N1  N2 N3    MB
    GRAV 1   3   32.2 0.0 0.0 -1.0
    """
    type = 'GRAV'

    def __init__(self, card=None, data=None):
        if card:
            ## Set identification number
            self.sid = card.field(1)
            ## Coordinate system identification number.
            self.cid = card.field(2, 0)
            ## scale factor
            self.scale = card.field(3)
            ## Acceleration vector components measured in coordinate system CID
            self.N = array(card.fields(4, 7, [0., 0., 0.]))
            ## Indicates whether the CID coordinate system is defined in the
            ## main Bulk Data Section (MB = -1) or the partitioned superelement
            ## Bulk Data Section (MB = 0). Coordinate systems referenced in the
            ## main Bulk Data Section are considered stationary with respect to
            ## the assembly basic coordinate system. See Remark 10.
            ## (Integer; Default = 0)
            self.mb = card.field(7, 0)
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
        fields = ['GRAV', self.sid, self.Cid(), self.scale] + N + [self.mb]
        return fields

    def reprFields(self):
        N = []
        for n in self.N:
            N.append(set_blank_if_default(n, 0.0))

        mb = set_blank_if_default(self.mb, 0)
        fields = ['GRAV', self.sid, self.Cid(), self.scale] + N + [mb]
        return fields


class ACCEL1(BaseCard):
    """
    Acceleration Load
    Defines static acceleration loads at individual GRID points.
    """
    type = 'ACCEL1'

    def __init__(self, card=None, data=None):
        ## Load set identification number (Integer>0)
        self.sid = card.field(1)

        ## Coordinate system identification number. (Integer>0: Default=0)
        self.cid = card.field(2, 0)

        ## Acceleration vector scale factor. (Real)
        self.scale = card.field(3)

        ## Components of the acceleration vector measured in coordinate system
        ## CID. (Real; at least one Ni != 0)
        self.N = array(card.fields(4, 7, [0., 0., 0.]))
        assert max(abs(self.N)) > 0.
        ## nodes to apply the acceleration to
        self.nodes = expand_thru_by(card.fields(9))

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
        ###
        assert 0 not in nodeIDs, 'nodeIDs = %s' % (nodeIDs)
        return nodeIDs

    def rawFields(self):
        fields = ['ACCEL1', self.sid, self.Cid(), self.scale,
                  self.N[0], self.N[1], self.N[2], None, None] + self.nodeIDs()
        return fields


class OneDeeLoad(Load):  # FORCE/MOMENT
    type = '1D_Load'

    def __init__(self, card, data):
        self.cid = None
        self.mag = None
        self.xyz = None
        Load.__init__(self, card, data)

    def getLoads(self):
        return [self]

    def transformLoad(self):
        #print("self.xyz = ",self.xyz)
        (xyz, matrix) = self.cid.transformToGlobal(self.xyz)
        if self.mag > 0.:
            #print("mag=%s xyz=%s" % (self.mag, xyz))
            return (True, self.node, self.mag * xyz)  # load
        return (False, self.node, xyz)  # enforced displacement

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


class Force(OneDeeLoad):
    """Generic class for all Forces"""
    type = '1D_Load'

    def __init__(self, card, data):
        OneDeeLoad.__init__(self, card, data)

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


class Moment(OneDeeLoad):
    """Generic class for all Moments"""
    type = 'Moment'

    def __init__(self, card, data):
        OneDeeLoad.__init__(self, card, data)

    def getLoads(self):
        return [self]

    def getReducedLoads(self):
        scaleFactors = [1.]
        loads = self.F()
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

#------------------------------------------------------------------------------


class FORCE(Force):
    type = 'FORCE'

    def __init__(self, card=None, data=None):
        """
        FORCE          3       1            100.      0.      0.      1.
        """
        Force.__init__(self, card, data)
        if card:
            self.sid = card.field(1)
            self.node = card.field(2)
            self.cid = card.field(3, 0)
            self.mag = card.field(4)
            xyz = card.fields(5, 8, [0., 0., 0.])
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
        """@todo cross reference and fix repr function"""
        self.cid = model.Coord(self.cid)

    def rawFields(self):
        fields = ['FORCE', self.sid, self.node, self.Cid(), self.mag
                 ] + list(self.xyz)
        return fields

    def reprFields(self):
        cid = set_blank_if_default(self.Cid(), 0)
        fields = ['FORCE', self.sid, self.node, cid, self.mag] + list(self.xyz)
        return fields


class FORCE1(Force):
    """
    Defines a static concentrated force at a grid point by specification of a
    magnitude and two grid points that determine the direction.
    """
    type = 'FORCE1'

    def __init__(self, card=None, data=None):
        Force.__init__(self, card, data)
        if card:
            self.sid = card.field(1)
            self.node = card.field(2)
            self.mag = card.field(3)
            self.g1 = card.field(4)
            self.g2 = card.field(5)
        else:
            self.sid = data[0]
            self.node = data[1]
            self.mag = data[2]
            self.g1 = data[3]
            self.g2 = data[4]
        ###

    def cross_reference(self, model):
        """@todo cross reference and fix repr function"""
        self.node = model.Node(self.node)
        self.g1 = model.Node(self.g1)
        self.g2 = model.Node(self.g2)
        self.xyz = self.g2.Position() - self.g1.Position()
        self.Normalize()

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
        (node, g1, g2) = self.nodeIDs([self.node, self.G1(), self.G2()])
        fields = ['FORCE1', self.sid, node, self.mag, g1, g2]
        return fields

    def reprFields(self):
        return self.rawFields()


class FORCE2(Force):
    """
    Defines a static concentrated force at a grid point by specification of a
    magnitude and four grid points that determine the direction.
    """
    type = 'FORCE2'

    def __init__(self, card=None, data=None):
        """
        FORCE2 SID G F G1 G2 G3 G4
        """
        Force.__init__(self, card, data)
        if card:
            self.sid = card.field(1)
            self.node = card.field(2)
            self.mag = card.field(3)
            self.g1 = card.field(4)
            self.g2 = card.field(5)
            self.g3 = card.field(5)
            self.g4 = card.field(5)
        else:
            self.sid = data[0]
            self.node = data[1]
            self.mag = data[2]
            self.g1 = data[3]
            self.g2 = data[4]
            self.g3 = data[5]
            self.g4 = data[6]
        ###

    def cross_reference(self, model):
        """@todo cross reference and fix repr function"""
        self.node = model.Node(self.node)

        v12 = model.Node(self.g2).Position() - model.Node(self.g1).Position()
        v34 = model.Node(self.g4).Position() - model.Node(self.g3).Position()
        v12 = v12 / norm(v12)
        v34 = v34 / norm(v34)
        self.xyz = cross(v12, v34)
        self.Normalize()

    def NodeID(self):
        if isinstance(self.node, int):
            return self.node
        return self.node.nid

    def rawFields(self):
        (node, g1, g2, g3, g4) = self.nodeIDs([self.node, self.g1, self.g2,
                                               self.g3, self.g4])
        fields = ['FORCE2', self.sid, node, self.mag, g1, g2, g3, g4]
        return fields

    def reprFields(self):
        return self.rawFields()


class MOMENT(Moment):
    type = 'MOMENT'

    def __init__(self, card=None, data=None):
        """
        Defines a static concentrated moment at a grid point by specifying a
        scale factor and a vector that determines the direction.

        MOMENT SID G CID M    N1  N2  N3
        MOMENT 2   5   6 2.9 0.0 1.0 0.0
        """
        Moment.__init__(self, card, data)
        self.sid = card.field(1)
        self.node = card.field(2)
        self.cid = card.field(3, 0)
        self.mag = card.field(4)

        xyz = card.fields(5, 8, [0., 0., 0.])
        assert len(xyz) == 3, 'xyz=%s' % (xyz)
        self.xyz = array(xyz)

    def Cid(self):
        if isinstance(self.cid, int):
            return self.cid
        return self.cid.cid

    def cross_reference(self, model):
        """@todo cross reference and fix repr function"""
        pass

    def rawFields(self):
        fields = ['MOMENT', self.sid, self.node, self.Cid(), self.mag
                 ] + list(self.xyz)
        return fields

    def reprFields(self):
        cid = set_blank_if_default(self.Cid(), 0)
        fields = ['MOMENT', self.sid, self.node, cid, self.mag
                 ] + list(self.xyz)
        return fields


class MOMENT1(Moment):
    type = 'MOMENT1'

    def __init__(self, card=None, data=None):
        """
        Defines a static concentrated moment at a grid point by specifying a
        magnitude and two grid points that determine the direction

        MOMENT1 SID G M G1 G2
        """
        Moment.__init__(self, card, data)
        if card:
            self.sid = card.field(1)
            self.node = card.field(2)
            self.mag = card.field(3)
            self.g1 = card.field(4)
            self.g2 = card.field(5)
            self.g3 = card.field(6)
            self.g4 = card.field(7)
            xyz = card.fields(5, 8, [0., 0., 0.])
        else:
            self.sid = data[0]
            self.node = data[1]
            self.mag = data[2]
            self.g1 = data[3]
            self.g2 = data[4]
            self.g3 = data[5]
            self.g4 = data[6]
            xyz = data[7:10]

        assert len(xyz) == 3, 'xyz=%s' % (xyz)
        self.xyz = array(xyz)

    def cross_reference(self, model):
        """@todo cross reference and fix repr function"""
        self.node = model.Node(self.node)
        self.xyz = model.Node(
            self.g2).Position() - model.Node(self.g1).Position()
        self.Normalize()

    def rawFields(self):
        (node, g1, g2) = self.nodeIDs([self.node, self.g1, self.g2])
        fields = ['MOMENT1', self.sid, node, self.mag, g1, g2]
        return fields

    def reprFields(self):
        return self.rawFields()


class MOMENT2(Moment):
    type = 'MOMENT2'

    def __init__(self, card=None, data=None):
        """
        Defines a static concentrated moment at a grid point by specification
        of a magnitude and four grid points that determine the direction.

        MOMENT2 SID G M G1 G2 G3 G4
        """
        Moment.__init__(self, card, data)
        if card:
            self.sid = card.field(1)
            self.node = card.field(2)
            self.mag = card.field(3)
            self.g1 = card.field(4)
            self.g2 = card.field(5)
            self.g3 = card.field(6)
            self.g4 = card.field(7)
            xyz = card.fields(5, 8, [0., 0., 0.])
        else:
            self.sid = data[0]
            self.node = data[1]
            self.mag = data[2]
            self.g1 = data[3]
            self.g2 = data[4]
            self.g3 = data[5]
            self.g4 = data[6]
            xyz = data[7:10]
        ###
        assert len(xyz) == 3, 'xyz=%s' % (xyz)
        self.xyz = array(xyz)

    def cross_reference(self, model):
        """@todo cross reference and fix repr function"""
        (self.g1, self.g2, self.g3, self.g4) = model.Nodes(self.g1, self.g2,
                                                           self.g3, self.g4)
        v12 = self.g2.Position() - self.g1.Position()
        v34 = self.g4.Position() - self.g3.Position()
        v12 = v12 / norm(v12)
        v34 = v34 / norm(v34)
        self.xyz = cross(v12, v34)

    def rawFields(self):
        (node, g1, g2, g3, g4) = self.nodeIDs([self.node, self.g1, self.g2,
                                                          self.g3, self.g4])
        fields = ['MOMENT2', self.sid, node, self.mag, g1, g2, g3, g4]
        return fields

    def reprFields(self):
        return self.rawFields()


class PLOAD(Load):
    type = 'PLOAD'

    def __init__(self, card=None, data=None):
        if card:
            self.sid = card.field(1)
            self.p = card.field(2)
            nodes = card.fields(3, 7)
            self.nodes = self._wipeEmptyFields(nodes)
        else:
            self.sid = data[0]
            self.p = data[1]
            self.nodes = data[2:]
            print("PLOAD = %s" % (data))
            raise NotImplementedError('PLOAD')
        assert len(self.nodes) in [3, 4], 'nodes=%s' % (self.nodes)

    def cross_reference(self, model):
        """@todo cross reference and fix repr function"""
        pass

    def getLoads(self):
        return [self]

    def rawFields(self):
        fields = ['PLOAD', self.sid, self.p] + self.nodeIDs()
        return fields

    def reprFields(self):
        return self.rawFields()


class PLOAD1(Load):
    type = 'PLOAD1'
    validTypes = ['FX', 'FY', 'FZ', 'FXE', 'FYE', 'FZE',
                  'MX', 'MY', 'MZ', 'MXE', 'MYE', 'MZE']
    validScales = ['LE', 'FR', 'LEPR', 'FRPR']

    def __init__(self, card=None, data=None):
        if card:
            self.sid = card.field(1)
            self.eid = card.field(2)
            self.Type = card.field(3)
            self.scale = card.field(4)
            self.x1 = card.field(5)
            self.p1 = card.field(6)
            self.x2 = card.field(7)
            self.p2 = card.field(8)
        else:
            self.sid = data[0]
            self.eid = data[1]
            self.Type = data[2]
            self.scale = data[3]
            self.x1 = data[4]
            self.p1 = data[5]
            self.x2 = data[6]
            self.p2 = data[7]
        ###
        assert self.Type in self.validTypes, '%s is an invalid type on the PLOAD1 card' % (self.Type)
        assert self.scale in self.validScales, '%s is an invalid scale on the PLOAD1 card' % (self.scale)

    def cross_reference(self, model):
        """@todo cross reference and fix repr function"""
        pass

    def getLoads(self):
        return [self]

    def rawFields(self):
        fields = ['PLOAD1', self.sid, self.eid, self.Type, self.scale, self.x1,
                  self.p1, self.x2, self.p2]
        return fields

    def reprFields(self):
        return self.rawFields()


class PLOAD2(Load):
    type = 'PLOAD2'

    def __init__(self, card=None, data=None):
        if card:
            self.sid = card.field(1)
            self.p = card.field(2)
            eids = card.fields(3, 9)

            if card.field(4) == 'THRU':
                eids = [i for i in xrange(eids[0], eids[2] + 1)]
            self.eids = eids
        else:
            self.sid = data[0]
            self.p = data[1]
            self.eids = list(data[2:])
            #print "PLOAD2 = ",data

    def cross_reference(self, model):
        """@todo cross reference and fix repr function"""
        pass

    def getLoads(self):
        return [self]

    def rawFields(self):
        fields = ['PLOAD2', self.sid, self.p]
        if len(self.eids) > 6:
            fields += [self.eids[0], 'THRU', self.eids[-1]]
        else:
            fields += self.eids
        return fields

    def reprFields(self):
        return self.rawFields()


class PLOAD4(Load):
    type = 'PLOAD4'

    def __init__(self, card=None, data=None):
        if card:
            self.sid = card.field(1)
            self.eid = card.field(2)
            p1 = card.field(3)
            p = card.fields(4, 7, [p1, p1, p1])  # [p1,p1,p1] are the defaults
            self.pressures = [p1] + p

            self.eids = [self.eid]
            if card.field(7) == 'THRU' and card.field(8):  # plates
                eid2 = card.field(8)
                if eid2:
                    self.eids = expand_thru([self.eid, 'THRU', eid2])

                self.g1 = None
                self.g34 = None
            else:
                ## used for CPENTA, CHEXA
                self.eids = [self.eid]
                ## used for solid element only
                self.g1 = card.field(7)
                ## g3/g4 - different depending on CHEXA/CPENTA or CTETRA
                self.g34 = card.field(8)

            ## Coordinate system identification number. See Remark 2.
            ## (Integer >= 0;Default=0)
            self.cid = card.field(9, 0)
            #print "PLOAD4 cid = ",self.cid
            self.NVector = card.fields(10, 13, [0., 0., 0.])
            self.sorl = card.field(13, 'SURF')
            self.ldir = card.field(14, 'NORM')
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
            #self.ldir    = data[8]
            #assert len(data)==8

            self.g1 = self.g1
            self.g34 = self.g34
            self.eids = [self.eid]
        ###

    def getLoads(self):
        return [self]

    def transformLoad(self):
        """
        @warning sorl=SURF is supported (not LINE)
        @warning ldir=NORM is supported (not X,Y,Z)
        """
        assert self.sorl == 'SURF', 'only surface loads are supported.  required_sorl=SURF.  actual=%s' % (self.sorl)
        assert self.ldir == 'NORM', 'only normal loads are supported.   required_ldir=NORM.  actual=%s' % (self.ldir)
        assert len(self.eids) == 1, 'only one load may be defined on each PLOAD4.  nLoads=%s\n%s' % (len(self.eids), str(self))

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
            ## @warning only supports normal pressures
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
        fields = ['PLOAD4', self.sid, eid, self.pressures[0], p2, p3, p4]

        #print "g3=|%s| g4=%s eids=|%s|" %(self.g3,self.g4,self.eids)
        if self.g1 is not None:  # is it a SOLID element
            (g1, g34) = self.nodeIDs([self.g1, self.g34])
            fields.append(g1)
            fields.append(g34)
        else:
            #print "eids = %s" %(self.eids)
            if len(self.eids) > 1:
                #print("self.eids = %s" %(self.eids))
                try:
                    fields.append('THRU')
                    eid = self.eids[-1]
                except:
                    print("g1  = %s" % (self.g1))
                    print("g34 = %s" % (self.g34))
                    print("self.eids = %s" % (self.eids))
                    raise
                fields.append(self.getElementIDs(eid))
            else:
                fields += [None, None]
        fields.append(cid)

        n1 = set_blank_if_default(self.NVector[0], 0.0)
        n2 = set_blank_if_default(self.NVector[1], 0.0)
        n3 = set_blank_if_default(self.NVector[2], 0.0)
        fields += [n1, n2, n3]
        fields.append(sorl)
        fields.append(ldir)
        return fields

    def reprFields(self):
        return self.rawFields()


class PLOADX1(Load):
    type = 'PLOADX1'

    def __init__(self, card=None, data=None):
        if card:
            self.sid = card.field(1)
            self.eid = card.field(2)
            self.pa = card.field(3)  # float
            self.pb = card.field(4, self.pa)  # float
            self.ga = card.field(5)  # int
            self.gb = card.field(6)  # int
            self.theta = card.field(7, 0.)
        else:
            self.sid = data[0]
            print("PLOADX1 = %s" % (data))
            raise NotImplementedError('PLOADX1')

    def cross_reference(self, model):
        #self.eid = model.Element(self.eid)
        #self.ga = model.Node(self.ga)
        #self.gb = model.Node(self.gb)
        pass

    def getLoads(self):
        return [self]

    def rawFields(self):
        fields = ['PLOADX1', self.sid, self.eid, self.pa, self.pb,
                  self.ga, self.gb, self.theta]
        return fields

    def reprFields(self):
        return self.rawFields()
