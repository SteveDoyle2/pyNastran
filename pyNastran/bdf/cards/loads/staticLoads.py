# pylint: disable=C0103,R0902,R0904,R0914
"""
All static loads are defined in this file.  This includes:

 * LOAD
 * GRAV
 * ACCEL
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
from six import integer_types
from six.moves import zip

from numpy import array, cross, allclose, unique, int32
from numpy.linalg import norm

from pyNastran.bdf.errors import CrossReferenceError
from pyNastran.bdf.cards.loads.loads import Load, LoadCombination
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.baseCard import BaseCard, expand_thru, expand_thru_by, range
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, string, string_or_blank,
    integer_or_string, fields, integer_string_or_blank, integer_or_double)
from pyNastran.bdf.field_writer_8 import print_card_8, print_float_8, set_string8_blank_if_default
from pyNastran.bdf.field_writer_16 import print_card_16, print_float_16, set_string16_blank_if_default
from pyNastran.bdf.field_writer_double import print_card_double, print_scientific_double


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
                #if isinstance(load, int):
                    #load_IDs += [load]

                if isinstance(load, LOAD):
                    lid = load.lid
                    if isinstance(lid, list):
                        load_IDs += load.lid
                    else:  # int
                        load_IDs += load.getLoadIDs()
                elif (isinstance(load, Force) or isinstance(load, Moment) or
                      isinstance(load, PLOAD4) or isinstance(load, GRAV)):
                    load_IDs += [load.sid]
                else:
                    msg = ('The getLoadIDs method doesnt support %s cards.\n'
                           '%s' % (load.__class__.__name__, str(load)))
                    raise NotImplementedError(msg)

        load_IDs = list(set(load_IDs))
        #print("load_IDs = ", load_IDs)
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
        #print("loadTypes = ", loadTypes)
        return loadTypes

    def write_calculix_GRAV(self, gx, gy, gz):
        msg = '*DLOAD\n'
        msg += 'AllElements,GRAV,%s,%s,%s\n' % (gx, gy, gz)
        return msg

    def write_code_aster_load(self, model, grid_word='node'):
        loadIDs = self.getLoadIDs()
        loadTypes = self.getLoadTypes()

        #msg = '# Loads\n'
        msg = ''
        (typesFound, forceLoads, momentLoads,
         forceConstraints, momentConstraints,
         gravityLoads) = self.organize_loads(model)

        nids = []
        for nid in forceLoads:
            nids.append(nid)
        for nid in momentLoads:
            nids.append(nid)

        if nids:
            msg += '# typesFound = %s\n' % (list(typesFound))
            msg += '# loadIDs    = %s\n' % (loadIDs)
            msg += "load_bc = AFFE_CHAR_MECA(MODELE=modmod,\n"
            #msg += "                        DDL_IMPO=(_F(GROUP_MA='Lleft',\n"
            msg += "                         FORCE_NODALE=(\n"

        #CHAR=AFFE_CHAR_MECA(MODELE=MODE,
        #             FORCE_NODALE=(
        #                     _F(NOEUD='N1',
        #                        FZ=-500.0),)

        spaces = "                           "
        for nid in sorted(nids):  # ,load in sorted(iteritems(forceLoads))
            msg += spaces + "_F(NOEUD='%s%s'," % (grid_word, nid)

            if nid in forceLoads:
                force = forceLoads[nid]
                if abs(force[0]) > 0.:
                    msg += " FX=%s," % force[0]
                if abs(force[1]) > 0.:
                    msg += " FY=%s," % force[1]
                if abs(force[2]) > 0.:
                    msg += " FZ=%s," % force[2]

            if nid in momentLoads:
                moment = momentLoads[nid]
                if abs(moment[0]) > 0.:
                    msg += " MX=%s," % moment[0]
                if abs(moment[1]) > 0.:
                    msg += " MY=%s," % moment[1]
                if abs(moment[2]) > 0.:
                    msg += " MZ=%s," % moment[2]
            #msg = msg[:-2]
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

        for gravity_load in gravityLoads:
            msg += 'CA_GRAVITY(%s);\n' % str(gravity_load)
        return msg, loadIDs, loadTypes

    def getReducedLoads(self):
        self.deprecated('self.getReducedLoads()', 'self.get_reduced_loads()', '0.8')
        return self.get_reduced_loads()

    def get_reduced_loads(self):
        """
        Get all load objects in a simplified form, which means all
        scale factors are already applied and only base objects
        (no LOAD cards) will be returned.

        .. todo:: lots more object types to support
        """
        scale_factors = []
        loads = []
        load_scale = self.scale # global
        for (loads_pack, i_scale) in zip(self.loadIDs, self.scaleFactors):
            scale = i_scale * load_scale # actual scale = global * local
            if isinstance(loads_pack, integer_types):
                raise RuntimeError('the load have not been cross-referenced')

            for load in loads_pack:
                if ['FORCE', 'FORCE1', 'FORCE2',
                    'MOMENT', 'MOMENT1', 'MOMENT2',
                    'PLOAD1', 'PLOAD2', 'PLOAD4',
                    'GRAV', 'ACCEL', 'ACCEL1']:
                    loads.append(load)
                    scale_factors.append(scale) # local
                elif isinstance(load, LOAD):
                    load_data = load.get_reduced_loads()
                    (reduced_scale_factors, reduced_loads) = load_data

                    loads += reduced_loads
                    scale_factors += [scale * j_scale
                                      for j_scale in reduced_scale_factors]
                else:
                    msg = ('%s isnt supported in get_reduced_loads method'
                           % load.__class__.__name__)
                    raise NotImplementedError(msg)

        return (scale_factors, loads)

    def organizeLoads(self, model):
        self.deprecated('organizeLoads(model)', 'organize_loads(model)', '0.8')
        return self.organize_loads(model)

    def organize_loads(self, model):
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
        (scale_factors, loads) = self.get_reduced_loads()

        for (scale_factor, load) in zip(scale_factors, loads):
            #print("*load = ",load)
            out = load.transformLoad()
            typesFound.add(load.__class__.__name__)
            if isinstance(load, Force):
                (isLoad, node, vector) = out
                if isLoad:  # load
                    if node not in forceLoads:
                        forceLoads[node] = vector * scale_factor
                    else:
                        forceLoads[node] += vector * scale_factor
                else:  # constraint
                    if node not in forceLoads:
                        forceConstraints[node] = vector * scale_factor
                    else:
                        forceConstraints[node] += vector * scale_factor

            elif isinstance(load, Moment):
                (isLoad, node, vector) = out
                if isLoad:  # load
                    if node not in momentLoads:
                        momentLoads[node] = vector * scale_factor
                    else:
                        momentLoads[node] += vector * scale_factor
                else:  # constraint
                    if node not in momentLoads:
                        momentConstraints[node] = vector * scale_factor
                    else:
                        momentConstraints[node] += vector * scale_factor

            elif isinstance(load, PLOAD4):
                (isLoad, nodes, vectors) = out
                for (nid, vector) in zip(nodes, vectors):
                    # not the same vector for all nodes
                    forceLoads[nid] = vector * scale_factor

            elif isinstance(load, GRAV):
                #(grav) = out
                gravityLoads.append(out * scale_factor)  # grav
            else:
                msg = '%s not supported' % (load.__class__.__name__)
                raise NotImplementedError(msg)

        return (typesFound, forceLoads, momentLoads, forceConstraints,
                momentConstraints, gravityLoads)

    def raw_fields(self):
        list_fields = ['LOAD', self.sid, self.scale]
        for (scale_factor, load_id) in zip(self.scaleFactors, self.loadIDs):
            list_fields += [scale_factor, self.LoadID(load_id)]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        else:
            return self.comment + print_card_16(card)

    def uncross_reference(self):
        ids = []
        for i, load_id in enumerate(self.loadIDs):
            idi = self.LoadID(load_id)
            ids.append(idi)
        self.loadIDs = ids


class GRAV(BaseCard):
    """
    Defines acceleration vectors for gravity or other acceleration loading.::

    +------+-----+-----+------+-----+-----+------+-----+
    | GRAV | SID | CID | A    | N1  | N2  | N3   |  MB |
    +------+-----+-----+------+-----+-----+------+-----+
    | GRAV | 1   | 3   | 32.2 | 0.0 | 0.0 | -1.0 |     |
    +------+-----+-----+------+-----+-----+------+-----+
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
                                                   'N=%s' % str(self.N))

    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

    def get_loads(self):
        return [self]

    def organizeLoads(self, model):
        self.deprecated('organizeLoads(model)', 'organize_loads(model)', '0.8')
        return self.organize_loads(model)

    def organize_loads(self, model):
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
        g2, matrix = self.cid.transformToGlobal(g)
        return g2

    #write_code_aster_load(self,mag):
        #p = self.GravityVector()
        #msg = 'GRAV([%s,%s,%s])' %(p)
        #return msg

    def cross_reference(self, model):
        msg = ' which is required by GRAV sid=%s' % self.sid
        self.cid = model.Coord(self.cid, msg=msg)

    def safe_cross_reference(self, model):
        # msg = "Couldn't find CORDx=%s which is required by GRAV sid=%s" % (self.cid, self.sid)
        msg = ' which is required by GRAV sid=%s' % self.sid
        self.cid = model.Coord(self.cid, msg=msg)

    def uncross_reference(self, model):
        self.cid = self.Cid()

    def Cid(self):
        if isinstance(self.cid, integer_types):
            return self.cid
        return self.cid.cid

    def GravityVector(self):
        """returns the gravity vector in absolute coordinates"""
        ## TODO: shouldn't be scaled by the ???
        p = self.cid.transform_vector_to_global(self.N)
        return self.scale * p

    def raw_fields(self):
        N = list(self.N)
        list_fields = ['GRAV', self.sid, self.Cid(), self.scale] + N + [
                       self.mb]
        return list_fields

    def repr_fields(self):
        N = []
        for n in self.N:
            N.append(set_blank_if_default(n, 0.0))

        mb = set_blank_if_default(self.mb, 0)
        list_fields = ['GRAV', self.sid, self.Cid(), self.scale] + N + [mb]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class ACCEL(BaseCard):
    """
    Acceleration Load

    Defines static acceleration loads, which may vary over a region of
    the structural model. The load variation is based upon the tabular
    input defined on this Bulk Data entry.

    +-------+------+------+------+------+-----+-----+--------+-----+
    |   1   |   2  |   3  |  4   |   5  |  6  |  7  | 8      |  9  |
    +-------+------+------+------+------+-----+-----+--------+-----+
    | ACCEL | SID  | CID  | N1   | N2   | N3  | DIR |        |     |
    +-------+------+------+------+------+-----+-----+--------+-----+
    |       | LOC1 | VAL1 | LOC2 | VAL2 | Continues in Groups of 2 |
    +-------+------+------+------+------+--------------------------+
    """
    type = 'ACCEL'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            #: Load set identification number (Integer>0)
            self.sid = integer(card, 1, 'sid')

            #: Coordinate system identification number. (Integer>0: Default=0)
            self.cid = integer_or_blank(card, 2, 'cid', 0)

            #: Components of the acceleration vector measured in coordinate system
            #: CID. (Real; at least one Ni != 0)
            self.N = array([double_or_blank(card, 3, 'N1', 0.0),
                            double_or_blank(card, 4, 'N2', 0.0),
                            double_or_blank(card, 5, 'N3', 0.0)])
            assert max(abs(self.N)) > 0.

            #: Acceleration vector scale factor. (Real)
            self.dir = string(card, 6, 'dir')
            assert self.dir in ['X', 'Y', 'Z'], 'dir=%r' % self.dir

            i = 9
            locs = []
            vals = []
            j = 0
            while i < len(card):
                loc = double(card, i, 'loc%i' % j)
                val = double(card, i, 'loc%i' % j)
                j += 1
                i += 2
            self.locs = array(locs)
            self.vals = array(vals)
        else:
            raise NotImplementedError(data)

    #def get_reduced_loads(self):
        #scale_factors = [1.]
        #loads = self.F()
        #return(scale_factors, loads)

    def cross_reference(self, model):
        msg = ' which is required by ACCEL sid=%s' % self.sid
        self.cid = model.Coord(self.cid, msg=msg)

    def uncross_reference(self, model):
        self.cid = self.Cid()

    def safe_cross_reference(self, model):
        msg = ' which is required by ACCEL sid=%s' % self.sid
        self.cid = model.Coord(self.cid, msg=msg)

    def Cid(self):
        if isinstance(self.cid, integer_types):
            return self.cid
        return self.cid.cid

    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

    def get_loads(self):
        return [self]

    def raw_fields(self):
        list_fields = [
            'ACCEL', self.sid, self.Cid(),
            self.N[0], self.N[1], self.N[2], self.dir, None, None,
        ]
        for loc, val in zip(self.locs, self.vals):
            list_fields += [loc, val]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)

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
        msg = ' which is required by ACCEL1 sid=%s' % self.sid
        self.cid = model.Coord(self.Cid(), msg=msg)
        self.nodes = model.Nodes(self.node_ids, allowEmptyNodes=True, msg=msg)

    def Cid(self):
        if isinstance(self.cid, integer_types):
            return self.cid
        return self.cid.cid

    @property
    def node_ids(self):
        #msg = ' which is required by ACCEL1 sid=%s' % self.sid
        #_node_ids(self.nodes, allowEmptyNodes=True, msg=msg)
        return self._nodeIDs()

    def _nodeIDs(self, nodes=None):  # this function comes from BaseCard.py
        """returns nodeIDs for repr functions"""
        if not nodes:
            nodes = self.nodes
        if isinstance(nodes[0], integer_types):
            nodeIDs = [node for node in nodes]
        else:
            nodeIDs = [node.nid for node in nodes]
        assert 0 not in nodeIDs, 'nodeIDs = %s' % (nodeIDs)
        return nodeIDs

    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

    def get_loads(self):
        return [self]

    def raw_fields(self):
        list_fields = [
            'ACCEL1', self.sid, self.Cid(), self.scale,
            self.N[0], self.N[1], self.N[2], None, None
            ] + self.node_ids
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


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
        (xyz, matrix) = self.cid.transformToGlobal(self.xyz)
        if self.mag > 0.:
            return (True, self.node, self.mag * xyz)  # load
        return (False, self.node, xyz)  # enforced displacement

    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

    def get_loads(self):
        return [self]

    def F(self):
        return self.xyz * self.mag

    def getReducedLoads(self):
        self.deprecated('getReducedLoads()', 'get_reduced_loads', '0.8')
        return self.get_reduced_loads()

    def get_reduced_loads(self):
        scale_factors = [1.]
        loads = self.F()
        return(scale_factors, loads)

    def organizeLoads(self, model):
        self.deprecated('organizeLoads(model)', 'organize_loads(model)', '0.8')
        return self.organize_loads(model)

    def organize_loads(self, model):
        (scale_factors, forceLoads) = self.get_reduced_loads()

        typesFound = [self.type]
        momentLoads = {}
        forceConstraints = {}
        momentConstraints = {}
        gravityLoads = []
        return (typesFound, forceLoads, momentLoads,
                forceConstraints, momentConstraints,
                gravityLoads)

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


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
            self.xyz /= normXYZ

    def transformLoad(self):
        #print("self.xyz = ",self.xyz)
        (xyz, matrix) = self.cid.transformToGlobal(self.xyz)
        if self.mag > 0.:
            #print("mag=%s xyz=%s" % (self.mag, xyz))
            return (True, self.node, self.mag * xyz)  # load
        return (False, self.node, xyz)  # enforced displacement

    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

    def get_loads(self):
        return [self]

    def getReducedLoads(self):
        self.deprecated('getReducedLoads()', 'get_reduced_loads()', '0.8')
        self.get_reduced_loads()

    def get_reduced_loads(self):
        scale_factors = [1.]
        loads = {self.node: self.M()}
        return(scale_factors, loads)

    def organizeLoads(self, model):
        self.deprecated('organizeLoads(model)', 'organize_loads(model)', '0.8')
        return self.organize_loads(model)

    def organize_loads(self, model):
        (scale_factors, momentLoads) = self.get_reduced_loads()

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

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


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
            self.xyz = array([double_or_blank(card, 5, 'X1', 0.0),
                              double_or_blank(card, 6, 'X2', 0.0),
                              double_or_blank(card, 7, 'X3', 0.0)])
            assert len(card) <= 8, 'len(FORCE card) = %i' % len(card)
        else:
            self.sid = data[0]
            self.node = data[1]
            self.cid = data[2]
            self.mag = data[3]
            xyz = data[4:7]
            assert len(xyz) == 3, 'xyz=%s' % str(xyz)
            self.xyz = array(xyz)

    @property
    def node_id(self):
        if isinstance(self.node, int):
            return self.node
        return self.node.nid

    def Cid(self):
        if isinstance(self.cid, integer_types):
            return self.cid
        return self.cid.cid

    def F(self):
        return {self.node_id : self.mag * self.xyz}

    #def nodeID(self):
        #return self.node

    def cross_reference(self, model):
        """
        .. todo:: cross reference and fix repr function
        """
        msg = ' which is required by FORCE sid=%s' % self.sid
        self.cid = model.Coord(self.cid, msg=msg)

    def uncross_reference(self):
        self.cid = self.Cid()

    def safe_cross_reference(self, model):
        msg = ' which is required by FORCE sid=%s' % self.sid
        # try:
        self.cid = model.Coord(self.cid, msg=msg)
        # except KeyError:
            # pass

    def raw_fields(self):
        list_fields = ['FORCE', self.sid, self.node, self.Cid(),
                       self.mag] + list(self.xyz)
        return list_fields

    def repr_fields(self):
        cid = set_blank_if_default(self.Cid(), 0)
        list_fields = ['FORCE', self.sid, self.node, cid,
                       self.mag] + list(self.xyz)
        return list_fields

    def write_card(self, size=8, is_double=False):
        if size == 8:
            cids = set_string8_blank_if_default(self.Cid(), 0)
            msg = 'FORCE   %8i%8i%8s%8s%8s%8s%8s\n' % (self.sid, self.node,
                cids, print_float_8(self.mag), print_float_8(self.xyz[0]),
                print_float_8(self.xyz[1]), print_float_8(self.xyz[2]))
        else:
            cids = set_string16_blank_if_default(self.Cid(), 0)
            if is_double:
                msg = ('FORCE*  %16i%16i%16s%s\n'
                       '*       %16s%16s%16s') % (self.sid, self.node,
                    cids, print_scientific_double(self.mag), print_scientific_double(self.xyz[0]),
                    print_scientific_double(self.xyz[1]), print_scientific_double(self.xyz[2]))
            else:
                msg = ('FORCE*  %16i%16i%16s%s\n'
                       '*       %16s%16s%16s') % (self.sid, self.node,
                    cids, print_float_16(self.mag), print_float_16(self.xyz[0]),
                    print_float_16(self.xyz[1]), print_float_16(self.xyz[2]))
        return self.comment + msg


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
        msg = ' which is required by FORCE1 sid=%s' % self.sid
        self.node = model.Node(self.node, msg=msg)
        self.g1 = model.Node(self.g1, msg=msg)
        self.g2 = model.Node(self.g2, msg=msg)
        self.xyz = self.g2.get_position() - self.g1.get_position()
        self.normalize()

    def uncross_reference(self):
        self.node =self.NodeID()
        self.g1 = self.G1()
        self.g2 = self.G2()

    def uncross_reference(self):
        self.node =self.NodeID()
        self.g1 = self.G1()
        self.g2 = self.G2()

    def safe_cross_reference(self, model):
        """
        .. todo:: cross reference and fix repr function
        """
        msg = ' which is required by FORCE1 sid=%s' % self.sid
        self.node = model.Node(self.node, msg=msg)
        self.g1 = model.Node(self.g1, msg=msg)
        self.g2 = model.Node(self.g2, msg=msg)
        self.xyz = self.g2.get_position() - self.g1.get_position()
        self.normalize()

    def G1(self):
        if isinstance(self.g1, integer_types) or isinstance(self.g1, float):
            return self.g1
        return self.g1.nid

    def G2(self):
        if isinstance(self.g2, integer_types) or isinstance(self.g1, float):
            return self.g2
        return self.g2.nid

    def NodeID(self):
        self.deprecated('self.NodeID()', 'self.node_id')
        return self.node_id

    @property
    def node_id(self):
        if isinstance(self.node, integer_types):
            return self.node
        return self.node.nid

    def raw_fields(self):
        (node, g1, g2) = self.nodeIDs([self.node, self.g1, self.g2])
        list_fields = ['FORCE1', self.sid, node, self.mag, g1, g2]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class FORCE2(Force):
    """
    Defines a static concentrated force at a grid point by specification of a
    magnitude and four grid points that determine the direction.
    """
    type = 'FORCE2'

    def __init__(self, card=None, data=None, comment=''):
        """
        +--------+-----+---+---+----+----+----+----+
        | FORCE2 | SID | G | F | G1 | G2 | G3 | G4 |
        +--------+-----+---+---+----+----+----+----+
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
        msg = ' which is required by FORCE2 sid=%s' % self.sid
        self.node = model.Node(self.node, msg=msg)
        self.g1 = model.Node(self.g1, msg=msg)
        self.g2 = model.Node(self.g2, msg=msg)
        self.g3 = model.Node(self.g3, msg=msg)
        self.g4 = model.Node(self.g4, msg=msg)

        v12 = self.g2.get_position() - self.g1.get_position()
        v34 = self.g4.get_position() - self.g3.get_position()
        try:
            v12 /= norm(v12)
        except FloatingPointError:
            msg = 'v12=%s norm(v12)=%s\n' % (v12, norm(v12))
            msg += 'g1.get_position()=%s\n' % self.g1.get_position()
            msg += 'g2.get_position()=%s' % self.g2.get_position()
            raise FloatingPointError(msg)

        try:
            v34 /= norm(v34)
        except FloatingPointError:
            msg = 'v34=%s norm(v34)=%s\n' % (v34, norm(v34))
            msg += 'g3.get_position()=%s\n' % self.g3.get_position()
            msg += 'g4.get_position()=%s' % self.g4.get_position()
            raise FloatingPointError(msg)
        self.xyz = cross(v12, v34)
        self.normalize()

    def safe_cross_reference(self, model):
        """
        .. todo:: cross reference and fix repr function
        """
        msg = ' which is required by FORCE2 sid=%s' % self.sid
        self.node = model.Node(self.node, msg=msg)
        self.g1 = model.Node(self.g1, msg=msg)
        self.g2 = model.Node(self.g2, msg=msg)
        self.g3 = model.Node(self.g3, msg=msg)
        self.g4 = model.Node(self.g4, msg=msg)

        v12 = self.g2.get_position() - self.g1.get_position()
        v34 = self.g4.get_position() - self.g3.get_position()
        try:
            v12 /= norm(v12)
        except FloatingPointError:
            msg = 'v12=%s norm(v12)=%s\n' % (v12, norm(v12))
            msg += 'g1.get_position()=%s\n' % self.g1.get_position()
            msg += 'g2.get_position()=%s' % self.g2.get_position()
            raise FloatingPointError(msg)

        try:
            v34 /= norm(v34)
        except FloatingPointError:
            msg = 'v34=%s norm(v34)=%s\n' % (v34, norm(v34))
            msg += 'g3.get_position()=%s\n' % self.g3.get_position()
            msg += 'g4.get_position()=%s' % self.g4.get_position()
            raise FloatingPointError(msg)
        self.xyz = cross(v12, v34)
        self.normalize()

    @property
    def node_id(self):
        if isinstance(self.node, integer_types):
            return self.node
        return self.node.nid

    def NodeID(self):
        self.deprecated('node_id', 'NodeID()', '0.8')
        return self.node_id

    def G1(self):
        if isinstance(self.g1, integer_types):
            return self.g1
        return self.g1.nid

    def G2(self):
        if isinstance(self.g2, integer_types):
            return self.g2
        return self.g2.nid

    def G3(self):
        if isinstance(self.g3, integer_types):
            return self.g3
        return self.g3.nid

    def G4(self):
        if isinstance(self.g4, integer_types):
            return self.g4
        return self.g4.nid

    def raw_fields(self):
        (node, g1, g2, g3, g4) = self.nodeIDs([self.node, self.g1, self.g2, self.g3, self.g4])
        list_fields = ['FORCE2', self.sid, node, self.mag, g1, g2, g3, g4]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class MOMENT(Moment):
    type = 'MOMENT'

    def __init__(self, card=None, data=None, comment=''):
        """
        Defines a static concentrated moment at a grid point by specifying a
        scale factor and a vector that determines the direction.::

        +--------+-----+---+-----+-----+-----+-----+-----+
        | MOMENT | SID | G | CID | M   |  N1 |  N2 |  N3 |
        +--------+-----+---+-----+-----+-----+-----+-----+
        | MOMENT | 2   | 5 |  6  | 2.9 | 0.0 | 1.0 | 0.0 |
        +--------+-----+---+-----+-----+-----+-----+-----+
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
            self.sid = data[0]
            self.node = data[1]
            self.cid = data[2]
            self.mag = data[3]
            xyz = data[4:7]
        assert len(xyz) == 3, 'xyz=%s' % str(xyz)
        self.xyz = xyz

    def Cid(self):
        if isinstance(self.cid, integer_types):
            return self.cid
        return self.cid.cid

    def cross_reference(self, model):
        """
        .. todo:: cross reference and fix repr function
        """
        #msg = ' which is required by MOMENT sid=%s' % self.sid
        pass

    @property
    def node_id(self):
        if isinstance(self.node, integer_types):
            return self.node
        return self.node.nid

    def safe_cross_reference(self, model):
        pass

    def raw_fields(self):
        list_fields = ['MOMENT', self.sid, self.node, self.Cid(),
                       self.mag] + list(self.xyz)
        return list_fields

    def repr_fields(self):
        cid = set_blank_if_default(self.Cid(), 0)
        list_fields = ['MOMENT', self.sid, self.node, cid,
                       self.mag] + list(self.xyz)
        return list_fields

    def write_card(self, size=8, is_double=False):
        if size == 8:
            cids = set_string8_blank_if_default(self.Cid(), 0)
            msg = 'MOMENT  %8i%8i%8s%8s%8s%8s%8s\n' % (self.sid, self.node,
                cids, print_float_8(self.mag), print_float_8(self.xyz[0]),
                print_float_8(self.xyz[1]), print_float_8(self.xyz[2]))
        else:
            cids = set_string16_blank_if_default(self.Cid(), 0)
            if is_double:
                msg = ('MOMENT* %16i%16i%16s%s\n'
                       '*       %16s%16s%16s') % (self.sid, self.node,
                    cids, print_scientific_double(self.mag), print_scientific_double(self.xyz[0]),
                    print_scientific_double(self.xyz[1]), print_scientific_double(self.xyz[2]))
            else:
                msg = ('MOMENT* %16i%16i%16s%s\n'
                       '*       %16s%16s%16s') % (self.sid, self.node,
                    cids, print_float_16(self.mag), print_float_16(self.xyz[0]),
                    print_float_16(self.xyz[1]), print_float_16(self.xyz[2]))
        return self.comment + msg


class MOMENT1(Moment):
    type = 'MOMENT1'

    def __init__(self, card=None, data=None, comment=''):
        """
        Defines a static concentrated moment at a grid point by specifying a
        magnitude and two grid points that determine the direction.::

        +---------+-----+---+---+----+----+
        | MOMENT1 | SID | G | M | G1 | G2 |
        +---------+-----+---+---+----+----+
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

        #assert len(xyz) == 3, 'xyz=%s' % str(xyz)
        #self.xyz = array(xyz)
        self.xyz = None

    def cross_reference(self, model):
        """
        .. todo:: cross reference and fix repr function
        """
        msg = ' which is required by MOMENT1 sid=%s' % self.sid
        self.node = model.Node(self.node, msg=msg)
        self.g1 = model.Node(self.g1, msg=msg)
        self.g2 = model.Node(self.g2, msg=msg)
        self.xyz = self.g2.get_position() - self.g1.get_position()
        self.normalize()

    def uncross_reference(self):
        self.node =self.NodeID()
        self.g1 = self.G1()
        self.g2 = self.G2()

    def safe_cross_reference(self, model):
        """
        .. todo:: cross reference and fix repr function
        """
        msg = ' which is required by MOMENT1 sid=%s' % self.sid
        self.node = model.Node(self.node, msg=msg)
        self.g1 = model.Node(self.g1, msg=msg)
        self.g2 = model.Node(self.g2, msg=msg)
        self.xyz = self.g2.get_position() - self.g1.get_position()
        self.normalize()

    def get_node_id(self):
        if isinstance(self.node, integer_types):
            return self.node
        return self.node.nid

    def G1(self):
        if isinstance(self.g1, integer_types):
            return self.g1
        return self.g1.nid

    def G2(self):
        if isinstance(self.g2, integer_types):
            return self.g2
        return self.g2.nid

    def raw_fields(self):
        node = self.get_node_id()
        g1 = self.G1()
        g2 = self.G2()
        list_fields = ['MOMENT1', self.sid, node, self.mag, g1, g2]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class MOMENT2(Moment):
    type = 'MOMENT2'

    def __init__(self, card=None, data=None, comment=''):
        """
        Defines a static concentrated moment at a grid point by specification
        of a magnitude and four grid points that determine the direction.::

        +---------+-----+---+---+----+----+----+----+
        | MOMENT2 | SID | G | M | G1 | G2 | G3 | G4 |
        +---------+-----+---+---+----+----+----+----+
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
            self.xyz = None
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
            assert len(xyz) == 3, 'xyz=%s' % str(xyz)

    def cross_reference(self, model):
        """
        .. todo:: cross reference and fix repr function
        """
        msg = ' which is required by MOMENT2 sid=%s' % self.sid
        self.node = model.Node(self.node, msg=msg)
        self.g1 = model.Node(self.g1, msg=msg)
        self.g2 = model.Node(self.g2, msg=msg)
        self.g3 = model.Node(self.g3, msg=msg)
        self.g4 = model.Node(self.g4, msg=msg)

        v12 = self.g2.get_position() - self.g1.get_position()
        v34 = self.g4.get_position() - self.g3.get_position()
        v12 = v12 / norm(v12)
        v34 = v34 / norm(v34)
        self.xyz = cross(v12, v34)

    def uncross_reference(self):
        self.node =self.NodeID()
        self.g1 = self.G1()
        self.g2 = self.G2()
        self.g3 = self.G3()
        self.g4 = self.G4()

    def safe_cross_reference(self, model):
        """
        .. todo:: cross reference and fix repr function
        """
        msg = ' which is required by MOMENT2 sid=%s' % self.sid
        self.node = model.Node(self.node, msg=msg)
        self.g1 = model.Node(self.g1, msg=msg)
        self.g2 = model.Node(self.g2, msg=msg)
        self.g3 = model.Node(self.g3, msg=msg)
        self.g4 = model.Node(self.g4, msg=msg)

        v12 = self.g2.get_position() - self.g1.get_position()
        v34 = self.g4.get_position() - self.g3.get_position()
        v12 = v12 / norm(v12)
        v34 = v34 / norm(v34)
        self.xyz = cross(v12, v34)

    def NodeID(self):
        if isinstance(self.node, integer_types):
            return self.node
        return self.node.nid

    def G1(self):
        if isinstance(self.g1, integer_types):
            return self.g1
        return self.g1.nid

    def G2(self):
        if isinstance(self.g2, integer_types):
            return self.g2
        return self.g2.nid

    def G3(self):
        if isinstance(self.g3, integer_types):
            return self.g3
        return self.g3.nid

    def G4(self):
        if isinstance(self.g4, integer_types):
            return self.g4
        return self.g4.nid

    def raw_fields(self):
        (node, g1, g2, g3, g4) = self.nodeIDs(nodes=[self.node, self.g1, self.g2,
                                                     self.g3, self.g4])
        assert isinstance(g1, integer_types), g1
        list_fields = ['MOMENT2', self.sid, node, self.mag, g1, g2, g3, g4]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class GMLOAD(Load):
    """
    Defines a static concentrated force at a grid point by specification of a
    magnitude and two grid points that determine the direction.
    """
    type = 'GMLOAD'

    def __init__(self, card=None, data=None, comment=''):
        Load.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.sid = integer(card, 1, 'sid')
            self.cid = integer_or_blank(card, 2, 'cid', 0)
            self.normal = array([
                double_or_blank(card, 3, 'N1', 0.),
                double_or_blank(card, 4, 'N2', 0.),
                double_or_blank(card, 5, 'N3', 1.),
            ])
            self.entity = string(card, 6, 'entity')
            self.entity_id = integer(card, 7, 'entity_id')
            self.method = string(card, 8, 'method')

            self.load_magnitudes = []
            for i in range(9, len(card)):
                ifield = i - 8
                load_mag = integer_or_double(card, i, 'load_magnitude_%s' % ifield)
                self.load_magnitudes.append(load_mag)
        else:
            raise NotImplemented()

    #def DEquation(self):
        #if isinstance(self.dequation, int):
            #return self.dequation
        #return self.dequation.equation_id

    def cross_reference(self, model):
        """
        .. todo:: cross reference and fix repr function
        """
        msg = ' which is required by GMLOAD sid=%s' % self.sid
        self.cid = model.Coord(self.Cid(), msg=msg)
        #self.node = model.Node(self.node, msg=msg)
        #self.g1 = model.Node(self.g1, msg=msg)
        #self.g2 = model.Node(self.g2, msg=msg)
        #self.xyz = self.g2.get_position() - self.g1.get_position()
        #self.normalize()

    #def G1(self):
        #if isinstance(self.g1, integer_types) or isinstance(self.g1, float):
            #return self.g1
        #return self.g1.nid

    #def G2(self):
        #if isinstance(self.g2, integer_types) or isinstance(self.g1, float):
            #return self.g2
        #return self.g2.nid

    #def NodeID(self):
        #if isinstance(self.node, integer_types):
            #return self.node
        #return self.node.nid

    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

    def get_loads(self):
        return [self]

    def raw_fields(self):
        list_fields = ['GMLOAD', self.sid, self.Cid()] + list(self.normal) + [
            self.entity, self.entity_id, self.method] + self.load_magnitudes
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


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
            n4 = integer_or_blank(card, 6, 'n4', 0)
            if n4:
                nodes.append(n4)
            #self.nodes = wipe_empty_fields(nodes)
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
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

    def get_loads(self):
        return [self]

    def raw_fields(self):
        list_fields = ['PLOAD', self.sid, self.p] + self.node_ids
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


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
            self.Type = string(card, 3, 'Type ("%s")' % '",  "'.join(self.validTypes))
            self.scale = string(card, 4, 'scale ("%s")' % '", "'.join(self.validScales))
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
            self.Type = self.validTypes[self.Type - 1]
            self.scale = self.validScales[self.scale - 1]

        if self.Type not in self.validTypes:
            msg = '%s is an invalid type on the PLOAD1 card; validTypes=[%s]' % (
                self.Type, ', '.join(self.validTypes).rstrip(', '))
            raise RuntimeError(msg)
        if self.scale not in self.validScales:
            msg = '%s is an invalid scale on the PLOAD1 card; validScales=[%s]' % (
                self.scale, ', '.join(self.validScales).rstrip(', '))
            raise RuntimeError(msg)

        assert 0.0 <= self.x1 <= self.x2, '0.0 <= x1 <= x2 -> x1=%s x2=%s' % (self.x1, self.x2)
        if self.scale in ['FR', 'FRPR']:
            assert self.x1 <= 1.0, 'x1=%r' % self.x1
            assert self.x2 <= 1.0, 'x2=%r' % self.x2
        assert self.scale in self.validScales, '%s is an invalid scale on the PLOAD1 card' % (self.scale)


    def cross_reference(self, model):
        """
        .. todo:: cross reference and fix repr function
        """
        msg = ' which is required by PLOAD1 sid=%s' % self.sid
        self.eid = model.Element(self.eid, msg=msg)

    def transformLoad(self):
        p1 = self.eid.ga.get_position()
        p2 = self.eid.gb.get_position()

        g0 = self.eid.g0_vector
        #if not isinstance(g0, ndarray):
        #    g0 = g0.get_position()

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
        self.deprecated('getReducedLoads()', 'get_reduced_loads()', '0.8')
        return self.get_reduced_loads()

    def get_reduced_loads(self):
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
        self.deprecated('organizeLoads(model)', 'organize_loads(model)', '0.8')
        return self.organize_loads(model)

    def organize_loads(self, model):
        """
        Figures out magnitudes of the loads to be applied to the various nodes.
        This includes figuring out scale factors.
        """
        forceLoads = {}  # spc enforced displacement (e.g. FORCE=0)
        momentLoads = {}
        forceConstraints = {}
        momentConstraints = {}
        gravityLoads = []

        typesFound = set()
        (scale_factors, loads) = self.get_reduced_loads()

        for scale_factor, load in zip(scale_factors, loads):
            out = load.transformLoad()
            typesFound.add(load.__class__.__name__)

            if isinstance(load, PLOAD1): # CBAR/CBEAM
                element = load.eid
                (ga, gb) = element.nodeIDs()
                load_type = load.Type

                scale = load.scale
                eType = element.type

                if load_type in ['FX', 'FY', 'FZ']:
                    p1 = element.ga.get_position()
                    p2 = element.gb.get_position()
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
                    msg = ('Point loads are not supported on...\n%s'
                           'Try setting x1=%r very close to x2=%r and\n'
                           'scaling p1=%r and p2=%r by x2-x1 (for "FR") and (x2-x1)/L (for "LE").'
                           % (repr(self), self.x1, self.x2, self.p1, self.p2))
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
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

    def get_loads(self):
        return [self]

    def Eid(self):
        if isinstance(self.eid, integer_types):
            return self.eid
        return self.eid.eid

    def raw_fields(self):
        list_fields = ['PLOAD1', self.sid, self.Eid(), self.Type, self.scale,
                       self.x1, self.p1, self.x2, self.p2]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class PLOAD2(Load):
    type = 'PLOAD2'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            self.sid = integer(card, 1, 'sid')
            self.pressure = double(card, 2, 'p')

            if integer_string_or_blank(card, 4, 'THRU') == 'THRU':
                e1 = integer(card, 3, 'Element1')
                e2 = integer(card, 5, 'Element1')
                eids = [i for i in range(e1, e2 + 1)]
                assert len(card) == 6, 'len(PLOAD2 card) = %i' % len(card)
            else:
                eids = fields(integer, card, 'eid', i=3, j=len(card))
            self.eids = eids
        else:
            self.sid = data[0]
            self.pressure = data[1]
            self.eids = list(data[2:])

    def cross_reference(self, model):
        """
        .. todo:: cross reference and fix repr function
        """
        pass

    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

    def get_loads(self):
        return [self]

    def raw_fields(self):
        list_fields = ['PLOAD2', self.sid, self.pressure]
        if len(self.eids) > 6:
            list_fields += [self.eids[0], 'THRU', self.eids[-1]]
        else:
            list_fields += self.eids
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


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
                    self.eids = list(unique(
                        expand_thru([self.eid, 'THRU', eid2], set_fields=False, sort_fields=False)
                    ))
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
            self.NVector = array([double_or_blank(card, 10, 'N1', 0.0),
                                  double_or_blank(card, 11, 'N2', 0.0),
                                  double_or_blank(card, 12, 'N3', 0.0)])
            self.sorl = string_or_blank(card, 13, 'sorl', 'SURF')
            self.ldir = string_or_blank(card, 14, 'ldir', 'NORM')
            assert len(card) <= 15, 'len(PLOAD4 card) = %i' % len(card)
        else:
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
            self.sorl = 'SURF'
            self.ldir = 'NORM'

        if self.sorl not in ['SURF', 'LINE']:
            raise RuntimeError(self.sorl)
        if self.ldir not in ['LINE', 'X', 'Y', 'Z', 'TANG', 'NORM']:
            raise RuntimeError(self.ldir)


    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

    def get_loads(self):
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
            (faceNodeIDs, Area) = self.eid.getFaceNodesAndArea(self, nid, nidOpposite)
        else:
            faceNodeIDs = self.eid.nodeIDs()
            Area = self.eid.Area()
        n = len(faceNodeIDs)

        vector = array(self.eid.Normal())
        vectors = []
        for (nid, p) in zip(faceNodeIDs, self.pressures):
            #: .. warning:: only supports normal pressures
            vectors.append(vector * p * Area / n)  # Force_i

        isLoad = None
        return (isLoad, faceNodeIDs, vectors)

    def Cid(self):
        if isinstance(self.cid, integer_types):
            return self.cid
        return self.cid.cid

    def cross_reference(self, model):
        msg = ' which is required by PLOAD4 sid=%s' % self.sid
        self.eid = model.Element(self.eid, msg=msg)
        self.cid = model.Coord(self.cid, msg=msg)
        if self.g1 is not None:
            self.g1 = model.Node(self.g1, msg=msg)
        if self.g34 is not None:
            self.g34 = model.Node(self.g34, msg=msg)
        if self.eids:
            self.eids = model.Elements(self.eids, msg=msg)

    def safe_cross_reference(self, model):
        msg = ' which is required by PLOAD4 sid=%s' % self.sid
        self.eid = model.Element(self.eid, msg=msg)
        self.cid = model.Coord(self.cid, msg=msg)
        if self.g1 is not None:
            self.g1 = model.Node(self.g1, msg=msg)
        if self.g34 is not None:
            self.g34 = model.Node(self.g34, msg=msg)
        if self.eids:
            self.eids = model.Elements(self.eids, msg=msg)

    def Eid(self, eid):
        if isinstance(eid, integer_types) or isinstance(eid, int32):
            return eid
        return eid.eid

    def nodeIDs(self, nodes=None):
        nodeIDs = [None, None]
        if isinstance(self.g1, integer_types):
            nodeIDs[0] = self.g1
        elif self.g1 is not None:
            nodeIDs[0] = self.g1.nid

        if isinstance(self.g34, integer_types):
            nodeIDs[1] = self.g34
        elif self.g34 is not None:
            nodeIDs[1] = self.g34.nid
        return nodeIDs

    def getElementIDs(self, eid=None):
        if eid:
            return self.Eid(eid)
        eids = []
        for element in self.eids:
            eids.append(self.Eid(element))
        return eids

    def raw_fields(self):
        eid = self.Eid(self.eid)
        cid = set_blank_if_default(self.Cid(), 0)
        sorl = set_blank_if_default(self.sorl, 'SURF')
        ldir = set_blank_if_default(self.ldir, 'NORM')
        p1 = self.pressures[0]
        p2 = set_blank_if_default(self.pressures[1], p1)
        p3 = set_blank_if_default(self.pressures[2], p1)
        p4 = set_blank_if_default(self.pressures[3], p1)
        list_fields = ['PLOAD4', self.sid, eid, self.pressures[0], p2, p3, p4]

        #print "g3=%r g4=%s eids=%r" %(self.g3, self.g4, self.eids)
        if self.g1 is not None:  # is it a SOLID element
            nodeIDs = self.nodeIDs([self.g1, self.g34])
            list_fields += nodeIDs
        else:
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

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


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
        #msg = ' which is required by PLOADX1 lid=%s' % self.sid
        #self.eid = model.Element(self.eid)
        #self.ga = model.Node(self.ga)
        #self.gb = model.Node(self.gb)
        pass

    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

    def get_loads(self):
        return [self]

    def raw_fields(self):
        list_fields = ['PLOADX1', self.sid, self.eid, self.pa, self.pb,
                  self.ga, self.gb, self.theta]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)
