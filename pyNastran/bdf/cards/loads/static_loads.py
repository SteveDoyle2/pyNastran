# pylint: disable=R0902,R0904,R0914
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
from six.moves import zip

from numpy import array, cross, allclose, unique
from numpy.linalg import norm

#from pyNastran.bdf.errors import CrossReferenceError
from pyNastran.utils import integer_types
from pyNastran.bdf.cards.loads.loads import Load, LoadCombination
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import BaseCard, expand_thru, expand_thru_by, range
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, string, string_or_blank,
    integer_or_string, fields, integer_string_or_blank, integer_or_double)
from pyNastran.bdf.field_writer_8 import print_card_8, print_float_8, set_string8_blank_if_default
from pyNastran.bdf.field_writer_16 import print_card_16, print_float_16, set_string16_blank_if_default
from pyNastran.bdf.field_writer_double import print_card_double, print_scientific_double


class LOAD(LoadCombination):
    type = 'LOAD'

    def __init__(self, sid, scale, scale_factors, load_ids, comment=''):
        LoadCombination.__init__(self, sid, scale, scale_factors, load_ids,
                                 comment=comment)

    def getLoadIDs(self):
        self.deprecated('getLoadIDs()', 'get_load_ids()', '0.8')
        return self.get_load_ids()

    def getReducedLoads(self):
        self.deprecated('self.getReducedLoads()', 'self.get_reduced_loads()', '0.8')
        return self.get_reduced_loads()

    def organizeLoads(self, model):
        self.deprecated('organizeLoads(model)', 'organize_loads(model)', '0.8')
        return self.organize_loads(model)

    def getLoadTypes(self):
        self.deprecated('getLoadTypes()', 'get_load_types()', '0.8')
        return self.get_load_types()

    def get_load_ids(self):
        """
        .. note:: requires a cross referenced load
        """
        load_ids = []
        for loads in self.load_ids:
            for load in loads:
                #if isinstance(load, int):
                    #load_ids += [load]

                if isinstance(load, LOAD):
                    lid = load.lid
                    if isinstance(lid, list):
                        load_ids += load.lid
                    else:  # int
                        load_ids += load.get_load_ids()
                elif isinstance(load, (Force, Moment, PLOAD4, GRAV)):
                    load_ids += [load.sid]
                else:
                    msg = ('The get_load_ids method doesnt support %s cards.\n'
                           '%s' % (load.__class__.__name__, str(load)))
                    raise NotImplementedError(msg)

        load_ids = list(set(load_ids))
        #print("load_ids = ", load_ids)
        return load_ids

    def get_load_types(self):
        """
        .. note:: requires a cross referenced load
        """
        load_types = []
        for loads in self.load_ids:
            for load in loads:
                if isinstance(load, LOAD):
                    lid = load.lid
                    if isinstance(lid, list):
                        load_types += load.type
                    else:  # int
                        load_types += [load.type] + load.getLoadTypes()
                elif isinstance(load, (Force, Moment, PLOAD4, GRAV)):
                    load_types += [load.type]
                else:
                    raise NotImplementedError(load)

        load_types = list(set(load_types))
        #print("load_types = ", load_types)
        return load_types

    def write_calculix_grav(self, gx, gy, gz):
        msg = '*DLOAD\n'
        msg += 'AllElements,GRAV,%s,%s,%s\n' % (gx, gy, gz)
        return msg

    def write_code_aster_load(self, model, grid_word='node'):
        load_ids = self.get_load_ids()
        load_types = self.get_load_types()

        #msg = '# Loads\n'
        msg = ''
        (types_found, force_loads, moment_loads,
         force_constraints, moment_constraints,
         gravity_loads) = self.organize_loads(model)

        nids = []
        for nid in force_loads:
            nids.append(nid)
        for nid in moment_loads:
            nids.append(nid)

        if nids:
            msg += '# types_found = %s\n' % (list(types_found))
            msg += '# load_ids    = %s\n' % (load_ids)
            msg += "load_bc = AFFE_CHAR_MECA(MODELE=modmod,\n"
            #msg += "                        DDL_IMPO=(_F(GROUP_MA='Lleft',\n"
            msg += "                         FORCE_NODALE=(\n"

        #CHAR=AFFE_CHAR_MECA(MODELE=MODE,
        #             FORCE_NODALE=(
        #                     _F(NOEUD='N1',
        #                        FZ=-500.0),)

        spaces = "                           "
        for nid in sorted(nids):  # ,load in sorted(iteritems(force_loads))
            msg += spaces + "_F(NOEUD='%s%s'," % (grid_word, nid)

            if nid in force_loads:
                force = force_loads[nid]
                if abs(force[0]) > 0.:
                    msg += " FX=%s," % force[0]
                if abs(force[1]) > 0.:
                    msg += " FY=%s," % force[1]
                if abs(force[2]) > 0.:
                    msg += " FZ=%s," % force[2]

            if nid in moment_loads:
                moment = moment_loads[nid]
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

        for gravity_load in gravity_loads:
            msg += 'CA_GRAVITY(%s);\n' % str(gravity_load)
        return msg, load_ids, load_types

    def get_reduced_loads(self):
        """
        Get all load objects in a simplified form, which means all
        scale factors are already applied and only base objects
        (no LOAD cards) will be returned.

        .. todo:: lots more object types to support
        """
        scale_factors = []
        loads = []
        simple_loads = [
            'FORCE', 'FORCE1', 'FORCE2',
            'MOMENT', 'MOMENT1', 'MOMENT2',
            'PLOAD1', 'PLOAD2', 'PLOAD4',
            'GRAV', 'ACCEL', 'ACCEL1']
        load_scale = self.scale # global
        for (loads_pack, i_scale) in zip(self.load_ids, self.scale_factors):
            scale = i_scale * load_scale # actual scale = global * local
            if isinstance(loads_pack, integer_types):
                raise RuntimeError('the load have not been cross-referenced')

            for load in loads_pack:
                if simple_loads:
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

    def organize_loads(self, model):
        """
        Figures out magnitudes of the loads to be applied to the various nodes.
        This includes figuring out scale factors.
        """
        force_loads = {}  # spc enforced displacement (e.g. FORCE=0)
        moment_loads = {}
        force_constraints = {}
        moment_constraints = {}
        gravity_loads = []

        types_found = set()
        (scale_factors, loads) = self.get_reduced_loads()

        for (scale_factor, load) in zip(scale_factors, loads):
            out = load.transform_load()
            types_found.add(load.__class__.__name__)
            if isinstance(load, Force):
                (is_load, node, vector) = out
                if is_load:  # load
                    if node not in force_loads:
                        force_loads[node] = vector * scale_factor
                    else:
                        force_loads[node] += vector * scale_factor
                else:  # constraint
                    if node not in force_loads:
                        force_constraints[node] = vector * scale_factor
                    else:
                        force_constraints[node] += vector * scale_factor

            elif isinstance(load, Moment):
                (is_load, node, vector) = out
                if is_load:  # load
                    if node not in moment_loads:
                        moment_loads[node] = vector * scale_factor
                    else:
                        moment_loads[node] += vector * scale_factor
                else:  # constraint
                    if node not in moment_loads:
                        moment_constraints[node] = vector * scale_factor
                    else:
                        moment_constraints[node] += vector * scale_factor

            elif isinstance(load, PLOAD4):
                (is_load, nodes, vectors) = out
                for (nid, vector) in zip(nodes, vectors):
                    # not the same vector for all nodes
                    force_loads[nid] = vector * scale_factor

            elif isinstance(load, GRAV):
                #(grav) = out
                gravity_loads.append(out * scale_factor)  # grav
            else:
                msg = '%s not supported' % (load.__class__.__name__)
                raise NotImplementedError(msg)

        return (types_found, force_loads, moment_loads, force_constraints,
                moment_constraints, gravity_loads)

    def raw_fields(self):
        list_fields = ['LOAD', self.sid, self.scale]
        for (scale_factor, load_id) in zip(self.scale_factors, self.load_ids):
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
        for i, load_id in enumerate(self.load_ids):
            idi = self.LoadID(load_id)
            ids.append(idi)
        self.load_ids = ids
        #assert ids == ['cat'], ids
        del self.load_ids_ref


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

    def __init__(self, sid, cid, scale, N, mb, comment=''):
        if comment:
            self._comment = comment

        #: Set identification number
        self.sid = sid

        #: Coordinate system identification number.
        self.cid = cid

        #: scale factor
        self.scale = scale

        #: Acceleration vector components measured in coordinate system CID
        self.N = N

        #: Indicates whether the CID coordinate system is defined in the
        #: main Bulk Data Section (MB = -1) or the partitioned superelement
        #: Bulk Data Section (MB = 0). Coordinate systems referenced in the
        #: main Bulk Data Section are considered stationary with respect to
        #: the assembly basic coordinate system. See Remark 10.
        #: (Integer; Default = 0)
        self.mb = mb

        assert not allclose(max(abs(self.N)), 0.), ('GRAV N is a zero vector, '
                                                    'N=%s' % str(self.N))

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        cid = integer_or_blank(card, 2, 'cid', 0)
        scale = double(card, 3, 'scale')
        N = array([double_or_blank(card, 4, 'N1', 0.0),
                   double_or_blank(card, 5, 'N2', 0.0),
                   double_or_blank(card, 6, 'N3', 0.0)])
        mb = integer_or_blank(card, 7, 'mb', 0)
        assert len(card) <= 8, 'len(GRAV card) = %i\ncard=%s' % (len(card), card)
        return GRAV(sid, cid, scale, N, mb, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        sid = data[0]
        cid = data[1]
        a = data[2]
        N = array(data[3:6])
        mb = data[6]
        scale = 1.
        assert len(data) == 7
        return GRAV(sid, cid, scale, N, mb, comment=comment)

    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

    def get_loads(self):
        return [self]

    def organizeLoads(self, model):
        self.deprecated('organizeLoads(model)', 'organize_loads(model)', '0.8')
        return self.organize_loads(model)

    def organize_loads(self, model):
        types_found = [self.type]
        force_loads = {}
        moment_loads = {}
        force_constraints = {}
        moment_constraints = {}
        gravity_load = self.transform_load()
        return (types_found, force_loads, moment_loads,
                force_constraints, moment_constraints,
                gravity_load)

    def transformLoad(self):
        self.deprecated('transformLoad()', 'transform_load()', '0.8')
        return self.transform_load()

    def transform_load(self):
        g = self.GravityVector()
        g2 = self.cid_ref.transform_node_to_global(g)
        return g2

    #write_code_aster_load(self,mag):
        #p = self.GravityVector()
        #msg = 'GRAV([%s,%s,%s])' %(p)
        #return msg

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by GRAV sid=%s' % self.sid
        self.cid = model.Coord(self.cid, msg=msg)
        self.cid_ref = self.cid

    def safe_cross_reference(self, model, debug=True):
        # msg = "Couldn't find CORDx=%s which is required by GRAV sid=%s" % (self.cid, self.sid)
        msg = ' which is required by GRAV sid=%s' % self.sid
        self.cid = model.Coord(self.cid, msg=msg)
        self.cid_ref = self.cid

    def uncross_reference(self):
        self.cid = self.Cid()
        del self.cid_ref

    def Cid(self):
        if isinstance(self.cid, integer_types):
            return self.cid
        return self.cid_ref.cid

    def GravityVector(self):
        """returns the gravity vector in absolute coordinates"""
        ## TODO: shouldn't be scaled by the ???
        p = self.cid_ref.transform_vector_to_global(self.N)
        return self.scale * p

    def raw_fields(self):
        N = list(self.N)
        list_fields = ['GRAV', self.sid, self.Cid(), self.scale] + N + [self.mb]
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
    |   1   |   2  |   3  |  4   |   5  |  6  |  7  |   8    |  9  |
    +=======+======+======+======+======+=====+=====+========+=====+
    | ACCEL | SID  | CID  | N1   | N2   | N3  | DIR |        |     |
    +-------+------+------+------+------+-----+-----+--------+-----+
    |       | LOC1 | VAL1 | LOC2 | VAL2 | Continues in Groups of 2 |
    +-------+------+------+------+------+--------------------------+
    """
    type = 'ACCEL'

    def __init__(self, sid, cid, N, direction, locs, vals, comment=''):
        if comment:
            self._comment = comment
        #: Load set identification number (Integer>0)
        self.sid = sid
        #: Coordinate system identification number. (Integer>0: Default=0)
        self.cid = cid

        #: Components of the acceleration vector measured in coordinate system
        #: CID. (Real; at least one Ni != 0)
        self.N = N

        #: Component direction of acceleration variation. (Character; one of X,Y or Z)
        self.direction = direction
        self.locs = array(locs, dtype='float64')
        self.vals = array(vals, dtype='float64')

        assert max(abs(self.N)) > 0.
        assert self.direction in ['X', 'Y', 'Z'], 'dir=%r' % self.direction

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        cid = integer_or_blank(card, 2, 'cid', 0)
        N = array([double_or_blank(card, 3, 'N1', 0.0),
                   double_or_blank(card, 4, 'N2', 0.0),
                   double_or_blank(card, 5, 'N3', 0.0)])
        direction = string(card, 6, 'dir')

        i = 9
        locs = []
        vals = []
        j = 0
        nfields = len(card)
        while i < nfields:
            #raise NotImplementedError('ACCEL-line 2')
            loc = double(card, i, 'loc%i' % j)
            val = double(card, i, 'loc%i' % j)
            #print('i=%s j=%s len=%s loc=%s val=%s' % (i, j, len(card), loc, val))
            locs.append(loc)
            vals.append(val)
            j += 1
            i += 2
        return ACCEL(sid, cid, N, direction, locs, vals, comment=comment)

    #def get_reduced_loads(self):
        #scale_factors = [1.]
        #loads = self.F()
        #return(scale_factors, loads)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by ACCEL sid=%s' % self.sid
        self.cid = model.Coord(self.cid, msg=msg)
        self.cid_ref = self.cid

    def uncross_reference(self):
        self.cid = self.Cid()
        del self.cid_ref

    def safe_cross_reference(self, model, debug=True):
        msg = ' which is required by ACCEL sid=%s' % self.sid
        self.cid = model.Coord(self.cid, msg=msg)

    def Cid(self):
        if isinstance(self.cid, integer_types):
            return self.cid
        return self.cid_ref.cid

    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

    def get_loads(self):
        return [self]

    def raw_fields(self):
        list_fields = [
            'ACCEL', self.sid, self.Cid(),
            self.N[0], self.N[1], self.N[2], self.direction, None, None,
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

    def __init__(self, sid, cid, scale, N, nodes, comment=''):
        if comment:
            self._comment = comment
        #: Load set identification number (Integer>0)
        self.sid = sid

        #: Coordinate system identification number. (Integer>0: Default=0)
        self.cid = cid

        #: Acceleration vector scale factor. (Real)
        self.scale = scale

        #: Components of the acceleration vector measured in coordinate system
        #: CID. (Real; at least one Ni != 0)
        self.N = N

        #: nodes to apply the acceleration to
        self.nodes = expand_thru_by(nodes)

        assert max(abs(self.N)) > 0.

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        cid = integer_or_blank(card, 2, 'cid', 0)
        scale = double(card, 3, 'scale')
        N = array([double_or_blank(card, 4, 'N1', 0.0),
                   double_or_blank(card, 5, 'N2', 0.0),
                   double_or_blank(card, 6, 'N3', 0.0)])

        nodes = fields(integer_or_string, card, 'node', i=9, j=len(card))
        return ACCEL1(sid, cid, scale, N, nodes, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by ACCEL1 sid=%s' % self.sid
        self.cid = model.Coord(self.Cid(), msg=msg)
        self.nodes = model.Nodes(self.node_ids, allow_empty_nodes=True, msg=msg)
        self.cid_ref = self.cid
        self.nodes_ref = self.nodes

    def safe_cross_reference(self, model):
        return self.cross_reference(model)

    def Cid(self):
        if isinstance(self.cid, integer_types):
            return self.cid
        return self.cid_ref.cid

    @property
    def node_ids(self):
        #msg = ' which is required by ACCEL1 sid=%s' % self.sid
        #_node_ids(self.nodes, allow_empty_nodes=True, msg=msg)
        return self._nodeIDs()

    def _nodeIDs(self, nodes=None):  # this function comes from BaseCard.py
        """returns node_ids for repr functions"""
        if not nodes:
            nodes = self.nodes
        if isinstance(nodes[0], integer_types):
            node_ids = [node for node in nodes]
        else:
            node_ids = [node.nid for node in nodes]
        assert 0 not in node_ids, 'node_ids = %s' % (node_ids)
        return node_ids

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

    def __init__(self):
        Load.__init__(self)

    def normalize(self):
        """
        adjust the vector to a unit length
        scale up the magnitude of the vector
        """
        if self.mag != 0.0:  # enforced displacement
            norm_xyz = norm(self.xyz)
            #mag = self.mag * norm_xyz
            self.mag *= norm_xyz
            try:
                self.xyz = self.xyz / norm_xyz
            except FloatingPointError:
                msg = 'xyz = %s\n' % self.xyz
                msg += 'norm_xyz = %s\n' % norm_xyz
                msg += 'card =\n%s' % str(self)
                raise FloatingPointError(msg)


    def transformLoad(self):
        self.deprecated('transformLoad()', 'transform_load()', '0.8')
        return self.transform_load()

    def transform_load(self):
        xyz = self.cid_ref.transform_node_to_global(self.xyz)
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
        (scale_factors, force_loads) = self.get_reduced_loads()

        types_found = [self.type]
        moment_loads = {}
        force_constraints = {}
        moment_constraints = {}
        gravity_loads = []
        return (types_found, force_loads, moment_loads,
                force_constraints, moment_constraints,
                gravity_loads)

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

    def __init__(self):
        Load.__init__(self)

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
        self.deprecated('transformLoad()', 'transform_load()', '0.8')
        return self.transform_load()

    def transform_load(self):
        #print("self.xyz = ",self.xyz)
        xyz = self.cid_ref.transform_node_to_global(self.xyz)
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
        loads = {
            self.node: self.M()
        }
        return(scale_factors, loads)

    def organizeLoads(self, model):
        self.deprecated('organizeLoads(model)', 'organize_loads(model)', '0.8')
        return self.organize_loads(model)

    def organize_loads(self, model):
        (scale_factors, moment_loads) = self.get_reduced_loads()

        types_found = [self.type]
        force_loads = {}
        force_constraints = {}
        moment_constraints = {}
        gravity_loads = []
        return (types_found, force_loads, moment_loads,
                force_constraints, moment_constraints,
                gravity_loads)

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

    def __init__(self, sid, node, cid, mag, xyz, comment=''):
        """
        +-------+-----+------+-------+------+------+------+------+
        |   1   |  2  |  3   |   4   |  5   |  6   |   7  |   8  |
        +-------+-----+------+-------+------+------+------+------+
        | FORCE | SID | NODE | CID   | MAG  |  FX  |  FY  |  FZ  |
        +-------+-----+------+-------+------+------+------+------+

        +-------+-----+------+-------+------+------+------+------+
        | FORCE |  3  |  1   |       | 100. |  0.  |  0.  |  1.  |
        +-------+-----+------+-------+------+------+------+------+
        """
        Force.__init__(self)
        if comment:
            self._comment = comment
        self.sid = sid
        self.node = node
        self.cid = cid
        self.mag = mag
        self.xyz = xyz

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        node = integer(card, 2, 'node')
        cid = integer_or_blank(card, 3, 'cid', 0)
        mag = double(card, 4, 'mag')
        xyz = array([double_or_blank(card, 5, 'X1', 0.0),
                     double_or_blank(card, 6, 'X2', 0.0),
                     double_or_blank(card, 7, 'X3', 0.0)])
        assert len(card) <= 8, 'len(FORCE card) = %i\ncard=%s' % (len(card), card)
        return FORCE(sid, node, cid, mag, xyz, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        sid = data[0]
        node = data[1]
        cid = data[2]
        mag = data[3]
        xyz = array(data[4:7])
        assert len(xyz) == 3, 'xyz=%s' % str(xyz)
        return FORCE(sid, node, cid, mag, xyz, comment=comment)

    @property
    def node_id(self):
        if isinstance(self.node, int):
            return self.node
        return self.node_ref.nid

    def Cid(self):
        if isinstance(self.cid, integer_types):
            return self.cid
        return self.cid_ref.cid

    def F(self):
        return {self.node_id : self.mag * self.xyz}

    #def nodeID(self):
        #return self.node

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by FORCE sid=%s' % self.sid
        self.cid = model.Coord(self.cid, msg=msg)
        self.cid_ref = self.cid

    def uncross_reference(self):
        self.cid = self.Cid()
        del self.cid_ref

    def safe_cross_reference(self, model, debug=True):
        msg = ' which is required by FORCE sid=%s' % self.sid
        # try:
        self.cid = model.Coord(self.cid, msg=msg)
        self.cid_ref = self.cid
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
            msg = 'FORCE   %8i%8i%8s%8s%8s%8s%8s\n' % (
                self.sid, self.node,
                cids, print_float_8(self.mag), print_float_8(self.xyz[0]),
                print_float_8(self.xyz[1]), print_float_8(self.xyz[2]))
        else:
            cids = set_string16_blank_if_default(self.Cid(), 0)
            if is_double:
                msg = ('FORCE*  %16i%16i%16s%s\n'
                       '*       %16s%16s%16s\n') % (
                           self.sid, self.node,
                           cids, print_scientific_double(self.mag),
                           print_scientific_double(self.xyz[0]),
                           print_scientific_double(self.xyz[1]),
                           print_scientific_double(self.xyz[2]))
            else:
                msg = ('FORCE*  %16i%16i%16s%s\n'
                       '*       %16s%16s%16s\n') % (
                           self.sid, self.node,
                           cids, print_float_16(self.mag), print_float_16(self.xyz[0]),
                           print_float_16(self.xyz[1]), print_float_16(self.xyz[2]))
        return self.comment + msg


class FORCE1(Force):
    """
    Defines a static concentrated force at a grid point by specification of a
    magnitude and two grid points that determine the direction.
    """
    type = 'FORCE1'

    def __init__(self, sid, node, mag, g1, g2, comment=''):
        Force.__init__(self)
        if comment:
            self._comment = comment
        self.sid = sid
        self.node = node
        self.mag = mag
        self.g1 = g1
        self.g2 = g2

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        node = integer(card, 2, 'node')
        mag = double(card, 3, 'mag')
        g1 = integer(card, 4, 'g1')
        g2 = integer(card, 5, 'g2')
        assert len(card) == 6, 'len(FORCE1 card) = %i\ncard=%s' % (len(card), card)
        return FORCE1(sid, node, mag, g1, g2, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        sid = data[0]
        node = data[1]
        mag = data[2]
        g1 = data[3]
        g2 = data[4]
        return FORCE1(sid, node, mag, g1, g2, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by FORCE1 sid=%s' % self.sid
        self.node = model.Node(self.node, msg=msg)
        self.node_ref = self.node

        self.g1 = model.Node(self.g1, msg=msg)
        self.g1_ref = self.g1

        self.g2 = model.Node(self.g2, msg=msg)
        self.g2_ref = self.g2

        self.xyz = self.g2_ref.get_position() - self.g1_ref.get_position()
        self.normalize()

    def uncross_reference(self):
        self.node = self.node_id
        self.g1 = self.G1()
        self.g2 = self.G2()
        del self.node_ref, self.g1_ref, self.g2_ref

    def safe_cross_reference(self, model, debug=True):
        """
        .. todo:: cross reference and fix repr function
        """
        return self.cross_reference(model)
        #msg = ' which is required by FORCE1 sid=%s' % self.sid
        #self.node = model.Node(self.node, msg=msg)
        #self.node_ref = self.node
        #self.g1 = model.Node(self.g1, msg=msg)
        #self.g1_ref = self.g1
        #self.g2 = model.Node(self.g2, msg=msg)
        #self.g2_ref = self.g2
        #self.xyz = self.g2.get_position() - self.g1.get_position()
        #self.normalize()

    def G1(self):
        if isinstance(self.g1, (integer_types, float)):
            return self.g1
        return self.g1_ref.nid

    def G2(self):
        if isinstance(self.g2, (integer_types, float)):
            return self.g2
        return self.g2_ref.nid

    @property
    def node_ids(self):
        return [self.node_id, self.G1(), self.G2()]

    @property
    def node_id(self):
        if isinstance(self.node, integer_types):
            return self.node
        return self.node_ref.nid

    def raw_fields(self):
        list_fields = ['FORCE1', self.sid, self.node_id, self.mag, self.G1(), self.G2()]
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

    def __init__(self, sid, node, mag, g1, g2, g3, g4, comment=''):
        """
        +--------+-----+---+---+----+----+----+----+
        | FORCE2 | SID | G | F | G1 | G2 | G3 | G4 |
        +--------+-----+---+---+----+----+----+----+
        """
        Force.__init__(self)
        if comment:
            self._comment = comment
        self.sid = sid
        self.node = node
        self.mag = mag
        self.g1 = g1
        self.g2 = g2
        self.g3 = g3
        self.g4 = g4

    def validate(self):
        assert self.g1 is not None, self.g1
        assert self.g2 is not None, self.g2
        assert self.g3 is not None, self.g3
        assert self.g4 is not None, self.g4
        #print('self.node_ids = ', self.node_ids)

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        node = integer(card, 2, 'node')
        mag = double(card, 3, 'mag')
        g1 = integer(card, 4, 'g1')
        g2 = integer(card, 5, 'g2')
        g3 = integer(card, 6, 'g3')
        g4 = integer(card, 7, 'g4')
        assert len(card) == 8, 'len(FORCE2 card) = %i\ncard=%s' % (len(card), card)
        return FORCE2(sid, node, mag, g1, g2, g3, g4, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        sid = data[0]
        node = data[1]
        mag = data[2]
        g1 = data[3]
        g2 = data[4]
        g3 = data[5]
        g4 = data[6]
        return FORCE2(sid, node, mag, g1, g2, g3, g4, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by FORCE2 sid=%s' % self.sid
        self.node = model.Node(self.node, msg=msg)
        self.g1 = model.Node(self.g1, msg=msg)
        self.g2 = model.Node(self.g2, msg=msg)
        self.g3 = model.Node(self.g3, msg=msg)
        self.g4 = model.Node(self.g4, msg=msg)

        self.node_ref = self.node
        self.g1_ref = self.g1
        self.g2_ref = self.g2
        self.g3_ref = self.g3
        self.g4_ref = self.g4

        v12 = self.g2_ref.get_position() - self.g1_ref.get_position()
        v34 = self.g4_ref.get_position() - self.g3_ref.get_position()
        try:
            v12 /= norm(v12)
        except FloatingPointError:
            msg = 'v12=%s norm(v12)=%s\n' % (v12, norm(v12))
            msg += 'g1.get_position()=%s\n' % self.g1_ref.get_position()
            msg += 'g2.get_position()=%s' % self.g2_ref.get_position()
            raise FloatingPointError(msg)

        try:
            v34 /= norm(v34)
        except FloatingPointError:
            msg = 'v34=%s norm(v34)=%s\n' % (v34, norm(v34))
            msg += 'g3.get_position()=%s\n' % self.g3_ref.get_position()
            msg += 'g4.get_position()=%s' % self.g4_ref.get_position()
            raise FloatingPointError(msg)
        self.xyz = cross(v12, v34)
        self.normalize()

    def safe_cross_reference(self, model, debug=True):
        """
        .. todo:: cross reference and fix repr function
        """
        msg = ' which is required by FORCE2 sid=%s' % self.sid
        self.node = model.Node(self.node, msg=msg)
        self.g1 = model.Node(self.g1, msg=msg)
        self.g2 = model.Node(self.g2, msg=msg)
        self.g3 = model.Node(self.g3, msg=msg)
        self.g4 = model.Node(self.g4, msg=msg)

        self.node_ref = self.node
        self.g1_ref = self.g1
        self.g2_ref = self.g2
        self.g3_ref = self.g3
        self.g4_ref = self.g4

        v12 = self.g2_ref.get_position() - self.g1_ref.get_position()
        v34 = self.g4_ref.get_position() - self.g3_ref.get_position()
        try:
            v12 /= norm(v12)
        except FloatingPointError:
            msg = 'v12=%s norm(v12)=%s\n' % (v12, norm(v12))
            msg += 'g1.get_position()=%s\n' % self.g1_ref.get_position()
            msg += 'g2.get_position()=%s' % self.g2_ref.get_position()
            raise FloatingPointError(msg)

        try:
            v34 /= norm(v34)
        except FloatingPointError:
            msg = 'v34=%s norm(v34)=%s\n' % (v34, norm(v34))
            msg += 'g3.get_position()=%s\n' % self.g3_ref.get_position()
            msg += 'g4.get_position()=%s' % self.g4_ref.get_position()
            raise FloatingPointError(msg)
        self.xyz = cross(v12, v34)
        self.normalize()

    def uncross_reference(self):
        self.node = self.node_id
        self.g1 = self.G1()
        self.g2 = self.G2()
        self.g3 = self.G3()
        self.g4 = self.G4()
        del self.node_ref, self.g1_ref, self.g2_ref, self.g3_ref, self.g4_ref

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
        return self.g1_ref.nid

    def G2(self):
        if isinstance(self.g2, integer_types):
            return self.g2
        return self.g2_ref.nid

    def G3(self):
        if isinstance(self.g3, integer_types):
            return self.g3
        return self.g3_ref.nid

    def G4(self):
        if isinstance(self.g4, integer_types):
            return self.g4
        return self.g4_ref.nid

    @property
    def node_ids(self):
        return [self.node_id, self.G1(), self.G2(), self.G3(), self.G4()]

    def raw_fields(self):
        (node, g1, g2, g3, g4) = self._nodeIDs([self.node, self.g1, self.g2, self.g3, self.g4])
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

    def __init__(self, sid, node, cid, mag, xyz, comment=''):
        """
        Defines a static concentrated moment at a grid point by specifying a
        scale factor and a vector that determines the direction.::

        +--------+-----+---+-----+-----+-----+-----+-----+
        | MOMENT | SID | G | CID |  M  |  N1 |  N2 |  N3 |
        +--------+-----+---+-----+-----+-----+-----+-----+
        | MOMENT |  2  | 5 |  6  | 2.9 | 0.0 | 1.0 | 0.0 |
        +--------+-----+---+-----+-----+-----+-----+-----+
        """
        Moment.__init__(self)
        if comment:
            self._comment = comment
        self.sid = sid
        self.node = node
        self.cid = cid
        self.mag = mag
        self.xyz = xyz

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        node = integer(card, 2, 'node')
        cid = integer_or_blank(card, 3, 'cid', 0)
        mag = double(card, 4, 'mag')

        xyz = array([double_or_blank(card, 5, 'X1', 0.0),
                     double_or_blank(card, 6, 'X2', 0.0),
                     double_or_blank(card, 7, 'X3', 0.0)])
        assert len(card) <= 8, 'len(MOMENT card) = %i\ncard=%s' % (len(card), card)
        assert len(xyz) == 3, 'xyz=%s' % str(xyz)
        return MOMENT(sid, node, cid, mag, xyz, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        sid = data[0]
        node = data[1]
        cid = data[2]
        mag = data[3]
        xyz = data[4:7]
        assert len(xyz) == 3, 'xyz=%s' % str(xyz)
        return MOMENT(sid, node, cid, mag, xyz, comment=comment)

    def Cid(self):
        if isinstance(self.cid, integer_types):
            return self.cid
        return self.cid_ref.cid

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by MOMENT sid=%s' % self.sid
        self.node = model.Node(self.node, msg=msg)
        self.node_ref = self.node
        self.cid = model.Coord(self.cid, msg=msg)
        self.cid_ref = self.cid

    def safe_cross_reference(self, model, debug=True):
        msg = ' which is required by MOMENT sid=%s' % self.sid
        self.node = model.Node(self.node, msg=msg)
        self.node_ref = self.node
        self.cid = model.Coord(self.cid, msg=msg)
        self.cid_ref = self.cid

    def uncross_reference(self):
        self.node = self.node_id
        self.cid = self.Cid()
        del self.node_ref, self.cid_ref

    @property
    def node_ids(self):
        """all the nodes referenced by the load"""
        return [self.node_id]

    @property
    def node_id(self):
        if isinstance(self.node, integer_types):
            return self.node
        return self.node_ref.nid

    def raw_fields(self):
        list_fields = ['MOMENT', self.sid, self.node, self.Cid(),
                       self.mag] + list(self.xyz)
        return list_fields

    def repr_fields(self):
        cid = set_blank_if_default(self.Cid(), 0)
        list_fields = ['MOMENT', self.sid, self.node_id, cid,
                       self.mag] + list(self.xyz)
        return list_fields

    def write_card(self, size=8, is_double=False):
        if size == 8:
            scid = set_string8_blank_if_default(self.Cid(), 0)
            msg = 'MOMENT  %8i%8i%8s%8s%8s%8s%8s\n' % (
                self.sid, self.node_id,
                scid, print_float_8(self.mag), print_float_8(self.xyz[0]),
                print_float_8(self.xyz[1]), print_float_8(self.xyz[2]))
        else:
            scid = set_string16_blank_if_default(self.Cid(), 0)
            if is_double:
                msg = ('MOMENT* %16i%16i%16s%s\n'
                       '*       %16s%16s%16s\n') % (
                           self.sid, self.node_id,
                           scid, print_scientific_double(self.mag),
                           print_scientific_double(self.xyz[0]),
                           print_scientific_double(self.xyz[1]),
                           print_scientific_double(self.xyz[2]))
            else:
                msg = ('MOMENT* %16i%16i%16s%s\n'
                       '*       %16s%16s%16s\n') % (
                           self.sid, self.node_id,
                           scid, print_float_16(self.mag), print_float_16(self.xyz[0]),
                           print_float_16(self.xyz[1]), print_float_16(self.xyz[2]))
        return self.comment + msg


class MOMENT1(Moment):
    type = 'MOMENT1'

    def __init__(self, sid, node, mag, g1, g2, comment=''):
        """
        Defines a static concentrated moment at a grid point by specifying a
        magnitude and two grid points that determine the direction.::

        +---------+-----+---+---+----+----+
        | MOMENT1 | SID | G | M | G1 | G2 |
        +---------+-----+---+---+----+----+
        """
        Moment.__init__(self)
        if comment:
            self._comment = comment
        self.xyz = None
        self.sid = sid
        self.node = node
        self.mag = mag
        self.g1 = g1
        self.g2 = g2

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        node = integer(card, 2, 'node')
        mag = double(card, 3, 'mag')
        g1 = integer(card, 4, 'g1')
        g2 = integer(card, 5, 'g2')
        assert len(card) == 6, 'len(MOMENT1 card) = %i\ncard=%s' % (len(card), card)
        return MOMENT1(sid, node, mag, g1, g2, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        sid = data[0]
        node = data[1]
        mag = data[2]
        g1 = data[3]
        g2 = data[4]
        return MOMENT1(sid, node, mag, g1, g2, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by MOMENT1 sid=%s' % self.sid
        self.node = model.Node(self.node, msg=msg)
        self.g1 = model.Node(self.g1, msg=msg)
        self.g2 = model.Node(self.g2, msg=msg)

        self.node_ref = self.node
        self.g1_ref = self.g1
        self.g2_ref = self.g2

        self.xyz = self.g2_ref.get_position() - self.g1_ref.get_position()
        self.normalize()

    def uncross_reference(self):
        self.node = self.node_id
        self.g1 = self.G1()
        self.g2 = self.G2()
        del self.node_ref, self.g1_ref, self.g2_ref

    def safe_cross_reference(self, model, debug=True):
        """
        .. todo:: cross reference and fix repr function
        """
        msg = ' which is required by MOMENT1 sid=%s' % self.sid
        self.node = model.Node(self.node, msg=msg)
        self.g1 = model.Node(self.g1, msg=msg)
        self.g2 = model.Node(self.g2, msg=msg)

        self.node_ref = self.node
        self.g1_ref = self.g1
        self.g2_ref = self.g2

        self.xyz = self.g2_ref.get_position() - self.g1_ref.get_position()
        self.normalize()

    @property
    def node_ids(self):
        """all the nodes referenced by the load"""
        return [self.node_id, self.G1(), self.G2()]

    @property
    def node_id(self):
        if isinstance(self.node, integer_types):
            return self.node
        return self.node_ref.nid

    def G1(self):
        if isinstance(self.g1, integer_types):
            return self.g1
        return self.g1_ref.nid

    def G2(self):
        if isinstance(self.g2, integer_types):
            return self.g2
        return self.g2_ref.nid

    def raw_fields(self):
        list_fields = ['MOMENT1', self.sid, self.node_id, self.mag, self.G1(), self.G2()]
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

    def __init__(self, sid, node, mag, g1, g2, g3, g4, xyz=None, comment=''):
        """
        Defines a static concentrated moment at a grid point by specification
        of a magnitude and four grid points that determine the direction.::

        +---------+-----+---+---+----+----+----+----+
        | MOMENT2 | SID | G | M | G1 | G2 | G3 | G4 |
        +---------+-----+---+---+----+----+----+----+
        """
        Moment.__init__(self)
        if comment:
            self._comment = comment
        self.sid = sid
        self.node = node
        self.mag = mag
        self.g1 = g1
        self.g2 = g2
        self.g3 = g3
        self.g4 = g4
        self.xyz = xyz

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        node = integer(card, 2, 'node')
        mag = double(card, 3, 'mag')
        g1 = integer(card, 4, 'g1')
        g2 = integer(card, 5, 'g2')
        g3 = integer(card, 6, 'g3')
        g4 = integer(card, 7, 'g4')
        xyz = None
        assert len(card) <= 8, 'len(MOMENT2 card) = %i\ncard=%s' % (len(card), card)
        return MOMENT2(sid, node, mag, g1, g2, g3, g4, xyz, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        sid = data[0]
        node = data[1]
        mag = data[2]
        g1 = data[3]
        g2 = data[4]
        g3 = data[5]
        g4 = data[6]
        assert len(data) == 7, data
        #xyz = array(data[7:10])
        #print('data =', data)
        #assert len(xyz) == 3, 'xyz=%s' % str(xyz)
        xyz = None
        return MOMENT2(sid, node, mag, g1, g2, g3, g4, xyz, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by MOMENT2 sid=%s' % self.sid
        self.node = model.Node(self.node, msg=msg)
        self.g1 = model.Node(self.g1, msg=msg)
        self.g2 = model.Node(self.g2, msg=msg)
        self.g3 = model.Node(self.g3, msg=msg)
        self.g4 = model.Node(self.g4, msg=msg)

        self.node_ref = self.node
        self.g1_ref = self.g1
        self.g2_ref = self.g2
        self.g3_ref = self.g3
        self.g4_ref = self.g4

        v12 = self.g2_ref.get_position() - self.g1_ref.get_position()
        v34 = self.g4_ref.get_position() - self.g3_ref.get_position()
        v12 = v12 / norm(v12)
        v34 = v34 / norm(v34)
        self.xyz = cross(v12, v34)

    def uncross_reference(self):
        self.node = self.NodeID()
        self.g1 = self.G1()
        self.g2 = self.G2()
        self.g3 = self.G3()
        self.g4 = self.G4()
        del self.node_ref, self.g1_ref, self.g2_ref, self.g3_ref, self.g4_ref

    def safe_cross_reference(self, model, debug=True):
        """
        .. todo:: cross reference and fix repr function
        """
        msg = ' which is required by MOMENT2 sid=%s' % self.sid
        self.node = model.Node(self.node, msg=msg)
        self.g1 = model.Node(self.g1, msg=msg)
        self.g2 = model.Node(self.g2, msg=msg)
        self.g3 = model.Node(self.g3, msg=msg)
        self.g4 = model.Node(self.g4, msg=msg)

        self.node_ref = self.node
        self.g1_ref = self.g1
        self.g2_ref = self.g2
        self.g3_ref = self.g3
        self.g4_ref = self.g4

        v12 = self.g2_ref.get_position() - self.g1_ref.get_position()
        v34 = self.g4_ref.get_position() - self.g3_ref.get_position()
        v12 = v12 / norm(v12)
        v34 = v34 / norm(v34)
        self.xyz = cross(v12, v34)

    def NodeID(self):
        if isinstance(self.node, integer_types):
            return self.node
        return self.node_ref.nid

    def G1(self):
        if isinstance(self.g1, integer_types):
            return self.g1
        return self.g1_ref.nid

    def G2(self):
        if isinstance(self.g2, integer_types):
            return self.g2
        return self.g2_ref.nid

    def G3(self):
        if isinstance(self.g3, integer_types):
            return self.g3
        return self.g3_ref.nid

    def G4(self):
        if isinstance(self.g4, integer_types):
            return self.g4
        return self.g4_ref.nid

    @property
    def node_ids(self):
        """all the nodes referenced by the load"""
        return [self.NodeID(), self.G1(), self.G2(), self.G3(), self.G4()]

    def raw_fields(self):
        list_fields = [
            'MOMENT2', self.sid, self.NodeID(), self.mag,
            self.G1(), self.G2(), self.G3(), self.G4()]
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

    def __init__(self, sid, cid, normal, entity, entity_id, method,
                 load_magnitudes, comment=''):
        Load.__init__(self)
        if comment:
            self._comment = comment
        self.sid = sid
        self.cid = cid
        self.normal = normal
        self.entity = entity
        self.entity_id = entity_id
        self.method = method
        self.load_magnitudes = load_magnitudes

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        cid = integer_or_blank(card, 2, 'cid', 0)
        normal = array([
            double_or_blank(card, 3, 'N1', 0.),
            double_or_blank(card, 4, 'N2', 0.),
            double_or_blank(card, 5, 'N3', 1.),
        ])
        entity = string(card, 6, 'entity')
        entity_id = integer(card, 7, 'entity_id')
        method = string(card, 8, 'method')

        load_magnitudes = []
        for i in range(9, len(card)):
            ifield = i - 8
            load_mag = integer_or_double(card, i, 'load_magnitude_%s' % ifield)
            load_magnitudes.append(load_mag)
        return GMLOAD(sid, cid, normal, entity, entity_id, method,
                      load_magnitudes, comment=comment)

    #def DEquation(self):
        #if isinstance(self.dequation, int):
            #return self.dequation
        #return self.dequation.equation_id

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by GMLOAD sid=%s' % self.sid
        self.cid = model.Coord(self.Cid(), msg=msg)
        self.cid_ref = self.cid
        #self.node = model.Node(self.node, msg=msg)
        #self.g1 = model.Node(self.g1, msg=msg)
        #self.g2 = model.Node(self.g2, msg=msg)
        #self.xyz = self.g2.get_position() - self.g1.get_position()
        #self.normalize()

    def safe_cross_reference(self, model):
        return self.cross_reference(model)

    def uncross_reference(self):
        self.cid = self.Cid()
        del self.cid_ref

    def Cid(self):
        if isinstance(self.cid, integer_types):
            return self.cid
        return self.cid_ref.cid

    #def G1(self):
        #if isinstance(self.g1, (integer_types, float)):
            #return self.g1
        #return self.g1_ref.nid

    #def G2(self):
        #if isinstance(self.g2, (integer_types, float)):
            #return self.g2
        #return self.g2_ref.nid

    #def NodeID(self):
        #if isinstance(self.node, integer_types):
            #return self.node
        #return self.node_ref.nid

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

    def __init__(self, sid, p, nodes, comment=''):
        if comment:
            self._comment = comment
        self.sid = sid
        self.p = p
        self.nodes = nodes
        assert len(self.nodes) in [3, 4], 'nodes=%s' % self.nodes

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        p = double(card, 2, 'p')
        nodes = [integer(card, 3, 'n1'),
                 integer(card, 4, 'n2'),
                 integer(card, 5, 'n3')]
        n4 = integer_or_blank(card, 6, 'n4', 0)
        if n4:
            nodes.append(n4)
        assert len(card) <= 7, 'len(PLOAD card) = %i\ncard=%s' % (len(card), card)
        return PLOAD(sid, p, nodes, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        sid = data[0]
        p = data[1]
        nodes = data[2:]
        return PLOAD(sid, p, nodes, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        pass

    def safe_cross_reference(self, model):
        return self.cross_reference(model)

    def uncross_reference(self):
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
    valid_types = ['FX', 'FY', 'FZ', 'FXE', 'FYE', 'FZE',
                   'MX', 'MY', 'MZ', 'MXE', 'MYE', 'MZE']
    valid_scales = ['LE', 'FR', 'LEPR', 'FRPR'] # LE: length-based; FR: fractional; PR:projected

    def __init__(self, sid, eid, Type, scale, x1, p1, x2, p2, comment=''):
        if comment:
            self._comment = comment
        self.sid = sid
        self.eid = eid
        self.Type = Type
        self.scale = scale
        self.x1 = x1
        self.p1 = p1
        self.x2 = x2
        self.p2 = p2
        self._validate_input()

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        eid = integer(card, 2, 'eid')
        Type = string(card, 3, 'Type ("%s")' % '",  "'.join(cls.valid_types))
        scale = string(card, 4, 'scale ("%s")' % '", "'.join(cls.valid_scales))
        x1 = double(card, 5, 'x1')
        p1 = double(card, 6, 'p1')
        x2 = double_or_blank(card, 7, 'x2', x1)
        p2 = double_or_blank(card, 8, 'p2', p1)
        assert len(card) <= 9, 'len(PLOAD1 card) = %i\ncard=%s' % (len(card), card)
        return PLOAD1(sid, eid, Type, scale, x1, p1, x2, p2, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        sid = data[0]
        eid = data[1]
        Type = data[2]
        scale = data[3]
        x1 = data[4]
        p1 = data[5]
        x2 = data[6]
        p2 = data[7]
        Type = cls.valid_types[Type - 1]
        scale = cls.valid_scales[scale - 1]
        return PLOAD1(sid, eid, Type, scale, x1, p1, x2, p2, comment=comment)

    def _validate_input(self):
        if self.Type not in self.valid_types:
            msg = '%s is an invalid type on the PLOAD1 card; valid_types=[%s]' % (
                self.Type, ', '.join(self.valid_types).rstrip(', '))
            raise RuntimeError(msg)
        if self.scale not in self.valid_scales:
            msg = '%s is an invalid scale on the PLOAD1 card; valid_scales=[%s]' % (
                self.scale, ', '.join(self.valid_scales).rstrip(', '))
            raise RuntimeError(msg)

        assert 0.0 <= self.x1 <= self.x2, '0.0 <= x1 <= x2 -> x1=%s x2=%s' % (self.x1, self.x2)
        if self.scale in ['FR', 'FRPR']:
            assert self.x1 <= 1.0, 'x1=%r' % self.x1
            assert self.x2 <= 1.0, 'x2=%r' % self.x2
        assert self.scale in self.valid_scales, '%s is an invalid scale on the PLOAD1 card' % (self.scale)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by PLOAD1 sid=%s' % self.sid
        self.eid = model.Element(self.eid, msg=msg)
        self.eid_ref = self.eid

    def safe_cross_reference(self, model):
        return self.cross_reference(model)

    def uncross_reference(self):
        self.eid = self.Eid()
        del self.eid_ref

    def transformLoad(self):
        self.deprecated('transformLoad()', 'transform_load()', '0.8')

    def transform_load(self):
        p1 = self.eid_ref.ga_ref.get_position()
        p2 = self.eid_ref.gb_ref.get_position()

        g0 = self.eid_ref.g0_vector
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
        force_loads = {}  # spc enforced displacement (e.g. FORCE=0)
        moment_loads = {}
        force_constraints = {}
        moment_constraints = {}
        gravity_loads = []

        types_found = set()
        (scale_factors, loads) = self.get_reduced_loads()

        for scale_factor, load in zip(scale_factors, loads):
            out = load.transformLoad()
            types_found.add(load.__class__.__name__)

            if isinstance(load, PLOAD1): # CBAR/CBEAM
                element = load.eid
                (ga, gb) = element.node_ids
                load_type = load.Type

                scale = load.scale
                eType = element.type

                if load_type in ['FX', 'FY', 'FZ']:
                    p1 = element.ga_ref.get_position()
                    p2 = element.gb_ref.get_position()
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
                        msg = 'p1=%r must be equal to p2=%r for x1=x2=%r'  % (
                            self.p1, self.p2, self.x1)
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
                        force_loads[ga] = Fv
                        force_loads[gb] = Fv
                        moment_loads[ga] = Mv
                        moment_loads[gb] = Mv
                    elif load_type in ['MX', 'MY', 'MZ']:
                        # these are really moments
                        Mv *= F
                        moment_loads[ga] = Mv
                        moment_loads[gb] = Mv
                    else:
                        raise NotImplementedError(load_type)
                else:
                    raise NotImplementedError(eType)
            else:
                msg = '%s not supported' % (load.__class__.__name__)
                raise NotImplementedError(msg)
        return (types_found, force_loads, moment_loads, force_constraints,
                moment_constraints, gravity_loads)

    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

    def get_loads(self):
        return [self]

    def Eid(self):
        if isinstance(self.eid, integer_types):
            return self.eid
        return self.eid_ref.eid

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

    def __init__(self, sid, pressure, eids, comment=''):
        if comment:
            self._comment = comment
        self.sid = sid
        self.pressure = pressure
        self.eids = eids

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        pressure = double(card, 2, 'p')

        if integer_string_or_blank(card, 4, 'THRU') == 'THRU':
            e1 = integer(card, 3, 'Element1')
            e2 = integer(card, 5, 'Element1')
            eids = [i for i in range(e1, e2 + 1)]
            assert len(card) == 6, 'len(PLOAD2 card) = %i\ncard=%s' % (len(card), card)
        else:
            eids = fields(integer, card, 'eid', i=3, j=len(card))
        return PLOAD2(sid, pressure, eids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        sid = data[0]
        pressure = data[1]
        eids = list(data[2:])
        return PLOAD2(sid, pressure, eids, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by PLOAD2 sid=%s' % self.sid
        self.eids = model.Elements(self.eids, msg=msg)
        self.eids_ref = self.eids

    def safe_cross_reference(self, model):
        return self.cross_reference(model)

    def uncross_reference(self):
        self.eids = self.element_ids
        del self.eids_ref

    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

    @property
    def element_ids(self):
        if isinstance(self.eids[0], int):
            eids = self.eids
        else:
            eids = [elem.Eid() for elem in self.eids]
        return eids

    def get_loads(self):
        return [self]

    def raw_fields(self):
        list_fields = ['PLOAD2', self.sid, self.pressure]
        eids = self.element_ids
        if len(eids) == 1:
            list_fields += eids
        else:
            eids.sort()
            delta_eid = eids[-1] - eids[0] + 1
            if delta_eid != len(eids):
                msg = 'eids=%s len(eids)=%s delta_eid=%s must be continuous' % (
                    eids, len(eids), delta_eid)
                raise RuntimeError(msg)
            #list_fields += eids
            list_fields += [eids[0], 'THRU', eids[-1]]
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

    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

    def transformLoad(self):
        self.deprecated('transformLoad()', 'transform_load()', '0.8')
        return self.transform_load()

    def __init__(self, sid, eids, pressures, g1, g34, cid, NVector, sorl, ldir, comment=''):
        if comment:
            self._comment = comment
        self.sid = sid

        # these can be greater than 1 if it's a shell (not a solid)
        self.eids = eids
        self.pressures = pressures

        #: used for solid element only
        self.g1 = g1
        #: g3/g4 - different depending on CHEXA/CPENTA or CTETRA
        self.g34 = g34

        #: Coordinate system identification number. See Remark 2.
        #: (Integer >= 0;Default=0)
        self.cid = cid
        self.NVector = NVector
        self.sorl = sorl
        self.ldir = ldir
        self._validate_input()

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        eid = integer(card, 2, 'eid')
        p1 = double_or_blank(card, 3, 'p1', 0.0)
        pressures = [
            p1,
            double_or_blank(card, 4, 'p2', p1),
            double_or_blank(card, 5, 'p3', p1),
            double_or_blank(card, 6, 'p4', p1)]

        eids = [eid]
        if (integer_string_or_blank(card, 7, 'g1/THRU') == 'THRU' and
            integer_or_blank(card, 8, 'eid2')):
            eid2 = integer(card, 8, 'eid2')
            if eid2:
                eids = list(unique(
                    expand_thru([eid, 'THRU', eid2], set_fields=False, sort_fields=False)
                ))
            g1 = None
            g34 = None
        else:
            eids = [eid]
            g1 = integer_or_blank(card, 7, 'g1')
            g34 = integer_or_blank(card, 8, 'g34')

        cid = integer_or_blank(card, 9, 'cid', 0)
        NVector = array([double_or_blank(card, 10, 'N1', 0.0),
                         double_or_blank(card, 11, 'N2', 0.0),
                         double_or_blank(card, 12, 'N3', 0.0)])
        sorl = string_or_blank(card, 13, 'sorl', 'SURF')
        ldir = string_or_blank(card, 14, 'ldir', 'NORM')
        assert len(card) <= 15, 'len(PLOAD4 card) = %i\ncard=%s' % (len(card), card)
        return PLOAD4(sid, eids, pressures, g1, g34, cid, NVector, sorl, ldir, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        sid = data[0]
        eid = data[1]
        pressures = data[2]

        g1 = data[3]
        g34 = data[4]
        if g1 == 0:
            g1 = None
        if g34 == 0:
            g34 = None
        cid = data[5]
        NVector = data[6]

        sorl = data[7]
        #self.ldir = data[8]
        #assert len(data)==8

        eids = [eid]
        sorl = 'SURF'
        ldir = 'NORM'
        return PLOAD4(sid, eids, pressures, g1, g34, cid, NVector, sorl, ldir, comment=comment)

    def _validate_input(self):
        if self.sorl not in ['SURF', 'LINE']:
            raise RuntimeError(self.sorl)
        if self.ldir not in ['LINE', 'X', 'Y', 'Z', 'TANG', 'NORM']:
            raise RuntimeError(self.ldir)
        assert self.g1 != 0, str(self)
        assert self.g34 != 0, str(self)

    def get_loads(self):
        return [self]

    def transform_load(self):
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
            msg = 'Only one load may be defined on each PLOAD4.  nLoads=%s\n%s' % (
                len(self.eids), str(self))
            raise RuntimeError(msg)

        elem = self.eids_ref[0]
        if self.g1 and self.g34:  # solid elements
            nid = self.g1_ref.nid
            nid_opposite = self.g34_ref.nid
            (face_node_ids, area) = elem.getFaceNodesAndArea(self, nid, nid_opposite)
        else:
            face_node_ids = elem.node_ids
            area = elem.Area()
        n = len(face_node_ids)

        elem = self.eids_ref[0]
        vector = array(elem.Normal())
        vectors = []
        for (nid, p) in zip(face_node_ids, self.pressures):
            #: .. warning:: only supports normal pressures
            vectors.append(vector * p * area / n)  # Force_i

        is_load = None
        return (is_load, face_node_ids, vectors)

    def Cid(self):
        if isinstance(self.cid, integer_types):
            return self.cid
        return self.cid_ref.cid

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by PLOAD4 sid=%s' % self.sid
        #self.eid = model.Element(self.eid, msg=msg)
        self.cid = model.Coord(self.cid, msg=msg)
        #self.eid_ref = self.eid
        self.cid_ref = self.cid
        if self.g1 is not None:
            self.g1 = model.Node(self.g1, msg=msg + '; g1')
            self.g1_ref = self.g1
        if self.g34 is not None:
            self.g34 = model.Node(self.g34, msg=msg + '; g34')
            self.g34_ref = self.g34
        if self.eids:
            self.eids = model.Elements(self.eids, msg=msg)
            self.eids_ref = self.eids

    def safe_cross_reference(self, model, debug=True):
        msg = ' which is required by PLOAD4 sid=%s' % self.sid
        #self.eid = model.Element(self.eid, msg=msg)
        self.cid = model.Coord(self.cid, msg=msg)
        #self.eid_ref = self.eid
        self.cid_ref = self.cid
        if self.g1 is not None:
            self.g1 = model.Node(self.g1, msg=msg)
            self.g1_ref = self.g1
        if self.g34 is not None:
            self.g34 = model.Node(self.g34, msg=msg)
            self.g34_ref = self.g34
        #if self.eids:
        self.eids = model.Elements(self.eids, msg=msg)
        self.eids_ref = self.eids

    def uncross_reference(self):
        #self.eid = self.Eid(self.eid)
        self.cid = self.Cid()
        if self.g1 is not None:
            self.g1 = self.G1()
            del self.g1_ref
        if self.g34 is not None:
            self.g34 = self.G34()
            del self.g34_ref
        self.eids = self.element_ids
        del self.cid_ref, self.eids_ref

    def G1(self):
        if isinstance(self.g1, (integer_types)):
            return self.g1
        return self.g1_ref.nid

    def G34(self):
        if isinstance(self.g34, (integer_types)):
            return self.g34
        return self.g34_ref.nid

    def Eid(self, eid):
        if isinstance(eid, (integer_types)):
            return eid
        return eid.eid

    @property
    def node_ids(self):
        node_ids = [None, None]
        if isinstance(self.g1, integer_types):
            node_ids[0] = self.g1
        elif self.g1 is not None:
            node_ids[0] = self.g1_ref.nid

        if isinstance(self.g34, integer_types):
            node_ids[1] = self.g34
        elif self.g34 is not None:
            node_ids[1] = self.g34_ref.nid
        return node_ids

    def get_element_ids(self, eid=None):
        if eid:
            eidi = self.Eid(eid)
            #assert len(str(eidi)) < 20, eidi
            return eidi
        eids = []
        for element in self.eids:
            eidi = self.Eid(element)
            #assert len(str(eidi)) < 20, eidi
            eids.append(eidi)
        return eids

    @property
    def element_ids(self):
        return self.get_element_ids()

    def raw_fields(self):
        eids = self.element_ids
        #print('eids = ', eids)
        eid = eids[0]
        #print('eid = ', eid)
        cid = set_blank_if_default(self.Cid(), 0)
        sorl = set_blank_if_default(self.sorl, 'SURF')
        ldir = set_blank_if_default(self.ldir, 'NORM')
        p1 = self.pressures[0]
        p2 = set_blank_if_default(self.pressures[1], p1)
        p3 = set_blank_if_default(self.pressures[2], p1)
        p4 = set_blank_if_default(self.pressures[3], p1)
        list_fields = ['PLOAD4', self.sid, eid, self.pressures[0], p2, p3, p4]

        if self.g1 is not None:  # is it a SOLID element
            node_ids = self.node_ids
            #node_ids = self.node_ids([self.g1, self.g34])
            list_fields += node_ids
        else:
            if len(eids) > 1:
                #print("self.eids = %s" %(self.eids))
                try:
                    list_fields.append('THRU')
                    eidi = eids[-1]
                    #print('eidi = ', eidi)
                except:
                    print("g1  = %s" % self.g1)
                    print("g34 = %s" % self.g34)
                    print("self.eids = %s" % self.eids)
                    raise
                list_fields.append(eidi)
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

    def __init__(self, sid, eid, pa, pb, ga, gb, theta=0., comment=''):
        if comment:
            self._comment = comment
        self.sid = sid
        self.eid = eid
        self.pa = pa
        self.pb = pb
        self.ga = ga
        self.gb = gb
        self.theta = theta

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        eid = integer(card, 2, 'eid')
        pa = double(card, 3, 'pa')
        pb = double_or_blank(card, 4, 'pb', pa)
        ga = integer(card, 5, 'ga')
        gb = integer(card, 6, 'gb')
        theta = double_or_blank(card, 7, 'theta', 0.)
        assert len(card) <= 8, 'len(PLOADX1 card) = %i\ncard=%s' % (len(card), card)
        return PLOADX1(sid, eid, pa, pb, ga, gb, theta, comment=comment)

    #@classmethod
    #def add_op2_data(cls, data, comment=''):
        #sid = data[0]
        #print("PLOADX1 = %s" % data)
        #raise NotImplementedError(data)
        #return PLOADX1(sid, eid, pa, pb, ga, gb, theta, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by PLOADX1 lid=%s' % self.sid
        self.eid = model.Element(self.eid, msg=msg)
        self.ga = model.Node(self.ga, msg=msg)
        self.gb = model.Node(self.gb, msg=msg)

        self.eid_ref = self.eid
        self.ga_ref = self.ga
        self.gb_ref = self.gb

    def safe_cross_reference(self, model):
        return self.cross_reference(model)

    def uncross_reference(self):
        self.eid = self.Eid()
        self.ga = self.Ga()
        self.gb = self.Gb()
        del self.eid_ref, self.ga_ref, self.gb_ref

    def Eid(self):
        if isinstance(self.eid, (integer_types)):
            return self.eid
        return self.eid_ref.eid

    def Ga(self):
        if isinstance(self.ga, (integer_types)):
            return self.ga
        return self.ga_ref.nid

    def Gb(self):
        if isinstance(self.gb, (integer_types)):
            return self.gb
        return self.gb_ref.nid

    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

    def get_loads(self):
        return [self]

    def raw_fields(self):
        list_fields = [
            'PLOADX1', self.sid, self.Eid(), self.pa, self.pb,
            self.Ga(), self.Gb(), self.theta]
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
