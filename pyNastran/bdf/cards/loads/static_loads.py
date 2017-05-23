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

import numpy as np
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
        """
        Creates a LOAD card

        Parameters
        ----------
        sid : int
            load id
        scale : float
            overall scale factor
        scale_factors : List[float]
            individual scale factors (corresponds to load_ids)
        load_ids : List[int]
            individual load_ids (corresponds to scale_factors)
        comment : str; default=''
            a comment for the card
        """
        LoadCombination.__init__(self, sid, scale, scale_factors, load_ids,
                                 comment=comment)

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
                        load_types += [load.type] + load.get_load_types()
                elif isinstance(load, (Force, Moment, PLOAD4, GRAV)):
                    load_types += [load.type]
                else:
                    raise NotImplementedError(load)

        load_types = list(set(load_types))
        #print("load_types = ", load_types)
        return load_types

    def get_reduced_loads(self, resolve_load_card=False, filter_zero_scale_factors=False):
        """
        Get all load objects in a simplified form, which means all
        scale factors are already applied and only base objects
        (no LOAD cards) will be returned.

        Parameters
        ----------
        resolve_load_card : bool; default=False
            Nastran requires that LOAD cards do not reference other load cards
            This feature can be enabled.
        filter_zero_scale_factors : bool; default=False
            Nastran does not filter loads with a 0.0 scale factor.  So, if you
            have a 0.0 load, but are missing load ids, Nastran will throw a
            fatal error.

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
            if scale == 0.0 and filter_zero_scale_factors:
                continue

            for load in loads_pack:
                if simple_loads:
                    loads.append(load)
                    scale_factors.append(scale) # local
                elif isinstance(load, LOAD):
                    if not resolve_load_card:
                        msg = 'A LOAD card cannot reference another LOAD card\n'
                        msg += 'current:\n%s\n' % str(self)
                        msg += 'new:\n%s' % str(load)
                        raise RuntimeError(msg)
                    load_data = load.get_reduced_loads(
                        resolve_load_card=True,
                        filter_zero_scale_factors=filter_zero_scale_factors)
                    (reduced_scale_factors, reduced_loads) = load_data

                    loads += reduced_loads
                    scale_factors += [scale * j_scale
                                      for j_scale in reduced_scale_factors]
                else:
                    msg = ('%s isnt supported in get_reduced_loads method'
                           % load.__class__.__name__)
                    raise NotImplementedError(msg)
        return (scale_factors, loads)

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
    |  1   |  2  |  3  |   4  |  5  |  6  |   7  |  8  |
    +======+=====+=====+======+=====+=====+======+=====+
    | GRAV | SID | CID |  A   | N1  | N2  |  N3  |  MB |
    +------+-----+-----+------+-----+-----+------+-----+
    | GRAV | 1   | 3   | 32.2 | 0.0 | 0.0 | -1.0 |     |
    +------+-----+-----+------+-----+-----+------+-----+
    """
    type = 'GRAV'

    def __init__(self, sid, scale, N, cid=0, mb=0, comment=''):
        """
        Creates an GRAV card

        Parameters
        ----------
        sid : int
            load id
        scale : float
            scale factor for load
        N : (3, ) float ndarray
            the acceleration vector in the cid frame
        cid : int; default=0
            the coordinate system for the load
        mb : int; default=0
            ???
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment

        #: Set identification number
        self.sid = sid

        #: Coordinate system identification number.
        self.cid = cid

        #: scale factor
        self.scale = scale

        #: Acceleration vector components measured in coordinate system CID
        self.N = np.asarray(N)

        #: Indicates whether the CID coordinate system is defined in the
        #: main Bulk Data Section (MB = -1) or the partitioned superelement
        #: Bulk Data Section (MB = 0). Coordinate systems referenced in the
        #: main Bulk Data Section are considered stationary with respect to
        #: the assembly basic coordinate system. See Remark 10.
        #: (Integer; Default = 0)
        self.mb = mb

        assert not allclose(max(abs(self.N)), 0.), ('GRAV N is a zero vector, '
                                                    'N=%s' % str(self.N))

    def validate(self):
        if not isinstance(self.scale, float):
            msg = 'scale=%s type=%s' % (self.scale, type(self.scale))
            raise TypeError(msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a GRAV card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        cid = integer_or_blank(card, 2, 'cid', 0)
        scale = double(card, 3, 'scale')
        N = array([double_or_blank(card, 4, 'N1', 0.0),
                   double_or_blank(card, 5, 'N2', 0.0),
                   double_or_blank(card, 6, 'N3', 0.0)])
        mb = integer_or_blank(card, 7, 'mb', 0)
        assert len(card) <= 8, 'len(GRAV card) = %i\ncard=%s' % (len(card), card)
        return GRAV(sid, scale, N, cid=cid, mb=mb, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a GRAV card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        sid = data[0]
        cid = data[1]
        a = data[2]
        N = array(data[3:6])
        mb = data[6]
        scale = 1.
        assert len(data) == 7
        return GRAV(sid, scale, N, cid=cid, mb=mb, comment=comment)

    def get_loads(self):
        return [self]

    def transform_load(self):
        g = self.GravityVector()
        g2 = self.cid_ref.transform_node_to_global(g)
        return g2

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

    def __init__(self, sid, N, direction, locs, vals, cid=0, comment=''):
        """
        Creates an ACCEL card

        Parameters
        ----------
        sid : int
            load id
        N : (3, ) float ndarray
            the acceleration vector in the cid frame
        direction : str
            Component direction of acceleration variation
            {X, Y, Z}
        locs : ???
            ???
        vals : ???
            ???
        cid : int; default=0
            the coordinate system for the load
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        #: Load set identification number (Integer>0)
        self.sid = sid
        #: Coordinate system identification number. (Integer>0: Default=0)
        self.cid = cid

        #: Components of the acceleration vector measured in coordinate system
        #: CID. (Real; at least one Ni != 0)
        self.N = np.asarray(N, dtype='float64')

        #: Component direction of acceleration variation. (Character; one of X,Y or Z)
        self.direction = direction
        self.locs = array(locs, dtype='float64')
        self.vals = array(vals, dtype='float64')

    def validate(self):
        assert max(abs(self.N)) > 0.
        assert self.direction in ['X', 'Y', 'Z'], 'dir=%r' % self.direction

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a ACCEL card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        cid = integer_or_blank(card, 2, 'cid', 0)
        N = [double_or_blank(card, 3, 'N1', 0.0),
             double_or_blank(card, 4, 'N2', 0.0),
             double_or_blank(card, 5, 'N3', 0.0)]
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
        return ACCEL(sid, N, direction, locs, vals, cid=cid, comment=comment)

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
        try:
            self.cid = model.Coord(self.cid, msg=msg)
            self.cid_ref = self.cid
        except KeyError:
            pass

    def Cid(self):
        if isinstance(self.cid, integer_types):
            return self.cid
        return self.cid_ref.cid

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

    def __init__(self, sid, scale, N, nodes, cid=0, comment=''):
        """
        Creates an ACCEL1 card

        Parameters
        ----------
        sid : int
            load id
        scale : float
            scale factor for load
        N : (3, ) float ndarray
            the acceleration vector in the cid frame
        direction : str
            Component direction of acceleration variation
            {X, Y, Z}
        nodes : List[int]
            the nodes to apply acceleration to
        cid : int; default=0
            the coordinate system for the load
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        #: Load set identification number (Integer>0)
        self.sid = sid

        #: Coordinate system identification number. (Integer>0: Default=0)
        self.cid = cid

        #: Acceleration vector scale factor. (Real)
        self.scale = scale

        #: Components of the acceleration vector measured in coordinate system
        #: CID. (Real; at least one Ni != 0)
        self.N = np.asarray(N)

        #: nodes to apply the acceleration to
        self.nodes = expand_thru_by(nodes)

        assert max(abs(self.N)) > 0.

    def validate(self):
        assert len(self.N) == 3, 'N=%r' % self.N
        assert isinstance(self.cid, integer_types), 'cid=%r' % self.cid
        assert isinstance(self.scale, float), 'cid=%r' % self.scale
        assert isinstance(self.nodes, list), 'nodes=%r' % self.nodes

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a ACCEL1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        cid = integer_or_blank(card, 2, 'cid', 0)
        scale = double(card, 3, 'scale')
        N = [double_or_blank(card, 4, 'N1', 0.0),
             double_or_blank(card, 5, 'N2', 0.0),
             double_or_blank(card, 6, 'N3', 0.0)]

        nodes = fields(integer_or_string, card, 'node', i=9, j=len(card))
        return ACCEL1(sid, scale, N, nodes, cid=cid, comment=comment)

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

    def normalize(self, msg=''):
        """
        adjust the vector to a unit length
        scale up the magnitude of the vector
        """
        assert abs(self.mag) > 0, 'mag=%s\n%s' % (self.mag, self)
        if abs(self.mag) != 0.0:  # enforced displacement
            norm_xyz = norm(self.xyz)
            if norm_xyz == 0.0:
                raise RuntimeError('xyz=%s norm_xyz=%s' % (self.xyz, norm_xyz))
            #mag = self.mag * norm_xyz
            self.mag *= norm_xyz
            try:
                self.xyz = self.xyz / norm_xyz
            except FloatingPointError:
                msgi = 'xyz = %s\n' % self.xyz
                msgi += 'norm_xyz = %s\n' % norm_xyz
                msgi += 'card =\n%s' % str(self)
                msgi += msg
                raise FloatingPointError(msgi)

    def transform_load(self):
        xyz = self.cid_ref.transform_node_to_global(self.xyz)
        if self.mag > 0.:
            return (True, self.node, self.mag * xyz)  # load
        return (False, self.node, xyz)  # enforced displacement

    def get_loads(self):
        return [self]

    def F(self):
        return self.xyz * self.mag

    def get_reduced_loads(self, resolve_load_card=False, filter_zero_scale_factors=False):
        scale_factors = [1.]
        loads = self.F()
        return(scale_factors, loads)

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
            norm_xyz = norm(self.xyz)
            #mag = self.mag*norm_xyz
            self.mag *= norm_xyz
            self.xyz /= norm_xyz

    def transform_load(self):
        #print("self.xyz = ",self.xyz)
        xyz = self.cid_ref.transform_node_to_global(self.xyz)
        if self.mag > 0.:
            #print("mag=%s xyz=%s" % (self.mag, xyz))
            return (True, self.node, self.mag * xyz)  # load
        return (False, self.node, xyz)  # enforced displacement

    def get_loads(self):
        return [self]

    def get_reduced_loads(self, resolve_load_card=False, filter_zero_scale_factors=False):
        scale_factors = [1.]
        loads = {
            self.node: self.M()
        }
        return(scale_factors, loads)

    def M(self):
        return {self.node_id : self.xyz * self.mag}

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class FORCE(Force):
    """
    +-------+-----+------+-------+------+------+------+------+
    |   1   |  2  |  3   |   4   |  5   |  6   |   7  |   8  |
    +=======+=====+======+=======+======+======+======+======+
    | FORCE | SID | NODE |  CID  | MAG  |  FX  |  FY  |  FZ  |
    +-------+-----+------+-------+------+------+------+------+

    +-------+-----+------+-------+------+------+------+------+
    | FORCE |  3  |  1   |       | 100. |  0.  |  0.  |  1.  |
    +-------+-----+------+-------+------+------+------+------+
    """
    type = 'FORCE'

    def __init__(self, sid, node, mag, xyz, cid=0, comment=''):
        """
        Creates a FORCE card

        Parameters
        ----------
        sid : int
            load id
        node : int
            the node to apply the load to
        mag : float
            the load's magnitude
        xyz : (3, ) float ndarray
            the load direction in the cid frame
        cid : int; default=0
            the coordinate system for the load
        comment : str; default=''
            a comment for the card
        """
        Force.__init__(self)
        if comment:
            self.comment = comment
        self.sid = sid
        self.node = node
        self.cid = cid
        self.mag = mag
        self.xyz = np.asarray(xyz, dtype='float64')
        assert self.xyz.size == 3, self.xyz.shape
        assert isinstance(self.cid, int), self.cid

    def validate(self):
        assert isinstance(self.cid, int), self.cid
        assert isinstance(self.mag, float), self.mag
        assert self.xyz.size == 3, self.xyz.shape

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a FORCE card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        node = integer(card, 2, 'node')
        cid = integer_or_blank(card, 3, 'cid', 0)
        mag = double(card, 4, 'mag')
        xyz = array([double_or_blank(card, 5, 'X1', 0.0),
                     double_or_blank(card, 6, 'X2', 0.0),
                     double_or_blank(card, 7, 'X3', 0.0)])
        assert len(card) <= 8, 'len(FORCE card) = %i\ncard=%s' % (len(card), card)
        return FORCE(sid, node, mag, xyz, cid=cid, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a FORCE card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        sid = data[0]
        node = data[1]
        cid = data[2]
        mag = data[3]
        xyz = array(data[4:7])
        return FORCE(sid, node, mag, xyz, cid=cid, comment=comment)

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
        """
        Creates a FORCE1 card

        Parameters
        ----------
        sid : int
            load id
        node : int
            the node to apply the load to
        mag : float
            the load's magnitude
        n1 / n2 : int / int
            defines the load direction
            n = n2 - n1
        comment : str; default=''
            a comment for the card
        """
        Force.__init__(self)
        if comment:
            self.comment = comment
        self.sid = sid
        self.node = node
        self.mag = mag
        self.g1 = g1
        self.g2 = g2

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a FORCE1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        node = integer(card, 2, 'node')
        mag = double(card, 3, 'mag')
        g1 = integer(card, 4, 'g1')
        g2 = integer(card, 5, 'g2')
        assert len(card) == 6, 'len(FORCE1 card) = %i\ncard=%s' % (len(card), card)
        return FORCE1(sid, node, mag, g1, g2, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a FORCE1 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
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

    +--------+-----+---+---+----+----+----+----+
    |   1    |  2  | 3 | 4 |  5 |  6 |  7 |  8 |
    +========+=====+===+===+====+====+====+====+
    | FORCE2 | SID | G | F | G1 | G2 | G3 | G4 |
    +--------+-----+---+---+----+----+----+----+
    """
    type = 'FORCE2'

    def __init__(self, sid, node, mag, g1, g2, g3, g4, comment=''):
        """
        Creates a FORCE2 card

        Parameters
        ----------
        sid : int
            load id
        node : int
            the node to apply the load to
        mag : float
            the load's magnitude
        g1 / g2 / g3 / g4 : int / int / int / int
            defines the load direction
            n = (g2 - g1) x (g4 - g3)
        comment : str; default=''
            a comment for the card
        """
        Force.__init__(self)
        if comment:
            self.comment = comment
        self.sid = sid
        self.node = node
        self.mag = mag
        self.g1 = g1
        self.g2 = g2
        self.g3 = g3
        self.g4 = g4

    def validate(self):
        assert isinstance(self.sid, integer_types), str(self)
        assert self.g1 is not None, self.g1
        assert self.g2 is not None, self.g2
        assert self.g3 is not None, self.g3
        assert self.g4 is not None, self.g4
        assert self.g1 != self.g2, 'g1=%s g2=%s' % (self.g1, self.g2)
        assert self.g3 != self.g4, 'g3=%s g4=%s' % (self.g3, self.g4)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a FORCE2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
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
        """
        Adds a FORCE2 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
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

        xyz1 = self.g1_ref.get_position()
        xyz2 = self.g2_ref.get_position()
        xyz3 = self.g3_ref.get_position()
        xyz4 = self.g4_ref.get_position()
        v21 = xyz2 - xyz1
        v43 = xyz4 - xyz3
        try:
            v21 /= norm(v21)
        except FloatingPointError:
            msg = 'v21=%s norm(v21)=%s\n' % (v21, norm(v21))
            msg += 'g1.get_position()=%s\n' % xyz1
            msg += 'g2.get_position()=%s' % xyz2
            raise FloatingPointError(msg)

        try:
            v43 /= norm(v43)
        except FloatingPointError:
            msg = 'v43=%s norm(v43)=%s\n' % (v43, norm(v43))
            msg += 'g3.get_position()=%s\n' % xyz3
            msg += 'g4.get_position()=%s' % xyz4
            raise FloatingPointError(msg)
        self.xyz = cross(v21, v43)

        msgi = 'xyz1=%s xyz2=%s xyz3=%s xyz4=%s\nv21=%s v43=%s\nxyz=%s' % (
            xyz1, xyz2, xyz3, xyz4, v21, v43, self.xyz)
        #print(msgi)
        self.normalize(msgi)
        #print(self.xyz)

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

        xyz1 = self.g1_ref.get_position()
        xyz2 = self.g2_ref.get_position()
        xyz3 = self.g3_ref.get_position()
        xyz4 = self.g4_ref.get_position()
        v21 = xyz2 - xyz1
        v43 = xyz4 - xyz3
        try:
            v21 /= norm(v21)
        except FloatingPointError:
            msg = 'v21=%s norm(v21)=%s\n' % (v21, norm(v21))
            msg += 'g1.get_position()=%s\n' % xyz1
            msg += 'g2.get_position()=%s' % xyz2
            raise FloatingPointError(msg)

        try:
            v43 /= norm(v43)
        except FloatingPointError:
            msg = 'v43=%s norm(v43)=%s\n' % (v43, norm(v43))
            msg += 'g3.get_position()=%s\n' % xyz3
            msg += 'g4.get_position()=%s' % xyz4
            raise FloatingPointError(msg)
        self.xyz = cross(v21, v43)
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
    """
    Defines a static concentrated moment at a grid point by specifying a
    scale factor and a vector that determines the direction.::

    +--------+-----+---+-----+-----+-----+-----+-----+
    |   1    |  2  | 3 |  4  |  5  |  6  |  7  |  8  |
    +========+=====+===+=====+=====+=====+=====+=====+
    | MOMENT | SID | G | CID |  M  |  N1 |  N2 |  N3 |
    +--------+-----+---+-----+-----+-----+-----+-----+
    | MOMENT |  2  | 5 |  6  | 2.9 | 0.0 | 1.0 | 0.0 |
    +--------+-----+---+-----+-----+-----+-----+-----+
    """
    type = 'MOMENT'

    def __init__(self, sid, node, mag, xyz, cid=0, comment=''):
        """
        Creates a MOMENT card

        Parameters
        ----------
        sid : int
            load id
        node : int
            the node to apply the load to
        mag : float
            the load's magnitude
        xyz : (3, ) float ndarray
            the load direction in the cid frame
        cid : int; default=0
            the coordinate system for the load
        comment : str; default=''
            a comment for the card
        """
        Moment.__init__(self)
        if comment:
            self.comment = comment
        self.sid = sid
        self.node = node
        self.cid = cid
        self.mag = mag
        self.xyz = np.asarray(xyz, dtype='float64')
        assert self.xyz.size == 3, self.xyz.shape

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a MOMENT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        node = integer(card, 2, 'node')
        cid = integer_or_blank(card, 3, 'cid', 0)
        mag = double(card, 4, 'mag')

        xyz = array([double_or_blank(card, 5, 'X1', 0.0),
                     double_or_blank(card, 6, 'X2', 0.0),
                     double_or_blank(card, 7, 'X3', 0.0)])
        assert len(card) <= 8, 'len(MOMENT card) = %i\ncard=%s' % (len(card), card)
        return MOMENT(sid, node, mag, xyz, cid=cid, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a MOMENT card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        sid = data[0]
        node = data[1]
        cid = data[2]
        mag = data[3]
        xyz = data[4:7]
        return MOMENT(sid, node, mag, xyz, cid=cid, comment=comment)

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
    """
    Defines a static concentrated moment at a grid point by specifying a
    magnitude and two grid points that determine the direction.::

    +---------+-----+---+---+----+----+
    |    1    |  2  | 3 | 4 | 5  | 6  |
    +=========+=====+===+===+====+====+
    | MOMENT1 | SID | G | M | G1 | G2 |
    +---------+-----+---+---+----+----+
    """
    type = 'MOMENT1'

    def __init__(self, sid, node, mag, g1, g2, comment=''):
        """
        Creates a MOMENT1 card

        Parameters
        ----------
        sid : int
            load id
        node : int
            the node to apply the load to
        mag : float
            the load's magnitude
        n1 / n2 : int / int
            defines the load direction
            n = n2 - n1
        comment : str; default=''
            a comment for the card
        """
        Moment.__init__(self)
        if comment:
            self.comment = comment
        self.xyz = None
        self.sid = sid
        self.node = node
        self.mag = mag
        self.g1 = g1
        self.g2 = g2

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a MOMENT1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        node = integer(card, 2, 'node')
        mag = double(card, 3, 'mag')
        g1 = integer(card, 4, 'g1')
        g2 = integer(card, 5, 'g2')
        assert len(card) == 6, 'len(MOMENT1 card) = %i\ncard=%s' % (len(card), card)
        return MOMENT1(sid, node, mag, g1, g2, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a MOMENT1 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
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
    """
    Defines a static concentrated moment at a grid point by specification
    of a magnitude and four grid points that determine the direction.::

    +---------+-----+---+---+----+----+----+----+
    |    1    |  2  | 3 | 4 |  5 |  6 |  7 |  8 |
    +=========+=====+===+===+====+====+====+====+
    | MOMENT2 | SID | G | M | G1 | G2 | G3 | G4 |
    +---------+-----+---+---+----+----+----+----+
    """
    type = 'MOMENT2'

    def __init__(self, sid, node, mag, g1, g2, g3, g4, xyz=None, comment=''):
        """
        Creates a MOMENT2 card

        Parameters
        ----------
        sid : int
            load id
        node : int
            the node to apply the load to
        mag : float
            the load's magnitude
        g1 / g2 / g3 / g4 : int / int / int / int
            defines the load direction
            n = (g2 - g1) x (g4 - g3)
        comment : str; default=''
            a comment for the card
        """
        Moment.__init__(self)
        if comment:
            self.comment = comment
        self.sid = sid
        self.node = node
        self.mag = mag
        self.g1 = g1
        self.g2 = g2
        self.g3 = g3
        self.g4 = g4
        self.xyz = xyz

    def validate(self):
        assert isinstance(self.sid, integer_types), str(self)
        assert self.g1 is not None, self.g1
        assert self.g2 is not None, self.g2
        assert self.g3 is not None, self.g3
        assert self.g4 is not None, self.g4
        assert self.g1 != self.g2, 'g1=%s g2=%s' % (self.g1, self.g2)
        assert self.g3 != self.g4, 'g3=%s g4=%s' % (self.g3, self.g4)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a MOMENT2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
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
        """
        Adds a MOMENT2 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
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
            self.comment = comment
        self.sid = sid
        self.cid = cid
        self.normal = normal
        self.entity = entity
        self.entity_id = entity_id
        self.method = method
        self.load_magnitudes = load_magnitudes

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a GMLOAD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
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

    def __init__(self, sid, pressure, nodes, comment=''):
        """
        Creates a PLOAD card, which defines a uniform pressure load on a
        shell/solid face or arbitrarily defined quad/tri face

        Parameters
        ----------
        sid : int
            load id
        pressure : float
            the pressure to apply
        nodes : List[int]
            The nodes that are used to define the normal are defined
            using the same method as the CTRIA3/CQUAD4 normal.
            n = 3 or 4
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        self.sid = sid
        self.pressure = pressure
        self.nodes = nodes
        assert len(self.nodes) in [3, 4], 'nodes=%s' % self.nodes

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PLOAD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        pressure = double(card, 2, 'pressure')
        nodes = [integer(card, 3, 'n1'),
                 integer(card, 4, 'n2'),
                 integer(card, 5, 'n3')]
        n4 = integer_or_blank(card, 6, 'n4', 0)
        if n4:
            nodes.append(n4)
        assert len(card) <= 7, 'len(PLOAD card) = %i\ncard=%s' % (len(card), card)
        return PLOAD(sid, pressure, nodes, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a PLOAD card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        sid = data[0]
        pressure = data[1]
        nodes = data[2:]
        return PLOAD(sid, pressure, nodes, comment=comment)

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

    def get_loads(self):
        return [self]

    def raw_fields(self):
        list_fields = ['PLOAD', self.sid, self.pressure] + self.node_ids
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

    def __init__(self, sid, eid, Type, scale, x1, p1, x2=None, p2=None, comment=''):
        """
        Creates a PLOAD1 card, which may be applied to a CBAR/CBEAM

        Parameters
        ----------
        sid : int
            load id
        eid : int
            element to apply the load to
        Type : str
            type of load that's applied
            valid_types = {FX, FY, FZ, FXE, FYE, FZE,
                           MX, MY, MZ, MXE, MYE, MZE}
        scale : float
            local pressure scaling factor
        x1 / x2 : float / float
            the starting/end position for the load application
            the default for x2 is x1
        p1 / p2 : float / float
            the magnitude of the load at x1 and x2
            the default for p2 is p1
        comment : str; default=''
            a comment for the card

        Point Load       : x1 == x2
        Distributed Load : x1 != x2
        """
        if comment:
            self.comment = comment
        if x2 is None:
            x2 = x1
        if p2 is None:
            p2 = p1
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
        """
        Adds a PLOAD1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
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
        """
        Adds a PLOAD1 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
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

    def transform_load(self):
        p1 = self.eid_ref.ga_ref.get_position()
        p2 = self.eid_ref.gb_ref.get_position()

        g0 = self.eid_ref.g0_vector
        #if not isinstance(g0, ndarray):
            #g0 = g0.get_position()

        x = p2 - p1
        y = p1 - g0
        z = cross(x, y)
        A = [x, y, z]
        #g = self.GravityVector()
        return A
        #(g2, matrix) = self.cid.transformToGlobal(A)
        #return (g2)

    def get_reduced_loads(self, resolve_load_card=False, filter_zero_scale_factors=False):
        """
        Get all load objects in a simplified form, which means all
        scale factors are already applied and only base objects
        (no LOAD cards) will be returned.

        .. todo:: lots more object types to support
        """
        scale_factors = [1.0]
        loads = [self]
        return scale_factors, loads

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
    """
    +--------+-----+------+------+------+------+------+------+------+
    |    1   |   2 |  3   |  4   |   5  |   6  |   7  |   8  |   9  |
    +========+=====+======+======+======+=============+======+======+
    | PLOAD2 | SID |  P   | EID1 | EID2 | EID3 | EID4 | EID5 | EID6 |
    +--------+-----+------+------+------+------+------+------+------+
    | PLOAD2 | 21  | -3.6 |  4   |  16  |  2   |      |      |      |
    +--------+-----+------+------+------+------+------+------+------+
    | PLOAD2 | SID |  P   | EID1 | THRU | EID2 |      |      |      |
    +--------+-----+------+------+------+------+------+------+------+
    """
    type = 'PLOAD2'

    def __init__(self, sid, pressure, eids, comment=''):
        """
        Creates a PLOAD2 card, which defines an applied load normal to the quad/tri face

        Parameters
        ----------
        sid : int
            load id
        pressure : float
            the pressure to apply to the elements
        eids : List[int]
            the elements to apply pressure to
            n < 6 or a continouus monotonic list of elements (e.g., [1, 2, ..., 1000])
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        if isinstance(eids, integer_types):
            self.eids = [eids]
        self.sid = sid
        self.pressure = pressure
        self.eids = eids

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PLOAD2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
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
        """
        Adds a PLOAD2 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
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

    @property
    def element_ids(self):
        if isinstance(self.eids[0], int):
            eids = self.eids
        else:
            eids = [elem.eid for elem in self.eids]
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
    """
    Standard Format
    ===============
    Defines a pressure load on a face of a CHEXA, CPENTA, CTETRA,
    CTRIA3, CTRIA6, CTRIAR, CQUAD4, CQUAD8, or CQUADR element.

    +--------+-----+-----+----+----+------+------+------+-------+
    |   1    |  2  |  3  |  4 |  5 |  6   |   7  |   8  |   9   |
    +========+=====+=====+====+====+======+======+======+=======+
    | PLOAD4 | SID | EID | P1 | P2 |  P3  |  P4  |  G1  | G3/G4 |
    +--------+-----+-----+----+----+------+------+------+-------+
    |        | CID | N1  | N2 | N3 | SORL | LDIR |      |       |
    +--------+-----+-----+----+----+------+------+------+-------+

    Alternate Format
    ================
    Defines a pressure load on a face of a CTRIA3, CTRIA6, CTRIAR,
    CQUAD4, CQUAD8, or CQUADR element.
    +--------+-----+-----+----+----+------+------+------+-------+
    |   1    |  2  |  3  |  4 |  5 |  6   |   7  |   8  |   9   |
    +========+=====+=====+====+====+======+======+======+=======+
    | PLOAD4 | SID | EID | P1 | P2 |  P3  |  P4  | THRU | EID2  |
    +--------+-----+-----+----+----+------+------+------+-------+
    |        | CID | N1  | N2 | N3 | SORL | LDIR |      |       |
    +--------+-----+-----+----+----+------+------+------+-------+

    .. warning:: NX does not support SORL and LDIR, MSC does
    """
    type = 'PLOAD4'

    def __init__(self, sid, eids, pressures,
                 g1=None, g34=None, cid=0, nvector=None, surf_or_line='SURF',
                 line_load_dir='NORM', comment=''):
        """
        Creates a PLOAD4 card

        Parameters
        ----------
        sid : int
            the load id
        eids : List[int, ...]
            shells : the range of element ids; must be sequential
            solids : must be length 1
        pressures : List[float, float, float, float]
            tri : must be length 4 (the last value should be the same as the 0th value)
            quad : must be length 4
        g1 : int/None
            only used for solid elements
        g34 : int / None
            only used for solid elements
        cid : int; default=0
            the coordinate system for ???
        nvector : (3, ) float ndarray
           blank : load acts normal to the face
           the local pressure vector
        surf_or_line : str; default='SURF'
           SURF : surface load
           LINE : line load    (only defined for QUADR, TRIAR)
           not supported
        line_load_dir : str; default='NORM'
           direction of the line load (see surf_or_line); {X, Y, Z, TANG, NORM}
           not supported
        comment : str; default=''
            a comment for the card

        TODO: fix the way "pressures" works
        """
        if nvector is None:
            nvector = np.zeros(3, dtype='float64')
        else:
            nvector = np.asarray(nvector, dtype='float64')

        if comment:
            self.comment = comment
        self.sid = sid

        # these can be greater than 1 if it's a shell (not a solid)
        self.eids = eids
        self.pressures = np.asarray(pressures)

        #: used for solid element only
        self.g1 = g1
        #: g3/g4 - different depending on CHEXA/CPENTA or CTETRA
        self.g34 = g34

        #: Coordinate system identification number. See Remark 2.
        #: (Integer >= 0;Default=0)
        self.cid = cid
        self.nvector = nvector

        # flag with values of SURF/LINE
        self.surf_or_line = surf_or_line

        # Line load direction
        #
        #   1. X, Y, Z : line load in x/y/z in the element coordinate
        #                system
        #   2. TANG    : line load is tangent to the edge pointing
        #                from G1 to G2
        #   3. NORM    : line load is in the mean plane, normal to the
        #                edge and pointing outwards from the element
        #
        #   if cid=N123 = 0: line_load_dir_default=NORM
        self.line_load_dir = line_load_dir

    def validate(self):
        if self.surf_or_line not in ['SURF', 'LINE']:
            raise RuntimeError(self.surf_or_line)
        if self.line_load_dir not in ['LINE', 'X', 'Y', 'Z', 'TANG', 'NORM']:
            raise RuntimeError(self.line_load_dir)
        assert self.g1 != 0, str(self)
        assert self.g34 != 0, str(self)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PLOAD4 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        eid = integer(card, 2, 'eid')
        p1 = double_or_blank(card, 3, 'p1', 0.0)
        pressures = [
            p1,
            double_or_blank(card, 4, 'p2', p1),
            double_or_blank(card, 5, 'p3', p1),
            double_or_blank(card, 6, 'p4', p1)]

        eids = [eid]
        g1_thru = integer_string_or_blank(card, 7, 'g1/THRU')
        if g1_thru == 'THRU' and integer_or_blank(card, 8, 'eid2'):
            # alternate form
            eid2 = integer(card, 8, 'eid2')
            if eid2:
                eids = list(unique(
                    expand_thru([eid, 'THRU', eid2], set_fields=False, sort_fields=False)
                ))
            g1 = None
            g34 = None
        else:
            # standard form
            eids = [eid]
            g1 = integer_or_blank(card, 7, 'g1')
            g34 = integer_or_blank(card, 8, 'g34')

        cid = integer_or_blank(card, 9, 'cid', 0)
        nvector = array([double_or_blank(card, 10, 'N1', 0.0),
                         double_or_blank(card, 11, 'N2', 0.0),
                         double_or_blank(card, 12, 'N3', 0.0)])
        surf_or_line = string_or_blank(card, 13, 'sorl', 'SURF')
        line_load_dir = string_or_blank(card, 14, 'ldir', 'NORM')
        assert len(card) <= 15, 'len(PLOAD4 card) = %i\ncard=%s' % (len(card), card)
        return PLOAD4(sid, eids, pressures, g1, g34, cid, nvector,
                      surf_or_line, line_load_dir, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a PLOAD4 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
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
        nvector = data[6]

        surf_or_line = data[7]
        #self.line_load_dir = data[8]
        #assert len(data)==8

        eids = [eid]
        surf_or_line = 'SURF'
        line_load_dir = 'NORM'
        return PLOAD4(sid, eids, pressures, g1, g34, cid, nvector,
                      surf_or_line, line_load_dir, comment=comment)

    def get_loads(self):
        return [self]

    def transform_load(self):
        """
        Considers single elememnts

        .. warning:: surf_or_line=SURF is supported (not LINE)
        .. warning:: line_load_dir=NORM is supported (not X,Y,Z)
        """
        if self.surf_or_line != 'SURF':
            msg = 'Only surface loads are supported.  required_surf_or_line=SURF.  actual=%r' % self.surf_or_line
            raise RuntimeError(msg)
        if self.line_load_dir != 'NORM':
            msg = 'Only normal loads are supported.   required_line_load_dir=NORM.  actual=%r' % self.line_load_dir
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

        if  self.surf_or_line != 'SURF':
            if norm(self.nvector) != 0.0 or self.cid != 0:
                vector = self.nvector / np.linalg.norm(self.nvector)
                assert self.Cid() == 0, 'cid=%r on a PLOAD4 is not supported\n%s' % (self.Cid(), str(self))
            else:
                # normal pressure
                assert len(self.eids_ref) == 1, 'only 1 element is supported by transform_load on PLOAD4\n%s' % (str(self))
                elem = self.eids_ref[0]
                vector = array(elem.Normal())
        else:
            raise NotImplementedError('surf_or_line=%r on PLOAD4 is not supported\n%s' % (self.surf_or_line, str(self)))

        vectors = []
        for (nid, p) in zip(face_node_ids, self.pressures):
            vectors.append(vector * p * area / n)  # Force_i
        is_load = None
        return (is_load, face_node_ids, vectors)

    def Cid(self):
        """gets the coordinate system object"""
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
        try:
            self.cid = model.Coord(self.cid, msg=msg)
            self.cid_ref = self.cid
        except KeyError:
            model.log.warning('Could not find cid=%s%s' % (self.cid, msg))

        #self.eid_ref = self.eid
        if self.g1 is not None:
            try:
                self.g1 = model.Node(self.g1, msg=msg)
                self.g1_ref = self.g1
            except KeyError:
                model.log.warning('Could not find g1=%s%s' % (self.g1, msg))

        if self.g34 is not None:
            try:
                self.g34 = model.Node(self.g34, msg=msg)
                self.g34_ref = self.g34
            except KeyError:
                model.log.warning('Could not find g34=%s%s' % (self.g34, msg))

        #if self.eids:
        msgia = 'Could not find element=%%s, %s\n' % msg
        self.eids, msgi = model.safe_get_elements(self.eids, msg=msgia)
        self.eids_ref = self.eids
        if msgi:
            model.log.warning(msgi.rstrip())

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
        if isinstance(self.g1, integer_types):
            return self.g1
        return self.g1_ref.nid

    def G34(self):
        if isinstance(self.g34, integer_types):
            return self.g34
        return self.g34_ref.nid

    def Eid(self, eid):
        if isinstance(eid, integer_types):
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
            return eidi
        eids = []
        for element in self.eids:
            eidi = self.Eid(element)
            eids.append(eidi)
        return eids

    @property
    def element_ids(self):
        return self.get_element_ids()

    def repr_fields(self):
        eids = self.element_ids
        eid = eids[0]
        p1 = self.pressures[0]
        p2 = set_blank_if_default(self.pressures[1], p1)
        p3 = set_blank_if_default(self.pressures[2], p1)
        p4 = set_blank_if_default(self.pressures[3], p1)
        list_fields = ['PLOAD4', self.sid, eid, self.pressures[0], p2, p3, p4]

        if self.g1 is not None:
            # is it a SOLID element
            node_ids = self.node_ids
            #node_ids = self.node_ids([self.g1, self.g34])
            list_fields += node_ids
        else:
            if len(eids) > 1:
                try:
                    list_fields.append('THRU')
                    eidi = eids[-1]
                except:
                    print("g1  = %s" % self.g1)
                    print("g34 = %s" % self.g34)
                    print("self.eids = %s" % self.eids)
                    raise
                list_fields.append(eidi)
            else:
                list_fields += [None, None]

        cid = self.Cid()
        if cid or norm(self.nvector) > 0.0:
            n1 = self.nvector[0]
            n2 = self.nvector[1]
            n3 = self.nvector[2]
            list_fields.append(cid)
            list_fields += [n1, n2, n3]
            surf_or_line = self.surf_or_line
            line_load_dir = self.line_load_dir
        else:
            list_fields += [None, None, None, None]
            surf_or_line = set_blank_if_default(self.surf_or_line, 'SURF')
            line_load_dir = set_blank_if_default(self.line_load_dir, 'NORM')
        list_fields.append(surf_or_line)
        if surf_or_line == 'LINE':
            list_fields.append(line_load_dir)
        return list_fields

    def raw_fields(self):
        eids = self.element_ids
        eid = eids[0]
        p1 = self.pressures[0]
        p2 = self.pressures[1]
        p3 = self.pressures[2]
        p4 = self.pressures[3]
        list_fields = ['PLOAD4', self.sid, eid, self.pressures[0], p2, p3, p4]

        if self.g1 is not None:
            # is it a SOLID element
            node_ids = self.node_ids
            #node_ids = self.node_ids([self.g1, self.g34])
            list_fields += node_ids
        else:
            if len(eids) > 1:
                try:
                    list_fields.append('THRU')
                    eidi = eids[-1]
                except:
                    print("g1  = %s" % self.g1)
                    print("g34 = %s" % self.g34)
                    print("self.eids = %s" % self.eids)
                    raise
                list_fields.append(eidi)
            else:
                list_fields += [None, None]

        cid = self.Cid()
        if cid or norm(self.nvector) > 0.0:
            n1 = self.nvector[0]
            n2 = self.nvector[1]
            n3 = self.nvector[2]
            list_fields.append(cid)
            list_fields += [n1, n2, n3]
        else:
            list_fields += [None, None, None, None]

        surf_or_line = self.surf_or_line
        line_load_dir = self.line_load_dir
        list_fields.append(surf_or_line)
        if surf_or_line == 'LINE':
            list_fields.append(line_load_dir)
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PLOADX1(Load):
    """
    Pressure Load on Axisymmetric Element

    Defines surface traction to be used with the CQUADX, CTRIAX, and CTRIAX6
    axisymmetric element.
    +---------+-----+-----+----+----+----+----+-------+
    |    1    |  2  |  3  |  4 |  5 |  6 |  7 |   8   |
    +=========+=====+=====+====+====+====+====+=======+
    | PLOADX1 | SID | EID | PA | PB | GA | GB | THETA |
    +---------+-----+-----+----+----+----+----+-------+
    """
    type = 'PLOADX1'

    def __init__(self, sid, eid, pa, nids, pb=None, theta=0., comment=''):
        """
        Creates a PLOADX1 card, which defines surface traction for
        axisymmetric elements.

        Parameters
        ----------
        sid : int
            load id
        eid : int
            element id (CQUADX, CTRIAX, or CTRIAX6)
        nids : List[int, int]
            Corner grid points.
            GA and GB are any two adjacent corner grid points of the element
        pa / pb : float / None
            Surface traction at grid point GA or GB
            pb : default is None -> pa
        theta : float; default=0.0
            Angle between surface traction and inward normal to the line
            segment.
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        if pb is None:
            pb = pa
        self.sid = sid
        self.eid = eid
        self.pa = pa
        self.pb = pb
        self.ga = nids[0]
        self.gb = nids[1]
        self.theta = theta

    def validate(self):
        assert isinstance(self.ga, integer_types), 'ga=%r' % self.ga
        assert isinstance(self.gb, integer_types), 'gb=%r' % self.gb
        assert isinstance(self.pa, float), 'pa=%r' % self.pa
        assert isinstance(self.pb, float), 'pb=%r' % self.pb

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PLOADX1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        eid = integer(card, 2, 'eid')
        pa = double(card, 3, 'pa')
        pb = double_or_blank(card, 4, 'pb', pa)
        ga = integer(card, 5, 'ga')
        gb = integer(card, 6, 'gb')
        theta = double_or_blank(card, 7, 'theta', 0.)
        assert len(card) <= 8, 'len(PLOADX1 card) = %i\ncard=%s' % (len(card), card)
        return PLOADX1(sid, eid, pa, ga, gb, pb, theta, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a PLOADX1 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        sid, eid, pa, pb, ga, gb, theta = data
        return PLOADX1(sid, eid, pa, ga, gb, pb, theta, comment=comment)

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

    @property
    def nodes(self):
        return [self.ga, self.gb]

    def safe_cross_reference(self, model):
        return self.cross_reference(model)

    def uncross_reference(self):
        self.eid = self.Eid()
        self.ga = self.Ga()
        self.gb = self.Gb()
        del self.eid_ref, self.ga_ref, self.gb_ref

    def Eid(self):
        if isinstance(self.eid, integer_types):
            return self.eid
        return self.eid_ref.eid

    def Ga(self):
        if isinstance(self.ga, integer_types):
            return self.ga
        return self.ga_ref.nid

    def Gb(self):
        if isinstance(self.gb, integer_types):
            return self.gb
        return self.gb_ref.nid

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
