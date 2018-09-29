# pylint: disable=R0902,R0904,R0914
"""
All static loads are defined in this file.  This includes:

 * LOAD
 * GRAV
 * ACCEL
 * ACCEL1
 * FORCE / MOMENT
 * FORCE1 / MOMENT1
 * FORCE2 / MOMENT2
 * MOMENT
 * PLOAD
 * PLOAD2
 * PLOAD4
 * PLOADX1

"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

import numpy as np
from numpy import array, cross, allclose, unique
from numpy.linalg import norm  # type: ignore

#from pyNastran.bdf.errors import CrossReferenceError
from pyNastran.utils.numpy_utils import integer_types, float_types
from pyNastran.bdf.cards.loads.loads import Load, LoadCombination
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import BaseCard, expand_thru, expand_thru_by #  _node_ids,
from pyNastran.bdf.cards.collpase_card import collapse_thru_by

from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, string, string_or_blank,
    integer_or_string, fields, integer_string_or_blank, integer_or_double)
from pyNastran.bdf.field_writer_8 import print_card_8, print_float_8, set_string8_blank_if_default
from pyNastran.bdf.field_writer_16 import (
    print_card_16, print_float_16, set_string16_blank_if_default)
from pyNastran.bdf.field_writer_double import print_card_double, print_scientific_double


class LOAD(LoadCombination):
    """
    +------+-----+------+------+----+-----+----+----+----+
    |   1  |  2  |  3   |  4   | 5  |  6  | 7  | 8  | 9  |
    +======+=====+======+======+====+=====+====+====+====+
    | LOAD | SID |  S   |  S1  | L1 | S2  | L2 | S3 | L3 |
    +------+-----+------+------+----+-----+----+----+----+
    |      | S4  |  L4  | etc. |    |     |    |    |    |
    +------+-----+------+------+----+-----+----+----+----+
    | LOAD | 101 | -0.5 | 1.0  | 3  | 6.2 | 4  |    |    |
    +------+-----+------+------+----+-----+----+----+----+
    """
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

        .. note::  MSC can handle self-referencing loads, NX cannot
        """
        LoadCombination.__init__(self, sid, scale, scale_factors, load_ids,
                                 comment=comment)

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
                        msg = (
                            'A LOAD card cannot reference another LOAD card\n'
                            'current:\n%s\n'
                            'new:\n%s' % (str(self), str(load))
                        )
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

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        load_ids2 = []
        msg = ', which is required by LOAD=%s' % (self.sid)
        for load_id in self.load_ids:
            if load_id == self.sid:
                msg = 'Type=%s sid=%s load_id=%s creates a recursion error' % (
                    self.type, self.sid, load_id)
                raise RuntimeError(msg)
            load_id2 = model.Load(load_id, consider_load_combinations=True, msg=msg)
            assert isinstance(load_id2, list), load_id2
            load_ids2.append(load_id2)
        self.load_ids_ref = load_ids2

    def safe_cross_reference(self, model, xref_errors, debug=True):
        load_ids2 = []
        msg = ', which is required by LOAD=%s' % (self.sid)
        for load_id in self.load_ids:
            try:
                load_id2 = model.Load(load_id, consider_load_combinations=True, msg=msg)
            except KeyError:
                if debug:
                    msg = 'Couldnt find load_id=%i, which is required by %s=%s' % (
                        load_id, self.type, self.sid)
                    print(msg)
                continue
            load_ids2.append(load_id2)
        self.load_ids_ref = load_ids2

    def raw_fields(self):
        list_fields = ['LOAD', self.sid, self.scale]
        load_ids = self.get_load_ids()
        for (scale_factor, load_id) in zip(self.scale_factors, load_ids):
            list_fields += [scale_factor, self.LoadID(load_id)]
        if len(load_ids) != len(self.scale_factors):
            msg = 'nload_ids=%s nscale_factors=%s and arent the same\n' % (
                len(load_ids), len(self.scale_factors))
            msg = 'load_ids=%s\n' % (load_ids)
            msg += 'scale_factors=%s\n' % (self.scale_factors)
            msg += print_card_8(list_fields)
            raise IndexError(msg)
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
        self.load_ids = self.get_load_ids()
        self.load_ids_ref = None


class GRAV(BaseCard):
    """
    Defines acceleration vectors for gravity or other acceleration loading.

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
        self.cid_ref = None

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
        unused_a = data[2]
        N = array(data[3:6])
        mb = data[6]
        scale = 1.
        assert len(data) == 7
        return GRAV(sid, scale, N, cid=cid, mb=mb, comment=comment)

    def get_loads(self):
        return [self]

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by GRAV sid=%s' % self.sid
        self.cid_ref = model.Coord(self.cid, msg=msg)

    def safe_cross_reference(self, model, xref_errors, debug=True):
        # msg = "Couldn't find CORDx=%s which is required by GRAV sid=%s" % (self.cid, self.sid)
        msg = ', which is required by GRAV sid=%s' % self.sid
        self.cid_ref = model.safe_coord(self.cid, self.sid, xref_errors, msg=msg)

    def uncross_reference(self):
        self.cid = self.Cid()
        self.cid_ref = None

    def Cid(self):
        if self.cid_ref is not None:
            return self.cid_ref.cid
        return self.cid

    def GravityVector(self):
        """returns the gravity vector in absolute coordinates"""
        if self.Cid() == 0:
            return self.N
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

    +-------+------+------+--------+------+-----+-----+--------+-----+
    |   1   |   2  |   3  |    4   |   5  |  6  |  7  |   8    |  9  |
    +=======+======+======+========+======+=====+=====+========+=====+
    | ACCEL | SID  | CID  |   N1   |  N2  | N3  | DIR |        |     |
    +-------+------+------+--------+------+-----+-----+--------+-----+
    |       | LOC1 | VAL1 |  LOC2  | VAL2 | Continues in Groups of 2 |
    +-------+------+------+--------+------+--------------------------+
    | ACCEL |  100 |   2  |   0.0  |  1.0 | 2.0 |  X  |        |     |
    +-------+------+------+--------+------+-----+-----+--------+-----+
    |       |  1.0 |  1.1 |   2.0  |  2.1 | 3.0 | 3.1 |  4.0   | 4.1 |
    +-------+------+------+--------+------+-----+-----+--------+-----+
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
        locs : List[float]
            Location along direction DIR in coordinate system CID for
            specification of a load scale factor.
        vals : List[float]
            The load scale factor associated with location LOCi
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
        self.cid_ref = None

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
        msg = ', which is required by ACCEL sid=%s' % self.sid
        self.cid_ref = model.Coord(self.cid, msg=msg)

    def uncross_reference(self):
        self.cid = self.Cid()
        self.cid_ref = None

    def safe_cross_reference(self, model, xref_errors, debug=True):
        msg = ', which is required by ACCEL sid=%s' % self.sid
        self.cid_ref = model.safe_coord(self.cid, self.sid, xref_errors, msg=msg)

    def Cid(self):
        if self.cid_ref is not None:
            return self.cid_ref.cid
        return self.cid

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

    +--------+---------+---------+-----+----+----+----+
    |    1   |    2    |    3    |  4  |  5 |  6 |  7 |
    +========+=========+=========+=====+====+====+====+
    | ACCEL1 |   SID   |   CID   |  A  | N1 | N2 | N3 |
    +--------+---------+---------+-----+----+----+----+
    |        | GRIDID1 | GRIDID2 | etc |    |    |    |
    +--------+---------+---------+-----+----+----+----+
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
        self.nodes_ref = None
        self.cid_ref = None

    def validate(self):
        assert len(self.N) == 3, 'N=%r' % self.N
        assert isinstance(self.cid, integer_types), 'cid=%r' % self.cid
        assert isinstance(self.scale, float_types), 'scale=%r' % self.scale
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
        msg = ', which is required by ACCEL1 sid=%s' % self.sid
        self.cid_ref = model.Coord(self.cid, msg=msg)
        self.nodes_ref = model.EmptyNodes(self.node_ids, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by ACCEL1 sid=%s' % self.sid
        self.cid_ref = model.safe_coord(self.cid, self.sid, xref_errors, msg=msg)
        self.nodes_ref = model.EmptyNodes(self.node_ids, msg=msg)

    def uncross_reference(self):
        self.cid = self.Cid()
        self.nodes = self.node_ids
        self.nodes_ref = None
        self.cid_ref = None

    def Cid(self):
        if self.cid_ref is not None:
            return self.cid_ref.cid
        return self.cid

    @property
    def node_ids(self):
        #msg = ', which is required by ACCEL1 sid=%s' % self.sid
        #_node_ids(self.nodes, allow_empty_nodes=True, msg=msg)
        return self._node_ids(nodes=self.nodes_ref)

    def _node_ids(self, nodes=None):  # this function comes from BaseCard.py
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
            ] + collapse_thru_by(self.node_ids)
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


#class Force(Load):
    #"""Generic class for all Forces"""
    #type = 'Force'

    #def __init__(self):
        #Load.__init__(self)

    #def get_loads(self):
        #return [self]

    #def F(self):
        #return self.xyz * self.mag

    #def get_reduced_loads(self, resolve_load_card=False, filter_zero_scale_factors=False):
        #scale_factors = [1.]
        #loads = self.F()
        #return(scale_factors, loads)

    #def write_card(self, size=8, is_double=False):
        #card = self.raw_fields()
        #if size == 8:
            #return self.comment + print_card_8(card)
        #if is_double:
            #return self.comment + print_card_double(card)
        #return self.comment + print_card_16(card)


#class Moment(Load):
    #"""Generic class for all Moments"""
    #type = 'Moment'

    #def __init__(self):
        #Load.__init__(self)

    #def get_loads(self):
        #return [self]

    #def get_reduced_loads(self, resolve_load_card=False, filter_zero_scale_factors=False):
        #scale_factors = [1.]
        #loads = {
            #self.node: self.M()
        #}
        #return(scale_factors, loads)

    #def write_card(self, size=8, is_double=False):
        #card = self.raw_fields()
        #if size == 8:
            #return self.comment + print_card_8(card)
        #if is_double:
            #return self.comment + print_card_double(card)
        #return self.comment + print_card_16(card)


class Load0(BaseCard):
    """common class for FORCE, MOMENT"""
    def __init__(self, sid, node, mag, xyz, cid=0, comment=''):
        """
        Creates a FORCE/MOMENT card

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
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.sid = sid
        self.node = node
        self.cid = cid
        self.mag = mag
        self.xyz = np.asarray(xyz, dtype='float64')
        assert self.xyz.size == 3, self.xyz.shape
        assert isinstance(self.cid, int), self.cid
        self.node_ref = None
        self.cid_ref = None

    def validate(self):
        assert isinstance(self.cid, int), self.cid
        assert isinstance(self.mag, float), self.mag
        assert self.xyz.size == 3, self.xyz.shape

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a FORCE/MOMENT card from ``BDF.add_card(...)``

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
        assert len(card) <= 8, 'len(%s card) = %i\ncard=%s' % (cls.type, len(card), card)
        return cls(sid, node, mag, xyz, cid=cid, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a FORCE/MOMENT card from the OP2

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
        return cls(sid, node, mag, xyz, cid=cid, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by %s sid=%s' % (self.type, self.sid)
        self.node_ref = model.Node(self.node, msg=msg)
        self.cid_ref = model.Coord(self.cid, msg=msg)

    def safe_cross_reference(self, model, xref_errors, debug=True):
        msg = ', which is required by %s sid=%s' % (self.type, self.sid)
        # try:
        self.node_ref = model.Node(self.node, msg=msg)
        self.cid_ref = model.safe_coord(self.cid, self.sid, xref_errors, msg=msg)

    def uncross_reference(self):
        self.cid = self.Cid()
        self.cid_ref = None

    def get_loads(self):
        return [self]

    @property
    def node_id(self):
        if self.node_ref is not None:
            return self.node_ref.nid
        return self.node

    def Cid(self):
        if self.cid_ref is not None:
            return self.cid_ref.cid
        return self.cid

    @property
    def scaled_vector(self):
        return self.xyz * self.mag

    def raw_fields(self):
        list_fields = [self.type, self.sid, self.node_id, self.Cid(),
                       self.mag] + list(self.xyz)
        return list_fields

    def repr_fields(self):
        cid = set_blank_if_default(self.Cid(), 0)
        list_fields = [self.type, self.sid, self.node_id, cid,
                       self.mag] + list(self.xyz)
        return list_fields


class FORCE(Load0):
    """
    Defines a static concentrated force at a grid point by specifying a
    scale factor and a vector that determines the direction.

    +-------+-----+------+-------+------+------+------+------+
    |   1   |  2  |  3   |   4   |  5   |  6   |   7  |   8  |
    +=======+=====+======+=======+======+======+======+======+
    | FORCE | SID | NODE |  CID  | MAG  |  FX  |  FY  |  FZ  |
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
        Load0.__init__(self, sid, node, mag, xyz, cid=cid, comment=comment)

    def write_card(self, size=8, is_double=False):
        if size == 8:
            cids = set_string8_blank_if_default(self.Cid(), 0)
            msg = 'FORCE   %8i%8i%8s%8s%8s%8s%8s\n' % (
                self.sid, self.node_id,
                cids, print_float_8(self.mag), print_float_8(self.xyz[0]),
                print_float_8(self.xyz[1]), print_float_8(self.xyz[2]))
        else:
            cids = set_string16_blank_if_default(self.Cid(), 0)
            if is_double:
                msg = ('FORCE*  %16i%16i%16s%s\n'
                       '*       %16s%16s%16s\n') % (
                           self.sid, self.node_id,
                           cids, print_scientific_double(self.mag),
                           print_scientific_double(self.xyz[0]),
                           print_scientific_double(self.xyz[1]),
                           print_scientific_double(self.xyz[2]))
            else:
                msg = ('FORCE*  %16i%16i%16s%s\n'
                       '*       %16s%16s%16s\n') % (
                           self.sid, self.node_id,
                           cids, print_float_16(self.mag), print_float_16(self.xyz[0]),
                           print_float_16(self.xyz[1]), print_float_16(self.xyz[2]))
        return self.comment + msg


class Load1(BaseCard):
    def __init__(self, sid, node, mag, g1, g2, comment=''):
        """
        Creates a FORCE1/MOMENT1 card

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
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.sid = sid
        self.node = node
        self.mag = mag
        self.g1 = g1
        self.g2 = g2
        self.node_ref = None
        self.g1_ref = None
        self.g2_ref = None
        self.xyz = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a FORCE1/MOMENT1 card from ``BDF.add_card(...)``

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
        assert len(card) == 6, 'len(%s card) = %i\ncard=%s' % (cls.type, len(card), card)
        return cls(sid, node, mag, g1, g2, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a FORCE1/MOMENT1 card from the OP2

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
        return cls(sid, node, mag, g1, g2, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by %s sid=%s' % (self.type, self.sid)
        self.node_ref = model.Node(self.node, msg=msg)
        self.g1_ref = model.Node(self.g1, msg=msg)
        self.g2_ref = model.Node(self.g2, msg=msg)

        self.xyz = self.g2_ref.get_position() - self.g1_ref.get_position()
        normalize(self)

    def uncross_reference(self):
        self.node = self.node_id
        self.g1 = self.G1()
        self.g2 = self.G2()
        self.node_ref = None
        self.g1_ref = None
        self.g2_ref = None

    def safe_cross_reference(self, model, safe_coord, debug=True):
        """
        .. todo:: cross reference and fix repr function
        """
        return self.cross_reference(model)
        #msg = ', which is required by FORCE1 sid=%s' % self.sid
        #self.node_ref = model.Node(self.node, msg=msg)
        #self.g1_ref = model.Node(self.g1, msg=msg)
        #self.g2_ref = model.Node(self.g2, msg=msg)
        #self.xyz = self.g2.get_position() - self.g1.get_position()
        #normalize(self)

    def get_loads(self):
        return [self]

    @property
    def scaled_vector(self):
        return self.xyz * self.mag

    @property
    def node_ids(self):
        return [self.node_id, self.G1(), self.G2()]

    def G1(self):
        if self.g1_ref is not None:
            return self.g1_ref.nid
        return self.g1

    def G2(self):
        if self.g2_ref is not None:
            return self.g2_ref.nid
        return self.g2

    @property
    def node_id(self):
        if self.node_ref is not None:
            return self.node_ref.nid
        return self.node

    def raw_fields(self):
        list_fields = [self.type, self.sid, self.node_id, self.mag, self.G1(), self.G2()]
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


class FORCE1(Load1):
    """
    Defines a static concentrated force at a grid point by specification of a
    magnitude and two grid points that determine the direction.

    +--------+-----+----+-------+----+----+
    |   1    |  2  | 3  |   4   | 5  | 6  |
    +========+=====+====+=======+====+====+
    | FORCE1 | SID | G  |   F   | G1 | G2 |
    +--------+-----+----+-------+----+----+
    | FORCE1 |  6  | 13 | -2.93 | 16 | 13 |
    +--------+-----+----+-------+----+----+
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
        Load1.__init__(self, sid, node, mag, g1, g2, comment)
        #Force.__init__(self)


class Load2(BaseCard):
    """common class for FORCE2, MOMENT2"""
    def __init__(self, sid, node, mag, g1, g2, g3, g4, comment=''):
        """
        Creates a FORCE2/MOMENT2 card

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
        if comment:
            self.comment = comment
        self.sid = sid
        self.node = node
        self.mag = mag
        self.g1 = g1
        self.g2 = g2
        self.g3 = g3
        self.g4 = g4
        self.node_ref = None
        self.g1_ref = None
        self.g2_ref = None
        self.g3_ref = None
        self.g4_ref = None
        self.xyz = None

    def validate(self):
        assert isinstance(self.sid, integer_types), str(self)
        assert self.g1 is not None, self.g1
        assert self.g2 is not None, self.g2
        assert self.g3 is not None, self.g3
        assert self.g1 != self.g2, 'g1=%s g2=%s' % (self.g1, self.g2)
        assert self.g3 != self.g4, 'g3=%s g4=%s' % (self.g3, self.g4)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a FORCE2/MOMENT2 card from ``BDF.add_card(...)``

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
        g4 = integer_or_blank(card, 7, 'g4')
        assert len(card) in [7, 8], 'len(%s card) = %i\ncard=%s' % (cls.type, len(card), card)
        return cls(sid, node, mag, g1, g2, g3, g4, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a FORCE2/MOMENT2 card from the OP2

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
        return cls(sid, node, mag, g1, g2, g3, g4, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by %s sid=%s' % (self.type, self.sid)
        self.node_ref = model.Node(self.node, msg=msg)
        self.g1_ref = model.Node(self.g1, msg=msg)
        self.g2_ref = model.Node(self.g2, msg=msg)
        self.g3_ref = model.Node(self.g3, msg=msg)

        xyz1 = self.g1_ref.get_position()
        xyz2 = self.g2_ref.get_position()
        xyz3 = self.g3_ref.get_position()
        v21 = xyz2 - xyz1

        try:
            v21 /= norm(v21)
        except FloatingPointError:
            msg = 'v1=v21=%s norm(v21)=%s\n' % (v21, norm(v21))
            msg += 'g1.get_position()=%s\n' % xyz1
            msg += 'g2.get_position()=%s' % xyz2
            raise FloatingPointError(msg)

        if self.g4 is None:
            xyz4 = None
            v2 = xyz3 - xyz1
            try:
                v2 /= norm(v2)
            except FloatingPointError:
                msg = 'v2=v31=%s norm(v31)=%s\n' % (v2, norm(v2))
                msg += 'g3.get_position()=%s\n' % xyz3
                msg += 'g1.get_position()=%s' % xyz1
                raise FloatingPointError(msg)
            xyz = cross(v21, v2)
        else:
            self.g4_ref = model.Node(self.g4, msg=msg)
            xyz4 = self.g4_ref.get_position()
            v2 = xyz4 - xyz3

            try:
                v2 /= norm(v2)
            except FloatingPointError:
                msg = 'v2=v43=%s norm(v43)=%s\n' % (v2, norm(v2))
                msg += 'g3.get_position()=%s\n' % xyz3
                msg += 'g4.get_position()=%s' % xyz4
                raise FloatingPointError(msg)
            xyz = cross(v21, v2)

        self.xyz = xyz

        msgi = 'xyz1=%s xyz2=%s xyz3=%s xyz4=%s\nv21=%s v43 (or v31)=%s\nxyz=%s' % (
            xyz1, xyz2, xyz3, xyz4, v21, v2, self.xyz)
        normalize(self, msgi)

    def safe_cross_reference(self, model, safe_coord, debug=True):
        """
        .. todo:: cross reference and fix repr function
        """
        msg = ', which is required by %s sid=%s' % (self.type, self.sid)
        is_failed = False
        try:
            self.node_ref = model.Node(self.node, msg=msg)
        except KeyError:
            is_failed = True
            model.log.warning('failed to cross-reference NODE=%i%s' % (self.node, msg))

        try:
            self.g1_ref = model.Node(self.g1, msg=msg)
            xyz1 = self.g1_ref.get_position()
        except KeyError:
            is_failed = True
            model.log.warning('failed to cross-reference G1=%i%s' % (self.g1, msg))

        try:
            self.g2_ref = model.Node(self.g2, msg=msg)
            xyz2 = self.g2_ref.get_position()
        except KeyError:
            is_failed = True
            model.log.warning('failed to cross-reference G2=%i%s' % (self.g2, msg))

        try:
            self.g3_ref = model.Node(self.g3, msg=msg)
            xyz3 = self.g3_ref.get_position()
        except KeyError:
            is_failed = True
            model.log.warning('failed to cross-reference G3=%i%s' % (self.g3, msg))

        if not is_failed:
            v21 = xyz2 - xyz1

        if self.g4 is not None:
            try:
                self.g4_ref = model.Node(self.g4, msg=msg)
            except KeyError:
                is_failed = True
            if not is_failed:
                xyz4 = self.g4_ref.get_position()
                model.log.warning('failed to cross-reference G4=%i%s' % (self.g4, msg))
        else:
            xyz3, xyz4 = xyz1, xyz3

        if not is_failed:
            v43 = xyz4 - xyz3
            v2 = v43
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

            #msgi = 'xyz1=%s xyz2=%s xyz3=%s xyz4=%s\nv21=%s v43 (or v31)=%s\nxyz=%s' % (
                #xyz1, xyz2, xyz3, xyz4, v21, v2, self.xyz)
            normalize(self, msg)

    @property
    def scaled_vector(self):
        return self.xyz * self.mag

    def uncross_reference(self):
        self.node = self.node_id
        self.g1 = self.G1()
        self.g2 = self.G2()
        self.g3 = self.G3()
        self.g4 = self.G4()
        self.node_ref = None
        self.g1_ref = None
        self.g2_ref = None
        self.g3_ref = None
        self.g4_ref = None
        self.xyz = None

    def get_loads(self):
        return [self]

    @property
    def node_id(self):
        if self.node_ref is not None:
            return self.node_ref.nid
        return self.node

    def G1(self):
        if self.g1_ref is not None:
            return self.g1_ref.nid
        return self.g1

    def G2(self):
        if self.g2_ref is not None:
            return self.g2_ref.nid
        return self.g2

    def G3(self):
        if self.g3_ref is not None:
            return self.g3_ref.nid
        return self.g3

    def G4(self):
        if self.g4_ref is not None:
            return self.g4_ref.nid
        return self.g4

    @property
    def node_ids(self):
        return [self.node_id, self.G1(), self.G2(), self.G3(), self.G4()]

    def _node_ids(self, nodes=None):
        """returns nodeIDs for repr functions"""
        if not nodes:
            nodes = self.nodes
        if isinstance(nodes[0], integer_types):
            return [node for node in nodes]
        else:
            return [node.nid for node in nodes]

    def raw_fields(self):
        (node, g1, g2, g3, g4) = self._node_ids([self.node, self.g1, self.g2, self.g3, self.g4])
        list_fields = [self.type, self.sid, node, self.mag, g1, g2, g3, g4]
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

    #def get_reduced_loads(self, resolve_load_card=False, filter_zero_scale_factors=False):
        #scale_factors = [1.]
        #loads = self.F()
        #return(scale_factors, loads)


class FORCE2(Load2):
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
        Load2.__init__(self, sid, node, mag, g1, g2, g3, g4, comment)



class MOMENT(Load0):
    """
    Defines a static concentrated moment at a grid point by specifying a
    scale factor and a vector that determines the direction.

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
        Load0.__init__(self, sid, node, mag, xyz, cid=cid, comment=comment)

    def uncross_reference(self):
        self.node = self.node_id
        self.cid = self.Cid()
        self.node_ref = None
        self.cid_ref = None

    @property
    def node_ids(self):
        """all the nodes referenced by the load"""
        return [self.node_id]

    @property
    def node_id(self):
        if self.node_ref is None:
            return self.node
        return self.node_ref.nid

    def raw_fields(self):
        list_fields = ['MOMENT', self.sid, self.node_id, self.Cid(),
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


class MOMENT1(Load1):
    """
    Defines a static concentrated moment at a grid point by specifying a
    magnitude and two grid points that determine the direction.

    +---------+-----+----+-------+----+----+
    |    1    |  2  | 3  |   4   | 5  | 6  |
    +=========+=====+====+=======+====+====+
    | MOMENT1 | SID | G  |   M   | G1 | G2 |
    +---------+-----+----+-------+----+----+
    | MOMENT1 |  6  | 13 | -2.93 | 16 | 13 |
    +---------+-----+----+-------+----+----+
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
        Load1.__init__(self, sid, node, mag, g1, g2, comment)
        #Moment.__init__(self)


class MOMENT2(Load2):
    """
    Defines a static concentrated moment at a grid point by specification
    of a magnitude and four grid points that determine the direction.

    +---------+-----+---+---+----+----+----+----+
    |    1    |  2  | 3 | 4 |  5 |  6 |  7 |  8 |
    +=========+=====+===+===+====+====+====+====+
    | MOMENT2 | SID | G | M | G1 | G2 | G3 | G4 |
    +---------+-----+---+---+----+----+----+----+
    """
    type = 'MOMENT2'
    def __init__(self, sid, node, mag, g1, g2, g3, g4, comment=''):
        Load2.__init__(self, sid, node, mag, g1, g2, g3, g4, comment)


class GMLOAD(Load):
    """
    Defines a static concentrated force at a grid point by specification of a
    magnitude and two grid points that determine the direction.
    """
    type = 'GMLOAD'

    def __init__(self, sid, normal, entity, entity_id, method,
                 load_magnitudes, cid=0, comment=''):
        """Creates a GMLOAD object"""
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
        self.cid_ref = None

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
        return GMLOAD(sid, normal, entity, entity_id, method,
                      load_magnitudes, cid=cid, comment=comment)

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
        msg = ', which is required by GMLOAD sid=%s' % self.sid
        self.cid_ref = model.Coord(self.Cid(), msg=msg)
        #self.node = model.Node(self.node, msg=msg)
        #self.g1 = model.Node(self.g1, msg=msg)
        #self.g2 = model.Node(self.g2, msg=msg)
        #self.xyz = self.g2.get_position() - self.g1.get_position()
        #normalize(self, msg)

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by GMLOAD sid=%s' % self.sid
        self.cid_ref = model.safe_coord(self.Cid(), self.sid, xref_errors, msg=msg)

    def uncross_reference(self):
        self.cid = self.Cid()
        self.cid_ref = None

    def Cid(self):
        if self.cid_ref is not None:
            return self.cid_ref.cid
        return self.cid

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
        # type: (int, bool) -> str
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)
        """
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class PLOAD(Load):
    """
    Static Pressure Load

    Defines a uniform static pressure load on a triangular or quadrilateral surface
    comprised of surface elements and/or the faces of solid elements.

    +-------+-----+------+----+----+----+----+
    |   1   |  2  |  3   | 4  | 5  | 6  | 7  |
    +=======+=====+======+====+====+====+====+
    | PLOAD | SID |  P   | G1 | G2 | G3 | G4 |
    +-------+-----+------+----+----+----+----+
    | PLOAD |  1  | -4.0 | 16 | 32 | 11 |    |
    +-------+-----+------+----+----+----+----+
    """
    type = 'PLOAD'

    def __init__(self, sid, pressure, nodes, comment=''):
        """
        Creates a PLOAD card, which defines a uniform pressure load on a
        shell/solid face or arbitrarily defined quad/tri face.

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
        if nodes[-1] == 0:
            nodes = list(nodes)
            nodes.pop()
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

    def safe_cross_reference(self, model, safe_coord):
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
        # type: (int, bool) -> str
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)
        """
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class PLOAD1(Load):
    """
    Applied Load on CBAR, CBEAM or CBEND Elements

    Defines concentrated, uniformly distributed, or linearly distributed
    applied loads to the CBAR or CBEAM elements at user-chosen points
    along the axis. For the CBEND element, only distributed loads over
    an entire length may be defined.

    +--------+-----+------+------+-------+-----+-------+-----+-------+
    |   1    |  2  |  3   |  4   |   5   |  6  |   7   |  8  |   9   |
    +========+=====+======+======+=======+=====+=======+=====+=======+
    | PLOAD1 | SID | EID  | TYPE | SCALE | X1  |  P1   |  X2 |  P2   |
    +--------+-----+------+------+-------+-----+-------+-----+-------+
    | PLOAD1 | 25  | 1065 |  MY  | FRPR  | 0.2 | 2.5E3 | 0.8 | 3.5E3 |
    +--------+-----+------+------+-------+-----+-------+-----+-------+
    """
    type = 'PLOAD1'
    valid_types = ['FX', 'FY', 'FZ', 'FXE', 'FYE', 'FZE',
                   'MX', 'MY', 'MZ', 'MXE', 'MYE', 'MZE']

    # LE: length-based; FR: fractional; PR:projected
    valid_scales = ['LE', 'FR', 'LEPR', 'FRPR']

    def __init__(self, sid, eid, load_type, scale, x1, p1, x2=None, p2=None, comment=''):
        """
        Creates a PLOAD1 card, which may be applied to a CBAR/CBEAM

        Parameters
        ----------
        sid : int
            load id
        eid : int
            element to apply the load to
        load_type : str
            type of load that's applied
            valid_types = {FX, FY, FZ, FXE, FYE, FZE,
                           MX, MY, MZ, MXE, MYE, MZE}
        scale : str
            Determines scale factor for X1, X2.
            {LE, FR, LEPR, FRPR}
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
        self.load_type = load_type
        self.scale = scale
        self.x1 = x1
        self.p1 = p1
        self.x2 = x2
        self.p2 = p2
        self.eid_ref = None

    @property
    def Type(self):
        return self.load_type

    @Type.setter
    def Type(self, load_type):
        self.load_type = load_type

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
        load_type = string(card, 3, 'Type ("%s")' % '",  "'.join(cls.valid_types))
        scale = string(card, 4, 'scale ("%s")' % '", "'.join(cls.valid_scales))
        x1 = double(card, 5, 'x1')
        p1 = double(card, 6, 'p1')
        x2 = double_or_blank(card, 7, 'x2', x1)
        p2 = double_or_blank(card, 8, 'p2', p1)
        assert len(card) <= 9, 'len(PLOAD1 card) = %i\ncard=%s' % (len(card), card)
        return PLOAD1(sid, eid, load_type, scale, x1, p1, x2, p2, comment=comment)

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
        load_type = data[2]
        scale = data[3]
        x1 = data[4]
        p1 = data[5]
        x2 = data[6]
        p2 = data[7]
        load_type = cls.valid_types[load_type - 1]
        scale = cls.valid_scales[scale - 1]
        return PLOAD1(sid, eid, load_type, scale, x1, p1, x2, p2, comment=comment)

    def validate(self):
        if self.load_type not in self.valid_types:
            msg = '%s is an invalid type on the PLOAD1 card; valid_types=[%s]' % (
                self.load_type, ', '.join(self.valid_types).rstrip(', '))
            raise RuntimeError(msg)
        if self.scale not in self.valid_scales:
            msg = '%s is an invalid scale on the PLOAD1 card; valid_scales=[%s]' % (
                self.scale, ', '.join(self.valid_scales).rstrip(', '))
            raise RuntimeError(msg)

        assert 0.0 <= self.x1 <= self.x2, '0.0 <= x1 <= x2 -> x1=%s x2=%s' % (self.x1, self.x2)
        if self.scale in ['FR', 'FRPR']:
            assert self.x1 <= 1.0, 'x1=%r' % self.x1
            assert self.x2 <= 1.0, 'x2=%r' % self.x2
        if self.scale not in self.valid_scales:
            msg = '%s is an invalid scale on the PLOAD1 card; valid_scales=[%s]' % (
                self.scale, ', '.join(self.valid_scales))
            raise RuntimeError(msg)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by PLOAD1 sid=%s' % self.sid
        self.eid_ref = model.Element(self.eid, msg=msg)

    def safe_cross_reference(self, model, safe_coord):
        return self.cross_reference(model)

    def uncross_reference(self):
        self.eid = self.Eid()
        self.eid_ref = None

    def get_loads(self):
        return [self]

    def Eid(self):
        if self.eid_ref is not None:
            return self.eid_ref.eid
        return self.eid

    def raw_fields(self):
        list_fields = ['PLOAD1', self.sid, self.Eid(), self.load_type, self.scale,
                       self.x1, self.p1, self.x2, self.p2]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        # type: (int, bool) -> str
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)
        """
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
        self.eids_ref = None

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
        msg = ', which is required by PLOAD2 sid=%s' % self.sid
        self.eids_ref = model.Elements(self.eids, msg=msg)

    def safe_cross_reference(self, model, safe_coord):
        return self.cross_reference(model)

    def uncross_reference(self):
        self.eids = self.element_ids
        self.eids_ref = None

    @property
    def element_ids(self):
        if self.eids_ref is not None:
            eids = [elem.eid for elem in self.eids_ref]
        else:
            eids = self.eids
        return self.eids

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
        # type: (int, bool) -> str
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)
        """
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)

#def PLOAD4_func(self, sid, eids, pressures,
                #g1=None, g34=None, cid=0, nvector=None, surf_or_line='SURF',
                #line_load_dir='NORM', comment=''):
    #"""
    #Creates a PLOAD4 card

    #Solid Format
    #============
    #Defines a pressure load on a face of a CHEXA, CPENTA, or CTETRA element.

    #+--------+-----+-----+----+----+------+------+------+-------+
    #|   1    |  2  |  3  |  4 |  5 |  6   |   7  |   8  |   9   |
    #+========+=====+=====+====+====+======+======+======+=======+
    #| PLOAD4 | SID | EID | P1 | P2 |  P3  |  P4  |  G1  | G3/G4 |
    #+--------+-----+-----+----+----+------+------+------+-------+
    #|        | CID | N1  | N2 | N3 | SORL | LDIR |      |       |
    #+--------+-----+-----+----+----+------+------+------+-------+

    #Shell Format
    #============
    #Defines a pressure load on a face of a CTRIA3, CTRIA6, CTRIAR,
    #CQUAD4, CQUAD8, or CQUADR element.
    #+--------+-----+-----+----+----+------+------+------+-------+
    #|   1    |  2  |  3  |  4 |  5 |  6   |   7  |   8  |   9   |
    #+========+=====+=====+====+====+======+======+======+=======+
    #| PLOAD4 | SID | EID | P1 | P2 |  P3  |  P4  | THRU | EID2  |
    #+--------+-----+-----+----+----+------+------+------+-------+
    #|        | CID | N1  | N2 | N3 | SORL | LDIR |      |       |
    #+--------+-----+-----+----+----+------+------+------+-------+

    #.. warning:: NX does not support SORL and LDIR, MSC does
    #"""
    #if g34 is None:
        #return PLOAD4Solid(
            #sid, eids, pressures,
            #g1=None, g34=None, cid=0, nvector=None, surf_or_line='SURF',
            #line_load_dir='NORM', comment='')
    #return PLOAD4Shell(
        #sid, eids, pressures, cid=0, nvector=None, surf_or_line='SURF',
        #line_load_dir='NORM', comment='')


#class PLOAD4Shell(PLOAD4):
    #def __init__(self, sid, eids, pressures, g1=None, g34=None, cid=0,
                 #nvector=None, surf_or_line='SURF',
                 #line_load_dir='NORM', comment=''):
        #PLOAD4.__init__(self, sid, eids, pressures, g1=None, g34=None,
                        #cid=0, nvector=None,
                        #surf_or_line='SURF',
                        #line_load_dir='NORM',
                        #comment='')
#class PLOAD4Shell(PLOAD4):
    #def __init__(self, sid, eids, pressures, g1=None, g34=None, cid=0,
                 #nvector=None, surf_or_line='SURF',
                 #line_load_dir='NORM', comment=''):
        #PLOAD4.__init__(self, sid, eids, pressures, g1=g1, g34=g34,
                        #cid=cid, nvector=nvector,
                        #surf_or_line=surf_or_line,
                        #line_load_dir=line_load_dir,
                        #comment=comment)

class PLOAD4(Load):
    """
    ``Solid Format``

    Defines a pressure load on a face of a CHEXA, CPENTA, or CTETRA element.

    +--------+-----+-----+----+----+------+------+------+-------+
    |   1    |  2  |  3  |  4 |  5 |  6   |   7  |   8  |   9   |
    +========+=====+=====+====+====+======+======+======+=======+
    | PLOAD4 | SID | EID | P1 | P2 |  P3  |  P4  |  G1  | G3/G4 |
    +--------+-----+-----+----+----+------+------+------+-------+
    |        | CID | N1  | N2 | N3 | SORL | LDIR |      |       |
    +--------+-----+-----+----+----+------+------+------+-------+

    ``Shell Format``

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

    def __init__(self, sid, eids, pressures, g1, g34,
                 cid=0, nvector=None, surf_or_line='SURF',
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
            the coordinate system for nvector
        nvector : (3, ) float ndarray
           blank : load acts normal to the face
           float : the local pressure vector
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
            nvector = np.full((3, ), np.nan, dtype='float64')
        else:
            nvector = np.asarray(nvector, dtype='float64')

        if comment:
            self.comment = comment
        if isinstance(eids, integer_types):
            eids = [eids]
        if isinstance(pressures, float_types):
            pressures = [pressures] * 4
        # TODO: handle default pressure as input

        self.sid = sid

        # these can be greater than 1 if it's a shell (not a solid)
        self.eids = eids
        self.pressures = np.asarray(pressures, dtype='float64')
        if surf_or_line == 'SURF':
            inan = np.isnan(self.pressures)
            self.pressures[inan] = pressures[0]

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
        #self.eid_ref = None
        self.g1_ref = None
        self.g34_ref = None
        self.cid_ref = None
        self.eids_ref = None

    def validate(self):
        if self.surf_or_line not in ['SURF', 'LINE']:
            raise RuntimeError('PLOAD4; sid=%s surf_or_line=%r' % (self.sid, self.surf_or_line))
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
            double_or_blank(card, 4, 'p2'),
            double_or_blank(card, 5, 'p3'),
            double_or_blank(card, 6, 'p4')]

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

        # If both (CID, N1, n2, N3) and LDIR are blank, then the default is
        # LDIR=NORM.
        cid = integer_or_blank(card, 9, 'cid')
        n1 = double_or_blank(card, 10, 'N1')
        n2 = double_or_blank(card, 11, 'N2')
        n3 = double_or_blank(card, 12, 'N3')
        is_not_none = [value for value in [cid, n1, n2, n3]
                       if value is not None]

        # set the nvector defaults only if one of them have a value
        # otherwise, they're all None
        if is_not_none:
            if cid is None:
                cid = 0
            if n1 is None:
                n1 = 0.
            if n2 is None:
                n2 = 0.
            if n3 is None:
                n3 = 0.
        nvector = array([n1, n2, n3])

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

        eids = [eid]
        if data[7] is None:
            surf_or_line = 'SURF'
            assert data[8] is None, data
            line_load_dir = 'NORM'
        else:
            surf_or_line = data[7]
            line_load_dir = data[8]
        pload4 = PLOAD4(sid, eids, pressures, g1, g34, cid, nvector,
                        surf_or_line, line_load_dir, comment=comment)
        assert sid < 10000000, pload4
        assert cid < 10000000, pload4
        return pload4

    def get_loads(self):
        return [self]

    def Cid(self):
        """gets the coordinate system object"""
        if self.cid_ref is not None:
            return self.cid_ref.cid
        return self.cid

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by PLOAD4 sid=%s' % self.sid
        if self.cid is not None:
            self.cid_ref = model.Coord(self.cid, msg=msg)
        if self.g1 is not None:
            self.g1_ref = model.Node(self.g1, msg=msg + '; g1')
        if self.g34 is not None:
            self.g34_ref = model.Node(self.g34, msg=msg + '; g34')
        if self.eids:
            self.eids_ref = model.Elements(self.eids, msg=msg)

    def safe_cross_reference(self, model, xref_errors, debug=True):
        msg = ', which is required by PLOAD4 sid=%s' % self.sid
        #self.eid = model.Element(self.eid, msg=msg)
        if self.cid is not None:
            self.cid_ref = model.safe_coord(self.cid, self.sid, xref_errors, msg=msg)

        #self.eid_ref = self.eid
        if self.g1 is not None:
            try:
                self.g1_ref = model.Node(self.g1, msg=msg)
            except KeyError:
                model.log.warning('Could not find g1=%s%s' % (self.g1, msg))

        if self.g34 is not None:
            try:
                self.g34_ref = model.Node(self.g34, msg=msg)
            except KeyError:
                model.log.warning('Could not find g34=%s%s' % (self.g34, msg))

        #if self.eids:
        msgia = 'Could not find element=%%s%s\n' % msg
        self.eids_ref, msgi = model.safe_get_elements(self.eids, msg=msgia)
        if msgi:
            model.log.warning(msgi.rstrip())

    def uncross_reference(self):
        self.cid = self.Cid()
        if self.g1 is not None:
            self.g1 = self.G1()
        if self.g34 is not None:
            self.g34 = self.G34()
        self.eids = self.element_ids
        self.g1_ref = None
        self.g34_ref = None
        self.cid_ref = None
        self.eids_ref = None

    def G1(self):
        if self.g1_ref is not None:
            return self.g1_ref.nid
        return self.g1

    def G34(self):
        if self.g34_ref is not None:
            return self.g34_ref.nid
        return self.g34

    @property
    def node_ids(self):
        node_ids = [self.G1(), self.G34()]
        return node_ids

    def get_element_ids(self, eid=None):
        if self.eids_ref is not None:
            try:
                eids = [eid_ref.eid for eid_ref in self.eids_ref]
            except AttributeError:
                eids = []
                for eid_ref in self.eids_ref:
                    if isinstance(eid_ref, integer_types):
                        # Nastran is NOT OK with elements that don't actually exist in the PLOAD4
                        # we do this for safe_cross_reference
                        eids.append(eid)
                    else:
                        eids.append(eid_ref.eid)
        else:
            eids = self.eids
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
        if cid is not None or not np.all(np.isnan(self.nvector)):
            n1, n2, n3 = self.nvector
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
        list_fields = ['PLOAD4', self.sid, eid, p1, p2, p3, p4]

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
        if cid is not None or not np.all(np.isnan(self.nvector)):
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
        # type: (int, bool) -> str
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)
        """
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

def update_pload4_vector_for_surf(pload4, normal, cid):
    """helper method"""
    if cid is None and np.all(np.isnan(pload4.nvector)):
        # element surface normal
        pass
    else:
        if norm(pload4.nvector) != 0.0 and cid == 0:
            normal = pload4.nvector / np.linalg.norm(pload4.nvector)
        else:
            raise NotImplementedError('cid=%r nvector=%s on a PLOAD4 is not supported\n%s' % (
                cid, pload4.nvector, str(pload4)))
    return normal


class PLOADX1(BaseCard):
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
        BaseCard.__init__(self)
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
        self.eid_ref = None
        self.ga_ref = None
        self.gb_ref = None

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
        nids = [ga, gb]
        return PLOADX1(sid, eid, pa, nids, pb=pb, theta=theta, comment=comment)

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
        nids = [ga, gb]
        return PLOADX1(sid, eid, pa, nids, pb=pb, theta=theta, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by PLOADX1 lid=%s' % self.sid
        self.eid_ref = model.Element(self.eid, msg=msg)
        self.ga_ref = model.Node(self.ga, msg=msg)
        self.gb_ref = model.Node(self.gb, msg=msg)

    @property
    def node_ids(self):
        return [self.Ga(), self.Gb()]

    @property
    def nodes(self):
        return [self.ga, self.gb]

    @property
    def nodes_ref(self):
        return [self.ga_ref, self.gb_ref]

    def safe_cross_reference(self, model, safe_coord):
        return self.cross_reference(model)

    def uncross_reference(self):
        self.eid = self.Eid()
        self.ga = self.Ga()
        self.gb = self.Gb()
        self.eid_ref = None
        self.ga_ref = None
        self.gb_ref = None

    def Eid(self):
        if self.eid_ref is not None:
            return self.eid_ref.eid
        return self.eid

    def Ga(self):
        if self.ga_ref is not None:
            return self.ga_ref.nid
        return self.ga

    def Gb(self):
        if self.gb_ref is not None:
            return self.gb_ref.nid
        return self.gb

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
        # type: (int, bool) -> str
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)
        """
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


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


#def normalize(self):
    #"""
    #adjust the vector to a unit length
    #scale up the magnitude of the vector
    #"""
    #if self.mag != 0.0:  # enforced displacement
        #norm_xyz = norm(self.xyz)
        ##mag = self.mag*norm_xyz
        #self.mag *= norm_xyz
        #self.xyz /= norm_xyz
