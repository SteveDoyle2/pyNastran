"""
All constraint cards are defined in this file.  This includes:

* Constraint
 * SUPORT
 * SUPORT1
 * SPC
 * SPC1
 * SPCAX
 * MPC
 * GMSPC
 * ConstraintADD
  * SPCADD
  * MPCADD

The ConstraintObject contain multiple constraints.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from itertools import count
from six import iteritems
from six.moves import zip, range

from pyNastran.utils import integer_types
from pyNastran.bdf.cards.base_card import BaseCard, _node_ids, expand_thru
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, components,
    components_or_blank, string)
from pyNastran.bdf.field_writer_8 import print_card_8, print_float_8
from pyNastran.bdf.field_writer_16 import print_float_16, print_card_16
from pyNastran.bdf.field_writer_double import print_scientific_double


class Constraint(BaseCard):
    def __init__(self):
        pass

    def raw_fields(self):
        fields = [self.type, self.conid]
        return fields

    def _nodeIDs(self, nodes=None, allow_empty_nodes=False, msg=''):
        """returns nodeIDs for repr functions"""
        return _node_ids(self, nodes, allow_empty_nodes, msg)


class SUPORT1(Constraint):
    """
    +---------+-----+-----+----+-----+----+-----+----+
    |    1    |  2  |  3  |  4 |  5  | 6  |  7  | 8  |
    +=========+=====+=====+====+=====+====+=====+====+
    | SUPORT1 | SID | ID1 | C1 | ID2 | C2 | ID3 | C3 |
    +---------+-----+-----+----+-----+----+-----+----+
    | SUPORT1 |  1  |  2  | 23 |  4  | 15 |  5  |  0 |
    +---------+-----+-----+----+-----+----+-----+----+
    """
    type = 'SUPORT1'

    def __init__(self, conid, IDs, Cs, comment=''):
        Constraint.__init__(self)
        if comment:
            self._comment = comment
        self.conid = conid
        self.IDs = IDs
        self.Cs = Cs
        assert len(self.IDs) > 0
        assert len(self.IDs) == len(self.Cs)

    @classmethod
    def add_card(cls, card, comment=''):
        conid = integer(card, 1, 'conid')  # really a support id sid

        nfields = len(card)
        assert len(card) > 2
        nterms = int((nfields - 1.) / 2.)
        n = 1
        IDs = []
        Cs = []
        for i in range(nterms):
            nstart = 2 + 2 * i
            ID = integer(card, nstart, 'ID%s' % n)
            C = components_or_blank(card, nstart + 1, 'component%s' % n, '0')
            IDs.append(ID)
            Cs.append(C)
            n += 1
        return SUPORT1(conid, IDs, Cs, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        conid = data[0]
        assert (len(data) - 1) % 2 == 0, data
        IDs = []
        Cs = []
        for i in range(1, len(data), 2):
            ID = data[i]
            C = data[i+1]
            IDs.append(ID)
            Cs.append(C)
        return SUPORT1(conid, IDs, Cs, comment=comment)

    def add_suport1_to_set(self, suport1):
        assert self.conid == suport1.conid, 'SUPORT1 conid=%s new_conid=%s; they must be the same' % (self.conid, suport1.conid)
        comment = self.comment + suport1.comment
        if comment:
            self._comment = comment
        self.IDs += suport1.IDs
        self.Cs += suport1.Cs

    @property
    def node_ids(self):
        msg = ', which is required by SUPORT1'
        return self._nodeIDs(nodes=self.IDs, allow_empty_nodes=True, msg=msg)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by SUPORT1'
        self.IDs = model.Nodes(self.IDs, allow_empty_nodes=True, msg=msg)
        self.IDs_ref = self.IDs

    def uncross_reference(self):
        self.IDs = self.node_ids
        del self.IDs_ref

    def raw_fields(self):
        fields = ['SUPORT1', self.conid]
        for ID, c in zip(self.node_ids, self.Cs):
            fields += [ID, c]
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class SUPORT(Constraint):
    """
    +---------+-----+-----+-----+-----+-----+-----+-----+----+
    |    1    |  2  |  3  |  4  |  5  |  6  |  7  |  8  | 9  |
    +=========+=====+=====+=====+=====+=====+=====+=====+====+
    | SUPORT  | ID1 | C1  | ID2 |  C2 | ID3 | C3  | ID4 | C4 |
    +---------+-----+-----+-----+-----+-----+-----+-----+----+
    """
    type = 'SUPORT'

    def __init__(self, IDs, Cs, comment=''):
        Constraint.__init__(self)
        if comment:
            self._comment = comment
        self.IDs = IDs ## TODO:  IDs reference nodes???
        self.Cs = Cs

    def validate(self):
        assert len(self.IDs) > 0
        assert len(self.IDs) == len(self.Cs)

    @classmethod
    def add_card(cls, card, comment=''):
        # TODO: remove fields...
        #fields = card.fields(1)

        nfields = len(card)
        assert len(card) > 1, card
        nterms = int(nfields / 2.)
        n = 1
        IDs = []
        Cs = []
        for i in range(nterms):
            nstart = 1 + 2 * i
            ID = integer(card, nstart, 'ID%s' % n)
            C = components_or_blank(card, nstart + 1, 'component%s' % n, '0')
            IDs.append(ID)
            Cs.append(C)
            n += 1
        return SUPORT(IDs, Cs, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        fields = data
        IDs = []
        Cs = []
        for i in range(0, len(fields), 2):
            IDs.append(fields[i])
            Cs.append(fields[i + 1])
        return SUPORT(IDs, Cs, comment=comment)

    @property
    def node_ids(self):
        msg = ', which is required by %s' % self.type
        return self._nodeIDs(nodes=self.IDs, allow_empty_nodes=True, msg=msg)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by %s' % self.type
        self.IDs = model.Nodes(self.IDs, allow_empty_nodes=True, msg=msg)
        self.IDs_ref = self.IDs

    def uncross_reference(self):
        self.IDs = self.node_ids
        del self.IDs_ref

    def raw_fields(self):
        fields = [self.type]
        for ID, c in zip(self.node_ids, self.Cs):
            fields += [ID, c]
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class SESUP(SUPORT):
    type = 'SESUP'

    def __init__(self):
        SUPORT.__init__(self)


class MPC(Constraint):
    """
    +-----+-----+----+----+-----+----+----+----+-----+
    |  1  |  2  |  3 |  4 |  5  |  6 |  7 |  8 |  9  |
    +=====+=====+====+====+=====+====+====+====+=====+
    | MPC | SID | G1 | C1 |  A1 | G2 | C2 | A2 |     |
    +-----+-----+----+----+-----+----+----+----+-----+
    |     |  G3 | C3 | A3 | ... |    |    |    |     |
    +-----+-----+----+----+-----+----+----+----+-----+
    """
    type = 'MPC'

    def __init__(self, conid, gids, constraints, enforced, comment=''):
        Constraint.__init__(self)
        if comment:
            self._comment = comment
        #: Set identification number. (Integer > 0)
        self.conid = conid
        #: Identification number of grid or scalar point. (Integer > 0)
        self.gids = gids
        #: Component number. (Any one of the Integers 1 through 6 for grid
        #: points; blank or zero for scalar points.)
        self.constraints = constraints
        #: Coefficient. (Real; Default = 0.0 except A1 must be nonzero.)
        self.enforced = enforced

    @classmethod
    def add_card(cls, card, comment=''):
        conid = integer(card, 1, 'conid')
        gids = []
        constraints = []
        enforced = []

        fields = card.fields(0)
        nfields = len(fields)

        i = 1
        for ifield in range(2, nfields, 8):
            grid = integer(card, ifield, 'G%i' % i)
            component = components_or_blank(card, ifield + 1, 'constraint%i' % i, 0)  # scalar point
            if i == 1:
                enforcedi = double(card, ifield + 2, 'enforced%i' % i)
                if enforcedi == 0.0:
                    raise RuntimeError('enforced1 must be nonzero; enforcedi=%r' % enforcedi)
            else:
                enforcedi = double_or_blank(card, ifield + 2, 'enforced%i' % i, 0.0)
            gids.append(grid)
            constraints.append(component)
            enforced.append(enforcedi)
            i += 1

            if ifield + 4 > nfields and i != 2:
                # if G2 is empty (it's ifield+4 because nfields is length based
                # and not loop friendly)
                break
            grid = integer(card, ifield + 3, 'G%i' % i)
            component = components_or_blank(card, ifield + 4, 'constraint%i' % i, 0)  # scalar point
            enforcedi = double_or_blank(card, ifield + 5, 'enforced%i' % i)
            gids.append(grid)
            constraints.append(component)
            enforced.append(enforcedi)
            i += 1
        return MPC(conid, gids, constraints, enforced, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        sid,
        msg = '%s has not implemented data parsing' % cls.type
        raise NotImplementedError(msg)
        return MPC(conid, gids, constraints, enforced, comment=comment)

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        msg = ', which is required by MPC=%s' % self.conid
        return self._nodeIDs(nodes=self.gids, allow_empty_nodes=True, msg=msg)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by MPC=%s' % self.conid
        self.gids = model.Nodes(self.gids, allow_empty_nodes=True, msg=msg)
        self.gids_ref = self.gids

    def uncross_reference(self):
        self.gids = self.node_ids
        del self.gids_ref

    def raw_fields(self):  # MPC
        fields = ['MPC', self.conid]
        for i, gid, constraint, enforced in zip(count(), self.node_ids, self.constraints, self.enforced):
            fields += [gid, constraint, enforced]
            if i % 2 == 1 and i > 0:
                fields.append(None)
                fields.append(None)
        return fields

    def write_card(self, size=8, is_double=False):
        """see BaseCard.write_card``"""
        if size == 8:
            return self.write_card_8()
        else:
            return self.write_card_16(is_double)

    def write_card_8(self):
        msg = 'MPC     %8s' % self.conid
        grids, constraints, enforceds = self.node_ids, self.constraints, self.enforced
        for i, grid, component, enforced in zip(count(), grids, constraints, enforceds):
            msg += '%8i%8s%8s' % (grid, component, print_float_8(enforced))
            if i % 2 == 1 and i > 0:
                msg += '\n%8s%8s' % ('', '')
        return self.comment + msg.rstrip() + '\n'

    def write_card_16(self, is_double=False):
        # TODO: we're sure MPCs support double precision?
        msg = 'MPC*    %16s' % self.conid
        grids, constraints, enforceds = self.node_ids, self.constraints, self.enforced
        if is_double:
            for i, grid, component, enforced in zip(count(), grids, constraints, enforceds):
                if i == 0:
                    msg += '%16i%16s%16s\n' % (
                        grid, component, print_scientific_double(enforced))
                elif i % 2 == 1:
                    msg += '%-8s%16i%16s%16s\n' % (
                        '*', grid, component, print_scientific_double(enforced))
                else:
                    msg += '%-8s%16s%16i%16s%16s\n' % (
                        '*', '', grid, component, print_scientific_double(enforced))
        else:
            for i, grid, component, enforced in zip(count(), grids, constraints, enforceds):
                if i == 0:
                    msg += '%16i%16s%16s\n' % (grid, component, print_float_16(enforced))
                elif i % 2 == 1:
                    msg += '%-8s%16i%16s%16s\n' % (
                        '*', grid, component, print_float_16(enforced))
                else:
                    msg += '%-8s%16s%16i%16s%16s\n' % (
                        '*', '', grid, component, print_float_16(enforced))
        if i % 2 == 0:
            msg += '*'
        return self.comment + msg.rstrip() + '\n'


class SPC(Constraint):
    """
    Defines enforced displacement/temperature (static analysis)
    velocity/acceleration (dynamic analysis).

     +-----+-----+----+----+------+----+----+----+
     | SPC | SID | G1 | C1 |  D1  | G2 | C2 | D2 |
     +-----+-----+----+----+------+----+----+----+
     | SPC |  2  | 32 | 3  | -2.6 |  5 |    |    |
     +-----+-----+----+----+------+----+----+----+
    """
    type = 'SPC'

    def __init__(self, conid, gids, constraints, enforced, comment=''):
        Constraint.__init__(self)
        if comment:
            self._comment = comment
        self.conid = conid
        self.gids = gids
        self.constraints = constraints
        self.enforced = enforced

    @classmethod
    def add_card(cls, card, comment=''):
        conid = integer(card, 1, 'sid')
        if card.field(5) in [None, '']:
            gids = [integer(card, 2, 'G1'),]
            constraints = [components_or_blank(card, 3, 'C1', 0)]
            enforced = [double_or_blank(card, 4, 'D1', 0.0)]
        else:
            gids = [
                integer(card, 2, 'G1'),
                integer_or_blank(card, 5, 'G2'),
            ]
            # :0 if scalar point 1-6 if grid
            constraints = [components_or_blank(card, 3, 'C1', 0),
                           components_or_blank(card, 6, 'C2', 0)]
            enforced = [double_or_blank(card, 4, 'D1', 0.0),
                        double_or_blank(card, 7, 'D2', 0.0)]
        return SPC(conid, gids, constraints, enforced, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        conid = data[0]
        gids = [data[1]]
        constraint = data[2]
        assert 0 <= constraint <= 123456, data
        enforced = [data[3]]
        assert conid > 0, data
        assert gids[0] > 0, data
        constraint_str = str(constraint)
        assert len(constraint_str) <= 6, data
        constraints = [constraint_str]
        #if constraints[0] == 0:
            #constraints[0] = 0
        #if constraints[0] == 16:
            #constraints[0] = '16'
        #else:
            #raise RuntimeError('SPC; constraints=%s data=%s' % (constraints, data))
        #assert 0 < constraints[0] > 1000, data
        return SPC(conid, gids, constraints, enforced, comment=comment)

    #def get_node_dofs(self, model):
        #pass

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        msg = ', which is required by %s=%s' % (self.type, self.conid)
        return self._nodeIDs(nodes=self.gids, allow_empty_nodes=True, msg=msg)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by %s=%s' % (self.type, self.conid)
        self.gids = model.Nodes(self.gids, allow_empty_nodes=True, msg=msg)
        self.gids_ref = self.gids

    def uncross_reference(self):
        self.gids = self.node_ids
        del self.gids_ref

    def raw_fields(self):
        fields = ['SPC', self.conid]
        for (gid, constraint, enforced) in zip(self.node_ids, self.constraints,
                                               self.enforced):
            fields += [gid, constraint, enforced]
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class GMSPC(Constraint):
    type = 'GMSPC'

    def __init__(self, conid, component, entity, entity_id, comment=''):
        Constraint.__init__(self)
        if comment:
            self._comment = comment
        self.conid = conid
        self.component = component
        self.entity = entity
        self.entity_id = entity_id

    @classmethod
    def add_card(cls, card, comment=''):
        conid = integer(card, 1, 'sid')
        component = components(card, 2, 'components')
        entity = string(card, 3, 'entity')
        entity_id = integer(card, 4, 'entity_id')
        return GMSPC(conid, component, entity, entity_id, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        raise NotImplementedError('GMSPC')

    def cross_reference(self, model):
        """TODO: xref"""
        #msg = ', which is required by %s=%s' % (self.type, self.conid)
        pass

    def uncross_reference(self):
        pass

    def raw_fields(self):
        fields = ['GMSPC', self.conid, self.component, self.entity, self.entity_id]
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class SPCAX(Constraint):
    """
    Defines a set of single-point constraints or enforced displacements
    for conical shell coordinates.

     +-------+-----+-----+-----+----+-----+
     | SPCAX | SID | RID | HID |  C |  D  |
     +-------+-----+-----+-----+----+-----+
     | SPCAX |  2  |  3  |  4  | 13 | 6.0 |
     +-------+-----+-----+-----+----+-----+
    """
    type = 'SPCAX'

    def __init__(self, conid, rid, hid, c, d, comment=''):
        # defines everything :) at least until cross-referencing methods are
        # implemented
        Constraint.__init__(self)
        if comment:
            self._comment = comment
        #: Identification number of a single-point constraint set.
        self.conid = conid
        #: Ring identification number. See RINGAX entry.
        self.rid = rid
        #: Harmonic identification number. (Integer >= 0)
        self.hid = hid
        #: Component identification number. (Any unique combination of the
        #: Integers 1 through 6.)
        self.c = c
        #: Enforced displacement value
        self.d = d

    @classmethod
    def add_card(cls, card, comment=''):
        conid = integer(card, 1, 'conid')
        rid = integer(card, 2, 'rid')
        hid = integer(card, 3, 'hid')
        c = components(card, 4, 'c')
        d = double(card, 5, 'd')
        return SPCAX(conid, rid, hid, c, d, comment=comment)

    #@classmethod
    #def add_op2_data(cls, data, comment=''):
        #msg = '%s has not implemented data parsing' % cls.type
        #raise NotImplementedError(msg)

    def cross_reference(self, model):
        pass
        #msg = ', which is required by %s=%s' % (self.type, self.conid)
        #self.rid = model.ring[self.rid]
        #self.hid = model.harmonic[self.hid]

    def uncross_reference(self):
        pass

    def raw_fields(self):
        fields = ['SPCAX', self.conid, self.rid, self.hid, self.c, self.d]
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class SPC1(Constraint):
    """
      +------+-----+----+----+--------+----+----+----+----+
      | SPC1 | SID | C  | G1 | G2     | G3 | G4 | G5 | G6 |
      +------+-----+----+----+--------+----+----+----+----+
      |      |  G7 | G8 | G9 | -etc.- |    |    |    |    |
      +------+-----+----+----+--------+----+----+----+----+

      +------+---+-----+--------+--------+--------+--------+--------+---+
      | SPC1 | 3 | 246 | 209075 | 209096 | 209512 | 209513 | 209516 |   |
      +------+---+-----+--------+--------+--------+--------+--------+---+
      | SPC1 | 3 |  2  |   1    |   3    |   10   |   9    |   6    | 5 |
      +------+---+-----+--------+--------+--------+--------+--------+---+
      |      | 2 |  8  |        |        |        |        |        |   |
      +------+---+-----+--------+--------+--------+--------+--------+---+

      +------+-----+-------+----+------+----+
      | SPC1 | SID | C     | G1 | THRU | G2 |
      +------+-----+-------+----+------+----+
      | SPC1 | 313 | 12456 | 6  | THRU | 32 |
      +------+-----+-------+----+------+----+
    """
    type = 'SPC1'

    def __init__(self, conid, constraints, nodes, comment=''):
        Constraint.__init__(self)
        if comment:
            self._comment = comment
        self.conid = conid
        self.constraints = constraints
        self.nodes = expand_thru(nodes)
        self.nodes.sort()

    @classmethod
    def add_card(cls, card, comment=''):
        conid = integer(card, 1, 'conid')
        constraints = components(card, 2, 'constraints')  # 246 = y; dx, dz dir
        nodes = card.fields(3)
        return SPC1(conid, constraints, nodes, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        conid = data[0]
        constraints = data[1]
        nodes = data[2:]
        assert conid > 0, data
        for nid in nodes:
            assert nid > 0, data
        return SPC1(conid, constraints, nodes, comment=comment)

    @property
    def node_ids(self):
        msg = ', which is required by SPC1; conid=%s' % self.conid
        return self._nodeIDs(self.nodes, allow_empty_nodes=True, msg=msg)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by SPC1; conid=%s' % self.conid
        self.nodes = model.Nodes(self.node_ids, allow_empty_nodes=True, msg=msg)
        self.nodes_ref = self.nodes

    def uncross_reference(self):
        self.nodes = self.node_ids
        del self.nodes_ref

    def raw_fields(self):
        fields = ['SPC1', self.conid, self.constraints] + self.node_ids
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class ConstraintADD(Constraint):
    def __init__(self):
        Constraint.__init__(self)


class SPCADD(ConstraintADD):
    """
    Defines a single-point constraint set as a union of single-point constraint
    sets defined on SPC or SPC1 entries.

    +--------+----+----+-----+
    |    1   | 2  |  3 |  4  |
    +========+====+====+=====+
    | SPCADD  | 2 | 1 |   3  |
    +---------+---+---+------+
    """
    type = 'SPCADD'

    def __init__(self, conid, sets, comment=''):
        ConstraintADD.__init__(self)
        if comment:
            self._comment = comment
        self.conid = conid
        self.sets = expand_thru(sets)
        self.sets.sort()

    @classmethod
    def add_card(cls, card, comment=''):
        conid = integer(card, 1, 'conid')
        sets = card.fields(2)
        return SPCADD(conid, sets, comment=comment)

    def organize_constraints(self, model):
        """
        Figures out magnitudes of the loads to be applied to the various nodes.
        This includes figuring out scale factors.
        """
        position_spcs = []
        types_found = ['SPCADD']
        (scale_factors, loads) = self.get_reduced_constraints()
        return (types_found, position_spcs)

    @property
    def spc_ids(self):
        spc_ids = []
        for spc in self.sets:
            if isinstance(spc, integer_types):
                spc_ids.append(spc)
            elif isinstance(spc, list):
                spc_ids.append(spc[0].conid)
            else:
                raise TypeError('type=%s; spc=\n%s' % (type(spc), spc))
        return spc_ids

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by SPCADD=%s' % self.conid
        for i, spc in enumerate(self.sets):
            self.sets[i] = model.SPC(spc, msg=msg)
        self.sets_ref = self.sets

    def uncross_reference(self):
        self.sets = self.spc_ids
        del self.sets_ref

    def raw_fields(self):
        fields = ['SPCADD', self.conid] + self.spc_ids
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)

    def write_card_16(self, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_16(card)


class MPCADD(ConstraintADD):
    r"""
    Defines a multipoint constraint equation of the form
    :math:`\Sigma_j A_j u_j =0` where :math:`u_j` represents
    degree-of-freedom :math:`C_j` at grid or scalar point :math:`G_j`.

    +--------+----+----+-----+
    |    1   | 2  |  3 |  4  |
    +========+====+====+=====+
    | MPCADD | 2  |  1 |  3  |
    +--------+----+----+-----+
    """
    type = 'MPCADD'

    def __init__(self, conid, sets, comment=''):
        ConstraintADD.__init__(self)
        if comment:
            self._comment = comment
        self.conid = conid
        self.sets = expand_thru(sets)
        self.sets.sort()

    @classmethod
    def add_card(cls, card, comment=''):
        conid = integer(card, 1, 'conid')
        sets = card.fields(2)
        return MPCADD(conid, sets, comment=comment)

    #@classmethod
    #def add_op2_data(cls, data, comment=''):
        #msg = '%s has not implemented data parsing' % self.type
        #raise NotImplementedError(msg)
        #return MPCADD(conid, sets, comment=comment)

    @property
    def mpc_ids(self):
        mpc_ids = []
        for mpc in self.sets:
            if isinstance(mpc, integer_types):
                mpc_ids.append(mpc)
            else:
                mpc_ids.append(mpc.conid)
        return mpc_ids

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by MPCADD=%s' % self.conid
        for i, mpc in enumerate(self.sets):
            self.sets[i] = model.MPC(mpc, msg=msg)
        self.sets_ref = self.sets

    def raw_fields(self):
        fields = ['MPCADD', self.conid]  # +self.sets
        for set_id in self.sets:
            if isinstance(set_id, integer_types):
                fields.append(set_id)
            else:
                fields.append(set_id[0].conid)
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)

    def write_card_16(self, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_16(card)
