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
from six import iteritems, integer_types
from six.moves import zip, range
from itertools import count

from pyNastran.bdf.cards.baseCard import BaseCard, _node_ids, expand_thru
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, components, components_or_blank, string)
from pyNastran.bdf.field_writer_8 import print_card_8, print_float_8
from pyNastran.bdf.field_writer_16 import print_float_16, print_card_16
from pyNastran.bdf.field_writer_double import print_scientific_double


class ConstraintObject(object):
    def __init__(self):
        self.constraints = {}  # SPC, SPC1, SPCAX, etc...
        self.add_constraints = {}  # SPCADD

    def Add(self, ADD_Constraint):
        """Bad name, but adds an SPCADD or MPCADD"""
        key = ADD_Constraint.conid
        if key in self.add_constraints:
            print("already has key...key=%s\n%s" % (key, str(ADD_Constraint)))
        else:
            self.add_constraints[key] = ADD_Constraint

    def append(self, constraint):
        key = constraint.conid
        if key in self.constraints:
            #print("already has key...key=%s\n%s" %(key,str(constraint)))
            self.constraints[key].append(constraint)
        else:
            self.constraints[key] = [constraint]

    def _spc(self, spcID):
        """could be an spcadd/mpcadd or spc/mpc,spc1,spcd,spcax,suport1"""
        if spcID in self.add_constraints:
            return self.add_constraints[spcID]
        return self.constraints[spcID]

    def cross_reference(self, model):
        #return
        #add_constraints2 = {}
        for key, add_constraint in sorted(iteritems(self.add_constraints)):
            for i, spc_id in enumerate(add_constraint.sets):
                add_constraint.sets[i] = self._spc(spc_id)
            self.add_constraints[key] = add_constraint
        #self.add_constraints = add_constraints2

        constraints2 = {}
        for key, constraints in sorted(iteritems(self.constraints)):
            constraints2[key] = []
            for constraint in constraints:
                constraints2[key].append(constraint.cross_reference(model))
        self.constraints = constraints2

    def ConstraintID(self):
        if isinstance(self.conid, integer_types):
            return self.conid
        return self.conid.conid

    def get_constraint_ids(self):
        ids = []
        for key, constraints in sorted(iteritems(self.add_constraints)):
            con_id = constraints.ConID()
            ids.append(con_id)

        for key, constraints in sorted(iteritems(self.constraints)):
            for constraint in constraints:
                con_id = constraint.ConID()
                ids.append(con_id)

        ids = list(set(ids))
        ids.sort()
        return ids

    def __repr__(self):
        msg = ''
        # write the SPCADD/MPCADD cards
        for key, add_constraint in sorted(iteritems(self.add_constraints)):
            msg += str(add_constraint)

        for key, constraints in sorted(iteritems(self.constraints)):
            for constraint in constraints:
                msg += str(constraint)
        return msg


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
    | SUPORT1 | SID | ID1 | C1 | ID2 | C2 | ID3 | C3 |
    +---------+-----+-----+----+-----+----+-----+----+
    | SUPORT1 |  1  |  2  | 23 |  4  | 15 |  5  |  0 |
    +---------+-----+-----+----+-----+----+-----+----+
    """
    type = 'SUPORT1'

    def __init__(self):
        Constraint.__init__(self)
        self.IDs = []
        self.Cs = []

    def add_card(self, card, comment=''):
        if comment:
            self._comment = comment
        self.conid = integer(card, 1, 'conid')  # really a support id sid

        nfields = len(card)
        assert len(card) > 2
        nterms = int((nfields - 1.) / 2.)
        n = 1
        for i in range(nterms):
            nstart = 2 + 2 * i
            ID = integer(card, nstart, 'ID%s' % n)
            C = components_or_blank(card, nstart + 1, 'component%s' % n, '0')
            self.IDs.append(ID)
            self.Cs.append(C)
            n += 1
        assert len(self.IDs) > 0
        assert len(self.IDs) == len(self.Cs)

    def add_op2_data(self, data, comment=''):
        if comment:
            self._comment = comment
        msg = '%s has not implemented data parsing' % self.type
        raise NotImplementedError(msg)
        #assert len(self.IDs) > 0
        #assert len(self.IDs) == len(self.Cs)

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
    | SUPORT  | ID1 | C1  | ID2 |  C2 | ID3 | C3  | ID4 | C4 |
    +---------+-----+-----+-----+-----+-----+-----+-----+----+
    """
    type = 'SUPORT'

    def __init__(self):
        Constraint.__init__(self)
        self.IDs = [] ## TODO:  IDs reference nodes???
        self.Cs = []

    def add_card(self, card, comment=''):
        if comment:
            self._comment = comment
        # TODO: remove fields...
        fields = card.fields(1)

        nfields = len(card)
        assert len(card) > 1, card
        nterms = int(nfields / 2.)
        n = 1
        for i in range(nterms):
            nstart = 1 + 2 * i
            ID = integer(card, nstart, 'ID%s' % n)
            C = components_or_blank(card, nstart + 1, 'component%s' % n, '0')
            self.IDs.append(ID)
            self.Cs.append(C)
            n += 1
        assert len(self.IDs) > 0
        assert len(self.IDs) == len(self.Cs)

    def add_op2_data(self, data):
        fields = data
        for i in range(0, len(fields), 2):
            self.IDs.append(fields[i])
            self.Cs.append(fields[i + 1])
        assert len(self.IDs) > 0
        assert len(self.IDs) == len(self.Cs)

    @property
    def node_ids(self):
        msg = ', which is required by %s' % self.type
        return self._nodeIDs(nodes=self.IDs, allow_empty_nodes=True, msg=msg)

    def cross_reference(self, model):
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
    type = 'MPC'

    def __init__(self):
        Constraint.__init__(self)

    def add_card(self, card, comment=''):
        if comment:
            self._comment = comment

        #: Set identification number. (Integer > 0)
        self.conid = integer(card, 1, 'conid')
        #: Identification number of grid or scalar point. (Integer > 0)
        self.gids = []
        #: Component number. (Any one of the Integers 1 through 6 for grid
        #: points; blank or zero for scalar points.)
        self.constraints = []
        #: Coefficient. (Real; Default = 0.0 except A1 must be nonzero.)
        self.enforced = []

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
            self.gids.append(grid)
            self.constraints.append(component)
            self.enforced.append(enforcedi)
            i += 1

            if ifield + 4 > nfields and i != 2:
                # if G2 is empty (it's ifield+4 because nfields is length based and not loop friendly)
                break
            grid = integer(card, ifield + 3, 'G%i' % i)
            component = components_or_blank(card, ifield + 4, 'constraint%i' % i, 0)  # scalar point
            enforcedi = double_or_blank(card, ifield + 5, 'enforced%i' % i)
            self.gids.append(grid)
            self.constraints.append(component)
            self.enforced.append(enforcedi)
            i += 1

    def add_op2_data(self, data, comment=''):
        if comment:
            self._comment = comment
        msg = '%s has not implemented data parsing' % self.type
        raise NotImplementedError(msg)

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        msg = ', which is required by %s=%s' % (self.type, self.conid)
        return self._nodeIDs(nodes=self.gids, allow_empty_nodes=True, msg=msg)

    def cross_reference(self, model):
        msg = ', which is required by %s=%s' % (self.type, self.conid)
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
                    msg += '%16i%16s%16s\n' % (grid, component, print_scientific_double(enforced))
                elif i % 2 == 1:
                    msg += '%-8s%16i%16s%16s\n' % ('*', grid, component, print_scientific_double(enforced))
                else:
                    msg += '%-8s%16s%16i%16s%16s\n' % ('*', '', grid, component, print_scientific_double(enforced))
        else:
            for i, grid, component, enforced in zip(count(), grids, constraints, enforceds):
                if i == 0:
                    msg += '%16i%16s%16s\n' % (grid, component, print_float_16(enforced))
                elif i % 2 == 1:
                    msg += '%-8s%16i%16s%16s\n' % ('*', grid, component, print_float_16(enforced))
                else:
                    msg += '%-8s%16s%16i%16s%16s\n' % ('*', '', grid, component, print_float_16(enforced))
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

    def __init__(self):
        Constraint.__init__(self)

    def add_card(self, card, comment=''):
        if comment:
            self._comment = comment

        self.conid = integer(card, 1, 'sid')
        if card.field(5) in [None, '']:
            self.gids = [integer(card, 2, 'G1'),]
            self.constraints = [components_or_blank(card, 3, 'C1', 0)]
            self.enforced = [double_or_blank(card, 4, 'D1', 0.0)]
        else:
            self.gids = [
                integer(card, 2, 'G1'),
                integer_or_blank(card, 5, 'G2'),
            ]
            # :0 if scalar point 1-6 if grid
            self.constraints = [components_or_blank(card, 3, 'C1', 0),
                                components_or_blank(card, 6, 'C2', 0)]
            self.enforced = [double_or_blank(card, 4, 'D1', 0.0),
                             double_or_blank(card, 7, 'D2', 0.0)]

    def add_op2_data(self, data, comment=''):
        if comment:
            self._comment = comment
        self.conid = data[0]
        self.gids = [data[1]]
        self.constraints = [data[2]]
        self.enforced = [data[3]]

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

    def __init__(self):
        Constraint.__init__(self)

    def add_card(self, card, comment=''):
        if comment:
            self._comment = comment
        self.conid = integer(card, 1, 'sid')
        self.components = components(card, 2, 'components')
        self.entity = string(card, 3, 'entity')
        self.entity_id = integer(card, 4, 'entity_id')

    def add_op2_data(self, data, comment=''):
        if comment:
            self._comment = comment
        raise NotImplementedError('GMSPC')

    def cross_reference(self, model):
        """TODO: xref"""
        msg = ', which is required by %s=%s' % (self.type, self.conid)

    def uncross_reference(self):
        pass

    def raw_fields(self):
        fields = ['GMSPC', self.conid, self.components, self.entity, self.entity_id]
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

    def __init__(self):
        # defines everything :) at least until cross-referencing methods are
        # implemented
        Constraint.__init__(self)

    def add_card(self, card, comment=''):
        if comment:
            self._comment = comment
        #: Identification number of a single-point constraint set.
        self.conid = integer(card, 1, 'conid')
        #: Ring identification number. See RINGAX entry.
        self.rid = integer(card, 2, 'rid')
        #: Harmonic identification number. (Integer >= 0)
        self.hid = integer(card, 3, 'hid')
        #: Component identification number. (Any unique combination of the
        #: Integers 1 through 6.)
        self.c = components(card, 4, 'c')
        #: Enforced displacement value
        self.d = double(card, 5, 'd')

    def add_op2_data(self, data, comment=''):
        if comment:
            self._comment = comment
        msg = '%s has not implemented data parsing' % self.type
        raise NotImplementedError(msg)

    def cross_reference(self, model):
        msg = ', which is required by %s=%s' % (self.type, self.conid)
        #self.rid = model.ring[self.rid]
        #self.hid = model.harmonic[self.hid]
        pass

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

    def __init__(self):
        Constraint.__init__(self)

    def add_card(self, card, comment=''):
        if comment:
            self._comment = comment
        self.conid = integer(card, 1, 'conid')
        self.constraints = components(card, 2, 'constraints')  # 246 = y; dx, dz dir
        nodes = card.fields(3)
        self.nodes = expand_thru(nodes)
        self.nodes.sort()

    def add_op2_data(self, data, comment=''):
        if comment:
            self._comment = comment
        self.conid = data[0]
        self.constraints = data[1]
        self.nodes = data[2:]

    @property
    def node_ids(self):
        msg = ', which is required by SPC1; conid=%s' % self.conid
        return self._nodeIDs(self.nodes, allow_empty_nodes=True, msg=msg)

    def cross_reference(self, model):
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

    +---------+---+---+------+
    | SPCADD  | 2 | 1 |   3  |
    +---------+---+---+------+
    """
    type = 'SPCADD'

    def __init__(self, card=None, data=None, comment=''):
        ConstraintADD.__init__(self)

    def add_card(self, card, comment=''):
        if comment:
            self._comment = comment
        self.conid = integer(card, 1, 'conid')
        sets = card.fields(2)
        self.sets = expand_thru(sets)
        self.sets.sort()

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
        msg = ', which is required by %s=%s' % (self.type, self.conid)
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
    +--------+----+----+-----+
    | MPCADD | 2  |  1 |  3  |
    +--------+----+----+-----+
    """
    type = 'MPCADD'

    def __init__(self):
        ConstraintADD.__init__(self)

    def add_card(self, card, comment=''):
        if comment:
            self._comment = comment
        self.conid = integer(card, 1, 'conid')
        sets = card.fields(2)
        self.sets = expand_thru(sets)
        self.sets.sort()

    def add_op2_data(self, data, comment=''):
        if comment:
            self._comment = comment
        msg = '%s has not implemented data parsing' % self.type
        raise NotImplementedError(msg)

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
        for i, mpc in enumerate(self.sets):
            self.sets[i] = model.MPC(mpc)
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
