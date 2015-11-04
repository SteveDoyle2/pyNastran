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
import warnings

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
            for i, spcID in enumerate(add_constraint.sets):
                add_constraint.sets[i] = self._spc(spcID)
            self.add_constraints[key] = add_constraint
        #self.add_constraints = add_constraints2

        constraints2 = {}
        for key, constraints in sorted(iteritems(self.constraints)):
            constraints2[key] = []
            for constraint in constraints:
                constraints2[key].append(constraint.cross_reference(model))

        self.constraints = constraints2

    def createConstraintsForID(self):
        """
        This function returns all the constraints with an given constraint ID.
        For example an MPCADD that references 2 MPCADDs which reference 4 MPCs
        should return 4 MPCs (or rather the IDs of those MPCs).

        .. todo:: This function *should* also find unassociated constraints.
         not really done yet, idea needs to be integrated/separated from
         cross-referencing.  no point in doing it twice
        """
        raise NotImplementedError()
        constraints2 = {}
        referencedConstraints = {}
        # some of the ADDConstraint keys are MPCADDs/SPCADDs, some are not
        for key, add_constraint in sorted(iteritems(self.add_constraints)):
            constraints = []
            for i, spc_id in enumerate(add_constraint.sets):
                constraintIDs = add_constraint.getConstraintIDs()
                constraints[spc_id] = constraintIDs
                constraints += constraintIDs
                #constraints.append(spc_id)
                constraints2[key] = [spc_id]

            constraints2[key] = constraints

        # not needed b/c there are no MPCADD/SPCADD
        #for key, constraints in sorted(iteritems(self.constraints)):
            #for constraint in constraints:
                #conID = constraint.ConID()
                #constraints2[conID]
        constraints3 = self.remapSPCs(constraints2)

    # def remapSPCs(self, constraints):
         #"""not really done yet"""
         ## takes the MPCADDs that reference MPCADDs and makes them
         ## reference MPCs
         #constraints2 = {}
         #key = constraints.keys()
         #nkeys = len(key) - 1
         #for i in range(nkeys):
             #key = keys[i]
             #constraints2[key]
             #for j in range(nkeys):
                 #if i > j:
                     #constraints2[key].append(constraints[key])
         #return constraints2

    def ConstraintID(self):
        if isinstance(self.conid, integer_types):
            return self.conid
        return self.conid.conid

    def getConstraintIDs(self):
        IDs = []
        for key, constraints in sorted(iteritems(self.add_constraints)):
            conID = constraints.ConID()
            IDs.append(conID)

        for key, constraints in sorted(iteritems(self.constraints)):
            for constraint in constraints:
                conID = constraint.ConID()
                IDs.append(conID)

        IDs = list(set(IDs))
        IDs.sort()
        return IDs

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
    def __init__(self, card, data):
        pass

    def raw_fields(self):
        fields = [self.type, self.conid]
        return fields

    def _nodeIDs(self, nodes=None, allowEmptyNodes=False, msg=''):
        """returns nodeIDs for repr functions"""
        return _node_ids(self, nodes, allowEmptyNodes, msg)


class SUPORT1(Constraint):
    """
    +---------+-----+-----+----+-----+----+-----+----+
    | SUPORT1 | SID | ID1 | C1 | ID2 | C2 | ID3 | C3 |
    +---------+-----+-----+----+-----+----+-----+----+
    | SUPORT1 |  1  |  2  | 23 |  4  | 15 |  5  |  0 |
    +---------+-----+-----+----+-----+----+-----+----+
    """
    type = 'SUPORT1'

    def __init__(self, card=None, data=None, comment=''):
        Constraint.__init__(self, card, data)

        self.IDs = []
        self.Cs = []
        if comment:
            self._comment = comment
        if card:
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
        else:
            msg = '%s has not implemented data parsing' % self.type
            raise NotImplementedError(msg)
        assert len(self.IDs) > 0
        assert len(self.IDs) == len(self.Cs)

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
        return self._nodeIDs(nodes=self.IDs, allowEmptyNodes=True, msg=msg)

    def cross_reference(self, model):
        msg = ', which is required by SUPORT1'
        self.IDs = model.Nodes(self.IDs, allowEmptyNodes=True, msg=msg)

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

    def __init__(self, card=None, data=None, comment=''):
        Constraint.__init__(self, card, data)
        if comment:
            self._comment = comment

        self.IDs = [] ## TODO:  IDs reference nodes???
        self.Cs = []
        if card:
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
        else:
            fields = data
            for i in range(0, len(fields), 2):
                self.IDs.append(fields[i])
                self.Cs.append(fields[i + 1])
        assert len(self.IDs) > 0
        assert len(self.IDs) == len(self.Cs)

    @property
    def node_ids(self):
        msg = ', which is required by %s' % self.type
        return self._nodeIDs(nodes=self.IDs, allowEmptyNodes=True, msg=msg)

    def cross_reference(self, model):
        msg = ', which is required by %s' % self.type
        self.IDs = model.Nodes(self.IDs, allowEmptyNodes=True, msg=msg)

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

    def __init__(self, card=None, data=None, comment=''):
        SUPORT.__init__(self, card, data)


class MPC(Constraint):
    type = 'MPC'

    def __init__(self, card=None, data=None, comment=''):
        Constraint.__init__(self, card, data)
        if comment:
            self._comment = comment

        if card:
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
        else:
            msg = '%s has not implemented data parsing' % self.type
            raise NotImplementedError(msg)

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        msg = ', which is required by %s=%s' % (self.type, self.conid)
        return self._nodeIDs(nodes=self.gids, allowEmptyNodes=True, msg=msg)

    def cross_reference(self, model):
        msg = ', which is required by %s=%s' % (self.type, self.conid)
        self.gids = model.Nodes(self.gids, allowEmptyNodes=True, msg=msg)

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

    def __init__(self, card=None, data=None, comment=''):
        Constraint.__init__(self, card, data)

        if comment:
            self._comment = comment
        if card:
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
        else:
            self.conid = data[0]
            self.gids = [data[1]]
            self.constraints = [data[2]]
            self.enforced = [data[3]]

    #def getNodeDOFs(self, model):
        #pass

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        msg = ', which is required by %s=%s' % (self.type, self.conid)
        return self._nodeIDs(nodes=self.gids, allowEmptyNodes=True, msg=msg)

    def cross_reference(self, model):
        msg = ', which is required by %s=%s' % (self.type, self.conid)
        self.gids = model.Nodes(self.gids, allowEmptyNodes=True, msg=msg)

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

    def __init__(self, card=None, data=None, comment=''):
        Constraint.__init__(self, card, data)

        if comment:
            self._comment = comment
        if card:
            self.conid = integer(card, 1, 'sid')
            self.components = components(card, 2, 'components')
            self.entity = string(card, 3, 'entity')
            self.entity_id = integer(card, 4, 'entity_id')
        else:
            raise NotImplementedError('GMSPC')

    def cross_reference(self, model):
        msg = ', which is required by %s=%s' % (self.type, self.conid)
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

    def __init__(self, card=None, data=None, comment=''):
        # defines everything :) at least until cross-referencing methods are
        # implemented
        Constraint.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
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
        else:
            msg = '%s has not implemented data parsing' % self.type
            raise NotImplementedError(msg)

    def cross_reference(self, model):
        msg = ', which is required by %s=%s' % (self.type, self.conid)
        #self.rid = model.ring[self.rid]
        #self.hid = model.harmonic[self.hid]
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

    def __init__(self, card=None, data=None, comment=''):
        Constraint.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.conid = integer(card, 1, 'conid')
            self.constraints = components(card, 2, 'constraints')  # 246 = y; dx, dz dir
            nodes = card.fields(3)
            self.nodes = expand_thru(nodes)
            self.nodes.sort()
        else:
            self.conid = data[0]
            self.constraints = data[1]
            self.nodes = data[2:]

    @property
    def node_ids(self):
        msg = ', which is required by SPC1; conid=%s' % self.conid
        return self._nodeIDs(self.nodes, allowEmptyNodes=True, msg=msg)

    def cross_reference(self, model):
        msg = ', which is required by SPC1; conid=%s' % self.conid
        self.nodes = model.Nodes(self.node_ids, allowEmptyNodes=True, msg=msg)

    def raw_fields(self):
        fields = ['SPC1', self.conid, self.constraints] + self.node_ids
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class ConstraintADD(Constraint):
    def __init__(self, card, data):
        Constraint.__init__(self, card, data)


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
        ConstraintADD.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.conid = integer(card, 1, 'conid')
            sets = card.fields(2)
            self.sets = expand_thru(sets)
            self.sets.sort()
        else:
            msg = '%s has not implemented data parsing' % self.type
            raise NotImplementedError(msg)

    def organizeConstraints(self, model):
        """
        Figures out magnitudes of the loads to be applied to the various nodes.
        This includes figuring out scale factors.
        """
        position_spcs = []
        types_found = ['SPCADD']
        (scale_factors, loads) = self.getReducedConstraints()
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

    def __init__(self, card=None, data=None, comment=''):
        ConstraintADD.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.conid = integer(card, 1, 'conid')
            sets = card.fields(2)
            self.sets = expand_thru(sets)
            self.sets.sort()
        else:
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
