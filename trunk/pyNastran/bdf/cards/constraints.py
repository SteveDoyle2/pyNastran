# pylint: disable=R0904,R0902
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from math import ceil
from itertools import izip, count

from pyNastran.bdf.cards.baseCard import BaseCard, expand_thru
from pyNastran.bdf.assign_type import (integer, integer_or_blank,
    double, double_or_blank,
    components, components_or_blank)

class constraintObject(object):
    def __init__(self):
        self.constraints = {}  # SPC, SPC1, SPCD, etc...
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
            #print "already has key...key=%s\n%s" %(key,str(constraint))
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
        for key, add_constraint in sorted(self.add_constraints.iteritems()):
            for i, spcID in enumerate(add_constraint.sets):
                add_constraint.sets[i] = self._spc(spcID)
            self.add_constraints[key] = add_constraint
        #self.add_constraints = add_constraints2

        constraints2 = {}
        for key, constraints in sorted(self.constraints.iteritems()):
            constraints2[key] = []
            for constraint in constraints:
                constraints2[key].append(constraint.cross_reference(model))

        self.constraints = constraints2

    def createConstraintsForID(self):
        """
        This function returns all the constraints with an given constraint ID.
        For example an MPCADD that references 2 MPCADDs which reference 4 MPCs
        should return 4 MPCs (or rather the IDs of those MPCs).

        @todo This function *should* also find unassociated constraints.
         not really done yet, idea needs to be integrated/separated from
         cross-referencing.  no point in doing it twice
        """
        constraints2 = {}
        referencedConstraints = {}
        # some of the ADDConstraint keys are MPCADDs/SPCADDs, some are not
        for key, add_constraint in sorted(self.add_constraints.iteritems()):
            constraints = []
            for i, spcID in enumerate(add_constraint.sets):
                constraintIDs = add_constraint.getConstraintIDs()
                constraints[spcID] = constraintIDs
                constraints += constraintIDs
                #constraints.append(spcID)
                constraints2[key] = [spcID]

            constraints2[key] = constraints

        ## not needed b/c there are no MPCADD/SPCADD
        #for key,constraints in sorted(self.constraints.iteritems()):
            #for constraint in constraints:
                #conID = constraint.ConID()
                #constraints2[conID]
        constraints3 = self.remapSPCs(constraints2)

    # def remapSPCs(self, constraints):
    #     """not really done yet"""
    #     ## takes the MPCADDs that reference MPCADDs and makes them
    #     ## reference MPCs
    #     constraints2 = {}
    #     key = constraints.keys()
    #     nkeys = len(key) - 1
    #     for i in xrange(nkeys):
    #         key = keys[i]
    #         constraints2[key]
    #         for j in xrange(nkeys):
    #             if i > j:
    #                 constraints2[key].append(constraints[key])
    #     return constraints2

    def ConstraintID(self):
        if isinstance(self.conid, int):
            return self.conid
        return self.conid.conid

    def getConstraintIDs(self):
        IDs = []
        for key, constraints in sorted(self.add_constraints.iteritems()):
            conID = constraints.ConID()
            IDs.append(conID)

        for key, constraints in sorted(self.constraints.iteritems()):
            for constraint in constraints:
                conID = constraint.ConID()
                IDs.append(conID)

        IDs = list(set(IDs))
        IDs.sort()
        return IDs

    def __repr__(self):
        msg = ''
        # write the SPCADD/MPCADD cards
        for key, add_constraint in sorted(self.add_constraints.iteritems()):
            msg += str(add_constraint)

        for key, constraints in sorted(self.constraints.iteritems()):
            for constraint in constraints:
                msg += str(constraint)
        return msg


class Constraint(BaseCard):
    def __init__(self, card, data):
        pass

    def rawFields(self):
        fields = [self.type, self.conid]
        return fields


class SUPORT1(Constraint):
    """
    @code
    SUPORT1 SID ID1 C1 ID2 C2 ID3 C3
    @endcode
    """
    type = 'SUPORT1'

    def __init__(self, card=None, data=None, comment=''):
        Constraint.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.conid = integer(card, 1, 'conid')  # really a support id sid

            self.IDs = []
            self.Cs = []
            nFields = len(card)
            nceil = int(ceil((nFields - 1.) / 2.))
            for i in xrange(2, nceil):
                ID = integer(card, 2 * i, 'ID%s' + str(i))
                C = components(card, i + 1, 'component%s' + str(i))
                self.IDs.append(ID)
                self.Cs.append(C)
        else:
            msg = '%s has not implemented data parsing' % self.type
            raise NotImplementedError(msg)

    def rawFields(self):
        fields = ['SUPORT1', self.conid]
        for ID, c in izip(self.IDs, self.Cs):
            fields += [ID, c]
        return fields


class SUPORT(Constraint):
    """
    @code
    SUPORT      ID1 C1 ID2 C2 ID3 C3 ID4 C4
    SUPORT1 SID ID1 C1 ID2 C2 ID3 C3
    @endcode
    """
    type = 'SUPORT'

    def __init__(self, card=None, data=None, comment=''):
        Constraint.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            fields = card.fields(1)
        else:
            fields = data

        self.IDs = []
        self.Cs = []
        for i in xrange(0, len(fields), 2):
            self.IDs.append(fields[i])
            self.Cs.append(fields[i + 1])

    def rawFields(self):
        fields = ['SUPORT']
        for ID, c in izip(self.IDs, self.Cs):
            fields += [ID, c]
        return fields


class MPC(Constraint):
    type = 'MPC'

    def __init__(self, card=None, data=None, comment=''):
        Constraint.__init__(self, card, data)
        if comment:
            self._comment = comment

        if card:
            ## Set identification number. (Integer > 0)
            self.conid = integer(card, 1, 'conid')
            ## Identification number of grid or scalar point. (Integer > 0)
            self.gids = []
            ## Component number. (Any one of the Integers 1 through 6 for grid
            ## points; blank or zero for scalar points.)
            self.constraints = []
            ## Coefficient. (Real; Default = 0.0 except A1 must be nonzero.)
            self.enforced = []

            fields = card.fields(0)
            nFields = len(fields) - 1
            for iField in xrange(2, nFields, 8):
                grid = integer(card, iField, 'gid'),
                component = components_or_blank(card, iField + 1, 'constraint', 0)  # scalar point
                value = double_or_blank(card, iField + 2, 'enforced', 0.0)
                self.gids.append(grid)
                self.constraints.append(component)
                self.enforced.append(value)

                if iField + 3 > nFields:
                    break
                grid = integer(card, iField + 3, 'gid')
                component = components_or_blank(card, iField + 4, 'constraint', 0)  # scalar point
                value = double_or_blank(card, iField + 5, 'enforced')
                self.gids.append(grid)
                self.constraints.append(component)
                self.enforced.append(value)

            # reduce the size if there are duplicate Nones
            #nConstraints = max(len(self.gids       ),
            #                   len(self.constraints),
            #                   len(self.enforced   ))
            #self.gids        = self.gids[       0:nConstraints]
            #self.constraints = self.constraints[0:nConstraints]
            #self.enforced    = self.enforced[   0:nConstraints]
        else:
            msg = '%s has not implemented data parsing' % self.type
            raise NotImplementedError(msg)

    def rawFields(self):  # MPC
        fields = ['MPC', self.conid]
        for (i, gid, constraint, enforced) in izip(count(), self.gids,
             self.constraints, self.enforced):
            #print [gid,constraint,enforced]
            fields += [gid, constraint, enforced]
            if i % 2 == 1 and i > 0:
                fields.append(None)
                fields.append(None)
        return fields


class SPC(Constraint):
    """
    Defines enforced displacement/temperature (static analysis)
    velocity/acceleration (dynamic analysis)

    @code
    SPC SID G1 C1 D1   G2 C2 D2
    SPC 2   32 3  -2.6  5
    @endcode
    """
    type = 'SPC'

    def __init__(self, card=None, data=None, comment=''):
        Constraint.__init__(self, card, data)

        if comment:
            self._comment = comment
        if card:
            self.conid = integer(card, 1, 'sid')
            self.gids = [integer(card, 2, 'G1'), integer_or_blank(card, 5, 'G2')]
            ## 0 if scalar point 1-6 if grid
            
            self.constraints = [components_or_blank(card, 3, 'C1', 0),
                                components_or_blank(card, 6, 'C2', 0)]
            self.enforced = [double_or_blank(card, 3, 'D1', 0.0),
                             double_or_blank(card, 7, 'D2', 0.0)]

            # reduce the size if there are duplicate Nones
            nConstraints = max(len(self.gids),
                               len(self.constraints),
                               len(self.enforced))
            self.gids = self.gids[0:nConstraints]
            self.constraints = self.constraints[0:nConstraints]
            self.enforced = self.enforced[0:nConstraints]
        else:
            self.conid = data[0]
            self.gids = [data[1]]
            self.constraints = [data[2]]
            self.enforced = [data[3]]

    def getNodeDOFs(self, model):
        pass
        #return conid,dofs

    def cross_reference(self, i, node):
        dofCount = 0
        self.gids[i] = node
        #for (i,constraint) in enumerate(self.constraints):
        #    if self.constraint is None:
        #        node = self.Node(self.gids[i])
        #        if not node.Is('GRID'): # SPOINT, EPOINT, DPOINT
        #            dofCount+=1
        #        else:
        #            dofCount+=6
        return dofCount

    def rawFields(self):
        fields = ['SPC', self.conid]
        for (gid, constraint, enforced) in izip(self.gids, self.constraints,
                self.enforced):
            fields += [gid, constraint, enforced]
        return fields


class SPCD(SPC):
    """
    Defines an enforced displacement value for static analysis and an enforced
    motion value (displacement, velocity or acceleration) in dynamic analysis.

    @code
    SPCD SID G1  C1   D1 G2 C2  D2
    SPCD 100 32 436 -2.6  5    2.9
    @endcode
    """
    type = 'SPCD'

    def __init__(self, card=None, data=None, comment=''):
        # defines everything :) at least until cross-referencing methods are
        # implemented
        SPC.__init__(self, card, data)
        if comment:
            self._comment = comment

    def rawFields(self):
        fields = ['SPCD', self.conid]
        for (gid, constraint, enforced) in izip(self.gids, self.constraints,
                self.enforced):
            fields += [gid, constraint, enforced]
        return fields


class SPCAX(Constraint):
    """
    Defines a set of single-point constraints or enforced displacements
    for conical shell coordinates.

    @code
    SPCAX SID RID HID C    D
    SPCAX 2   3     4 13 6.0
    @endcode
    """
    type = 'SPCAX'

    def __init__(self, card=None, data=None, comment=''):
        # defines everything :) at least until cross-referencing methods are
        # implemented
        SPC.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            ## Identification number of a single-point constraint set.
            self.conid = integer(card, 1, 'conid')
            ## Ring identification number. See RINGAX entry.
            self.rid = integer(card, 2, 'rid')
            ## Harmonic identification number. (Integer >= 0)
            self.hid = integer(card, 3, 'hid')
            ## Component identification number. (Any unique combination of the
            ## Integers 1 through 6.)
            self.c = components(card, 4, 'c')
            ## Enforced displacement value
            self.d = double(card, 5, 'd')
        else:
            msg = '%s has not implemented data parsing' % self.type
            raise NotImplementedError(msg)

    def cross_reference(self, i, node):
        pass

    def rawFields(self):
        fields = ['SPCAX', self.conid, self.rid, self.hid, self.c, self.d]
        return fields


class SPC1(Constraint):
    """
    @code
    SPC1 SID C G1 G2 G3 G4 G5 G6
    G7 G8 G9 -etc.-

    SPC1     3       246     209075  209096  209512  209513  209516
    SPC1 3 2 1 3 10 9 6 5
    2 8

    SPC1 SID C    G1 THRU G2
    SPC1 313 12456 6 THRU 32
    @endcode
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
            msg = '%s has not implemented data parsing' % self.type
            raise NotImplementedError(msg)

    def cross_reference(self, i, node):
        self.nodes[i] = node

    def rawFields(self):  # SPC1
        #test = [i for i in xrange(self.nodes[0],self.nodes[-1]+1)]
        #print "self.nodes = ",self.nodes
        #print "test       = ",test
        #if self.nodes==test:
        #    nodes = [self.nodes[0],'THRU',self.nodes[-1]]
        #else:
        #print "SPC1 self.nodes = ",self.nodes
        nodes = [int(i) for i in self.nodes]  # SPC1
        fields = ['SPC1', self.conid, self.constraints] + nodes
        return fields

#class ADDConstraint(Constraint):
#    def __init__(self,card,data):
#        self.__init__(Constraint)


class ConstraintADD(Constraint):
    def __init__(self, card, data):
        Constraint.__init__(self, card, data)

    def _reprSpcMpcAdd(self, fields):
        outSPCs = ''
        fieldSets = []

        #print "repr---mpcadd"
        #print self.sets
        for spcsets in sorted(self.sets):
            #print "*spcsets",spcsets
            if isinstance(spcsets, int):  # spcset wasnt found
                #print "int unfound...%s" %(spcsets)
                #outSPCs += str(spcsets)
                fields.append(spcsets)
            elif isinstance(spcsets, list):
                #print 'list'
                for spcset in sorted(spcsets):
                    fieldSets.append(spcset.conid)
                    outSPCs += str(spcset)
            else:
                #print 'dict'
                #outSPCs += str(spcsets.conid)
                for (key, spcset) in sorted(spcsets.iteritems()):
                    fieldSets.append(spcsets.conid)

        # SPCADD
        return self.print_card(fields + list(set(fieldSets))) + outSPCs


class SPCADD(ConstraintADD):
    """
    Defines a single-point constraint set as a union of single-point constraint
    sets defined on SPC or SPC1 entries.

    @code
    SPCADD   2       1       3
    @endcode
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
        positionSPCs = []
        typesFound = ['SPCADD']
        (scaleFactors, loads) = self.getReducedConstraints()
        return (typesFound, positionSPCs)

    def cross_reference(self, i, node):
        #dofCount = 0
        self.sets.sort()
        self.sets[i] = node

    def rawFields(self):
        fields = ['SPCADD', self.conid]  # +self.sets
        for setID in self.sets:
            #print "setID = ",setID
            fields.append(setID)
        return fields
        #return self._reprSpcMpcAdd(fields)


class MPCADD(ConstraintADD):
    r"""
    Defines a multipoint constraint equation of the form
    \f$ \Sigma_j A_j u_j =0 \f$ where \f$ u_j \f$ represents
    degree-of-freedom \f$ C_j \f$ at grid or scalar point \f$ G_j \f$.
    mPCADD   2       1       3
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

    def cross_reference(self, i, node):
        #dofCount = 0
        self.sets.sort()
        self.sets[i] = node

    def rawFields(self):
        #outSPCs = ''
        fields = ['MPCADD', self.conid]  # +self.sets
        for setID in self.sets:
            fields.append(setID)
        return fields
        #return self._reprSpcMpcAdd(fields)