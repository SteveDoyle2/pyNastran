# pylint: disable=R0904,R0902
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from itertools import izip, count

from pyNastran.bdf.cards.baseCard import BaseCard, expand_thru


class constraintObject2(object):
    def __init__(self):
        self.constraints = {}  # SPC, SPC1, SPCD, etc...
        self.addConstraints = {}  # SPCADD

    def Add(self, ADD_Constraint):
        """bad name, but adds an SPCADD or MPCADD"""
        key = ADD_Constraint.conid
        if key in self.addConstraints:
            print("already has key...key=%s\n%s" % (key, str(ADD_Constraint)))
        else:
            self.addConstraints[key] = ADD_Constraint

    def append(self, constraint):
        key = constraint.conid
        if key in self.constraints:
            #print "already has key...key=%s\n%s" %(key,str(constraint))
            self.constraints[key].append(constraint)
        else:
            self.constraints[key] = [constraint]
        ###

    def _spc(self, spcID):
        """could be an spcadd/mpcadd or spc/mpc,spc1,spcd,spcax,suport1"""
        if spcID in self.addConstraints:
            return self.addConstraints[spcID]
        return self.constraints[spcID]

    def cross_reference(self, model):
        #return
        #addConstraints2 = {}
        for key, addConstraint in sorted(self.addConstraints.iteritems()):
            for i, spcID in enumerate(addConstraint.sets):
                addConstraint.sets[i] = self._spc(spcID)
            self.addConstraints[key] = addConstraint
        #self.addConstraints = addConstraints2

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
        for key, addConstraint in sorted(self.addConstraints.iteritems()):
            constraints = []
            for i, spcID in enumerate(addConstraint.sets):
                constraintIDs = addConstraint.getConstraintIDs()
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

    def remapSPCs(self, constraints):
        """not really done yet"""
        ## takes the MPCADDs that reference MPCADDs and makes them
        ## reference MPCs

        constraints2 = {}
        Keys = constraints.keys()
        nKeys = len(Keys) - 1
        for i in xrange(nKeys):
            Key = Keys[i]
            constraints2[Key]
            for j in xrange(nKeys):
                if i > j:
                    constraints2[Key].append(constraints[Key])

        return constraints2

    def ConstraintID(self):
        if isinstance(self.conid, int):
            return self.conid
        return self.conid.conid

    def getConstraintIDs(self):
        IDs = []
        for key, constraints in sorted(self.addConstraints.iteritems()):
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
        for key, addConstraint in sorted(self.addConstraints.iteritems()):
            msg += str(addConstraint)

        for key, constraints in sorted(self.constraints.iteritems()):
            for constraint in constraints:
                msg += str(constraint)
        return msg


class constraintObject(object):

    def __init__(self):
        self.constraints = {}  # SPC, SPC1, SPCD, etc...
        self.addConstraints = {}  # SPCADD
        self.resolvedConstraints = []

    def add(self, constraint):
        conid = constraint.conid
        assert conid not in self.addConstraints
        self.addConstraints[conid] = constraint

    def append(self, constraint):
        key = constraint.conid
        if key in self.constraints:
            self.constraints[key].append(constraint)
        else:
            self.constraints[key] = [constraint]

    def getConstraintIDs(self):
        IDs = []
        for key, constraints in sorted(self.constraints.iteritems()):
            for constraint in constraints:
                conID = constraint.ConID()
                IDs.append(conID)
        return IDs

    def ConstraintID(self):
        if isinstance(self.conid, int):
            return self.conid
        return self.conid.conid

    def crossReference(self, model):
        #print "xref spcadds..."
        #print "spcadds = ",self.addConstraints
        if self.addConstraints:
            for (key, addConstraint) in sorted(self.addConstraints.iteritems()):  # SPCADDs
                self.crossReference_AddConstraint(key, addConstraint)
            #print "spcadds2 = ",self.addConstraints
        else:
            pass  # not done, no spcsets

        # xrefs nodes...not done...
        #print "xref spc/spc1/spcd..."
        return
        for key, constraints in sorted(self.constraints.iteritems()):  # SPC, SPC1, SPCD
            for constraint in constraints:
                #constraint.crossR
                pass
                #if constraint.type=='SPCADD'
                #else:
                #(conid,nodeDOFs) = constraint.getNodeDOFs(model) # the constrained nodes
                #for nodeDOF in nodeDOFs:
                #    self.addConstraint(conid,nodeDOF)

    def crossReference_AddConstraint(self, key, addConstraint):
        """
        cross references MPCSETs and SPCSETs
        """
        #print "add key=%s" %(key)
        #sets = type(addConstraint)
        spcsets = addConstraint.sets
        #sys.stdout.flush()
        #print str(addConstraint.sets)
        #sys.stdout.flush()
        #print "spcsets = ",spcsets
        for (i, conid) in enumerate(spcsets):
            #print "conid = ",conid
            #conid = spcset.conid
            #print "self.addConstraints[conid] = ",self.getConstraint(conid)
            constraint = self.getConstraint(conid)
            #print 'newSlot = ',self.addConstraints[key].gids[i]
            self.addConstraints[key].cross_reference(i, constraint)
            #self.addConstraints[key].gids[i] = self.getConstraint(conid)
            #print "spcadds* = ",self.addConstraints

    #def popConstraint(self, conid):  # need to not throw away constraints...
    #    if conid in self.addConstraints:
    #        return self.addConstraints[conid]
    #    elif conid in self.constraints:
    #        return self.constraints[conid]
    #    else:
    #        return conid
    #    ###

    def getConstraint(self, conid):
        #print "sid=%s" %(conid)
        if conid in self.addConstraints:
            return self.addConstraints[conid]
        elif conid in self.constraints:
            return self.constraints[conid]
        else:
            return conid

    def addConstraint(self, conid, nodeDOF):
        (nid, dofs) = nodeDOF
        for dof in dofs:
            self.resolvedConstraints.append((nid, dof))
        ###

    def __repr__(self):
        msg = ''
        for addID, spcadd in sorted(self.addConstraints.iteritems()):
            msg += str(spcadd)  # this writes the SPC cards as well
        return msg

        msg = ''
        #print "repr %s" %(self.addConstraints)
        if self.addConstraints:
            for addID, spcadd in sorted(self.addConstraints.iteritems()):
                msg += str(spcadd)  # this writes the SPC cards as well
        else:
            for key, constraintSets in sorted(self.constraints.iteritems()):
                for constraint in constraintSets:
                    msg += str(constraint)

        #print msg
        return msg

        # works for spc, spc1, spcd
        #for key,constraintSets in sorted(self.constraints.iteritems()):
        #    for constraint in constraintSets:
        #        msg += str(constraint)
        return msg


class Constraint(BaseCard):
    def __init__(self, card, data):
        pass

    def rawFields(self):
        fields = [self.type, self.conid]
        return fields


class SUPORT1(Constraint):
    """
    SUPORT1 SID ID1 C1 ID2 C2 ID3 C3
    """
    type = 'SUPORT1'

    def __init__(self, card=None, data=None):
        Constraint.__init__(self, card, data)
        self.conid = card.field(1)  # really a support id sid
        fields = card.fields(2)

        self.IDs = []
        self.Cs = []
        for i in xrange(0, len(fields), 2):
            self.IDs.append(fields[i])
            self.Cs.append(fields[i + 1])

    def rawFields(self):
        fields = ['SUPORT1', self.conid]
        for ID, c in izip(self.IDs, self.Cs):
            fields += [ID, c]
        return fields


class SUPORT(Constraint):
    """
    SUPORT      ID1 C1 ID2 C2 ID3 C3 ID4 C4
    SUPORT1 SID ID1 C1 ID2 C2 ID3 C3
    """
    type = 'SUPORT'

    def __init__(self, card=None, data=None):
        Constraint.__init__(self, card, data)
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

    def __init__(self, card=None, data=None):
        Constraint.__init__(self, card, data)
        self.conid = card.field(1)
        #self.gids        = [card.field(2),card.field(5,None)]
        #self.constraints = [card.field(3),card.field(6,None)] # 0 if scalar point 1-6 if grid
        #self.enforced    = [card.field(4),card.field(7,None)]

        self.gids = []
        self.constraints = []  # 0 if scalar point 1-6 if grid
        self.enforced = []
        #print "-----------"

        fields = card.fields(0)
        nFields = len(fields) - 1
        #print "fields = ",fields
        #print "nFields = ",nFields
        for iField in xrange(2, nFields, 8):
            pack1 = [card.field(iField),
                     card.field(iField + 1, 0),
                     card.field(iField + 2, 0.)]
            #print "pack1 = ",pack1
            self.gids.append(pack1[0])
            self.constraints.append(pack1[1])  # default=0 scalar point
            self.enforced.append(pack1[2])  # default=0.0

            pack2 = [card.field(iField + 3),
                     card.field(iField + 4, 0),
                     card.field(iField + 5, 0.)]
            if pack2 != [None, 0, 0.]:
                #print "pack2 = ",pack2
                #print "adding pack2"
                self.gids.append(pack2[0])
                self.constraints.append(pack2[1])  # default=0 scalar point
                self.enforced.append(pack2[2])  # default=0.0

        # reduce the size if there are duplicate Nones
        #nConstraints = max(len(self.gids       ),
        #                   len(self.constraints),
        #                   len(self.enforced   ))
        #self.gids        = self.gids[       0:nConstraints]
        #self.constraints = self.constraints[0:nConstraints]
        #self.enforced    = self.enforced[   0:nConstraints]

    def rawFields(self):  # MPC
        fields = ['MPC', self.conid]
        for (i, gid, constraint, enforced) in izip(count(), self.gids, self.constraints, self.enforced):
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
    SPC SID G1 C1 D1   G2 C2 D2
    SPC 2   32 3  -2.6  5
    """
    type = 'SPC'

    def __init__(self, card=None, data=None):
        Constraint.__init__(self, card, data)

        if card:
            self.conid = card.field(1)
            self.gids = [card.field(2), card.field(5, None)]
            ## 0 if scalar point 1-6 if grid
            self.constraints = [card.field(3), card.field(6, None)]
            self.enforced = [card.field(4), card.field(7, None)]

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
        for (gid, constraint, enforced) in izip(self.gids, self.constraints, self.enforced):
            fields += [gid, constraint, enforced]
        return fields


class SPCD(SPC):
    """
    Defines an enforced displacement value for static analysis and an enforced
    motion value (displacement, velocity or acceleration) in dynamic analysis.
    SPCD SID G1  C1   D1 G2 C2  D2
    SPCD 100 32 436 -2.6  5    2.9
    """
    type = 'SPCD'

    def __init__(self, card=None, data=None):
        # defines everything :) at least until cross-referencing methods are
        # implemented
        SPC.__init__(self, card, data)

    def rawFields(self):
        fields = ['SPCD', self.conid]
        for (gid, constraint, enforced) in izip(self.gids, self.constraints, self.enforced):
            fields += [gid, constraint, enforced]
        return fields


class SPCAX(Constraint):
    """
    Defines a set of single-point constraints or enforced displacements
    for conical shell coordinates.
    SPCAX SID RID HID C    D
    SPCAX 2   3     4 13 6.0
    """
    type = 'SPCAX'

    def __init__(self, card=None, data=None):
        # defines everything :) at least until cross-referencing methods are
        # implemented
        SPC.__init__(self, card, data)

        ## Identification number of a single-point constraint set.
        self.conid = card.field(1)
        ## Ring identification number. See RINGAX entry.
        self.rid = card.field(2)
        ## Harmonic identification number. (Integer >= 0)
        self.hid = card.field(3)
        ## Component identification number. (Any unique combination of the
        ## Integers 1 through 6.)
        self.c = card.field(4)
        ## Enforced displacement value
        self.d = card.field(5)

    def cross_reference(self, i, node):
        pass

    def rawFields(self):
        fields = ['SPCAX', self.conid, self.rid, self.hid, self.c, self.d]
        return fields


class SPC1(Constraint):
    """
    SPC1 SID C G1 G2 G3 G4 G5 G6
    G7 G8 G9 -etc.-

    SPC1     3       246     209075  209096  209512  209513  209516
    SPC1 3 2 1 3 10 9 6 5
    2 8

    SPC1 SID C    G1 THRU G2
    SPC1 313 12456 6 THRU 32
    """
    type = 'SPC1'

    def __init__(self, card=None, data=None):
        Constraint.__init__(self, card, data)
        self.conid = card.field(1)
        self.constraints = str(card.field(2, ''))  # 246 = y; dx, dz dir
        nodes = card.fields(3)
        self.nodes = expand_thru(nodes)
        self.nodes.sort()

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
        return self.printCard(fields + list(set(fieldSets))) + outSPCs  # SPCADD


class SPCADD(ConstraintADD):
    """
    Defines a single-point constraint set as a union of single-point constraint
    sets defined on SPC or SPC1 entries.
    SPCADD   2       1       3
    """
    type = 'SPCADD'

    def __init__(self, card=None, data=None):
        ConstraintADD.__init__(self, card, data)
        self.conid = card.field(1)
        sets = card.fields(2)
        self.sets = expand_thru(sets)
        self.sets.sort()

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

    def __init__(self, card=None, data=None):
        ConstraintADD.__init__(self, card, data)
        self.conid = card.field(1)
        sets = card.fields(2)
        self.sets = expand_thru(sets)
        self.sets.sort()

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
