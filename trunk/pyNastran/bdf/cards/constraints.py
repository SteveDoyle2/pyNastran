# my code
import sys
from baseCard import BaseCard


class constraintObject(object):
    """
    todo rename constraint id from cid b/c thats the coordinate system id
    """
    def __init__(self):
        self.constraints    = {} # SPC, SPC1, SPCD, etc...
        self.addConstraints = {} # SPCADD
        self.resolvedConstraints = []

    def add(self,constraint):
        cid = constraint.cid
        assert cid not in self.addConstraints
        self.addConstraints[cid] = constraint

    def append(self,constraint):
        key = constraint.cid
        if self.constraints.has_key(key):
            self.constraints[key].append(constraint)
        else:
            self.constraints[key] = [constraint]

    def crossReference(self,mesh):
        #print "xref spcadds..."
        #print "spcadds = ",self.addConstraints
        if self.addConstraints:
            for (key,addConstraint) in sorted(self.addConstraints.items()):  # SPCADDs
                self.crossReference_AddConstraint(key,addConstraint)
            #print "spcadds2 = ",self.addConstraints
        else:
            pass  # not done, no spcsets
        ###

        # xrefs nodes...not done...
        #print "xref spc/spc1/spcd..."
        for key,constraints in sorted(self.constraints.items()): # SPC, SPC1, SPCD
            for constraint in constraints:
                #constraint.crossR
                pass
                #if constraint.type=='SPCADD'
                #else:
                #(cid,nodeDOFs) = constraint.getNodeDOFs(mesh) # the constrained nodes
                #for nodeDOF in nodeDOFs:
                #    self.addConstraint(cid,nodeDOF)
                ###
            ###
        ###
    ###

    def crossReference_AddConstraint(self,key,addConstraint):
        """
        cross references MPCSETs and SPCSETs
        """
        #print "add key=%s" %(key)
        #sets = type(addConstraint)
        spcsets = addConstraint.sets
        #sys.stdout.flush()
        #print str(addConstraint.sets)
        #sys.stdout.flush()
        #sys.exit('xxx1---constraints.py')
        #print "spcsets = ",spcsets
        for i,cid in enumerate(spcsets):
            #print "cid = ",cid
            #cid = spcset.cid
            #print "self.addConstraints[cid] = ",self.getConstraint(cid)
            constraint = self.getConstraint(cid)
            #print 'newSlot = ',self.addConstraints[key].gids[i]
            self.addConstraints[key].crossReference(i,constraint)
            #self.addConstraints[key].gids[i] = self.getConstraint(cid)
            #print "spcadds* = ",self.addConstraints
        ###
    ###

    def popConstraint(self,cid):  # need to not throw away constraints...
        if cid in self.addConstraints:
            return self.addConstraints[cid]
        elif cid in self.constraints:
            return self.constraints[cid]
        else:
            return cid
        ###

    def getConstraint(self,cid):
        if cid in self.addConstraints:
            return self.addConstraints[cid]
        elif cid in self.constraints:
            return self.constraints[cid]
        else:
            return cid
        ###
    def addConstraint(self,cid,nodeDOF):
        (nid,dofs) = nodeDOF
        for dof in dofs:
            self.resolvedConstraints.append( (nid,dof) )
        ###

    def __repr__(self):
        msg = ''
        #print "repr %s" %(self.addConstraints)
        if self.addConstraints:
            for addID,spcadd in sorted(self.addConstraints.items()):
                msg += str(spcadd)  # deceptively this writes the SPC cards as well
        else:
            for key,constraintSets in sorted(self.constraints.items()):
                for constraint in constraintSets:
                    msg += str(constraint)
                ###
            ###
        ###
        #print msg
        #sys.exit('asd')
        return msg
        
        # works for spc, spc1, spcd
        #for key,constraintSets in sorted(self.constraints.items()):
        #    for constraint in constraintSets:
        #        msg += str(constraint)
        return msg

class Constraint(BaseCard):
    def __init__(self,card):
        self.cid  = card.field(1)
    
    def __repr__(self):
        fields = [self.type,self.cid]
        return self.printCard(fields)

class SUPORT1(Constraint):
    """
    SUPORT1 SID ID1 C1 ID2 C2 ID3 C3
    """
    type = 'SUPORT1'
    def __init__(self,card):
        Constraint.__init__(self,card)
        self.cid = card.field(1)  # really a support id sid
        fields   = card.fields(2)
        
        self.IDs = []
        self.Cs  = []
        #print "fields = ",fields
        for i in range(0,len(fields),2):
            #print "i = ",i
            self.IDs.append(fields[i  ])
            self.Cs.append( fields[i+1])
            if fields[i+1]==None:
                break
            ###
        ###

    def __repr__(self):
        fields = ['SUPORT1',self.cid]
        for ID,c in zip(self.IDs,self.Cs):
            fields += [ID,c]
        return self.printCard(fields)

class MPC(Constraint):
    type = 'MPC'
    def __init__(self,card):
        Constraint.__init__(self,card)    # defines self.cid
        self.gids        = [card.field(2),card.field(5,None)]
        self.constraints = [card.field(3),card.field(6,None)] # 0 if scalar point 1-6 if grid
        self.enforced    = [card.field(4),card.field(7,None)]

        # reduce the size if there are duplicate Nones
        nConstraints = max(len(self.gids       ),
                           len(self.constraints),
                           len(self.enforced   ))
        self.gids        = self.gids[       0:nConstraints]
        self.constraints = self.constraints[0:nConstraints]
        self.enforced    = self.enforced[   0:nConstraints]

    def __repr__(self): # MPC
        fields = ['MPC',self.cid]
        for (gid,constraint,enforced) in zip(self.gids,self.constraints,self.enforced):
            fields += [gid,constraint,enforced]
        return self.printCard(fields)
    
class SPC(Constraint):
    """
    Defines enforced displacement/temperature (static analysis)
    velocity/acceleration (dynamic analysis)
    SPC SID G1 C1 D1   G2 C2 D2
    SPC 2   32 3  -2.6  5
    """
    type = 'SPC'
    def __init__(self,card):
        Constraint.__init__(self,card)    # defines self.cid
        self.gids        = [card.field(2),card.field(5,None)]
        self.constraints = [card.field(3),card.field(6,None)] # 0 if scalar point 1-6 if grid
        self.enforced    = [card.field(4),card.field(7,None)]
        
        # reduce the size if there are duplicate Nones
        nConstraints = max(len(self.gids       ),
                           len(self.constraints),
                           len(self.enforced   ))
        self.gids        = self.gids[       0:nConstraints]
        self.constraints = self.constraints[0:nConstraints]
        self.enforced    = self.enforced[   0:nConstraints]

    def getNodeDOFs(self,mesh):
        pass
        return cid,dofs

    def crossReference(self,i,node):
        dofCount = 0
        self.gids[i] = node
        #for (i,constraint) in enumerate(self.constraints):
        #    if self.constraint is None:
        #        node = self.Node(self.gids[i])
        #        if not node.Is('GRID'): # SPOINT, EPOINT, DPOINT
        #            dofCount+=1
        #        else:
        #            dofCount+=6
        #        ###
        #    ###
        return dofCount

    def __repr__(self): # SPC
        fields = ['SPC',self.cid]
        for (gid,constraint,enforced) in zip(self.gids,self.constraints,self.enforced):
            fields += [gid,constraint,enforced]
        return self.printCard(fields)

class SPCD(Constraint):
    """
    Defines an enforced displacement value for static analysis and an enforced motion
    value (displacement, velocity or acceleration) in dynamic analysis.
    SPCD SID G1  C1   D1 G2 C2  D2
    SPCD 100 32 436 -2.6  5    2.9
    """
    type = 'SPCD'
    def __init__(self,card):
        SPC.__init__(self,card)  # defines everything :) at least until cross-referencing methods are implemented

class SPCAX(Constraint):
    """
    Defines a set of single-point constraints or enforced displacements
    for conical shell coordinates.
    SPCAX SID RID HID C    D
    SPCAX 2   3     4 13 6.0
    """
    type = 'SPCAX'
    def __init__(self,card):
        SPC.__init__(self,card)  # defines everything :) at least until cross-referencing methods are implemented

        ## Identification number of a single-point constraint set.
        self.cid = card.field(1)
        ## Ring identification number. See RINGAX entry.
        self.rid = card.field(2)
        ## Harmonic identification number. (Integer >= 0)
        self.hid = card.field(3)
        ## Component identification number. (Any unique combination of the
        ## Integers 1 through 6.)
        self.c   = card.field(4)
        ## Enforced displacement value
        self.d   = card.field(5)

    def crossReference(self,i,node):
        pass

    def __repr__(self): # SPCAX
        fields = ['SPCAX',self.cid,self.rid,self.hid,self.c,self.d]
        return self.printCard(fields)

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
    def __init__(self,card):
        Constraint.__init__(self,card)
        self.constraints = card.field(2)  # 246 = y; dx, dz dir
        nodes = card.fields(3)
        self.nodes = self.expandThru(nodes)
        #print "nodes = ",nodes

    def crossReference(self,i,node):
        self.nodes[i] = node

    def __repr__(self): # SPC1
        #test = [i for i in range(self.nodes[0],self.nodes[-1]+1)]
        #print "self.nodes = ",self.nodes
        #print "test       = ",test
        #if self.nodes==test:
        #    nodes = [self.nodes[0],'THRU',self.nodes[-1]]
        #else:
        #print "SPC1 self.nodes = ",self.nodes
        nodes = [int(i) for i in self.nodes] # SPC1
        fields = ['SPC1',self.cid,self.constraints]+nodes
        return self.printCard(fields)

class SPCADD(Constraint):
    """
    Defines a single-point constraint set as a union of single-point constraint
    sets defined on SPC or SPC1 entries.
    SPCADD   2       1       3
    """
    type = 'SPCADD'
    def __init__(self,card):
        Constraint.__init__(self,card)
        sets = card.fields(2)
        
        self.sets = self.expandThru(sets)
        #print "self.nodes = ",self.nodes

    def crossReference(self,i,node):
        dofCount = 0
        self.sets[i] = node

    def __repr__(self):
        fields = ['SPCADD',self.cid] #+self.sets
        return self._reprSpcMpcAdd(fields)

    def _reprSpcMpcAdd(self,fields):
        outSPCs = ''
        fieldSets = []

        #print "repr---mpcadd"
        #print self.sets
        for spcsets in self.sets:
            #print "*spcsets",spcsets
            if isinstance(spcsets,int):  # spcset wasnt found
                #print "int unfound...%s" %(spcsets)
                #outSPCs += str(spcsets)
                fields.append(spcsets)
            elif isinstance(spcsets,list):
                #print 'list'
                for spcset in spcsets:
                    fieldSets.append(spcset.cid)
                    outSPCs += str(spcset)
                
            else:
                #print 'dict'
                #outSPCs += str(spcsets.cid)
                for key,spcset in spcsets.items():
                    fieldSets.append(spcsets.cid)

        return self.printCard(fields+list(set(fieldSets)))+outSPCs  # SPCADD

class MPCADD(SPCADD):
    """
    Defines a multipoint constraint equation of the form \f$ \Sigma_j A_j u_j =0 \f$
    where \f$ u_j \f$ represents degree-of-freedom \f$ C_j \f$ at grid or scalar point \f$ G_j \f$.
    mPCADD   2       1       3
    """
    type = 'MPCADD'
    def __init__(self,card):
        Constraint.__init__(self,card)
        sets = card.fields(2)
        
        self.sets = self.expandThru(sets)
        #print "self.nodes = ",self.nodes

    def crossReference(self,i,node):
        dofCount = 0
        self.sets[i] = node

    def __repr__(self):
        outSPCs = ''
        fields = ['MPCADD',self.cid] #+self.sets
        return self._reprSpcMpcAdd(fields)

