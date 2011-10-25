# my code
from baseCard import BaseCard


class constraintObject(object):
    def __init__(self):
        self.constraints    = {} # SPC, SPC1, SPCD, etc...
        self.addConstraints = {} # SPCADD
        self.resolvedConstraints = []

    def add(self,constraint):
        self.addConstraints[key] = [constraint]

    def append(self,constraint):
        key = constraint.cid
        if self.constraints.has_key(key):
            self.constraints[key].append(constraint)
        else:
            self.constraints[key] = [constraint]

    def crossReference(self,mesh):
        for key,addConstraint in sorted(self.addConstraints):  # SPCADDs
            nodes = addConstraint.nodes
            for i,node in enumerate(nodes):
                nodes[i] = self.constraints[node]

        for key,constraints in sorted(self.constraints.items()): # SPC, SPC1, SPCD
            for constraint in constraints:
                #if constraint.type=='SPCADD'
                #else:
                (cid,nodeDOFs) = constraint.getNodeDOFs(mesh) # the constrained nodes
                for nodeDOF in nodeDOFs:
                    self.addConstraint(cid,nodeDOF)
                ###
            ###
        ###
    ###

    def addConstraint(self,cid,nodeDOF):
        (nid,dofs) = nodeDOF
        for dof in dofs:
            self.resolvedConstraints.append( (nid,dof) )
        ###
    ###

class Constraint(BaseCard):
    def __init__(self,card):
        self.cid  = card.field(1)

    def cleanNodes(self,nodes):
        """
        nodes are cleaned to get rid of blank fields...which shouldnt be there...
        """
        nodes2 = []
        for node in nodes:
            if node=="":
                pass
            else:
                nodes2.append(node)
            ###
        ###
        #print "nodes2 = ",nodes2
        nodes = nodes2
        if len(nodes2)>1 and nodes2[1]=='THRU':
            nodes = [int(i) for i in range(nodes2[0],nodes2[2]+1)]

        self.nodes = nodes
        #print "*nodes = ",nodes
        #return nodes2

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

class MPC(Constraint): # not done...
    type = 'MPC'
    def __init__(self,card):
        Constraint.__init__(self,card)    # defines self.cid

    def __repr__(self): # MPC
        fields = ['MPC',self.cid]
        for (gid,constraint,enforced) in zip(self.gids,self.constraints,self.enforced):
            fields += [gid,constriant,enforced]
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

    #def getNodeDOFs():
    #    pass

    def crossReference(self,mesh):
        dofCount = 0
        for (i,constraint) in enumerate(self.constraints):
            if self.constraint is None:
                node = self.Node(self.gids[i])
                if not node.Is('GRID'): # SPOINT, EPOINT, DPOINT
                    dofCount+=1
                else:
                    dofCount+=6
                ###
            ###
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
        self.cleanNodes(nodes)
        #print "nodes = ",nodes

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
    Defines a single-point constraint set as a union of single-point constraint sets defined
    on SPC or SPC1 entries.
    
    SPCADD   2       1       3
    
    Defines a single-point constraint set as a union of single-point constraint sets defined
    on SPC or SPC1 entries.
    """
    type = 'SPCADD'
    def __init__(self,card):
        Constraint.__init__(self,card)
        nodes = card.fields(2)
        
        self.cleanNodes(nodes)
        #print "self.nodes = ",self.nodes

    def __repr__(self):
        fields = ['SPCADD',self.cid]+self.nodes
        return self.printCard(fields)
