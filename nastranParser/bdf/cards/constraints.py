# my code
from baseCard import BaseCard

class Constraint(BaseCard):
    def __init__(self,card):
        #self.type = card[0]
        self.lid  = card[1]

    def cleanNodes(self,nodes):
        """
        why do the nodes need to be cleaned?
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

class SPC1(Constraint):
    #SPC1     3       246     209075  209096  209512  209513  209516
    type = 'SPC1'
    def __init__(self,card):
        Constraint.__init__(self,card)
        self.constrained = card[2]  # 246 = y; dx, dz dir
        nodes = card[3:]
        self.cleanNodes(nodes)
        #print "nodes = ",nodes

    def __repr__(self):
        #test = [i for i in range(self.nodes[0],self.nodes[-1]+1)]
        #print "self.nodes = ",self.nodes
        #print "test       = ",test
        #if self.nodes==test:
        #    nodes = [self.nodes[0],'THRU',self.nodes[-1]]
        #else:
        nodes = [int(i) for i in self.nodes]
        fields = [self.type,self.lid,self.constrained]+nodes
        return self.printCard(fields)

class SPCADD(Constraint):
#SPCADD   2       1       3
    type = 'SPCADD'
    def __init__(self,card):
        Constraint.__init__(self,card)
        nodes = card[2:]
        self.cleanNodes(nodes)
        #print "self.nodes = ",self.nodes

    def __repr__(self):
        fields = [self.type,self.lid]+self.nodes
        return self.printCard(fields)
