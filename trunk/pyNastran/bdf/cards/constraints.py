# my code
from baseCard import BaseCard

class Constraint(BaseCard):
    def __init__(self,card):
        #self.type = card.field(0)
        self.id  = card.field(1)

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
        fields = [self.type,self.id]
        return self.printCard(fields)

class SUPORT1(Constraint):
    """
    #SUPORT1 SID ID1 C1 ID2 C2 ID3 C3
    """
    type = 'SUPORT1'
    def __init__(self,card):
        Constraint.__init__(self,card)
        self.sid = card.field(1)
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
        fields = [self.type,self.sid]
        for ID,c in zip(self.IDs,self.Cs):
            fields += [ID,c]
        return self.printCard(fields)

class SPC1(Constraint):
    #SPC1     3       246     209075  209096  209512  209513  209516
    type = 'SPC1'
    def __init__(self,card):
        Constraint.__init__(self,card)
        self.constrained = card.field(2)  # 246 = y; dx, dz dir
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
        print "SPC1 self.nodes = ",self.nodes
        nodes = [int(i) for i in self.nodes] # SPC1
        fields = [self.type,self.id,self.constrained]+nodes
        return self.printCard(fields)

class SPCADD(Constraint):
    #SPCADD   2       1       3
    type = 'SPCADD'
    def __init__(self,card):
        Constraint.__init__(self,card)
        nodes = card.fields(2)
        self.cleanNodes(nodes)
        #print "self.nodes = ",self.nodes

    def __repr__(self):
        fields = [self.type,self.id]+self.nodes
        return self.printCard(fields)
