#import sys

from pyNastran.bdf.cards.baseCard import BaseCard

class Load(BaseCard):
    """defines the DefaultLoad class"""
    type = 'DefLoad'
    def __init__(self,card,data):
        pass
        #self.type = card[0]

    #def normalize(self,v):
    #    #print "v = ",v
    #    return v/norm(v)

    def Cid(self):
        if isinstance(self.cid,int):
            return self.cid
        else:
            return self.cid.cid
        ###

    def nodeIDs(self,nodes=None):
        """returns nodeIDs for repr functions"""
        if not nodes:
           nodes = self.nodes
        if isinstance(nodes[0],int):
            #print 'if'
            return [node     for node in nodes]
        else:
            #print 'else'
            return [node.nid for node in nodes]
        ###

    def rawFields(self):
        fields = [self.type,self.lid]
        return fields

    def reprFields(self):
        return self.rawFields()


class LSEQ(BaseCard): # Requires LOADSET in case control deck
    """
    Defines a sequence of static load sets
    @todo how does this work...
    """
    type = 'LSEQ'
    def __init__(self,card=None,data=None):
        self.sid  = card.field(1)
        self.exciteID = card.field(2)
        self.lid = card.field(3)
        self.tid = card.field(4)

    def nodeIDs(self,nodes=None):
        """returns nodeIDs for repr functions"""
        if not nodes:
           nodes = self.nodes
        if isinstance(nodes[0],int):
            #print 'if'
            return [node     for node in nodes]
        else:
            #print 'else'
            return [node.nid for node in nodes]
        ###

    def crossReference(self,model):
        self.lid = model.Load(self.lid)
        if self.tid:
            self.tid = model.Load(self.tid)
        ###
    
    def Lid(self):
        if isinstance(self.lid,int):
            return self.lid
        return self.lid.lid
        
    def Tid(self):
        if self.tid is None:
            return None
        if isinstance(self.tid,int):
            return self.tid
        return self.tid.tid

    def rawFields(self):
        fields = ['LSEQ',self.sid,self.exciteID,self.Lid(),self.Tid()]
        return fields

    def reprFields(self):
        return self.rawFields()

class DLOAD(Load):
    type = 'DLOAD'
    def __init__(self,card=None,data=None):
        ## load ID
        self.lid   = card.field(1)
        self.scale = card.field(2)

        fields = card.fields(3)
        n = len(fields)//2
        if len(fields)%2==1:
            n+=1
            raise Exception('missing last magnitude on DLOAD card=%s' %(card.fields()) )

        self.sids = []
        self.mags = []
        for i in range(n):
            j = 2*i
            self.mags.append(fields[j  ])
            self.sids.append(fields[j+1])  # RLOADx,TLOADx,ACSRC
        ###

    def crossReference(self,model):
        for (i,sid) in enumerate(self.sids):
            self.sids[i] = model.Load(sid)
        ###

    def Sid(self,sid):
        if isinstance(sid,int):
            return sid
        return sid.lid

    def rawFields(self):
        fields = ['DLOAD',self.lid,self.scale]
        for (mag,sid) in zip(self.mags,self.sids):
            fields += [mag,self.Sid(sid)]
        return fields

    def reprFields(self):
        return self.rawFields()

class SLOAD(Load):
    type = 'SLOAD'
    def __init__(self,card=None,data=None):
        ## load ID
        self.lid = card.field(1)
        
        fields = card.fields(2)
        n = len(fields)//2
        if len(fields)%2==1:
            n+=1
            raise Exception('missing last magnitude on SLOAD card=%s' %(card.fields()) )

        self.sids = []
        self.mags = []
        for i in range(n):
            j = 2*i
            self.sids.append(fields[j  ])
            self.mags.append(fields[j+1])
        ###

    def crossReference(self,model):
        for (i,sid) in enumerate(self.sids):
            self.sids[i] = model.Load(sid)
        ###

    def Sid(self,sid):
        if isinstance(sid,int):
            return sid
        return sid.lid

    def rawFields(self):
        fields = ['SLOAD',self.lid]
        for sid,mag in zip(self.sids,self.mags):
            fields += [self.Sid(sid),mag]
        return fields

    def reprFields(self):
        return self.rawFields()

class RANDPS(BaseCard):
    """
    Power Spectral Density Specification
    Defines load set power spectral density factors for use in random analysis having the
    frequency dependent form
    \f[ S_{jk}(F) = (X+iY)G(F) \f]
    """
    type = 'RANDPS'
    def __init__(self,card=None,data=None):
        if card:
            ## Random analysis set identification number. (Integer > 0)
            ## Defined by RANDOM in the Case Control Deck.
            self.lid = card.field(2)
            ## Subcase identification number of the excited load set. (Integer > 0)
            self.j   = card.field(3)
            ## Subcase identification number of the applied load set. (Integer >= 0; K >= J)
            self.k   = card.field(4)
            ## Components of the complex number. (Real)
            self.x   = card.field(5)
            self.y   = card.field(6)
            ## Identification number of a TABRNDi entry that defines G(F).
            self.tid = card.field(7)

    def crossReference(self,model):
        self.tid = model.Table(self.tid)
    
    def Tid(self):
        if isinstance(self.tid,int):
            return self.tid
        return self.tid.tid

    def rawFields(self):
        fields = [self.lid,self.j,self.k,self.x,self.y,self.Tid()]
        return fields

