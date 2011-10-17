class CardInstantiationError(RuntimeError):
    pass

from elements import Element

class CROD(Element):
    type = 'CROD'
    def __init__(self,card):
        Element.__init__(self,card)
        self.id  = card.field(1)
        self.pid = card.field(2,self.id)

        nids = card.fields(3,5)
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==2

    def __repr__(self):
        fields = [self.type,self.eid,self.pid]+self.nodes
        return self.printCard(fields)

class CTUBE(CROD):
    type = 'CTUBE'
    def __init__(self,card):
        CROD.__init__(self,card)
    ###
###


class CBAR(Element):
    """
    CBAR EID PID GA GB X1 X2 X3 OFFT
    PA PB W1A W2A W3A W1B W2B W3B
    or
    CBAR EID PID GA GB G0 OFFT
    PA PB W1A W2A W3A W1B W2B W3B

    CBAR       22062       4   21798   21799   0.0+0   0.0+0     -1.
                               0.0+0   0.0+0     -9.   0.0+0   0.0+0     -9.
    """
    type = 'CBAR'
    def __init__(self,card):
        Element.__init__(self,card)
        self.pid = card.field(2)
        self.ga = card.field(3)
        self.gb = card.field(4)
        self.initX_G0()

        self.offt = card.field(8,'GGG')
        assert self.offt[0] in ['G','B','O'],'invalid offt parameter of CBEAM...offt=%s' %(self.offt)
        assert self.offt[1] in ['G','B','O'],'invalid offt parameter of CBEAM...offt=%s' %(self.offt)
        assert self.offt[2] in ['G','B','O'],'invalid offt parameter of CBEAM...offt=%s' %(self.offt)

        self.pa = card.field(9)
        self.pb = card.field(10)
        
        self.w1a = card.field(11,0.0)
        self.w2a = card.field(12,0.0)
        self.w3a = card.field(13,0.0)

        self.w1b = card.field(14,0.0)
        self.w2b = card.field(15,0.0)
        self.w3b = card.field(16,0.0)

    def initX_G0(self,card):
        field5 = card.field(5)
        if isinstance(field5,int):
            self.g0 = field5
            self.x1 = None
            self.x2 = None
            self.x3 = None
        elif isinstance(field5,float):
            self.g0 = None
            self.x1 = card.field(5,0.0)
            self.x2 = card.field(6,0.0)
            self.x3 = card.field(7,0.0)
        else:
            msg = 'field5 on cbar is the wrong type...field5=%s' %(field5)
            raise RuntimeError(msg)

    def crossReference(self,nodes):
        """
        set g0-ga to x1,x2,x3
        """
        if self.g0:
            v = nodes[self.g0].Position()-nodes[self.ga].Position()
            self.x1 = v[0]
            self.x2 = v[1]
            self.x3 = v[2]
        ###

    def updateNodes(self,nodes):
        """@todo maybe improve"""
        self.crossReference(self,nodes)

    def getX_G0_defaults(self):
        if self.g0:
            x1 = self.setBlankIfDefault(self.x1,0.0)
            x2 = self.setBlankIfDefault(self.x2,0.0)
            x3 = self.setBlankIfDefault(self.x3,0.0)
            return (x1,x2,x3)
        else:
            return (self.g0,None,None)
        ###

    def __repr__(self):
        w1a = self.setBlankIfDefault(self.w1a,0.0)
        w2a = self.setBlankIfDefault(self.w2a,0.0)
        w3a = self.setBlankIfDefault(self.w3a,0.0)
        w1b = self.setBlankIfDefault(self.w1b,0.0)
        w2b = self.setBlankIfDefault(self.w2b,0.0)
        w3b = self.setBlankIfDefault(self.w3b,0.0)
        (x1,x2,x3) = self.getX_G0_defaults()
        offt = self.setBlankIfDefault(self.offt,'GGG')
        fields = [self.type,self.eid,self.pid,self.ga,self.gb,x1,x2,x3,offt,
                  self.pa,self.pb,w1a,w2a,w3a,w1b,w2b,w3b]

        return self.printCard(fields)

class CBEAM(CBAR):
    """
    CBEAM EID PID GA GB X1 X2 X3 OFFT/BIT
    PA PB W1A W2A W3A W1B W2B W3B
    SA SB
    or
    CBEAM EID PID GA GB G0 - - OFFT/BIT
    PA PB W1A W2A W3A W1B W2B W3B
    SA SB
    """
    type = 'CBEAM'
    def __init__(self,card):
        Element.__init__(self,card)
        self.pid = card.field(2,self.id)
        self.ga = card.field(3)
        self.gb = card.field(4)

        self.initX_G0(card)
        self.initOfftBit(card)
        self.pa = card.field(9)
        self.pb = card.field(10)
        
        self.w1a = card.field(11,0.0)
        self.w2a = card.field(12,0.0)
        self.w3a = card.field(13,0.0)

        self.w1b = card.field(14,0.0)
        self.w2b = card.field(15,0.0)
        self.w3b = card.field(16,0.0)

        self.sa = card.field(17)
        self.sb = card.field(18)
    
    def initOfftBit(self,card):
        field8 = card.field(8)
        if isinstance(field8,float):
            self.bit = field8
            self.offt = None
        elif field8 is None:
            self.offt = 'GGG' # default
            self.bit = None
        elif isinstance(field8,str):
            self.offt = field8
            self.bit = None
            assert self.offt[0] in ['G','B','O'],'invalid offt parameter of CBEAM...offt=%s' %(self.offt)
            assert self.offt[1] in ['G','B','O'],'invalid offt parameter of CBEAM...offt=%s' %(self.offt)
            assert self.offt[2] in ['G','B','O'],'invalid offt parameter of CBEAM...offt=%s' %(self.offt)
        else:
            msg = 'field8 on CBEAM card is not a string(offt) or bit (float)...field8=%s\n' %(field8)
            raise CardInstantiationError(msg)
        ###

    def getOfft_Bit_defaults(self):
        if self.isOfft:
            field8 = self.offt
        else:
            field8 = self.setBlankIfDefault(self.bit,0.0)
        return field8

    def __repr__(self):
        w1a = self.setBlankIfDefault(self.w1a,0.0)
        w2a = self.setBlankIfDefault(self.w2a,0.0)
        w3a = self.setBlankIfDefault(self.w3a,0.0)
        w1b = self.setBlankIfDefault(self.w1b,0.0)
        w2b = self.setBlankIfDefault(self.w2b,0.0)
        w3b = self.setBlankIfDefault(self.w3b,0.0)
        (x1,x2,x3) = self.getX_G0_defaults()
        offt = self.getOfft_Bit_defaults()
        fields = [self.type,self.eid,self.pid,self.ga,self.gb,x1,x2,x3,offt,
                  self.pa,self.pb,w1a,w2a,w3a,w1b,w2b,w3b,
                  self.sa,self.sb]
        return self.printCard(fields)
