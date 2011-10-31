from baseCard import BaseCard

class RLOAD1(BaseCard): # not integrated
    """
    Defines a frequency-dependent dynamic load of the form
    for use in frequency response problems.
    RLOAD1 SID EXCITEID DELAY DPHASE TC TD TYPE
    \f[ \large \left{ P(f)  \right}  = \left{A\right} [ C(f)+iD(f)] e^{  i \left{\theta - 2 \pi f \tau \right} }
    RLOAD1 5   3                     1
    """
    type = 'RLOAD1'
    def __init__(self,card):
        self.sid      = card.field(1)
        self.exciteID = card.field(2)
        self.delay    = card.field(3)
        self.dphase   = card.field(4)
        self.tc    = card.field(5)
        self.td    = card.field(6)
        self.Type  = card.field(7,0)

        if   self.Type in [0,'L','LO'.'LOA','LOAD']: self.Type = 'LOAD'
        elif self.Type in [1,'D','DI','DIS','DIPS']: self.Type = 'DISP'
        elif self.Type in [2,'V','VE','VEL','VELO']: self.Type = 'VELO'
        elif self.Type in [3,'A','AC','ACC','ACCE']: self.Type = 'ACCE'
        else: raise RuntimeError('invalid RLOADi type  Type=|%s|' %(self.Type))

    def __repr__(self):
        Type = self.setBlankIfDefault(self.Type,'LOAD')
        fields = ['RLOAD1',self.sid,self.exciteID,self.delay,self.dphase,self.tc,self.td,Type]
        return self.printCard(fields)

class RLOAD2(BaseCard): # not integrated
    """
    Defines a frequency-dependent dynamic load of the form
    for use in frequency response problems.
    
    \f[ \large \left{ P(f)  \right}  = \left{A\right} * B(f) e^{  i \left{ \phi(f) + \theta - 2 \pi f \tau \right} }
    RLOAD2 SID EXCITEID DELAY DPHASE TB TP TYPE
    RLOAD2 5   3                     1
    """
    type = 'RLOAD2'
    def __init__(self,card):
        self.sid      = card.field(1)
        self.exciteID = card.field(2)
        self.delay    = card.field(3)
        self.dphase   = card.field(4)
        self.tb    = card.field(5)
        self.tp    = card.field(6)
        self.Type  = card.field(7,0)

        if   self.Type in [0,'L','LO'.'LOA','LOAD']: self.Type = 'LOAD'
        elif self.Type in [1,'D','DI','DIS','DIPS']: self.Type = 'DISP'
        elif self.Type in [2,'V','VE','VEL','VELO']: self.Type = 'VELO'
        elif self.Type in [3,'A','AC','ACC','ACCE']: self.Type = 'ACCE'
        else: raise RuntimeError('invalid RLOADi type  Type=|%s|' %(self.Type))

    def __repr__(self):
        Type = self.setBlankIfDefault(self.Type,0.0)
        fields = ['RLOAD2',self.sid,self.exciteID,self.delay,self.dphase,self.tb,self.tp,Type]
        return self.printCard(fields)
