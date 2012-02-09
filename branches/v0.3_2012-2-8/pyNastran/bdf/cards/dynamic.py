from baseCard import BaseCard

class FREQ(BaseCard):
    """
    Defines a set of frequencies to be used in the solution of frequency response problems.
    FREQ SID F1 F2 F3 F4 F5 F6 F7
    F8 F9 F10
    """
    type = 'FREQ'
    def __init__(self,card=None,data=None):
        self.sid = card.field(1)
        self.freqs = card.fields(2)
        self.cleanFreqs()
    
    def cleanFreqs(self):
        self.freqs = list(set(self.freqs))
        self.freqs.sort()

    def getFreqs(self):
        return self.freqs
    
    def addFrequencies(self,freqs):
        """
        Combines the frequencies from 1 FREQx object with another.
        All FREQi entries with the same frequency set identification numbers will be used. Duplicate
        frequencies will be ignored.
        @param self the object pointer
        @param freqs the frequencies for a FREQx object
        """
        #print "self.freqs = ",self.freqs
        #print "freqs = ",freqs
        self.freqs += freqs
        self.cleanFreqs()

    def addFrequencyObject(self,freq):
        """
        @see addFrequencies
        @param self the object pointer
        @param freq a FREQx object
        """
        self.addFrequencies(freq.freqs)

    def rawFields(self):
        fields = ['FREQ',self.sid]+self.freqs
        return fields

class FREQ1(FREQ):
    """
    Defines a set of frequencies to be used in the solution of frequency response problems
    by specification of a starting frequency, frequency increment, and the number of
    increments desired.
    FREQ1 SID F1 DF NDF
    @note this card rewrites as a FREQ card
    """
    type = 'FREQ1'
    def __init__(self,card=None,data=None):
        self.sid = card.field(1)
        f1  = card.field(2,0.0)
        df  = card.field(3)
        ndf = card.field(4,1)
        
        self.freqs = []
        for i in range(ndf):
            self.freqs.append(f1+i*df)
        ###
        self.cleanFreqs()

class FREQ2(FREQ):
    """
    Defines a set of frequencies to be used in the solution of frequency response problems
    by specification of a starting frequency, final frequency, and the number of
    logarithmic increments desired.
    FREQ2 SID F1 F2 NDF
    @note this card rewrites as a FREQ card
    """
    type = 'FREQ2'
    def __init__(self,card=None,data=None):
        self.sid = card.field(1)
        f1 = card.field(2,0.0)
        f2 = card.field(3)
        nf = card.field(4,1)
        
        d = 1./nf*log(f2/f1)
        self.freqs = []
        for i in range(nf):
            self.freqs.append(f1*exp(i*d)) # 0 based index
        ###
        self.cleanFreqs()

#class FREQ3(FREQ):

class FREQ4(FREQ):
    """
    Defines a set of frequencies used in the solution of modal frequency-response
    problems by specifying the amount of 'spread' around each natural frequency and
    the number of equally spaced excitation frequencies within the spread.
    FREQ4 SID F1 F2 FSPD NFM
    @note this card rewrites as a FREQ card
    @todo not done...
    """
    type = 'FREQ4'
    def __init__(self,card=None,data=None):
        self.sid = card.field(1)
        f1   = card.field(2,0.0)
        f2   = card.field(3,1.e20)
        fspd = card.field(4,0.1)
        nfm  = card.field(5,3)

#class FREQ5(FREQ):

class RLOAD1(BaseCard): # not integrated
    """
    Defines a frequency-dependent dynamic load of the form
    for use in frequency response problems.
    RLOAD1 SID EXCITEID DELAY DPHASE TC TD TYPE
    \f[ \large \left{ P(f)  \right}  = \left{A\right} [ C(f)+iD(f)] e^{  i \left{\theta - 2 \pi f \tau \right} } \f]
    RLOAD1 5   3                     1
    """
    type = 'RLOAD1'
    def __init__(self,card=None,data=None):
        self.sid      = card.field(1)
        self.exciteID = card.field(2)
        self.delay    = card.field(3)
        self.dphase   = card.field(4)
        self.tc    = card.field(5)
        self.td    = card.field(6)
        self.Type  = card.field(7,0)

        if   self.Type in [0,'L','LO','LOA','LOAD']: self.Type = 'LOAD'
        elif self.Type in [1,'D','DI','DIS','DIPS']: self.Type = 'DISP'
        elif self.Type in [2,'V','VE','VEL','VELO']: self.Type = 'VELO'
        elif self.Type in [3,'A','AC','ACC','ACCE']: self.Type = 'ACCE'
        else: raise RuntimeError('invalid RLOADi type  Type=|%s|' %(self.Type))

    def rawFields(self):
        fields = ['RLOAD1',self.sid,self.exciteID,self.delay,self.dphase,self.tc,self.td,self.Type]
        return fields

    def reprFields(self):
        Type = self.setBlankIfDefault(self.Type,'LOAD')
        fields = ['RLOAD1',self.sid,self.exciteID,self.delay,self.dphase,self.tc,self.td,Type]
        return fields

class RLOAD2(BaseCard): # not integrated
    """
    Defines a frequency-dependent dynamic load of the form
    for use in frequency response problems.
    
    \f[ \large \left{ P(f)  \right}  = \left{A\right} * B(f) e^{  i \left{ \phi(f) + \theta - 2 \pi f \tau \right} } \f]
    RLOAD2 SID EXCITEID DELAY DPHASE TB TP TYPE
    RLOAD2 5   3                     1
    """
    type = 'RLOAD2'
    def __init__(self,card=None,data=None):
        self.sid      = card.field(1)
        self.exciteID = card.field(2)
        self.delay    = card.field(3)
        self.dphase   = card.field(4)
        self.tb    = card.field(5)
        self.tp    = card.field(6)
        self.Type  = card.field(7,0)

        if   self.Type in [0,'L','LO','LOA','LOAD']: self.Type = 'LOAD'
        elif self.Type in [1,'D','DI','DIS','DIPS']: self.Type = 'DISP'
        elif self.Type in [2,'V','VE','VEL','VELO']: self.Type = 'VELO'
        elif self.Type in [3,'A','AC','ACC','ACCE']: self.Type = 'ACCE'
        else: raise RuntimeError('invalid RLOADi type  Type=|%s|' %(self.Type))

    def rawFields(self):
        fields = ['RLOAD2',self.sid,self.exciteID,self.delay,self.dphase,self.tb,self.tp,self.Type]
        return fields

    def reprFields(self):
        Type = self.setBlankIfDefault(self.Type,0.0)
        fields = ['RLOAD2',self.sid,self.exciteID,self.delay,self.dphase,self.tb,self.tp,Type]
        return fields

class NLPARM(BaseCard):
    """
    Defines a set of parameters for nonlinear static analysis iteration strategy
    NLPARM ID NINC DT KMETHOD KSTEP MAXITER CONV INTOUT
    EPSU EPSP EPSW MAXDIV MAXQN MAXLS FSTRESS LSTOL
    MAXBIS MAXR RTOLB
    """
    type = 'NLPARM'
    def __init__(self,card=None,data=None):
        self.nid       = card.field(1)
        self.ninc      = card.field(2,10)
        self.dt        = card.field(3,0.0)
        self.kMethod   = card.field(4,'AUTO')
        self.kStep     = card.field(5,5)
        self.maxIter   = card.field(6,25)
        self.conv      = card.field(7,'PW')
        self.intOut    = card.field(7,'NO')
        self.epsU      = card.field(9,1e-2)
        self.epsP      = card.field(10,1e-2)
        self.epsW      = card.field(11,1e-2)
        self.maxDiv    = card.field(12,3)
        self.maxQn     = card.field(13,self.maxIter)
        self.maxLs     = card.field(14,4)
        self.fStress   = card.field(15,0.2)
        self.lsTol     = card.field(16,0.5)
        self.maxBisect = card.field(17,5)
        self.maxR      = card.field(21,20.)
        self.rTolB     = card.field(23,20.)

    def rawFields(self):
        fields = ['NLPARM',self.nid,self.ninc,self.dt,self.kMethod,self.kStep,self.maxIter,self.conv,self.intOut,
                           self.epsU,self.epsP,self.epsW,self.maxDiv,self.maxQn,self.maxLs,self.fStress,self.lsTol,
                           self.maxBisect,None,None,None,self.maxR,None,self.rTolB]
        return fields

    def reprFields(self):
        ninc      = self.setBlankIfDefault(self.ninc,10)
        dt        = self.setBlankIfDefault(self.dt,0.0)
        kMethod   = self.setBlankIfDefault(self.kMethod,'AUTO')
        kStep     = self.setBlankIfDefault(self.kStep,5)
        maxIter   = self.setBlankIfDefault(self.maxIter,25)
        conv      = self.setBlankIfDefault(self.conv,'PW')
        intOut    = self.setBlankIfDefault(self.intOut,'NO')
        epsU      = self.setBlankIfDefault(self.epsU,1e-2)
        epsP      = self.setBlankIfDefault(self.epsP,1e-2)
        epsW      = self.setBlankIfDefault(self.epsW,1e-2)
        maxDiv    = self.setBlankIfDefault(self.maxDiv,3)
        maxQn     = self.setBlankIfDefault(self.maxQn,self.maxIter)
        maxLs     = self.setBlankIfDefault(self.maxLs,4)
        fStress   = self.setBlankIfDefault(self.fStress,0.2)
        lsTol     = self.setBlankIfDefault(self.lsTol,0.5)
        maxBisect = self.setBlankIfDefault(self.maxBisect,5)
        maxR      = self.setBlankIfDefault(self.maxR,20.)
        rTolB     = self.setBlankIfDefault(self.rTolB,20.)

        fields = ['NLPARM',self.nid,ninc,dt,kMethod,kStep,maxIter,conv,intOut,
                           epsU,epsP,epsW,maxDiv,maxQn,maxLs,fStress,lsTol,maxBisect,None,None,None,
                           maxR,None,rTolB]
        return fields



