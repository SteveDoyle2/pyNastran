# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from math import log, exp
from itertools import izip

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import BaseCard


class FREQ(BaseCard):
    """
    Defines a set of frequencies to be used in the solution of frequency
    response problems.
    FREQ SID F1 F2 F3 F4 F5 F6 F7
    F8 F9 F10
    """
    type = 'FREQ'

    def __init__(self, card=None, data=None):
        self.sid = card.field(1)
        self.freqs = card.fields(2)
        self.cleanFreqs()

    def cleanFreqs(self):
        self.freqs = list(set(self.freqs))
        self.freqs.sort()

    def getFreqs(self):
        return self.freqs

    def addFrequencies(self, freqs):
        """
        Combines the frequencies from 1 FREQx object with another.
        All FREQi entries with the same frequency set identification numbers
        will be used. Duplicate frequencies will be ignored.
        @param self the object pointer
        @param freqs the frequencies for a FREQx object
        """
        #print "self.freqs = ",self.freqs
        #print "freqs = ",freqs
        self.freqs += freqs
        self.cleanFreqs()

    def addFrequencyObject(self, freq):
        """
        @see addFrequencies
        @param self the object pointer
        @param freq a FREQx object
        """
        self.addFrequencies(freq.freqs)

    def rawFields(self):
        fields = ['FREQ', self.sid] + self.freqs
        return fields


class FREQ1(FREQ):
    """
    Defines a set of frequencies to be used in the solution of frequency
    response problems by specification of a starting frequency, frequency
    increment, and the number of increments desired.
    FREQ1 SID F1 DF NDF
    @note this card rewrites as a FREQ card
    """
    type = 'FREQ1'

    def __init__(self, card=None, data=None):
        self.sid = card.field(1)
        f1 = card.field(2, 0.0)
        df = card.field(3)
        ndf = card.field(4, 1)

        self.freqs = []
        for i in xrange(ndf):
            self.freqs.append(f1 + i * df)
        ###
        self.cleanFreqs()


class FREQ2(FREQ):
    """
    Defines a set of frequencies to be used in the solution of frequency
    response problems by specification of a starting frequency, final
    frequency, and the number of logarithmic increments desired.
    FREQ2 SID F1 F2 NDF
    @note this card rewrites as a FREQ card
    """
    type = 'FREQ2'

    def __init__(self, card=None, data=None):
        self.sid = card.field(1)
        f1 = card.field(2, 0.0)
        f2 = card.field(3)
        nf = card.field(4, 1)

        d = 1. / nf * log(f2 / f1)
        self.freqs = []
        for i in xrange(nf):
            self.freqs.append(f1 * exp(i * d))  # 0 based index
        self.cleanFreqs()

#class FREQ3(FREQ):


class FREQ4(FREQ):
    """
    Defines a set of frequencies used in the solution of modal frequency
    response problems by specifying the amount of 'spread' around each natural
    frequency and the number of equally spaced excitation frequencies within
    the spread.
    FREQ4 SID F1 F2 FSPD NFM
    @note this card rewrites as a FREQ card
    @todo not done...
    """
    type = 'FREQ4'

    def __init__(self, card=None, data=None):
        self.sid = card.field(1)
        self.f1 = card.field(2, 0.0)
        self.f2 = card.field(3, 1.e20)
        self.fspd = card.field(4, 0.1)
        self.nfm = card.field(5, 3)

    def rawFields(self):
        fields = ['FREQ4', self.sid, self.f1, self.f2, self.fspd, self.nfm]
        return fields

    def reprFields(self):
        return self.rawFields()


#class FREQ5(FREQ):

class TSTEP(BaseCard):
    """
    Transient Time Step
    Defines time step intervals at which a solution will be generated and
    output in transient analysis.
    """
    type = 'TSTEP'

    def __init__(self, card=None, data=None):
        self.sid = card.field(1)
        self.N = []
        self.DT = []
        self.NO = []
        fields = card.fields(1)

        i = 1
        nFields = len(fields)
        while i < nFields:
            self.N.append(card.field(i + 1, 1))
            self.DT.append(card.field(i + 2, 0.))
            self.NO.append(card.field(i + 3, 1))
            i += 8

    def rawFields(self):
        fields = ['TSTEP', self.sid]
        for (n, dt, no) in izip(self.N, self.DT, self.NO):
            fields += [n, dt, no, None, None, None, None, None]
        return fields

    def reprFields(self):
        return self.rawFields()


class TSTEPNL(BaseCard):
    """
    Defines parametric controls and data for nonlinear transient structural or
    heat transfer analysis. TSTEPNL is intended for SOLs 129, 159, and 600.
    Parameters for Nonlinear Transient Analysis
    """
    type = 'TSTEPNL'

    def __init__(self, card=None, data=None):
        if card:
            self.sid = card.field(1)
            self.ndt = card.field(2)
            self.dt = card.field(3)
            self.no = card.field(4, 1)
            ## @note not listed in all QRGs
            self.method = card.field(5, 'ADAPT')
            if self.method == 'ADAPT':
                self.kStep = card.field(6, 2)
            elif self.method == 'ITER':
                self.kStep = card.field(6, 10)
            elif self.method in ['AUTO', 'TSTEP']:
                self.kStep = card.field(6)
            else:
                msg = 'invalid TSTEPNL Method.  method=|%s|' % (self.method)
                raise RuntimeError(msg)
            self.maxIter = card.field(7, 10)
            self.conv = card.field(8, 'PW')

            # line 2
            self.epsU = card.field(9, 1.E-2)
            self.epsP = card.field(10, 1.E-3)
            self.epsW = card.field(11, 1.E-6)
            self.maxDiv = card.field(12, 2)
            self.maxQn = card.field(13, 10)
            self.MaxLs = card.field(14, 2)
            self.fStress = card.field(15, 0.2)

            # line 3
            self.maxBisect = card.field(17, 5)
            self.adjust = card.field(18, 5)
            self.mStep = card.field(19)
            self.rb = card.field(20, 0.6)
            self.maxR = card.field(21, 32.)
            self.uTol = card.field(22, 0.1)
            self.rTolB = card.field(23, 20.)
            self.minIter = card.field(24)  # not listed in all QRGs
        else:
            (sid, ndt, dt, no, method, kStep, maxIter, conv, epsU, epsP, epsW,
             maxDiv, maxQn, maxLs, fStress, maxBisect,
             adjust, mStep, rb, maxR, uTol, rTolB) = data
            self.sid = sid
            self.ndt = ndt
            self.dt = dt
            self.no = no
            self.method = method
            self.kStep = kStep
            self.maxIter = maxIter
            self.conv = conv

            # line 2
            self.epsU = epsU
            self.epsP = epsP
            self.epsW = epsW
            self.maxDiv = maxDiv
            self.maxQn = maxQn
            self.MaxLs = maxLs
            self.fStress = fStress

            # line 3
            self.maxBisect = maxBisect
            self.adjust = adjust
            self.mStep = mStep
            self.rb = rb
            self.maxR = maxR
            self.uTol = uTol
            self.rTolB = rTolB
            self.minIter = None  # not listed in DMAP 2005
        ###

    def rawFields(self):
        fields = ['TSTEPNL', self.sid, self.ndt, self.dt, self.no, self.method, self.kStep, self.maxIter, self.conv,
                  self.epsU, self.epsP, self.epsW, self.maxDiv, self.maxQn, self.MaxLs, self.fStress, None,
                  self.maxBisect, self.adjust, self.mStep, self.rb, self.maxR, self.uTol, self.rTolB, self.minIter]
        return fields

    def reprFields(self):
        #no      = set_blank_if_default(self.no,1)
        no = self.no
        method = set_blank_if_default(self.method, 'ADAPT')

        kStep = self.kStep
        #if self.method=='ADAPT':
            #kStep = set_blank_if_default(self.kStep, 2)
        #elif self.method=='ITER':
            #kStep = set_blank_if_default(self.kStep, 10)
        #else:
            #msg = 'invalid TSTEPNL Method.  method=|%s|' %(self.method)
            #raise RuntimeError(msg)

        #maxIter = set_blank_if_default(self.maxIter, 10)
        conv = set_blank_if_default(self.conv, 'PW')

        epsU = set_blank_if_default(self.epsU, 1e-2)
        epsP = set_blank_if_default(self.epsP, 1e-3)
        epsW = set_blank_if_default(self.epsW, 1e-6)
        maxDiv = set_blank_if_default(self.maxDiv, 2)
        maxQn = set_blank_if_default(self.maxQn, 10)
        MaxLs = set_blank_if_default(self.MaxLs, 2)
        fStress = set_blank_if_default(self.fStress, 0.2)

        maxBisect = set_blank_if_default(self.maxBisect, 5)
        adjust = set_blank_if_default(self.adjust, 5)
        rb = set_blank_if_default(self.rb, 0.6)
        maxR = set_blank_if_default(self.maxR, 32.)
        uTol = set_blank_if_default(self.uTol, 0.1)
        rTolB = set_blank_if_default(self.rTolB, 20.)

        fields = ['TSTEPNL', self.sid, self.ndt, self.dt, no, method, kStep, self.maxIter, conv,
                  epsU, epsP, epsW, maxDiv, maxQn, MaxLs, fStress, None,
                  maxBisect, adjust, self.mStep, rb, maxR, uTol, rTolB, self.minIter]
        return fields


class NLPARM(BaseCard):
    """
    Defines a set of parameters for nonlinear static analysis iteration
    strategy.
    NLPARM ID NINC DT KMETHOD KSTEP MAXITER CONV INTOUT
    EPSU EPSP EPSW MAXDIV MAXQN MAXLS FSTRESS LSTOL
    MAXBIS MAXR RTOLB
    """
    type = 'NLPARM'

    def __init__(self, card=None, data=None):
        if card:
            self.nid = card.field(1)
            self.ninc = card.field(2, 10)
            self.dt = card.field(3, 0.0)
            self.kMethod = card.field(4, 'AUTO')
            self.kStep = card.field(5, 5)
            self.maxIter = card.field(6, 25)
            self.conv = card.field(7, 'PW')
            self.intOut = card.field(8, 'NO')

            # line 2
            self.epsU = card.field(9, 0.01)
            self.epsP = card.field(10, 0.01)
            self.epsW = card.field(11, 0.01)
            self.maxDiv = card.field(12, 3)
            self.maxQn = card.field(13, self.maxIter)
            self.maxLs = card.field(14, 4)
            self.fStress = card.field(15, 0.2)
            self.lsTol = card.field(16, 0.5)

            # line 3
            self.maxBisect = card.field(17, 5)
            self.maxR = card.field(21, 20.)
            self.rTolB = card.field(23, 20.)
        else:
            (sid, ninc, dt, kMethod, kStep, maxIter, conv, intOut, epsU, epsP, epsW,
             maxDiv, maxQn, maxLs, fStress, lsTol, maxBisect, maxR, rTolB) = data
            self.nid = sid
            self.ninc = ninc
            self.dt = dt
            self.kMethod = kMethod
            self.kStep = kStep
            self.maxIter = maxIter
            self.conv = conv
            self.intOut = intOut

            # line 2
            self.epsU = epsU
            self.epsP = epsP
            self.epsW = epsW
            self.maxDiv = maxDiv
            self.maxQn = maxQn
            self.maxLs = maxLs
            self.fStress = fStress
            self.lsTol = lsTol

            # line 3
            self.maxBisect = maxBisect
            self.maxR = maxR
            self.rTolB = rTolB
        ###

    def rawFields(self):
        fields = ['NLPARM', self.nid, self.ninc, self.dt, self.kMethod, self.kStep, self.maxIter, self.conv, self.intOut,
                  self.epsU, self.epsP, self.epsW, self.maxDiv, self.maxQn, self.maxLs, self.fStress, self.lsTol,
                  self.maxBisect, None, None, None, self.maxR, None, self.rTolB]
        return fields

    def reprFields(self):
        ninc = set_blank_if_default(self.ninc, 10)
        dt = set_blank_if_default(self.dt, 0.0)
        kMethod = set_blank_if_default(self.kMethod, 'AUTO')
        kStep = set_blank_if_default(self.kStep, 5)
        maxIter = set_blank_if_default(self.maxIter, 25)
        conv = set_blank_if_default(self.conv, 'PW')
        intOut = set_blank_if_default(self.intOut, 'NO')
        epsU = set_blank_if_default(self.epsU, 0.01)
        epsP = set_blank_if_default(self.epsP, 0.01)
        epsW = set_blank_if_default(self.epsW, 0.01)
        maxDiv = set_blank_if_default(self.maxDiv, 3)
        maxQn = set_blank_if_default(self.maxQn, self.maxIter)
        maxLs = set_blank_if_default(self.maxLs, 4)
        fStress = set_blank_if_default(self.fStress, 0.2)
        lsTol = set_blank_if_default(self.lsTol, 0.5)
        maxBisect = set_blank_if_default(self.maxBisect, 5)
        maxR = set_blank_if_default(self.maxR, 20.)
        rTolB = set_blank_if_default(self.rTolB, 20.)

        fields = ['NLPARM', self.nid, ninc, dt, kMethod, kStep, maxIter, conv, intOut,
                  epsU, epsP, epsW, maxDiv, maxQn, maxLs, fStress, lsTol,
                  maxBisect, None, None, None, maxR, None, rTolB]
        return fields
