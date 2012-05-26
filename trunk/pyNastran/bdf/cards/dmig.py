import os
import sys

from pyNastran.bdf.cards.baseCard import BaseCard

from math import sin,sinh,cos,cosh,tan,tanh,sqrt,atan,atan2 #,acosh,acos,asin,asinh,atanh #,atanh2
from numpy import average

def ssq(*listA):
    out = 0.
    for x in listA:
        out += x*x
    return out

def sum2(*listA):
    return sum(listA)

def mod(x,y):
    return x%y

def logx(x,y):
    log(y,x)

def dim(x,y):
    return x-min(x,y)

def db(p,pref):
    """
    sound pressure in decibels
    would capitalize it, but you wouldnt be able to call the function...
    """
    return 20.*log(p/pref)

class DEQATN(BaseCard):# needs work...
    type = 'DEQATN'
    def __init__(self,card=None,data=None):
        newCard = ''
        foundNone = False
        for field in card.card:
            if foundNone is False and field is not None:
                newCard += field+','
                foundNone = True
            elif foundNone is True and field is not None:
                newCard += field

        #if len(card.card)>1:
        #    print "card.card = ",card.card
        #    line0 = ','.join(card.card)
        #else:
        #    line0 = ''.join(card.card)
        line0 = newCard
        self.eqID = line0[8:16]
        
        assert len(self.eqID)==8,'len(eqID)==%s' %(len(self.eqID))
        eq = line0[16:]
        eq = eq.replace(' ','').lower()
        (self.name,self.eq) = eq.split('=')
        print "EQ = ",self.eq
        
    def evaluate(args):
        #eqLow = self.eq.lower()
        #eval(self.eq)
        pass

    def __repr__(self):
        eq = self.name+'='+self.eq
        eqLine = eq[0:56]
        eq = eq[56:]
        fields = ['DEQATN  ','%8s'%(self.eqID),eqLine]

        if len(eq):
            eqLine = eq[0:72]
            eq = eq[72:]
            fields += ['        '+eqLine]
        return ''.join(fields)

class DMIG(BaseCard):  # not done...
    """
    Defines direct input matrices related to grid, extra, and/or scalar points. The matrix
    is defined by a single header entry and one or more column entries. A column entry
    is required for each column with nonzero elements.
    """
    type = 'DMIG'
    def addColumn(self,card=None,data=None):
        #print "hi column"
        Gj = card.field(2)
        Cj = card.field(3)
        #
        
        nFields = card.nFields()
        nLoops = (nFields-5)//4
        minLoop = nLoops-1
        if minLoop==0:
            minLoop = 1
        #assert nFields <= 8,'nFields=%s' %(nFields)

        print "minLoop = ",minLoop
        print "nLoops  = ",nLoops
        for i in range(nLoops):
            self.GCj.append((Gj,Cj))

        for i in range(nLoops):
            n = 5+4*i
            Gi = card.field(n)
            Ci = card.field(n+1)
            self.GCi.append((Gi,Ci))

            self.Real.append(card.field(n+2))
            self.Complex.append(card.field(n+3))
        
        assert len(self.GCj)==len(self.GCi),'(len(GCj)=%s len(GCi)=%s' %(len(self.GCj),len(self.GCi))
        #if isComplex():
            #self.Complex(card.field(v)
    def buildMatrices(self):
        pass
    #def renameMatrix(self,nameNew):
        #pass
    def isComplex(self):
        if self.tin in [3,4]:
            return True
        return False

    def __init__(self,card=None,data=None):
        self.name = card.field(1)
        #zero
        self.ifo   = card.field(3)
        self.tin   = card.field(4)
        self.tout  = card.field(5,0)
        
        self.polar = card.field(6)
        self.ncol  = card.field(8)

        #self.name = card.field(1)
        self.GCj = []
        self.GCi = []
        #self.Gi = []
        #self.Ci = []
        self.Real = []
        self.Complex = []

    def __repr__(self):
        fields = ['DMIG',self.name,0,self.ifo,self.tin,self.tout,self.polar,None,self.ncol]
        return fields
