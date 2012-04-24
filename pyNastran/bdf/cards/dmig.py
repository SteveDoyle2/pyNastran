## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
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
    def __init__(self,card=None,data=None):
        self.name = card.field(1)
        #zero
        self.ifo   = card.field(3)
        self.tin   = card.field(4)
        self.tout  = card.field(5,0)
        self.polar = card.field(6)
        self.ncol  = card.field(8)

        #self.name = card.field(1)
        #self.Gj = []
        #self.Cj = []
        #self.Gi = []
        #self.Ci = []
        #self.Ai = []
        #self.Bi = []
        
