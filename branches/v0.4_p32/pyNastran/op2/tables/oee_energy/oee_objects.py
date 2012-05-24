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
import sys
from math import isnan

from pyNastran.op2.resultObjects.op2_Objects import scalarObject
from pyNastran.op2.op2Errors import *

class StrainEnergyObject(scalarObject):
    """
                               E L E M E N T   S T R A I N   E N E R G I E S
 
    ELEMENT-TYPE = QUAD4      * TOTAL ENERGY OF ALL ELEMENTS IN PROBLEM     =   9.817708E+08
    SUBCASE               1   * TOTAL ENERGY OF ALL ELEMENTS IN SET       1 =   4.192036E+08
 
       ELEMENT-ID   STRAIN-ENERGY  PERCENT OF TOTAL  STRAIN-ENERGY-DENSITY
               12   2.291087E+07        2.3336            2.291087E+02
               13   1.582968E+07        1.6124            1.055312E+02
               14   6.576075E+07        6.6982            3.288037E+02
    """
    def __init__(self,dataCode,iSubcase,dt=None):
        scalarObject.__init__(self,dataCode,iSubcase)
        self.energy  = {}
        self.percent = {}
        self.density = {}
        #print self.dataCode
        #print "numWide = %s %s"  %(self.dataCode['numWide'],type(self.dataCode['numWide']))
        self.updateNumWide()

        if dt is not None:
            self.dt = dt
            #print "*****dt = ",self.dt
            self.addNewTransient()
            self.isTransient = True
            self.add = self.addTransient
            #asdf
        ###

    def updateNumWide(self):
        #print "***dataCode = ",self.dataCode
        if self.dataCode['numWide']==4:
            self.getLength    = self.getLength4
            self.add          = self.add4
            self.addTransient = self.addTransient4
        elif self.dataCode['numWide']==5:
            self.getLength    = self.getLength5
            self.add          = self.add5
            self.addTransient = self.addTransient5
        else:
            raise RuntimeError('invalid numWide=%s' %(self.numWide))
        ###
        #print "self.dt = ",self.dt
        if self.dt is not None:
            self.add = self.addTransient

    def getLength4(self):
        return(16,'ifff')

    def getLength5(self):
        return(20,'ssssssssfff')

    def updateDt(self,dataCode,dt):
        """
        this method is called if the object
        already exits and a new time step is found
        """
        self.dataCode = dataCode
        self.applyDataCode()
        #assert dt>=0.
        self.log.debug("updating %s...%s=%s  iSubcase=%s" %(self.name,self.name,dt,self.iSubcase))
        #print "dataCode = ",self.dataCode
        if dt is not None:
            self.dt = dt
            self.addNewTransient()
        ###
        self.updateNumWide()

    def deleteTransient(self,dt):
        del self.energy[dt]
        del self.percent[dt]
        del self.density[dt]

    def getTransients(self):
        k = self.energy.keys()
        k.sort()
        return k

    def addNewTransient(self):
        """
        initializes the transient variables
        """
        if self.dt not in self.energy:
            self.energy[self.dt]  = {}
            self.percent[self.dt] = {}
            self.density[self.dt] = {}
        
    def add4(self,out):
        #print "add4"
        (eid,energy,percent,density) = out
        eid = (eid-self.deviceCode)//10
        #print "energyGridIDs = %s" %(self.energy.keys())
        #assert grid not in self.energy,'grid=%s out=%s' %(grid,out)
        if eid<=0:
            raise InvalidGridID_Error('eid=%s' %(eid))
        self.energy[eid]  = energy
        self.percent[eid] = percent
        self.density[eid] = density

    def add5(self,out):
        #print "add5"
        (a,b,c,d,e,f,g,h,energy,percent,density) = out
        #print "out = ",out
        word = a+b+c+d+e+f+g+h
        word = word.strip()
        #assert word not in self.energy,'%s in energy...' %(word)
        #if grid<=0:
        #    raise InvalidGridID_Error('grid=%s' %(grid))
        self.energy[word]  = energy
        self.percent[word] = percent
        self.density[word] = density

    def addTransient4(self,out):
        #print "addTransient4"
        dt = self.dt
        (eid,energy,percent,density) = out
        eid = (eid-self.deviceCode)//10
        #print str(self)
        #assert grid not in self.energy[dt],'grid=%s dt=%s energy=%s percent=%s density=%s' %(grid,dt,energy,percent,density)
        if eid<=0:
            raise InvalidGridID_Error('eid=%s' %(eid))

        self.energy[dt][eid]  = energy
        #print "self.energy = ",self.energy
        self.percent[dt][eid] = percent
        self.density[dt][eid] = density
    
    def addTransient5(self,out):
        #print "addTransient5"
        dt = self.dt
        (a,b,c,d,energy,percent,density) = out
        word = a+b+c+d
        assert word not in self.energy[dt]
        #if grid<=0:
        #    raise InvalidGridID_Error('grid=%s' %(grid))
        self.energy[dt][word]  = energy
        self.percent[dt][word] = percent
        self.density[dt][word] = density

    def __reprTransient__(self):
        msg  = '---Transient Strain Energy Object---\n'
        for dt,Energy in sorted(self.energy.items()):
            msg += "%s = %s\n" %(self.dataCode['name'],dt)
            msg += "%-10s %-14s% -14s% -14s\n" %('EID','Energy','PercentTotal','Density')
            #print "dt=%s Energy=%s" %(dt,Energy)
            for eid,energy in sorted(Energy.items()):
                percent = self.percent[dt][eid]
                density = self.density[dt][eid]
                if isnan(density):
                    density2 = '%-14s\n' %('-----')
                else:
                    density2 = '%-14g\n' %(density)
                
                #print "eid = ",eid
                #print "energy = ",energy
                #print "percent = ",percent
                #print "density = ",density
                #print "density2 = ",density2
                #sys.stdout.flush()
                msg += "%-10s %-14g %-14g %s" %(eid,energy,percent,density2)
            ###
        ###
        return msg

    def __repr__(self):
        if self.isTransient:
            return self.__reprTransient__()

        msg  = '---Strain Energy Object---\n'
        msg += "%-10s %-14s% -14s% -14s\n" %('EID','Energy','PercentTotal','Density')
        for eid,energy in sorted(self.energy.items()):
            percent = self.percent[eid]
            density = self.density[eid]
            if isnan(density):
                density2 = '%-14s\n' %('-----')
            else:
                density2 = '%-14g\n' %(density)
            msg += "%-10s %-14g %-14g %s" %(eid,energy,percent,density2)
        return msg
