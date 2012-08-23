from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import sys
from math import isnan

from pyNastran.op2.resultObjects.op2_Objects import scalarObject


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
    def __init__(self, dataCode, isSort1, iSubcase, dt=None):
        scalarObject.__init__(self, dataCode, iSubcase)
        self.energy = {}
        self.percent = {}
        self.density = {}
        #print self.dataCode
        #print "numWide = %s %s"  %(self.dataCode['numWide'],type(self.dataCode['numWide']))

        self.dt = dt
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
        else:
            assert dt is not None
            self.add = self.addSort2

    def updateDt(self, dataCode, dt):
        """
        this method is called if the object
        already exits and a new time step is found
        """
        self.dataCode = dataCode
        self.applyDataCode()
        #assert dt>=0.
        self.log.debug("updating %s...%s=%s  iSubcase=%s" % (
            self.name, self.name, dt, self.iSubcase))
        #print "dataCode = ",self.dataCode
        if dt is not None:
            self.dt = dt
            self.addNewTransient(dt)
        ###
        self.updateNumWide()

    def deleteTransient(self, dt):
        del self.energy[dt]
        del self.percent[dt]
        del self.density[dt]

    def getTransients(self):
        k = self.energy.keys()
        k.sort()
        return k

    def addNewTransient(self, dt):
        """
        initializes the transient variables
        """
        self.energy[dt] = {}
        self.percent[dt] = {}
        self.density[dt] = {}

    def add(self, dt, out):
        #print "add4"
        (eid, energy, percent, density) = out
        #print "energyGridIDs = %s" %(self.energy.keys())
        #assert grid not in self.energy,'grid=%s out=%s' %(grid,out)
        if isinstance(eid, int) and eid <= 0:
            raise ValueError("Invalid Grid ID: eid=%s" % (eid))
        self.energy[eid] = energy
        self.percent[eid] = percent
        self.density[eid] = density

    def addSort1(self, dt, out):
        if dt not in self.energy:
            self.addNewTransient(dt)

        (eid, energy, percent, density) = out
        #print str(self)
        #assert grid not in self.energy[dt],'grid=%s dt=%s energy=%s percent=%s density=%s' %(grid,dt,energy,percent,density)
        if eid <= 0:
            raise ValueError("Invalid Grid ID: eid=%s" % (eid))

        self.energy[dt][eid] = energy
        #print "self.energy = ",self.energy
        self.percent[dt][eid] = percent
        self.density[dt][eid] = density

    def __reprTransient__(self):
        msg = '---Transient Strain Energy Object---\n'
        for dt, Energy in sorted(self.energy.iteritems()):
            msg += "%s = %s\n" % (self.dataCode['name'], dt)
            msg += "%-10s %-14s% -14s% -14s\n" % (
                'EID', 'Energy', 'PercentTotal', 'Density')
            #print "dt=%s Energy=%s" %(dt,Energy)
            for eid, energy in sorted(Energy.iteritems()):
                percent = self.percent[dt][eid]
                density = self.density[dt][eid]
                if isnan(density):
                    density2 = '%-14s\n' % ('-----')
                else:
                    density2 = '%-14g\n' % (density)

                #print "eid = ",eid
                #print "energy = ",energy
                #print "percent = ",percent
                #print "density = ",density
                #print "density2 = ",density2
                #sys.stdout.flush()
                msg += "%-10s %-14g %-14g %s" % (
                    eid, energy, percent, density2)
            ###
        ###
        return msg

    def __repr__(self):
        if self.nonlinearFactor is not None:
            return self.__reprTransient__()

        msg = '---Strain Energy Object---\n'
        msg += "%-10s %-14s% -14s% -14s\n" % (
            'EID', 'Energy', 'PercentTotal', 'Density')
        for eid, energy in sorted(self.energy.iteritems()):
            percent = self.percent[eid]
            density = self.density[eid]
            if isnan(density):
                density2 = '%-14s\n' % ('-----')
            else:
                density2 = '%-14g\n' % (density)
            msg += "%-10s %-14g %-14g %s" % (eid, energy, percent, density2)
        return msg
