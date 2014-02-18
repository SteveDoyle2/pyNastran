from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from math import isnan

from pyNastran.op2.resultObjects.op2_Objects import ScalarObject


class StrainEnergyObject(ScalarObject):
    """
    ::

                                 E L E M E N T   S T R A I N   E N E R G I E S

      ELEMENT-TYPE = QUAD4      * TOTAL ENERGY OF ALL ELEMENTS IN PROBLEM     =   9.817708E+08
      SUBCASE               1   * TOTAL ENERGY OF ALL ELEMENTS IN SET       1 =   4.192036E+08

         ELEMENT-ID   STRAIN-ENERGY  PERCENT OF TOTAL  STRAIN-ENERGY-DENSITY
                 12   2.291087E+07        2.3336            2.291087E+02
                 13   1.582968E+07        1.6124            1.055312E+02
                 14   6.576075E+07        6.6982            3.288037E+02
    """
    def __init__(self, data_code, is_sort1, isubcase, dt=None, read_mode=0):
        ScalarObject.__init__(self, data_code, isubcase, read_mode)

        self.energy = {}
        self.percent = {}
        self.density = {}
        #if read_mode == 0:
            #return
        #print self.data_code
        #print "num_wide = %s %s"  %(self.data_code['num_wide'],type(self.data_code['num_wide']))

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.addSort2

    def get_stats(self):
        msg = self._get_data_code()
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.energy)
            time0 = self.energy.keys()[0]
            nelements = len(self.energy[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.energy)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  energy, percent, density\n')
        return msg

    def update_dt(self, data_code, dt):
        """
        this method is called if the object
        already exits and a new time step is found
        """
        self.data_code = data_code
        self.apply_data_code()
        #assert dt>=0.
        self.log.debug("updating %s...%s=%s  isubcase=%s" % (
            self.name, self.name, dt, self.isubcase))
        #print "data_code = ",self.data_code
        if dt is not None:
            self.dt = dt
            self.add_new_transient(dt)
        self.updateNumWide()

    def delete_transient(self, dt):
        del self.energy[dt]
        del self.percent[dt]
        del self.density[dt]

    def get_transients(self):
        k = self.energy.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
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

    def add_sort1(self, dt, out):
        if dt not in self.energy:
            self.add_new_transient(dt)

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
            msg += "%s = %s\n" % (self.data_code['name'], dt)
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
        return msg
