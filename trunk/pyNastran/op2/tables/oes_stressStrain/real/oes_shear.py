from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from .oes_objects import StressObject, StrainObject


class RealShearStressObject(StressObject):
    """
    ::

      # format_code=1 sort_code=0 stressCode=0
                                     S T R E S S E S   I N   S H E A R   P A N E L S      ( C S H E A R )
      ELEMENT            MAX            AVG        SAFETY         ELEMENT            MAX            AVG        SAFETY
        ID.             SHEAR          SHEAR       MARGIN           ID.             SHEAR          SHEAR       MARGIN
          328        1.721350E+03   1.570314E+03   7.2E+01
    """
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        StressObject.__init__(self, data_code, isubcase)
        self.eType = 'CSHEAR'

        self.code = [self.format_code, self.sort_code, self.s_code]
        self.maxShear = {}
        self.avgShear = {}
        self.margin = {}

        self.getLength = self.getLength
        self.isImaginary = False
        #if dt is not None:
        #    self.add_new_transient = self.add_new_transient
        #    self.add_new_eid       = self.addNewEidTransient
        #else:
        #    self.add_new_eid = self.add_new_eid

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            self.add = self.addSort2
            self.add_new_eid = self.add_new_eid_sort2

    def get_stats(self):
        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.maxShear)
            s0 = self.maxShear.keys()[0]
            nelements = len(self.maxShear[s0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.maxShear)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, maxShear, avgShear, margin\n')
        return msg

    def getLength(self):
        return (16, 'fff')

    def delete_transient(self, dt):
        del self.maxShear[dt]
        del self.avgShear[dt]
        del self.margin[dt]

    def get_transients(self):
        k = self.maxShear.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        self.dt = dt
        self.maxShear[dt] = {}
        self.avgShear[dt] = {}
        self.margin[dt] = {}

    def add_f06_data(self, data, dt):
        if dt:
            if dt not in self.maxShear:
                self.maxShear[dt] = {}
                self.avgShear[dt] = {}
                self.margin[dt] = {}
            for datai in data:
                (eid, max_shear, avg_shear, margin) = datai
                self.maxShear[dt][eid] = max_shear
                self.avgShear[dt][eid] = avg_shear
                self.margin[dt][eid] = margin
                return

        for datai in data:
            (eid, max_shear, avg_shear, margin) = datai
            self.maxShear[eid] = max_shear
            self.avgShear[eid] = avg_shear
            self.margin[eid] = margin

    def add_new_eid(self, dt, eid, out):
        #print "Rod Stress add..."
        (maxShear, avgShear, margin) = out
        assert isinstance(eid, int)
        self.maxShear = {}
        self.avgShear = {}
        self.margin = {}
        self.maxShear[eid] = maxShear
        self.avgShear[eid] = avgShear
        self.margin[eid] = margin

    def add_new_eid_sort1(self, dt, eid, out):
        (maxShear, avgShear, margin) = out
        if dt not in self.maxShear:
            self.add_new_transient(dt)
        assert isinstance(eid, int)
        assert eid >= 0
        self.maxShear[dt][eid] = maxShear
        self.avgShear[dt][eid] = avgShear
        self.margin[dt][eid] = margin


class RealShearStrainObject(StrainObject):

    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        StrainObject.__init__(self, data_code, isubcase)
        self.eType = 'CSHEAR'
        #raise Exception('not supported...CSHEAR strain')
        self.code = [self.format_code, self.sort_code, self.s_code]
        self.maxShear = {}
        self.avgShear = {}
        self.margin = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            self.add = self.addSort2
            self.add_new_eid = self.add_new_eid_sort2

    def get_stats(self):
        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.maxShear)
            s0 = self.maxShear.keys()[0]
            nelements = len(self.maxShear[s0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.maxShear)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, maxShear, avgShear, margin\n')
        return msg

    def getLength(self):
        return (16, 'fff')

    def delete_transient(self, dt):
        del self.maxShear[dt]
        del self.avgShear[dt]
        del self.margin[dt]

    def get_transients(self):
        k = self.maxShear.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        .. note:: make sure you set self.dt first
        """
        self.dt = dt
        self.maxShear[dt] = {}
        self.avgShear[dt] = {}
        self.margin[dt] = {}

    def add_f06_data(self, data, dt):
        if dt:
            if dt not in self.maxShear:
                self.maxShear[dt] = {}
                self.avgShear[dt] = {}
                self.margin[dt] = {}
            for datai in data:
                (eid, max_shear, avg_shear, margin) = datai
                self.maxShear[dt][eid] = max_shear
                self.avgShear[dt][eid] = avg_shear
                self.margin[dt][eid] = margin
                return

        for datai in data:
            (eid, max_shear, avg_shear, margin) = datai
            self.maxShear[eid] = max_shear
            self.avgShear[eid] = avg_shear
            self.margin[eid] = margin

    def add_new_eid(self, dt, eid, out):
        raise NotImplementedError()
        (axial, SMa, torsion, SMt) = out
        #print "Rod Strain add..."
        assert eid >= 0
        #self.eType = self.eType
        self.maxShearl[eid] = axial
        self.avgShear[eid] = SMa
        self.margin[eid] = torsion

    def add_new_eid_sort1(self, dt, eid, out):
        (maxShear, avgShear, margin) = out
        if dt not in self.maxShear:
            self.add_new_transient(dt)
        assert eid >= 0

        #self.eType[eid] = self.element_type
        self.maxShear[dt][eid] = maxShear
        self.avgShear[dt][eid] = avgShear
        self.margin[dt][eid] = margin

    def add_new_eid_sort2(self, eid, dt, out):
        (maxShear, avgShear, margin) = out
        if dt not in self.maxShear:
            self.add_new_transient(dt)
        assert eid >= 0

        #self.eType[eid] = self.element_type
        self.maxShear[dt][eid] = maxShear
        self.avgShear[dt][eid] = avgShear
        self.margin[dt][eid] = margin