from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
from numpy import argsort

from pyNastran.op2.resultObjects.op2_Objects import ScalarObject


class OES_Object(ScalarObject):
    def __init__(self, data_code, isubcase, apply_data_code=True):
        self.element_type = None
        self.element_name = None
        self.nonlinear_factor = None
        self._times = None
        ScalarObject.__init__(self, data_code, isubcase, apply_data_code=apply_data_code)
        #self.log.debug("starting OES...element_name=%-6s isubcase=%s" % (self.element_name, self.isubcase))
        #print self.data_code

    def isCurvatureOld(self):
        if self.stress_bits[2] == 0:
            return True
        return False

    def isCurvature(self):
        if self.s_code in [0, 1, 14, 15, 16, 17, 27, 30, 31]:  # fiber distance
            return False
        elif self.s_code in [10, 11, 26, ]:  # fiber curvature
            return True
        raise NotImplementedError('add s_code=%s' % self.s_code)

    def isFiberDistance(self):
        return not(self.isCurvature())

    def isVonMises(self):
        #print self.stress_bits
        #iMs = not(self.isMaxShear())
        #print 'isVonMises = ',iMs
        return not(self.isMaxShear())

    def isMaxShear(self):
        #print self.stress_bits
        if self.stress_bits[4] == 0:
            #print 'isMaxShear = True'
            return True
        #print 'isMaxShear = False'
        return False

    def getOrderedETypes(self, valid_types):
        """
        Groups element IDs by element type

        :param valid_types: list of valid element types
                           e.g. ['CTRIA3', 'CTRIA6', 'CQUAD4', 'CQUAD8']
        :returns TypesOut:      the ordered list of types
        :returns orderedETypes: dictionary of Type-IDs to write
        """
        orderedETypes = {}

        #valid_types = ['CTRIA3', 'CTRIA6', 'CQUAD4', 'CQUAD8']
        for eType in valid_types:
            orderedETypes[eType] = []
        for eid, eType in sorted(iteritems(self.eType)):
            assert eType in valid_types, 'unsupported eType=%r; valid_type=%s' % (eType, str(['%r' % str(t) for t in valid_types]))
            orderedETypes[eType].append(eid)

        minVals = []
        for eType in valid_types:
            vals = orderedETypes[eType]
            if len(vals) == 0:
                minVals.append(-1)
            else:
                minVals.append(min(vals))

        #print "minVals = ",minVals
        argList = argsort(minVals)

        TypesOut = []
        for i in argList:
            TypesOut.append(valid_types[i])
        return (TypesOut, orderedETypes)


class StressObject(OES_Object):
    def __init__(self, data_code, isubcase):
        OES_Object.__init__(self, data_code, isubcase)

    def update_dt(self, data_code, dt):
        self.data_code = data_code
        self.apply_data_code()
        #assert dt>=0.
        #print "data_code=",self.data_code
        self.element_name = self.data_code['element_name']
        if dt is not None:
            #print("updating stress...%s=%s element_name=%s" %
            #     (self.data_code['name'], dt, self.element_name))
            self.dt = dt
            self.add_new_transient(dt)

    def isStrain(self):
        return self.is_stress()

    def isStress(self):
        return self.is_strain()

    def is_strain(self):
        return True

    def is_stress(self):
        return False


class StrainObject(OES_Object):
    def __init__(self, data_code, isubcase):
        OES_Object.__init__(self, data_code, isubcase)

    def update_dt(self, data_code, dt):
        self.data_code = data_code
        self.apply_data_code()
        self.element_name = self.data_code['element_name']
        if dt is not None:
            #print("updating strain...%s=%s element_name=%s" %
            #     (self.data_code['name'], dt, self.element_name))
            self.dt = dt
            self.add_new_transient(dt)

    def isStrain(self):
        return self.is_stress()

    def isStress(self):
        return self.is_strain()

    def is_strain(self):
        return False

    def is_stress(self):
        return True
