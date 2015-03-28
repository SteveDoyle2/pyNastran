from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
from numpy import argsort

from pyNastran.op2.resultObjects.op2_Objects import ScalarObject

class OES_Object_Deprecated(object):
    def __init__(self):
        pass

    def isVonMises(self):
        return self.is_von_mises()

    def isMaxShear(self):
        return self.is_max_shear()

    def isFiberDistance(self):
        return self.is_fiber_distance()

    def isCurvatureOld(self):
        if self.stress_bits[2] == 0:
            return True
        return False

    def isCurvature(self):
        return self.is_curvature()

class OES_Object(OES_Object_Deprecated, ScalarObject):
    def __init__(self, data_code, isubcase, apply_data_code=True):
        self.element_type = None
        self.element_name = None
        self.nonlinear_factor = None
        self._times = None
        OES_Object_Deprecated.__init__(self)
        ScalarObject.__init__(self, data_code, isubcase, apply_data_code=apply_data_code)
        #self.log.debug("starting OES...element_name=%-6s isubcase=%s" % (self.element_name, self.isubcase))
        #print self.data_code

    def is_curvature(self):
        if self.s_code in [0, 1, 14, 15, 16, 17, 27, 30, 31]:  # fiber distance
            return False
        elif self.s_code in [10, 11, 26, ]:  # fiber curvature
            return True
        raise NotImplementedError('add s_code=%s' % self.s_code)

    def is_fiber_distance(self):
        return not self.is_curvature()

    def is_von_mises(self):
        #print self.stress_bits
        #iMs = not(self.is_max_shear())
        #print 'is_von_mises = ', iMs
        return not self.is_max_shear()

    def is_max_shear(self):
        #print self.stress_bits
        if self.stress_bits[4] == 0:
            #print 'is_max_shear = True'
            return True
        #print 'is_max_shear = False'
        return False

    def getOrderedETypes(self, valid_types):
        """
        Groups element IDs by element type

        :param valid_types: list of valid element types
                           e.g. ['CTRIA3', 'CTRIA6', 'CQUAD4', 'CQUAD8']
        :returns types_out:      the ordered list of types
        :returns ordered_etypes: dictionary of Type-IDs to write
        """
        ordered_etypes = {}

        #valid_types = ['CTRIA3', 'CTRIA6', 'CQUAD4', 'CQUAD8']
        for etype in valid_types:
            ordered_etypes[etype] = []
        for eid, etype in sorted(iteritems(self.eType)):
            assert etype in valid_types, 'unsupported eType=%r; valid_type=%s' % (etype, str(['%r' % str(t) for t in valid_types]))
            ordered_etypes[etype].append(eid)

        min_vals = []
        for etype in valid_types:
            vals = ordered_etypes[etype]
            if len(vals) == 0:
                min_vals.append(-1)
            else:
                min_vals.append(min(vals))
        arg_list = argsort(min_vals)

        types_out = []
        for i in arg_list:
            types_out.append(valid_types[i])
        return (types_out, ordered_etypes)


class StressObject_Deprecated(object):
    def __init__(self):
        pass

    def isStrain(self):
        return self.is_stress()

    def isStress(self):
        return self.is_strain()


class StressObject(StressObject_Deprecated, OES_Object):
    def __init__(self, data_code, isubcase):
        StressObject_Deprecated.__init__(self)
        OES_Object.__init__(self, data_code, isubcase)

    def update_dt(self, data_code, dt):
        self.data_code = data_code
        self.apply_data_code()
        #assert dt >= 0.
        self.element_name = self.data_code['element_name']
        if dt is not None:
            #print("updating stress...%s=%s element_name=%s" %
            #     (self.data_code['name'], dt, self.element_name))
            self.dt = dt
            self.add_new_transient(dt)

    def is_strain(self):
        return True

    def is_stress(self):
        return False

class StrainObject_Deprecated(object):
    def __init__(self):
        pass

    def isStrain(self):
        return self.is_stress()

    def isStress(self):
        return self.is_strain()

class StrainObject(StrainObject_Deprecated, OES_Object):
    def __init__(self, data_code, isubcase):
        StrainObject_Deprecated.__init__(self)
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

    def is_strain(self):
        return False

    def is_stress(self):
        return True
