from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import numpy as np
from pyNastran.op2.result_objects.op2_objects import ScalarObject

SORT2_TABLE_NAME_MAP = {
    'OES2' : 'OES1',
    'OES2C' : 'OES1C',
    'OESATO2' : 'OESATO1',
    'OESCRM2' : 'OESCRM1',
    'OESNO2' : 'OESNO1',
    'OESPSD2' : 'OESPSD1',
    'OESRMS2' : 'OESRMS1',

    'OSTR2' : 'OSTR1',
    'OSTR2C' : 'OSTR1C',
    'OSTRATO2' : 'OSTRATO1',
    'OSTRCRM2' : 'OSTRCRM1',
    'OSTRNO2' : 'OSTRNO1',
    'OSTRPSD2' : 'OSTRPSD1',
    'OSTRRMS2' : 'OSTRRMS1',
}

class OES_Object(ScalarObject):
    def __init__(self, data_code, isubcase, apply_data_code=True):
        self.element_type = None
        self.element_name = None
        self.nonlinear_factor = np.nan
        self._times = None
        ScalarObject.__init__(self, data_code, isubcase, apply_data_code=apply_data_code)
        #self.log.debug("starting OES...element_name=%-6s isubcase=%s" % (self.element_name, self.isubcase))
        #print self.data_code

    def finalize(self):
        """it's required that the object be in SORT1"""
        self.set_as_sort1()

    def set_as_sort1(self):
        """the data is in SORT1, but the flags are wrong"""
        if self.is_sort1:
            return
        self.table_name = SORT2_TABLE_NAME_MAP[self.table_name]
        self.sort_bits[1] = 0 # sort1
        self.sort_method = 1
        assert self.is_sort1 is True, self.is_sort1

    @property
    def is_curvature(self):
        if self.is_stress:
            curvature_flag = False
        else:
            # strain only
            curvature_flag = True if self.stress_bits[2] == 0 else False
        if self.s_code in [10, 11, 20, 27]:
            assert curvature_flag, curvature_flag
            return True
        assert not curvature_flag, curvature_flag
        return False

    @property
    def is_fiber_distance(self):
        return not self.is_curvature

    @property
    def is_von_mises(self):
        return not self.is_max_shear

    @property
    def is_max_shear(self):
        return True if self.stress_bits[4] == 0 else False

    def _get_headers(self):
        raise NotImplementedError('overwrite this')

    @property
    def is_stress(self):
        raise NotImplementedError('overwrite this')

class StressObject(OES_Object):
    def __init__(self, data_code, isubcase):
        OES_Object.__init__(self, data_code, isubcase)
        assert self.is_stress, self.code_information()
        assert not self.is_strain, self.code_information()

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

    @property
    def is_strain(self):
        class_name = self.__class__.__name__
        assert self.stress_bits[1] == self.stress_bits[3], 'class_name=%s scode=%s stress_bits=%s' % (class_name, self.s_code, self.stress_bits)
        assert self.stress_bits[1] == 0, 'class_name=%s scode=%s stress_bits=%s' % (class_name, self.s_code, self.stress_bits)
        return False

    @property
    def is_stress(self):
        class_name = self.__class__.__name__
        assert self.stress_bits[1] == self.stress_bits[3], 'class_name=%s scode=%s stress_bits=%s' % (class_name, self.s_code, self.stress_bits)
        assert self.stress_bits[1] == 0, 'class_name=%s scode=%s stress_bits=%s' % (class_name, self.s_code, self.stress_bits)
        return True


class StrainObject(OES_Object):
    def __init__(self, data_code, isubcase):
        OES_Object.__init__(self, data_code, isubcase)
        assert self.is_strain, self.code_information()
        assert not self.is_stress, self.code_information()

    def update_dt(self, data_code, dt):
        self.data_code = data_code
        self.apply_data_code()
        self.element_name = self.data_code['element_name']
        if dt is not None:
            #print("updating strain...%s=%s element_name=%s" %
            #     (self.data_code['name'], dt, self.element_name))
            self.dt = dt
            self.add_new_transient(dt)

    @property
    def is_strain(self):
        class_name = self.__class__.__name__
        assert self.stress_bits[1] == self.stress_bits[3], 'class_name=%s scode=%s stress_bits=%s; table_name=%r' % (class_name, self.s_code, self.stress_bits, self.table_name)
        assert self.stress_bits[1] == 1, 'class_name=%s scode=%s stress_bits=%s; table_name=%r' % (class_name, self.s_code, self.stress_bits, self.table_name)
        return True

    @property
    def is_stress(self):
        class_name = self.__class__.__name__
        assert self.stress_bits[1] == self.stress_bits[3], 'class_name=%s scode=%s stress_bits=%s; table_name=%r' % (class_name, self.s_code, self.stress_bits, self.table_name)
        assert self.stress_bits[1] == 1, 'class_name=%s is_stress=False scode=%s stress_bits=%s; element_type=%s element_name=%s; table_name=%r' % (class_name, self.s_code, self.stress_bits, self.element_type, self.element_name, self.table_name)
        return False
