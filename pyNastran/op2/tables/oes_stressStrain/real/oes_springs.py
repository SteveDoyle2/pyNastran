from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from .oes_objects import StressObject, StrainObject


class CelasStressObject(StressObject):
    """
    ::

                                S T R E S S E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )
          TIME         STRESS              TIME         STRESS              TIME         STRESS              TIME         STRESS
      0.0            0.0               5.000000E-02   0.0               1.000000E-01   0.0               1.500000E-01   0.0
      2.000000E-01   0.0               2.500000E-01   0.0               3.000000E-01   0.0               3.500000E-01   0.0
    """
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        StressObject.__init__(self, data_code, isubcase, read_mode)
        self.eType = {}
        self.element_name = self.data_code['element_name']

        self.code = [self.format_code, self.sort_code, self.s_code]
        self.stress = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                #self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            #self.add = self.addSort2
            self.add_new_eid = self.add_new_eid_sort2

    def get_stats(self):
        nelements = len(self.eType)

        msg = self.get_data_code()
        eTypes = list(set(self.eType.values()))
        if self.dt is not None:  # transient
            ntimes = len(self.stress)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, stress\n')
        msg.append('  eTypes = %s\n' %(', '.join(eTypes)))
        return msg

    def getLength(self):
        return (8, 'f')

    def get_transients(self):
        k = self.stress.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """initializes the transient variables"""
        self.element_name = self.data_code['element_name']
        self.dt = dt
        self.stress[dt] = {}

    def add_new_eid(self, dt, eid, out):
        (stress,) = out
        self.eType[eid] = self.element_name
        self.stress[eid] = stress

    def add_new_eid_sort1(self, dt, eid, out):
        if dt not in self.stress:
            self.add_new_transient(dt)
        (stress,) = out
        self.eType[eid] = self.element_name
        self.stress[dt][eid] = stress

    def add_new_eid_sort2(self, eid, dt, out):
        if dt not in self.stress:
            self.add_new_transient(dt)
        (stress,) = out
        self.eType[eid] = self.element_name
        self.stress[dt][eid] = stress

    def __repr__(self):
        return self.get_stats()


class CelasStrainObject(StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        StrainObject.__init__(self, data_code, isubcase, read_mode)
        self.eType = {}
        self.element_name = self.data_code['element_name']

        self.code = [self.format_code, self.sort_code, self.s_code]

        self.isTransient = False
        self.strain = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                #self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            #self.add = self.addSort2
            self.add_new_eid = self.add_new_eid_sort2

    def get_stats(self):
        nelements = len(self.eType)
        eTypes = list(set(self.eType.values()))

        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.stress)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, strain\n')
        msg.append('  eTypes = %s\n' %(', '.join(eTypes)))
        return msg

    def getLength(self):
        return (8, 'f')

    def get_transients(self):
        k = self.strain.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        self.strain[dt] = {}

    def add_new_eid(self, dt, eid, out):
        (strain,) = out
        assert eid >= 0
        self.eType[eid] = self.element_name
        self.strain[eid] = strain

    def add_new_eid_sort1(self, dt, eid, out):
        (strain,) = out
        assert eid >= 0

        self.eType[eid] = self.element_type
        self.strain[dt][eid] = strain

    def __repr__(self):
        return self.get_stats()


class NonlinearSpringStressObject(StressObject):
    """
    """
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        StressObject.__init__(self, data_code, isubcase, read_mode)
        self.eType = {}
        self.element_name = self.data_code['element_name']

        self.code = [self.format_code, self.sort_code, self.s_code]
        self.force = {}
        self.stress = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            self.add_new_eid = self.add_new_eid_sort2

    def get_stats(self):
        nelements = len(self.eType)
        eTypes = list(set(self.eType.values()))

        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.stress)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, force, stress\n')
        msg.append('  eTypes = %s\n' %(', '.join(eTypes)))
        return msg

    def get_transients(self):
        k = self.stress.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """initializes the transient variables"""
        self.dt = dt
        self.force[dt] = {}
        self.stress[dt] = {}

    def add_new_eid(self, eType, dt, eid, force, stress):
        self.eType[eid] = eType
        self.force[eid] = force
        self.stress[eid] = stress

    def add_new_eid_sort1(self, eType, dt, eid, force, stress):
        if dt not in self.stress:
            self.add_new_transient(dt)
        self.eType[eid] = eType
        self.force[dt][eid] = force
        self.stress[dt][eid] = stress

    def add_new_eid_sort2(self, eType, eid, dt, force, stress):
        if dt not in self.stress:
            self.add_new_transient(dt)
        self.eType[eid] = eType
        self.force[dt][eid] = force
        self.stress[dt][eid] = stress

    def __repr__(self):
        return self.get_stats()