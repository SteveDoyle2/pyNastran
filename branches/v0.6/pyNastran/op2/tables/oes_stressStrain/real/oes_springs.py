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
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        StressObject.__init__(self, data_code, isubcase)
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

    def delete_transient(self, dt):
        del self.stress[dt]

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

    def add_f06_data(self, data, dt):
        if dt is not None:
            for datai in data:
                (eid, stressi) = datai
                self.stress[dt][eid] = stressi
            return

        for datai in data:
            (eid, stressi) = datai
            self.stress[eid] = stressi

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f)
        msg = header + ['                              S T R E S S E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )\n',
                        '      ELEMENT         STRESS           ELEMENT         STRESS           ELEMENT         STRESS           ELEMENT         STRESS\n',
                        '        ID.                              ID.                              ID.                              ID.\n',
                        ]
        stresses = []
        #elements = []
        line = '   '
        for eid, stress in sorted(self.stress.iteritems()):
            #elements.append(eid)
            stresses.append(stress)
            line += '%10s  %10.6E     ' % (eid, stress)
            if len(stresses) == 3:
                stresses = []
                msg.append(line.rstrip() + '\n')

        if stresses:
            msg.append(line.rstrip() + '\n')
        msg.append(pageStamp % pageNum + '\n')

        f.write(''.join(msg))
        return pageNum


class CelasStrainObject(StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        StrainObject.__init__(self, data_code, isubcase)
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

    def delete_transient(self, dt):
        del self.strain[dt]

    def get_transients(self):
        k = self.strain.keys()
        k.sort()
        return k

    def add_f06_data(self, data, dt):
        if dt is not None:
            for datai in data:
                (eid, straini) = datai
                self.strain[dt][eid] = straini
            return

        for datai in data:
            (eid, straini) = datai
            self.strain[eid] = straini

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f)
        msg = header + ['                               S T R A I N S    I N   S C A L A R   S P R I N G S        ( C E L A S 2 )\n',
                        '      ELEMENT         STRAIN           ELEMENT         STRAIN           ELEMENT         STRAIN           ELEMENT         STRAIN\n',
                        '        ID.                              ID.                              ID.                              ID.\n',
                        ]
        strains = []
        #elements = []
        line = '   '
        for eid, strain in sorted(self.strain.iteritems()):
            #elements.append(eid)
            strains.append(strain)
            line += '%10s  %10.6E     ' % (eid, strain)
            if len(strains) == 3:
                strains = []
                msg.append(line.rstrip() + '\n')

        if strains:
            msg.append(line.rstrip() + '\n')
        msg.append(pageStamp % pageNum + '\n')

        f.write(''.join(msg))
        return pageNum


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


class NonlinearSpringStressObject(StressObject):
    """
    """
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        StressObject.__init__(self, data_code, isubcase)
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

    def delete_transient(self, dt):
        del self.force[dt]
        del self.stress[dt]

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