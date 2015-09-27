from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
from six.moves import zip

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject
from pyNastran.f06.f06_formatting import writeFloats13E


class RealCelasStress(StressObject):
    """
    ::

                                S T R E S S E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )
          TIME         STRESS              TIME         STRESS              TIME         STRESS              TIME         STRESS
      0.0            0.0               5.000000E-02   0.0               1.000000E-01   0.0               1.500000E-01   0.0
      2.000000E-01   0.0               2.500000E-01   0.0               3.000000E-01   0.0               3.500000E-01   0.0
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
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

    def is_real(self):
        return True

    def is_complex(self):
        return False

    def get_stats(self):
        nelements = len(self.eType)

        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.stress)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, stress\n')
        msg.append('  %s\n' % self.element_name)
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

    def add_f06_data(self, data, transient):
        if transient is not None:
            dt = transient[1]
            assert not isinstance(dt, list)
            if dt not in self.stress:
                self.stress[dt] = {}
            for datai in data:
                (eid, stressi) = datai
                self.stress[dt][eid] = stressi
            return

        for datai in data:
            (eid, stressi) = datai
            self.stress[eid] = stressi

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        msg = header + ['                              S T R E S S E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )\n',
                        '      ELEMENT         STRESS           ELEMENT         STRESS           ELEMENT         STRESS           ELEMENT         STRESS\n',
                        '        ID.                              ID.                              ID.                              ID.\n',
                        ]
        f.write(''.join(msg))

        _write_f06_springs(f, self.stress)
        #_write_f06_springs_transient(f, self.stress)
        f.write(page_stamp % page_num)
        return page_num

def _write_f06_springs_transient(f, stress, header, words, name):
    for dt, datai in sorted(iteritems(data)):
        header[1] = ' %s = %10.4E\n' % (name, dt)
        msg += header + words
        f.write(''.join(msg))

        eids = []
        stresses = []
        for eid, stress in sorted(iteritems(datai)):
            eids.append(eid)
            stresses.append(stress)
            if len(stresses) == 4:
                stresses, is_all_zeros = writeFloats13E(stresses)
                f.write('    %10i  %13s    %10i  %13s    %10i  %13s    %10i  %13s\n' % (
                    eids[0], stresses[0],
                    eids[1], stresses[1],
                    eids[2], stresses[2],
                    eids[3], stresses[3]))
                eids = []
                stresses = []

        if stresses:
            line = '    '
            stresses, is_all_zeros = writeFloats13E(stresses)
            for eid, stress in zip(eids, stresses):
                line += '%10i  %13s    ' % (eid, stress)
            f.write(line.rstrip() + '\n')

            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            msg = ['']
        page_num += 1

    return page_num - 1


def _write_f06_springs(f, data):
    eids = []
    stresses = []
    for eid, stress in sorted(iteritems(data)):
        eids.append(eid)
        stresses.append(stress)
        if len(stresses) == 4:
            stresses, is_all_zeros = writeFloats13E(stresses)
            f.write('    %10i  %13s    %10i  %13s    %10i  %13s    %10i  %13s\n' % (
                eids[0], stresses[0],
                eids[1], stresses[1],
                eids[2], stresses[2],
                eids[3], stresses[3]))
            eids = []
            stresses = []

    if stresses:
        line = '    '
        stresses, is_all_zeros = writeFloats13E(stresses)
        for eid, stress in zip(eids, stresses):
            line += '%10i  %13s    ' % (eid, stress)
        f.write(line.rstrip() + '\n')


class RealCelasStrain(StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
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
            #self.add_new_eid = self.add_new_eid_sort2

    def get_stats(self):
        nelements = len(self.eType)
        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.stress)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, strain\n')
        msg.append('  %s\n' % self.element_name)
        return msg

    def getLength(self):
        return (8, 'f')

    def delete_transient(self, dt):
        del self.strain[dt]

    def get_transients(self):
        k = self.strain.keys()
        k.sort()
        return k

    def add_f06_data(self, data, transient):
        if transient is not None:
            dt = transient[1]
            if dt not in self.strain:
                self.strain[dt] = {}
            for datai in data:
                (eid, straini) = datai
                self.strain[dt][eid] = straini
            return

        for datai in data:
            (eid, straini) = datai
            self.strain[eid] = straini

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        msg = header + ['                               S T R A I N S    I N   S C A L A R   S P R I N G S        ( C E L A S 2 )\n',
                        '      ELEMENT         STRAIN           ELEMENT         STRAIN           ELEMENT         STRAIN           ELEMENT         STRAIN\n',
                        '        ID.                              ID.                              ID.                              ID.\n',
                        ]
        f.write(''.join(msg))
        _write_f06_springs(f, self.strain)
        f.write(page_stamp % page_num)
        return page_num


    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        self.strain[dt] = {}

    def add_new_eid(self, dt, eid, out):
        (strain,) = out
        assert eid >= 0, eid
        self.eType[eid] = self.element_name
        self.strain[eid] = strain

    def add_new_eid_sort1(self, dt, eid, out):
        (strain,) = out
        assert eid >= 0, eid

        self.eType[eid] = self.element_type
        self.strain[dt][eid] = strain


class NonlinearSpringStress(StressObject):
    """
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
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
        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.stress)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, force, stress\n')
        msg.append('  %s\n' % self.element_name)
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
