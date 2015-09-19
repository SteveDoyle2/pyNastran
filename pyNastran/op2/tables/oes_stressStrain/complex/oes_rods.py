from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
from six.moves import range
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject
from pyNastran.f06.f06_formatting import writeImagFloats13E, get_key0

class ComplexRodDamper(StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StressObject.__init__(self, data_code, isubcase)
        self.eType = 'CBUSH'

        self.code = [self.format_code, self.sort_code, self.s_code]
        self.axial = {}
        self.torsion = {}


class ComplexRodStress(StressObject):

    def __init__(self, data_code, is_sort1, isubcase, dt):
        StressObject.__init__(self, data_code, isubcase)
        self.eType = 'CROD'

        self.code = [self.format_code, self.sort_code, self.s_code]
        self.axial = {}
        self.torsion = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                #self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            #self.add = self.addSort2
            self.add_new_eid = self.add_new_eid_sort2

    def getLength(self):
        raise RuntimeError('is this used???')
        return (20, '4f')

    def get_stats(self):
        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.axial)
            s0 = get_key0(self.axial)
            nelements = len(self.axial[s0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.axial)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  axial, torsion\n')
        return msg

    def add_f06_data(self, data, transient):
        if transient is None:
            for line in data:
                (eid, axial, torsion) = line
                self.axial[eid] = axial
                self.torsion[eid] = torsion
            return

        (dtName, dt) = transient
        self.dt = dt
        self.data_code['name'] = dtName
        if dt not in self.axial:
            self.update_dt(self.data_code, dt)

        for line in data:
            (eid, axial, torsion) = line
            self.axial[dt][eid] = axial
            self.torsion[dt][eid] = torsion

    def delete_transient(self, dt):
        del self.axial[dt]
        del self.torsion[dt]

    def get_transients(self):
        k = self.axial.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        self.dt = dt
        self.axial[dt] = {}
        self.torsion[dt] = {}

    def add_new_eid(self, dt, eid, axial, torsion):
        #print "Rod Stress add..."
        #(axial, torsion) = out
        assert isinstance(eid, int)
        self.axial[eid] = axial
        self.torsion[eid] = torsion

    def add_new_eid_sort1(self, dt, eid, axial, torsion):
        if dt not in self.axial:
            self.add_new_transient(dt)
        self.axial[dt][eid] = axial
        self.torsion[dt][eid] = torsion

    def add_new_eid_sort2(self, eid, dt, axial, torsion):
        if dt not in self.axial:
            self.add_new_transient(dt)
        self.axial[dt][eid] = axial
        self.torsion[dt][eid] = torsion

    def _write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        raise RuntimeError('this should never happen')

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.element_type == 1:
            element_header = '                           C O M P L E X   S T R E S S E S   I N   R O D   E L E M E N T S   ( C R O D )\n'
        elif self.element_type == 3:
            element_header = '                          C O M P L E X   S T R E S S E S   I N   R O D   E L E M E N T S   ( C T U B E )\n'
        elif self.element_type == 10:
            element_header = '                         C O M P L E X   S T R E S S E S   I N   R O D   E L E M E N T S   ( C O N R O D )\n'
        else:
            raise NotImplementedError('element_name=%r element_type=%s' % (self.element_name, self.element_type))

        if is_mag_phase:
            mag_phase = '                                                          (MAG/PHASE)\n'  # not tested
        else:
            mag_phase = '                                                          (REAL/IMAGINARY)\n'

        words = [element_header,
                    mag_phase,
                    ' \n',
                    '                 ELEMENT                             AXIAL                                         TORQUE\n',
                    '                   ID.                               FORCE\n',]
                    #'                       1                 -2.459512E+05 /  3.377728E+04                  0.0          /  0.0\n',]

        msg = []
        for dt, axials in sorted(iteritems(self.axial)):
            dtLine = '%14s = %12.5E\n' % (self.data_code['name'], dt)
            header[2] = dtLine
            msg += header + words
            out = []
            for eid in sorted(axials):
                axial = self.axial[dt][eid]
                torsion = self.torsion[dt][eid]
                ([axialr, torsionr, axiali, torsioni], is_all_zeros) = writeImagFloats13E([axial, torsion], is_mag_phase)
                f.write('                %8i                 %-13s / %-13s                 %-13s / %s\n' % (eid, axialr, axiali, torsionr, torsioni))

            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            page_num += 1
        return page_num - 1


class ComplexRodStrain(StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StrainObject.__init__(self, data_code, isubcase)
        self.eType = 'CROD'  # {} # 'CROD/CONROD/CTUBE'

        self.code = [self.format_code, self.sort_code, self.s_code]

        self.axial = {}
        self.torsion = {}

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
        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.axial)
            s0 = get_key0(self.axial)
            nelements = len(self.axial[s0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.axial)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  axial, torsion\n')
        return msg

    def add_f06_data(self, data, transient):
        if transient is None:
            for line in data:
                (eid, axial, torsion) = line
                self.axial[eid] = axial
                self.torsion[eid] = torsion
            return

        (dtName, dt) = transient
        self.dt = dt
        self.data_code['name'] = dtName
        if dt not in self.axial:
            self.update_dt(self.data_code, dt)

        for line in data:
            (eid, axial, torsion) = line
            self.axial[dt][eid] = axial
            self.torsion[dt][eid] = torsion

    def delete_transient(self, dt):
        del self.axial[dt]
        del self.torsion[dt]

    def get_transients(self):
        k = self.axial.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        self.dt = dt
        self.axial[self.dt] = {}
        self.torsion[self.dt] = {}

    def add_new_eid(self, dt, eid, axial, torsion):
        assert eid >= 0, eid
        self.axial[eid] = axial
        self.torsion[eid] = torsion

    def add_new_eid_sort1(self, dt, eid, axial, torsion):
        assert eid >= 0, eid
        if dt not in self.axial:
            self.add_new_transient(dt)
        self.axial[dt][eid] = axial
        self.torsion[dt][eid] = torsion

    def add_new_eid_sort2(self, eid, dt, out):
        (axial, torsion) = out
        assert eid >= 0, eid
        if dt not in self.axial:
            self.add_new_transient(dt)
        self.axial[dt][eid] = axial
        self.torsion[dt][eid] = torsion

    def _write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.dt is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        raise RuntimeError('this should never happen')

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.element_type == 1:
            element_header = '                           C O M P L E X    S T R A I N S    I N   R O D   E L E M E N T S   ( C R O D )\n'
        elif self.element_type == 3:
            element_header = '                          C O M P L E X    S T R A I N S    I N   R O D   E L E M E N T S   ( C T U B E )\n'
        elif self.element_type == 10:
            element_header = '                         C O M P L E X    S T R A I N S    I N   R O D   E L E M E N T S   ( C O N R O D )\n'
        else:
            raise NotImplementedError('element_name=%r element_type=%s' % (self.element_name, self.element_type))

        if is_mag_phase:
            mag_phase = '                                                          (MAG/PHASE)\n'  # not tested
        else:
            mag_phase = '                                                          (REAL/IMAGINARY)\n'

        words = [element_header,
                    mag_phase,
                    ' \n',
                    '                 ELEMENT                             AXIAL                                         TORQUE\n',
                    '                   ID.                               FORCE\n',
                    '                       1                 -2.459512E+05 /  3.377728E+04                  0.0          /  0.0\n',]
        words = ['                                       S T R A I N S   I N   R O D   E L E M E N T S      ( C R O D )\n',
                 '       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                 '         ID.        STRAIN       MARGIN        STRAIN      MARGIN         ID.        STRAIN       MARGIN        STRAIN      MARGIN\n']
        msg = []
        for dt, axials in sorted(iteritems(self.axial)):
            dtLine = '%14s = %12.5E\n' % (self.data_code['name'], dt)
            header[2] = dtLine
            msg += header + words
            out = []
            for eid in sorted(axials):
                axial = self.axial[dt][eid]
                torsion = self.torsion[dt][eid]
                ([axialr, torsionr, axiali, torsioni], is_all_zeros) = writeImagFloats13E([axial, torsion], is_mag_phase)
                f.write('                %8i                 %-13s / %-13s                 %-13s / %s\n' % (eid, axialr, axiali, torsionr, torsioni))
            msg.append(page_stamp % page_num)
            page_num += 1
            f.write(''.join(msg))
        return page_num - 1
