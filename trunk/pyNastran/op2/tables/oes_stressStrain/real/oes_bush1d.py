from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from ..real.oes_objects import StressObject
from pyNastran.f06.f06_formatting import writeFloats13E


class Bush1DStressObject(StressObject):
    """
    # s_code=0
                           C O M P L E X   S T R E S S E S   I N   B A R   E L E M E N T S   ( C B A R )
                                                         (MAGNITUDE/PHASE)

            ELEMENT                    LOCATION       LOCATION       LOCATION       LOCATION             AVERAGE
              ID.                          1              2              3              4             AXIAL STRESS

                  1     ENDA          9.331276E+04   9.331276E+04   9.331276E+04   9.331276E+04        0.0
                                      180.0000         0.0            0.0          180.0000              0.0
    """
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        StressObject.__init__(self, data_code, isubcase)
        self.eType = {}
        self.code = [self.format_code, self.sort_code, self.s_code]

        self.element_force = {}
        self.axial_displacement = {}
        self.axial_velocity = {}
        self.axial_stress = {}
        self.axial_strain = {}
        self.plastic_strain = {}
        self.is_failed = {}

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
        if self.dt is not None:  # transient
            ntimes = len(self.element_force)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  imaginary type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, element_force, axial_displacement, axial_velocity, axial_stress, axial_strain, plastic_strain, is_failed\n')
        return msg

    def add_f06_data(self, data, transient):
        raise NotImplementedError('CBUSH1D')
        if transient is None:
            for line in data:
                (eType, eid, fe, ue, ve, ao, ae, ep, fail) = line
                self.eType[eid] = 'CBUSH1D'
                self.element_force[eid] = fe
                self.axial_displacement[eid] = ue
                self.axial_velocity[eid] = ve
                self.axial_stress[eid] = ao
                self.axial_strain[eid] = ae
                self.plastic_strain[eid] = ep
                self.is_failed[eid] = fail
            return

        (dtName, dt) = transient
        self.data_code['name'] = dtName

        if dt not in self.element_force:
            self.update_dt(self.data_code, dt)

        for line in data:
            (eType, eid, fe, ue, ao, ae) = line
            self.eType[eid] = 'CBUSH1D'
            self.element_force[dt][eid] = fe
            self.axial_displacement[dt][eid] = ue
            self.axial_velocity[dt][eid] = ve
            self.axial_stress[dt][eid] = ao
            self.axial_strain[dt][eid] = ae
            self.plastic_strain[dt][eid] = ep
            self.is_failed[dt][eid] = fail

    def delete_transient(self, dt):
        del self.element_force[dt]
        del self.axial_displacement[dt]
        del self.axial_velocity[dt]
        del self.axial_stress[dt]
        del self.axial_strain[dt]
        del self.plastic_strain[dt]
        del self.is_failed[dt]

    def get_transients(self):
        k = self.element_force.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        self.dt = dt
        self.element_force[dt] = {}
        self.axial_displacement[dt] = {}
        self.axial_velocity[dt] = {}
        self.axial_stress[dt] = {}
        self.axial_strain[dt] = {}
        self.plastic_strain[dt] = {}
        self.is_failed[dt] = {}

    def add_new_eid(self, eType, dt, eid, fe, ue, ve, ao, ae, ep, fail):
        self.eType[eid] = eType
        self.element_force[eid] = fe
        self.axial_displacement[eid] = ue
        self.axial_velocity[eid] = ve
        self.axial_stress[eid] = ao
        self.axial_strain[eid] = ae
        self.plastic_strain[eid] = ep
        self.is_failed[eid] = fail

    def add_new_eid_sort1(self, eType, dt, eid, fe, ue, ve, ao, ae, ep, fail):
        if dt not in self.element_force:
            self.add_new_transient(dt)
        self.eType[eid] = eType
        self.element_force[dt][eid] = fe
        self.axial_displacement[dt][eid] = ue
        self.axial_velocity[dt][eid] = ve
        self.axial_stress[dt][eid] = ao
        self.axial_strain[dt][eid] = ae
        self.plastic_strain[dt][eid] = ep
        self.is_failed[dt][eid] = fail

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f, is_mag_phase)

        raise NotImplementedError('CBUSH1D')
        msg = header + [
            '                                 S T R E S S E S   I N   B A R   E L E M E N T S          ( C B A R )\n',
            '  ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T\n',
            '    ID.          SB1            SB2            SB3            SB4           STRESS         SB-MAX         SB-MIN     M.S.-C\n',
        ]

        for eid, S1s in sorted(self.s1.iteritems()):
            element_force = self.element_force[eid]
            axial_displacement = self.axial_displacement[eid]
            axial_velocity = self.axial_velocity[eid]
            axial_stress = self.axial_stress[eid]
            axial_strain = self.axial_strain[eid]
            plastic_strain = self.plastic_strain[eid]
            is_failed = self.is_failed[eid]
            #eType = self.eType[eid]
            #axial = self.axial[eid]
            #s1 = self.s1[eid]
            #s2 = self.s2[eid]
            #s3 = self.s3[eid]
            #s4 = self.s4[eid]

            vals = [element_force, axial_displacement, axial_velocity, axial_stress, axial_strain, plastic_strain, is_failed]
            (vals2, isAllZeros) = self.writeImagFloats13E(vals, is_mag_phase)
            [element_force, axial_displacement, axial_velocity, axial_stress, axial_strain, plastic_strain, is_failed] = vals2
            msg.append('0%8i   %-13s  %-13s  %-13s  %-13s  %s\n' % (eid, element_force, axial_displacement, axial_velocity, axial_stress, axial_strain, plastic_strain, is_failed))
            msg.append(' %8s   %-13s  %-13s  %-13s  %-13s  %s\n' % ('', s1ai, s2ai, s3ai, s4ai, axiali))

            msg.append(' %8s   %-13s  %-13s  %-13s  %s\n' % ('', s1br, s2br, s3br, s4br))
            msg.append(' %8s   %-13s  %-13s  %-13s  %s\n' % ('', s1bi, s2bi, s3bi, s4bi))

        msg.append(pageStamp % pageNum)
        f.write(''.join(msg))
        return pageNum

    def _write_f06_transient(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        #raise NotImplementedError('CBUSH1D')
        words = [
#      ELEMENT-ID =    7001
                    '                    S T R E S S E S   ( F O R C E S )   I N   B U S H 1 D   E L E M E N T S   ( C B U S H 1 D )',
                    ' \n',
                    '                        AXIAL          AXIAL          AXIAL       AXIAL         AXIAL         PLASTIC\n',
                    '        TIME            FORCE       DISPLACEMENT    VELOCITY      STRESS        STRAIN        STRAIN        STATUS\n',
        ]
        msg = []
        for dt, ElementForce in sorted(self.element_force.iteritems()):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for eid, element_force in sorted(ElementForce.iteritems()):
                #eType = self.eType[eid]
                #element_force = self.element_force[dt][eid]
                axial_displacement = self.axial_displacement[dt][eid]
                axial_velocity = self.axial_velocity[dt][eid]
                axial_stress = self.axial_stress[dt][eid]
                axial_strain = self.axial_strain[dt][eid]
                plastic_strain = self.plastic_strain[dt][eid]
                is_failed = self.is_failed[dt][eid]
                vals = [element_force, axial_displacement, axial_velocity, axial_stress, axial_strain, plastic_strain, is_failed]
                (vals2, isAllZeros) = writeFloats13E(vals)
                [element_force, axial_displacement, axial_velocity, axial_stress, axial_strain, plastic_strain, is_failed] = vals2

                msg.append(' %-13s   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (eid,
                    element_force, axial_displacement, axial_velocity, axial_stress, axial_strain, plastic_strain, is_failed))

            msg.append(pageStamp % pageNum)
            f.write(''.join(msg))
            msg = ['']
            pageNum += 1
        return pageNum - 1
