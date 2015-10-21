from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
from numpy import zeros

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, OES_Object
from pyNastran.f06.f06_formatting import writeFloats13E


class RealBush1DStressArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        self.eType = {}
        #self.code = [self.format_code, self.sort_code, self.s_code]
        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.ielement = 0
        self.nelements = 0  # result specific

    def is_real(self):
        return True

    def is_complex(self):
        return False

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def _get_msgs(self):
        raise NotImplementedError('%s needs to implement _get_msgs' % self.__class__.__name__)

    def get_headers(self):
        headers = ['element_force', 'axial_displacement', 'axial_velocity', 'axial_stress', 'axial_strain', 'plastic_strain', 'slipV', 'is_failed']
        return headers

    def build(self):
        if self.is_built:
            return
        #print("self.ielement =", self.ielement)
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal

        if self.element_type == 34:
            nnodes_per_element = 1
        else:
            raise NotImplementedError(self.element_type)

        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("***name=%s type=%s nnodes_per_element=%s ntimes=%s nelements=%s ntotal=%s" % (
            #self.element_name, self.element_type, nnodes_per_element, self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, int):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.elements = zeros(self.ntotal, dtype='int32')

        #self.element_force = {}
        #self.axial_displacement = {}
        #self.axial_velocity = {}
        #self.axial_stress = {}
        #self.axial_strain = {}
        #self.plastic_strain = {}
        #self.is_failed = {}
        # [element_force, axial_displacement, axial_velocity, axial_stress, axial_strain, plastic_strain, is_failed]
        self.data = zeros((self.ntimes, self.ntotal, 7), dtype='float32')

    def add_sort1(self, element_force, axial_displacement, axial_velocity, axial_stress, axial_strain, plastic_strain, is_failed):
        assert isinstance(eid, int)
        self._times[self.itime] = dt
        self.elements[self.itotal] = eid
        self.data[self.itime, self.itotal, :] = [element_force, axial_displacement, axial_velocity, axial_stress, axial_strain, plastic_strain, is_failed]
        self.itotal += 1
        self.ielement += 1

    def get_stats(self):
        if not self.is_built:
            return ['<%s>\n' % self.__class__.__name__,
                    '  ntimes: %i\n' % self.ntimes,
                    '  ntotal: %i\n' % self.ntotal,
                    ]

        nelements = self.ntotal
        ntimes = self.ntimes
        ntotal = self.ntotal
        nelements = self.ntotal

        msg = []
        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i\n'
                       % (self.__class__.__name__, ntimes, nelements))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i\n'
                       % (self.__class__.__name__, nelements))
            ntimes_word = 1
        headers = self.get_headers()

        n = len(headers)
        assert n == self.data.shape[2], 'nheaders=%s shape=%s' % (n, str(self.data.shape))
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.elements)  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        ind = ravel([searchsorted(self.elements == eid) for eid in eids])
        return ind

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        msg = self._get_msgs()
        (ntimes, ntotal) = self.data.shape[:2]
        eids = self.elements
        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg))

            #[element_force, axial_displacement, axial_velocity, axial_stress, axial_strain, plastic_strain, is_failed]
            element_force = self.data[itime, :, 0]
            axial_displacement = self.data[itime, :, 1]
            axial_velocity = self.data[itime, :, 2]
            axial_stress = self.data[itime, :, 3]
            axial_strain = self.data[itime, :, 4]
            plastic_strain = self.data[itime, :, 5]
            is_failed = self.data[itime, :, 6]

            # loop over all the elements
            for (i, eid, element_forcei, axial_displacementi, axial_velocityi, axial_stressi, axial_straini, plastic_straini, is_failedi) in zip(
                count(), eids, element_force, axial_displacement, axial_velocity, axial_stress, axial_strain, plastic_strain, is_failed):

                vals = [element_forcei, axial_displacementi, axial_velocityi, axial_stressi, axial_straini, plastic_straini, is_failedi]
                (vals2, is_all_zeros) = writeFloats13E(vals)
                [element_forcei, axial_displacementi, axial_velocityi, axial_stressi, axial_straini, plastic_straini, is_failedi] = vals2
                f.write('0%8i   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %%s\n'
                        % (eid, element_forcei, axial_displacementi, axial_velocityi, axial_stressi, axial_straini, plastic_straini, is_failedi))
            f.write(page_stamp % page_num)
            page_num += 1
        if self.nonlinear_factor is None:
            page_num -= 1
        return page_num


class RealBush1DStress(StressObject):
    """
    # s_code=0
                           C O M P L E X   S T R E S S E S   I N   B A R   E L E M E N T S   ( C B A R )
                                                         (MAGNITUDE/PHASE)

            ELEMENT                    LOCATION       LOCATION       LOCATION       LOCATION             AVERAGE
              ID.                          1              2              3              4             AXIAL STRESS

                  1     ENDA          9.331276E+04   9.331276E+04   9.331276E+04   9.331276E+04        0.0
                                      180.0000         0.0            0.0          180.0000              0.0
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
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
            #self.add_new_eid = self.add_new_eid_sort2

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

    def _write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase, is_sort1=is_sort1)

        raise NotImplementedError('CBUSH1D')
        msg = header + ['    S T R A I N S   I N   B U S H   E L E M E N T S        ( C B U S H 1 D)',
                        '',
                        '    ELEMENT-ID        STRAIN-TX     STRAIN-TY     STRAIN-TZ    STRAIN-RX     STRAIN-RY     STRAIN-RZ ',
                        #'0                          1      1.000000E-06  0.0           0.0           0.0           0.0           0.0',
        ]
        for eid, S1s in sorted(iteritems(self.s1)):
            element_force = self.element_force[eid]
            axial_displacement = self.axial_displacement[eid]
            axial_velocity = self.axial_velocity[eid]
            axial_stress = self.axial_stress[eid]
            axial_strain = self.axial_strain[eid]
            plastic_strain = self.plastic_strain[eid]
            is_failed = self.is_failed[eid]

            vals = [element_force, axial_displacement, axial_velocity, axial_stress, axial_strain, plastic_strain, is_failed]
            (vals2, is_all_zeros) = self.writeImagFloats13E(vals, is_mag_phase)
            [element_force, axial_displacement, axial_velocity, axial_stress, axial_strain, plastic_strain, is_failed] = vals2
            msg.append('0%8i   %-13s  %-13s  %-13s  %-13s  %s\n' % (eid, element_force, axial_displacement, axial_velocity, axial_stress, axial_strain, plastic_strain, is_failed))

        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        #raise NotImplementedError('CBUSH1D')
        words = [
#      ELEMENT-ID =    7001
                    '                    S T R E S S E S   ( F O R C E S )   I N   B U S H 1 D   E L E M E N T S   ( C B U S H 1 D )',
                    ' \n',
                    '                        AXIAL          AXIAL          AXIAL       AXIAL         AXIAL         PLASTIC\n',
                    '        TIME            FORCE       DISPLACEMENT    VELOCITY      STRESS        STRAIN        STRAIN        STATUS\n',
        ]
        msg = []
        for dt, ElementForce in sorted(iteritems(self.element_force)):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for eid, element_force in sorted(iteritems(ElementForce)):
                #eType = self.eType[eid]
                #element_force = self.element_force[dt][eid]
                axial_displacement = self.axial_displacement[dt][eid]
                axial_velocity = self.axial_velocity[dt][eid]
                axial_stress = self.axial_stress[dt][eid]
                axial_strain = self.axial_strain[dt][eid]
                plastic_strain = self.plastic_strain[dt][eid]
                is_failed = self.is_failed[dt][eid]
                vals = [element_force, axial_displacement, axial_velocity, axial_stress, axial_strain, plastic_strain, is_failed]
                (vals2, is_all_zeros) = writeFloats13E(vals)
                [element_force, axial_displacement, axial_velocity, axial_stress, axial_strain, plastic_strain, is_failed] = vals2

                msg.append(' %-13s   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (eid,
                    element_force, axial_displacement, axial_velocity, axial_stress, axial_strain, plastic_strain, is_failed))

            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            msg = ['']
            page_num += 1
        return page_num - 1
