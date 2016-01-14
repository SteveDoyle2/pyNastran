from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, OES_Object
from pyNastran.f06.f06_formatting import write_imag_floats_13e
try:
    import pandas as pd
except ImportError:
    pass


class ComplexCBush1DArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        #self.eType = {}
        #self.code = [self.format_code, self.sort_code, self.s_code]

        self.nelements = 0  # result specific

    def is_real(self):
        return False

    def is_complex(self):
        return True

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def _get_msgs(self):
        raise NotImplementedError()

    def get_headers(self):
        raise NotImplementedError()

    def build(self):
        if self.is_built:
            return

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        self.nelements //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, int):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.nelements, dtype='int32')

        #[tx, ty, tz, rx, ry, rz]
        self.data = zeros((self.ntimes, self.nelements, 6), dtype='complex64')

    def build_dataframe(self):
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
        self.data_frame.columns.names = column_names
        self.data_frame.index.names = ['ElementID', 'Item']

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()
        assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (self.table_name, table.table_name)
        assert self.approach_code == table.approach_code
        if not array_equal(self.element, table.element):
            assert self.element.shape == table.element.shape, 'shape=%s element.shape=%s' % (self.element.shape, table.element.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            for eid, eid2 in zip(self.element, table.element):
                msg += '%s, %s\n' % (eid, eid2)
            print(msg)
            raise ValueError(msg)
        if not array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            i = 0
            if self.is_sort1():
                for itime in range(ntimes):
                    for ieid, eid in enumerate(self.element):
                        t1 = self.data[itime, ieid, :]
                        t2 = table.data[itime, ieid, :]
                        (tx1, ty1, tz1, rx1, ry1, rz1) = t1
                        (tx2, ty2, tz2, rx2, ry2, rz2) = t2
                        d = t1 - t2
                        if not allclose([tx1.real, tx1.imag, ty1.real, ty1.imag],
                                        [tx2.real, tx2.imag, ty2.real, ty2.imag], atol=0.0001):
                        #if not array_equal(t1, t2):
                            msg += '%-4s  (%s, %sj, %s, %sj)\n      (%s, %sj, %s, %sj)\n  dt12=(%s, %sj, %s, %sj)\n' % (
                                eid,
                                tx1.real, tx1.imag, ty1.real, ty1.imag,
                                tx2.real, tx2.imag, ty2.real, ty2.imag,
                                d[0].real, d[0].imag, d[1].real, d[1].imag,)
                            i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
            else:
                raise NotImplementedError(self.is_sort2())
            if i > 0:
                print(msg)
                raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, tx, ty, tz, rx, ry, rz):
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [tx, ty, tz, rx, ry, rz]
        self.ielement += 1

    def get_stats(self):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        ntimes, nelements, _ = self.data.shape
        assert self.ntimes == ntimes, 'ntimes=%s expected=%s' % (self.ntimes, ntimes)
        assert self.nelements == nelements, 'nelements=%s expected=%s' % (self.nelements, nelements)

        msg = []
        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i\n'
                       % (self.__class__.__name__, ntimes, nelements))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i\n'
                       % (self.__class__.__name__, nelements))
            ntimes_word = 1
        msg.append('  eType\n')
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_element_index(self, eids):
        itot = searchsorted(eids, self.element)
        return itot

    def eid_to_element_node_index(self, eids):
        ind = searchsorted(eids, self.element)
        return ind

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        msg_temp = self.get_f06_header(is_mag_phase)

        if self.is_sort1():
            page_num = self._write_sort1_as_sort1(header, page_stamp, page_num, f, msg_temp, is_mag_phase)
        else:
            raise NotImplementedError()
        return page_num

    def _write_sort1_as_sort1(self, header, page_stamp, page_num, f, msg_temp, is_mag_phase):
        ntimes = self.data.shape[0]

        eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            tx = self.data[itime, :, 0]
            ty = self.data[itime, :, 1]
            tz = self.data[itime, :, 2]
            rx = self.data[itime, :, 3]
            ry = self.data[itime, :, 4]
            rz = self.data[itime, :, 5]
            for eid, itx, ity, itz, irx, iry, irz in zip(eids, tx, ty, tz, rx, ry, rz):
                [txr, tyr, tzr, rxr, ryr, rzr,
                 txi, tyi, tzi, rxi, ryi, rzi] = write_imag_floats_13e([itx, ity, itz, irx, iry, irz], is_mag_phase)
                #'0               1.000000E-01      0.0           2.912573E+00  0.0           0.0           0.0           0.0'
                #'                                    0.0         179.9942        0.0           0.0           0.0           0.0'
                f.write('                %8i    %-13s %-13s %-13s %-13s %-13s %s\n'
                        '                %8s    %-13s %-13s %-13s %-13s %-13s %s\n' % (
                            eid, txr, tyr, tzr, rxr, ryr, rzr,
                            '', txi, tyi, tzi, rxi, ryi, rzi))

            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class ComplexCBush1DStressArray(ComplexCBush1DArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexCBush1DArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        headers = ['tx', 'ty', 'tz', 'rx', 'ry', 'rz']
        return headers

    def get_f06_header(self, is_mag_phase=True):
        #'                         C O M P L E X   S T R E S S E S   I N   B U S H   E L E M E N T S   ( C B U S H ) '
        #' '
        #'                  FREQUENCY         STRESS-TX     STRESS-TY     STRESS-TZ    STRESS-RX     STRESS-RY     STRESS-RZ '
        #'0               1.000000E-01      0.0           2.912573E+00  0.0           0.0           0.0           0.0'
        #'                                    0.0         179.9942        0.0           0.0           0.0           0.0'

        if self.element_type == 102:
            element_header = '                         C O M P L E X   S T R E S S E S   I N   B U S H   E L E M E N T S   ( C B U S H ) \n'
            ''
            #' '
            #
            #'0               1.000000E-01      0.0           2.912573E+00  0.0           0.0           0.0           0.0'
            #'                                    0.0         179.9942        0.0           0.0           0.0           0.0'
        else:
            raise NotImplementedError('element_name=%r element_type=%s' % (self.element_name, self.element_type))

        if is_mag_phase:
            mag_phase = '                                                         (MAGNITUDE/PHASE)\n \n'
        else:
            mag_phase = '                                                          (REAL/IMAGINARY)\n \n'  # not tested

        words = [
            element_header,
            mag_phase,
            '                  FREQUENCY         STRESS-TX     STRESS-TY     STRESS-TZ    STRESS-RX     STRESS-RY     STRESS-RZ \n',
            '                   ID.                               FORCE\n',]
        return words


#class ComplexCBush1dStrainArray(ComplexCBushArray, StrainObject):
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #ComplexCBush1dArray.__init__(self, data_code, is_sort1, isubcase, dt)
        #StrainObject.__init__(self, data_code, isubcase)

    #def get_headers(self):
        #headers = ['tx', 'ty', 'tz', 'rx', 'ry', 'rz'] # tx, ty, tz, rx, ry, rz
        #return headers

    #def get_f06_header(self, is_mag_phase=True):
        ##if self.element_type == 1:
            ##element_header = '                           C O M P L E X    S T R A I N S    I N   R O D   E L E M E N T S   ( C R O D )\n'
        ##elif self.element_type == 3:
            ##element_header = '                          C O M P L E X    S T R A I N S    I N   R O D   E L E M E N T S   ( C T U B E )\n'
        ##elif self.element_type == 10:
            ##element_header = '                         C O M P L E X    S T R A I N S    I N   R O D   E L E M E N T S   ( C O N R O D )\n'
        ##else:
        #raise NotImplementedError('element_name=%r element_type=%s' % (self.element_name, self.element_type))

        #if is_mag_phase:
            #mag_phase = '                                                          (MAG/PHASE)\n'  # not tested
        #else:
            #mag_phase = '                                                          (REAL/IMAGINARY)\n'

        #words = [
            #element_header,
            #mag_phase,
            #' \n',
            #'                 ELEMENT                             AXIAL                                         TORQUE\n',
            #'                   ID.                               FORCE\n',
           ##'                       1                 -2.459512E+05 /  3.377728E+04                  0.0          /  0.0\n',
        #]
        #return words


class ComplexBush1DStress(StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StressObject.__init__(self, data_code, isubcase)
        self.eType = {}

        self.code = [self.format_code, self.sort_code, self.s_code]

        self.element_force = {}
        self.axial_displacement = {}
        self.axial_stress = {}
        self.axial_strain = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                #self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            #self.add = self.add_sort2
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
        msg.append('  eType, element_force, axial_displacement, axial_velocity, axial_stress, axial_strain\n')
        return msg

    def add_f06_data(self, data, transient):
        raise NotImplementedError('CBUSH1D')
        if transient is None:
            for line in data:
                (eType, eid, fe, ue, ao, ae) = line
                self.eType[eid] = 'CBUSH1D'
                self.element_force[eid] = fe
                self.axial_displacement[eid] = ue
                self.axial_stress[eid] = ao
                self.axial_strain[eid] = ae
            return

        (dt_name, dt) = transient
        self.data_code['name'] = dt_name

        if dt not in self.element_force:
            self.update_dt(self.data_code, dt)

        for line in data:
            (eType, eid, fe, ue, ao, ae) = line
            self.eType[eid] = 'CBUSH1D'
            self.element_force[dt][eid] = fe
            self.axial_displacement[dt][eid] = ue
            self.axial_stress[dt][eid] = ao
            self.axial_strain[dt][eid] = ae

    def delete_transient(self, dt):
        del self.element_force[dt]
        del self.axial_displacement[dt]
        del self.axial_stress[dt]
        del self.axial_strain[dt]

    def get_transients(self):
        k = self.element_force.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        #print "****add new transient****"
        self.dt = dt
        self.element_force[dt] = {}
        self.axial_displacement[dt] = {}
        self.axial_stress[dt] = {}
        self.axial_strain[dt] = {}

    def add_new_eid(self, eType, dt, eid, fe, ue, ao, ae):
        #print "Bush1d Stress add..."
        self.eType[eid] = eType

        self.element_force[eid] = fe
        self.axial_displacement[eid] = ue
        self.axial_stress[eid] = ao
        self.axial_strain[eid] = ae

    def add_new_eid_sort1(self, eType, dt, eid, fe, ue, ao, ae):
        #msg = "dt=%s eid=%s fe=%s" % (dt, eid, fe)
        #print msg
        if dt not in self.element_force:
            self.add_new_transient(dt)
        self.eType[eid] = eType
        self.element_force[dt][eid] = fe
        self.axial_displacement[dt][eid] = ue
        self.axial_stress[dt][eid] = ao
        self.axial_strain[dt][eid] = ae

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        raise NotImplementedError('CBUSH1D')
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        msg = header + [
            '                                 S T R E S S E S   I N   B A R   E L E M E N T S          ( C B A R )\n',
            '  ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T\n',
            '    ID.          SB1            SB2            SB3            SB4           STRESS         SB-MAX         SB-MIN     M.S.-C\n',
        ]

        for eid, S1s in sorted(iteritems(self.s1)):
            eType = self.eType[eid]
            axial = self.axial[eid]
            s1 = self.s1[eid]
            s2 = self.s2[eid]
            s3 = self.s3[eid]
            s4 = self.s4[eid]

            vals = [s1[0], s2[0], s3[0], s4[0], axial,
                    s1[1], s2[1], s3[1], s4[1], ]
            vals2 = write_imag_floats_13e(vals, is_mag_phase)
            [s1ar, s2ar, s3ar, s4ar, axialr,
             s1br, s2br, s3br, s4br,
             s1ai, s2ai, s3ai, s4ai, axiali,
             s1bi, s2bi, s3bi, s4bi, ] = vals2
            msg.append('0%8i   %-13s  %-13s  %-13s  %-13s  %s\n' % (eid, s1ar, s2ar, s3ar, s4ar, axialr))
            msg.append(' %8s   %-13s  %-13s  %-13s  %-13s  %s\n' % ('', s1ai, s2ai, s3ai, s4ai, axiali))

            msg.append(' %8s   %-13s  %-13s  %-13s  %s\n' % ('', s1br, s2br, s3br, s4br))
            msg.append(' %8s   %-13s  %-13s  %-13s  %s\n' % ('', s1bi, s2bi, s3bi, s4bi))

        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num

    def __write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        raise NotImplementedError('CBUSH1D')
        words = [
            '                                 S T R E S S E S   I N   B A R   E L E M E N T S          ( C B A R )\n',
            '  ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T\n',
            '    ID.          SB1            SB2            SB3            SB4           STRESS         SB-MAX         SB-MIN     M.S.-C\n',
        ]
        msg = []
        for dt, S1s in sorted(iteritems(self.s1)):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for eid, S1 in sorted(iteritems(S1s)):
                eType = self.eType[eid]
                axial = self.axial[dt][eid]
                s1 = self.s1[dt][eid]
                s2 = self.s2[dt][eid]
                s3 = self.s3[dt][eid]
                s4 = self.s4[dt][eid]
                vals = [s1[0], s2[0], s3[0], s4[0], axial,
                        s1[1], s2[1], s3[1], s4[1], ]
                vals2 = write_imag_floats_13e(vals, is_mag_phase)
                [s1ar, s2ar, s3ar, s4ar, axialr,
                 s1br, s2br, s3br, s4br,
                 s1ai, s2ai, s3ai, s4ai, axiali,
                 s1bi, s2bi, s3bi, s4bi, ] = vals2
                msg.append('0%8i   %-13s  %-13s  %-13s  %-13s  %s\n' % (eid, s1ar, s2ar, s3ar, s4ar, axialr))
                msg.append(' %8s   %-13s  %-13s  %-13s  %-13s  %s\n' % ('', s1ai, s2ai, s3ai, s4ai, axiali))

                msg.append(' %8s   %-13s  %-13s  %-13s  %s\n' % ('', s1br, s2br, s3br, s4br))
                msg.append(' %8s   %-13s  %-13s  %-13s  %s\n' % ('', s1bi, s2bi, s3bi, s4bi))

            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            msg = ['']
            page_num += 1
        return page_num - 1