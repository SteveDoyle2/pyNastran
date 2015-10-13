from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
from six.moves import zip, range
from numpy import zeros, searchsorted

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import writeFloats13E, _eigenvalue_header, get_key0


# there is a bug for mode_solid_shell_bar.op2 for multiple times
class RealRodArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        self.eType = {}
        #self.code = [self.format_code, self.sort_code, self.s_code]

        self.nelements = 0  # result specific

        if is_sort1:
            #sort1
            self.add_new_eid = self.add_new_eid_sort1
        else:
            raise NotImplementedError('SORT2')

    def is_real(self):
        return True

    def is_complex(self):
        return False

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def _get_msgs(self):
        raise NotImplementedError()

    def get_headers(self):
        raise NotImplementedError()

    def build(self):
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        self.nelements //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, int):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.nelements, dtype='int32')

        #[axial, torsion, SMa, SMt]
        self.data = zeros((self.ntimes, self.nelements, 4), dtype='float32')

    def add_new_eid_sort1(self, dt, eid, axial, SMa, torsion, SMt):
        self._times[self.itime] = dt
        #if self.itime == 0:
        #print('itime=%s eid=%s' % (self.itime, eid))
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [axial, SMa, torsion, SMt]
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

    def get_f06_header(self, is_mag_phase=True):
        crod_msg, conrod_msg, ctube_msg = self._get_msgs()
        if 'CROD' in self.element_name:
            msg = crod_msg
        elif 'CONROD' in self.element_name:
            msg = conrod_msg
        elif 'CTUBE' in self.element_name:
            msg = ctube_msg
        else:
            raise NotImplementedError(self.element_name)
        return self.element_name, msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element)  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.element_node[:, 0] == eid) for eid in eids])
        ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        (elem_name, msg_temp) = self.get_f06_header(is_mag_phase)

        if self.is_sort1():
            page_num = self._write_sort1_as_sort1(header, page_stamp, page_num, f, msg_temp)
        return page_num

    def _write_sort1_as_sort1(self, header, page_stamp, page_num, f, msg_temp):
        ntimes = self.data.shape[0]

        eids = self.element
        is_odd = False
        nwrite = len(eids)
        if len(eids) % 2 == 1:
            nwrite -= 1
            is_odd = True

        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            # TODO: can I get this without a reshape?
            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            axial = self.data[itime, :, 0]
            SMa = self.data[itime, :, 1]
            torsion = self.data[itime, :, 2]
            SMt = self.data[itime, :, 3]

            # loop over all the elements
            out = []
            for eid, axiali, SMai, torsioni, SMti in zip(eids, axial, SMa, torsion, SMt):
                ([axiali, torsioni, SMai, SMti],
                is_all_zeros) = writeFloats13E([axiali, torsioni, SMai, SMti])
                out.append([eid, axiali, SMai, torsioni, SMti])

            for i in range(0, nwrite, 2):
                out_line = '      %8i %-13s  %-13s %-13s  %-13s %-8i   %-13s  %-13s %-13s  %-s\n' % (tuple(out[i] + out[i + 1]))
                f.write(out_line)
            if is_odd:
                out_line = '      %8i %-13s  %-13s %-13s  %13s\n' % (tuple(out[-1]))
                f.write(out_line)
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealBushStressArray(RealRodArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealRodArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        headers = ['axial', 'SMa', 'torsion', 'SMt']
        return headers

    def _get_msgs(self):
        base_msg = ['       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                    '         ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN\n']
        crod_msg   = ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )\n', ]
        conrod_msg = ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C O N R O D )\n', ]
        ctube_msg  = ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C T U B E )\n', ]
        #cbush_msg  = ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C B U S H )\n', ]
        crod_msg += base_msg
        conrod_msg += base_msg
        ctube_msg += base_msg
        #cbush_msg += base_msg
        return crod_msg, conrod_msg, ctube_msg


class RealRodStressArray(RealRodArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealRodArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        headers = ['axial', 'SMa', 'torsion', 'SMt']
        return headers

    def _get_msgs(self):
        base_msg = ['       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                    '         ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN\n']
        crod_msg   = ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )\n', ]
        conrod_msg = ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C O N R O D )\n', ]
        ctube_msg  = ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C T U B E )\n', ]
        crod_msg += base_msg
        conrod_msg += base_msg
        ctube_msg += base_msg
        return crod_msg, conrod_msg, ctube_msg

class RealRodStrainArray(RealRodArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealRodArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        headers = ['axial', 'SMa', 'torsion', 'SMt']
        return headers

    def _get_msgs(self):
        base_msg = ['       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                    '         ID.        STRAIN       MARGIN        STRAIN      MARGIN         ID.        STRAIN       MARGIN        STRAIN      MARGIN\n']
        crod_msg   = ['                                       S T R A I N S   I N   R O D   E L E M E N T S      ( C R O D )\n', ]
        conrod_msg = ['                                       S T R A I N S   I N   R O D   E L E M E N T S      ( C O N R O D )\n', ]
        ctube_msg  = ['                                       S T R A I N S   I N   R O D   E L E M E N T S      ( C T U B E )\n', ]
        crod_msg += base_msg
        conrod_msg += base_msg
        ctube_msg += base_msg
        return crod_msg, conrod_msg, ctube_msg


class RodDamper(StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StressObject.__init__(self, data_code, isubcase)
        self.eType = 'CBUSH'

        self.code = [self.format_code, self.sort_code, self.s_code]
        self.axial = {}
        self.torsion = {}

    def get_stats(self):
        msg = self.get_data_code()
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.axial)
            a0 = get_key0(self.stress)
            nelements = len(self.axial[a0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.axial)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, axial, torsion\n')
        msg.append('  eTypes = %s\n' % self.element_name)
        return msg


class RealRodStress(StressObject):
    """
    ::

      # format_code=1 stressCode=0
                                       S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )
         ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY
           ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN
               1    5.000000E+03              0.0                               2    0.0                       0.0

      # format_code=1 stressCode=0
                                       S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )
        ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY
          ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StressObject.__init__(self, data_code, isubcase)
        self.eType = 'CROD'

        self.code = [self.format_code, self.sort_code, self.s_code]
        self.axial = {}
        self.torsion = {}

        self.MS_axial = {}
        self.MS_torsion = {}
        self.isImaginary = False

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
            a0 = get_key0(self.axial)
            nelements = len(self.axial[a0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.axial)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, axial, torsion, MS_axial, MS_torsion\n')
        return msg

    def getLength(self):
        return (20, '4f')

    def add_f06_data(self, data, transient):
        if transient is None:
            for line in data:
                (eid, axial, MSa, torsion, MSt) = line
                if MSa is None:
                    MSa = 0.
                if MSt is None:
                    MSt = 0.
                self.axial[eid] = axial
                self.MS_axial[eid] = MSa
                self.torsion[eid] = torsion
                self.MS_torsion[eid] = MSt
            return

        (dtName, dt) = transient
        self.dt = dt
        self.data_code['name'] = dtName
        if dt not in self.axial:
            self.update_dt(self.data_code, dt)

        for line in data:
            (eid, axial, MSa, torsion, MSt) = line
            if MSa is None:
                MSa = 0.
            if MSt is None:
                MSt = 0.
            self.axial[dt][eid] = axial
            self.MS_axial[dt][eid] = MSa
            self.torsion[dt][eid] = torsion
            self.MS_torsion[dt][eid] = MSt

    def delete_transient(self, dt):
        del self.axial[dt]
        del self.torsion[dt]
        del self.MS_axial[dt]
        del self.MS_torsion[dt]

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
        self.MS_axial[dt] = {}
        self.torsion[dt] = {}
        self.MS_torsion[dt] = {}

    def add_new_eid(self, dt, eid, axial, SMa, torsion, SMt):
        #print "Rod Stress add..."
        assert isinstance(eid, int)
        self.axial[eid] = axial
        self.MS_axial[eid] = SMa
        self.torsion[eid] = torsion
        self.MS_torsion[eid] = SMt

    def add_new_eid_sort1(self, dt, eid, axial, SMa, torsion, SMt):
        if dt not in self.axial:
            self.add_new_transient(dt)
        self.axial[dt][eid] = axial
        self.MS_axial[dt][eid] = SMa
        self.torsion[dt][eid] = torsion
        self.MS_torsion[dt][eid] = SMt

    def add_new_eid_sort2(self, eid, dt, axial, SMa, torsion, SMt):
        if dt not in self.axial:
            self.add_new_transient(dt)
        self.axial[dt][eid] = axial
        self.MS_axial[dt][eid] = SMa
        self.torsion[dt][eid] = torsion
        self.MS_torsion[dt][eid] = SMt

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        words = header + ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )\n',
                        '       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                        '         ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN\n']
        return self._write_f06(words, page_stamp, page_num, f)

    def _write_f06(self, msg, page_stamp, page_num=1, f=None):
        out = []
        for eid in sorted(self.axial):
            axial = self.axial[eid]
            MSa = self.MS_axial[eid]
            torsion = self.torsion[eid]
            MSt = self.MS_torsion[eid]
            (vals2, is_all_zeros) = writeFloats13E([axial, torsion])
            (axial, torsion) = vals2
            out.append([eid, axial, MSa, torsion, MSt])

        nOut = len(out)
        nWrite = nOut
        if nOut % 2 == 1:
            nWrite = nOut - 1
        for i in range(0, nWrite, 2):
            #print i,out[i:]
            out_line = '      %8i   %13s  %10.4E %13s  %10.4E   %8i   %13s  %10.4E %13s  %10.4E\n' % (tuple(out[i] + out[i + 1]))
            msg.append(out_line)

        if nOut % 2 == 1:
            out_line = '      %8i   %13s  %10.4E %13s  %10.4E\n' % (
                tuple(out[-1]))
            msg.append(out_line)
        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )\n',
                 '       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                 '         ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN\n']
        msg = []
        itime = 0
        ntimes = len(self.axial)
        for dt, axials in sorted(iteritems(self.axial)):
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            msg += header + words
            out = []
            for eid in sorted(axials):
                axial = self.axial[dt][eid]
                MSa = self.MS_axial[dt][eid]
                torsion = self.torsion[dt][eid]
                MSt = self.MS_torsion[dt][eid]

                (vals2, is_all_zeros) = writeFloats13E([axial, torsion])
                (axial, torsion) = vals2
                out.append([eid, axial, MSa, torsion, MSt])

            nOut = len(out)
            nWrite = nOut
            if nOut % 2 == 1:
                nWrite = nOut - 1
            for i in range(0, nWrite, 2):
                out_line = '      %8i   %13s  %10.4E %13s  %10.4E   %8i   %13s  %10.4E %13s  %10.4E\n' % (tuple(out[i] + out[i + 1]))
                msg.append(out_line)

            if nOut % 2 == 1:
                out_line = '      %8i   %13s  %10.4E %13s  %10.4E\n' % (tuple(out[-1]))
                msg.append(out_line)
            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            page_num += 1
            itime += 1
        return page_num - 1


class RealRodStrain(StrainObject):
    """
    ::

      # s_code=1
                                       S T R A I N S   I N   R O D   E L E M E N T S      ( C R O D )
      ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY
        ID.        STRAIN       MARGIN        STRAIN      MARGIN

      # s_code=10
                                         S T R A I N S   I N   R O D   E L E M E N T S      ( C O N R O D )
      ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY
        ID.        STRAIN       MARGIN        STRAIN      MARGIN         ID.        STRAIN       MARGIN        STRAIN      MARGIN
         1001    1.000000E+00   1.0E+00    1.250000E+00   3.0E+00         1007    1.000000E+00   1.0E+00    1.250000E+00   3.0E+00
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        assert not isinstance(data_code, bool)
        StrainObject.__init__(self, data_code, isubcase)
        self.eType = 'CROD'  # {} # 'CROD/CONROD/CTUBE'

        self.code = [self.format_code, self.sort_code, self.s_code]

        self.axial = {}
        self.torsion = {}

        self.MS_axial = {}
        self.MS_torsion = {}
        self.isImaginary = False

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
            a0 = get_key0(self.axial)
            nelements = len(self.axial[a0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.axial)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, axial, torsion, MS_axial, MS_torsion\n')
        return msg

    def add_f06_data(self, data, transient):
        if transient is None:
            for line in data:
                (eid, axial, MSa, torsion, MSt) = line
                if MSa is None:
                    MSa = 0.
                if MSt is None:
                    MSt = 0.
                self.axial[eid] = axial
                self.MS_axial[eid] = MSa
                self.torsion[eid] = torsion
                self.MS_torsion[eid] = MSt
            return

        (dtName, dt) = transient
        self.dt = dt
        self.data_code['name'] = dtName
        if dt not in self.axial:
            self.update_dt(self.data_code, dt)

        for line in data:
            (eid, axial, MSa, torsion, MSt) = line
            if MSa is None:
                MSa = 0.
            if MSt is None:
                MSt = 0.
            self.axial[dt][eid] = axial
            self.MS_axial[dt][eid] = MSa
            self.torsion[dt][eid] = torsion
            self.MS_torsion[dt][eid] = MSt

    def getLength(self):
        return (20, '4f')

    def delete_transient(self, dt):
        del self.axial[dt]
        del self.torsion[dt]
        del self.MS_axial[dt]
        del self.MS_torsion[dt]

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
        self.MS_axial[self.dt] = {}
        self.torsion[self.dt] = {}
        self.MS_torsion[self.dt] = {}

    def add_new_eid(self, dt, eid, axial, SMa, torsion, SMt):
        assert eid >= 0, eid
        self.axial[eid] = axial
        self.MS_axial[eid] = SMa
        self.torsion[eid] = torsion
        self.MS_torsion[eid] = SMt

    def add_new_eid_sort1(self, dt, eid, axial, SMa, torsion, SMt):
        assert eid >= 0, eid
        if dt not in self.axial:
            self.add_new_transient(dt)
        self.axial[dt][eid] = axial
        self.MS_axial[dt][eid] = SMa
        self.torsion[dt][eid] = torsion
        self.MS_torsion[dt][eid] = SMt

    def add_new_eid_sort2(self, eid, dt, axial, SMa, torsion, SMt):
        assert eid >= 0, eid
        #self.eType[eid] = self.element_type
        if dt not in self.axial:
            self.add_new_transient(dt)
        self.axial[dt][eid] = axial
        self.MS_axial[dt][eid] = SMa
        self.torsion[dt][eid] = torsion
        self.MS_torsion[dt][eid] = SMt

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.dt is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        words = header + ['                                       S T R A I N S   I N   R O D   E L E M E N T S      ( C R O D )\n',
                        '       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                        '         ID.        STRAIN       MARGIN        STRAIN      MARGIN         ID.        STRAIN       MARGIN        STRAIN      MARGIN\n']
        return self._write_f06(words, page_stamp, page_num, f)

    def _write_f06(self, msg, page_stamp, page_num, f):
        out = []
        for eid in sorted(self.axial):
            axial = self.axial[eid]
            MSa = self.MS_axial[eid]
            torsion = self.torsion[eid]
            MSt = self.MS_torsion[eid]
            (vals2, is_all_zeros) = writeFloats13E([axial, torsion])
            (axial, torsion) = vals2
            out.append([eid, axial, MSa, torsion, MSt])

        nOut = len(out)
        nWrite = nOut
        if nOut % 2 == 1:
            nWrite = nOut - 1
        for i in range(0, nWrite, 2):
            out_line = '      %8i   %13s  %10.4E %13s  %10.4E   %8i   %13s  %10.4E %13s  %10.4E\n' % (tuple(out[i] + out[i + 1]))
            msg.append(out_line)

        if nOut % 2 == 1:
            out_line = '      %8i   %13s  %10.4E %13s  %10.4E\n' % (tuple(out[-1]))
            msg.append(out_line)
        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                                       S T R A I N S   I N   R O D   E L E M E N T S      ( C R O D )\n',
                 '       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                 '         ID.        STRAIN       MARGIN        STRAIN      MARGIN         ID.        STRAIN       MARGIN        STRAIN      MARGIN\n']
        msg = []
        itime = 0
        ntimes = len(self.axial)
        for dt, axials in sorted(iteritems(self.axial)):
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            msg += header + words
            out = []
            for eid in sorted(axials):
                axial = self.axial[dt][eid]
                MSa = self.MS_axial[dt][eid]
                torsion = self.torsion[dt][eid]
                MSt = self.MS_torsion[dt][eid]

                out.append([eid, axial, MSa, torsion, MSt])

            nOut = len(out)
            nWrite = nOut
            if nOut % 2 == 1:
                nWrite = nOut - 1
            for i in range(0, nWrite, 2):
                out_line = '      %8i   %13.6E  %10.4E %13.6E  %10.4E   %8i   %13.6E  %10.4E %13.6E  %10.4E\n' % (tuple(out[i] + out[i + 1]))
                msg.append(out_line)

            if nOut % 2 == 1:
                out_line = '      %8i   %13.6E  %10.4E %13.6E  %10.4E\n' % (
                    tuple(out[-1]))
                msg.append(out_line)
            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            page_num += 1
            itime += 1
        return page_num - 1


class ConrodStress(RealRodStress):
    eType = 'CONROD'
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealRodStress.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        words = header + ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C O N R O D )\n',
                        '       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                        '         ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN\n']
        return self._write_f06(words, page_stamp, page_num, f)


class CtubeStress(RealRodStress):
    eType = 'CTUBE'
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealRodStress.__init__(self, is_sort1, data_code, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        words = header + ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C T U B E )\n',
                        '       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                        '         ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN\n']
        return self._write_f06(words, page_stamp, page_num, f)


class ConrodStrain(RealRodStrain):
    eType = 'CONROD'
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealRodStrain.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        words = ['                                       S T R A I N S   I N   R O D   E L E M E N T S      ( C O N R O D )\n',
                 '       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                 '         ID.        STRAIN       MARGIN        STRAIN      MARGIN         ID.        STRAIN       MARGIN        STRAIN      MARGIN\n']
        return self._write_f06(words, page_stamp, page_num, f)

class CtubeStrain(RealRodStrain):
    eType = 'CTUBE'
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealRodStrain.__init__(self, is_sort1, data_code, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        words = ['                                       S T R A I N S   I N   R O D   E L E M E N T S      ( C T U B E )\n',
                 '       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                 '         ID.        STRAIN       MARGIN        STRAIN      MARGIN         ID.        STRAIN       MARGIN        STRAIN      MARGIN\n']
        return self._write_f06(words, page_stamp, page_num, f)
