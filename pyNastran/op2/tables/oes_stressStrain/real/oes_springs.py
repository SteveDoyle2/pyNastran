from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
from six.moves import zip
from numpy import zeros

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import writeFloats13E, write_float_13E, _eigenvalue_header

class RealSpringArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        #self.eType = {}

        self.nelements = 0  # result specific
        if is_sort1:
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

        #[stress]
        self.data = zeros((self.ntimes, self.nelements, 1), dtype='float32')

    def add_new_eid_sort1(self, dt, eid, stress):
        self._times[self.itime] = dt
        #if self.itime == 0:
        #print('itime=%s eid=%s' % (self.itime, eid))
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [stress]
        self.ielement += 1

    def get_stats(self):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        #print(self.data.shape[:1])
        #ntimes, nelements = self.data.shape[:1]
        ntimes = self.data.shape[0]
        nelements = self.data.shape[1]
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
        msg_temp = self.get_f06_header(is_mag_phase)

        if self.is_sort1():
            page_num = self._write_sort1_as_sort1(header, page_stamp, page_num, f, msg_temp)
        else:
            raise NotImplementedError(self.code_information())
            #page_num = self._write_sort2_as_sort2(header, page_stamp, page_num, f, msg_temp)
        return page_num

    def _write_sort1_as_sort1(self, header, page_stamp, page_num, f, msg_temp):
        ntimes = self.data.shape[0]

        eids = self.element
        is_odd = False
        nwrite = len(eids)
        nrows = nwrite // 4
        nleftover = nwrite - nrows * 4

        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))
            stress = self.data[itime, :, 0]

            out = []
            for eid, stressi in zip(eids, stress):
                out.append([eid, write_float_13E(stressi)])

            for i in range(0, nrows * 4, 4):
                f.write('    %10i  %13s    %10i  %13s    %10i  %13s    %10i  %13s\n' % (
                    tuple(out[i] + out[i + 1] + out[i + 2] + out[i + 3])))

            i = nrows * 4
            if nleftover == 3:
                f.write('    %10i  %13s    %10i  %13s    %10i  %13s\n' % (
                    tuple(out[i] + out[i + 1] + out[i + 2])))
            elif nleftover == 2:
                f.write('    %10i  %13s    %10i  %13s\n' % (
                    tuple(out[i] + out[i + 1])))
            elif nleftover == 1:
                f.write('    %10i  %13s\n' % tuple(out[i]))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealSpringStressArray(RealSpringArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealSpringArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        headers = ['spring_stress']
        return headers

    def get_f06_header(self, is_mag_phase=True):
        if self.element_type == 11:  # CELAS1
            msg = ['                              S T R E S S E S   I N   S C A L A R   S P R I N G S        ( C E L A S 1 )\n']
        elif self.element_type == 12:  # CELAS2
            msg = ['                              S T R E S S E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )\n']
        elif self.element_type == 13:  # CELAS3
            msg = ['                              S T R E S S E S   I N   S C A L A R   S P R I N G S        ( C E L A S 3 )\n']
        elif self.element_type == 14:  # CELAS4
            msg = ['                              S T R E S S E S   I N   S C A L A R   S P R I N G S        ( C E L A S 4 )\n']
        else:
            raise NotImplementedError('element_name=%s element_type=%s' % (self.element_name, self.element_type))

        msg += [
            '      ELEMENT         STRESS           ELEMENT         STRESS           ELEMENT         STRESS           ELEMENT         STRESS\n'
            '        ID.                              ID.                              ID.                              ID.\n'
        ]
        return msg



class RealSpringStrainArray(RealSpringArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealSpringArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        headers = ['spring_strain']
        return headers

    def get_f06_header(self, is_mag_phase=True):
        if self.element_type == 11:  # CELAS1
            msg = ['                               S T R A I N S    I N   S C A L A R   S P R I N G S        ( C E L A S 1 )\n']
        elif self.element_type == 12:  # CELAS2
            msg = ['                               S T R A I N S    I N   S C A L A R   S P R I N G S        ( C E L A S 2 )\n']
        elif self.element_type == 13:  # CELAS3
            msg = ['                               S T R A I N S    I N   S C A L A R   S P R I N G S        ( C E L A S 3 )\n']
        elif self.element_type == 14:  # CELAS4
            msg = ['                               S T R A I N S    I N   S C A L A R   S P R I N G S        ( C E L A S 4 )\n']
        else:
            raise NotImplementedError('element_name=%s element_type=%s' % (self.element_name, self.element_type))

        msg += [
            '      ELEMENT         STRAIN           ELEMENT         STRAIN           ELEMENT         STRAIN           ELEMENT         STRAIN\n'
            '        ID.                              ID.                              ID.                              ID.\n'
        ]
        return msg


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
