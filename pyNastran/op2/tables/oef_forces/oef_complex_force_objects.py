from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from pyNastran.op2.result_objects.op2_objects import ScalarObject
from pyNastran.f06.f06_formatting import write_imag_floats_13e, get_key0, write_float_12E
from pyNastran.f06.f06_formatting import _eigenvalue_header
import numpy as np
from numpy import zeros, array_equal, searchsorted, allclose
try:
    import pandas as pd
except ImportError:
    pass

class ComplexRodForceArray(ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        if is_sort1:
            pass
        else:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        headers = ['axial_force', 'torque']
        return headers

    def is_real(self):
        return True

    def is_complex(self):
        return False

    #def get_headers(self):
        #headers = ['axial', 'torque']
        #return headers

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

        #[axial_force, torque]
        self.data = zeros((self.ntimes, self.ntotal, 2), dtype='complex64')

    def build_dataframe(self):
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
        self.data_frame.columns.names = column_names
        self.data_frame.index.names = ['ElementID', 'Item']

    def __eq__(self, table):
        self._eq_header(table)
        assert self.is_sort1() == table.is_sort1()
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, e in enumerate(self.element):
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (axial1, torque1) = t1
                    (axial2, torque2) = t2

                    if not allclose(t1, t2):
                        msg += '(%s, %s)    (%s, %s)  (%s, %s)\n' % (
                            eid, nid,
                            axial1, torque1,
                            axial2, torque2)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, axial, torque):
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [axial, torque]
        self.ielement += 1

    def get_stats(self):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        #ntotal = self.ntotal

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
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True, is_sort1=True):
        if self.element_type == 1: # CROD
            msg = ['                             C O M P L E X   F O R C E S   I N   R O D   E L E M E N T S   ( C R O D )\n']
        elif self.element_type == 10: # CONROD
            msg = ['                           C O M P L E X   F O R C E S   I N   R O D   E L E M E N T S   ( C O N R O D )\n']
        elif self.element_type == 3: # CTUBE
            msg = ['                            C O M P L E X   F O R C E S   I N   R O D   E L E M E N T S   ( C T U B E )\n']
            #pass
        else:
            raise NotImplementedError('element_name=%s element_type=%s' % (self.element_name, self.element_type))

        if is_mag_phase:
            msg += ['                                                          (MAGNITUDE/PHASE)\n']
        else:
            msg += ['                                                          (REAL/IMAGINARY)\n']

        if is_sort1:
            msg += [
                ' \n'
                '                 ELEMENT                             AXIAL                                       TORSIONAL\n'
                '                   ID.                              STRAIN                                         STRAIN\n'
            ]
            #'                      14                  0.0          /  0.0                           0.0          /  0.0'
        else:
            raise NotImplementedError('sort2')


        return self.element_name, msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element)  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.element == eid) for eid in eids])
        ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        (elem_name, msg_temp) = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        # write the f06
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element
        is_odd = False
        nwrite = len(eids)
        if len(eids) % 2 == 1:
            nwrite -= 1
            is_odd = True

        #print('len(eids)=%s nwrite=%s is_odd=%s' % (len(eids), nwrite, is_odd))
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            axial = self.data[itime, :, 0]
            torsion = self.data[itime, :, 1]

            for eid, axiali, torsioni in zip(eids, axial, torsion):
                [raxial, rtorsion, iaxial, itorsion] = write_imag_floats_13e([axiali, torsioni], is_mag_phase)
                #ELEMENT                             AXIAL                                       TORSIONAL
                    #ID.                              STRESS                                         STRESS
                    #14                  0.0          /  0.0                           0.0          /  0.0

                f.write('      %8i   %-13s / %-13s   %-13s / %s\n' % (eid, raxial, iaxial, rtorsion, itorsion))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class ComplexCShearForceArray(ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        if is_sort1:
            pass
        else:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        headers = [
            'force41', 'force14', 'force21', 'force12', 'force32', 'force23',
            'force43', 'force34', 'kickForce1', 'kickForce2', 'kickForce3',
            'kickForce4', 'shear12', 'shear23', 'shear34', 'shear41'
        ]
        return headers

    #def get_headers(self):
        #headers = ['axial', 'torque']
        #return headers

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

        #[force41, force14, force21, force12, force32, force23, force43, force34,
        #kick_force1, kick_force2, kick_force3, kick_force4,
        #shear12, shear23, shear34, shear41]
        self.data = zeros((self.ntimes, self.ntotal, 16), dtype='complex64')

    def build_dataframe(self):
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
        self.data_frame.columns.names = column_names
        self.data_frame.index.names = ['ElementID', 'Item']

    def __eq__(self, table):
        self._eq_header(table)
        assert self.is_sort1() == table.is_sort1()
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, eid in enumerate(self.element):
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (force41a, force14a, force21a, force12a, force32a, force23a, force43a, force34a, kick_force1a, kick_force2a, kick_force3a, kick_force4a, shear12a, shear23a, shear34a, shear41a) = t1
                    (force41b, force14b, force21b, force12b, force32b, force23b, force43b, force34b, kick_force1b, kick_force2b, kick_force3b, kick_force4b, shear12b, shear23b, shear34b, shear41b) = t2

                    if not allclose(t1, t2):
                        msg += (
                            '%s   (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)\n'
                            '     (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                                eid,
                                force41a, force14a, force21a, force12a, force32a, force23a, force43a, force34a, kick_force1a, kick_force2a, kick_force3a, kick_force4a, shear12a, shear23a, shear34a, shear41a,
                                force41b, force14b, force21b, force12b, force32b, force23b, force43b, force34b, kick_force1b, kick_force2b, kick_force3b, kick_force4b, shear12b, shear23b, shear34b, shear41b
                            ))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid,
                  force41, force14, force21, force12, force32, force23, force43, force34,
                  kick_force1, kick_force2, kick_force3, kick_force4,
                  shear12, shear23, shear34, shear41):
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [
            force41, force14, force21, force12, force32, force23, force43, force34,
            kick_force1, kick_force2, kick_force3, kick_force4,
            shear12, shear23, shear34, shear41]
        self.ielement += 1

    def get_stats(self):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        #ntotal = self.ntotal

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
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True, is_sort1=True):
        msg = ['                C O M P L E X   F O R C E S   A C T I N G   O N   S H E A R   P A N E L   E L E M E N T S   (CSHEAR)\n']

        if is_mag_phase:
            msg += ['                                                          (MAGNITUDE/PHASE)\n \n']
        else:
            msg += ['                                                          (REAL/IMAGINARY)\n \n']

        if is_sort1:
            msg += [
                '                  ====== POINT  1 ======      ====== POINT  2 ======      ====== POINT  3 ======      ====== POINT  4 ======'
                ' ELEMENT          F-FROM-4      F-FROM-2      F-FROM-1      F-FROM-3      F-FROM-2      F-FROM-4      F-FROM-3      F-FROM-1'
                '         ID               KICK-1       SHEAR-12       KICK-2       SHEAR-23       KICK-3       SHEAR-34       KICK-4       SHEAR-41'
            ]
        else:
            raise NotImplementedError('sort2')
        return msg

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        # write the f06
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element

        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))

            ## TODO: I'm sure this ordering is wrong...
            force41 = self.data[itime, :, 0]
            force14 = self.data[itime, :, 1]
            force21 = self.data[itime, :, 2]
            force12 = self.data[itime, :, 3]
            force32 = self.data[itime, :, 4]
            force23 = self.data[itime, :, 5]
            force43 = self.data[itime, :, 6]
            force34 = self.data[itime, :, 7]
            kick_force1 = self.data[itime, :, 8]
            kick_force2 = self.data[itime, :, 9]
            kick_force3 = self.data[itime, :, 10]
            kick_force4 = self.data[itime, :, 11]
            shear12 = self.data[itime, :, 12]
            shear23 = self.data[itime, :, 13]
            shear34 = self.data[itime, :, 14]
            shear41 = self.data[itime, :, 15]

            for (eid, iforce41, force14i, iforce21, iforce12, iforce32, iforce23, iforce43, iforce34,
                 ikick_force1, ikick_force2, ikick_force3, ikick_force4,
                 ishear12, ishear23, ishear34, ishear41) in zip(
                     eids, force41, force14, force21, force12, force32, force23, force43, force34,
                     kick_force1, kick_force2, kick_force3, kick_force4,
                     shear12, shear23, shear34, shear41):

                vals2 = write_imag_floats_13e([
                     iforce41, force14i, iforce21, iforce12, iforce32, iforce23, iforce43, iforce34,
                     ikick_force1, ikick_force2, ikick_force3, ikick_force4,
                     ishear12, ishear23, ishear34, ishear41], is_mag_phase)

                [
                    force41r, force14r, force21i, force12r, force32r, force23r, force43r, force34r,
                    kick_force1r, kick_force2r, kick_force3r, kick_force4r,
                    shear12r, shear23r, shear34r, shear41r,

                    force41i, force14i, force21i, force12i, force32i, force23i, force43i, force34i,
                    kick_force1i, kick_force2i, kick_force3i, kick_force4i,
                    shear12i, shear23i, shear34i, shear41i
                ] = vals2
                for val in vals2:
                    assert val == ' 0.0', vals2
                #complex_cshear_force_f06
                #'                  ====== POINT  1 ======      ====== POINT  2 ======      ====== POINT  3 ======      ====== POINT  4 ======'
                #' ELEMENT          F-FROM-4      F-FROM-2      F-FROM-1      F-FROM-3      F-FROM-2      F-FROM-4      F-FROM-3      F-FROM-1'
                #'         ID               KICK-1       SHEAR-12       KICK-2       SHEAR-23       KICK-3       SHEAR-34       KICK-4       SHEAR-41'
                #'            25  0.0           0.0           0.0           0.0           0.0           0.0           0.0           0.0'
                #'                0.0           0.0           0.0           0.0           0.0           0.0           0.0           0.0'

                #f.write('      %8i   %-13s %-13s %-13s / %s\n' % (eid,
                                                                  #force41r, 0force42r, force21r, kick
                                                                  #force41i, ))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class ComplexSpringDamperForceArray(ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        if is_sort1:
            pass
        else:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        headers = ['spring_force']
        return headers

    #def get_headers(self):
        #headers = ['axial', 'torque']
        #return headers

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

        #[axial_force, torque]
        self.data = zeros((self.ntimes, self.ntotal, 1), dtype='complex64')

    def build_dataframe(self):
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
        self.data_frame.columns.names = column_names
        self.data_frame.index.names = ['ElementID', 'Item']

    def __eq__(self, table):
        self._eq_header(table)
        assert self.is_sort1() == table.is_sort1()
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, eid in enumerate(self.element):
                    t1 = self.data[itime, ie, 0]
                    t2 = table.data[itime, ie, 0]

                    if not allclose([t1.real, t1.imag], [t2.real, t2.imag], atol=0.0001):
                        msg += '%s    (%s, %s)  (%s, %s)\n' % (
                            eid,
                            t1.real, t1.imag,
                            t2.real, t2.imag)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, force):
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, 0] = force
        self.ielement += 1

    def get_stats(self):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        #ntotal = self.ntotal

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
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True, is_sort1=True):
        # 11-CELAS1, 12-CELAS2, 13-CELAS3, 14-CELAS4
        if self.element_type == 11:
            msg = ['                         C O M P L E X   F O R C E S   I N   S C A L A R   S P R I N G S   ( C E L A S 1 )\n']
        elif self.element_type == 12:
            msg = ['                         C O M P L E X   F O R C E S   I N   S C A L A R   S P R I N G S   ( C E L A S 2 )\n']
        elif self.element_type == 13:
            msg = ['                         C O M P L E X   F O R C E S   I N   S C A L A R   S P R I N G S   ( C E L A S 3 )\n']
        elif self.element_type == 14:
            msg = ['                         C O M P L E X   F O R C E S   I N   S C A L A R   S P R I N G S   ( C E L A S 4 )\n']
        elif self.element_type == 20: # CDAMP1
            msg = ['                         C O M P L E X   F O R C E S   I N   S C A L A R   D A M P E R S   ( C D A M P 1 )\n']
        elif self.element_type == 21: # CDAMP2
            msg = ['                         C O M P L E X   F O R C E S   I N   S C A L A R   D A M P E R S   ( C D A M P 2 )\n']
        else:
            raise NotImplementedError('element_name=%s element_type=%s' % (self.element_name, self.element_type))

        if is_mag_phase:
            msg += ['                                                          (MAGNITUDE/PHASE)\n \n']
        else:
            msg += ['                                                          (REAL/IMAGINARY)\n \n']

        if is_sort1:
            msg += [
                '                ELEMENT                                                   ELEMENT\n'
                '                  ID.                    FORCE                              ID.                    FORCE\n'
            ]
            #'                      14                  0.0          /  0.0                           0.0          /  0.0'
        else:
            msg += ['            FREQUENCY                    FORCE                        FREQUENCY                    FORCE\n']


        return msg

    #def get_element_index(self, eids):
        ## elements are always sorted; nodes are not
        #itot = searchsorted(eids, self.element)  #[0]
        #return itot

    #def eid_to_element_node_index(self, eids):
        ##ind = ravel([searchsorted(self.element == eid) for eid in eids])
        #ind = searchsorted(eids, self.element)
        ##ind = ind.reshape(ind.size)
        ##ind.sort()
        #return ind

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        # write the f06
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element
        is_odd = False
        nwrite = len(eids)
        if len(eids) % 2 == 1:
            nwrite -= 1
            is_odd = True

        #print('len(eids)=%s nwrite=%s is_odd=%s' % (len(eids), nwrite, is_odd))
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            spring_force = self.data[itime, :, 0]

            for eid, spring_forcei in zip(eids, spring_force):
                [rspring, ispring] = write_imag_floats_13e([spring_forcei], is_mag_phase)
                #ELEMENT                             AXIAL                                       TORSIONAL
                    #ID.                              STRESS                                         STRESS
                    #14                  0.0          /  0.0                           0.0          /  0.0


                f.write('      %8i   %-13s / %-13s\n' % (eid, rspring, ispring))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

class ComplexSpringForceArray(ComplexSpringDamperForceArray):  # 11-CELAS1,12-CELAS2,13-CELAS3, 14-CELAS4
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexSpringDamperForceArray.__init__(self, data_code, is_sort1, isubcase, dt)

class ComplexDamperForceArray(ComplexSpringDamperForceArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexSpringDamperForceArray.__init__(self, data_code, is_sort1, isubcase, dt)


class ComplexViscForceArray(ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        if is_sort1:
            pass
        else:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        headers = ['axial_force', 'torque']
        return headers

    #def get_headers(self):
        #headers = ['axial', 'torque']
        #return headers

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

        #[axial_force, torque]
        self.data = zeros((self.ntimes, self.ntotal, 2), dtype='complex64')

    def __eq__(self, table):
        self._eq_header(table)
        assert self.is_sort1() == table.is_sort1()
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, eid in enumerate(self.element):
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (axial1, torque1) = t1
                    (axial2, torque2) = t2

                    if not allclose(t1, t2):
                        msg += '(%s, %s)    (%s, %s)  (%s, %s)\n' % (
                            eid, nid,
                            axial1, torque1,
                            axial2, torque2)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, axial, torque):
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [axial, torque]
        self.ielement += 1

    def get_stats(self):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        #ntotal = self.ntotal

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
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True, is_sort1=True):
        #if self.element_type == 1: # CROD
            #msg = ['                             C O M P L E X   F O R C E S   I N   R O D   E L E M E N T S   ( C R O D )\n']
        #elif self.element_type == 10: # CONROD
            #msg = ['                           C O M P L E X   F O R C E S   I N   R O D   E L E M E N T S   ( C O N R O D )\n']
        #elif self.element_type == 3: # CTUBE
            #msg = ['                            C O M P L E X   F O R C E S   I N   R O D   E L E M E N T S   ( C T U B E )\n']
            ##pass
        if self.element_type == 24:
            msg = ['                           C O M P L E X   F O R C E S   I N   V I S C   E L E M E N T S   ( C V I S C )\n']
        else:
            raise NotImplementedError('element_name=%s element_type=%s' % (self.element_name, self.element_type))
        if is_mag_phase:
            msg += ['                                                          (MAGNITUDE/PHASE)\n']
        else:
            msg += ['                                                          (REAL/IMAGINARY)\n']

        if is_sort1:
            msg += [
                ' \n'
                '                 ELEMENT                             AXIAL                                       TORSIONAL\n'
                '                   ID.                              STRAIN                                         STRAIN\n'
            ]
            #'                      14                  0.0          /  0.0                           0.0          /  0.0'
        else:
            raise NotImplementedError('sort2')


        return self.element_name, msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element)  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.element == eid) for eid in eids])
        ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        (elem_name, msg_temp) = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        # write the f06
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element
        is_odd = False
        nwrite = len(eids)
        if len(eids) % 2 == 1:
            nwrite -= 1
            is_odd = True

        #print('len(eids)=%s nwrite=%s is_odd=%s' % (len(eids), nwrite, is_odd))
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            axial = self.data[itime, :, 0]
            torsion = self.data[itime, :, 1]

            for eid, axiali, torsioni in zip(eids, axial, torsion):
                [raxial, rtorsion, iaxial, itorsion] = write_imag_floats_13e([axiali, torsioni], is_mag_phase)
                #ELEMENT                             AXIAL                                       TORSIONAL
                    #ID.                              STRESS                                         STRESS
                    #14                  0.0          /  0.0                           0.0          /  0.0
                f.write('      %8i   %-13s / %-13s   %-13s / %s\n' % (eid, raxial, iaxial, rtorsion, itorsion))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class ComplexPlateForceArray(ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        if is_sort1:
            pass
        else:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        headers = ['mx', 'my', 'mxy', 'bmx', 'bmy', 'bmxy', 'tx', 'ty']
        return headers

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

        #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
        self.data = zeros((self.ntimes, self.ntotal, 8), dtype='complex64')

    def build_dataframe(self):
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
        self.data_frame.columns.names = column_names
        self.data_frame.index.names = ['ElementID', 'Item']

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()
        self._eq_header(table)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, eid in enumerate(self.element):
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (mx1, my1, mxy1, bmx1, bmy1, bmxy1, tx1, ty1) = t1
                    (mx2, my2, mxy2, bmx2, bmy2, bmxy2, tx2, ty2) = t2

                    if not allclose(t1, t2):
                    #if not np.array_equal(t1.real, t2.real):
                        msg += ('%-8s (%s, %s, %s, %s, %s, %s, %s, %s)\n'
                                '%-8s (%s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                            eid,
                            #mx1.real, my1.real, mxy1.real, bmx1.real, bmy1.real, bmxy1.real, tx1.real, ty1.real,
                            mx1, my1, mxy1, bmx1, bmy1, bmxy1, tx1, ty1,
                            '',
                            mx2, my2, mxy2, bmx2, bmy2, bmxy2, tx2, ty2,
                            #mx2.real, my2.real, mxy2.real, bmx2.real, bmy2.real, bmxy2.real, tx2.real, ty2.real,
                        ))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [mx, my, mxy, bmx, bmy, bmxy, tx, ty]
        self.ielement += 1
        self.itotal += 1

    def get_stats(self):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        #ntotal = self.ntotal

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
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True, is_sort1=True):
        loads = ['    ELEMENT                - MEMBRANE  FORCES -                        - BENDING MOMENTS -               - TRANSVERSE SHEAR FORCES -\n'
                 '      ID              FX            FY            FXY             MX            MY            MXY             QX            QY\n',]
        if is_mag_phase:
            mag_real = ['                                                         (MAGNITUDE/PHASE)\n \n']
        else:
            mag_real = ['                                                          (REAL/IMAGINARY)\n \n']

        cquad4_bilinear = ['                  C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n']  # good
        cquad4_linear = ['                  C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n']  # good
        ctria3 = ['                     C O M P L E X   F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n']  # good
        cquad8 = ['                  C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n']
        cquadr = ['                  C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n']
        ctria6 = ['                     C O M P L E X   F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n']
        ctriar = ['                     C O M P L E X   F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n']

        is_bilinear = False
        if self.element_type == 144: # CQUAD4
            msg = cquad4_linear + mag_real + loads
        elif self.element_type == 33: # CQUAD4
            msg = cquad4_bilinear + mag_real + loads
        elif self.element_type == 64:  #CQUAD8
            msg = cquad8 + mag_real + loads
        elif self.element_type == 82:  # CQUADR
            msg = cquadr + mag_real + loads
        elif self.element_type == 74: # CTRIA3
            msg = ctria3 + mag_real + loads
        elif self.element_type == 75:  # CTRIA6
            msg = ctria6 + mag_real + loads
        elif self.element_type == 70:  # CTRIAR
            msg = ctriar + mag_real + loads
        else:
            raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))
        return msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element)  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.element == eid) for eid in eids])
        ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        # write the f06
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            mx = self.data[itime, :, 0]
            my = self.data[itime, :, 1]
            mxy = self.data[itime, :, 2]
            bmx = self.data[itime, :, 3]
            bmy = self.data[itime, :, 4]
            bmxy = self.data[itime, :, 5]
            tx = self.data[itime, :, 6]
            ty = self.data[itime, :, 7]

            for eid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi in zip(eids, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
                [smxr, smyr, smxyr, sbmxr, sbmyr, sbmxyr, stxr, styr,
                 smxi, smyi, smxyi, sbmxi, sbmyi, sbmxyi, stxi, styi] = write_imag_floats_13e(
                      [mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi], is_mag_phase)
                #"""
                    #ELEMENT                - MEMBRANE  FORCES -                        - BENDING MOMENTS -               - TRANSVERSE SHEAR FORCES -
                      #ID              FX            FY            FXY             MX            MY            MXY             QX            QY
                #0       564       1.543439E+03  7.311177E+02  1.322702E+02    1.080178E+00  1.699104E+00  2.618547E-01    3.877034E+01  4.518554E+00
                                  #358.3129      358.0245      177.5593        177.5292      178.2112        0.0907        358.1465      179.4567
                #"""
                #                fx     fy     fxy     mx     my     mxy    qx      qy
                f.write('0  %8i   %-13s  %-13s  %-13s   %-13s  %-13s  %-13s   %-13s  %s\n'
                        '   %8s   %-13s  %-13s  %-13s   %-13s  %-13s  %-13s   %-13s  %s\n' % (
                        eid, smxr, smyr, smxyr, sbmxr, sbmyr, sbmxyr, stxr, styr,
                        '', smxi, smyi, smxyi, sbmxi, sbmyi, sbmxyi, stxi, styi))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class ComplexPlate2ForceArray(ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        if is_sort1:
            pass
        else:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        headers = ['mx', 'my', 'mxy', 'bmx', 'bmy', 'bmxy', 'tx', 'ty']
        return headers

    def build(self):
        if self.is_built:
            return

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        self.nelements //= self.ntimes
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
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
        self.element_node = zeros((self.ntotal, 2), dtype='int32')

        #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
        self.data = zeros((self.ntimes, self.ntotal, 8), dtype='complex64')

    def build_dataframe(self):
        headers = self.get_headers()
        assert 0 not in self.element
        #print(self.element_node)
        element_node = [self.element_node[:, 0], self.element_node[:, 1]]
        assert 0 not in self.element_node[:, 0]
        if self.nonlinear_factor is not None:
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=element_node, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
        else:
            self.data_frame = pd.Panel(self.data, major_axis=element_node, minor_axis=headers).to_frame()
            self.data_frame.columns.names = ['Static']
        self.data_frame.index.names = ['ElementID', 'NodeID', 'Item']

    def __eq__(self, table):
        self._eq_header(table)
        assert self.is_sort1() == table.is_sort1()
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, e in enumerate(self.element_node):
                    (eid, nid) = e
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (mx1, my1, mxy1, bmx1, bmy1, bmxy1, tx1, ty1) = t1
                    (mx2, my2, mxy2, bmx2, bmy2, bmxy2, tx2, ty2) = t2

                    if not allclose(t1, t2):
                        base1 = '(%s, %s)   ' % (eid, nid)
                        base2 = ' ' * len(base1)
                        msg += '%s (%s, %s, %s, %s, %s, %s, %s, %s)\n%s(%s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                            base1,
                            mx1, my1, mxy1, bmx1, bmy1, bmxy1, tx1, ty1,
                            base2,
                            mx2, my2, mxy2, bmx2, bmy2, bmxy2, tx2, ty2)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_new_element_sort1(self, dt, eid, term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.element_node[self.itotal, :] = [eid, nid]
        self.data[self.itime, self.itotal, :] = [mx, my, mxy, bmx, bmy, bmxy, tx, ty]
        self.itotal += 1
        self.ielement += 1

    def add_sort1(self, dt, eid, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        self._times[self.itime] = dt
        #assert self.element[self.ielement - 1] == eid, eid
        self.element_node[self.itotal, :] = [eid, nid]
        self.data[self.itime, self.itotal, :] = [mx, my, mxy, bmx, bmy, bmxy, tx, ty]
        self.itotal += 1

    def get_stats(self):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        #ntotal = self.ntotal

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
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True, is_sort1=True):
        loads = ['    ELEMENT                - MEMBRANE  FORCES -                        - BENDING MOMENTS -               - TRANSVERSE SHEAR FORCES -\n'
                 '      ID              FX            FY            FXY             MX            MY            MXY             QX            QY\n',]
        if is_mag_phase:
            mag_real = ['                                                         (MAGNITUDE/PHASE)\n \n']
        else:
            mag_real = ['                                                          (REAL/IMAGINARY)\n \n']

        cquad4_bilinear = ['                  C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n']  # good
        cquad4_linear = ['                  C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n']  # good
        ctria3 = ['                     C O M P L E X   F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n']  # good
        cquad8 = ['                  C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n']
        cquadr = ['                  C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n']
        ctria6 = ['                     C O M P L E X   F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n']
        ctriar = ['                     C O M P L E X   F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n']

        is_bilinear = False
        if self.element_type == 144: # CQUAD4
            msg = cquad4_linear + mag_real + loads
        elif self.element_type == 33: # CQUAD4
            msg = cquad4_bilinear + mag_real + loads
        elif self.element_type == 64:  #CQUAD8
            msg = cquad8 + mag_real + loads
        elif self.element_type == 82:  # CQUADR
            msg = cquadr + mag_real + loads
        elif self.element_type == 74: # CTRIA3
            msg = ctria3 + mag_real + loads
        elif self.element_type == 75:  # CTRIA6
            msg = ctria6 + mag_real + loads
        elif self.element_type == 70:  # CTRIAR
            msg = ctriar + mag_real + loads
        else:
            raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))
        return msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element)  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.element == eid) for eid in eids])
        ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        # write the f06
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            mx = self.data[itime, :, 0]
            my = self.data[itime, :, 1]
            mxy = self.data[itime, :, 2]
            bmx = self.data[itime, :, 3]
            bmy = self.data[itime, :, 4]
            bmxy = self.data[itime, :, 5]
            tx = self.data[itime, :, 6]
            ty = self.data[itime, :, 7]

            for eid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi in zip(eids, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
                [smxr, smyr, smxyr, sbmxr, sbmyr, sbmxyr, stxr, styr,
                 smxi, smyi, smxyi, sbmxi, sbmyi, sbmxyi, stxi, styi] = write_imag_floats_13e(
                      [mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi], is_mag_phase)
                #"""
                    #ELEMENT                - MEMBRANE  FORCES -                        - BENDING MOMENTS -               - TRANSVERSE SHEAR FORCES -
                      #ID              FX            FY            FXY             MX            MY            MXY             QX            QY
                #0       564       1.543439E+03  7.311177E+02  1.322702E+02    1.080178E+00  1.699104E+00  2.618547E-01    3.877034E+01  4.518554E+00
                                  #358.3129      358.0245      177.5593        177.5292      178.2112        0.0907        358.1465      179.4567
                #"""
                #                fx     fy     fxy     mx     my     mxy    qx      qy
                f.write('0  %8i   %-13s  %-13s  %-13s   %-13s  %-13s  %-13s   %-13s  %s\n'
                        '   %8s   %-13s  %-13s  %-13s   %-13s  %-13s  %-13s   %-13s  %s\n' % (
                        eid, smxr, smyr, smxyr, sbmxr, sbmyr, sbmxyr, stxr, styr,
                        '', smxi, smyi, smxyi, sbmxi, sbmyi, sbmxyi, stxi, styi))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class ComplexCBarForceArray(ScalarObject):
    def get_headers(self):
        headers = ['bending_moment_1a', 'bending_moment_2a',
                   'bending_moment_1b', 'bending_moment_2b',
                   'shear1', 'shear2', 'axial', 'torque', ]
        return headers

    def __init__(self, data_code, is_sort1, isubcase, dt):
        #ForceObject.__init__(self, data_code, isubcase)
        ScalarObject.__init__(self, data_code, isubcase)

        #self.eType = {}
        self.result_flag = 0
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.itime = 0
        self.nelements = 0  # result specific
        self.element_type = 'CBAR'
        #self.cid = {}  # gridGauss

        if is_sort1:
            #sort1
            pass
        else:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def is_real(self):
        return False

    def is_complex(self):
        return True

    def build(self):
        #print('ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (self.ntimes, self.nelements, self.ntotal, self.subtitle))
        if self.is_built:
            return
        nnodes = 1

        #self.names = []
        #self.nelements //= nnodes
        self.nelements //= self.ntimes
        #self.ntotal //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        self.is_built = True
        #print('ntotal=%s ntimes=%s nelements=%s' % (self.ntotal, self.ntimes, self.nelements))

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        self._times = zeros(self.ntimes, 'float32')
        self.element = zeros(self.ntotal, 'int32')

        # the number is messed up because of the offset for the element's properties

        if not self.nelements * nnodes == self.ntotal:
            msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (self.ntimes,
                                                                           self.nelements, nnodes,
                                                                           self.nelements * nnodes,
                                                                           self.ntotal)
            raise RuntimeError(msg)
        #[bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq]
        self.data = zeros((self.ntimes, self.ntotal, 8), 'complex64')


    def build_dataframe(self):
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
        self.data_frame.columns.names = column_names
        self.data_frame.index.names = ['ElementID', 'Item']

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()
        self._eq_header(table)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            i = 0
            if self.is_sort1():
                for itime in range(ntimes):
                    for ieid, eid in enumerate(self.element):
                        t1 = self.data[itime, ieid, :]
                        t2 = table.data[itime, ieid, :]
                        (s1a1, s2a1, s3a1, s4a1, axial1, s2a1, s2b1, s2c1, s2d1) = t1
                        (s1a2, s2a2, s3a2, s4a2, axial2, s2a2, s2b2, s2c2, s2d2) = t2
                        #d = t1 - t2
                        if not allclose([s1a1.real, s2a1.real, s3a1.real, s4a1.real, axial1.real, s2a1.real, s2b1.real, s2c1.real, s2d1.real],
                                        [s1a2.real, s2a2.real, s3a2.real, s4a2.real, axial2.real, s2a2.real, s2b2.real, s2c2.real, s2d2.real], atol=0.0001):
                        #if not np.array_equal(t1, t2):
                            msg += '%-4s  (%s, %s, %s, %s, %s, %s, %s, %s, %s)\n      (%s, %s, %s, %s, %s, %s, %s, %s, %s)\m' % (
                                eid,
                                s1a1.real, s2a1.real, s3a1.real, s4a1.real, axial1.real, s2a1.real, s2b1.real, s2c1.real, s2d1.real,
                                s1a2.real, s2a2.real, s3a2.real, s4a2.real, axial2.real, s2a2.real, s2b2.real, s2c2.real, s2d2.real,
                                )
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

    def add_sort1(self, dt, eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq):
        self._times[self.itime] = dt
        self.data[self.itime, self.itotal, :] = [bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq]
        self.element[self.itotal] = eid
        self.itotal += 1

    def get_stats(self):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        msg = []

        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%i\n' % (self.__class__.__name__, nelements))
        msg.append('  eType, cid\n')
        msg.append('  data: [ntimes, nelements, 8] where 8=[%s]\n' % str(', '.join(self.get_headers())))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  is_sort1=%s is_sort2=%s\n' % (self.is_sort1(), self.is_sort2()))
        msg.append('  CBAR\n  ')
        msg += self.get_data_code()
        return msg

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        #msg_temp, nnodes = get_f06_header(self, is_mag_phase, is_sort1)


        #is_sort1 = False
        if is_mag_phase:
            mag_phase = '                                                          (MAGNITUDE/PHASE)\n \n'
        else:
            mag_phase = '                                                          (REAL/IMAGINARY)\n \n'


        name = self.data_code['name']
        if name == 'freq':
            name = 'FREQUENCY'
        #else: # mode
            #raise RuntimeError(name)

        if is_sort1:
            line1 = '0    ELEMENT         BEND-MOMENT-END-A            BEND-MOMENT-END-B                  SHEAR\n'
            line2 = '       ID.         PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE\n'
        else:
            line1 = '                    BEND-MOMENT-END-A            BEND-MOMENT-END-B                  SHEAR\n'
            line2 = '   %16s       PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE\n' % name

        # force
        msg_temp = header + [
            '                             C O M P L E X   F O R C E S   I N   B A R   E L E M E N T S   ( C B A R )\n',
            mag_phase,
            ' ',
            line1,
            line2,
        ]
        if self.is_sort1():
            assert self.is_sort1() == True, str(self)
            if is_sort1:
                page_num = self._write_sort1_as_sort1(f, page_num, page_stamp, header, msg_temp, is_mag_phase)
            else:
                self._write_sort1_as_sort2(f, page_num, page_stamp, header, msg_temp, is_mag_phase)
        else:
            assert self.is_sort1() == True, str(self)
        return page_num - 1

    def _write_sort1_as_sort1(self, f, page_num, page_stamp, header, msg_temp, is_mag_phase):
        eids = self.element
        times = self._times
        ntimes = self.data.shape[0]
        for itime in range(ntimes):
            dt = self._times[itime]
            dt_line = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
            header[1] = dt_line
            msg = header + msg_temp
            f.write(''.join(msg))

            #bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq
            assert self.is_sort1() == True, str(self)
            bm1a = self.data[itime, :, 0]
            bm2a = self.data[itime, :, 1]
            bm1b = self.data[itime, :, 2]
            bm2b = self.data[itime, :, 3]
            ts1 = self.data[itime, :, 4]
            ts2 = self.data[itime, :, 5]
            af = self.data[itime, :, 6]
            trq = self.data[itime, :, 7]

            for eid, bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi in zip(eids, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq):
                vals = (bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi)
                vals2 = write_imag_floats_13e(vals, is_mag_phase)
                (bm1air, bm2air, bm1bir, bm2bir, ts1ir, ts2ir, afir, trqir,
                 bm1aii, bm2aii, bm1bii, bm2bii, ts1ii, ts2ii, afii, trqii) = vals2

                f.write('0%16i   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                        ' %14s   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                            eid, bm1air, bm2air, bm1bir, bm2bir, ts1ir, ts2ir, afir, trqir,
                            '', bm1aii, bm2aii, bm1bii, bm2bii, ts1ii, ts2ii, afii, trqii))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num

    def _write_sort1_as_sort2(self, f, page_num, page_stamp, header, msg_temp, is_mag_phase):
        eids = self.element
        times = self._times
        ntimes = self.data.shape[0]
        for ieid, eid in enumerate(eids):
            eid_line = ' ELEMENT-ID = %s' % (eid)
            header[1] = eid_line
            msg = header + msg_temp
            f.write(''.join(msg))

            #bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq
            bm1a = self.data[:, ieid, 0]
            bm2a = self.data[:, ieid, 1]
            bm1b = self.data[:, ieid, 2]
            bm2b = self.data[:, ieid, 3]
            ts1 = self.data[:, ieid, 4]
            ts2 = self.data[:, ieid, 5]
            af = self.data[:, ieid, 6]
            trq = self.data[:, ieid, 7]

            for dt, bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi in zip(times, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq):
                vals = (bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi)
                vals2 = write_imag_floats_13e(vals, is_mag_phase)
                (bm1air, bm2air, bm1bir, bm2bir, ts1ir, ts2ir, afir, trqir,
                 bm1aii, bm2aii, bm1bii, bm2bii, ts1ii, ts2ii, afii, trqii) = vals2

                f.write('0%16s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                        ' %15s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                            write_float_12E(dt),
                            bm1air, bm2air, bm1bir, bm2bir, ts1ir, ts2ir, afir, trqir,
                            '', bm1aii, bm2aii, bm1bii, bm2bii, ts1ii, ts2ii, afii, trqii))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num


class ComplexCBeamForceArray(ScalarObject):
    def get_headers(self):
        headers = [
            'sd', 'bending_moment1', 'bending_moment2', 'shear1', 'shear2',
            'axial_force', 'total_torque', 'warping_torque', ]
        return headers

    def __init__(self, data_code, is_sort1, isubcase, dt):
        #ForceObject.__init__(self, data_code, isubcase)
        ScalarObject.__init__(self, data_code, isubcase)

        self.result_flag = 0
        self.itime = 0
        self.nelements = 0  # result specific
        self.element_type = 'CBEAM'

        if is_sort1:
            #sort1
            pass
        else:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def is_real(self):
        return False

    def is_complex(self):
        return True

    def build(self):
        #print('ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (self.ntimes, self.nelements, self.ntotal, self.subtitle))
        if self.is_built:
            return
        nnodes = 11

        #self.names = []
        #self.nelements //= nnodes
        self.nelements //= self.ntimes
        #self.ntotal //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        self.is_built = True
        #print('ntotal=%s ntimes=%s nelements=%s' % (self.ntotal, self.ntimes, self.nelements))

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        self._times = zeros(self.ntimes, 'float32')
        self.element = zeros(self.ntotal, 'int32')
        self.element_node = zeros((self.ntotal, 2), 'int32')

        # the number is messed up because of the offset for the element's properties

        if not self.nelements * nnodes == self.ntotal:
            msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (self.ntimes,
                                                                           self.nelements, nnodes,
                                                                           self.nelements * nnodes,
                                                                           self.ntotal)
            raise RuntimeError(msg)
        #[sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq]
        self.data = zeros((self.ntimes, self.ntotal, 8), 'complex64')

    def finalize(self):
        sd = self.data[0, :, 0].real
        i_sd_zero = np.where(sd != 0.0)[0]
        i_node_zero = np.where(self.element_node[:, 1] != 0)[0]
        assert i_node_zero.max() > 0, 'CBEAM element_node hasnt been filled'
        i = np.union1d(i_sd_zero, i_node_zero)
        self.element = self.element[i]
        self.element_node = self.element_node[i, :]
        self.data = self.data[:, i, :]

    def build_dataframe(self):
        headers = self.get_headers()[1:]
        column_names, column_values = self._build_dataframe_transient_header()
        element_location = [
            self.element_node[:, 0],
            self.data[0, :, 0].real,
        ]
        self.data_frame = pd.Panel(self.data[:, :, 1:], items=column_values, major_axis=element_location, minor_axis=headers).to_frame()
        self.data_frame.columns.names = column_names
        self.data_frame.index.names = ['ElementID', 'Location', 'Item']
        #print(self.data_frame)

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()
        self._eq_header(table)
        if not np.array_equal(self.data, table.data):
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
                        #if not np.array_equal(t1, t2):
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

    #def add_new_element_sort1(self, dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq):
        #return self.add_sort1(dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq)

    def add_sort1(self, dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq):
        self._times[self.itime] = dt
        self.data[self.itime, self.itotal, :] = [sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq]
        self.element[self.itotal] = eid
        self.element_node[self.itotal, :] = [eid, nid]
        self.itotal += 1

    def get_stats(self):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        msg = []

        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%i\n' % (self.__class__.__name__, nelements))
        #msg.append('  eType, cid\n')
        msg.append('  data: [ntimes, nelements, 8] where 8=[%s]\n' % str(', '.join(self.get_headers())))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  is_sort1=%s is_sort2=%s\n' % (self.is_sort1(), self.is_sort2()))
        msg.append('  CBEAM\n  ')
        msg += self.get_data_code()
        return msg

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        # option B
        #'                           C O M P L E X   F O R C E S   I N   B E A M   E L E M E N T S   ( C B E A M ) '
        #'                                                          (REAL/IMAGINARY)'
        #'                    STAT DIST/   - BENDING MOMENTS -            - WEB  SHEARS -           AXIAL          TOTAL          WARPING'
        #'   ELEMENT-ID  GRID   LENGTH    PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE         TORQUE'
        #'0        20'
        #'0                11   0.000    0.0           0.0            0.0           0.0            0.0            0.0            0.0'
        #'                               0.0           0.0            0.0           0.0            0.0            0.0            0.0'
        #'0                12   1.000    0.0           0.0            0.0           0.0            0.0            0.0            0.0'
        #'                               0.0           0.0            0.0           0.0            0.0            0.0            0.0'

        #msg_temp, nnodes = get_f06_header(self, is_mag_phase, is_sort1)
        #print('write_f06 not implemented for ComplexCBeamForceArray')
        #return page_num
        #asdf

        #is_sort1 = False
        if is_mag_phase:
            mag_phase = '                                                          (MAGNITUDE/PHASE)\n \n'
        else:
            mag_phase = '                                                          (REAL/IMAGINARY)\n \n'


        #name = self.data_code['name']
        #if name == 'freq':
            #name = 'FREQUENCY'
        #else: # mode
            #raise RuntimeError(name)

        if is_sort1:
            line1 = '0    ELEMENT         BEND-MOMENT-END-A            BEND-MOMENT-END-B                  SHEAR\n'
            line2 = '       ID.         PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE\n'
        else:
            line1 = '                    BEND-MOMENT-END-A            BEND-MOMENT-END-B                  SHEAR\n'
            line2 = '   %16s       PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE\n' % name

        # force
        msg_temp = header + [
            '                             C O M P L E X   F O R C E S   I N   B A R   E L E M E N T S   ( C B E A M )\n',
            mag_phase,
            ' ',
            line1,
            line2,
        ]
        if self.is_sort1():
            assert self.is_sort1() == True, str(self)
            #if is_sort1:
            page_num = self._write_sort1_as_sort1(f, page_num, page_stamp, header, msg_temp, is_mag_phase)
            #else:
                #self._write_sort1_as_sort2(f, page_num, page_stamp, header, msg_temp, is_mag_phase)
        else:
            assert self.is_sort1() == True, str(self)
        return page_num - 1

    def _write_sort1_as_sort1(self, f, page_num, page_stamp, header, msg_temp, is_mag_phase):
        eids = self.element
        times = self._times
        ntimes = self.data.shape[0]
        for itime in range(ntimes):
            dt = self._times[itime]
            dt_line = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
            header[1] = dt_line
            msg = header + msg_temp
            f.write(''.join(msg))

            #bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq
            assert self.is_sort1() == True, str(self)
            #sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq
            sd = self.data[itime, :, 0]
            bm1 = self.data[itime, :, 1]
            bm2 = self.data[itime, :, 2]
            ts1 = self.data[itime, :, 3]
            ts2 = self.data[itime, :, 4]
            af = self.data[itime, :, 5]
            ttrq = self.data[itime, :, 6]
            wtrq = self.data[itime, :, 7]

            for eid, sdi, bm1i, bm2i, ts1i, ts2i, afi, ttrqi, wtrqi in zip(eids, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq):
                vals = (sdi, bm1i, bm2i, ts1i, ts2i, afi, ttrqi, wtrqi)
                vals2 = write_imag_floats_13e(vals, is_mag_phase)
                (sdir, bm1ir, bm2ir, ts1ir, ts2ir, afir, ttrqir, wtrqir,
                 sdii, bm1ii, bm2ii, ts1ii, ts2ii, afii, ttrqii, wtrqii) = vals2

                f.write('0%16i   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                        ' %14s   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                            eid, sdir, bm1ir, bm2ir, ts1ir, ts2ir, afir, ttrqir, wtrqir,
                            '',  sdii, bm1ii, bm2ii, ts1ii, ts2ii, afii, ttrqii, wtrqii))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num


class ComplexCBendForceArray(ScalarObject):  # 69-CBEND
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        if is_sort1:
            pass
        else:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        headers = [
            'bending_moment_1a', 'bending_moment_2a', 'shear_1a', 'shear_2a', 'axial_a', 'torque_a',
            'bending_moment_1b', 'bending_moment_2b', 'shear_1b', 'shear_2b', 'axial_b', 'torque_b',
        ]
        return headers

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
        self.element_nodes = zeros((self.nelements, 3), dtype='int32')

        #[bending_moment_1a, bending_moment_2a, shear_1a, shear_2a, axial_a, torque_a
        # bending_moment_1b, bending_moment_2b, shear_1b, shear_2b, axial_b, torque_b]
        self.data = zeros((self.ntimes, self.nelements, 12), dtype='complex64')

    def build_dataframe(self):
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        element = self.element_nodes[:, 0]
        self.data_frame = pd.Panel(self.data, items=column_values, major_axis=element, minor_axis=headers).to_frame()
        self.data_frame.columns.names = column_names
        self.data_frame.index.names = ['ElementID', 'Item']

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()
        if not np.array_equal(self.element_nodes, table.element_nodes):
            assert self.element_nodes.shape == table.element_nodes.shape, 'element_nodes shape=%s table.shape=%s' % (self.element_nodes.shape, table.element_nodes.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += 'Eid, Nid_A, Nid_B\n'
            for (eid1, nida1, nidb1), (eid2, nida2, nidb2) in zip(self.element_nodes, table.element_nodes):
                msg += '(%s, %s, %s), (%s, %s, %s)\n' % (eid1, nida1, nidb1, eid2, nida2, nidb2)
            print(msg)
            raise ValueError(msg)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            eids = self.element_nodes[:, 0]
            for itime in range(self.ntimes):
                for ie, eid in enumerate(eids):
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (bending_moment_1a1, bending_moment_2a1, shear_1a1, shear_2a1, axial_a1, torque_a1,
                     bending_moment_1b1, bending_moment_2b1, shear_1b1, shear_2b1, axial_b1, torque_b1) = t1
                    (bending_moment_1a2, bending_moment_2a2, shear_1a2, shear_2a2, axial_a2, torque_a2,
                     bending_moment_1b2, bending_moment_2b2, shear_1b2, shear_2b2, axial_b2, torque_b2) = t2

                    if not allclose(t1, t2):
                        msg += '(%s)    (%s, %s)  (%s, %s)\n' % (
                            eid,
                            bending_moment_1a1.real,
                            bending_moment_1b1.real,

                            bending_moment_1a2.real,
                            bending_moment_1b2.real, )
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)

                    #if not allclose(t1, t2):
                        #msg += '(%s)    (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)  (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                            #eid,
                            #bending_moment_1a1, bending_moment_2a1, shear_1a1, shear_2a1, axial_a1, torque_a1,
                            #bending_moment_1b1, bending_moment_2b1, shear_1b1, shear_2b1, axial_b1, torque_b1,

                            #bending_moment_1a2, bending_moment_2a2, shear_1a2, shear_2a2, axial_a2, torque_a2,
                            #bending_moment_1b2, bending_moment_2b2, shear_1b2, shear_2b2, axial_b2, torque_b2)
                        #i += 1
                        #if i > 10:
                            #print(msg)
                            #raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid,
                  nid_a, bending_moment_1a, bending_moment_2a, shear_1a, shear_2a, axial_a, torque_a,
                  nid_b, bending_moment_1b, bending_moment_2b, shear_1b, shear_2b, axial_b, torque_b):
        bending_moment_1a, bending_moment_2a, shear_1a, shear_2a, axial_a, torque_a,
        bending_moment_1b, bending_moment_2b, shear_1b, shear_2b, axial_b, torque_b

        self._times[self.itime] = dt
        self.element_nodes[self.ielement] = [eid, nid_a, nid_b]
        self.data[self.itime, self.ielement, :] = [
            bending_moment_1a, bending_moment_2a, shear_1a, shear_2a, axial_a, torque_a,
            bending_moment_1b, bending_moment_2b, shear_1b, shear_2b, axial_b, torque_b
        ]
        self.ielement += 1
        if self.ielement == self.nelements:
            self.ielement = 0

    def get_stats(self):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        #ntotal = self.ntotal

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
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True, is_sort1=True):
        msg = ['                           C O M P L E X   F O R C E S   I N   B E N D    E L E M E N T S   ( C B E N D )\n']

        if is_mag_phase:
            msg += ['                                                          (MAGNITUDE/PHASE)\n']
        else:
            msg += ['                                                          (REAL/IMAGINARY)\n']

        if is_sort1:
            msg += [
                '                                 - BENDING MOMENTS -            -   SHEARS   -            AXIAL'
                '   ELEMENT-ID  GRID    END      PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE'
            ]
        else:
            raise NotImplementedError('sort2')
        return msg

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        #'                           C O M P L E X   F O R C E S   I N   B E N D    E L E M E N T S   ( C B E N D )'
        #'                                                          (REAL/IMAGINARY)'
        #'                                 - BENDING MOMENTS -            -   SHEARS   -            AXIAL'
        #'   ELEMENT-ID  GRID    END      PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE'
        #'0        27      21      A     0.0           0.0            0.0           0.0            0.0            0.0'
        #'                               0.0           0.0            0.0           0.0            0.0            0.0'
        #'0                22      B     0.0           0.0            0.0           0.0            0.0            0.0'
        #'                               0.0           0.0            0.0           0.0            0.0            0.0'
        msg_temp = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        # write the f06
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element_nodes[:, 0]
        nid_a = self.element_nodes[:, 1]
        nid_b = self.element_nodes[:, 2]
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            bending_moment_1a = self.data[itime, :, 0]
            bending_moment_2a = self.data[itime, :, 1]
            shear_1a = self.data[itime, :, 2]
            shear_2a = self.data[itime, :, 3]
            axial_a = self.data[itime, :, 4]
            torque_a = self.data[itime, :, 5]
            bending_moment_1b = self.data[itime, :, 6]
            bending_moment_2b = self.data[itime, :, 7]
            shear_1b = self.data[itime, :, 8]
            shear_2b = self.data[itime, :, 9]
            axial_b = self.data[itime, :, 10]
            torque_b = self.data[itime, :, 11]

            for (eid,
                 nid_ai, bending_moment_1ai, bending_moment_2ai, shear_1ai, shear_2ai, axial_ai, torque_ai,
                 nid_bi, bending_moment_1bi, bending_moment_2bi, shear_1bi, shear_2bi, axial_bi, torque_bi) in zip(eids,
                                                                                                                   nid_a, bending_moment_1a, bending_moment_2a, shear_1a, shear_2a, axial_a, torque_a,
                                                                                                                   nid_b, bending_moment_1b, bending_moment_2b, shear_1b, shear_2b, axial_b, torque_b):
                [bending_moment_1ari, bending_moment_2ari, shear_1ari, shear_2ari, axial_ari, torque_ari,
                 bending_moment_1bri, bending_moment_2bri, shear_1bri, shear_2bri, axial_bri, torque_bri,
                 bending_moment_1aii, bending_moment_2aii, shear_1aii, shear_2aii, axial_aii, torque_aii,
                 bending_moment_1bii, bending_moment_2bii, shear_1bii, shear_2bii, axial_bii, torque_bii,
                 ] = write_imag_floats_13e(
                     [bending_moment_1ai, bending_moment_2ai, shear_1ai, shear_2ai, axial_ai, torque_ai,
                      bending_moment_1bi, bending_moment_2bi, shear_1bi, shear_2bi, axial_bi, torque_bi],
                     is_mag_phase)
                f.write(
                    '0  %8s%8s      A    %13s %13s  %13s %13s  %13s  %s\n'
                    '                              %13s %13s  %13s %13s  %13s  %s\n'
                    '0  %8s%8s      B    %13s %13s  %13s %13s  %13s  %s\n'
                    '                              %13s %13s  %13s %13s  %13s  %s\n'
                    % (
                    eid, nid_ai,
                    bending_moment_1ari, bending_moment_2ari, shear_1ari, shear_2ari, axial_ari, torque_ari,
                    bending_moment_1aii, bending_moment_2aii, shear_1aii, shear_2aii, axial_aii, torque_aii,
                    '', nid_bi,
                    bending_moment_1bri, bending_moment_2bri, shear_1bri, shear_2bri, axial_bri, torque_bri,
                    bending_moment_1bii, bending_moment_2bii, shear_1bii, shear_2bii, axial_bii, torque_bii,))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class ComplexSolidPressureForceArray(ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        if is_sort1:
            pass
        else:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        headers = ['ax', 'ay', 'az', 'vx', 'vy', 'vz', 'pressure']
        return headers

    #def get_headers(self):
        #headers = ['axial', 'torque']
        #return headers

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

        #[ax, ay, az, vx, vy, vz, pressure]
        self.data = zeros((self.ntimes, self.ntotal, 7), dtype='complex64')

    def build_dataframe(self):
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
        self.data_frame.columns.names = column_names
        self.data_frame.index.names = ['ElementID', 'Item']

    def __eq__(self, table):
        self._eq_header(table)
        assert self.is_sort1() == table.is_sort1()
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, eid in enumerate(self.element):
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (ax1, ay1, az1, vx1, vy1, vz1, pressure1) = t1
                    (ax2, ay2, az2, vx2, vy2, vz2, pressure2) = t2
                    #rpressure1 = pressure1.real
                    #rpressure2 = pressure2.real
                    if not allclose([ax1, ay1, az1, vx1, vy1, vz1],
                                    [ax2, ay2, az2, vx2, vy2, vz2]):
                        msg += '%s    (%s, %s)  (%s, %s)\n' % (
                            eid,
                            ax1.real, t1.imag,
                            ax2.real, t2.imag)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, ename, ax, ay, az, vx, vy, vz, pressure):
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [ax, ay, az, vx, vy, vz, pressure]
        self.ielement += 1

    def get_stats(self):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        #ntotal = self.ntotal

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
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True, is_sort1=True):
        msg = [
            '                                                 ( R O O T   M E A N   S Q U A R E )      \n'
            '              C O M P L E X   A C C E L E R A T I O N S   V E L O C I T I E S   A N D   P R E S S U R E   L E V E L S\n'
            #'                                                          (REAL/IMAGINARY)'
            #'     ELE-ID   EL-TYPE    X-ACCELERATION  Y-ACCELERATION  Z-ACCELERATION     X-VELOCITY     Y-VELOCITY     Z-VELOCITY   PRESSURE (DB)'
            #'       2000    PENPR      6.883253E+06    1.066544E+07   -6.883253E+06     7.288279E+05  -3.134843E+04  -7.288279E+05   1.162309E+02'
            #'                          1.831744E+07   -7.878719E+05   -1.831744E+07    -2.738759E+05  -4.243642E+05   2.738759E+05'
            #''
        ]
        #msg = ['                                   C O M P L E X   A C O U S T I C   P R E S S U R E   R E S U L T S']
        #'                                   C O M P L E X   A C O U S T I C   P R E S S U R E   R E S U L T S'
        #'                                                         (MAGNITUDE/PHASE)'
        #' '
        #'      POINT ID.               TYPE                    P                   P(RMS)              DB                  DB(A)'
        #'0           57                 S                   7.339671E+05        5.189931E+05        1.173135E+02        3.011353E+01'
        #'                                                       249.9102            249.9102            249.9102            249.9102'

        if is_mag_phase:
            msg += ['                                                          (MAGNITUDE/PHASE)\n \n']
        else:
            msg += ['                                                          (REAL/IMAGINARY)\n \n']

        if is_sort1:
            msg += ['     ELE-ID   EL-TYPE    X-ACCELERATION  Y-ACCELERATION  Z-ACCELERATION     X-VELOCITY     Y-VELOCITY     Z-VELOCITY   PRESSURE (DB)\n']
            #msg += [
                #'      POINT ID.               TYPE                    P                   P(RMS)              DB                  DB(A)\n'
            #]
            #'                      14                  0.0          /  0.0                           0.0          /  0.0'
        else:
            sort2


        return msg

    #def get_element_index(self, eids):
        ## elements are always sorted; nodes are not
        #itot = searchsorted(eids, self.element)  #[0]
        #return itot

    #def eid_to_element_node_index(self, eids):
        ##ind = ravel([searchsorted(self.element == eid) for eid in eids])
        #ind = searchsorted(eids, self.element)
        ##ind = ind.reshape(ind.size)
        ##ind.sort()
        #return ind

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        # write the f06
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element

        #print('len(eids)=%s nwrite=%s is_odd=%s' % (len(eids), nwrite, is_odd))
        etypei = self.element_type
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            ax = self.data[itime, :, 0]
            ay = self.data[itime, :, 0]
            az = self.data[itime, :, 0]
            vx = self.data[itime, :, 0]
            vy = self.data[itime, :, 0]
            vz = self.data[itime, :, 0]
            pressure = self.data[itime, :, 0]

            for eid, axi, ayi, azi, vxi, vyi, vzi, pressurei in zip(eids, ax, ay, az, vx, vy, vz, pressure):
                [saxr, sayr, sazr, svxr, svyr, svzr, spressurer,
                 saxi, sayi, sazi, svxi, svyi, svzi, spressurei] = write_imag_floats_13e([axi, ayi, azi, vxi, vyi, vzi, pressurei], is_mag_phase)
                #'       1000    HEXPR      1.582050E-08    5.505425E+06    2.598164E-09    -8.884337E-10  -4.806934E+04   1.046571E-10   9.968034E+01'
                #'                         -1.116439E-08   -6.040572E+05    1.315160E-09    -1.258955E-09  -4.381078E+05  -2.067553E-10'

                f.write('      %8i %8s %-13s %-13s %-13s %-13s %-13s %-13s %s\n'
                        '      %8s %8s %-13s %-13s %-13s %-13s %-13s %s\n\n'
                        % (eid, etypei, saxr, sayr, sazr, svxr, svyr, svzr, spressurer,
                           '', '',      saxi, sayi, sazi, svxi, svyi, svzi))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class ComplexCBushForceArray(ScalarObject):
    def get_headers(self):
        headers = ['fx', 'fy', 'fz', 'mx', 'my', 'mz']
        return headers

    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)

        self.result_flag = 0
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.itime = 0
        self.nelements = 0  # result specific
        self.element_type = 'CBUSH'

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def build(self):
        #print('ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (self.ntimes, self.nelements, self.ntotal, self.subtitle))
        if self.is_built:
            return
        nnodes = 1

        #self.names = []
        #self.nelements //= nnodes
        self.nelements /= self.ntimes
        #self.ntotal //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        self.is_built = True
        #print('ntotal=%s ntimes=%s nelements=%s' % (self.ntotal, self.ntimes, self.nelements))

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        self._times = zeros(self.ntimes, 'float32')
        self.element = zeros(self.ntotal, 'int32')

        # the number is messed up because of the offset for the element's properties

        if not self.nelements * nnodes == self.ntotal:
            msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (self.ntimes,
                                                                           self.nelements, nnodes,
                                                                           self.nelements * nnodes,
                                                                           self.ntotal)
            raise RuntimeError(msg)
        #[fx, fy, fz, mx, my, mz]
        self.data = zeros((self.ntimes, self.ntotal, 6), 'complex64')

    def build_dataframe(self):
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
        self.data_frame.columns.names = column_names
        self.data_frame.index.names = ['ElementID', 'Item']

    def __eq__(self, table):
        self._eq_header(table)
        assert self.is_sort1() == table.is_sort1()
        if not np.array_equal(self.data, table.data):
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
                        #if not np.array_equal(t1, t2):
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

    def add_sort1(self, dt, eid, fx, fy, fz, mx, my, mz):
        #[fx, fy, fz, mx, my, mz]
        self._times[self.itime] = dt
        self.data[self.itime, self.itotal, :] = [fx, fy, fz, mx, my, mz]
        self.element[self.itotal] = eid
        self.itotal += 1

    def get_stats(self):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        #ntotal = self.ntotal
        msg = []

        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%i\n' % (self.__class__.__name__, nelements))
        msg.append('  eType, cid\n')
        msg.append('  data: [ntimes, nelements, 6] where 6=[%s]\n' % str(', '.join(self.get_headers())))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        # msg.append('  is_sort1=%s is_sort2=%s\n' % (self.is_sort1(), self.is_sort2()))
        msg.append('  CBUSH\n  ')
        msg += self.get_data_code()
        return msg

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        #msg_temp, nnodes = get_f06_header(self, is_mag_phase, is_sort1)

        # write the f06

        #is_sort1 = False
        if is_mag_phase:
            mag_phase = '                                                          (MAGNITUDE/PHASE)\n\n'
        else:
            mag_phase = '                                                          (REAL/IMAGINARY)\n\n'


        name = self.data_code['name']
        if name == 'freq':
            name = 'FREQUENCY'
        else:
            raise RuntimeError(name)

        # is_sort1 = True
        if is_sort1:
            line2 = '       ID.     FORCE-X       FORCE-Y       FORCE-Z      MOMENT-X      MOMENT-Y      MOMENT-Z  \n'
        else:
            line2 = '   %26s        FORCE-X       FORCE-Y       FORCE-Z      MOMENT-X      MOMENT-Y      MOMENT-Z  \n' % name

        # force
        msg_temp = header + [
            '                         C O M P L E X   F O R C E S   I N   B U S H   E L E M E N T S   ( C B U S H ) \n',
            mag_phase,
            ' ',
            # line1,
            line2,
        ]
        if self.is_sort1():
            if is_sort1:
                page_num = self._write_sort1_as_sort1(f, page_num, page_stamp, header, msg_temp, is_mag_phase)
            else:
                page_num = self._write_sort1_as_sort2(f, page_num, page_stamp, header, msg_temp, is_mag_phase)
        else:
            assert self.is_sort1() == True, str(self)
        return page_num - 1

    def _write_sort1_as_sort1(self, f, page_num, page_stamp, header, msg_temp, is_mag_phase):
        ntimes = self.data.shape[0]
        eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]
            dt_line = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
            header[1] = dt_line
            msg = header + msg_temp
            f.write(''.join(msg))

            #fx, fy, fz, mx, my, mz
            if self.is_sort1():
                fx = self.data[itime, :, 0]
                fy = self.data[itime, :, 1]
                fz = self.data[itime, :, 2]
                mx = self.data[itime, :, 3]
                my = self.data[itime, :, 4]
                mz = self.data[itime, :, 5]
            else:
                fx = self.data[:, itime, 0]
                fy = self.data[:, itime, 1]
                fz = self.data[:, itime, 2]
                mx = self.data[:, itime, 3]
                my = self.data[:, itime, 4]
                mz = self.data[:, itime, 5]

            for eid, fxi, fyi, fzi, mxi, myi, mzi in zip(eids, fx, fy, fz, mx, my, mz):
                vals = (fxi, fyi, fzi, mxi, myi, mzi)
                vals2 = write_imag_floats_13e(vals, is_mag_phase)
                (fxir, fyir, fzir, mxir, myir, mzir,
                 fxii, fyii, fzii, mxii, myii, mzii) = vals2

                f.write('0%26i   %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                        ' %26s   %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                            eid, fxir, fyir, fzir, mxir, myir, mzir,
                            '', fxii, fyii, fzii, mxii, myii, mzii))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num

    def _write_sort1_as_sort2(self, f, page_num, page_stamp, header, msg_temp, is_mag_phase):
        eids = self.element
        times = self._times
        for ieid, eid in enumerate(eids):
            eid_line = ' ELEMENT-ID = %s' % (eid)
            header[1] = eid_line
            msg = header + msg_temp
            f.write(''.join(msg))

            if self.is_sort1():
                fx = self.data[:, ieid, 0]
                fy = self.data[:, ieid, 1]
                fz = self.data[:, ieid, 2]
                mx = self.data[:, ieid, 3]
                my = self.data[:, ieid, 4]
                mz = self.data[:, ieid, 5]
            else:
                raise RuntimeError()

            for dt, fxi, fyi, fzi, mxi, myi, mzi in zip(times, fx, fy, fz, mx, my, mz):
                vals = (fxi, fyi, fzi, mxi, myi, mzi)
                vals2 = write_imag_floats_13e(vals, is_mag_phase)
                (fxir, fyir, fzir, mxir, myir, mzir,
                 fxii, fyii, fzii, mxii, myii, mzii) = vals2

                f.write('0%26s   %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                        ' %26s   %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                            write_float_12E(dt),
                            fxir, fyir, fzir, mxir, myir, mzir,
                            '', fxii, fyii, fzii, mxii, myii, mzii))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num


class ComplexForce_VU(ScalarObject):  # 191-VUBEAM
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
        self.parent = {}
        self.coord = {}
        self.icord = {}

        self.forceX = {}
        self.shearY = {}
        self.shearZ = {}
        self.torsion = {}
        self.bendingY = {}
        self.bendingZ = {}

        # TODO if dt=None, handle SORT1 case
        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.add_sort2

    def get_stats(self):
        msg = ['  '] + self.get_data_code()
        nelements = len(self.coord)
        if self.dt is not None:  # transient
            ntimes = len(self.forceX)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  parent, coord, icord, forceX, shearY, shearZ, torsion, '
                   'bendingY, bendingZ\n')
        return msg

    def add_new_transient(self, dt):
        self.dt = dt
        self.forceX[dt] = {}
        self.shearY[dt] = {}
        self.shearZ[dt] = {}
        self.torsion[dt] = {}
        self.bendingY[dt] = {}
        self.bendingZ[dt] = {}

    def add(self, nnodes, dt, data):
        [eid, parent, coord, icord, forces] = data
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord

        self.forceX[eid] = {}
        self.shearY[eid] = {}
        self.shearZ[eid] = {}
        self.torsion[eid] = {}
        self.bendingY[eid] = {}
        self.bendingZ[eid] = {}

        for force in forces:
            [nid, posit, forceX, shearY, shearZ, torsion,
                bendingY, bendingZ] = force
            self.forceX[eid][nid] = forceX
            self.shearY[eid][nid] = shearY
            self.shearZ[eid][nid] = shearZ
            self.torsion[eid][nid] = torsion
            self.bendingY[eid][nid] = bendingY
            self.bendingZ[eid][nid] = bendingZ

    def add_sort1(self, nnodes, dt, data):
        [eid, parent, coord, icord, forces] = data
        if dt not in self.forceX:
            self.add_new_transient(dt)
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord

        self.forceX[dt][eid] = {}
        self.shearY[dt][eid] = {}
        self.shearZ[dt][eid] = {}
        self.torsion[dt][eid] = {}
        self.bendingY[dt][eid] = {}
        self.bendingZ[dt][eid] = {}

        for force in forces:
            [nid, posit, forceX, shearY, shearZ, torsion,
                bendingY, bendingZ] = force
            self.forceX[dt][eid][nid] = forceX
            self.shearY[dt][eid][nid] = shearY
            self.shearZ[dt][eid][nid] = shearZ
            self.torsion[dt][eid][nid] = torsion
            self.bendingY[dt][eid][nid] = bendingY
            self.bendingZ[dt][eid][nid] = bendingZ

    def add_sort2(self, nnodes, eid, data):
        [dt, parent, coord, icord, forces] = data
        if dt not in self.forceX:
            self.add_new_transient(dt)
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord

        self.forceX[dt][eid] = {}
        self.shearY[dt][eid] = {}
        self.shearZ[dt][eid] = {}
        self.torsion[dt][eid] = {}
        self.bendingY[dt][eid] = {}
        self.bendingZ[dt][eid] = {}
        for force in forces:
            [nid, posit, forceX, shearY, shearZ, torsion,
                bendingY, bendingZ] = force
            self.forceX[dt][eid][nid] = forceX
            self.shearY[dt][eid][nid] = shearY
            self.shearZ[dt][eid][nid] = shearZ
            self.torsion[dt][eid][nid] = torsion
            self.bendingY[dt][eid][nid] = bendingY
            self.bendingZ[dt][eid][nid] = bendingZ


class ComplexForce_VU_2D(ScalarObject):  # 189-VUQUAD,190-VUTRIA
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
        self.parent = {}
        self.coord = {}
        self.icord = {}
        self.theta = {}

        self.membraneX = {}
        self.membraneY = {}
        self.membraneXY = {}
        self.bendingX = {}
        self.bendingY = {}
        self.bendingXY = {}
        self.shearYZ = {}
        self.shearXZ = {}

        # TODO if dt=None, handle SORT1 case
        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.add_sort2

    def get_stats(self):
        msg = ['  '] + self.get_data_code()
        nelements = len(self.coord)
        if self.dt is not None:  # transient
            ntimes = len(self.membraneX)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  parent, coord, icord, theta, membraneX, membraneY, '
                   'membraneXY, bendingX, bendingY, bendingXY, '
                   'shearYZ, shearXZ\n')
        return msg

    def add_new_transient(self, dt):
        self.membraneX[dt] = {}
        self.membraneY[dt] = {}
        self.membraneXY[dt] = {}
        self.bendingX[dt] = {}
        self.bendingY[dt] = {}
        self.bendingXY[dt] = {}
        self.shearYZ[dt] = {}
        self.shearXZ[dt] = {}

    def add(self, nnodes, dt, data):
        [eid, parent, coord, icord, theta, forces] = data
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        self.theta[eid] = theta

        self.membraneX[eid] = {}
        self.membraneY[eid] = {}
        self.membraneXY[eid] = {}
        self.bendingX[eid] = {}
        self.bendingY[eid] = {}
        self.bendingXY[eid] = {}
        self.shearYZ[eid] = {}
        self.shearXZ[eid] = {}

        for force in forces:
            [nid, membraneX, membraneY, membraneXY, bendingX,
                bendingY, bendingXY, shearYZ, shearXZ] = force
            self.membraneX[eid][nid] = membraneX
            self.membraneY[eid][nid] = membraneY
            self.membraneXY[eid][nid] = membraneXY
            self.bendingX[eid][nid] = bendingX
            self.bendingY[eid][nid] = bendingY
            self.bendingXY[eid][nid] = bendingXY
            self.shearYZ[eid][nid] = shearYZ
            self.shearXZ[eid][nid] = shearXZ

    def add_sort1(self, nnodes, dt, data):
        [eid, parent, coord, icord, theta, forces] = data
        self._fill_object(dt, eid, parent, coord, icord, theta, forces)

    def add_sort2(self, nnodes, eid, data):
        [dt, parent, coord, icord, theta, forces] = data
        self._fill_object(dt, eid, parent, coord, icord, theta, forces)

    def _fill_object(self, dt, eid, parent, coord, icord, theta, forces):
        if dt not in self.membraneX:
            self.add_new_transient(dt)
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        self.theta[eid] = theta

        self.membraneX[dt][eid] = {}
        self.membraneY[dt][eid] = {}
        self.membraneXY[dt][eid] = {}
        self.bendingX[dt][eid] = {}
        self.bendingY[dt][eid] = {}
        self.bendingXY[dt][eid] = {}
        self.shearYZ[dt][eid] = {}
        self.shearXZ[dt][eid] = {}

        for force in forces:
            [nid, membraneX, membraneY, membraneXY, bendingX,
                bendingY, bendingXY, shearYZ, shearXZ] = force
            self.membraneX[dt][eid][nid] = membraneX
            self.membraneY[dt][eid][nid] = membraneY
            self.membraneXY[dt][eid][nid] = membraneXY
            self.bendingX[dt][eid][nid] = bendingX
            self.bendingY[dt][eid][nid] = bendingY
            self.bendingXY[dt][eid][nid] = bendingXY
            self.shearYZ[dt][eid][nid] = shearYZ
            self.shearXZ[dt][eid][nid] = shearXZ


'                  C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )'
'                                                          (REAL/IMAGINARY)'
' '
'    ELEMENT                    - MEMBRANE  FORCES -                      - BENDING   MOMENTS -            - TRANSVERSE SHEAR FORCES -'
'      ID       GRID-ID     FX            FY            FXY           MX            MY            MXY           QX            QY'
'0       100    CEN/8  0.0           0.0           0.0           0.0           0.0           0.0          -3.492460E-10 -1.368206E-09'
'                      0.0           0.0           0.0           0.0           0.0           0.0           2.910383E-11  5.088840E-10'
''
