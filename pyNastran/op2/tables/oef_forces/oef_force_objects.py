#pylint disable=C0301
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
from six.moves import zip, range
from itertools import cycle
import numpy as np
from numpy import zeros, searchsorted, array_equal, allclose, sqrt, pi

from pyNastran.op2.result_objects.op2_objects import ScalarObject
from pyNastran.f06.f06_formatting import write_floats_13e, writeFloats12E, _eigenvalue_header, get_key0, write_float_13e
try:
    import pandas as pd
except ImportError:
    pass

class RealForceObject(ScalarObject):
    def __init__(self, data_code, isubcase, apply_data_code=True):
        self.element_type = None
        self.element_name = None
        self.nonlinear_factor = None
        self.element = None
        self._times = None
        ScalarObject.__init__(self, data_code, isubcase, apply_data_code=apply_data_code)

    def is_real(self):
        return True

    def is_complex(self):
        return False

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        raise NotImplementedError()

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


class FailureIndices(RealForceObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealForceObject.__init__(self, data_code, isubcase)
        self.nelements = 0  # result specific
        if is_sort1:
            self.add = self.add_sort1
        else:
            raise NotImplementedError('SORT2')

    def build(self):
        if self.is_built:
            return

        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
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

        #[force]
        self.data = zeros((self.ntimes, self.nelements, 1), dtype='float32')

    def build_dataframe(self):
        headers = self.get_headers()
        if self.nonlinear_factor is not None:
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
        else:
            self.data_frame = pd.Panel(self.data, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = ['Static']
        self.data_frame.index.names = ['ElementID', 'Item']

    def add(self, dt, eid, force):
        self.add_sort1(dt, eid, force)

    def add_sort1(self, dt, eid, force):
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [force]
        self.ielement += 1

    def get_stats(self):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

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
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True, is_sort1=True):
        raise NotImplementedError('this should be overwritten by %s' % (self.__class__.__name__))

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        asd
        if self.is_sort1():
            page_num = self._write_sort1_as_sort1(header, page_stamp, page_num, f, msg_temp)
        else:
            raise NotImplementedError(self.code_information())
            #page_num = self._write_sort2_as_sort2(header, page_stamp, page_num, f, msg_temp)

        '          F A I L U R E   I N D I C E S   F O R   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )\n'
        '   ELEMENT  FAILURE        PLY   FP=FAILURE INDEX FOR PLY    FB=FAILURE INDEX FOR BONDING   FAILURE INDEX FOR ELEMENT      FLAG\n'
        '     ID      THEORY         ID  (DIRECT STRESSES/STRAINS)     (INTER-LAMINAR STRESSES)      MAX OF FP,FB FOR ALL PLIES\n'
        '         1   HOFFMAN       101      6.987186E-02      \n'
        '                                                                     1.687182E-02                                              \n'
        '                           102      9.048269E-02      \n'
        '                                                                     1.721401E-02                                               \n'
        return page_num


class RealSpringDamperForceArray(RealForceObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealForceObject.__init__(self, data_code, isubcase)
        self.nelements = 0  # result specific
        if is_sort1:
            self.add = self.add_sort1
        else:
            raise NotImplementedError('SORT2')

    def build(self):
        if self.is_built:
            return

        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
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

        #[force]
        self.data = zeros((self.ntimes, self.nelements, 1), dtype='float32')

    def build_dataframe(self):
        headers = self.get_headers()
        if self.nonlinear_factor is not None:
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
        else:
            self.data_frame = pd.Panel(self.data, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = ['Static']
        self.data_frame.index.names = ['ElementID', 'Item']

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()
        assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (self.table_name, table.table_name)
        assert self.approach_code == table.approach_code
        if self.nonlinear_factor is not None:
            assert np.array_equal(self._times, table._times), 'ename=%s-%s times=%s table.times=%s' % (
                self.element_name, self.element_type, self._times, table._times)
        if not np.array_equal(self.element, table.element):
            assert self.element.shape == table.element.shape, 'shape=%s element.shape=%s' % (self.element.shape, table.element.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\nEid1, Eid2\n' % str(self.code_information())
            for eid, eid2 in zip(self.element, table.element):
                msg += '%s, %s\n' % (eid, eid2)
            print(msg)
            raise ValueError(msg)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            i = 0
            if self.is_sort1():
                for itime in range(ntimes):
                    for ieid, eid, in enumerate(self.element):
                        t1 = self.data[itime, inid, :]
                        t2 = table.data[itime, inid, :]
                        (force1, stress1) = t1
                        (force2, stress2) = t2
                        if not allclose(t1, t2):
                        #if not np.array_equal(t1, t2):
                            msg += '%s\n  (%s, %s, %s, %s, %s, %s)\n  (%s, %s, %s, %s, %s, %s)\n' % (
                                eid,
                                force1, stress1,
                                force2, stress2)
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

    def add(self, dt, eid, force):
        self.add_sort1(dt, eid, force)

    def add_sort1(self, dt, eid, force):
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [force]
        self.ielement += 1

    def get_stats(self):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

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
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True, is_sort1=True):
        raise NotImplementedError('this should be overwritten by %s' % (self.__class__.__name__))

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)

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
                out.append([eid, write_float_13e(stressi)])

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

class RealSpringForceArray(RealSpringDamperForceArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealSpringDamperForceArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def get_headers(self):
        headers = ['spring_force']
        return headers

    def get_f06_header(self, is_mag_phase=True, is_sort1=True):
        if self.element_type == 11:  # CELAS1
            msg = ['                              F O R C E S   I N   S C A L A R   S P R I N G S        ( C E L A S 1 )\n']
        elif self.element_type == 12:  # CELAS2
            msg = ['                              F O R C E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )\n']
        elif self.element_type == 13:  # CELAS3
            msg = ['                              F O R C E S   I N   S C A L A R   S P R I N G S        ( C E L A S 3 )\n']
        elif self.element_type == 14:  # CELAS4
            msg = ['                              F O R C E S   I N   S C A L A R   S P R I N G S        ( C E L A S 4 )\n']
        else:
            raise NotImplementedError('element_name=%s element_type=%s' % (self.element_name, self.element_type))

        msg += [
            '      ELEMENT         FORCE            ELEMENT         FORCE            ELEMENT         FORCE            ELEMENT         FORCE\n'
            '        ID.                              ID.                              ID.                              ID.\n'
        ]
        return msg

class RealDamperForceArray(RealSpringDamperForceArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealSpringDamperForceArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def get_headers(self):
        headers = ['damper_force']
        return headers

    def get_f06_header(self, is_mag_phase=True, is_sort1=True):
        if self.element_type == 20:  # CDAMP1
            msg = ['                              F O R C E S   I N   S C A L A R   D A M P E R S        ( C D A M P 1 )\n']
        elif self.element_type == 21:  # CDAMP2
            msg = ['                              S T R E S S E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )\n']
        #elif self.element_type == 13:  # CDAMP3
            #msg = ['                              S T R E S S E S   I N   S C A L A R   S P R I N G S        ( C E L A S 3 )\n']
        #elif self.element_type == 14:  # CDAMP4
            #msg = ['                              S T R E S S E S   I N   S C A L A R   S P R I N G S        ( C E L A S 4 )\n']
        else:
            raise NotImplementedError('element_name=%s element_type=%s' % (self.element_name, self.element_type))

        if is_sort1:
            msg += [
                '      ELEMENT         FORCE            ELEMENT         FORCE            ELEMENT         FORCE            ELEMENT         FORCE\n'
                '        ID.                              ID.                              ID.                              ID.\n'
            ]
        else:
            msg += [
                '                         AXIAL                                                       AXIAL\n'
                '       TIME              FORCE         TORQUE                      TIME              FORCE         TORQUE\n'
        ]
        return msg


class RealRodForceArray(RealForceObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealForceObject.__init__(self, data_code, isubcase)
        self.nelements = 0  # result specific

        if is_sort1:
            self.add = self.add_sort1
        else:
            raise NotImplementedError('SORT2')

    def get_headers(self):
        headers = ['axial', 'torsion']
        return headers

    #def get_headers(self):
        #headers = ['axial', 'torque']
        #return headers

    def _get_msgs(self):
        base_msg = ['       ELEMENT       AXIAL       TORSIONAL     ELEMENT       AXIAL       TORSIONAL\n',
                    '         ID.         FORCE        MOMENT        ID.          FORCE        MOMENT\n']
        crod_msg = ['                                     F O R C E S   I N   R O D   E L E M E N T S      ( C R O D )\n', ]
        conrod_msg = ['                                     F O R C E S   I N   R O D   E L E M E N T S      ( C O N R O D )\n', ]
        ctube_msg = ['                                     F O R C E S   I N   R O D   E L E M E N T S      ( C T U B E )\n', ]
        crod_msg += base_msg
        conrod_msg += base_msg
        ctube_msg += base_msg
        return crod_msg, conrod_msg, ctube_msg

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

        #[axial_force, torque]
        self.data = zeros((self.ntimes, self.nelements, 2), dtype='float32')

    def build_dataframe(self):
        headers = self.get_headers()
        if self.nonlinear_factor is not None:
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
        else:
            self.data_frame = pd.Panel(self.data, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = ['Static']
        self.data_frame.index.names = ['ElementID', 'Item']

    def add(self, dt, eid, axial, torque):
        self.add_sort1(dt, eid, axial, torque)

    def add_sort1(self, dt, eid, axial, torque):
        # print('ielement=%s' % self.ielement)
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [axial, torque]
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
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n  ' % self.element_name)
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
        return self.element_name, msg

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        (elem_name, msg_temp) = self.get_f06_header(is_mag_phase)

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

            out = []
            for eid, axiali, torsioni in zip(eids, axial, torsion):
                [axiali, torsioni] = write_floats_13e([axiali, torsioni])
                out.append([eid, axiali, torsioni])

            for i in range(0, nwrite, 2):
                out_line = '      %8i   %-13s  %-13s  %8i   %-13s  %s\n' % tuple(out[i] + out[i + 1])
                f.write(out_line)
            if is_odd:
                out_line = '      %8i   %-13s  %s\n' % tuple(out[-1])
                f.write(out_line)
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()
        assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (self.table_name, table.table_name)
        assert self.approach_code == table.approach_code
        if self.nonlinear_factor is not None:
            assert np.array_equal(self._times, table._times), 'ename=%s-%s times=%s table.times=%s' % (
                self.element_name, self.element_type, self._times, table._times)
        if not np.array_equal(self.element, table.element):
            assert self.element.shape == table.element.shape, 'element shape=%s table.shape=%s' % (self.element.shape, table.element.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += 'Eid\n'
            for eid, eid2 in zip(self.element, table.element):
                msg += '%s, %s\n' % (eid, eid2)
            print(msg)
            raise ValueError(msg)
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

                    if not np.array_equal(t1, t2):
                        msg += '(%s)    (%s, %s)  (%s, %s)\n' % (
                            eid,
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


class RealCBeamForceArray(ScalarObject):
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
        return True

    def is_complex(self):
        return False

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
        self.data = zeros((self.ntimes, self.ntotal, 8), 'float32')

    def finalize(self):
        sd = self.data[0, :, 0]
        i_sd_zero = np.where(sd != 0.0)[0]
        i_node_zero = np.where(self.element_node[:, 1] != 0)[0]
        assert i_node_zero.max() > 0, 'CBEAM element_node hasnt been filled'
        i = np.union1d(i_sd_zero, i_node_zero)

        #self.nelements = len(self.element) // 11
        self.element = self.element[i]
        self.element_node = self.element_node[i, :]
        self.data = self.data[:, i, :]

    def build_dataframe(self):
        headers = self.get_headers()
        element_location = [
            self.element_node[:, 0],
            self.data[0, :, 0],
        ]
        if self.nonlinear_factor is not None:
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data[:, :, 1:], items=column_values, major_axis=element_location, minor_axis=headers[1:]).to_frame()
            self.data_frame.columns.names = column_names
            self.data_frame.index.names = ['ElementID', 'Location', 'Item']
        else:
            df1 = pd.DataFrame(element_location).T
            df1.columns = ['ElementID', 'Location']
            df2 = pd.DataFrame(self.data[0])
            df2.columns = headers
            self.data_frame = df1.join([df2])
        #self.data_frame = self.data_frame.reset_index().replace({'NodeID': {0:'CEN'}}).set_index(['ElementID', 'NodeID'])
        #print(self.data_frame)

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
                    for ieid, eid, in enumerate(self.element):
                        t1 = self.data[itime, inid, :]
                        t2 = table.data[itime, inid, :]
                        (axial_stress1, equiv_stress1, total_strain1, effective_plastic_creep_strain1, effective_creep_strain1, linear_torsional_stress1) = t1
                        (axial_stress2, equiv_stress2, total_strain2, effective_plastic_creep_strain2, effective_creep_strain2, linear_torsional_stress2) = t2
                        if not np.allclose(t1, t2):
                        #if not np.array_equal(t1, t2):
                            msg += '%s\n  (%s, %s, %s, %s, %s, %s)\n  (%s, %s, %s, %s, %s, %s)\n' % (
                                eid,
                                axial_stress1, equiv_stress1, total_strain1, effective_plastic_creep_strain1, effective_creep_strain1, linear_torsional_stress1,
                                axial_stress2, equiv_stress2, total_strain2, effective_plastic_creep_strain2, effective_creep_strain2, linear_torsional_stress2)
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

    def add_new_element_sort1(self, dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq):
        return self.add_sort1(dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq)

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
        #name = self.data_code['name']
        #if name == 'freq':
            #name = 'FREQUENCY'
        #else: # mode
            #raise RuntimeError(name)
        if is_sort1:
            msg_temp = [
                '                                 F O R C E S   I N   B E A M   E L E M E N T S        ( C B E A M )\n',
                '                    STAT DIST/   - BENDING MOMENTS -            - WEB  SHEARS -           AXIAL          TOTAL          WARPING\n',
                '   ELEMENT-ID  GRID   LENGTH    PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE         TORQUE\n']
        else:
            raise NotImplementedError('CBEAM-SORT2')

        if self.is_sort1():
            assert self.is_sort1() == True, str(self)
            #if is_sort1:
            page_num = self._write_sort1_as_sort1(f, page_num, page_stamp, header, msg_temp)
            #else:
                #self._write_sort1_as_sort2(f, page_num, page_stamp, header, msg_temp)
        else:
            assert self.is_sort1() == True, str(self)
        return page_num - 1

    def get_headers(self):
        headers = [
            'sd', 'bending_moment1', 'bending_moment2', 'shear1', 'shear2',
            'axial_force', 'total_torque', 'warping_torque', ]
        return headers

    def _write_sort1_as_sort1(self, f, page_num, page_stamp, header, msg_temp):
        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]
        long_form = False
        if nids.min() == 0:
            #header = [
            #    '                         S T R E S S   D I S T R I B U T I O N   I N   B A R   E L E M E N T S       ( C B A R )\n'
            #    '0    ELEMENT  STATION    SXC           SXD           SXE           SXF            AXIAL          S-MAX         S-MIN         M.S.-T\n'
            #    '       ID.     (PCT)                                                                                                         M.S.-C\n'
            #]
            msg = header + [
                 '                                 F O R C E S   I N   B E A M   E L E M E N T S        ( C B E A M )\n',
                 '                    STAT DIST/   - BENDING MOMENTS -            - WEB  SHEARS -           AXIAL          TOTAL          WARPING\n',
                 '   ELEMENT-ID  GRID   LENGTH    PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE         TORQUE\n']

            long_form = True

        times = self._times
        ntimes = self.data.shape[0]
        for itime in range(ntimes):
            if self.nonlinear_factor is not None:
                dt = self._times[itime]
                dt_line = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
                header[1] = dt_line
            msg = header + msg_temp
            f.write(''.join(msg))

            #sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq
            assert self.is_sort1() == True, str(self)
            sd = self.data[itime, :, 0]
            bm1 = self.data[itime, :, 1]
            bm2 = self.data[itime, :, 2]
            ts1 = self.data[itime, :, 3]
            ts2 = self.data[itime, :, 4]
            af = self.data[itime, :, 5]
            ttrq = self.data[itime, :, 6]
            wtrq = self.data[itime, :, 7]

            for eid, nid, sdi, bm1i, bm2i, ts1i, ts2i, afi, ttrqi, wtrqi in zip(eids, nids, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq):
                vals = (bm1i, bm2i, ts1i, ts2i, afi, ttrqi, wtrqi)
                vals2 = write_floats_13e(vals)
                (sbm1i, sbm2i, sts1i, sts2i, safi, sttrqi, swtrq) = vals2

                if long_form:
                    f.write('           %8i   %.3f   %-13s %-13s  %-13s %-13s  %-13s  %-13s  %s\n' % (
                        eid, sdi, sbm1i, sbm2i, sts1i, sts2i, safi, sttrqi, swtrq))
                else:
                    if sdi == 0.:
                        f.write('0  %8i\n' % eid)
                    f.write('           %8i   %.3f   %-13s %-13s  %-13s %-13s  %-13s  %-13s  %s\n' % (
                        nid, sdi, sbm1i, sbm2i, sts1i, sts2i, safi, sttrqi, swtrq))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num


class RealCShearForceArray(ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        if is_sort1:
            self.add = self.add_sort1
        else:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        headers = [
            'force41', 'force21', 'force12', 'force32', 'force23', 'force43',
            'force34', 'force14',
            'kick_force1', 'shear12',  'kick_force2', 'shear23',
            'kick_force3', 'shear34', 'kick_force4', 'shear41',
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

        #[force41, force21, force12, force32, force23, force43,
        # force34, force14,
        # kick_force1, shear12, kick_force2, shear23,
        # kick_force3, shear34, kick_force4, shear41]
        self.data = zeros((self.ntimes, self.ntotal, 16), dtype='float32')

    def build_dataframe(self):
        headers = self.get_headers()
        if self.nonlinear_factor is not None:
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
            self.data_frame.index.names = ['ElementID', 'Item']
        else:
            self.data_frame = pd.Panel(self.data, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = ['Static']
            self.data_frame.index.names = ['ElementID', 'Item']

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()
        assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (self.table_name, table.table_name)
        assert self.approach_code == table.approach_code
        if self.nonlinear_factor is not None:
            assert np.array_equal(self._times, table._times), 'ename=%s-%s times=%s table.times=%s' % (
                self.element_name, self.element_type, self._times, table._times)
        if not np.array_equal(self.element, table.element):
            assert self.element.shape == table.element.shape, 'element shape=%s table.shape=%s' % (self.element.shape, table.element.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += 'Eid\n'
            for eid1, eid2 in zip(self.element, table.element):
                msg += '%s, %s\n' % (eid1, eid2)
            print(msg)
            raise ValueError(msg)
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

                    if not np.array_equal(t1, t2):
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

    def add(self, dt, eid,
            force41, force14, force21, force12, force32, force23, force43, force34,
            kick_force1, kick_force2, kick_force3, kick_force4,
            shear12, shear23, shear34, shear41):
        self.add_sort1(dt, eid,
                       force41, force14, force21, force12, force32, force23, force43, force34,
                       kick_force1, kick_force2, kick_force3, kick_force4,
                       shear12, shear23, shear34, shear41)

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
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = [
            '                           F O R C E S   A C T I N G   O N   S H E A R   P A N E L   E L E M E N T S   (CSHEAR)\n'
            ' \n'
            '                  ====== POINT  1 ======      ====== POINT  2 ======      ====== POINT  3 ======      ====== POINT  4 ======\n'
            '   ELEMENT        F-FROM-4      F-FROM-2      F-FROM-1      F-FROM-3      F-FROM-2      F-FROM-4      F-FROM-3      F-FROM-1\n'
            '         ID               KICK-1       SHEAR-12       KICK-2       SHEAR-23       KICK-3       SHEAR-34       KICK-4       SHEAR-41\n'
        ]

        #(elem_name, msg_temp) = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element

        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            f14 = self.data[itime, :, 0]
            f12 = self.data[itime, :, 1]
            f21 = self.data[itime, :, 2]
            f23 = self.data[itime, :, 3]
            f32 = self.data[itime, :, 4]
            f34 = self.data[itime, :, 5]
            f43 = self.data[itime, :, 6]
            f41 = self.data[itime, :, 7]

            kick1 = self.data[itime, :, 8]
            tau12 = self.data[itime, :, 9]
            kick2 = self.data[itime, :, 10]
            tau23 = self.data[itime, :, 11]
            kick3 = self.data[itime, :, 12]
            tau34 = self.data[itime, :, 13]
            kick4 = self.data[itime, :, 14]
            tau41 = self.data[itime, :, 15]

            zip_in = [
                f14, f12, f21, f23, f32, f34, f43, f41,
                kick1, tau12, kick2, tau23, kick3, tau34, kick4, tau41,
            ]
            for (eid, f14i, f12i, f21i, f23i, f32i, f34i, f43i, f41i, kick1i, tau12i, kick2i, tau23i, kick3i, tau34i, kick4i, tau41i) in zip(
                eids, f14, f12, f21, f23, f32, f34, f43, f41, kick1, tau12, kick2, tau23, kick3, tau34, kick4, tau41):
                (vals2, is_all_zeros) = writeFloats12E(
                    [f14i, f12i, f21i, f23i, f32i, f34i, f43i, f41i,
                     kick1i, tau12i, kick2i, tau23i, kick3i, tau34i, kick4i, tau41i])
                [
                    f14i, f12i,
                    f21i, f23i,
                    f32i, f34i,
                    f43i, f41i,
                    kick1i, tau12i, kick2i, tau23i,
                    kick3i, tau34i, kick4i, tau41i
                ] = vals2
                f.write('0%13i%-13s %-13s %-13s %-13s %-13s %-13s %-13s %s\n'
                        '                     %-13s %-13s %-13s %-13s %-13s %-13s %-13s %s\n' % (
                            eid, f14, f12, f21, f23, f32, f34, f43, f41,
                            kick1, tau12, kick2, tau23, kick3, tau34, kick4, tau41))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealViscForceArray(RealForceObject):  # 24-CVISC
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealForceObject.__init__(self, data_code, isubcase)

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        if is_sort1:
            self.add = self.add_sort1
        else:
            raise NotImplementedError('SORT2')

    def get_headers(self):
        headers = ['axial', 'torsion']
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
        self.data = zeros((self.ntimes, self.ntotal, 2), dtype='float32')

    def build_dataframe(self):
        headers = self.get_headers()
        if self.nonlinear_factor is not None:
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
        else:
            self.data_frame = pd.Panel(self.data, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = ['Static']
        self.data_frame.index.names = ['ElementID', 'Item']

    def add(self, dt, eid, axial, torque):
        self.add_sort1(dt, eid, axial, torque)

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
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True, is_sort1=True):
        if is_sort1:
            msg = [
                '                                     F O R C E S   I N   R O D   E L E M E N T S      ( C R O D )\n'
                ' \n'
                '       ELEMENT       AXIAL       TORSIONAL     ELEMENT       AXIAL       TORSIONAL\n'
                '         ID.         FORCE        MOMENT        ID.          FORCE        MOMENT\n'
            ]
        else:
            msg = [
                '                                   F O R C E S   I N   V I S C   E L E M E N T S   ( C V I S C )'
                ' \n'
                '                         AXIAL                                                       AXIAL\n'
                '       TIME              FORCE         TORQUE                      TIME              FORCE         TORQUE\n'
                #'   0.0                0.0            0.0                       1.000000E+00      -5.642718E-04   0.0\n'
                #'   2.000000E+00      -1.905584E-06   0.0                       3.000000E+00       9.472010E-07   0.0\n'
            ]
        return msg

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = self.get_f06_header(is_mag_phase)

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

            out = []
            for eid, axiali, torsioni in zip(eids, axial, torsion):
                [axiali, torsioni] = write_floats_13e([axiali, torsioni])
                out.append([eid, axiali, torsioni])

            for i in range(0, nwrite, 2):
                out_line = '      %8i   %-13s  %-13s  %8i   %-13s  %s\n' % tuple(out[i] + out[i + 1])
                f.write(out_line)
            if is_odd:
                out_line = '      %8i   %-13s  %s\n' % tuple(out[-1])
                f.write(out_line)
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()
        assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (self.table_name, table.table_name)
        assert self.approach_code == table.approach_code
        if not np.array_equal(self.element, table.element):
            assert self.element.shape == table.element.shape, 'element shape=%s table.shape=%s' % (self.element.shape, table.element.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += 'Eid\n'
            for eid1, eid2 in zip(self.element, table.element):
                msg += '%s, %s\n' % (eid1, eid2)
            print(msg)
            raise ValueError(msg)
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

                    if not np.array_equal(t1, t2):
                        msg += '(%s)    (%s, %s)  (%s, %s)\n' % (
                            eid,
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


class RealPlateForceArray(RealForceObject):  # 33-CQUAD4, 74-CTRIA3
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealForceObject.__init__(self, data_code, isubcase)

        self.dt = dt
        self.nelements = 0

        #if is_sort1:
            #if dt is not None:
                #self.add = self.add_sort1
        #else:
            #assert dt is not None
            #self.add = self.add_sort2

    def _get_msgs(self):
        raise NotImplementedError()

    def get_headers(self):
        return ['mx', 'my', 'mxy', 'bmx', 'bmy', 'bmxy', 'tx', 'ty']

    def build(self):
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        #self.nelements //= self.ntimes
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
        self.element = zeros(self.ntotal, dtype='int32')

        #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
        self.data = zeros((self.ntimes, self.ntotal, 8), dtype='float32')

    def build_dataframe(self):
        headers = self.get_headers()
        assert 0 not in self.element
        if self.nonlinear_factor is not None:
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
        else:
            self.data_frame = pd.Panel(self.data, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = ['Static']
        self.data_frame.index.names = ['ElementID', 'Item']

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()

        self._eq_header(table)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, e in enumerate(self.element_node):
                    (eid, nid) = e
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (fiber_dist1, oxx1, oyy1, txy1, angle1, majorP1, minorP1, ovm1) = t1
                    (fiber_dist2, oxx2, oyy2, txy2, angle2, majorP2, minorP2, ovm2) = t2

                    # vm stress can be NaN for some reason...
                    if not np.array_equal(t1[:-1], t2[:-1]):
                        msg += '(%s, %s)    (%s, %s, %s, %s, %s, %s, %s, %s)  (%s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                            eid, nid,
                            fiber_dist1, oxx1, oyy1, txy1, angle1, majorP1, minorP1, ovm1,
                            fiber_dist2, oxx2, oyy2, txy2, angle2, majorP2, minorP2, ovm2)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    #def add_new_eid_sort1(self, dt, eid, axial, SMa, torsion, SMt):
        #self._times[self.itime] = dt
        #self.element[self.ielement] = eid
        #self.data[self.itime, self.ielement, :] = [axial, SMa, torsion, SMt]
        #self.ielement += 1

    #def add(self, dt, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        #self.add_sort1(dt, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)

    def add_sort1(self, dt, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        self._times[self.itime] = dt
        self.element[self.itotal] = eid
        self.data[self.itime, self.itotal, :] = [mx, my, mxy, bmx, bmy, bmxy, tx, ty]
        self.itotal += 1

    def add_sort2(self, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        raise NotImplementedError('SORT2')
        #if dt not in self.mx:
            #self.add_new_transient(dt)
        self.data[self.itime, self.itotal, :] = [mx, my, mxy, bmx, bmy, bmxy, tx, ty]

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
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True):
        if 'CTRIA3' in self.element_name:
            msg = [
                '                             F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'
                ' \n'
                '    ELEMENT                - MEMBRANE  FORCES -                        - BENDING  MOMENTS -              - TRANSVERSE SHEAR FORCES -\n'
                '      ID              FX            FY            FXY             MX            MY            MXY             QX            QY\n'
            ]
            nnodes = 3
        elif 'CQUAD4' in self.element_name:
            msg = [
                '                          F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'
                ' \n'
                '    ELEMENT                    - MEMBRANE  FORCES -                      - BENDING   MOMENTS -            - TRANSVERSE SHEAR FORCES -\n'
                '      ID       GRID-ID     FX            FY            FXY           MX            MY            MXY           QX            QY\n'
            ]
            nnodes = 4
        else:
            raise NotImplementedError(self.element_name)
        return self.element_name, nnodes, msg

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        (elem_name, nnodes, msg_temp) = self.get_f06_header(is_mag_phase)

        # write the f06
        ntimes = self.data.shape[0]

        eids = self.element
        cen_word = 'CEN/%i' % nnodes
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
            mx = self.data[itime, :, 0]
            my = self.data[itime, :, 1]
            mxy = self.data[itime, :, 2]
            bmx = self.data[itime, :, 3]
            bmy = self.data[itime, :, 4]
            bmxy = self.data[itime, :, 5]
            tx = self.data[itime, :, 6]
            ty = self.data[itime, :, 7]

            if self.element_type == 74:
                for eid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi in zip(eids, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
                    [mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi] = write_floats_13e(
                     [mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi])
                    # ctria3
                    #          8      -7.954568E+01  2.560061E+03 -4.476376E+01    1.925648E+00  1.914048E+00  3.593237E-01    8.491534E+00  5.596094E-01  #
                    f.write('   %8i %18s %13s %13s   %13s %13s %13s   %13s %s\n' % (eid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi))

            elif self.element_type == 33:
                for eid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi in zip(eids, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
                    [mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi] = write_floats_13e(
                        [mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi])
                    # cquad4
                    #0         6    CEN/4  1.072685E+01  2.504399E+03 -2.455727E+01 -5.017930E+00 -2.081427E+01 -5.902618E-01 -9.126162E+00  4.194400E+01#
                    #Fmt = '% 8i   ' + '%27.20E   ' * 8 + '\n'
                    #f.write(Fmt % (eid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi))
                    #
                    f.write('0 %8i %8s %13s %13s %13s %13s %13s %13s %13s %s\n' % (eid, cen_word, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi))
            else:
                raise NotImplementedError(self.element_type)
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealPlateBilinearForceArray(RealForceObject):  # 144-CQUAD4
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealForceObject.__init__(self, data_code, isubcase)

        self.dt = dt
        self.nelements = 0

        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.add_sort2

    def _get_msgs(self):
        raise NotImplementedError()

    def get_headers(self):
        return ['mx', 'my', 'mxy', 'bmx', 'bmy', 'bmxy', 'tx', 'ty']

    def build(self):
        # print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        #self.nelements //= self.ntimes
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
        self.element_node = zeros((self.ntotal, 2), dtype='int32')

        # -MEMBRANE FORCES-   -BENDING MOMENTS- -TRANSVERSE SHEAR FORCES -
        #     FX FY FXY           MX MY MXY            QX QY
        #[fx, fy, fxy,  mx,  my,  mxy, qx, qy]
        #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
        self.data = zeros((self.ntimes, self.ntotal, 8), dtype='float32')

    def build_dataframe(self):
        headers = self.get_headers()
        element_node = [self.element_node[:, 0], self.element_node[:, 1]]
        if self.nonlinear_factor is not None:
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=element_node, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
            self.data_frame.index.names = ['ElementID', 'NodeID', 'Item']
        else:
            df1 = pd.DataFrame(element_node).T
            df1.columns = ['ElementID', 'NodeID']
            df2 = pd.DataFrame(self.data[0])
            df2.columns = headers
            self.data_frame = df1.join(df2)
        self.data_frame = self.data_frame.reset_index().replace({'NodeID': {0:'CEN'}}).set_index(['ElementID', 'NodeID'])
        #print(self.data_frame)

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

                    if not np.array_equal(t1, t2):
                        msg += '(%s, %s)    (%s, %s, %s, %s, %s, %s, %s, %s)  (%s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                            eid, nid,
                            mx1, my1, mxy1, bmx1, bmy1, bmxy1, tx1, ty1,
                            mx2, my2, mxy2, bmx2, bmy2, bmxy2, tx2, ty2)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add(self, dt, eid, term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        self.add_sort1(dt, eid, term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)

    def add_sort1(self, dt, eid, term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        self._times[self.itime] = dt
        self.element_node[self.itotal] = [eid, nid]
        self.data[self.itime, self.itotal, :] = [mx, my, mxy, bmx, bmy, bmxy, tx, ty]
        self.itotal += 1

    def add_sort2(self, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        raise NotImplementedError('SORT2')
        #if dt not in self.mx:
            #self.add_new_transient(dt)
        self.data[self.itime, self.itotal, :] = [mx, my, mxy, bmx, bmy, bmxy, tx, ty]

    @property
    def nnodes_per_element(self):
        if self.element_type == 144:  # CQUAD4
            nnodes_element = 5
        elif self.element_type == 64:  # CQUAD8
            nnodes_element = 5
        elif self.element_type == 82:  # CQUADR
            nnodes_element = 5
        elif self.element_type == 75:  # CTRIA6
            nnodes_element = 4
        elif self.element_type == 70:  # CTRIAR
            nnodes_element = 4
        else:
            raise NotImplementedError('element_type=%s element_name=%s' % (self.element_type, self.element_name))
        return nnodes_element

    def get_stats(self):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        ntotal = self.ntotal

        msg = []
        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i ntotal=%i nnodes/element=%i\n'
                       % (self.__class__.__name__, ntimes, nelements, ntotal, self.nnodes_per_element))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i ntotal=%i nnodes/element=%i\n'
                       % (self.__class__.__name__, nelements, ntotal, self.nnodes_per_element))
            ntimes_word = 1

        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True):
        # if 'CTRIA3' in self.element_name:
            # msg = [
                # '                             F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'
                # ' \n'
                # '    ELEMENT                - MEMBRANE  FORCES -                        - BENDING  MOMENTS -              - TRANSVERSE SHEAR FORCES -\n'
                # '      ID              FX            FY            FXY             MX            MY            MXY             QX            QY\n'
            # ]
            # nnodes = 3

        if self.element_type == 70:
            # CQUAD4
            element_name = 'CTRIAR'
            msg = [
                '                             F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'
                ' \n'
                '    ELEMENT                    - MEMBRANE  FORCES -                      - BENDING   MOMENTS -            - TRANSVERSE SHEAR FORCES -\n'
                '      ID       GRID-ID     FX            FY            FXY           MX            MY            MXY           QX            QY\n'
            ]
            nnodes = 6
        elif self.element_type == 75:
            # CQUAD4
            element_name = 'CTRIA6'
            msg = [
                '                             F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n'
                ' \n'
                '    ELEMENT                    - MEMBRANE  FORCES -                      - BENDING   MOMENTS -            - TRANSVERSE SHEAR FORCES -\n'
                '      ID       GRID-ID     FX            FY            FXY           MX            MY            MXY           QX            QY\n'
            ]
            nnodes = 6
        elif self.element_type == 64:
            # CQUAD4
            element_name = 'CQUAD8'
            msg = [
                '                          F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n'
                ' \n'
                '    ELEMENT                    - MEMBRANE  FORCES -                      - BENDING   MOMENTS -            - TRANSVERSE SHEAR FORCES -\n'
                '      ID       GRID-ID     FX            FY            FXY           MX            MY            MXY           QX            QY\n'
            ]
            nnodes = 8
        elif self.element_type == 82:
            # CQUAD4
            element_name = 'CQUADR'
            msg = [
                '                          F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n'
                ' \n'
                '    ELEMENT                    - MEMBRANE  FORCES -                      - BENDING   MOMENTS -            - TRANSVERSE SHEAR FORCES -\n'
                '      ID       GRID-ID     FX            FY            FXY           MX            MY            MXY           QX            QY\n'
            ]
            nnodes = 4
        elif self.element_type == 144:
            # CQUAD4
            element_name = 'CQUAD4'
            msg = [
                '                          F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n'
                ' \n'
                '    ELEMENT                    - MEMBRANE  FORCES -                      - BENDING   MOMENTS -            - TRANSVERSE SHEAR FORCES -\n'
                '      ID       GRID-ID     FX            FY            FXY           MX            MY            MXY           QX            QY\n'
            ]
            nnodes = 4
        else:
            raise NotImplementedError('element_name=%s element_type=%s' % (self.element_name, self.element_type))
        return element_name, nnodes, msg

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        (elem_name, nnodes, msg_temp) = self.get_f06_header(is_mag_phase)

        # write the f06
        ntimes = self.data.shape[0]

        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]
        cen_word = 'CEN/%i' % nnodes
        if self.element_type  in [64, 82, 144]: # CQUAD8, CQUADR, CQUAD4
            cyc = cycle([0, 1, 2, 3, 4])
        elif self.element_type  in [70, 75]: # CTRIAR, CTRIA6
            cyc = cycle([0, 1, 2, 3])
        else:
            raise NotImplementedError(self.element_type)

        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
            mx = self.data[itime, :, 0]
            my = self.data[itime, :, 1]
            mxy = self.data[itime, :, 2]
            bmx = self.data[itime, :, 3]
            bmy = self.data[itime, :, 4]
            bmxy = self.data[itime, :, 5]
            tx = self.data[itime, :, 6]
            ty = self.data[itime, :, 7]

            for i, eid, nid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi in zip(cyc, eids, nids, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
                [mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi] = write_floats_13e(
                    [mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi])
                # ctria3
                #          8      -7.954568E+01  2.560061E+03 -4.476376E+01    1.925648E+00  1.914048E+00  3.593237E-01    8.491534E+00  5.596094E-01  #
                if i == 0:
                    f.write('0  %8i    CEN/4 %-13s %-13s %-13s %-13s %-13s %-13s %-13s %s\n' % (eid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi))
                else:
                    f.write('            %8i %-13s %-13s %-13s %-13s %-13s %-13s %-13s %s\n' % (nid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi))
            # else:
                # raise NotImplementedError(self.element_type)
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealCBarForceArray(ScalarObject):  # 34-CBAR
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        if is_sort1:
            self.add = self.add_sort1
        else:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def is_real(self):
        return True

    def is_complex(self):
        return False

    def get_headers(self):
        headers = [
            'bending_moment_a1', 'bending_moment_a2',
            'bending_moment_b1', 'bending_moment_b2',
            'shear1', 'shear2',
            'axial', 'torque']
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

        #[bending_moment_a1, bending_moment_a2, bending_moment_b1, bending_moment_b2, shear1, shear2, axial, torque]
        self.data = zeros((self.ntimes, self.ntotal, 8), dtype='float32')

    def build_dataframe(self):
        headers = self.get_headers()
        if self.nonlinear_factor is not None:
            column_names, column_values = self._build_dataframe_transient_header()
            # Create a 3D Panel
            #column_values = [modes, freq]
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
            #self.data_frame = self.data_frame.to_frame()

            # Define names for column labels
            self.data_frame.columns.names = column_names
        else:
            self.data_frame = pd.Panel(self.data, major_axis=self.element, minor_axis=headers).to_frame()
            #cbar_forces = cbar_forces.to_frame()
            self.data_frame.columns.names = ['Static']
        # Define names for the row labels
        self.data_frame.index.names = ['ElementID', 'Item']

    def add(self, dt, data):
        self.add_sort1(dt, data)

    def add_sort1(self, dt, data):
        data = [eid, bending_moment_a1, bending_moment_a2,
                bending_moment_b1, bending_moment_b2, shear1, shear2, axial, torque] = data
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [
            bending_moment_a1, bending_moment_a2,
            bending_moment_b1, bending_moment_b2,
            shear1, shear2, axial, torque]
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
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def eid_to_element_node_index(self, eids):
        ind = searchsorted(eids, self.element)
        return ind

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        words = ['                                 F O R C E S   I N   B A R   E L E M E N T S         ( C B A R )\n',
                 '0    ELEMENT         BEND-MOMENT END-A            BEND-MOMENT END-B                - SHEAR -               AXIAL\n',
                 '       ID.         PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2         FORCE         TORQUE\n']
        msg = []
        #header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
        eids = self.element
        #f.write(''.join(words))

        ntimes = self.data.shape[0]
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + words))

            bm1a = self.data[itime, :, 0]
            bm2a = self.data[itime, :, 1]
            bm1b = self.data[itime, :, 2]
            bm2b = self.data[itime, :, 3]
            ts1 = self.data[itime, :, 4]
            ts2 = self.data[itime, :, 5]
            af = self.data[itime, :, 6]
            trq = self.data[itime, :, 7]
            for eid, bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi in zip(eids, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq):
                [bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi] = write_floats_13e([
                 bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi])
            f.write('     %8i    %-13s %-13s  %-13s %-13s  %-13s %-13s  %-13s  %s\n' % (eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq))
            f.write(page_stamp % page_num)
        return page_num

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
                    (bm1a1, bm2a1, bm1b1, bm2b1, ts11, ts21, af1, trq1) = t1
                    (bm1a2, bm2a2, bm1b2, bm2b2, ts12, ts22, af2, trq2) = t2

                    if not np.array_equal(t1, t2):
                        msg += '(%s)    (%s, %s, %s, %s, %s, %s, %s, %s)  (%s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                            eid,
                            bm1a1, bm2a1, bm1b1, bm2b1, ts11, ts21, af1, trq1,
                            bm1a2, bm2a2, bm1b2, bm2b2, ts12, ts22, af2, trq2)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True


class RealConeAxForceArray(ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        if is_sort1:
            self.add = self.add_sort1
        else:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        headers = [
            'hopa', 'bmu', 'bmv', 'tm', 'su', 'sv'
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

        #[hopa, bmu, bmv, tm, su, sv]
        self.data = zeros((self.ntimes, self.ntotal, 6), dtype='float32')

    def build_dataframe(self):
        headers = self.get_headers()
        if self.nonlinear_factor is not None:
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
            self.data_frame.index.names = ['ElementID', 'Item']
        else:
            df1 = pd.DataFrame(self.element)
            df1.columns = ['ElementID']
            df2 = pd.DataFrame(self.data[0])
            df2.columns = headers
            self.data_frame = df1.join([df2])
        #print(self.data_frame)

    def __eq__(self, table):
        self._eq_header(table)
        assert self.is_sort1() == table.is_sort1()

        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, e in enumerate(self.element):
                    eid = e
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (hopa1, bmu1, bmv1, tm1, su1, sv1) = t1
                    (hopa2, bmu2, bmv2, tm2, su2, sv2) = t2

                    if not np.array_equal(t1, t2):
                        msg += (
                            '%s   (%s, %s, %s, %s, %s, %s)\n'
                            '     (%s, %s, %s, %s, %s, %s)\n' % (
                                eid,
                                hopa1, bmu1, bmv1, tm1, su1, sv1,
                                hopa2, bmu2, bmv2, tm2, su2, sv2,
                            ))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add(self, dt, eid, hopa, bmu, bmv, tm, su, sv):
        self.add_sort1(dt, eid, hopa, bmu, bmv, tm, su, sv)

    def add_sort1(self, dt, eid, hopa, bmu, bmv, tm, su, sv):
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [hopa, bmu, bmv, tm, su, sv]
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
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = [
            '             F O R C E S   I N   A X I S - S Y M M E T R I C   C O N I C A L   S H E L L   E L E M E N T S   (CCONEAX)\n'
            ' \n'
            '  ELEMENT     HARMONIC    POINT           BEND-MOMENT       BEND-MOMENT      TWIST-MOMENT           SHEAR            SHEAR\n'
            '   ID.         NUMBER     ANGLE               V                 U                                     V                U\n'
            #'      101        0                       5.864739E-09      1.759422E-09      0.0               0.0               0.0'
            #'      101                0.0000          5.864739E-09      1.759422E-09      0.0               0.0               0.0'
        ]

        #(elem_name, msg_temp) = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            hopa = self.data[itime, :, 0]
            bmu = self.data[itime, :, 1]
            bmv = self.data[itime, :, 2]
            tm = self.data[itime, :, 3]
            su = self.data[itime, :, 4]
            sv = self.data[itime, :, 5]

            for (eid, hopai, bmui, bmvi, tmi, sui, svi) in zip(
                eids, hopa, bmu, bmv, tm, su, sv):
                if hopai > 0.1:
                    raise NotImplementedError(hopai)

                vals2 = write_floats_13e(
                    [hopai, bmui, bmvi, tmi, sui, svi])
                [hopai, bmui, bmvi, tmi, sui, svi] = vals2

                # TODO: hopa is probably the wrong type
                              # hopa        # Mu       Mv      twist       Vy       Vu
                f.write(' %8i  %-13s  %-13s %-13s     %-13s     %-13s     %-13s     %s\n' % (
                    eid, 0.0, hopai, bmui, bmvi, tmi, sui, svi))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealCBar100ForceArray(RealForceObject):  # 100-CBAR
    """
    CBAR-34s are converted to CBAR-100s when you have PLOAD1s
    (distributed bar loads).  The number of stations by default is 2,
    but with a CBARAO, you can change this (max of 8 points; 6 internal
    points).

    If you use a CBARO without PLOAD1s, you wil turn CBAR-34s into
    CBAR-100s as well.
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealForceObject.__init__(self, data_code, isubcase)

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        if is_sort1:
            self.add = self.add_sort1
        else:
            raise NotImplementedError('SORT2')

    def get_headers(self):
        headers = [
            'station', 'bending_moment1', 'bending_moment2', 'shear1', 'shear2', 'axial', 'torque'
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
        self.element = zeros(self.nelements, dtype='int32')

        # [station, bending_moment1, bending_moment2, shear1, shear2, axial, torque]
        self.data = zeros((self.ntimes, self.ntotal, 7), dtype='float32')

    #def finalize(self):
        #sd = self.data[0, :, 0]
        #i_sd_zero = np.where(sd != 0.0)[0]
        #i_node_zero = np.where(self.element_node[:, 1] != 0)[0]
        #assert i_node_zero.max() > 0, 'CBAR element_node hasnt been filled'
        #i = np.union1d(i_sd_zero, i_node_zero)
        #self.element = self.element[i]
        #self.element_node = self.element_node[i, :]
        #self.data = self.data[:, i, :]

    def build_dataframe(self):
        headers = self.get_headers()
        element_location = [
            self.element,
            self.data[0, :, 0],
        ]
        if self.nonlinear_factor is not None:
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data[:, :, 1:], items=column_values, major_axis=element_location, minor_axis=headers[1:]).to_frame()
            self.data_frame.columns.names = column_names
            self.data_frame.index.names = ['ElementID', 'Location', 'Item']
        else:
            df1 = pd.DataFrame(element_location).T
            df1.columns = ['ElementID', 'Location']
            df2 = pd.DataFrame(self.data[0])
            df2.columns = headers
            self.data_frame = df1.join([df2])
        #self.data_frame = self.data_frame.reset_index().replace({'NodeID': {0:'CEN'}}).set_index(['ElementID', 'NodeID'])
        #print(self.data_frame)

    def add(self, eid, sd, bm1, bm2, ts1, ts2, af, trq):
        self.add_sort1(dt, eid, sd, bm1, bm2, ts1, ts2, af, trq)

    def add_sort1(self, dt, eid, sd, bm1, bm2, ts1, ts2, af, trq):
        self._times[self.itime] = dt
        self.element[self.ielement] = eid

        # station, bending_moment1, bending_moment2, shear1, shear2, axial, torque
        self.data[self.itime, self.ielement, :] = [sd, bm1, bm2, ts1, ts2, af, trq]
        self.ielement += 1

    def get_stats(self):
        if not self.is_built:
            return ['<%s>\n' % self.__class__.__name__,
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
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def eid_to_element_node_index(self, eids):
        ind = searchsorted(eids, self.element)
        return ind

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        #'                         F O R C E   D I S T R I B U T I O N   I N   B A R   E L E M E N T S          ( C B A R )'
        #'0    ELEMENT  STATION         BEND-MOMENT                      SHEAR FORCE                     AXIAL'
        #'       ID.     (PCT)     PLANE 1        PLANE 2           PLANE 1        PLANE 2               FORCE              TORQUE'
        #'           10   0.000   0.0           -5.982597E+06      0.0           -7.851454E+03          0.0                0.0'
        #'           10   1.000   0.0            1.868857E+06      0.0           -7.851454E+03          0.0                0.0'
        #'           11   0.000   0.0            1.868857E+06      0.0           -7.851454E+03          0.0                0.0'
        #'           11   0.050   0.0            2.261429E+06      0.0           -7.851454E+03          0.0                0.0'
        #'           11   0.100   0.0            2.654002E+06      0.0           -7.851454E+03          0.0                0.0'
        words = [
            '                         F O R C E   D I S T R I B U T I O N   I N   B A R   E L E M E N T S          ( C B A R )\n'
            '0    ELEMENT  STATION         BEND-MOMENT                      SHEAR FORCE                     AXIAL\n'
            '       ID.     (PCT)     PLANE 1        PLANE 2           PLANE 1        PLANE 2               FORCE              TORQUE\n']
            # '        15893   0.000   1.998833E+02   9.004551E+01      2.316835E+00   1.461960E+00         -2.662207E+03       9.795244E-02'
        msg = []
        #header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
        eids = self.element
        #f.write(''.join(words))

        ntimes = self.data.shape[0]
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + words))

            # sd, bm1, bm2, ts1, ts2, af, trq
            sd = self.data[itime, :, 0]
            bm1 = self.data[itime, :, 1]
            bm2 = self.data[itime, :, 2]
            ts1 = self.data[itime, :, 3]
            ts2 = self.data[itime, :, 4]
            af = self.data[itime, :, 5]
            trq = self.data[itime, :, 6]
            for eid, sdi, bm1i, bm2i, ts1i, ts2i, afi, trqi in zip(eids, sd, bm1, bm2, ts1, ts2, af, trq):
                [bm1i, bm2i, ts1i, ts2i, afi, trqi] = write_floats_13e([
                 bm1i, bm2i, ts1i, ts2i, afi, trqi])
            #f.write('     %8i  %4.3f  %-13s  %-13s     %-13s  %-13s     %-13s     %s\n' % (
                #eid, sd, bm1i, bm2i, ts1i, ts2i, afi, trqi))
            f.write('     %8i   %4.3f  %-13s  %-13s     %-13s  %-13s         %-13s      %s\n' % (
                eid, sdi, bm1i, bm2i, ts1i, ts2i, afi, trqi))
            f.write(page_stamp % page_num)
        return page_num

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
                    (sd1, bm11, bm21, ts11, ts21, af1, trq1) = t1
                    (sd2, bm12, bm22, ts12, ts22, af2, trq2) = t2

                    if not np.array_equal(t1, t2):
                        msg += '(%s)    (%s, %s, %s, %s, %s, %s, %s)  (%s, %s, %s, %s, %s, %s, %s)\n' % (
                            eid,
                            sd1, bm11, bm21, ts11, ts21, af1, trq1,
                            sd2, bm12, bm22, ts12, ts22, af2, trq2)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True


class RealCGapForceArray(ScalarObject):  # 38-CGAP
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        if is_sort1:
            self.add = self.add_sort1
        else:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        headers = [
            'fx', 'sfy', 'sfz', 'u', 'v', 'w', 'sv', 'sw'
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
        self.element = zeros(self.nelements, dtype='int32')

        # [fx, sfy, sfz, u, v, w, sv, sw]
        self.data = zeros((self.ntimes, self.ntotal, 8), dtype='float32')

    def build_dataframe(self):
        headers = self.get_headers()
        if self.nonlinear_factor is not None:
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
            self.data_frame.index.names = ['ElementID', 'Item']
        else:
            self.data_frame = pd.Panel(self.data, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = ['Static']
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
                    eid = e
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (fx1, sfy1, sfz1, u1, v1, w1, sv1, sw1) = t1
                    (fx2, sfy2, sfz2, u2, v2, w2, sv2, sw2) = t2

                    if not np.array_equal(t1, t2):
                        msg += (
                            '%s   (%s, %s, %s, %s, %s, %s, %s, %s)\n'
                            '     (%s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                                eid,
                                fx1, sfy1, sfz1, u1, v1, w1, sv1, sw1,
                                fx2, sfy2, sfz2, u2, v2, w2, sv2, sw2,
                            ))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add(self, dt, eid, fx, sfy, sfz, u, v, w, sv, sw):
        self.add_sort1(dt, eid, fx, sfy, sfz, u, v, w, sv, sw)

    def add_sort1(self, dt, eid, fx, sfy, sfz, u, v, w, sv, sw):
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [fx, sfy, sfz, u, v, w, sv, sw]
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
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = [
            '                                    F O R C E S   I N   G A P   E L E M E N T S   ( C G A P )\n'
            '  ELEMENT   - F O R C E S  I N  E L E M  S Y S T E M -         - D I S P L A C E M E N T S  I N  E L E M  S Y S T E M -\n'
            '    ID       COMP-X         SHEAR-Y        SHEAR-Z        AXIAL-U        TOTAL-V        TOTAL-W        SLIP-V         SLIP-W\n'
            #'     101   3.333333E+04   0.0            0.0            1.149425E-04   0.0            0.0            0.0            0.0\n'
        ]

        ##(elem_name, msg_temp) = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        #ntimes, ntotal = self.data.shape[:1]
        ntimes = self.data.shape[0]

        eids = self.element

        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            # [fx, sfy, sfz, u, v, w, sv, sw]
            fx = self.data[itime, :, 0]
            sfy = self.data[itime, :, 1]
            sfz = self.data[itime, :, 2]
            u = self.data[itime, :, 3]
            v = self.data[itime, :, 4]
            w = self.data[itime, :, 5]
            sv = self.data[itime, :, 6]
            sw = self.data[itime, :, 7]

            for (eid, fxi, sfyi, sfzi, ui, vi, wi, svi, swi) in zip(eids, fx, sfy, sfz, u, v, w, sv, sw):
                (vals2, is_all_zeros) = writeFloats12E(
                    [fxi, sfyi, sfzi, ui, vi, wi, svi, swi])
                [fxi, sfyi, sfzi, ui, vi, wi, svi, swi] = vals2
                f.write('0%13i%-13s %-13s %-13s %-13s %-13s %-13s %-13s %s\n' % (
                            eid, fxi, sfyi, sfzi, ui, vi, wi, svi, swi))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealBendForceArray(RealForceObject):  # 69-CBEND
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealForceObject.__init__(self, data_code, isubcase)
        self.nelements = 0  # result specific

    def get_headers(self):
        headers = [
            'bending_moment_1a', 'bending_moment_2a', 'shear_1a', 'shear_2a', 'axial_a', 'torque_a',
            'bending_moment_1b', 'bending_moment_2b', 'shear_1b', 'shear_2b', 'axial_b', 'torque_b',
        ]
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
        self.element_nodes = zeros((self.nelements, 3), dtype='int32')

        #[bending_moment_1a, bending_moment_2a, shear_1a, shear_2a, axial_a, torque_a
        # bending_moment_1b, bending_moment_2b, shear_1b, shear_2b, axial_b, torque_b]
        self.data = zeros((self.ntimes, self.nelements, 12), dtype='float32')

    def build_dataframe(self):
        element = self.element_nodes[:, 0]
        headers = self.get_headers()
        if self.nonlinear_factor is not None:
            # TODO: add NodeA, NodeB
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
            self.data_frame.index.names = ['ElementID', 'Item']
            #print(self.data_frame)
        else:
            df1 = pd.DataFrame(self.element_nodes)
            df1.columns = ['ElementID', 'NodeA', 'NodeB']
            df2 = pd.DataFrame(self.data[0])
            df2.columns = headers
            self.data_frame = df1.join(df2)

    #def add(self, dt, eid, axial, torque):
        #self.add_sort1(dt, eid, axial, torque)

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
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = [
            '                                 F O R C E S   I N   B E N D   E L E M E N T S        ( C B E N D )\n'
            '                                 - BENDING MOMENTS -            -   SHEARS   -            AXIAL\n'
            '   ELEMENT-ID  GRID    END      PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE\n'
            #'0      6901    6901      A     0.0           0.0            0.0           0.0            0.0           -6.305720E-16'
            #'               6902      B    -5.000000E-01  5.000000E-01   1.000000E+00 -1.000000E+00  -5.000000E-07  -1.666537E-07'
        ]
        # write the f06
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element_nodes[:, 0]
        nid_a = self.element_nodes[:, 1]
        nid_b = self.element_nodes[:, 2]

        #print('len(eids)=%s nwrite=%s is_odd=%s' % (len(eids), nwrite, is_odd))
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
                [bending_moment_1ai, bending_moment_2ai, shear_1ai, shear_2ai, axial_ai, torque_ai,
                 bending_moment_1bi, bending_moment_2bi, shear_1bi, shear_2bi, axial_bi, torque_bi] = write_floats_13e(
                     [bending_moment_1ai, bending_moment_2ai, shear_1ai, shear_2ai, axial_ai, torque_ai,
                      bending_moment_1bi, bending_moment_2bi, shear_1bi, shear_2bi, axial_bi, torque_bi])
                f.write(
                    '0  %8i%8i      A    %13s %13s  %13s %13s  %13s  %13s\n'
                    '           %8i      B    %13s %13s  %13s %13s  %13s  %13s\n' % (
                        eid, nid_ai, bending_moment_1ai, bending_moment_2ai, shear_1ai, shear_2ai, axial_ai, torque_ai,
                        nid_bi, bending_moment_1bi, bending_moment_2bi, shear_1bi, shear_2bi, axial_bi, torque_bi
                    ))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def __eq__(self, table):
        self._eq_header(table)
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


                    if not np.array_equal(t1, t2):
                        msg += '(%s)    (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)  (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                            eid,
                            bending_moment_1a1, bending_moment_2a1, shear_1a1, shear_2a1, axial_a1, torque_a1,
                            bending_moment_1b1, bending_moment_2b1, shear_1b1, shear_2b1, axial_b1, torque_b1,

                            bending_moment_1a2, bending_moment_2a2, shear_1a2, shear_2a2, axial_a2, torque_a2,
                            bending_moment_1b2, bending_moment_2b2, shear_1b2, shear_2b2, axial_b2, torque_b2)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True


class RealSolidPressureForceArray(ScalarObject):  # 77-PENTA_PR,78-TETRA_PR
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        if is_sort1:
            self.add = self.add_sort1
        else:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        headers = [
            'ax', 'ay', 'az', 'vx', 'vy', 'vz', 'pressure'
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

        #[ax, ay, az, vx, vy, vz, pressure]
        self.data = zeros((self.ntimes, self.ntotal, 7), dtype='float32')

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

                    if not np.array_equal(t1, t2):
                        msg += (
                            '%s   (%s, %s, %s, %s, %s, %s, %s)\n'
                            '     (%s, %s, %s, %s, %s, %s, %s)\n' % (
                                eid,
                                ax1, ay1, az1, vx1, vy1, vz1, pressure1,
                                ax2, ay2, az2, vx2, vy2, vz2, pressure2,
                            ))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add(self, dt, eid, etype, ax, ay, az, vx, vy, vz, pressure):
        self.add_sort1(dt, eid, etype, ax, ay, az, vx, vy, vz, pressure)

    def add_sort1(self, dt, eid, etype, ax, ay, az, vx, vy, vz, pressure):
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
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        #(elem_name, msg_temp) = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        #(ntimes, ntotal, two) = self.data.shape

        if self.is_sort1():
            page_num = self._write_sort1_as_sort1(header, page_stamp, page_num, f)
        else:
            raise NotImplementedError('sort2')
        return page_num

    def _write_sort2_as_sort1(self, header, page_stamp, page_num=1, f=None):
        msg_temp = ['                                   P E A K   A C C E L E R A T I O N S   A N D   P R E S S U R E S\n',
                    ' \n',
                    '    TIME         EL-TYPE             X-ACCELERATION            Y-ACCELERATION            Z-ACCELERATION            PRESSURE (DB)\n']
        ntimes = self.data.shape[0]

        eids = self.element
        etype = self.element_name
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            vx = self.data[itime, :, 0]
            vy = self.data[itime, :, 1]
            vz = self.data[itime, :, 2]
            ax = self.data[itime, :, 3]
            ay = self.data[itime, :, 4]
            az = self.data[itime, :, 5]
            pressure = self.data[itime, :, 5]

            for (eid, vxi, vyi, vzi, axi, ayi, azi, pressurei) in zip(
                eids, vx, vy, vz, ax, ay, az, pressure):

                vals2 = write_floats_13e([axi, ayi, azi, pressurei])
                [sax, say, saz, spressure] = vals2

                #etype = 'PENPR'
                f.write('0%13s    %5s               %-13s             %-13s             %-13s             %s\n' % (eid, etype, sax, say, saz, spressure))

            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def _write_sort1_as_sort1(self, header, page_stamp, page_num=1, f=None):
        msg_temp = ['                                   P E A K   A C C E L E R A T I O N S   A N D   P R E S S U R E S\n',
                    ' \n',
                    '    ELEMENT-ID   EL-TYPE             X-ACCELERATION            Y-ACCELERATION            Z-ACCELERATION            PRESSURE (DB)\n']  # TODO: bad line...
        ntimes = self.data.shape[0]

        eids = self.element
        etype = self.element_name
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            vx = self.data[itime, :, 0]
            vy = self.data[itime, :, 1]
            vz = self.data[itime, :, 2]
            ax = self.data[itime, :, 3]
            ay = self.data[itime, :, 4]
            az = self.data[itime, :, 5]
            pressure = self.data[itime, :, 5]

            for (eid, vxi, vyi, vzi, axi, ayi, azi, pressurei) in zip(
                eids, vx, vy, vz, ax, ay, az, pressure):

                vals2 = write_floats_13e([axi, ayi, azi, pressurei])
                [sax, say, saz, spressure] = vals2

                #etype = 'PENPR'
                f.write('0%13s    %5s               %-13s             %-13s             %-13s             %s\n' % (eid, etype, sax, say, saz, spressure))

            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealCBushForceArray(ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        if is_sort1:
            self.add = self.add_sort1
        else:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        headers = ['fx', 'fy', 'fz', 'mx', 'my', 'mz']
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

        #[fx, fy, fz, mx, my, mz]
        self.data = zeros((self.ntimes, self.nelements, 6), dtype='float32')

    def build_dataframe(self):
        headers = self.get_headers()
        if self.nonlinear_factor is not None:
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
            self.data_frame.index.names = ['ElementID', 'Item']
        else:
            self.data_frame = pd.Panel(self.data, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = ['Static']
            self.data_frame.index.names = ['ElementID', 'Item']

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()
        assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (self.table_name, table.table_name)
        assert self.approach_code == table.approach_code
        if self.nonlinear_factor is not None:
            assert np.array_equal(self._times, table._times), 'ename=%s-%s times=%s table.times=%s' % (
                self.element_name, self.element_type, self._times, table._times)
        if not np.array_equal(self.element, table.element):
            assert self.element.shape == table.element.shape, 'shape=%s element.shape=%s' % (self.element.shape, table.element.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            for eid, eid2 in zip(self.element, table.element):
                msg += '%s, %s\n' % (eid, eid2)
            print(msg)
            raise ValueError(msg)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            i = 0
            if self.is_sort1():
                for itime in range(ntimes):
                    for ieid, eid, in enumerate(self.element):
                        t1 = self.data[itime, ieid, :]
                        t2 = table.data[itime, ieid, :]
                        (fx1, fy1, fz1, mx1, my1, mz1) = t1
                        (fx2, fy2, fz2, mx2, my2, mz2) = t2
                        if not allclose(t1, t2):
                        #if not np.array_equal(t1, t2):
                            msg += '%s\n  (%s, %s, %s, %s, %s, %s)\n  (%s, %s, %s, %s, %s, %s)\n' % (
                                eid,
                                fx1, fy1, fz1, mx1, my1, mz1,
                                fx2, fy2, fz2, mx2, my2, mz2)
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

    def add(self, dt, eid, fx, fy, fz, mx, my, mz):
        self.add_sort1(dt, eid, fx, fy, fz, mx, my, mz)

    def add_sort1(self, dt, eid, fx, fy, fz, mx, my, mz):
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [fx, fy, fz, mx, my, mz]
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
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self):
        msg = ['                                 F O R C E S   I N   B U S H   E L E M E N T S        ( C B U S H )\n'
        ' \n'
        '                  ELEMENT-ID        FORCE-X       FORCE-Y       FORCE-Z      MOMENT-X      MOMENT-Y      MOMENT-Z  \n']
        #'0                        599      0.0           2.000000E+00  3.421458E-14  1.367133E-13 -3.752247E-15  1.000000E+00\n']
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
        msg_temp = self.get_f06_header()

        # write the f06
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            fx = self.data[itime, :, 0]
            fy = self.data[itime, :, 1]
            fz = self.data[itime, :, 2]
            mx = self.data[itime, :, 3]
            my = self.data[itime, :, 4]
            mz = self.data[itime, :, 5]

            for eid, fxi, fyi, fzi, mxi, myi, mzi in zip(eids, fx, fy, fz, mx, my, mz):
                [fxi, fyi, fzi, mxi, myi, mzi] = write_floats_13e([fxi, fyi, fzi, mxi, myi, mzi])
                f.write('                    %8i     %-13s %-13s %-13s %-13s %-13s %s\n' % (eid, fxi, fyi, fzi, mxi, myi, mzi))
                #'0                        599      0.0           2.000000E+00  3.421458E-14  1.367133E-13 -3.752247E-15  1.000000E+00\n']
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealForce_VU(ScalarObject):  # 191-VUBEAM
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
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

        # handle SORT1 case
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


class RealForce_VU_2D(ScalarObject):  # 190-VUTRIA # 189-VUQUAD
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
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

        # handle SORT1 case
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
