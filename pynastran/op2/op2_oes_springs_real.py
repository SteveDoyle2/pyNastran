from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
from six.moves import zip
import numpy as np
from numpy import zeros, array_equal
from itertools import count

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import write_floats_13e, write_float_13e, _eigenvalue_header
try:
    import pandas as pd
except ImportError:
    pass


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
        self._eq_header(table)
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

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
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


def _write_f06_springs_transient(f, stress, header, words, name):
    raise RuntimeError('is this used?')
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
                stresses = write_floats_13e(stresses)
                f.write('    %10i  %13s    %10i  %13s    %10i  %13s    %10i  %13s\n' % (
                    eids[0], stresses[0],
                    eids[1], stresses[1],
                    eids[2], stresses[2],
                    eids[3], stresses[3]))
                eids = []
                stresses = []

        if stresses:
            line = '    '
            stresses = write_floats_13e(stresses)
            for eid, stress in zip(eids, stresses):
                line += '%10i  %13s    ' % (eid, stress)
            f.write(line.rstrip() + '\n')

            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            msg = ['']
        page_num += 1

    return page_num - 1


def _write_f06_springs(f, data):
    raise RuntimeError('is this used?')
    eids = []
    stresses = []
    for eid, stress in sorted(iteritems(data)):
        eids.append(eid)
        stresses.append(stress)
        if len(stresses) == 4:
            stresses = write_floats_13e(stresses)
            f.write('    %10i  %13s    %10i  %13s    %10i  %13s    %10i  %13s\n' % (
                eids[0], stresses[0],
                eids[1], stresses[1],
                eids[2], stresses[2],
                eids[3], stresses[3]))
            eids = []
            stresses = []

    if stresses:
        line = '    '
        stresses = write_floats_13e(stresses)
        for eid, stress in zip(eids, stresses):
            line += '%10i  %13s    ' % (eid, stress)
        f.write(line.rstrip() + '\n')


class RealNonlinearSpringStressArray(OES_Object):
    """
    ::

      #ELEMENT-ID =     102
                               #N O N L I N E A R   S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )
        #TIME          AXIAL STRESS         EQUIVALENT         TOTAL STRAIN       EFF. STRAIN          EFF. CREEP        LIN. TORSIONAL
                                             #STRESS                             PLASTIC/NLELAST          STRAIN              STRESS
      #2.000E-02        1.941367E+01        1.941367E+01        1.941367E-04        0.0                 0.0                 0.0
      #3.000E-02        1.941367E+01        1.941367E+01        1.941367E-04        0.0                 0.0                 0.0
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=True)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        self.nelements = 0  # result specific

        if is_sort1:
            pass
        else:
            raise NotImplementedError('SORT2')

    def is_real(self):
        return True

    def is_complex(self):
        return False

    def is_stress(self):
        return True

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def _get_msgs(self):
        raise NotImplementedError()

    def get_headers(self):
        headers = ['force', 'stress']
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

        #[force, stress]
        self.data = zeros((self.ntimes, self.nelements, 2), dtype='float32')

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

    def add_sort1(self, dt, eid, force, stress):
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [force, stress]
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

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        if self.is_sort1():
            if self.element_type == 224:
                nspring = 1 # CELAS1
            elif self.element_type == 225:
                nspring = 3 # CELAS3
            else:
                raise NotImplementedError('type=%s name=%s' % (self.element_type, self.element_name))

            msg = [
                '          N O N L I N E A R   F O R C E S  A N D  S T R E S S E S  I N   S C A L A R   S P R I N G S    ( C E L A S %i )\n'
                ' \n'
                 '         ELEMENT-ID          FORCE         STRESS                    ELEMENT-ID          FORCE         STRESS\n' % nspring
                #'         5.000000E-02     2.000000E+01   1.000000E+01                1.000000E-01     4.000000E+01   2.000000E+01'
            ]
            page_num = self._write_sort1_as_sort1(header, page_stamp, page_num, f, msg)
        else:
            msg = [
                '          N O N L I N E A R   F O R C E S  A N D  S T R E S S E S  I N   S C A L A R   S P R I N G S    ( C E L A S %i )\n'
                ' \n'
                '             STEP            FORCE         STRESS                        STEP            FORCE         STRESS\n' % nspring
                #'         5.000000E-02     2.000000E+01   1.000000E+01                1.000000E-01     4.000000E+01   2.000000E+01'
            ]
            raise NotImplementedError('RealNonlinearSpringStressArray-sort2')
        return page_num

    def _write_sort1_as_sort1(self, header, page_stamp, page_num, f, msg_temp):
        """
        ::
              ELEMENT-ID =      20
                  N O N L I N E A R   F O R C E S  A N D  S T R E S S E S  I N   S C A L A R   S P R I N G S    ( C E L A S 1 )

                     STEP            FORCE         STRESS                        STEP            FORCE         STRESS
                 5.000000E-02     2.000000E+01   1.000000E+01                1.000000E-01     4.000000E+01   2.000000E+01
                 1.500000E-01     6.000000E+01   3.000000E+01                2.000000E-01     8.000000E+01   4.000000E+01
        """
        ntimes = self.data.shape[0]

        eids = self.element
        neids = len(eids)
        is_odd = neids % 2 == 1
        if is_odd:
            neids -= 1

        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))
            force = self.data[itime, :, 0]
            stress = self.data[itime, :, 1]
            for i, eid, forcei, stressi, in zip(count(step=2), eids[:neids:2], force[:neids:2], stress[:neids:2]):
                f.write('         %-13i   %-13s  %-13s                %-13s   %-13s  %s\n' % (
                    eid,
                    write_float_13e(forcei),
                    write_float_13e(stressi),
                    eids[i + 1],
                    write_float_13e(force[i + 1]),
                    write_float_13e(stress[i + 1])
                ))
            if is_odd:
                f.write('         %-13i   %-13s  %s\n' % (
                    eids[neids],
                    write_float_13e(force[neids]),
                    write_float_13e(stress[neids])
                ))

            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1
