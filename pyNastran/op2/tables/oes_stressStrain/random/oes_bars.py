from typing import List

import numpy as np
from numpy import zeros, searchsorted, ravel

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import (
    StressObject, StrainObject, OES_Object)
from pyNastran.f06.f06_formatting import write_floats_13e, _eigenvalue_header
from pyNastran.utils.numpy_utils import integer_types

#oxx = 0. # max from bending and axial
#txz = 1. # from transverse shear; txz=Vz/(Kz*A)
#txy = 1. # from transverse shear; txy=Vz/(Ky*A)
#t = 2. # from torsional stress; t=T*C/J
#ovm = (oxx**2 + 3 * (txy**2 + txz**2 + t**2))**0.5

class RandomBarArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.ielement = 0
        self.nelements = 0  # result specific

        #if not is_sort1:
            #raise NotImplementedError('SORT2')
            #assert dt is not None
            #self.add = self.add_sort2
            #self.add_new_eid = self.add_new_eid_sort2
            #self.addNewNode = self.addNewNodeSort2

    @property
    def is_real(self):
        return True

    @property
    def is_complex(self):
        return False

    def _reset_indices(self):
        self.itotal = 0
        if self.table_name not in ['OESRMS2', 'OESNO2', 'OSTRRMS2', 'OSTRNO2']:
            self.ielement = 0

    def _get_msgs(self):
        raise NotImplementedError('%s needs to implement _get_msgs' % self.__class__.__name__)

    def get_headers(self):
        raise NotImplementedError('%s needs to implement get_headers' % self.__class__.__name__)

    def build(self):
        """sizes the vectorized attributes of the RealBarArray"""
        #print("self.ielement =", self.ielement)
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal

        # buggy MSC 2005 (was this ever fixed?)
        # NX doesn't have this bug
        if self.table_name in ['OESRMS2', 'OESNO2', 'OSTRRMS2', 'OSTRNO2']:
            self.ntotal = self.nelements

        #if self.element_type == 34:
            #nnodes_per_element = 1
        #else:
            #raise NotImplementedError(self.element_type)

        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("***name=%s type=%s nnodes_per_element=%s ntimes=%s nelements=%s ntotal=%s" % (
            #self.element_name, self.element_type, nnodes_per_element, self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.ntotal, dtype='int32')

        #[s1a, s2a, s3a, s4a, axial,
        # s1b, s2b, s3b, s4b]
        self.data = zeros((self.ntimes, self.ntotal, 9), dtype='float32')

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()
        if self.nonlinear_factor not in (None, np.nan):
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
            self.data_frame.index.names = ['ElementID', 'Item']
        else:
            self.data_frame = pd.Panel(self.data, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = ['Static']
            self.data_frame.index.names = ['ElementID', 'Item']

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            i = 0
            if self.is_sort1:
                for itime in range(ntimes):
                    for ieid, eid, in enumerate(self.element):
                        t1 = self.data[itime, ieid, :]
                        t2 = table.data[itime, ieid, :]
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
                raise NotImplementedError(self.is_sort2)
            if i > 0:
                print(msg)
                raise ValueError(msg)
        return True

    def add_new_eid_sort1(self, dt, eid,
                          s1a, s2a, s3a, s4a, axial,
                          s1b, s2b, s3b, s4b):

        assert isinstance(eid, integer_types)
        assert eid > 0, eid
        self._times[self.itime] = dt
        self.element[self.itotal] = eid
        self.data[self.itime, self.itotal, :] = [s1a, s2a, s3a, s4a, axial,
                                                 s1b, s2b, s3b, s4b]
        self.itotal += 1
        self.ielement += 1

    #def add_sort1(self, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm):
        #assert eid is not None
        #msg = "i=%s dt=%s eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g ovmShear=%g" % (
            #self.itotal, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm)
        ##print(msg)
        #if isinstance(nodeID, str):
            #nodeID = 0
        ##assert isinstance(nodeID, integer_types), nodeID
        #self.element_node[self.itotal, :] = [eid, nodeID]
        #self.data[self.itime, self.itotal, :] = [fd, oxx, oyy, txy, angle, majorP, minorP, ovm]
        #self.itotal += 1

    def get_stats(self, short=False) -> List[str]:
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.ntotal
        ntimes = self.ntimes
        unused_ntotal = self.ntotal
        nelements = self.ntotal

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i\n'
                       % (self.__class__.__name__, ntimes, nelements))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i\n'
                       % (self.__class__.__name__, nelements))
            ntimes_word = '1'
        headers = self.get_headers()

        n = len(headers)
        assert n == self.data.shape[2], 'nheaders=%s shape=%s' % (n, str(self.data.shape))
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element.shape = %s\n' % str(self.element.shape).replace('L', ''))
        msg.append('  element type: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element)  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        ind = ravel([searchsorted(self.element == eid) for eid in eids])
        #ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg = self._get_msgs()
        ntimes = self.data.shape[0]
        eids = self.element
        #print('CBAR ntimes=%s ntotal=%s' % (ntimes, ntotal))
        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg))

            s1a = self.data[itime, :, 0]
            s2a = self.data[itime, :, 1]
            s3a = self.data[itime, :, 2]
            s4a = self.data[itime, :, 3]

            axial = self.data[itime, :, 4]
            s1b = self.data[itime, :, 5]
            s2b = self.data[itime, :, 6]
            s3b = self.data[itime, :, 7]
            s4b = self.data[itime, :, 8]

            for (eid, s1ai, s2ai, s3ai, s4ai, axiali, s1bi, s2bi, s3bi, s4bi) in zip(
                eids, s1a, s2a, s3a, s4a, axial, s1b, s2b, s3b, s4b):

                vals = [s1ai, s2ai, s3ai, s4ai, axiali,
                        s1bi, s2bi, s3bi, s4bi,]
                vals2 = write_floats_13e(vals)
                [s1ai, s2ai, s3ai, s4ai, axiali,
                 s1bi, s2bi, s3bi, s4bi,] = vals2

                f06_file.write('\n      %-13s     ENDA         %-13s  %-13s  %-13s  %-13s       %s\n'
                               '0     %-13s     ENDB         %-13s  %-13s  %-13s  %s\n'
                               % (eid, s1ai, s2ai, s3ai, s4ai, axiali,
                                  '', s1bi, s2bi, s3bi, s4bi))
                #f06_file.write('0%8i   %-13s  %-13s  %-13s  %-13s  %s\n'
                               #' %8s   %-13s  %-13s  %-13s  %s\n'
                               #% (eid, s1ai, s2ai, s3ai, s4ai, axiali,
                                  #'', s1bi, s2bi, s3bi, s4bi))

            f06_file.write(page_stamp % page_num)
            page_num += 1

        if self.nonlinear_factor in (None, np.nan):
            page_num -= 1
        return page_num


class RandomBarStressArray(RandomBarArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RandomBarArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        headers = ['s1a', 's2a', 's3a', 's4a', 'axial',
                   's1b', 's2b', 's3b', 's4b']
        return headers

    def _get_msgs(self):
        if self.element_type == 34:
            pass
        else:
            raise NotImplementedError(self.element_type)

        msg = [
            '                                 S T R E S S E S   I N   B A R   E L E M E N T S          ( C B A R )\n',
        ]
        if self.table_name in ['OESATO1', 'OESATO2']:
            msg += ['                                                 ( AUTO-CORRELATION FUNCTION )\n']
        elif self.table_name in ['OESPSD1', 'OESPSD2']:
            msg += ['                                             ( POWER SPECTRAL DENSITY FUNCTION )\n']
        elif self.table_name in ['OESRMS1', 'OESRMS2', 'OESXRMS1']:
            msg += ['                                                     ( ROOT MEAN SQUARE )\n']
        elif self.table_name in ['OESCRM1', 'OESCRM2']:
            msg += ['                                               ( CUMULATIVE ROOT MEAN SQUARE )\n']
        elif self.table_name in ['OESNO1', 'OESNO2']:
            msg += ['                                                 ( NUMBER OF ZERO CROSSINGS )\n']
        else:
            raise NotImplementedError(self.table_name)
        msg += [
            '\n  ELEMENT        SA1            SA2            SA3            SA4           AXIAL\n',
            '    ID.          SB1            SB2            SB3            SB4           STRESS\n',
        ]
        return msg


class RandomBarStrainArray(RandomBarArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RandomBarArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        headers = ['e1a', 'e2a', 'e3a', 'e4a', 'axial',
                   'e1b', 'e2b', 'e3b', 'e4b',]
        return headers

    def _get_msgs(self):
        if self.element_type == 34:
            pass
        else:
            raise NotImplementedError(self.element_type)

        msg = [
            '                                  S T R A I N S    I N   B A R   E L E M E N T S          ( C B A R )\n',
        ]
        if self.table_name in ['OSTRATO1', 'OSTRATO2']:
            msg += ['                                                 ( AUTO-CORRELATION FUNCTION )\n']
        elif self.table_name in ['OSTRPSD1', 'OSTRPSD2']:
            msg += ['                                             ( POWER SPECTRAL DENSITY FUNCTION )\n']
        elif self.table_name in ['OSTRRMS1', 'OSTRRMS2']:
            msg += ['                                                     ( ROOT MEAN SQUARE )\n']
        elif self.table_name in ['OSTRCRM1', 'OSTRCRM2']:
            msg += ['                                               ( CUMULATIVE ROOT MEAN SQUARE )\n']
        elif self.table_name in ['OSTRNO1', 'OSTRNO2']:
            msg += ['                                                 ( NUMBER OF ZERO CROSSINGS )\n']
        else:
            raise NotImplementedError(self.table_name)
        msg += [
            '\n  ELEMENT        SA1            SA2            SA3            SA4           AXIAL\n',
            '    ID.          SB1            SB2            SB3            SB4           STRAIN\n',
        ]

        return msg

