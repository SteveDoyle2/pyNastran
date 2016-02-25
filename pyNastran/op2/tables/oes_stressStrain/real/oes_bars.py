from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
from itertools import count
import numpy as np
from numpy import zeros, searchsorted, ravel

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import write_floats_13e, _eigenvalue_header
try:
    import pandas as pd
except ImportError:
    pass


#oxx = 0. # max from bending and axial
#txz = 1. # from transverse shear; txz=Vz/(Kz*A)
#txy = 1. # from transverse shear; txy=Vz/(Ky*A)
#t = 2. # from torsional stress; t=T*C/J
#ovm = (oxx**2 + 3 * (txy**2 + txz**2 + t**2))**0.5

class RealBarArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.ielement = 0
        self.nelements = 0  # result specific

        if is_sort1:
            if dt is not None:
                #self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
                #self.addNewNode = self.addNewNodeSort1
        else:
            raise NotImplementedError('SORT2')
            #assert dt is not None
            #self.add = self.add_sort2
            #self.add_new_eid = self.add_new_eid_sort2
            #self.addNewNode = self.addNewNodeSort2

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
        raise NotImplementedError('%s needs to implement get_headers' % self.__class__.__name__)

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
        self.element = zeros(self.ntotal, dtype='int32')

        #[s1a, s2a, s3a, s4a, axial, smaxa, smina, MS_tension,
        # s1b, s2b, s3b, s4b,        sminb, sminb, MS_compression]
        self.data = zeros((self.ntimes, self.ntotal, 15), dtype='float32')

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

    def add_new_eid(self, dt, eid, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,
                    s1b, s2b, s3b, s4b, smaxb, sminb, MSc):
        self.add_new_eid_sort1(dt, eid,
                               s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,
                               s1b, s2b, s3b, s4b, smaxb, sminb, MSc)

    def add_new_eid_sort1(self, dt, eid,
                          s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,
                          s1b, s2b, s3b, s4b, smaxb, sminb, MSc):

        assert isinstance(eid, int)
        self._times[self.itime] = dt
        self.element[self.itotal] = eid
        self.data[self.itime, self.itotal, :] = [s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,
                                                 s1b, s2b, s3b, s4b,        smaxb, sminb, MSc]
        self.itotal += 1
        self.ielement += 1

    #def add_sort1(self, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm):
        #assert eid is not None
        #msg = "i=%s dt=%s eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g ovmShear=%g" % (
            #self.itotal, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm)
        ##print(msg)
        #if isinstance(nodeID, string_types):
            #nodeID = 0
        ##assert isinstance(nodeID, int), nodeID
        #self.element_node[self.itotal, :] = [eid, nodeID]
        #self.data[self.itime, self.itotal, :] = [fd, oxx, oyy, txy, angle, majorP, minorP, ovm]
        #self.itotal += 1

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
        itot = searchsorted(eids, self.element)  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        ind = ravel([searchsorted(self.element == eid) for eid in eids])
        #ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg = self._get_msgs()
        (ntimes, ntotal) = self.data.shape[:2]
        eids = self.element
        #print('CBAR ntimes=%s ntotal=%s' % (ntimes, ntotal))
        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg))

            s1a = self.data[itime, :, 0]
            s2a = self.data[itime, :, 1]
            s3a = self.data[itime, :, 2]
            s4a = self.data[itime, :, 3]

            axial = self.data[itime, :, 4]
            smaxa = self.data[itime, :, 5]
            smina = self.data[itime, :, 6]
            MSt = self.data[itime, :, 7]

            s1b = self.data[itime, :, 8]
            s2b = self.data[itime, :, 9]
            s3b = self.data[itime, :, 10]
            s4b = self.data[itime, :, 11]

            smaxb = self.data[itime, :, 12]
            sminb = self.data[itime, :, 13]
            MSc = self.data[itime, :, 14]

            for (i, eid, s1ai, s2ai, s3ai, s4ai, axiali, smaxai, sminai, MSti,
                         s1bi, s2bi, s3bi, s4bi,         smaxbi, sminbi, MSci) in zip(
                count(), eids, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,
                               s1b, s2b, s3b, s4b,        smaxb, sminb, MSc):

                vals = [s1ai, s2ai, s3ai, s4ai, axiali, smaxai, sminai, MSti,
                        s1bi, s2bi, s3bi, s4bi,         smaxbi, sminbi, MSci]
                vals2 = write_floats_13e(vals)
                [s1ai, s2ai, s3ai, s4ai, axiali, smaxai, sminai, MSti,
                 s1bi, s2bi, s3bi, s4bi,         smaxbi, sminbi, MSci] = vals2
                f.write('0%8i   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s %s\n'
                        ' %8s   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s %s\n'
                        % (eid, s1ai, s2ai, s3ai, s4ai, axiali, smaxai, sminai, MSti,
                            '', s1bi, s2bi, s3bi, s4bi, '',     smaxbi, sminbi, MSci))

            f.write(page_stamp % page_num)
            page_num += 1

        if self.nonlinear_factor is None:
            page_num -= 1
        return page_num


class RealBarStressArray(RealBarArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealBarArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        headers = ['s1a', 's2a', 's3a', 's4a', 'axial', 'smaxa', 'smina', 'MS_tension',
                   's1b', 's2b', 's3b', 's4b',          'smaxb', 'sminb', 'MS_compression']
        return headers

    def _get_msgs(self):
        if self.element_type == 34:
            pass
        else:
            raise NotImplementedError(self.element_type)

        msg = [
            '                                 S T R E S S E S   I N   B A R   E L E M E N T S          ( C B A R )\n',
            '  ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T\n',
            '    ID.          SB1            SB2            SB3            SB4           STRESS         SB-MAX         SB-MIN     M.S.-C\n',
        ]
        return msg


class RealBarStrainArray(RealBarArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealBarArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        headers = ['e1a', 'e2a', 'e3a', 'e4a', 'axial', 'emaxa', 'emina', 'MS_tension',
                   'e1b', 'e2b', 'e3b', 'e4b',          'emaxb', 'eminb', 'MS_compression']
        return headers

    def _get_msgs(self):
        if self.element_type == 34:
            pass
        else:
            raise NotImplementedError(self.element_type)

        msg = [
            '                                  S T R A I N S    I N   B A R   E L E M E N T S          ( C B A R )\n',
            '  ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T\n',
            '    ID.          SB1            SB2            SB3            SB4           STRAIN         SB-MAX         SB-MIN     M.S.-C\n',
        ]
        return msg
