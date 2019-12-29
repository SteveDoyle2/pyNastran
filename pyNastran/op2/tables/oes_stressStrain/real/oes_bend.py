from itertools import cycle
from typing import List

import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import (
    StressObject, StrainObject, OES_Object)
from pyNastran.f06.f06_formatting import write_floats_13e, write_floats_8p1e


class RealBendArray(OES_Object):
    """
    Common class for:
     - RealBendStressArray
     - RealBendStrainArray
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)   ## why???
        self.element_node = None
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        #self.itime = 0
        self.nelements = 0  # result specific

        #if is_sort1:
            #pass
        #else:
            #raise NotImplementedError('SORT2')

    @property
    def is_real(self):
        return False

    @property
    def is_complex(self):
        return True

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    #def get_nnodes(self):
        #return get_nnodes(self)

    def build(self):
        """sizes the vectorized attributes of the ComplexBendArray"""
        if not hasattr(self, 'subtitle'):
            self.subtitle = self.data_code['subtitle']
        #print('ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (
            #self.ntimes, self.nelements, self.ntotal, self.subtitle))
        nnodes = 1

        #self.names = []
        #self.nelements //= nnodes
        self.nelements //= self.ntimes
        self.ntotal = self.nelements * nnodes * 2
        #self.ntotal
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        self.is_built = True
        #print('ntotal=%s ntimes=%s nelements=%s' % (self.ntotal, self.ntimes, self.nelements))

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        self._times = np.zeros(self.ntimes, 'float32')
        #self.ntotal = self.nelements * nnodes

        self.element_node = np.zeros((self.ntotal, 2), 'int32')

        # the number is messed up because of the offset for the element's properties
        if not self.nelements * nnodes * 2 == self.ntotal:
            msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (
                self.ntimes, self.nelements, nnodes, self.nelements * nnodes,
                self.ntotal)
            raise RuntimeError(msg)

        # [angle, sc, sd, se, sf, omax, omin, mst, msc]
        self.data = np.zeros((self.ntimes, self.ntotal, 9), 'float32')

    #def build_dataframe(self):
        #"""creates a pandas dataframe"""
        #import pandas as pd
        #print(self.data_code)
        #headers = self.headers
        #column_names, column_values = self._build_dataframe_transient_header()
        #self.data_frame = pd.Panel(self.data, items=column_values,
                                   #major_axis=self.element_node, minor_axis=headers).to_frame()
        #self.data_frame.columns.names = column_names
        #self.data_frame.index.names = ['ElementID', 'Item']

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
                    for ieid, (eid, nid) in enumerate(self.element_node):
                        t1 = self.data[itime, ieid, :]
                        t2 = table.data[itime, ieid, :]
                        (angle1, sc1, sd1, se1, sf1, omax1, omin1, mst1, msc1) = t1
                        (angle2, sc2, sd2, se2, sf2, omax2, omin2, mst2, msc2) = t2
                        delta = t1 - t2
                        if not np.allclose([angle1, sc1, sd1, se1, sf1],
                                           [angle2, sc2, sd2, se2, sf2], atol=0.0001):
                        #if not np.array_equal(t1, t2):
                            msg += '%-4s  (%s, %s, %s, %s)\n      (%s, %s, %s, %s)\n  dt12=(%s, %s, %s, %s)\n' % (
                                eid,
                                sc1, sd1, se1, sf1,
                                sc2, sd2, se2, sf2,
                                delta[0], delta[1], delta[2], delta[3],)
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

    def add_sort1(self, dt, eid, grid, angle, sc, sd, se, sf, omax, omin, mst, msc):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.data[self.itime, self.itotal] = [angle, sc, sd, se, sf, omax, omin, mst, msc]
        self.element_node[self.itotal] = [eid, grid]
        #self.ielement += 1
        self.itotal += 1

    def get_stats(self, short=False) -> List[str]:
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        nnodes = self.element_node.shape[0]
        #ntotal = self.ntotal
        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i nnodes=%i\n'
                       % (self.__class__.__name__, ntimes, nelements, nnodes))
        else:
            msg.append('  type=%s nelements=%i nnodes=%i\n' % (self.__class__.__name__, nelements, nnodes))
        msg.append('  data: [ntimes, nnodes, 5] where 5=[%s]\n' % str(', '.join(self._get_headers())))
        msg.append('  element_node.shape = %s\n' % str(self.element_node.shape).replace('L', ''))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    @property
    def headers(self):
        return self._get_headers()

    def get_headers(self) -> List[str]:
        return self.headers

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        """
        '                                   S T R A I N S    I N   B E N D   E L E M E N T S        ( C B E N D )'
        '                        CIRC.'
        '   ELEMENT-ID  GRID END  ANG.   SXC           SXD           SXE           SXF           S-MAX         S-MIN         M.S.-T   M.S.-C'
        '0      6901    6901   A    0   4.372282E-16 -5.960465E-15  0.0           0.0           4.372282E-16 -5.960465E-15          '
        '               6902   B    0  -6.533992E-15  5.000001E-07 -5.000000E-13 -5.000000E-13  5.000001E-07 -5.000000E-13 -6.0E-01  6.0E+05'
        """
        #raise NotImplementedError('CBEND.stress/strain.real write_f06')
        msg_temp = _get_cbend_msg(self.is_stress, is_mag_phase, is_sort1)
        ntimes = self.data.shape[0]
        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]
        counter = cycle([0, 1])
        if self.is_sort1:
            if is_sort1:
                for itime in range(ntimes):
                    dt = self._times[itime]
                    if self.nonlinear_factor in (None, np.nan):
                        dt_line = ''
                    else:

                        dt_line = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
                    header[1] = dt_line
                    msg = header + msg_temp
                    f06_file.write('\n'.join(msg))

                    # [angle, sc, sd, se, sf]
                    angles = self.data[itime, :, 0]
                    scs = self.data[itime, :, 1]
                    sds = self.data[itime, :, 2]
                    ses = self.data[itime, :, 3]
                    sfs = self.data[itime, :, 4]

                    maxs = self.data[itime, :, 5]
                    mins = self.data[itime, :, 6]
                    msts = self.data[itime, :, 7]
                    mscs = self.data[itime, :, 8]
                    assert len(eids) == len(angles)
                    assert len(angles) > 0, angles
                    for i, eid, nid, anglei, sci, sdi, sei, sfi, maxi, mini, msti, msci in zip(counter, eids, nids, angles, scs, sds, ses, sfs, maxs, mins, msts, mscs):
                        assert isinstance(eid, integer_types), 'eid=%s type=%s' % (eid, type(eid))
                        [angle, sc, sd, se, sf, omax, omin] = write_floats_13e([anglei, sci, sdi, sei, sfi, maxi, mini])
                        [mst, msc] = write_floats_8p1e([msti, msci])

                        #f.write('                      28                  0.0          /  0.0                           0.0          /  0.0\n')

                        #'      ELEMENT-ID =    6901'
                        #'                         C O M P L E X   S T R E S S E S   I N   B E N D   E L E M E N T S   ( C B E N D ) '
                        #'                                                          (REAL/IMAGINARY)'
                        #'                        CIRC.      LOCATION         LOCATION         LOCATION         LOCATION'
                        #'   FREQUENCY   GRID END  ANG.         C                D                E                F'
                        #'0 0.0          6901   A    0     1.384767E+01     6.258920E-01    -1.217803E+01     1.043753E+00'
                        #'                                -4.615430E-01    -2.086098E-02     4.058937E-01    -3.478828E-02'
                        if i == 0:
                            f06_file.write(
                                '0 %9i%8i   A  %.2f %-13s %-13s %-13s %-13s %-13s %-13s %s %s\n'% (
                                    eid, nid, anglei,
                                    sc, sd, se, sf, omax, omin, mst, msc,
                                ))
                        else:
                            f06_file.write(
                                '  %9s%8i   B  %.2f %-13s %-13s %-13s %-13s %-13s %-13s %s %s\n'% (
                                    '', nid, anglei,
                                    sc, sd, se, sf, omax, omin, mst, msc,
                                ))
                    f06_file.write(page_stamp % page_num)
                    page_num += 1
            else:
                raise NotImplementedError('RealBendArray-sort2')
        else:
            raise NotImplementedError('RealBendArray-sort2')
        return page_num - 1

def _get_cbend_msg(is_stress, is_mag_phase, is_sort1):
    """get the header for the CBEND result"""

    if is_stress:
        stress_strain = '                                  S T R E S S E S   I N   B E N D   E L E M E N T S        ( C B E N D )'
    else:
        stress_strain = '                                   S T R A I N S    I N   B E N D   E L E M E N T S        ( C B E N D )'

    assert is_sort1 is True
    sort1 = '   ELEMENT-ID  GRID END  ANG.   SXC           SXD           SXE           SXF           S-MAX         S-MIN         M.S.-T   M.S.-C\n'
    msg = [
        stress_strain,
        '                        CIRC.',
        sort1,
    ]
    #'0      6901    6901   A    0   4.372282E-16 -5.960465E-15  0.0           0.0           4.372282E-16 -5.960465E-15          '
    #'               6902   B    0  -6.533992E-15  5.000001E-07 -5.000000E-13 -5.000000E-13  5.000001E-07 -5.000000E-13 -6.0E-01  6.0E+05'

    #if is_sort1:
        #msg.append('   ELEMENT-ID  GRID END  ANG.         C                D                E                F\n')
    #else:
        #msg.append('   FREQUENCY   GRID END  ANG.         C                D                E                F\n')
    return msg

class RealBendStressArray(RealBendArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealBendArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def _get_headers(self):
        return ['angle', 'sc', 'sd', 'se', 'sf', 'omax', 'omin', 'mst', 'msc']

class RealBendStrainArray(RealBendArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealBendArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)
        assert self.is_strain, self.stress_bits

    def _get_headers(self):
        return ['angle', 'sc', 'sd', 'se', 'sf', 'emax', 'emin', 'mst', 'msc']
