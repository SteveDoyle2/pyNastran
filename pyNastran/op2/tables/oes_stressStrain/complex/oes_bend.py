from itertools import cycle
from typing import List

import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import (
    StressObject, StrainObject, OES_Object)
from pyNastran.f06.f06_formatting import write_imag_floats_13e


class ComplexBendArray(OES_Object):
    """
    Common class for:
     - ComplexBendStressArray
     - ComplexBendStrainArray
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
        if self.is_built:
            return
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

        # [angle, sc, sd, se, sf]
        self.data = np.zeros((self.ntimes, self.ntotal, 5), 'complex64')

    def build_dataframe(self):
        """creates a pandas dataframe"""
        #Freq                                   0.0                  2.5
        #ElementID NodeID Item
        #6901      6901   angle   0.000000+0.000000j   0.000000+0.000000j
        #                 sc     13.847674-0.461543j  13.855294-0.462052j
        #                 sd      0.625892-0.020861j   0.623742-0.020717j
        #                 se    -12.178029+0.405894j -12.185331+0.406381j
        #                 sf      1.043753-0.034788j   1.046222-0.034953j
        #          6904   angle   0.000000+0.000000j   0.000000+0.000000j
        #                 sc     -1.660571-0.416504j  -1.663256-0.416978j
        #                 sd     -2.790551+0.024178j  -2.789738+0.024356j
        #                 se      0.627616+0.450933j   0.629571+0.451455j
        #                 sf      1.757596+0.010251j   1.756053+0.010121j
        #6902      6901   angle   0.000000+0.000000j   0.000000+0.000000j
        headers = self.headers
        column_names, column_values = self._build_dataframe_transient_header()
        self.data_frame = self._build_pandas_transient_element_node(
            column_values, column_names,
            headers, self.element_node, self.data)

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
                        (angle1, sc1, sd1, se1, sf1) = t1
                        (angle2, sc2, sd2, se2, sf2) = t2
                        d = t1 - t2
                        if not np.allclose([angle1.real, sc1.real, sc1.imag, sd1.real, sd1.imag, se1.real, se1.imag, sf1.real, sf1.imag, ],
                                           [angle2.real, sc2.real, sc2.imag, sd2.real, sd2.imag, se2.real, se2.imag, sf2.real, sf2.imag, ], atol=0.0001):
                        #if not np.array_equal(t1, t2):
                            msg += '%-4s  (%s, %sj, %s, %sj)\n      (%s, %sj, %s, %sj)\n  dt12=(%s, %sj, %s, %sj)\n' % (
                                eid,
                                sc1.real, sc1.imag, sd1.real, sd1.imag,
                                sc2.real, sc2.imag, sd2.real, sd2.imag,
                                d[0].real, d[0].imag, d[1].real, d[1].imag,)
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

    def add_sort1(self, dt, eid, grid, angle, sc, sd, se, sf):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.data[self.itime, self.itotal] = [angle, sc, sd, se, sf]
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
            msg.append('  type=%s ntimes=%i nelements=%i nnodes=%i; table_name=%r\n' % (
                self.__class__.__name__, ntimes, nelements, nnodes, self.table_name))
        else:
            msg.append('  type=%s nelements=%i nnodes=%i; table_name=%r\n' % (
                self.__class__.__name__, nelements, nnodes, self.table_name))
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
              ELEMENT-ID =    6901'
                                 C O M P L E X   S T R E S S E S   I N   B E N D   E L E M E N T S   ( C B E N D )
                                                                  (REAL/IMAGINARY)
                                CIRC.      LOCATION         LOCATION         LOCATION         LOCATION
           FREQUENCY   GRID END  ANG.         C                D                E                F
        '0 0.0          6901   A    0     1.384767E+01     6.258920E-01    -1.217803E+01     1.043753E+00
        '                                -4.615430E-01    -2.086098E-02     4.058937E-01    -3.478828E-02
        """
        msg_temp = _get_cbend_msg(is_mag_phase, is_sort1)
        ntimes = self.data.shape[0]
        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]
        counter = cycle([0, 1])
        if self.is_sort1:
            if is_sort1:
                for itime in range(ntimes):
                    dt = self._times[itime]

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
                    assert len(eids) == len(angles)
                    assert len(angles) > 0, angles
                    for i, eid, nid, anglei, sci, sdi, sei, sfi in zip(counter, eids, nids, angles, scs, sds, ses, sfs):
                        assert isinstance(eid, integer_types), 'eid=%s type=%s' % (eid, type(eid))
                        [sc_real, sc_imag,
                         sd_real, sd_imag,
                         se_real, se_imag,
                         sf_real, sf_imag,] = write_imag_floats_13e([sci, sdi, sei, sfi], is_mag_phase)

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
                                '0%12i %8i  A  %.2f %-13s    %-13s    %-13s    %s\n'
                                ' %12s %8s          %-13s    %-13s    %-13s    %s\n'% (
                                    eid, nid, anglei.real,
                                    sc_real, sd_real, se_real, sf_real,
                                    '', '',
                                    sc_imag, sd_imag, se_imag, sf_imag,
                                ))
                        else:
                            f06_file.write(
                                '0%12s %8i  B  %.2f %-13s    %-13s    %-13s    %s\n'
                                ' %12s %8s          %-13s    %-13s    %-13s    %s\n'% (
                                    '', nid, anglei.real,
                                    sc_real, sd_real, se_real, sf_real,
                                    '', '',
                                    sc_imag, sd_imag, se_imag, sf_imag,
                                ))
                    f06_file.write(page_stamp % page_num)
                    page_num += 1
            else:
                raise NotImplementedError('ComplexBendArray-sort2')
        else:
            raise NotImplementedError('ComplexBendArray-sort2')
        return page_num - 1

def _get_cbend_msg(is_mag_phase, is_sort1):
    """get the header for the CBEND result"""
    if is_mag_phase:
        raise NotImplementedError()
    else:
        realimag_magphase = '                                                          (REAL/IMAGINARY)'

    msg = [
        '                         C O M P L E X   S T R E S S E S   I N   B E N D   E L E M E N T S   ( C B E N D ) ',
        realimag_magphase,
        '                        CIRC.      LOCATION         LOCATION         LOCATION         LOCATION',
    ]
    if is_sort1:
        msg.append('   ELEMENT-ID  GRID END  ANG.         C                D                E                F\n')
    else:
        msg.append('   FREQUENCY   GRID END  ANG.         C                D                E                F\n')
    return msg

class ComplexBendStressArray(ComplexBendArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexBendArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def _get_headers(self):
        return ['angle', 'sc', 'sd', 'se', 'sf']

class ComplexBendStrainArray(ComplexBendArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexBendArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)
        assert self.is_strain, self.stress_bits

    def _get_headers(self):
        return ['angle', 'sc', 'sd', 'se', 'sf']
