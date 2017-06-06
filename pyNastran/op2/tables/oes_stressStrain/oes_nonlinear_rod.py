from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import integer_types
from six.moves import range
import numpy as np
from numpy import zeros

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import OES_Object
from pyNastran.f06.f06_formatting import write_floats_13e, _eigenvalue_header
try:
    import pandas as pd
except ImportError:
    pass


class RealNonlinearRodArray(OES_Object): # 89-CRODNL, 92-CONRODNL
    """
    ::

      ELEMENT-ID =     102
                               N O N L I N E A R   S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )
        TIME          AXIAL STRESS         EQUIVALENT         TOTAL STRAIN       EFF. STRAIN          EFF. CREEP        LIN. TORSIONAL
                                             STRESS                             PLASTIC/NLELAST          STRAIN              STRESS
      2.000E-02        1.941367E+01        1.941367E+01        1.941367E-04        0.0                 0.0                 0.0
      3.000E-02        1.941367E+01        1.941367E+01        1.941367E-04        0.0                 0.0                 0.0
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=True)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        self.nelements = 0  # result specific

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
        headers = ['axial_stress', 'equiv_stress', 'total_strain',
                   'effective_plastic_creep_strain', 'effective_creep_strain',
                   'linear_torsional_stress']
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
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.nelements, dtype='int32')

        #[axial_stress, equiv_stress, total_strain, effective_plastic_creep_strain,
        # effective_creep_strain, linear_torsional_stress]
        self.data = zeros((self.ntimes, self.nelements, 6), dtype='float32')

    def build_dataframe(self):
        headers = self.get_headers()
        if self.nonlinear_factor is not None:
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
            self.data_frame.index.names = ['ElementID', 'Item']
        else:
            df1 = pd.DataFrame(self.element).T
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
            ntimes = self.data.shape[0]

            i = 0
            if self.is_sort1():
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
                raise NotImplementedError(self.is_sort2())
            if i > 0:
                print(msg)
                raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, axial_stress, equiv_stress, total_strain,
                  effective_plastic_creep_strain, effective_creep_strain, linear_torsional_stress):
        """unvectorized method for adding SORT1 transient data"""
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [
            axial_stress, equiv_stress, total_strain, effective_plastic_creep_strain,
            effective_creep_strain, linear_torsional_stress
        ]
        self.ielement += 1

    def get_stats(self, short=False):
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
            ntimes_word = '1'
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
        if is_sort1:
            msg = [
                '                         N O N L I N E A R   S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )\n',
                ' \n',
                '    ELEMENT-ID    AXIAL STRESS         EQUIVALENT         TOTAL STRAIN       EFF. STRAIN          EFF. CREEP        LIN. TORSIONAL\n',
                '                                         STRESS                             PLASTIC/NLELAST          STRAIN              STRESS\n'
            ]
        else:
            msg = [
                '                         N O N L I N E A R   S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )\n',
                ' \n',
                '    TIME          AXIAL STRESS         EQUIVALENT         TOTAL STRAIN       EFF. STRAIN          EFF. CREEP        LIN. TORSIONAL\n',
                '                                         STRESS                             PLASTIC/NLELAST          STRAIN              STRESS\n'
            ]

        if self.is_sort1():
            page_num = self._write_sort1_as_sort1(header, page_stamp, page_num, f, msg)
        else:
            raise NotImplementedError('RealNonlinearRodArray')
        return page_num

    def _write_sort1_as_sort1(self, header, page_stamp, page_num, f, msg_temp):
        ntimes = self.data.shape[0]

        eids = self.element
        #is_odd = False
        #nwrite = len(eids)

        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            axial = self.data[itime, :, 0]
            eqs = self.data[itime, :, 1]
            total = self.data[itime, :, 2]
            epcs = self.data[itime, :, 3]
            ecs = self.data[itime, :, 4]
            lts = self.data[itime, :, 5]

            #print "dt=%s axials=%s eqs=%s ts=%s epcs=%s ecs=%s lts=%s" %(dt,axial,eqs,ts,epcs,ecs,lts)
            #msgE[eid] = '      ELEMENT-ID = %8i\n' % (eid)
            #if eid not in msgT:
                #msgT[eid] = []
            #msgT[eid].append('  %9.3E       %13.6E       %13.6E       %13.6E       %13.6E       %13.6E       %13.6E\n' % (dt, axial, eqs, ts, epcs, ecs, lts))

            for eid, axiali, eqsi, totali, epcsi, ecsi, ltsi in zip(eids, axial, eqs, total, epcs, ecs, lts):
                ([saxial, seqs, stotal, sepcs, secs, slts]) = write_floats_13e(
                    [axiali, eqsi, totali, epcsi, ecsi, ltsi])
                f.write('  %8i       %-13s       %-13s       %-13s       %-13s       %-13s       %s\n' % (
                    eid, saxial, seqs, stotal, sepcs, secs, slts))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1
