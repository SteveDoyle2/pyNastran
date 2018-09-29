from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import integer_types
from itertools import count
import numpy as np
from numpy import zeros, searchsorted, ravel
ints = (int, np.int32)

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, OES_Object
from pyNastran.f06.f06_formatting import write_floats_13e, _eigenvalue_header
try:
    import pandas as pd  # type: ignore
except ImportError:
    pass


class RealBush1DStressArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=True)
        #self.code = [self.format_code, self.sort_code, self.s_code]
        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.ielement = 0
        self.nelements = 0  # result specific

    @property
    def is_stress(self):
        return True

    @property
    def is_real(self):
        return True

    @property
    def is_complex(self):
        return False

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def _get_msgs(self):
        words = [
            #'      ELEMENT-ID =     104'
            '                    S T R E S S E S   ( F O R C E S )   I N   B U S H 1 D   E L E M E N T S   ( C B U S H 1 D )\n',
            ' \n',
            '                        AXIAL          AXIAL          AXIAL       AXIAL         AXIAL         PLASTIC\n',
            '        TIME            FORCE       DISPLACEMENT    VELOCITY      STRESS        STRAIN        STRAIN        STATUS\n',
            #'    2.000000E-02      1.960396E+01  1.960396E-04  1.940792E-02  1.960396E+01  1.960396E-04  0.000000E+00    \n',
        ]
        return words
        # raise NotImplementedError('%s needs to implement _get_msgs' % self.__class__.__name__)

    def get_headers(self):
        headers = ['element_force', 'axial_displacement', 'axial_velocity',
                   'axial_stress', 'axial_strain', 'plastic_strain']
        return headers

    def build(self):
        """sizes the vectorized attributes of the RealBush1DStressArray"""
        if self.is_built:
            return
        #print("self.ielement =", self.ielement)
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal

        if self.element_type == 40:
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
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.ntotal, dtype='int32')
        self.is_failed = zeros((self.ntimes, self.ntotal, 1), dtype='int32')

        # [element_force, axial_displacement, axial_velocity, axial_stress, axial_strain, plastic_strain, is_failed]
        self.data = zeros((self.ntimes, self.ntotal, 6), dtype='float32')

    def build_dataframe(self):
        """creates a pandas dataframe"""
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
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)

        i = 0
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            if self.is_sort1:
                for itime in range(ntimes):
                    for ieid, eid, in enumerate(self.element):
                        t1 = self.data[itime, ieid, :]
                        t2 = table.data[itime, ieid, :]
                        #i_not_nan = np.isnp.where(t1 != np.nan)[0]
                        i_not_nan = np.isfinite(t1)
                        (axial_stress1, equiv_stress1, total_strain1, effective_plastic_creep_strain1, effective_creep_strain1, linear_torsional_stress1) = t1
                        (axial_stress2, equiv_stress2, total_strain2, effective_plastic_creep_strain2, effective_creep_strain2, linear_torsional_stress2) = t2
                        if not np.allclose(t1[i_not_nan], t2[i_not_nan]):
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

    def add_sort1(self, dt, eid, element_force, axial_displacement, axial_velocity,
                  axial_stress, axial_strain, plastic_strain, is_failed):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, ints)
        assert isinstance(eid, int) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        # pyNastran_examples\move_tpl\ar29scb1.op2
        #print('dt=%s eid=%s force=%s' % (dt, eid, element_force))
        #print('element.shape=%s' % self.element.shape)
        #print('data.shape=%s' % str(self.data.shape))
        #print('times.shape=%s' % self._times.shape)
        #print('itime=%s ielement=%s itotal=%s' % (self.itime, self.itotal, self.ielement))
        self._times[self.itime] = dt
        self.element[self.itotal] = eid
        self.is_failed[self.itime, self.itotal, 0] = is_failed
        self.data[self.itime, self.itotal, :] = [
            element_force, axial_displacement, axial_velocity,
            axial_stress, axial_strain, plastic_strain]
        self.itotal += 1
        self.ielement += 1

    def get_stats(self, short=False):
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
            ntimes_word = '1'
        headers = self.get_headers()

        n = len(headers)
        assert n == self.data.shape[2], 'nheaders=%s shape=%s' % (n, str(self.data.shape))
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  element.shape = %s\n' % str(self.element.shape).replace('L', ''))
        msg.append('  is_failed.shape = %s\n' % str(self.is_failed.shape).replace('L', ''))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element)  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        ind = ravel([searchsorted(self.element == eid) for eid in eids])
        return ind

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg = self._get_msgs()
        ntimes = self.data.shape[0]
        eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg))

            #[element_force, axial_displacement, axial_velocity, axial_stress, axial_strain, plastic_strain, is_failed]
            element_force = self.data[itime, :, 0]
            axial_displacement = self.data[itime, :, 1]
            axial_velocity = self.data[itime, :, 2]
            axial_stress = self.data[itime, :, 3]
            axial_strain = self.data[itime, :, 4]
            plastic_strain = self.data[itime, :, 5]
            is_failed = self.is_failed[itime, :, 0]

            for (i, eid, element_forcei, axial_displacementi, axial_velocityi, axial_stressi, axial_straini, plastic_straini, is_failedi) in zip(
                count(), eids, element_force, axial_displacement, axial_velocity, axial_stress, axial_strain, plastic_strain, is_failed):

                vals = [element_forcei, axial_displacementi, axial_velocityi, axial_stressi, axial_straini, plastic_straini, is_failedi]
                vals2 = write_floats_13e(vals)
                [element_forcei, axial_displacementi, axial_velocityi, axial_stressi, axial_straini, plastic_straini, is_failedi] = vals2
                f06_file.write(
                    '0%8i   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                    % (eid, element_forcei, axial_displacementi, axial_velocityi, axial_stressi,
                       axial_straini, plastic_straini, is_failedi))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        if self.nonlinear_factor is None:
            page_num -= 1
        return page_num
