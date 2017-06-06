from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import itervalues, integer_types
from math import isnan
from collections import OrderedDict
from numpy import zeros, empty, array_equal
import numpy as np

from pyNastran.op2.result_objects.op2_objects import ScalarObject
from pyNastran.f06.f06_formatting import get_key0, _eigenvalue_header, write_float_13e
try:
    import pandas as pd
except ImportError:
    pass


class RealStrainEnergyArray(ScalarObject):
    """
    ::

                                E L E M E N T   S T R A I N   E N E R G I E S

      ELEMENT-TYPE = QUAD4     * TOTAL ENERGY OF ALL ELEMENTS IN PROBLEM   = 9.817708E+08
      SUBCASE               1  * TOTAL ENERGY OF ALL ELEMENTS IN SET     1 = 4.192036E+08

         ELEMENT-ID   STRAIN-ENERGY  PERCENT OF TOTAL  STRAIN-ENERGY-DENSITY
                 12   2.291087E+07        2.3336            2.291087E+02
                 13   1.582968E+07        1.6124            1.055312E+02
                 14   6.576075E+07        6.6982            3.288037E+02
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific
        self.itime = None
        self.itotal2 = 0
        #self.element_name_count = OrderedDict()
        self.dt_temp = None

        if is_sort1:
            pass
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
        headers = [
            'strain_energy', 'percent', 'strain_energy_density'
        ]
        return headers

    def build(self):
        if self.is_built:
            return
        del self.dt_temp

        #print(self._ntotals)

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        self.ntotal = max(self._ntotals)
        #if max(self._ntotals) != min(self._ntotals):
            #raise RuntimeError('variable length in RealStrainEnergyArray')

        #self.names = []
        #self.nelements = self.ntotal // self.ntimes
        self.nelements = self.ntotal
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.itotal2 = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'
        self.build_data(dtype)

    def build_data(self, dtype):
        """actually performs the build step"""
        self._times = zeros(self.ntimes, dtype=dtype)
        #self.element = zeros(self.nelements, dtype='int32')
        self.element = zeros((self.ntimes, self.nelements), dtype='int32')
        #self.element_data_type = empty(self.nelements, dtype='|U8')

        #[energy, percent, density]
        assert isinstance(self.ntimes, integer_types), self.ntimes
        assert isinstance(self.ntotal, integer_types), self.ntotal
        self.data = zeros((self.ntimes, self.nelements, 3), dtype='float32')

    def build_dataframe(self):
        """
        major-axis - the axis

        mode              1     2   3
        freq              1.0   2.0 3.0
        ElementID Item
        1         T1
                  T2
                  ...

        major_axis / top = [
            [1, 2, 3],
            [1.0, 2.0, 3.0]
        ]
        minor_axis / headers = [ese, %, sed]
        name = mode
        """
        headers = self.get_headers()
        ntimes = self.element.shape[0]
        nelements = self.element.shape[1]
        if ntimes == 1:
            column_names, column_values = self._build_dataframe_transient_header()
            element = self.element.ravel()
            self.data_frame = pd.Panel(self.data, items=column_values,
                                       major_axis=element,
                                       minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
        else:
            nvalues = ntimes * nelements
            element = self.element.ravel()
            if self.nonlinear_factor is not None:
                column_names, column_values = self._build_dataframe_transient_header()
                #column_names = column_names[0]
                #column_values = column_values[0]

                column_values2 = []
                for value in column_values:
                    values2 = []
                    for valuei in value:
                        values = np.ones(nelements) * valuei
                        values2.append(values)
                    values3 = np.vstack(values2).ravel()
                    column_values2.append(values3)
                df1 = pd.DataFrame(column_values2).T
                df1.columns = column_names

                df2 = pd.DataFrame(element)
                df2.columns = ['ElementID']

                dfs = [df2]
                for i, header in enumerate(headers):
                    df = pd.DataFrame(self.data[:, :, i].ravel())
                    df.columns = [header]
                    dfs.append(df)
                self.data_frame = df1.join(dfs)
                #self.data_frame.columns.names = column_names

            # remove empty rows
            self.data_frame = self.data_frame[self.data_frame.ElementID != 0]

    def __eq__(self, table):
        return self.assert_equal(table)

    def assert_equal(self, table, rtol=1.e-5, atol=1.e-8):
        self._eq_header(table)
        assert self.is_sort1() == table.is_sort1()

        if not np.array_equal(self.element, table.element):
            assert self.element.shape == table.element.shape, 'element shape=%s table.shape=%s' % (self.element.shape, table.element.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += 'itime: eid1 eid2\n'

            i = 0
            for itime in range(self.ntimes):
                for eid1, eid2 in zip(self.element[itime, :], table.element[itime, :]):
                    msg += '%s: %s %s\n' % (itime, eid1, eid2)
                    if eid1 != eid2 and np.isnan(eid1):
                        i += 1
                        if i > 10:
                            print(msg)
                        raise ValueError(msg)
            if i > 0:
                raise ValueError(msg)


        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, eid in enumerate(self.element[itime, :]):
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (energyi1, percenti1, densityi1) = t1
                    (energyi2, percenti2, densityi2) = t2

                    if np.isnan(densityi1) or not np.isfinite(densityi1):
                        if not np.array_equal(t1[:2], t2[:2]):
                            msg += (
                                '%s (%s, %s)\n'
                                '%s (%s, %s)\n' % (
                                    eid, energyi1, percenti1,
                                    ' ' * len(str(eid)),
                                    energyi2, percenti2,
                                ))
                            i += 1
                            if i > 10:
                                print(msg)
                                raise ValueError(msg)
                    elif not np.array_equal(t1, t2):
                        msg += (
                            '%s (%s, %s, %s)\n'
                            '%s (%s, %s, %s)\n' % (
                                eid, energyi1, percenti1, densityi1,
                                ' ' * len(str(eid)),
                                energyi2, percenti2, densityi2,
                            ))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, energyi, percenti, densityi):
        """unvectorized method for adding SORT1 transient data"""
        #itime = self.itime // self.nelement_types
        itime = self.itime
        self._times[itime] = dt
        try:
            self.element[itime, self.ielement] = eid
            #self.element_data_type[self.ielement] = etype
            self.data[itime, self.ielement, :] = [energyi, percenti, densityi]
        except IndexError:
            print('RealStrainEnergyArray', dt, eid, energyi, percenti, densityi)
            raise
        self.ielement += 1
        self.itotal += 1

    def get_stats(self, short=False):
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
            ntimes_word = '1'
        headers = self.get_headers()
        n = len(headers)
        msg.append('  element: [%s, nelements]; eid=100000000 -> total\n' % (ntimes_word))
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        #msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        """
        '      EIGENVALUE =  2.005177E+05'
        '      CYCLES =  7.126832E+01'
        '                                           E L E M E N T   S T R A I N   E N E R G I E S'
        ' '
        '                ELEMENT-TYPE = TETRA               * TOTAL ENERGY OF ALL ELEMENTS IN PROBLEM     =   1.002589E+05'
        '                   MODE               1            * TOTAL ENERGY OF ALL ELEMENTS IN SET      -1 =   1.002589E+05'
        '0'
        '                                    ELEMENT-ID          STRAIN-ENERGY           PERCENT OF TOTAL    STRAIN-ENERGY-DENSITY'
        '                                             4          3.247409E+00                 0.0032              1.948445E+01'
        '                                             5          3.977916E+00                 0.0040              2.386749E+01'
        ''
        '                        TYPE = TETRA    SUBTOTAL        7.225325E+00                 0.0072'
        """
        msg_temp = (
            '                                           E L E M E N T   S T R A I N   E N E R G I E S\n'
            ' \n'
            '                ELEMENT-TYPE = %s               * TOTAL ENERGY OF ALL ELEMENTS IN PROBLEM     =   %s\n'
            '                   MODE        %8i            * TOTAL ENERGY OF ALL ELEMENTS IN SET      -1 =   %s\n'
            '0\n'
            '                                    ELEMENT-ID          STRAIN-ENERGY           PERCENT OF TOTAL    STRAIN-ENERGY-DENSITY\n'
        )
        ntimes = self.data.shape[0]

        #etype = self.element_data_type
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            total_energy = 0.
            total_set_energy = 0.

            eids = self.element[itime, :]
            # energy, percent, density
            energy = self.data[itime, :, 0]
            percent = self.data[itime, :, 1]
            density = self.data[itime, :, 2]

            itotal = np.where(eids == 100000000)[0][0]
            #total_energy = self.data[:, :, 0].sum()
            #total_set_energy = energy.sum()
            #total_set_energy = energy[itotal]
            #total_percent = percent.sum()


            msg_temp2 = [msg_temp % (self.element_name, total_energy, itime + 1, total_set_energy)]
            f.write(''.join(header + msg_temp2))


            fmt1 = ' ' * 36 + '%10i         %-13s                 %.4f             %s\n'
            fmt2 = '\n                        TYPE = %-8s SUBTOTAL       %13s                 %.4f\n'

            for (eid, energyi, percenti, densityi) in zip(eids, energy, percent, density):
                senergyi = write_float_13e(energyi)
                sdensityi = write_float_13e(densityi)
                # ELEMENT-ID    STRAIN-ENERGY   PERCENT OF TOTAL  STRAIN-ENERGY-DENSITY
                #          1   -8.307121E-12         0.0052           -2.886861E-12
                if eid == 100000000:
                    f.write(fmt2 % (self.element_name, senergyi, percenti))
                    break
                f.write(fmt1 % (
                    eid, senergyi, percenti, sdensityi))
            f.write(page_stamp % page_num)
            page_num += 1
            break
        return page_num - 1


class ComplexStrainEnergyArray(ScalarObject):
    """
    ::

            FREQUENCY =  2.000000E+03
                                     E L E M E N T   S T R A I N   E N E R G I E S   ( A V E R A G E )

                      ELEMENT-TYPE = QUAD4               * TOTAL ENERGY OF ALL ELEMENTS IN PROBLEM     =   1.611784E-08
                      SUBCASE               1            * TOTAL ENERGY OF ALL ELEMENTS IN SET      -1 =   1.611784E-08
      0
                             ELEMENT-ID       STRAIN-ENERGY (MAG/PHASE)               PERCENT OF TOTAL    STRAIN-ENERGY-DENSITY
                                      5       2.027844E-10 /   0.0                         1.2581              2.027844E-09
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific
        self.itime = None
        self.itotal2 = 0
        #self.element_name_count = OrderedDict()
        self.dt_temp = None

        if is_sort1:
            pass
        else:
            raise NotImplementedError('SORT2')

    def is_real(self):
        return False

    def is_complex(self):
        return True

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        headers = [
            'strain_energy', 'percent', 'strain_energy_density'
        ]
        return headers

    def build(self):
        if self.is_built:
            return
        del self.dt_temp

        #print(self._ntotals)

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        self.ntotal = max(self._ntotals)
        #if max(self._ntotals) != min(self._ntotals):
            #raise RuntimeError('variable length in RealStrainEnergyArray')

        #self.names = []
        #self.nelements = self.ntotal // self.ntimes
        self.nelements = self.ntotal
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.itotal2 = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'
        self.build_data(dtype)

    def build_data(self, dtype):
        """actually performs the build step"""
        self._times = zeros(self.ntimes, dtype=dtype)
        #self.element = zeros(self.nelements, dtype='int32')
        self.element = zeros((self.ntimes, self.nelements), dtype='int32')
        #self.element_data_type = empty(self.nelements, dtype='|U8')

        #[energy, percent, density]
        assert isinstance(self.ntimes, integer_types), self.ntimes
        assert isinstance(self.ntotal, integer_types), self.ntotal
        self.data = zeros((self.ntimes, self.nelements, 4), dtype='float32')

    #def build_dataframe(self):
        #"""
        #major-axis - the axis

        #mode              1     2   3
        #freq              1.0   2.0 3.0
        #ElementID Item
        #1         T1
                  #T2
                  #...

        #major_axis / top = [
            #[1, 2, 3],
            #[1.0, 2.0, 3.0]
        #]
        #minor_axis / headers = [ese, %, sed]
        #name = mode
        #"""
        #headers = self.get_headers()
        #ntimes = self.element.shape[0]
        #nelements = self.element.shape[1]
        #if ntimes == 1:
            #column_names, column_values = self._build_dataframe_transient_header()
            #element = self.element.ravel()
            #self.data_frame = pd.Panel(self.data, items=column_values,
                                       #major_axis=element,
                                       #minor_axis=headers).to_frame()
            #self.data_frame.columns.names = column_names
        #else:
            #nvalues = ntimes * nelements
            #element = self.element.ravel()
            #if self.nonlinear_factor is not None:
                #column_names, column_values = self._build_dataframe_transient_header()
                ##column_names = column_names[0]
                ##column_values = column_values[0]

                #column_values2 = []
                #for value in column_values:
                    #values2 = []
                    #for valuei in value:
                        #values = np.ones(nelements) * valuei
                        #values2.append(values)
                    #values3 = np.vstack(values2).ravel()
                    #column_values2.append(values3)
                #df1 = pd.DataFrame(column_values2).T
                #df1.columns = column_names

                #df2 = pd.DataFrame(element)
                #df2.columns = ['ElementID']

                #dfs = [df2]
                #for i, header in enumerate(headers):
                    #df = pd.DataFrame(self.data[:, :, i].ravel())
                    #df.columns = [header]
                    #dfs.append(df)
                #self.data_frame = df1.join(dfs)
                ##self.data_frame.columns.names = column_names

            ## remove empty rows
            #self.data_frame = self.data_frame[self.data_frame.ElementID != 0]

    def __eq__(self, table):
        return self.assert_equal(table)

    def assert_equal(self, table, rtol=1.e-5, atol=1.e-8):
        self._eq_header(table)
        assert self.is_sort1() == table.is_sort1()

        if not np.array_equal(self.element, table.element):
            assert self.element.shape == table.element.shape, 'element shape=%s table.shape=%s' % (self.element.shape, table.element.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += 'itime: eid1 eid2\n'

            i = 0
            for itime in range(self.ntimes):
                for eid1, eid2 in zip(self.element[itime, :], table.element[itime, :]):
                    msg += '%s: %s %s\n' % (itime, eid1, eid2)
                    if eid1 != eid2 and np.isnan(eid1):
                        i += 1
                        if i > 10:
                            print(msg)
                        raise ValueError(msg)
            if i > 0:
                raise ValueError(msg)


        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, eid in enumerate(self.element[itime, :]):
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (energyi1r, engery1i, percenti1, densityi1) = t1
                    (energyi2r, engery2i, percenti2, densityi2) = t2
                    print(t1, t2)
                    if np.isnan(densityi1) or not np.isfinite(densityi1):
                        if not np.array_equal(t1[:2], t2[:2]):
                            msg += (
                                '%s (%s+%si, %s)\n'
                                '%s (%s+%si, %s)\n' % (
                                    eid, energyi1r, engery1i, percenti1,
                                    ' ' * len(str(eid)),
                                    energyi2r, engery2i, percenti2,
                                ))
                            i += 1
                            if i > 10:
                                print(msg)
                                raise ValueError(msg)
                    elif not np.array_equal(t1, t2):
                        msg += (
                            '%s (%s+%si, %s, %s)\n'
                            '%s (%s+%si, %s, %s)\n' % (
                                eid, energyi1r, engery1i, percenti1, densityi1,
                                ' ' * len(str(eid)),
                                energyi2r, engery2i, percenti2, densityi2,
                            ))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, energyr, energyi, percenti, densityi):
        """unvectorized method for adding SORT1 transient data"""
        #itime = self.itime // self.nelement_types
        itime = self.itime
        self._times[itime] = dt
        try:
            self.element[itime, self.ielement] = eid
            #self.element_data_type[self.ielement] = etype
            self.data[itime, self.ielement, :] = [energyr, energyi, percenti, densityi]
        except IndexError:
            print('ComplexStrainEnergyArray', dt, eid, energyr, energyi, percenti, densityi)
            raise
        self.ielement += 1
        self.itotal += 1

    def get_stats(self, short=False):
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
            ntimes_word = '1'
        headers = self.get_headers()
        n = len(headers)
        msg.append('  element: [%s, nelements]; eid=100000000 -> total\n' % (ntimes_word))
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        #msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = (
            '                               E L E M E N T   S T R A I N   E N E R G I E S   ( A V E R A G E )                 \n'
            ' \n'
            '                ELEMENT-TYPE = %-5s               * TOTAL ENERGY OF ALL ELEMENTS IN PROBLEM     =   %s\n'
            '                SUBCASE               1            * TOTAL ENERGY OF ALL ELEMENTS IN SET      -1 =   %s\n'
            '0\n'
            '                       ELEMENT-ID       STRAIN-ENERGY (MAG/PHASE)               PERCENT OF TOTAL    STRAIN-ENERGY-DENSITY\n'
            #'                                5       2.027844E-10 /   0.0                         1.2581              2.027844E-09'
        )
        ntimes = self.data.shape[0]

        #etype = self.element_data_type
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            total_energy = 0.
            total_set_energy = 0.

            eids = self.element[itime, :]
            # energyr, energyi, percent, density
            energyr = self.data[itime, :, 0]
            energyi = self.data[itime, :, 1]
            percent = self.data[itime, :, 2]
            density = self.data[itime, :, 3]

            #total_energy = self.data[:, :, 0].sum()
            #total_set_energy = energy.sum()
            #total_set_energy = energy[itotal]
            #total_percent = percent.sum()


            msg_temp2 = [msg_temp % (self.element_name, total_energy, total_set_energy)]
            f.write(''.join(header + msg_temp2))


            fmt1 = ' ' * 23 + '%10i      %-13s /  %-13s               %7.4f             %s\n'
            fmt2 = '\n                        TYPE = %-8s SUBTOTAL       %13s                 %.4f\n'

            for (eid, energyri, energyii, percenti, densityi) in zip(eids, energyr, energyi, percent, density):
                senergyr = write_float_13e(energyri)
                senergyi = write_float_13e(energyii)
                sdensityi = write_float_13e(densityi)
                # ELEMENT-ID    STRAIN-ENERGY   PERCENT OF TOTAL  STRAIN-ENERGY-DENSITY
                #          1   -8.307121E-12         0.0052           -2.886861E-12
                #if eid == 100000000:
                    #f.write(fmt2 % (self.element_name, senergyi, percenti))
                    #break
                f.write(fmt1 % (
                    eid, senergyr, senergyi, percenti, sdensityi))
            f.write(page_stamp % page_num)
            page_num += 1
            #break
        return page_num - 1
