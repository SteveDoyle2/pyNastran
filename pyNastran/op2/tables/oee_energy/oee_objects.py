from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import string_types, binary_type
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.result_objects.op2_objects import BaseElement
from pyNastran.f06.f06_formatting import _eigenvalue_header, write_float_13e
from pyNastran.op2.op2_interface.write_utils import set_table3_field

SORT2_TABLE_NAME_MAP = {
    'ONRGY2' : 'ONRGY1',
}
class RealStrainEnergyArray(BaseElement):
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
        BaseElement.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific
        self.itime = None
        self.itotal2 = 0
        #self.element_name_count = OrderedDict()
        self.dt_temp = None

        #if is_sort1:
            #pass
        #else:
            #raise NotImplementedError('SORT2')

    @property
    def is_real(self):
        return True

    @property
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
        """sizes the vectorized attributes of the RealStrainEnergyArray"""
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
        self._times = np.zeros(self.ntimes, dtype=dtype)
        #self.element = zeros(self.nelements, dtype='int32')
        #if dtype in 'DMIG':
        #print(self.element_name, self.element_type)
        if self.element_name == 'DMIG':
            self.element = np.zeros((self.ntimes, self.nelements), dtype='|U8')
        else:
            self.element = np.zeros((self.ntimes, self.nelements), dtype='int32')
        #self.element_data_type = empty(self.nelements, dtype='|U8')

        #[energy, percent, density]
        assert isinstance(self.ntimes, integer_types), self.ntimes
        assert isinstance(self.ntotal, integer_types), self.ntotal
        self.data = np.zeros((self.ntimes, self.nelements, 3), dtype='float32')

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
        import pandas as pd
        #print(''.join(self.get_stats()))
        #print(self.element)
        #print(self.data)
        headers = self.get_headers()
        ntimes = self.element.shape[0]
        nelements = self.element.shape[1]

        element = self.element.ravel()
        if element.dtype is np.dtype(np.int32):
            compare = 0
        else:
            # unicode
            #value = value.tolist()

            element = np.asarray(element, dtype='|U8')
            compare = ''

        #print('ntimes=%s' % ntimes)
        if ntimes == 1:
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values,
                                       major_axis=element,
                                       minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
        else:
            # we can get into this in a linear case
            # F:\work\pyNastran\examples\Dropbox\move_tpl\setp04.op2

            nvalues = ntimes * nelements

            #if self.nonlinear_factor not in (None, np.nan):
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
            try:
                self.data_frame.columns.names = column_names
            except ValueError:
                #print('headers =', headers)
                print('self.cannot apply column_names=%s to RealStrainEnergyArray: %r' % (
                    column_names, self.element_name))

            # remove empty rows
            assert self.data_frame is not None
            self.data_frame = self.data_frame[self.data_frame.ElementID != compare]

    def __eq__(self, table):
        return self.assert_equal(table)

    def assert_equal(self, table, rtol=1.e-5, atol=1.e-8):
        self._eq_header(table)
        assert self.is_sort1 == table.is_sort1

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
        assert (isinstance(eid, int) and eid > 0) or isinstance(eid, binary_type), 'dt=%s eid=%s' % (dt, eid)
        itime = self.itime
        self._times[itime] = dt
        self.element[itime, self.ielement] = eid

        if self.element_name == 'DMIG':
            if not np.isnan(densityi):
                raise RuntimeError(
                    'RealStrainEnergyArray: itime=%s ielement=%s; '
                    'dt=%s eid=%s energyi=%s percenti=%s densityi=%s' % (
                        self.itime, self.ielement, dt, eid, energyi, percenti, densityi))
            self.data[itime, self.ielement, :] = [energyi, percenti, np.nan]
        else:
            try:
                #self.element_data_type[self.ielement] = etype
                self.data[itime, self.ielement, :] = [energyi, percenti, densityi]
            except (ValueError, IndexError):
                print('RealStrainEnergyArray: itime=%s ielement=%s; '
                      'dt=%s eid=%s energyi=%s percenti=%s densityi=%s' % (
                    self.itime, self.ielement, dt, eid, energyi, percenti, densityi))
                raise
        self.ielement += 1
        self.itotal += 1

    def finalize(self):
        self.set_as_sort1()

    def set_as_sort1(self):
        """changes the table into SORT1"""
        if self.is_sort1:
            return
        try:
            analysis_method = self.analysis_method
        except AttributeError:
            print(self.code_information())
            raise
        #print(self.get_stats())
        #print(self.node_gridtype)
        #print(self.data.shape)
        #aaa
        self.sort_method = 1
        self.sort_bits[1] = 0
        bit0, bit1, bit2 = self.sort_bits
        self.table_name = SORT2_TABLE_NAME_MAP[self.table_name]
        self.sort_code = bit0 + 2*bit1 + 4*bit2
        #print(self.code_information())
        assert self.is_sort1
        if analysis_method != 'N/A':
            self.data_names[0] = analysis_method
            #print(self.table_name_str, analysis_method, self._times)
            setattr(self, self.analysis_method + 's', self._times)
        del self.analysis_method

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
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s element_name=%r ntimes=%i nelements=%i\n'
                       % (self.__class__.__name__, self.element_name, ntimes, nelements))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s element_name=%r nelements=%i\n'
                       % (self.__class__.__name__, self.element_name, nelements))
            ntimes_word = '1'
        headers = self.get_headers()
        n = len(headers)
        msg.append('  element: [%s, nelements]; eid=100000000 -> total\n' % (ntimes_word))
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        #msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        # '      EIGENVALUE =  2.005177E+05'
        # '      CYCLES =  7.126832E+01'
        # '                                           E L E M E N T   S T R A I N   E N E R G I E S'
        # ' '
        # '                ELEMENT-TYPE = TETRA               * TOTAL ENERGY OF ALL ELEMENTS IN PROBLEM     =   1.002589E+05'
        # '                   MODE               1            * TOTAL ENERGY OF ALL ELEMENTS IN SET      -1 =   1.002589E+05'
        # '0'
        # '                                    ELEMENT-ID          STRAIN-ENERGY           PERCENT OF TOTAL    STRAIN-ENERGY-DENSITY'
        # '                                             4          3.247409E+00                 0.0032              1.948445E+01'
        # '                                             5          3.977916E+00                 0.0040              2.386749E+01'
        # ''
        # '                        TYPE = TETRA    SUBTOTAL        7.225325E+00                 0.0072'

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

            #itotal = np.where(eids == 100000000)[0][0]
            #total_energy = self.data[:, :, 0].sum()
            #total_set_energy = energy.sum()
            #total_set_energy = energy[itotal]
            #total_percent = percent.sum()


            msg_temp2 = [msg_temp % (self.element_name, total_energy, itime + 1, total_set_energy)]
            f06_file.write(''.join(header + msg_temp2))


            fmt1 = ' ' * 36 + '%10s         %-13s                 %.4f             %s\n'
            fmt1_nan = ' ' * 36 + '%10s         %-13s                 %.4f             %s\n'
            fmt2 = '\n                        TYPE = %-8s SUBTOTAL       %13s                 %.4f\n'

            for (eid, energyi, percenti, densityi) in zip(eids, energy, percent, density):
                senergyi = write_float_13e(energyi)
                sdensityi = write_float_13e(densityi)
                # ELEMENT-ID    STRAIN-ENERGY   PERCENT OF TOTAL  STRAIN-ENERGY-DENSITY
                #          1   -8.307121E-12         0.0052           -2.886861E-12
                if eid == 100000000:
                    f06_file.write(fmt2 % (self.element_name, senergyi, percenti))
                    break
                try:
                    f06_file.write(fmt1 % (eid, senergyi, percenti, sdensityi))
                except TypeError:
                    #print('eid = %r; type=%s' % (eid, type(eid)))
                    #print('senergyi = %r; type=%s' % (senergyi, type(senergyi)))
                    #print('percenti = %r; type=%s' % (percenti, type(percenti)))
                    #print('sdensityi = %r; type=%s' % (sdensityi, type(sdensityi)))
                    assert np.isnan(sdensityi), 'eid=%s sdensityi=%s' % (eid, sdensityi)
                    f06_file.write(fmt1_nan % (eid, senergyi, percenti, ''))
                    #if 0:
                        #print('senergyi = %r; type=%s' % (senergyi, type(senergyi)))
                        #print('percenti = %r; type=%s' % (percenti, type(percenti)))
                        #print('sdensityi = %r; type=%s' % (sdensityi, type(sdensityi)))
                        #msg = fmt1 % (eid, senergyi, percenti, sdensityi)
                        #raise TypeError(msg)
                    #raise RuntimeError(msg)

            f06_file.write(page_stamp % page_num)
            page_num += 1
            break
        return page_num - 1

    def write_op2(self, op2, op2_ascii, itable, new_result, date,
                  is_mag_phase=False, endian='>'):
        """writes an OP2"""
        import inspect
        from struct import Struct, pack
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write('%s.write_op2: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        if itable == -1:
            self._write_table_header(op2, op2_ascii, date)
            itable = -3

        ntotali = self.num_wide

        device_code = self.device_code
        op2_ascii.write('  ntimes = %s\n' % self.ntimes)

        if self.is_sort1:
            struct1 = Struct(endian + b'i 3f')
        else:
            raise NotImplementedError('SORT2')

        for itime in range(self.ntimes):
            eids = self.element[itime, :]
            nelements = len(eids)
            ntotal = ntotali * nelements
            op2_ascii.write('nelements=%i\n' % nelements)

            eids_device = eids * 10 + self.device_code
            self._write_table_3(op2, op2_ascii, new_result, itable, itime)

            # record 4
            itable -= 1
            header = [4, itable, 4,
                      4, 1, 4,
                      4, 0, 4,
                      4, ntotal, 4,
                      4 * ntotal]
            op2.write(pack('%ii' % len(header), *header))
            op2_ascii.write('r4 [4, 0, 4]\n')
            op2_ascii.write('r4 [4, %s, 4]\n' % (itable))
            op2_ascii.write('r4 [4, %i, 4]\n' % (4 * ntotal))

            energy = self.data[itime, :, 0]
            percent = self.data[itime, :, 1]
            density = self.data[itime, :, 2]
            #print(eids_device)
            for (eid, eid_device, energyi, percenti, densityi) in zip(eids, eids_device, energy, percent, density):
                data = [eid_device, energyi, percenti, densityi]
                #print(data)

                #vals = (fxi, fyi, fzi, mxi, myi, mzi)
                #vals2 = write_imag_floats_13e(vals, is_mag_phase)
                #(fxir, fyir, fzir, mxir, myir, mzir,
                 #fxii, fyii, fzii, mxii, myii, mzii) = vals2
                #op2_ascii.write('0%26i   %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                               #' %26s   %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                                   #eid, fxir, fyir, fzir, mxir, myir, mzir,
                                   #'', fxii, fyii, fzii, mxii, myii, mzii))
                op2.write(struct1.pack(*data))

            itable -= 1
            header = [4 * ntotal,]
            op2.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable

    def _write_table_3(self, op2, op2_ascii, new_result, itable, itime): #itable=-3, itime=0):
        import inspect
        from struct import pack
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write('%s.write_table_3: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        #print('new_result=%s itable=%s' % (new_result, itable))
        if new_result and itable != -3:
            header = [
                4, 146, 4,
            ]
        else:
            header = [
                4, itable, 4,
                4, 1, 4,
                4, 0, 4,
                4, 146, 4,
            ]
        op2.write(pack(b'%ii' % len(header), *header))
        op2_ascii.write('table_3_header = %s\n' % header)

        approach_code = self.approach_code
        table_code = self.table_code
        isubcase = self.isubcase
        element_name = ('%-8s' % self.element_name).encode('ascii')
        #element_type = self.element_type
        #assert isinstance(self.element_type, int), self.element_type

        #[
            #'aCode', 'tCode', 'element_type', 'isubcase',
            #'???', '???', '???', 'load_set'
            #'format_code', 'num_wide', 's_code', '???',
            #'???', '???', '???', '???',
            #'???', '???', '???', '???',
            #'???', '???', '???', '???',
            #'???', 'Title', 'subtitle', 'label']
        #random_code = self.random_code
        format_code = self.format_code
        s_code = 0 # self.s_code
        num_wide = self.num_wide
        acoustic_flag = 0
        thermal = 0
        title = b'%-128s' % self.title.encode('ascii')
        subtitle = b'%-128s' % self.subtitle.encode('ascii')
        label = b'%-128s' % self.label.encode('ascii')
        ftable3 = b'50i 128s 128s 128s'
        oCode = 0
        load_set = 0
        #print(self.code_information())

        ftable3 = b'i' * 50 + b'128s 128s 128s'
        field6 = 0
        field7 = 0
        if 0:
            pass
        elif self.analysis_code == 1:
            field5 = self.lsdvmns[itime]
        elif self.analysis_code == 2:
            field5 = self.modes[itime]
            #field6 = self.eigns[itime]
            #field7 = self.cycles[itime]
            field5 = int(field5)
            assert isinstance(field5, int), type(field5)
            #assert isinstance(field6, float), type(field6)
            #assert isinstance(field7, float), type(field7)
            #ftable3 = set_table3_field(ftable3, 6, b'f') # field 6
            #ftable3 = set_table3_field(ftable3, 7, b'f') # field 7
        elif self.analysis_code == 5:
            #try:
            #print(self)
            field5 = self.freq2s[itime]
            #except AttributeError:  # pragma: no cover
                #print(self)
                #raise
            ftable3 = set_table3_field(ftable3, 5, b'f') # field 5
        elif self.analysis_code == 6:
            #if hasattr(self, 'times'):
            try:
                field5 = self.times[itime]
            ##elif hasattr(self, 'dts'):
                ##field5 = self.times[itime]
            #else:  # pragma: no cover
            except:
                print(self.get_stats())
                raise NotImplementedError('cant find times or dts on analysis_code=8')
            ftable3 = set_table3_field(ftable3, 5, b'f') # field 5
        #elif self.analysis_code == 7:  # pre-buckling
            #field5 = self.loadIDs[itime] # load set number
        #elif self.analysis_code == 8:  # post-buckling
            #field5 = self.lsdvmns[itime] # load set number
            #if hasattr(self, 'eigns'):
                #field6 = self.eigns[itime]
            #elif hasattr(self, 'eigrs'):
                #field6 = self.eigrs[itime]
            #else:  # pragma: no cover
                #print(self.get_stats())
                #raise NotImplementedError('cant find eigns or eigrs on analysis_code=8')
            #assert isinstance(field6, float_types), type(field6)
            #ftable3 = set_table3_field(ftable3, 6, b'f') # field 6
        elif self.analysis_code == 9:  # complex eigenvalues
            field5 = self.modes[itime]
            #if hasattr(self, 'eigns'):
                #field6 = self.eigns[itime]
            #elif hasattr(self, 'eigrs'):
                #field6 = self.eigrs[itime]
            #else:  # pragma: no cover
                #print(self.get_stats())
                #raise NotImplementedError('cant find eigns or eigrs on analysis_code=8')
            #ftable3 = set_table3_field(ftable3, 6, b'f') # field 6
            #field7 = self.eigis[itime]
            #ftable3 = set_table3_field(ftable3, 7, b'f') # field 7
            assert isinstance(field5, int), type(field5)
        elif self.analysis_code == 10:  # nonlinear statics
            field5 = self.loadFactors[itime]
            ftable3 = set_table3_field(ftable3, 5, b'f') # field 5; load step
        #elif self.analysis_code == 11:  # old geometric nonlinear statics
            #field5 = self.loadIDs[itime] # load set number
        else:
            raise NotImplementedError(self.analysis_code)
        # we put these out of order, so we can set the element_name field
        # (that spans 2 fields) in the right order
        ftable3 = set_table3_field(ftable3, 7, b'4s') # field 7
        ftable3 = set_table3_field(ftable3, 6, b'4s') # field 7
        element_name0 = element_name[:4]
        element_name1 = element_name[4:]
        print(str(self))
        table3 = [
            approach_code, table_code, 0, isubcase, field5,
            element_name0, element_name1, load_set, format_code, num_wide,
            s_code, acoustic_flag, 0, 0, 0,
            0, 0, 0, 0, 0,
            0, 0, thermal, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0,
            title, subtitle, label,
        ]
        assert table3[22] == thermal

        n = 0
        for v in table3:
            if isinstance(v, (int, float, np.float32)):
                n += 4
            elif isinstance(v, string_types):
                #print('%i %r' % (len(v), v))
                n += len(v)
            else:
                #print('write_table_3', v)
                n += len(v)
        assert n == 584, n
        data = [584] + table3 + [584]
        fmt = b'i' + ftable3 + b'i'
        #f.write(pack(fascii, '%s header 3c' % self.table_name, fmt, data))
        op2_ascii.write('%s header 3c = %s\n' % (self.table_name, data))

        op2.write(pack(fmt, *data))


class ComplexStrainEnergyArray(BaseElement):
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
        BaseElement.__init__(self, data_code, isubcase)
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

    @property
    def is_real(self):
        return False

    @property
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
        """sizes the vectorized attributes of the ComplexStrainEnergyArray"""
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
        self._times = np.zeros(self.ntimes, dtype=dtype)
        #self.element = np.zeros(self.nelements, dtype='int32')
        self.element = np.zeros((self.ntimes, self.nelements), dtype='int32')
        #self.element_data_type = empty(self.nelements, dtype='|U8')

        #[energy, percent, density]
        assert isinstance(self.ntimes, integer_types), self.ntimes
        assert isinstance(self.ntotal, integer_types), self.ntotal
        self.data = np.zeros((self.ntimes, self.nelements, 4), dtype='float32')

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
        #import pandas as pd
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
            #if self.nonlinear_factor not in (None, np.nan):
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
        assert self.is_sort1 == table.is_sort1

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
                    #print(t1, t2)
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
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
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
        msg.append('  element: [%s, nelements]; eid=100000000 -> total\n' % (ntimes_word))
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        #msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
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
            f06_file.write(''.join(header + msg_temp2))


            fmt1 = ' ' * 23 + '%10i      %-13s /  %-13s               %7.4f             %s\n'
            fmt2 = '\n                        TYPE = %-8s SUBTOTAL       %13s                 %.4f\n'

            for (eid, energyri, energyii, percenti, densityi) in zip(eids, energyr, energyi, percent, density):
                senergyr = write_float_13e(energyri)
                senergyi = write_float_13e(energyii)
                sdensityi = write_float_13e(densityi)
                # ELEMENT-ID    STRAIN-ENERGY   PERCENT OF TOTAL  STRAIN-ENERGY-DENSITY
                #          1   -8.307121E-12         0.0052           -2.886861E-12
                #if eid == 100000000:
                    #f06_file.write(fmt2 % (self.element_name, senergyi, percenti))
                    #break
                f06_file.write(fmt1 % (
                    eid, senergyr, senergyi, percenti, sdensityi))
            f06_file.write(page_stamp % page_num)
            page_num += 1
            #break
        return page_num - 1
