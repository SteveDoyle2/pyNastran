from typing import TextIO
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.result_objects.op2_objects import BaseElement, get_times_dtype
from pyNastran.op2.result_objects.utils_pandas import build_dataframe_transient_header
from pyNastran.f06.f06_formatting import (
    write_float_13e, write_floats_13e_long, _eigenvalue_header)
from pyNastran.op2.op2_interface.write_utils import set_table3_field
from pyNastran.op2.writer.utils import fix_table3_types
from pyNastran.op2.result_objects.op2_objects import set_as_sort1


TABLE_NAME_TO_TABLE_CODE = {
    'ONRGY1': 18,
}
element_name_to_element_type = {
    'CELAS1': 11,
    'CELAS2': 12,
    'CELAS3': 13,
    'CELAS4': 14,
}


class RealEnergyArray(BaseElement):
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
    def __init__(self, data_code, is_sort1: bool, isubcase: int, dt):
        self.element_type = None
        self.element_name = None
        BaseElement.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific
        self.itime = None
        self.itotal2 = 0
        #self.element_name_count = {}
        self.dt_temp = None

        #if is_sort1:
            #pass
        #else:
            #raise NotImplementedError('SORT2')

    @property
    def is_real(self) -> bool:
        return True

    @property
    def is_complex(self) -> bool:
        return False

    def _reset_indices(self) -> None:
        self.itotal = 0
        self.ielement = 0

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

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size, self.analysis_fmt)
        self.build_data(dtype, idtype, fdtype)

    def build_data(self, dtype, idtype, fdtype):
        """actually performs the build step"""
        self._times = np.zeros(self.ntimes, dtype=self.analysis_fmt)
        #self.element = zeros(self.nelements, dtype='int32')
        #if dtype in 'DMIG':
        #print(self.element_name, self.element_type)
        if self.element_name == 'DMIG':
            self.element = np.zeros((self.ntimes, self.nelements), dtype='|U8')
        else:
            self.element = np.zeros((self.ntimes, self.nelements), dtype=idtype)
        #self.element_data_type = empty(self.nelements, dtype='|U8')

        #[energy, percent, density]
        assert isinstance(self.ntimes, integer_types), self.ntimes
        assert isinstance(self.ntotal, integer_types), self.ntotal
        self.data = np.zeros((self.ntimes, self.nelements, 3), dtype=fdtype)

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
            column_names, column_values = build_dataframe_transient_header(self)
            # Static     strain_energy    percent  strain_energy_density
            # ElementID
            # 6               0.375997   0.878471               1.503990
            # 7               0.462838   1.081362               1.851350
            # 16              0.590421   1.379446               0.590421
            # 17              1.399199   3.269053               0.932799
            # 23              8.086797  18.893791               8.086797
            # 100000000      10.915252  25.502123                    NaN
            self.data_frame = pd.DataFrame(self.data[0], columns=headers, index=element)
            self.data_frame.index.name = 'ElementID'
            self.data_frame.columns.names = ['Static']
        else:
            # we can get into this in a linear case
            # F:\work\pyNastran\examples\Dropbox\move_tpl\setp04.op2

            #nvalues = ntimes * nelements

            #if self.nonlinear_factor not in (None, np.nan):
            column_names, column_values = build_dataframe_transient_header(self)
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
                print(f'self.cannot apply column_names={column_names} '
                      f'to {self.class_name}: {self.element_name}')

            # remove empty rows
            assert self.data_frame is not None
            self.data_frame = self.data_frame[self.data_frame.ElementID != compare]

    @classmethod
    def add_static_case(cls, table_name, element_name, element, data, isubcase,
                        is_sort1=True, is_random=False, is_msc=True,
                        random_code=0, title='', subtitle='', label=''):
        assert len(element.shape) == 1, element.shape
        assert data.shape[0] == 1, data.shape
        assert data.shape[2] == 1, data.shape

        analysis_code = 1  # static
        data_code = oee_data_code(table_name, analysis_code,
                                  is_sort1=is_sort1, is_random=is_random,
                                  random_code=random_code,
                                  title=title, subtitle=subtitle, label=label,
                                  is_msc=is_msc)
        #data_code['loadIDs'] = [0] # TODO: ???
        nmodes = data.shape[0]
        etot = data.sum(axis=0)
        assert data.ndim == 3, data.shape
        assert len(etot) == nmodes, etot
        assert nmodes == 1, data.shape
        data_code['etotpos'] = etot  # TODO: is this the right axis?
        data_code['etotneg'] = etot * 0.
        data_code['lsdvmns'] = [0]  # TODO: ???
        data_code['data_names'] = []

        # I'm only sure about the 1s in the strains and the
        # corresponding 0s in the stresses.
        #if is_stress:
            #data_code['stress_bits'] = [0, 0, 0, 0]
            #data_code['s_code'] = 0
        #else:
            #data_code['stress_bits'] = [0, 1, 0, 1]
            #data_code['s_code'] = 1  # strain?
        element_type = element_name_to_element_type[element_name]
        data_code['element_name'] = element_name
        data_code['element_type'] = element_type
        #data_code['load_set'] = 1

        ntimes = data.shape[0]
        nnodes = data.shape[1]
        dt = None
        obj = cls(data_code, is_sort1, isubcase, dt)
        nelements = len(element)
        obj.element = element.reshape(1, nelements)

        # [energy, percent, density]
        ntimes = 1
        data += 1
        obj.data = np.full((ntimes, nelements+1, 3), np.nan, dtype=data.dtype)
        totals = data.sum(axis=0)
        #print(obj.data[0, :, :])
        #print('totals =', totals)
        obj.data[:, :-1, 0] = data
        percent = data / totals
        obj.data[:, :-1, 1] = percent
        #obj.data[:, :-1, 2] = data / totals  # density???
        obj.data[:, -1, 0] = np.nansum(data, axis=0)
        obj.data[:, -1, 1] = np.nansum(percent, axis=0)
        #print(obj.data[0, :, :])

        ntotals = len(totals)
        assert ntimes == ntotals, ntotals

        obj.ntimes = ntimes
        obj.ntotal = nnodes
        obj._times = [None]
        obj.is_built = True
        return obj

    def __eq__(self, table):  # pragma: no cover
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
                for ielem, eid in enumerate(self.element[itime, :]):
                    t1 = self.data[itime, ielem, :]
                    t2 = table.data[itime, ielem, :]
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
        assert self.sort_method == 1, self
        #itime = self.itime // self.nelement_types
        assert (isinstance(eid, int) and eid > 0) or isinstance(eid, bytes), 'dt=%s eid=%s' % (dt, eid)
        itime = self.itime
        self._times[itime] = dt
        self.element[itime, self.ielement] = eid

        if self.element_name == 'DMIG':
            if not np.isnan(densityi):
                raise RuntimeError(
                    f'{self.class_name}: itime={self.itime} ielement={self.ielement}; '
                    f'dt={dt} eid={eid} energyi={energyi} percenti={percenti} densityi={densityi}')
            self.data[itime, self.ielement, :] = [energyi, percenti, np.nan]
        else:
            try:
                #self.element_data_type[self.ielement] = etype
                self.data[itime, self.ielement, :] = [energyi, percenti, densityi]
            except (ValueError, IndexError):
                print(f'{self.class_name}: itime={self.itime} ielement={self.ielement}; '
                    f'dt={dt} eid={eid} energyi={energyi} percenti={percenti} densityi={densityi}')
                raise
        self.ielement += 1
        self.itotal += 1

    def finalize(self):
        self.set_as_sort1()

    def set_as_sort1(self):
        """changes the table into SORT1"""
        set_as_sort1(self)

    def get_stats(self, short: bool=False) -> list[str]:
        if not self.is_built:
            return [
                f'<{self.__class__.__name__}>; table_name={self.table_name!r}\n',
                f'  ntimes: {self.ntimes:d}\n',
                f'  ntotal: {self.ntotal:d}\n',
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
        msg.append(f'  etotpos={self.etotpos}; etotneg={self.etotneg}; thresh={self.thresh}\n')
        msg.append(f'  element: [{ntimes_word:s}, nelements]; eid=100000000 -> total\n')
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append(f'  data.shape = {self.data.shape}\n')
        #msg.append('  element type: %s\n' % self.element_type)
        #msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_csv(self, csv_file: TextIO,
                  is_exponent_format: bool=False,
                  is_mag_phase: bool=False, is_sort1: bool=True,
                  write_header: bool=True):
        """
        Strain Energy Table
        -------------------
        Flag, SubcaseID, iTime,       Eid, BLANK, BLANK,  StrainEnergy,     Percent, StrainEnergyDensity, BLANK, BLANK, BLANK
        14,           1,     0,        14,     0,     0,  1.861520E-01,  3.7039E+01,        1.861519E+00,     0,     0,     0
        14,           1,     0,        15,     0,     0,  1.285691E-01,  2.5383E+01,        1.285691E+00,     0,     0,     0
        14,           1,     0, 100000000,     0,     0,  3.147211E-01,  6.2422E+01,                 NAN,     0,     0,     0
        14,           1,     1,        14,     0,     0,  1.944612E-03,  3.8224E-01,        1.944612E-02,     0,     0,     0
        14,           1,     1,        15,     0,     0,  8.724089E-03,  1.7818E+00,        8.724089E-02,     0,     0,     0
        14,           1,     1, 100000000,     0,     0,  1.066870E-02,  2.1740E+00,                 NAN,     0,     0,     0
        14,           1,     2,        14,     0,     0,  2.729307E-02,  5.4615E+00,        2.729307E-01,     0,     0,     0
        14,           1,     2,        15,     0,     0,  2.109067E-01,  4.2134E+01,        2.109067E+00,     0,     0,     0

        """
        name = str(self.__class__.__name__)
        if write_header:
            csv_file.write('# %s\n' % name)
            #'                                    ELEMENT-ID                  STRAIN-ENERGY     PERCENT OF TOTAL    STRAIN-ENERGY-DENSITY\n'
            headers = ['Flag', 'SubcaseID', 'iTime', 'Eid', 'BLANK', 'BLANK', 'StrainEnergy',
                       'Percent', 'StrainEnergyDensity', 'BLANK', 'BLANK', 'BLANK']
            csv_file.write('# ' + ','.join(headers) + '\n')

        flag = 14  # strain energy
        isubcase = self.isubcase
        #times = self._times

        # write the csv
        ntimes = self.data.shape[0]
        eids = self.element
        eid_len = len(str(eids.max()))

        #zero = ' 0.000000E+00'

        ntimes = self.data.shape[0]

        #etype = self.element_data_type
        for itime in range(ntimes):
            #dt = self._times[itime]  # TODO: rename this...
            #header = _eigenvalue_header(self, header, itime, ntimes, dt)
            #total_energy = 0.
            #total_set_energy = 0.

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

            #msg_temp2 = [msg_temp % (self.element_name, total_energy, itime + 1, total_set_energy)]
            #csv_file.write(''.join(header + msg_temp2))

            #fmt1 = ' ' * 36 + '%10s         %-13s                 %.4f             %s\n'
            #fmt1_nan = ' ' * 36 + '%10s         %-13s                 %.4f             %s\n'
            fmt = f'{flag}, {isubcase}, {itime}, %{eid_len}d, 0, 0, %s, %s, %s, 0, 0, 0\n'

            # headers = ['Flag', 'SubcaseID', 'iTime', 'Eid', 'BLANK', 'BLANK', 'StrainEnergy',
            #            'Percent', 'StrainEnergyDensity', 'BLANK', 'BLANK', 'BLANK']
            for (eid, energyi, percenti, densityi) in zip(eids, energy, percent, density):
                if is_exponent_format:
                    [senergyi, percenti, sdensityi] = write_floats_13e_long([energyi, percenti, densityi])
                # ELEMENT-ID    STRAIN-ENERGY   PERCENT OF TOTAL  STRAIN-ENERGY-DENSITY
                #          1   -8.307121E-12         0.0052           -2.886861E-12
                #if eid == 100000000:
                    #csv_file.write(fmt2 % (self.element_name, senergyi, percenti))
                    #break
                #try:
                csv_file.write(fmt % (eid, senergyi, percenti, sdensityi))
                #except TypeError:
                    #print('eid = %r; type=%s' % (eid, type(eid)))
                    #print('senergyi = %r; type=%s' % (senergyi, type(senergyi)))
                    #print('percenti = %r; type=%s' % (percenti, type(percenti)))
                    #print('sdensityi = %r; type=%s' % (sdensityi, type(sdensityi)))
                    #assert np.isnan(sdensityi), 'eid=%s sdensityi=%s' % (eid, sdensityi)
                    #csv_file.write(fmt1_nan % (eid, senergyi, percenti, ''))
                    #if 0:
                        #print('senergyi = %r; type=%s' % (senergyi, type(senergyi)))
                        #print('percenti = %r; type=%s' % (percenti, type(percenti)))
                        #print('sdensityi = %r; type=%s' % (sdensityi, type(sdensityi)))
                        #msg = fmt1 % (eid, senergyi, percenti, sdensityi)
                        #raise TypeError(msg)
                    #raise RuntimeError(msg)
        return

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
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

        if self.analysis_code in {1, 5, 6, 10}:  # statics, frequency, transient, nonlinear statics
            subcase_mode = 'SUBCASE'
        elif self.analysis_code in {2, 9}:  # modes, complex modes
            subcase_mode = 'MODE'
        else:  # pragma: no cover
            raise NotImplementedError(self.code_information())

        msg_temp = (
            '                                           E L E M E N T   S T R A I N   E N E R G I E S\n'
            ' \n'
            '                ELEMENT-TYPE = %-4s                * TOTAL ENERGY OF ALL ELEMENTS IN PROBLEM     =  %s\n'
            '                %-8s       %8i            * TOTAL ENERGY OF ALL ELEMENTS IN SET      -1 =  %s\n'
            '0\n'
            '                                    ELEMENT-ID          STRAIN-ENERGY           PERCENT OF TOTAL    STRAIN-ENERGY-DENSITY\n'
        )
        ntimes = self.data.shape[0]

        #etype = self.element_data_type
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            total_energy = write_float_13e(self.etotpos)
            total_set_energy = write_float_13e(self.etotpos)

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

            msg_temp2 = [msg_temp % (self.element_name, total_energy,
                                     subcase_mode, itime + 1, total_set_energy)]
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

    def write_op2(self, op2_file, op2_ascii, itable: int,
                  new_result, date,
                  is_mag_phase: bool=False, endian: bytes=b'>'):
        """writes an OP2"""
        import inspect
        from struct import Struct, pack
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write(f'{self.__class__.__name__}.write_op2: {call_frame[1][3]}\n')

        if itable == -1:
            self._write_table_header(op2_file, op2_ascii, date)
            itable = -3

        ntotali = self.num_wide

        device_code = self.device_code
        op2_ascii.write(f'  ntimes = {self.ntimes}\n')

        if not self.is_sort1:
            raise NotImplementedError('SORT2')

        apply_device_code = True
        eids_dtype_kind = self.element.dtype.kind
        if eids_dtype_kind == 'i':
            struct1 = Struct(endian + b'i 3f')
        elif eids_dtype_kind == 'U':
            apply_device_code = False
            struct1 = Struct(endian + b'4s 3f')
        else:  # pragma: no cover
            raise NotImplementedError(f'{self.class_name}: element.dtype.kind = {eids_dtype_kind!r}')

        for itime in range(self.ntimes):
            eids = self.element[itime, :]
            nelements = len(eids)
            ntotal = ntotali * nelements
            op2_ascii.write(f'nelements={nelements:d}\n')

            if apply_device_code:
                eids_device = eids * 10 + self.device_code
            else:
                eids_device = eids.view('|S4')
            self._write_table_3(op2_file, op2_ascii, new_result, itable, itime)

            # record 4
            itable -= 1
            header = [4, itable, 4,
                      4, 1, 4,
                      4, 0, 4,
                      4, ntotal, 4,
                      4 * ntotal]
            op2_file.write(pack('%ii' % len(header), *header))
            op2_ascii.write('r4 [4, 0, 4]\n')
            op2_ascii.write(f'r4 [4, {itable:d}, 4]\n')
            op2_ascii.write(f'r4 [4, {4 * ntotal:d}, 4]\n')

            energy = self.data[itime, :, 0]
            percent = self.data[itime, :, 1]
            density = self.data[itime, :, 2]
            for (eid, eid_device, energyi, percenti, densityi) in zip(eids, eids_device, energy, percent, density):
                data = [eid_device, energyi, percenti, densityi]

                #vals = (fxi, fyi, fzi, mxi, myi, mzi)
                #vals2 = write_imag_floats_13e(vals, is_mag_phase)
                #(fxir, fyir, fzir, mxir, myir, mzir,
                 #fxii, fyii, fzii, mxii, myii, mzii) = vals2
                #op2_ascii.write('0%26i   %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                               #' %26s   %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                                   #eid, fxir, fyir, fzir, mxir, myir, mzir,
                                   #'', fxii, fyii, fzii, mxii, myii, mzii))
                op2_file.write(struct1.pack(*data))

            itable -= 1
            header = [4 * ntotal,]
            op2_file.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable

    def _write_table_3(self, op2_file, op2_ascii, new_result,
                       itable, itime):
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
        op2_file.write(pack(b'%ii' % len(header), *header))
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
        s_code = 0  # self.s_code
        num_wide = self.num_wide
        acoustic_flag = 0
        thermal = 0
        title = b'%-128s' % self.title.encode('ascii')
        subtitle = b'%-128s' % self.subtitle.encode('ascii')
        label = b'%-128s' % self.label.encode('ascii')
        #ftable3 = b'50i 128s 128s 128s'
        #oCode = 0
        load_set = 0
        #print(self.code_information())

        ftable3 = b'i' * 50 + b'128s 128s 128s'
        #field6 = 0
        #field7 = 0
        if 0:
            pass
        elif self.analysis_code == 1:
            field5 = self.lsdvmns[itime]
        elif self.analysis_code == 2:
            field5 = self.modes[itime]
            #field6 = self.eigns[itime]
            #field7 = self.cycles[itime]
            field5 = int(field5)
            assert isinstance(field5, integer_types), type(field5)
            #assert isinstance(field6, float), type(field6)
            #assert isinstance(field7, float), type(field7)
            #ftable3 = set_table3_field(ftable3, 6, b'f') # field 6
            #ftable3 = set_table3_field(ftable3, 7, b'f') # field 7
        elif self.analysis_code == 5:
            #try:
            #print(self)
            field5 = self.freqs[itime]
            #except AttributeError:  # pragma: no cover
                #print(self)
                #raise
            ftable3 = set_table3_field(ftable3, 5, b'f')  # field 5
        elif self.analysis_code == 6:
            #if hasattr(self, 'times'):
            try:
                field5 = self.times[itime]
            ##elif hasattr(self, 'dts'):
                ##field5 = self.times[itime]
            #else:  # pragma: no cover
            except Exception:
                print(self.get_stats())
                raise NotImplementedError('cant find times or dts on analysis_code=8')
            ftable3 = set_table3_field(ftable3, 5, b'f')  # field 5
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
            assert isinstance(field5, integer_types), type(field5)
        elif self.analysis_code == 10:  # nonlinear statics
            field5 = self.load_factors[itime]
            ftable3 = set_table3_field(ftable3, 5, b'f')  # field 5; load step
        #elif self.analysis_code == 11:  # old geometric nonlinear statics
            #field5 = self.loadIDs[itime] # load set number
        else:
            raise NotImplementedError(self.analysis_code)
        # we put these out of order, so we can set the element_name field
        # (that spans 2 fields) in the right order
        ftable3 = set_table3_field(ftable3, 7, b'4s')  # field 7
        ftable3 = set_table3_field(ftable3, 6, b'4s')  # field 6
        element_name0 = element_name[:4]
        element_name1 = element_name[4:]
        #str(self)
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

        table3 = fix_table3_types(table3, size=4)
        data = [584] + table3 + [584]
        fmt = b'i' + ftable3 + b'i'
        #f.write(pack(fascii, '%s header 3c' % self.table_name, fmt, data))
        op2_ascii.write('%s header 3c = %s\n' % (self.table_name, data))
        op2_file.write(pack(fmt, *data))


class RealStrainEnergyArray(RealEnergyArray):
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
    @property
    def headers(self) -> list[str]:
        headers = [
            'strain_energy', 'percent', 'strain_energy_density',
        ]
        return headers


class RealKineticEnergyArray(RealEnergyArray):
    @property
    def headers(self) -> list[str]:
        headers = [
            'kinetic_energy', 'percent', 'kinetic_energy_density',
        ]
        return headers


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
        #self.element_name_count = {}
        self.dt_temp = None

        if not is_sort1:
            raise NotImplementedError('SORT2')

    @property
    def is_real(self) -> bool:
        return False

    @property
    def is_complex(self) -> bool:
        return True

    def _reset_indices(self) -> None:
        self.itotal = 0
        self.ielement = 0

    @property
    def headers(self) -> list[str]:
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

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size, self.analysis_fmt)
        self.build_data(dtype)

    def build_data(self, dtype):
        """actually performs the build step"""
        self._times = np.zeros(self.ntimes, dtype=self.analysis_fmt)
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
            #column_names, column_values = build_dataframe_transient_header(self)
            #element = self.element.ravel()
            #self.data_frame = pd.Panel(self.data, items=column_values,
                                       #major_axis=element,
                                       #minor_axis=headers).to_frame()
            #self.data_frame.columns.names = column_names
        #else:
            #nvalues = ntimes * nelements
            #element = self.element.ravel()
            #if self.nonlinear_factor not in (None, np.nan):
                #column_names, column_values = build_dataframe_transient_header(self)
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

    def __eq__(self, table):  # pragma: no cover
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

        if not np.array_equal(self.data, table.data):  # pragma: no cover
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ielem, eid in enumerate(self.element[itime, :]):
                    t1 = self.data[itime, ielem, :]
                    t2 = table.data[itime, ielem, :]
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
        assert self.sort_method == 1, self
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

    def get_stats(self, short: bool=False) -> list[str]:
        if not self.is_built:
            return [
                f'<{self.__class__.__name__}>; table_name={self.table_name!r}\n',
                f'  ntimes: {self.ntimes:d}\n',
                f'  ntotal: {self.ntotal:d}\n',
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
        msg.append(f'  element: [{ntimes_word:s}, nelements]; eid=100000000 -> total\n')
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append(f'  data.shape = {self.data.shape}\n')
        #msg.append('  element type: %s\n' % self.element_type)
        #msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f06_file, header=None, page_stamp: str='PAGE %s', page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        if header is None:
            header = []

        if self.class_name == 'ComplexStrainEnergyArray':
            msg_temp = (
                '                               E L E M E N T   S T R A I N   E N E R G I E S   ( A V E R A G E )                 \n'
                ' \n'
                '                ELEMENT-TYPE = %-5s               * TOTAL ENERGY OF ALL ELEMENTS IN PROBLEM     =   %s\n'
                '                SUBCASE               1            * TOTAL ENERGY OF ALL ELEMENTS IN SET      -1 =   %s\n'
                '0\n'
                '                       ELEMENT-ID       STRAIN-ENERGY (MAG/PHASE)               PERCENT OF TOTAL    STRAIN-ENERGY-DENSITY\n'
                #'                                5       2.027844E-10 /   0.0                         1.2581              2.027844E-09'
            )
        else:
            msg_temp = (
                '                               E L E M E N T   K I N E T I C   E N E R G I E S   ( A V E R A G E )                 \n'
                ' \n'
                '                ELEMENT-TYPE = %-5s               * TOTAL ENERGY OF ALL ELEMENTS IN PROBLEM     =   %s\n'
                '                SUBCASE               1            * TOTAL ENERGY OF ALL ELEMENTS IN SET      -1 =   %s\n'
                '0\n'
                '                       ELEMENT-ID       KINETIC-ENERGY (MAG/PHASE)               PERCENT OF TOTAL    KINETIC-ENERGY-DENSITY\n'
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
            #fmt2 = '\n                        TYPE = %-8s SUBTOTAL       %13s                 %.4f\n'

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


def oee_data_code(table_name, analysis_code,
                  is_sort1=True, is_random=False,
                  random_code=0, title='', subtitle='', label='', is_msc=True):
    """helper for result creation writing"""
    sort1_sort_bit = 0 if is_sort1 else 1
    random_sort_bit = 1 if is_random else 0
    sort_method = 1 if is_sort1 else 2
    assert analysis_code != 0, analysis_code
    #if format_code == 1:
        #format_word = "Real"
    #elif format_code == 2:
        #format_word = "Real/Imaginary"
    #elif format_code == 3:
        #format_word = "Magnitude/Phase"
    #DEVICE_CODE_MAP = {
        #1 : "Print",
        #2 : "Plot",
        #3 : "Print and Plot",
        #4 : "Punch",
        #5 : "Print and Punch",
        #6 : "Plot and Punch",
        #7 : "Print, Plot, and Punch",
    #}

    table_code = TABLE_NAME_TO_TABLE_CODE[table_name]
    sort_code = 1  # TODO: what should this be???

    #table_code = tCode % 1000
    #sort_code = tCode // 1000
    tCode = table_code * 1000 + sort_code

    device_code = 2  # Plot
    approach_code = analysis_code * 10 + device_code
    #print(f'approach_code={approach_code} analysis_code={analysis_code} device_code={device_code}')
    data_code = {
        'nonlinear_factor': None,
        'approach_code': approach_code,
        'analysis_code': analysis_code,
        'sort_bits': [0, sort1_sort_bit, random_sort_bit],  # real, sort1, random
        'sort_method': sort_method,
        'is_msc': is_msc,
        'format_code': 1,  # real
        'table_code': table_code,
        'tCode': tCode,
        'table_name': table_name,  ## TODO: should this be a string?
        'device_code': device_code,
        'random_code': random_code,
        'thermal': 0,
        'title': title,
        'subtitle': subtitle,
        'label': label,
        'num_wide': 8,  # displacement-style table
    }
    return data_code
