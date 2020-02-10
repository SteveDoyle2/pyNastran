from math import sqrt
from struct import pack
from typing import List

import numpy as np

from pyNastran.op2.result_objects.op2_objects import BaseScalarObject
from pyNastran.op2.op2_interface.write_utils import write_table_header # set_table3_field,
from pyNastran.f06.f06_formatting import write_floats_13e


class RealEigenvalues(BaseScalarObject):
    """
    cycle = sqrt(abs(eigenvalue)) / (2. * pi)
    radians = sqrt(abs(eigenvalue))
    """
    def __init__(self, title, table_name: str, nmodes=0):
        #self.modeNumber = []
        BaseScalarObject.__init__(self)
        self.title = title
        self.table_name = table_name
        self.mode = np.zeros(nmodes, dtype='int32')
        self.extraction_order = np.zeros(nmodes, dtype='int32')
        self.eigenvalues = np.zeros(nmodes, dtype='float32')
        self.radians = np.zeros(nmodes, dtype='float32')
        self.cycles = np.zeros(nmodes, dtype='float32')
        self.generalized_mass = np.zeros(nmodes, dtype='float32')
        self.generalized_stiffness = np.zeros(nmodes, dtype='float32')
        self.data_frame = None

    def __eq__(self, table):  # pragma: no cover
        return True

    def get_stats(self, short=False) -> List[str]:
        msg = []
        neigenvalues = len(self.extraction_order)
        msg.append('  type=%s neigenvalues=%s\n' % (self.__class__.__name__,
                                                    neigenvalues))
        msg.append('  title, extraction_order, eigenvalues, radians, '
                   'cycles, generalized_mass, generalized_stiffness\n')
        return msg

    @property
    def is_real(self):
        return True

    @property
    def is_complex(self):
        return False

    def add_f06_line(self, data, imode):
        (mode_num, extract_order, eigenvalue, radian, cycle, gen_mass, gen_stiffness) = data
        self.mode[imode] = mode_num
        self.extraction_order[imode] = extract_order
        self.eigenvalues[imode] = eigenvalue
        self.radians[imode] = radian
        #cyclei = sqrt(abs(eigenvalue)) / (2. * pi)
        #if not allclose(cycle, cyclei):
            #print('cycle=%s cyclei=%s' % (cycle, cyclei))
        self.cycles[imode] = cycle
        self.generalized_mass[imode] = gen_mass
        self.generalized_stiffness[imode] = gen_stiffness

    def get_headers(self) -> List[str]:
        headers = ['eigenvalue', 'radians', 'cycle', 'generalized_mass', 'generalized_stiffness']
        return headers

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()
        #cycle = sqrt(abs(eigenvalue)) / (2. * pi)
        data = np.vstack([self.eigenvalues, self.radians, self.cycles,
                          self.generalized_mass, self.generalized_stiffness]).T
        modes_extraction_order = np.vstack([self.mode, self.extraction_order]).T

        df1 = pd.DataFrame(modes_extraction_order)
        df1.columns = ['Mode', 'ExtractionOrder']
        df2 = pd.DataFrame(data)
        df2.columns = headers
        self.data_frame = df1.join(df2)

    def add_f06_data(self, data):
        for i, line in enumerate(data):
            self.add_f06_line(line, i)

    def write_f06(self, f06_file, header, page_stamp, page_num=1):
        title = ''
        if self.title is not None:
            title = '%s' % str(self.title).center(124).rstrip() + '\n'
        msg = header + ['                                              R E A L   E I G E N V A L U E S\n', title,
                        '   MODE    EXTRACTION      EIGENVALUE            RADIANS             CYCLES            GENERALIZED         GENERALIZED\n',
                        '    NO.       ORDER                                                                       MASS              STIFFNESS\n']
        for (imode, mode_num) in enumerate(self.mode):
            order = self.extraction_order[imode]
            eigenvalue = self.eigenvalues[imode]
            #cycle = sqrt(abs(eigenvalue)) / (2. * pi)

            omega = self.radians[imode]
            freq = self.cycles[imode]
            mass = self.generalized_mass[imode]
            stiff = self.generalized_stiffness[imode]
            [eigen, omega, freq, mass, stiff] = write_floats_13e(
                [eigenvalue, omega, freq, mass, stiff])
            msg.append(' %8s  %8s       %-13s       %-13s       %-13s       %-13s       %s\n' % (
                mode_num, order, eigen, omega, freq, mass, stiff))
        msg.append(page_stamp % page_num)
        f06_file.write(''.join(msg))
        return page_num

    def write_op2(self, op2, op2_ascii, itable, new_result, date,
                  is_mag_phase=False, endian='>'):
        """writes an OP2"""
        import inspect
        from struct import Struct, pack
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write('%s.write_op2: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        if itable == -1:
            _write_table_header(self.table_name, op2, op2_ascii, date)
            itable = -3

        #if isinstance(self.nonlinear_factor, float):
            #op2_format = '%sif' % (7 * self.ntimes)
            #raise NotImplementedError()
        #else:
            #op2_format = 'i21f'
        #s = Struct(op2_format)

        # table 4 info
        #ntimes = self.data.shape[0]
        #nnodes = self.data.shape[1]

        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)

        #ntotali = self.num_wide
        #ntotal = ntotali * nelements

        #print('shape = %s' % str(self.data.shape))
        #assert self.ntimes == 1, self.ntimes

        #fmt = '%2i %6f'
        #print('ntotal=%s' % (ntotal))
        #assert ntotal == 193, ntotal
        nmodes = len(self.mode)
        ntotal = nmodes * 7

        structi = Struct(endian + b'ii5f')

        self._write_table_3(op2, op2_ascii, new_result, itable, 0)

        # record 4
        #print('stress itable = %s' % itable)
        itable -= 1
        #print('4, %s' % itable)
        header = [4, itable, 4,
                  4, 1, 4,
                  4, 0, 4,
                  4, ntotal, 4,
                  4 * ntotal]
        op2.write(pack('%ii' % len(header), *header))
        op2_ascii.write('r4 [4, 0, 4]\n')
        op2_ascii.write('r4 [4, %s, 4]\n' % (itable))
        op2_ascii.write('r4 [4, %i, 4]\n' % (4 * ntotal))

        for (imode, mode_num) in enumerate(self.mode):
            extract_order = self.extraction_order[imode]
            eigenvalue = self.eigenvalues[imode]
            #cycle = sqrt(abs(eigenvalue)) / (2. * pi)

            omega = self.radians[imode]
            freq = self.cycles[imode]
            gen_mass = self.generalized_mass[imode]
            gen_stiffness = self.generalized_stiffness[imode]
            #(mode_num, extract_order, eigenvalue, radian, cycle, gen_mass, gen_stiffness) = data
            data = [mode_num, extract_order, eigenvalue, omega, freq, gen_mass, gen_stiffness]

            [eigen, omega, freq, gen_mass, gen_stiffness] = write_floats_13e(
                [eigenvalue, omega, freq, gen_mass, gen_stiffness])
            op2_ascii.write(' %8s  %8s       %-13s       %-13s       %-13s       %-13s       %s\n' % (
                mode_num, extract_order, eigen, omega, freq, gen_mass, gen_stiffness))

            op2.write(structi.pack(*data))

        itable -= 1
        header = [4 * ntotal,]
        op2.write(pack('i', *header))
        op2_ascii.write('footer = %s\n' % header)
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

        #approach_code = self.approach_code
        approach_code = 0
        #table_code = self.table_code
        table_code = 0
        #isubcase = self.isubcase
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
        #s_code = 0 # self.s_code
        #num_wide = self.num_wide
        num_wide = 7
        #acoustic_flag = 0
        #thermal = 0
        title = b'%-128s' % self.title.encode('ascii')
        subtitle = b' '*128
        label = b' '*128
        assert len(title) == 128
        assert len(subtitle) == 128
        assert len(label) == 128
        #subtitle = b'%-128s' % self.subtitle.encode('ascii')
        #label = b'%-128s' % self.label.encode('ascii')
        ftable3 = b'50i 128s 128s 128s'
        #oCode = 0
        #load_set = 0
        #print(self.code_information())

        #print(title, len(title))
        ftable3 = b'i' * 50 + b'128s 128s 128s'
        #field6 = 0
        #field7 = 0

        #self.seven = self.add_data_parameter(data, 'seven', b'i', 10, False)  # seven
        #self.residual_flag = self.add_data_parameter(data, 'residual_flag', b'i', 11, False)
        #self.fluid_flag = self.add_data_parameter(data, 'fluid_flag', b'i', 12, False)

        #seven = 1
        fluid_flag = 1
        residual_flag = 1
        table3 = [
            approach_code, table_code, 0, 0, 0,
            0, 0, 0, 0, num_wide,
            residual_flag, fluid_flag, 0, 0, 0,
            0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0,
            title, subtitle, label,
        ]
        assert len(table3) == 53, len(table3)
        assert table3[12-1] == fluid_flag, fluid_flag
        assert table3[11-1] == residual_flag, residual_flag

        n = 0
        for i, v in enumerate(table3):
            if isinstance(v, (int, float)):
                n += 4
            elif isinstance(v, str):
                #print(len(v), v)
                n += len(v)
            else:
                print('write_table_3', i, v)
                n += len(v)
        assert n == 584, n
        data = [584] + table3 + [584]
        fmt = b'i' + ftable3 + b'i'
        #print(fmt)
        #print(data)
        #f.write(pack(fascii, '%s header 3c' % self.table_name, fmt, data))
        op2_ascii.write('%s header 3c = %s\n' % (self.table_name, data))
        op2.write(pack(fmt, *data))

    def __repr__(self):
        if self.data_frame is not None:
            return str(self.data_frame)

        msg = '%-7s %15s %15s %10s %10s %10s %15s\n' % (
            'ModeNum', 'ExtractionOrder', 'Eigenvalue', 'Radians', 'Cycles', 'GenMass', 'GenStiffness')
        for imode, mode_num in enumerate(self.mode):
            extract_order = self.extraction_order[imode]
            eigenvalue = self.eigenvalues[imode]
            radian = self.radians[imode]

            cycle = sqrt(abs(eigenvalue)) / (2. * np.pi)
            #cycle = self.cycles[imode]
            gen_m = self.generalized_mass[imode]
            gen_k = self.generalized_stiffness[imode]
            msg += '%-7s %15s %15s %10s %10s %10s %15s\n' % (
                mode_num, extract_order, eigenvalue, radian, cycle, gen_m, gen_k)
        return msg


def _write_table_header(table_name: str, op2_file, fascii, date):
    table_name = '%-8s' % table_name # 'BOUGV1  '
    fascii.write('%s._write_table_header\n' % table_name)
    #get_nmarkers- [4, 0, 4]
    #marker = [4, 2, 4]
    #table_header = [8, 'BOUGV1  ', 8]
    write_table_header(op2_file, fascii, table_name)


    #read_markers -> [4, -1, 4]
    #get_nmarkers- [4, 0, 4]
    #read_record - marker = [4, 7, 4]
    #read_record - record = [28, recordi, 28]

    #write_markers(op2_file, fascii, '  %s header1a' % self.table_name, [-1, 0, 7])
    data_a = [4, -1, 4,]
    #data_a = []
    #data_b = [4, -1, 4,]
    data_c = [4, 7, 4,]
    data = data_a + data_c
    blank = ' ' * len(table_name)
    fascii.write('%s header1a_i = %s\n' % (table_name, data_a))
    #fascii.write('%s            = %s\n' % (blank, data_b))
    fascii.write('%s            = %s\n' % (blank, data_c))
    op2_file.write(pack('<6i', *data))

    table1_fmt = b'<9i'
    table1 = [
        28,
        102, 0, 0, 0, 512, 0, 0,
        28,
    ]
    fascii.write('%s header1b = %s\n' % (table_name, table1))
    op2_file.write(pack(table1_fmt, *table1))

    #recordi = [subtable_name, month, day, year, 0, 1]

    data = [
        4, -2, 4,
        4, 1, 4,
        4, 0, 4,
        4, 7, 4,
    ]
    fascii.write('%s header2a = %s\n' % (table_name, data))
    op2_file.write(pack(b'<12i', *data))

    month, day, year = date

    subtable_name = b'OUG1    '
    table2 = [
        28,  # 4i -> 13i
        # subtable,todays date 3/6/2014, 0, 1  ( year=year-2000)
        b'%-8s' % subtable_name, month, day, year - 2000, 0, 1,
        28,
        ]
    table2_format = 'i8s6i'
    fascii.write('%s header2b = %s\n' % (table_name, table2))
    op2_file.write(pack(table2_format, *table2))


class ComplexEigenvalues(BaseScalarObject):
    """
    cycle = freq = eigi / (2*pi)
    radians = eigi
    damping = atan2(eigi, eigr) * 2
    """
    def __init__(self, title, table_name, nmodes):
        BaseScalarObject.__init__(self)
        self.title = title
        self.table_name = table_name

        self.mode = np.zeros(nmodes, dtype='int32')
        self.extraction_order = np.zeros(nmodes, dtype='int32')
        self.eigenvalues = np.zeros(nmodes, dtype='complex64')
        self.cycles = np.zeros(nmodes, dtype='float32')
        self.damping = np.zeros(nmodes, dtype='float32')

        self.data_frame = None

    def __eq__(self, table):  # pragma: no cover
        return True

    def get_stats(self, short=False) -> List[str]:
        neigenvalues = len(self.extraction_order)
        msg = []
        msg.append('  type=%s neigenvalues=%s\n' % (self.__class__.__name__, neigenvalues))
        msg.append('  isubcase, extraction_order, eigenvalues, '
                   'cycles, damping\n')
        return msg

    @property
    def is_real(self):
        return False

    @property
    def is_complex(self):
        return True

    def add_op2_line(self, data, i):
        (root_num, extract_order, eigr, eigi, cycle, damping) = data
        self.mode[i] = root_num
        self.extraction_order[i] = extract_order
        self.eigenvalues[i] = complex(eigr, eigi)
        self.cycles[i] = cycle
        self.damping[i] = damping

    def add_op2_data(self, data):
        for imode, line in enumerate(data):
            self.add_op2_line(line, imode)

    def get_headers(self) -> List[str]:
        headers = ['eigenvalue', 'frequency', 'damping']
        return headers

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()

        cdata = self.eigenvalues
        fdata = np.vstack([self.cycles, self.damping]).T
        modes_extraction_order = np.vstack([self.mode, self.extraction_order]).T

        df1 = pd.DataFrame(modes_extraction_order)
        df1.columns = ['Mode', 'ExtractionOrder']
        df2 = pd.DataFrame(cdata)
        df2.columns = [headers[0]]
        df3 = pd.DataFrame(fdata)
        df3.columns = headers[1:]
        self.data_frame = df1.join([df2, df3])
        #print(self.data_frame)

    def write_f06(self, f06_file, header, page_stamp, page_num=1):  # not proper msg start
        title = ''
        if self.title is not None:
            title = '%s' % str(self.title).center(124).rstrip() + '\n'
        msg = header + ['                                        C O M P L E X   E I G E N V A L U E   S U M M A R Y\n', title,
                        '0                ROOT     EXTRACTION                  EIGENVALUE                     FREQUENCY              DAMPING\n',
                        '                  NO.        ORDER             (REAL)           (IMAG)                (CYCLES)            COEFFICIENT\n']

        for (imode, mode) in enumerate(self.mode):
            extract_order = self.extraction_order[imode]
            eigr = self.eigenvalues[imode].real
            eigi = self.eigenvalues[imode].imag

            freq = self.cycles[imode]
            damping = self.damping[imode]
            [eigr, eigi, freq, damping] = write_floats_13e([eigr, eigi, freq, damping])
            #            imode order      eigr     eigi          freq        damping
            msg.append(' %22s  %10s         %-15s  %-13s         %-13s         %s\n' % (
                mode, extract_order, eigr, eigi, freq, damping))

        msg.append(page_stamp % page_num)
        f06_file.write(''.join(msg))
        return page_num

    def write_op2(self, op2, op2_ascii, itable, new_result, date,
                  is_mag_phase=False, endian='>'):
        """writes an OP2"""
        import inspect
        from struct import Struct, pack
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write('%s.write_op2: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        if itable == -1:
            _write_table_header(self.table_name, op2, op2_ascii, date)
            itable = -3

        #if isinstance(self.nonlinear_factor, float):
            #op2_format = '%sif' % (7 * self.ntimes)
            #raise NotImplementedError()
        #else:
            #op2_format = 'i21f'
        #s = Struct(op2_format)

        # table 4 info
        #ntimes = self.data.shape[0]
        #nnodes = self.data.shape[1]

        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)

        #ntotali = self.num_wide
        #ntotal = ntotali * nelements

        #print('shape = %s' % str(self.data.shape))
        #assert self.ntimes == 1, self.ntimes

        #fmt = '%2i %6f'
        #print('ntotal=%s' % (ntotal))
        #assert ntotal == 193, ntotal
        nmodes = len(self.mode)
        ntotal = nmodes * 7

        structi = Struct(endian + b'ii4f')

        self._write_table_3(op2, op2_ascii, new_result, itable, 0)

        # record 4
        #print('stress itable = %s' % itable)
        itable -= 1
        #print('4, %s' % itable)
        header = [4, itable, 4,
                  4, 1, 4,
                  4, 0, 4,
                  4, ntotal, 4,
                  4 * ntotal]
        op2.write(pack('%ii' % len(header), *header))
        op2_ascii.write('r4 [4, 0, 4]\n')
        op2_ascii.write('r4 [4, %s, 4]\n' % (itable))
        op2_ascii.write('r4 [4, %i, 4]\n' % (4 * ntotal))

        for (imode, mode_num) in enumerate(self.mode):
            extract_order = self.extraction_order[imode]
            eigr = self.eigenvalues[imode].real
            eigi = self.eigenvalues[imode].imag

            freq = self.cycles[imode]
            damping = self.damping[imode]
            data = [mode_num, extract_order, eigr, eigi, freq, damping]


            [eigr, eigi, freq, damping] = write_floats_13e([eigr, eigi, freq, damping])
            #            imode order      eigr     eigi          freq        damping
            op2_ascii.write(' %22s  %10s         %-15s  %-13s         %-13s         %s\n' % (
                mode_num, extract_order, eigr, eigi, freq, damping))

            op2.write(structi.pack(*data))

        itable -= 1
        header = [4 * ntotal,]
        op2.write(pack('i', *header))
        op2_ascii.write('footer = %s\n' % header)
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

        #approach_code = self.approach_code
        approach_code = 0
        #table_code = self.table_code
        table_code = 0
        #isubcase = self.isubcase
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
        #s_code = 0 # self.s_code
        #num_wide = self.num_wide
        num_wide = 7
        #acoustic_flag = 0
        #thermal = 0
        title = b'%-128s' % self.title.encode('ascii')
        subtitle = b' '*128
        label = b' '*128
        assert len(title) == 128
        assert len(subtitle) == 128
        assert len(label) == 128
        #subtitle = b'%-128s' % self.subtitle.encode('ascii')
        #label = b'%-128s' % self.label.encode('ascii')
        ftable3 = b'50i 128s 128s 128s'
        #oCode = 0
        #load_set = 0
        #print(self.code_information())

        #print(title, len(title))
        ftable3 = b'i' * 50 + b'128s 128s 128s'
        #field6 = 0
        #field7 = 0

        #self.seven = self.add_data_parameter(data, 'seven', b'i', 10, False)  # seven
        #self.residual_flag = self.add_data_parameter(data, 'residual_flag', b'i', 11, False)
        #self.fluid_flag = self.add_data_parameter(data, 'fluid_flag', b'i', 12, False)

        #seven = 1
        fluid_flag = 1
        residual_flag = 1
        table3 = [
            approach_code, table_code, 0, 0, 0,
            0, 0, 0, 0, num_wide,
            residual_flag, fluid_flag, 0, 0, 0,
            0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0,
            title, subtitle, label,
        ]
        assert len(table3) == 53, len(table3)
        assert table3[12-1] == fluid_flag, fluid_flag
        assert table3[11-1] == residual_flag, residual_flag

        n = 0
        for i, v in enumerate(table3):
            if isinstance(v, (int, float)):
                n += 4
            elif isinstance(v, str):
                #print(len(v), v)
                n += len(v)
            else:
                print('write_table_3', i, v)
                n += len(v)
        assert n == 584, n
        data = [584] + table3 + [584]
        fmt = b'i' + ftable3 + b'i'
        #print(fmt)
        #print(data)
        #f.write(pack(fascii, '%s header 3c' % self.table_name, fmt, data))
        op2_ascii.write('%s header 3c = %s\n' % (self.table_name, data))
        op2.write(pack(fmt, *data))

    def __repr__(self):
        msg = '%-7s %15s %15s %10s %10s %10s\n' % (
            'RootNum', 'ExtractionOrder', 'Eigenvalue', '', 'Cycles', 'Damping')
        msg += '%-7s %15s %15s %10s\n' % ('', '', 'Real', 'Imaginary')
        for imode, unused_mode in enumerate(self.mode):
            extract_order = self.extraction_order[imode]
            eigenvalue = self.eigenvalues[imode]
            cycle = self.cycles[imode]
            damping = self.damping[imode]
            msg += '%-7s %15s %15s %10s %10s %10s\n' % (
                imode, extract_order,
                eigenvalue.real, eigenvalue.imag, cycle, damping)
        return msg


class BucklingEigenvalues(BaseScalarObject):
    def __init__(self, title, table_name, nmodes=0):
        BaseScalarObject.__init__(self)
        self.title = title
        self.table_name = table_name
        self.mode = np.zeros(nmodes, dtype='int32')
        self.extraction_order = np.zeros(nmodes, dtype='int32')
        self.eigenvalues = np.zeros(nmodes, dtype='float32')
        self.freqs = np.zeros(nmodes, dtype='float32')
        self.omegas = np.zeros(nmodes, dtype='float32')
        self.generalized_mass = np.zeros(nmodes, dtype='float32')
        self.generalized_stiffness = np.zeros(nmodes, dtype='float32')
        self.data_frame = None

    def __eq__(self, table):  # pragma: no cover
        return True

    def get_stats(self, short=False) -> List[str]:
        neigenvalues = len(self.extraction_order)
        msg = []
        msg.append('  type=%s neigenvalues=%s\n' % (self.__class__.__name__, neigenvalues))
        msg.append('  imode, extraction_order, eigenvalues, '
                   'radians, cycles, generalized_mass, generalized_stiffness\n')
        return msg

    @property
    def is_real(self):
        return False

    @property
    def is_complex(self):
        return False

    def is_buckling(self):
        return True

    def add_op2_line(self, data, imode):
        (root_num, extract_order, eigr, omega, freq, mass, stiff) = data
        self.mode[imode] = root_num
        self.extraction_order[imode] = extract_order
        self.eigenvalues[imode] = eigr
        self.freqs[imode] = freq
        self.omegas[imode] = omega
        self.generalized_mass[imode] = mass
        self.generalized_stiffness[imode] = stiff

    #def add_op2_data(self, data):
        #for i, line in enumerate(data):
            #self.add_op2_line(line, i)

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()
        nmodes = len(self.eigenvalues)

        modes_extraction_order = np.zeros((nmodes, 2), dtype='float32')
        fdata = np.zeros((nmodes, 5), dtype='float32')

        imodei = 0
        for (imode, unused_mode) in enumerate(self.mode):
            eigi = self.eigenvalues[imode]
            extraction_order = self.extraction_order[imode]
            freq = self.freqs[imode]
            omega = self.omegas[imode]
            gen_m = self.generalized_mass[imode]
            gen_k = self.generalized_stiffness[imode]

            fdata[imodei, :] = [eigi, freq, omega, gen_m, gen_k]
            modes_extraction_order[imodei, :] = [imode, extraction_order]
            imodei += 1
        df1 = pd.DataFrame(modes_extraction_order)
        df1.columns = ['Mode', 'ExtractionOrder']
        df2 = pd.DataFrame(fdata)
        df2.columns = headers
        self.data_frame = df1.join([df2])
        #print(self.data_frame)

    def get_headers(self) -> List[str]:
        headers = ['eigenvalue', 'radians', 'cycles', 'generalized_mass', 'generalized_stiffness']
        return headers

    def write_f06(self, f06_file, header, page_stamp, page_num=1):  # not proper msg start
        title = ''
        if self.title is not None:
            title = '%s' % str(self.title).center(124).rstrip() + '\n'
        msg = header + ['                                              R E A L   E I G E N V A L U E S\n', title,
                        '   MODE    EXTRACTION      EIGENVALUE            RADIANS             CYCLES            GENERALIZED         GENERALIZED\n',
                        '    NO.       ORDER                                                                       MASS              STIFFNESS\n']
        f06_file.write(''.join(msg))

        for (imode, unused_mode) in enumerate(self.mode):
            order = self.extraction_order[imode]
            eigr = self.eigenvalues[imode]
            freq = self.freqs[imode]
            omega = self.omegas[imode]
            mass = self.generalized_mass[imode]
            stiff = self.generalized_stiffness[imode]
            [eigr, freq, omega, mass, stiff] = write_floats_13e([eigr, freq, omega, mass, stiff])
            #            i  ord eig ome f   m          k
            f06_file.write(' %8s%10s%20s%20s%20s%20s       %s\n' % (
                imode, order, eigr, omega, freq, mass, stiff))
        f06_file.write(page_stamp % page_num)
        return page_num
