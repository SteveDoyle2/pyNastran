from __future__ import print_function
from math import sqrt

import numpy as np

from pyNastran.op2.result_objects.op2_objects import BaseScalarObject
from pyNastran.f06.f06_formatting import write_floats_13e
try:
    import pandas as pd  # type: ignore
except ImportError:
    pass


class RealEigenvalues(BaseScalarObject):
    """
    cycle = sqrt(abs(eigenvalue)) / (2. * pi)
    radians = sqrt(abs(eigenvalue))
    """
    def __init__(self, title, nmodes=0):
        #self.modeNumber = []
        BaseScalarObject.__init__(self)
        self.title = title
        self.mode = np.zeros(nmodes, dtype='int32')
        self.extraction_order = np.zeros(nmodes, dtype='int32')
        self.eigenvalues = np.zeros(nmodes, dtype='float32')
        self.radians = np.zeros(nmodes, dtype='float32')
        self.cycles = np.zeros(nmodes, dtype='float32')
        self.generalized_mass = np.zeros(nmodes, dtype='float32')
        self.generalized_stiffness = np.zeros(nmodes, dtype='float32')
        self.data_frame = None

    def __eq__(self, table):
        return True

    def get_stats(self, short=False):
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

    def get_headers(self):
        headers = ['eigenvalue', 'radians', 'cycle', 'generalized_mass', 'generalized_stiffness']
        return headers

    def build_dataframe(self):
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


class ComplexEigenvalues(BaseScalarObject):
    """
    cycle = freq = eigi / (2*pi)
    radians = eigi
    damping = atan2(eigi, eigr) * 2
    """
    def __init__(self, title, nmodes):
        BaseScalarObject.__init__(self)
        #self.rootNumber = []
        self.title = title
        #self.extraction_order = {}
        #self.eigenvalues = {}
        #self.cycles = {}
        #self.damping = {}

        self.mode = np.zeros(nmodes, dtype='int32')
        self.extraction_order = np.zeros(nmodes, dtype='int32')
        self.eigenvalues = np.zeros(nmodes, dtype='complex64')
        self.cycles = np.zeros(nmodes, dtype='float32')
        self.damping = np.zeros(nmodes, dtype='float32')

        self.data_frame = None

    def __eq__(self, table):
        return True

    def get_stats(self, short=False):
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

    def add_f06_line(self, data, i):
        (root_num, extract_order, eigr, eigi, cycle, damping) = data
        self.mode[i] = root_num
        self.extraction_order[i] = extract_order
        self.eigenvalues[i] = complex(eigr, eigi)
        self.cycles[i] = cycle
        self.damping[i] = damping

    def add_f06_data(self, data):
        for imode, line in enumerate(data):
            self.add_f06_line(line, imode)

    def get_headers(self):
        headers = ['eigenvalue', 'frequency', 'damping']
        return headers

    def build_dataframe(self):
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
    def __init__(self, title, nmodes=0):
        BaseScalarObject.__init__(self)
        self.title = title
        self.mode = np.zeros(nmodes, dtype='int32')
        self.extraction_order = np.zeros(nmodes, dtype='int32')
        self.eigenvalues = np.zeros(nmodes, dtype='float32')
        self.freqs = np.zeros(nmodes, dtype='float32')
        self.omegas = np.zeros(nmodes, dtype='float32')
        self.generalized_mass = np.zeros(nmodes, dtype='float32')
        self.generalized_stiffness = np.zeros(nmodes, dtype='float32')
        self.data_frame = None

    def __eq__(self, table):
        return True

    def get_stats(self, short=False):
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

    def add_f06_line(self, data, imode):
        (root_num, extract_order, eigr, omega, freq, mass, stiff) = data
        self.mode[imode] = root_num
        self.extraction_order[imode] = extract_order
        self.eigenvalues[imode] = eigr
        self.freqs[imode] = freq
        self.omegas[imode] = omega
        self.generalized_mass[imode] = mass
        self.generalized_stiffness[imode] = stiff

      #def add_f06_data(self, data):
        #for i, line in enumerate(data):
            #self.add_f06_line(line, i)

    def build_dataframe(self):
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

    def get_headers(self):
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

    #def __repr__(self):
        #msg = '%-7s %15s %15s %10s %10s %10s\n' % (
            #'RootNum', 'ExtractionOrder', 'Eigenvalue', '', 'Cycles', 'Damping')
        #msg += '%-7s %15s %15s %10s\n' % ('', '', 'Real', 'Imaginary')
        #for root_num, extract_order in sorted(iteritems(self.extraction_order)):
            #eigenvalue = self.eigenvalues[root_num]
            #cycle = self.cycles[root_num]
            #damping = self.damping[root_num]
            #msg += '%-7s %15s %15s %10s %10s %10s\n' % (
                #root_num, extract_order, eigenvalue[0], eigenvalue[1], cycle, damping)
        #return msg
