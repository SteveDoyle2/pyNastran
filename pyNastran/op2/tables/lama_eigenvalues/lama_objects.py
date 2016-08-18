from __future__ import print_function
from six import iteritems
from math import sqrt

import numpy as np

from pyNastran.op2.result_objects.op2_objects import BaseScalarObject
from pyNastran.f06.f06_formatting import write_floats_13e
try:
    import pandas as pd
except ImportError:
    pass


class RealEigenvalues(BaseScalarObject):
    """
    cycle = sqrt(abs(eigenvalue)) / (2. * pi)
    radians = sqrt(abs(eigenvalue))
    """
    def __init__(self, title):
        #self.modeNumber = []
        BaseScalarObject.__init__(self)
        self.title = title
        self.extraction_order = {}
        self.eigenvalues = {}
        self.radians = {}
        self.cycles = {}
        self.generalized_mass = {}
        self.generalized_stiffness = {}

    def __eq__(self, table):
        return True

    def get_stats(self):
        msg = []
        neigenvalues = len(self.extraction_order)
        msg.append('  type=%s neigenvalues=%s\n' % (self.__class__.__name__,
                                                    neigenvalues))
        msg.append('  title, extraction_order, eigenvalues, radians, '
                   'cycles, generalized_mass, generalized_stiffness\n')
        return msg

    def is_real(self):
        return True

    def is_complex(self):
        return False

    def add_f06_line(self, data):
        (mode_num, extract_order, eigenvalue, radian, cycle, gen_mass, gen_stiffness) = data
        self.extraction_order[mode_num] = extract_order
        self.eigenvalues[mode_num] = eigenvalue
        self.radians[mode_num] = radian
        #cyclei = sqrt(abs(eigenvalue)) / (2. * pi)
        #if not allclose(cycle, cyclei):
            #print('cycle=%s cyclei=%s' % (cycle, cyclei))
        self.cycles[mode_num] = cycle
        self.generalized_mass[mode_num] = gen_mass
        self.generalized_stiffness[mode_num] = gen_stiffness

    def get_headers(self):
        headers = ['eigenvalue', 'radians', 'cycle', 'generalized_mass', 'generalized_stiffness']
        return headers

    def build_dataframe(self):
        headers = self.get_headers()
        nmodes = len(self.eigenvalues)

        modes_extraction_order = np.zeros((nmodes, 2), dtype='float32')
        data = np.zeros((nmodes, 5), dtype='float32')

        imodei = 0
        for (imode, eigi) in sorted(iteritems(self.eigenvalues)):
            #cycle = sqrt(abs(eigenvalue)) / (2. * pi)
            extraction_order = self.extraction_order[imode]
            omega = self.radians[imode]
            freq = self.cycles[imode]
            mass = self.generalized_mass[imode]
            stiff = self.generalized_stiffness[imode]
            data[imodei, :] = [eigi, omega, freq, mass, stiff]
            modes_extraction_order[imodei, :] = [imode, extraction_order]
            imodei += 1

        df1 = pd.DataFrame(modes_extraction_order)
        df1.columns = ['Mode', 'ExtractionOrder']
        df2 = pd.DataFrame(data)
        df2.columns = headers
        self.data_frame = df1.join(df2)

    def add_f06_data(self, data):
        for line in data:
            self.add_f06_line(line)

    def write_f06(self, f, header, page_stamp, page_num=1):
        title = ''
        if self.title is not None:
            title = '%s' % str(self.title).center(124).rstrip() + '\n'
        msg = header + ['                                              R E A L   E I G E N V A L U E S\n', title,
                        '   MODE    EXTRACTION      EIGENVALUE            RADIANS             CYCLES            GENERALIZED         GENERALIZED\n',
                        '    NO.       ORDER                                                                       MASS              STIFFNESS\n']
        for (imode, order) in sorted(iteritems(self.extraction_order)):
            eigenvalue = self.eigenvalues[imode]
            #cycle = sqrt(abs(eigenvalue)) / (2. * pi)

            omega = self.radians[imode]
            freq = self.cycles[imode]
            mass = self.generalized_mass[imode]
            stiff = self.generalized_stiffness[imode]
            [eigen, omega, freq, mass, stiff] = write_floats_13e([eigenvalue, omega, freq, mass, stiff])
            msg.append(' %8s  %8s       %-13s       %-13s       %-13s       %-13s       %s\n' % (
                imode, order, eigen, omega, freq, mass, stiff))
        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num

    def __repr__(self):
        msg = '%-7s %15s %15s %10s %10s %10s %15s\n' % (
            'ModeNum', 'ExtractionOrder', 'Eigenvalue', 'Radians', 'Cycles', 'GenMass', 'GenStiffness')
        for mode_num, extract_order in sorted(iteritems(self.extraction_order)):
            eigenvalue = self.eigenvalues[mode_num]
            radian = self.radians[mode_num]

            cycle = sqrt(abs(eigenvalue)) / (2. * np.pi)
            #cycle = self.cycles[mode_num]
            gen_m = self.generalized_mass[mode_num]
            gen_k = self.generalized_stiffness[mode_num]
            msg += '%-7s %15s %15s %10s %10s %10s %15s\n' % (
                mode_num, extract_order, eigenvalue, radian, cycle, gen_m, gen_k)
        return msg


class ComplexEigenvalues(BaseScalarObject):
    """
    cycle = freq = eigi / (2*pi)
    radians = eigi
    damping = atan2(eigi, eigr) * 2
    """
    def __init__(self, title):
        BaseScalarObject.__init__(self)
        #self.rootNumber = []
        self.title = title
        self.extraction_order = {}
        self.eigenvalues = {}
        self.cycles = {}
        self.damping = {}

    def __eq__(self, table):
        return True

    def get_stats(self):
        neigenvalues = len(self.extraction_order)
        msg = []
        msg.append('  type=%s neigenvalues=%s\n' % (self.__class__.__name__, neigenvalues))
        msg.append('  isubcase, extraction_order, eigenvalues, '
                   'cycles, damping\n')
        return msg

    def is_real(self):
        return False

    def is_complex(self):
        return True

    def add_f06_line(self, data):
        (root_num, extract_order, eigr, eigi, cycle, damping) = data
        self.extraction_order[root_num] = extract_order
        self.eigenvalues[root_num] = complex(eigr, eigi)
        self.cycles[root_num] = cycle
        self.damping[root_num] = damping

    def add_f06_data(self, data):
        for line in data:
            self.add_f06_line(line)

    def get_headers(self):
        headers = ['eigenvalue', 'frequency', 'damping']
        return headers

    def build_dataframe(self):
        headers = self.get_headers()
        nmodes = len(self.eigenvalues)

        modes_extraction_order = np.zeros((nmodes, 2), dtype='float32')
        cdata = np.zeros(nmodes, dtype='complex64')
        fdata = np.zeros((nmodes, 2), dtype='float32')

        imodei = 0
        for (imode, eigi) in sorted(iteritems(self.eigenvalues)):
            extraction_order = self.extraction_order[imode]
            freq = self.cycles[imode]
            damping = self.damping[imode]
            cdata[imodei] = eigi
            fdata[imodei, :] = [freq, damping]
            modes_extraction_order[imodei, :] = [imode, extraction_order]
            imodei += 1
        df1 = pd.DataFrame(modes_extraction_order)
        df1.columns = ['Mode', 'ExtractionOrder']
        df2 = pd.DataFrame(cdata)
        df2.columns = [headers[0]]
        df3 = pd.DataFrame(fdata)
        df3.columns = headers[1:]
        self.data_frame = df1.join([df2, df3])
        #print(self.data_frame)

    def write_f06(self, f, header, page_stamp, page_num=1):  # not proper msg start
        title = ''
        if self.title is not None:
            title = '%s' % str(self.title).center(124).rstrip() + '\n'
        msg = header + ['                                        C O M P L E X   E I G E N V A L U E   S U M M A R Y\n', title,
                        '0                ROOT     EXTRACTION                  EIGENVALUE                     FREQUENCY              DAMPING\n',
                        '                  NO.        ORDER             (REAL)           (IMAG)                (CYCLES)            COEFFICIENT\n']

        for (imode, order) in sorted(iteritems(self.extraction_order)):
            eigr = self.eigenvalues[imode].real
            eigi = self.eigenvalues[imode].imag

            freq = self.cycles[imode]
            damping = self.damping[imode]
            [eigr, eigi, freq, damping] = write_floats_13e([eigr, eigi, freq, damping])
            #            imode order      eigr     eigi          freq        damping
            msg.append(' %22s  %10s         %-15s  %-13s         %-13s         %s\n' % (
                imode, order, eigr, eigi, freq, damping))

        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num

    def __repr__(self):
        msg = '%-7s %15s %15s %10s %10s %10s\n' % ('RootNum', 'ExtractionOrder', 'Eigenvalue', '', 'Cycles', 'Damping')
        msg += '%-7s %15s %15s %10s\n' % ('', '', 'Real', 'Imaginary')
        for imode, extract_order in sorted(iteritems(self.extraction_order)):
            eigenvalue = self.eigenvalues[imode]
            cycle = self.cycles[imode]
            damping = self.damping[imode]
            msg += '%-7s %15s %15s %10s %10s %10s\n' % (
                imode, extract_order,
                eigenvalue[0], eigenvalue[1], cycle, damping)
        return msg


class BucklingEigenvalues(BaseScalarObject):
    def __init__(self, title):
        BaseScalarObject.__init__(self)
        self.title = title
        self.extraction_order = {}
        self.eigenvalues = {}
        self.freqs = {}
        self.omegas = {}
        self.generalized_mass = {}
        self.generalized_stiffness = {}

    def __eq__(self, table):
        return True

    def get_stats(self):
        neigenvalues = len(self.extraction_order)
        msg = []
        msg.append('  type=%s neigenvalues=%s\n' % (self.__class__.__name__, neigenvalues))
        msg.append('  imode, extraction_order, eigenvalues, '
                   'radians, cycles, generalized_mass, generalized_stiffness\n')
        return msg

    def is_real(self):
        return False

    def is_complex(self):
        return False

    def is_buckling(self):
        return True

    def add_f06_line(self, data):
        (root_num, extract_order, eigr, omega, freq, mass, stiff) = data
        self.extraction_order[root_num] = extract_order
        self.eigenvalues[root_num] = eigr
        self.freqs[root_num] = freq
        self.omegas[root_num] = omega
        self.generalized_mass[root_num] = mass
        self.generalized_stiffness[root_num] = stiff

      #def add_f06_data(self, data):
        #for line in data:
            #self.add_f06_line(line)

    def build_dataframe(self):
        headers = self.get_headers()
        nmodes = len(self.eigenvalues)

        modes_extraction_order = np.zeros((nmodes, 2), dtype='float32')
        fdata = np.zeros((nmodes, 5), dtype='float32')

        imodei = 0
        for (imode, eigi) in sorted(iteritems(self.eigenvalues)):
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

    def write_f06(self, f, header, page_stamp, page_num=1):  # not proper msg start
        title = ''
        if self.title is not None:
            title = '%s' % str(self.title).center(124).rstrip() + '\n'
        msg = header + ['                                              R E A L   E I G E N V A L U E S\n', title,
                        '   MODE    EXTRACTION      EIGENVALUE            RADIANS             CYCLES            GENERALIZED         GENERALIZED\n',
                        '    NO.       ORDER                                                                       MASS              STIFFNESS\n']

        for (imode, order) in sorted(iteritems(self.extraction_order)):
            eigr = self.eigenvalues[imode]
            freq = self.freqs[imode]
            omega = self.omegas[imode]
            mass = self.generalized_mass[imode]
            stiff = self.generalized_stiffness[imode]
            [eigr, freq, omega, mass, stiff] = write_floats_13e([eigr, freq, omega, mass, stiff])
            #            i  ord eig ome f   m          k
            msg.append(' %8s%10s%20s%20s%20s%20s       %s\n' % (
                imode, order, eigr, omega, freq, mass, stiff))
        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num

    #def __repr__(self):
        #msg = '%-7s %15s %15s %10s %10s %10s\n' % ('RootNum', 'ExtractionOrder', 'Eigenvalue', '', 'Cycles', 'Damping')
        #msg += '%-7s %15s %15s %10s\n' % ('', '', 'Real', 'Imaginary')
        #for root_num, extract_order in sorted(iteritems(self.extraction_order)):
            #eigenvalue = self.eigenvalues[root_num]
            #cycle = self.cycles[root_num]
            #damping = self.damping[root_num]
            #msg += '%-7s %15s %15s %10s %10s %10s\n' % (root_num, extract_order,
                                                        #eigenvalue[0], eigenvalue[1], cycle, damping)
        #return msg
