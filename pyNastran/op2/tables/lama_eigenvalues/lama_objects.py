from __future__ import print_function
from six import iteritems
from math import sqrt

from numpy import array, pi

from pyNastran.op2.resultObjects.op2_Objects import BaseScalarObject
from pyNastran.f06.f06_formatting import write_floats_13e
#from pyNastran.op2.resultObjects.op2_Objects import scalarObject,array


class RealEigenvalues(BaseScalarObject):

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

    def get_stats(self):
        msg = []
        neigenvalues = len(self.extraction_order)
        msg.append('  type=%s neigenvalues=%s\n' % (self.__class__.__name__,
                                                    neigenvalues))
        msg.append('  title, extraction_order, eigenvalues, radians, '
                   'cycles, generalized_mass, generalized_stiffness\n')
        return msg

    def isReal(self):
        return True

    def isComplex(self):
        return False

    def addF06Line(self, data):
        (modeNum, extract_order, eigenvalue, radian, cycle, gen_mass, gen_stiffness) = data
        #print('data =', data)
        self.extraction_order[modeNum] = extract_order
        self.eigenvalues[modeNum] = eigenvalue
        self.radians[modeNum] = radian
        #cyclei = sqrt(abs(eigenvalue)) / (2. * pi)
        #if not allclose(cycle, cyclei):
            #print('cycle=%s cyclei=%s' % (cycle, cyclei))
        self.cycles[modeNum] = cycle
        self.generalized_mass[modeNum] = gen_mass
        self.generalized_stiffness[modeNum] = gen_stiffness

    def add_f06_data(self, data):
        for line in data:
            self.addF06Line(line)

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

            cycle = sqrt(abs(eigenvalue)) / (2. * pi)
            #cycle = self.cycles[mode_num]
            genM = self.generalized_mass[mode_num]
            genK = self.generalized_stiffness[mode_num]
            msg += '%-7s %15s %15s %10s %10s %10s %15s\n' % (
                mode_num, extract_order, eigenvalue, radian, cycle, genM, genK)
        return msg


class ComplexEigenvalues(BaseScalarObject):
    def __init__(self, title):
        BaseScalarObject.__init__(self)
        #self.rootNumber = []
        self.title = title
        self.extraction_order = {}
        self.eigenvalues = {}
        self.cycles = {}
        self.damping = {}

    def get_stats(self):
        neigenvalues = len(self.extraction_order)
        msg = []
        msg.append('  type=%s neigenvalues=%s\n' % (self.__class__.__name__, neigenvalues))
        msg.append('  isubcase, extraction_order, eigenvalues, '
                   'cycles, damping\n')
        return msg

    def isReal(self):
        return False

    def isComplex(self):
        return True

    def addF06Line(self, data):
        (root_num, extract_order, eigr, eigi, cycle, damping) = data
        self.extraction_order[root_num] = extract_order
        self.eigenvalues[root_num] = array([eigr, eigi])
        self.cycles[root_num] = cycle
        self.damping[root_num] = damping

    def add_f06_data(self, data):
        for line in data:
            self.addF06Line(line)

    def write_f06(self, f, header, page_stamp, page_num=1):  # not proper msg start
        title = ''
        if self.title is not None:
            title = '%s' % str(self.title).center(124).rstrip() + '\n'
        msg = header + ['                                        C O M P L E X   E I G E N V A L U E   S U M M A R Y\n', title,
                        '0                ROOT     EXTRACTION                  EIGENVALUE                     FREQUENCY              DAMPING\n',
                        '                  NO.        ORDER             (REAL)           (IMAG)                (CYCLES)            COEFFICIENT\n']

        for (imode, order) in sorted(iteritems(self.extraction_order)):
            eigr = self.eigenvalues[imode][0]
            eigi = self.eigenvalues[imode][1]

            freq = self.cycles[imode]
            damping = self.damping[imode]
            [eigr, eigi, freq, damping] = write_floats_13e([eigr, eigi, freq, damping])
            #            imode order      eigr     eigi          freq        damping
            msg.append(' %22s  %10s         %-15s  %-13s         %-13s         %s\n' % (imode, order, eigr, eigi, freq, damping))

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
            msg += '%-7s %15s %15s %10s %10s %10s\n' % (imode, extract_order,
                                                        eigenvalue[0], eigenvalue[1], cycle, damping)
        return msg


class BucklingEigenvalues(BaseScalarObject):
    def __init__(self, title):
        BaseScalarObject.__init__(self)
        #self.rootNumber = []
        self.title = title
        self.extraction_order = {}
        self.eigenvalues = {}
        self.freqs = {}
        self.omegas = {}
        self.generalized_mass = {}
        self.generalized_stiffness = {}
        #self.cycles = {}
        #self.damping = {}

    def get_stats(self):
        neigenvalues = len(self.extraction_order)
        msg = []
        msg.append('  type=%s neigenvalues=%s\n' % (self.__class__.__name__, neigenvalues))
        msg.append('  isubcase, extraction_order, eigenvalues, '
                   'cycles, damping\n')
        return msg

    def isReal(self):
        return False

    def isComplex(self):
        return False

    def isBuckling(self):
        return True

    def addF06Line(self, data):
        #print('data =', data)
         #(iMode, order, eigen, omega, freq, mass, stiff)
        (root_num, extract_order, eigr, omega, freq, mass, stiff) = data
        self.extraction_order[root_num] = extract_order
        self.eigenvalues[root_num] = eigr
        self.freqs[root_num] = freq
        self.omegas[root_num] = omega
        self.generalized_mass[root_num] = mass
        self.generalized_stiffness[root_num] = stiff

      #def add_f06_data(self, data):
        #for line in data:
            #self.addF06Line(line)

    def write_f06(self, f, header, page_stamp, page_num=1):  # not proper msg start
        title = ''
        if self.title is not None:
            title = '%s' % str(self.title).center(124).rstrip() + '\n'
        msg = header + ['                                              R E A L   E I G E N V A L U E S\n', title,
                        '   MODE    EXTRACTION      EIGENVALUE            RADIANS             CYCLES            GENERALIZED         GENERALIZED\n',
                        '    NO.       ORDER                                                                       MASS              STIFFNESS\n']

        for (imode, order) in sorted(iteritems(self.extraction_order)):
            eigr = self.eigenvalues[imode]

            omega = self.omegas[imode]
            freq = self.freqs[imode]
            mass = self.generalized_mass[imode]
            stiff = self.generalized_stiffness[imode]
            [eigr, freq, omega, mass, stiff] = write_floats_13e([eigr, freq, omega, mass, stiff])
            #            i  ord eig ome f   m          k
            msg.append(' %8s%10s%20s%20s%20s%20s       %s\n' % (imode, order, eigr, omega, freq, mass, stiff))
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
