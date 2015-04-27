from __future__ import print_function
from six import iteritems
from math import sqrt

from numpy import array, pi, allclose

from pyNastran.op2.resultObjects.op2_Objects import BaseScalarObject
from pyNastran.f06.f06_formatting import writeFloats13E
#from pyNastran.op2.resultObjects.op2_Objects import scalarObject,array


class RealEigenvalues(BaseScalarObject):

    def __init__(self, title):
        #self.modeNumber = []
        BaseScalarObject.__init__(self)
        self.title = title
        self.extractionOrder = {}
        self.eigenvalues = {}
        self.radians = {}
        self.cycles = {}
        self.generalizedMass = {}
        self.generalizedStiffness = {}

    def get_stats(self):
        msg = []
        neigenvalues = len(self.extractionOrder)
        msg.append('  type=%s neigenvalues=%s\n' % (self.__class__.__name__,
                                                 neigenvalues))
        msg.append('  title, extractionOrder, eigenvalues, radians, '
                   'cycles, generalizedMass, generalizedStiffness\n')
        return msg

    def isReal(self):
        return True

    def isComplex(self):
        return False

    def addF06Line(self, data):
        (modeNum, extractOrder, eigenvalue, radian, cycle, genM, genK) = data
        #print('data =', data)
        self.extractionOrder[modeNum] = extractOrder
        self.eigenvalues[modeNum] = eigenvalue
        self.radians[modeNum] = radian
        #cyclei = sqrt(abs(eigenvalue)) / (2. * pi)
        #if not allclose(cycle, cyclei):
            #print('cycle=%s cyclei=%s' % (cycle, cyclei))
        self.cycles[modeNum] = cycle
        self.generalizedMass[modeNum] = genM
        self.generalizedStiffness[modeNum] = genK

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
        for (iMode, order) in sorted(iteritems(self.extractionOrder)):
            eigenvalue = self.eigenvalues[iMode]
            #cycle = sqrt(abs(eigenvalue)) / (2. * pi)

            omega = self.radians[iMode]
            freq = self.cycles[iMode]
            mass = self.generalizedMass[iMode]
            stiff = self.generalizedStiffness[iMode]
            ([eigen, omega, freq, mass, stiff], is_all_zeros) = writeFloats13E([eigenvalue, omega, freq, mass, stiff])
            msg.append(' %8s  %8s       %-13s       %-13s       %-13s       %-13s       %s\n' % (iMode, order, eigen, omega, freq, mass, stiff))
        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num

    def __repr__(self):
        msg = '%-7s %15s %15s %10s %10s %10s %15s\n' % ('ModeNum', 'ExtractionOrder', 'Eigenvalue', 'Radians', 'Cycles', 'GenMass', 'GenStiffness')
        for modeNum, extractOrder in sorted(iteritems(self.extractionOrder)):
            eigenvalue = self.eigenvalues[modeNum]
            radian = self.radians[modeNum]

            cycle = sqrt(abs(eigenvalue)) / (2. * pi)
            #cycle = self.cycles[modeNum]
            genM = self.generalizedMass[modeNum]
            genK = self.generalizedStiffness[modeNum]
            msg += '%-7s %15s %15s %10s %10s %10s %15s\n' % (modeNum, extractOrder, eigenvalue, radian, cycle, genM, genK)
        return msg


class ComplexEigenvalues(BaseScalarObject):
    def __init__(self, title):
        BaseScalarObject.__init__(self)
        #self.rootNumber = []
        self.title = title
        self.extractionOrder = {}
        self.eigenvalues = {}
        self.cycles = {}
        self.damping = {}

    def get_stats(self):
        neigenvalues = len(self.extractionOrder)
        msg = []
        msg.append('  type=%s neigenvalues=%s\n' % (self.__class__.__name__, neigenvalues))
        msg.append('  isubcase, extractionOrder, eigenvalues, '
                   'cycles, damping\n')
        return msg

    def isReal(self):
        return False

    def isComplex(self):
        return True

    def addF06Line(self, data):
        (root_num, extract_order, eigr, eigi, cycle, damping) = data
        self.extractionOrder[root_num] = extract_order
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

        for (imode, order) in sorted(iteritems(self.extractionOrder)):
            eigr = self.eigenvalues[imode][0]
            eigi = self.eigenvalues[imode][1]

            freq = self.cycles[imode]
            damping = self.damping[imode]
            ([eigr, eigi, freq, damping], is_all_zeros) = writeFloats13E([eigr, eigi, freq, damping])
            #            imode order      eigr     eigi          freq        damping
            msg.append(' %22s  %10s         %-15s  %-13s         %-13s         %s\n' % (imode, order, eigr, eigi, freq, damping))

        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num

    def __repr__(self):
        msg = '%-7s %15s %15s %10s %10s %10s\n' % ('RootNum', 'ExtractionOrder', 'Eigenvalue', '', 'Cycles', 'Damping')
        msg += '%-7s %15s %15s %10s\n' % ('', '', 'Real', 'Imaginary')
        for imode, extractOrder in sorted(iteritems(self.extractionOrder)):
            eigenvalue = self.eigenvalues[imode]
            cycle = self.cycles[imode]
            damping = self.damping[imode]
            msg += '%-7s %15s %15s %10s %10s %10s\n' % (imode, extractOrder,
                                                        eigenvalue[0], eigenvalue[1], cycle, damping)
        return msg


class BucklingEigenvalues(BaseScalarObject):
    def __init__(self, title):
        BaseScalarObject.__init__(self)
        #self.rootNumber = []
        self.title = title
        self.extractionOrder = {}
        self.eigenvalues = {}
        self.freqs = {}
        self.omegas = {}
        self.generalized_mass = {}
        self.generalized_stiffness = {}
        #self.cycles = {}
        #self.damping = {}

    def get_stats(self):
        neigenvalues = len(self.extractionOrder)
        msg = []
        msg.append('  type=%s neigenvalues=%s\n' % (self.__class__.__name__, neigenvalues))
        msg.append('  isubcase, extractionOrder, eigenvalues, '
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
        (rootNum, extractOrder, eigr, omega, freq, mass, stiff) = data
        self.extractionOrder[rootNum] = extractOrder
        self.eigenvalues[rootNum] = eigr
        self.freqs[rootNum] = freq
        self.omegas[rootNum] = omega
        self.generalized_mass[rootNum] = mass
        self.generalized_stiffness[rootNum] = stiff

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

        for (iMode, order) in sorted(iteritems(self.extractionOrder)):
            eigr = self.eigenvalues[iMode]

            omega = self.omegas[iMode]
            freq = self.freqs[iMode]
            mass = self.generalized_mass[iMode]
            stiff = self.generalized_stiffness[iMode]
            ([eigr, freq, omega, mass, stiff], is_all_zeros) = writeFloats13E([eigr, freq, omega, mass, stiff])
            #            i  ord eig ome f   m          k
            msg.append(' %8s%10s%20s%20s%20s%20s       %s\n' % (iMode, order, eigr, omega, freq, mass, stiff))
        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num

    #def __repr__(self):
        #msg = '%-7s %15s %15s %10s %10s %10s\n' % ('RootNum', 'ExtractionOrder', 'Eigenvalue', '', 'Cycles', 'Damping')
        #msg += '%-7s %15s %15s %10s\n' % ('', '', 'Real', 'Imaginary')
        #for rootNum, extractOrder in sorted(iteritems(self.extractionOrder)):
            #eigenvalue = self.eigenvalues[rootNum]
            #cycle = self.cycles[rootNum]
            #damping = self.damping[rootNum]
            #msg += '%-7s %15s %15s %10s %10s %10s\n' % (rootNum, extractOrder,
                                                        #eigenvalue[0], eigenvalue[1], cycle, damping)
        #return msg
