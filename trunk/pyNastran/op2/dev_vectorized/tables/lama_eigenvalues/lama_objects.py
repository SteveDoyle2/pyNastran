from numpy import array
from pyNastran.op2.resultObjects.op2_Objects import BaseScalarObject
from pyNastran.f06.f06_formatting import writeFloats13E
#from pyNastran.op2.resultObjects.op2_Objects import scalarObject,array


class RealEigenvalues(BaseScalarObject):

    def __init__(self, isubcase):
        #self.modeNumber = []
        BaseScalarObject.__init__(self)
        self. subcase = isubcase
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
        msg.append('  isubcase, extractionOrder, eigenvalues, radians, '
                   'cycles, generalizedMass, generalizedStiffness\n')
        return msg

    def isReal(self):
        return True

    def isComplex(self):
        return False

    def addF06Line(self, data):
        (modeNum, extractOrder, eigenvalue, radian, cycle, genM, genK) = data
        self.extractionOrder[modeNum] = extractOrder
        self.eigenvalues[modeNum] = eigenvalue
        self.radians[modeNum] = radian
        self.cycles[modeNum] = cycle
        self.generalizedMass[modeNum] = genM
        self.generalizedStiffness[modeNum] = genK

    def add_f06_data(self, data):
        for line in data:
            self.addF06Line(line)

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        msg = header + ['                                              R E A L   E I G E N V A L U E S\n',
                        '   MODE    EXTRACTION      EIGENVALUE            RADIANS             CYCLES            GENERALIZED         GENERALIZED\n',
                        '    NO.       ORDER                                                                       MASS              STIFFNESS\n']
        for (iMode, order) in sorted(self.extractionOrder.iteritems()):
            eigen = self.eigenvalues[iMode]
            omega = self.radians[iMode]
            freq = self.cycles[iMode]
            mass = self.generalizedMass[iMode]
            stiff = self.generalizedStiffness[iMode]
            ([eigen, omega, freq, mass, stiff], isAllZeros) = writeFloats13E([eigen, omega, freq, mass, stiff])
            msg.append(' %8s  %8s       %13s       %13s       %13s       %13s       %13s\n' % (iMode, order, eigen, omega, freq, mass, stiff))

        msg.append(pageStamp + str(pageNum) + '\n')
        f.write(''.join(msg))
        return pageNum

    def __repr__(self):
        msg = '%-7s %15s %15s %10s %10s %10s %15s\n' % ('ModeNum', 'ExtractionOrder', 'Eigenvalue', 'Radians', 'Cycles', 'GenMass', 'GenStiffness')
        for modeNum, extractOrder in sorted(self.extractionOrder.iteritems()):
            eigenvalue = self.eigenvalues[modeNum]
            radian = self.radians[modeNum]
            cycle = self.cycles[modeNum]
            genM = self.generalizedMass[modeNum]
            genK = self.generalizedStiffness[modeNum]
            msg += '%-7s %15s %15s %10s %10s %10s %15s\n' % (modeNum, extractOrder, eigenvalue, radian, cycle, genM, genK)
        return msg


class ComplexEigenvalues(BaseScalarObject):
    def __init__(self, isubcase):
        BaseScalarObject.__init__(self)
        #self.rootNumber = []
        self. subcase = isubcase
        self.extractionOrder = {}
        self.eigenvalues = {}
        self.cycles = {}
        self.damping = {}

    def get_stats(self):
        neigenvalues = len(self.extractionOrder)
        msg = []
        msg.append('  type=%s neigenvalues=%s\n' % (self.__class__.__name__,
                                                 neigenvalues))
        msg.append('  isubcase, extractionOrder, eigenvalues, '
                   'cycles, damping\n')
        return msg

    def isReal(self):
        return False

    def isComplex(self):
        return True

    def addF06Line(self, data):
        (rootNum, extractOrder, eigr, eigi, cycle, damping) = data
        self.extractionOrder[rootNum] = extractOrder
        self.eigenvalues[rootNum] = array([eigr, eigi])
        self.cycles[rootNum] = cycle
        self.damping[rootNum] = damping

    def add_f06_data(self, data):
        for line in data:
            self.addF06Line(line)

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):  # not proper msg start
        msg = header + ['                                        C O M P L E X   E I G E N V A L U E S\n',
                        '   MODE    EXTRACTION      EIGENVALUE            CYCLES            DAMPING\n',
                        '    NO.       ORDER\n']
        #raise NotImplementedError()
        for (iMode, order) in sorted(self.extractionOrder.iteritems()):
            eigen = self.eigenvalues[iMode]
            freq = self.cycles[iMode]
            damping = self.damping[iMode]
            ([eigen, freq, damping], isAllZeros) = writeFloats13E([eigen, freq, damping])
            msg.append(' %8s  %8s       %13s       %13s       %13s       %13s       %13s\n' % (iMode, order, eigen, freq, damping))

        msg.append(pageStamp + str(pageNum) + '\n')
        f.write(''.join(msg))
        return pageNum

    def __repr__(self):
        msg = '%-7s %15s %15s %10s %10s %10s\n' % ('RootNum',
                                                   'ExtractionOrder', 'Eigenvalue', '', 'Cycles', 'Damping')
        msg += '%-7s %15s %15s %10s\n' % ('', '', 'Real', 'Imaginary')
        for rootNum, extractOrder in sorted(self.extractionOrder.iteritems()):
            eigenvalue = self.eigenvalues[rootNum]
            cycle = self.cycles[rootNum]
            damping = self.damping[rootNum]
            msg += '%-7s %15s %15s %10s %10s %10s\n' % (rootNum, extractOrder,
                                                        eigenvalue[0], eigenvalue[1], cycle, damping)
        return msg
