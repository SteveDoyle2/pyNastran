import sys
from numpy import array
from pyNastran.op2.resultObjects.op2_Objects import baseScalarObject
#from pyNastran.op2.resultObjects.op2_Objects import scalarObject,array


class RealEigenvalues(baseScalarObject):
    def __init__(self, iSubcase):
        #self.modeNumber = []
        baseScalarObject.__init__(self)
        self.iSubcase = iSubcase
        self.extractionOrder = {}
        self.eigenvalues = {}
        self.radians = {}
        self.cycles = {}
        self.generalizedMass = {}
        self.generalizedStiffness = {}

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

    def addF06Data(self, data):
        for line in data:
            self.addF06Line(line)

    def writeMatlab(self, iSubcase, f=None, isMagPhase=False):
        iModesMsg = 'fem.eigenvalues(%i).iModes    = [' % (iSubcase)
        modesMsg = 'fem.eigenvalues(%i).modes     = [' % (iSubcase)
        orderMsg = 'fem.eigenvalues(%i).order     = [' % (iSubcase)
        omegaMsg = 'fem.eigenvalues(%i).radians   = [' % (iSubcase)
        cyclesMsg = 'fem.eigenvalues(%i).cycles    = [' % (iSubcase)
        massMsg = 'fem.eigenvalues(%i).mass      = [' % (iSubcase)
        stiffMsg = 'fem.eigenvalues(%i).stiffness = [' % (iSubcase)

        for (iMode, order) in sorted(self.extractionOrder.iteritems()):
            iModesMsg += '%s,' % (iMode)
            orderMsg += '%s,' % (order)
            modesMsg += '%s,' % (self.eigenvalues[iMode])
            omegaMsg += '%s,' % (self.radians[iMode])
            cyclesMsg += '%s,' % (self.cycles[iMode])
            massMsg += '%s,' % (self.generalizedMass[iMode])
            stiffMsg += '%s,' % (self.generalizedStiffness[iMode])
        f.write(iModesMsg + '];\n')
        f.write(orderMsg + '];\n')
        f.write(modesMsg + '];\n')
        f.write(omegaMsg + '];\n')
        f.write(cyclesMsg + '];\n')
        f.write(massMsg + '];\n')
        f.write(stiffMsg + '];\n')

    def writeF06(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        msg = header + ['                                              R E A L   E I G E N V A L U E S\n',
                        '   MODE    EXTRACTION      EIGENVALUE            RADIANS             CYCLES            GENERALIZED         GENERALIZED\n',
                        '    NO.       ORDER                                                                       MASS              STIFFNESS\n']
        for (iMode, order) in sorted(self.extractionOrder.iteritems()):
            eigen = self.eigenvalues[iMode]
            omega = self.radians[iMode]
            freq = self.cycles[iMode]
            mass = self.generalizedMass[iMode]
            stiff = self.generalizedStiffness[iMode]
            ([eigen, omega, freq, mass, stiff], isAllZeros) = self.writeFloats13E([eigen, omega, freq, mass, stiff])
            msg.append(' %8s  %8s       %13s       %13s       %13s       %13s       %13s\n' % (iMode, order, eigen, omega, freq, mass, stiff))
        ###
        msg.append(pageStamp + str(pageNum) + '\n')
        return (''.join(msg), pageNum)

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


class ComplexEigenvalues(baseScalarObject):
    def __init__(self, iSubcase):
        baseScalarObject.__init__(self)
        #self.rootNumber = []
        self.iSubcase = iSubcase
        self.extractionOrder = {}
        self.eigenvalues = {}
        self.cycles = {}
        self.damping = {}

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

    def addF06Data(self, data):
        for line in data:
            self.addF06Line(line)

    def writeF06(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):  # not proper msg start
        msg = header + ['                                        C O M P L E X   E I G E N V A L U E S\n',
                        '   MODE    EXTRACTION      EIGENVALUE            CYCLES            DAMPING\n',
                        '    NO.       ORDER\n']
        #raise NotImplementedError()
        for (iMode, order) in sorted(self.extractionOrder.iteritems()):
            eigen = self.eigenvalues[iMode]
            freq = self.cycles[iMode]
            damping = self.damping[iMode]
            ([eigen, freq, damping], isAllZeros) = self.writeFloats13E(
                [eigen, freq, damping])
            msg.append(' %8s  %8s       %13s       %13s       %13s       %13s       %13s\n' % (iMode, order, eigen, freq, damping))
        ###
        msg.append(pageStamp + str(pageNum) + '\n')
        return (''.join(msg), pageNum)

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
