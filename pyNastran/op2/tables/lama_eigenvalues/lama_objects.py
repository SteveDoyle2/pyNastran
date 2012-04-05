class RealEigenvalues(object):
    def __init__(self,iSubcase):
        #self.modeNumber = []
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

    def addF06Data(self,data):
        for line in data:
            (modeNum,extractOrder,eigenvalue,radian,cycle,genM,genK) = line
            self.extractionOrder[modeNum] = extractOrder
            self.eigenvalues[modeNum] = eigenvalue
            self.radians[modeNum] = radian
            self.cycles[modeNum] = cycle
            self.generalizedMass[modeNum] = genM
            self.generalizedStiffness[modeNum] = genK

    def __repr__(self):
        msg = '%-7s %15s %15s %10s %10s %10s %15s\n' %('ModeNum','ExtractionOrder','Eigenvalue','Radians','Cycles','GenMass','GenStiffness')
        for modeNum,extractOrder in sorted(self.extractionOrder.items()):
            eigenvalue = self.eigenvalues[modeNum]
            radian = self.radians[modeNum]
            cycle = self.cycles[modeNum]
            genM = self.generalizedMass[modeNum]
            genK = self.generalizedStiffness[modeNum]
            msg += '%-7s %15s %15s %10s %10s %10s %15s\n' %(modeNum,extractOrder,eigenvalue,radian,cycle,genM,genK)
        return msg

class ComplexEigenvalues(object):
    def __init__(self,iSubcase):
        #self.rootNumber = []
        self.iSubcase = iSubcase
        self.extractionOrder = {}
        self.eigenvalues = {}
        self.cycles  = {}
        self.damping = {}

    def isReal(self):
        return False
    def isComplex(self):
        return True

    def addF06Data(self,data):
        for line in data:
            (rootNum,extractOrder,eigr,eigi,cycle,damping) = line
            self.extractionOrder[rootNum] = extractOrder
            self.eigenvalues[rootNum] = array([eigr,eigi])
            self.cycles[rootNum]      = cycle
            self.damping[rootNum]     = damping

    def __repr__(self):
        msg  = '%-7s %15s %15s %10s %10s %10s\n' %('RootNum','ExtractionOrder','Eigenvalue','','Cycles','Damping')
        msg += '%-7s %15s %15s %10s\n' %('','','Real','Imaginary')
        for rootNum,extractOrder in sorted(self.extractionOrder.items()):
            eigenvalue = self.eigenvalues[rootNum]
            cycle      = self.cycles[rootNum]
            damping    = self.damping[rootNum]
            msg += '%-7s %15s %15s %10s %10s %10s\n' %(rootNum,extractOrder,eigenvalue[0],eigenvalue[1],cycle,damping)
        return msg
