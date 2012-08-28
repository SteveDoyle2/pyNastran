from numpy import array

from pyNastran.op2.resultObjects.op2_Objects import scalarObject
from pyNastran.op2.resultObjects.tableObject import TableObject, ComplexTableObject


class EigenVectorObject(TableObject):  # approachCode=2, sortCode=0, thermal=0
    """
    @code
    EIGENVALUE =  6.158494E+07
        CYCLES =  1.248985E+03         R E A L   E I G E N V E C T O R   N O .          1

    POINT ID.   TYPE          T1             T2             T3             R1             R2             R3
           1      G      2.547245E-17  -6.388945E-16   2.292728E+00  -1.076928E-15   2.579163E-17   0.0
        2002      G     -6.382321E-17  -1.556607E-15   3.242408E+00  -6.530917E-16   1.747180E-17   0.0
        2003      G      0.0            0.0            0.0            0.0            0.0            0.0
    @endcode
    """
    def __init__(self, dataCode, isSort1, iSubcase, iMode):
        TableObject.__init__(self, dataCode, isSort1, iSubcase, iMode)
        #self.caseVal = mode
        self.updateDt = self.updateMode
        #print "mode = %s" %(mode)
        #print "dataCode = ",self.dataCode
        self.setDataMembers()

        #assert mode>=0.
        self.gridTypes = {}
        #self.translations = {iMode: {}}
        #self.rotations    = {iMode: {}}

    def readF06Data(self, dataCode, data):
        iMode = dataCode['mode']
        if iMode not in self.translations:
            self.updateMode(dataCode, iMode)

        for line in data:
            (nid, gridType, t1, t2, t3, r1, r2, r3) = line
            self.gridTypes[nid] = gridType
            self.translations[iMode][nid] = array([t1, t2, t3])
            self.rotations[iMode][nid] = array([r1, r2, r3])
        ###

    def updateMode(self, dataCode, iMode):
        """
        this method is called if the object
        already exits and a new time step is found
        """
        #assert mode>=0.
        self.dataCode = dataCode
        self.applyDataCode()
        #self.caseVal = iMode
        #print "mode = %s" %(str(mode))
        self.addNewMode(iMode)
        self.setDataMembers()

    def addNewMode(self, iMode):
        self.translations[iMode] = {}
        self.rotations[iMode] = {}

    def eigenvalues(self):
        return self.eigrs

    def writeMatlab(self, iSubcase, f=None, isMagPhase=False):
        name = 'eigenvectors'
        return self._writeMatlabTransient(name, iSubcase, f, isMagPhase)

    def writeF06(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        """
        @code
        EIGENVALUE =  6.158494E+07
            CYCLES =  1.248985E+03         R E A L   E I G E N V E C T O R   N O .          1

        POINT ID.   TYPE          T1             T2             T3             R1             R2             R3
               1      G      2.547245E-17  -6.388945E-16   2.292728E+00  -1.076928E-15   2.579163E-17   0.0
            2002      G     -6.382321E-17  -1.556607E-15   3.242408E+00  -6.530917E-16   1.747180E-17   0.0
            2003      G      0.0            0.0            0.0            0.0            0.0            0.0
        @endcode
        """
        msg = []
        hasCycle = hasattr(self, 'modeCycle')

        for i, (iMode, eigVals) in enumerate(sorted(self.translations.iteritems())):
            msg += header
            freq = self.eigrs[i]
            msg.append('%16s = %13E\n' % ('EIGENVALUE', freq))

            if hasCycle:
                msg.append('%16s = %13E         R E A L   E I G E N V E C T O R   N O . %10i\n \n' % ('CYCLES', self.modeCycle, iMode))
            else:
                msg.append('                                         R E A L   E I G E N V E C T O R   N O . %10i\n \n' % (iMode))

            msg.append('      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n')
            for nodeID, displacement in sorted(eigVals.iteritems()):
                rotation = self.rotations[iMode][nodeID]
                gridType = self.gridTypes[nodeID]
                (dx, dy, dz) = displacement
                (rx, ry, rz) = rotation

                vals = [dx, dy, dz, rx, ry, rz]
                (vals2, isAllZeros) = self.writeFloats13E(vals)
                [dx, dy, dz, rx, ry, rz] = vals2
                msg.append('%14i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' % (nodeID, gridType, dx, dy, dz, rx, ry, rz.rstrip()))
            ###
            msg.append(pageStamp + str(pageNum) + '\n')
            if f is not None:
                f.write(''.join(msg))
                msg = ['']
            pageNum += 1
        ###
        return (''.join(msg), pageNum - 1)

    def __repr__(self):
        msg = '---EIGENVECTORS---\n'
        msg += self.printDataMembers()
        name = self.dataCode['name']

        headers = ['Tx', 'Ty', 'Tz', 'Rx', 'Ry', 'Rz']
        headerLine = '%-8s %8s ' % ('nodeID', 'GridType',)
        for header in headers:
            headerLine += '%10s ' % (header)
        headerLine += '\n'

        for i, (iMode, eigVals) in enumerate(sorted(self.translations.iteritems())):
            freq = self.eigrs[i]
            msg += '%s = %g\n' % (name, iMode)
            msg += 'eigenvalueReal = %g\n' % (freq)
            msg += headerLine
            for nodeID, displacement in sorted(eigVals.iteritems()):
                rotation = self.rotations[iMode][nodeID]
                gridType = self.gridTypes[nodeID]
                (dx, dy, dz) = displacement
                (rx, ry, rz) = rotation

                msg += '%-8i %8s ' % (nodeID, gridType)
                vals = [dx, dy, dz, rx, ry, rz]
                for val in vals:
                    if abs(val) < 1e-6:
                        msg += '%10s ' % (0)
                    else:
                        msg += '%10.3g ' % (val)
                    ###
                msg += '\n'
            ###
            msg += '\n'
            #print msg
        return msg


class realEigenVectorObject(scalarObject):  # approachCode=2, sortCode=0, thermal=0
    """
    @code
                                         R E A L   E I G E N V E C T O R   N O .          1
      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3
             1      G      0.0            0.0            0.0            0.0            1.260264E-01   0.0
             7      G      9.999849E-01   0.0            6.728968E-03   0.0            8.021386E-03   0.0
    @endcode
    """
    def __init__(self, dataCode, iSubcase, iMode):
        scalarObject.__init__(self, dataCode, iSubcase)
        #self.caseVal = mode
        #print "mode = %s" %(iMode)
        self.caseVal = self.getUnsteadyValue()
        self.setDataMembers()

        #assert mode>=0.
        self.gridTypes = {}
        self.translations = {iMode: {}}
        self.rotations = {iMode: {}}

    def addNewMode(self, iMode):
        self.translations[iMode] = {}
        self.rotations[iMode] = {}

    def updateDt(self, dataCode, dt):
        #print " self.dataCode = ",self.dataCode
        self.dataCode = dataCode
        self.applyDataCode()
        self.setDataMembers()
        self.caseVal = dt

        #print "*self.dataCode = ",self.dataCode
        self.translations[self.caseVal] = {}
        self.rotations[self.caseVal] = {}

    def deleteTransient(self, dt):
        del self.translations[dt]
        del self.rotations[dt]

    def getTransients(self):
        k = self.translations.keys()
        k.sort()
        return k

    def add(self, nodeID, gridType, v1, v2, v3, v4, v5, v6):
        msg = "nodeID=%s v1=%s v2=%s v3=%s" % (nodeID, v1, v2, v3)
        msg += "           v4=%s v5=%s v6=%s" % (v4, v5, v6)
        #print msg
        assert 0 < nodeID < 1000000000, msg
        #assert nodeID not in self.translations

        if gridType == 1:
            Type = 'G'
        elif gridType == 2:
            Type = 'S'
        elif gridType == 7:
            Type = 'L'
        else:
            raise ValueError('invalid grid type...gridType=%s' % (gridType))

        self.gridTypes[nodeID] = Type
        #print 'self.caseVal = %s' %(self.caseVal),type(self.caseVal)
        #print "d = ",self.translations
        self.translations[self.caseVal][nodeID] = [v1, v2, v3]
        self.rotations[self.caseVal][nodeID] = [v4, v5, v6]
    ###

    def modes(self):
        return sorted(self.translations.keys())

    def eigenvalues(self):
        return self.eigrs

    def writeMatlab(self, iSubcase, f=None, isMagPhase=False):
        name = 'eigenvectors'
        return self._writeMatlabTransient(name, iSubcase, f, isMagPhase)

    def writeF06(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        """
        @code
        EIGENVALUE =  6.158494E+07
                                           R E A L   E I G E N V E C T O R   N O .          1

        POINT ID.   TYPE          T1             T2             T3             R1             R2             R3
               1      G      2.547245E-17  -6.388945E-16   2.292728E+00  -1.076928E-15   2.579163E-17   0.0
            2002      G     -6.382321E-17  -1.556607E-15   3.242408E+00  -6.530917E-16   1.747180E-17   0.0
            2003      G      0.0            0.0            0.0            0.0            0.0            0.0
        @endcode
        """
        msg = []
        #print self.dataCode
        for i, (iMode, eigVals) in enumerate(sorted(self.translations.iteritems())):
            msg += header
            freq = self.eigrs[i]
            msg.append('%16s = %12E\n' % ('EIGENVALUE', freq))
            msg.append('                                         R E A L   E I G E N V E C T O R   N O . %10i\n \n' % (iMode))
            msg.append('      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n')
            for nodeID, translation in sorted(eigVals.iteritems()):
                rotation = self.rotations[iMode][nodeID]
                gridType = self.gridTypes[nodeID]
                (dx, dy, dz) = translation
                (rx, ry, rz) = rotation

                vals = [dx, dy, dz, rx, ry, rz]
                (vals2, isAllZeros) = self.writeFloats13E(vals)
                [dx, dy, dz, rx, ry, rz] = vals2
                msg.append('%14i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' % (nodeID, gridType, dx, dy, dz, rx, ry, rz.rstrip()))
            ###
            msg.append(pageStamp + str(pageNum) + '\n')
            pageNum += 1
        ###
        return (''.join(msg), pageNum - 1)

    def __repr__(self):
        msg = '---REAL EIGENVECTORS---\n'
        msg += self.printDataMembers()
        name = self.dataCode['name']

        headers = ['T']
        headerLine = '%-8s %8s ' % ('nodeID', 'GridType',)
        for header in headers:
            headerLine += '%10s ' % (header)
        headerLine += '\n'

        for iMode, eigVals in sorted(self.translations.iteritems()):
            msg += '%s = %s\n' % (name, iMode)
            msg += headerLine
            for nodeID, translation in sorted(eigVals.iteritems()):
                Type = self.gridTypes[nodeID]

                rotation = self.rotations[iMode][nodeID]
                (dx, dy, dz) = translation
                (rx, ry, rz) = rotation

                vals = [dx, dy, dz, rx, ry, rz]
                msg += '%-8i %8s ' % (nodeID, Type)
                for v in vals:
                    if abs(v) < 1e-6:
                        msg += '%10s ' % (0)
                    else:
                        msg += '%10.3f ' % (v)
                msg += '\n'
            msg += '\n'
        return msg


class ComplexEigenVectorObject(ComplexTableObject):  # approachCode=2, sortCode=0, thermal=0
    def __init__(self, dataCode, isSort1, iSubcase, iMode):
        ComplexTableObject.__init__(self, dataCode, isSort1, iSubcase, iMode)
        self.caseVal = iMode
        self.updateDt = self.updateMode
        self.setDataMembers()

        #print "mode = %s" %(mode)

        #assert mode>=0.
        #self.gridTypes = {}
        #self.translations = {iMode: {}}
        #self.rotations    = {iMode: {}}

    def updateMode(self, dataCode, iMode):
        """
        this method is called if the object
        already exits and a new time step is found
        """
        #assert mode>=0.
        self.caseVal = iMode
        self.dataCode = dataCode
        self.applyDataCode()
        self.setDataMembers()
        #print "mode = %s" %(str(mode))
        self.addNewMode(iMode)

    def addNewMode(self, iMode):
        self.translations[iMode] = {}
        self.rotations[iMode] = {}

    def eigenvalues(self):
        return sorted(self.translations.keys())

    def writeF06(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        msg = []
        #print self.dataCode
        hasCycle = hasattr(self, 'modeCycle')
        for i, (iMode, eigVals) in enumerate(sorted(self.translations.iteritems())):
            msg += header
            freq = self.eigrs[i]
            #freq = 0.0
            msg.append('%16s = %12E\n' % ('EIGENVALUE', freq))
            if hasCycle:
                msg.append('%16s = %12E          C O M P L E X   E I G E N V E C T O R   N O . %10i\n \n' % ('CYCLES', self.modeCycle, iMode))
            else:
                msg.append('                                         C O M P L E X   E I G E N V E C T O R   N O . %10i\n \n' % (iMode))
            msg.append('      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n')
            for nodeID, displacement in sorted(eigVals.iteritems()):
                rotation = self.rotations[iMode][nodeID]
                gridType = self.gridTypes[nodeID]
                (dx, dy, dz) = displacement
                (rx, ry, rz) = rotation

                vals = [dx, dy, dz, rx, ry, rz]
                (vals2, isAllZeros) = self.writeImagFloats13E(vals, isMagPhase)
                [dxr, dyr, dzr, rxr, ryr, rzr, dxi, dyi,
                    dzi, rxi, ryi, rzi] = vals2

                msg.append('%14i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' % (nodeID, gridType, dxr, dyr, dzr, rxr, ryr, rzr.rstrip()))
                msg.append('%14s %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' % ('', '', dxi, dyi, dzi, rxi, ryi, rzi.rstrip()))
            ###
            msg.append(pageStamp + str(pageNum) + '\n')
            if f is not None:
                f.write(''.join(msg))
                msg = ['']
            pageNum += 1
        ###
        return (''.join(msg), pageNum - 1)

    def __repr__(self):
        msg = '---EIGENVECTORS---\n'
        msg += self.printDataMembers()

        headers = ['T1', 'T2', 'T3', 'R1', 'R2', 'R3']
        headerLine = '%-8s %8s ' % ('nodeID', 'GridType',)
        for header in headers:
            headerLine += '%10s ' % (header)
        headerLine += '\n'
        name = self.dataCode['name']

        for i, (iMode, eigVals) in enumerate(sorted(self.translations.iteritems())):
            msg += '%s = %g\n' % (name, iMode)
            msg += headerLine
            for nodeID, translation in sorted(eigVals.iteritems()):
                rotation = self.rotations[iMode][nodeID]
                Type = self.gridTypes[nodeID]
                #(dx,dy,dz) = displacement
                #(rx,ry,rz) = rotation

                msg += '%-8i %8s ' % (nodeID, Type)
                vals = translation + rotation
                for val in vals:
                    if abs(val) < 1e-6:
                        msg += '%10s ' % (0)
                    else:
                        msg += '%10s ' % (str(val))
                        #msg += '%10.3g ' %(val)
                msg += '\n'
            msg += '\n'
        return msg
