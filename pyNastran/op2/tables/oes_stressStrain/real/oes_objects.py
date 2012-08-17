from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from numpy import argsort

from pyNastran.op2.resultObjects.op2_Objects import scalarObject


class OES_Object(scalarObject):
    def __init__(self, dataCode, iSubcase):
        scalarObject.__init__(self, dataCode, iSubcase)
        self.log.debug("starting OES...elementName=%s iSubcase=%s" %
                       (self.elementName, self.iSubcase))
        #print self.dataCode

    def isCurvatureOld(self):
        if self.stressBits[2] == 0:
            return True
        return False

    def isCurvature(self):
        if self.sCode in [0, 1, 14, 15, 16, 17, 27, 30, 31]:  # fiber distance
            return False
        elif self.sCode in [10, 11, 26, ]:  # fiber curvature
            return True
        raise NotImplementedError('add sCode=%s' % (self.sCode))

    def isFiberDistance(self):
        return not(self.isCurvature())

    def isVonMises(self):
        #print self.stressBits
        #iMs = not(self.isMaxShear())
        #print 'isVonMises = ',iMs
        return not(self.isMaxShear())

    def isMaxShear(self):
        #print self.stressBits
        if self.stressBits[4] == 0:
            #print 'isMaxShear = True'
            return True
        #print 'isMaxShear = False'
        return False

    def getOrderedETypes(self, validTypes):
        """
        @param validTypes list of valid element types
               e.g. ['CTRIA3','CTRIA6','CQUAD4','CQUAD8']

        @retval TypesOut the ordered list of types
        @retval orderedETypes dictionary of Type-IDs to write
        """
        orderedETypes = {}

        #validTypes = ['CTRIA3','CTRIA6','CQUAD4','CQUAD8']
        for eType in validTypes:
            orderedETypes[eType] = []
        for eid, eType in sorted(self.eType.items()):
            #print "eType = ",eType
            assert eType in validTypes, 'unsupported eType=%s' % (eType)
            orderedETypes[eType].append(eid)

        minVals = []
        for eType in validTypes:
            vals = orderedETypes[eType]
            #print "len(%s) = %s" %(eType,len(vals))
            if len(vals) == 0:
                minVals.append(-1)
            else:
                minVals.append(min(vals))

        #print "minVals = ",minVals
        argList = argsort(minVals)

        TypesOut = []
        for i in argList:
            TypesOut.append(validTypes[i])
        #print "validTypes = %s" %(validTypes)
        #print "minVals    = %s" %(minVals)
        #print "argList    = %s" %(argList)
        #print "TypesOut   = %s" %(TypesOut)
        #print "orderedETypes.keys = %s" %(orderedETypes.keys())
        return (TypesOut, orderedETypes)


class stressObject(OES_Object):
    def __init__(self, dataCode, iSubcase):
        OES_Object.__init__(self, dataCode, iSubcase)

    def updateDt(self, dataCode, dt):
        self.dataCode = dataCode
        self.applyDataCode()
        #assert dt>=0.
        #print "dataCode=",self.dataCode
        self.elementName = self.dataCode['elementName']
        if dt is not None:
            self.log.debug("updating stress...%s=%s elementName=%s" %
                           (self.dataCode['name'], dt, self.elementName))
            self.dt = dt
            self.addNewTransient(dt)
        ###

    def isStrain(self):
        return True

    def isStress(self):
        return False


class strainObject(OES_Object):
    def __init__(self, dataCode, iSubcase):
        OES_Object.__init__(self, dataCode, iSubcase)

    def updateDt(self, dataCode, dt):
        self.dataCode = dataCode
        self.applyDataCode()
        #print "dataCode=",self.dataCode
        self.elementName = self.dataCode['elementName']
        #assert dt>=0.
        if dt is not None:
            self.log.debug("updating strain...%s=%s elementName=%s" %
                           (self.dataCode['name'], dt, self.elementName))
            self.dt = dt
            self.addNewTransient()

    def isStress(self):
        return False

    def isStrain(self):
        return True
