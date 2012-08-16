from pyNastran.op2.tables.oqg_constraintForces.oqg_spcForces import SPCForcesObject  # ,ComplexSPCForcesObject
from pyNastran.op2.tables.oqg_constraintForces.oqg_mpcForces import MPCForcesObject  # ,ComplexMPCForcesObject


class OQG(object):
    def __init__(self):
        self.spcForces = {}
        self.mpcForces = {}

    def getSpcForces(self):
        (subcaseName, iSubcase, transient, dt, analysisCode,
            isSort1) = self.readSubcaseNameID()
        headers = self.skip(2)
        #print "headers = %s" %(headers)

        dataTypes = [int, str, float, float, float, float, float, float]
        data = self.readTable(dataTypes)

        dataCode = {'log': self.log, 'analysisCode': analysisCode,
                    'deviceCode': 1, 'tableCode': 3, 'sortCode': 0,
                    'sortBits': [0, 0, 0], 'numWide': 8, 'tableName': 'OQG',
                    'nonlinearFactor': dt, }

        if iSubcase in self.spcForces:
            self.spcForces[iSubcase].addF06Data(data, transient)
        else:
            isSort1 = True
            spc = SPCForcesObject(dataCode, isSort1, iSubcase)
            spc.addF06Data(data, transient)
            self.spcForces[iSubcase] = spc
        self.iSubcases.append(iSubcase)

    def getMpcForces(self):
        (subcaseName, iSubcase, transient, dt, analysisCode,
            isSort1) = self.readSubcaseNameID()
        headers = self.skip(2)
        #print "headers = %s" %(headers)

        dataTypes = [int, str, float, float, float, float, float, float]
        data = self.readTable(dataTypes)

        dataCode = {'log': self.log, 'analysisCode': analysisCode,
                    'deviceCode': 1, 'tableCode': 39,
                    'sortCode': 0, 'sortBits': [0, 0, 0], 'numWide': 8,
                    'tableName': 'OQG', 'nonlinearFactor': dt, }

        if iSubcase in self.mpcForces:
            self.mpcForces[iSubcase].addF06Data(data, transient)
        else:
            isSort1 = True
            mpc = MPCForcesObject(dataCode, isSort1, iSubcase)
            mpc.addF06Data(data, transient)
            self.mpcForces[iSubcase] = mpc
        self.iSubcases.append(iSubcase)
