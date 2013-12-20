from pyNastran.op2.tables.oqg_constraintForces.oqg_spcForces import SPCForcesObject  # ,ComplexSPCForcesObject
from pyNastran.op2.tables.oqg_constraintForces.oqg_mpcForces import MPCForcesObject  # ,ComplexMPCForcesObject


class OQG(object):
    def __init__(self):
        self.spcForces = {}
        self.mpcForces = {}

    def getSpcForces(self):
        (subcaseName, isubcase, transient, dt, analysis_code,
            is_sort1) = self.readSubcaseNameID()
        headers = self.skip(2)
        #print "headers = %s" %(headers)

        dataTypes = [int, str, float, float, float, float, float, float]
        data = self.readTable(dataTypes)

        data_code = {'log': self.log, 'analysis_code': analysis_code,
                    'device_code': 1, 'table_code': 3, 'sort_code': 0,
                    'sort_bits': [0, 0, 0], 'num_wide': 8, 'table_name': 'OQG',
                    'nonlinear_factor': dt, }

        if isubcase in self.spcForces:
            self.spcForces[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            spc = SPCForcesObject(data_code, is_sort1, isubcase)
            spc.add_f06_data(data, transient)
            self.spcForces[isubcase] = spc
        self.iSubcases.append(isubcase)

    def getMpcForces(self):
        (subcaseName, isubcase, transient, dt, analysis_code,
            is_sort1) = self.readSubcaseNameID()
        headers = self.skip(2)
        #print "headers = %s" %(headers)

        dataTypes = [int, str, float, float, float, float, float, float]
        data = self.readTable(dataTypes)

        data_code = {'log': self.log, 'analysis_code': analysis_code,
                    'device_code': 1, 'table_code': 39,
                    'sort_code': 0, 'sort_bits': [0, 0, 0], 'num_wide': 8,
                    'table_name': 'OQG', 'nonlinear_factor': dt, }

        if isubcase in self.mpcForces:
            self.mpcForces[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            mpc = MPCForcesObject(data_code, is_sort1, isubcase)
            mpc.add_f06_data(data, transient)
            self.mpcForces[isubcase] = mpc
        self.iSubcases.append(isubcase)
