#pylint: disable=C0301,C0111
from pyNastran.op2.tables.oqg_constraintForces.oqg_spcForces import RealSPCForces  # ,ComplexSPCForcesObject
from pyNastran.op2.tables.oqg_constraintForces.oqg_mpcForces import RealMPCForces  # ,ComplexMPCForcesObject
from pyNastran.op2.tables.opg_appliedLoads.opg_loadVector import RealLoadVector

class OQG(object):
    def _read_f06_table(self, data_types, debug=False):
        pass
    def __init__(self):
        self.iSubcases = []
        self.i = 0

    def _load_vector(self):
        (subcase_name, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
        headers = self.skip(2)
        #print "headers = %s" % (headers)

        data = self._real_f06_real_table_data(allow_blanks=False)

        data_code = {
            'analysis_code': analysis_code,
            'device_code': 1,

            'table_code': 1, # ???
            'table_name': 'OPG',
            'format_code': 1, # ???

            'sort_code': 0,
            'sort_bits': [0, 0, 0], 'num_wide': 8,
            'nonlinear_factor': dt,
            'dataNames':['lsdvmn'],
            'lsdvmn': 1,
            }

        if isubcase in self.load_vectors:
            self.load_vectors[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            spc = RealLoadVector(data_code, is_sort1, isubcase, dt)
            spc.add_f06_data(data, transient)
            self.load_vectors[isubcase] = spc
        self.iSubcases.append(isubcase)

    def _forces_of_single_point_constraints(self):
        (subcase_name, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
        headers = self.skip(2)
        #print "headers = %s" %(headers)

        data = self._real_f06_real_table_data(allow_blanks=False)

        data_code = {
            'analysis_code': analysis_code,
            'device_code': 1, 'table_code': 3, 'sort_code': 0,
            'sort_bits': [0, 0, 0], 'num_wide': 8, 'table_name': 'OQG',
            'nonlinear_factor': dt,
            'format_code': 3,  # ???
            'dataNames':['lsdvmn'],
            'lsdvmn': 1,
        }

        if isubcase in self.spc_forces:
            self.spc_forces[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            spc = RealSPCForces(data_code, is_sort1, isubcase, dt)
            spc.add_f06_data(data, transient)
            self.spc_forces[isubcase] = spc
        self.iSubcases.append(isubcase)

    def _forces_of_multi_point_constraints(self):
        (subcase_name, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
        headers = self.skip(2)
        #print "headers = %s" %(headers)

        data = self._real_f06_real_table_data(allow_blanks=False)

        data_code = {
            'analysis_code': analysis_code,
            'device_code': 1, 'table_code': 39,
            'sort_code': 0, 'sort_bits': [0, 0, 0], 'num_wide': 8,
            'table_name': 'OQG', 'nonlinear_factor': dt,
            'dataNames':['lsdvmn']
        }

        if isubcase in self.mpc_forces:
            self.mpc_forces[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            mpc = RealMPCForces(data_code, is_sort1, isubcase, dt)
            mpc.add_f06_data(data, transient)
            self.mpc_forces[isubcase] = mpc
        self.iSubcases.append(isubcase)
