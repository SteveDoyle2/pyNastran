from struct import pack
from pyNastran.op2.resultObjects.tableObject import RealTableVector, ComplexTableVector, TableObject, ComplexTableObject


class RealVelocityVector(RealTableVector):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableVector.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        words = ['                                                   V E L O C I T Y   V E C T O R\n', ]
                 #' \n',
                 #'      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']

        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(words, header, pageStamp, pageNum, f)
        #words += self.get_table_marker()
        return self._write_f06_block(words, header, pageStamp, pageNum, f)


class ComplexVelocityVector(ComplexTableVector):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableVector.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        words = ['                                       C O M P L E X   V E L O C I T Y   V E C T O R\n']
        return self._write_f06_transient_block(words, header, pageStamp, pageNum, f, is_mag_phase)


class RealVelocityObject(TableObject):  # approach_code=10, thermal=0
    def __init__(self, data_code, is_sort1, isubcase, dt):
        TableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        words = ['                                                   V E L O C I T Y   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        words += self.get_table_marker()
        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(words, header, pageStamp, pageNum, f)
        return self._write_f06_block(words, header, pageStamp, pageNum, f)


class ComplexVelocityObject(ComplexTableObject):  # table_code=10, approach_code=???
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        words = ['                                       C O M P L E X   V E L O C I T Y   V E C T O R\n']
        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(words, header, pageStamp, pageNum, f, is_mag_phase)
        return self._write_f06_block(words, header, pageStamp, pageNum, f, is_mag_phase)
