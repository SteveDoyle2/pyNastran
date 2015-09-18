from struct import pack
from pyNastran.op2.resultObjects.tableObject import RealTableArray, ComplexTableArray, RealTableObject, ComplexTableObject


class RealVelocityArray(RealTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                                                   V E L O C I T Y   V E C T O R\n', ]
                 #' \n',
                 #'      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']

        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(words, header, page_stamp, page_num, f, is_sort1)
        #words += self.get_table_marker()
        return self._write_f06_block(words, header, page_stamp, page_num, f)


class ComplexVelocityArray(ComplexTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                                       C O M P L E X   V E L O C I T Y   V E C T O R\n']
        return self._write_f06_transient_block(words, header, page_stamp, page_num, f, is_mag_phase, is_sort1)


class RealVelocity(RealTableObject):  # approach_code=10, thermal=0
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                                                   V E L O C I T Y   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        words += self.get_table_marker()
        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(words, header, page_stamp, page_num, f, is_sort1)
        return self._write_f06_block(words, header, page_stamp, page_num, f)


class ComplexVelocity(ComplexTableObject):  # table_code=10, approach_code=???
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                                       C O M P L E X   V E L O C I T Y   V E C T O R\n']
        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(words, header, page_stamp, page_num, f, is_mag_phase, is_sort1)
        return self._write_f06_block(words, header, page_stamp, page_num, f, is_mag_phase)
