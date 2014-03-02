from struct import pack
from pyNastran.op2.resultObjects.tableObject import RealTableVector, ComplexTableVector, RealTableObject, ComplexTableObject


class RealAccelerationVector(RealTableVector):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableVector.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
        words = ['                                             A C C E L E R A T I O N   V E C T O R\n', ]
                 #' \n',
                 #'      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #words += self.get_table_marker()
        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(words, header, pageStamp, page_num, f)
        return self._write_f06_block(words, header, pageStamp, page_num, f)


class ComplexAccelerationVector(ComplexTableVector):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableVector.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
        words = ['                                       C O M P L E X   A C C E L E R A T I O N   V E C T O R\n']
        #words += self.get_table_marker()
        return self._write_f06_transient_block(words, header, pageStamp, page_num, f, is_mag_phase)


class RealAcceleration(RealTableObject):  # approach_code=11, thermal=0
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, page_num, f)
        words = ['                                             A C C E L E R A T I O N   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        words += self.get_table_marker()
        return self._write_f06_block(words, header, pageStamp, page_num, f)

    def _write_f06_transient(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
        words = ['                                             A C C E L E R A T I O N   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        words += self.get_table_marker()
        return self._write_f06_transient_block(words, header, pageStamp, page_num, f)


class ComplexAcceleration(ComplexTableObject):  # table_code=11, approach_code=???
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        ComplexTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, page_num, f, is_mag_phase)

        words = ['                                       C O M P L E X   A C C E L E R A T I O N   V E C T O R\n']
        return self._write_f06_block(words, header, pageStamp, page_num, f, is_mag_phase)

    def _write_f06_transient(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
        words = ['                                       C O M P L E X   A C C E L E R A T I O N   V E C T O R\n']
        return self._write_f06_transient_block(words, header, pageStamp, page_num, f, is_mag_phase)
