from pyNastran.op2.resultObjects.tableObject import RealTableArray, ComplexTableArray


class RealAccelerationArray(RealTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                                             A C C E L E R A T I O N   V E C T O R\n', ]
                 #' \n',
                 #'      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #words += self.get_table_marker()
        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(words, header, page_stamp, page_num, f, write_words=True,
                                                   is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        return self._write_f06_block(words, header, page_stamp, page_num, f, write_words=True)


class ComplexAccelerationArray(ComplexTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                                       C O M P L E X   A C C E L E R A T I O N   V E C T O R\n']
        #words += self.get_table_marker()
        return self._write_f06_transient_block(words, header, page_stamp, page_num, f,
                                               is_mag_phase=is_mag_phase, is_sort1=is_sort1)
