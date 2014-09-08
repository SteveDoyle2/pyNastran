from struct import pack
from pyNastran.op2.resultObjects.tableObject import RealTableArray, ComplexTableArray, RealTableObject, ComplexTableObject


def make_pack_form(data):
    N = 0
    n = 0
    Form = ''
    old = None
    for d in data:
        if isinstance(d, str):
            n = len(d)
            f = 's'
        elif isinstance(d, int):
            n = 4
            f = 'i'
        elif isinstance(d, float):
            n = 4
            f = 'f'
        else:
            raise NotImplementedError(type(d))
        if old and f != fold:
            form = str(N) + fold
            Form += form
            N = n
        else:
            N += n
        old = d
        fold = f
    if N:
        form = str(N) + f
        Form += form
    return form


class RealDisplacementArray(RealTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
        words = ['                                             D I S P L A C E M E N T   V E C T O R\n', ]
                 #' \n',
                 #'      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #words += self.get_table_marker()
        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(words, header, pageStamp, page_num, f)
        return self._write_f06_block(words, header, pageStamp, page_num, f)


class ComplexDisplacementArray(ComplexTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
        words = ['                                       C O M P L E X   D I S P L A C E M E N T   V E C T O R\n']
        return self._write_f06_transient_block(words, header, pageStamp, page_num, f, is_mag_phase)


class RealDisplacement(RealTableObject):  # approach_code=1, thermal=0
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
        words = ['                                             D I S P L A C E M E N T   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        words += self.get_table_marker()

        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(words, header, pageStamp, page_num, f)
        return self._write_f06_block(words, header, pageStamp, page_num, f)


class ComplexDisplacement(ComplexTableObject):  # approach_code=1, sort_code=0, thermal=0
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
        words = ['                                       C O M P L E X   D I S P L A C E M E N T   V E C T O R\n']
        return self._write_f06_transient_block(words, header, pageStamp, page_num, f, is_mag_phase)
