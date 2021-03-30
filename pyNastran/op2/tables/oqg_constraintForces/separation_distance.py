import numpy as np
from pyNastran.op2.result_objects.scalar1_table_object import RealScalarTableArray

class SeparationDistanceArray(RealScalarTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealScalarTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        '                    I N I T I A L  C O N T A C T  S E P A R A T I O N  D I S T A N C E'
        ' '
        '      POINT ID.   TYPE       DISTANCE'
        '             1      G      1.000000E-01'

        if header is None:
            header = []
        #' \n',
        if self.table_name == 'OSPDSI1': # initial
            words = ['                    I N I T I A L  C O N T A C T  S E P A R A T I O N  D I S T A N C E\n']
        elif self.table_name == 'OSPDS1': # deformed
            words = ['                    D E F O R M E D  C O N T A C T  S E P A R A T I O N  D I S T A N C E\n']
        else:
            raise NotImplementedError(self.table_name)
        words += [
            ' \n'
            '      POINT ID.   TYPE       DISTANCE\n'
        ]
        #words += self.get_table_marker()
        write_words = True
        if self.nonlinear_factor not in (None, np.nan):
            return self._write_f06_transient_block(
                words, header, page_stamp, page_num, f06_file, write_words,
                is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        return self._write_f06_block(
            words, header, page_stamp, page_num, f06_file, write_words,
            is_mag_phase=is_mag_phase, is_sort1=is_sort1)
