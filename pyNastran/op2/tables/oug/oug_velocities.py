import numpy as np
from pyNastran.op2.result_objects.table_object import RealTableArray, ComplexTableArray


class RealVelocityArray(RealTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        words = ['                                                   V E L O C I T Y   V E C T O R\n', ]
        #' \n',
        #'      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        if self.table_name in ['OUGV1', 'OUGV2', 'OUPV1', 'BOUGV1']:
            pass
        elif self.table_name in ['OUXY1', 'OUXY2']:
            words = ['                                           V E L O C I T Y   V E C T O R   (SOLUTION SET)']
        elif self.table_name in ['ROUGV1', 'ROUGV2']:
            words += ['                                                (RELATIVE TO ENFORCED MOTION INPUT)']
        elif self.table_name in ['OVGATO1', 'OVGATO2']:
            words += ['                                                 ( AUTO-CORRELATION FUNCTION )']
        elif self.table_name in ['OUGPSD1', 'OUGPSD2', 'OVGPSD1', 'OVGPSD2']:
            words += ['                                             ( POWER SPECTRAL DENSITY FUNCTION )']
        elif self.table_name in ['OUGRMS1', 'OUGRMS2', 'OVGRMS1', 'OVGRMS2']:
            words += ['                                                     ( ROOT MEAN SQUARE )']
        elif self.table_name in ['OUGCRM1', 'OUGCRM2', 'OVGCRM1', 'OVGCRM2']:
            words += ['                                               ( CUMULATIVE ROOT MEAN SQUARE )']
        elif self.table_name in ['OUGNO1', 'OUGNO2', 'OVGNO1', 'OVGNO2']:  # , 'OVGNO1', 'OVGNO2'
            words += ['                                                 ( NUMBER OF ZERO CROSSINGS )']
        else:
            raise NotImplementedError(self.table_name)

        write_words = True
        if self.nonlinear_factor not in (None, np.nan):
            return self._write_f06_transient_block(words, header, page_stamp, page_num, f06_file, write_words,
                                                   is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        #words += self.get_table_marker()
        return self._write_f06_block(words, header, page_stamp, page_num, f06_file, write_words,
                                     is_mag_phase=is_mag_phase, is_sort1=is_sort1)


class ComplexVelocityArray(ComplexTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        words = ['                                       C O M P L E X   V E L O C I T Y   V E C T O R\n']
        return self._write_f06_transient_block(words, header, page_stamp, page_num, f06_file,
                                               is_mag_phase=is_mag_phase, is_sort1=is_sort1)
