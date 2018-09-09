import numpy as np
from pyNastran.op2.result_objects.table_object import RealTableArray

class RealTemperatureGradientAndFluxArray(RealTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        words = header + [
            '                   F I N I T E   E L E M E N T   T E M P E R A T U R E   G R A D I E N T S   A N D   F L U X E S\n',
            ' \n',
            '   ELEMENT-ID   EL-TYPE        X-GRADIENT       Y-GRADIENT       Z-GRADIENT        X-FLUX           Y-FLUX           Z-FLUX\n']
        if self.table_name in ['OQG1', 'OQG2']:
            pass
        elif self.table_name in ['OQGATO1', 'OQGATO2']:
            words += ['                                                 ( AUTO-CORRELATION FUNCTION )']
        elif self.table_name in ['OQGPSD1', 'OQGPSD2']:
            words += ['                                             ( POWER SPECTRAL DENSITY FUNCTION )']
        elif self.table_name in ['OQGRMS1', 'OQGRMS2']:
            words += ['                                                     ( ROOT MEAN SQUARE )']
        elif self.table_name in ['OQGCRM1', 'OQGCRM2']:
            words += ['                                               ( CUMULATIVE ROOT MEAN SQUARE )']
        elif self.table_name in ['OQGNO1', 'OQGNO2']:
            words += ['                                                 ( NUMBER OF ZERO CROSSINGS )']
        else:
            raise NotImplementedError('table_name=%r' % self.table_name)
        #words += self.get_table_marker()
        write_words = False
        if self.nonlinear_factor not in (None, np.nan):
            return self._write_f06_transient_block(words, header, page_stamp, page_num, f06_file, write_words)
        return self._write_f06_block(words, header, page_stamp, page_num, f06_file, write_words)
