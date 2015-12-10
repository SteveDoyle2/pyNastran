from six import iteritems

from pyNastran.f06.f06_formatting import write_floats_13e
from pyNastran.op2.resultObjects.tableObject import RealTableArray, ComplexTableArray, RealTableObject, ComplexTableObject

class RealTemperatureGradientAndFluxArray(RealTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        msg = header + ['                   F I N I T E   E L E M E N T   T E M P E R A T U R E   G R A D I E N T S   A N D   F L U X E S\n',
                        ' \n',
                        '   ELEMENT-ID   EL-TYPE        X-GRADIENT       Y-GRADIENT       Z-GRADIENT        X-FLUX           Y-FLUX           Z-FLUX\n']
        #words += self.get_table_marker()
        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(words, header, page_stamp, page_num, f)
        return self._write_f06_block(words, header, page_stamp, page_num, f)


class RealTemperatureGradientAndFlux(RealTableObject):

    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        msg = header + ['                   F I N I T E   E L E M E N T   T E M P E R A T U R E   G R A D I E N T S   A N D   F L U X E S\n',
                        ' \n',
                        '   ELEMENT-ID   EL-TYPE        X-GRADIENT       Y-GRADIENT       Z-GRADIENT        X-FLUX           Y-FLUX           Z-FLUX\n']
        f.write(''.join(msg))
        for nodeID, translation in sorted(iteritems(self.translations)):
            rotation = self.rotations[nodeID]
            grid_type = self.gridTypes[nodeID]

            (dx, dy, dz) = translation
            (rx, ry, rz) = rotation
            vals = [dx, dy, dz, rx, ry, rz]
            vals2 = write_floats_13e(vals)
            #if not is_all_zeros:
            [dx, dy, dz, rx, ry, rz] = vals2
            f.write('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (nodeID, grid_type, dx, dy, dz, rx, ry, rz))
        f.write(page_stamp % page_num)
        return page_num

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                   F I N I T E   E L E M E N T   T E M P E R A T U R E   G R A D I E N T S   A N D   F L U X E S\n',
                 ' \n',
                 '   ELEMENT-ID   EL-TYPE        X-GRADIENT       Y-GRADIENT       Z-GRADIENT        X-FLUX           Y-FLUX           Z-FLUX\n']
        return self._write_f06_transient_block(words, header, page_stamp, page_num, f)


class ComplexTemperatureGradientAndFlux(ComplexTableObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        raise NotImplementedError()
        ComplexTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        msg = header + ['                               F O R C E S   O F   S I N G L E - P O I N T   C O N S T R A I N T\n',
                        ' \n',
                        '   ELEMENT-ID   EL-TYPE        X-GRADIENT       Y-GRADIENT       Z-GRADIENT        X-FLUX           Y-FLUX           Z-FLUX\n']
        raise RuntimeError('is this valid...')
        for nodeID, translation in sorted(iteritems(self.translations)):
            rotation = self.rotations[nodeID]
            grid_type = self.gridTypes[nodeID]

            (dx, dy, dz) = translation
            #dxr=dx.real; dyr=dy.real; dzr=dz.real;
            #dxi=dx.imag; dyi=dy.imag; dzi=dz.imag

            (rx, ry, rz) = rotation
            #rxr=rx.real; ryr=ry.real; rzr=rz.real
            #rxi=rx.imag; ryi=ry.imag; rzi=rz.imag

            #vals = [dxr,dyr,dzr,rxr,ryr,rzr,dxi,dyi,dzi,rxi,ryi,rzi]
            vals = list(translation) + list(rotation)
            vals2 = write_floats_13e(vals)
            #if not is_all_zeros:
            [dx, dy, dz, rx, ry, rz] = vals2
            msg.append('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (nodeID, grid_type, dx, dy, dz, rx, ry, rz))
        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                         C O M P L E X   F O R C E S   O F   S I N G L E   P O I N T   C O N S T R A I N T\n']
        return self._write_f06_transient_block(words, header, page_stamp, page_num, f, is_mag_phase)

