from pyNastran.op2.resultObjects.tableObject import RealTableVector, ComplexTableVector, RealTableObject, ComplexTableObject
from pyNastran.f06.f06_formatting import writeFloats13E


class RealLoadVectorVector(RealTableVector):  # table_code=2, sort_code=0, thermal=0

    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        RealTableVector.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        words = ['                                                     L O A D   V E C T O R\n', ]
        #words += self.get_table_marker()
        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(words, header, pageStamp, pageNum, f)
        return self._write_f06_block(words, header, pageStamp, pageNum, f)


class ComplexLoadVectorVector(ComplexTableVector):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableVector.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        words = ['                                               C O M P L E X   L O A D   V E C T O R\n', ]
        return self._write_f06_transient_block(words, header, pageStamp, pageNum, f, is_mag_phase)


class RealLoadVectorObject(RealTableObject):  # table_code=2, sort_code=0, thermal=0

    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        RealTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f)
        msg = header + ['                                                     L O A D   V E C T O R\n',
                        ' \n',
                        '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        f.write(''.join(msg))
        for nodeID, translation in sorted(self.translations.iteritems()):
            rotation = self.rotations[nodeID]
            grid_type = self.gridTypes[nodeID]

            (dx, dy, dz) = translation
            (rx, ry, rz) = rotation
            vals = [dx, dy, dz, rx, ry, rz]
            (vals2, is_all_zeros) = writeFloats13E(vals)
            #if not is_all_zeros:
            [dx, dy, dz, rx, ry, rz] = vals2
            f.write('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (nodeID, grid_type, dx, dy, dz, rx, ry, rz))
        f.write(pageStamp % pageNum)
        return pageNum

    def _write_f06_transient(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        words = ['                                                     L O A D   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']

        for dt, translations in sorted(self.translations.iteritems()):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            f.write(''.join(header + words))
            for nodeID, translation in sorted(translations.iteritems()):
                rotation = self.rotations[dt][nodeID]
                grid_type = self.gridTypes[nodeID]

                (dx, dy, dz) = translation
                (rx, ry, rz) = rotation

                vals = [dx, dy, dz, rx, ry, rz]
                (vals2, is_all_zeros) = writeFloats13E(vals)
                #if not is_all_zeros:
                [dx, dy, dz, rx, ry, rz] = vals2
                f.write('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (nodeID, grid_type, dx, dy, dz, rx, ry, rz))
            f.write(pageStamp % pageNum)
            pageNum += 1
        return pageNum - 1


class ComplexLoadVectorObject(ComplexTableObject):  # table_code=11, approach_code=???
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f, is_mag_phase)
        msg = header + ['                                               C O M P L E X   L O A D   V E C T O R\n',
                        '                                                          (REAL/IMAGINARY)\n',
                        ' \n',
                        '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        f.write(''.join(msg))
        for nodeID, translation in sorted(self.translations.iteritems()):
            rotation = self.rotations[nodeID]
            grid_type = self.gridTypes[nodeID]

            (dx, dy, dz) = translation
            (rx, ry, rz) = rotation

            if is_mag_phase:
                dxr = abs(dx)
                dxi = angle(dx, deg=True)
                dyr = abs(dy)
                dyi = angle(dy, deg=True)
                dzr = abs(dz)
                dzi = angle(dz, deg=True)

                rxr = abs(rx)
                rxi = angle(rx, deg=True)
                ryr = abs(ry)
                ryi = angle(ry, deg=True)
                rzr = abs(rz)
                rzi = angle(rz, deg=True)
            else:
                dxr = dx.real
                dyr = dy.real
                dzr = dz.real
                dxi = dx.imag
                dyi = dy.imag
                dzi = dz.imag

                rxr = rx.real
                ryr = ry.real
                rzr = rz.real
                rxi = rx.imag
                ryi = ry.imag
                rzi = rz.imag

            vals = [dxr, dyr, dzr, rxr, ryr, rzr, dxi, dyi, dzi, rxi, ryi, rzi]
            (vals2, is_all_zeros) = writeFloats13E(vals)
            [dxr, dyr, dzr, rxr, ryr, rzr,
             dxi, dyi, dzi, rxi, ryi, rzi] = vals2
            f.write('0 %12i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                    '  %12s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                      nodeID, grid_type, dxr, dyr, dzr, rxr, ryr, rzr,
                      '', '',            dxi, dyi, dzi, rxi, ryi, rzi))
        f.write(pageStamp % pageNum)
        return pageNum

    def _write_f06_transient(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        words = ['                                               C O M P L E X   L O A D   V E C T O R\n',
                 '                                                          (REAL/IMAGINARY)\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #return self.__write_f06_transient_block(words,header,pageStamp,pageNum,f,is_mag_phase)

        for dt, translations in sorted(self.translations.iteritems()):
            header[2] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            f.write(''.join(header + words))
            for nodeID, translation in sorted(translations.iteritems()):
                rotation = self.rotations[dt][nodeID]
                grid_type = self.gridTypes[nodeID]

                (dx, dy, dz) = translation
                (rx, ry, rz) = rotation

                if is_mag_phase:
                    dxr = abs(dx)
                    dxi = angle(dx, deg=True)
                    dyr = abs(dy)
                    dyi = angle(dy, deg=True)
                    dzr = abs(dz)
                    dzi = angle(dz, deg=True)

                    rxr = abs(rx)
                    rxi = angle(rx, deg=True)
                    ryr = abs(ry)
                    ryi = angle(ry, deg=True)
                    rzr = abs(rz)
                    rzi = angle(rz, deg=True)
                else:
                    dxr = dx.real
                    dyr = dy.real
                    dzr = dz.real
                    dxi = dx.imag
                    dyi = dy.imag
                    dzi = dz.imag

                    rxr = rx.real
                    ryr = ry.real
                    rzr = rz.real
                    rxi = rx.imag
                    ryi = ry.imag
                    rzi = rz.imag

                vals = [dxr, dyr, dzr, rxr, ryr, rzr,
                        dxi, dyi, dzi, rxi, ryi, rzi]
                (vals2, is_all_zeros) = writeFloats13E(vals)
                [dxr, dyr, dzr, rxr, ryr, rzr,
                 dxi, dyi, dzi, rxi, ryi, rzi] = vals2
                #if not is_all_zeros:
                f.write('0 %12i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                        '  %12s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                            nodeID, grid_type, dxr, dyr, dzr, rxr, ryr, rzr,
                            '', '',            dxi, dyi, dzi, rxi, ryi, rzi))

            f.write(pageStamp % pageNum)
            pageNum += 1
        return pageNum - 1


class RealThermalVectorObject(RealTableObject):
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        RealTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f)
        msg = header + ['                                              T E M P E R A T U R E   V E C T O R\n',
                        ' \n',
                        '      POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE\n']
        f.write(''.join(msg))
        for nodeID, translation in sorted(self.translations.iteritems()):
            rotation = self.rotations[nodeID]
            grid_type = self.gridTypes[nodeID]

            (dx, dy, dz) = translation
            (rx, ry, rz) = rotation
            vals = [dx, dy, dz, rx, ry, rz]
            (vals2, is_all_zeros) = writeFloats13E(vals)
            #if not is_all_zeros:
            [dx, dy, dz, rx, ry, rz] = vals2
            f.write('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (nodeID, grid_type, dx, dy, dz, rx, ry, rz))

        f.write(pageStamp % pageNum)
        return pageNum

    def _write_f06_transient(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        words = ['                                              T E M P E R A T U R E   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE\n']

        for dt, translations in sorted(self.translations.iteritems()):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            f.write(''.join(header + words))
            for nodeID, translation in sorted(translations.iteritems()):
                rotation = self.rotations[dt][nodeID]
                grid_type = self.gridTypes[nodeID]

                (dx, dy, dz) = translation
                (rx, ry, rz) = rotation

                vals = [dx, dy, dz, rx, ry, rz]
                (vals2, is_all_zeros) = writeFloats13E(vals)
                if not is_all_zeros:
                    [dx, dy, dz, rx, ry, rz] = vals2
                    f.write('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (nodeID, grid_type, dx, dy, dz, rx, ry, rz))

            f.write(pageStamp % pageNum)
            pageNum += 1
        return pageNum - 1


class ThermalLoadVectorObject(RealThermalVectorObject):     # table_code=2, thermal=1
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        RealThermalVectorObject.__init__(self, data_code, is_sort1, isubcase, dt)


class ThermalVelocityVectorObject(RealThermalVectorObject):  # table_code=10, thermal=1
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        RealThermalVectorObject.__init__(self, data_code, is_sort1, isubcase, dt)
