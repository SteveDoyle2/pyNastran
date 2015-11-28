from six import iteritems
from pyNastran.op2.resultObjects.tableObject import RealTableArray, ComplexTableArray, RealTableObject, ComplexTableObject
from pyNastran.f06.f06_formatting import writeFloats13E


class RealLoadVectorArray(RealTableArray):  # table_code=2, sort_code=0, thermal=0

    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                                                     L O A D   V E C T O R\n', ]
        #words += self.get_table_marker()
        write_words = True
        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(words, header, page_stamp, page_num, f, write_words,
                                                   is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        return self._write_f06_block(words, header, page_stamp, page_num, f, write_words,
                         is_mag_phase=False, is_sort1=True)


class ComplexLoadVectorArray(ComplexTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                                               C O M P L E X   L O A D   V E C T O R\n', ]
        return self._write_f06_transient_block(words, header, page_stamp, page_num, f,
                                               is_mag_phase, is_sort1)


class RealLoadVector(RealTableObject):  # table_code=2, sort_code=0, thermal=0

    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_sort1)
        msg = header + ['                                                     L O A D   V E C T O R\n',
                        ' \n',
                        '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        f.write(''.join(msg))
        for nodeID, translation in sorted(iteritems(self.translations)):
            rotation = self.rotations[nodeID]
            grid_type = self.gridTypes[nodeID]

            (dx, dy, dz) = translation
            (rx, ry, rz) = rotation
            vals = [dx, dy, dz, rx, ry, rz]
            (vals2, is_all_zeros) = writeFloats13E(vals)
            #if not is_all_zeros:
            [dx, dy, dz, rx, ry, rz] = vals2
            f.write('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (nodeID, grid_type, dx, dy, dz, rx, ry, rz))
        f.write(page_stamp % page_num)
        return page_num

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                                                     L O A D   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']

        for dt, translations in sorted(iteritems(self.translations)):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            f.write(''.join(header + words))
            for nodeID, translation in sorted(iteritems(translations)):
                rotation = self.rotations[dt][nodeID]
                grid_type = self.gridTypes[nodeID]

                (dx, dy, dz) = translation
                (rx, ry, rz) = rotation

                vals = [dx, dy, dz, rx, ry, rz]
                (vals2, is_all_zeros) = writeFloats13E(vals)
                #if not is_all_zeros:
                [dx, dy, dz, rx, ry, rz] = vals2
                f.write('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (nodeID, grid_type, dx, dy, dz, rx, ry, rz))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class ComplexLoadVector(ComplexTableObject):  # table_code=11, approach_code=???
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        msg = header + ['                                               C O M P L E X   L O A D   V E C T O R\n',
                        '                                                          (REAL/IMAGINARY)\n',
                        ' \n',
                        '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        f.write(''.join(msg))
        for nodeID, translation in sorted(iteritems(self.translations)):
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
        f.write(page_stamp % page_num)
        return page_num

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None,
                             is_mag_phase=False, is_sort1=True):
        words = ['                                               C O M P L E X   L O A D   V E C T O R\n',
                 '                                                          (REAL/IMAGINARY)\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #return self.__write_f06_transient_block(words,header,page_stamp,page_num,f,is_mag_phase)

        for dt, translations in sorted(iteritems(self.translations)):
            header[2] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            f.write(''.join(header + words))
            for nodeID, translation in sorted(iteritems(translations)):
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

            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealTemperatureVectorArray(RealTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = [
            '                                              T E M P E R A T U R E   V E C T O R\n',
            ' \n',
            '      POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE\n'
        ]
        #words += self.get_table_marker()
        write_words = False
        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(words, header, page_stamp, page_num, f, write_words,
                                                   is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        return self._write_f06_block(words, header, page_stamp, page_num, f, write_words,
                                         is_mag_phase=is_mag_phase, is_sort1=is_sort1)

class RealThermalVector(RealTableObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None,
                  is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f,
                                             is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        msg = header + ['                                              T E M P E R A T U R E   V E C T O R\n',
                        ' \n',
                        '      POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE\n']
        f.write(''.join(msg))
        for nodeID, translation in sorted(iteritems(self.translations)):
            rotation = self.rotations[nodeID]
            grid_type = self.gridTypes[nodeID]

            (dx, dy, dz) = translation
            (rx, ry, rz) = rotation
            vals = [dx, dy, dz, rx, ry, rz]
            (vals2, is_all_zeros) = writeFloats13E(vals)
            #if not is_all_zeros:
            [dx, dy, dz, rx, ry, rz] = vals2
            f.write('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (nodeID, grid_type, dx, dy, dz, rx, ry, rz))

        f.write(page_stamp % page_num)
        return page_num

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None,
                             is_mag_phase=False, is_sort1=True):
        words = ['                                              T E M P E R A T U R E   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE\n']

        for dt, translations in sorted(iteritems(self.translations)):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            f.write(''.join(header + words))
            for nodeID, translation in sorted(iteritems(translations)):
                rotation = self.rotations[dt][nodeID]
                grid_type = self.gridTypes[nodeID]

                (dx, dy, dz) = translation
                (rx, ry, rz) = rotation

                vals = [dx, dy, dz, rx, ry, rz]
                (vals2, is_all_zeros) = writeFloats13E(vals)
                if not is_all_zeros:
                    [dx, dy, dz, rx, ry, rz] = vals2
                    f.write('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (nodeID, grid_type, dx, dy, dz, rx, ry, rz))

            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


#class RealThermalLoadVector(RealThermalVector):     # table_code=2, thermal=1
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #RealThermalVector.__init__(self, data_code, is_sort1, isubcase, dt)


#class RealThermalVelocityVector(RealThermalVector):  # table_code=10, thermal=1
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #RealThermalVector.__init__(self, data_code, is_sort1, isubcase, dt)

class RealThermalVelocityVectorArray(RealTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = [
            '                                              THERMAL VELOCITY   V E C T O R\n',
            ' \n',
            '      POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE\n'
        ]
        #words += self.get_table_marker()
        write_words = False
        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(words, header, page_stamp, page_num, f, write_words,
                                                   is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        return self._write_f06_block(words, header, page_stamp, page_num, f, write_words,
                                         is_mag_phase=is_mag_phase, is_sort1=is_sort1)

