from pyNastran.op2.resultObjects.tableObject import TableObject, ComplexTableObject
from pyNastran.f06.f06_formatting import writeFloats13E


class ForceVectorObject(TableObject):  # table_code=12, sort_code=0, thermal=0

    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        TableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_matlab(self, isubcase, f=None, is_mag_phase=False):
        name = 'forceVector'
        if self.nonlinear_factor is None:
            return self._write_matlab(name, isubcase, f)
        else:
            return self._write_matlab_transient(name, isubcase, f)

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f)

        msg = header + ['                                         N O N - L I N E A R - F O R C E   V E C T O R\n'
                        ' \n',
                        '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        for nodeID, translation in sorted(self.translations.iteritems()):
            rotation = self.rotations[nodeID]
            gridType = self.gridTypes[nodeID]

            (dx, dy, dz) = translation
            (rx, ry, rz) = rotation
            vals = [dx, dy, dz, rx, ry, rz]
            (vals2, isAllZeros) = writeFloats13E(vals)
            if not isAllZeros:
                [dx, dy, dz, rx, ry, rz] = vals2
                msg.append('%14i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' % (nodeID, gridType, dx, dy, dz, rx, ry, rz.rstrip()))

        msg.append(pageStamp % pageNum)
        f.write(''.join(msg))
        return pageNum

    def _write_f06_transient(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        msg = []
        words = ['                                         N O N - L I N E A R - F O R C E   V E C T O R\n'
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']

        for dt, translations in sorted(self.translations.iteritems()):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for nodeID, translation in sorted(translations.iteritems()):
                rotation = self.rotations[dt][nodeID]
                gridType = self.gridTypes[nodeID]

                (dx, dy, dz) = translation
                (rx, ry, rz) = rotation

                vals = [dx, dy, dz, rx, ry, rz]
                (vals2, isAllZeros) = writeFloats13E(vals)
                if not isAllZeros:
                    [dx, dy, dz, rx, ry, rz] = vals2
                    msg.append('%14i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' % (nodeID, gridType, dx, dy, dz, rx, ry, rz.rstrip()))

            msg.append(pageStamp % pageNum)
            f.write(''.join(msg))
            msg = ['']
            pageNum += 1
        return pageNum - 1


class ComplexForceVectorObject(ComplexTableObject):  # table_code=12, approach_code=???
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_matlab(self, isubcase, f=None, is_mag_phase=False):
        name = 'forceVector'
        if self.nonlinear_factor is None:
            return self._write_matlab(name, isubcase, f, is_mag_phase)
        else:
            return self._write_matlab_transient(name, isubcase, f, is_mag_phase)

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f, is_mag_phase=False)
        msg = header + ['                                       C O M P L E X   FORCE   V E C T O R\n',
                        '                                                          (REAL/IMAGINARY)\n',
                        ' \n',
                        '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        for nodeID, translation in sorted(self.translations.iteritems()):
            rotation = self.rotations[nodeID]
            gridType = self.gridTypes[nodeID]

            (dx, dy, dz) = translation
            dxr = dx.real
            dyr = dy.real
            dzr = dz.real
            dxi = dx.imag
            dyi = dy.imag
            dzi = dz.imag

            (rx, ry, rz) = rotation
            rxr = rx.real
            ryr = ry.real
            rzr = rz.real
            rxi = rx.imag
            ryi = ry.imag
            rzi = rz.imag

            vals = [dxr, dyr, dzr, rxr, ryr, rzr, dxi, dyi, dzi, rxi, ryi, rzi]
            (vals2, isAllZeros) = writeFloats13E(vals)
            [dxr, dyr, dzr, rxr, ryr, rzr,
             dxi, dyi, dzi, rxi, ryi, rzi] = vals2
            msg.append('0 %12i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' % (nodeID, gridType, dxr, dyr, dzr, rxr, ryr, rzr.rstrip()))
            msg.append('  %12s %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' % ('', '',           dxi, dyi, dzi, rxi, ryi, rzi.rstrip()))

        msg.append(pageStamp % pageNum)
        f.write(''.join(msg))
        return pageNum

    def _write_f06_transient(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        words = ['                                       C O M P L E X   F O R C E   V E C T O R\n',
                 '                                                          (REAL/IMAGINARY)\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        msg = []
        for dt, translations in sorted(self.translations.iteritems()):
            header[2] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for nodeID, translation in sorted(translations.iteritems()):
                rotation = self.rotations[dt][nodeID]
                gridType = self.gridTypes[nodeID]

                (dx, dy, dz) = translation
                dxr = dx.real
                dyr = dy.real
                dzr = dz.real
                dxi = dx.imag
                dyi = dy.imag
                dzi = dz.imag

                (rx, ry, rz) = rotation
                rxr = rx.real
                ryr = ry.real
                rzr = rz.real
                rxi = rx.imag
                ryi = ry.imag
                rzi = rz.imag

                vals = [dxr, dyr, dzr, rxr, ryr, rzr, dxi,
                        dyi, dzi, rxi, ryi, rzi]
                (vals2, isAllZeros) = writeFloats13E(vals)
                [dxr, dyr, dzr, rxr, ryr, rzr,
                 dxi, dyi, dzi, rxi, ryi, rzi] = vals2
                msg.append('0 %12i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' % (nodeID, gridType, dxr, dyr, dzr, rxr, ryr, rzr.rstrip()))
                msg.append('  %12s %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' % ('', '',           dxi, dyi, dzi, rxi, ryi, rzi.rstrip()))
            msg.append(pageStamp % pageNum)
            f.write(''.join(msg))
            msg = ['']
            pageNum += 1
        return pageNum - 1