from six import iteritems
from pyNastran.op2.resultObjects.tableObject import RealTableArray, ComplexTableArray, RealTableObject, ComplexTableObject
from pyNastran.f06.f06_formatting import write_floats_13e


class RealForceVectorArray(RealTableArray):  # table_code=2, sort_code=0, thermal=0

    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                                         N O N - L I N E A R - F O R C E   V E C T O R\n']
        #words += self.get_table_marker()
        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(words, header, page_stamp, page_num, f, write_words=True,
                                                   is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        return self._write_f06_block(words, header, page_stamp, page_num, f, write_words=True)


class ComplexForceVectorArray(ComplexTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                                       C O M P L E X   F O R C E   V E C T O R\n',]
        return self._write_f06_transient_block(words, header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)


class RealForceVector(RealTableObject):  # table_code=12, sort_code=0, thermal=0

    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_sort1)

        msg = header + ['                                         N O N - L I N E A R - F O R C E   V E C T O R\n'
                        ' \n',
                        '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        for nodeID, translation in sorted(iteritems(self.translations)):
            rotation = self.rotations[nodeID]
            grid_type = self.gridTypes[nodeID]

            (dx, dy, dz) = translation
            (rx, ry, rz) = rotation
            vals = [dx, dy, dz, rx, ry, rz]
            vals2 = write_floats_13e(vals)
            if not is_all_zeros:
                [dx, dy, dz, rx, ry, rz] = vals2
                msg.append('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (nodeID, grid_type, dx, dy, dz, rx, ry, rz))

        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        msg = []
        words = ['                                         N O N - L I N E A R - F O R C E   V E C T O R\n'
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']

        for dt, translations in sorted(iteritems(self.translations)):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for nodeID, translation in sorted(iteritems(translations)):
                rotation = self.rotations[dt][nodeID]
                grid_type = self.gridTypes[nodeID]

                (dx, dy, dz) = translation
                (rx, ry, rz) = rotation

                vals = [dx, dy, dz, rx, ry, rz]
                vals2 = write_floats_13e(vals)
                if not is_all_zeros:
                    [dx, dy, dz, rx, ry, rz] = vals2
                    msg.append('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (nodeID, grid_type, dx, dy, dz, rx, ry, rz))

            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            msg = ['']
            page_num += 1
        return page_num - 1


class ComplexForceVector(ComplexTableObject):  # table_code=12, approach_code=???
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=False, is_sort1=is_sort1)
        msg = header + ['                                       C O M P L E X   F O R C E   V E C T O R\n',
                        '                                                          (REAL/IMAGINARY)\n',
                        ' \n',
                        '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        for nodeID, translation in sorted(iteritems(self.translations)):
            rotation = self.rotations[nodeID]
            grid_type = self.gridTypes[nodeID]

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
            vals2 = write_floats_13e(vals)
            [dxr, dyr, dzr, rxr, ryr, rzr,
             dxi, dyi, dzi, rxi, ryi, rzi] = vals2
            msg.append('0 %12i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (nodeID, grid_type, dxr, dyr, dzr, rxr, ryr, rzr))
            msg.append('  %12s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % ('', '',            dxi, dyi, dzi, rxi, ryi, rzi))

        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                                       C O M P L E X   F O R C E   V E C T O R\n',
                 '                                                          (REAL/IMAGINARY)\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        msg = []
        for dt, translations in sorted(iteritems(self.translations)):
            header[2] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for nodeID, translation in sorted(iteritems(translations)):
                rotation = self.rotations[dt][nodeID]
                grid_type = self.gridTypes[nodeID]

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
                vals2 = write_floats_13e(vals)
                [dxr, dyr, dzr, rxr, ryr, rzr,
                 dxi, dyi, dzi, rxi, ryi, rzi] = vals2
                msg.append('0 %12i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (nodeID, grid_type, dxr, dyr, dzr, rxr, ryr, rzr))
                msg.append('  %12s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % ('', '',           dxi, dyi, dzi, rxi, ryi, rzi))
            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            msg = ['']
            page_num += 1
        return page_num - 1