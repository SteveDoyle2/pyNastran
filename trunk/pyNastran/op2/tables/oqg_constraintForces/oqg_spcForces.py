from pyNastran.op2.resultObjects.tableObject import TableObject, ComplexTableObject
from pyNastran.f06.f06_formatting import writeFloats13E


class SPCForcesObject(TableObject):

    def __init__(self, data_code, is_sort1, isubcase, dt=None, read_mode=0):
        TableObject.__init__(self, data_code, is_sort1, isubcase, dt, read_mode)

    def write_matlab(self, isubcase, f, is_mag_phase=False):
        name = 'spcForces'
        if self.nonlinear_factor is None:
            return self._write_matlab(name, isubcase, f)
        else:
            return self._write_matlab_transient(name, isubcase, f)

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, f, pageNum)
        msg = header + ['                               F O R C E S   O F   S I N G L E - P O I N T   C O N S T R A I N T\n',
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

        msg.append(pageStamp + str(pageNum) + '\n')
        f.write(''.join(msg))
        return pageNum

    def _write_f06_transient(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        words = ['                               F O R C E S   O F   S I N G L E - P O I N T   C O N S T R A I N T\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        return self._write_f06_transient_block(words, header, pageStamp, f, pageNum)

    def __reprTransient__(self):
        msg = '---SPC FORCES---\n'
        if self.nonlinear_factor is not None:
            msg += 'dt = %g\n' % self.dt

        headers = ['T1', 'T2', 'T3', 'R1', 'R2', 'R3']
        msg += '%-8s ' % ('GRID')
        for header in headers:
            msg += '%10s ' % header
        msg += '\n'

        for dt, translations in sorted(self.translations.iteritems()):
            msg += 'dt = %s' % (dt)
            for nodeID, translation in sorted(translations.iteritems()):
                rotation = self.rotations[dt][nodeID]
                (Fx, Fy, Fz) = translation
                (Mx, My, Mz) = rotation

                msg += '%-8i ' % nodeID
                vals = [Fx, Fy, Fz, Mx, My, Mx]
                for val in vals:
                    if abs(val) < 1e-6:
                        msg += '%10s ' % 0
                    else:
                        msg += '%10.2f ' % val
                msg += '\n'
        return msg

    def __repr__(self):
        return self.write_f06(['', ''], 'PAGE ', 1)[0]
        if self.nonlinear_factor is not None:
            return self.__reprTransient__()

        msg = '---SPC FORCES---\n'
        if self.dt is not None:
            msg += 'dt = %g\n' % (self.dt)

        headers = ['T1', 'T2', 'T3', 'R1', 'R2', 'R3']
        msg += '%-8s ' % 'GRID'
        for header in headers:
            msg += '%10s ' % header
        msg += '\n'

        for nodeID, translation in sorted(self.translations.iteritems()):
            rotation = self.rotations[nodeID]
            (Fx, Fy, Fz) = translation
            (Mx, My, Mz) = rotation

            msg += '%-8i ' % (nodeID)
            vals = [Fx, Fy, Fz, Mx, My, Mx]
            for val in vals:
                if abs(val) < 1e-6:
                    msg += '%10s ' % 0
                else:
                    msg += '%10.2f ' % val
            msg += '\n'
        return msg


class ComplexSPCForcesObject(ComplexTableObject):
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        ComplexTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_matlab(self, isubcase, f, is_mag_phase=False):
        name = 'spcForces'
        if self.nonlinear_factor is None:
            return self._write_matlab(name, isubcase, f)
        else:
            return self._write_matlab_transient(name, isubcase, f)

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, f, pageNum, is_mag_phase)
        msg = header + ['                               F O R C E S   O F   S I N G L E - P O I N T   C O N S T R A I N T\n',
                        ' \n',
                        '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        raise RuntimeError('is this valid...')
        for nodeID, translation in sorted(self.translations.iteritems()):
            rotation = self.rotations[nodeID]
            gridType = self.gridTypes[nodeID]

            (dx, dy, dz) = translation
            #dxr=dx.real; dyr=dy.real; dzr=dz.real;
            #dxi=dx.imag; dyi=dy.imag; dzi=dz.imag

            (rx, ry, rz) = rotation
            #rxr=rx.real; ryr=ry.real; rzr=rz.real
            #rxi=rx.imag; ryi=ry.imag; rzi=rz.imag

            #vals = [dxr,dyr,dzr,rxr,ryr,rzr,dxi,dyi,dzi,rxi,ryi,rzi]
            vals = list(translation) + list(rotation)
            (vals2, isAllZeros) = writeFloats13E(vals)
            if not isAllZeros:
                [dx, dy, dz, rx, ry, rz] = vals2
                msg.append('%14i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' % (nodeID, gridType, dx, dy, dz, rx, ry, rz.rstrip()))
        msg.append(pageStamp + str(pageNum) + '\n')
        f.write(''.join(msg))
        return pageNum

    def _write_f06_transient(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        words = ['                         C O M P L E X   F O R C E S   O F   S I N G L E   P O I N T   C O N S T R A I N T\n']
        return self._write_f06_transient_block(words, header, pageStamp, f, pageNum, is_mag_phase)

    def __reprTransient__(self):
        msg = '---COMPLEX SPC FORCES---\n'
        if self.nonlinear_factor is not None:
            msg += 'dt = %g\n' % self.dt

        raise RuntimeError('is this valid...')
        headers = ['T1', 'T2', 'T3', 'R1', 'R2', 'R3']
        msg += '%-8s ' % ('GRID')
        for header in headers:
            msg += '%10s ' % header
        msg += '\n'

        for dt, translations in sorted(self.translations.iteritems()):
            msg += 'dt = %s' % (dt)
            for nodeID, translation in sorted(translations.iteritems()):
                rotation = self.rotations[dt][nodeID]
                msg += '%-8i ' % (nodeID)
                vals = translation + rotation
                for val in vals:
                    if abs(val) < 1e-6:
                        msg += '%10s ' % 0
                    else:
                        msg += '%10.2f ' % val
                msg += '\n'
        return msg

    def __repr__(self):
        return self.write_f06(['', ''], 'PAGE ', 1)[0]
        if self.nonlinear_factor is not None:
            return self.__reprTransient__()

        msg = '---COMPLEX SPC FORCES---\n'
        if self.nonlinear_factor is not None:
            msg += 'dt = %g\n' % self.dt

        raise RuntimeError('is this valid...')
        headers = ['T1', 'T2', 'T3', 'R1', 'R2', 'R3']
        msg += '%-8s ' % ('GRID')
        for header in headers:
            msg += '%10s ' % header
        msg += '\n'

        for nodeID, translation in sorted(self.translations.iteritems()):
            rotation = self.rotations[nodeID]
            (Fx, Fy, Fz) = translation
            (Mx, My, Mz) = rotation

            msg += '%-8i ' % nodeID
            vals = [Fx, Fy, Fz, Mx, My, Mx]
            for val in vals:
                if abs(val) < 1e-6:
                    msg += '%10s ' % 0
                else:
                    msg += '%10.2f ' % val
            msg += '\n'
        return msg
