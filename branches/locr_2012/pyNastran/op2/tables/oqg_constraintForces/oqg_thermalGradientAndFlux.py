from pyNastran.op2.resultObjects.tableObject import TableObject, ComplexTableObject, writeFloats13E


class TemperatureGradientAndFluxObject(TableObject):
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        TableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_matlab(self, isubcase, f=None, isMagPhase=False):
        name = 'spcForces'
        if self.nonlinear_factor is None:
            return self._write_matlab(name, isubcase, f)
        else:
            return self._write_matlab_transient(name, isubcase, f)

    def write_f06(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f)
        msg = header + ['                   F I N I T E   E L E M E N T   T E M P E R A T U R E   G R A D I E N T S   A N D   F L U X E S\n',
                        ' \n',
                        '   ELEMENT-ID   EL-TYPE        X-GRADIENT       Y-GRADIENT       Z-GRADIENT        X-FLUX           Y-FLUX           Z-FLUX\n']
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
        if f is not None:
            f.write(''.join(msg))
            msg = ['']
        return (''.join(msg), pageNum)

    def _write_f06_transient(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        words = ['                   F I N I T E   E L E M E N T   T E M P E R A T U R E   G R A D I E N T S   A N D   F L U X E S\n',
                 ' \n',
                 '   ELEMENT-ID   EL-TYPE        X-GRADIENT       Y-GRADIENT       Z-GRADIENT        X-FLUX           Y-FLUX           Z-FLUX\n']
        return self._write_f06_transient_block(words, header, pageStamp, pageNum, f)

    def __reprTransient__(self):
        msg = '---SPC FORCES---\n'
        if self.nonlinear_factor is not None:
            msg += 'dt = %g\n' % (self.dt)

        headers = ['T1', 'T2', 'T3', 'R1', 'R2', 'R3']
        msg += '%-8s ' % ('GRID')
        for header in headers:
            msg += '%10s ' % (header)
        msg += '\n'

        for dt, translations in sorted(self.translations.iteritems()):
            msg += 'dt = %s' % (dt)
            for nodeID, translation in sorted(translations.iteritems()):
                rotation = self.rotations[dt][nodeID]
                (Fx, Fy, Fz) = translation
                (Mx, My, Mz) = rotation

                msg += '%-8i ' % (nodeID)
                vals = [Fx, Fy, Fz, Mx, My, Mx]
                for val in vals:
                    if abs(val) < 1e-6:
                        msg += '%10s ' % (0)
                    else:
                        msg += '%10.2f ' % (val)
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
        msg += '%-8s ' % ('GRID')
        for header in headers:
            msg += '%10s ' % (header)
        msg += '\n'

        for nodeID, translation in sorted(self.translations.iteritems()):
            rotation = self.rotations[nodeID]
            (Fx, Fy, Fz) = translation
            (Mx, My, Mz) = rotation

            msg += '%-8i ' % (nodeID)
            vals = [Fx, Fy, Fz, Mx, My, Mx]
            for val in vals:
                if abs(val) < 1e-6:
                    msg += '%10s ' % (0)
                else:
                    msg += '%10.2f ' % (val)
            msg += '\n'
        return msg


class ComplexTemperatureGradientAndFluxObject(ComplexTableObject):
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        asdf
        ComplexTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_matlab(self, isubcase, f=None, isMagPhase=False):
        name = 'spcForces'
        if self.nonlinear_factor is None:
            return self._write_matlab(name, isubcase, f)
        else:
            return self._write_matlab_transient(name, isubcase, f)

    def write_f06(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f, isMagPhase)
        msg = header + ['                               F O R C E S   O F   S I N G L E - P O I N T   C O N S T R A I N T\n',
                        ' \n',
                        '   ELEMENT-ID   EL-TYPE        X-GRADIENT       Y-GRADIENT       Z-GRADIENT        X-FLUX           Y-FLUX           Z-FLUX\n']
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
        if f is not None:
            f.write(''.join(msg))
            msg = ['']
        return (''.join(msg), pageNum)

    def _write_f06_transient(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        words = ['                         C O M P L E X   F O R C E S   O F   S I N G L E   P O I N T   C O N S T R A I N T\n']
        return self._write_f06_transient_block(words, header, pageStamp, pageNum, f, isMagPhase)

    def __reprTransient__(self):
        msg = '---COMPLEX SPC FORCES---\n'
        if self.nonlinear_factor is not None:
            msg += 'dt = %g\n' % (self.dt)

        raise RuntimeError('is this valid...')
        headers = ['T1', 'T2', 'T3', 'R1', 'R2', 'R3']
        msg += '%-8s ' % ('GRID')
        for header in headers:
            msg += '%10s ' % (header)
        msg += '\n'

        for dt, translations in sorted(self.translations.iteritems()):
            msg += 'dt = %s' % (dt)
            for nodeID, translation in sorted(translations.iteritems()):
                rotation = self.rotations[dt][nodeID]
                msg += '%-8i ' % (nodeID)
                vals = translation + rotation
                for val in vals:
                    if abs(val) < 1e-6:
                        msg += '%10s ' % (0)
                    else:
                        msg += '%10.2f ' % (val)
                msg += '\n'
        return msg

    def __repr__(self):
        return self.write_f06(['', ''], 'PAGE ', 1)[0]
        if self.nonlinear_factor is not None:
            return self.__reprTransient__()

        msg = '---COMPLEX SPC FORCES---\n'
        if self.nonlinear_factor is not None:
            msg += 'dt = %g\n' % (self.dt)

        raise RuntimeError('is this valid...')
        headers = ['T1', 'T2', 'T3', 'R1', 'R2', 'R3']
        msg += '%-8s ' % ('GRID')
        for header in headers:
            msg += '%10s ' % (header)
        msg += '\n'

        for nodeID, translation in sorted(self.translations.iteritems()):
            rotation = self.rotations[nodeID]
            (Fx, Fy, Fz) = translation
            (Mx, My, Mz) = rotation

            msg += '%-8i ' % (nodeID)
            vals = [Fx, Fy, Fz, Mx, My, Mx]
            for val in vals:
                if abs(val) < 1e-6:
                    msg += '%10s ' % (0)
                else:
                    msg += '%10.2f ' % (val)
            msg += '\n'
        return msg
