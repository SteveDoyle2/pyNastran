from pyNastran.op2.resultObjects.tableObject import TableObject, ComplexTableObject


class VelocityObject(TableObject):  # approach_code=10, thermal=0
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        TableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def writeMatlab(self, isubcase, f=None, isMagPhase=False):
        name = 'velocities'
        if self.nonlinear_factor is None:
            return self._writeMatlab(name, isubcase, f)
        else:
            return self._writeMatlabTransient(name, isubcase, f)

    def write_f06(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        if self.nonlinear_factor is not None:
            return self.writeF06Transient(header, pageStamp, pageNum, f)
        words = ['                                                   V E L O C I T Y   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        words += self.getTableMarker()
        return self._writeF06Block(words, header, pageStamp, pageNum, f)

    def writeF06Transient(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        words = ['                                                   V E L O C I T Y   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        words += self.getTableMarker()
        return self._writeF06TransientBlock(words, header, pageStamp, pageNum, f)

    def __repr__(self):
        if self.nonlinear_factor is not None:
            return self.__reprTransient__()

        msg = '---VELOCITIES---\n'
        msg += self.writeHeader()

        for nodeID, translation in sorted(self.translations.iteritems()):
            rotation = self.rotations[nodeID]
            gridType = self.gridTypes[nodeID]

            (dx, dy, dz) = translation
            (rx, ry, rz) = rotation

            msg += '%-10i %-8s ' % (nodeID, gridType)
            vals = [dx, dy, dz, rx, ry, rz]
            for val in vals:
                if abs(val) < 1e-6:
                    msg += '%10s ' % (0)
                else:
                    msg += '%10.3e ' % (val)
            msg += '\n'
        return msg

    def __reprTransient__(self):
        msg = '---TRANSIENT VELOCITY---\n'
        msg += self.writeHeader()

        for dt, translations in sorted(self.translations.iteritems()):
            msg += '%s = %g\n' % (self.data_code['name'], dt)
            for nodeID, translation in sorted(translations.iteritems()):
                rotation = self.rotations[dt][nodeID]
                gridType = self.gridTypes[nodeID]
                (dx, dy, dz) = translation
                (rx, ry, rz) = rotation

                msg += '%-10i %8s ' % (nodeID, gridType)
                vals = [dx, dy, dz, rx, ry, rz]
                for val in vals:
                    if abs(val) < 1e-6:
                        msg += '%10s ' % (0)
                    else:
                        msg += '%10.3e ' % (val)
                msg += '\n'
        return msg


class ComplexVelocityObject(ComplexTableObject):  # table_code=10, approach_code=???
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        ComplexTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def writeMatlab(self, isubcase, f=None, isMagPhase=False):
        name = 'velocities'
        if self.nonlinear_factor is None:
            return self._writeMatlab(name, isubcase, f, isMagPhase)
        else:
            return self._writeMatlabTransient(name, isubcase, f, isMagPhase)

    def write_f06(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        if self.nonlinear_factor is not None:
            return self.writeF06Transient(header, pageStamp, pageNum, f, isMagPhase)

        words = ['                                       C O M P L E X   V E L O C I T Y   V E C T O R\n']
        return self._writeF06Block(words, header, pageStamp, pageNum, f, isMagPhase)

    def writeF06Transient(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        words = ['                                       C O M P L E X   V E L O C I T Y   V E C T O R\n']
        return self._writeF06TransientBlock(words, header, pageStamp, pageNum, f, isMagPhase)

    def __repr__(self):
        return self.write_f06(['', '', ''], 'PAGE ', 1)[0]

        msg = '---COMPLEX VELOCITIES---\n'
        #if self.dt is not None:
        #    msg += '%s = %g\n' %(self.data_code['name'],self.dt)
        headers = ['DxReal', 'DxImag', 'DyReal', 'DyImag', 'DzReal', 'DyImag', 'RxReal', 'RxImag', 'RyReal', 'RyImag', 'RzReal', 'RzImag']
        msg += '%-10s ' % ('nodeID')
        for header in headers:
            msg += '%10s ' % (header)
        msg += '\n'

        for freq, translations in sorted(self.translations.iteritems()):
            msg += '%s = %g\n' % (self.data_code['name'], freq)

            for nodeID, translation in sorted(translations.iteritems()):
                rotation = self.rotations[freq][nodeID]

                msg += '%-10i ' % (nodeID)
                vals = translation + rotation
                for val in vals:
                    if abs(val) < 1e-6:
                        msg += '%10s ' % (0)
                    else:
                        msg += '%10.3e ' % (val)
                msg += '\n'
        return msg
