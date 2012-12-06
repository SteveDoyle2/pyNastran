import sys

# pyNastran
from pyNastran.op2.resultObjects.tableObject import (TableObject,
                                                     ComplexTableObject)


class DisplacementObject(TableObject):  # approach_code=1, thermal=0
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        TableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def writeMatlab(self, isubcase, f=None, isMagPhase=False):
        name = 'displacements'
        if self.nonlinear_factor is None:
            return self._writeMatlab(name, isubcase, f)
        else:
            return self._writeMatlabTransient(name, isubcase, f)

    def write_f06(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        if self.nonlinear_factor is not None:
            return self.writeF06Transient(header, pageStamp, pageNum, f)
        words = ['                                             D I S P L A C E M E N T   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        words += self.getTableMarker()
        return self._writeF06Block(words, header, pageStamp, pageNum, f)

    def writeF06Transient(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        words = ['                                             D I S P L A C E M E N T   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        words += self.getTableMarker()
        return self._writeF06TransientBlock(words, header, pageStamp, pageNum, f)

    def __repr__(self):
        if self.nonlinear_factor is not None:
            return self.__reprTransient__()

        msg = ['---DISPLACEMENTS---\n']
        msg.append(self.writeHeader())

        for nodeID, translation in sorted(self.translations.iteritems()):
            rotation = self.rotations[nodeID]
            gridType = self.gridTypes[nodeID]

            (dx, dy, dz) = translation
            (rx, ry, rz) = rotation

            msg2 = '%-10i %-8s ' % (nodeID, gridType)
            vals = [dx, dy, dz, rx, ry, rz]
            for val in vals:
                if abs(val) < 1e-6:
                    msg2 += '%10s ' % (0)
                else:
                    msg2 += '%10.3e ' % (val)
            msg2 = '\n'
            msg.append(msg2)
        return ''.join(msg)

    def __reprTransient__(self):
        msg = ['---TRANSIENT DISPLACEMENTS---\n']
        msg.append(self.writeHeader())

        for dt, translations in sorted(self.translations.iteritems()):
            msg2 = '%s = %g\n' % (self.data_code['name'], dt)
            for nodeID, translation in sorted(translations.iteritems()):
                rotation = self.rotations[dt][nodeID]
                gridType = self.gridTypes[nodeID]
                (dx, dy, dz) = translation
                (rx, ry, rz) = rotation

                msg2 += '%-10i %8s ' % (nodeID, gridType)
                vals = [dx, dy, dz, rx, ry, rz]
                for val in vals:
                    if abs(val) < 1e-6:
                        msg2 += '%10s ' % (0)
                    else:
                        msg2 += '%10.3e ' % (val)
                msg2 += '\n'
                msg.append(msg2)
        return ''.join(msg)


class ComplexDisplacementObject(ComplexTableObject):  # approach_code=1, sort_code=0, thermal=0
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        ComplexTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def writeMatlab(self, isubcase, f=None, isMagPhase=False):
        name = 'displacements'
        if self.nonlinear_factor is None:
            return self._writeMatlab(name, isubcase, f, isMagPhase)
        else:
            return self._writeMatlabTransient(name, isubcase, f, isMagPhase)

    def write_f06(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        if self.nonlinear_factor is not None:
            return self.writeF06Transient(header, pageStamp, pageNum, f, isMagPhase)

        words = ['                                       C O M P L E X   D I S P L A C E M E N T   V E C T O R\n']
        return self._writeF06Block(words, header, pageStamp, pageNum, f, isMagPhase)

    def writeF06Transient(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        words = ['                                       C O M P L E X   D I S P L A C E M E N T   V E C T O R\n']
        return self._writeF06TransientBlock(words, header, pageStamp, pageNum, f, isMagPhase)

    def __repr__(self):
        return self.write_f06(['', '', ''], 'PAGE ', 1)[0]

        msg = '---COMPLEX DISPLACEMENTS---\n'
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
